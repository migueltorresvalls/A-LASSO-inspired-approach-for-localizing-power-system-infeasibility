%% download and add matpower
addpath "C:\Program Files\MATLAB\R2024a\bin\matpower7.1"
%% load case data
case_name = 'case14.m';
M_in = loadcase(case_name); %load a feasible case
M_in = runpf(M_in); % run power flow (update the casedata with power flow solution)
[Nbus ,~] = size(M_in.bus);

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ... 
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
MU_PMAX , MU_PMIN , MU_QMAX , MU_QMIN , PC1 , PC2 , QC1MIN , QC1MAX , ... 
QC2MIN , QC2MAX , RAMP_AGC , RAMP_10 , RAMP_30 , RAMP_Q , APF] = idx_gen;

%% manipulate case to create infeasible system 
M_in.bus(:,PD) = 4.2*M_in.bus(:,PD); %increase Pload
M_in.bus(:,QD) = 4.2*M_in.bus(:,QD); %increase Qload 

% Get the original voltage magnitudes
original_Vm = M_in.bus(:, VM); % VM is the column index for voltage magnitudes from idx_bus

%% Prepare for constraints Ibus = Ybus dot V + n
[Ybus, ~, ~] = makeYbus(M_in); %Ybus is a complex sparse matrix
Ybus_real = real(Ybus);
Ybus_imag = imag(Ybus);

% Get Pbus, Qbus (which would be used to write the nonlinear I(v)function
Sbus = makeSbus(M_in.baseMVA, M_in.bus, M_in.gen);
Pbus = real(Sbus); 
Qbus = imag(Sbus);

%% use optimization tool (ipopt)
c = 0;
num_iterations = 100;
best_sparsity = -inf; % Start with the lowest possible sparsity best_nr = [];
best_ni = [];
best_Vr = []; 
best_Vi = [];

threshold = 0.8; % Elements with absolute value less than this are considered sparse

for i = 1:num_iterations
  % Solve the optimization problem
  [Vr_opt , Vi_opt , nr_opt , ni_opt] = optimizeInfeasibilityCurrents(Ybus_real, Ybus_imag, Pbus, Qbus, c, Nbus, original_Vm); 
  nsol = nr_opt+1j*ni_opt;

  % Calculate sparsity based on the defined threshold
  current_sparsity = sum(abs(nsol) > threshold);

  % Check if the current result has higher sparsity
  if current_sparsity > best_sparsity 
      best_sparsity = current_sparsity;
      best_Vr = Vr_opt; best_Vi = Vi_opt;
      best_nr = nr_opt; best_ni = ni_opt;
  end 
end


fprintf('Highest sparsity achieved (number of elements < %.2f): %d\n', threshold , best_sparsity);
fprintf('Real Voltage with Highest Sparsity:\n'); disp(best_Vr);
fprintf('Imaginary Voltage with Highest Sparsity:\n'); disp(best_Vi);
fprintf('Total Voltage with Highest Sparsity:\n');
disp(abs(best_Vr + 1j*best_Vi));
fprintf('Real Infeasibility Currents with Highest Sparsity:\n');
disp(best_nr);
fprintf('Imaginary Infeasibility Currents with Highest Sparsity:\n'); disp(best_ni);
fprintf('Total Infeasibility Currents with Highest Sparsity:\n'); disp(abs(best_nr + 1j*best_ni));

%Since the system is set up based on ECF and includes inputting current %instead of inputting power, putting voltage back to the system will of %course not converge.
%We check convergences by KCL:
V_mag2 = best_Vr.^2 + best_Vi.^2; % Square of the magnitude of the voltage at each bus
I_re_calc = (Pbus .* best_Vr + Qbus .* best_Vi) ./ V_mag2; % Real part of injected current
I_im_calc = -(Qbus .* best_Vr - Pbus .* best_Vi) ./ V_mag2; % Imaginary part of injected current


% Constructing the net current injections into the network, andincluding slack variable adjustments
ceq_current = [I_re_calc; I_im_calc] - ...
    [Ybus_real, -Ybus_imag, eye(Nbus), zeros(Nbus, Nbus);
    Ybus_imag, Ybus_real, zeros(Nbus, Nbus), eye(Nbus)] * [best_Vr; best_Vi; best_nr; best_ni];

% Checking if the Kirchhoff's Current Law (KCL) is satisfied to within a numerical tolerance
tolerance = 1e-6; % Define a small tolerance
if all(abs(ceq_current) < tolerance)
    fprintf('KCL is satisfied within the specified tolerance.\n');
else
    fprintf('KCL is not satisfied. Maximum deviation: %g\n', max(abs( ceq_current)));
end


function [Vr_opt , Vi_opt , nr_opt , ni_opt] = optimizeInfeasibilityCurrents(Ybus_real, Ybus_imag, Pbus, Qbus, c, numBuses , original_Vm)
    % Initial guess
    x0 = rand(4*numBuses , 1);

    % Set new VMAX and VMIN based on 1.2 and 0.8 times the original voltages
    VMAX = 1.2 * original_Vm; VMIN = 0.8 * original_Vm;

    % Define objective function
    obj = @(x) objectiveFunction(x, c, numBuses);

    % Set options for Ipopt via the OPTI Toolbox
    opts = optiset('solver', 'ipopt', 'display', 'iter'); 
    
    % Create OPTI Object with nonlinear constraints
    nlcon = @(x) nonlinearConstraints(x, Ybus_real, Ybus_imag, Pbus, Qbus, numBuses, VMAX, VMIN);
    nlrhs = zeros(4 * numBuses, 1); % All constraints are zero-bound ( inequalities are rewritten to <=0 form)
    nle = [zeros(2 * numBuses, 1); % Equality constraints for current injections
        ones(numBuses , 1); % Inequality constraints for VMAX( <=0)
        ones(numBuses , 1)]; % Inequality constraints for VMIN( <=0)
    
    Opt = opti('fun', obj, 'nl', nlcon, nlrhs, nle, 'x0', x0, 'options', opts);

    % Solve the optimization problem using Ipopt
    [x_opt , ~] = solve(Opt);

    % Extract optimized values
    Vr_opt = x_opt(1:numBuses);
    Vi_opt = x_opt(numBuses+1:2*numBuses);
    nr_opt = x_opt(2*numBuses+1:3*numBuses); ni_opt = x_opt(3*numBuses+1:end);
end

function f = objectiveFunction(x, c, numBuses) nr = x(2*numBuses+1:3*numBuses);
    ni = x(3*numBuses+1:end);
    n_mag = sqrt(nr.^2 + ni.^2);
    f = 0.5 * sum(n_mag.^2) + c * sum(n_mag);
end

function ceq = nonlinearConstraints(x, Ybus_real, Ybus_imag, Pbus, Qbus, numBuses, VMAX, VMIN)
    Vr = x(1:numBuses);
    Vi = x(numBuses+1:2*numBuses);
    V_mag2 = Vr.^2 + Vi.^2;

    % Compute the real and imaginary parts of current injections
    I_re_calc = (Pbus .* Vr + Qbus .* Vi) ./ V_mag2;
    I_im_calc = -(Qbus .* Vr - Pbus .* Vi) ./ V_mag2;

    % Current injection constraints (equality constraints)
    ceq_current = [I_re_calc; I_im_calc] - ...
    [Ybus_real , -Ybus_imag , eye(numBuses), zeros(numBuses , numBuses);Ybus_imag , Ybus_real , zeros(numBuses , numBuses), eye(numBuses)] * x;

    % Voltage magnitude constraints (inequality constraints) 
    ceq_voltage = [V_mag2 - VMAX.^2; % upper bound
    VMIN.^2 - V_mag2]; % lower bound (fixed sign)
    
    ceq = [ceq_current; ceq_voltage];
end
