% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025
% Sript to generate data for optimal operating conditions of acetone
% cracking in a PFR
% -------------------------------------------------------------------
% Script to solve the PFR boundary value problem using a non-linear 
% finite volume scheme.
% Uses the CUT loop to solve the coupled transport problems:
%   C - Concentration equation
%   U - U momentum equation
%   T - Temperature equation
% which is iteratively solved until convergence
% -------------------------------------------------------------------

close all; clear variables; clc;

% cA0 and Ta values to test
c_vec = linspace(5,100,50);
T_vec = linspace(1050,1300,51);

% Preallocate conversion, ketene molar flow, and power consumption 
% arrays
exit_X = zeros(length(c_vec),length(T_vec));
F = zeros(length(c_vec),length(T_vec));
W = zeros(length(c_vec),length(T_vec));

% ----- INPUTS -----

% Initial temperature, T0 (K)
T0 = 1035;

% Reactor length, L (m)
L = 5;

% PFR velocity, U (m/s)
U0 = 2;

% Stoichiometric coefficients [a,b,c,d]
stoich = [-1;1;1;0];

% Molecular weights (kg/mol), 
% components arranged in same order as stoich vector
MW = [58.08;42.04;16.04;28.01]/1000;

% Heat capacity coefficients: Cp = a1 + a2*T + a3T^2 (J/mol*K), 
% components arranged in same order as stoich vector

a1 = [26.63;20.04;13.39;6.25];
a2 = [0.183;0.0945;0.077;8.78*1e-3];
a3 = [-45.86;-30.95;-18.71;-0.021]*1e-6;

% Arrhenius frequency factor, k0, (1/s)
k0 = exp(34.34);

% Activation energy, Ea, (J/mol)
Ea = 284500;

% Gas constant, R, (J/mol*K)
R = 8.314;

% Reference temperature, Tref (K)
Tref = 298;

% Enthalpy of reaction at Tref, Href (J/mol*K)
H_f_298 = [-216670;-61090;-71840;0];
Href = stoich'*H_f_298;

% Mixture average thermal conductivity, kappa (J/m*K)
kappa = 2e-1;

% Mixture average diffusivity, D (m^2/s)
D = 1e-6;

% Effective heat transfer coefficient, h = UA/V where U is the overall
% heat transfer coefficient (W/m^2*K), A is reactor surface area (m^2), 
% and V is the reactor volume (m^3)
h = 16500;

% Mixture average viscocity, mu (kg/m*s)
mu = 30*1e-6;

% Steam heat capacity Cp = A + BT + DT^-2
AA = 3.470;
BB = 1.450e-3;
DD = 0.121e5;
Cp = @(T) R*(AA + BB*T + DD*T.^-2);

% Boiling point of water
Tb = 373.15;


% -------------------------------------------------------------------
% ----- SOLVER -----

% --- SOLVER PARAMETERS ---
% Number of nodes to solve for
nodes = 250;

% Error tolerance for relative residuals, generally, anything below 
% 1e-6 is not recommended as it leads to oscillations in the iterations
% without further convergence
rel_tol = 1e-8;

% Max number of iterations for each time step
max_iter = 25;

% Max number of time steps
max_time_iter = 500;    

% Below this tolerance, the solver will switch to a Newton algorithm,
% set to 0 if only fixed point iterations are desired. The Newton
% algorithm is efficient for problems below 10000 unknowns, but quickly
% becomes too slow for any bigger problems.
% Generally, do not set above 1e-3 since the Newton solver will converge 
% to unphysical solutions without a good initial guess
newton_tol = 0;

% If true, averaging will be performed on the outputs to make the 
% functions look visually smoother; introduces a small amount of error.
% However, by the mean value theorem, this error vanishes as dz -> 0.
smoothed_output = false;     

% Relaxation factors, 0.6 is good for many cases with a 1e-6 rel_tol,
% but they can be any number between 0 and 1
RFCA = 0.3;           % CA relaxation factor
RFT = RFCA;           % theta relaxation factor
RFU = RFCA;           % U relaxation factor

% Time step sizes
dt = Inf;

% Generate a uniform computational mesh
z = linspace(0,1,nodes);

% Spatial finite difference
dz = z(2) - z(1);

tic
for ii = 1:length(c_vec)
for jj = 1:length(T_vec)

% Initial concentration of A, cA0 (mol/m^3)
cA0 = c_vec(ii);
% Initial concentrations of other species (mol/m^3)
% i = 1        2       3        4
%     Acetone  Ketene  Methane  Nitrogen (inert)
c0 = [cA0;0;0;0];
% Heat exchanger temperature, Ta (K)
Ta = T_vec(jj);
  
% -------------------------------------------------------------------
% ----- PRE-SOLVE SETUP -----

% --- CALCULATED QUANTITIES ---
rho_0 = sum(c0.*MW);            % Initial mixture density
THETA = c0/cA0;                 % Normalized initial concentrations
% Initial mixture heat capacity
K = cA0*(THETA'*a1+THETA'*a2*T0+THETA'*a3*T0^2);           
p0 = sum(c0)*R*T0;              % Initial ideal gas pressure of mixture

% Calulate dimensionless numbers
Pe_M = U0*L/D;                   % Mass Peclet number
Pe_T = K*U0*L/kappa;             % Thermal Peclet number
Da = k0*L/U0;                    % First order Damkohler number
Re = rho_0*U0*L/mu;              % Reynolds number
Eu = 1/(rho_0*U0^2);             % Euler number

% Calculate recurring constants
beta = -Ea/(R*T0);               % Dimensionless activation energy
TCp_ref = stoich'*a1*Tref+(1/2)*stoich'*a2*Tref^2+(1/3)*stoich'*a3*Tref^3;

Gamma = cA0*(Href - TCp_ref)/(K*T0);  % Dimensionless reference enthalpy
N = h*L/(K*U0);                       % Number of transfer units
epsilon = -sum(stoich)*cA0/(sum(c0)*stoich(1));  % Change in moles of rxn

% Store parameters in structure for easy usage in functions
params = struct('cA0',cA0,'T0',T0,'Ta',Ta,'stoich',stoich, ...
                'rho_0',rho_0,'THETA',THETA,'K',K,'p0',p0,'MW',MW, ...
                'Pe_M',Pe_M,'Pe_T',Pe_T,'Da',Da,'Re',Re,'Eu',Eu, ...
                'beta',beta,'Gamma',Gamma,'N',N,'epsilon',epsilon, ...
                'a1',a1,'a2',a2,'a3',a3,'dz',dz,'dt',dt,'nodes',nodes);

% --- SOLUTION INITIALIZATION ---
% Maximum number of time steps
t = 1:max_time_iter;

% Initial guess for CA, assume an exponential curve
CA = (1-exp(-z))';

% Initial guess for theta, assume constant
theta = ones(nodes,1);

% Initial guess for U, assume constant
U = ones(nodes,1);

% Count the number of linear and non-linear iterations used
Linear_Solves = 0;

% Preallocate residual histories
err_CA = zeros(max_time_iter,1);
err_U = zeros(max_time_iter,1);
err_T = zeros(max_time_iter,1);


% ----------------------------------------------------------------
% ----- ITERATIVE SOLVER -----
for i = t
        
    % --- CA SOLVER ---
    % Save previous iteration to current
    CA_old = CA;

    % Gather CA stencil coefficients using values from the 
    % previous iteration, perform a linear solve for the updated 
    % CA
    [A,f] = stencil_coefficients_CA(CA,U,theta,params);
    %CA = A\f;
    [l,u] = ilu(A);
    CA = gmres(A,f,[],rel_tol,[],l,u,CA);

    % Update CA using underrelaxation
    CA = (1-RFCA)*CA_old + RFCA*CA; 

    % --- U SOLVER ---
    % Predict X, P, rho
    X = 1 - U.*CA;
    % Make sure X is bounded between 0 and 1, turn on if having
    % convergence issues in X
    % X(X<0) = 0; X(X>1) = 1;  

    P = pressure(X,U,theta,params);
    rho = density(X,theta,P,params);
    % Make sure density is positive, turn on if having convergence
    % issues
    % rho(rho<0) = 1e-14;
    
    % Save previous iteration to current
    U_old = U;
    % Gather U stencil coefficients using values from the previous
    % iteration, perform a linear solve for the updated U
    [A,f] = stencil_coefficients_U(U,rho,P,params);
    %U = A\f;
    [l,u] = ilu(A);
    U = gmres(A,f,[],rel_tol,[],l,u,U);

    % Update U using underrelaxation
    U = (1-RFU)*U_old + RFU*U;
    
    % --- THETA SOLVER ---
    % Update X,P field using newly calculated U 
    X = 1 - U.*CA;
    %X(X<0) = 1e-14; X(X>1) = 1;
    P = pressure(X,U,theta,params);
    % Save previous iteration to current
    theta_old = theta;

    % Gather theta stencil coefficients using values from the previous
    % iteration, perform a linear solve for the updated theta
    [A,f] = stencil_coefficients_theta(X,U,P,theta,params);
    %theta = A\f;
    [l,u] = ilu(A);
    theta = gmres(A,f,[],rel_tol,[],l,u,theta);

    % Update theta using underrelaxation
    theta = (1-RFT)*theta_old + RFT*theta;

    % If theta drops below 0 anywhere, reset it to a small number. 
    % Required to prevent the divergence of the first few iterations 
    % with a poor initial guess.
    theta(theta<0) = 1e-14;

    % Increment the number of linear solves (3 per iteration)
    Linear_Solves = Linear_Solves + 3;

    % Calculate L2 relative residuals
    res_CA = norm(CA-CA_old,2)/norm(CA,2)
    res_U = norm(U-U_old,2)/norm(U,2)
    res_T = norm(theta-theta_old,2)/norm(theta,2)

    % Store residual history
    err_CA(i) = res_CA; err_U(i) = res_U; err_T(i) = res_T;

    % Convergence criterion: all residuals less than the relative 
    % tolerance
    if res_CA < rel_tol && res_U < rel_tol && res_T < rel_tol
        break
    end

end

    H_rxn = Href - TCp_ref + ....
        stoich'*a1*T0+(1/2)*stoich'*a2*T0^2+(1/3)*stoich'*a3*T0^3;
    X = 1 - U(end)*CA(end);
    FB = X*cA0*U0;
    W_tot = U0*R*T0*cA0.*log(cA0*R*T0/101325) + ...
            cA0.*(2*U0*H_rxn./((Ta-T0).*Cp(Ta)))*R.* ...
            (AA*(Ta-Tb)+0.5*BB*(Ta^2-Tb^2)-DD*(Ta^-1-Tb^-1));
    exit_X(ii,jj) = X;
    F(ii,jj) = FB;
    W(ii,jj) = W_tot;

end
end
toc

% Save cost function surface
save('PFR_optimization','exit_X','W','F','c_vec','T_vec')