% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025
% Test script plot error convergence for PFR CUT scheme. MUST run
% 'PFR_FV_test_CUT_err.m' (previous script) first to generate the 
% data.
% -------------------------------------------------------------------
close all; clc;

% Mimimum/maximum exponent of the refined/coarsest mesh (h = 2^k_max)
k_max = 20;
k_min = 4;

% Load maximum refined mesh and store for comparison
temp = load(['sol_' num2str(2^k_max) '.mat'],'X','U','theta','dz');
X_star = temp.X;
U_star = temp.U;
theta_star = temp.theta;
h_star = temp.dz;

% Preallocate error vectors
e_L2_X = zeros(k_max-k_min-1,1);
e_inf_X = zeros(k_max-k_min-1,1);

e_L2_U = zeros(k_max-k_min-1,1);
e_inf_U = zeros(k_max-k_min-1,1);

e_L2_theta = zeros(k_max-k_min-1,1);
e_inf_theta = zeros(k_max-k_min-1,1);

% Loop through meshes, starting with the most refined and moving to
% the coarsest mesh
i = 1;
N = 2.^(k_max:-1:k_min);
h = zeros(k_max-k_min-1,1);

for n = N
    % Load the coarser mesh
    temp = load(['sol_' num2str(n) '.mat']);
    X = temp.X;
    U = temp.U;
    theta = temp.theta;
    dz = temp.dz;
    h(i) = dz;
    
    % Calculate the error by subtracting the refined mesh evaluated at
    % the coarse mesh grid points. Since the refinement halves each time
    % this means evaluating the u_star at every 'step'. 
    step = 2^(i-1);
    e_X = X_star(1:step:end)-X;
    e_U = U_star(1:step:end)-U;
    e_theta = theta_star(1:step:end)-theta;

    % Calculate normalized L2 error over the domain
    e_L2_X(i) = norm(e_X)*sqrt(dz);
    % Calculate the absolute maximum error
    e_inf_X(i) = max(abs(e_X),[],'all');

    % Calculate normalized L2 error over the domain
    e_L2_U(i) = norm(e_U)*sqrt(dz);
    % Calculate the absolute maximum error
    e_inf_U(i) = max(abs(e_U),[],'all');

    % Calculate normalized L2 error over the domain
    e_L2_theta(i) = norm(e_theta)*sqrt(dz);
    % Calculate the absolute maximum error
    e_inf_theta(i) = max(abs(e_theta),[],'all');

    % Increment loop counter
    i = i + 1;

end

% Plot error on a log-log plot, compare to O(h) and O(h^2) growth
figure

% L2 error
loglog(h, e_L2_X,marker='o',color='b',linewidth=2)
hold on
loglog(h, e_L2_U,marker='o',color=[1,180,40]/256,linewidth=2)
loglog(h, e_L2_theta,marker='o',color='r',linewidth=2)

% O(h)
loglog(h,2.^-(k_max:-1:k_min),linestyle=':',color='k',linewidth=2)
% O(h^2)
f = h.^2;
loglog(h,f,linestyle='--',color='k',linewidth=2)

% Pretty plot parameters
grid on
box on
axis padded
yticks(10.^(-8:2:0))
ylim([10^-8,1])
xlabel('$\Delta z$',interpreter='latex')
ylabel('$||e||_2$',interpreter='latex')
legend('$X$', '$U$', '$\theta$', '$O(\Delta z)$', ...
        '$O(\Delta z^2)$', interpreter='latex', location='northwest')
fontname('Serif'); fontsize(16,'points');

figure
% Maximum error
loglog(h, e_inf_X,marker='o',color='b',linewidth=2)
hold on
loglog(h, e_inf_U,marker='o',color=[1,180,40]/256,linewidth=2)
loglog(h, e_inf_theta,marker='o',color='r',linewidth=2)

% O(h)
loglog(h,h,linestyle=':',color='k',linewidth=2)
% O(h^2)

loglog(h,f,linestyle='--',color='k',linewidth=2)

% Pretty plot parameters
grid on
box on
axis padded
ylim([10^-8,1])
xticks(10.^(-8:2:0))
yticks(10.^(-8:2:0))
xlabel('$\Delta z$',interpreter='latex')
ylabel('$||e||_{\infty}$',interpreter='latex')

legend('$X$', '$U$', '$\theta$', '$O(\Delta z)$', ...
        '$O(\Delta z^2)$', interpreter='latex', location='northwest')
fontname('Serif'); fontsize(16,'points');