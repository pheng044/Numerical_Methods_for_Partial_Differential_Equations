% -------------------------------------------------------------------
% Patrick Heng
% 23 Feb. 2025
% Test script plot error convergence from our homework. MUST run
% 'HW_1_Err_Conv.m' (previous script) first to generate the data.
% -------------------------------------------------------------------
close all; clc;

% Maximum exponent of the refined mesh (h = 2^k_max)
k_max = 10;

% Load maximum refined mesh and store for comparison
temp = load('u_mesh_1024_1024.mat','u');
u_star = temp.u;

% Preallocate error vectors
e_L2 = zeros(k_max-1,1);
e_inf = zeros(k_max-1,1);

% Loop through meshes, starting with the most refined and moving to
% the coarsest mesh
i = 1;
N = 2.^(k_max:-1:2);

for n = N
    % Load the coarser mesh
    temp = load(['u_mesh_' num2str(n) '_' num2str(n) '.mat']);
    u = temp.u;
    
    % Calculate the error by subtracting the refined mesh evaluated at
    % the coarse mesh grid points. Since the refinement halves each time
    % this means evaluating the u_star at every 'step'. 
    % (I think something might be wrong here)
    step = 2^(i-1);
    e = u_star(1:step:end,1:step:end)-u;
    e = e(1:end,:);
    
    % Calculate normalized L2 error over the domain
    e_L2(i) = norm(e,'fro')*(1/(n-1));
    % Calculate the absolute maximum error
    e_inf(i) = max(e,[],'all');
    
    % Increment loop counter
    i = i + 1;

end

% Plot error on a log-log plot, compare to O(h) and O(h^2) growth
figure
% L2 error
loglog((N-1).^-1, e_L2,marker='o',color='b',linewidth=2)

hold on
% Maximum error
loglog((N-1).^-1, e_inf,marker='o',color='r',linewidth=2)

% O(h)
loglog((N-1).^-1,(N-1).^-1,linestyle=':',color='k',linewidth=2)
% O(h^2)
loglog((N-1).^-1,(N-1).^-2,linestyle='--',color='k',linewidth=2)

% Pretty plot parameters
grid on
box on
axis padded
yticks(10.^(-6:1:0))
xlabel('$h$',interpreter='latex')
ylabel('$||e||$',interpreter='latex')

legend('$||e||_{2}$', '$||e||_{\infty}$', '$O(h)$', '$O(h^2)$', interpreter='latex', location='best')
fontname('Serif'); fontsize(16,'points');
