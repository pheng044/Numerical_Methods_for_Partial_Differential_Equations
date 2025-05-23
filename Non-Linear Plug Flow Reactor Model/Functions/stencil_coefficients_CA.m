% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025
% 1D stencil coefficients for non-dimensionalized PFR problem
%  -> concentration equation
% -------------------------------------------------------------------

function [A,f] = stencil_coefficients_CA(CA,U,theta,params)

    % Load parameters needed to create coefficients
    Pe_M = params.Pe_M;
    Da = params.Da;
    beta = params.beta;
    dt = params.dt;
    dz = params.dz;
    nodes = params.nodes;

    % Compute stencil coefficients
    W = - U/(2*dz) - (1/Pe_M)/dz^2;
    C =  1/dt + (1/Pe_M)*2/dz^2 + Da*exp(beta./theta);
    E = U/(2*dz) - (1/Pe_M)/dz^2;
    
    % Return coefficients in convenient form for using spdiags
    diag = [W,C,E];
    
    % Generate A at internal nodes
    A = spdiags(diag,[-1,0,1],nodes,nodes+1);
    A = A(:,1:nodes);

    % Clear boundary rows
    A(1,:) = 0; 
    A(nodes,:) = 0;
    
    % Left flux boundary
    A(1,1) = U(1) - 2/(Pe_M*dz);
    % Right flux boundary
    A(nodes,nodes) = 1/dz;
    A(nodes,nodes-1) = -1/dz;
  
    % Forcing function with appropriate boundary conditions
    f = CA/dt;                  % Internal
    f(1) = 1 - 2/(Pe_M*dz);     % Right
    f(nodes) = 0;               % Left
    
end