% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025
% 1D stencil coefficients for non-dimensionalized PFR problem
%  -> momentum equation
% -------------------------------------------------------------------

function [A,f] = stencil_coefficients_U(U,rho,P,params)
    
    % Load parameters needed to create coefficients
    Re = params.Re;
    Eu = params.Eu;
    dt = params.dt;
    dz = params.dz;
    nodes = params.nodes;

    % Compute stencil coefficients
    W = -rho.*U/(2*dz) - (1/Re)/dz^2;
    C = rho/dt + (1/Re)*2/dz^2;
    E = rho.*U/(2*dz) - (1/Re)/dz^2;
    
    % Return coefficients in convenient form for using spdiags
    diag = [W,C,E];

    % Generate A at internal nodes
    A = spdiags(diag,[-1,0,1],nodes,nodes+1);
    A = A(:,1:nodes);

    % Clear boundary rows
    A(1,:) = 0; A(nodes,:) = 0;

    % Left flux boundary
    A(1,1) = rho(1)*U(1) - 2/(Re*dz);
    % Right flux boundary
    A(nodes,nodes) = 1/dz;
    A(nodes,nodes-1) = -1/dz;
    
    % Forcing function with appropriate boundary conditions
    f = U/dt + Eu*([P(2:nodes);0]-[0;P(1:nodes-1)])/(2*dz); %Internal
    f(1) = 1 - 2/(Re*dz);           % Left
    f(nodes) = 0;                   % Right

end