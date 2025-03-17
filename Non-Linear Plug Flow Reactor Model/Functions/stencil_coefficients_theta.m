% -------------------------------------------------------------------
% Patrick Heng
% 3 Mar 2025
% 1D stencil coefficients for non-dimensionalized PFR equation
% -------------------------------------------------------------------

function [A,f] = stencil_coefficients_theta(X,U,P,theta,params)
    
    % Load parameters needed to create coefficients
    cA0 = params.cA0;
    stoich = params.stoich;
    THETA = params.THETA;
    p0 = params.p0;
    Pe_T = params.Pe_T;
    Da = params.Da;
    T0 = params.T0;
    Ta = params.Ta;
    K = params.K;
    beta = params.beta;
    Gamma = params.Gamma;
    N = params.N;
    a1 = params.a1;
    a2 = params.a2;
    a3 = params.a3;
    dt = params.dt;
    dz = params.dz;
    nodes = params.nodes;

    % Calculate enthalpy coefficient such that H*theta = total enthalpy
    AA = THETA'*a1 + THETA'*a2*T0*theta + THETA'*a3*T0^2*theta.^2  ...
            - (stoich'*a1+(stoich'*a2*T0)*theta+(stoich'*a3*T0^2) ...
                *theta.^2).*(X/stoich(1));
    H = (cA0/K)*AA./U;

    % Compute stencil coefficients
    W = -U.*H/(2*dz) - (1/Pe_T)*(1/dz^2);
    C = H/dt + (1/Pe_T)*(2/dz^2) + Da*exp(beta./theta).*(1-X)...
        .*(cA0*(stoich'*a1+(stoich'*a2*T0/2)...
        *theta+((1/3)*stoich'*a3*T0^2)*theta.^2)/K)./U + N;
    E = U.*H/(2*dz) - (1/Pe_T)*(1/dz^2);
    
    % Return coefficients in convenient form for using spdiags
    diag = [W,C,E];

    % Generate A at internal nodes
    A = spdiags(diag,[-1,0,1],nodes,nodes+1);
    A = A(:,1:nodes);

    % Clear boundary rows
    A(1,:) = 0; A(nodes,:) = 0;

    % Left flux boundary
    A(1,1) = U(1)*H(1) + 3/(2*Pe_T*dz); A(1,2) = -4/(2*Pe_T*dz); 
    A(1,3) = 1/(2*Pe_T*dz);

    % Right flux boundary
    A(nodes,nodes) = 3/(2*Pe_T*dz); A(nodes,nodes-1) = -4/(2*Pe_T*dz); 
    A(nodes,nodes-2) = 1/(2*Pe_T*dz);

    % Forcing function with appropriate boundary conditions
    % Internal
    f = H.*theta/dt - Da*exp(beta./theta).*(1-X).*Gamma./U + N*Ta/T0 ...
            + (p0/K/T0)*U.*([P(2:nodes);0]-[0;P(1:nodes-1)])/(2*dz);
    f(1) = 1;               % Left
    f(nodes) = 0;           % Right

end