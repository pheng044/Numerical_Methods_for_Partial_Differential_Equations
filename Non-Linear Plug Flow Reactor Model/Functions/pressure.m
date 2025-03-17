% -------------------------------------------------------------------
% Patrick Heng
% 3 Mar 2025
% Calculate pressure field for non-dimensionalized PFR equation
% -------------------------------------------------------------------

function P = pressure(X,U,theta,params)

    epsilon = params.epsilon;
    P = (1 + epsilon*X).*theta./U;

end