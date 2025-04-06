% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025
% Calculate pressure field for non-dimensionalized PFR equation
% -------------------------------------------------------------------

function P = pressure(X,U,theta,params)

    epsilon = params.epsilon;
    P = (1 + epsilon*X).*theta./U;

end