% -------------------------------------------------------------------
% Patrick Heng
% 7 Feb 2025
% 2D normal vectors for finite differences over a uniform computational 
% domain. Returns normal vector to lines of constant xi given mapping
% derivatives and Jacobian.
% -------------------------------------------------------------------

function [N_xi_xi, N_xi_eta] = normals_xi(x_xi,x_eta,y_xi,y_eta,J)
    % Normalization factor for normal vector
    mag_xi = 1./sqrt(x_eta.^2+y_eta.^2);
    
    n_xi(:,1) = y_eta.*mag_xi;      % x component
    n_xi(:,2) = -x_eta.*mag_xi;     % y component

    % Normal derivative coefficients, notation defined in report
    N_xi_xi = (y_eta.*n_xi(:,1) - x_eta.*n_xi(:,2))./J;
    N_xi_eta = (-y_xi.*n_xi(:,1) + x_xi.*n_xi(:,2))./J;

end
