% -------------------------------------------------------------------
% Patrick Heng
% 7 Feb 2025
% 2D normal vectors for finite differences over a uniform computational 
% domain. Returns normal vector to lines of constant eta given mapping
% derivatives and Jacobian.
% -------------------------------------------------------------------

function [N_eta_xi, N_eta_eta] = normals_eta(x_xi,x_eta,y_xi,y_eta,J)
    % Normalization factor for normal vector
    mag_eta = 1./sqrt(x_xi.^2+y_xi.^2);
    
    n_eta(:,1) = -y_xi.*mag_eta;    % x component
    n_eta(:,2) = x_xi.*mag_eta;     % y component
    
    % Normal derivative coefficients, notation defined in report
    N_eta_xi = (y_eta.*n_eta(:,1) - x_eta.*n_eta(:,2))./J;
    N_eta_eta = (-y_xi.*n_eta(:,1) + x_xi.*n_eta(:,2))./J;

end
