% -------------------------------------------------------------------
% Patrick Heng
% 7 Feb 2025
% 2D transformation coefficients for Laplacian and derivatives for a 
% bilinear map. Used for finite differences.
% -------------------------------------------------------------------

function [x_xi,x_eta,y_xi,y_eta,x_xi_eta,...
            y_xi_eta,J,a,b,c,d,e,alpha,beta] = ...
            bilinear_map_derivatives(Xi,Eta,coeffs_x,coeffs_y)

    % Use coefficients as defined in class for better readability
    Rx = coeffs_x(1); Sx = coeffs_x(2); Tx = coeffs_x(3);
    Ry = coeffs_y(1); Sy = coeffs_y(2); Ty = coeffs_y(3);
    
    % Compute first derivatives of x,y over the entire mesh
    x_xi = Rx*Eta + Sx;
    x_eta = Rx*Xi + Tx;
    
    y_xi = Ry*Eta + Sy;
    y_eta = Ry*Xi + Ty;
    
    % Compute mixed derivatives of x,y over the entire mesh
    x_xi_eta = Rx;
    y_xi_eta = Ry;

    % Compute Jacobian determinant
    J = x_xi.*y_eta - x_eta.*y_xi;
    
    % Compute mapping coefficients for the Laplacian over the 
    % entire domain (notation as defined in class)
    a = x_eta.^2 + y_eta.^2;
    b = x_eta.*x_xi + y_eta.*y_xi;
    c = x_xi.^2 + y_xi.^2;
    
    alpha = -2*b.*x_xi_eta;
    beta = -2*b.*y_xi_eta;
    
    e = (x_eta.*beta - y_eta.*alpha)./J;
    d = (y_xi.*alpha - x_xi.*beta)./J;

end
