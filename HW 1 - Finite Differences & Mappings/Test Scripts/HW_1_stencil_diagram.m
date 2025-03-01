% -------------------------------------------------------------------
% Patrick Heng
% 21 Feb. 2025
% Script to generate diagram for 9-point stencil coefficients.
% -------------------------------------------------------------------
close all; clear all; clc;

% Generate 3x3 grid
x = 1:3;
y = 1:3;
h = 0.1;
[X,Y] = meshgrid(x,y);

% Plot nodal points and gridlines
hold on
plot(X,Y,marker='.',markersize=30,color='k')
affine_grid(X,Y);

% Add labels to the nodal points
% Top
text(1+h,3-h,texlabel('$u_{i-1,\ j+1}: NW_{i,j}$'),interpreter='latex')
text(2+h,3-h,texlabel('$u_{i,\ j+1}: N_{i,j}$'),interpreter='latex')
text(3+h,3-h,texlabel('$u_{i+1,\ j+1}: NE_{i,j}$'),interpreter='latex')

% Center
text(1+h,2-h,texlabel('$u_{i-1,\ j}: W_{i,j}$'),interpreter='latex')
text(2+h,2-h,texlabel('$u_{i,j}: P_{i,j}$'),interpreter='latex')
text(3+h,2-h,texlabel('$u_{i+1,\ j}: E_{i,j}$'),interpreter='latex')

% Bottom
text(1+h,1-h,texlabel('$u_{i-1,\ j-1}: SW_{i,j}$'),interpreter='latex')
text(2+h,1-h,texlabel('$u_{i,\ j-1}: S_{i,j}$'),interpreter='latex')
text(3+h,1-h,texlabel('$u_{i+1,\ j-1}: SE_{i,j}$'),interpreter='latex')

% Pretty plot parameters
pos = get(gca, 'Position');
pos(1) = 0.05; pos(3) = 0.6;
set(gca, 'Position', pos)
fontsize(16,'points')

xlabel(''); ylabel(''); xticks([]); yticks([])

% Plot gridlines for transformed coordinates
function affine_grid(X,Y)

    % Horizontal gridlines
    for i = 1:numel(X(1,:))
        line([X(1,1),X(1,end)],[Y(i,1),Y(i,end)],color='k')
    end
    
    % Vertical gridlines
    for i = 1:numel(Y(1,:))
        line([X(1,i),X(end,i)],[Y(1,i),Y(end,i)],color='k')
    end

end
