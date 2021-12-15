function [X, Y, MAP] = Regrid_map(x,y,map,nx,ny,xPol)
%SCALE_EFF_MAP scales 2D characteristic maps and the axis vector to a 
%desired grid size. Designed to scale electric motor efficiency maps, but
%can be used for any 2D map.

%Author: Sebastian Krapf
%05/12/2017

%% Explanation
% INPUT:
%     x: vector of axis of abscissa; SIZE: AAx1
%     y: vector of axis of ordinates; SIZE: 1xBB
%     map: matrix that contains data; SIZE: AAxBB
%     nx: number of desired data points abscissa
%     ny: number of desired data points ordinate
% OUTPUT:
%     X: scaled vector of axis of abscissa
%     Y: scaled vector of axis of ordinates
%     MAP: scaled matrix that contains data

%% Surpress warning
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId'); %Turn off interpolation warning

%% Data treatment
% create unisized vectors from x,y and map
[xData,yData] = ndgrid(x, y);

% create meshgrid with desired resolution (needed f. (extra-)interpolation)
X = min(x) : (max(x)-min(x))/(nx-1) : max(x); %grid size depends on the input vector. if x(1) ~= 0, grid size from n x is 
Y = min(y) : (max(y)-min(y))/(ny-1) : max(y);
[xq, yq] = meshgrid(X, Y);

% map (extra-)interpolation
F = griddedInterpolant(xData,yData,map','linear',xPol); %change of extrapolation method can lead to different results!
MAP = F( xq, yq);
X = X';

%% Plot - comment out, if not wished
% figure; mesh(x,y,map); 
% title('Raw data'); colorbar; 
% figure; mesh(X,Y,MAP); 
% title('Interpolated data'); colorbar; 

%% Turn on warning again
warning('on','MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
end

