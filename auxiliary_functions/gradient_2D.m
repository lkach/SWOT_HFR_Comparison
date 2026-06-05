% GRAD = gradient_2D(F,X,Y)
%   
% Calculates gradient "GRAD" from a 2D scalar field F using a central
% finite difference method. This function uses matrices for X and Y,
% allowing non-uniformly spaced or angled grids (though they still must be
% rectangular matrices).
% 
% IN:   F = scalar field whose gradient is desired (size Ny-by-Nx)
% IN:   X = x field (same size as F)
% IN:   Y = y field (same size as F)
% 
% OUT:  GRAD = matrix, gradient of F, size(GRAD) = [Ny Nx 2]
%              GRAD(:,:,1) = dF/dx
%              GRAD(:,:,2) = dF/dy

function GRAD = gradient_2D(F,X,Y)

if [sum(size(X) == size(F)) == 2] && [sum(size(Y) == size(F)) == 2]
else
    error('Inputs "F", "X", and "Y" must be the same size and shape.')
end

Mat_nan = nan(size(F,1) + 2, size(F,2) + 2);
F_nan = Mat_nan;
F_nan(2:(end-1),2:(end-1)) = F;
F = F_nan;
clear F_nan
X_nan = Mat_nan;
X_nan(2:(end-1),2:(end-1)) = X;
X = X_nan;
clear X_nan
Y_nan = Mat_nan;
Y_nan(2:(end-1),2:(end-1)) = Y;
Y = Y_nan;
clear Mat_nan Y_nan

% Angle displacement from Cartesian (in radians), size Ny by Nx
% x-hat vector in the input reference frame, rotated to the true reference frame:
Angle_x = [angle([X(3:end,2:[end-1]) - X(1:[end-2],2:[end-1])] + 1i*[Y(3:end,2:[end-1]) - Y(1:[end-2],2:[end-1])]) - pi/2];
% y-hat vector in the input reference frame, rotated to the true reference frame:
Angle_y =  angle([X(2:[end-1],3:end) - X(2:[end-1],1:[end-2])] + 1i*[Y(2:[end-1],3:end) - Y(2:[end-1],1:[end-2])]);

Angle_x = mod(Angle_x+pi/2, pi) - pi/2; Angle_y = mod(Angle_y+pi/2, pi) - pi/2;
% [median(Angle_x(:),'omitnan')    median(Angle_y(:),'omitnan')]*180/pi

% First index = y, Second index = x
GRAD = nan([[size(F)-[2 2]] 2]); % GRAD(:,:,1), GRAD(:,:,2)
% In the reference frame along the grid of X,Y:
GRADx = [F(2:[end-1],3:end) - F(2:[end-1],1:[end-2])] ./ [X(2:[end-1],3:end) - X(2:[end-1],1:[end-2])];
GRADy = [F(3:end,2:[end-1]) - F(1:[end-2],2:[end-1])] ./ [Y(3:end,2:[end-1]) - Y(1:[end-2],2:[end-1])];
% In the north/south Cartesian reference frame:
GRAD(:,:,1) = GRADx.*cos(Angle_x ) - GRADy.*sin(Angle_y );
GRAD(:,:,2) = GRADx.*sin(Angle_x ) + GRADy.*cos(Angle_y );

end

% https://scicomp.stackexchange.com/questions/21915/discrete-definitions-of-curl-nabla-times-f
% https://en.wikipedia.org/wiki/Finite_difference_method
% https://en.wikipedia.org/wiki/Rotation_matrix

%% DEMONSTRATION

% close all
% 
% [X,Y] = meshgrid(-10:0.5:10,-10:0.5:10);
% % Rotate the grid:
% THETA1 = pi/8;
% THETA2 = -pi/16;
% X_ = X*cos(THETA1) - Y*sin(THETA1);
% Y_ = X*sin(THETA2) + Y*cos(THETA2);
% X = X_ ; Y = Y_;
% 
% % F = (X-5).*(Y+7);
% % F = (X).*(Y);
% F = sin(X./3).*(Y);
% 
% GRAD = gradient_2D(F,X,Y);
% % %
% figure
% 
% subplot(121)
% pcolor_centered(X,Y,F); shading flat
% caxis([-1 1]*abs(max(F(:)))); colormap(bwr)
% colorbar
% title('Field')
% 
% subplot(122)
% pcolor_centered(X,Y, sqrt( (GRAD(:,:,1)).^2 + (GRAD(:,:,2)).^2) ); hold on; shading flat
% quiver(X,Y, -GRAD(:,:,1), -GRAD(:,:,2),'color',[1 1 1]*0);
% CB = colorbar; CB.Label.String = 'Abs(GRAD)';
% title('-Gradient of Field')
% 
% set(gcf,'Position',[255 378 1143 420])
% 
% figure
% pcolor_centered(X,Y, angle(GRAD(:,:,1) + 1i*GRAD(:,:,2)) ); hold on; shading flat
% quiver(X,Y, GRAD(:,:,1), GRAD(:,:,2),'color',[1 1 1]*0); colormap(hsv)
% CB = colorbar; CB.Label.String = 'Angle(GRAD)';
% title('Angle of the Gradient of the Field')
% 
% set(gcf,'Position',[255 378 1143 420])