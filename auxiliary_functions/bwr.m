function n = bwr(m)
%   BWR   Linear blue-white-red color map
%   BWR(M) returns an M-by-3 matrix containing a BWR colormap.
%   BWR, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB uses the length of the
%   default colormap.
%
%   To reset the colormap of the current figure:
%
%             colormap(bwr)
%


if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

l_m = round(m/2);
r_m = m - l_m;

% % Linear color gradient. This choice narrows the range of values that are
% % plotted nearly white.
r = [linspace(0,255,r_m)';255*ones(l_m,1)];
g = [linspace(0,255,r_m)';linspace(255,0,r_m)'];
b = [255*ones(l_m,1);linspace(255,0,r_m)'];

% % SQRT of the linear gradient above. This choice broadens the range of
% % values that are plotted nearly white.
% r = [sqrt(linspace(0,255,r_m)')*sqrt(255);255*ones(l_m,1)];
% g = [sqrt(linspace(0,255,r_m)')*sqrt(255);sqrt(linspace(255,0,r_m)')*sqrt(255)];
% b = [255*ones(l_m,1);sqrt(linspace(255,0,r_m)')*sqrt(255)];

% % Quadratic color gradient. This choice strongly narrows the range of
% % values that are plotted nearly white.
% r = [(linspace(0,255,r_m)'.^2)/255;255*ones(l_m,1)];
% g = [(linspace(0,255,r_m)'.^2)/255;(linspace(255,0,r_m)'.^2)/255];
% b = [255*ones(l_m,1);(linspace(255,0,r_m)'.^2)/255];

% % Consider this as an input option:
%   BWR(M,polyn) returns an M-by-3 matrix containing a BWR colormap, with
%       the gradient following a power law given by exponent "polyn"
%       (default 1). For example, polyn = 2 gives a wide range of
%       blue/red and a narrow range for white, polyn = 0.5 gives the
%       opposite.

n = [r g b];
n = n/255;

