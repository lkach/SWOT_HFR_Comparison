% 6-term 2D polynomial fit
% 
% IN:
% xx = 2D x-grid (NxM)
% yy = 2D y-grid
% dd = 2D data grid
% nn = number in either direction from a point to fit (e.g. nn = 2 means a
%      5x5 grid)
% 
% OUT:
% AA = coefficients corresponding to the terms, at each point (NxMx6):
%      [1, x, y, x^2, y^2, xy]
% dd_fit = 2D smoothed data grid
% 
% See:
% Yann-Treden Tranchant, Benoit Legresy, Annie Foppert, et al. SWOT reveals
% fine-scale balanced motions and dispersion properties in the Antarctic
% Circumpolar Current. ESS Open Archive . January 11, 2025.
% DOI: 10.22541/essoar.173655552.25945463/v1

function [AA,dd_fit] = fit_6term_2d(xx,yy,dd,nn)
xx.*yy.*dd; % size check
Collimate = @(IN) IN(:);
dd_fit = nan(size(dd));
AA = nan(size(dd,1), size(dd,2), 6);
II = size(xx,2);
JJ = size(xx,1);
for ii = [nn+1]:[II - nn]
    for jj = [nn+1]:[JJ - nn]
        HH = [Collimate(ones(size(xx([jj-nn]:[jj+nn],[ii-nn]:[ii+nn])))) ... 
              Collimate(xx([jj-nn]:[jj+nn],[ii-nn]:[ii+nn])) ...
              Collimate(yy([jj-nn]:[jj+nn],[ii-nn]:[ii+nn])) ...
              Collimate(xx([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]).^2) ...
              Collimate(yy([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]).^2) ...
              Collimate(xx([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]).*yy([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]))];
        HH_ = [1 ... 
               xx(jj,ii) ...
               yy(jj,ii) ...
               xx(jj,ii).^2 ...
               yy(jj,ii).^2 ...
               xx(jj,ii).*yy(jj,ii)];
        dd_ij = Collimate(dd([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]));
        if sum(isfinite(dd_ij)) == 0
            dd_fit(jj,ii) = nan;
            AA(jj,ii,:) = nan;
        else
            aa = HH(isfinite(dd_ij),:)\dd_ij(isfinite(dd_ij));
            dd_fit(jj,ii) = HH_*aa;
            AA(jj,ii,:) = aa;
        end
        % aa = HH\dd_ij;
        % dd_fit(jj,ii) = HH_*aa;
    end
end
dd_fit = dd_fit.*isfinite(dd)./isfinite(dd);
end
