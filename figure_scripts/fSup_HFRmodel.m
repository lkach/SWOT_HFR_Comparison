% This script presents a fake eddy that is observed by idealized HFR
% antennas, which observe radial velocity. Then, those velocities are
% mapped to the u,v components, in order to see how good or bad the
% reconstruction is. The two major effects that the HFR array imparts, and
% the eddy shift due to meanflow, are all presented here.

% Establish coordinate system:
dx = 100; % m
xx = [(-100000):dx:0]';
yy = [( -50000):dx:(50000)]';
[XX,YY] = meshgrid(xx,yy);

% Surface height:
SIGMA = 15000; % 20000; % Width of the Gaussian SSH feature
HH0 = -0.1; % SSH at center of eddy
XX0 = -50000; YY0 = 0; % Center of the eddy

% Gaussian:
HH = HH0*exp(-[[XX - XX0].^2 + [YY - YY0].^2]/[2 * SIGMA^2]);
% HH = 2*HH; HH(abs(HH) > 0.1) = -0.1;

% % Yukey (flat top):
% HH = zeros(size(XX));
% TUKEY = tukeywin(length(xx),1);
% rr = xx - XX0;
% for ii = 1:size(XX,1)
%     for jj = 1:size(XX,2)
%         HH(ii,jj) = HH0*interp1(rr,TUKEY,sqrt([XX(ii,jj) - XX0].^2 + [YY(ii,jj) - YY0].^2));
%     end
% end
% HH(~isfinite(HH)) = 0;

% % Add an SSH tilt that induces a southward coastal current:
% HH = HH - abs(HH0)*tanh([XX - XX0]/SIGMA);


% Geostrophic velocity
[HH_x,HH_y] = gradient(HH,dx,dx);
% sideral day = 23.9345 hours = 86164.2 seconds
LAT = 40; % degrees
ff = 2*[2*pi/(86164.2)]*sind(LAT);
gg = 9.81;
UUg = [-gg/ff]*HH_y; VVg = [+gg/ff]*HH_x;

% % % Add ageostrophic current that does not appear in SSH:
UU = UUg; VV = VVg; % Only geostrophic
% UU = UUg + [0]; VV = VVg + [-0.2]; % Eddy + southward current
% UU = UUg + [0]; VV = VVg - max(abs(UUg(:) + 1i*VVg(:))); % Eddy + southward current
% UU = UUg + [-0.25]; VV = VVg + [0]; % Offshore current (e.g. Ekman)
% UU = UUg + 0.05*randn(size(UUg)); VV = VVg + 0.05*randn(size(VVg)); % Add a noise term to velocity

% %% Define locations of HFR antennas

% x = real, y = imaginary
ANTENNAS = [ ...
    ...0 + 1i*50000; ...
    ...0 + 1i*30000; ...
    0 + 1i*15000; ...
    ...0 + 1i*0; ...
    0 - 1i*15000; ...
    ...0 - 1i*30000; ...
    ...0 - 1i*50000; ...
    ];

% Antenna observation locations
% Template centered at an antenna at {0,0}
dR = 3000; % radial distance between observation points (meters)
RR = [dR:dR:120000]';
dTHETA = 360/200;
THETA = [0:dTHETA:[360 - dTHETA]];
% THETA = [[90 + dTHETA]:dTHETA:[270 - dTHETA]];
Radial_Template = RR*exp(1i*THETA*[pi/180]);

% Define each antenna's observations points
Radial_Locations = {};
for ii = 1:length(ANTENNAS)
    OBS = ANTENNAS(ii) + Radial_Template(:);
    % Trim points where there are no data:
    OBS = OBS([real(OBS) < max(XX(:))] & ...
              [real(OBS) > min(XX(:))] & ...
              [imag(OBS) < max(YY(:))] & ...
              [imag(OBS) > min(YY(:))] );
    Radial_Locations{ii} = OBS;
end

% %% Define radial velocities

% Define each antenna's observations points
Vel_Radial = {};
for ii = 1:length(ANTENNAS)
    OBS = ANTENNAS(ii) + Radial_Template(:);
    % Vq = interp2(X,Y,V,Xq,Yq)
    UU_radinterp = interp2(XX,YY,UU,real(OBS),imag(OBS));
    VV_radinterp = interp2(XX,YY,VV,real(OBS),imag(OBS));
    UU_radinterp = UU_radinterp( [real(OBS) < max(XX(:))] & ...
                                 [real(OBS) > min(XX(:))] & ...
                                 [imag(OBS) < max(YY(:))] & ...
                                 [imag(OBS) > min(YY(:))] );
    VV_radinterp = VV_radinterp( [real(OBS) < max(XX(:))] & ...
                                 [real(OBS) > min(XX(:))] & ...
                                 [imag(OBS) < max(YY(:))] & ...
                                 [imag(OBS) > min(YY(:))] );
    rHat = Radial_Locations{ii} - ANTENNAS(ii);
    rHat = rHat./abs(rHat);
    Vel_Radial{ii} = [UU_radinterp.*real(rHat) + VV_radinterp.*imag(rHat)].*rHat;
end

% %% Make the 6-km grid

% Establish coordinate system:
dx_grid = 6000; % m
xx_grid = -1*flip([0:dx_grid:(100000)]');
yy_grid = [(-48000):dx_grid:(48000)]';
[XX_grid,YY_grid] = meshgrid(xx_grid,yy_grid);

% %% Visualize before the fit

% close all
% 
% figure
% pcolor_centered(XX/1000,YY/1000,HH); shading flat; hold on
% quiver(XX(1:10:end,1:10:end)/1000,...
%        YY(1:10:end,1:10:end)/1000,...
%        UU(1:10:end,1:10:end),...
%        VV(1:10:end,1:10:end),1,"Color","w");
% plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'k*')
% plot(XX_grid(:)/1000, YY_grid(:)/1000, 'r.')
% 
% % Add radial information to verifythat it is correct:
% a_i = 2;
% plot(  real(Radial_Locations{a_i}(:)/1000), ...
%        imag(Radial_Locations{a_i}(:)/1000), 'c.')
% quiver(real(Radial_Locations{a_i}(:)/1000), ...
%        imag(Radial_Locations{a_i}(:)/1000), ...
%        real(Vel_Radial{      a_i}(:)), ...
%        imag(Vel_Radial{      a_i}(:)),1,"Color","c")
% % scaling looks wrong conpared to the black arrows created above, but
% % inspecting the values confirms that they are correct, MATLAB's vector
% % size scaling is just bad.
% 
% plot(xx([1 end])/1000,YY0*[1 1]'/1000,'w')
% plot(XX0*[1 1]'/1000,yy([1 end])/1000,'w')
% 
% xlabel('X (km)'); ylabel('Y (km)');
% axis equal; axis tight
% colormap("turbo")
% CB = colorbar; CB.Label.String = 'SSH (m)';
% clim([-1 1]*max(abs(HH(:))))
% 
% figure
% subplot(211)
% plot(xx/1000,HH(round(size(HH,1)/2),:),'k'); hold on
% plot(xx/1000,VV(round(size(HH,1)/2),:),'r')
% plot(xx/1000,0*xx,'k--')
% legend('\eta (m)','v (m/s)')
% xlabel('X (km)'); % ylabel('\eta (m)');
% subplot(212)
% plot(yy/1000,HH(:,round(size(HH,2)/2)),'k'); hold on
% plot(yy/1000,UU(:,round(size(HH,2)/2)),'b')
% plot(yy/1000,0*yy,'k--')
% legend('\eta (m)','u (m/s)')
% xlabel('Y (km)'); % ylabel('\eta (m)');

% %% Obtain {U_grid, V_grid} from radials using the script that Mark Otero shared with Luke Kachelein in 2021

tic
SEARCH_RADIUS = 10000; % search radius for fitting (m)
UU_grid_default = nan(size(XX_grid));
VV_grid_default = nan(size(XX_grid));
DOPX_grid_default = nan(size(XX_grid));
DOPY_grid_default = nan(size(XX_grid));
HDOP_grid_default = nan(size(XX_grid));
UU_grid_planar = nan(size(XX_grid));
VV_grid_planar = nan(size(XX_grid));
DOPX_grid_planar = nan(size(XX_grid));
DOPY_grid_planar = nan(size(XX_grid));
HDOP_grid_planar = nan(size(XX_grid));
UU_grid_quadratic = nan(size(XX_grid));
VV_grid_quadratic = nan(size(XX_grid));
DOPX_grid_quadratic = nan(size(XX_grid));
DOPY_grid_quadratic = nan(size(XX_grid));
HDOP_grid_quadratic = nan(size(XX_grid));
for ii = 1:size(XX_grid,1)
    for jj = 1:size(XX_grid,2)
        FitCenter = XX_grid(ii,jj) + 1i*YY_grid(ii,jj);
        XY_Radial_ij = [];
        VelRadial_ij = [];
        for kk = 1:length(Vel_Radial)
            IND = abs([Radial_Locations{kk} - FitCenter]) < SEARCH_RADIUS;
            if sum(IND) == 0
            else
                XY_Radial_ij = [XY_Radial_ij ; Radial_Locations{kk}(IND)];
                % ^ This is actually unused by the algorithm, as it is
                % unweighted by position (I might not do that if I designed
                % the system, but maybe that would just reduce the amount
                % of usable data?).
                VelRadial_ij = [VelRadial_ij ; Vel_Radial{kk}(IND)];
            end
        end
        % % Default method, mean in search area:
        [uGrid,vGrid,dopx,dopy,hdop,~] = uwlsTotal( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi );
        UU_grid_default(ii,jj) = uGrid;
        VV_grid_default(ii,jj) = vGrid;
        DOPX_grid_default(ii,jj) = dopx;
        DOPY_grid_default(ii,jj) = dopy;
        HDOP_grid_default(ii,jj) = hdop;

        % Alternative method, planar fit in search area:
        [uGrid,vGrid,dopx,dopy,hdop] = uwlsTotal_planar( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi , real(XY_Radial_ij), imag(XY_Radial_ij) , real(FitCenter) , imag(FitCenter) , SEARCH_RADIUS );
        UU_grid_planar(ii,jj) = uGrid;
        VV_grid_planar(ii,jj) = vGrid;
        DOPX_grid_planar(ii,jj) = dopx;
        DOPY_grid_planar(ii,jj) = dopy;
        HDOP_grid_planar(ii,jj) = hdop;

        % Alternative method, planar fit in search area:
        [uGrid,vGrid,dopx,dopy,hdop] = uwlsTotal_quadratic( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi , real(XY_Radial_ij), imag(XY_Radial_ij) , real(FitCenter) , imag(FitCenter) , SEARCH_RADIUS );
        UU_grid_quadratic(ii,jj) = uGrid;
        VV_grid_quadratic(ii,jj) = vGrid;
        DOPX_grid_quadratic(ii,jj) = dopx;
        DOPY_grid_quadratic(ii,jj) = dopy;
        HDOP_grid_quadratic(ii,jj) = hdop;
        

    end
end
toc



%% Visualize with the fit:

close all

figure % 1
subplot(2,3,[1 2 4 5])
pcolor_centered(XX/1000,YY/1000,HH); shading flat; hold on

% % % arrows and streamlines for TRUE CURRENT
% quiver(XX(1:10:end,1:10:end)/1000,...
%        YY(1:10:end,1:10:end)/1000,...
%        UU(1:10:end,1:10:end),...
%        VV(1:10:end,1:10:end),1,"Color","w");
% SL = streamline(XX(1:10:end,1:10:end)/1000,...
%                 YY(1:10:end,1:10:end)/1000,...
%                 UU(1:10:end,1:10:end),...
%                 VV(1:10:end,1:10:end),...
%                 xx(1:10:end)/1000, 0 + 0*xx(1:10:end)/1000);
%                 % [-37000,-36000,-35000]/1000, [0 0 0] + 0.01);
%      set(SL,'Color','w')

plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ko','MarkerFaceColor','k')

% quiver(XX_grid(:)/1000,YY_grid(:)/1000,...
%        UU_grid(:),VV_grid(:),...
%        1,"Color","r");
% plot(XX_grid(:)/1000, YY_grid(:)/1000, 'r.')

plot(xx([1 end])/1000,YY0*[1 1]'/1000,'w')
plot(XX0*[1 1]'/1000,yy([1 end])/1000,'w')

xlabel('X (km)'); ylabel('Y (km)');
axis equal; axis tight
colormap("turbo")
CB = colorbar; CB.Label.String = 'SSH (m)';
clim([-1 1]*max(abs(HH(:))))

% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid,XX,YY,'spine');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid,XX,YY,'spine');

% % % arrows and streamlines for INTERPOLATED GRIDDED HFR
quiver(XX(1:10:end,1:10:end)/1000,...
       YY(1:10:end,1:10:end)/1000,...
       UU_grid_interp(1:10:end,1:10:end),...
       VV_grid_interp(1:10:end,1:10:end),...
       1,"Color","r");
SL = streamline(XX(1:10:end,1:10:end)/1000,...
                YY(1:10:end,1:10:end)/1000,...
                UU_grid_interp(1:10:end,1:10:end),...
                VV_grid_interp(1:10:end,1:10:end),...
                ...-20 + 0*yy(1:10:end)/1000, yy(1:10:end)/1000);
                [-9:0]*60/10, zeros(1,10));
     set(SL,'Color','r')

subplot(2,3,3)
pcolor_centered(XX(1:10:end,1:10:end)/1000,...
                YY(1:10:end,1:10:end)/1000,...
                sqrt(UU(1:10:end,1:10:end).^2 + VV(1:10:end,1:10:end).^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'wo','MarkerFaceColor','w')
axis equal; axis tight
xlabel('X (km)'); ylabel('Y (km)');
clim([0 1]*max(sqrt([UU(:).^2 + VV(:).^2 ; UU_grid_interp(:).^2 + VV_grid_interp(:).^2])))
CB = colorbar; CB.Label.String = '|u| (m)';

subplot(2,3,6)
pcolor_centered(XX(1:10:end,1:10:end)/1000,...
                YY(1:10:end,1:10:end)/1000,...
                sqrt(UU_grid_interp(1:10:end,1:10:end).^2 + VV_grid_interp(1:10:end,1:10:end).^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'wo','MarkerFaceColor','w')
axis equal; axis tight
clim([0 1]*max(sqrt([UU(:).^2 + VV(:).^2 ; UU_grid_interp(:).^2 + VV_grid_interp(:).^2])))
CB = colorbar; CB.Label.String = '|u_{gridded}| (m)';



figure % 2
contour(XX/1000,YY/1000,abs(UU + 1i*VV),"Color","k"); hold on
contour(XX/1000,YY/1000,abs(UU_grid_interp + 1i*VV_grid_interp),"Color","r"); 
% set(gca,'Color','k')
legend('u + iv','u_{gridded} + iv_{gridded}')
axis equal; axis tight



figure % 3
% pcolor_centered(XX/1000,YY/1000,abs([UU_grid_interp + 1i*VV_grid_interp] - [UU + 1i*VV]) ./ ...
%                                 abs(UU + 1i*VV)); clim([0 1]) % ratio
pcolor_centered(XX/1000,YY/1000,abs([UU_grid_interp + 1i*VV_grid_interp] - [UU + 1i*VV]) ); % value
shading flat; hold on
quiver(XX(             1:20:end, ...
                       1:20:end)/1000,...
       YY(             1:20:end, ...
                       1:20:end)/1000,...
       UU_grid_interp( 1:20:end, ...
                       1:20:end) - ...
       UU(             1:20:end, ...
                       1:20:end),...
       VV_grid_interp( 1:20:end, ...
                       1:20:end) - ...
       VV(             1:20:end, ...
                       1:20:end),1,"Color","w");
plot(xx([1 end])/1000,YY0*[1 1]'/1000,'w')
plot(XX0*[1 1]'/1000,yy([1 end])/1000,'w')
colormap("turbo")
CB = colorbar; CB.Label.String = '$|\mathbf{u}_\mathrm{gridded} - \mathbf{u}|$';
CB.Label.Interpreter = 'LaTeX';
axis equal; axis tight


figure % 4
VelocityMagnitudeRatio = abs(UU_grid_interp + 1i*VV_grid_interp)./abs(UU + 1i*VV);

% MASK = VelocityMagnitudeRatio < 2; MASK_CONDITION = '$|\mathbf{u}_\mathrm{gridded}|/|\mathbf{u}|<$ 2';
MASK = abs(UU_grid_interp + 1i*VV_grid_interp) > 0.1; MASK_CONDITION = '$|\mathbf{u}_\mathrm{gridded}|>$0.1';

MASK = MASK./MASK;
pcolor_centered(XX/1000,YY/1000,... % ratio
                VelocityMagnitudeRatio.*[MASK]); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ko','MarkerFaceColor','k')
CB = colorbar; CB.Label.String = ...
    ['$|\mathbf{u}_\mathrm{gridded}|/|\mathbf{u}|$ where ' , MASK_CONDITION];
CB.Label.Interpreter = 'LaTeX';
CB.Label.FontSize = 20;
clim([0 2])
colormap("turbo")
axis equal; axis tight


figure % 5
subplot(211)
plot(xx/1000,HH(round(size(HH,1)/2),:),'k'); hold on
plot(xx/1000,VV(round(size(HH,1)/2),:),'r')
plot(xx_grid/1000,VV_grid(9,:),'.--m')
plot(xx/1000,0*xx/1000,'k--')
legend('\eta (m)','v (m/s)','v gridded')
xlabel('X (km)'); % ylabel('\eta (m)');
subplot(212)
plot(yy/1000,HH(:,round(size(HH,2)/2)),'k'); hold on
plot(yy/1000,UU(:,round(size(HH,2)/2)),'b')
plot(yy_grid/1000,UU_grid(:,9),'.--c')
plot(yy/1000,0*yy/1000,'k--')
legend('\eta (m)','u (m/s)','u gridded')
xlabel('Y (km)'); % ylabel('\eta (m)');

%% Supplemental Figure for Publication, compare fitting methods

close all

% UU = UUg; VV = VVg; % Only geostrophic
% UU = UUg; VV = VVg - 0.1; % Geostrophic + Southward

figure('Color','w') % 1
tiledlayout(2,2,"TileSpacing","tight")

nexttile(1)
pcolor_centered(XX/1000,...
                YY/1000,...
                sqrt(UU.^2 + VV.^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-50,0,'+','Color','w')
axis equal; axis tight
xticks([-100 -75 -50 -25 0]); yticks([-50 -25 0 25 50])
% clim([0 1]*max(sqrt([UU(:).^2 + VV(:).^2 ; UU_grid_interp(:).^2 + VV_grid_interp(:).^2])))
CLIM = [0 1]*max(sqrt([UU(:).^2 + VV(:).^2]));
% CLIM = [0 1]*0.45;
clim(CLIM)
CB = colorbar; CB.Label.String = '$|\vec{u}|$ (m s$^{-1}$)';
CB.Label.Interpreter = "latex";
CB.Label.FontSize = 18; CB.Location = "northoutside";
text(-100,50,'$\quad|\vec{u}_g|$','VerticalAlignment','top','Color','w','Interpreter','latex','FontSize',20)
set(gca,'FontSize',14)

nexttile(2)
% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid_default,XX,YY,'spine');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid_default,XX,YY,'spine');
pcolor_centered(XX/1000,...
                YY/1000,...
                sqrt(UU_grid_interp.^2 + VV_grid_interp.^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-50,0,'+','Color','w')
axis equal; axis tight
xticks([-100 -75 -50 -25 0]); yticks([-50 -25 0 25 50])
% clim([0 1]*max(sqrt([UU(:).^2 + VV(:).^2 ; UU_grid_interp(:).^2 + VV_grid_interp(:).^2])))
clim(CLIM)
% CB = colorbar; CB.Label.String = '|u_{gridded,default}| (m)';
% CB.Label.FontSize = 12;
text(-100,50,'$\quad|\vec{u}_\mathrm{default}|$','VerticalAlignment','top','Color','w','Interpreter','latex','FontSize',20)
set(gca,'FontSize',14)

nexttile(3)
% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid_planar,XX,YY,'spine');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid_planar,XX,YY,'spine');
pcolor_centered(XX/1000,...
                YY/1000,...
                sqrt(UU_grid_interp.^2 + VV_grid_interp.^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
xlabel('X (km)'); ylabel('Y (km)');
text(-50,0,'+','Color','w')
axis equal; axis tight
xticks([-100 -75 -50 -25 0]); yticks([-50 -25 0 25 50])
% clim([0 1]*max(sqrt([UU(:).^2 + VV(:).^2 ; UU_grid_interp(:).^2 + VV_grid_interp(:).^2])))
clim(CLIM)
% CB = colorbar; CB.Label.String = '|u_{gridded,planar}| (m)';
% CB.Label.FontSize = 12;
text(-100,50,'$\quad|\vec{u}_\mathrm{planar}|$','VerticalAlignment','top','Color','w','Interpreter','latex','FontSize',20)
set(gca,'FontSize',14)

nexttile(4)
% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid_quadratic,XX,YY,'spine');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid_quadratic,XX,YY,'spine');
pcolor_centered(XX/1000,...
                YY/1000,...
                sqrt(UU_grid_interp.^2 + VV_grid_interp.^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-50,0,'+','Color','w')
axis equal; axis tight
xticks([-100 -75 -50 -25 0]); yticks([-50 -25 0 25 50])
% clim([0 1]*max(sqrt([UU(:).^2 + VV(:).^2 ; UU_grid_interp(:).^2 + VV_grid_interp(:).^2])))
clim(CLIM)
% CB = colorbar; CB.Label.String = '|u_{gridded,quadratic}| (m)';
% CB.Label.FontSize = 12;
colormap('turbo')
text(-100,50,'$\quad|\vec{u}_\mathrm{quadratic}|$','VerticalAlignment','top','Color','w','Interpreter','latex','FontSize',20)
set(gca,'FontSize',14)
set(gcf,'Position',[708    67   733   730])

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

figure('Color','w') % 2
tiledlayout(2,2,"TileSpacing","tight")

% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid_default,XX,YY,'spine');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid_default,XX,YY,'spine');

nexttile(1)
pcolor_centered(XX/1000,...
                YY/1000,...
                0*sqrt(UU.^2 + VV.^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-50,0,'+','Color','w')
axis equal; axis tight
xticks([-100 -75 -50 -25 0]); yticks([-50 -25 0 25 50])
Collimate = @(IN) IN(:);
CLIM = [0 1]*max(Collimate(sqrt([UU_grid_interp - UU].^2 + ...
                                [VV_grid_interp - VV].^2)));
clim(CLIM)
CB = colorbar; CB.Label.String = '$|\Delta\vec{u}|$ (m s$^{-1}$)';
CB.Label.Interpreter = "latex";
CB.Label.FontSize = 18; CB.Location = "northoutside";
text(-100,50,'$\quad|\vec{u}_g - \vec{u}_g|$','VerticalAlignment','top','Color','w','Interpreter','latex','FontSize',20)
set(gca,'FontSize',14)

nexttile(2)
pcolor_centered(XX/1000,...
                YY/1000,...
                sqrt([UU_grid_interp - UU].^2 + ...
                     [VV_grid_interp - VV].^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-50,0,'+','Color','w')
axis equal; axis tight
xticks([-100 -75 -50 -25 0]); yticks([-50 -25 0 25 50])
% clim([0 1]*max(sqrt([UU(:).^2 + VV(:).^2 ; UU_grid_interp(:).^2 + VV_grid_interp(:).^2])))
clim(CLIM)
% CB = colorbar; CB.Label.String = '|u_{gridded,default}| (m)';
% CB.Label.FontSize = 12;
text(-100,50,'$\quad|\vec{u}_\mathrm{default} - \vec{u}_g|$','VerticalAlignment','top','Color','w','Interpreter','latex','FontSize',20)
set(gca,'FontSize',14)

nexttile(3)
% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid_planar,XX,YY,'spine');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid_planar,XX,YY,'spine');
pcolor_centered(XX/1000,...
                YY/1000,...
                sqrt([UU_grid_interp - UU].^2 + ...
                     [VV_grid_interp - VV].^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
xlabel('X (km)'); ylabel('Y (km)');
text(-50,0,'+','Color','w')
axis equal; axis tight
xticks([-100 -75 -50 -25 0]); yticks([-50 -25 0 25 50])
% clim([0 1]*max(sqrt([UU(:).^2 + VV(:).^2 ; UU_grid_interp(:).^2 + VV_grid_interp(:).^2])))
clim(CLIM)
% CB = colorbar; CB.Label.String = '|u_{gridded,planar}| (m)';
% CB.Label.FontSize = 12;
text(-100,50,'$\quad|\vec{u}_\mathrm{planar} - \vec{u}_g|$','VerticalAlignment','top','Color','w','Interpreter','latex','FontSize',20)
set(gca,'FontSize',14)

nexttile(4)
% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid_quadratic,XX,YY,'spine');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid_quadratic,XX,YY,'spine');
pcolor_centered(XX/1000,...
                YY/1000,...
                sqrt([UU_grid_interp - UU].^2 + ...
                     [VV_grid_interp - VV].^2)); hold on
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-50,0,'+','Color','w')
axis equal; axis tight
xticks([-100 -75 -50 -25 0]); yticks([-50 -25 0 25 50])
% clim([0 1]*max(sqrt([UU(:).^2 + VV(:).^2 ; UU_grid_interp(:).^2 + VV_grid_interp(:).^2])))
clim(CLIM)
% CB = colorbar; CB.Label.String = '|u_{gridded,quadratic}| (m)';
% CB.Label.FontSize = 12;
colormap('turbo')
text(-100,50,'$\quad|\vec{u}_\mathrm{quadratic} - \vec{u}_g|$','VerticalAlignment','top','Color','w','Interpreter','latex','FontSize',20)
set(gca,'FontSize',14)
set(gcf,'Position',[708    67   733   730])

%%
% figure(1)
% exportgraphics(gcf,...
% ['/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/FSup_HFRgriddingmethods.pdf'],...
% 'BackgroundColor','none','ContentType','image')

%% Supplemental Figure for Publication, examine shortcomings of gridding in general and effect of meanflow

DIFFERENT_SOUTHWARD_CURRENTS = false;

close all

figure('Color','w')
tiledlayout(3,2,'TileSpacing','tight')
FONTSIZE = 16;
BG_white = [1 1 1 0.5];
BG_black = [0 0 0 0.2];
% BG_white = [1 1 1 0];
% BG_black = [0 0 0 0];

AX1 = nexttile(1); % Geostrophic cyclonic eddy
pcolor_centered(XX/1000,YY/1000,HH); hold on
Stream_points = {[-50 -50 -50 -50 -50], ...
                 [  5  10  15  20  25]};
    UU = UUg; VV = VVg; % Only geostrophic
SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU(1:30:end,1:30:end),...
       VV(1:30:end,1:30:end),1,"filled",'Color','k')
text(-50,-43,'$\vec{u} = \vec{u}_\mathrm{g}$',...
     'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','BackgroundColor',BG_white)
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
clim([-1 1]*abs(HH0))
% xlabel('x (km)');
ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)
CB = colorbar; CB.Label.String = '$\eta$ (m)';
CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
CB.Ticks = [-1 0 1]*abs(HH0);
CB.Location = "northoutside";


AX2 = nexttile(2); % Geostrophic cyclonic eddy, velocity as the color quantity
UU = UUg; VV = VVg; % Only geostrophic
pcolor_centered(XX/1000,YY/1000,abs(UU + 1i*VV)); hold on
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU(1:30:end,1:30:end),...
       VV(1:30:end,1:30:end),1,"filled",'Color','w')
text(-50,-43,'$\vec{u} = \vec{u}_\mathrm{g}$',...
     'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','Color','w','BackgroundColor',BG_black)
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
clim([0 1]*max(abs(UU(:) + 1i*VV(:))))
% xlabel('x (km)');
% ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)
CB = colorbar; CB.Label.String = '$|\vec{u}|$ (m/s)';
CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
CB.Ticks = [0:0.1:0.4];
CB.Location = "northoutside";


AX3 = nexttile(3); % Geostrophic cyclonic eddy shifted by strong meanflow
Stream_points = {[-50 -50 -50 -50 -50], ...
                 [  5  10  15  20  25]};
pcolor_centered(XX/1000,YY/1000,HH); hold on
SL = streamline(XX/1000,YY/1000,UU_grid_interp,VV_grid_interp,  Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
SL = streamline(XX/1000,YY/1000,-UU_grid_interp,-VV_grid_interp,Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU_grid_interp(1:30:end,1:30:end),...
       VV_grid_interp(1:30:end,1:30:end),1,"filled",'Color','k')
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-50,-43,'$\vec{u} = $ gridded $\vec{u}_\mathrm{g}$',...
     'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','Color','k','BackgroundColor',BG_white)
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
clim([-1 1]*abs(HH0))
% xlabel('x (km)');
ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)


AX4 = nexttile(4); % Geostrophic cyclonic eddy, offshore diminishment from gridding
UU = UUg; VV = VVg; % Only geostrophic
% % % 
SEARCH_RADIUS = 10000; % search radius for fitting (m)
UU_grid = nan(size(XX_grid));
VV_grid = nan(size(XX_grid));
DOPX_grid = nan(size(XX_grid));
DOPY_grid = nan(size(XX_grid));
HDOP_grid = nan(size(XX_grid));
for ii = 1:size(XX_grid,1)
    for jj = 1:size(XX_grid,2)
        FitCenter = XX_grid(ii,jj) + 1i*YY_grid(ii,jj);
        XY_Radial_ij = [];
        VelRadial_ij = [];
        for kk = 1:length(Vel_Radial)
            IND = abs([Radial_Locations{kk} - FitCenter]) < SEARCH_RADIUS;
            if sum(IND) == 0
            else
                XY_Radial_ij = [XY_Radial_ij ; Radial_Locations{kk}(IND)];
                % ^ This is actually unused by the algorithm, as it is
                % unweighted by position (I might not do that if I designed
                % the system, but maybe that would just reduce the amount
                % of usable data?).
                VelRadial_ij = [VelRadial_ij ; Vel_Radial{kk}(IND)];
            end
        end
        % % Default method, mean in search area:
        [uGrid,vGrid,dopx,dopy,hdop] = uwlsTotal( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi );

        % % Alternative method, planat or quadratic in search area:
        % [uGrid,vGrid,dopx,dopy,hdop] = ...
        %     uwlsTotal_modified( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi , real(XY_Radial_ij), imag(XY_Radial_ij) , real(FitCenter) , imag(FitCenter) , SEARCH_RADIUS, 'PLANAR' );
        
        UU_grid(ii,jj) = uGrid;
        VV_grid(ii,jj) = vGrid;
        DOPX_grid(ii,jj) = dopx;
        DOPY_grid(ii,jj) = dopy;
        HDOP_grid(ii,jj) = hdop;
    end
end
% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid,XX,YY,'spine');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid,XX,YY,'spine');
% % % 
pcolor_centered(XX/1000,YY/1000,abs(UU_grid_interp + 1i*VV_grid_interp)); hold on
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU_grid_interp(1:30:end,1:30:end),...
       VV_grid_interp(1:30:end,1:30:end),1,"filled",'Color','w')
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-50,-43,'$\vec{u} = $ gridded $\vec{u}_\mathrm{g}$',...
     'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','Color','w','BackgroundColor',BG_black)
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
clim([0 1]*max(abs(UU(:) + 1i*VV(:))))
% xlabel('x (km)');
% ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)


if DIFFERENT_SOUTHWARD_CURRENTS
    AX5 = nexttile(5); % Geostrophic cyclonic eddy shifted by weak meanflow
    pcolor_centered(XX/1000,YY/1000,HH); hold on
    Stream_points = {[-57 -62 -67 -72 -77 -12], ...
                     [  0   0   0   0   0   0]}; 
        VVa = - 0.1;
        UU = UUg; VV = VVg + VVa; % Add a weak southward current to geostrophic flow
    SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
         set(SL,'Color','w','LineWidth',1);
    SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
         set(SL,'Color','w','LineWidth',1);
    quiver(XX(1:30:end,1:30:end)/1000,...
           YY(1:30:end,1:30:end)/1000,...
           UU(1:30:end,1:30:end),...
           VV(1:30:end,1:30:end),1,"filled",'Color','k')
    text(-50,-43,['$\vec{u} = \vec{u}_\mathrm{g}$ - (' num2str(abs(VVa)) ' m/s)$\hat{y}$'],...
         'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','BackgroundColor',BG_white)
    % annotation('textbox','String','12345','FitBoxToText','on');
    plot(-50,0,'w+','LineWidth',1.5)
    axis equal; axis tight
    colormap('turbo')
    clim([-1 1]*abs(HH0))
    xlabel('x (km)');
    ylabel('y (km)')
    xticks(-100:25:0); yticks(-50:25:50)
    
    AX6 = nexttile(6); % Geostrophic cyclonic eddy shifted by strong meanflow
    pcolor_centered(XX/1000,YY/1000,HH); hold on
    Stream_points = {[-50 -50 -10 -45 -35 -28], ...
                     [  5  15  50  50   0   7]};
        VVa = - 0.4;
        % VVa = - max(abs(UU(:) + 1i*VV(:)));
        UU = UUg; VV = VVg + VVa; % Add a weak southward current to geostrophic flow
    SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
         set(SL,'Color','w','LineWidth',1);
    SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
         set(SL,'Color','w','LineWidth',1);
    quiver(XX(1:30:end,1:30:end)/1000,...
           YY(1:30:end,1:30:end)/1000,...
           UU(1:30:end,1:30:end),...
           VV(1:30:end,1:30:end),1,"filled",'Color','k')
    text(-50,-43,['$\vec{u} = \vec{u}_\mathrm{g}$ - (' num2str(abs(VVa)) ' m/s)$\hat{y}$'],...
         'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','BackgroundColor',BG_white)
    plot(-50,0,'w+','LineWidth',1.5)
    axis equal; axis tight
    colormap('turbo')
    clim([-1 1]*abs(HH0))
    xlabel('x (km)');
    % ylabel('y (km)')
    xticks(-100:25:0); yticks(-50:25:50)
else

    AX5 = nexttile(5); % Geostrophic cyclonic eddy shifted by strong meanflow
    pcolor_centered(XX/1000,YY/1000,HH); hold on
    Stream_points = {[-50 -50 -10 -45 -35 -28], ...
                     [  5  15  50  50   0   7]};
        VVa = - 0.2; % VVa = - 0.4;
        % VVa = - max(abs(UU(:) + 1i*VV(:)));
        UU = UUg; VV = VVg + VVa; % Add a weak southward current to geostrophic flow
    SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
         set(SL,'Color','w','LineWidth',1);
    SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
         set(SL,'Color','w','LineWidth',1);
    quiver(XX(1:30:end,1:30:end)/1000,...
           YY(1:30:end,1:30:end)/1000,...
           UU(1:30:end,1:30:end),...
           VV(1:30:end,1:30:end),1,"filled",'Color','k')
    text(-50,-43,['$\vec{u} = \vec{u}_\mathrm{g}$ - (' num2str(abs(VVa)) ' m/s)$\hat{y}$'],...
         'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','BackgroundColor',BG_white)
    plot(-50,0,'w+','LineWidth',1.5)
    axis equal; axis tight
    colormap('turbo')
    clim([-1 1]*abs(HH0))
    xlabel('x (km)');
    % ylabel('y (km)')
    xticks(-100:25:0); yticks(-50:25:50)

    AX6 = nexttile(6); % Geostrophic cyclonic eddy shifted by strong meanflow
    Stream_points = {[-50 -50 -10 -45 -35 -28], ...
                     [  5  15  50  50   0   7]};
        % VVa = - 0.4;
        % % VVa = - max(abs(UU(:) + 1i*VV(:)));
        % UU = UUg; VV = VVg + VVa; % Add a weak southward current to geostrophic flow
    pcolor_centered(XX/1000,YY/1000,abs(UU + 1i*VV)); hold on
    SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
         set(SL,'Color','w','LineWidth',1);
    SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
         set(SL,'Color','w','LineWidth',1);
    quiver(XX(1:30:end,1:30:end)/1000,...
           YY(1:30:end,1:30:end)/1000,...
           UU(1:30:end,1:30:end),...
           VV(1:30:end,1:30:end),1,"filled",'Color','k')
    text(-50,-43,['$\vec{u} = \vec{u}_\mathrm{g}$ - (' num2str(abs(VVa)) ' m/s)$\hat{y}$'],...
         'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','BackgroundColor',BG_white)
    plot(-50,0,'w+','LineWidth',1.5)
    axis equal; axis tight
    colormap('turbo')
    clim([0 1]*max(abs(UU(:) + 1i*VV(:))))
    xlabel('x (km)');
    % ylabel('y (km)')
    xticks(-100:25:0); yticks(-50:25:50)
end







text(AX1,-90,40,'(a)','FontSize',20,'Interpreter','LaTeX','Color','k','HorizontalAlignment','center','BackgroundColor',BG_white)
text(AX2,-90,40,'(b)','FontSize',20,'Interpreter','LaTeX','Color','w','HorizontalAlignment','center','BackgroundColor',BG_black)
text(AX3,-90,40,'(c)','FontSize',20,'Interpreter','LaTeX','Color','k','HorizontalAlignment','center','BackgroundColor',BG_white)
text(AX4,-90,40,'(d)','FontSize',20,'Interpreter','LaTeX','Color','w','HorizontalAlignment','center','BackgroundColor',BG_black)
text(AX5,-90,40,'(e)','FontSize',20,'Interpreter','LaTeX','Color','k','HorizontalAlignment','center','BackgroundColor',BG_white)
text(AX6,-90,40,'(f)','FontSize',20,'Interpreter','LaTeX','Color','k','HorizontalAlignment','center','BackgroundColor',BG_white)
% text(AX5,-90,40,'(g)','FontSize',20,'Interpreter','LaTeX','Color','k','HorizontalAlignment','center','BackgroundColor',BG_white)
% text(AX6,-90,40,'(h)','FontSize',20,'Interpreter','LaTeX','Color','k','HorizontalAlignment','center','BackgroundColor',BG_white)

% set(gcf,'Position',[388          63        1053         734])
set(gcf,'Position',[951    63   490   734])

%%
% % % figure(1)
% % % exportgraphics(gcf,...
% % % ['/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/FSup_HFRmodel_vector.pdf'],...
% % % 'BackgroundColor','white','ContentType','vector')

% figure(1)
% exportgraphics(gcf,...
% ['/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/FSup_HFRmodel_raster.pdf'],...
% 'BackgroundColor','none','ContentType','image')

%% Divergence Supplemental Figure for Publication

close all

figure('Color','w')
tiledlayout(2,3,'TileSpacing','tight')
FONTSIZE = 16;

nexttile % Geostrophic cyclonic eddy
Stream_points = {[-50 -50 -50 -50 -50], ...
                 [  5  10  15  20  25]};
    UU = UUg; VV = VVg; % Only geostrophic
    DIV = divergence(XX,YY,UU,VV);
pcolor_centered(XX/1000,YY/1000,DIV/ff); hold on
SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU(1:30:end,1:30:end),...
       VV(1:30:end,1:30:end),1,"filled",'Color','k')
text(-30,-43,'$\vec{u} = \vec{u}_\mathrm{g}$','Interpreter','latex','FontSize',FONTSIZE)
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
% clim([-1 1]*abs(HH0))
xlabel('x (km)'); ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)
CB = colorbar; CB.Label.String = '$\nabla\cdot\vec{u}/f$'; % (s$^{-1}$)
CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% CB.Ticks = [-1 0 1]*abs(HH0);
CB.Location = "northoutside";

nexttile % Geostrophic cyclonic eddy shifted by weak meanflow
Stream_points = {[-50 -50 -50 -50 -10 -45], ...
                 [  5  10  15  20  50  50]}; 
    VVa = - 0.1;
    UU = UUg; VV = VVg + VVa; % Add a weak southward current to geostrophic flow
    DIV = divergence(XX,YY,UU,VV);
pcolor_centered(XX/1000,YY/1000,DIV/ff); hold on
SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU(1:30:end,1:30:end),...
       VV(1:30:end,1:30:end),1,"filled",'Color','k')
text(-50,40,['$\vec{u} = \vec{u}_\mathrm{g}$ - (' num2str(abs(VVa)) ' m/s)$\hat{y}$'],...
     'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center')
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
% clim([-1 1]*abs(HH0))
xlabel('x (km)'); ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)
CB = colorbar; CB.Label.String = '$\nabla\cdot\vec{u}/f$'; % (s$^{-1}$)
CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% CB.Ticks = [-1 0 1]*abs(HH0);
CB.Location = "northoutside";

nexttile % Geostrophic cyclonic eddy shifted by strong meanflow
Stream_points = {[-50 -50 -10 -45 -35 -28], ...
                 [  5  15  50  50   0   7]};
    VVa = - 0.4;
    % VVa = - max(abs(UU(:) + 1i*VV(:)));
    UU = UUg; VV = VVg + VVa; % Add a weak southward current to geostrophic flow
    DIV = divergence(XX,YY,UU,VV);
pcolor_centered(XX/1000,YY/1000,DIV/ff); hold on
SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU(1:30:end,1:30:end),...
       VV(1:30:end,1:30:end),1,"filled",'Color','k')
text(-50,40,['$\vec{u} = \vec{u}_\mathrm{g}$ - (' num2str(abs(VVa)) ' m/s)$\hat{y}$'],...
     'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center')
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
% clim([-1 1]*abs(HH0))
xlabel('x (km)'); ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)
CB = colorbar; CB.Label.String = '$\nabla\cdot\vec{u}/f$'; % (s$^{-1}$)
CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% CB.Ticks = [-1 0 1]*abs(HH0);
CB.Location = "northoutside";

% % % % % % % % % % %

nexttile % Geostrophic cyclonic eddy, velocity as the color quantity
% UU = UUg; VV = VVg; % Only geostrophic
% pcolor_centered(XX/1000,YY/1000,abs(UU + 1i*VV)); hold on
% quiver(XX(1:30:end,1:30:end)/1000,...
%        YY(1:30:end,1:30:end)/1000,...
%        UU(1:30:end,1:30:end),...
%        VV(1:30:end,1:30:end),1,"filled",'Color','w')
% text(-30,-43,'$\vec{u} = \vec{u}_\mathrm{g}$','Interpreter','latex','FontSize',FONTSIZE,'Color','w')
% plot(-50,0,'w+','LineWidth',1.5)
% axis equal; axis tight
% colormap('turbo')
% clim([0 1]*max(abs(UU(:) + 1i*VV(:))))
% xlabel('x (km)'); ylabel('y (km)')
% xticks(-100:25:0); yticks(-50:25:50)
% CB = colorbar; CB.Label.String = '$|\vec{u}|$ (m/s)';
% CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% CB.Ticks = [0:0.1:0.4];
% CB.Location = "southoutside";


nexttile % Geostrophic cyclonic eddy, offshore diminishment from gridding
UU = UUg; VV = VVg; % Only geostrophic
% % % 
SEARCH_RADIUS = 10000; % search radius for fitting (m)
UU_grid = nan(size(XX_grid));
VV_grid = nan(size(XX_grid));
DOPX_grid = nan(size(XX_grid));
DOPY_grid = nan(size(XX_grid));
HDOP_grid = nan(size(XX_grid));
for ii = 1:size(XX_grid,1)
    for jj = 1:size(XX_grid,2)
        FitCenter = XX_grid(ii,jj) + 1i*YY_grid(ii,jj);
        XY_Radial_ij = [];
        VelRadial_ij = [];
        for kk = 1:length(Vel_Radial)
            IND = abs([Radial_Locations{kk} - FitCenter]) < SEARCH_RADIUS;
            if sum(IND) == 0
            else
                XY_Radial_ij = [XY_Radial_ij ; Radial_Locations{kk}(IND)];
                % ^ This is actually unused by the algorithm, as it is
                % unweighted by position (I might not do that if I designed
                % the system, but maybe that would just reduce the amount
                % of usable data?).
                VelRadial_ij = [VelRadial_ij ; Vel_Radial{kk}(IND)];
            end
        end
        % % Default method, mean in search area:
        % [uGrid,vGrid,dopx,dopy,hdop] = uwlsTotal( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi );

        % Alternative method, planat or quadratic in search area:
        [uGrid,vGrid,dopx,dopy,hdop] = uwlsTotal_modified( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi , real(XY_Radial_ij), imag(XY_Radial_ij) , real(FitCenter) , imag(FitCenter) , SEARCH_RADIUS );

        UU_grid(ii,jj) = uGrid;
        VV_grid(ii,jj) = vGrid;
        DOPX_grid(ii,jj) = dopx;
        DOPY_grid(ii,jj) = dopy;
        HDOP_grid(ii,jj) = hdop;
    end
end
% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid,XX,YY,'spline');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid,XX,YY,'spline');
DIV = divergence(XX,YY,UU_grid_interp,VV_grid_interp);
% % % 
pcolor_centered(XX/1000,YY/1000,DIV/ff); hold on
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU_grid_interp(1:30:end,1:30:end),...
       VV_grid_interp(1:30:end,1:30:end),1,"filled",'Color','w')
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-75,-43,'$\vec{u} = $ gridded $\vec{u}_\mathrm{g}$','Interpreter','latex','FontSize',FONTSIZE,'Color','w')
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
% clim([0 1]*max(abs(UU(:) + 1i*VV(:))))
xlabel('x (km)'); ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)
CB = colorbar; CB.Label.String = '$\nabla\cdot\vec{u}/f$'; % (s$^{-1}$)
CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% CB.Ticks = [-1 0 1]*abs(HH0);
CB.Location = "northoutside";


nexttile % Geostrophic cyclonic eddy falsely shifted by antenna geometry
% pcolor_centered(XX/1000,YY/1000,HH); hold on
% SL = streamline(XX/1000,YY/1000,UU_grid_interp,VV_grid_interp,[-50 -50 -50 -50 -10 ],[5 10 15 20 50 ]);
%      set(SL,'Color','w','LineWidth',1);
% SL = streamline(XX/1000,YY/1000,-UU_grid_interp,-VV_grid_interp,[-50 -50 -50 -50 -10 ],[5 10 15 20 50 ]);
%      set(SL,'Color','w','LineWidth',1);
% quiver(XX(1:30:end,1:30:end)/1000,...
%        YY(1:30:end,1:30:end)/1000,...
%        UU_grid_interp(1:30:end,1:30:end),...
%        VV_grid_interp(1:30:end,1:30:end),1,"filled",'Color','k')
% plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
% text(-75,-43,'$\vec{u} = $ gridded $\vec{u}_\mathrm{g}$','Interpreter','latex','FontSize',FONTSIZE,'Color','k')
% plot(-50,0,'w+','LineWidth',1.5)
% axis equal; axis tight
% colormap('turbo')
% clim([-1 1]*abs(HH0))
% xlabel('x (km)'); ylabel('y (km)')
% xticks(-100:25:0); yticks(-50:25:50)


%% Curl Supplemental Figure for Publication

close all

figure('Color','w')
tiledlayout(1,3,'TileSpacing','tight')
nexttile % Geostrophic cyclonic eddy
Stream_points = {[-50 -50 -50 -50 -50], ...
                 [  5  10  15  20  25]};
    UU = UUg; VV = VVg; % Only geostrophic
    CURL = curl(XX,YY,UU,VV);
    CURL_g = CURL;
    maxCURL = max(abs(CURL(:)));
pcolor_centered(XX/1000,YY/1000,CURL); hold on
SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU(1:30:end,1:30:end),...
       VV(1:30:end,1:30:end),1,"filled",'Color','k')
text(-30,-43,'$\vec{u} = \vec{u}_\mathrm{g}$','Interpreter','latex','FontSize',FONTSIZE,'Color','k')
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
clim([-1 1]*maxCURL)
xlabel('x (km)'); ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)
CB = colorbar; CB.Label.String = '$\nabla\times\vec{u}_\mathrm{g}$ (s$^{-1}$)';
CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% CB.Ticks = [-1 0 1]*abs(HH0);
CB.Location = "northoutside";


nexttile % Geostrophic cyclonic eddy, offshore diminishment from gridding
UU = UUg; VV = VVg; % Only geostrophic
% % %
% Define radial velocities
% Define each antenna's observations points
Vel_Radial = {};
for ii = 1:length(ANTENNAS)
    OBS = ANTENNAS(ii) + Radial_Template(:);
    % Vq = interp2(X,Y,V,Xq,Yq)
    UU_radinterp = interp2(XX,YY, UUg + 0, real(OBS),imag(OBS));
    VV_radinterp = interp2(XX,YY, VVg + 0*[-0.4], real(OBS),imag(OBS));
    UU_radinterp = UU_radinterp( [real(OBS) < max(XX(:))] & ...
                                 [real(OBS) > min(XX(:))] & ...
                                 [imag(OBS) < max(YY(:))] & ...
                                 [imag(OBS) > min(YY(:))] );
    VV_radinterp = VV_radinterp( [real(OBS) < max(XX(:))] & ...
                                 [real(OBS) > min(XX(:))] & ...
                                 [imag(OBS) < max(YY(:))] & ...
                                 [imag(OBS) > min(YY(:))] );
    rHat = Radial_Locations{ii} - ANTENNAS(ii);
    rHat = rHat./abs(rHat);
    Vel_Radial{ii} = [UU_radinterp.*real(rHat) + VV_radinterp.*imag(rHat)].*rHat;
end
SEARCH_RADIUS = 10000; % search radius for fitting (m)
UU_grid = nan(size(XX_grid));
VV_grid = nan(size(XX_grid));
DOPX_grid = nan(size(XX_grid));
DOPY_grid = nan(size(XX_grid));
HDOP_grid = nan(size(XX_grid));
for ii = 1:size(XX_grid,1)
    for jj = 1:size(XX_grid,2)
        FitCenter = XX_grid(ii,jj) + 1i*YY_grid(ii,jj);
        XY_Radial_ij = [];
        VelRadial_ij = [];
        for kk = 1:length(Vel_Radial)
            IND = abs([Radial_Locations{kk} - FitCenter]) < SEARCH_RADIUS;
            if sum(IND) == 0
            else
                XY_Radial_ij = [XY_Radial_ij ; Radial_Locations{kk}(IND)];
                % ^ This is actually unused by the algorithm, as it is
                % unweighted by position (I might not do that if I designed
                % the system, but maybe that would just reduce the amount
                % of usable data?).
                VelRadial_ij = [VelRadial_ij ; Vel_Radial{kk}(IND)];
            end
        end
        % % Default method, mean in search area:
        % [uGrid,vGrid,dopx,dopy,hdop] = uwlsTotal( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi );

        % Alternative method, planat or quadratic in search area:
        [uGrid,vGrid,dopx,dopy,hdop] = uwlsTotal_modified( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi , real(XY_Radial_ij), imag(XY_Radial_ij) , real(FitCenter) , imag(FitCenter) , SEARCH_RADIUS );

        UU_grid(ii,jj) = uGrid;
        VV_grid(ii,jj) = vGrid;
        DOPX_grid(ii,jj) = dopx;
        DOPY_grid(ii,jj) = dopy;
        HDOP_grid(ii,jj) = hdop;
    end
end
% Interpolate gridded velocities to the original 0.1 km grid to compare directly
UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid,XX,YY,'spline');
VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid,XX,YY,'spline');
CURL = curl(XX,YY,UU_grid_interp,VV_grid_interp);
    CURL_grid = CURL;
% % % 
pcolor_centered(XX/1000,YY/1000,CURL); hold on
SL = streamline(XX/1000,YY/1000,UU_grid_interp,VV_grid_interp,  Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
SL = streamline(XX/1000,YY/1000,-UU_grid_interp,-VV_grid_interp,Stream_points{1},Stream_points{2});
     set(SL,'Color','w','LineWidth',1);
quiver(XX(1:30:end,1:30:end)/1000,...
       YY(1:30:end,1:30:end)/1000,...
       UU_grid_interp(1:30:end,1:30:end),...
       VV_grid_interp(1:30:end,1:30:end),1,"filled",'Color','k')
plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
text(-75,-43,'$\vec{u} = $ gridded $\vec{u}_\mathrm{g}$','Interpreter','latex','FontSize',FONTSIZE,'Color','k')
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
clim([-1 1]*maxCURL)
xlabel('x (km)'); ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)
CB = colorbar; CB.Label.String = '$\nabla\times\vec{u}_\mathrm{gridded}$ (s$^{-1}$)';
CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% CB.Ticks = [-1 0 1]*abs(HH0);
CB.Location = "northoutside";

AX3 = nexttile; % Geostrophic cyclonic eddy
pcolor_centered(XX/1000,YY/1000,CURL_grid./CURL_g .* ...
                [[abs(CURL_g) > 5*10^-6]./[abs(CURL_g) > 5*10^-6]]); hold on
                % [[abs(CURL_grid./CURL_g) < 2]./[abs(CURL_grid./CURL_g) < 2]]); hold on
plot(-50,0,'w+','LineWidth',1.5)
axis equal; axis tight
colormap('turbo')
clim([-1 1]*2)
xlabel('x (km)'); ylabel('y (km)')
xticks(-100:25:0); yticks(-50:25:50)
CB = colorbar; CB.Label.String = '$\nabla\times\vec{u}_\mathrm{gridded}$ / $\nabla\times\vec{u}_\mathrm{g}$';
CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% CB.Ticks = [-1 0 1]*abs(HH0);
CB.Location = "northoutside";

set(AX3,'Color',[1 1 1]*0.7)
colormap(AX3,'bwr')

figure
imagesc(CURL_grid./CURL_g .* [[abs(CURL_g) > 5*10^-6]./[abs(CURL_g) > 5*10^-6]])
axis equal; axis tight; colormap('bwr'); clim([-1 1]*2)

% %%
% 
% 
% 
% figure('Color','w')
% tiledlayout(2,3,'TileSpacing','tight')
% FONTSIZE = 16;
% 
% nexttile % Geostrophic cyclonic eddy
% Stream_points = {[-50 -50 -50 -50 -50], ...
%                  [  5  10  15  20  25]};
%     UU = UUg; VV = VVg; % Only geostrophic
%     CURL = curl(XX,YY,UU,VV);
%     maxCURL = max(abs(CURL(:)));
% pcolor_centered(XX/1000,YY/1000,CURL); hold on
% SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
%      set(SL,'Color','w','LineWidth',1);
% SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
%      set(SL,'Color','w','LineWidth',1);
% quiver(XX(1:30:end,1:30:end)/1000,...
%        YY(1:30:end,1:30:end)/1000,...
%        UU(1:30:end,1:30:end),...
%        VV(1:30:end,1:30:end),1,"filled",'Color','k')
% text(-30,-43,'$\vec{u} = \vec{u}_\mathrm{g}$','Interpreter','latex','FontSize',FONTSIZE,'Color','k')
% plot(-50,0,'w+','LineWidth',1.5)
% axis equal; axis tight
% colormap('turbo')
% clim([-1 1]*maxCURL)
% xlabel('x (km)'); ylabel('y (km)')
% xticks(-100:25:0); yticks(-50:25:50)
% CB = colorbar; CB.Label.String = '$\nabla\times\vec{u}$ (s$^{-1}$)';
% CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% % CB.Ticks = [-1 0 1]*abs(HH0);
% CB.Location = "northoutside";
% 
% nexttile % Geostrophic cyclonic eddy shifted by weak meanflow
% Stream_points = {[-50 -50 -50 -50 -10 -45], ...
%                  [  5  10  15  20  50  50]}; 
%     VVa = - 0.1;
%     UU = UUg; VV = VVg + VVa; % Add a weak southward current to geostrophic flow
%     CURL = curl(XX,YY,UU,VV);
% pcolor_centered(XX/1000,YY/1000,CURL); hold on
% SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
%      set(SL,'Color','w','LineWidth',1);
% SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
%      set(SL,'Color','w','LineWidth',1);
% quiver(XX(1:30:end,1:30:end)/1000,...
%        YY(1:30:end,1:30:end)/1000,...
%        UU(1:30:end,1:30:end),...
%        VV(1:30:end,1:30:end),1,"filled",'Color','k')
% text(-50,40,['$\vec{u} = \vec{u}_\mathrm{g}$ - (' num2str(abs(VVa)) ' m/s)$\hat{y}$'],...
%      'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','Color','k')
% plot(-50,0,'w+','LineWidth',1.5)
% axis equal; axis tight
% colormap('turbo')
% clim([-1 1]*maxCURL)
% xlabel('x (km)'); ylabel('y (km)')
% xticks(-100:25:0); yticks(-50:25:50)
% CB = colorbar; CB.Label.String = '$\nabla\times\vec{u}$ (s$^{-1}$)';
% CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% % CB.Ticks = [-1 0 1]*abs(HH0);
% CB.Location = "northoutside";
% 
% nexttile % Geostrophic cyclonic eddy shifted by strong meanflow
% Stream_points = {[-50 -50 -10 -45 -35 -28], ...
%                  [  5  15  50  50   0   7]};
%     VVa = - 0.4;
%     % VVa = - max(abs(UU(:) + 1i*VV(:)));
%     UU = UUg; VV = VVg + VVa; % Add a weak southward current to geostrophic flow
%     CURL = curl(XX,YY,UU,VV);
% pcolor_centered(XX/1000,YY/1000,CURL); hold on
% SL = streamline(XX/1000,YY/1000,UU,VV,  Stream_points{1},Stream_points{2});
%      set(SL,'Color','w','LineWidth',1);
% SL = streamline(XX/1000,YY/1000,-UU,-VV,Stream_points{1},Stream_points{2});
%      set(SL,'Color','w','LineWidth',1);
% quiver(XX(1:30:end,1:30:end)/1000,...
%        YY(1:30:end,1:30:end)/1000,...
%        UU(1:30:end,1:30:end),...
%        VV(1:30:end,1:30:end),1,"filled",'Color','k')
% text(-50,40,['$\vec{u} = \vec{u}_\mathrm{g}$ - (' num2str(abs(VVa)) ' m/s)$\hat{y}$'],...
%      'Interpreter','latex','FontSize',FONTSIZE,'HorizontalAlignment','center','Color','k')
% plot(-50,0,'w+','LineWidth',1.5)
% axis equal; axis tight
% colormap('turbo')
% clim([-1 1]*maxCURL)
% xlabel('x (km)'); ylabel('y (km)')
% xticks(-100:25:0); yticks(-50:25:50)
% CB = colorbar; CB.Label.String = '$\nabla\times\vec{u}$ (s$^{-1}$)';
% CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% % CB.Ticks = [-1 0 1]*abs(HH0);
% CB.Location = "northoutside";
% 
% % % % % % % % % % % %
% 
% nexttile % Geostrophic cyclonic eddy, velocity as the color quantity
% 
% 
% nexttile % Geostrophic cyclonic eddy, offshore diminishment from gridding
% UU = UUg; VV = VVg; % Only geostrophic
% % % % 
% SEARCH_RADIUS = 10000; % search radius for fitting (m)
% UU_grid = nan(size(XX_grid));
% VV_grid = nan(size(XX_grid));
% DOPX_grid = nan(size(XX_grid));
% DOPY_grid = nan(size(XX_grid));
% HDOP_grid = nan(size(XX_grid));
% for ii = 1:size(XX_grid,1)
%     for jj = 1:size(XX_grid,2)
%         FitCenter = XX_grid(ii,jj) + 1i*YY_grid(ii,jj);
%         XY_Radial_ij = [];
%         VelRadial_ij = [];
%         for kk = 1:length(Vel_Radial)
%             IND = abs([Radial_Locations{kk} - FitCenter]) < SEARCH_RADIUS;
%             if sum(IND) == 0
%             else
%                 XY_Radial_ij = [XY_Radial_ij ; Radial_Locations{kk}(IND)];
%                 % ^ This is actually unused by the algorithm, as it is
%                 % unweighted by position (I might not do that if I designed
%                 % the system, but maybe that would just reduce the amount
%                 % of usable data?).
%                 VelRadial_ij = [VelRadial_ij ; Vel_Radial{kk}(IND)];
%             end
%         end
%         [uGrid,vGrid,dopx,dopy,hdop] = ...
%             uwlsTotal( abs(VelRadial_ij), angle(VelRadial_ij)*180/pi );
%         UU_grid(ii,jj) = uGrid;
%         VV_grid(ii,jj) = vGrid;
%         DOPX_grid(ii,jj) = dopx;
%         DOPY_grid(ii,jj) = dopy;
%         HDOP_grid(ii,jj) = hdop;
%     end
% end
% % Interpolate gridded velocities to the original 0.1 km grid to compare directly
% UU_grid_interp = interp2(XX_grid,YY_grid,UU_grid,XX,YY,'spline');
% VV_grid_interp = interp2(XX_grid,YY_grid,VV_grid,XX,YY,'spline');
% CURL = curl(XX,YY,UU_grid_interp,VV_grid_interp);
% % % % 
% pcolor_centered(XX/1000,YY/1000,CURL); hold on
% quiver(XX(1:30:end,1:30:end)/1000,...
%        YY(1:30:end,1:30:end)/1000,...
%        UU_grid_interp(1:30:end,1:30:end),...
%        VV_grid_interp(1:30:end,1:30:end),1,"filled",'Color','w')
% plot(real(ANTENNAS)/1000, imag(ANTENNAS)/1000, 'ro','MarkerFaceColor','r')
% text(-75,-43,'$\vec{u} = $ gridded $\vec{u}_\mathrm{g}$','Interpreter','latex','FontSize',FONTSIZE,'Color','k')
% plot(-50,0,'w+','LineWidth',1.5)
% axis equal; axis tight
% colormap('turbo')
% clim([-1 1]*maxCURL)
% xlabel('x (km)'); ylabel('y (km)')
% xticks(-100:25:0); yticks(-50:25:50)
% CB = colorbar; CB.Label.String = '$\nabla\times\vec{u}$ (s$^{-1}$)';
% CB.Label.Interpreter = 'LaTeX'; CB.Label.FontSize = 16;
% % CB.Ticks = [-1 0 1]*abs(HH0);
% CB.Location = "northoutside";
% 
% 
% nexttile % Geostrophic cyclonic eddy falsely shifted by antenna geometry



%%
%%
%%
%% Assuming quadratic approximation of the Gaussian eddy
%  and that the ageostrophic current on top of the geostrophic flow is
%  constant and southward only (i.e., u_a = 0, v_g = -|v_g|):

% Delta_x_center = x_center = x_0
% where x_0 = center of the geostrophic eddy.

Va = -0.1;
Va = [-0.25:0.01:0.25];
% SIGMA_ = [15000:100:25000];
% HH0_ = [-0.2:0.01:0];

% Delta_x_center = [(-gg*HH0 + sqrt((gg*HH0)^2 + (ff*SIGMA*Va)^2))/(ff * Va) ; ...
%                   (-gg*HH0 - sqrt((gg*HH0)^2 + (ff*SIGMA*Va)^2))/(ff * Va)   ];
Delta_x_center = (-gg*HH0 - sqrt((gg*HH0).^2 + (ff*SIGMA*Va).^2))./(ff * Va);
% Delta_x_center = (-gg*HH0 - sqrt((gg*HH0).^2 + (ff*SIGMA_*Va).^2))./(ff * Va);

figure
plot(Va,Delta_x_center,'.-')
% plot(SIGMA_,Delta_x_center,'.-')
% plot(HH0_,Delta_x_center,'.-')

%% Analytic form of the offset (using the product log function, requires "Symbolic Math Toolbox")

Va = [-0.25:0.01:0.25]';

Delta_x_center_analytic = (1i/sqrt(2)) * SIGMA * sqrt(lambertw(-([ff.*SIGMA.*Va]./[2*gg*HH0]).^2)) .* [-1 1];

figure
plot(Va,Delta_x_center,'k.-'); hold on
plot(Va,Delta_x_center_analytic,'.-')
legend('Approximation','Analytical')


%%
%%
%%
%%
%%
function [ u, v, dopx, dopy, hdop, varargout ] = uwlsTotal( rSpeed, rHeading )
% UWLSTOTAL Computes a total velocity from radials using least-squares
%
%    [ u, v, dopx, dopy, hdop ] = uwlsTotal( rSpeed, rHeading ) returns the
%    eastward (u) and northward (v) components of the total velocity
%    resulting from the unweighted least squares fit of radial velocity
%    measurements provided by rSpeed and rHeading.  The horizontal dilution
%    of precision (hdop) and corresponding eastward (dopx) and northward
%    (dopy) components are also returned as a metric for assessing the
%    quality of the solution.
%
%    INPUT:
%        rSpeed   - Column vector (n* x 1) of radial velocity magnitude
%        rHeading - Column vector (n* x 1) of radial velocity heading in
%                   degrees counterclockwise from +x (east)
%
%        *n >= 2
%
%    OUTPUT:
%        u    - Eastward total velocity
%        v    - Northward total velocity
%        dopx - Eastward dilution of precision
%        dopy - Northward dilution of precision
%        hdop - Horizontal dilution of precision
%
%    All output values are scalars and the units of eastward and northward
%    velocity are the same as the units for rSpeed.
%
%    Radial velocities are related to the total velocity by projection of
%    the eastward and northward components of the total velocity onto the
%    radial heading:
%
%        rSpeed = u*cosd(rHeading) + v*sind(rHeading)
%
%    and in matrix form as
%
%        rSpeed = X*b
%
%    where
%
%        X = [ cosd(rHeading) sind(rHeading) ]
% b = [ u; v ] %
%    The unweighted least squares solution for b is given by
%
%        b = inv( X'*X ) * X' * rSpeed
%
%    The geometric portion of the covariance matrix for b is
%
%        C = inv( X'*X )
%
%    which provides the dilution of precision along the diagonal where
%
%        dopx = sqrt( C(1,1) )
%        dopy = sqrt( C(2,2) )
%        hdop = sqrt( C(1,1) + C(2,2) )
%
%    References
%
%        Least-squares methods for the extraction of surface currents from
%        CODAR crossed-loop data: Application at ARSLOE, Lipa, B. J., D. E.
%        Barrick, IEEE J. Ocean. Eng., 1983.
%
%        Shipborne measurement of surface current fields by HF radar,
%        Gurgel, K. -W., L'Onde Electr., 1994.
%
%        On the accuracy of HF radar suface current measurements:
%        Intercomparisons with ship-based sensors, Chapman, R. D., L. K.
%        Shay, H. C. Graber, J. B. Edson, A. Karachintsev, C. L. Trump, D.
%        B. Ross, J. Geophys. Res., 1997.
%
%        Guide to GPS Positioning, Wells, D. E., et al., Can. GPS Assoc.,
%        Fredericton, N. B., Canada, 1986.
%
%    Acknowledgement
%
%         This function was adapted from the TOTCALC function (ver. 3.0, Apr
%         2004) included in the HFRadarmap toolbox (ver. 4.1, Sept 2004) written
%         by Mike Cook while at the Naval Postgraduate School Department of
%         Oceanography.
% 
%     HF-Radar Network
%     Scripps Institution of Oceanography
%     Coastal Observing Research and Development Center

X = zeros( numel( rHeading ) , 2 );
X(:,1) = cosd( rHeading );
X(:,2) = sind( rHeading );
C = inv( X' * X );
dopx = sqrt( C(1,1) );
dopy = sqrt( C(2,2) );
hdop = sqrt( C(1,1) + C(2,2) );
b = C * X' * rSpeed;
u = b(1);
v = b(2);

% Added for examining the regressor matrix:
if nargout == 6
    varargout{1} = X;
else
end

end
% Shared with Luke Kachelein by Mark Otero, May 3, 2021
% [ u, v, dopx, dopy, hdop ] = uwlsTotal([1 1]',[45 135]') % simple example
% that passes a sanity check


%% "uwlsTotal" BUT MODIFIED BY LUKE KACHELEIN FOR EXPERIMENTATION:

function [ u, v, dopx, dopy, hdop , varargout] = uwlsTotal_modified( rSpeed, rHeading , xx , yy , xx0 , yy0 , SEARCH_RADIUS, METHOD )

% % % Decide which to use:

% METHOD = 'PLANAR';
% METHOD = 'QUADRATIC';

if strcmp(METHOD,'PLANAR')
    % Planar model in x and y:
    X = zeros( numel( rHeading ) , 6 );
    X(:,1)  = cosd( rHeading );
    X(:,2)  = cosd( rHeading ) .* [xx - xx0]/SEARCH_RADIUS;
    X(:,3)  = cosd( rHeading ) .* [yy - yy0]/SEARCH_RADIUS;
    X(:,4)  = sind( rHeading );
    X(:,5)  = sind( rHeading ) .* [xx - xx0]/SEARCH_RADIUS;
    X(:,6)  = sind( rHeading ) .* [yy - yy0]/SEARCH_RADIUS;
    C = inv( X' * X );
    dopx = sqrt( C(1,1) );
    dopy = sqrt( C(4,4) );
    hdop = sqrt( C(1,1) + C(4,4) );
    b = C * X' * rSpeed;
    u = b(1);
    v = b(4);
elseif strcmp(METHOD,'QUADRATIC')
    % Quadratic model in x and y:
    X = zeros( numel( rHeading ) , 12 );
    X(:,1)  = cosd( rHeading );
    X(:,2)  = cosd( rHeading ) .* [xx - xx0]/SEARCH_RADIUS;
    X(:,3)  = cosd( rHeading ) .* [yy - yy0]/SEARCH_RADIUS;
    X(:,4)  = cosd( rHeading ) .* [xx - xx0].^2 /SEARCH_RADIUS^2;
    X(:,5)  = cosd( rHeading ) .* [yy - yy0].^2 /SEARCH_RADIUS^2;
    X(:,6)  = cosd( rHeading ) .* [xx - xx0].*[yy - yy0] /SEARCH_RADIUS^2;
    X(:,7)  = sind( rHeading );
    X(:,8)  = sind( rHeading ) .* [xx - xx0]/SEARCH_RADIUS;
    X(:,9)  = sind( rHeading ) .* [yy - yy0]/SEARCH_RADIUS;
    X(:,10) = sind( rHeading ) .* [xx - xx0].^2 /SEARCH_RADIUS^2;
    X(:,11) = sind( rHeading ) .* [yy - yy0].^2 /SEARCH_RADIUS^2;
    X(:,12) = sind( rHeading ) .* [xx - xx0].*[yy - yy0] /SEARCH_RADIUS^2;
    C = inv( X' * X );
    dopx = sqrt( C(1,1) );
    dopy = sqrt( C(7,7) );
    hdop = sqrt( C(1,1) + C(7,7) );
    b = C * X' * rSpeed;
    u = b(1);
    v = b(7);
else
end
% Added for examining the regressor matrix:
if nargout == 6
    varargout{1} = X;
else
end

end

%% Planar and quadratic functions, shorthand

function [ u, v, dopx, dopy, hdop] = uwlsTotal_planar( rSpeed, rHeading , xx , yy , xx0 , yy0 , SEARCH_RADIUS)
[ u, v, dopx, dopy, hdop ] = ...
    uwlsTotal_modified( rSpeed, rHeading , xx , yy , xx0 , yy0 , SEARCH_RADIUS, 'PLANAR' );
end

function [ u, v, dopx, dopy, hdop ] = uwlsTotal_quadratic( rSpeed, rHeading , xx , yy , xx0 , yy0 , SEARCH_RADIUS)
[ u, v, dopx, dopy, hdop ] = ...
    uwlsTotal_modified( rSpeed, rHeading , xx , yy , xx0 , yy0 , SEARCH_RADIUS, 'QUADRATIC' );
end

%% PCOLOR_CENTERED (first 3 basic pcolor arguments)
% 
% Makes pcolor actually plot the matrix you give it (e.g. give it a 10x10
% matrix, plot a 10x10 matrix, not a 9x9 matrix).

% function PC = pcolor_centered(X,Y,D)
function PC = pcolor_centered(varargin)
if nargin == 3
X = varargin{1};
Y = varargin{2};
D = varargin{3};
dx = X(1,2) - X(1,1);
X_ = [X , X(:,end) + dx ; ...
      X(end,:) , X(end,end) + dx];
dy = Y(2,1) - Y(1,1);
Y_ = [Y , Y(:,end) ; ...
      Y(end,:) + dy , Y(end,end) + dy];
elseif nargin == 1
D = varargin{1};
[X,Y] = meshgrid(1:size(D,2), 1:size(D,1));
dx = 1;
X_ = [X , X(:,end) + dx ; ...
      X(end,:) , X(end,end) + dx];
dy = 1;
Y_ = [Y , Y(:,end) ; ...
      Y(end,:) + dy , Y(end,end) + dy];
else
    error('1 or 3 inputs expected')
end
D_ = [D , nan(size(D,1),1) ; ...
      nan(1,size(D,2)) , nan ];
PC = pcolor(X_ - dx/2, Y_ - dy/2, D_);
shading flat
end
