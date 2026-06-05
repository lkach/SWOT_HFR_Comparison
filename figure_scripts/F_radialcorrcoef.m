% Examine radial data directly

% Determine antenna curator and name at:
% https://cordc.ucsd.edu/projects/hfrnet/
% after selecting Overlays -> Station Placemarks

% Download data from:
% https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:IOOS-HFRadarRadial
% HTTPS:
% https://www.ncei.noaa.gov/data/oceans/ndbc/hfradar/radial/
% The radial that sees the southern eddy is BML/PAFS (May given below):
% https://www.ncei.noaa.gov/data/oceans/ndbc/hfradar/radial/2023/202305/BML/PAFS/
% A single data file is structured like this:
% https://www.ncei.noaa.gov/data/oceans/ndbc/hfradar/radial/2023/202305/BML/PAFS/RDL_i_BML_PAFS_2023_05_09_0200.ruv

% Point Arena Field Station, CA
% Station ID: PAFS
% Affiliation: BML
% Coords: 38.9284, -123.7278
% N 38 55.7060, W 123 43.6660
% 
% Ctr Freq: 4.460 MHz
% Time: 2025-05-05 20:21:44 UTC
% Age: 1:16 (HH:MM)
% Format: ruv
% Station Diagnostics

% Give where this script is regardles of how it is run,
% credit to <https://www.mathworks.com/matlabcentral/
% answers/81148-get-path-from-running-script>.
mfile = mfilename('fullpath');
if contains(mfile,'LiveEditorEvaluationHelper')
    mfile = matlab.desktop.editor.getActiveFilename;
end
mfileDashInd = mfile == '/';
mfolder = mfile(1:max(mfileDashInd.*[1:length(mfile)]));
cd(mfolder); cd ..;

%% Load SWOT/Radial CorrCoef Data

TRIN_latlon = 1i*41.0736 - 124.1578;
GCVE_latlon = 1i*38.5672 - 123.3315;
SMOA_latlon = 1i*40.7688 - 124.2188;
SHEL_latlon = 1i*40.0334 - 124.0789;
BRAG_latlon = 1i*39.4380 - 123.8161;
PAFS_latlon = 1i*38.9284 - 123.7278;

PATH = './data/SWOT_HFR_CorrCoef/';
% DIR = dir([PATH 'Radial_vel_comparison_*_quadfit*pix.mat']);
DIR = dir([PATH 'Radial_vel_comparison_*_quadfit*pix_weighted.mat']);

for ii = 1:length(DIR)
    eval(['SWOTRadCC.' DIR(ii).name(23:26) '.' replace(DIR(ii).name,{'Radial_vel_comparison_','quadfit','.mat'},{'','',''}) ...
          ' = load(''' DIR(ii).folder '/' DIR(ii).name ''');']);
end

% % % % % % % Example:
% SWOTRadCC.GCVE.GCVE_2pix
%                  Antenna_ID: 'GCVE'
%                   KM_search: 3
%                      PIXELS: 2
%             Radial_cg_corrP: [285×69 double]
%             Radial_cg_corrR: [285×69 double]
%     Radial_cg_lowpass_corrP: [285×69 double]
%     Radial_cg_lowpass_corrR: [285×69 double]
%      Radial_cg_lowpass_rmsd: [285×69 double]
%              Radial_cg_rmsd: [285×69 double]
%              Radial_g_corrP: [285×69 double]
%              Radial_g_corrR: [285×69 double]
%      Radial_g_lowpass_corrP: [285×69 double]
%      Radial_g_lowpass_corrR: [285×69 double]
%       Radial_g_lowpass_rmsd: [285×69 double]
%               Radial_g_rmsd: [285×69 double]

%% Load SWOT

% % % First and last SWOT flyover average times:
% '08-Apr-2023 07:18:35'
% '09-Jul-2023 16:47:16'
NORCAL.SWOT = load('./data/NORCAL_SWOTdata_CCS.mat');
NORCAL_SWOT_cell = fields(NORCAL.SWOT);
NORCAL.SWOT.mean_time = nan(length(NORCAL.SWOT.time),1);
for ii = 1:length(NORCAL.SWOT.mean_time)
    NORCAL.SWOT.mean_time(ii) = mean(NORCAL.SWOT.time{ii}(:),'omitnan');
end
for jj = 1:length(NORCAL_SWOT_cell)
    if strcmp(NORCAL_SWOT_cell{jj},'mean_time')
    else
        for ii = 1:length(NORCAL.SWOT.mean_time)
            FIELD = NORCAL_SWOT_cell{jj};
            eval(['NORCAL.SWOT.' FIELD '{ii} = flip(NORCAL.SWOT.' FIELD '{ii}'');']);
        end
    end
end

% Remove unused variables (comment lines if they are needed later):
NORCAL.SWOT = rmfield(NORCAL.SWOT,'time');
NORCAL.SWOT = rmfield(NORCAL.SWOT,'mean_sea_surface_cnescls');

whos

Collimate = @(IN) IN(:);
gg = 9.807; % m s^-2

SWOT_polygon = [flip(NORCAL.SWOT.lon{1}(:,1)) ; NORCAL.SWOT.lon{1}(1,:)' ; NORCAL.SWOT.lon{1}(:,end) ; flip(NORCAL.SWOT.lon{1}(end,:))'] + ...
               [flip(NORCAL.SWOT.lat{1}(:,1)) ; NORCAL.SWOT.lat{1}(1,:)' ; NORCAL.SWOT.lat{1}(:,end) ; flip(NORCAL.SWOT.lat{1}(end,:))']*1i;

% figure
% scatter(real(SWOT_polygon),imag(SWOT_polygon),10,1:length(SWOT_polygon),'filled');colormap(turbo)

%% Calculate SWOT geostrophic velocity

if true
    warning(['This section is included for completeness ' ...
             '(to refer to how geostrophic velocity is ' ...
             'calculated), but does not need to be run'])
else


% Recall that at the surface, geostrophic velocity is:
% u_g = -(g/f)(d\eta/dy)
% v_g = +(g/f)(d\eta/dx)
% where both derivatives are partial derivatives.
% Simple textbook citation:
% https://uw.pressbooks.pub/ocean285/chapter/geostrophic-balance/

close all

SWOT_geostr_vel_calculated = true;
gg = 9.807; % m s^-2
GaussianKernelSTD_meters = 4000;
PIXELS = 4;
% Number of pixels in one direction for the quadratic fit, e.g.:
% PIXELS = 6 -> 13x13 box, i.e. [26 km]^2

U_geostr = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
V_geostr = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% U_geostr_ = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% V_geostr_ = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
for ti_swot = 1:length(NORCAL.SWOT.ssha_karin_2)
    
    % % % Calculate ssha to obtain geostrophic velocity:
    SSHA = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
           [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
           [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
    

    % % % % % Smoothing is necessary to reduce high-k noise that can really throw
    % % % % % off the velocity when taking the first derivative.


    % % % Less good Gaussian smoothing procedures:
    % SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');
    % SSHA = smoothdata2(SSHA,'movmedian',7,'omitnan');
    % SSHA = smoothdata2_interpnadirgap(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, SSHA, 2000,GaussianKernelSTD_meters);
    

    % % % Isotropic Gaussian blur (best Gaussian blur version so far):
    % SSHA = bettergaussiansmooth2(SSHA,2000,GaussianKernelSTD_meters);
    % SSHA = SSHA .* ...
    %    [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
    %    [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
    % % % Built-in gradientm function (requires mapping toolbox, better):
    % [~,~,gradSSHA_y, gradSSHA_x] = gradientm(NORCAL.SWOT.lat{ti_swot},NORCAL.SWOT.lon{ti_swot},SSHA);
    % SWOT_ug = -gg*gradSSHA_y./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
    %     SWOT_ug = SWOT_ug .* [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}]./[~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}];
    % SWOT_vg =  gg*gradSSHA_x./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
    %     SWOT_vg = SWOT_vg .* [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}]./[~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}];
    % 
    % U_geostr(:,:,ti_swot) = SWOT_ug;
    % V_geostr(:,:,ti_swot) = SWOT_vg;


    % % % My custom gradient function (not as good as gradientm):
    % gradSSHA = gradient_2D(SSHA,NORCAL.SWOT.lon{ti_swot}.*cosd(NORCAL.SWOT.lat{ti_swot}),...
    %                             NORCAL.SWOT.lat{ti_swot})*[1/111000];
    % gradSSHA_x = gradSSHA(:,:,1); gradSSHA_y = gradSSHA(:,:,2);
    % SWOT_ug = -gg*gradSSHA_y./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
    %     SWOT_ug = SWOT_ug .* [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}]./[~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}];
    % SWOT_vg =  gg*gradSSHA_x./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
    %     SWOT_vg = SWOT_vg .* [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}]./[~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}];
    % 
    % U_geostr_(:,:,ti_swot) = SWOT_ug;
    % V_geostr_(:,:,ti_swot) = SWOT_vg;

    % % % Obtain geostrophic velocity directly using the quadratic fit
    % % % method of Tranchant et al. (2025):
    [AA,~] = fit_6term_2d(NORCAL.SWOT.lon{ti_swot}, ...
                          NORCAL.SWOT.lat{ti_swot}, ...
                          SSHA, PIXELS, true); % 6 -> 13x13 box, i.e. [26 km]^2
                          % PIXELS -> [2*PIXELS + 1]x[2*PIXELS + 1] box
    % To get geostrophic velocity from inputting surface height:
    U_geostr(:,:,ti_swot) = -(gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,3))/111000;
    V_geostr(:,:,ti_swot) =  (gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,2))./[111000*cosd(NORCAL.SWOT.lat{ti_swot})];


    disp(100*ti_swot/length(NORCAL.SWOT.ssha_karin_2))
end

end

%% Apply the iterative method of Penven et al. (2014) to get the Cyclogeostrophic current

if true
    warning(['This section is included for completeness ' ...
             '(to refer to how cyclogeostrophic velocity is ' ...
             'calculated), but does not need to be run'])
else

% https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2013JC009528

% For \vec{u} = u + iv:
% \vec{u}^{(n+1)} = \vec{u}_g + \frac{\hat{k}}{f}\times(\vec{u}^{(n)}\cdot\nabla\vec{u}^{(n)})

% Recall that {0,0,1}x{a,b,c} = {-b,a,0}

% Vetor equation broken down into 2 scalar equations:
% u^{(n+1)} = u_g + \frac{1}{f}(-u^{(n)}\partial_x v^{(n)} - v^{(n)}\partial_y v^{(n)})
% v^{(n+1)} = v_g + \frac{1}{f}( u^{(n)}\partial_x u^{(n)} + v^{(n)}\partial_y u^{(n)})

% Because each iteration shrinks the borders, just do one step.

one_over_f = [1./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})];

U_cyclogeostr_1it = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
V_cyclogeostr_1it = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
for ti_swot = 1:size(U_geostr,3)
    U_n = U_geostr(:,:,ti_swot);
    V_n = V_geostr(:,:,ti_swot);
    for CG_Iteration = 1:1
        dxU_n = divergence(NORCAL.SWOT.lon{ti_swot}*111000,NORCAL.SWOT.lat{ti_swot}*111000,...
                           U_n,0*U_n);
        dyU_n = divergence(NORCAL.SWOT.lon{ti_swot}*111000,NORCAL.SWOT.lat{ti_swot}*111000,...
                           0*U_n,U_n);
        dxV_n = divergence(NORCAL.SWOT.lon{ti_swot}*111000,NORCAL.SWOT.lat{ti_swot}*111000,...
                           V_n,0*V_n);
        dyV_n = divergence(NORCAL.SWOT.lon{ti_swot}*111000,NORCAL.SWOT.lat{ti_swot}*111000,...
                           0*V_n,V_n);
    
        U_n1 = U_geostr(:,:,ti_swot) - one_over_f.*(U_n.*dxV_n + V_n.*dyV_n);
        V_n1 = V_geostr(:,:,ti_swot) + one_over_f.*(U_n.*dxU_n + V_n.*dyU_n);
    
        U_n = U_n1;
        V_n = V_n1;
    end
    U_cyclogeostr_1it(:,:,ti_swot) = U_n1;
    V_cyclogeostr_1it(:,:,ti_swot) = V_n1;
    disp(100*ti_swot/size(U_geostr,3))
end

U_cyclogeostr_2it = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
V_cyclogeostr_2it = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
for ti_swot = 1:size(U_geostr,3)
    U_n = U_geostr(:,:,ti_swot);
    V_n = V_geostr(:,:,ti_swot);
    for CG_Iteration = 1:2
        dxU_n = divergence(NORCAL.SWOT.lon{ti_swot}*111000,NORCAL.SWOT.lat{ti_swot}*111000,...
                           U_n,0*U_n);
        dyU_n = divergence(NORCAL.SWOT.lon{ti_swot}*111000,NORCAL.SWOT.lat{ti_swot}*111000,...
                           0*U_n,U_n);
        dxV_n = divergence(NORCAL.SWOT.lon{ti_swot}*111000,NORCAL.SWOT.lat{ti_swot}*111000,...
                           V_n,0*V_n);
        dyV_n = divergence(NORCAL.SWOT.lon{ti_swot}*111000,NORCAL.SWOT.lat{ti_swot}*111000,...
                           0*V_n,V_n);
    
        U_n1 = U_geostr(:,:,ti_swot) - one_over_f.*(U_n.*dxV_n + V_n.*dyV_n);
        V_n1 = V_geostr(:,:,ti_swot) + one_over_f.*(U_n.*dxU_n + V_n.*dyU_n);
    
        U_n = U_n1;
        V_n = V_n1;
    end
    U_cyclogeostr_2it(:,:,ti_swot) = U_n1;
    V_cyclogeostr_2it(:,:,ti_swot) = V_n1;
    disp(100*ti_swot/size(U_geostr,3))
end

CG_Iteration = 1;
U_cyclogeostr_Nit = U_cyclogeostr_1it;
V_cyclogeostr_Nit = V_cyclogeostr_1it;

% %% Compare G, 1-iteration, and 2-iterations:
% figure;
% subplot(131); imagesc(sqrt(U_geostr(:,:,ti_swot).^2 + V_geostr(:,:,ti_swot).^2)); axis equal; colorbar; clim([0 2]); colormap(turbo)
% subplot(132); imagesc(sqrt(U_cyclogeostr_1it(:,:,ti_swot).^2 + V_cyclogeostr_1it(:,:,ti_swot).^2)); axis equal; colorbar; clim([0 2])
% subplot(133); imagesc(sqrt(U_cyclogeostr_2it(:,:,ti_swot).^2 + V_cyclogeostr_2it(:,:,ti_swot).^2)); axis equal; colorbar; clim([0 2])

end

error('Forced to stop by user.')

%% Map RMSD

% SWOTRadCC.GCVE.GCVE_2pix
%                  Antenna_ID: 'GCVE'
%                   KM_search: 3
%                      PIXELS: 2
%             Radial_cg_corrP: [285×69 double]
%             Radial_cg_corrR: [285×69 double]
%     Radial_cg_lowpass_corrP: [285×69 double]
%     Radial_cg_lowpass_corrR: [285×69 double]
%      Radial_cg_lowpass_rmsd: [285×69 double]
%              Radial_cg_rmsd: [285×69 double]
%              Radial_g_corrP: [285×69 double]
%              Radial_g_corrR: [285×69 double]
%      Radial_g_lowpass_corrP: [285×69 double]
%      Radial_g_lowpass_corrR: [285×69 double]
%       Radial_g_lowpass_rmsd: [285×69 double]
%               Radial_g_rmsd: [285×69 double]

Antenna_list = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};
Antenna_ID = Antenna_list{4};
PIXELS_list = [2 4 6 8 10 12];
PLOT_MATRICES = {'Radial_g_rmsd','Radial_g_lowpass_rmsd','Radial_cg_rmsd','Radial_cg_lowpass_rmsd'};
HIST_COLOR = {'r','g','b','k'};
P_threshold = 0.05;

close all
for jj = 1:4
    figure('Color','w')
    tiledlayout(2,3,'TileSpacing','tight')
end
figure('Color','w')

for jj = 1:4
    PLOT_MATRIX = PLOT_MATRICES{jj};
    figure(jj)
    for ii = 1:length(PIXELS_list)
        eval(['AX' num2str(ii) ' = nexttile;'])
        M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                     'longitudes',[-126 -122.] + [-.1 .1],...
                     'latitudes',[37 42] + [-.1 .1]);
        set(gcf,'color','w'); hold on
        % m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
        %                   eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.' PLOT_MATRIX]));
        m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                          eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.' PLOT_MATRIX]));
        COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
        m_grid('box','fancy', 'backgroundcolor','none');
        m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',15)
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1,imag(eval([Antenna_ID '_latlon'])),Antenna_ID)
        if ii == 1
            title({replace(PLOT_MATRIX,'_','\_');[num2str([2*PIXELS_list(ii) + 1]*2) ' km fitting radius']})
        else
            title({[num2str([2*PIXELS_list(ii) + 1]*2) ' km fitting radius']})
        end
        set(gca,'FontSize',16)
        CB = colorbar; clim([0 1]*1)
        CB.Label.String = 'RMSD (m/s)'; CB.Label.FontSize = 16;
    end
    colormap(turbo)
    linkaxes([AX1 AX2 AX3 AX4 AX5 AX6],'xy')

    figure(5)
    % histogram(eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.' PLOT_MATRIX]), ...
    %           [0:0.01:1], 'Normalization','pdf', 'DisplayStyle', 'stairs', 'EdgeColor',HIST_COLOR{jj}, 'LineWidth',0.5*3); hold on
    histogram(eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.' PLOT_MATRIX]), ...
              [0:0.01:1], 'Normalization','pdf', 'DisplayStyle', 'stairs', 'EdgeColor',HIST_COLOR{jj}, 'LineWidth',0.5*3); hold on
end



%% Map correlation (geostrophic vs. full HFR)


Antenna_list = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};
PIXELS_list = [2 4 6 8 10 12];
Antenna_ID = Antenna_list{5};


close all
figure('Color','w')
tiledlayout(2,3,'TileSpacing','tight')
% tiledlayout(1,4,'TileSpacing','tight')

P_threshold = 0.05;

for ii = 1:length(PIXELS_list)
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color',[1 1 1]*0.7); hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.Radial_g_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.Radial_g_corrP']);
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',15)
    title({Antenna_ID ; [num2str([2*PIXELS_list(ii) + 1]*2) ' km fitting radius']})
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1]*1)
    CB.Label.String = ['CC, SWOT_g vs HFR (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
end

colormap(bwr)

%% Map correlation (cyclogeostrophic vs. full HFR)

Antenna_list = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};
PIXELS_list = [2 4 6 8 10 12];
HIST_COLOR = {'r','g','b','m','y','c'};

Antenna_ID = Antenna_list{4};


close all
figure('Color','w')
tiledlayout(2,3,'TileSpacing','tight')
% tiledlayout(1,4,'TileSpacing','tight')

P_threshold = 0.05;

figure('Color','w')
HIST_BINS = [-1:0.05:1];

disp(' ')

for ii = 1:length(PIXELS_list)
    figure(1)
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color',[1 1 1]*0.7); hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.Radial_cg_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.Radial_cg_corrP']);
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',15)
    title({Antenna_ID ; [num2str([2*PIXELS_list(ii) + 1]*2) ' km fitting radius']})
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1]*1)
    CB.Label.String = ['CC, SWOT_{cg} vs HFR (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

    figure(2)
    % disp([num2str(PIXELS_list(ii)) ' pixels median C.C.: ' num2str( prctile(Collimate(CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]),50) )])
    % disp([num2str(PIXELS_list(ii)) ' pixels median C.C.: ' num2str( prctile(Collimate(CorrCoef_R),50) )])
    disp([num2str(PIXELS_list(ii)) ' pixels median C.C.: ' num2str( prctile(Collimate(CorrCoef_R(CorrCoef_R>0.4)),50) )])
    histogram(CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold],HIST_BINS,...
              'Normalization','pdf', 'DisplayStyle', 'stairs', 'EdgeColor',HIST_COLOR{ii}, 'LineWidth',0.5*3); hold on
end

figure(2)
set(gca,'Color',[1 1 1]*0.7)
figure(1)
colormap(bwr)

%% LOW PASS
%% Map correlation (geostrophic vs. lowpass HFR)


Antenna_list = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};
PIXELS_list = [2 4 6 8 10 12];
Antenna_ID = Antenna_list{2};


close all
figure('Color','w')
tiledlayout(2,3,'TileSpacing','tight')
% tiledlayout(1,4,'TileSpacing','tight')

P_threshold = 0.05;

for ii = 1:length(PIXELS_list)
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color',[1 1 1]*0.7); hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.Radial_g_lowpass_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.Radial_g_lowpass_corrP']);
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',15)
    title({Antenna_ID ; [num2str([2*PIXELS_list(ii) + 1]*2) ' km fitting radius']})
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1]*1)
    CB.Label.String = ['CC, SWOT_g vs HFR_{lp} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
end

colormap(bwr)

%% Map correlation (cyclogeostrophic vs. lowpass HFR)


Antenna_list = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};
PIXELS_list = [2 4 6 8 10 12];
Antenna_ID = Antenna_list{5};


close all
figure('Color','w')
tiledlayout(2,3,'TileSpacing','tight')
% tiledlayout(1,4,'TileSpacing','tight')

P_threshold = 0.05;

for ii = 1:length(PIXELS_list)
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color',[1 1 1]*0.7); hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.Radial_cg_lowpass_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix.Radial_cg_lowpass_corrP']);
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',15)
    title({Antenna_ID ; [num2str([2*PIXELS_list(ii) + 1]*2) ' km fitting radius']})
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1]*1)
    CB.Label.String = ['CC, SWOT_{cg} vs HFR_{lp} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
end

colormap(bwr)

%% Calculate highest correlations and lowest RMSD

% Antenna_list = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};

% Remove GCVE because it is too far from the the swath:
Antenna_list = {'TRIN','SMOA','SHEL','BRAG','PAFS'};

for ii = 1:length(PIXELS_list)
    % % % % % % % % % % % CorrCoef
    % % % % % % % GEOSTROPHIC
    % % % Full HFR
    eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(ii)) ' = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(Antenna_list));']); jj = 1;
    eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrR .* [SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold]./[SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrR .* [SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold]./[SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrR .* [SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold]./[SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrR .* [SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold]./[SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrR .* [SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold]./[SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold];']); jj = jj + 1;
    % eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrR .* [SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold]./[SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(ii)) ' = max(maxCorrCoef_g.pix' num2str(PIXELS_list(ii)) ',[],3,''omitnan'');']);
    % % % Filtered HFR
    eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(ii)) ' = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(Antenna_list));']); jj = 1;
    eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrR .* [SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold]./[SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrR .* [SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold]./[SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrR .* [SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold]./[SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrR .* [SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold]./[SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrR .* [SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold]./[SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold];']); jj = jj + 1;
    % eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrR .* [SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold]./[SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(ii)) ' = max(maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(ii)) ',[],3,''omitnan'');'])
    % % % % % % % CYCLOGEOSTROPHIC
    % % % Full HFR
    eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(ii)) ' = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(Antenna_list));']); jj = 1;
    eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrR .* [SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold]./[SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrR .* [SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold]./[SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrR .* [SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold]./[SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrR .* [SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold]./[SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrR .* [SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold]./[SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold];']); jj = jj + 1;
    % eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrR .* [SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold]./[SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(ii)) ' = max(maxCorrCoef_cg.pix' num2str(PIXELS_list(ii)) ',[],3,''omitnan'');'])
    % % % Filtered HFR
    eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(ii)) ' = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(Antenna_list));']); jj = 1;
    eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrR .* [SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold]./[SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrR .* [SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold]./[SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrR .* [SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold]./[SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrR .* [SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold]./[SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrR .* [SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold]./[SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold];']); jj = jj + 1;
    % eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrR .* [SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold]./[SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_corrP<P_threshold];']); jj = jj + 1;
    eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(ii)) ' = max(maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(ii)) ',[],3,''omitnan'');'])

    % % % % % % % % % % % RMSD
    % % % % % % % GEOSTROPHIC
    % % % Full HFR
    eval(['minRMSD_g.pix' num2str(PIXELS_list(ii)) ' = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(Antenna_list));']); jj = 1;
    eval(['minRMSD_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_rmsd;']); jj = jj + 1;
    % eval(['minRMSD_g.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g.pix' num2str(PIXELS_list(ii)) ' = min(minRMSD_g.pix' num2str(PIXELS_list(ii)) ',[],3,''omitnan'');']);
    % % % Filtered HFR
    eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(ii)) ' = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(Antenna_list));']); jj = 1;
    eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_rmsd;']); jj = jj + 1;
    % eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_g_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(ii)) ' = min(minRMSD_g_lowpass.pix' num2str(PIXELS_list(ii)) ',[],3,''omitnan'');'])
    % % % % % % % CYCLOGEOSTROPHIC
    % % % Full HFR
    eval(['minRMSD_cg.pix' num2str(PIXELS_list(ii)) ' = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(Antenna_list));']); jj = 1;
    eval(['minRMSD_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_rmsd;']); jj = jj + 1;
    % eval(['minRMSD_cg.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg.pix' num2str(PIXELS_list(ii)) ' = min(minRMSD_cg.pix' num2str(PIXELS_list(ii)) ',[],3,''omitnan'');'])
    % % % Filtered HFR
    eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(ii)) ' = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(Antenna_list));']); jj = 1;
    eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.TRIN.TRIN_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SMOA.SMOA_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.SHEL.SHEL_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.BRAG.BRAG_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.PAFS.PAFS_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_rmsd;']); jj = jj + 1;
    % eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(ii)) '(:,:,jj) = SWOTRadCC.GCVE.GCVE_' num2str(PIXELS_list(ii)) 'pix_weighted.Radial_cg_lowpass_rmsd;']); jj = jj + 1;
    eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(ii)) ' = min(minRMSD_cg_lowpass.pix' num2str(PIXELS_list(ii)) ',[],3,''omitnan'');'])

    disp(ii)
end


%% Map all antennas' highest correlations together

CMAP = bwr(100).^0.667;
% CMAP(:,2) = [linspace(0,0.5,50) linspace(0.5,0,50)]';

close all

for pixel_i = 1 % 1:length(PIXELS_list)
    figure('Color','w')
    tiledlayout(1,2,'TileSpacing','tight')
    P_threshold = 0.05;
    
    AX1 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(pixel_i))]) );
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    m_text(-122.75, 41.75, '(a)','Interpreter','latex','FontSize',32)
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1])
    % CB.Label.String = ['$u_\mathrm{radial}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    % CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT,HFR),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    % CB.Label.String = ['C(u$^\mathrm{rad}_\mathrm{g}$,u$^\mathrm{rad}_\mathrm{HFR}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';
    title(['CorrCoef(u$^\mathrm{rad}_\mathrm{g}$,u$^\mathrm{rad}_\mathrm{HFR}$),  P$<$' num2str(P_threshold) ''],'Interpreter','latex')
    
    AX2 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(pixel_i))]) );
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    m_text(-122.75, 41.75, '(b)','Interpreter','latex','FontSize',32)
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1])
    % CB.Label.String = ['$u_\mathrm{radial, lowpass}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    % CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT,HFR$_\mathrm{LP}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    % CB.Label.String = ['C(u$^\mathrm{rad}_\mathrm{g}$,u$^\mathrm{rad}_\mathrm{HFR, low pass}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';
    % title([num2str([2*PIXELS_list(pixel_i) + 1]*2) ' km fitting box'])
    title(['CorrCoef(u$^\mathrm{rad}_\mathrm{g}$,u$^\mathrm{rad}_\mathrm{HFR, low pass}$),  P$<$' num2str(P_threshold) ''],'Interpreter','latex')
    colormap(CMAP)
    % set(gcf,'Position',[-919   238   882   540])
    set(gcf,'Position',[-1546 240 950 684])
    linkaxes([AX1 AX2],'xy')
end


for pixel_i = 1 % 1:length(PIXELS_list)
    figure('Color','w')
    tiledlayout(1,2,'TileSpacing','tight')
    P_threshold = 0.05;
    
    AX1 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(pixel_i))]) );
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    m_text(-122.75, 41.75, '(a)','Interpreter','latex','FontSize',32)
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1])
    % CB.Label.String = ['$u_\mathrm{radial}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    % CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT,HFR),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    % CB.Label.String = ['C(u$^\mathrm{rad}_\mathrm{cg}$,u$^\mathrm{rad}_\mathrm{HFR}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';
    title(['CorrCoef(u$^\mathrm{rad}_\mathrm{cg}$,u$^\mathrm{rad}_\mathrm{HFR}$),  P$<$' num2str(P_threshold) ''],'Interpreter','latex')
    
    AX2 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(pixel_i))]) );
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    m_text(-122.75, 41.75, '(b)','Interpreter','latex','FontSize',32)
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1])
    % CB.Label.String = ['$u_\mathrm{radial, lowpass}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    % CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT,HFR$_\mathrm{LP}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    % CB.Label.String = ['C(u$^\mathrm{rad}_\mathrm{cg}$,u$^\mathrm{rad}_\mathrm{HFR, low pass}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';
    % title([num2str([2*PIXELS_list(pixel_i) + 1]*2) ' km fitting box'])
    title(['CorrCoef(u$^\mathrm{rad}_\mathrm{cg}$,u$^\mathrm{rad}_\mathrm{HFR, low pass}$),  P$<$' num2str(P_threshold) ''],'Interpreter','latex')
    colormap(CMAP)
    % set(gcf,'Position',[-919   238   882   540])
    set(gcf,'Position',[-1546 240 950 684])
    linkaxes([AX1 AX2],'xy')
end

% %%
% figure(1)
% exportgraphics(gcf,...
% ['/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/' ...
%  'radialcorrcoef_SWOTg_' num2str([2*PIXELS_list(pixel_i) + 1]*2) 'km_HFR_HFRlp_2panel.pdf'],...
% 'BackgroundColor','none','ContentType','vector')

%% Map the improvement in CorrCoef as the fitting window expands

CMAP = bwr(100).^0.667;
% CMAP(:,2) = [linspace(0,0.5,50) linspace(0.5,0,50)]';

close all

for pixel_i = 2:length(PIXELS_list)

    % Comparison_index = [2*PIXELS_list(pixel_i) + 1]*2;
    Comparison_index = 1;

    figure('Color','w')
    tiledlayout(1,2,'TileSpacing','tight')
    P_threshold = 0.05;
    
    AX1 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(pixel_i         ))]) - ...
                      eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(Comparison_index))]) );
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1]*0.5)
    % CB.Label.String = ['$u_\mathrm{radial}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT,HFR),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';
    title({['Difference between '];...
           [num2str([2*PIXELS_list(pixel_i         ) + 1]*2) ' km and ' ...
            num2str([2*PIXELS_list(Comparison_index) + 1]*2) ' fitting boxes']})
    
    AX2 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(pixel_i         ))]) - ...
                      eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(Comparison_index))]) );
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1]*0.5)
    % CB.Label.String = ['$u_\mathrm{radial, lowpass}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    CB.Label.String = ['$\Delta$CorrCoef$_\mathrm{radial}$(SWOT,HFR$_\mathrm{LP}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';
    colormap(CMAP)
    set(gcf,'Position',[-919   238   882   540])
    linkaxes([AX1 AX2],'xy')
end

%% CorrCoef distribution as a function of fitting window size (PUBLICATION FIGURE)

close all

prctile_range = [10 25 50 75 90];

maxCorrCoef_g.PRCTILE = [];
maxCorrCoef_cg.PRCTILE = [];
maxCorrCoef_g_lowpass.PRCTILE = [];
maxCorrCoef_cg_lowpass.PRCTILE = [];

for pixel_i = 1:length(PIXELS_list)
    
    PRCTILE_i = prctile(eval(['maxCorrCoef_g.pix' num2str(PIXELS_list(pixel_i)) '(:)']) , prctile_range)';
    maxCorrCoef_g.PRCTILE(:,pixel_i) = PRCTILE_i;

    PRCTILE_i = prctile(eval(['maxCorrCoef_cg.pix' num2str(PIXELS_list(pixel_i)) '(:)']) , prctile_range)';
    maxCorrCoef_cg.PRCTILE(:,pixel_i) = PRCTILE_i;
    
    PRCTILE_i = prctile(eval(['maxCorrCoef_g_lowpass.pix' num2str(PIXELS_list(pixel_i)) '(:)']) , prctile_range)';
    maxCorrCoef_g_lowpass.PRCTILE(:,pixel_i) = PRCTILE_i;

    PRCTILE_i = prctile(eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(pixel_i)) '(:)']) , prctile_range)';
    maxCorrCoef_cg_lowpass.PRCTILE(:,pixel_i) = PRCTILE_i;

end

figure('Color','w')
TL = tiledlayout(1,2,'TileSpacing','tight');
P_threshold = 0.05;

% Geostr. vs. unfiltered HFR
nexttile
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g.PRCTILE(3,:), 'k.-', 'LineWidth',0.5*3 ,'MarkerSize',6*3 ); hold on
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g.PRCTILE(2,:), 'k.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g.PRCTILE(4,:), 'k.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g.PRCTILE(1,:), 'k.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g.PRCTILE(5,:), 'k.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% cyclogeostr. vs. unfiltered HFR
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg.PRCTILE(3,:), 'r.-', 'LineWidth',0.5*3 ,'MarkerSize',6*3 ); hold on
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg.PRCTILE(2,:), 'r.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg.PRCTILE(4,:), 'r.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg.PRCTILE(1,:), 'r.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg.PRCTILE(5,:), 'r.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
text(44,maxCorrCoef_cg.PRCTILE(1,end) - 0.03, [num2str(prctile_range(1)) '\%ile'],'Interpreter','latex')
text(44,maxCorrCoef_cg.PRCTILE(2,end) - 0.03, [num2str(prctile_range(2)) '\%ile'],'Interpreter','latex')
text(44,maxCorrCoef_cg.PRCTILE(3,end) - 0.03, [num2str(prctile_range(3)) '\%ile'],'Interpreter','latex')
text(44,maxCorrCoef_cg.PRCTILE(4,end) - 0.03, [num2str(prctile_range(4)) '\%ile'],'Interpreter','latex')
text(44,maxCorrCoef_cg.PRCTILE(5,end) - 0.03, [num2str(prctile_range(5)) '\%ile'],'Interpreter','latex')
xlabel(TL,'Fitting Box Side Length (km)','Interpreter','latex','FontSize',16)
ylabel('CorrCoef','Interpreter','latex','FontSize',16)
set(gca,'ylim',[0 1])
LEG = legend('C(u$^\mathrm{rad}_\mathrm{geo}$, u$^\mathrm{rad}_\mathrm{HFR}$)','C(u$^\mathrm{rad}_\mathrm{cyclogeo}$, u$^\mathrm{rad}_\mathrm{HFR}$)');
LEG.Interpreter = 'latex';
LEG.FontSize = 16;

% Geostr. vs. filtered HFR
nexttile
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g_lowpass.PRCTILE(3,:), 'k.-', 'LineWidth',0.5*3 ,'MarkerSize',6*3 ); hold on
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g_lowpass.PRCTILE(2,:), 'k.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g_lowpass.PRCTILE(4,:), 'k.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g_lowpass.PRCTILE(1,:), 'k.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_g_lowpass.PRCTILE(5,:), 'k.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% cyclogeostr. vs. filtered HFR
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg_lowpass.PRCTILE(3,:), 'r.-', 'LineWidth',0.5*3 ,'MarkerSize',6*3 ); hold on
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg_lowpass.PRCTILE(2,:), 'r.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg_lowpass.PRCTILE(4,:), 'r.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg_lowpass.PRCTILE(1,:), 'r.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
plot([2*PIXELS_list + 1]*2, maxCorrCoef_cg_lowpass.PRCTILE(5,:), 'r.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
text(44,maxCorrCoef_cg_lowpass.PRCTILE(1,end) - 0.03, [num2str(prctile_range(1)) '\%ile'],'Interpreter','latex')
text(44,maxCorrCoef_cg_lowpass.PRCTILE(2,end) - 0.03, [num2str(prctile_range(2)) '\%ile'],'Interpreter','latex')
text(44,maxCorrCoef_cg_lowpass.PRCTILE(3,end) - 0.03, [num2str(prctile_range(3)) '\%ile'],'Interpreter','latex')
text(44,maxCorrCoef_cg_lowpass.PRCTILE(4,end) - 0.03, [num2str(prctile_range(4)) '\%ile'],'Interpreter','latex')
text(44,maxCorrCoef_cg_lowpass.PRCTILE(5,end) - 0.03, [num2str(prctile_range(5)) '\%ile'],'Interpreter','latex')
% xlabel('Fitting Box Side Length (km)','Interpreter','latex','FontSize',16)
% ylabel('CorrCoef','Interpreter','latex','FontSize',16)
set(gca,'ylim',[0 1])
LEG = legend('C(u$^\mathrm{rad}_\mathrm{geo}$, u$^\mathrm{rad}_\mathrm{HFR, low pass}$)','C(u$^\mathrm{rad}_\mathrm{cyclogeo}$, u$^\mathrm{rad}_\mathrm{HFR, low pass}$)');
LEG.Interpreter = 'latex';
LEG.FontSize = 16;

set(gcf,'Position',[-1078         366         599         445])

% %%
% figure(1)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/radial_SWOT_corrcoef_vs_fitwindow.pdf',...
% 'BackgroundColor','none','ContentType','vector')


%% Map all antennas' minimum RMSD together

% CMAP = bwr(100).^0.667;
CMAP = "turbo";
CLIM = [0 1]*0.5;

close all

for pixel_i = 1:length(PIXELS_list)
    figure('Color','w')
    tiledlayout(2,2,'TileSpacing','tight')
    P_threshold = 0.05;
    
    AX1 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['minRMSD_g.pix' num2str(PIXELS_list(pixel_i))]) );
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['minRMSD_g.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['minRMSD_g.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    set(gca,'FontSize',16)
    CB = colorbar; clim(CLIM)
    % CB.Label.String = ['$u_\mathrm{radial}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    CB.Label.String = ['RMSD$_\mathrm{radial}$(SWOT$_\mathrm{g}$,HFR)']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';
    title([num2str([2*PIXELS_list(pixel_i) + 1]*2) ' km fitting box'])
    
    AX2 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(pixel_i))]) );
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['minRMSD_g_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    set(gca,'FontSize',16)
    CB = colorbar; clim(CLIM)
    % CB.Label.String = ['$u_\mathrm{radial, lowpass}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    CB.Label.String = ['RMSD$_\mathrm{radial}$(SWOT$_\mathrm{g}$,HFR$_\mathrm{LP}$)']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';

    AX3 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['minRMSD_cg.pix' num2str(PIXELS_list(pixel_i))]) );
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['minRMSD_cg.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['minRMSD_cg.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    set(gca,'FontSize',16)
    CB = colorbar; clim(CLIM)
    % CB.Label.String = ['$u_\mathrm{radial}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT$_\mathrm{cg}$,HFR),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';
    title([num2str([2*PIXELS_list(pixel_i) + 1]*2) ' km fitting box'])
    
    AX4 = nexttile;
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(pixel_i))]) );
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
    % m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
    %                   eval(['minRMSD_cg_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    for ii = 1:length(Antenna_list) % unfiltered
        Antenna_ID = Antenna_list{ii};
        m_plot(real(eval([Antenna_ID '_latlon'])), ...
               imag(eval([Antenna_ID '_latlon'])),'o', ...
               'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
        m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
               imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
    end
    set(gca,'FontSize',16)
    CB = colorbar; clim(CLIM)
    % CB.Label.String = ['$u_\mathrm{radial, lowpass}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT$_\mathrm{cg}$,HFR$_\mathrm{LP}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
    CB.Label.Interpreter = 'latex';

    colormap(CMAP)
    set(gcf,'Position',[-776    82   739   895])
    linkaxes([AX1 AX2],'xy')
end

%% Better Radial Correlation Maps and RMSD Maps

close all

Antenna_list = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};
PIXELS_list = [2 4 6 8 10 12];
i_pixel = 4;

close all
figure('Color','w')
% tiledlayout(1,length(Antenna_list),'TileSpacing','tight')
tiledlayout(2,3,'TileSpacing','tight')


P_threshold = 0.05;

for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    % set(gcf,'color',[1 1 1]*0.7);
    hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrP']);
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.5 0.5],'k','LineWidth',1)
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.75 0.75],'k')
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',15)
    for jj = 1:length(Antenna_list)
        m_plot(real(eval([Antenna_list{jj} '_latlon'])),imag(eval([Antenna_list{jj} '_latlon'])),'ok','MarkerSize',5)
    end
    % title({Antenna_ID ; [num2str([2*PIXELS_list(i_pixel) + 1]*2) ' km fitting radius']})
    title({Antenna_ID})
    set(gca,'FontSize',16)
    clim([-1 1]*1)
end
CB = colorbar;
% CB.Label.String = ['CC, SWOT_g vs HFR_{lp} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.String = ['Radial C.C._{ug,hfr_{LP}} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

colormap(turbo)

set(gcf,'Position',[-922    83   923   894])

Axes_list = '';
for ii = 1:length(Antenna_list)
    Axes_list = [Axes_list ' AX' num2str(ii)];
end
eval(['linkaxes([' Axes_list '],''xy'')'])


% SWOTRadCC.SHEL.SHEL_8pix_weighted
%                  Antenna_ID: 'SHEL'
%                   KM_search: 3
%                      PIXELS: 8
%             Radial_cg_corrP: [285×69 double]
%             Radial_cg_corrR: [285×69 double]
%     Radial_cg_lowpass_corrP: [285×69 double]
%     Radial_cg_lowpass_corrR: [285×69 double]
%      Radial_cg_lowpass_rmsd: [285×69 double]
%              Radial_cg_rmsd: [285×69 double]
%              Radial_g_corrP: [285×69 double]
%              Radial_g_corrR: [285×69 double]
%      Radial_g_lowpass_corrP: [285×69 double]
%      Radial_g_lowpass_corrR: [285×69 double]
%       Radial_g_lowpass_rmsd: [285×69 double]
%               Radial_g_rmsd: [285×69 double]

figure('Color','w')
% tiledlayout(1,length(Antenna_list),'TileSpacing','tight')
tiledlayout(2,3,'TileSpacing','tight')
for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    % set(gcf,'color',[1 1 1]*0.7);
    hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrP']);
    CorrCoef_R =    eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P =    eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);
    RMSD_ant_swot = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      RMSD_ant_swot.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.5 0.5],'k','LineWidth',1)
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.75 0.75],'k')
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',15)
    for jj = 1:length(Antenna_list)
        m_plot(real(eval([Antenna_list{jj} '_latlon'])),imag(eval([Antenna_list{jj} '_latlon'])),'ok','MarkerSize',5)
    end
    % title({Antenna_ID ; [num2str([2*PIXELS_list(i_pixel) + 1]*2) ' km fitting radius']})
    title({Antenna_ID})
    set(gca,'FontSize',16)
    clim([0 1]*0.5)
end
CB = colorbar;
% CB.Label.String = ['CC, SWOT_g vs HFR_{lp} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.String = ['RMSD_{ug,hfr_{LP}} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

colormap(turbo)

set(gcf,'Position',[-922    83   923   894])

Axes_list = '';
for ii = 1:length(Antenna_list)
    Axes_list = [Axes_list ' AX' num2str(ii)];
end
eval(['linkaxes([' Axes_list '],''xy'')'])


%% 3-panel Better Radial Correlation Maps (formerly the publication figure)

disp('%%%%%%%%%%%%%%%%%%%%%')

close all

Antenna_list = {'SHEL','BRAG','PAFS'};
Antenna_list_full = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};
PIXELS_list = [2 4 6 8 10 12];
i_pixel = 4;
Panel_label = {'(a)','(b)','(c)'};

close all
figure('Color','w')
tiledlayout(1,length(Antenna_list),'TileSpacing','tight')
% tiledlayout(1,3,'TileSpacing','tight')

P_threshold = 0.05;

for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-125.6 -123.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    % set(gcf,'color',[1 1 1]*0.7);
    hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrP']);
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.5 0.5],'k','LineWidth',1)
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.75 0.75],'k')
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',25)
    m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, imag(eval([Antenna_ID '_latlon'])), Antenna_ID,'FontSize',16)
    m_text(-125.25, 41.5, Panel_label{ii},'Interpreter','latex','FontSize',30,'HorizontalAlignment','center')
    for jj = 1:length(Antenna_list_full)
        m_plot(real(eval([Antenna_list_full{jj} '_latlon'])), imag(eval([Antenna_list_full{jj} '_latlon'])),'ok','MarkerSize',8)
    end
    % title({Antenna_ID ; [num2str([2*PIXELS_list(i_pixel) + 1]*2) ' km fitting radius']})
    % title({Antenna_ID})
    set(gca,'FontSize',16)
    clim([-1 1]*1)

    disp('%%%%%%%%%%%%%%%%%%%%%')
    disp(prctile(Collimate(CorrCoef_R),[0 25 50 75 90 100]))
    disp(prctile(Collimate(CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]),[0 25 50 75 90 100]))
end
CB = colorbar;
% CB.Label.String = ['CC, SWOT_g vs HFR_{lp} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.Interpreter = 'latex';
CB.Label.String = ['Radial C.C. $u_g$ vs. $u_\mathrm{hfr,LP}$ (P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

colormap(turbo)

set(gcf,'Position',[-922    83   923   894])

Axes_list = '';
for ii = 1:length(Antenna_list)
    Axes_list = [Axes_list ' AX' num2str(ii)];
end
eval(['linkaxes([' Axes_list '],''xy'')'])
set(gcf,'Position',[-1073          83        1074         894])

% % % % % % % % % % %
% % % % % % % % % % %
% % % % % % % % % % %

figure('Color','w')
tiledlayout(1,length(Antenna_list),'TileSpacing','tight')

P_threshold = 0.05;

for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-125.6 -123.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    % set(gcf,'color',[1 1 1]*0.7);
    hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrP']);
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);
    RMSD_ant_swot = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      RMSD_ant_swot.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.5 0.5],'k','LineWidth',1)
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.75 0.75],'k')
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',25)
    m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, imag(eval([Antenna_ID '_latlon'])), Antenna_ID,'FontSize',16)
    m_text(-125.25, 41.5, Panel_label{ii},'Interpreter','latex','FontSize',30,'HorizontalAlignment','center')
    for jj = 1:length(Antenna_list_full)
        m_plot(real(eval([Antenna_list_full{jj} '_latlon'])), imag(eval([Antenna_list_full{jj} '_latlon'])),'ok','MarkerSize',8)
    end
    % title({Antenna_ID ; [num2str([2*PIXELS_list(i_pixel) + 1]*2) ' km fitting radius']})
    % title({Antenna_ID})
    set(gca,'FontSize',16)
    clim([0 1]*0.5)

    disp('%%%%%%%%%%%%%%%%%%%%%')
    disp(prctile(Collimate(CorrCoef_R),[0 25 50 75 90 100]))
    disp(prctile(Collimate(CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]),[0 25 50 75 90 100]))
end
CB = colorbar;
% CB.Label.String = ['CC, SWOT_g vs HFR_{lp} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.Interpreter = 'latex';
CB.Label.String = ['RMSD $u_g$ vs. $u_\mathrm{hfr,LP}$ (P$<$' num2str(P_threshold) ') (m s$^{-1}$)']; CB.Label.FontSize = 16;

colormap(turbo)

set(gcf,'Position',[-922    83   923   894])

Axes_list = '';
for ii = 1:length(Antenna_list)
    Axes_list = [Axes_list ' AX' num2str(ii)];
end
eval(['linkaxes([' Axes_list '],''xy'')'])
set(gcf,'Position',[-1073          83        1074         894])

% %%
% figure(1)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/radial_SWOT_corrcoef_3Antennas.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% 
% figure(2)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/radial_SWOT_rmsd_3Antennas.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% 6-panel Better Radial Correlation Maps, only good antennas (PUBLICATION FIGURE)

disp('%%%%%%%%%%%%%%%%%%%%%')

close all

Antenna_list = {'SHEL','BRAG','PAFS'};
Antenna_list_full = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};
PIXELS_list = [2 4 6 8 10 12];
i_pixel = 4;
% Panel_label = {'(a)','(b)','(c)'};
Panel_label = {'(a)','(b)','(c)','(d)','(e)','(f)'};

close all
figure('Color','w')
tiledlayout(2,length(Antenna_list),'TileSpacing','tight')

P_threshold = 0.05;

for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-125.6 -123.] + [-.1 .1],...
                 'latitudes',[37 41.25] + [-.1 .1]);
    % set(gcf,'color',[1 1 1]*0.7);
    hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrP']);
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.5 0.5],'k','LineWidth',1)
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.75 0.75],'k')
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',25)
    m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, imag(eval([Antenna_ID '_latlon'])), Antenna_ID,'FontSize',16)
    m_text(-125.25, 40.5, Panel_label{ii},'Interpreter','latex','FontSize',30,'HorizontalAlignment','center')
    for jj = 1:length(Antenna_list_full)
        m_plot(real(eval([Antenna_list_full{jj} '_latlon'])), imag(eval([Antenna_list_full{jj} '_latlon'])),'ok','MarkerSize',8)
    end
    % title({Antenna_ID ; [num2str([2*PIXELS_list(i_pixel) + 1]*2) ' km fitting radius']})
    % title({Antenna_ID})
    set(gca,'FontSize',16)
    clim([-1 1]*1)

    disp('%%%%%%%%%%%%%%%%%%%%%')
    disp(prctile(Collimate(CorrCoef_R),[0 25 50 75 90 100]))
    disp(prctile(Collimate(CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]),[0 25 50 75 90 100]))
end
II = ii;
CB = colorbar;
% CB.Label.String = ['CC, SWOT_g vs HFR_{lp} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.Interpreter = 'latex';
CB.Label.String = ['Radial C.C. $u_g^{\mathrm{rad}}$ vs. $u_\mathrm{hfr,LP}^{\mathrm{rad}}$ where P$<$' num2str(P_threshold)]; CB.Label.FontSize = 16;
colormap(turbo)

set(gcf,'Position',[-922    83   923   894])

Axes_list = '';
for ii = 1:length(Antenna_list)
    Axes_list = [Axes_list ' AX' num2str(ii)];
end

for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii+II) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-125.6 -123.] + [-.1 .1],...
                 'latitudes',[37 41.25] + [-.1 .1]);
    % set(gcf,'color',[1 1 1]*0.7);
    hold on

    % CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrR']);
    % CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrP']);
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);
    RMSD_ant_swot = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      RMSD_ant_swot.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.5 0.5],'k','LineWidth',1)
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R,[0.75 0.75],'k')
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',25)
    m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, imag(eval([Antenna_ID '_latlon'])), Antenna_ID,'FontSize',16)
    m_text(-125.25, 40.5, Panel_label{ii+II},'Interpreter','latex','FontSize',30,'HorizontalAlignment','center')
    for jj = 1:length(Antenna_list_full)
        m_plot(real(eval([Antenna_list_full{jj} '_latlon'])), imag(eval([Antenna_list_full{jj} '_latlon'])),'ok','MarkerSize',8)
    end
    % title({Antenna_ID ; [num2str([2*PIXELS_list(i_pixel) + 1]*2) ' km fitting radius']})
    % title({Antenna_ID})
    set(gca,'FontSize',16)
    clim([0 1]*0.5)

    disp('%%%%%%%%%%%%%%%%%%%%%')
    disp(prctile(Collimate(CorrCoef_R),[0 25 50 75 90 100]))
    disp(prctile(Collimate(CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold]),[0 25 50 75 90 100]))
end
CB = colorbar;
% CB.Label.String = ['CC, SWOT_g vs HFR_{lp} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.Interpreter = 'latex';
CB.Label.String = ['RMS($u_g^{\mathrm{rad}}$ - $u_\mathrm{hfr,LP}^{\mathrm{rad}}$) where P$<$' num2str(P_threshold) ' (m s$^{-1}$)']; CB.Label.FontSize = 16;

colormap(turbo)

for ii = 1:length(Antenna_list)
    Axes_list = [Axes_list ' AX' num2str(ii+II)];
end
eval(['linkaxes([' Axes_list '],''xy'')'])
set(gcf,'Position',[-801    90   802   882])

disp('%%%%%%%%%%%%%%%%%%%%%')
disp(['SHEL minimum RMSD:'])
disp(min(Collimate(eval(['SWOTRadCC.' 'SHEL' '.' 'SHEL' '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd']))))
disp(['BRAG minimum RMSD:'])
disp(min(Collimate(eval(['SWOTRadCC.' 'BRAG' '.' 'BRAG' '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd']))))
disp(['PAFS minimum RMSD:'])
disp(min(Collimate(eval(['SWOTRadCC.' 'PAFS' '.' 'PAFS' '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd']))))

disp('%%%%%%%%%%%%%%%%%%%%%')
disp(['SHEL [1 5 10]%ile RMSD:'])
disp(prctile(Collimate(eval(['SWOTRadCC.' 'SHEL' '.' 'SHEL' '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd'])),[1 5 10]))
disp(['BRAG [1 5 10]%ile RMSD:'])
disp(prctile(Collimate(eval(['SWOTRadCC.' 'BRAG' '.' 'BRAG' '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd'])),[1 5 10]))
disp(['PAFS [1 5 10]%ile RMSD:'])
disp(prctile(Collimate(eval(['SWOTRadCC.' 'PAFS' '.' 'PAFS' '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd'])),[1 5 10]))

disp('%%%%%%%%%%%%%%%%%%%%%')
for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);
    RMSD_ant_swot = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd']);
    disp([Antenna_ID ' minimum RMSD:'])
end


% %%
% figure(1)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/radial_SWOT_corrcoef_RMSD_3Antennas.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% Hist2 CorrCoef and RMSD

close all

figure('Color','w')
tiledlayout(1,length(Antenna_list),'TileSpacing','tight')
for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii) ' = nexttile;'])
    
    hold on
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);
    RMSD_ant_swot = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_rmsd']);

    histogram2(CorrCoef_R   ,....*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold],...
               RMSD_ant_swot,....*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold],...
               'DisplayStyle','tile','ShowEmptyBins','off','Normalization','pdf');
    xlabel('CC'); ylabel('RMSD')
    title(Antenna_ID)
end

figure('Color','w')
for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    
    hold on
    CorrCoef_R = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);

    histogram(CorrCoef_R.*[CorrCoef_P<P_threshold]./[CorrCoef_P<P_threshold],...
              'DisplayStyle','stairs','Normalization','pdf','LineWidth',2);
    xlabel('CC'); ylabel('RMSD')
    title(Antenna_ID)
end

%% Difference between CG and G correlation

close all

Antenna_list = {'SHEL','BRAG','PAFS'};
Antenna_list_full = {'TRIN','SMOA','SHEL','BRAG','PAFS','GCVE'};
PIXELS_list = [2 4 6 8 10 12];
i_pixel = 4;
Panel_label = {'(a)','(b)','(c)'};

close all
figure('Color','w')
tiledlayout(1,length(Antenna_list),'TileSpacing','tight')

P_threshold = 0.05;

for ii = 1:length(Antenna_list)
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
    hold on

    % CorrCoef_R_g  = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrR']);
    % CorrCoef_P_g  = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_g_lowpass_corrP']);
    % CorrCoef_R_cg = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_cg_lowpass_corrR']);
    % CorrCoef_P_cg = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix.Radial_cg_lowpass_corrP']);
    CorrCoef_R_g  = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrR']);
    CorrCoef_P_g  = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_g_lowpass_corrP']);
    CorrCoef_R_cg = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_cg_lowpass_corrR']);
    CorrCoef_P_cg = eval(['SWOTRadCC.' Antenna_ID '.' Antenna_ID '_' num2str(PIXELS_list(i_pixel)) 'pix_weighted.Radial_cg_lowpass_corrP']);

    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      [CorrCoef_R_g.*[CorrCoef_P_g<P_threshold]./[CorrCoef_P_g<P_threshold]] - ...
                      [CorrCoef_R_cg.*[CorrCoef_P_cg<P_threshold]./[CorrCoef_P_cg<P_threshold]]);
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R_g,[0.5 0.5],'k','LineWidth',1)
    m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
              CorrCoef_R_g,[0.75 0.75],'k')
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    m_plot(real(eval([Antenna_ID '_latlon'])),imag(eval([Antenna_ID '_latlon'])),'.k','MarkerSize',25)
    m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, imag(eval([Antenna_ID '_latlon'])), Antenna_ID,'FontSize',16)
    m_text(-122.5, 41.5, Panel_label{ii},'Interpreter','latex','FontSize',24,'HorizontalAlignment','center')
    for jj = 1:length(Antenna_list_full)
        m_plot(real(eval([Antenna_list_full{jj} '_latlon'])), imag(eval([Antenna_list_full{jj} '_latlon'])),'ok','MarkerSize',8)
    end
    set(gca,'FontSize',16)
    clim([-1 1]*0.1)

    disp('%%%%%%%%%%%%%%%%%%%%%')
    disp(prctile(Collimate(CorrCoef_R_g),[0 25 50 75 90 100]))
    disp(prctile(Collimate(CorrCoef_R_g.*[CorrCoef_P_g<P_threshold]./[CorrCoef_P_g<P_threshold]),[0 25 50 75 90 100]))
end
CB = colorbar;
% CB.Label.String = ['CC, SWOT_g vs HFR_{lp} (P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.Interpreter = 'latex';
CB.Label.String = {['$\Delta$ Radial C.C. $u_g$ vs. $u_\mathrm{hfr,LP}$ (P$<$' num2str(P_threshold) ')'] ; ...
                    'C(g,hfr) - C(cg,hfr)'};
CB.Label.FontSize = 16;

colormap(turbo)

set(gcf,'Position',[-922    83   923   894])

Axes_list = '';
for ii = 1:length(Antenna_list)
    Axes_list = [Axes_list ' AX' num2str(ii)];
end
eval(['linkaxes([' Axes_list '],''xy'')'])
set(gcf,'Position',[-1073          83        1074         894])

%% Time series comparison

close all
figure('Color','w')
tiledlayout(2,4,'TileSpacing','tight')

nexttile([2 1])
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
                 'longitudes',[-126 -122.] + [-.1 .1],...
                 'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(pixel_i))]) );
m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.5 0.5],'k','LineWidth',1);
m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  eval(['maxCorrCoef_cg_lowpass.pix' num2str(PIXELS_list(pixel_i))]) ,[0.75 0.75],'k');
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
for ii = 1:length(Antenna_list) % unfiltered
    Antenna_ID = Antenna_list{ii};
    m_plot(real(eval([Antenna_ID '_latlon'])), ...
           imag(eval([Antenna_ID '_latlon'])),'o', ...
           'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
    m_text(real(eval([Antenna_ID '_latlon'])) + 0.1, ...
           imag(eval([Antenna_ID '_latlon'])), Antenna_ID)
end
m_text(-122.75, 41.75, '(b)','Interpreter','latex','FontSize',32)
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
CB.Label.Interpreter = 'latex';
title(['CorrCoef(u$^\mathrm{rad}_\mathrm{cg}$,u$^\mathrm{rad}_\mathrm{HFR, low pass}$),  P$<$' num2str(P_threshold) ''],'Interpreter','latex')
colormap(CMAP)

[x_click,y_click] = m_ginput(1);
[yi_click,xi_click] = find(abs([NORCAL.SWOT.lon{1}    + 1i*NORCAL.SWOT.lat{1}]    - [x_click + 1i*y_click]) == ...
                       min(abs([NORCAL.SWOT.lon{1}(:) + 1i*NORCAL.SWOT.lat{1}(:)] - [x_click + 1i*y_click])));

m_plot(x_click,y_click,'+w')

% nexttile([1 3])
% plot()


%%
%%
%%
%% $
%% Difference in CorrCoef after time-filtering HFR

% % % % Full HFR
% maxCorrCoef = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(antenna_IDs));
% maxCorrCoef(:,:,1) = TRIN.Radial_vel_corrR;
% maxCorrCoef(:,:,2) = GCVE.Radial_vel_corrR;
% maxCorrCoef(:,:,3) = SHEL.Radial_vel_corrR;
% maxCorrCoef(:,:,4) = BRAG.Radial_vel_corrR;
% maxCorrCoef(:,:,5) = PAFS.Radial_vel_corrR;
% maxCorrCoef = max(maxCorrCoef,[],3,'omitnan');
% 
% % % % Filtered HFR
% maxCorrCoef_lowpass = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(antenna_IDs));
% maxCorrCoef_lowpass(:,:,1) = TRIN.Radial_vel_lowpass_corrR;
% maxCorrCoef_lowpass(:,:,2) = GCVE.Radial_vel_lowpass_corrR;
% maxCorrCoef_lowpass(:,:,3) = SHEL.Radial_vel_lowpass_corrR;
% maxCorrCoef_lowpass(:,:,4) = BRAG.Radial_vel_lowpass_corrR;
% maxCorrCoef_lowpass(:,:,5) = PAFS.Radial_vel_lowpass_corrR;
% maxCorrCoef_lowpass = max(maxCorrCoef_lowpass,[],3,'omitnan');

close all
figure('Color','w')
P_threshold = 0.05;

M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  maxCorrCoef_g_lowpass - maxCorrCoef_g);
% m_scatter(real(unique(Radial.LON + 1i*Radial.LAT)),...
%           imag(unique(Radial.LON + 1i*Radial.LAT)),1,'k','filled')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
for ii = 1:length(Antenna_ID) % unfiltered
    Antenna_ID = Antenna_list{ii};
    m_plot(real(eval([Antenna_ID '.Origin_antenna'])), ...
           imag(eval([Antenna_ID '.Origin_antenna'])), 'o', ...
           'MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    m_text(real(eval([Antenna_ID '.Origin_antenna'])) + 0.1, ...
           imag(eval([Antenna_ID '.Origin_antenna'])), Antenna_ID)
end
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1]*.5)
% CB.Label.String = ['CorrCoef($u_\mathrm{radial, lowpass}$) - CorrCoef($u_\mathrm{radial}$) (P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT,HFR) - CorrCoef$_\mathrm{radial}$(SWOT,HFR$_\mathrm{LP}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
CB.Label.Interpreter = 'latex';
colormap([1/5]*round(5*bwr(100)))
% colormap(turbo)

set(gcf,'Position',[-704   238   667   655])
%%
figure('Color','w')
histogram(maxCorrCoef_g_lowpass - maxCorrCoef_g,[-1:0.01:1],'Normalization','pdf'); hold on
text(-1,1,['Mean = ' num2str(mean(maxCorrCoef_g_lowpass(:) - maxCorrCoef_g(:),'omitnan'))])
text(-1,0.9,['Median = ' num2str(median(maxCorrCoef_g_lowpass(:) - maxCorrCoef_g(:),'omitnan'))])
text(-1,0.8,['STDev = ' num2str(std(maxCorrCoef_g_lowpass(:) - maxCorrCoef_g(:),'omitnan'))])


% %%
% figure(1)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/radial_SWOT_corrcoef_diff.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% All antennas

figure('Color','w')
tiledlayout(2,length(Antenna_ID),'TileSpacing','tight')

for ii = 1:length(Antenna_ID) % unfiltered
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[-126 -122.] + [-.1 .1],...
        'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval([Antenna_ID '.Radial_vel_corrR .* [' ...
                            Antenna_ID '.Radial_vel_corrP<P_threshold]./[' ...
                            Antenna_ID '.Radial_vel_corrP<P_threshold]']));
    m_plot(real(eval([Antenna_ID '.Origin_antenna'])), imag(eval([Antenna_ID '.Origin_antenna'])), 'k.','MarkerSize',25)
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1])
    title(['BML: ' Antenna_ID])
end
CB.Label.String = ['Vel_{radial}: corr. coef. (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
for ii = 1:length(Antenna_ID) % filtered
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(length(Antenna_ID) + ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[-126 -122.] + [-.1 .1],...
        'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval([Antenna_ID '.Radial_vel_lowpass_corrR .* [' ...
                            Antenna_ID '.Radial_vel_lowpass_corrP<P_threshold]./[' ...
                            Antenna_ID '.Radial_vel_lowpass_corrP<P_threshold]']));
    m_plot(real(eval([Antenna_ID '.Origin_antenna'])), imag(eval([Antenna_ID '.Origin_antenna'])), 'k.','MarkerSize',25)
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1])
    title(['BML: ' Antenna_ID])
end
CB.Label.String = ['Vel_{radial, lowpass}: corr. coef. (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
colormap(CMAP)
linkaxes([AX1 AX2 AX3 ...
          AX4 AX5 AX6],'xy')
set(gcf,'Position',[-1439         271        1205         540])
% %%
% figure(2)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/radial_SWOT_corrcoef_8panel.pdf',...
% 'BackgroundColor','none','ContentType','vector')

% % % % % % % % % % % 


figure('Color','w')
% tiledlayout(2,2,'TileSpacing','tight')
tiledlayout(2,4,'TileSpacing','tight')

% nexttile([2 1]);
% M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
%              'longitudes',[-126 -122.] + [-.1 .1],...
%              'latitudes',[37 42] + [-.1 .1]);
% set(gcf,'color','w'); hold on
% m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
%                   maxCorrCoef_lowpass - maxCorrCoef);
% COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
% m_grid('box','fancy', 'backgroundcolor','none');
% set(gca,'FontSize',16)
% CB = colorbar; clim([-1 1])
% CB.Label.String = ['CorrCoef_{lowpass} - CorrCoef_{full}']; CB.Label.FontSize = 16;
% colormap(bwr)

for ii = 1:length(Antenna_ID)
    Antenna_ID = Antenna_list{ii};
    eval(['AX' num2str(ii) ' = nexttile([2,1]);'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[-126 -122.] + [-.1 .1],...
        'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval([Antenna_ID '.Radial_vel_lowpass_corrR' ...
                    ' .* [' Antenna_ID '.Radial_vel_lowpass_corrP<P_threshold]./[' ... dot comment for low-P areas
                            Antenna_ID '.Radial_vel_lowpass_corrP<P_threshold]' ... dot comment for low-P areas
                            ' - ' ...
                            Antenna_ID '.Radial_vel_corrR' ...
                    ' .* [' Antenna_ID '.Radial_vel_corrP<P_threshold]./[' ... dot comment for low-P areas
                            Antenna_ID '.Radial_vel_corrP<P_threshold]' ... dot comment for low-P areas
                            ]));
    m_plot(real(eval([Antenna_ID '.Origin_antenna'])), imag(eval([Antenna_ID '.Origin_antenna'])), 'k.','MarkerSize',25)
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1]*0.5)
    CB.Label.String = ['U_{radial, lowpass}: corr. coef. (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    title(['BML: ' Antenna_ID])
end
colormap(bwr)
linkaxes([AX1 AX2 AX3],'xy')

nexttile;
histogram(maxCorrCoef_g_lowpass - maxCorrCoef_g,[-1:0.01:1],'Normalization','pdf')
xlabel('CorrCoef_{lowpass} - CorrCoef_{full}')
xlabel('PDF')
title({['Mean = ' num2str(mean(maxCorrCoef_g_lowpass(:) - maxCorrCoef_g(:),'omitnan'))] ; ...
       ['Median = ' num2str(median(maxCorrCoef_g_lowpass(:) - maxCorrCoef_g(:),'omitnan'))] ; ...
       ['\sigma = ' num2str(std(maxCorrCoef_g_lowpass(:) - maxCorrCoef_g(:),'omitnan'))]})

nexttile;
histogram([SHEL.Radial_vel_lowpass_corrR .* [SHEL.Radial_vel_lowpass_corrP<P_threshold]./[SHEL.Radial_vel_lowpass_corrP<P_threshold]] - ...
          [SHEL.Radial_vel_corrR .* [SHEL.Radial_vel_corrP<P_threshold]./[SHEL.Radial_vel_corrP<P_threshold]], ...
          [-1:0.01:1],'Normalization','pdf'); hold on
histogram([BRAG.Radial_vel_lowpass_corrR .* [BRAG.Radial_vel_lowpass_corrP<P_threshold]./[BRAG.Radial_vel_lowpass_corrP<P_threshold]] - ...
          [BRAG.Radial_vel_corrR .* [BRAG.Radial_vel_corrP<P_threshold]./[BRAG.Radial_vel_corrP<P_threshold]], ...
          [-1:0.01:1],'Normalization','pdf');
histogram([PAFS.Radial_vel_lowpass_corrR .* [PAFS.Radial_vel_lowpass_corrP<P_threshold]./[PAFS.Radial_vel_lowpass_corrP<P_threshold]] - ...
          [PAFS.Radial_vel_corrR .* [PAFS.Radial_vel_corrP<P_threshold]./[PAFS.Radial_vel_corrP<P_threshold]], ...
          [-1:0.01:1],'Normalization','pdf');
xlabel('CorrCoef_{lowpass} - CorrCoef_{full}')
xlabel('PDF')

SHEL.diff_mean = mean([SHEL.Radial_vel_lowpass_corrR(:) .* [SHEL.Radial_vel_lowpass_corrP(:)<P_threshold]./[SHEL.Radial_vel_lowpass_corrP(:)<P_threshold]] - ...
          [SHEL.Radial_vel_corrR(:) .* [SHEL.Radial_vel_corrP(:)<P_threshold]./[SHEL.Radial_vel_corrP(:)<P_threshold]],'omitnan');
SHEL.diff_std = std([SHEL.Radial_vel_lowpass_corrR(:) .* [SHEL.Radial_vel_lowpass_corrP(:)<P_threshold]./[SHEL.Radial_vel_lowpass_corrP(:)<P_threshold]] - ...
          [SHEL.Radial_vel_corrR(:) .* [SHEL.Radial_vel_corrP(:)<P_threshold]./[SHEL.Radial_vel_corrP(:)<P_threshold]],'omitnan');
BRAG.diff_mean = mean([BRAG.Radial_vel_lowpass_corrR(:) .* [BRAG.Radial_vel_lowpass_corrP(:)<P_threshold]./[BRAG.Radial_vel_lowpass_corrP(:)<P_threshold]] - ...
          [BRAG.Radial_vel_corrR(:) .* [BRAG.Radial_vel_corrP(:)<P_threshold]./[BRAG.Radial_vel_corrP(:)<P_threshold]],'omitnan');
BRAG.diff_std = std([BRAG.Radial_vel_lowpass_corrR(:) .* [BRAG.Radial_vel_lowpass_corrP(:)<P_threshold]./[BRAG.Radial_vel_lowpass_corrP(:)<P_threshold]] - ...
          [BRAG.Radial_vel_corrR(:) .* [BRAG.Radial_vel_corrP(:)<P_threshold]./[BRAG.Radial_vel_corrP(:)<P_threshold]],'omitnan');
PAFS.diff_mean = mean([PAFS.Radial_vel_lowpass_corrR(:) .* [PAFS.Radial_vel_lowpass_corrP(:)<P_threshold]./[PAFS.Radial_vel_lowpass_corrP(:)<P_threshold]] - ...
          [PAFS.Radial_vel_corrR(:) .* [PAFS.Radial_vel_corrP(:)<P_threshold]./[PAFS.Radial_vel_corrP(:)<P_threshold]],'omitnan');
PAFS.diff_std = std([PAFS.Radial_vel_lowpass_corrR(:) .* [PAFS.Radial_vel_lowpass_corrP(:)<P_threshold]./[PAFS.Radial_vel_lowpass_corrP(:)<P_threshold]] - ...
          [PAFS.Radial_vel_corrR(:) .* [PAFS.Radial_vel_corrP(:)<P_threshold]./[PAFS.Radial_vel_corrP(:)<P_threshold]],'omitnan');
legend(['SHEL - (\mu,\sigma) = (' num2str(SHEL.diff_mean) ', ' num2str(SHEL.diff_std) ')'],...
       ['BRAG - (\mu,\sigma) = (' num2str(BRAG.diff_mean) ', ' num2str(BRAG.diff_std) ')'],...
       ['PAFS - (\mu,\sigma) = (' num2str(PAFS.diff_mean) ', ' num2str(PAFS.diff_std) ')'])

% % % % % % % % % % % 

figure('Color','w')
tiledlayout(1,3,'TileSpacing','tight')

% histogram([SHEL.Radial_vel_lowpass_corrR .* [SHEL.Radial_vel_lowpass_corrP<P_threshold]./[SHEL.Radial_vel_lowpass_corrP<P_threshold]], ...
%           [-1:0.01:1],'Normalization','pdf'); hold on
% histogram([SHEL.Radial_vel_corrR .* [SHEL.Radial_vel_corrP<P_threshold]./[SHEL.Radial_vel_corrP<P_threshold]], ...
%           [-1:0.01:1],'Normalization','pdf');

nexttile
histogram(SHEL.Radial_vel_lowpass_corrR, [-1:0.05:1],'Normalization','pdf'); hold on
histogram(SHEL.Radial_vel_corrR,         [-1:0.05:1],'Normalization','pdf');
title(['BML: SHEL'])

nexttile
histogram(BRAG.Radial_vel_lowpass_corrR, [-1:0.05:1],'Normalization','pdf'); hold on
histogram(BRAG.Radial_vel_corrR,         [-1:0.05:1],'Normalization','pdf');
title(['BML: BRAG'])

nexttile
histogram(PAFS.Radial_vel_lowpass_corrR, [-1:0.05:1],'Normalization','pdf'); hold on
histogram(PAFS.Radial_vel_corrR,         [-1:0.05:1],'Normalization','pdf');
title(['BML: PAFS'])
legend('correlation coefficient (lowpass)','correlation coefficient (full)')


xlabel('CorrCoef_{lowpass} - CorrCoef_{full}')
xlabel('PDF')


%%
%%
%% Frequency Spectra evolving over time at some point (similar to wavelet)
%% This does not provide much information; the eddy may not persist enough
%  for a clear signal of its vorticity modulating NIOs (or no strong NIOs
%  occurred at this time)
%%
%%
%%

User_chosen_date = datenum('2023-05-09 02:00:00');
ti_swot = dsearchn(NORCAL.SWOT.mean_time, Radial.TIME(dsearchn(Radial.TIME,User_chosen_date)));
SSHA = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
% SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');

% SSHA snapshot:
close all
figure('Color','w')
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(Radial.LON) max(Radial.LON)] + [-1 1]/1,...
             'latitudes', [min(Radial.LAT) max(Radial.LAT)] + [-1 1]/1);
set(gcf,'color','w'); hold on
VecScale = 0.001;
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                  NORCAL.SWOT.lat{ti_swot},...
                  SSHA);
m_plot(real(Radial.Origin_antenna), imag(Radial.Origin_antenna), 'ko')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo)
CB = colorbar;
CB.Label.String = 'SWOT: \eta (m)'; CB.Label.FontSize = 20; clim([-1 1]*0.15)
set(gca,'FontSize',16)


[x_click,y_click] = m_ginput(1);
[yi_click,xi_click] = find(abs([NORCAL.SWOT.lon{ti_swot}    + 1i*NORCAL.SWOT.lat{ti_swot}]    - [x_click + 1i*y_click]) == ...
                       min(abs([NORCAL.SWOT.lon{ti_swot}(:) + 1i*NORCAL.SWOT.lat{ti_swot}(:)] - [x_click + 1i*y_click])));


m_plot(NORCAL.SWOT.lon{ti_swot}(yi_click,xi_click),...
       NORCAL.SWOT.lat{ti_swot}(yi_click,xi_click),'*w')
%%
KM_search = 0.4;
KM_search = 2.5;
KM_search = 3.0;
% KM_search = 6.0;
[Distance_from_click,~,~] = m_idist(x_click,y_click,   Radial.LON,Radial.LAT);
IND = [Distance_from_click/1000 < KM_search];

T_HFR_at_point = unique(Radial.TIME(IND));
Ur_HFR_at_point = nan(length(T_HFR_at_point),1);
UrSTD_HFR_at_point = nan(length(T_HFR_at_point),1);
tic
for ii = 1:length(T_HFR_at_point)
    Ur_HFR_at_point(ii)    = mean(Radial.VEL(IND & Radial.TIME==T_HFR_at_point(ii))/100,"omitmissing");
    UrSTD_HFR_at_point(ii) = std( Radial.VEL(IND & Radial.TIME==T_HFR_at_point(ii))/100,"omitmissing");
    % T_HFR_at_point(ii)     = mean(DATA(IND & Radial.TIME==unique_Time_HFR(ii),21));
    % Ur_HFR_at_point(ii)    = mean(DATA(IND & Radial.TIME==unique_Time_HFR(ii),18)/100,"omitmissing");
    % UrSTD_HFR_at_point(ii) = std( DATA(IND & Radial.TIME==unique_Time_HFR(ii),18)/100,"omitmissing");
end
toc

%% nunanspec

Spectrum_time_evolving = [];
T_center = [];
N_hours = 14*24; IND = 1:N_hours;
TIME_UNITS = 'day';
warning('requires "nunanspectrum" function')
II = 0;
% Initialize frequency and perform first spectral estimate
    [Spec_jj,Freq,~] = ...
        nunanspectrum(Ur_HFR_at_point(II + IND),T_HFR_at_point(II + IND),'day',...
        'Segments',1, 'Plot',false, 'Window','rectwin');
    Spectrum_time_evolving = [Spectrum_time_evolving, Spec_jj];
    T_center = [T_center, mean(T_HFR_at_point(II + IND),'omitnan')];
    II = II + N_hours;
while [II + N_hours] < length(T_HFR_at_point)
    [Spec_jj,~,~] = ...
        nunanspectrum(Ur_HFR_at_point(II + IND),T_HFR_at_point(II + IND),TIME_UNITS,...
        'Segments',1, 'Plot',false, 'Window','rectwin', 'Freq',Freq);
    Spectrum_time_evolving = [Spectrum_time_evolving, Spec_jj];
    T_center = [T_center, mean(T_HFR_at_point(II + IND),'omitnan')];
    II = II + N_hours/3;
end

[TT,FF] = meshgrid(T_center,Freq);

figure('Color','w')
pcolor_centered(TT,FF,log10(Spectrum_time_evolving)); shading flat
% loglog(Freq,Spectrum_time_evolving); shading flat
datetick('x')
ylabel([TIME_UNITS '^{-1}'])


%%
%%
%%
%%
%%
%% Auxiliary functions
%%
% Interpolate as average withing some radius:
function D_out = mean_in_new_grid(X_in,Y_in,D_in,X_out,Y_out,ThreshholdRadius)
D_out = nan(size(X_out));
for ii=1:size(X_out,1) % make sure the grids are how you expect them, i.e.
                       % imagesc(X_out) looks like it should from space
    for jj=1:size(X_out,2)
        IND = sqrt([X_in - X_out(ii,jj)].^2 + [Y_in - Y_out(ii,jj)].^2) < ThreshholdRadius;
        D_out(ii,jj) = mean((D_in(IND)),'omitnan');
    end
end
end
%%
% Interpolate over the nadir gap before smoothing (only works for SWOT
% passes where N-S is the vertical direction):
function BB = smoothdata2_interpnadirgap(XX,YY,AA,DELTA,SIGMA)
    % Error will occur if any of the inputs are different sizes:
    XX.*YY.*AA;
    for ii = 1:size(AA,1)
        AA(ii,:) = interp_over_nan(XX(ii,:),AA(ii,:),'linear');
    end
    % BB = smoothdata2(AA,'gaussian',nn,'omitnan');
    BB = bettergaussiansmooth2(AA,DELTA,SIGMA);
end
%%
% Coriolis frequency in s^-1, given input of degrees latitude:
function FCOR = fcor_degrees_cps(LAT)
    % sideral day = 23.9345 hours = 86164.2 seconds
    FCOR = 2*[2*pi/(86164.2)]*sind(LAT);
end
%% Improved smoothing function:
% This function was not written by any of the authors (see attribution
% below):
function c = nanconv(a, k, varargin)
% NANCONV Convolution in 1D or 2D ignoring NaNs.
%   C = NANCONV(A, K) convolves A and K, correcting for any NaN values
%   in the input vector A. The result is the same size as A (as though you
%   called 'conv' or 'conv2' with the 'same' shape).
%
%   C = NANCONV(A, K, 'param1', 'param2', ...) specifies one or more of the following:
%     'edge'     - Apply edge correction to the output.
%     'noedge'   - Do not apply edge correction to the output (default).
%     'nanout'   - The result C should have NaNs in the same places as A.
%     'nonanout' - The result C should have ignored NaNs removed (default).
%                  Even with this option, C will have NaN values where the
%                  number of consecutive NaNs is too large to ignore.
%     '2d'       - Treat the input vectors as 2D matrices (default).
%     '1d'       - Treat the input vectors as 1D vectors.
%                  This option only matters if 'a' or 'k' is a row vector,
%                  and the other is a column vector. Otherwise, this
%                  option has no effect.
%
%   NANCONV works by running 'conv2' either two or three times. The first
%   time is run on the original input signals A and K, except all the
%   NaN values in A are replaced with zeros. The 'same' input argument is
%   used so the output is the same size as A. The second convolution is
%   done between a matrix the same size as A, except with zeros wherever
%   there is a NaN value in A, and ones everywhere else. The output from
%   the first convolution is normalized by the output from the second 
%   convolution. This corrects for missing (NaN) values in A, but it has
%   the side effect of correcting for edge effects due to the assumption of
%   zero padding during convolution. When the optional 'noedge' parameter
%   is included, the convolution is run a third time, this time on a matrix
%   of all ones the same size as A. The output from this third convolution
%   is used to restore the edge effects. The 'noedge' parameter is enabled
%   by default so that the output from 'nanconv' is identical to the output
%   from 'conv2' when the input argument A has no NaN values.
%
% See also conv, conv2
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2013, Benjamin Kraus
% $Id: nanconv.m 4861 2013-05-27 03:16:22Z bkraus $
% Process input arguments
for arg = 1:nargin-2
    switch lower(varargin{arg})
        case 'edge'; edge = true; % Apply edge correction
        case 'noedge'; edge = false; % Do not apply edge correction
        case {'same','full','valid'}; shape = varargin{arg}; % Specify shape
        case 'nanout'; nanout = true; % Include original NaNs in the output.
        case 'nonanout'; nanout = false; % Do not include NaNs in the output.
        case {'2d','is2d'}; is1D = false; % Treat the input as 2D
        case {'1d','is1d'}; is1D = true; % Treat the input as 1D
    end
end
% Apply default options when necessary.
if(exist('edge','var')~=1); edge = false; end
if(exist('nanout','var')~=1); nanout = false; end
if(exist('is1D','var')~=1); is1D = false; end
if(exist('shape','var')~=1); shape = 'same';
elseif(~strcmp(shape,'same'))
    error([mfilename ':NotImplemented'],'Shape ''%s'' not implemented',shape);
end
% Get the size of 'a' for use later.
sza = size(a);
% If 1D, then convert them both to columns.
% This modification only matters if 'a' or 'k' is a row vector, and the
% other is a column vector. Otherwise, this argument has no effect.
if(is1D);
    if(~isvector(a) || ~isvector(k))
        error('MATLAB:conv:AorBNotVector','A and B must be vectors.');
    end
    a = a(:); k = k(:);
end
% Flat function for comparison.
o = ones(size(a));
% Flat function with NaNs for comparison.
on = ones(size(a));
% Find all the NaNs in the input.
n = isnan(a);
% Replace NaNs with zero, both in 'a' and 'on'.
a(n) = 0;
on(n) = 0;
% Check that the filter does not have NaNs.
if(any(isnan(k)));
    error([mfilename ':NaNinFilter'],'Filter (k) contains NaN values.');
end
% Calculate what a 'flat' function looks like after convolution.
if(any(n(:)) || edge)
    flat = conv2(on,k,shape);
else flat = o;
end
% The line above will automatically include a correction for edge effects,
% so remove that correction if the user does not want it.
if(any(n(:)) && ~edge); flat = flat./conv2(o,k,shape); end
% Do the actual convolution
c = conv2(a,k,shape)./flat;
% If requested, replace output values with NaNs corresponding to input.
if(nanout); c(n) = NaN; end
% If 1D, convert back to the original shape.
if(is1D && sza(1) == 1); c = c.'; end
end

% This one was written by LK:
function D_smoothed = bettergaussiansmooth2(DD,DELTA,SIGMA)
% DD = input matrix for smoothing (must be on a grid with equal-sized
%   spacing, i.e. dx = dy = DELTA)
% Delta = length of a grid step of DELTA in physical units (e.g. km)
% SIGMA = standard deviation of the desired Gaussian smoothing
%   window/kernel (in the same physical units as DELTA)
NN = round(4*SIGMA/DELTA); % pixels along the edge of the square window
[XX,YY] = meshgrid(DELTA*[(0:[NN-1]) - NN/2], ...
                   DELTA*[(0:[NN-1]) - NN/2]); % window grid coordinates
WW = exp(-[XX.^2 + YY.^2]/[SIGMA^2]);
% WW(WW < exp(-[(NN*DELTA)^2]/[(2*SIGMA)^2])) = 0;
WW(WW < exp(-4)) = 0; % make isotropic
WW = WW/[sum(WW(:))];
D_smoothed = nanconv(DD,WW,'same','edge');
end

%% 6-term 2D polynomial fit
 % 6-term 2D polynomial fit
 % 
 % IN:
 % xx = 2D x-grid (NxM)
 % yy = 2D y-grid
 % dd = 2D data grid
 % nn = number in either direction from a point to fit (e.g. nn = 2 means a
 %      5x5 grid)
 % wq = optional input, can be anything; if given, this enables a weighted
 %      least squares fit, where the center is weighted as 1 and increasing
 %      distance results in 1/r^2 dropoff.
 % 
 % OUT:
 % AA = coefficients corresponding to the terms, at each point (NxMx6):
 %      [1, x, y, x^2, y^2, xy]
 % dd_fit = 2D smoothed data grid
 % 
 % To get geostrophic velocity from inputting surface height:
 %      u_g = -(g/f) * squeeze(AA(:,:,3))
 %      v_g = +(g/f) * squeeze(AA(:,:,2))
 % 
 % See:
 % Yann-Treden Tranchant, Benoit Legresy, Annie Foppert, et al. SWOT reveals
 % fine-scale balanced motions and dispersion properties in the Antarctic
 % Circumpolar Current. ESS Open Archive . January 11, 2025.
 % DOI: 10.22541/essoar.173655552.25945463/v1

function [AA,dd_fit] = fit_6term_2d(xx,yy,dd,nn,varargin)
xx.*yy.*dd; % size check
Collimate = @(IN) IN(:);
dd_fit = nan(size(dd));
AA = nan(size(dd,1), size(dd,2), 6);
II = size(xx,2);
JJ = size(xx,1);

if nargin < 5
    isWeighted = false;
    WW = speye(length(xx));
else
    isWeighted = true;
end

dxy = sqrt(median(diff(xx(1,:)),'omitnan').^2 + ...
           median(diff(yy(:,1)),'omitnan').^2);

for ii = [nn+1]:[II - nn]
    for jj = [nn+1]:[JJ - nn]

        % % % % % Centered at each point:
        HH = [Collimate(ones(size(xx([jj-nn]:[jj+nn],[ii-nn]:[ii+nn])))) ... 
              Collimate( xx([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]) - xx(jj,ii) ) ...
              Collimate( yy([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]) - yy(jj,ii) ) ...
              Collimate([xx([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]) - xx(jj,ii)].^2) ...
              Collimate([yy([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]) - yy(jj,ii)].^2) ...
              Collimate([xx([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]) - xx(jj,ii)] .* ...
                        [yy([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]) - yy(jj,ii)] )];
        HH_ = [1 0 0 0 0 0];



        if isWeighted
            % disp('IS WEIGHTED')
            XXYY2 = HH(:,4) + HH(:,5);
            DIAG = exp(XXYY2/[([2*nn + 1]*dxy/6).^2]); % Gaussian, +-3Sigma at edges
            RR = speye(length(DIAG));
            RR = spdiags(DIAG,0,RR);
            WW = inv(RR);

            dd_ij = Collimate(dd([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]));
            if sum(isfinite(dd_ij)) == 0
                dd_fit(jj,ii) = nan;
                AA(jj,ii,:) = nan;
            else
                aa = inv(HH(isfinite(dd_ij),:)'*WW(isfinite(dd_ij),isfinite(dd_ij))*HH(isfinite(dd_ij),:))*...
                     HH(isfinite(dd_ij),:)'*WW(isfinite(dd_ij),isfinite(dd_ij))*dd_ij(isfinite(dd_ij));
                dd_fit(jj,ii) = HH_*aa;
                AA(jj,ii,:) = aa;
            end
            %$# figure; plot(HH(:,2),'.-');hold on;plot(HH(:,3),'.-');error
            %$# figure; plot(diag(WW),'.-');error
            %$# figure; imagesc(exp(-[[xx([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]) - xx(jj,ii)].^2 + [yy([jj-nn]:[jj+nn],[ii-nn]:[ii+nn]) - yy(jj,ii)].^2]/([2*nn + 1]*dxy/6).^2));error
            % aa = HH\dd_ij;
            % dd_fit(jj,ii) = HH_*aa;
        else
            % disp('IS NOT WEIGHTED')
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
end
AA = AA.*repmat(isfinite(dd),1,1,size(AA,3))./repmat(isfinite(dd),1,1,size(AA,3));
dd_fit = dd_fit.*isfinite(dd)./isfinite(dd);
end

%% Discarded material

% figure('Color','w')
% tiledlayout(1,4,'TileSpacing','tight')
% P_threshold = 0.05;
% 
% % Geostr. vs. unfiltered HFR
% nexttile
% plot(PIXELS_list, maxCorrCoef_g.PRCTILE(3,:), 'k.-', 'LineWidth',0.5*3 ,'MarkerSize',6*3 ); hold on
% plot(PIXELS_list, maxCorrCoef_g.PRCTILE(2,:), 'k.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_g.PRCTILE(4,:), 'k.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_g.PRCTILE(1,:), 'k.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_g.PRCTILE(5,:), 'k.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% set(gca,'ylim',[0 1])
% title('C(u$_\mathrm{geo}$, u$_\mathrm{HFR}$)','Interpreter','latex')
% 
% % cyclogeostr. vs. unfiltered HFR
% nexttile
% plot(PIXELS_list, maxCorrCoef_cg.PRCTILE(3,:), 'r.-', 'LineWidth',0.5*3 ,'MarkerSize',6*3 ); hold on
% plot(PIXELS_list, maxCorrCoef_cg.PRCTILE(2,:), 'r.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_cg.PRCTILE(4,:), 'r.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_cg.PRCTILE(1,:), 'r.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_cg.PRCTILE(5,:), 'r.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% set(gca,'ylim',[0 1])
% title('C(u$_\mathrm{cyclogeo}$, u$_\mathrm{HFR}$)','Interpreter','latex')
% 
% % Geostr. vs. filtered HFR
% nexttile
% plot(PIXELS_list, maxCorrCoef_g_lowpass.PRCTILE(3,:), 'k.-', 'LineWidth',0.5*3 ,'MarkerSize',6*3 ); hold on
% plot(PIXELS_list, maxCorrCoef_g_lowpass.PRCTILE(2,:), 'k.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_g_lowpass.PRCTILE(4,:), 'k.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_g_lowpass.PRCTILE(1,:), 'k.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_g_lowpass.PRCTILE(5,:), 'k.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% set(gca,'ylim',[0 1])
% title('C(u$_\mathrm{geo}$, u$_\mathrm{HFR, low pass}$)','Interpreter','latex')
% 
% % cyclogeostr. vs. filtered HFR
% nexttile
% plot(PIXELS_list, maxCorrCoef_cg_lowpass.PRCTILE(3,:), 'r.-', 'LineWidth',0.5*3 ,'MarkerSize',6*3 ); hold on
% plot(PIXELS_list, maxCorrCoef_cg_lowpass.PRCTILE(2,:), 'r.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_cg_lowpass.PRCTILE(4,:), 'r.-', 'LineWidth',0.5*2 ,'MarkerSize',6*2, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_cg_lowpass.PRCTILE(1,:), 'r.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% plot(PIXELS_list, maxCorrCoef_cg_lowpass.PRCTILE(5,:), 'r.-', 'LineWidth',0.5*1 ,'MarkerSize',6*1, 'HandleVisibility','off' );
% set(gca,'ylim',[0 1])
% title('C(u$_\mathrm{cyclogeo}$, u$_\mathrm{HFR, low pass}$)','Interpreter','latex')
