%% Data for visualization:

% HFR:
% https://cordc.ucsd.edu/projects/hfrnet/#

% SWOT:
% https://swot-calval.oceandatalab.com/?date=2023-07-10T12:00:00&timespan=2w&extent=-14945849.77604_4333355.5476237_-13184740.644594_5193119.2416557&center=-14065295.210317_4763237.3946397&zoom=7&products=3857_SWOT_L3_roughness!3857_SWOT_L3_ssha_unedited_nolr!3857_SWOT_L3_ssha_v1_0&opacity=100_100_100&stackLevel=10.04_10.05_120.27&selection=001

%% Load Northern California Data (all data together is ~2.2 GB)

% wget "https://hfrnet-tds.ucsd.edu/thredds/ncss/HFR/USWC/6km/hourly/RTV/HFRADAR_US_West_Coast_6km_Resolution_Hourly_RTV_best.ncd?var=hdop&var=number_of_radials&var=number_of_sites&var=u&var=v&north=42&west=-125.4&east=-122&south=37&disableProjSubset=on&horizStride=1&time_start=2023-04-01T00%3A00%3A00Z&time_end=2023-07-10T21%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf" -O "./HFR_NorCal_CalVal.nc"

% Give where this script is regardles of how it is run,
% credit to <https://www.mathworks.com/matlabcentral/
% answers/81148-get-path-from-running-script>.
mfile = mfilename('fullpath');
if contains(mfile,'LiveEditorEvaluationHelper')
    mfile = matlab.desktop.editor.getActiveFilename;
end
mfileDashInd = mfile == '/';
mfolder = mfile(1:max(mfileDashInd.*[1:length(mfile)]));
cd(mfolder); cd ../data;

tic

if ~exist('NORCAL','var')
T0 = datenum('2012-01-01 00:00:00'); % HFR 0-time
SWOT_calval_t0 = '08-Apr-2023 07:00:00';
SWOT_calval_tf = '09-Jul-2023 17:00:00';
% NORCAL.HFR = ncreadall('./HFR_NorCal_CalVal.nc');
NORCAL.HFR = ncreadall('./HFR_NorCal_CalVal_2023.nc');
calval_ind = [dsearchn(T0 + NORCAL.HFR.time/24, datenum(SWOT_calval_t0)):...
              dsearchn(T0 + NORCAL.HFR.time/24, datenum(SWOT_calval_tf))];
NORCAL.HFR.lat = flip(NORCAL.HFR.lat,1);
NORCAL.HFR.u = flip(permute(NORCAL.HFR.u,[2 1 3]),1);
    NORCAL.HFR.u_mean_2023 = mean(NORCAL.HFR.u,3,'omitmissing');
    % NORCAL.HFR.u = NORCAL.HFR.u - mean(NORCAL.HFR.u,3,'omitmissing');
    NORCAL.HFR.u =   NORCAL.HFR.u(:,:,calval_ind) ...
                   - NORCAL.HFR.u_mean_2023;
NORCAL.HFR.v = flip(permute(NORCAL.HFR.v,[2 1 3]),1);
    NORCAL.HFR.v_mean_2023 = mean(NORCAL.HFR.v,3,'omitmissing');
    % NORCAL.HFR.v = NORCAL.HFR.v - mean(NORCAL.HFR.v,3,'omitmissing');
    NORCAL.HFR.v =   NORCAL.HFR.v(:,:,calval_ind) ...
                   - NORCAL.HFR.v_mean_2023;
NORCAL.HFR.hdop = flip(permute(NORCAL.HFR.hdop,[2 1 3]),1);
NORCAL.HFR.hdop = NORCAL.HFR.hdop(:,:,calval_ind);
NORCAL.HFR.time_run = NORCAL.HFR.time_run(calval_ind);
NORCAL.HFR.time = NORCAL.HFR.time(calval_ind);
NORCAL.HFR.number_of_radials = flip(permute(NORCAL.HFR.number_of_radials,[2 1 3]),1);
NORCAL.HFR.number_of_sites = flip(permute(NORCAL.HFR.number_of_sites,[2 1 3]),1);
NORCAL.HFR = rmfield(NORCAL.HFR,{'number_of_radials','number_of_sites'}); % eliminate until otherwise needed
[NORCAL.HFR.LON,NORCAL.HFR.LAT] = meshgrid(NORCAL.HFR.lon,NORCAL.HFR.lat);

% % % First and last SWOT flyover average times:
% '08-Apr-2023 07:18:35'
% '09-Jul-2023 16:47:16'
NORCAL.SWOT = load('./NORCAL_SWOTdata_CCS.mat');
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
for ti_swot = 1:length(NORCAL.SWOT.ssha_karin_2)
    SSHA(:,:,ii) = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
                   [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
                   [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
end
NORCAL.SWOT.mean_SSHA = mean(SSHA,3,'omitnan'); clear SSHA

DIR_SST = dir('./SST/viirs_*_calval_mendocino/*.nc');
for file_i = 1:length(DIR_SST)
    FOLDER = DIR_SST(file_i).folder;
    FILE = DIR_SST(file_i).name;
    NORCAL.SST.longitude{file_i} = flip(flip(ncread([FOLDER '/' FILE],'/navigation_data/longitude'),2),1)';
    NORCAL.SST.latitude{file_i}  = flip(flip(ncread([FOLDER '/' FILE],'/navigation_data/latitude'),2),1)';
    NORCAL.SST.sst{file_i}       = flip(flip(ncread([FOLDER '/' FILE],'/geophysical_data/sst'),2),1)';
    NORCAL.SST.qual_sst{file_i}  = flip(flip(ncread([FOLDER '/' FILE],'/geophysical_data/qual_sst'),2),1)';

    YEAR = ncread([FOLDER '/' FILE],'/scan_line_attributes/year');
    DAY  = ncread([FOLDER '/' FILE],'/scan_line_attributes/day');
    MSEC = ncread([FOLDER '/' FILE],'/scan_line_attributes/msec');

    NORCAL.SST.datetime{file_i} = datetime(YEAR, 1, DAY, 0, 0, 0, MSEC);
end

NORCAL.SST.mean_datenum = nan(length(NORCAL.SST.datetime),1);
for ii = 1:length(NORCAL.SST.datetime)
    NORCAL.SST.mean_datenum(ii) = mean(datenum(NORCAL.SST.datetime{ii}),'omitnan');
end
previous_SST_ti = -1; % initialize this variable to check if recalculated grad(T) is needed

% Remove unused variables (comment lines if they are needed later):
NORCAL.SWOT = rmfield(NORCAL.SWOT,'time');
NORCAL.SWOT = rmfield(NORCAL.SWOT,'mean_sea_surface_cnescls');
% NORCAL.SST = rmfield(NORCAL.SST,'datetime');

whos

T0 = datenum('2012-01-01 00:00:00');
Collimate = @(IN) IN(:);

else
    warning('Variables were already loaded into memory.')
end

toc

% % % Determine whis SST passes have the most good data:
SST_is_good = nan(length(NORCAL.SST.sst),1);
SST_indices = [1:length(SST_is_good)]';
for ti_SST = 1:length(NORCAL.SST.sst)
    SST_mat = NORCAL.SST.sst{ti_SST}.*[[NORCAL.SST.qual_sst{ti_SST}<3]./[NORCAL.SST.qual_sst{ti_SST}<3]];
    SST_is_good(ti_SST) = sum(isfinite(SST_mat(:)))/numel(SST_mat);
end
figure
P3 = plot3(NORCAL.SST.mean_datenum,SST_is_good,SST_indices,'.-');
datetick
P3.Parent.CameraPosition = 1.0e+05 * ...
   [7.390235000000000   0.000003500000000   0.013465573435631];
[~,Low_to_High] = sort(SST_is_good);

close all

%% Calculate SWOT geostrophic velocity

% Recall that at the surface, geostrophic velocity is:
% u_g = -(g/f)(d\eta/dy)
% v_g = +(g/f)(d\eta/dx)
% where both derivatives are partial derivatives.
% Simple textbook citation:
% https://uw.pressbooks.pub/ocean285/chapter/geostrophic-balance/

close all

PIXEL_LIST = [2,4,6,8,10,12];

for pixel_i = 1:length(PIXEL_LIST)

% PIXELS = 4;
PIXELS = PIXEL_LIST(pixel_i);

SWOT_geostr_vel_calculated = true;
gg = 9.807; % m s^-2
GaussianKernelSTD_meters = 4000;

U_geostr = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
V_geostr = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% U_geostr_ = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% V_geostr_ = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
warning ('off','all')
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
                          SSHA, PIXELS, ... % 6 -> 13x13 box, i.e. [26 km]^2
                          true); % <--- final input forces weighted LSF (Gaussian decorr. with 3Sigma at edges)
    % To get geostrophic velocity from inputting surface height:
    U_geostr(:,:,ti_swot) = -(gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,3))/111000;
    V_geostr(:,:,ti_swot) =  (gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,2))./[111000*cosd(NORCAL.SWOT.lat{ti_swot})];


    disp(100*ti_swot/length(NORCAL.SWOT.ssha_karin_2))
end
warning ('on','all')

%% Example to verify that u_g is correct:

ti_swot = 70;

figure
% tiledlayout(1,2,'TileSpacing','compact')
% 
SSHA = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];

[~,SSHA] = fit_6term_2d(NORCAL.SWOT.lon{ti_swot}, ...
                      NORCAL.SWOT.lat{ti_swot}, ...
                      SSHA, PIXELS); % 6 -> 13x13 box, i.e. [26 km]^2
                      % PIXELS -> [2*PIXELS + 1]x[2*PIXELS + 1] box


pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                NORCAL.SWOT.lat{ti_swot}, ...
                SSHA); hold on; shading flat
quiver(NORCAL.SWOT.lon{ti_swot},...
       NORCAL.SWOT.lat{ti_swot},...
       U_geostr(:,:,ti_swot),...
       V_geostr(:,:,ti_swot),   5,'k')
set(gcf,'Position',[-1439          42         470         769])
title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))])
axis equal
colormap(bwr)
clim([-1 1]*0.2)

%% Apply the iterative method of Penven et al. (2014) to get the Cyclogeostrophic current
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

%% Interpolate the SWOT GEOSTROPHIC AND CYCLOGEOSTROPHIC velocities to the HFR grid to compare in the next step:

% if exist('HFR_scatter_times','var')
%     warning('You already calculated these variables, which take a long time to calculate.')
% else


close all

LON_LAT_lims = [-180 180 -90 90]; % all data (the entire world)
% LON_LAT_lims = [-125 -123.5 38 39.5]; % generous limits of where the southern eddy is
% LON_LAT_lims = [-124.55 -123.75 38.5 39.1]; % stricter limits of where the southern eddy is
% LON_LAT_lims = [-125 -124.25 40 40.5]; % stricter limits of where the northern eddy is
% LON_LAT_lims = [-125 -123.5 38.5 40.5]; % encompassing both eddies and between

T_start = datenum('2023-04-01 00:00:00'); % All Cal/Val
T_end   = datenum('2023-07-11 00:00:00'); % All Cal/Val
% T_start = datenum('2023-05-01 00:00:00'); % All May
% T_end   = datenum('2023-05-31 00:00:00'); % All May
% T_start = datenum('2023-05-09 00:00:00'); % Southern eddy duration
% T_end   = datenum('2023-05-24 00:00:00'); % Southern eddy duration

SWOT_scatter_times = [dsearchn(NORCAL.SWOT.mean_time,T_start):...
                      dsearchn(NORCAL.SWOT.mean_time,T_end  )];
HFR_scatter_times = [];

% SWOT_ug_all = []; % = U_geostr
% SWOT_vg_all = []; % = V_geostr
Ucg_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
Vcg_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
Anglecg_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
tic
for ti = 1:length(SWOT_scatter_times)
    ti_swot = SWOT_scatter_times(ti);
    ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));
    Delta_t(ti) = [T0 + NORCAL.HFR.time(ti_hfr)/24] - NORCAL.SWOT.mean_time(ti_swot);
    U_cg_slice = U_cyclogeostr_Nit(:,:,ti_swot);
    V_cg_slice = V_cyclogeostr_Nit(:,:,ti_swot);

    Ucg_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                                 NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)),...
                                        NORCAL.SWOT.lat{ti_swot}(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                                 NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)),...
                                        U_cg_slice(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                       NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)), ...
                                        NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
    Vcg_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                                 NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)),...
                                        NORCAL.SWOT.lat{ti_swot}(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                                 NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)),...
                                        V_cg_slice(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                       NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)), ...
                                        NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
    Anglecg_SWOT_HFRgrid = angle(Ucg_SWOT_HFRgrid + 1i*Vcg_SWOT_HFRgrid);

    % Ucg_SWOT_HFRgrid_all = [Ucg_SWOT_HFRgrid_all; Ucg_SWOT_HFRgrid(NORCAL.HFR.LON>LON_LAT_lims(1) & NORCAL.HFR.LON<LON_LAT_lims(2) & NORCAL.HFR.LAT>LON_LAT_lims(3) & NORCAL.HFR.LAT<LON_LAT_lims(4))];
    % Vcg_SWOT_HFRgrid_all = [Vcg_SWOT_HFRgrid_all; Vcg_SWOT_HFRgrid(NORCAL.HFR.LON>LON_LAT_lims(1) & NORCAL.HFR.LON<LON_LAT_lims(2) & NORCAL.HFR.LAT>LON_LAT_lims(3) & NORCAL.HFR.LAT<LON_LAT_lims(4))];
    % Anglecg_SWOT_HFRgrid_all = [Anglecg_SWOT_HFRgrid_all; Anglecg_SWOT_HFRgrid(NORCAL.HFR.LON>LON_LAT_lims(1) & NORCAL.HFR.LON<LON_LAT_lims(2) & NORCAL.HFR.LAT>LON_LAT_lims(3) & NORCAL.HFR.LAT<LON_LAT_lims(4))];
    Ucg_SWOT_HFRgrid_all(:,:,ti) = Ucg_SWOT_HFRgrid;
    Vcg_SWOT_HFRgrid_all(:,:,ti) = Vcg_SWOT_HFRgrid;
    Anglecg_SWOT_HFRgrid_all(:,:,ti) = Anglecg_SWOT_HFRgrid;

    HFR_scatter_times = [HFR_scatter_times; ti_hfr];

    disp(100*ti/length(SWOT_scatter_times))
end

Ug_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
Vg_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
Angleg_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
tic
for ti = 1:length(SWOT_scatter_times)
    ti_swot = SWOT_scatter_times(ti);
    ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));
    Delta_t(ti) = [T0 + NORCAL.HFR.time(ti_hfr)/24] - NORCAL.SWOT.mean_time(ti_swot);
    U_geostr_slice = U_geostr(:,:,ti_swot);
    V_geostr_slice = V_geostr(:,:,ti_swot);

    Ug_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                                NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)),...
                                       NORCAL.SWOT.lat{ti_swot}(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                                NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)),...
                                       U_geostr_slice(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                      NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)), ...
                                       NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
    Vg_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                                NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)),...
                                       NORCAL.SWOT.lat{ti_swot}(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                                NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)),...
                                       V_geostr_slice(NORCAL.SWOT.lon{ti_swot}>LON_LAT_lims(1) & NORCAL.SWOT.lon{ti_swot}<LON_LAT_lims(2) & ...
                                                      NORCAL.SWOT.lat{ti_swot}>LON_LAT_lims(3) & NORCAL.SWOT.lat{ti_swot}<LON_LAT_lims(4)), ...
                                       NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
    Angleg_SWOT_HFRgrid = angle(Ug_SWOT_HFRgrid + 1i*Vg_SWOT_HFRgrid);

    % Ug_SWOT_HFRgrid_all = [Ug_SWOT_HFRgrid_all; Ug_SWOT_HFRgrid(NORCAL.HFR.LON>LON_LAT_lims(1) & NORCAL.HFR.LON<LON_LAT_lims(2) & NORCAL.HFR.LAT>LON_LAT_lims(3) & NORCAL.HFR.LAT<LON_LAT_lims(4))];
    % Vg_SWOT_HFRgrid_all = [Vg_SWOT_HFRgrid_all; Vg_SWOT_HFRgrid(NORCAL.HFR.LON>LON_LAT_lims(1) & NORCAL.HFR.LON<LON_LAT_lims(2) & NORCAL.HFR.LAT>LON_LAT_lims(3) & NORCAL.HFR.LAT<LON_LAT_lims(4))];
    % Angleg_SWOT_HFRgrid_all = [Angleg_SWOT_HFRgrid_all; Angleg_SWOT_HFRgrid(NORCAL.HFR.LON>LON_LAT_lims(1) & NORCAL.HFR.LON<LON_LAT_lims(2) & NORCAL.HFR.LAT>LON_LAT_lims(3) & NORCAL.HFR.LAT<LON_LAT_lims(4))];
    Ug_SWOT_HFRgrid_all(:,:,ti) = Ug_SWOT_HFRgrid;
    Vg_SWOT_HFRgrid_all(:,:,ti) = Vg_SWOT_HFRgrid;
    Angleg_SWOT_HFRgrid_all(:,:,ti) = Angleg_SWOT_HFRgrid;

    % HFR_scatter_times = [HFR_scatter_times; ti_hfr];

    disp(100*ti/length(SWOT_scatter_times))
end

toc

% end

%% Define the filtered part of the HFR current to compare to SWOT geostrophic velocity

close all

if exist('WINDOW','var')
    warning('You already calculated these variables, which take a long time to calculate.')
else


WINDOW = hanning(36); WINDOW = WINDOW/sum(WINDOW); % very low pass
% WINDOW = hanning(24); WINDOW = WINDOW/sum(WINDOW); % low pass
% WINDOW = hanning(6); WINDOW = WINDOW/sum(WINDOW); % tides pass
% ConvWithWindow = @(IN) conv(squeeze(IN),WINDOW,'same');
ConvWithWindow = @(IN) nanconv(squeeze(IN),WINDOW,'same','edge');
% ConvWithWindow = @(IN) IN; % control to test lack of convolution

% Low-pass only:
L_U = nan(size(NORCAL.HFR.u));
L_V = nan(size(NORCAL.HFR.v));
tic
for i1 = 1:size(NORCAL.HFR.LAT,1)
    for i2 = 1:size(NORCAL.HFR.LAT,2)
        % if sum(isfinite(squeeze(U_rot_unfiltered(i1,i2,3:end))))/size(squeeze(U_rot_unfiltered(i1,i2,3:end)),3) > 0.5
            L_U(i1,i2,:) = ConvWithWindow(NORCAL.HFR.u(i1,i2,:) - 1.0*mean(NORCAL.HFR.u(i1,i2,:),3,'omitnan'));
            L_V(i1,i2,:) = ConvWithWindow(NORCAL.HFR.v(i1,i2,:) - 1.0*mean(NORCAL.HFR.v(i1,i2,:),3,'omitnan'));
        % else
        % end
    end
    disp(100*i1/size(NORCAL.HFR.u,1))
end
toc
end

%% Helmholtz Decomposition at all times

if exist('U_rot_unfiltered','var')
    warning('You already calculated these variables, which take a long time to calculate.')
else


U_rot_filtered = [];
V_rot_filtered = [];
U_div_filtered = [];
V_div_filtered = [];
U_rot_unfiltered = [];
V_rot_unfiltered = [];
U_div_unfiltered = [];
V_div_unfiltered = [];
tic
for ti_hfr = 1:size(NORCAL.HFR.u,3)
    [~, ~, udiv, vdiv, urot, vrot] = helmholtz_Cgrid(...
        L_U(:,:,ti_hfr), ...
        L_V(:,:,ti_hfr),...
        isfinite(L_U(:,:,ti_hfr)),...
        6000,... *ones(size(NORCAL.HFR.u(:,:,ti_hfr)))
        6000,...
        6000,...
        6000, [false false]); % , 'VERBOSE'
    urot = urot.*urot./urot;
    vrot = vrot.*vrot./vrot;
    udiv = udiv.*udiv./udiv;
    vdiv = vdiv.*vdiv./vdiv;

    U_rot_filtered(:,:,ti_hfr) = urot;
    V_rot_filtered(:,:,ti_hfr) = vrot;
    U_div_filtered(:,:,ti_hfr) = udiv;
    V_div_filtered(:,:,ti_hfr) = vdiv;
    % disp(100*ti_hfr/size(NORCAL.HFR.u,3))

    [~, ~, udiv, vdiv, urot, vrot] = helmholtz_Cgrid(...
        NORCAL.HFR.u(:,:,ti_hfr) - 0*mean(NORCAL.HFR.u,3,'omitnan'), ...
        NORCAL.HFR.v(:,:,ti_hfr) - 0*mean(NORCAL.HFR.v,3,'omitnan'),...
        isfinite(NORCAL.HFR.u(:,:,ti_hfr)),...
        6000,... *ones(size(NORCAL.HFR.u(:,:,ti_hfr)))
        6000,...
        6000,...
        6000, [false false]); % , 'VERBOSE'
    urot = urot.*urot./urot;
    vrot = vrot.*vrot./vrot;
    udiv = udiv.*udiv./udiv;
    vdiv = vdiv.*vdiv./vdiv;

    U_rot_unfiltered(:,:,ti_hfr) = urot;
    V_rot_unfiltered(:,:,ti_hfr) = vrot;
    U_div_unfiltered(:,:,ti_hfr) = udiv;
    V_div_unfiltered(:,:,ti_hfr) = vdiv;

    disp(100*ti_hfr/size(NORCAL.HFR.u,3))
end
toc

end

clear urot vrot udiv vdiv

% error('Forced stop by user')

%% % Save time series of SWOT g and cg velocity, and if already defined,
% % % the equivalent velocities from HFR at SWOT times
% error
if exist('U_rot_unfiltered','var')
    % pwd = .../SWOT_HFR_Analysis/data
    % HFR at/near SWOT times
    U_rot_filtered_SWOTtimes = nan(size(Ug_SWOT_HFRgrid_all));
    V_rot_filtered_SWOTtimes = nan(size(Ug_SWOT_HFRgrid_all));
    U_rot_unfiltered_SWOTtimes = nan(size(Ug_SWOT_HFRgrid_all));
    V_rot_unfiltered_SWOTtimes = nan(size(Ug_SWOT_HFRgrid_all));
    tic
    for ti = 1:length(NORCAL.SWOT.mean_time)
        ti_HFR_closest = dsearchn(NORCAL.HFR.time/24 + T0,NORCAL.SWOT.mean_time(ti));
        U_rot_filtered_SWOTtimes(:,:,ti) = U_rot_filtered(:,:,ti_HFR_closest);
        V_rot_filtered_SWOTtimes(:,:,ti) = V_rot_filtered(:,:,ti_HFR_closest);
        U_rot_unfiltered_SWOTtimes(:,:,ti) = U_rot_unfiltered(:,:,ti_HFR_closest);
        V_rot_unfiltered_SWOTtimes(:,:,ti) = V_rot_unfiltered(:,:,ti_HFR_closest);
    end
    toc
    SWOT_time = NORCAL.SWOT.mean_time;
    tic
    % Save the relevant variables; note that the files
    % <SWOT_and_HFR_velocities_*pix.mat> are also generated and saved with
    % this script, but the user will need to re-run this script without the
    % weighting option in "fit_6term_2d", and then save the variables below
    % but without the "_weighted" string in the file name.
    save(['./SWOT_and_HFR_velocities_' num2str(PIXELS) 'pix_weighted.mat'],... save(['./SWOT_and_HFR_velocities_' num2str(PIXELS) 'pix.mat'],...
         'SWOT_time','WINDOW','ConvWithWindow',...
         'Ug_SWOT_HFRgrid_all','Vg_SWOT_HFRgrid_all',...
         'Ucg_SWOT_HFRgrid_all','Vcg_SWOT_HFRgrid_all',...
         'PIXELS','U_geostr','V_geostr',...
         'CG_Iteration','U_cyclogeostr_Nit','V_cyclogeostr_Nit',...
         'U_rot_filtered_SWOTtimes','V_rot_filtered_SWOTtimes',...
         'U_rot_unfiltered_SWOTtimes','V_rot_unfiltered_SWOTtimes')
    toc
    % error('Forced stop')
else
    warning('No time for this step yet')
end

end % pixel_i loop

error('Forced stop')

%% Map {U,V}_geostr interpolated to HFR grid alongside HFR data

% Chose whether to plot 2km or 6km (smoothed to HFR specifications) u,v_g
PLOT_GRIDDED_GEOSTROPHIC_VELOCITY = true;

Antenna_Locations = ...
    [1i*42.836667 - 124.5655; ...
     1i*42.828883 - 70.8145; ...
     1i*42.58625 - 70.650556; ...
     1i*42.159917 - 70.704883; ...
     1i*42.078833 - 70.220333; ...
     1i*41.844217 - 69.947233; ...
     1i*41.7845 - 124.254817; ...
     1i*41.5088 - 71.0687; ...
     1i*41.447989 - 71.431856; ...
     1i*41.3484 - 70.640167; ...
     1i*41.322933 - 71.804217; ...
     1i*41.25045 - 69.974067; ...
     1i*41.241933 - 70.106903; ...
     1i*41.1528 - 71.551183; ...
     1i*41.152583 - 71.551833; ...
     1i*41.073567 - 124.157783; ...
     1i*41.070967 - 71.856733; ...
     1i*40.98255 - 73.623783; ...
     1i*40.969333 - 72.1237; ...
     1i*40.90875 - 73.587283; ...
     1i*40.788183 - 72.74525; ...
     1i*40.768783 - 124.218767; ...
     1i*40.586783 - 73.589817; ...
     1i*40.543633 - 74.124533; ...
     1i*40.462117 - 74.2535; ...
     1i*40.4413 - 74.099333; ...
     1i*40.4301 - 73.983583; ...
     1i*40.366817 - 73.973517; ...
     1i*40.195583 - 74.00885; ...
     1i*40.033367 - 124.078867; ...
     1i*39.9325 - 74.07255; ...
     1i*39.736233 - 74.11705; ...
     1i*39.532033 - 74.262133; ...
     1i*39.438017 - 123.816117; ...
     1i*39.378417 - 74.398517; ...
     1i*39.378367 - 74.399017; ...
     1i*39.198867 - 74.652967; ...
     1i*39.073767 - 74.912367; ...
     1i*38.99095 - 74.799; ...
     1i*38.94215 - 74.886117; ...
     1i*38.931267 - 74.961333; ...
     1i*38.928433 - 123.727767; ...
     1i*38.794033 - 75.089183; ...
     1i*38.781783 - 75.12825; ...
     1i*38.567183 - 123.33155; ...
     1i*38.319533 - 123.073617; ...
     1i*38.317383 - 123.072317; ...
     1i*38.205017 - 75.1529; ...
     1i*38.04715 - 122.989133; ...
     1i*37.8899 - 122.445933; ...
     1i*37.872433 - 122.597717; ...
     1i*37.853317 - 122.4192; ...
     1i*37.84395 - 122.4772; ...
     1i*37.815483 - 122.5299; ...
     1i*37.806367 - 122.465967; ...
     1i*37.8029 - 122.3971; ...
     1i*37.712267 - 122.501317; ...
     1i*37.533717 - 122.519217; ...
     1i*37.49675 - 122.499367; ...
     1i*37.1379 - 75.972133; ...
     1i*37.088517 - 122.2734];
LatRange = [37 42];
LonRange = [min(NORCAL.HFR.lon(:)) -122];
Antenna_Locations = Antenna_Locations(real(Antenna_Locations) > min(LonRange) & ...
                                      real(Antenna_Locations) < max(LonRange) & ...
                                      imag(Antenna_Locations) > min(LatRange) & ...
                                      imag(Antenna_Locations) < max(LatRange));

USER_PICKED_DATE = '2023-05-03 02:00:00';
USER_PICKED_DATE = '2023-05-09 02:00:00';
USER_PICKED_DATE = '2023-05-10 02:00:00';
% USER_PICKED_DATE = '2023-05-15 02:00:00';
% USER_PICKED_DATE = '2023-05-21 02:00:00';
% USER_PICKED_DATE = '2023-05-26 02:00:00';
% USER_PICKED_DATE = '2023-06-01 22:00:00';
% USER_PICKED_DATE = '2023-06-04 22:00:00';
% USER_PICKED_DATE = '2023-06-16 00:00:00';
% USER_PICKED_DATE = '2023-06-26 19:00:00';
% USER_PICKED_DATE = '2023-07-05 17:00:00';
% USER_PICKED_DATE = '2023-07-09 17:00:00';
ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum(USER_PICKED_DATE));
ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));

Ug_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot},...
                                   U_geostr(:,:,ti_swot), ...
                                   NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
Vg_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot},...
                                   V_geostr(:,:,ti_swot), ...
                                   NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);

Ucg_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot},...
                                    U_cyclogeostr_Nit(:,:,ti_swot), ...
                                    NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
Vcg_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot},...
                                    V_cyclogeostr_Nit(:,:,ti_swot), ...
                                    NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);

% figure;
% subplot(121); imagesc(U_geostr(:,:,ti_swot));colorbar;clim([-1 1])
% subplot(122); imagesc(Ug_SWOT_HFRgrid);colorbar;clim([-1 1])

close all

figure('Color','w')
tiledlayout(1,2,'TileSpacing','tight')
VecScale = 0.12;
AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                  NORCAL.SWOT.lat{ti_swot},...
                  sqrt(U_geostr(:,:,ti_swot).^2 + V_geostr(:,:,ti_swot).^2));
% m_quiver(NORCAL.SWOT.lon{ti_swot},...
%          NORCAL.SWOT.lat{ti_swot},...
%          VecScale*U_geostr(:,:,ti_swot),VecScale*V_geostr(:,:,ti_swot),0,'k')
m_quiver(NORCAL.SWOT.lon{ti_swot}(1:3:end,1:3:end),...
         NORCAL.SWOT.lat{ti_swot}(1:3:end,1:3:end),...
         VecScale*U_geostr(1:3:end,1:3:end,ti_swot),VecScale*V_geostr(1:3:end,1:3:end,ti_swot),0,'k')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
m_text(-122.2568, 41.8015, '(a)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center');
% CB = colorbar; CB.Label.String = 'm s$^{-1}$'; CB.Label.FontSize = 20; CB.Label.Interpreter = 'LaTeX';
CLIM = [0 1]*1;
clim(CLIM);
xlabel('SWOT: $\bf{u}_\mathrm{g}$','Interpreter','latex','Position',[2.4813e-08 0.0515 -1.0000e+16])
set(gca,'FontSize',16)
ylabel(['SWOT time: ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'Interpreter','latex','VerticalAlignment','top')

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                  NORCAL.SWOT.lat{ti_swot},...
                  sqrt(U_cyclogeostr_Nit(:,:,ti_swot).^2 + V_cyclogeostr_Nit(:,:,ti_swot).^2));
% m_quiver(NORCAL.SWOT.lon{ti_swot},...
%          NORCAL.SWOT.lat{ti_swot},...
%          VecScale*U_cyclogeostr_Nit(:,:,ti_swot),VecScale*V_cyclogeostr_Nit(:,:,ti_swot),0,'k')
m_quiver(NORCAL.SWOT.lon{ti_swot}(1:3:end,1:3:end),...
         NORCAL.SWOT.lat{ti_swot}(1:3:end,1:3:end),...
         VecScale*U_cyclogeostr_Nit(1:3:end,1:3:end,ti_swot),VecScale*V_cyclogeostr_Nit(1:3:end,1:3:end,ti_swot),0,'k')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
m_text(-122.2568, 41.8015, '(b)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center');
CB = colorbar; CB.Label.String = 'm s$^{-1}$'; CB.Label.FontSize = 20; CB.Label.Interpreter = 'LaTeX';
clim(CLIM);
xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
set(gca,'FontSize',16)
% title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'Interpreter','latex')
linkaxes([AX1 AX2],'xy')
colormap(turbo)

set(gcf,'Position',[-1336 273 711 540])

% %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

figure('Color','w')
% tiledlayout(1,3,'TileSpacing','compact')
tiledlayout(2,2,'TileSpacing','tight')
% Small function to plot the edges of a matrix (recall that <Collimate = @(IN) IN(:);>):
Edges_Only = @(IN) [Collimate(IN(1,:)); Collimate(IN(:,end)); flip(Collimate(IN(end,:))); flip(Collimate(IN(:,1)))];

VecScale = 0.12;

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
if PLOT_GRIDDED_GEOSTROPHIC_VELOCITY
    m_pcolor_centered(NORCAL.HFR.LON,...
                      NORCAL.HFR.LAT,...
                      sqrt(Ug_SWOT_HFRgrid.^2 + Vg_SWOT_HFRgrid.^2));
    m_quiver(NORCAL.HFR.LON,...
             NORCAL.HFR.LAT,...
             VecScale*Ug_SWOT_HFRgrid,VecScale*Vg_SWOT_HFRgrid,0,'k')
    % m_contour(NORCAL.HFR.LON,...
    %           NORCAL.HFR.LAT,...
    %           sqrt(abs(Ug_SWOT_HFRgrid + 1i*Vg_SWOT_HFRgrid)),[0:0.1:1],'k')
    % M_PLOT = m_plot(real(Antenna_Locations),imag(Antenna_Locations),'k.','MarkerSize',25,'MarkerFaceColor','k');
    xlabel('SWOT: $\bf{u}_\mathrm{g}$ on HFR grid','Interpreter','latex','Position',[2.4813e-08 0.0515 -1.0000e+16])
else
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                      NORCAL.SWOT.lat{ti_swot},...
                      sqrt(U_geostr(:,:,ti_swot).^2 + V_geostr(:,:,ti_swot).^2));
    m_quiver(NORCAL.SWOT.lon{ti_swot}(1:3:end,1:3:end),...
             NORCAL.SWOT.lat{ti_swot}(1:3:end,1:3:end),...
             VecScale*U_geostr(1:3:end,1:3:end,ti_swot),VecScale*V_geostr(1:3:end,1:3:end,ti_swot),0,'k')
    xlabel('SWOT: $\bf{u}_\mathrm{g}$','Interpreter','latex','Position',[2.4813e-08 0.0515 -1.0000e+16])
end
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
m_text(-122.2568, 41.8015, '(a)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center');
% CB = colorbar; CB.Label.String = 'm s$^{-1}$'; CB.Label.FontSize = 20; CB.Label.Interpreter = 'LaTeX';
CLIM = [0 1]*1;
clim(CLIM);
set(gca,'FontSize',16)
ylabel(['SWOT time: ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'Interpreter','latex','VerticalAlignment','top')

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
if PLOT_GRIDDED_GEOSTROPHIC_VELOCITY
    m_pcolor_centered(NORCAL.HFR.LON,...
                      NORCAL.HFR.LAT,...
                      sqrt(Ucg_SWOT_HFRgrid.^2 + Vcg_SWOT_HFRgrid.^2));
    m_quiver(NORCAL.HFR.LON,...
             NORCAL.HFR.LAT,...
             VecScale*Ucg_SWOT_HFRgrid,VecScale*Vcg_SWOT_HFRgrid,0,'k')
    % m_contour(NORCAL.HFR.LON,...
    %           NORCAL.HFR.LAT,...
    %           sqrt(abs(Ucg_SWOT_HFRgrid + 1i*Vcg_SWOT_HFRgrid)),[0:0.1:1],'k')
    % M_PLOT = m_plot(real(Antenna_Locations),imag(Antenna_Locations),'k.','MarkerSize',25,'MarkerFaceColor','k');
    xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
else
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                      NORCAL.SWOT.lat{ti_swot},...
                      sqrt(U_cyclogeostr_Nit(:,:,ti_swot).^2 + V_cyclogeostr_Nit(:,:,ti_swot).^2));
    m_quiver(NORCAL.SWOT.lon{ti_swot}(1:3:end,1:3:end),...
             NORCAL.SWOT.lat{ti_swot}(1:3:end,1:3:end),...
             VecScale*U_cyclogeostr_Nit(1:3:end,1:3:end,ti_swot),VecScale*V_cyclogeostr_Nit(1:3:end,1:3:end,ti_swot),0,'k')
    xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
end
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
m_text(-122.2568, 41.8015, '(b)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center');
% CB = colorbar; CB.Label.String = 'm s$^{-1}$'; CB.Label.FontSize = 20; CB.Label.Interpreter = 'LaTeX';
clim(CLIM);
set(gca,'FontSize',16)
% title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'Interpreter','latex')

AX3 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,...
                  NORCAL.HFR.LAT,...
                  sqrt(U_rot_unfiltered(:,:,ti_hfr).^2 + V_rot_unfiltered(:,:,ti_hfr).^2));
m_quiver(NORCAL.HFR.LON,...
         NORCAL.HFR.LAT,...
         VecScale*U_rot_unfiltered(:,:,ti_hfr),VecScale*V_rot_unfiltered(:,:,ti_hfr),0,'k')
M_PLOT = m_plot(real(Antenna_Locations),imag(Antenna_Locations),'k.','MarkerSize',25,'MarkerFaceColor','k');
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'w','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'w--','LineWidth',1)
xlabel('HFR: $\bf{u}_\mathrm{rot}$','Interpreter','latex','Position',[2.4813e-08 0.0515 -1.0000e+16])
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
m_text(-122.2568, 41.8015, '(c)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center');
% CB = colorbar; CB.Label.String = 'm s$^{-1}$'; CB.Label.FontSize = 20; CB.Label.Interpreter = 'LaTeX';
clim(CLIM);
set(gca,'FontSize',16)
ylabel(['HFR time: ' datestr(T0 + NORCAL.HFR.time(ti_hfr)/24)],'Interpreter','latex','VerticalAlignment','top')


AX4 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,...
                  NORCAL.HFR.LAT,...
                  sqrt(U_rot_filtered(:,:,ti_hfr).^2 + V_rot_filtered(:,:,ti_hfr).^2));
m_quiver(NORCAL.HFR.LON,...
         NORCAL.HFR.LAT,...
         VecScale*U_rot_filtered(:,:,ti_hfr),VecScale*V_rot_filtered(:,:,ti_hfr),0,'k')
M_PLOT = m_plot(real(Antenna_Locations),imag(Antenna_Locations),'k.','MarkerSize',25,'MarkerFaceColor','k');
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'w','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'w--','LineWidth',1)
xlabel('HFR: $\bf{u}_\mathrm{LP,rot}$','Interpreter','latex','Position',[2.4813e-08 0.0515 -1.0000e+16])
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
m_text(-122.2568, 41.8015, '(d)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center');
CB = colorbar; CB.Label.String = 'm s$^{-1}$'; CB.Label.FontSize = 20; CB.Label.Interpreter = 'LaTeX';
clim(CLIM);
set(gca,'FontSize',16)

colormap('turbo')

linkaxes([AX1 AX2 AX3 AX4],'xy')
% set([AX1 AX2 AX3 AX4], 'XAxisLocation', 'top')

% set(gcf,'Position',[-1708         117        1321         757])
% set(gcf,'Position',[113          62        1328         735])
% set(gcf,'Position',[113          62        1321         757])
% set(gcf,'Position',[389          62        1045         735])
% set(gcf,'Position',[509    62   925   735])
% set(gcf,'Position',[-646     1   647   976])
set(gcf,'Position',[-597    83   598   894])

% %%
% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/SWOT_ug_ucg.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% %%
% figure(2)
% exportgraphics(gcf,...
% '../figures/draft/SWOT_ug_ucg_HFR_rotfiltered.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% Scatter GEOSTR. vs. CYCLO.

close all
uv_bins = [-2:0.05:2];
speed_bins = [0:0.05:3];

figure
tiledlayout(1,4,'TileSpacing','compact')

nexttile
HIST_u = ...
histogram2(U_geostr(:),U_cyclogeostr_Nit(:),...
... histogram2(round(U_geostr(:)),U_cyclogeostr_Nit(:),...
           uv_bins,uv_bins,...
           'DisplayStyle','tile','Normalization','pdf');
axis equal; xlabel('u_g'); ylabel('u_{cg}')

nexttile
HIST_v = ...
histogram2(V_geostr(:),V_cyclogeostr_Nit(:),...
           uv_bins,uv_bins,...
           'DisplayStyle','tile','Normalization','pdf');
axis equal; xlabel('v_g'); ylabel('v_{cg}')

nexttile
HIST_speed = ...
histogram2(abs(U_geostr(:) + 1i*V_geostr(:)),...
           abs(U_cyclogeostr_Nit(:) + 1i*V_cyclogeostr_Nit(:)),...
           speed_bins,speed_bins,...
           'DisplayStyle','tile','Normalization','pdf');
axis equal; xlabel('|\bf{u}_g|'); ylabel('|\bf{u}_{cg}|')

nexttile
HIST_angle = ...
histogram2([180/pi]*angle(U_geostr(:) + 1i*V_geostr(:)), ...
           [180/pi]*angle(U_cyclogeostr_Nit(:) + 1i*V_cyclogeostr_Nit(:)),...
           [-180:4.5:180],[-180:4.5:180],'DisplayStyle','tile','Normalization','pdf');
axis equal; xlabel('\theta_g'); ylabel('\theta_{cg}')

colormap('turbo')


CLIM = [-3.5 1];

figure('Color','w')
tiledlayout(1,4,'TileSpacing','compact')

AX1 = nexttile;
XBinEdges = [HIST_u.XBinEdges(2:end) + HIST_u.XBinEdges(1:[end-1])]/2;
YBinEdges = [HIST_u.YBinEdges(2:end) + HIST_u.YBinEdges(1:[end-1])]/2;
[XBinEdges,YBinEdges] = meshgrid(XBinEdges,YBinEdges);
pcolor_centered(XBinEdges,YBinEdges,log10(HIST_u.Values).*HIST_u.Values./HIST_u.Values); hold on
% for r_times_f = [5000 20000 Inf]*fcor_degrees_cps(39.5) % [0.1 1 4 Inf]
%     plot(   [0:0.01:2],    2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
%     plot(-1*[0:0.01:2], -1*2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
% end
% text(-1.75,1.75,'(a)','Interpreter','latex','FontSize',20,'Color','w','VerticalAlignment','middle','HorizontalAlignment','center')
% TTT1=text(-1.4768,-0.5008,'$r$ = 5 km','Interpreter','latex','FontSize',12,'Color',    'w','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',12);
% TTT2=text(-1.4551,-0.8180,'$r$ = 20 km','Interpreter','latex','FontSize',12,'Color',   'w','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',26);
% TTT3=text(-1.6659,-1.5051,'$r$ = $\infty$','Interpreter','latex','FontSize',12,'Color','w','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',45);
clim(CLIM)
axis equal
ylabel('$u_\mathrm{g}$ (m/s)','Interpreter','latex','FontSize',20);
xlabel('$u_\mathrm{cg}$','Interpreter','latex','FontSize',20);
CB = colorbar;
CB.Label.String = 'log$_{10}$(PDF)'; CB.Label.Interpreter = 'latex'; CB.Label.FontSize = 14;
CB.Location = "northoutside";
pbaspect([1 1 1]);

AX2 = nexttile;
XBinEdges = [HIST_v.XBinEdges(2:end) + HIST_v.XBinEdges(1:[end-1])]/2;
YBinEdges = [HIST_v.YBinEdges(2:end) + HIST_v.YBinEdges(1:[end-1])]/2;
[XBinEdges,YBinEdges] = meshgrid(XBinEdges,YBinEdges);
pcolor_centered(XBinEdges,YBinEdges,log10(HIST_v.Values).*HIST_v.Values./HIST_v.Values); hold on
% for r_times_f = [5000 20000 Inf]*fcor_degrees_cps(39.5)
%     plot(   [0:0.01:2],    2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
%     plot(-1*[0:0.01:2], -1*2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
% end
clim(CLIM)
text(-1.75,1.75,'(b)','Interpreter','latex','FontSize',20,'Color','w','VerticalAlignment','middle','HorizontalAlignment','center')
axis equal
ylabel('$v_\mathrm{g}$ (m/s)','Interpreter','latex','FontSize',20);
xlabel('$v_\mathrm{cg}$','Interpreter','latex','FontSize',20);
CB = colorbar;
CB.Label.String = ''; % 'log$_{10}$(PDF)';
CB.Label.Interpreter = 'latex'; CB.Label.FontSize = 14; CB.Location = "northoutside";
pbaspect([1 1 1]);

AX3 = nexttile;
XBinEdges = [HIST_speed.XBinEdges(2:end) + HIST_speed.XBinEdges(1:[end-1])]/2;
YBinEdges = [HIST_speed.YBinEdges(2:end) + HIST_speed.YBinEdges(1:[end-1])]/2;
[XBinEdges,YBinEdges] = meshgrid(XBinEdges,YBinEdges);
pcolor_centered(XBinEdges,YBinEdges,log10(HIST_speed.Values).*HIST_speed.Values./HIST_speed.Values); hold on
% for r_times_f = [5000 20000 Inf]*fcor_degrees_cps(39.5)
%     plot(   [0:0.01:2],    2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
%     plot(-1*[0:0.01:2], -1*2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
% end
clim(CLIM)
text(0.225,2.825,'(c)','Interpreter','latex','FontSize',20,'Color','w','VerticalAlignment','middle','HorizontalAlignment','center')
axis equal
ylabel('$\mathbf{u}_\mathrm{g}$ (m/s)','Interpreter','latex','FontSize',20);
xlabel('$\mathbf{u}_\mathrm{cg}$','Interpreter','latex','FontSize',20);
CB = colorbar;
CB.Label.String = ''; % 'log$_{10}$(PDF)';
CB.Label.Interpreter = 'latex'; CB.Label.FontSize = 14; CB.Location = "northoutside";
pbaspect([1 1 1]);

AX4 = nexttile;
XBinEdges = [HIST_angle.XBinEdges(2:end) + HIST_angle.XBinEdges(1:[end-1])]/2;
YBinEdges = [HIST_angle.YBinEdges(2:end) + HIST_angle.YBinEdges(1:[end-1])]/2;
[XBinEdges,YBinEdges] = meshgrid(XBinEdges,YBinEdges);
pcolor_centered(XBinEdges,YBinEdges,log10(HIST_angle.Values).*HIST_angle.Values./HIST_angle.Values);
clim([-7.5 -4])
text(-157.5,157.5,'(d)','Interpreter','latex','FontSize',20,'Color','k','VerticalAlignment','middle','HorizontalAlignment','center')
axis equal
ylabel('$\theta_\mathrm{g}$ ($^\circ$)','Interpreter','latex','FontSize',20);
xlabel('$\theta_\mathrm{cg}$','Interpreter','latex','FontSize',20);
set(gca,'XTick',[-180:90:180],'YTick',[-180:90:180])
CB = colorbar;
CB.Label.String = ''; % 'log$_{10}$(PDF)';
CB.Label.Interpreter = 'latex'; CB.Label.FontSize = 14; CB.Location = "northoutside";
CB.Ticks = [-7:-4];
pbaspect([1 1 1]);
colormap('turbo')

% sgtitle('THE AXES ARE FLIPPED, DO NOT USE UNTIL THIS IS SORTED OUT','Color','r')

linkaxes([AX1 AX2],'xy')
AX1.Color = [1 1 1]*0.5;
AX2.Color = [1 1 1]*0.5;
AX3.Color = [1 1 1]*0.5;
AX4.Color = [1 1 1]*0.5;
set(gcf,'Position',[-1439         271         942         540])

% %%
% figure(2)
% exportgraphics(gcf,...
% '../figures/draft/F_UVAngle_GvsCG.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp(['Image saved'])

%% Plot time series side-by-side

close all
figure('Color','w')
% tiledlayout(1,2,'TileSpacing','compact')

VecScale = 0.12;
Plot_uvMagnitude_Cuttoff = 0.0; % For every velocity magnitude below this, u,v,angle will not be plotted

AX0 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-1 1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,...
                  NORCAL.HFR.LAT,...
                  sqrt(Ug_SWOT_HFRgrid.^2 + Vg_SWOT_HFRgrid.^2));
m_contour(NORCAL.HFR.LON, NORCAL.HFR.LAT, sum(isfinite(NORCAL.HFR.u),3)./size(NORCAL.HFR.u,3),[0:0.1:1],'k','LineWidth',2)
m_quiver(NORCAL.HFR.LON,...
         NORCAL.HFR.LAT,...
         VecScale*Ug_SWOT_HFRgrid,VecScale*Vg_SWOT_HFRgrid,0,'k')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
CB = colorbar; CB.Label.String = 'SWOT: |u_g| at HFR grid (m s^{-1})'; CB.Label.FontSize = 20; clim([0 1]*1)
set(gca,'FontSize',16)
title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))])

[x_click,y_click] = ginput(1);
[x_click,y_click] = m_xy2ll(x_click,y_click);
[yi_click,xi_click] = find(abs([NORCAL.HFR.LON    + 1i*NORCAL.HFR.LAT]    - [x_click + 1i*y_click]) == ...
                       min(abs([NORCAL.HFR.LON(:) + 1i*NORCAL.HFR.LAT(:)] - [x_click + 1i*y_click])));

m_plot(NORCAL.HFR.LON(yi_click,xi_click),...
       NORCAL.HFR.LAT(yi_click,xi_click),'*w')

% % Before further interpolating all SWOT geostrophic currents to the HFR
% grid, first just compare the average of u_v,v_g within some radius:
warning(['This only works because each pass is identicial; ' ...
         'this script will not work for multi-pass analysis.'])
[yi_click_swot,xi_click_swot] = ...
    find(abs([NORCAL.SWOT.lon{1}    + 1i*NORCAL.SWOT.lat{1}]    - [x_click + 1i*y_click]) == ...
     min(abs([NORCAL.SWOT.lon{1}(:) + 1i*NORCAL.SWOT.lat{1}(:)] - [x_click + 1i*y_click])));

figure('Color','w')
tiledlayout(4,1,'TileSpacing','compact')
NORCAL_HFR_uv_magnitude = squeeze(sqrt([NORCAL.HFR.u(yi_click,xi_click,:)].^2 + [NORCAL.HFR.v(yi_click,xi_click,:)].^2));

AX1 = nexttile;
% plot(NORCAL.HFR.time/24 + T0 , squeeze(NORCAL.HFR.u(yi_click,xi_click,:)) , '.-'); hold on
plot(NORCAL.HFR.time/24 + T0 , squeeze(NORCAL.HFR.u(yi_click,xi_click,:))./[NORCAL_HFR_uv_magnitude>Plot_uvMagnitude_Cuttoff] , '.-'); hold on
plot(NORCAL.HFR.time/24 + T0 , squeeze(U_rot_filtered(yi_click,xi_click,:)) , '.-r'); hold on
% plot(NORCAL.HFR.time/24 + T0 , squeeze(U_rot_filtered(yi_click,xi_click,:) + U_div_filtered(yi_click,xi_click,:)) , '.-r'); hold on%$
plot(NORCAL.SWOT.mean_time , squeeze(U_geostr(yi_click_swot,xi_click_swot,:)) , '.-k')
plot(NORCAL.SWOT.mean_time , squeeze(U_cyclogeostr_Nit(yi_click_swot,xi_click_swot,:)) , '.-c')
ylabel('u (m/s)'); datetick
legend('Unfiltered HFR','Filtered HFR','SWOT geostrophic',['SWOT cyclogeostrophic (' num2str(CG_Iteration) ' iteration)'])
title({[num2str(NORCAL.HFR.LAT(yi_click,xi_click)) '\circ lat' ] ; ...
       [num2str(NORCAL.HFR.LON(yi_click,xi_click)) '\circ lon' ]})

AX2 = nexttile;
% plot(NORCAL.HFR.time/24 + T0 , squeeze(NORCAL.HFR.v(yi_click,xi_click,:)) , '.-','HandleVisibility','off'); hold on
plot(NORCAL.HFR.time/24 + T0 , squeeze(NORCAL.HFR.v(yi_click,xi_click,:))./[NORCAL_HFR_uv_magnitude>Plot_uvMagnitude_Cuttoff] , '.-','HandleVisibility','off'); hold on
plot(NORCAL.HFR.time/24 + T0 , squeeze(V_rot_filtered(yi_click,xi_click,:)) , '.-r'); hold on
plot(NORCAL.SWOT.mean_time , squeeze(V_geostr(yi_click_swot,xi_click_swot,:)) , '.-k')
plot(NORCAL.SWOT.mean_time , squeeze(V_cyclogeostr_Nit(yi_click_swot,xi_click_swot,:)) , '.-c')
ylabel('v'); datetick

AX3 = nexttile;
plot(NORCAL.HFR.time/24 + T0 , NORCAL_HFR_uv_magnitude , '.-','HandleVisibility','off'); hold on
plot(NORCAL.HFR.time/24 + T0 , squeeze(sqrt([U_rot_filtered(yi_click,xi_click,:)].^2 + ...
                                            [V_rot_filtered(yi_click,xi_click,:)].^2)) , '.-r'); hold on
plot(NORCAL.SWOT.mean_time , squeeze(sqrt([U_geostr(yi_click_swot,xi_click_swot,:)].^2 + ...
                                          [V_geostr(yi_click_swot,xi_click_swot,:)].^2)) , '.-k')
plot(NORCAL.SWOT.mean_time , squeeze(sqrt([U_cyclogeostr_Nit(yi_click_swot,xi_click_swot,:)].^2 + ...
                                          [V_cyclogeostr_Nit(yi_click_swot,xi_click_swot,:)].^2)) , '.-c')
ylabel('$\sqrt(u^2 + v^2)$','Interpreter','LaTeX'); datetick

AX4 = nexttile;
% plot(NORCAL.HFR.time/24 + T0 , squeeze(angle([NORCAL.HFR.u(yi_click,xi_click,:)] + ...
%                                              [NORCAL.HFR.v(yi_click,xi_click,:)]*1i))*180/pi , '.','HandleVisibility','off'); hold on
plot(NORCAL.HFR.time/24 + T0 , squeeze(angle([NORCAL.HFR.u(yi_click,xi_click,:)] + ...
                                             [NORCAL.HFR.v(yi_click,xi_click,:)]*1i))*[180/pi]./[NORCAL_HFR_uv_magnitude>Plot_uvMagnitude_Cuttoff] , '.','HandleVisibility','off'); hold on
plot(NORCAL.HFR.time/24 + T0 , squeeze(angle([U_rot_filtered(yi_click,xi_click,:)] + ...
                                             [V_rot_filtered(yi_click,xi_click,:)]*1i))*180/pi , '.r'); hold on
plot(NORCAL.SWOT.mean_time , squeeze(angle([U_geostr(yi_click_swot,xi_click_swot,:)] + ...
                                           [V_geostr(yi_click_swot,xi_click_swot,:)]*1i))*180/pi , '.-k')
plot(NORCAL.SWOT.mean_time , squeeze(angle([U_cyclogeostr_Nit(yi_click_swot,xi_click_swot,:)] + ...
                                           [V_cyclogeostr_Nit(yi_click_swot,xi_click_swot,:)]*1i))*180/pi , '.-c')
% AX4 = nexttile; % UNWRAPPED VERSION:
% % plot(NORCAL.HFR.time/24 + T0 , squeeze(unwrap(angle([NORCAL.HFR.u(yi_click,xi_click,:)] + ...
% %                                                     [NORCAL.HFR.v(yi_click,xi_click,:)]*1i)))*180/pi , '.-'); hold on
% plot(NORCAL.HFR.time/24 + T0 , squeeze(unwrap(angle([U_rot_filtered(yi_click,xi_click,:)] + ...
%                                                     [V_rot_filtered(yi_click,xi_click,:)]*1i)))*180/pi , '.-'); hold on
% plot(NORCAL.SWOT.mean_time , squeeze(unwrap(angle([U_geostr(yi_click_swot,xi_click_swot,:)] + ...
%                                                   [V_geostr(yi_click_swot,xi_click_swot,:)]*1i)))*180/pi , '.-')

yticks(-180:45:180)
% set(gca,'ylim',[-180 180])
ylabel('angle(u + iv)'); datetick % xtickformat('MM-dd')

set(gcf,'Position',[-1515         118        1243         850])
set([AX1 AX2 AX3 AX4],'Color',[1 1 1]*0.7)

linkaxes([AX1 AX2 AX3 AX4],'x')
linkaxes([AX1 AX2],'y')
set(gca,'ylim',[-180 180]) % for some reason this needs to be after the axis link

AX1.FontSize = 14;
AX2.FontSize = 14;
AX3.FontSize = 14;
AX4.FontSize = 14;

% %% Plot filtered and unfiltered spectra (OPTIONAL)
% 
% figure
% subplot(211)
% [S_uu,f_vec] = nanspectrum1(NORCAL.HFR.u(yi_click,xi_click,nn:end), 1/24, 'day', SEGS, '.-', true, 0, 'hanning'); hold on
% [S_lulu,~]   = nanspectrum1(L_U(yi_click,xi_click,nn:end), 1/24, 'day', SEGS, '.-', true, 0, 'hanning');
% legend('S_{u}','','S_{LP(u)}')
% subplot(212)
% semilogx(f_vec,S_lulu./S_uu,'k.-')
% ylabel('S_{LP(u)}/S_{u}')

%% RMSD and correlation between geostrophic velocity to HFR velocity

% HFR at/near SWOT times
U_rot_filtered_SWOTtimes = nan(size(Ug_SWOT_HFRgrid_all));
V_rot_filtered_SWOTtimes = nan(size(Ug_SWOT_HFRgrid_all));
tic
for ti = 1:length(NORCAL.SWOT.mean_time)
    ti_HFR_closest = dsearchn(NORCAL.HFR.time/24 + T0,NORCAL.SWOT.mean_time(ti));
    U_rot_filtered_SWOTtimes(:,:,ti) = U_rot_filtered(:,:,ti_HFR_closest);
    V_rot_filtered_SWOTtimes(:,:,ti) = V_rot_filtered(:,:,ti_HFR_closest);
end
toc

%
% Correlation
CUTOFF = 0.8; % minimum fraction of good data required to perform calculation
U_SWOTg_HFR_Corr_R = nan(size(NORCAL.HFR.LON));
U_SWOTg_HFR_Corr_P = nan(size(NORCAL.HFR.LON));
V_SWOTg_HFR_Corr_R = nan(size(NORCAL.HFR.LON));
V_SWOTg_HFR_Corr_P = nan(size(NORCAL.HFR.LON));
RMSD_velmag_SWOTg_HFR = nan(size(NORCAL.HFR.LON));
U_SWOTcg_HFR_Corr_R = nan(size(NORCAL.HFR.LON));
U_SWOTcg_HFR_Corr_P = nan(size(NORCAL.HFR.LON));
V_SWOTcg_HFR_Corr_R = nan(size(NORCAL.HFR.LON));
V_SWOTcg_HFR_Corr_P = nan(size(NORCAL.HFR.LON));
RMSD_velmag_SWOTcg_HFR = nan(size(NORCAL.HFR.LON));
W_SWOTg_HFR_Corr_R = nan(size(NORCAL.HFR.LON)); % Complex (no p-values possible with them)
W_SWOTg_HFR_Corr_angle = nan(size(NORCAL.HFR.LON));
W_SWOTcg_HFR_Corr_R = nan(size(NORCAL.HFR.LON)); % Complex (no p-values possible with them)
W_SWOTcg_HFR_Corr_angle = nan(size(NORCAL.HFR.LON));
tic
for ii = 1:size(Ug_SWOT_HFRgrid_all,1)
    for jj = 1:size(Ug_SWOT_HFRgrid_all,2)
        % Geostrophy
        if sum(isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)))/numel(Ug_SWOT_HFRgrid_all(ii,jj,:)) > CUTOFF & ...
           sum(isfinite(U_rot_filtered_SWOTtimes(ii,jj,:)))/numel(U_rot_filtered_SWOTtimes(ii,jj,:)) > CUTOFF
            [RR,PP] = corrcoef(Ug_SWOT_HFRgrid_all(     ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:))), ...
                               U_rot_filtered_SWOTtimes(ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:))));
            U_SWOTg_HFR_Corr_R(ii,jj) = RR(1,2);
            U_SWOTg_HFR_Corr_P(ii,jj) = PP(1,2);
            [RR,PP] = corrcoef(Vg_SWOT_HFRgrid_all(     ii,jj,isfinite(Vg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))), ...
                               V_rot_filtered_SWOTtimes(ii,jj,isfinite(Vg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))));
            V_SWOTg_HFR_Corr_R(ii,jj) = RR(1,2);
            V_SWOTg_HFR_Corr_P(ii,jj) = PP(1,2);
            [RR] = corrcoef(   Ug_SWOT_HFRgrid_all(  ii,jj,isfinite(   Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))) + ...
                            1i*Vg_SWOT_HFRgrid_all(  ii,jj,isfinite(   Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))), ...
                               U_rot_filtered_SWOTtimes(ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))) + ...
                            1i*V_rot_filtered_SWOTtimes(ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))));
            W_SWOTg_HFR_Corr_R(ii,jj) = abs(RR(1,2));
            W_SWOTg_HFR_Corr_angle(ii,jj) = angle(RR(1,2))*180/pi;
            RMSD_velmag_SWOTg_HFR(ii,jj) = rms( [Ug_SWOT_HFRgrid_all(ii,jj,:) + 1i*Vg_SWOT_HFRgrid_all(ii,jj,:)] - ...
                                                [U_rot_filtered_SWOTtimes(ii,jj,:) + 1i*V_rot_filtered_SWOTtimes(ii,jj,:)],'omitnan');
        else
        end
        % Cyclogeostrophy
        if sum(isfinite(Ucg_SWOT_HFRgrid_all(ii,jj,:)))/numel(Ucg_SWOT_HFRgrid_all(ii,jj,:)) > CUTOFF & ...
           sum(isfinite(U_rot_filtered_SWOTtimes(ii,jj,:)))/numel(U_rot_filtered_SWOTtimes(ii,jj,:)) > CUTOFF
            [RR,PP] = corrcoef(Ucg_SWOT_HFRgrid_all(     ii,jj,isfinite(Ucg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:))), ...
                               U_rot_filtered_SWOTtimes( ii,jj,isfinite(Ucg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:))));
            U_SWOTcg_HFR_Corr_R(ii,jj) = RR(1,2);
            U_SWOTcg_HFR_Corr_P(ii,jj) = PP(1,2);
            [RR,PP] = corrcoef(Vcg_SWOT_HFRgrid_all(     ii,jj,isfinite(Vcg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))), ...
                               V_rot_filtered_SWOTtimes( ii,jj,isfinite(Vcg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))));
            V_SWOTcg_HFR_Corr_R(ii,jj) = RR(1,2);
            V_SWOTcg_HFR_Corr_P(ii,jj) = PP(1,2);
            [RR] = corrcoef(   Ucg_SWOT_HFRgrid_all(    ii,jj,isfinite( Ucg_SWOT_HFRgrid_all(ii,jj,:)+Vcg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))) + ...
                            1i*Vcg_SWOT_HFRgrid_all(    ii,jj,isfinite( Ucg_SWOT_HFRgrid_all(ii,jj,:)+Vcg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))), ...
                               U_rot_filtered_SWOTtimes(ii,jj,isfinite( Ucg_SWOT_HFRgrid_all(ii,jj,:)+Vcg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))) + ...
                            1i*V_rot_filtered_SWOTtimes(ii,jj,isfinite( Ucg_SWOT_HFRgrid_all(ii,jj,:)+Vcg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))));
            W_SWOTcg_HFR_Corr_R(ii,jj) = abs(RR(1,2));
            W_SWOTcg_HFR_Corr_angle(ii,jj) = angle(RR(1,2))*180/pi;
            RMSD_velmag_SWOTcg_HFR(ii,jj) = rms([Ucg_SWOT_HFRgrid_all(    ii,jj,:) + 1i*Vcg_SWOT_HFRgrid_all(    ii,jj,:)] - ...
                                                [U_rot_filtered_SWOTtimes(ii,jj,:) + 1i*V_rot_filtered_SWOTtimes(ii,jj,:)],'omitnan');
        else
        end
    end
    disp(100*ii/size(Ug_SWOT_HFRgrid_all,1))
end
toc

%% Save key variables to compare to the results of other approaches

error('Already done.')

save(['./CorrCoef_SWOT_HFRcartesian_quadfit_' num2str([2*PIXELS + 1]*2) 'kmBox.mat'],...
    'U_SWOTg_HFR_Corr_R',...
    'U_SWOTg_HFR_Corr_P',...
    'V_SWOTg_HFR_Corr_R',...
    'V_SWOTg_HFR_Corr_P',...
    'RMSD_velmag_SWOTg_HFR',...
    'U_SWOTcg_HFR_Corr_R',...
    'U_SWOTcg_HFR_Corr_P',...
    'V_SWOTcg_HFR_Corr_R',...
    'V_SWOTcg_HFR_Corr_P',...
    'RMSD_velmag_SWOTcg_HFR',...
    'PIXELS')

%% Map correlation

close all
figure('Color','w')
tiledlayout(2,3,'TileSpacing','tight')

P_threshold = 0.05;
% P_threshold = 0.1;

CorrCoefColorLim = [-1 1];

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  RMSD_velmag_SWOTg_HFR);
m_text(-122.5,41.5,'(a)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
% CB = colorbar;
clim([0 1]*0.5)
% TITLE_TEXT = 'RMSE (Geostr. vs. HFR) (m/s)';
TITLE_TEXT = 'RMSE (m s$^{-1}$)';
% CB.Label.FontSize = 16;
title(TITLE_TEXT,'Interpreter','latex')
ylabel('Geostrophic vs. HFR','Interpreter','latex')

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  U_SWOTg_HFR_Corr_R .* [U_SWOTg_HFR_Corr_P<P_threshold]./[U_SWOTg_HFR_Corr_P<P_threshold]);
m_text(-122.5,41.5,'(b)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
% CB = colorbar;
clim(CorrCoefColorLim)
% TITLE_TEXT = {['U (Geostr. vs. HFR): corr. coef. R'];['(where P<' num2str(P_threshold) ')']};
TITLE_TEXT = {['$C_{uu}$ where $P<$ ' num2str(P_threshold) '']};
% CB.Label.FontSize = 16;
title(TITLE_TEXT,'Interpreter','latex')

AX3 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  V_SWOTg_HFR_Corr_R .* [V_SWOTg_HFR_Corr_P<P_threshold]./[V_SWOTg_HFR_Corr_P<P_threshold]);
m_text(-122.5,41.5,'(c)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
% CB = colorbar;
clim(CorrCoefColorLim)
% TITLE_TEXT = {['V (Geostr. vs. HFR): corr. coef. R'];['(where P<' num2str(P_threshold) ')']};
TITLE_TEXT = {['$C_{vv}$ where $P<$ ' num2str(P_threshold) '']};
% CB.Label.FontSize = 16;
title(TITLE_TEXT,'Interpreter','latex')

% % % % % % % % % % % 

AX4 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  RMSD_velmag_SWOTcg_HFR);
m_text(-122.5,41.5,'(d)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([0 1]*0.5)
% TITLE_TEXT = 'RMSE (CG vs. HFR) (m/s)';
CB.Label.String = 'm s$^{-1}$';
CB.Label.Interpreter = 'latex';
CB.Label.FontSize = 16; CB.Location = 'southoutside';
% title(TITLE_TEXT)
ylabel('Cyclogeostrophic vs. HFR','Interpreter','latex')

AX5 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  U_SWOTcg_HFR_Corr_R .* [U_SWOTcg_HFR_Corr_P<P_threshold]./[U_SWOTcg_HFR_Corr_P<P_threshold]);
m_text(-122.5,41.5,'(e)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim(CorrCoefColorLim)
% TITLE_TEXT = {['U (CG vs. HFR): corr. coef. R'];['(where P<' num2str(P_threshold) ')']};
CB.Label.String = 'Correlation Coefficient';
CB.Label.Interpreter = 'latex';
CB.Label.FontSize = 16; CB.Location = 'southoutside';
% title(TITLE_TEXT)

AX6 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  V_SWOTcg_HFR_Corr_R .* [V_SWOTcg_HFR_Corr_P<P_threshold]./[V_SWOTcg_HFR_Corr_P<P_threshold]);
m_text(-122.5,41.5,'(f)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim(CorrCoefColorLim)
% TITLE_TEXT = {['V (CG vs. HFR): corr. coef. R'];['(where P<' num2str(P_threshold) ')']};
CB.Label.String = 'Correlation Coefficient';
CB.Label.Interpreter = 'latex';
CB.Label.FontSize = 16; CB.Location = 'southoutside';
% title(TITLE_TEXT)

linkaxes([AX1 AX2 AX3 AX4 AX5 AX6],'xy')

colormap(turbo)
colormap([AX1],flip(pink))
colormap([AX4],flip(pink))
set([AX1 AX2 AX3 AX4 AX5 AX6],'Color',[1 1 1]*0.7)
% set(gcf,'Position',[-1158 81 1159 896])
set(gcf,'Position',[-951 83 952 894])


% % % COMPLEX CORRELATION AND ANGLE
figure('Color','w')
tiledlayout(2,3,'TileSpacing','tight')

P_threshold = 0.05;
% P_threshold = 0.1;

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  RMSD_velmag_SWOTg_HFR);
m_text(-122.5,41.63,'(a)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
% CB = colorbar;
clim([0 1]*0.5)
% TITLE_TEXT = 'RMSE (Geostr. vs. HFR) (m/s)';
TITLE_TEXT = 'RMSE (m s$^{-1}$)';
% CB.Label.FontSize = 16;
title(TITLE_TEXT,'Interpreter','latex')
ylabel('Geostrophic vs. HFR','Interpreter','latex')

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  W_SWOTg_HFR_Corr_R);% .* [U_SWOTg_HFR_Corr_P<P_threshold]./[U_SWOTg_HFR_Corr_P<P_threshold]);
m_text(-122.5,41.63,'(b)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
% CB = colorbar;
clim([0 1])
% TITLE_TEXT = {['U (Geostr. vs. HFR): corr. coef. R'];['(where P<' num2str(P_threshold) ')']};
TITLE_TEXT = {['$C_{ww}$']};
% CB.Label.FontSize = 16;
title(TITLE_TEXT,'Interpreter','latex')

AX3 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  W_SWOTg_HFR_Corr_angle.*[W_SWOTg_HFR_Corr_R>=0.5]./[W_SWOTg_HFR_Corr_R>=0.5]);% .* [V_SWOTg_HFR_Corr_P<P_threshold]./[V_SWOTg_HFR_Corr_P<P_threshold]);
m_text(-122.5,41.63,'(c)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
% CB = colorbar;
clim([-1 1]*180)
% TITLE_TEXT = {['V (Geostr. vs. HFR): corr. coef. R'];['(where P<' num2str(P_threshold) ')']};
TITLE_TEXT = {['$\theta_{ww}$']};
% CB.Label.FontSize = 16;
title(TITLE_TEXT,'Interpreter','latex')

% % % % % % % % % % % 

AX4 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  RMSD_velmag_SWOTcg_HFR);
m_text(-122.5,41.63,'(d)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([0 1]*0.5)
% TITLE_TEXT = 'RMSE (CG vs. HFR) (m/s)';
CB.Label.String = 'm s$^{-1}$';
CB.Label.Interpreter = 'latex';
CB.Label.FontSize = 16; CB.Location = 'southoutside';
% title(TITLE_TEXT)
ylabel('Cyclogeostrophic vs. HFR','Interpreter','latex')

AX5 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  W_SWOTcg_HFR_Corr_R);% .* [U_SWOTcg_HFR_Corr_P<P_threshold]./[U_SWOTcg_HFR_Corr_P<P_threshold]);
m_text(-122.5,41.63,'(e)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([0 1])
% TITLE_TEXT = {['U (CG vs. HFR): corr. coef. R'];['(where P<' num2str(P_threshold) ')']};
CB.Label.String = 'Correlation Coefficient';
CB.Label.Interpreter = 'latex';
CB.Label.FontSize = 16; CB.Location = 'southoutside';
% title(TITLE_TEXT)

AX6 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  W_SWOTcg_HFR_Corr_angle.*[W_SWOTcg_HFR_Corr_R>=0.5]./[W_SWOTcg_HFR_Corr_R>=0.5]);% .* [V_SWOTcg_HFR_Corr_P<P_threshold]./[V_SWOTcg_HFR_Corr_P<P_threshold]);
m_text(-122.5,41.63,'(f)','Interpreter','latex','HorizontalAlignment','center','FontSize',20)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar;
clim([-1 1]*180)
% TITLE_TEXT = {['V (CG vs. HFR): corr. coef. R'];['(where P<' num2str(P_threshold) ')']};
CB.Label.String = 'Angle (degrees)';
CB.Label.Interpreter = 'latex';
CB.Label.FontSize = 16; CB.Location = 'southoutside';
CB.Ticks = [-180 -90 0 90 180];
% title(TITLE_TEXT)

linkaxes([AX1 AX2 AX3 AX4 AX5 AX6],'xy')

colormap(turbo)
colormap([AX1],flip(pink))
colormap([AX4],flip(pink))
colormap([AX3],cmocean('phase'))
colormap([AX6],cmocean('phase'))
set([AX1 AX2 AX3 AX4 AX5 AX6],'Color',[1 1 1]*0.7)
% set(gcf,'Position',[-1158 81 1159 896])
set(gcf,'Position',[-951 83 952 894])

INSET_POSITION = [0.2216    0.7926    0.1051    0.0717];
axes('Position',INSET_POSITION)
histogram(RMSD_velmag_SWOTg_HFR,[0:0.06:0.6],'Normalization','pdf','FaceColor',[1 1 1]*0.2); yticklabels([])
axes('Position',INSET_POSITION + 1*[0.2815 0 0 0])
histogram(W_SWOTg_HFR_Corr_R,[0:0.05:1],'Normalization','pdf','FaceColor',[1 1 1]*0.2); yticklabels([])
axes('Position',INSET_POSITION + 2*[0.2815 0 0 0])
histogram(W_SWOTg_HFR_Corr_angle.*[W_SWOTg_HFR_Corr_R>=0.5]./[W_SWOTg_HFR_Corr_R>=0.5],[-180:10:180],'Normalization','pdf','FaceColor',[1 1 1]*0.2); yticklabels([])
xticks([-180 -90 0 90 180])
axes('Position',INSET_POSITION + [0 -0.4204 0 0])
histogram(RMSD_velmag_SWOTcg_HFR,[0:0.06:0.6],'Normalization','pdf','FaceColor',[1 1 1]*0.2); yticklabels([])
axes('Position',INSET_POSITION + 1*[0.2815 0 0 0] + [0 -0.4204 0 0])
histogram(W_SWOTcg_HFR_Corr_R,[0:0.05:1],'Normalization','pdf','FaceColor',[1 1 1]*0.2); yticklabels([])
axes('Position',INSET_POSITION + 2*[0.2815 0 0 0] + [0 -0.4204 0 0])
histogram(W_SWOTcg_HFR_Corr_angle.*[W_SWOTcg_HFR_Corr_R>=0.5]./[W_SWOTcg_HFR_Corr_R>=0.5],[-180:10:180],'Normalization','pdf','FaceColor',[1 1 1]*0.2); yticklabels([])
xticks([-180 -90 0 90 180])
% % % ^^^^^^^ COMPLEX CORRELATION AND ANGLE


% figure('Color','w')
% plot(squeeze(U_rot_filtered_SWOTtimes(yi_click,xi_click,:)),...
%      squeeze(V_rot_filtered_SWOTtimes(yi_click,xi_click,:)), '.-'); hold on
% plot(squeeze(Ug_SWOT_HFRgrid_all(yi_click,xi_click,:)),...
%      squeeze(Vg_SWOT_HFRgrid_all(yi_click,xi_click,:)), '.-'); hold on

% error
% % % % % % % % % % % 
% % % % % % % % % % % 
% % % % % % % % % % % 

figure('Color','w')
tiledlayout(1,3,'TileSpacing','compact')

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  RMSD_velmag_SWOTg_HFR - RMSD_velmag_SWOTcg_HFR);
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1]*0.2)
CB.Label.String = {'\Delta RMSE (Geostr. vs. HFR) (m/s)';...
                   '\Delta<0 = Geostr. is better  ---  \Delta>0 = CycloGeostr. is better'}; CB.Label.FontSize = 16;

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  [U_SWOTg_HFR_Corr_R .* [U_SWOTg_HFR_Corr_P<P_threshold]./[U_SWOTg_HFR_Corr_P<P_threshold]] - ...
                  [U_SWOTcg_HFR_Corr_R .* [U_SWOTcg_HFR_Corr_P<P_threshold]./[U_SWOTcg_HFR_Corr_P<P_threshold]]);
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1]*0.2)
CB.Label.String = {['\Delta U (Geostr. vs. HFR): corr. coef. R (where P<' num2str(P_threshold) ')'];...
                   '\Delta>0 = Geostr. is better  ---  \Delta<0 = CycloGeostr. is better'}; CB.Label.FontSize = 16;

AX3 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  [V_SWOTg_HFR_Corr_R .* [V_SWOTg_HFR_Corr_P<P_threshold]./[V_SWOTg_HFR_Corr_P<P_threshold]] - ...
                  [V_SWOTcg_HFR_Corr_R .* [V_SWOTcg_HFR_Corr_P<P_threshold]./[V_SWOTcg_HFR_Corr_P<P_threshold]]);
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1]*0.2)
CB.Label.String = {['\Delta V (Geostr. vs. HFR): corr. coef. R (where P<' num2str(P_threshold) ')'];...
                   '\Delta>0 = Geostr. is better  ---  \Delta<0 = CycloGeostr. is better'}; CB.Label.FontSize = 16;
linkaxes([AX1 AX2 AX3],'xy')

colormap(bwr)
set([AX1 AX2 AX3],'Color',[1 1 1]*0.7)

figure(3)

[x_click,y_click] = ginput(1);
[x_click,y_click] = m_xy2ll(x_click,y_click);
[yi_click,xi_click] = find(abs([NORCAL.HFR.LON    + 1i*NORCAL.HFR.LAT]    - [x_click + 1i*y_click]) == ...
                       min(abs([NORCAL.HFR.LON(:) + 1i*NORCAL.HFR.LAT(:)] - [x_click + 1i*y_click])));


AX1; m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'*w');
AX2; m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'*w');
AX3; m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'*w');



% % % % % % % % % % % 
% % % % % % % % % % % 
% % % % % % % % % % % 

figure('Color','w')
subplot(211)
plot(NORCAL.SWOT.mean_time,squeeze(U_rot_filtered_SWOTtimes(yi_click,xi_click,:)), 'k.-'); hold on
plot(NORCAL.SWOT.mean_time,squeeze(Ug_SWOT_HFRgrid_all(yi_click,xi_click,:)), 'r.-')
plot(NORCAL.SWOT.mean_time,squeeze(Ucg_SWOT_HFRgrid_all(yi_click,xi_click,:)), 'c.-')
title({['CorrCoef(U_{HFR} , U_g) = ' num2str(   U_SWOTg_HFR_Corr_R( yi_click,xi_click))];...
       ['CorrCoef(U_{HFR} , U_{cg}) = ' num2str(U_SWOTcg_HFR_Corr_R(yi_click,xi_click))]})
legend('HFR U rotational filtered at SWOT times', ...
       'SWOT U geostr. SWOT on HFR grid', ...
       'SWOT U cyclogeostr. SWOT on HFR grid')
datetick
subplot(212)
plot(NORCAL.SWOT.mean_time,squeeze(V_rot_filtered_SWOTtimes(yi_click,xi_click,:)), 'k.-'); hold on
plot(NORCAL.SWOT.mean_time,squeeze(Vg_SWOT_HFRgrid_all(yi_click,xi_click,:)), 'r.-')
plot(NORCAL.SWOT.mean_time,squeeze(Vcg_SWOT_HFRgrid_all(yi_click,xi_click,:)), 'c.-')
title({['CorrCoef(V_{HFR} , V_g) = ' num2str(   V_SWOTg_HFR_Corr_R( yi_click,xi_click))];...
       ['CorrCoef(V_{HFR} , V_{cg}) = ' num2str(V_SWOTcg_HFR_Corr_R(yi_click,xi_click))]})
legend('HFR V rotational filtered at SWOT times', ...
       'SWOT V geostr. SWOT on HFR grid', ...
       'SWOT V cyclogeostr. SWOT on HFR grid')
datetick



%% Histogram of CorrCoef:

close all

figure('Color','w')
HIST_BINS = [-1:0.01:1];
subplot(311)
histogram(U_SWOTg_HFR_Corr_R  .* [U_SWOTg_HFR_Corr_P<P_threshold] ./[U_SWOTg_HFR_Corr_P<P_threshold], HIST_BINS);hold on
histogram(U_SWOTcg_HFR_Corr_R .* [U_SWOTcg_HFR_Corr_P<P_threshold]./[U_SWOTcg_HFR_Corr_P<P_threshold],HIST_BINS);
subplot(312)
histogram(V_SWOTg_HFR_Corr_R  .* [V_SWOTg_HFR_Corr_P<P_threshold] ./[V_SWOTg_HFR_Corr_P<P_threshold], HIST_BINS);hold on
histogram(V_SWOTcg_HFR_Corr_R .* [V_SWOTcg_HFR_Corr_P<P_threshold]./[V_SWOTcg_HFR_Corr_P<P_threshold],HIST_BINS);
subplot(313)
histogram(RMSD_velmag_SWOTg_HFR .* [U_SWOTg_HFR_Corr_P<P_threshold] ./[U_SWOTg_HFR_Corr_P<P_threshold] ...
                                .* [V_SWOTg_HFR_Corr_P<P_threshold] ./[V_SWOTg_HFR_Corr_P<P_threshold], [0:0.01:1]);hold on
histogram(RMSD_velmag_SWOTcg_HFR.* [U_SWOTg_HFR_Corr_P<P_threshold] ./[U_SWOTg_HFR_Corr_P<P_threshold] ...
                                .* [V_SWOTg_HFR_Corr_P<P_threshold] ./[V_SWOTg_HFR_Corr_P<P_threshold], [0:0.01:1]);
legend('RMSD: G vs. HFR','RMSD: CG vs. HFR')



figure('Color','w')
HIST_BINS = [-1:0.01:1];
subplot(311)
histogram(U_SWOTg_HFR_Corr_R  , HIST_BINS);hold on
histogram(U_SWOTcg_HFR_Corr_R , HIST_BINS);
legend('C_{uu}: G vs. HFR','C_{uu}: CG vs. HFR')
subplot(312)
histogram(V_SWOTg_HFR_Corr_R  , HIST_BINS);hold on
histogram(V_SWOTcg_HFR_Corr_R , HIST_BINS);
legend('C_{vv}: G vs. HFR','C_{vv}: CG vs. HFR')
subplot(313)
histogram(RMSD_velmag_SWOTg_HFR  , [0:0.01:1]);hold on
histogram(RMSD_velmag_SWOTcg_HFR , [0:0.01:1]);
legend('RMSD: G vs. HFR','RMSD: CG vs. HFR')

%%

close all

figure('Color','w')
HIST_BINS = [0:0.01:0.5];
subplot(211)
histogram(rms(Ug_SWOT_HFRgrid_all + 1i*Vg_SWOT_HFRgrid_all,3,'omitnan'),   HIST_BINS, 'Normalization','pdf'); hold on
histogram(rms(U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes,3,'omitnan'), HIST_BINS, 'Normalization','pdf')
histogram(rms([Ug_SWOT_HFRgrid_all + 1i*Vg_SWOT_HFRgrid_all] - ...
              [U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes],3,'omitnan'), HIST_BINS, 'Normalization','pdf')
hold on
% histogram(RMSD_velmag_SWOTg_HFR, HIST_BINS, 'Normalization','pdf')
% legend('|u+iv|_{g}','|u+iv|_{HFR}','|difference|')
legend('RMS (u+iv)_{g}','RMS (u+iv)_{HFR}','RMS (difference)')
title({'RMS';'Geostrophy'})

subplot(212)
histogram(rms(Ucg_SWOT_HFRgrid_all + 1i*Vcg_SWOT_HFRgrid_all,3,'omitnan'), HIST_BINS, 'Normalization','pdf'); hold on
histogram(rms(U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes,3,'omitnan'), HIST_BINS, 'Normalization','pdf')
histogram(rms([Ucg_SWOT_HFRgrid_all + 1i*Vcg_SWOT_HFRgrid_all] - ...
              [U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes],3,'omitnan'), HIST_BINS, 'Normalization','pdf')
% legend('|u+iv|_{cg}','|u+iv|_{HFR}','|difference|')
legend('RMS (u+iv)_{cg}','RMS (u+iv)_{HFR}','RMS (difference)')
title('Cyclogeostrophy')



figure('Color','w')
subplot(211)
histogram(std(Ug_SWOT_HFRgrid_all + 1i*Vg_SWOT_HFRgrid_all,0,3,'omitnan'),   HIST_BINS, 'Normalization','pdf'); hold on
histogram(std(U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes,0,3,'omitnan'), HIST_BINS, 'Normalization','pdf')
histogram(std([Ug_SWOT_HFRgrid_all + 1i*Vg_SWOT_HFRgrid_all] - ...
              [U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes],0,3,'omitnan'), HIST_BINS, 'Normalization','pdf')
hold on
% legend('|u+iv|_{g}','|u+iv|_{HFR}','|difference|')
legend('STD (u+iv)_{g}','STD (u+iv)_{HFR}','STD (difference)')
title({'STD';'Geostrophy'})

subplot(212)
histogram(std(Ucg_SWOT_HFRgrid_all + 1i*Vcg_SWOT_HFRgrid_all,0,3,'omitnan'), HIST_BINS, 'Normalization','pdf'); hold on
histogram(std(U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes,0,3,'omitnan'), HIST_BINS, 'Normalization','pdf')
histogram(std([Ucg_SWOT_HFRgrid_all + 1i*Vcg_SWOT_HFRgrid_all] - ...
              [U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes],0,3,'omitnan'), HIST_BINS, 'Normalization','pdf')
% legend('|u+iv|_{cg}','|u+iv|_{HFR}','|difference|')
legend('STD (u+iv)_{cg}','STD (u+iv)_{HFR}','STD (difference)')
title('Cyclogeostrophy')

% [Ucg_SWOT_HFRgrid_all(    ii,jj,:) + 1i*Vcg_SWOT_HFRgrid_all(    ii,jj,:)] - ...
% [U_rot_filtered_SWOTtimes(ii,jj,:) + 1i*V_rot_filtered_SWOTtimes(ii,jj,:)]

figure('Color','y')
subplot(211)
histogram(rms([Ug_SWOT_HFRgrid_all + 1i*Vg_SWOT_HFRgrid_all] - ...
              [U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes],3,'omitnan'), HIST_BINS, 'Normalization','pdf');hold on
histogram(RMSD_velmag_SWOTg_HFR, HIST_BINS, 'Normalization','pdf')
subplot(212)
histogram(rms([Ucg_SWOT_HFRgrid_all + 1i*Vcg_SWOT_HFRgrid_all] - ...
              [U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes],3,'omitnan'), HIST_BINS, 'Normalization','pdf');hold on
histogram(RMSD_velmag_SWOTcg_HFR, HIST_BINS, 'Normalization','pdf')

%% Percentiles:
disp(['%%%%%%%%%%%%%%%%%%%%%%% CORRELATION COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%'])
disp(['        5          10           25          50          75          90           95%ile'])
disp(['Ug:    ' num2str(prctile(U_SWOTg_HFR_Corr_R(:), [5    10   25   50   75   90   95]))])
disp(['Ucg:   ' num2str(prctile(U_SWOTcg_HFR_Corr_R(:),[5    10   25   50   75   90   95]))])
disp(['Vg:    ' num2str(prctile(V_SWOTg_HFR_Corr_R(:), [5    10   25   50   75   90   95]))])
disp(['Vcg:   ' num2str(prctile(V_SWOTcg_HFR_Corr_R(:),[5    10   25   50   75   90   95]))])

disp(['%%%%%%%%%%%%%%%%%%%%%%% STD OF TOTAL CURRENTS %%%%%%%%%%%%%%%%%%%%%%%%'])
disp(['        5          10           25          50          75          90           95%ile'])
disp(['U+iVg: ' num2str(prctile(Collimate(std(Ug_SWOT_HFRgrid_all      + 1i*Vg_SWOT_HFRgrid_all     ,0,3,'omitnan')),[5    10   25   50   75   90   95]))])
disp(['U+iVcg:' num2str(prctile(Collimate(std(Ucg_SWOT_HFRgrid_all     + 1i*Vcg_SWOT_HFRgrid_all    ,0,3,'omitnan')),[5    10   25   50   75   90   95]))])
disp(['HFR:   ' num2str(prctile(Collimate(std(U_rot_filtered_SWOTtimes + 1i*V_rot_filtered_SWOTtimes,0,3,'omitnan')),[5    10   25   50   75   90   95]))])

%% Map vorticity

USER_PICKED_DATE = '2023-05-03 02:00:00';
USER_PICKED_DATE = '2023-05-09 02:00:00';
% USER_PICKED_DATE = '2023-05-10 02:00:00';
% USER_PICKED_DATE = '2023-05-11 02:00:00';
% USER_PICKED_DATE = '2023-05-15 02:00:00';
% USER_PICKED_DATE = '2023-05-21 02:00:00';
% USER_PICKED_DATE = '2023-05-26 02:00:00';
% USER_PICKED_DATE = '2023-06-01 22:00:00';
% USER_PICKED_DATE = '2023-06-04 22:00:00';
% USER_PICKED_DATE = '2023-06-16 00:00:00';
% USER_PICKED_DATE = '2023-06-26 19:00:00';
ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum(USER_PICKED_DATE));
ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));


% % % pixel average uv_g then curl (more analogous to HFR):
VORT_g = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
              mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                          U_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000),...
              mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                          V_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000)...
                          )*[1/111000]./...
         fcor_degrees_cps(NORCAL.HFR.LAT);
VORT_cg = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
               mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                           U_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000),...
               mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                           V_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000)...
                           )*[1/111000]./...
          fcor_degrees_cps(NORCAL.HFR.LAT);

% VORT_g = curl(NORCAL.SWOT.lon{ti_swot}.*cosd(NORCAL.SWOT.lat{ti_swot}), ...
%               NORCAL.SWOT.lat{ti_swot},...
%               U_geostr(:,:,ti_swot), V_geostr(:,:,ti_swot))*[1/111000]./...
%          fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
% 
% VORT_cg = curl(NORCAL.SWOT.lon{ti_swot}.*cosd(NORCAL.SWOT.lat{ti_swot}), ...
%                NORCAL.SWOT.lat{ti_swot},...
%                U_cyclogeostr_1it(:,:,ti_swot), V_cyclogeostr_1it(:,:,ti_swot))*[1/111000]./...
%           fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});

% VORT_hfr = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,... HFR filtered
%                 U_rot_filtered(:,:,ti_hfr), U_rot_filtered(:,:,ti_hfr))*[1/111000]./...
%            fcor_degrees_cps(NORCAL.HFR.LAT);

VORT_hfr = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,... HFR unfiltered
                U_rot_unfiltered(:,:,ti_hfr), V_rot_unfiltered(:,:,ti_hfr))*[1/111000]./...
           fcor_degrees_cps(NORCAL.HFR.LAT);

Vort_CLim = [-1 1]*0.8;
QUIV_STEP = 2; % number of steps between plotted arrows (1 is too crowded)
VEC_SCALE = 0.08;

close all
figure('Color','w')
tiledlayout(1,3,"TileSpacing","compact")

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
set(gcf,'color','w')
COAST = m_gshhs_i('patch',0.8*[1 1 1]);
m_grid('box','fancy', 'backgroundcolor','none'); hold on
set(gcf,'Position',[262    69   962   728])
% m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, VORT_g);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_g);
% m_quiver(NORCAL.SWOT.lon{ti_swot}(1:QUIV_STEP:end,1:QUIV_STEP:end), ...
%          NORCAL.SWOT.lat{ti_swot}(1:QUIV_STEP:end,1:QUIV_STEP:end), ...
%          squeeze(U_geostr(1:QUIV_STEP:end,1:QUIV_STEP:end,ti_swot))*VEC_SCALE, ...
%          squeeze(V_geostr(1:QUIV_STEP:end,1:QUIV_STEP:end,ti_swot))*VEC_SCALE,0,'k');
% CB = colorbar; CB.Label.String = '$\zeta_\mathrm{g}/f$'; CB.Label.Interpreter = 'latex';
% CB.FontSize = 12; CB.Label.FontSize = 20;
xlabel('$\zeta_\mathrm{g}/f$','Interpreter','latex','FontSize',32)
clim(Vort_CLim)
MT = m_text(-123.7568, 40.8015, '(a)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center');
title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'FontSize',14)

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
set(gcf,'color','w')
COAST = m_gshhs_i('patch',0.8*[1 1 1]);
m_grid('box','fancy', 'backgroundcolor','none'); hold on
set(gcf,'Position',[262    69   962   728])
% m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, VORT_cg);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_cg);
% m_quiver(NORCAL.SWOT.lon{ti_swot}(1:QUIV_STEP:end,1:QUIV_STEP:end), ...
%          NORCAL.SWOT.lat{ti_swot}(1:QUIV_STEP:end,1:QUIV_STEP:end), ...
%          squeeze(U_cyclogeostr_Nit(1:QUIV_STEP:end,1:QUIV_STEP:end,ti_swot))*VEC_SCALE, ...
%          squeeze(V_cyclogeostr_Nit(1:QUIV_STEP:end,1:QUIV_STEP:end,ti_swot))*VEC_SCALE,0,'k');
% CB = colorbar; CB.Label.String = '$\zeta_\mathrm{cg}/f$'; CB.Label.Interpreter = 'latex';
% CB.FontSize = 12; CB.Label.FontSize = 20;
xlabel('$\zeta_\mathrm{cg}/f$','Interpreter','latex','FontSize',32)
clim(Vort_CLim)
m_text(-123.7568, 40.8015, '(b)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center')
% title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'FontSize',14)

AX3 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
set(gcf,'color','w')
COAST = m_gshhs_i('patch',0.8*[1 1 1]);
m_grid('box','fancy', 'backgroundcolor','none'); hold on
set(gcf,'Position',[262    69   962   728])
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_hfr);
% m_quiver(NORCAL.HFR.LON, NORCAL.HFR.LAT,... HFR unfiltered
%          U_rot_unfiltered(:,:,ti_hfr)*VEC_SCALE, V_rot_unfiltered(:,:,ti_hfr)*VEC_SCALE,0,'k')
CB = colorbar; CB.Label.String = '$\zeta/f$'; CB.Label.Interpreter = 'latex';
CB.FontSize = 12; CB.Label.FontSize = 20;
xlabel('$\zeta_\mathrm{HFR}/f$','Interpreter','latex','FontSize',32)
clim(Vort_CLim)
m_text(-123.7568, 40.8015, '(c)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center')
title(['HFR time = ' datestr(T0 + NORCAL.HFR.time(ti_hfr)/24)],'FontSize',14)

colormap('bwr')

set(gcf,'Position',[1          63        1440         734])

% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/F_vorticity_maps_g_cg_hfr.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp(['Image saved'])

%% Scatter Vorticity (preparation)

if exist('VORT_GEOSTR','var')
    error('You already calculated these variables, which take a long time to calculate.')
else
end

HFR_choice = 2;
% 1 for u, 2 for u_rot, 3 for u_rot_filtered

USER_PICKED_DATE = '2023-04-01 00:00:00'; % full time range
% USER_PICKED_DATE = '2023-05-01 00:00:00'; % May only
ti_swot_1 = dsearchn(NORCAL.SWOT.mean_time,datenum(USER_PICKED_DATE));

USER_PICKED_DATE = '2023-07-31 23:00:00'; % full time range
% USER_PICKED_DATE = '2023-05-31 23:00:00'; % May only
ti_swot_end = dsearchn(NORCAL.SWOT.mean_time,datenum(USER_PICKED_DATE));

VORT_GEOSTR = [];
VORT_CYCLOG = [];
VORT_HFR = [];

for ti_swot = ti_swot_1:ti_swot_end
    SWOT_t = NORCAL.SWOT.mean_time(ti_swot);
    if isfinite(SWOT_t)
        ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, SWOT_t);
        
        % tiledlayout(1,2,'TileSpacing','compact')

        % % % pixel average uv_g then curl (more analogous to HFR):
        VORT_g = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
                      mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                                  U_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000),...
                      mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                                  V_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000)...
                                  )*[1/111000]./...
                 fcor_degrees_cps(NORCAL.HFR.LAT);
        VORT_cg = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
                       mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                                   U_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000),...
                       mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                                   V_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000)...
                                   )*[1/111000]./...
                  fcor_degrees_cps(NORCAL.HFR.LAT);

        VORT_GEOSTR = [VORT_GEOSTR; Collimate(VORT_g) ];
        VORT_CYCLOG = [VORT_CYCLOG; Collimate(VORT_cg) ];
        
        if HFR_choice == 1 % HFR u
            Vort_hfr = curl(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                            NORCAL.HFR.u(:,:,ti_hfr), NORCAL.HFR.v(:,:,ti_hfr))*[1/111000];
        elseif HFR_choice == 2 % U_rot (unfiltered)
            Vort_hfr = curl(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                            U_rot_unfiltered(:,:,ti_hfr), V_rot_unfiltered(:,:,ti_hfr))*[1/111000];
        elseif HFR_choice == 3 % U_rot_filtered
            Vort_hfr = curl(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                            U_rot_filtered(:,:,ti_hfr), V_rot_filtered(:,:,ti_hfr))*[1/111000];
        else
        end

        VORT_HFR = [VORT_HFR ; Collimate(Vort_hfr./fcor_degrees_cps(NORCAL.HFR.LAT))];
    end
    disp([ti_swot, ti_swot_end])
end

%% Scatter Plot

% LON_LAT_lims = [-180 180 -90 90]; % all data (the entire world)
LON_LAT_lims = [-125 -123.5 38 39.5]; % generous limits of where the southern eddy is
% LON_LAT_lims = [-124.55 -123.75 38.5 39.1]; % stricter limits of where the southern eddy is
% LON_LAT_lims = [-125 -124.25 40 40.5]; % stricter limits of where the northern eddy is
% LON_LAT_lims = [-125 -123.5 38.5 40.5]; % encompassing both eddies and between

% TIME_lims = datenum(['2023-04-01 00:00:00';'2023-07-31 00:00:00']); % full time range
% TIME_lims = datenum(['2023-04-01 00:00:00';'2023-05-01 00:00:00']); % April only
TIME_lims = datenum(['2023-05-01 00:00:00';'2023-06-01 00:00:00']); % May only
% TIME_lims = datenum(['2023-06-01 00:00:00';'2023-07-01 00:00:00']); % June only
% TIME_lims = datenum(['2023-07-01 00:00:00';'2023-08-01 00:00:00']); % July only

SpaceScreen = [];
TimeScreen = [];
for ti_swot = ti_swot_1:ti_swot_end
    SWOT_t = NORCAL.SWOT.mean_time(ti_swot);
    if isfinite(SWOT_t)
        SpaceScreen_i = NORCAL.HFR.LON(:) > LON_LAT_lims(1) & ...
                        NORCAL.HFR.LON(:) < LON_LAT_lims(2) & ...
                        NORCAL.HFR.LAT(:) > LON_LAT_lims(3) & ...
                        NORCAL.HFR.LAT(:) < LON_LAT_lims(4);
        SpaceScreen = [SpaceScreen; SpaceScreen_i(:)];

        if SWOT_t > TIME_lims(1) && SWOT_t < TIME_lims(2)
            TimeScreen = [TimeScreen; true(size(SpaceScreen_i))];
        else
            TimeScreen = [TimeScreen; false(size(SpaceScreen_i))];
        end
    end
end

HBINS = [-2:0.01:2];
FONTSIZE = 12;

% XLIM_a = [-0.4 0.83]; YLIM_a = [-0.6 0.6];
% XLIM_b = [-0.4 0.4]; YLIM_b = [-0.6 .83];
% XLIM_c = XLIM_b; YLIM_c = YLIM_b;

XLIM_a = [-0.6 1]; YLIM_a = [-0.8 0.8];
XLIM_b = [-0.4 0.4]; YLIM_b = [-0.8 0.9];
XLIM_c = XLIM_b; YLIM_c = YLIM_b;

% XLIM_a = [-1 1]; YLIM_a = [-1 1];
% % XLIM_a = [-0.4 0.83]; YLIM_a = [-0.6 0.6];
% XLIM_b = XLIM_a; YLIM_b = YLIM_a;
% XLIM_c = XLIM_b; YLIM_c = YLIM_b;

XTICK = [-1:0.2:1];
YTICK = [-1:0.2:1];

close all
figure('Color','w')
tiledlayout(2,2,"TileSpacing","compact")

AX1 = nexttile([1 2]);
histogram2(VORT_GEOSTR(isfinite(VORT_HFR + VORT_GEOSTR + VORT_CYCLOG) & SpaceScreen & TimeScreen), ...
           VORT_CYCLOG(isfinite(VORT_HFR + VORT_GEOSTR + VORT_CYCLOG) & SpaceScreen & TimeScreen),...
           HBINS,HBINS,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none'); hold on
for r_times_f = [5000 20000 Inf]*fcor_degrees_cps(39.5) % [0.1 1 4 Inf]
    plot(   [0:0.01:2],    2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',1)
    % plot(-1*[0:0.01:2], -1*2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
end
plot([-2 2], [-2 2], 'w','LineWidth',0.67)
text(-0.5,0.7,'(a)','Color','w','FontSize',20,'Interpreter','latex','HorizontalAlignment','center')
TTT1=text(0.7002,0.3334,'$r$ = 5 km','Interpreter','latex','FontSize',12,'Color',    'w','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',15);
TTT2=text(0.6782,0.4730,'$r$ = 20 km','Interpreter','latex','FontSize',12,'Color',   'w','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',30);
TTT3=text(0.4800,0.5290,'$r$ = $\infty$','Interpreter','latex','FontSize',12,'Color','w','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',45);
axis equal
% clim([0 40])
% set(gca,'FontSize',16,'XLim',[-0.5 0.5],'YLim',[-0.5 0.8],'XTick',[-0.4:0.2:0.4],'YTick',[-0.4:0.2:0.8])
% set(gca,'FontSize',16,'XLim',[-0.5 0.5],'YLim',[-0.7 0.8],'XTick',[-0.4:0.2:0.4],'YTick',[-0.8:0.2:0.8])
set(gca,'FontSize',FONTSIZE,'XLim',XLIM_a,'YLim',YLIM_a,'XTick',XTICK,'YTick',YTICK)
% set(gca,'FontSize',FONTSIZE,'XLim',XLIM_a,'YLim',YLIM_a)
xlabel('$\zeta_\mathrm{g}/f$','Interpreter','latex','FontSize',FONTSIZE+8)
ylabel('$\zeta_\mathrm{cg}/f$','Interpreter','latex','FontSize',FONTSIZE+8)

AX2 = nexttile;
histogram2(VORT_HFR(   isfinite(VORT_HFR + VORT_GEOSTR + VORT_CYCLOG) & SpaceScreen & TimeScreen), ...
           VORT_GEOSTR(isfinite(VORT_HFR + VORT_GEOSTR + VORT_CYCLOG) & SpaceScreen & TimeScreen),...
           HBINS,HBINS,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none'); hold on
plot([-0.4 0.4] , [-0.4 0.4], 'w','LineWidth',1)
text(-0.3,0.8,'(b)','Color','w','FontSize',20,'Interpreter','latex','HorizontalAlignment','center')
axis equal
clim([0 40])
% set(gca,'FontSize',16,'XLim',[-0.5 0.5],'YLim',[-0.5 0.8],'XTick',[-0.4:0.2:0.4],'YTick',[-0.4:0.2:0.8])
% set(gca,'FontSize',16,'XLim',[-0.5 0.5],'YLim',[-0.7 0.8],'XTick',[-0.4:0.2:0.4],'YTick',[-0.8:0.2:0.8])
set(gca,'FontSize',FONTSIZE,'XLim',XLIM_b,'YLim',YLIM_b,'XTick',XTICK,'YTick',YTICK)
% set(gca,'FontSize',FONTSIZE,'XLim',XLIM_b,'YLim',YLIM_b)
xlabel('$\zeta_\mathrm{HFR}/f$','Interpreter','latex','FontSize',FONTSIZE+8)
ylabel('$\zeta_\mathrm{g}/f$','Interpreter','latex','FontSize',FONTSIZE+8)

AX3 = nexttile;
histogram2(VORT_HFR(   isfinite(VORT_HFR + VORT_GEOSTR + VORT_CYCLOG) & SpaceScreen & TimeScreen),...
           VORT_CYCLOG(isfinite(VORT_HFR + VORT_GEOSTR + VORT_CYCLOG) & SpaceScreen & TimeScreen),...
           HBINS,HBINS,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none'); hold on
plot([-0.4 0.4] , [-0.4 0.4], 'w','LineWidth',1)
text(-0.3,0.8,'(c)','Color','w','FontSize',20,'Interpreter','latex','HorizontalAlignment','center')
axis equal
clim([0 40])
CB = colorbar;
CB.Label.String = 'PDF'; CB.Label.Interpreter = 'latex'; CB.Label.FontSize = 14;
% set(gca,'FontSize',16,'XLim',[-0.5 0.5],'YLim',[-0.5 0.8],'XTick',[-0.4:0.2:0.4],'YTick',[-0.4:0.2:0.8])
% set(gca,'FontSize',16,'XLim',[-0.5 0.5],'YLim',[-0.7 0.8],'XTick',[-0.4:0.2:0.4],'YTick',[-0.8:0.2:0.8])
set(gca,'FontSize',FONTSIZE,'XLim',XLIM_c,'YLim',YLIM_c,'XTick',XTICK,'YTick',YTICK)
% set(gca,'FontSize',FONTSIZE,'XLim',XLIM_c,'YLim',YLIM_c)
xlabel('$\zeta_\mathrm{HFR}/f$','Interpreter','latex','FontSize',FONTSIZE+8)
ylabel('$\zeta_\mathrm{cg}/f$','Interpreter','latex','FontSize',FONTSIZE+8)

colormap(turbo)

AX1.Color = [1 1 1]*0.5;
AX2.Color = AX1.Color;
AX3.Color = AX1.Color;

% set(gcf,'Position',[-1039         368         739         443])
% set(gcf,'Position',[1    64   541   733])
% set(gcf,'Position',[45    48   465   733])
set(gcf,'Position',[45    154   411   733])

% corrcoef(VORT_HFR(isfinite(VORT_HFR) & isfinite(VORT_GEOSTR)),VORT_GEOSTR(isfinite(VORT_HFR) & isfinite(VORT_GEOSTR)))
% corrcoef(VORT_HFR(isfinite(VORT_HFR) & isfinite(VORT_CYCLOG)),VORT_CYCLOG(isfinite(VORT_HFR) & isfinite(VORT_CYCLOG)))

% [ones(size(VORT_HFR(isfinite(VORT_HFR(:) + VORT_CYCLOG(:))))) VORT_HFR(isfinite(VORT_HFR(:) + VORT_CYCLOG(:)))]\VORT_CYCLOG(isfinite(VORT_HFR(:) + VORT_CYCLOG(:)))
% [ones(size(VORT_HFR(isfinite(VORT_HFR(:) + VORT_GEOSTR(:))))) VORT_HFR(isfinite(VORT_HFR(:) + VORT_GEOSTR(:)))]\VORT_GEOSTR(isfinite(VORT_HFR(:) + VORT_GEOSTR(:)))

% %%
% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/F_vorticity_g_cg_HFR.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')

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
            DIAG = exp(XXYY2/[2*([2*nn + 1]*dxy/4).^2]); % Gaussian, +-2Sigma at edges
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
