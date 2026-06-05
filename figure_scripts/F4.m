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

Multiple_CG_iterations = input('Do you want to run the CG calculation for 2 and 3 iterations? 1 (yes) or 0 (no): ');
Load_Existing_Geostrophic_Variables = input('Do you want to load [ SWOT_and_HFR_velocities_8pix_weighted.mat ]? 1 (yes) or 0 (no): ');

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


% % % Load geostrophic velocities if you don't want/need to calculate them
% % % again:
if Load_Existing_Geostrophic_Variables
    GEOSTR_FILE_TO_LOAD = './SWOT_and_HFR_velocities_8pix_weighted.mat';
    load(GEOSTR_FILE_TO_LOAD)
    disp(['Loading ' GEOSTR_FILE_TO_LOAD])
else
end

% Load DUACS (both products), gridded to HFR and SWOT grids:
DUACS_8_struct = load('./DUACS_8_gridded_to_HFR_and_SWOT_grids.mat');
DUACS_4_struct = load('./DUACS_4_gridded_to_HFR_and_SWOT_grids.mat');

close all

%% Calculate SWOT geostrophic velocity

% Recall that at the surface, geostrophic velocity is:
% u_g = -(g/f)(d\eta/dy)
% v_g = +(g/f)(d\eta/dx)
% where both derivatives are partial derivatives.
% Simple textbook citation:
% https://uw.pressbooks.pub/ocean285/chapter/geostrophic-balance/

close all

if exist('U_geostr','var')
    warning('Skipping the *_geostr calculations.')
else

SWOT_geostr_vel_calculated = true;
gg = 9.807; % m s^-2
GaussianKernelSTD_meters = 4000;

U_geostr = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
V_geostr = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% U_geostr_ = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% V_geostr_ = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);

SSHA_matrix = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);

warning('off','all')



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
                          SSHA, 8, ... % 6 -> 13x13 box, i.e. [26 km]^2
                          true); % <--- final input forces weighted LSF (Gaussian decorr. with 3Sigma at edges)
    % To get geostrophic velocity from inputting surface height:
    U_geostr(:,:,ti_swot) = -(gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,3))/111000;
    V_geostr(:,:,ti_swot) =  (gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,2))./[111000*cosd(NORCAL.SWOT.lat{ti_swot})];

    SSHA_matrix(:,:,ti_swot) = SSHA;

    disp(100*ti_swot/length(NORCAL.SWOT.ssha_karin_2))
end

warning('on','all')

end


%% Show the effect of the spatial smoothing function
% % close all
% SSHA_ = zeros(size(SSHA));
% SSHA_(143,35) = NaN;
% % SSHA_ = bettergaussiansmooth2(SSHA_,2000,GaussianKernelSTD_meters*3);
% % SSHA_ = smoothdata2(SSHA_,'gaussian',15,'omitnan');
% % SSHA_ = filter2(fspecial("gaussian",15,2),SSHA_);
% SSHA_ = nanfilter2(fspecial("gaussian",15,2),SSHA_);
% figure;imagesc(SSHA_)
% axis equal

%% Example to verify that u_g is correct:

ti_swot = 40;

figure
% tiledlayout(1,2,'TileSpacing','compact')
% 
SSHA = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
% 
% nexttile
% pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
%                 NORCAL.SWOT.lat{ti_swot}, ...
%                 SSHA); hold on; shading flat
% contour(NORCAL.SWOT.lon{ti_swot},...
%         NORCAL.SWOT.lat{ti_swot},...
%         SSHA,[-0.3:0.025:0.3],'k')
% set(gcf,'Position',[-1439          42         470         769])
% title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))])
% axis equal
% colormap(bwr)
% clim([-1 1]*0.2)
% 
% nexttile
% 
% SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');


% SSHA = bettergaussiansmooth2(SSHA,2000,GaussianKernelSTD_meters);
% SSHA = SSHA .* ...
%        [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
%        [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
warning('off','all')
[~,SSHA] = fit_6term_2d(NORCAL.SWOT.lon{ti_swot}, ...
                      NORCAL.SWOT.lat{ti_swot}, ...
                      SSHA, 8, true); % 6 -> 13x13 box, i.e. [26 km]^2
warning('on','all')

pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                NORCAL.SWOT.lat{ti_swot}, ...
                SSHA); hold on; shading flat
% pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
%                 NORCAL.SWOT.lat{ti_swot}, sqrt(...
%                 U_geostr(:,:,ti_swot).^2 + ...
%                 V_geostr(:,:,ti_swot).^2)); hold on; shading flat
quiver(NORCAL.SWOT.lon{ti_swot},...
       NORCAL.SWOT.lat{ti_swot},...
       U_geostr(:,:,ti_swot),...
       V_geostr(:,:,ti_swot),   5,'k')
% contour(NORCAL.SWOT.lon{ti_swot},...
%         NORCAL.SWOT.lat{ti_swot},...
%         SSHA,[-0.3:0.025:0.3],'k')
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

% if exist('U_cyclogeostr_Nit','var')
%     warning('Skipping the *_cyclogeostr_*it calculations.')
% else

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


U_cyclogeostr_3it = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
V_cyclogeostr_3it = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
for ti_swot = 1:size(U_geostr,3)
    U_n = U_geostr(:,:,ti_swot);
    V_n = V_geostr(:,:,ti_swot);
    for CG_Iteration = 1:3
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
    U_cyclogeostr_3it(:,:,ti_swot) = U_n1;
    V_cyclogeostr_3it(:,:,ti_swot) = V_n1;
    disp(100*ti_swot/size(U_geostr,3))
end

U_cyclogeostr_4it = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
V_cyclogeostr_4it = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
for ti_swot = 1:size(U_geostr,3)
    U_n = U_geostr(:,:,ti_swot);
    V_n = V_geostr(:,:,ti_swot);
    for CG_Iteration = 1:4
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
    U_cyclogeostr_4it(:,:,ti_swot) = U_n1;
    V_cyclogeostr_4it(:,:,ti_swot) = V_n1;
    disp(100*ti_swot/size(U_geostr,3))
end


CG_Iteration = 1;
U_cyclogeostr_Nit = U_cyclogeostr_1it;
V_cyclogeostr_Nit = V_cyclogeostr_1it;

% end

%% Interpolate the SWOT GEOSTROPHIC AND CYCLOGEOSTROPHIC velocities to the HFR grid to compare in the next step:

if exist('HFR_scatter_times','var')
    error('You already calculated these variables, which take a long time to calculate.')
else
end

% if exist('Ucg_SWOT_HFRgrid_all','var')
%     warning('Skipping the *cg_SWOT_HFRgrid_all calculations.')
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
%%
if Multiple_CG_iterations
    warning('Performing the 2 and 3 iteration CG gridding calculations.')

    Ucg2_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
    Vcg2_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
    Anglecg2_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
    tic
    for ti = 1:length(SWOT_scatter_times)
        ti_swot = SWOT_scatter_times(ti);
        ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));
        Delta_t(ti) = [T0 + NORCAL.HFR.time(ti_hfr)/24] - NORCAL.SWOT.mean_time(ti_swot);
        U_cg_slice = U_cyclogeostr_2it(:,:,ti_swot);
        V_cg_slice = V_cyclogeostr_2it(:,:,ti_swot);

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
        Ucg2_SWOT_HFRgrid_all(:,:,ti) = Ucg_SWOT_HFRgrid;
        Vcg2_SWOT_HFRgrid_all(:,:,ti) = Vcg_SWOT_HFRgrid;
        Anglecg2_SWOT_HFRgrid_all(:,:,ti) = Anglecg_SWOT_HFRgrid;

        disp(100*ti/length(SWOT_scatter_times))
    end

    Ucg3_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
    Vcg3_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
    Anglecg3_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
    tic
    for ti = 1:length(SWOT_scatter_times)
        ti_swot = SWOT_scatter_times(ti);
        ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));
        Delta_t(ti) = [T0 + NORCAL.HFR.time(ti_hfr)/24] - NORCAL.SWOT.mean_time(ti_swot);
        U_cg_slice = U_cyclogeostr_3it(:,:,ti_swot);
        V_cg_slice = V_cyclogeostr_3it(:,:,ti_swot);

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
        Ucg3_SWOT_HFRgrid_all(:,:,ti) = Ucg_SWOT_HFRgrid;
        Vcg3_SWOT_HFRgrid_all(:,:,ti) = Vcg_SWOT_HFRgrid;
        Anglecg3_SWOT_HFRgrid_all(:,:,ti) = Anglecg_SWOT_HFRgrid;

        disp(100*ti/length(SWOT_scatter_times))
    end

    Ucg4_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
    Vcg4_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
    Anglecg4_SWOT_HFRgrid_all = nan(size(NORCAL.HFR.u,1),size(NORCAL.HFR.u,2),length(SWOT_scatter_times));
    tic
    for ti = 1:length(SWOT_scatter_times)
        ti_swot = SWOT_scatter_times(ti);
        ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));
        Delta_t(ti) = [T0 + NORCAL.HFR.time(ti_hfr)/24] - NORCAL.SWOT.mean_time(ti_swot);
        U_cg_slice = U_cyclogeostr_4it(:,:,ti_swot);
        V_cg_slice = V_cyclogeostr_4it(:,:,ti_swot);

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
        Ucg4_SWOT_HFRgrid_all(:,:,ti) = Ucg_SWOT_HFRgrid;
        Vcg4_SWOT_HFRgrid_all(:,:,ti) = Vcg_SWOT_HFRgrid;
        Anglecg4_SWOT_HFRgrid_all(:,:,ti) = Anglecg_SWOT_HFRgrid;

        disp(100*ti/length(SWOT_scatter_times))
    end

else
    warning('Skipping the 2 and 3 iteration CG gridding calculations.')
end

toc

% end

%% Define the filtered part of the HFR current to compare to SWOT geostrophic velocity

close all

% if exist('WINDOW','var')
%     warning('Skipping the L_* calculations.')
% else

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
            L_U(i1,i2,:) = ConvWithWindow(NORCAL.HFR.u(i1,i2,:) - 0.0*mean(NORCAL.HFR.u(i1,i2,:),3,'omitnan'));
            L_V(i1,i2,:) = ConvWithWindow(NORCAL.HFR.v(i1,i2,:) - 0.0*mean(NORCAL.HFR.v(i1,i2,:),3,'omitnan'));
        % else
        % end
    end
    disp(100*i1/size(NORCAL.HFR.u,1))
end
toc

% end

%% Helmholtz Decomposition at all times

if exist('U_rot_filtered','var')
    % error('You already calculated these variables, which take a long time to calculate.')
    warning('Skipping the *cg_SWOT_HFRgrid calculations.')
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

clear urot vrot udiv vdiv

end


%% RMSD and correlation between geostrophic velocity to HFR velocity

if exist('U_rot_filtered_SWOTtimes','var')
    warning('Skipping the L_* calculations.')
else

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

end

%
% Correlation
CUTOFF = 0.8; % minimum fraction of good data required to perform calculation
U_SWOTg_HFR_Corr_R = nan(size(NORCAL.HFR.LON));
U_SWOTg_HFR_Corr_P = nan(size(NORCAL.HFR.LON));
V_SWOTg_HFR_Corr_R = nan(size(NORCAL.HFR.LON));
V_SWOTg_HFR_Corr_P = nan(size(NORCAL.HFR.LON));
RMSD_velmag_SWOTg_HFR = nan(size(NORCAL.HFR.LON));

U_SWOTg_HFRunfiltered_Corr_R = nan(size(NORCAL.HFR.LON));
U_SWOTg_HFRunfiltered_Corr_P = nan(size(NORCAL.HFR.LON));
V_SWOTg_HFRunfiltered_Corr_R = nan(size(NORCAL.HFR.LON));
V_SWOTg_HFRunfiltered_Corr_P = nan(size(NORCAL.HFR.LON));
RMSD_velmag_SWOTg_HFRunfiltered = nan(size(NORCAL.HFR.LON));

U_SWOTcg_HFR_Corr_R = nan(size(NORCAL.HFR.LON));
U_SWOTcg_HFR_Corr_P = nan(size(NORCAL.HFR.LON));
V_SWOTcg_HFR_Corr_R = nan(size(NORCAL.HFR.LON));
V_SWOTcg_HFR_Corr_P = nan(size(NORCAL.HFR.LON));
RMSD_velmag_SWOTcg_HFR = nan(size(NORCAL.HFR.LON));

U_SWOTcg_HFRunfiltered_Corr_R = nan(size(NORCAL.HFR.LON));
U_SWOTcg_HFRunfiltered_Corr_P = nan(size(NORCAL.HFR.LON));
V_SWOTcg_HFRunfiltered_Corr_R = nan(size(NORCAL.HFR.LON));
V_SWOTcg_HFRunfiltered_Corr_P = nan(size(NORCAL.HFR.LON));
RMSD_velmag_SWOTcg_HFRunfiltered = nan(size(NORCAL.HFR.LON));

W_SWOTg_HFR_Corr_R = nan(size(NORCAL.HFR.LON)); % Complex (no p-values possible with them)
W_SWOTg_HFR_Corr_angle = nan(size(NORCAL.HFR.LON));
W_SWOTcg_HFR_Corr_R = nan(size(NORCAL.HFR.LON)); % Complex (no p-values possible with them)
W_SWOTcg_HFR_Corr_angle = nan(size(NORCAL.HFR.LON));

W_SWOTg_HFRunfiltered_Corr_R = nan(size(NORCAL.HFR.LON)); % Complex (no p-values possible with them)
W_SWOTg_HFRunfiltered_Corr_angle = nan(size(NORCAL.HFR.LON));
W_SWOTcg_HFRunfiltered_Corr_R = nan(size(NORCAL.HFR.LON)); % Complex (no p-values possible with them)
W_SWOTcg_HFRunfiltered_Corr_angle = nan(size(NORCAL.HFR.LON));

tic
for ii = 1:size(Ug_SWOT_HFRgrid_all,1)
    for jj = 1:size(Ug_SWOT_HFRgrid_all,2)
        % Geostrophy
        if sum(isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)))/numel(Ug_SWOT_HFRgrid_all(ii,jj,:)) > CUTOFF & ...
           sum(isfinite(U_rot_filtered_SWOTtimes(ii,jj,:)))/numel(U_rot_filtered_SWOTtimes(ii,jj,:)) > CUTOFF
            % Ug vs. filtered HFR
            [RR,PP] = corrcoef(Ug_SWOT_HFRgrid_all(     ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:))), ...
                               U_rot_filtered_SWOTtimes(ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:))));
            U_SWOTg_HFR_Corr_R(ii,jj) = RR(1,2);
            U_SWOTg_HFR_Corr_P(ii,jj) = PP(1,2);

            % Ug vs. unfiltered HFR
            [RR,PP] = corrcoef(Ug_SWOT_HFRgrid_all(     ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:))), ...
                               U_rot_unfiltered_SWOTtimes(ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:))));
            U_SWOTg_HFRunfiltered_Corr_R(ii,jj) = RR(1,2);
            U_SWOTg_HFRunfiltered_Corr_P(ii,jj) = PP(1,2);

            % Vg vs. filtered HFR
            [RR,PP] = corrcoef(Vg_SWOT_HFRgrid_all(     ii,jj,isfinite(Vg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))), ...
                               V_rot_filtered_SWOTtimes(ii,jj,isfinite(Vg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))));
            V_SWOTg_HFR_Corr_R(ii,jj) = RR(1,2);
            V_SWOTg_HFR_Corr_P(ii,jj) = PP(1,2);
            
            % Vg vs. unfiltered HFR
            [RR,PP] = corrcoef(Vg_SWOT_HFRgrid_all(     ii,jj,isfinite(Vg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))), ...
                               V_rot_unfiltered_SWOTtimes(ii,jj,isfinite(Vg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))));
            V_SWOTg_HFRunfiltered_Corr_R(ii,jj) = RR(1,2);
            V_SWOTg_HFRunfiltered_Corr_P(ii,jj) = PP(1,2);

            % Ug+iVg vs. filtered HFR
            [RR] = corrcoef(   Ug_SWOT_HFRgrid_all(  ii,jj,isfinite(   Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))) + ...
                            1i*Vg_SWOT_HFRgrid_all(  ii,jj,isfinite(   Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))), ...
                               U_rot_filtered_SWOTtimes(ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))) + ...
                            1i*V_rot_filtered_SWOTtimes(ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_filtered_SWOTtimes(ii,jj,:)+V_rot_filtered_SWOTtimes(ii,jj,:))));
            W_SWOTg_HFR_Corr_R(ii,jj) = abs(RR(1,2));
            W_SWOTg_HFR_Corr_angle(ii,jj) = angle(RR(1,2))*180/pi;
            RMSD_velmag_SWOTg_HFR(ii,jj) = rms( [Ug_SWOT_HFRgrid_all(ii,jj,:) + 1i*Vg_SWOT_HFRgrid_all(ii,jj,:)] - ...
                                                [U_rot_filtered_SWOTtimes(ii,jj,:) + 1i*U_rot_filtered_SWOTtimes(ii,jj,:)],'omitnan');

            % Ug+iVg vs. unfiltered HFR
            [RR] = corrcoef(   Ug_SWOT_HFRgrid_all(  ii,jj,isfinite(   Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))) + ...
                            1i*Vg_SWOT_HFRgrid_all(  ii,jj,isfinite(   Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))), ...
                               U_rot_unfiltered_SWOTtimes(ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))) + ...
                            1i*V_rot_unfiltered_SWOTtimes(ii,jj,isfinite(Ug_SWOT_HFRgrid_all(ii,jj,:)+Vg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))));
            W_SWOTg_HFRunfiltered_Corr_R(ii,jj) = abs(RR(1,2));
            W_SWOTg_HFRunfiltered_Corr_angle(ii,jj) = angle(RR(1,2))*180/pi;
            RMSD_velmag_SWOTg_HFRunfiltered(ii,jj) = rms( [Ug_SWOT_HFRgrid_all(ii,jj,:) + 1i*Vg_SWOT_HFRgrid_all(ii,jj,:)] - ...
                                                          [U_rot_unfiltered_SWOTtimes(ii,jj,:) + 1i*U_rot_unfiltered_SWOTtimes(ii,jj,:)],'omitnan');
        else
        end
        % Cyclogeostrophy
        if sum(isfinite(Ucg_SWOT_HFRgrid_all(ii,jj,:)))/numel(Ucg_SWOT_HFRgrid_all(ii,jj,:)) > CUTOFF & ...
           sum(isfinite(U_rot_filtered_SWOTtimes(ii,jj,:)))/numel(U_rot_filtered_SWOTtimes(ii,jj,:)) > CUTOFF
            % CG vs. filtered HFR
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
                                                [U_rot_filtered_SWOTtimes(ii,jj,:) + 1i*U_rot_filtered_SWOTtimes(ii,jj,:)],'omitnan');

            % CG vs. unfiltered HFR
            [RR,PP] = corrcoef(Ucg_SWOT_HFRgrid_all(       ii,jj,isfinite(Ucg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:))), ...
                               U_rot_unfiltered_SWOTtimes( ii,jj,isfinite(Ucg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:))));
            U_SWOTcg_HFRunfiltered_Corr_R(ii,jj) = RR(1,2);
            U_SWOTcg_HFRunfiltered_Corr_P(ii,jj) = PP(1,2);
            [RR,PP] = corrcoef(Vcg_SWOT_HFRgrid_all(       ii,jj,isfinite(Vcg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))), ...
                               V_rot_unfiltered_SWOTtimes( ii,jj,isfinite(Vcg_SWOT_HFRgrid_all(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))));
            V_SWOTcg_HFRunfiltered_Corr_R(ii,jj) = RR(1,2);
            V_SWOTcg_HFRunfiltered_Corr_P(ii,jj) = PP(1,2);
            [RR] = corrcoef(   Ucg_SWOT_HFRgrid_all(      ii,jj,isfinite( Ucg_SWOT_HFRgrid_all(ii,jj,:)+Vcg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))) + ...
                            1i*Vcg_SWOT_HFRgrid_all(      ii,jj,isfinite( Ucg_SWOT_HFRgrid_all(ii,jj,:)+Vcg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))), ...
                               U_rot_unfiltered_SWOTtimes(ii,jj,isfinite( Ucg_SWOT_HFRgrid_all(ii,jj,:)+Vcg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))) + ...
                            1i*V_rot_unfiltered_SWOTtimes(ii,jj,isfinite( Ucg_SWOT_HFRgrid_all(ii,jj,:)+Vcg_SWOT_HFRgrid_all(ii,jj,:)+U_rot_unfiltered_SWOTtimes(ii,jj,:)+V_rot_unfiltered_SWOTtimes(ii,jj,:))));
            W_SWOTcg_HFRunfiltered_Corr_R(ii,jj) = abs(RR(1,2));
            W_SWOTcg_HFRunfiltered_Corr_angle(ii,jj) = angle(RR(1,2))*180/pi;
            RMSD_velmag_SWOTcg_HFRunfiltered(ii,jj) = rms([Ucg_SWOT_HFRgrid_all(      ii,jj,:) + 1i*Vcg_SWOT_HFRgrid_all(      ii,jj,:)] - ...
                                                          [U_rot_unfiltered_SWOTtimes(ii,jj,:) + 1i*U_rot_unfiltered_SWOTtimes(ii,jj,:)],'omitnan');
        else
        end
    end
    disp(100*ii/size(Ug_SWOT_HFRgrid_all,1))
end
toc

%% Calculate all vorticity

% if exist('VORT_g_MAT','var')
%     error('You already calculated these variables, which take a long time to calculate.')
% else
% end

tic
VORT_g_MAT   = nan(size(NORCAL.HFR.LON,1), size(NORCAL.HFR.LON,2), length(NORCAL.SWOT.mean_time));
VORT_cg_MAT  = nan(size(NORCAL.HFR.LON,1), size(NORCAL.HFR.LON,2), length(NORCAL.SWOT.mean_time));
if Multiple_CG_iterations
    VORT_cg2_MAT  = nan(size(NORCAL.HFR.LON,1), size(NORCAL.HFR.LON,2), length(NORCAL.SWOT.mean_time));
    VORT_cg3_MAT  = nan(size(NORCAL.HFR.LON,1), size(NORCAL.HFR.LON,2), length(NORCAL.SWOT.mean_time));
    VORT_cg4_MAT  = nan(size(NORCAL.HFR.LON,1), size(NORCAL.HFR.LON,2), length(NORCAL.SWOT.mean_time));
else
end
VORT_hfr_MAT = nan(size(NORCAL.HFR.LON,1), size(NORCAL.HFR.LON,2), length(NORCAL.SWOT.mean_time));
VORT_hfrunfiltered_MAT = nan(size(NORCAL.HFR.LON,1), size(NORCAL.HFR.LON,2), length(NORCAL.SWOT.mean_time));
UV_rot_unfiltered_SWOTtimes = nan(size(NORCAL.HFR.LON,1), size(NORCAL.HFR.LON,2), length(NORCAL.SWOT.mean_time));
for ti_swot = 1:length(NORCAL.SWOT.mean_time)
    ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));

    % % % pixel average uv_g then curl (more analogous to HFR):
    VORT_g = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
                  Ug_SWOT_HFRgrid_all(:,:,ti_swot),...
                  ...mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                  ...            U_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000),...
                  Vg_SWOT_HFRgrid_all(:,:,ti_swot) ...
                  ...mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                  ...            V_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000)...
                              )*[1/111000]./...
             fcor_degrees_cps(NORCAL.HFR.LAT);
    VORT_cg = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
                   Ucg_SWOT_HFRgrid_all(:,:,ti_swot), ...
                   ...mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                   ...            U_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000),...
                   Vcg_SWOT_HFRgrid_all(:,:,ti_swot) ...
                   ...mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
                   ...            V_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000)...
                               )*[1/111000]./...
              fcor_degrees_cps(NORCAL.HFR.LAT);

    if Multiple_CG_iterations
        VORT_cg2 = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
                       Ucg2_SWOT_HFRgrid_all(:,:,ti_swot), ...
                       Vcg2_SWOT_HFRgrid_all(:,:,ti_swot) ...
                                   )*[1/111000]./...
                  fcor_degrees_cps(NORCAL.HFR.LAT);
        VORT_cg3 = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
                       Ucg3_SWOT_HFRgrid_all(:,:,ti_swot), ...
                       Vcg3_SWOT_HFRgrid_all(:,:,ti_swot) ...
                                   )*[1/111000]./...
                  fcor_degrees_cps(NORCAL.HFR.LAT);
        VORT_cg4 = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
                       Ucg4_SWOT_HFRgrid_all(:,:,ti_swot), ...
                       Vcg4_SWOT_HFRgrid_all(:,:,ti_swot) ...
                                   )*[1/111000]./...
                  fcor_degrees_cps(NORCAL.HFR.LAT);
    else
    end
    
    % % % curl then pixel average (not as analogous to HFR, but a better
    % % % representation of SWOT vorticity):
    % VORT_g =   curl(NORCAL.SWOT.lon{ti_swot}.*cosd(NORCAL.SWOT.lat{ti_swot}), ...
    %                 NORCAL.SWOT.lat{ti_swot},...
    %                 U_geostr(:,:,ti_swot), V_geostr(:,:,ti_swot))*[1/111000]./...
    %                 fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
    % VORT_g = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
    %                           VORT_g,NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
    % 
    % VORT_cg =  curl(NORCAL.SWOT.lon{ti_swot}.*cosd(NORCAL.SWOT.lat{ti_swot}), ...
    %                 NORCAL.SWOT.lat{ti_swot},...
    %                 U_cyclogeostr_1it(:,:,ti_swot), V_cyclogeostr_1it(:,:,ti_swot))*[1/111000]./...
    %                 fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
    % VORT_cg = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
    %                            VORT_cg,NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
    
    if Load_Existing_Geostrophic_Variables
        VORT_hfr = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,... HFR filtered
                        U_rot_filtered_SWOTtimes(:,:,ti_swot), V_rot_filtered_SWOTtimes(:,:,ti_swot))*[1/111000]./...
                        fcor_degrees_cps(NORCAL.HFR.LAT);
        VORT_hfrunfiltered = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,... HFR filtered
                                  U_rot_unfiltered_SWOTtimes(:,:,ti_swot), V_rot_unfiltered_SWOTtimes(:,:,ti_swot))*[1/111000]./...
                                  fcor_degrees_cps(NORCAL.HFR.LAT);
    else
        VORT_hfr = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,... HFR filtered
                        U_rot_filtered(:,:,ti_hfr), V_rot_filtered(:,:,ti_hfr))*[1/111000]./...
                        fcor_degrees_cps(NORCAL.HFR.LAT);
        VORT_hfrunfiltered = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,... HFR filtered
                                  U_rot_unfiltered(:,:,ti_hfr), V_rot_unfiltered(:,:,ti_hfr))*[1/111000]./...
                                  fcor_degrees_cps(NORCAL.HFR.LAT);
    end
    
    % VORT_hfr = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,... HFR unfiltered
    %                 U_rot_unfiltered(:,:,ti_hfr), V_rot_unfiltered(:,:,ti_hfr))*[1/111000]./...
    %                 fcor_degrees_cps(NORCAL.HFR.LAT);

    VORT_g_MAT(:,:,ti_swot) = VORT_g;
    VORT_cg_MAT(:,:,ti_swot) = VORT_cg;
    if Multiple_CG_iterations
        VORT_cg2_MAT(:,:,ti_swot) = VORT_cg2;
        VORT_cg3_MAT(:,:,ti_swot) = VORT_cg3;
        VORT_cg4_MAT(:,:,ti_swot) = VORT_cg4;
    else
    end
    VORT_hfr_MAT(:,:,ti_swot) = VORT_hfr;
    VORT_hfrunfiltered_MAT(:,:,ti_swot) = VORT_hfrunfiltered;
    if Load_Existing_Geostrophic_Variables
        UV_rot_unfiltered_SWOTtimes(:,:,ti_swot) = U_rot_unfiltered_SWOTtimes(:,:,ti_swot) + 1i*V_rot_unfiltered_SWOTtimes(:,:,ti_swot);
    else
        UV_rot_unfiltered_SWOTtimes(:,:,ti_swot) = U_rot_unfiltered(:,:,ti_hfr) + 1i*V_rot_unfiltered(:,:,ti_hfr);
    end

    disp([ti_swot, length(NORCAL.SWOT.mean_time)])
end
toc

tic
VORT_duacs4_hfrgrid_MAT = nan(size(DUACS_4_struct.LAT_HFR,1),size(DUACS_4_struct.LAT_HFR,2),length(DUACS_4_struct.time));
VORT_duacs8_hfrgrid_MAT = nan(size(DUACS_8_struct.LAT_HFR,1),size(DUACS_8_struct.LAT_HFR,2),length(DUACS_8_struct.time));
VORT_duacs4_swotgrid_MAT = nan(size(DUACS_4_struct.LAT_SWOT,1),size(DUACS_4_struct.LAT_SWOT,2),length(DUACS_4_struct.time));
VORT_duacs8_swotgrid_MAT = nan(size(DUACS_8_struct.LAT_SWOT,1),size(DUACS_8_struct.LAT_SWOT,2),length(DUACS_8_struct.time));
for ti_duacs = 1:length(DUACS_4_struct.time)
    % % % SWOT grid
    VORT_duacs4s = curl(DUACS_4_struct.LON_SWOT.*cosd(...
                        DUACS_4_struct.LAT_SWOT), ...
                        DUACS_4_struct.LAT_SWOT,...
                        DUACS_4_struct.U_geostr(:,:,ti_duacs),...
                        DUACS_4_struct.V_geostr(:,:,ti_duacs))*[1/111000]./fcor_degrees_cps(...
                        DUACS_4_struct.LAT_SWOT);
    VORT_duacs8s = curl(DUACS_8_struct.LON_SWOT.*cosd(...
                        DUACS_8_struct.LAT_SWOT), ...
                        DUACS_8_struct.LAT_SWOT,...
                        DUACS_8_struct.U_geostr(:,:,ti_duacs),...
                        DUACS_8_struct.V_geostr(:,:,ti_duacs))*[1/111000]./fcor_degrees_cps(...
                        DUACS_8_struct.LAT_SWOT);
    % % % HFR grid
    VORT_duacs4h = curl(DUACS_4_struct.LON_HFR.*cosd(...
                        DUACS_4_struct.LAT_HFR), ...
                        DUACS_4_struct.LAT_HFR,...
                        DUACS_4_struct.U_geostr_HFRgrid(:,:,ti_duacs),...
                        DUACS_4_struct.V_geostr_HFRgrid(:,:,ti_duacs))*[1/111000]./fcor_degrees_cps(...
                        DUACS_4_struct.LAT_HFR);
    VORT_duacs8h = curl(DUACS_8_struct.LON_HFR.*cosd(...
                        DUACS_8_struct.LAT_HFR), ...
                        DUACS_8_struct.LAT_HFR,...
                        DUACS_8_struct.U_geostr_HFRgrid(:,:,ti_duacs),...
                        DUACS_8_struct.V_geostr_HFRgrid(:,:,ti_duacs))*[1/111000]./fcor_degrees_cps(...
                        DUACS_8_struct.LAT_HFR);

    VORT_duacs4_swotgrid_MAT(:,:,ti_duacs) = VORT_duacs4s;
    VORT_duacs8_swotgrid_MAT(:,:,ti_duacs) = VORT_duacs8s;
    
    VORT_duacs4_hfrgrid_MAT(:,:,ti_duacs) = VORT_duacs4h;
    VORT_duacs8_hfrgrid_MAT(:,:,ti_duacs) = VORT_duacs8h;
end
toc

% figure;
% subplot(1,2,1);pcolor_centered(DUACS_4_struct.LON_HFR,DUACS_4_struct.LAT_HFR,VORT_duacs8_hfrgrid_MAT(:,:,50));axis equal;axis tight;
% subplot(1,2,2);pcolor_centered(DUACS_4_struct.LON_SWOT,DUACS_4_struct.LAT_SWOT,VORT_duacs8_swotgrid_MAT(:,:,50));axis equal;axis tight

%% HFR vorticity at all HFR times, and SWOT vorticities at full 2km resolution

VORT_hfr_MAT_alltimes = nan(size(U_rot_unfiltered));
tic
for ii = 1:size(U_rot_unfiltered,3)
    % VORT_hfr_MAT_alltimes(:,:,ii) = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
    %                                      U_rot_unfiltered(:,:,ii), V_rot_unfiltered(:,:,ii))*[1/111000]./...
    %                                      fcor_degrees_cps(NORCAL.HFR.LAT);
    
    VORT_hfr_MAT_alltimes(:,:,ii) = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
                                         U_rot_filtered(:,:,ii), V_rot_filtered(:,:,ii))*[1/111000]./...
                                         fcor_degrees_cps(NORCAL.HFR.LAT);
    if ~rem(ii,20)
        disp([num2str(100*ii/size(U_rot_unfiltered,3))])
    end
end
% toc

VORT_g_MAT_fullres = nan(size(U_geostr));
VORT_cg_MAT_fullres = nan(size(U_geostr));
VORT_cg2_MAT_fullres = nan(size(U_geostr));
VORT_cg3_MAT_fullres = nan(size(U_geostr));
VORT_cg4_MAT_fullres = nan(size(U_geostr));
% tic
for ii = 1:size(U_geostr,3)
    VORT_g_MAT_fullres(:,:,ii) = curl(NORCAL.SWOT.lon{ii}.*cosd(NORCAL.SWOT.lat{ii}), NORCAL.SWOT.lat{ii},...
                                       U_geostr(:,:,ii), V_geostr(:,:,ii)*[1/111000]./...
                                       NORCAL.SWOT.lat{ii});

    VORT_cg_MAT_fullres(:,:,ii) = curl(NORCAL.SWOT.lon{ii}.*cosd(NORCAL.SWOT.lat{ii}), NORCAL.SWOT.lat{ii},...
                                        U_cyclogeostr_Nit(:,:,ii), V_cyclogeostr_Nit(:,:,ii)*[1/111000]./...
                                        NORCAL.SWOT.lat{ii});

    if Multiple_CG_iterations
        VORT_cg2_MAT_fullres(:,:,ii) = curl(NORCAL.SWOT.lon{ii}.*cosd(NORCAL.SWOT.lat{ii}), NORCAL.SWOT.lat{ii},...
                                            U_cyclogeostr_2it(:,:,ii), V_cyclogeostr_2it(:,:,ii)*[1/111000]./...
                                            NORCAL.SWOT.lat{ii});
        VORT_cg3_MAT_fullres(:,:,ii) = curl(NORCAL.SWOT.lon{ii}.*cosd(NORCAL.SWOT.lat{ii}), NORCAL.SWOT.lat{ii},...
                                            U_cyclogeostr_3it(:,:,ii), V_cyclogeostr_3it(:,:,ii)*[1/111000]./...
                                            NORCAL.SWOT.lat{ii});
        VORT_cg4_MAT_fullres(:,:,ii) = curl(NORCAL.SWOT.lon{ii}.*cosd(NORCAL.SWOT.lat{ii}), NORCAL.SWOT.lat{ii},...
                                            U_cyclogeostr_4it(:,:,ii), V_cyclogeostr_4it(:,:,ii)*[1/111000]./...
                                            NORCAL.SWOT.lat{ii});
    else
    end

        disp([num2str(100*ii/size(U_geostr,3))])
end
toc


%%
error('Forced stop by user')
%% Plot a time series of the HFR coverage to determine good dates to compare to HFR

close all

figure
plot(T0 + NORCAL.HFR.time/24, ...
     squeeze(sum(sum(isfinite(U_rot_filtered),1),2))/(size(U_rot_filtered,1)*size(U_rot_filtered,2)), '.-')
title('Fraction of valid filtered, rotational HFR data')
datetick

%% Power spectra comparison

close all

figure('Color','w')
VecScale = 0.12;
Plot_uvMagnitude_Cuttoff = 0.0; % For every velocity magnitude below this, u,v,angle will not be plotted
AX0 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-1 1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,sqrt(Ug_SWOT_HFRgrid.^2 + Vg_SWOT_HFRgrid.^2));
m_contour(NORCAL.HFR.LON, NORCAL.HFR.LAT, sum(isfinite(NORCAL.HFR.u),3)./size(NORCAL.HFR.u,3),[0:0.1:1],'k','LineWidth',2)
m_quiver(NORCAL.HFR.LON,NORCAL.HFR.LAT,VecScale*Ug_SWOT_HFRgrid,VecScale*Vg_SWOT_HFRgrid,0,'k')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
CB = colorbar; CB.Label.String = 'SWOT: |u_g| at HFR grid (m s^{-1})'; CB.Label.FontSize = 20; clim([0 1]*1)
set(gca,'FontSize',16)
title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))])
[x_click,y_click] = ginput(1);
[x_click,y_click] = m_xy2ll(x_click,y_click);
[yi_click,xi_click] = find(abs([NORCAL.HFR.LON    + 1i*NORCAL.HFR.LAT]    - [x_click + 1i*y_click]) == ...
                       min(abs([NORCAL.HFR.LON(:) + 1i*NORCAL.HFR.LAT(:)] - [x_click + 1i*y_click])));
m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'*w')

figure
tiledlayout(3,1,"TileSpacing","compact")
nexttile
plot(NORCAL.HFR.time/24 + T0 , squeeze(NORCAL.HFR.u(yi_click,xi_click,:))   ); hold on
plot(NORCAL.HFR.time/24 + T0 , squeeze(U_rot_filtered(yi_click,xi_click,:)) )
datetick
AX1 = nexttile;
[SPEC_unfiltered, Freq_u, ~] = nunanspectrum(squeeze(NORCAL.HFR.u(yi_click,xi_click,:)) + ...
                                             squeeze(NORCAL.HFR.v(yi_click,xi_click,:))*1i,...
                                             [NORCAL.HFR.time - NORCAL.HFR.time(1)]/24,...
                                             'day','Segments',3,'Window','hanning','Method','nufft');
[SPEC_filtered, Freq_f, ~] = nunanspectrum(squeeze(U_rot_filtered(yi_click,xi_click,:)) + ...
                                           squeeze(V_rot_filtered(yi_click,xi_click,:))*1i,...
                                           [NORCAL.HFR.time - NORCAL.HFR.time(1)]/24,...
                                           'day','Segments',3,'Window','hanning','Method','nufft');
loglog(Freq_u,SPEC_unfiltered(:,1),'k.-'); hold on
loglog(Freq_f,SPEC_filtered(:,1),'g.-');
xlabel('freq (/day)')
ylabel('PSD')
legend('CCW Unfiltered','CCW Filtered')

AX2 = nexttile;
loglog(Freq_u,SPEC_unfiltered(:,2),'k.-'); hold on
loglog(Freq_f,SPEC_filtered(:,2),'g.-');
xlabel('freq (/day)')
ylabel('PSD')
legend('CW Unfiltered','CW Filtered')

set(gcf,'Position',[-1439         122         481         855])

linkaxes([AX1, AX2],'xy')

%% Map {U,V}_geostr interpolated to HFR grid alongside HFR data

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
% tiledlayout(1,3,'TileSpacing','compact')
tiledlayout(2,2,'TileSpacing','tight')

VecScale = 0.12;

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,...
                  NORCAL.HFR.LAT,...
                  sqrt(Ug_SWOT_HFRgrid.^2 + Vg_SWOT_HFRgrid.^2));
m_quiver(NORCAL.HFR.LON,...
         NORCAL.HFR.LAT,...
         VecScale*Ug_SWOT_HFRgrid,VecScale*Vg_SWOT_HFRgrid,0,'k')
% m_contour(NORCAL.HFR.LON,...
%           NORCAL.HFR.LAT,...
%           sqrt(abs(Ug_SWOT_HFRgrid + 1i*Vg_SWOT_HFRgrid)),[0:0.1:1],'k')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
m_text(-122.2568, 41.8015, '(a)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center');
% CB = colorbar; CB.Label.String = 'm s$^{-1}$'; CB.Label.FontSize = 20; CB.Label.Interpreter = 'LaTeX';
CLIM = [0 1]*1;
clim(CLIM);
xlabel('SWOT: $\bf{u}_\mathrm{g}$ on HFR grid','Interpreter','latex','Position',[2.4813e-08 0.0515 -1.0000e+16])
set(gca,'FontSize',16)
ylabel(['SWOT time: ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'Interpreter','latex','VerticalAlignment','top')

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -122.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,...
                  NORCAL.HFR.LAT,...
                  sqrt(Ucg_SWOT_HFRgrid.^2 + Vcg_SWOT_HFRgrid.^2));
m_quiver(NORCAL.HFR.LON,...
         NORCAL.HFR.LAT,...
         VecScale*Ucg_SWOT_HFRgrid,VecScale*Vcg_SWOT_HFRgrid,0,'k')
% m_contour(NORCAL.HFR.LON,...
%           NORCAL.HFR.LAT,...
%           sqrt(abs(Ucg_SWOT_HFRgrid + 1i*Vcg_SWOT_HFRgrid)),[0:0.1:1],'k')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
m_text(-122.2568, 41.8015, '(b)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center');
% CB = colorbar; CB.Label.String = 'm s$^{-1}$'; CB.Label.FontSize = 20; CB.Label.Interpreter = 'LaTeX';
clim(CLIM);
xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
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
% '../figures/draft/SWOT_ug_ucg_HFR_rotfiltered.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% Scatter GEOSTR. vs. CYCLO.

close all
uv_bins = [-2:0.05:2];

figure
tiledlayout(1,3,'TileSpacing','compact')

nexttile
HIST_u = ...
histogram2(U_geostr(:),U_cyclogeostr_Nit(:),... round(U_geostr(:))
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
HIST_angle = ...
histogram2([180/pi]*angle(U_geostr(:) + 1i*V_geostr(:)), ...
           [180/pi]*angle(U_cyclogeostr_Nit(:) + 1i*V_cyclogeostr_Nit(:)),...
           [-180:4.5:180],[-180:4.5:180],'DisplayStyle','tile','Normalization','pdf');
axis equal; xlabel('\theta_g'); ylabel('\theta_{cg}')

colormap('turbo')


CLIM = [-3.5 1];

figure('Color','w')
tiledlayout(1,3,'TileSpacing','compact')
AX1 = nexttile;
XBinEdges = [HIST_u.XBinEdges(2:end) + HIST_u.XBinEdges(1:[end-1])]/2;
YBinEdges = [HIST_u.YBinEdges(2:end) + HIST_u.YBinEdges(1:[end-1])]/2;
[XBinEdges,YBinEdges] = meshgrid(XBinEdges,YBinEdges);
pcolor_centered(XBinEdges,YBinEdges,log10(HIST_u.Values).*HIST_u.Values./HIST_u.Values); hold on
for r_times_f = [5000 20000 Inf]*fcor_degrees_cps(39.5) % [0.1 1 4 Inf]
    plot(   [0:0.01:2],    2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
    plot(-1*[0:0.01:2], -1*2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
end
clim(CLIM)
text(-1.75,1.75,'(a)','Interpreter','latex','FontSize',20,'Color','w','VerticalAlignment','middle','HorizontalAlignment','center')
TTT1=text(-1.4768,-0.5008,'$r$ = 5 km','Interpreter','latex','FontSize',12,'Color',    'w','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',12);
TTT2=text(-1.4551,-0.8180,'$r$ = 20 km','Interpreter','latex','FontSize',12,'Color',   'w','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',26);
TTT3=text(-1.6659,-1.5051,'$r$ = $\infty$','Interpreter','latex','FontSize',12,'Color','w','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',45);
axis equal
xlabel('$u_\mathrm{g}$ (m/s)','Interpreter','latex','FontSize',20);
ylabel('$u_\mathrm{cg}$','Interpreter','latex','FontSize',20);
CB = colorbar;
CB.Label.String = 'log$_{10}$(PDF)'; CB.Label.Interpreter = 'latex'; CB.Label.FontSize = 14;
CB.Location = "northoutside";
pbaspect([1 1 1]);
AX2 = nexttile;
XBinEdges = [HIST_v.XBinEdges(2:end) + HIST_v.XBinEdges(1:[end-1])]/2;
YBinEdges = [HIST_v.YBinEdges(2:end) + HIST_v.YBinEdges(1:[end-1])]/2;
[XBinEdges,YBinEdges] = meshgrid(XBinEdges,YBinEdges);
pcolor_centered(XBinEdges,YBinEdges,log10(HIST_v.Values).*HIST_v.Values./HIST_v.Values); hold on
for r_times_f = [5000 20000 Inf]*fcor_degrees_cps(39.5)
    plot(   [0:0.01:2],    2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
    plot(-1*[0:0.01:2], -1*2*[0:0.01:2]./[1 + sqrt(1 + 4*[0:0.01:2]./[r_times_f])],'w','LineWidth',0.67)
end
clim(CLIM)
text(-1.75,1.75,'(b)','Interpreter','latex','FontSize',20,'Color','w','VerticalAlignment','middle','HorizontalAlignment','center')
axis equal
xlabel('$v_\mathrm{g}$ (m/s)','Interpreter','latex','FontSize',20);
ylabel('$v_\mathrm{cg}$','Interpreter','latex','FontSize',20);
CB = colorbar;
CB.Label.String = ''; % 'log$_{10}$(PDF)';
CB.Label.Interpreter = 'latex'; CB.Label.FontSize = 14; CB.Location = "northoutside";
pbaspect([1 1 1]);
AX3 = nexttile;
XBinEdges = [HIST_angle.XBinEdges(2:end) + HIST_angle.XBinEdges(1:[end-1])]/2;
YBinEdges = [HIST_angle.YBinEdges(2:end) + HIST_angle.YBinEdges(1:[end-1])]/2;
[XBinEdges,YBinEdges] = meshgrid(XBinEdges,YBinEdges);
pcolor_centered(XBinEdges,YBinEdges,log10(HIST_angle.Values).*HIST_angle.Values./HIST_angle.Values);
clim([-7.5 -4])
text(-157.5,157.5,'(c)','Interpreter','latex','FontSize',20,'Color','k','VerticalAlignment','middle','HorizontalAlignment','center')
axis equal
xlabel('$\theta_\mathrm{g}$ ($^\circ$)','Interpreter','latex','FontSize',20);
ylabel('$\theta_\mathrm{cg}$','Interpreter','latex','FontSize',20);
set(gca,'XTick',[-180:90:180],'YTick',[-180:90:180])
CB = colorbar;
CB.Label.String = ''; % 'log$_{10}$(PDF)';
CB.Label.Interpreter = 'latex'; CB.Label.FontSize = 14; CB.Location = "northoutside";
CB.Ticks = [-7:-4];
pbaspect([1 1 1]);
colormap('turbo')

sgtitle('THE AXES ARE FLIPPED, DO NOT USE UNTIL THIS IS SORTED OUT','Color','r')

linkaxes([AX1 AX2],'xy')
AX1.Color = [1 1 1]*0.5;
AX2.Color = [1 1 1]*0.5;
AX3.Color = [1 1 1]*0.5;
set(gcf,'Position',[-1439         271         942         540])


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


%% Map correlation

close all
figure('Color','w')
tiledlayout(2,3,'TileSpacing','compact')

P_threshold = 0.05;
% P_threshold = 0.1;

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  RMSD_velmag_SWOTg_HFR);
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([0 1]*0.5)
CB.Label.String = 'RMSE (Geostr. vs. HFR) (m/s)'; CB.Label.FontSize = 16;

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  U_SWOTg_HFR_Corr_R .* [U_SWOTg_HFR_Corr_P<P_threshold]./[U_SWOTg_HFR_Corr_P<P_threshold]);
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
CB.Label.String = ['U (Geostr. vs. HFR): corr. coef. R (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

AX3 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  V_SWOTg_HFR_Corr_R .* [V_SWOTg_HFR_Corr_P<P_threshold]./[V_SWOTg_HFR_Corr_P<P_threshold]);
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
CB.Label.String = ['V (Geostr. vs. HFR): corr. coef. R (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

% % % % % % % % % % % 

AX4 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  RMSD_velmag_SWOTcg_HFR);
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([0 1]*0.5)
CB.Label.String = 'RMSE (CG vs. HFR) (m/s)'; CB.Label.FontSize = 16;

AX5 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  U_SWOTcg_HFR_Corr_R .* [U_SWOTcg_HFR_Corr_P<P_threshold]./[U_SWOTcg_HFR_Corr_P<P_threshold]);
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
CB.Label.String = ['U (CG vs. HFR): corr. coef. R (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

AX6 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  V_SWOTcg_HFR_Corr_R .* [V_SWOTcg_HFR_Corr_P<P_threshold]./[V_SWOTcg_HFR_Corr_P<P_threshold]);
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
CB.Label.String = ['V (CG vs. HFR): corr. coef. R (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
linkaxes([AX1 AX2 AX3 AX4 AX5 AX6],'xy')

colormap(turbo)
set([AX1 AX2 AX3 AX4 AX5 AX6],'Color',[1 1 1]*0.7)


% figure('Color','w')
% plot(squeeze(U_rot_filtered_SWOTtimes(yi_click,xi_click,:)),...
%      squeeze(V_rot_filtered_SWOTtimes(yi_click,xi_click,:)), '.-'); hold on
% plot(squeeze(Ug_SWOT_HFRgrid_all(yi_click,xi_click,:)),...
%      squeeze(Vg_SWOT_HFRgrid_all(yi_click,xi_click,:)), '.-'); hold on


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
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
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
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
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
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1]*0.2)
CB.Label.String = {['\Delta V (Geostr. vs. HFR): corr. coef. R (where P<' num2str(P_threshold) ')'];...
                   '\Delta>0 = Geostr. is better  ---  \Delta<0 = CycloGeostr. is better'}; CB.Label.FontSize = 16;
linkaxes([AX1 AX2 AX3],'xy')

colormap(bwr)
set([AX1 AX2 AX3],'Color',[1 1 1]*0.7)

figure(2)

[x_click,y_click] = ginput(1);
[x_click,y_click] = m_xy2ll(x_click,y_click);
[yi_click,xi_click] = find(abs([NORCAL.HFR.LON    + 1i*NORCAL.HFR.LAT]    - [x_click + 1i*y_click]) == ...
                       min(abs([NORCAL.HFR.LON(:) + 1i*NORCAL.HFR.LAT(:)] - [x_click + 1i*y_click])));


AX1; m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'*w');
AX2; m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'*w');
AX3; m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'*w');


%% Plot the difference in iteration magnitude to investigate convergence
% For more information, see Section 2.4 of Penven et al. (2014)
% https://doi.org/10.1002/2013JC009528

close all

figure
histogram(abs([Ug_SWOT_HFRgrid_all   + 1i*Vg_SWOT_HFRgrid_all] - ...
              [Ucg_SWOT_HFRgrid_all  + 1i*Vcg_SWOT_HFRgrid_all])); hold on

histogram(abs([Ucg_SWOT_HFRgrid_all  + 1i*Vcg_SWOT_HFRgrid_all] - ...
              [Ucg2_SWOT_HFRgrid_all + 1i*Vcg2_SWOT_HFRgrid_all]));

histogram(abs([Ucg2_SWOT_HFRgrid_all + 1i*Vcg2_SWOT_HFRgrid_all] - ...
              [Ucg3_SWOT_HFRgrid_all + 1i*Vcg3_SWOT_HFRgrid_all]));

histogram(abs([Ucg3_SWOT_HFRgrid_all + 1i*Vcg3_SWOT_HFRgrid_all] - ...
              [Ucg4_SWOT_HFRgrid_all + 1i*Vcg4_SWOT_HFRgrid_all]));

figure
subplot(1,4,1)
pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
    mean(abs([Ug_SWOT_HFRgrid_all   + 1i*Vg_SWOT_HFRgrid_all] - ...
             [Ucg_SWOT_HFRgrid_all  + 1i*Vcg_SWOT_HFRgrid_all]),3,'omitnan'));
clim([0 1]*0.1)
subplot(1,4,2)
pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
    mean(abs([Ucg_SWOT_HFRgrid_all  + 1i*Vcg_SWOT_HFRgrid_all] - ...
             [Ucg2_SWOT_HFRgrid_all + 1i*Vcg2_SWOT_HFRgrid_all]),3,'omitnan'));
clim([0 1]*0.1)
subplot(1,4,3)
pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
    mean(abs([Ucg2_SWOT_HFRgrid_all + 1i*Vcg2_SWOT_HFRgrid_all] - ...
             [Ucg3_SWOT_HFRgrid_all + 1i*Vcg3_SWOT_HFRgrid_all]),3,'omitnan'));
clim([0 1]*0.1)
subplot(1,4,4)
pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
    mean(abs([Ucg3_SWOT_HFRgrid_all + 1i*Vcg3_SWOT_HFRgrid_all] - ...
             [Ucg4_SWOT_HFRgrid_all + 1i*Vcg4_SWOT_HFRgrid_all]),3,'omitnan'));
clim([0 1]*0.1)
colorbar
colormap(turbo)



%% Map vorticity and vorticity time series (PUBLICATION FIGURE)

USER_PICKED_DATE = '2023-05-03 02:00:00';
USER_PICKED_DATE = '2023-05-06 02:00:00';
USER_PICKED_DATE = '2023-05-09 02:00:00';
% USER_PICKED_DATE = '2023-05-10 02:00:00';
% USER_PICKED_DATE = '2023-05-15 02:00:00';
% USER_PICKED_DATE = '2023-05-21 02:00:00';
% USER_PICKED_DATE = '2023-05-26 02:00:00';
% USER_PICKED_DATE = '2023-06-01 22:00:00';
% USER_PICKED_DATE = '2023-06-04 22:00:00';
% USER_PICKED_DATE = '2023-06-16 00:00:00';
% USER_PICKED_DATE = '2023-06-26 19:00:00';
ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum(USER_PICKED_DATE));
ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));

VORT_g = VORT_g_MAT(:,:,ti_swot);
VORT_cg = VORT_cg_MAT(:,:,ti_swot);
if Multiple_CG_iterations
    VORT_cg2 = VORT_cg2_MAT(:,:,ti_swot);
    VORT_cg3 = VORT_cg3_MAT(:,:,ti_swot);
    VORT_cg4 = VORT_cg4_MAT(:,:,ti_swot);
else
end
VORT_hfr = VORT_hfr_MAT(:,:,ti_swot);
VORT_hfrunfiltered = VORT_hfrunfiltered_MAT(:,:,ti_swot);


% % % % pixel average uv_g then curl (more analogous to HFR):
% VORT_g = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
%               mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
%                           U_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000),...
%               mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
%                           V_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000)...
%                           )*[1/111000]./...
%          fcor_degrees_cps(NORCAL.HFR.LAT);
% VORT_cg = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,...
%                mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
%                            U_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000),...
%                mean_in_new_grid(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot},...
%                            V_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000)...
%                            )*[1/111000]./...
%           fcor_degrees_cps(NORCAL.HFR.LAT);
% 
% % VORT_g = curl(NORCAL.SWOT.lon{ti_swot}.*cosd(NORCAL.SWOT.lat{ti_swot}), ...
% %               NORCAL.SWOT.lat{ti_swot},...
% %               U_geostr(:,:,ti_swot), V_geostr(:,:,ti_swot))*[1/111000]./...
% %          fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
% % 
% % VORT_cg = curl(NORCAL.SWOT.lon{ti_swot}.*cosd(NORCAL.SWOT.lat{ti_swot}), ...
% %                NORCAL.SWOT.lat{ti_swot},...
% %                U_cyclogeostr_1it(:,:,ti_swot), V_cyclogeostr_1it(:,:,ti_swot))*[1/111000]./...
% %           fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
% 
% % VORT_hfr = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,... HFR filtered
% %                 U_rot_filtered(:,:,ti_hfr), U_rot_filtered(:,:,ti_hfr))*[1/111000]./...
% %            fcor_degrees_cps(NORCAL.HFR.LAT);
% 
% VORT_hfr = curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT,... HFR unfiltered
%                 U_rot_unfiltered(:,:,ti_hfr), V_rot_unfiltered(:,:,ti_hfr))*[1/111000]./...
%            fcor_degrees_cps(NORCAL.HFR.LAT);

Vort_CLim = [-1 1]*0.601;
MapBackgroundColor = [1 1 1];
MapLandColor = [1 1 1]*0.8;
QUIV_STEP = 2; % number of steps between plotted arrows (1 is too crowded)
VEC_SCALE = 0.08;

% Small function to plot the edges of a matrix (recall that <Collimate = @(IN) IN(:);>):
Edges_Only = @(IN) [Collimate(IN(1,:)); Collimate(IN(:,end)); flip(Collimate(IN(end,:))); flip(Collimate(IN(:,1)))];

close all
figure('Color','w')
% tiledlayout(1,3,"TileSpacing","compact")
tiledlayout(2,2,"TileSpacing","compact")

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
% m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, VORT_g);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_g); hold on
% m_quiver(NORCAL.SWOT.lon{ti_swot}(1:QUIV_STEP:end,1:QUIV_STEP:end), ...
%          NORCAL.SWOT.lat{ti_swot}(1:QUIV_STEP:end,1:QUIV_STEP:end), ...
%          squeeze(U_geostr(1:QUIV_STEP:end,1:QUIV_STEP:end,ti_swot))*VEC_SCALE, ...
%          squeeze(V_geostr(1:QUIV_STEP:end,1:QUIV_STEP:end,ti_swot))*VEC_SCALE,0,'k');
% CB = colorbar; CB.Label.String = '$\zeta_\mathrm{g}/f$'; CB.Label.Interpreter = 'latex';
% CB.FontSize = 12; CB.Label.FontSize = 20;
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'k','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'k--','LineWidth',1)
COAST = m_gshhs_i('patch',MapLandColor);
m_grid('box','fancy', 'backgroundcolor',MapBackgroundColor);
% xlabel('$\zeta_\mathrm{g}/f$','Interpreter','latex','FontSize',32)
m_text(-123.2568, 40.5, '$\frac{\zeta_\mathrm{g}}{f}$', 'fontsize', 40,'Interpreter','latex','HorizontalAlignment','center');
clim(Vort_CLim)
% m_text(-123.7568, 40.8015, '(a)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center');
m_text(-125.25, 40.8015, '(a)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center')
title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'FontSize',14)

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
% m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, VORT_cg);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_cg); hold on
% m_quiver(NORCAL.SWOT.lon{ti_swot}(1:QUIV_STEP:end,1:QUIV_STEP:end), ...
%          NORCAL.SWOT.lat{ti_swot}(1:QUIV_STEP:end,1:QUIV_STEP:end), ...
%          squeeze(U_cyclogeostr_Nit(1:QUIV_STEP:end,1:QUIV_STEP:end,ti_swot))*VEC_SCALE, ...
%          squeeze(V_cyclogeostr_Nit(1:QUIV_STEP:end,1:QUIV_STEP:end,ti_swot))*VEC_SCALE,0,'k');
% CB = colorbar; CB.Label.String = '$\zeta_\mathrm{cg}/f$'; CB.Label.Interpreter = 'latex';
% CB.FontSize = 12; CB.Label.FontSize = 20;
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'k','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'k--','LineWidth',1)
COAST = m_gshhs_i('patch',MapLandColor);
m_grid('box','fancy', 'backgroundcolor',MapBackgroundColor); hold on
% xlabel('$\zeta_\mathrm{cg}/f$','Interpreter','latex','FontSize',32)
m_text(-123.2568, 40.5, '$\frac{\zeta^{(1)}_\mathrm{cg}}{f}$', 'fontsize', 40,'Interpreter','latex','HorizontalAlignment','center')
clim(Vort_CLim)
% m_text(-123.7568, 40.8015, '(b)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center')
m_text(-125.25, 40.8015, '(b)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center')
% title(['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'FontSize',14)

AX3 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_hfrunfiltered); hold on
% m_quiver(NORCAL.HFR.LON, NORCAL.HFR.LAT,... HFR unfiltered
%          U_rot_unfiltered(:,:,ti_hfr)*VEC_SCALE, V_rot_unfiltered(:,:,ti_hfr)*VEC_SCALE,0,'k')
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'k','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'k--','LineWidth',1)
COAST = m_gshhs_i('patch',MapLandColor);
m_grid('box','fancy', 'backgroundcolor',MapBackgroundColor);
% CB = colorbar; CB.Label.String = '$\zeta/f$'; CB.Label.Interpreter = 'latex';
% CB.FontSize = 12; CB.Label.FontSize = 20;
% xlabel('$\zeta_{\mathrm{HFR}_\mathrm{unfiltered}}/f$','Interpreter','latex','FontSize',32)
m_text(-123.2568, 40.5, '$\frac{\zeta^{\mathrm{unfiltered}}_{\mathrm{HFR}}}{f}$', 'fontsize', 40,'Interpreter','latex','HorizontalAlignment','center')
clim(Vort_CLim)
% m_text(-123.7568, 40.8015, '(c)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center')
m_text(-125.25, 40.8015, '(c)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center')
title(['HFR time = ' datestr(T0 + NORCAL.HFR.time(ti_hfr)/24)],'FontSize',14)

AX4 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_hfr); hold on
% m_quiver(NORCAL.HFR.LON, NORCAL.HFR.LAT,... HFR unfiltered
%          U_rot_unfiltered(:,:,ti_hfr)*VEC_SCALE, V_rot_unfiltered(:,:,ti_hfr)*VEC_SCALE,0,'k')
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'k','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'k--','LineWidth',1)
COAST = m_gshhs_i('patch',MapLandColor);
m_grid('box','fancy', 'backgroundcolor',MapBackgroundColor);
CB = colorbar; CB.Label.String = '$\zeta/f$'; CB.Label.Interpreter = 'latex';
CB.FontSize = 12; CB.Label.FontSize = 20;
% xlabel('$\zeta_{\mathrm{HFR}_\mathrm{L.P.}}/f$','Interpreter','latex','FontSize',32)
m_text(-123.2568, 40.5, '$\frac{\zeta^{\mathrm{Low-pass}}_{\mathrm{HFR}}}{f}$', 'fontsize', 40,'Interpreter','latex','HorizontalAlignment','center')
clim(Vort_CLim)
% m_text(-123.7568, 40.8015, '(d)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center')
m_text(-125.25, 40.8015, '(d)', 'fontsize', 24,'Interpreter','latex','HorizontalAlignment','center')
% title(['HFR time = ' datestr(T0 + NORCAL.HFR.time(ti_hfr)/24)],'FontSize',14)

% colormap(bwr*0.999)
colormap(cmocean('balance'))

set(gcf,'color','w') % set(gcf,'Position',[262    69   962   728])
% set(gcf,'Position',[1          63        1440         734])
set(gcf,'Position',[-747    90   748   882])


% %%
figure(2)
set(gcf,'Color','w')

HIST2 = histogram2( VORT_g_MAT(isfinite(VORT_g_MAT + VORT_cg_MAT)),...
                   VORT_cg_MAT(isfinite(VORT_g_MAT + VORT_cg_MAT)),'DisplayStyle','tile','ShowEmptyBins','off','EdgeColor','none','Normalization','pdf');
% HIST2 = histogram2(VORT_g_MAT(:,:,ti_swot + [-10:10]),VORT_cg_MAT(:,:,ti_swot + [-10:10]),'DisplayStyle','tile','ShowEmptyBins','off','EdgeColor','none');
set(gca,'ColorScale','log')
colormap('turbo')
clim([min(HIST2.Values(HIST2.Values > 0)) max(HIST2.Values(HIST2.Values > 0))])
CB = colorbar;
CB.Location = "east";
CB.AxisLocation = "in";
CB.Ticks = 10.^[-1 0 1 2];
CB.Position = [0.6540    0.1278    0.0222    0.4278]; % 0.6832
CB.Label.String = 'PDF';
CB.Label.Interpreter = 'latex';
set(gca,'FontSize',16)
xlabel('$\zeta_\mathrm{g}/f$', 'Interpreter','latex','FontSize',20)
xticks([-0.5 0 0.5])
ylabel('$\zeta_\mathrm{cg}/f$','Interpreter','latex','FontSize',20)
yticks([-0.5 0 0.5])
hold on

vel_vec = [0.1:0.01:max(HIST2.XBinLimits)];
for r_times_f = [2500 5000 10000 Inf]*fcor_degrees_cps(39.5)
    plot(   vel_vec,    2*vel_vec./[1 + sqrt(1 + 4*vel_vec./[r_times_f])],'k','LineWidth',1)
    % plot(   vel_vec,    2*vel_vec./[1 + sqrt(1 + 4*vel_vec./[r_times_f])],'w--','LineWidth',1)
    TTT1 = text(max(HIST2.XBinLimits) + 0.02, 2*vel_vec(end)./[1 + sqrt(1 + 4*vel_vec(end)./[r_times_f])], ...
                replace(['$r$ = ' num2str([r_times_f/fcor_degrees_cps(39.5)]/1000) ' km'],'Inf km','$\infty$'), ...
                'Interpreter','latex','FontSize',16,'Color','k','VerticalAlignment','middle','HorizontalAlignment','left','Rotation',0*15);
end
plot(   -vel_vec,    2*-vel_vec./[1 + sqrt(1 + 4*-vel_vec./[r_times_f])],'k','LineWidth',1)
% plot(   -vel_vec,    2*-vel_vec./[1 + sqrt(1 + 4*-vel_vec./[r_times_f])],'w--','LineWidth',1)

% TTT1=text(max(HIST2.XBinLimits) + 0.02, 0.3523,'$r$ = 5 km',    'Interpreter','latex','FontSize',16,'Color','k','VerticalAlignment','middle','HorizontalAlignment','left','Rotation',0*15);
% TTT2=text(max(HIST2.XBinLimits) + 0.02, 0.4252,'$r$ = 10 km',   'Interpreter','latex','FontSize',16,'Color','k','VerticalAlignment','middle','HorizontalAlignment','left','Rotation',0*30);
% TTT3=text(max(HIST2.XBinLimits) + 0.02, 0.6200,'$r$ = $\infty$','Interpreter','latex','FontSize',16,'Color','k','VerticalAlignment','middle','HorizontalAlignment','left','Rotation',0*45);

axis equal
axis tight

% % % % % % % Another, Supplemental figure to show the higher-iteration CG estimates:
figure('Color','w')
TL = tiledlayout(2,2,"TileSpacing","compact");

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_cg); hold on
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'k','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'k--','LineWidth',1)
COAST = m_gshhs_i('patch',MapLandColor);
m_grid('box','fancy', 'backgroundcolor',MapBackgroundColor); hold on
xlabel('$\zeta^{(1)}_\mathrm{cg}/f$','Interpreter','latex','FontSize',24)
clim(Vort_CLim)
m_text(-123.7568, 40.8015, '(a)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center')
title(TL,['SWOT time = ' datestr(NORCAL.SWOT.mean_time(ti_swot))],'FontSize',14)

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_cg2); hold on
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'k','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'k--','LineWidth',1)
COAST = m_gshhs_i('patch',MapLandColor);
m_grid('box','fancy', 'backgroundcolor',MapBackgroundColor); hold on
xlabel('$\zeta^{(2)}_\mathrm{cg}/f$','Interpreter','latex','FontSize',24)
clim(Vort_CLim)
m_text(-123.7568, 40.8015, '(b)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center')

AX3 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_cg3); hold on
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'k','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'k--','LineWidth',1)
COAST = m_gshhs_i('patch',MapLandColor);
m_grid('box','fancy', 'backgroundcolor',MapBackgroundColor); hold on
% CB = colorbar; CB.Label.String = '$\zeta/f$'; CB.Label.Interpreter = 'latex';
% CB.FontSize = 12; CB.Label.FontSize = 20;
xlabel('$\zeta^{(3)}_\mathrm{cg}/f$','Interpreter','latex','FontSize',24)
clim(Vort_CLim)
m_text(-123.7568, 40.8015, '(c)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center')

AX4 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[min(NORCAL.HFR.lon(:)) -122.5] + [-.1 .1],...
        'latitudes',[38 41] + [-.1 .1]);
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_cg4); hold on
m_plot(Edges_Only(NORCAL.SWOT.lon{ti_swot}), Edges_Only(NORCAL.SWOT.lat{ti_swot}), 'k','LineWidth',1) % SWOT track lines
m_plot(NORCAL.SWOT.lon{ti_swot}(:,35), NORCAL.SWOT.lat{ti_swot}(:,35), 'k--','LineWidth',1)
COAST = m_gshhs_i('patch',MapLandColor);
m_grid('box','fancy', 'backgroundcolor',MapBackgroundColor); hold on
CB = colorbar; CB.Label.String = '$\zeta/f$'; CB.Label.Interpreter = 'latex';
CB.FontSize = 12; CB.Label.FontSize = 20;
xlabel('$\zeta^{(4)}_\mathrm{cg}/f$','Interpreter','latex','FontSize',24)
clim(Vort_CLim)
m_text(-123.7568, 40.8015, '(d)', 'fontsize', 32,'Interpreter','latex','HorizontalAlignment','center')


% colormap(bwr*0.999)
colormap(cmocean('balance'))

set(gcf,'color','w') % set(gcf,'Position',[262    69   962   728])
% set(gcf,'Position',[1          63        1440         734]) % if 1x3
set(gcf,'Position',[830    74   611   718]) % if 2x2


% figure('Color','w')
% histogram(VORT_g_MAT)


% %%
% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/F_vorticity_maps_g_cg_hfr_4panel.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp(['Image saved'])
% 
% figure(2)
% exportgraphics(gcf,...
% '../figures/draft/F_vorticity_scatter_g_vs_cg.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp(['Image saved'])
% 
% figure(3)
% exportgraphics(gcf,...
% '../figures/draft/FS_vorticity_cg1234.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp(['Image saved'])

%% % % % % % % % % % 
% % % % % % % % % % % TIME SERIES (PUBLICATION FIGURE)
% % % % % % % % % % % 

CLICK_POINT = false; % set to false for a pre-set location for comparing time series.

HFR_ROT_MARKER = 'k.-';
SWOTg_ROT_MARKER = 'r.-';
SWOTc_ROT_MARKER = 'c.-';
DUACS_PLOT_COLOR = [0.4 0.7 1]; % [0 0.5 1];
DATE_LINE_ALPHA = 0.25;

% % % % % % % Time for the snapshot map:
USER_PICKED_DATE = '2023-05-03 02:00:00';
USER_PICKED_DATE = '2023-05-09 02:00:00';
% USER_PICKED_DATE = '2023-05-05 02:00:00';
% USER_PICKED_DATE = '2023-05-10 02:00:00';
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
Ug_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot},U_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
Vg_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot},V_geostr(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
Ucg_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot},U_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);
Vcg_SWOT_HFRgrid = mean_in_new_grid(NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot},V_cyclogeostr_Nit(:,:,ti_swot),NORCAL.HFR.LON,NORCAL.HFR.LAT,10000/111000);

close all

CLIM = [0 1];
VecScale = 0.12;

figure('Color','w') % Just to pick
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -123.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                  sqrt(Ug_SWOT_HFRgrid.^2 + Vg_SWOT_HFRgrid.^2) ); TITLE = '|u+iv|';
                  ...U_SWOTg_HFR_Corr_R .* [U_SWOTg_HFR_Corr_P<P_threshold]./[U_SWOTg_HFR_Corr_P<P_threshold]); TITLE = 'C.C._{ug,hfr}';
                  ...V_SWOTg_HFR_Corr_R .* [V_SWOTg_HFR_Corr_P<P_threshold]./[V_SWOTg_HFR_Corr_P<P_threshold]); TITLE = 'C.C._{vg,hfr}';
                  ...sqrt(U_SWOTg_HFR_Corr_R.^2 + V_SWOTg_HFR_Corr_R.^2)); TITLE = 'sqrt(C.C._{ug,hfr}^2 + C.C._{vg,hfr}^2)';
                  ...VORT_hfr); TITLE = '\zeta_{hfr}';
                  ...VORT_g); TITLE = '\zeta_{g}';
                  ...mean(VORT_g_MAT,3,'omitnan')); TITLE = '<\zeta_{g}>';
                  ...std(VORT_g_MAT,0,3,'omitnan')); TITLE = '\sigma(\zeta_{g})';
                  ...RMSD_velmag_SWOTg_HFR);
axis equal; axis tight
title(TITLE)
set(gca,'FontSize',16)
colormap(turbo)

if CLICK_POINT
    [x_click,y_click] = ginput(1);
    [x_click,y_click] = m_xy2ll(x_click,y_click);
    [yi_click,xi_click] = find(abs([NORCAL.HFR.LON    + 1i*NORCAL.HFR.LAT]    - [x_click + 1i*y_click]) == ...
                           min(abs([NORCAL.HFR.LON(:) + 1i*NORCAL.HFR.LAT(:)] - [x_click + 1i*y_click])));
else
    % % % Eddy:
    % xi_click = 23; yi_click = 59;
    % xi_click = 22; yi_click = 58;
    % xi_click = 21; yi_click = 58; % <--- high C_uu
    % xi_click = 21; yi_click = 59;
    xi_click = 21; yi_click = 61; % <--- Main figure; CG has more realistic \zeta than G
    % xi_click = 21; yi_click = 62;
    % xi_click = 23; yi_click = 63;
    % xi_click = 20; yi_click = 62;
    % xi_click = 20; yi_click = 63;
    % % % Coast:
    % xi_click = 19; yi_click = 20;
    % xi_click = 19; yi_click = 21; % <- potentially interesting
    % xi_click = 25; yi_click = 50; % <- potentially interesting
    % xi_click = 23; yi_click = 51;
    % % % Supplemental figures:
    % xi_click = 18; yi_click = 23; % f_timeseries_comparison_supp1.pdf
    % xi_click = 24; yi_click = 51; % f_timeseries_comparison_supp2.pdf
    xi_click = 14; yi_click = 36; % f_timeseries_comparison_supp3.pdf
    
end

m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'ow','MarkerSize',10,'LineWidth',2)
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');


figure('Color','w')
tiledlayout(3,4,'TileSpacing','compact')

nexttile([3 1])
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON,...
                  NORCAL.HFR.LAT,...
                  sqrt(Ug_SWOT_HFRgrid.^2 + Vg_SWOT_HFRgrid.^2));
m_quiver(NORCAL.HFR.LON,...
         NORCAL.HFR.LAT,...
         VecScale*Ug_SWOT_HFRgrid,VecScale*Vg_SWOT_HFRgrid,0,'k')
m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'+k','MarkerSize',10,'LineWidth',2)%,'MarkerEdgeColor','w')
disp([NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click)])
m_text(-125.25,41.75,'(a)','FontSize',20,'Interpreter','latex','HorizontalAlignment','center')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim(CLIM);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
set(gca,'FontSize',16)
title({'SWOT Pass at';[datestr(NORCAL.SWOT.mean_time(ti_swot)) ' UTC']},'Interpreter','latex','FontSize',16)
CB = colorbar; clim([0 1])
CB.Location = 'southoutside'; CB.Label.String = '$|\mathbf{u}_\mathrm{g}|$ (m s$^{-1}$)';
CB.Label.Interpreter = 'LaTeX';
colormap('turbo')



nexttile([1 3]);
M_SIZE = 10; L_WIDTH = 1.5;
plot(NORCAL.SWOT.mean_time, squeeze(VORT_hfr_MAT(yi_click,xi_click,:)), HFR_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'HandleVisibility','on'); hold on
% plot(NORCAL.SWOT.mean_time, squeeze(VORT_hfrunfiltered_MAT(yi_click,xi_click,:)), '.-','MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'Color',[1 1 1]*0.7,'HandleVisibility','off');
plot(NORCAL.SWOT.mean_time, squeeze(VORT_g_MAT(yi_click,xi_click,:)),   SWOTg_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'HandleVisibility','on')
plot(NORCAL.SWOT.mean_time, squeeze(VORT_cg_MAT(yi_click,xi_click,:)),  SWOTc_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)

plot(NORCAL.SWOT.mean_time, squeeze(VORT_cg2_MAT(yi_click,xi_click,:)), '--c', 'MarkerSize',M_SIZE-0*3,'LineWidth',L_WIDTH); % 'Color', [0.2 1 2/3]
plot(NORCAL.SWOT.mean_time, squeeze(VORT_cg3_MAT(yi_click,xi_click,:)), '-.c' , 'MarkerSize',M_SIZE-0*6,'LineWidth',L_WIDTH+0.5); % 'Color', [0.2 1 1/3]
plot(NORCAL.SWOT.mean_time, squeeze(VORT_cg4_MAT(yi_click,xi_click,:)), ':c' , 'MarkerSize',M_SIZE-0*6,'LineWidth',L_WIDTH+0.5); % 'Color', [0.2 1 1/3]

% plot(DUACS_4_struct.time, squeeze(VORT_duacs4_hfrgrid_MAT(yi_click,xi_click,:)), ':','LineWidth',L_WIDTH,'Color',DUACS_PLOT_COLOR,'HandleVisibility','off')
plot(DUACS_8_struct.time, squeeze(VORT_duacs8_hfrgrid_MAT(yi_click,xi_click,:)), '--' ,'LineWidth',L_WIDTH + 0.5,'Color',DUACS_PLOT_COLOR,'HandleVisibility','on')
YLIM = [-0.25 0.5]; % ylim;
plot([1 1]*NORCAL.SWOT.mean_time(ti_swot), YLIM, '-','Color',[0 0 0 DATE_LINE_ALPHA],'LineWidth',L_WIDTH,'HandleVisibility','off')

LEG2 = legend({'$\zeta_{\mathrm{HFR}_\mathrm{LP}}$',...
               '$\zeta_\mathrm{g}$',...
               '$\zeta^{(1)}_\mathrm{cg}$',...
               '$\zeta^{(2)}_\mathrm{cg}$',...
               '$\zeta^{(3)}_\mathrm{cg}$',...
               '$\zeta^{(4)}_\mathrm{cg}$',...
               '$\zeta_{\mathrm{DUACS\ 1/8}^\circ}$'},...
               'NumColumns',4,'Interpreter','LaTeX','FontSize',18,'EdgeColor','none','BackgroundAlpha',0,...
               'Position',[0.6196    0.8356    0.2848    0.0897]); % 0.6240

AA = VORT_hfr_MAT(yi_click,xi_click,:);
AA_uf = VORT_hfrunfiltered_MAT(yi_click,xi_click,:);
BB = VORT_g_MAT(yi_click,xi_click,:);
VORT_SWOTg_HFR_Corr_R = corrcoef(AA(isfinite(AA) & isfinite(BB)),...
                                 BB(isfinite(AA) & isfinite(BB)));
VORT_SWOTg_HFRunfiltered_Corr_R = corrcoef(AA_uf(isfinite(AA_uf) & isfinite(BB)),...
                                           BB(   isfinite(AA_uf) & isfinite(BB)));
BB = VORT_cg_MAT(yi_click,xi_click,:);
VORT_SWOTcg_HFR_Corr_R = corrcoef(AA(isfinite(AA) & isfinite(BB)),...
                                  BB(isfinite(AA) & isfinite(BB)));
VORT_SWOTcg_HFRunfiltered_Corr_R = corrcoef(AA_uf(isfinite(AA_uf) & isfinite(BB)),...
                                            BB(   isfinite(AA_uf) & isfinite(BB)));
BB = VORT_cg2_MAT(yi_click,xi_click,:);
VORT_SWOTcg2_HFR_Corr_R = corrcoef(AA(isfinite(AA) & isfinite(BB)),...
                                   BB(isfinite(AA) & isfinite(BB)));
BB = VORT_cg3_MAT(yi_click,xi_click,:);
VORT_SWOTcg3_HFR_Corr_R = corrcoef(AA(isfinite(AA) & isfinite(BB)),...
                                   BB(isfinite(AA) & isfinite(BB)));
BB = VORT_cg4_MAT(yi_click,xi_click,:);
VORT_SWOTcg4_HFR_Corr_R = corrcoef(AA(isfinite(AA) & isfinite(BB)),...
                                   BB(isfinite(AA) & isfinite(BB)));

text(datenum('10-04-2023','dd-mm-yyyy'),YLIM(2) - diff(YLIM)*0.25,'(b)','FontSize',20,'Interpreter','latex')
title({['$C(\zeta_\mathrm{g},\zeta_{\mathrm{HFR}_\mathrm{LP}})$ = '  num2str(  VORT_SWOTg_HFR_Corr_R(1,2), '%.2f') '\qquad\qquad' ...
        '$C(\zeta^{(1)}_\mathrm{cg},\zeta_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  VORT_SWOTcg_HFR_Corr_R(1,2), '%.2f') '\qquad\qquad' ...
        '$C(\zeta^{(4)}_\mathrm{cg},\zeta_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  VORT_SWOTcg4_HFR_Corr_R(1,2), '%.2f')] ; ...
       ...['$C(\zeta_\mathrm{g},\zeta_{\mathrm{HFR}_\mathrm{full}})$ = '  num2str(  VORT_SWOTg_HFRunfiltered_Corr_R(1,2), '%.2f') '\qquad\qquad' ...
       ... '$C(\zeta_\mathrm{cg},\zeta_{\mathrm{HFR}_\mathrm{full}})$ = ' num2str(  VORT_SWOTcg_HFRunfiltered_Corr_R(1,2), '%.2f')] ...
       ...['$C(\zeta^{(2)}_\mathrm{cg},\zeta_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  VORT_SWOTcg2_HFR_Corr_R(1,2), '%.2f') '\qquad\qquad' ...
       ... '$C(\zeta^{(3)}_\mathrm{cg},\zeta_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  VORT_SWOTcg3_HFR_Corr_R(1,2), '%.2f')] ; ...
       },...
       'Interpreter','latex','FontSize',16)
ylabel(['$\zeta/f$'],'Interpreter','latex','FontSize',16)
set(gca,'YAxisLocation','right')
% LEG = legend('HFR \zeta rotational filtered at SWOT times', ...
%              'SWOT \zeta geostr. on HFR grid', ...
%              'SWOT \zeta cyclogeostr. on HFR grid');
% LEG.Location = 'southeast';
datetick
xlim([datenum('2023-04-07') datenum('2023-07-10')]);
ylim(YLIM)



nexttile([1 3])

AA = U_rot_filtered_SWOTtimes(yi_click,xi_click,:);
BB = Ucg4_SWOT_HFRgrid_all(yi_click,xi_click,:);
U_SWOTcg4_HFR_Corr_R = corrcoef(AA(isfinite(AA) & isfinite(BB)),...
                                BB(isfinite(AA) & isfinite(BB)));
M_SIZE = 10; L_WIDTH = 1.5;
plot(NORCAL.SWOT.mean_time,squeeze(U_rot_filtered_SWOTtimes(yi_click,xi_click,:)), HFR_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH); hold on
% plot(NORCAL.SWOT.mean_time,squeeze(U_rot_unfiltered_SWOTtimes(yi_click,xi_click,:)), '.-','MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'Color',[1 1 1]*0.7);%,'HandleVisibility','off');
plot(NORCAL.SWOT.mean_time,squeeze(Ug_SWOT_HFRgrid_all(yi_click,xi_click,:)),      SWOTg_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)
plot(NORCAL.SWOT.mean_time,squeeze(Ucg_SWOT_HFRgrid_all(yi_click,xi_click,:)),     SWOTc_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)
plot(NORCAL.SWOT.mean_time,squeeze(Ucg4_SWOT_HFRgrid_all(yi_click,xi_click,:)),     'c:','MarkerSize',M_SIZE,'LineWidth',L_WIDTH+0.5)
% plot(DUACS_4_struct.time, squeeze(DUACS_4_struct.U_geostr_HFRgrid(yi_click,xi_click,:)), ':','LineWidth',L_WIDTH,'Color',DUACS_PLOT_COLOR)
plot(DUACS_8_struct.time, squeeze(DUACS_8_struct.U_geostr_HFRgrid(yi_click,xi_click,:)), '--' ,'LineWidth',L_WIDTH + 0.5,'Color',DUACS_PLOT_COLOR)
YLIM = [-0.52 0.67]; % ylim;
text(datenum('10-04-2023','dd-mm-yyyy'),YLIM(2) - diff(YLIM)*0.25,'(c)','FontSize',20,'Interpreter','latex')
title({['$C(u_\mathrm{g},u_{\mathrm{HFR}_\mathrm{LP}})$ = '  num2str(  U_SWOTg_HFR_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
        '$C(u^{(1)}_\mathrm{cg},u_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  U_SWOTcg_HFR_Corr_R(yi_click,xi_click), '%.2f') '\qquad\qquad' ...
        '$C(u^{(4)}_\mathrm{cg},u_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  U_SWOTcg4_HFR_Corr_R(1,2), '%.2f')] ; ...
       ...['$C(u_\mathrm{g},u_{\mathrm{HFR}_\mathrm{full}})$ = '  num2str(U_SWOTg_HFRunfiltered_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
       ... '$C(u_\mathrm{cg},u_{\mathrm{HFR}_\mathrm{full}})$ = ' num2str(U_SWOTcg_HFRunfiltered_Corr_R(yi_click,xi_click), '%.2f')] ...
        },...
        'Interpreter','latex','FontSize',16)
ylabel('$u$ (m/s)','Interpreter','latex','FontSize',16)
set(gca,'YAxisLocation','right')
% LEG1 = legend({'HFR rotational filtered at SWOT times', ...
%               ...'HFR rotational unfiltered at SWOT times', ...
%               'SWOT geostr. on HFR grid', ...
%               'SWOT cyclogeostr. on HFR grid', ...
%               ...'DUACS 1/4^\circ', ...
%               'DUACS 1/8^\circ'},'FontSize',11);
% LEG1.NumColumns = length(LEG1.String); LEG1.EdgeColor = 'w';
% LEG1.Location = 'southeast'; LEG1.Position = [0.2954    0.0628    0.6205    0.0498];
datetick
xlim([datenum('2023-04-07') datenum('2023-07-10')]);
plot([1 1]*NORCAL.SWOT.mean_time(ti_swot), YLIM, '-','Color',[0 0 0 DATE_LINE_ALPHA],'LineWidth',L_WIDTH,'HandleVisibility','off')


nexttile([1 3])

AA = V_rot_filtered_SWOTtimes(yi_click,xi_click,:);
BB = Vcg4_SWOT_HFRgrid_all(yi_click,xi_click,:);
V_SWOTcg4_HFR_Corr_R = corrcoef(AA(isfinite(AA) & isfinite(BB)),...
                                BB(isfinite(AA) & isfinite(BB)));

plot(NORCAL.SWOT.mean_time,squeeze(V_rot_filtered_SWOTtimes(yi_click,xi_click,:)), HFR_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH); hold on
% plot(NORCAL.SWOT.mean_time,squeeze(V_rot_unfiltered_SWOTtimes(yi_click,xi_click,:)), '.-','MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'Color',[1 1 1]*0.7,'HandleVisibility','off');
plot(NORCAL.SWOT.mean_time,squeeze(Vg_SWOT_HFRgrid_all(yi_click,xi_click,:)),      SWOTg_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)
plot(NORCAL.SWOT.mean_time,squeeze(Vcg_SWOT_HFRgrid_all(yi_click,xi_click,:)),     SWOTc_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)
plot(NORCAL.SWOT.mean_time,squeeze(Vcg4_SWOT_HFRgrid_all(yi_click,xi_click,:)),     'c:','MarkerSize',M_SIZE,'LineWidth',L_WIDTH+0.5)
% plot(DUACS_4_struct.time, squeeze(DUACS_4_struct.V_geostr_HFRgrid(yi_click,xi_click,:)), ':','LineWidth',L_WIDTH,'Color',DUACS_PLOT_COLOR)
plot(DUACS_8_struct.time, squeeze(DUACS_8_struct.V_geostr_HFRgrid(yi_click,xi_click,:)), '--' ,'LineWidth',L_WIDTH + 0.5,'Color',DUACS_PLOT_COLOR)
YLIM = ylim;
plot([1 1]*NORCAL.SWOT.mean_time(ti_swot), YLIM, '-','Color',[0 0 0 DATE_LINE_ALPHA],'LineWidth',L_WIDTH,'HandleVisibility','off')

YLIM = [-0.52 0.67]; % ylim;
text(datenum('10-04-2023','dd-mm-yyyy'),YLIM(2) - diff(YLIM)*0.25,'(d)','FontSize',20,'Interpreter','latex')
% title(['$C(v_\mathrm{g},v_\mathrm{HFR})$ = '  num2str(  V_SWOTg_HFR_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
%        '$C(v_\mathrm{cg},v_\mathrm{HFR})$ = ' num2str(  V_SWOTcg_HFR_Corr_R(yi_click,xi_click), '%.2f')],...
%        'Interpreter','latex','FontSize',16)
title({['$C(v_\mathrm{g},v_{\mathrm{HFR}_\mathrm{LP}})$ = '  num2str(  V_SWOTg_HFR_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
        '$C(v^{(1)}_\mathrm{cg},v_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  V_SWOTcg_HFR_Corr_R(yi_click,xi_click), '%.2f') '\qquad\qquad' ...
        '$C(v^{(4)}_\mathrm{cg},v_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  V_SWOTcg4_HFR_Corr_R(1,2), '%.2f')] ; ...
       ...['$C(v_\mathrm{g},v_{\mathrm{HFR}_\mathrm{full}})$ = '  num2str(V_SWOTg_HFRunfiltered_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
       ... '$C(v_\mathrm{cg},v_{\mathrm{HFR}_\mathrm{full}})$ = ' num2str(V_SWOTcg_HFRunfiltered_Corr_R(yi_click,xi_click), '%.2f')] ...
        },...
        'Interpreter','latex','FontSize',16)
ylabel('$v$ (m/s)','Interpreter','latex','FontSize',16)
set(gca,'YAxisLocation','right')
% LEG = legend('HFR V rotational filtered at SWOT times', ...
%              'SWOT V geostr. on HFR grid', ...
%              'SWOT V cyclogeostr. on HFR grid');
% LEG.Location = 'southeast';
datetick
xlim([datenum('2023-04-07') datenum('2023-07-10')]);


set(gcf,'Position',[-1504         142        1141         669])

%% OLD VERSION (Commented)

% figure('Color','w')
% tiledlayout(3,4,'TileSpacing','compact')
% 
% nexttile([3 1])
% M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
%              'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
%              'latitudes',[37 42] + [-.1 .1]);
% set(gcf,'color','w'); hold on
% m_pcolor_centered(NORCAL.HFR.LON,...
%                   NORCAL.HFR.LAT,...
%                   sqrt(Ug_SWOT_HFRgrid.^2 + Vg_SWOT_HFRgrid.^2));
% m_quiver(NORCAL.HFR.LON,...
%          NORCAL.HFR.LAT,...
%          VecScale*Ug_SWOT_HFRgrid,VecScale*Vg_SWOT_HFRgrid,0,'k')
% m_plot(NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click),'+k','MarkerSize',10,'LineWidth',2)%,'MarkerEdgeColor','w')
% disp([NORCAL.HFR.LON(yi_click,xi_click),NORCAL.HFR.LAT(yi_click,xi_click)])
% m_text(-125.25,41.75,'(a)','FontSize',20,'Interpreter','latex','HorizontalAlignment','center')
% COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
% m_grid('box','fancy', 'backgroundcolor','none');
% clim(CLIM);
% % xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
% set(gca,'FontSize',16)
% title({'SWOT Pass at';[datestr(NORCAL.SWOT.mean_time(ti_swot)) ' UTC']},'Interpreter','latex','FontSize',16)
% CB = colorbar; clim([0 1])
% CB.Location = 'southoutside'; CB.Label.String = '$|\mathbf{u}_\mathrm{g}|$ (m s$^{-1}$)';
% CB.Label.Interpreter = 'LaTeX';
% colormap('turbo')
% 
% % subplot(211)
% nexttile([1 3])
% M_SIZE = 10; L_WIDTH = 1.5;
% plot(NORCAL.SWOT.mean_time,squeeze(U_rot_filtered_SWOTtimes(yi_click,xi_click,:)), HFR_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH); hold on
% plot(NORCAL.SWOT.mean_time,squeeze(U_rot_unfiltered_SWOTtimes(yi_click,xi_click,:)), '.-','MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'Color',[1 1 1]*0.7);%,'HandleVisibility','off');
% plot(NORCAL.SWOT.mean_time,squeeze(Ug_SWOT_HFRgrid_all(yi_click,xi_click,:)),      SWOTg_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)
% plot(NORCAL.SWOT.mean_time,squeeze(Ucg_SWOT_HFRgrid_all(yi_click,xi_click,:)),     SWOTc_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)
% plot(DUACS_4_struct.time, squeeze(DUACS_4_struct.U_geostr_HFRgrid(yi_click,xi_click,:)), '--','LineWidth',L_WIDTH,'Color',DUACS_PLOT_COLOR)
% plot(DUACS_8_struct.time, squeeze(DUACS_8_struct.U_geostr_HFRgrid(yi_click,xi_click,:)), ':b' ,'LineWidth',L_WIDTH + 0.5,'Color',DUACS_PLOT_COLOR)
% YLIM = ylim;
% text(datenum('10-04-2023','dd-mm-yyyy'),YLIM(2) - diff(YLIM)*0.25,'(b)','FontSize',20,'Interpreter','latex')
% title({['$C(u_\mathrm{g},u_{\mathrm{HFR}_\mathrm{LP}})$ = '  num2str(  U_SWOTg_HFR_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
%         '$C(u_\mathrm{cg},u_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  U_SWOTcg_HFR_Corr_R(yi_click,xi_click), '%.2f')] ; ...
%        ['$C(u_\mathrm{g},u_{\mathrm{HFR}_\mathrm{full}})$ = '  num2str(U_SWOTg_HFRunfiltered_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
%         '$C(u_\mathrm{cg},u_{\mathrm{HFR}_\mathrm{full}})$ = ' num2str(U_SWOTcg_HFRunfiltered_Corr_R(yi_click,xi_click), '%.2f')]},...
%         'Interpreter','latex','FontSize',16)
% ylabel('$u$ (m/s)','Interpreter','latex','FontSize',16)
% set(gca,'YAxisLocation','right')
% LEG = legend({'HFR rotational filtered at SWOT times', ...
%               'HFR rotational unfiltered at SWOT times', ...
%               'SWOT geostr. on HFR grid', ...
%               'SWOT cyclogeostr. on HFR grid', ...
%               'DUACS 1/4^\circ', ...
%               'DUACS 1/8^\circ'},'NumColumns',2);
% LEG.Location = 'southeast';
% LEG.Position = [0.5855    0.7299    0.3269    0.0643];
% datetick
% xlim([datenum('2023-04-07') datenum('2023-07-10')]);
% YLIM = ylim;
% plot([1 1]*NORCAL.SWOT.mean_time(ti_swot), YLIM, '-','Color',[0 0 0 DATE_LINE_ALPHA],'LineWidth',L_WIDTH,'HandleVisibility','off')
% 
% % subplot(212)
% nexttile([1 3])
% plot(NORCAL.SWOT.mean_time,squeeze(V_rot_filtered_SWOTtimes(yi_click,xi_click,:)), HFR_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH); hold on
% plot(NORCAL.SWOT.mean_time,squeeze(V_rot_unfiltered_SWOTtimes(yi_click,xi_click,:)), '.-','MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'Color',[1 1 1]*0.7,'HandleVisibility','off');
% plot(NORCAL.SWOT.mean_time,squeeze(Vg_SWOT_HFRgrid_all(yi_click,xi_click,:)),      SWOTg_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)
% plot(NORCAL.SWOT.mean_time,squeeze(Vcg_SWOT_HFRgrid_all(yi_click,xi_click,:)),     SWOTc_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)
% plot(DUACS_4_struct.time, squeeze(DUACS_4_struct.V_geostr_HFRgrid(yi_click,xi_click,:)), '--','LineWidth',L_WIDTH,'Color',DUACS_PLOT_COLOR)
% plot(DUACS_8_struct.time, squeeze(DUACS_8_struct.V_geostr_HFRgrid(yi_click,xi_click,:)), ':b' ,'LineWidth',L_WIDTH + 0.5,'Color',DUACS_PLOT_COLOR)
% YLIM = ylim;
% plot([1 1]*NORCAL.SWOT.mean_time(ti_swot), YLIM, '-','Color',[0 0 0 DATE_LINE_ALPHA],'LineWidth',L_WIDTH,'HandleVisibility','off')
% 
% YLIM = ylim;
% text(datenum('10-04-2023','dd-mm-yyyy'),YLIM(2) - diff(YLIM)*0.25,'(c)','FontSize',20,'Interpreter','latex')
% % title(['$C(v_\mathrm{g},v_\mathrm{HFR})$ = '  num2str(  V_SWOTg_HFR_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
% %        '$C(v_\mathrm{cg},v_\mathrm{HFR})$ = ' num2str(  V_SWOTcg_HFR_Corr_R(yi_click,xi_click), '%.2f')],...
% %        'Interpreter','latex','FontSize',16)
% title({['$C(v_\mathrm{g},v_{\mathrm{HFR}_\mathrm{LP}})$ = '  num2str(  V_SWOTg_HFR_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
%         '$C(v_\mathrm{cg},v_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  V_SWOTcg_HFR_Corr_R(yi_click,xi_click), '%.2f')] ; ...
%        ['$C(v_\mathrm{g},v_{\mathrm{HFR}_\mathrm{full}})$ = '  num2str(V_SWOTg_HFRunfiltered_Corr_R( yi_click,xi_click), '%.2f') '\qquad\qquad' ...
%         '$C(v_\mathrm{cg},v_{\mathrm{HFR}_\mathrm{full}})$ = ' num2str(V_SWOTcg_HFRunfiltered_Corr_R(yi_click,xi_click), '%.2f')]},...
%         'Interpreter','latex','FontSize',16)
% ylabel('$v$ (m/s)','Interpreter','latex','FontSize',16)
% set(gca,'YAxisLocation','right')
% % LEG = legend('HFR V rotational filtered at SWOT times', ...
% %              'SWOT V geostr. on HFR grid', ...
% %              'SWOT V cyclogeostr. on HFR grid');
% % LEG.Location = 'southeast';
% datetick
% xlim([datenum('2023-04-07') datenum('2023-07-10')]);
% 
% % Added later: vorticity
% nexttile([1 3]);
% plot(NORCAL.SWOT.mean_time, squeeze(VORT_hfr_MAT(yi_click,xi_click,:)), HFR_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'HandleVisibility','off'); hold on
% plot(NORCAL.SWOT.mean_time, squeeze(VORT_hfrunfiltered_MAT(yi_click,xi_click,:)), '.-','MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'Color',[1 1 1]*0.7,'HandleVisibility','off');
% plot(NORCAL.SWOT.mean_time, squeeze(VORT_g_MAT(yi_click,xi_click,:)),   SWOTg_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH,'HandleVisibility','off')
% plot(NORCAL.SWOT.mean_time, squeeze(VORT_cg_MAT(yi_click,xi_click,:)),  SWOTc_ROT_MARKER,'MarkerSize',M_SIZE,'LineWidth',L_WIDTH)
% 
% plot(NORCAL.SWOT.mean_time, squeeze(VORT_cg2_MAT(yi_click,xi_click,:)), '--c', 'MarkerSize',M_SIZE-0*3,'LineWidth',L_WIDTH); % 'Color', [0.2 1 2/3]
% plot(NORCAL.SWOT.mean_time, squeeze(VORT_cg3_MAT(yi_click,xi_click,:)), '-.c' , 'MarkerSize',M_SIZE-0*6,'LineWidth',L_WIDTH+0.5); % 'Color', [0.2 1 1/3]
% plot(NORCAL.SWOT.mean_time, squeeze(VORT_cg4_MAT(yi_click,xi_click,:)), ':c' , 'MarkerSize',M_SIZE-0*6,'LineWidth',L_WIDTH+0.5); % 'Color', [0.2 1 1/3]
% 
% plot(DUACS_4_struct.time, squeeze(VORT_duacs4_hfrgrid_MAT(yi_click,xi_click,:)), '--','LineWidth',L_WIDTH,'Color',DUACS_PLOT_COLOR,'HandleVisibility','off')
% plot(DUACS_8_struct.time, squeeze(VORT_duacs8_hfrgrid_MAT(yi_click,xi_click,:)), ':b' ,'LineWidth',L_WIDTH + 0.5,'Color',DUACS_PLOT_COLOR,'HandleVisibility','off')
% YLIM = ylim;
% plot([1 1]*NORCAL.SWOT.mean_time(ti_swot), YLIM, '-','Color',[0 0 0 DATE_LINE_ALPHA],'LineWidth',L_WIDTH,'HandleVisibility','off')
% 
% LEG = legend({'$\zeta^{(1)}_\mathrm{cg}$',...
%               '$\zeta^{(2)}_\mathrm{cg}$',...
%               '$\zeta^{(3)}_\mathrm{cg}$',...
%               '$\zeta^{(4)}_\mathrm{cg}$'},...
%               'NumColumns',3,'Interpreter','LaTeX','FontSize',18,'EdgeColor','w');
% 
% AA = VORT_hfr_MAT(yi_click,xi_click,:);
% AA_uf = VORT_hfrunfiltered_MAT(yi_click,xi_click,:);
% BB = VORT_g_MAT(yi_click,xi_click,:);
% VORT_SWOTg_HFR_Corr_R = corrcoef(AA(isfinite(AA) & isfinite(BB)),...
%                                  BB(isfinite(AA) & isfinite(BB)));
% VORT_SWOTg_HFRunfiltered_Corr_R = corrcoef(AA_uf(isfinite(AA_uf) & isfinite(BB)),...
%                                            BB(   isfinite(AA_uf) & isfinite(BB)));
% BB = VORT_cg_MAT(yi_click,xi_click,:);
% VORT_SWOTcg_HFR_Corr_R = corrcoef(AA(isfinite(AA) & isfinite(BB)),...
%                                   BB(isfinite(AA) & isfinite(BB)));
% VORT_SWOTcg_HFRunfiltered_Corr_R = corrcoef(AA_uf(isfinite(AA_uf) & isfinite(BB)),...
%                                             BB(   isfinite(AA_uf) & isfinite(BB)));
% YLIM = ylim;
% text(datenum('10-04-2023','dd-mm-yyyy'),YLIM(2) - diff(YLIM)*0.25,'(d)','FontSize',20,'Interpreter','latex')
% title({['$C(\zeta_\mathrm{g},\zeta_{\mathrm{HFR}_\mathrm{LP}})$ = '  num2str(  VORT_SWOTg_HFR_Corr_R(1,2), '%.2f') '\qquad\qquad' ...
%         '$C(\zeta_\mathrm{cg},\zeta_{\mathrm{HFR}_\mathrm{LP}})$ = ' num2str(  VORT_SWOTcg_HFR_Corr_R(1,2), '%.2f')] ; ...
%        ['$C(\zeta_\mathrm{g},\zeta_{\mathrm{HFR}_\mathrm{full}})$ = '  num2str(  VORT_SWOTg_HFRunfiltered_Corr_R(1,2), '%.2f') '\qquad\qquad' ...
%        '$C(\zeta_\mathrm{cg},\zeta_{\mathrm{HFR}_\mathrm{full}})$ = ' num2str(  VORT_SWOTcg_HFRunfiltered_Corr_R(1,2), '%.2f')]},...
%        'Interpreter','latex','FontSize',16)
% ylabel(['$\zeta/f$'],'Interpreter','latex','FontSize',16)
% set(gca,'YAxisLocation','right')
% % LEG = legend('HFR \zeta rotational filtered at SWOT times', ...
% %              'SWOT \zeta geostr. on HFR grid', ...
% %              'SWOT \zeta cyclogeostr. on HFR grid');
% % LEG.Location = 'southeast';
% datetick
% xlim([datenum('2023-04-07') datenum('2023-07-10')]);
% ylim(YLIM)
% 
% % set(gcf,'Position',[-1504         271        1025         540])
% % set(gcf,'Position',[-1504         142        1025         669])
% set(gcf,'Position',[-1504         142        1141         669])

%%
% figure(2); exportgraphics(gcf,...
% '../figures/draft/f_timeseries_comparison_v2.pdf',...
% 'BackgroundColor','none','ContentType','vector')

% figure(2); exportgraphics(gcf,...
% '../figures/draft/f_timeseries_comparison_supp3.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% VORTICITY TOTAL HISTOGRAM

HIST_BINS = -20:0.1:25;

close all
figure('Color','w')
tiledlayout(2,1,"TileSpacing","tight")

nexttile
H1 = histogram(VORT_g_MAT_fullres(isfinite(VORT_g_MAT_fullres)), HIST_BINS, ...
               'Normalization','pdf', ...
               'DisplayStyle','stairs', ...
               'EdgeColor',hex2rgb('#1F77B4'), ...
               'LineWidth',3); hold on
H2 = histogram(VORT_cg_MAT_fullres(isfinite(VORT_cg_MAT_fullres)), HIST_BINS, ...
               'Normalization','pdf', ...
               'DisplayStyle','stairs', ...
               'EdgeColor',hex2rgb('#FF7F0E'), ...
               'LineWidth',3);
H2_2 = histogram(VORT_cg2_MAT_fullres(isfinite(VORT_cg2_MAT_fullres)), HIST_BINS, ...
                 'Normalization','pdf', ...
                 'DisplayStyle','stairs', ...
                 'EdgeColor',hex2rgb('#CF660A'), ...
                 'LineWidth',2);
H2_3 = histogram(VORT_cg3_MAT_fullres(isfinite(VORT_cg3_MAT_fullres)), HIST_BINS, ...
                 'Normalization','pdf', ...
                 'DisplayStyle','stairs', ...
                 'EdgeColor',hex2rgb('#A65208'), ...
                 'LineWidth',1);
H3 = histogram(VORT_hfr_MAT_alltimes(isfinite(VORT_hfr_MAT_alltimes) & mean(NORCAL.HFR.hdop,3,'omitnan') < 0.6), HIST_BINS, ...
               'Normalization','pdf', ...
               'DisplayStyle','stairs', ...
               'EdgeColor',hex2rgb('#2CA02C'), ...
               'LineWidth',3);

AX = gca;
AX.FontSize = 16;
AX.LineWidth = 1.2;
AX.TickDir = 'out';
AX.Box = 'on';
AX.YScale = 'log';
AX.XLim = [-5 8];
AX.YLim = 10.^[-5 1];
AX.YMinorTick = 'on';
AX.XMinorTick = 'on';
% xlabel('Relative Vorticity $\zeta/f$','Interpreter','latex','FontSize',16)
xlabel('','Interpreter','latex','FontSize',16)
    AX.XTickLabels = {};
ylabel('Probability Density','Interpreter','latex','FontSize',16)
grid on
AX.GridAlpha = 0.15;
AX.MinorGridAlpha = 0.08;

% Legend
LEG = legend(AX,[H1 H2 H2_2 H2_3 H3], ...
             {['$\zeta_\mathrm{g}/f$, 2 km'], ...
              ['$\zeta_\mathrm{cg}/f$, 2 km'], ...
              ['$\zeta^{(2)}_\mathrm{cg}/f$, 2 km'], ...
              ['$\zeta^{(3)}_\mathrm{cg}/f$, 2 km'], ...
              ['$\zeta_\mathrm{HFR}/f$, hourly']}, ...
              'Interpreter','latex', ...
              'FontSize',20);%,...
% LEG.Position = [0.6957    0.6906    0.1596    0.2127];

text(-4.8252,1.5172,{'(a)'},'Interpreter','latex','FontSize',20,'VerticalAlignment','bottom')
SKEWNESS_CHART = {'Skewness' ; ...
    [num2str(skewness(H1.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H1.Data)),'%.3f')] ; 
    [num2str(skewness(H2.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H2.Data)),'%.3f')] ; 
    [num2str(skewness(H2_2.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H2_2.Data)),'%.3f')] ; 
    [num2str(skewness(H2_3.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H2_3.Data)),'%.3f')] ; 
    [num2str(skewness(H3.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H3.Data)),'%.3f')] };
text(5.5,12*10^-4,SKEWNESS_CHART,'Interpreter','latex','FontSize',16,'VerticalAlignment','bottom')
text(5.5,12*10^-4,{'';'g: ';'cg: ';'HFR: '},'Interpreter','latex','FontSize',16,...
     'HorizontalAlignment','right','VerticalAlignment','bottom')
% set(gcf,'Position',[-1122         449         643         362])

nexttile
H4 = histogram(VORT_g_MAT(isfinite(VORT_g_MAT)), HIST_BINS, ...
               'Normalization','pdf', ...
               'DisplayStyle','stairs', ...
               'EdgeColor',hex2rgb('#1F77B4'), ...
               'LineWidth',3); hold on
H5 = histogram(VORT_cg_MAT(isfinite(VORT_cg_MAT)), HIST_BINS, ...
               'Normalization','pdf', ...
               'DisplayStyle','stairs', ...
               'EdgeColor',hex2rgb('#FF7F0E'), ...
               'LineWidth',3);
H5_2 = histogram(VORT_cg2_MAT(isfinite(VORT_cg2_MAT)), HIST_BINS, ...
                 'Normalization','pdf', ...
                 'DisplayStyle','stairs', ...
                 'EdgeColor',hex2rgb('#CF660A'), ...
                 'LineWidth',2);
H5_3 = histogram(VORT_cg3_MAT(isfinite(VORT_cg3_MAT)), HIST_BINS, ...
                 'Normalization','pdf', ...
                 'DisplayStyle','stairs', ...
                 'EdgeColor',hex2rgb('#A65208'), ...
                 'LineWidth',1);
H6 = histogram(VORT_hfr_MAT(isfinite(VORT_hfr_MAT)), HIST_BINS, ...
               'Normalization','pdf', ...
               'DisplayStyle','stairs', ...
               'EdgeColor',hex2rgb('#2CA02C'), ...
               'LineWidth',3);

AX2 = gca;
AX2.FontSize = 16;
AX2.LineWidth = 1.2;
AX2.TickDir = 'out';
AX2.Box = 'on';
AX2.YScale = 'log';
AX2.XLim = [-5 8];
AX2.YLim = AX.YLim;
AX2.YMinorTick = 'on';
AX2.XMinorTick = 'on';
xlabel('Relative Vorticity $\zeta/f$','Interpreter','latex','FontSize',16)
ylabel('Probability Density','Interpreter','latex','FontSize',16)
grid on
AX2.GridAlpha = 0.15;
AX2.MinorGridAlpha = 0.08;

% Legend
LEG = legend(AX2,[H4 H5 H5_2 H5_3 H6], ...
             {['$\zeta_\mathrm{g}/f$, on HFR grid'], ...
              ['$\zeta_\mathrm{cg}/f$, on HFR grid'], ...
              ['$\zeta^{(2)}_\mathrm{cg}/f$, on HFR grid'], ...
              ['$\zeta^{(3)}_\mathrm{cg}/f$, on HFR grid'], ...
              ['$\zeta_\mathrm{HFR}/f$, at SWOT times']}, ...
              'Interpreter','latex', ...
              'FontSize',20);%,...
% LEG.Position = [0.6957    0.6906    0.1596    0.2127];

text(-4.8252,1.5172,{'(b)'},'Interpreter','latex','FontSize',20,'VerticalAlignment','bottom')
SKEWNESS_CHART = {'Skewness' ; ...
    [num2str(skewness(H4.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H4.Data)),'%.2f')] ; 
    [num2str(skewness(H5.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H5.Data)),'%.2f')] ; 
    [num2str(skewness(H5_2.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H5_2.Data)),'%.2f')] ; 
    [num2str(skewness(H5_3.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H5_3.Data)),'%.2f')] ; 
    [num2str(skewness(H6.Data,0),'%.2f') ' $\pm$ ' num2str(sqrt(6/length(H6.Data)),'%.2f')] };
text(5.5,12*10^-4,SKEWNESS_CHART,'Interpreter','latex','FontSize',16,'VerticalAlignment','bottom')
text(5.5,12*10^-4,{'';'g: ';'cg: ';'HFR: '},'Interpreter','latex','FontSize',16,...
     'HorizontalAlignment','right','VerticalAlignment','bottom')


set(gcf,'Position',[-1148         215         669         642])

% %%
% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/F_vorticity_histograms.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp(['Image saved'])

%%
ii = 10;
% [1 89]
close all
get_rid_of_nans = @(IN) IN(isfinite(IN));
figure
set(gcf,'Position',[-1535 292 1313 540])


for ii = 1:size(VORT_g_MAT_fullres,3)
    subplot(131)
    M_P = m_proj('Lambert',...
                 'longitudes',[-127 -122.5] + [-.1 .1],...
                 'latitudes',[35 42] + [-.1 .1]);
    % m_imagesc(VORT_g_MAT_alltimes(:,:,ii))
    m_pcolor_centered(NORCAL.SWOT.lon{ii},NORCAL.SWOT.lat{ii},...
                      VORT_g_MAT_fullres(:,:,ii));
    axis equal; axis tight; colorbar
    clim([-2 2]);
    colormap(cmocean('balance'))
    POSITION1 = get(gca,'Position');
    
    subplot(132)
    M_P = m_proj('Lambert',...
                 'longitudes',[-127 -122.5] + [-.1 .1],...
                 'latitudes',[35 42] + [-.1 .1]);
    % m_imagesc(VORT_g_MAT(:,:,ii))
    m_pcolor_centered(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
                      VORT_g_MAT(:,:,ii));
    axis equal; axis tight; colorbar
    clim([-2 2]*std(VORT_g_MAT(:),'omitnan')/std(VORT_g_MAT_fullres(:),'omitnan'));
    colormap(cmocean('balance'))
    POSITION2 = get(gca,'Position');
    set(gca,'Position',[POSITION2(1) POSITION1(2:4)]);
    
    subplot(133)
    hold off
    histogram(get_rid_of_nans(VORT_g_MAT_fullres(:,:,ii)),'Normalization','pdf'); hold on
    histogram(get_rid_of_nans(VORT_g_MAT(:,:,ii)),'Normalization','pdf')
    legend('\zeta/f at full resolution','\zeta/f at 6km resolution')
    xlim([-1 1]*6)

    disp(ii)
    pause(0.1)
end


%% % % % % % % % % %  VORTICITY
% % % % % % % % % % % TIME SERIES
% % % % % % % % % % % 

close all

figure('Color','w')
tiledlayout(2,4,'TileSpacing','compact')

nexttile([2 1])
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, VORT_hfr);
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim(Vort_CLim*6/8);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
set(gca,'FontSize',16)
colormap('bwr')

[x_click,y_click] = ginput(2);
[x_click,y_click] = m_xy2ll(x_click,y_click);
[yi_click(1),xi_click(1)] = find(abs([NORCAL.HFR.LON    + 1i*NORCAL.HFR.LAT]    - [x_click(1) + 1i*y_click(1)]) == ...
                             min(abs([NORCAL.HFR.LON(:) + 1i*NORCAL.HFR.LAT(:)] - [x_click(1) + 1i*y_click(1)])));
[yi_click(2),xi_click(2)] = find(abs([NORCAL.HFR.LON    + 1i*NORCAL.HFR.LAT]    - [x_click(2) + 1i*y_click(2)]) == ...
                             min(abs([NORCAL.HFR.LON(:) + 1i*NORCAL.HFR.LAT(:)] - [x_click(2) + 1i*y_click(2)])));

m_plot(NORCAL.HFR.LON(yi_click(1),xi_click(1)),NORCAL.HFR.LAT(yi_click(1),xi_click(1)),'ok','MarkerSize',10,'LineWidth',2)
m_plot(NORCAL.HFR.LON(yi_click(2),xi_click(2)),NORCAL.HFR.LAT(yi_click(2),xi_click(2)),'+k','MarkerSize',10,'LineWidth',2)
CB = colorbar; CB.Label.String = '$\zeta_{\mathrm{HFR}}/f$'; CB.Label.Interpreter = 'latex';
CB.FontSize = 12; CB.Label.FontSize = 20;


nexttile([1 3])
plot(NORCAL.SWOT.mean_time,squeeze(VORT_hfr_MAT(yi_click(1),xi_click(1),:)), 'k.-','LineWidth',1.5,'MarkerSize',12); hold on
plot(NORCAL.SWOT.mean_time,squeeze(VORT_g_MAT(  yi_click(1),xi_click(1),:)), 'r.-','LineWidth',1.5,'MarkerSize',12)
plot(NORCAL.SWOT.mean_time,squeeze(VORT_cg_MAT( yi_click(1),xi_click(1),:)), 'c.-','LineWidth',1.5,'MarkerSize',12)
ylim([-1 1]*0.6)
% title({['CorrCoef(U_g , U_{HFR} ) = ' num2str(   U_SWOTg_HFR_Corr_R( yi_click,xi_click))];...
%        ['CorrCoef(U_{cg} , U_{HFR} ) = ' num2str(U_SWOTcg_HFR_Corr_R(yi_click,xi_click))]})
% legend('HFR U rotational filtered at SWOT times', ...
%        'SWOT U geostr. SWOT on HFR grid', ...
%        'SWOT U cyclogeostr. SWOT on HFR grid')
LEG = legend('HFR \zeta/f', 'SWOT \zeta_g/f', 'SWOT \zeta_{cg}/f');
LEG.FontSize = 16;
ylabel('$\zeta/f$','Interpreter','latex')
datetick


nexttile([1 3])
plot(NORCAL.SWOT.mean_time,squeeze(VORT_hfr_MAT(yi_click(2),xi_click(2),:)), 'k.-','LineWidth',1.5,'MarkerSize',12); hold on
plot(NORCAL.SWOT.mean_time,squeeze(VORT_g_MAT(  yi_click(2),xi_click(2),:)), 'r.-','LineWidth',1.5,'MarkerSize',12)
plot(NORCAL.SWOT.mean_time,squeeze(VORT_cg_MAT( yi_click(2),xi_click(2),:)), 'c.-','LineWidth',1.5,'MarkerSize',12)
ylim([-1 1]*0.6)
% title({['CorrCoef(U_g , U_{HFR} ) = ' num2str(   U_SWOTg_HFR_Corr_R( yi_click,xi_click))];...
%        ['CorrCoef(U_{cg} , U_{HFR} ) = ' num2str(U_SWOTcg_HFR_Corr_R(yi_click,xi_click))]})
% legend('HFR U rotational filtered at SWOT times', ...
%        'SWOT U geostr. SWOT on HFR grid', ...
%        'SWOT U cyclogeostr. SWOT on HFR grid')
ylabel('$\zeta/f$','Interpreter','latex')
datetick

set(gcf,'Position',[-1504         271        1295         540])

%% Scatter Plot

HBINS = [-2:0.01:2];
FONTSIZE = 12;

% % % Filter to only display vorticity that fulfills certain qualifications:
IS_GOOD = isfinite(VORT_hfr_MAT + VORT_g_MAT + VORT_cg_MAT);
% IS_GOOD = [isfinite(VORT_hfr_MAT + VORT_g_MAT + VORT_cg_MAT)] & ...
%           [mean(NORCAL.HFR.hdop,3,'omitnan') < 0.25] & ...
%           [abs(UV_rot_unfiltered_SWOTtimes) > 0.5];
%           % [abs(Ug_SWOT_HFRgrid_all + 1i*Vg_SWOT_HFRgrid_all) > 0.2];

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
histogram2(VORT_g_MAT(IS_GOOD), ...
           VORT_cg_MAT(IS_GOOD),...
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
histogram2(VORT_hfr_MAT(   IS_GOOD), ...
           VORT_g_MAT(IS_GOOD),...
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
histogram2(VORT_hfr_MAT(   IS_GOOD),...
           VORT_cg_MAT(IS_GOOD),...
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

% corrcoef(VORT_hfr_MAT(isfinite(VORT_hfr_MAT) & isfinite(VORT_g_MAT)),VORT_g_MAT(isfinite(VORT_hfr_MAT) & isfinite(VORT_g_MAT)))
% corrcoef(VORT_hfr_MAT(isfinite(VORT_hfr_MAT) & isfinite(VORT_cg_MAT)),VORT_cg_MAT(isfinite(VORT_hfr_MAT) & isfinite(VORT_cg_MAT)))

% [ones(size(VORT_hfr_MAT(isfinite(VORT_hfr_MAT(:) + VORT_cg_MAT(:))))) VORT_hfr_MAT(isfinite(VORT_hfr_MAT(:) + VORT_cg_MAT(:)))]\VORT_cg_MAT(isfinite(VORT_hfr_MAT(:) + VORT_cg_MAT(:)))
% [ones(size(VORT_hfr_MAT(isfinite(VORT_hfr_MAT(:) + VORT_g_MAT(:))))) VORT_hfr_MAT(isfinite(VORT_hfr_MAT(:) + VORT_g_MAT(:)))]\VORT_g_MAT(isfinite(VORT_hfr_MAT(:) + VORT_g_MAT(:)))

% %%
% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/F_vorticity_g_cg_HFR.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')

%% Scatter Plot of change in RMSD(vel) versus zeta/f for effect of CG on DeltaRMSD

Collimate = @(IN) IN(:);
vortBINS = [-1:0.01:1];
rmsdBINS = [-0.5:0.01:0.75];
[XX,YY] = meshgrid([vortBINS(2:end) + vortBINS(1:[end-1])]/2,...
                   [rmsdBINS(2:end) + rmsdBINS(1:[end-1])]/2);

close all

figure
subplot(1,2,1)
HIST = ...
histogram2(Collimate(VORT_g_MAT),...
           abs(Collimate([Ucg_SWOT_HFRgrid_all     + 1i*Vcg_SWOT_HFRgrid_all] - [U_rot_filtered_SWOTtimes + 1i*U_rot_filtered_SWOTtimes])) - ...
           abs(Collimate([Ug_SWOT_HFRgrid_all      + 1i*Vcg_SWOT_HFRgrid_all] - [U_rot_filtered_SWOTtimes + 1i*U_rot_filtered_SWOTtimes])),...
           vortBINS,rmsdBINS,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
axis equal
xlabel(['\zeta_g/f'])
ylabel({['|w_{cg} - w_{HFR,r,f}| - |w_{g} - w_{HFR,r,f}|'] ; ...
        ['Positive: G is better -- Negative: CG is better']})
% scatter(Collimate(VORT_g_MAT),...
%            Collimate([Ucg_SWOT_HFRgrid_all     + 1i*Vcg_SWOT_HFRgrid_all] - ...
%                      [U_rot_filtered_SWOTtimes + 1i*U_rot_filtered_SWOTtimes]),...
%        '.')
subplot(1,2,2)
HIST_Values = HIST.Values'./~~HIST.Values';
pcolor_centered(XX,YY,log10(HIST_Values));
% axis equal;
grid on;
hold on
xlabel(['\zeta_g/f'])
ylabel({['|w_{cg} - w_{HFR,r,f}| - |w_{g} - w_{HFR,r,f}|'] ; ...
        ['Positive: G is better -- Negative: CG is better']})

% Evenly-space
DX_eb = 0.1; % Vorticity bin for averaging
X_eb = -0.6:DX_eb:0.6; Y_eb = nan(size(X_eb)); YNEG_eb = Y_eb; YPOS_eb = Y_eb;
XNEG_eb = ones(size(X_eb))*DX_eb/2;
XPOS_eb = ones(size(X_eb))*DX_eb/2;

% Manually-chosen spacing:
X_eb =    [-0.3 0.0 0.3]; Y_eb = nan(size(X_eb)); YNEG_eb = Y_eb; YPOS_eb = Y_eb;
XNEG_eb = [ 0.7 0.2 0.1];
XPOS_eb = [ 0.1 0.2 0.7];

for ii = 1:length(X_eb)
    % Mean and STD:
    % Y_eb(ii) = mean(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
    %                           HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),'omitnan');
    % YNEG_eb = std(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
    %                         HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),'omitnan');

    % Median and interquartile range:
    Y_eb(ii) = median(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
                                HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),'omitnan');
    YNEG_eb(ii) = prctile(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
                                    HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),25);
    YPOS_eb(ii) = prctile(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
                                    HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),75);
end
% YPOS_eb = YNEG_eb; % <--- Mean and STD
YPOS_eb = YPOS_eb - Y_eb; YNEG_eb = Y_eb - YNEG_eb; % <--- Median and interquartile range

errorbar(X_eb,Y_eb,YNEG_eb,YPOS_eb,XNEG_eb,XPOS_eb,'.r')

% %%
% % This produces incorrect results because of NaNs:
% disp(['Fraction of points with increased RMSD for G versus CG (zeta/f < -0.2):   ' ...
%       num2str(sum(HIST.Data(HIST.Data(:,1)<-0.2,2) >= 0 )/sum(HIST.Data(:,1)<-0.2))])
% disp(['Fraction of points with increased RMSD for G versus CG (zeta/f > -0.2 and < 0.2):   ' ...
%       num2str(sum(HIST.Data(HIST.Data(:,1)>-0.2 & HIST.Data(:,1)<0.2,2) >= 0 )/sum(HIST.Data(:,1)>-0.2 & HIST.Data(:,1)<0.2))])
% disp(['Fraction of points with increased RMSD for G versus CG (zeta/f > +0.2):   ' ...
%       num2str(sum(HIST.Data(HIST.Data(:,1)>0.2,2) >= 0 )/sum(HIST.Data(:,1)>0.2))])


% This produces correct results after eliminating NaNs:
HIST_Data = HIST.Data(isfinite(HIST.Data(:,1)) & isfinite(HIST.Data(:,2)),:);
disp(['Fraction of points with increased RMSD for G versus CG (zeta/f < -0.2):   ' ...
      num2str(sum(HIST_Data(HIST_Data(:,1)<-0.2,2) >= 0 )/sum(HIST_Data(:,1)<-0.2))])
inv_prctile(HIST_Data(HIST_Data(:,1)<-0.2,2),0)

disp(['Fraction of points with increased RMSD for G versus CG (zeta/f > -0.2 and < 0.2):   ' ...
      num2str(sum(HIST_Data(HIST_Data(:,1)>-0.2 & HIST_Data(:,1)<0.2,2) >= 0 )/sum(HIST_Data(:,1)>-0.2 & HIST_Data(:,1)<0.2))])
inv_prctile(HIST_Data(HIST_Data(:,1)>-0.2 & HIST_Data(:,1)<0.2,2),0)

disp(['Fraction of points with increased RMSD for G versus CG (zeta/f > +0.2):   ' ...
      num2str(sum(HIST_Data(HIST_Data(:,1)>0.2,2) >= 0 )/sum(HIST_Data(:,1)>0.2))])
inv_prctile(HIST_Data(HIST_Data(:,1)>0.2,2),0)

%% Scatter Plot of change in RMSD(vort) versus zeta/f for effect of CG on DeltaRMSD
% This finding may be referenced in the manuscript in section 3.2

Collimate = @(IN) IN(:);
vortBINS = [-1:0.01:1];
rmsdBINS = [-0.5:0.01:0.75];
[XX,YY] = meshgrid([vortBINS(2:end) + vortBINS(1:[end-1])]/2,...
                   [rmsdBINS(2:end) + rmsdBINS(1:[end-1])]/2);

close all

figure
subplot(1,2,1)
HIST = ...
histogram2(Collimate(VORT_g_MAT),... Collimate([VORT_g_MAT + VORT_cg_MAT]/2)
           abs(Collimate([VORT_cg_MAT] - [VORT_hfr_MAT])) - ...
           abs(Collimate([VORT_g_MAT]  - [VORT_hfr_MAT])),...
           vortBINS,rmsdBINS,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
axis equal
xlabel(['\zeta_g/f'])
ylabel({['|\zeta_{cg} - \zeta_{HFR,r,f}| - |\zeta_{g} - \zeta_{HFR,r,f}|'] ; ...
        ['Positive: G is better -- Negative: CG is better']})
% scatter(Collimate(VORT_g_MAT),...
%            Collimate([VORT_cg_MAT] - ...
%                      [VORT_hfr_MAT]),...
%        '.')
subplot(1,2,2)
HIST_Values = HIST.Values'./~~HIST.Values';
pcolor_centered(XX,YY,log10(HIST_Values));
% axis equal;
grid on;
hold on
xlabel(['\zeta_g/f'])
ylabel({['|\zeta_{cg} - \zeta_{HFR,r,f}| - |\zeta_{g} - \zeta_{HFR,r,f}|'] ; ...
        ['Positive: G is better -- Negative: CG is better']})

% Evenly-space
DX_eb = 0.1; % Vorticity bin for averaging
X_eb = -0.6:DX_eb:0.6; Y_eb = nan(size(X_eb)); YNEG_eb = Y_eb; YPOS_eb = Y_eb;
XNEG_eb = ones(size(X_eb))*DX_eb/2;
XPOS_eb = ones(size(X_eb))*DX_eb/2;

% Manually-chosen spacing:
X_eb =    [-0.3 0.0 0.3]; Y_eb = nan(size(X_eb)); YNEG_eb = Y_eb; YPOS_eb = Y_eb;
XNEG_eb = [ 0.7 0.2 0.1];
XPOS_eb = [ 0.1 0.2 0.7];

for ii = 1:length(X_eb)
    % Mean and STD:
    Y_eb(ii) = mean(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
                              HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),'omitnan');
    YNEG_eb = std(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
                            HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),'omitnan');

    % % Median and interquartile range:
    % Y_eb(ii) = median(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
    %                             HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),'omitnan');
    % YNEG_eb(ii) = prctile(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
    %                                 HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),25);
    % YPOS_eb(ii) = prctile(HIST.Data(HIST.Data(:,1) > X_eb(ii)-XNEG_eb(ii) & ...
    %                                 HIST.Data(:,1) < X_eb(ii)+XPOS_eb(ii),2),75);
end
YPOS_eb = YNEG_eb; % <--- Mean and STD
% YPOS_eb = YPOS_eb - Y_eb; YNEG_eb = Y_eb - YNEG_eb; % <--- Median and interquartile range

errorbar(X_eb,Y_eb,YNEG_eb,YPOS_eb,XNEG_eb,XPOS_eb,'.r')

% %%
% % This produces incorrect results because of NaNs:
% disp(['Fraction of points with increased RMSD for G versus CG (zeta/f < -0.2):   ' ...
%       num2str(sum(HIST.Data(HIST.Data(:,1)<-0.2,2) >= 0 )/sum(HIST.Data(:,1)<-0.2))])
% disp(['Fraction of points with increased RMSD for G versus CG (zeta/f > -0.2 and < 0.2):   ' ...
%       num2str(sum(HIST.Data(HIST.Data(:,1)>-0.2 & HIST.Data(:,1)<0.2,2) >= 0 )/sum(HIST.Data(:,1)>-0.2 & HIST.Data(:,1)<0.2))])
% disp(['Fraction of points with increased RMSD for G versus CG (zeta/f > +0.2):   ' ...
%       num2str(sum(HIST.Data(HIST.Data(:,1)>0.2,2) >= 0 )/sum(HIST.Data(:,1)>0.2))])


% This produces correct results after eliminating NaNs:
HIST_Data = HIST.Data(isfinite(HIST.Data(:,1)) & isfinite(HIST.Data(:,2)),:);
disp(['Fraction of points with increased RMSD for G versus CG (zeta/f < -0.2):   ' ...
      num2str(sum(HIST_Data(HIST_Data(:,1)<-0.2,2) >= 0 )/sum(HIST_Data(:,1)<-0.2))])
inv_prctile(HIST_Data(HIST_Data(:,1)<-0.2,2),0)

disp(['Fraction of points with increased RMSD for G versus CG (zeta/f > -0.2 and < 0.2):   ' ...
      num2str(sum(HIST_Data(HIST_Data(:,1)>-0.2 & HIST_Data(:,1)<0.2,2) >= 0 )/sum(HIST_Data(:,1)>-0.2 & HIST_Data(:,1)<0.2))])
inv_prctile(HIST_Data(HIST_Data(:,1)>-0.2 & HIST_Data(:,1)<0.2,2),0)

disp(['Fraction of points with increased RMSD for G versus CG (zeta/f > +0.2):   ' ...
      num2str(sum(HIST_Data(HIST_Data(:,1)>0.2,2) >= 0 )/sum(HIST_Data(:,1)>0.2))])
inv_prctile(HIST_Data(HIST_Data(:,1)>0.2,2),0)

%% Evolution of CG estimate after multiple iterations:

LINPLOT = false;

% % % % % % % Time for the snapshot map:
USER_PICKED_DATE = '2023-05-03 02:00:00';
USER_PICKED_DATE = '2023-05-05 02:00:00';
USER_PICKED_DATE = '2023-05-09 02:00:00';
% USER_PICKED_DATE = '2023-05-10 02:00:00';
% USER_PICKED_DATE = '2023-05-15 02:00:00';
% USER_PICKED_DATE = '2023-05-21 02:00:00';
% USER_PICKED_DATE = '2023-05-26 02:00:00';
USER_PICKED_DATE = '2023-06-01 22:00:00';
% USER_PICKED_DATE = '2023-06-04 22:00:00';
% USER_PICKED_DATE = '2023-06-16 00:00:00';
USER_PICKED_DATE = '2023-06-26 19:00:00';
% USER_PICKED_DATE = '2023-07-05 17:00:00';
% USER_PICKED_DATE = '2023-07-09 17:00:00';
ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum(USER_PICKED_DATE));

if LINPLOT
    CLIM = [0 0.25];
else
    CLIM = [-5 0];
end
CMAP = 'turbo';

close all

figure('Color','w')
tiledlayout(1,3,"TileSpacing","compact")

% G -> CG(1)
nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
if LINPLOT
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                  (abs(U_cyclogeostr_1it(:,:,ti_swot) - U_geostr(:,:,ti_swot))));
                  ...log10(abs(U_cyclogeostr_1it(:,:,ti_swot) - U_geostr(:,:,ti_swot))));
else
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                  ...(abs(U_cyclogeostr_1it(:,:,ti_swot) - U_geostr(:,:,ti_swot))));
                  log10(abs(U_cyclogeostr_1it(:,:,ti_swot) - U_geostr(:,:,ti_swot))));
end
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim(CLIM);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
title(['$|u^{(1)}_\mathrm{cg} - u_\mathrm{g}|$'],'Interpreter','latex')
set(gca,'FontSize',16)
colormap(CMAP)

% CG(1) -> CG(2)
nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
if LINPLOT
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                      (abs(U_cyclogeostr_2it(:,:,ti_swot) - U_cyclogeostr_1it(:,:,ti_swot))));
                      ...log10(abs(U_cyclogeostr_2it(:,:,ti_swot) - U_cyclogeostr_1it(:,:,ti_swot))));
else
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                      ...(abs(U_cyclogeostr_2it(:,:,ti_swot) - U_cyclogeostr_1it(:,:,ti_swot))));
                      log10(abs(U_cyclogeostr_2it(:,:,ti_swot) - U_cyclogeostr_1it(:,:,ti_swot))));
end
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim(CLIM);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
title(['$|u^{(2)}_\mathrm{cg} - u^{(1)}_\mathrm{cg}|$'],'Interpreter','latex')
set(gca,'FontSize',16)
colormap(CMAP)

% CG(2) -> CG(3)
nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
if LINPLOT
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                      (abs(U_cyclogeostr_3it(:,:,ti_swot) - U_cyclogeostr_2it(:,:,ti_swot))));
                      ...log10(abs(U_cyclogeostr_3it(:,:,ti_swot) - U_cyclogeostr_2it(:,:,ti_swot))));
else
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                      ...(abs(U_cyclogeostr_3it(:,:,ti_swot) - U_cyclogeostr_2it(:,:,ti_swot))));
                      log10(abs(U_cyclogeostr_3it(:,:,ti_swot) - U_cyclogeostr_2it(:,:,ti_swot))));
end
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim(CLIM);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
title(['$|u^{(3)}_\mathrm{cg} - u^{(2)}_\mathrm{cg}|$'],'Interpreter','latex')
set(gca,'FontSize',16)
colormap(CMAP)


figure('Color','w')
tiledlayout(1,3,"TileSpacing","compact")

% G -> CG(1)
nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
if LINPLOT
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                  (abs(V_cyclogeostr_1it(:,:,ti_swot) - V_geostr(:,:,ti_swot))));
                  ...log10(abs(V_cyclogeostr_1it(:,:,ti_swot) - V_geostr(:,:,ti_swot))));
else
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                  ...(abs(V_cyclogeostr_1it(:,:,ti_swot) - V_geostr(:,:,ti_swot))));
                  log10(abs(V_cyclogeostr_1it(:,:,ti_swot) - V_geostr(:,:,ti_swot))));
end
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim(CLIM);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
title(['$|v^{(1)}_\mathrm{cg} - v_\mathrm{g}|$'],'Interpreter','latex')
set(gca,'FontSize',16)
colormap(CMAP)

% CG(1) -> CG(2)
nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
if LINPLOT
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                      (abs(V_cyclogeostr_2it(:,:,ti_swot) - V_cyclogeostr_1it(:,:,ti_swot))));
                      ...log10(abs(V_cyclogeostr_2it(:,:,ti_swot) - V_cyclogeostr_1it(:,:,ti_swot))));
else
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                      ...(abs(V_cyclogeostr_2it(:,:,ti_swot) - V_cyclogeostr_1it(:,:,ti_swot))));
                      log10(abs(V_cyclogeostr_2it(:,:,ti_swot) - V_cyclogeostr_1it(:,:,ti_swot))));
end
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim(CLIM);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
title(['$|v^{(2)}_\mathrm{cg} - v^{(1)}_\mathrm{cg}|$'],'Interpreter','latex')
set(gca,'FontSize',16)
colormap(CMAP)

% CG(2) -> CG(3)
nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
if LINPLOT
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                      (abs(V_cyclogeostr_3it(:,:,ti_swot) - V_cyclogeostr_2it(:,:,ti_swot))));
                      ...log10(abs(V_cyclogeostr_3it(:,:,ti_swot) - V_cyclogeostr_2it(:,:,ti_swot))));
else
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                      ...(abs(V_cyclogeostr_3it(:,:,ti_swot) - V_cyclogeostr_2it(:,:,ti_swot))));
                      log10(abs(V_cyclogeostr_3it(:,:,ti_swot) - V_cyclogeostr_2it(:,:,ti_swot))));
end
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim(CLIM);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
title(['$|v^{(3)}_\mathrm{cg} - v^{(2)}_\mathrm{cg}|$'],'Interpreter','latex')
set(gca,'FontSize',16)
colormap(CMAP)

% % % % % % % 
% % % % % % % 
% % % % % % % 

% %% Histograms (to do)
% 
% figure('Color','w')
% tiledlayout(1,2,"TileSpacing","compact")
% nexttile
% histogram(abs(U_cyclogeostr_1it(:,:,ti_swot) - U_geostr(:,:,ti_swot)))

% % % % % % % 
% % % % % % % 
% % % % % % % 

% %% Change in step size (third versus first)

figure('Color','w')
tiledlayout(1,2,"TileSpacing","compact")

nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
              (abs(U_cyclogeostr_3it(:,:,ti_swot) - U_cyclogeostr_2it(:,:,ti_swot))) - ...
              (abs(U_cyclogeostr_1it(:,:,ti_swot) - U_geostr(:,:,ti_swot))));
% m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
%               log10((abs(U_cyclogeostr_3it(:,:,ti_swot) - U_cyclogeostr_2it(:,:,ti_swot))) - ...
%                      abs(U_cyclogeostr_1it(:,:,ti_swot) - U_geostr(:,:,ti_swot))));
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim([-1 1]*0.1);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
title(['$|u^{(3)}_\mathrm{cg} - u^{(2)}_\mathrm{cg}| - |u^{(1)}_\mathrm{cg} - u_\mathrm{g}|$'],'Interpreter','latex')
set(gca,'FontSize',16)
colormap(bwr)
CB = colorbar; CB.Label.String = 'Positive means more iterations are worse';

nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-.15 .15],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
              (abs(V_cyclogeostr_3it(:,:,ti_swot) - V_cyclogeostr_2it(:,:,ti_swot))) - ...
              (abs(V_cyclogeostr_1it(:,:,ti_swot) - V_geostr(:,:,ti_swot))));
% m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
%               log10((abs(V_cyclogeostr_3it(:,:,ti_swot) - V_cyclogeostr_2it(:,:,ti_swot))) - ...
%                      abs(V_cyclogeostr_1it(:,:,ti_swot) - V_geostr(:,:,ti_swot))));
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
clim([-1 1]*0.1);
% xlabel(['SWOT: $\bf{u}^{(' num2str(CG_Iteration) ')}_\mathrm{cg}$ on HFR grid'],'Interpreter','latex','Position',[2.4813e-08 0.0525 -1.0000e+16])
title(['$|v^{(3)}_\mathrm{cg} - v^{(2)}_\mathrm{cg}| - |v^{(1)}_\mathrm{cg} - v_\mathrm{g}|$'],'Interpreter','latex')
set(gca,'FontSize',16)
colormap(bwr)
CB = colorbar; CB.Label.String = 'Positive means more iterations are worse';




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
%%
% PERCENTILE = inv_prctile(DATA,VALUE)
%
% Basically the opposite of MATLAB's defaul "prctile". "inv_prctile" will
% give the approximate percentile of a value you give it compared to the
% given data set:
%
% IN:   DATA  = Vector of data with an arbitrary distribution. Note that,
%               unlike the built-in function "prctile", this function has
%               no option to format "DATA" as a matrix. If "DATA" is a
%               matrix and you want to use all elements of it, you will
%               have to reshape it.
% IN:   VALUE = Some value whose percentile in "DATA" is sought. If "VALUE"
%               is greater than max(DATA), then "PERCENTILE" = 100. If it's
%               less than min(DATA), "PERCENTILE" = 0.
%
% OUT:  PERCENTILE = The exact percentile in "DATA" which "VALUE" has if it
%                    is an element of "DATA", or the percentile it would
%                    have if it were in data.
%
%
%
% NOTE: If there are multiple occurances of "VALUE" in "DATA", then the
% value of "PERCENTILE" given will be the "most middle" possible, e.g.:
%
% inv_prctile([1 2 3 3 3],3) gives 70. See the folowing example for
% clarification for why this isn't 60:
%
% [[10:10:100]; prctile([1 2 3 3 3],[10:10:100])]
%
% Basically (I think), it's the average of 60 (because 60% are below the
% middle "3") and 80 (because 100%-80% = 20% are above the middle "3")
% 
% For short "DATA" vectors, the output might not make a lot of sense; feel
% free to tweak this function as necessary.

function PERCENTILE = inv_prctile(DATA,VALUE)

if VALUE < min(DATA)
    PERCENTILE = 0;
elseif VALUE > max(DATA)
    PERCENTILE = 100;
else
    
    if isrow(DATA)
        DATA = DATA';
    else
    end
    
    DATA = sort(DATA);
    
    if sum(DATA == VALUE) == 1 % i.e. if VALUE is in DATA once
        PERCENTILE = 100*(dsearchn(DATA,VALUE)/length(DATA)) - 50/length(DATA);
    elseif sum(DATA == VALUE) > 1 % i.e. if VALUE is in DATA more than once
        % This is janky, but it should work:
        N_occurances = sum(DATA == VALUE);
        dDATA = abs(diff(DATA));
        dDATA(dDATA == 0) = nan;
        min_diff = min(dDATA); % to prevent overshoot
        min_prc = 100*(dsearchn(DATA,VALUE)/length(DATA));
        for ii = 1:(N_occurances-1) % slightly shrink all but one occurance of VALUE
            DATA(dsearchn(DATA,VALUE)) = VALUE - 0.5*min_diff;
        end
        max_prc = 100*(dsearchn(DATA,VALUE)/length(DATA));
        PERCENTILE = mean([min_prc max_prc]) - 50/length(DATA);
    else % i.e. if VALUE is not in DATA
        PERCENTILE = 100*(dsearchn(DATA,VALUE)/length(DATA)) - 50/length(DATA);
    end
    
end
end