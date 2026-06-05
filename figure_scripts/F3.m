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
    % Optional: clear some variables that prevent later necessary
    % calculations:
    % clear U_rot_unfiltered 
else
end

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

if exist('U_cyclogeostr_Nit','var')
    warning('Skipping the *_cyclogeostr_*it calculations.')
else

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

CG_Iteration = 1;
U_cyclogeostr_Nit = U_cyclogeostr_1it;
V_cyclogeostr_Nit = V_cyclogeostr_1it;

end

%% Interpolate the SWOT GEOSTROPHIC AND CYCLOGEOSTROPHIC velocities to the HFR grid to compare in the next step:

if exist('HFR_scatter_times','var')
    error('You already calculated these variables, which take a long time to calculate.')
else
end


T_start = datenum('2023-04-01 00:00:00'); % All Cal/Val
T_end   = datenum('2023-07-11 00:00:00'); % All Cal/Val
SWOT_scatter_times = [dsearchn(NORCAL.SWOT.mean_time,T_start):...
                      dsearchn(NORCAL.SWOT.mean_time,T_end  )];
HFR_scatter_times = [];
for ti = 1:length(SWOT_scatter_times)
    ti_swot = SWOT_scatter_times(ti);
    ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, NORCAL.SWOT.mean_time(ti_swot));
    HFR_scatter_times = [HFR_scatter_times; ti_hfr];
end



if exist('Ucg_SWOT_HFRgrid_all','var')
    warning('Skipping the *cg_SWOT_HFRgrid_all calculations.')
else

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
else
    warning('Skipping the 2 and 3 iteration CG gridding calculations.')
end

toc

end

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

if exist('U_rot_unfiltered','var')
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

error('Forced stop by user')

%% G: 4-panel Histograms (GEOSTROPHIC SWOT velocity gridded to HFR grid)


% LON_LAT_lims = [-180 180 -90 90]; % all data (the entire world)
% LON_LAT_lims = [-125 -123.5 38 39.5]; % generous limits of where the southern eddy is
LON_LAT_lims = [-124.55 -123.75 38.5 39.1]; % stricter limits of where the southern eddy is
% LON_LAT_lims = [-125 -124.25 40 40.5]; % stricter limits of where the northern eddy is
% LON_LAT_lims = [-125 -123.5 38.5 40.5]; % encompassing both eddies and between

% Histogram_Time_Range = dsearchn(NORCAL.SWOT.mean_time,datenum( ...
%                            '2023-04-01 00:00:00','yyyy-mm-dd HH:MM:SS')):...
%                        dsearchn(NORCAL.SWOT.mean_time,datenum( ...
%                            '2023-07-31 00:00:00','yyyy-mm-dd HH:MM:SS'));
Histogram_Time_Range = dsearchn(NORCAL.SWOT.mean_time,datenum( ...
                           '2023-05-01 00:00:00','yyyy-mm-dd HH:MM:SS')):...
                       dsearchn(NORCAL.SWOT.mean_time,datenum( ...
                           '2023-06-01 00:00:00','yyyy-mm-dd HH:MM:SS'));

NN = 70; % bin number
HFR_Bins =     [-1:[2/NN]:1]*1.0;
SWOT_Bins =    [-1:[2/NN]:1]*1.0;
CLIM =         [0 10]; % PDF
Center = @(IN) [IN(2:end) + IN(1:[end-1])]/2;
close all
% MinMag = 0.75; % <--- This is the minimum magnitide of SWOT velocity that will be considered

clear Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g Median_Ug_over_Uhfr Median_Ug_over_Vhfr; II = 1;
for MinMag = 0 % [0:0.05:0.75] % <--- This is the minimum magnitide of SWOT velocity that will be considered
MagnitudeLimit = abs(Ug_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4), ...
                                          NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2), ...
                                          SWOT_scatter_times(Histogram_Time_Range)) + ...
                     Vg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4), ...
                                          NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2), ...
                                          SWOT_scatter_times(Histogram_Time_Range))*1i ) >= ...
                 MinMag; % <--- This is the minimum magnitide of SWOT velocity that will be considered
MagnitudeLimit = MagnitudeLimit./MagnitudeLimit;




% [R_uu,P_uu] = ...
%     corrcoef(Ug_SWOT_HFRgrid_all,Collimate(U_rot_filtered(:,:,HFR_scatter_times)),'Rows','complete');
% [R_vv,P_vv] = ...
%     corrcoef(Vg_SWOT_HFRgrid_all,Collimate(V_rot_filtered(:,:,HFR_scatter_times)),'Rows','complete');
% [R_ang,P_ang] = ...
%     corrcoef(Angleg_SWOT_HFRgrid_all,Collimate(angle(U_rot_filtered(:,:,HFR_scatter_times) + ...
%                                                   1i*V_rot_filtered(:,:,HFR_scatter_times))),'Rows','complete');
% [R_speed,P_speed] = ...
%     corrcoef(sqrt(Ug_SWOT_HFRgrid_all.^2 + Vg_SWOT_HFRgrid_all.^2),...
%              sqrt(Collimate(U_rot_filtered(:,:,HFR_scatter_times)).^2 + ...
%                   Collimate(V_rot_filtered(:,:,HFR_scatter_times)).^2),'Rows','complete');
[R_uu,P_uu,RL_uu,RH_uu] = ...
    corrcoef(Collimate(Ug_SWOT_HFRgrid_all(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_vv,P_vv,RL_vv,RH_vv] = ...
    corrcoef(Collimate(Vg_SWOT_HFRgrid_all(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_ww] = ...
    corrcoef(Collimate(Ug_SWOT_HFRgrid_all(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit) +  ...
             Collimate(1i*Vg_SWOT_HFRgrid_all(   NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit) + ...
             Collimate(   1i*V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_ang,P_ang,RL_ang,RH_ang] = ...
    corrcoef(Collimate(Angleg_SWOT_HFRgrid_all(  NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(angle(U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)) + ...
                          1i*V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range))).*MagnitudeLimit),'Rows','complete');
[R_speed,P_speed,RL_speed,RH_speed] = ...
    corrcoef(sqrt(Collimate(Ug_SWOT_HFRgrid_all( NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2 + ...
                  Collimate(Vg_SWOT_HFRgrid_all( NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2),...
             sqrt(Collimate(      U_rot_filtered(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2 + ...
                  Collimate(      V_rot_filtered(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2),'Rows','complete');

Median_Ug_over_Uhfr(II) = median([Collimate(Ug_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                                     NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)]./...
                                 [Collimate(      U_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                                     NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)],"omitmissing");
Median_Vg_over_Vhfr(II) = median([Collimate(Vg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                                     NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)]./...
                                 [Collimate(      V_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                                     NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)],"omitmissing");

% %
figure('Color','w')
tiledlayout(2,2,'TileSpacing','tight')
nexttile
% HIST_uu    = histogram2(Ug_SWOT_HFRgrid_all,Collimate(U_rot_filtered(:,:,HFR_scatter_times)),...
%                         HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
HIST_uu    = histogram2(Collimate(Ug_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        Collimate(U_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_uu = HIST_uu.Values; [X_uu,Y_uu] = meshgrid(Center(HIST_uu.XBinEdges),Center(HIST_uu.YBinEdges));
text(-0.7*max(HFR_Bins),0.7*max(HFR_Bins),['Corr. Coef. = ' num2str(R_uu(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_uu(1,2) R_uu(1,2) RH_uu(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex')
ylabel('HFR $u$ (m/s)','Interpreter','latex');
xlabel('SWOT $u_\mathrm{g}$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';

nexttile
% HIST_vv    = histogram2(Vg_SWOT_HFRgrid_all,Collimate(V_rot_filtered(:,:,HFR_scatter_times)),...
%                         HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
HIST_vv    = histogram2(Collimate(Vg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        Collimate(V_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_vv = HIST_vv.Values; [X_vv,Y_vv] = meshgrid(Center(HIST_vv.XBinEdges),Center(HIST_vv.YBinEdges));
text(-0.7*max(HFR_Bins),0.7*max(HFR_Bins),['Corr. Coef. = ' num2str(R_vv(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_vv(1,2) R_vv(1,2) RH_vv(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex')
ylabel('HFR $v$ (m/s)','Interpreter','latex');
xlabel('SWOT $v_\mathrm{g}$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';

nexttile
% HIST_aa    = histogram2(Angleg_SWOT_HFRgrid_all*180/pi,Collimate(angle(U_rot_filtered(:,:,HFR_scatter_times) + ...
%                                                                     1i*V_rot_filtered(:,:,HFR_scatter_times)))*180/pi,...
%                         [-180:5:180], [-180:5:180],'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
HIST_aa    = histogram2(Collimate(Angleg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)*180/pi,...
                        Collimate(angle(U_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)) + ...
                                     1i*V_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range))).*MagnitudeLimit)*180/pi,...
                        [-180:5:180], [-180:5:180],'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_aa = HIST_aa.Values; [X_aa,Y_aa] = meshgrid(Center(HIST_aa.XBinEdges),Center(HIST_aa.YBinEdges));
ylabel('$\theta_\mathrm{HFR}$ ($^\circ$)','Interpreter','latex');
xlabel('$\theta_\mathrm{g}$ ($^\circ$)','Interpreter','latex')
xticks([-180:90:180]); yticks([-180:90:180]);
% clim([-7.5 -4.5]); % xlim([-180 180]); ylim([-180 180])
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';

nexttile
% HIST_speed = histogram2(sqrt(Ug_SWOT_HFRgrid_all.^2 + Vg_SWOT_HFRgrid_all.^2),...
%                         Collimate(sqrt(U_rot_filtered(:,:,SWOT_scatter_times).^2 + V_rot_filtered(:,:,SWOT_scatter_times).^2)),...
%                         HFR_Bins(HFR_Bins>=0), SWOT_Bins(SWOT_Bins>=0),'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
HIST_speed = histogram2(Collimate(sqrt(Ug_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                            NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).^2 + ...
                                       Vg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                            NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).^2).*MagnitudeLimit),...
                        Collimate(sqrt(U_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                            NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).^2 + ...
                                       V_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                            NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).^2).*MagnitudeLimit),...
                        HFR_Bins(HFR_Bins>=0), SWOT_Bins(SWOT_Bins>=0),'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_speed = HIST_speed.Values; [X_speed,Y_speed] = meshgrid(Center(HIST_speed.XBinEdges),Center(HIST_speed.YBinEdges));
text(0.15*max(HFR_Bins),0.85*max(HFR_Bins),['Corr. Coef. = ' num2str(R_speed(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_speed(1,2) R_speed(1,2) RH_speed(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex')
ylabel('HFR $|\bf\mathrm{u}|$ (m/s)','Interpreter','latex');
xlabel('SWOT $|\bf{u}_\mathrm{g}|$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';
% set(gcf,'Position',[-1516         495        1176         354])
set(gcf,'Position',[-1516          84         856         765])
colormap('turbo')


figure('Color','w')
ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum('2023-05-09'));
SSHA = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
% SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-1 1]/10,...
             'latitudes',[38 42] + [-.1 .1]);
% COAST = m_gshhs_i('patch',0.8*[1 1 1]); hold on
set(gcf,'color','w')
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                  NORCAL.SWOT.lat{ti_swot},...
                  SSHA); hold on
m_plot(LON_LAT_lims([1 2 2 1 1]), LON_LAT_lims([3 3 4 4 3]),'k-','LineWidth',2)
COAST = m_gshhs_i('patch',0.8*[1 1 1],'FaceAlpha',0.5); hold on
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo); clim([-1 1]*0.67*max(abs(SSHA(:))))
set(gcf,'Position',[-798   271   319   540])

Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(II,:) = [MinMag, R_uu(1,2) + 1i*R_vv(1,2), RL_uu(1,2) + 1i*RL_vv(1,2), RH_uu(1,2) + 1i*RH_vv(1,2), R_ww(1,2)];

if II > 1
    close all
end

II = II + 1;
disp(MinMag)

end

% %%
% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/geovel_jointpdf_May2023.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')
% 
% % figure(1)
% % exportgraphics(gcf,...
% % '../figures/draft/geovel_jointpdf_May2023_demeanedHFR.pdf',...
% % 'BackgroundColor','none','ContentType','vector')
% % disp('Image saved')
% 
% figure(2)
% exportgraphics(gcf,...
% '../figures/draft/geovel_jointpdf_May2023_mapoverlay.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')

%% CG: 4-panel Histograms (CYCLOGEOSTROPHIC SWOT velocity gridded to HFR grid)


% LON_LAT_lims = [-180 180 -90 90]; % all data (the entire world)
% LON_LAT_lims = [-125 -123.5 38 39.5]; % generous limits of where the southern eddy is
LON_LAT_lims = [-124.55 -123.75 38.5 39.1]; % stricter limits of where the southern eddy is
% LON_LAT_lims = [-125 -124.25 40 40.5]; % stricter limits of where the northern eddy is
% LON_LAT_lims = [-125 -123.5 38.5 40.5]; % encompassing both eddies and between

% Histogram_Time_Range = dsearchn(NORCAL.SWOT.mean_time,datenum( ...
%                            '2023-04-01 00:00:00','yyyy-mm-dd HH:MM:SS')):...
%                        dsearchn(NORCAL.SWOT.mean_time,datenum( ...
%                            '2023-07-31 00:00:00','yyyy-mm-dd HH:MM:SS'));
Histogram_Time_Range = dsearchn(NORCAL.SWOT.mean_time,datenum( ...
                           '2023-05-01 00:00:00','yyyy-mm-dd HH:MM:SS')):...
                       dsearchn(NORCAL.SWOT.mean_time,datenum( ...
                           '2023-06-01 00:00:00','yyyy-mm-dd HH:MM:SS'));

NN = 70; % bin number
HFR_Bins =     [-1:[2/NN]:1]*1.0;
SWOT_Bins =    [-1:[2/NN]:1]*1.0;
CLIM =         [0 10]; % PDF
Center = @(IN) [IN(2:end) + IN(1:[end-1])]/2;
close all
% MinMag = 0.75; % <--- This is the minimum magnitide of SWOT velocity that will be considered

clear Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg Median_Ucg_over_Uhfr Median_Ucg_over_Vhfr; II = 1;
for MinMag = 0 % [0:0.05:0.75] % <--- This is the minimum magnitide of SWOT velocity that will be considered
MagnitudeLimit = abs(Ucg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4), ...
                                          NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2), ...
                                          SWOT_scatter_times(Histogram_Time_Range)) + ...
                     Vcg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4), ...
                                          NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2), ...
                                          SWOT_scatter_times(Histogram_Time_Range))*1i ) >= ...
                 MinMag; % <--- This is the minimum magnitide of SWOT velocity that will be considered
MagnitudeLimit = MagnitudeLimit./MagnitudeLimit;




% [R_uu,P_uu] = ...
%     corrcoef(Ucg_SWOT_HFRgrid_all,Collimate(U_rot_filtered(:,:,HFR_scatter_times)),'Rows','complete');
% [R_vv,P_vv] = ...
%     corrcoef(Vcg_SWOT_HFRgrid_all,Collimate(V_rot_filtered(:,:,HFR_scatter_times)),'Rows','complete');
% [R_ang,P_ang] = ...
%     corrcoef(Anglecg_SWOT_HFRgrid_all,Collimate(angle(U_rot_filtered(:,:,HFR_scatter_times) + ...
%                                                   1i*V_rot_filtered(:,:,HFR_scatter_times))),'Rows','complete');
% [R_speed,P_speed] = ...
%     corrcoef(sqrt(Ucg_SWOT_HFRgrid_all.^2 + Vcg_SWOT_HFRgrid_all.^2),...
%              sqrt(Collimate(U_rot_filtered(:,:,HFR_scatter_times)).^2 + ...
%                   Collimate(V_rot_filtered(:,:,HFR_scatter_times)).^2),'Rows','complete');
[R_uu,P_uu,RL_uu,RH_uu] = ...
    corrcoef(Collimate(Ucg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_vv,P_vv,RL_vv,RH_vv] = ...
    corrcoef(Collimate(Vcg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_ww] = ...
    corrcoef(Collimate(Ucg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit) +  ...
             Collimate(1i*Vcg_SWOT_HFRgrid_all(  NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit) + ...
             Collimate(   1i*V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_ang,P_ang,RL_ang,RH_ang] = ...
    corrcoef(Collimate(Anglecg_SWOT_HFRgrid_all( NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(angle(U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)) + ...
                          1i*V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range))).*MagnitudeLimit),'Rows','complete');
[R_speed,P_speed,RL_speed,RH_speed] = ...
    corrcoef(sqrt(Collimate(Ucg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2 + ...
                  Collimate(Vcg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2),...
             sqrt(Collimate(      U_rot_filtered(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2 + ...
                  Collimate(      V_rot_filtered(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2),'Rows','complete');

Median_Ucg_over_Uhfr(II) = median([Collimate(Ucg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                      NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)]./...
                                  [Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                      NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)],"omitmissing");
Median_Vcg_over_Vhfr(II) = median([Collimate(Vcg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                      NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)]./...
                                  [Collimate(      V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                      NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)],"omitmissing");


% %
figure('Color','w')
tiledlayout(2,2,'TileSpacing','tight')
nexttile
% HIST_uu    = histogram2(Ucg_SWOT_HFRgrid_all,Collimate(U_rot_filtered(:,:,HFR_scatter_times)),...
%                         HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
HIST_uu    = histogram2(Collimate(Ucg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        Collimate(U_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_uu = HIST_uu.Values; [X_uu,Y_uu] = meshgrid(Center(HIST_uu.XBinEdges),Center(HIST_uu.YBinEdges));
text(-0.7*max(HFR_Bins),0.7*max(HFR_Bins),['Corr. Coef. = ' num2str(R_uu(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_uu(1,2) R_uu(1,2) RH_uu(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex')
ylabel('HFR $u$ (m/s)','Interpreter','latex');
xlabel('SWOT $u_\mathrm{cg}$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';

nexttile
% HIST_vv    = histogram2(Vcg_SWOT_HFRgrid_all,Collimate(V_rot_filtered(:,:,HFR_scatter_times)),...
%                         HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
HIST_vv    = histogram2(Collimate(Vcg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        Collimate(V_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_vv = HIST_vv.Values; [X_vv,Y_vv] = meshgrid(Center(HIST_vv.XBinEdges),Center(HIST_vv.YBinEdges));
text(-0.7*max(HFR_Bins),0.7*max(HFR_Bins),['Corr. Coef. = ' num2str(R_vv(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_vv(1,2) R_vv(1,2) RH_vv(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex')
ylabel('HFR $v$ (m/s)','Interpreter','latex');
xlabel('SWOT $v_\mathrm{cg}$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';

nexttile
% HIST_aa    = histogram2(Anglecg_SWOT_HFRgrid_all*180/pi,Collimate(angle(U_rot_filtered(:,:,HFR_scatter_times) + ...
%                                                                     1i*V_rot_filtered(:,:,HFR_scatter_times)))*180/pi,...
%                         [-180:5:180], [-180:5:180],'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
HIST_aa    = histogram2(Collimate(Anglecg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)*180/pi,...
                        Collimate(angle(U_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)) + ...
                                     1i*V_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range))).*MagnitudeLimit)*180/pi,...
                        [-180:5:180], [-180:5:180],'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_aa = HIST_aa.Values; [X_aa,Y_aa] = meshgrid(Center(HIST_aa.XBinEdges),Center(HIST_aa.YBinEdges));
ylabel('$\theta_\mathrm{HFR}$ ($^\circ$)','Interpreter','latex');
xlabel('$\theta_\mathrm{cg}$ ($^\circ$)','Interpreter','latex')
xticks([-180:90:180]); yticks([-180:90:180]);
% clim([-7.5 -4.5]); % xlim([-180 180]); ylim([-180 180])
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';

nexttile
% HIST_speed = histogram2(sqrt(Ucg_SWOT_HFRgrid_all.^2 + Vcg_SWOT_HFRgrid_all.^2),...
%                         Collimate(sqrt(U_rot_filtered(:,:,SWOT_scatter_times).^2 + V_rot_filtered(:,:,SWOT_scatter_times).^2)),...
%                         HFR_Bins(HFR_Bins>=0), SWOT_Bins(SWOT_Bins>=0),'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
HIST_speed = histogram2(Collimate(sqrt(Ucg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                            NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).^2 + ...
                                       Vcg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                            NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).^2).*MagnitudeLimit),...
                        Collimate(sqrt(U_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                            NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).^2 + ...
                                       V_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                            NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).^2).*MagnitudeLimit),...
                        HFR_Bins(HFR_Bins>=0), SWOT_Bins(SWOT_Bins>=0),'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_speed = HIST_speed.Values; [X_speed,Y_speed] = meshgrid(Center(HIST_speed.XBinEdges),Center(HIST_speed.YBinEdges));
text(0.15*max(HFR_Bins),0.85*max(HFR_Bins),['Corr. Coef. = ' num2str(R_speed(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_speed(1,2) R_speed(1,2) RH_speed(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex')
ylabel('HFR $|\bf\mathrm{u}|$ (m/s)','Interpreter','latex');
xlabel('SWOT $|\bf{u}_\mathrm{cg}|$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';
% set(gcf,'Position',[-1516         495        1176         354])
set(gcf,'Position',[-1516          84         856         765])
colormap('turbo')


figure('Color','w')
ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum('2023-05-09'));
SSHA = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
% SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-1 1]/10,...
             'latitudes',[38 42] + [-.1 .1]);
% COAST = m_gshhs_i('patch',0.8*[1 1 1]); hold on
set(gcf,'color','w')
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                  NORCAL.SWOT.lat{ti_swot},...
                  SSHA); hold on
m_plot(LON_LAT_lims([1 2 2 1 1]), LON_LAT_lims([3 3 4 4 3]),'k-','LineWidth',2)
COAST = m_gshhs_i('patch',0.8*[1 1 1],'FaceAlpha',0.5); hold on
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo); clim([-1 1]*0.67*max(abs(SSHA(:))))
set(gcf,'Position',[-798   271   319   540])

Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(II,:) = [MinMag, R_uu(1,2) + 1i*R_vv(1,2), RL_uu(1,2) + 1i*RL_vv(1,2), RH_uu(1,2) + 1i*RH_vv(1,2), R_ww(1,2)];

if II > 1
    close all
end

II = II + 1;
disp(MinMag)

end

% %%
% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/cyclogeovel_jointpdf_May2023.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')
% 
% % figure(1)
% % exportgraphics(gcf,...
% % '../figures/draft/cyclogeovel_jointpdf_May2023_demeanedHFR.pdf',...
% % 'BackgroundColor','none','ContentType','vector')
% % disp('Image saved')
% 
% figure(2)
% exportgraphics(gcf,...
% '../figures/draft/cyclogeovel_jointpdf_May2023_mapoverlay.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')



%% G: 3-panel Histograms (GEOSTROPHIC SWOT velocity gridded to HFR grid)


LON_LAT_lims = [-124.55 -123.75 38.5 39.1]; % stricter limits of where the southern eddy is
Histogram_Time_Range = dsearchn(NORCAL.SWOT.mean_time,datenum( ...
                           '2023-05-01 00:00:00','yyyy-mm-dd HH:MM:SS')):...
                       dsearchn(NORCAL.SWOT.mean_time,datenum( ...
                           '2023-06-01 00:00:00','yyyy-mm-dd HH:MM:SS'));

NN = 70; % bin number
HFR_Bins =     [-1:[2/NN]:1]*1.0;
SWOT_Bins =    [-1:[2/NN]:1]*1.0;
CLIM =         [0 10]; % PDF
Center = @(IN) [IN(2:end) + IN(1:[end-1])]/2;
close all
% MinMag = 0.75; % <--- This is the minimum magnitide of SWOT velocity that will be considered

clear Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g Median_Ug_over_Uhfr Median_Ug_over_Vhfr; II = 1;
for MinMag = 0 % [0:0.05:0.75] % <--- This is the minimum magnitide of SWOT velocity that will be considered
MagnitudeLimit = abs(Ug_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4), ...
                                          NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2), ...
                                          SWOT_scatter_times(Histogram_Time_Range)) + ...
                     Vg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4), ...
                                          NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2), ...
                                          SWOT_scatter_times(Histogram_Time_Range))*1i ) >= ...
                 MinMag; % <--- This is the minimum magnitide of SWOT velocity that will be considered
MagnitudeLimit = MagnitudeLimit./MagnitudeLimit;




[R_uu,P_uu,RL_uu,RH_uu] = ...
    corrcoef(Collimate(Ug_SWOT_HFRgrid_all(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_vv,P_vv,RL_vv,RH_vv] = ...
    corrcoef(Collimate(Vg_SWOT_HFRgrid_all(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_ww] = ...
    corrcoef(Collimate(Ug_SWOT_HFRgrid_all(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit) +  ...
             Collimate(1i*Vg_SWOT_HFRgrid_all(   NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit) + ...
             Collimate(   1i*V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_ang,P_ang,RL_ang,RH_ang] = ...
    corrcoef(Collimate(Angleg_SWOT_HFRgrid_all(  NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(angle(U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)) + ...
                          1i*V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range))).*MagnitudeLimit),'Rows','complete');
[R_speed,P_speed,RL_speed,RH_speed] = ...
    corrcoef(sqrt(Collimate(Ug_SWOT_HFRgrid_all( NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2 + ...
                  Collimate(Vg_SWOT_HFRgrid_all( NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2),...
             sqrt(Collimate(      U_rot_filtered(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2 + ...
                  Collimate(      V_rot_filtered(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2),'Rows','complete');

Median_Ug_over_Uhfr(II) = median([Collimate(Ug_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                                     NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)]./...
                                 [Collimate(      U_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                                     NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)],"omitmissing");
Median_Vg_over_Vhfr(II) = median([Collimate(Vg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                                     NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)]./...
                                 [Collimate(      V_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                                     NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)],"omitmissing");

% %
figure('Color','w')
tiledlayout(2,2,'TileSpacing','tight')
nexttile(3)
HIST_uu    = histogram2(Collimate(Ug_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        Collimate(U_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none'); hold on
plot([-1 0.5],[-1 0.5],'--','LineWidth',1.5,'Color',[1 1 1]*0.6);
PDF_uu = HIST_uu.Values; [X_uu,Y_uu] = meshgrid(Center(HIST_uu.XBinEdges),Center(HIST_uu.YBinEdges));
text(0,0.7*max(HFR_Bins),['Corr. Coef. = ' num2str(R_uu(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_uu(1,2) R_uu(1,2) RH_uu(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex','HorizontalAlignment','center')
ylabel('HFR $u$ (m/s)','Interpreter','latex');
xlabel('SWOT $u_\mathrm{g}$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';
axis equal; axis tight

nexttile(4)
HIST_vv    = histogram2(Collimate(Vg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        Collimate(V_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_vv = HIST_vv.Values; [X_vv,Y_vv] = meshgrid(Center(HIST_vv.XBinEdges),Center(HIST_vv.YBinEdges)); hold on
plot([-1 0.5],[-1 0.5],'--','LineWidth',1.5,'Color',[1 1 1]*0.6);
text(0,0.7*max(HFR_Bins),['Corr. Coef. = ' num2str(R_vv(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_vv(1,2) R_vv(1,2) RH_vv(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex','HorizontalAlignment','center')
ylabel('HFR $v$ (m/s)','Interpreter','latex');
xlabel('SWOT $v_\mathrm{g}$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';
axis equal; axis tight

nexttile(2)
HIST_aa    = histogram2(Collimate(Angleg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)*180/pi,...
                        Collimate(angle(U_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)) + ...
                                     1i*V_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range))).*MagnitudeLimit)*180/pi,...
                        [-180:5:180], [-180:5:180],'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_aa = HIST_aa.Values; [X_aa,Y_aa] = meshgrid(Center(HIST_aa.XBinEdges),Center(HIST_aa.YBinEdges));
ylabel('$\theta_\mathrm{HFR}$ ($^\circ$)','Interpreter','latex');
xlabel('$\theta_\mathrm{g}$ ($^\circ$)','Interpreter','latex')
xticks([-180:90:180]); yticks([-180:90:180]);
% clim([-7.5 -4.5]); % xlim([-180 180]); ylim([-180 180])
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';
axis equal; axis tight

colormap('turbo')

nexttile(1)
ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum('2023-05-09'));
SSHA = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
% SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-1 1]/10,...
             'latitudes',[38 40] + [-.1 .1]);
% COAST = m_gshhs_i('patch',0.8*[1 1 1]); hold on
set(gcf,'color','w')
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                  NORCAL.SWOT.lat{ti_swot},...
                  SSHA); hold on
m_plot(LON_LAT_lims([1 2 2 1 1]), LON_LAT_lims([3 3 4 4 3]),'k-','LineWidth',2)
COAST = m_gshhs_i('patch',0.8*[1 1 1],'FaceAlpha',0.5); hold on
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo); clim([-1 1]*0.67*max(abs(SSHA(:))))

set(gcf,'Position',[-816    88   817   735])

Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(II,:) = [MinMag, R_uu(1,2) + 1i*R_vv(1,2), RL_uu(1,2) + 1i*RL_vv(1,2), RH_uu(1,2) + 1i*RH_vv(1,2), R_ww(1,2)];

if II > 1
    close all
end

II = II + 1;
disp(MinMag)

end

% %%
% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/geovel_jointpdf_May2023.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')
% 
% % figure(1)
% % exportgraphics(gcf,...
% % '../figures/draft/geovel_jointpdf_May2023_demeanedHFR.pdf',...
% % 'BackgroundColor','none','ContentType','vector')
% % disp('Image saved')
% 
% figure(2)
% exportgraphics(gcf,...
% '../figures/draft/geovel_jointpdf_May2023_mapoverlay.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')

%% CG: 3-panel Histograms (CYCLOGEOSTROPHIC SWOT velocity gridded to HFR grid)


LON_LAT_lims = [-124.55 -123.75 38.5 39.1]; % stricter limits of where the southern eddy is
Histogram_Time_Range = dsearchn(NORCAL.SWOT.mean_time,datenum( ...
                           '2023-05-01 00:00:00','yyyy-mm-dd HH:MM:SS')):...
                       dsearchn(NORCAL.SWOT.mean_time,datenum( ...
                           '2023-06-01 00:00:00','yyyy-mm-dd HH:MM:SS'));

NN = 70; % bin number
HFR_Bins =     [-1:[2/NN]:1]*1.0;
SWOT_Bins =    [-1:[2/NN]:1]*1.0;
CLIM =         [0 10]; % PDF
Center = @(IN) [IN(2:end) + IN(1:[end-1])]/2;
close all
% MinMag = 0.75; % <--- This is the minimum magnitide of SWOT velocity that will be considered

clear Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg Median_Ucg_over_Uhfr Median_Ucg_over_Vhfr; II = 1;
for MinMag = 0 % [0:0.05:0.75] % <--- This is the minimum magnitide of SWOT velocity that will be considered
MagnitudeLimit = abs(Ucg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4), ...
                                          NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2), ...
                                          SWOT_scatter_times(Histogram_Time_Range)) + ...
                     Vcg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4), ...
                                          NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2), ...
                                          SWOT_scatter_times(Histogram_Time_Range))*1i ) >= ...
                 MinMag; % <--- This is the minimum magnitide of SWOT velocity that will be considered
MagnitudeLimit = MagnitudeLimit./MagnitudeLimit;




[R_uu,P_uu,RL_uu,RH_uu] = ...
    corrcoef(Collimate(Ucg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_vv,P_vv,RL_vv,RH_vv] = ...
    corrcoef(Collimate(Vcg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_ww] = ...
    corrcoef(Collimate(Ucg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit) +  ...
             Collimate(1i*Vcg_SWOT_HFRgrid_all(  NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit) + ...
             Collimate(   1i*V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),'Rows','complete');
[R_ang,P_ang,RL_ang,RH_ang] = ...
    corrcoef(Collimate(Anglecg_SWOT_HFRgrid_all( NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
             Collimate(angle(U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)) + ...
                          1i*V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range))).*MagnitudeLimit),'Rows','complete');
[R_speed,P_speed,RL_speed,RH_speed] = ...
    corrcoef(sqrt(Collimate(Ucg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2 + ...
                  Collimate(Vcg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2),...
             sqrt(Collimate(      U_rot_filtered(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2 + ...
                  Collimate(      V_rot_filtered(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                 NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit).^2),'Rows','complete');

Median_Ucg_over_Uhfr(II) = median([Collimate(Ucg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                      NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)]./...
                                  [Collimate(      U_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                      NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)],"omitmissing");
Median_Vcg_over_Vhfr(II) = median([Collimate(Vcg_SWOT_HFRgrid_all(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                      NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)]./...
                                  [Collimate(      V_rot_filtered(     NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                      NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)],"omitmissing");


% %
figure('Color','w')
tiledlayout(2,2,'TileSpacing','tight')
nexttile(3)
HIST_uu    = histogram2(Collimate(Ucg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        Collimate(U_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none'); hold on
plot([-1 0.5],[-1 0.5],'--','LineWidth',1.5,'Color',[1 1 1]*0.6);
PDF_uu = HIST_uu.Values; [X_uu,Y_uu] = meshgrid(Center(HIST_uu.XBinEdges),Center(HIST_uu.YBinEdges));
text(0,0.7*max(HFR_Bins),['Corr. Coef. = ' num2str(R_uu(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_uu(1,2) R_uu(1,2) RH_uu(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex','HorizontalAlignment','center')
ylabel('HFR $u$ (m/s)','Interpreter','latex');
xlabel('SWOT $u_\mathrm{cg}$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';

nexttile(4)
HIST_vv    = histogram2(Collimate(Vcg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        Collimate(V_rot_filtered(      NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                       NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)).*MagnitudeLimit),...
                        HFR_Bins, SWOT_Bins,'Normalization','pdf','DisplayStyle','tile','EdgeColor','none'); hold on
plot([-1 0.5],[-1 0.5],'--','LineWidth',1.5,'Color',[1 1 1]*0.6);
PDF_vv = HIST_vv.Values; [X_vv,Y_vv] = meshgrid(Center(HIST_vv.XBinEdges),Center(HIST_vv.YBinEdges));
text(0,0.7*max(HFR_Bins),['Corr. Coef. = ' num2str(R_vv(1,2),'%.3f') ' $\pm$ ' num2str(mean(diff([RL_vv(1,2) R_vv(1,2) RH_vv(1,2)])),'%.3f')],'FontSize',18,'Interpreter','latex','HorizontalAlignment','center')
ylabel('HFR $v$ (m/s)','Interpreter','latex');
xlabel('SWOT $v_\mathrm{cg}$ (m/s)','Interpreter','latex')
clim(CLIM)
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';

nexttile(2)
HIST_aa    = histogram2(Collimate(Anglecg_SWOT_HFRgrid_all(NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),SWOT_scatter_times(Histogram_Time_Range)).*MagnitudeLimit)*180/pi,...
                        Collimate(angle(U_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range)) + ...
                                     1i*V_rot_filtered(    NORCAL.HFR.lat>LON_LAT_lims(3) & NORCAL.HFR.lat<LON_LAT_lims(4),...
                                                           NORCAL.HFR.lon>LON_LAT_lims(1) & NORCAL.HFR.lon<LON_LAT_lims(2),HFR_scatter_times(Histogram_Time_Range))).*MagnitudeLimit)*180/pi,...
                        [-180:5:180], [-180:5:180],'Normalization','pdf','DisplayStyle','tile','EdgeColor','none');
PDF_aa = HIST_aa.Values; [X_aa,Y_aa] = meshgrid(Center(HIST_aa.XBinEdges),Center(HIST_aa.YBinEdges));
ylabel('$\theta_\mathrm{HFR}$ ($^\circ$)','Interpreter','latex');
xlabel('$\theta_\mathrm{cg}$ ($^\circ$)','Interpreter','latex')
xticks([-180:90:180]); yticks([-180:90:180]);
% clim([-7.5 -4.5]); % xlim([-180 180]); ylim([-180 180])
CB = colorbar;
CB.FontSize = 14; set(gca,'FontSize',14);
% CB.Label.String = 'PDF';

colormap('turbo')

nexttile(1)
ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum('2023-05-09'));
SSHA = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
       [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
% SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(NORCAL.HFR.lon(:)) -123.] + [-1 1]/10,...
             'latitudes',[38 40] + [-.1 .1]);
% COAST = m_gshhs_i('patch',0.8*[1 1 1]); hold on
set(gcf,'color','w')
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                  NORCAL.SWOT.lat{ti_swot},...
                  SSHA); hold on
m_plot(LON_LAT_lims([1 2 2 1 1]), LON_LAT_lims([3 3 4 4 3]),'k-','LineWidth',2)
COAST = m_gshhs_i('patch',0.8*[1 1 1],'FaceAlpha',0.5); hold on
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo); clim([-1 1]*0.67*max(abs(SSHA(:))))

set(gcf,'Position',[-816    88   817   735])

Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(II,:) = [MinMag, R_uu(1,2) + 1i*R_vv(1,2), RL_uu(1,2) + 1i*RL_vv(1,2), RH_uu(1,2) + 1i*RH_vv(1,2), R_ww(1,2)];

if II > 1
    close all
end

II = II + 1;
disp(MinMag)

end

% %%
% figure(1)
% exportgraphics(gcf,...
% '../figures/draft/cyclogeovel_jointpdf_May2023.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')
% 
% % figure(1)
% % exportgraphics(gcf,...
% % '../figures/draft/cyclogeovel_jointpdf_May2023_demeanedHFR.pdf',...
% % 'BackgroundColor','none','ContentType','vector')
% % disp('Image saved')
% 
% figure(2)
% exportgraphics(gcf,...
% '../figures/draft/cyclogeovel_jointpdf_May2023_mapoverlay.pdf',...
% 'BackgroundColor','none','ContentType','vector')
% disp('Image saved')




%% Plot CorrCoef as a function of the minimum SWOT velocity magnitude we examine:

close all

% figure("Color","w")
% tiledlayout(1,2,"TileSpacing","tight")
% nexttile
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,2)),'b.-'); hold on
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,3)),'b--');
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,4)),'b--');
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,2)),'r.-');
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,3)),'r--');
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,4)),'r--');
% ylim([0 1])
% title('Geostrophic')
% nexttile
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,2)),'b.-'); hold on
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,3)),'b--');
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,4)),'b--');
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,2)),'r.-');
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,3)),'r--');
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,4)),'r--');
% ylim([0 1])
% title('CycloGeostrophic')

% figure("Color","w")
% tiledlayout(1,2,"TileSpacing","tight")
% nexttile
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,1),Median_Ug_over_Uhfr','.-b');hold on
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_g(:,1),Median_Vg_over_Vhfr','.-r');
% ylim([0 4])
% title('Geostrophic')
% nexttile
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,1),Median_Ucg_over_Uhfr','.-b');hold on
% plot(Mag_Ruuvv_RLuuvv_RHuuvv_cg(:,1),Median_Vcg_over_Vhfr','.-r');
% ylim([0 4])
% title('CycloGeostrophic')


figure("Color","w")
tiledlayout(1,2,"TileSpacing","tight")
nexttile
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,2)),'b.-'); hold on
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,3)),'b--','HandleVisibility','off');
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,4)),'b--','HandleVisibility','off');
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,2)),'g.-');
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,3)),'g--','HandleVisibility','off');
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),real(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,4)),'g--','HandleVisibility','off');
legend('G vs. HFR','CG vs. HFR')
ylim([0 1])
title('U')
nexttile
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,2)),'r.-'); hold on
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,3)),'r--','HandleVisibility','off');
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,4)),'r--','HandleVisibility','off');
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,2)),'m.-');
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,3)),'m--','HandleVisibility','off');
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),imag(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,4)),'m--','HandleVisibility','off');
legend('G vs. HFR','CG vs. HFR')
ylim([0 1])
title('V')

figure("Color","w")
tiledlayout(1,2,"TileSpacing","tight")
nexttile
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),Median_Ug_over_Uhfr','.-b');hold on
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),Median_Ucg_over_Uhfr','.-g')
legend('G/HFR','CG/HFR')
ylim([0 4])
title('U')
nexttile
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),Median_Vg_over_Vhfr','.-r');;hold on
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),Median_Vcg_over_Vhfr','.-m');
legend('G/HFR','CG/HFR')
ylim([0 4])
title('V')

figure("Color","w")
tiledlayout(1,2,"TileSpacing","tight")
nexttile
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),abs(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,5)),'b.-'); hold on
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),abs(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,5)),'g.-')
title('Full Corr.Coef.')
legend('G','CG')
nexttile
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,1),angle(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_g(:,5))*180/pi,'b.-'); hold on
plot(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,1),angle(Mag_Ruuvv_RLuuvv_RHuuvv_Rww_cg(:,5))*180/pi,'g.-')
title('Angle')
legend('G','CG')


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