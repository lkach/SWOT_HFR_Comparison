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

%% Load Radial Data

% Antenna_ID = 'TRIN';
% Antenna_ID = 'GCVE';
% Antenna_ID = 'SMOA';
% Antenna_ID = 'SHEL';
Antenna_ID = 'BRAG';
% Antenna_ID = 'PAFS';

PATH = './data/HFR_radials/';
% DIR = dir([PATH '*.ruv']);
DIR = dir([PATH 'Radial_data_allcalvaltimes_' Antenna_ID '.mat']);
% cd(PATH)

% DATA = [];
% TimeStamp = [];
% for ii = 1:length(DIR)
%     DATA_struct = importdata(DIR(ii).name,' ',55);
%     TimeStamp_ii = datenum( ...
%         replace( ...
%                 replace(DATA_struct.textdata{7},{'%TimeStamp: ','  '},{'','$'}), ...
%                 {' ','$'},{'-',' '}) , 'yyyy-mm-dd HH');
%     DATA_ii = [DATA_struct.data, TimeStamp_ii*ones(size(DATA_struct.data,1),1)];
%     DATA = [DATA ; DATA_ii];
%     TimeStamp = [TimeStamp; TimeStamp_ii];
%     disp(ii)
% end
% Variable_names = [DATA_struct.textdata{end-1} '   Time(UTC)'];
% Origin_antenna = str2num(replace(replace(DATA_struct.textdata{10},{'%Origin:  '},{'1i*'}),{' '},{' + '}));

Radial = load([DIR.folder '/' DIR.name],'LON','LAT','VEL','TIME',...
              'Origin_antenna','Variable_names','Unique_LONLAT','mean_VEL');
DATA = []; % Placeholder

if strcmp(Antenna_ID,'GCVE')
    % For some reason, GCVE didn't save the correct origin, even though
    % the data seem to be in the correct locations, so instead save it
    % manually:
    Radial.Origin_antenna = 1i*38.5671833 - 123.3315500;
else
end

% Unique locations that the radials sample:
Radial.Unique_LONLAT = Radial.LON + Radial.LAT*1i;
Radial.Unique_LONLAT = unique(Radial.Unique_LONLAT);

% %% Mean radial velocity for HFR calculated
% % error % takes a long time to run (287.975086 seconds ~ almost 5 minutes)
% Radial.mean_VEL = nan(size(Radial.Unique_LONLAT));
% tic
% for ii = 1:length(Radial.Unique_LONLAT)
%     Radial.mean_VEL(ii) = mean(Radial.VEL([Radial.LON + 1i*Radial.LAT] == Radial.Unique_LONLAT(ii)),'omitnan');
% end
% toc
% Unique_LONLAT = Radial.Unique_LONLAT; mean_VEL = Radial.mean_VEL;
% save([DIR.folder '/' DIR.name],'Unique_LONLAT','mean_VEL','-append')
% %%

% direcitonal vectors:
clear rHat;
rHat(:,1) = [real(Radial.Origin_antenna - [Radial.LON + 1i*Radial.LAT])]./...
              abs(Radial.Origin_antenna - [Radial.LON + 1i*Radial.LAT]);
rHat(:,2) = [imag(Radial.Origin_antenna - [Radial.LON + 1i*Radial.LAT])]./...
              abs(Radial.Origin_antenna - [Radial.LON + 1i*Radial.LAT]);

warning('Uncomment these if you have already calculated these quantities:')
load([PATH 'Dist_to_nearest_radial_' Antenna_ID '.mat'])
load([PATH 'Radial_vel_comparison_' Antenna_ID '.mat'])

% %% Load gridded HFR for histogram comparison only

T0 = datenum('2012-01-01 00:00:00'); % HFR 0-time
SWOT_calval_t0 = '01-Apr-2023 07:00:00';
SWOT_calval_tf = '10-Jul-2023 17:00:00';
NORCAL.HFR.lat = ncread('./data/HFR_NorCal_CalVal_2023.nc','lat');
NORCAL.HFR.lon = ncread('./data/HFR_NorCal_CalVal_2023.nc','lon');
NORCAL.HFR.u = ncread('./data/HFR_NorCal_CalVal_2023.nc','u');
NORCAL.HFR.v = ncread('./data/HFR_NorCal_CalVal_2023.nc','v');
NORCAL.HFR.time = ncread('./data/HFR_NorCal_CalVal_2023.nc','time');
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
NORCAL.HFR = rmfield(NORCAL.HFR,'v_mean_2023');
NORCAL.HFR = rmfield(NORCAL.HFR,'u_mean_2023');


Radial.ANGLE  = angle([Radial.LON - real(Radial.Origin_antenna)] + 1i*[Radial.LAT - imag(Radial.Origin_antenna)])*180/pi;
Radial.R_dist = abs([Radial.LON - real(Radial.Origin_antenna)] + 1i*[Radial.LAT - imag(Radial.Origin_antenna)]);



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
                          SSHA, 6); % 6 -> 13x13 box, i.e. [26 km]^2
    % To get geostrophic velocity from inputting surface height:
    U_geostr(:,:,ti_swot) = -(gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,3))/111000;
    V_geostr(:,:,ti_swot) =  (gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,2))./[111000*cosd(NORCAL.SWOT.lat{ti_swot})];


    disp(100*ti_swot/length(NORCAL.SWOT.ssha_karin_2))
end

% % % BROKEN FOR SOME REASON:


% % % % Recall that at the surface, geostrophic velocity is:
% % % % u_g = -(g/f)(d\eta/dy)
% % % % v_g = +(g/f)(d\eta/dx)
% % % % where both derivatives are partial derivatives.
% % % % Simple textbook citation:
% % % % https://uw.pressbooks.pub/ocean285/chapter/geostrophic-balance/
% % % 
% % % close all
% % % 
% % % SWOT_geostr_vel_calculated = true;
% % % gg = 9.807; % m s^-2
% % % GaussianKernelSTD_meters = 4000;
% % % PIXELS = 4;
% % % % Number of pixels in one direction for the quadratic fit, e.g.:
% % % % PIXELS = 6 -> 13x13 box, i.e. [26 km]^2
% % % 
% % % U_geostr = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% % % V_geostr = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% % % % U_geostr_ = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% % % % V_geostr_ = nan([size(NORCAL.SWOT.ssha_karin_2{1}) length(NORCAL.SWOT.ssha_karin_2)]);
% % % for ti_swot = 1:length(NORCAL.SWOT.ssha_karin_2)
% % % 
% % %     % % % Calculate ssha to obtain geostrophic velocity:
% % %     SSHA = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
% % %            [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
% % %            [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
% % % 
% % % 
% % %     % % % % % Smoothing is necessary to reduce high-k noise that can really throw
% % %     % % % % % off the velocity when taking the first derivative.
% % % 
% % % 
% % %     % % % Less good Gaussian smoothing procedures:
% % %     % SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');
% % %     % SSHA = smoothdata2(SSHA,'movmedian',7,'omitnan');
% % %     % SSHA = smoothdata2_interpnadirgap(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, SSHA, 2000,GaussianKernelSTD_meters);
% % % 
% % % 
% % %     % % % Isotropic Gaussian blur (best Gaussian blur version so far):
% % %     % SSHA = bettergaussiansmooth2(SSHA,2000,GaussianKernelSTD_meters);
% % %     % SSHA = SSHA .* ...
% % %     %    [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
% % %     %    [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
% % %     % % % Built-in gradientm function (requires mapping toolbox, better):
% % %     % [~,~,gradSSHA_y, gradSSHA_x] = gradientm(NORCAL.SWOT.lat{ti_swot},NORCAL.SWOT.lon{ti_swot},SSHA);
% % %     % SWOT_ug = -gg*gradSSHA_y./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
% % %     %     SWOT_ug = SWOT_ug .* [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}]./[~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}];
% % %     % SWOT_vg =  gg*gradSSHA_x./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
% % %     %     SWOT_vg = SWOT_vg .* [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}]./[~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}];
% % %     % 
% % %     % U_geostr(:,:,ti_swot) = SWOT_ug;
% % %     % V_geostr(:,:,ti_swot) = SWOT_vg;
% % % 
% % % 
% % %     % % % My custom gradient function (not as good as gradientm):
% % %     % gradSSHA = gradient_2D(SSHA,NORCAL.SWOT.lon{ti_swot}.*cosd(NORCAL.SWOT.lat{ti_swot}),...
% % %     %                             NORCAL.SWOT.lat{ti_swot})*[1/111000];
% % %     % gradSSHA_x = gradSSHA(:,:,1); gradSSHA_y = gradSSHA(:,:,2);
% % %     % SWOT_ug = -gg*gradSSHA_y./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
% % %     %     SWOT_ug = SWOT_ug .* [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}]./[~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}];
% % %     % SWOT_vg =  gg*gradSSHA_x./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot});
% % %     %     SWOT_vg = SWOT_vg .* [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}]./[~NORCAL.SWOT.ssha_karin_2_qual{ti_swot}];
% % %     % 
% % %     % U_geostr_(:,:,ti_swot) = SWOT_ug;
% % %     % V_geostr_(:,:,ti_swot) = SWOT_vg;
% % % 
% % %     % % % Obtain geostrophic velocity directly using the quadratic fit
% % %     % % % method of Tranchant et al. (2025):
% % %     [AA,~] = fit_6term_2d(NORCAL.SWOT.lon{ti_swot}, ...
% % %                           NORCAL.SWOT.lat{ti_swot}, ...
% % %                           SSHA, PIXELS); % 6 -> 13x13 box, i.e. [26 km]^2
% % %                           % PIXELS -> [2*PIXELS + 1]x[2*PIXELS + 1] box
% % %     % To get geostrophic velocity from inputting surface height:
% % %     U_geostr(:,:,ti_swot) = -(gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,3))/111000;
% % %     V_geostr(:,:,ti_swot) =  (gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot})) .* squeeze(AA(:,:,2))./[111000*cosd(NORCAL.SWOT.lat{ti_swot})];
% % % 
% % % 
% % %     disp(100*ti_swot/length(NORCAL.SWOT.ssha_karin_2))
% % % end

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

error('Forced stop by user')

%% Radial unit vectors in direction from SWOT grid points toward antenna

SWOT_radial_unit_complex = ...
    [   -[NORCAL.SWOT.lon{1} + 1i*NORCAL.SWOT.lat{1}] + [Radial.Origin_antenna]]./...
    abs(-[NORCAL.SWOT.lon{1} + 1i*NORCAL.SWOT.lat{1}] + [Radial.Origin_antenna]);
SWOT_radial_unit_x = real(SWOT_radial_unit_complex);
SWOT_radial_unit_y = imag(SWOT_radial_unit_complex);
clear SWOT_radial_unit_complex


figure('Color','w')
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(Radial.LON) max(Radial.LON)] + [-1 1]/1,...
             'latitudes', [min(Radial.LAT) max(Radial.LAT)] + [-1 1]/1);
% set(gcf,'color',[1 1 1]*0.5);
hold on
VecScale = 0.03;
m_quiver(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},VecScale*SWOT_radial_unit_x,VecScale*SWOT_radial_unit_y,0,'k');
m_plot(real(Radial.Origin_antenna), imag(Radial.Origin_antenna), 'r*')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
CB = colorbar;
CB.Label.String = '{\bf r}'; CB.Label.FontSize = 20; clim([-1 1]*1)
set(gca,'FontSize',16)



%% visualize a snapshot

User_chosen_date = datenum('2023-05-09 02:00:00');
User_chosen_date = datenum('2023-05-31 02:00:00');
time_IND = Radial.TIME == Radial.TIME(dsearchn(Radial.TIME,User_chosen_date));
% time_IND = time_IND & [DATA(:,6) == 0 | DATA(:,7) == 0];

close all
figure('Color','w')

M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(Radial.LON) max(Radial.LON)] + [-1 1]/-10,...
             'latitudes', [min(Radial.LAT) max(Radial.LAT)] + [-1 1]/-10);
set(gcf,'color',[1 1 1]*0.5); hold on
VecScale = 0.001;

% m_quiver(Radial.LON(time_IND),...
%          Radial.LAT(time_IND),...
%          VecScale*DATA(time_IND,3), ...
%          VecScale*DATA(time_IND,4),0,'k')
% m_scatter(Radial.LON(time_IND),Radial.LAT(time_IND),10,... magnitude
          % 0.01*sqrt(DATA(time_IND,3).^2 + DATA(time_IND,4).^2),'filled'); colormap(turbo)

m_scatter(Radial.LON(time_IND),Radial.LAT(time_IND),10,... magnitude and direction
          0.01*Radial.VEL(time_IND),'filled'); colormap(bwr)
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
CB = colorbar;
CB.Label.String = 'HFR: |u| (m s^{-1})'; CB.Label.FontSize = 20; clim([-1 1]*1)
set(gca,'FontSize',16)
title(['HFR time = ' datestr(Radial.TIME(dsearchn(Radial.TIME,User_chosen_date)))])


%% Overlay radials over SWOT SSHA

User_chosen_date = datenum('2023-04-01 02:00:00');
% User_chosen_date = datenum('2023-05-08 02:00:00');
% User_chosen_date = datenum('2023-05-09 02:00:00');
User_chosen_date = datenum('2023-05-29 02:00:00');
% User_chosen_date = datenum('2023-05-31 02:00:00');
% User_chosen_date = datenum('2023-06-09 02:00:00');
% User_chosen_date = datenum('2023-06-15 02:00:00');
time_IND = Radial.TIME == Radial.TIME(dsearchn(Radial.TIME,User_chosen_date));
ti_swot = dsearchn(NORCAL.SWOT.mean_time, Radial.TIME(dsearchn(Radial.TIME,User_chosen_date)));

    SSHA = [ NORCAL.SWOT.ssha_karin_2{ti_swot}      +   NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
           [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
           [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
    % SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');
    % SSHA = smoothdata2(SSHA,'gaussian',3,'omitnan');
    % SSHA = smoothdata2(SSHA,'movmedian',7,'omitnan');
    % SSHA = smoothdata2_interpnadirgap(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, SSHA, 7);
    SSHA = bettergaussiansmooth2(SSHA,2000,GaussianKernelSTD_meters);

    % Smoothing is necessary to reduce high-k noise that can really throw
    % off the velocity when taking the first derivative via "gradient".


% close all
figure('Color','w')

M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(Radial.LON) max(Radial.LON)] + [-1 1]/1/2,...
             'latitudes', [min(Radial.LAT) max(Radial.LAT)] + [-1 1]/1/2);
set(gcf,'color','w'); hold on
VecScale = 0.001;
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                  NORCAL.SWOT.lat{ti_swot},...
                  SSHA);
% m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
%                   NORCAL.SWOT.lat{ti_swot},...
%                   NORCAL.SWOT.ocean_tide_fes{ti_swot});
m_contour(NORCAL.SWOT.lon{ti_swot},...
          NORCAL.SWOT.lat{ti_swot},...
          SSHA,[-0.15:0.01:0.15],'w')
m_quiver(Radial.LON(time_IND),...
         Radial.LAT(time_IND),...
         VecScale*Radial.VEL(time_IND).*[rHat(time_IND,1)*1 + rHat(time_IND,2)*0], ...
         VecScale*Radial.VEL(time_IND).*[rHat(time_IND,1)*0 + rHat(time_IND,2)*1],0,'k')
         % VecScale*DATA(time_IND,3), ...
         % VecScale*DATA(time_IND,4),0,'k')
% m_scatter(Radial.LON(time_IND),Radial.LAT(time_IND),10,...
%          sqrt(VecScale*DATA(time_IND,3).^2 + VecScale*DATA(time_IND,4).^2),'filled')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo)
CB = colorbar;
CB.Label.String = 'SWOT: \eta (m)'; CB.Label.FontSize = 20; clim([-1 1]*0.15)
set(gca,'FontSize',16)
title(['HFR time = ' datestr(Radial.TIME(dsearchn(Radial.TIME,User_chosen_date)))])

% % % % % % % % % % % 
%% Gradient*g/f
close all
gradSSHA = gradient_2D(SSHA,NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot})*[1/111000];
figure('Color','w')
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(Radial.LON) max(Radial.LON)] + [-1 1]/-1,...
             'latitudes', [min(Radial.LAT) max(Radial.LAT)] + [-1 1]/-1/2);
set(gcf,'color','w'); hold on
VecScale = 0.001;
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
                  NORCAL.SWOT.lat{ti_swot},...
                  sqrt(gradSSHA(:,:,1).^2 + gradSSHA(:,:,2).^2)*gg./fcor_degrees_cps(NORCAL.SWOT.lat{ti_swot}));
m_contour(NORCAL.SWOT.lon{ti_swot},...
          NORCAL.SWOT.lat{ti_swot},...
          SSHA,[-0.15:0.01:0.15],'w')
m_quiver(Radial.LON(time_IND),...
         Radial.LAT(time_IND),...
         VecScale*Radial.VEL(time_IND).*[rHat(time_IND,1)*1 + rHat(time_IND,2)*0], ...
         VecScale*Radial.VEL(time_IND).*[rHat(time_IND,1)*0 + rHat(time_IND,2)*1],0,'k')
% m_scatter(Radial.LON(time_IND),Radial.LAT(time_IND),10,...
%          sqrt(VecScale*DATA(time_IND,3).^2 + VecScale*DATA(time_IND,4).^2),'filled')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo)
CB = colorbar;
CB.Label.String = 'SWOT: $\frac{g}{f}|\nabla\eta|$ (m/s)'; CB.Label.FontSize = 20;
CB.Label.Interpreter = 'LaTeX';
% clim([0 1]*0.00001)
clim([0 1]*1)
set(gca,'FontSize',16)
title(['HFR time = ' datestr(Radial.TIME(dsearchn(Radial.TIME,User_chosen_date)))])




%% Overlay radial velocities over SWOT "radial" velocities relative to antenna

User_chosen_date = datenum('2023-05-07 02:00:00');
User_chosen_date = datenum('2023-05-09 02:00:00');
% User_chosen_date = datenum('2023-05-16 02:00:00');
% User_chosen_date = datenum('2023-05-29 02:00:00');
% User_chosen_date = datenum('2023-05-31 02:00:00');
time_IND = Radial.TIME == Radial.TIME(dsearchn(Radial.TIME,User_chosen_date));
ti_swot = dsearchn(NORCAL.SWOT.mean_time, Radial.TIME(dsearchn(Radial.TIME,User_chosen_date)));

    SSHA = [ NORCAL.SWOT.ssha_karin_2{     ti_swot} +   NORCAL.SWOT.height_cor_xover{     ti_swot}] .* ...
           [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
           [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
    SSHA = smoothdata2(SSHA,'gaussian',7,'omitnan');
    % SSHA = smoothdata2(SSHA,'gaussian',3,'omitnan');
    % SSHA = smoothdata2(SSHA,'movmedian',7,'omitnan');
    % SSHA = smoothdata2_interpnadirgap(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, SSHA, 7);
    % Smoothing is necessary to reduce high-k noise that can really throw
    % off the velocity when taking the first derivative via "gradient".


close all
figure('Color','w')

M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(Radial.LON) max(Radial.LON)] + [-1 1]/1/2,...
             'latitudes', [min(Radial.LAT) max(Radial.LAT)] + [-1 1]/1/2);
set(gcf,'color',[1 1 1]*0.5); hold on
VecScale = 0.001;
% m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},... |u_g + i*v_g|
%                   NORCAL.SWOT.lat{ti_swot},...
%                   sqrt([U_geostr(:,:,ti_swot)].^2 + [V_geostr(:,:,ti_swot)].^2))
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},... (u_g,v_g).rHat
                  NORCAL.SWOT.lat{ti_swot},...
                  squeeze(U_geostr(:,:,ti_swot)).*SWOT_radial_unit_x + ...
                  squeeze(V_geostr(:,:,ti_swot)).*SWOT_radial_unit_y);
% m_contour(NORCAL.SWOT.lon{ti_swot},...
%           NORCAL.SWOT.lat{ti_swot},...
%           SSHA,[-0.15:0.015:0.15],'k')
% m_quiver(Radial.LON(time_IND),...
%          Radial.LAT(time_IND),...
%          VecScale*DATA(time_IND,3), ...
%          VecScale*DATA(time_IND,4),0,'k')
m_scatter(Radial.LON(time_IND),Radial.LAT(time_IND),10,... magnitude and direction
          0.01*Radial.VEL(time_IND),'filled');
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
% colormap(turbo)
colormap(bwr)
CB = colorbar;
CB.Label.String = '{\bfu}\cdot{\bfr} (m/s)'; CB.Label.FontSize = 20; clim([-1 1]*1.0)
set(gca,'FontSize',16)
title(['HFR time = ' datestr(Radial.TIME(dsearchn(Radial.TIME,User_chosen_date)))])




%% Histogram of radial velocity magnitudes

% Antenna_ID_test = 'TRIN';
% Antenna_ID_test = 'GCVE';
Antenna_ID_test = 'SMOA';
% Antenna_ID_test = 'SHEL';
% Antenna_ID_test = 'BRAG';
% Antenna_ID_test = 'PAFS';

% % Display the standard deviation of radial velocities from different
% % antennas;
% for Antenna_ID_ii = {'TRIN','GCVE','SMOA','SHEL','BRAG','PAFS'}
%     DIR = dir([PATH 'Radial_data_allcalvaltimes_' Antenna_ID_ii{1} '.mat']);
%     Radial_VEL = load([DIR.folder '/' DIR.name],'VEL');
%     Radial_VEL = Radial_VEL.VEL;
%     disp(std(0.01*Radial_VEL(:),'omitnan'));
% end
% % % % OUTPUT:
% % 0.2505
% % 0.2327
% % 0.2886
% % 0.2561
% % 0.2733
% % 0.2832


DIR = dir([PATH 'Radial_data_allcalvaltimes_' Antenna_ID_test '.mat']);

Radial_VEL = load([DIR.folder '/' DIR.name],'VEL');
Radial_VEL = Radial_VEL.VEL; Radial_VEL = Radial_VEL(isfinite(Radial_VEL));

histogram_nan = @(IN) histogram(IN(isfinite(IN)),round( sqrt(sum(isfinite(IN(:))))/2 ),'Normalization','pdf');
histogram_nan_lims = @(IN,LIMS) histogram(IN(isfinite(IN)),LIMS,'Normalization','pdf');

close all
figure('Color','w')
AX1 = subplot(211);
% histogram_nan(squeeze(U_geostr(:,:,ti_swot)).*SWOT_radial_unit_x + ...
%               squeeze(V_geostr(:,:,ti_swot)).*SWOT_radial_unit_y); hold on
histogram_nan_lims(U_geostr.*repmat(SWOT_radial_unit_x,1,1,size(U_geostr,3)) + ...
                   V_geostr.*repmat(SWOT_radial_unit_y,1,1,size(V_geostr,3)),[-2:0.05:2]); hold on
histogram_nan_lims(0.01*Radial_VEL,[-2:0.05:2]);
% histogram_nan(0.01*Radial.VEL(time_IND));
% histogram_nan(0.01*Radial.VEL(time_IND & inpolygon(Radial.LON,Radial.LAT,real(SWOT_polygon),imag(SWOT_polygon)) ));
xlim([-1 1]*2)
% set(gca,'YScale','log')
LEG = legend(['SWOT: {\bfu_g\cdotr},   \sigma = ' ...
              num2str(std(Collimate(U_geostr.*repmat(SWOT_radial_unit_x,1,1,size(U_geostr,3)) + ...
                                    V_geostr.*repmat(SWOT_radial_unit_y,1,1,size(U_geostr,3)) ),'omitnan'))],...
             ['u_{radial HFR, ' Antenna_ID_test '},   \sigma = ' num2str(std(0.01*Radial_VEL,'omitnan'))]);
% LEG = legend(['SWOT: {\bfu_g\cdotr},   \sigma = ' num2str(std(Collimate(squeeze(U_geostr).*SWOT_radial_unit_x + squeeze(V_geostr).*SWOT_radial_unit_y),'omitnan'))],...
%              ['u_{radial HFR},   \sigma = ' num2str(std(0.01*Radial.VEL(time_IND),'omitnan'))]);
LEG.FontSize = 16;

AX2 = subplot(212);
% histogram_nan(NORCAL.HFR.u)
% histogram_nan(NORCAL.HFR.v)
histogram(nan,'HandleVisibility','off'); hold on
histogram(nan,'HandleVisibility','off')
histogram(NORCAL.HFR.u(isfinite(NORCAL.HFR.u)),[-1:0.05:1],'Normalization','pdf')
histogram(NORCAL.HFR.v(isfinite(NORCAL.HFR.v)),[-1:0.05:1],'Normalization','pdf')
xlim([-1 1]*2)
% set(gca,'YScale','log')
LEG = legend(['u_{gridded HFR},   \sigma = ' num2str(std(NORCAL.HFR.u(isfinite(NORCAL.HFR.u)),'omitnan'))], ...
             ['v_{gridded HFR},   \sigma = ' num2str(std(NORCAL.HFR.v(isfinite(NORCAL.HFR.v)),'omitnan'))]);
LEG.FontSize = 16;

linkaxes([AX1 AX2],'xy')

% std(0.01*Radial.VEL(time_IND),'omitnan')
% std(Collimate(U_geostr.*SWOT_radial_unit_x + ...
%               V_geostr.*SWOT_radial_unit_y),'omitnan')
% ans =
%     0.2502
% ans =
%     0.2280

figure('Color','w')
Velg_radial_isfinite = abs(U_geostr.*SWOT_radial_unit_x + V_geostr.*SWOT_radial_unit_y);
    Velg_radial_isfinite = Velg_radial_isfinite(isfinite(Velg_radial_isfinite));
histogram(Velg_radial_isfinite,[0:0.05:2],'Normalization','pdf','DisplayStyle','stairs','LineWidth',1); hold on
histogram(abs(0.01*Radial_VEL),[0:0.05:2],'Normalization','pdf','DisplayStyle','stairs','LineWidth',1);
histogram(abs(NORCAL.HFR.u(isfinite(NORCAL.HFR.u) & isfinite(NORCAL.HFR.v)) + ...
           1i*NORCAL.HFR.v(isfinite(NORCAL.HFR.u) & isfinite(NORCAL.HFR.v)) ),...
           [0:0.05:2],'Normalization','pdf','DisplayStyle','stairs','LineWidth',1)
legend('|(u_g + iv_g)\cdotr|','|u_r|','|u_{HFR} + iv_{HFR}|')

%% Click a point and plot HFR velocity and parallel SWOT geostrophic velocities side-by-side

% % SSHA snapshot:
% close all
% figure('Color','w')
% M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
%              'longitudes',[min(Radial.LON) max(Radial.LON)] + [-1 1]/1,...
%              'latitudes', [min(Radial.LAT) max(Radial.LAT)] + [-1 1]/1);
% set(gcf,'color','w'); hold on
% VecScale = 0.001;
% m_pcolor_centered(NORCAL.SWOT.lon{ti_swot},...
%                   NORCAL.SWOT.lat{ti_swot},...
%                   SSHA);
% m_plot(real(Radial.Origin_antenna), imag(Radial.Origin_antenna), 'ko')
% COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
% m_grid('box','fancy', 'backgroundcolor','none');
% colormap(turbo)
% CB = colorbar;
% CB.Label.String = 'SWOT: \eta (m)'; CB.Label.FontSize = 20; clim([-1 1]*0.15)
% set(gca,'FontSize',16)

% Pre-calculated correlation map
P_threshold = 0.05;
close all
figure('Color','w')
AX1 = subplot(121);
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
PC = ...
m_pcolor(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  ... SSHA);
                  Radial_vel_corrR);
PC.ZData = PC.CData;
m_scatter(real(unique(Radial.LON(isfinite(Radial.VEL)) + 1i*Radial.LAT(isfinite(Radial.VEL)))),...
          imag(unique(Radial.LON(isfinite(Radial.VEL)) + 1i*Radial.LAT(isfinite(Radial.VEL)))),1,'k','filled')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo)
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
AX2 = subplot(122);
CB.Label.String = ['U_{radial}: corr. coef. (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
PC = ...
m_pcolor(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  ...SSHA);
                  Radial_vel_lowpass_corrR);
PC.ZData = PC.CData;
m_scatter(real(unique(Radial.LON + 1i*Radial.LAT)),...
          imag(unique(Radial.LON + 1i*Radial.LAT)),1,'k','filled')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo)
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
CB.Label.String = ['U_{radial, lowpass}: corr. coef. (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
linkaxes([AX1 AX2],'xy')

[x_click,y_click] = m_ginput(1);
[yi_click,xi_click] = find(abs([NORCAL.SWOT.lon{ti_swot}    + 1i*NORCAL.SWOT.lat{ti_swot}]    - [x_click + 1i*y_click]) == ...
                       min(abs([NORCAL.SWOT.lon{ti_swot}(:) + 1i*NORCAL.SWOT.lat{ti_swot}(:)] - [x_click + 1i*y_click])));


m_plot(NORCAL.SWOT.lon{ti_swot}(yi_click,xi_click),...
       NORCAL.SWOT.lat{ti_swot}(yi_click,xi_click),'ow','MarkerSize',10,'LineWidth',2)

subplot(121)
m_plot(NORCAL.SWOT.lon{ti_swot}(yi_click,xi_click),...
       NORCAL.SWOT.lat{ti_swot}(yi_click,xi_click),'ow','MarkerSize',10,'LineWidth',2)

% IND = [Radial.LON - x_click == min(Radial.LON - x_click)] & [Radial.LAT - y_click == min(Radial.LAT - y_click)];
% %%
KM_search = 0.4;
KM_search = 2.5;
KM_search = 3.0;
% KM_search = 6.0;
% IND = [abs(Radial.LON - x_click) < KM_search/111] & [abs(Radial.LAT - y_click) < KM_search/111];
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
%%

% close all

Vel_g_radialdir = squeeze(U_geostr(yi_click,xi_click,:)).*SWOT_radial_unit_x(yi_click,xi_click) + ...
                  squeeze(V_geostr(yi_click,xi_click,:)).*SWOT_radial_unit_y(yi_click,xi_click);
Ur_HFR_at_point_atSWOTtimes = interp1(T_HFR_at_point,Ur_HFR_at_point,NORCAL.SWOT.mean_time,'linear');

WINDOW = hanning(24); WINDOW = WINDOW/sum(WINDOW); % low pass
ConvWithWindow = @(IN) conv(squeeze(IN),WINDOW,'same');
Ur_lowpass = ConvWithWindow(Ur_HFR_at_point);
Ur_HFR_at_point_atSWOTtimes_lowpass = interp1(T_HFR_at_point,Ur_lowpass,NORCAL.SWOT.mean_time,'linear');


figure('Color','w')
AX1 = gca;
plot(NORCAL.SWOT.mean_time, Vel_g_radialdir , 'k.-','MarkerSize',18,'LineWidth',1.5); hold on
% plot(Radial.TIME(IND), Radial.VEL(IND)/100, 'c.-') % {time,velocity mag (m/s)}
% % % 
plot(T_HFR_at_point,Ur_HFR_at_point, 'r.-') % {time,velocity mag (m/s)}
errorbar(T_HFR_at_point,Ur_HFR_at_point,UrSTD_HFR_at_point,'r','HandleVisibility','off')
% % % 
plot(T_HFR_at_point,Ur_lowpass, 'g.-')
legend('SWOT',['HFR at <' num2str(KM_search) ' km'],['HFR at <' num2str(KM_search) ' km low passed'])
datetick('x')

[R_vel_full,P_vel_full] = corrcoef(Vel_g_radialdir(isfinite(Vel_g_radialdir) & isfinite(Ur_HFR_at_point_atSWOTtimes)),...
          Ur_HFR_at_point_atSWOTtimes(isfinite(Vel_g_radialdir) & isfinite(Ur_HFR_at_point_atSWOTtimes)));
[R_vel_lowp,P_vel_lowp] = corrcoef(Vel_g_radialdir(        isfinite(Vel_g_radialdir) & isfinite(Ur_HFR_at_point_atSWOTtimes_lowpass)),...
          Ur_HFR_at_point_atSWOTtimes_lowpass(isfinite(Vel_g_radialdir) & isfinite(Ur_HFR_at_point_atSWOTtimes_lowpass)));
title({['CorrCoef(full signal)    = ' num2str(R_vel_full(1,2)) '  (P = ' num2str(P_vel_full(1,2)) ')'] ; ...
       ['CorrCoef(lowpass signal) = ' num2str(R_vel_lowp(1,2)) '  (P = ' num2str(P_vel_lowp(1,2)) ')'] ; ' '})


figure('Color','w')
AX2 = gca;
plot(NORCAL.SWOT.mean_time,Vel_g_radialdir , 'k.-','MarkerSize',12,'LineWidth',1); hold on
% % % 
plot(NORCAL.SWOT.mean_time,Ur_HFR_at_point_atSWOTtimes, 'r.-','MarkerSize',12,'LineWidth',1) % {time,velocity mag (m/s)}
% % % 
plot(NORCAL.SWOT.mean_time,Ur_HFR_at_point_atSWOTtimes_lowpass, 'g.-','MarkerSize',12,'LineWidth',1)
legend('SWOT',['HFR at <' num2str(KM_search) ' km'],['HFR at <' num2str(KM_search) ' km low passed'])
datetick('x')
set(gca,'FontSize',14)
ylabel('m/s')
title({['CorrCoef(full signal)    = ' num2str(R_vel_full(1,2)) '  (P = ' num2str(P_vel_full(1,2)) ')'] ; ...
       ['CorrCoef(lowpass signal) = ' num2str(R_vel_lowp(1,2)) '  (P = ' num2str(P_vel_lowp(1,2)) ')'] ; ' '})

linkaxes([AX1 AX2],'xy')

% figure('Color','w')
% plot(NORCAL.SWOT.mean_time, ...
%      squeeze(U_geostr(yi_click,xi_click,:)).*SWOT_radial_unit_x(yi_click,xi_click) + ...
%      squeeze(V_geostr(yi_click,xi_click,:)).*SWOT_radial_unit_y(yi_click,xi_click) - ...
%      mean(squeeze(U_geostr(yi_click,xi_click,:)).*SWOT_radial_unit_x(yi_click,xi_click) + ...
%           squeeze(V_geostr(yi_click,xi_click,:)).*SWOT_radial_unit_y(yi_click,xi_click),'omitnan'), 'k.-'); hold on
% plot(    T_HFR_at_point,Ur_HFR_at_point - mean(Ur_HFR_at_point,'omitnan'), 'r.-') % {time,velocity mag (m/s)}
% errorbar(T_HFR_at_point,Ur_HFR_at_point - mean(Ur_HFR_at_point,'omitnan'),UrSTD_HFR_at_point,'r')
% title('both de-meaned')
% legend('SWOT',['HFR at <' num2str(KM_search) ' km'])
% datetick('x')

% Calculate power spectrum of the HFR radial velocity:
figure
warning('requires "nunanspectrum" function')
nunanspectrum(Ur_HFR_at_point(2:end),                    T_HFR_at_point(2:end),       'day','Segments',9,'Plot',true,'Plot_option','b.-');
nunanspectrum(Ur_HFR_at_point_atSWOTtimes(2:end),        NORCAL.SWOT.mean_time(2:end),'day','Segments',4,'Plot',true,'Plot_option','b.--');
nunanspectrum(Ur_lowpass(2:end),                         T_HFR_at_point(2:end),       'day','Segments',9,'Plot',true,'Plot_option','r.-');
nunanspectrum(Ur_HFR_at_point_atSWOTtimes_lowpass(2:end),NORCAL.SWOT.mean_time(2:end),'day','Segments',4,'Plot',true,'Plot_option','r.--');
nunanspectrum(Vel_g_radialdir(2:end),                    NORCAL.SWOT.mean_time(2:end),'day','Segments',4,'Plot',true,'Plot_option','k.--');
title('Power Spectrum (using NUFFT)')
legend(['HFR at <' num2str(KM_search) ' km'],...
       ['HFR at <' num2str(KM_search) ' km (SWOT times)'],...
       ['HFR at <' num2str(KM_search) ' km low passed'],...
       ['HFR at <' num2str(KM_search) ' km low passed (SWOT times)'],...
       ['SWOT u_g \cdot r (SWOT times, obviously)'])

%% Tidal fit
%% Fit at a single location to test

% Radial.Unique_LONLAT = Radial.LON + Radial.LAT*1i;
% Radial.Unique_LONLAT = unique(Radial.Unique_LONLAT);

% I could use my pre-existing red_tide package, and may eventually, but for
% now I will craft the least-squares fit from the ground up:

% Frequencies:
df = 1/[24*[Radial.TIME(end) - Radial.TIME(1)]]; % hr^-1
NIO_Frequencies = ([2*[7.2921*10^-5]*sind(min(Radial.LAT))*60*60]/[2*pi]):df:...
                  ([2*[7.2921*10^-5]*sind(max(Radial.LAT))*60*60]/[2*pi]); % hr^-1
NIO_Frequencies = 2*pi*NIO_Frequencies; % Redundant step to convert to rad/hr
Tidal_frequencies = 1./[25.81933871   24   12.65834751   12.4206012   12]; % in hours
% Tidal_frequencies = 1./[25.81933871   24   ]; % in hours
%     Tidal_frequencies = [Tidal_frequencies, [1/13]:df:[1/11.5]];
Omega = 2*pi.*Tidal_frequencies;
Omega = unique([Omega , NIO_Frequencies]); % rad/hr

% Test with the click-selected point from earlier:
TT = [T_HFR_at_point - T_HFR_at_point(1)]*24; % hours after t(1)
TT = TT(isfinite(TT));
Ur = Ur_HFR_at_point(isfinite(TT)) - 0*Ur_lowpass;
% Ur = sin(TT*2*pi/12.4206012 + 2*pi*rand);
HH = [ones(size(TT)) sin(TT*Omega) cos(TT*Omega)];
Ur_detided = Ur - HH*[HH\Ur] + 0*Ur_lowpass;
Ur_tidal   = HH*[HH\Ur] + 0*Ur_lowpass;

close all

figure
plot(TT, Ur, '.-k'); hold on
plot(TT, Ur_detided, '.-r')
% plot(TT, Ur_tidal, '.-r')
legend('HFR','HFR detided')

% Calculate power spectrum of the HFR radial velocity, and compare that to
% the p.s. of the detided data:
figure
warning('requires "nunanspectrum" function')
nunanspectrum(Ur_HFR_at_point(2:end),T_HFR_at_point(2:end),'day','Segments',9,'Plot',true);
nunanspectrum(Ur_detided(2:end),TT(2:end)/24,'day','Segments',9,'Plot',true,'Plot_option','r.-');
loglog(Collimate([Omega;Omega;nan(size(Omega))])*24/[2*pi],...
       Collimate([ones(size(Omega))*10^-4; ones(size(Omega))*10^-1; nan(size(Omega))]), 'k')
title('Power Spectrum (using NUFFT)')
legend('HFR','HFR detided')

%% Fit at each unique location

% Unnecessary, as tidal fitting is not as advantageous to this analysis as
% originally anticipated.


%% Mean radial velocity for HFR and SWOT plotted

close all
figure('Color','w')

AX1 = subplot(121);
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(Radial.LON) max(Radial.LON)] + [-1 1]/1,...
             'latitudes', [min(Radial.LAT) max(Radial.LAT)] + [-1 1]/1);
set(gcf,'color',[1 1 1]*0.5); hold on
VecScale = 0.001;
m_pcolor_centered(NORCAL.SWOT.lon{1},...
                  NORCAL.SWOT.lat{1},...
                  mean(U_geostr.*SWOT_radial_unit_x + ...
                       V_geostr.*SWOT_radial_unit_y , 3, 'omitnan'));
m_plot(real(Radial.Origin_antenna), imag(Radial.Origin_antenna), 'ko')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
colormap(turbo)
CB = colorbar;
CB.Label.String = 'SWOT: |{\bf u}_g\cdot{\bf r}| (m)'; CB.Label.FontSize = 20; clim([-1 1]*1)
set(gca,'FontSize',16)

AX2 = subplot(122);
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[min(Radial.LON) max(Radial.LON)] + [-1 1]/1,...
             'latitudes', [min(Radial.LAT) max(Radial.LAT)] + [-1 1]/1);
% set(gcf,'color','w');
hold on
VecScale = 0.001;
M_S = m_scatter(real(Radial.Unique_LONLAT), imag(Radial.Unique_LONLAT),10,... magnitude and direction
                0.01*Radial.mean_VEL,'filled');
m_plot(real(Radial.Origin_antenna), imag(Radial.Origin_antenna), 'ko')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
colormap(bwr)
CB = colorbar;
CB.Label.String = 'HFR: |radial velocity| (m)'; CB.Label.FontSize = 20; clim([-1 1]*1)
set(gca,'FontSize',16)
title('Change dot size by redefining M\_S.SizeData')

linkaxes([AX1 AX2],'xy')

%%
%%
%% The following needs to be run on the server because it will take a very long time to run

% Perform on server using the separate scripts
% process_workspace_temp_****.m after saving full workspace and manually
% transferring it to the server:

error('Already saved.')

save(['./data/workspace_temp_' Antenna_ID '.mat'])

error('Force stop by user, run the rest of this script on a server.')

%%
%% Calculate the difference in velocity at each radar point over time,
%  and the correlation coefficient between the two.
%% First figure out how far the closest radial is from each SWOT point,
% which will be used in the next step.
% This should take ~30 minutes to complete:
if exist('Dist_to_nearest_radial','var')
    error('You already calculated these variables, which take a long time to calculate.')
else
end
Dist_to_nearest_radial = nan(size(NORCAL.SWOT.lat{1}));
unique_LonLAT_HFR = unique(Radial.LON+1i*Radial.LAT);
for ii = 1:size(NORCAL.SWOT.lat{1},1)
    for jj = 1:size(NORCAL.SWOT.lat{1},2)
        Dist_to_nearest_radial(ii,jj) = ...
            min(m_idist(NORCAL.SWOT.lon{1}(ii,jj),NORCAL.SWOT.lat{1}(ii,jj), ...
                        real(unique_LonLAT_HFR) , imag(unique_LonLAT_HFR) ));
    end
    disp(100*ii/size(NORCAL.SWOT.lat{1},1))
end

mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
CurrentDir = dir(mfilePath);
cd(CurrentDir.folder)
cd ../data
save(['./HFR_radials/Dist_to_nearest_radial_' Antenna_ID '.mat'],'Dist_to_nearest_radial','Antenna_ID')
disp('done')
%% Now compare the SWOT and HFR radial velocities head-to-head
% Make this on the SWOT grid, with NaN where there are not enough nearby
% radar points to give useful information. I have to do it as a loop; maybe
% there is a more efficient way, but I only need an answer.
if exist('Radial_vel_rmsd','var')
    error('You already calculated these variables, which take a long time to calculate.')
else
end
KM_search = 3.0; % search radius for averaging radial velocities to compare to a SWOT grid point
SWOT_mean_time = NORCAL.SWOT.mean_time; % For naming simplicity
Radial_vel_rmsd = nan(size(NORCAL.SWOT.lat{1}));
Radial_vel_corrR = nan(size(NORCAL.SWOT.lat{1}));
Radial_vel_corrP = nan(size(NORCAL.SWOT.lat{1}));
Radial_vel_lowpass_rmsd = nan(size(NORCAL.SWOT.lat{1}));
Radial_vel_lowpass_corrR = nan(size(NORCAL.SWOT.lat{1}));
Radial_vel_lowpass_corrP = nan(size(NORCAL.SWOT.lat{1}));

WINDOW = hanning(24); WINDOW = WINDOW/sum(WINDOW); % low pass
ConvWithWindow = @(IN) conv(squeeze(IN),WINDOW,'same');

for ii = 1:size(NORCAL.SWOT.lat{1},1)
    for jj = 1:size(NORCAL.SWOT.lat{1},2)
        % ii = 124; jj = 50; %$
        if sum(isfinite(U_geostr(ii,jj,:)),3)/size(U_geostr,3) >= 0.9 & Dist_to_nearest_radial(ii,jj) < 5*KM_search*1000
            
            SWOT_radial_vel_ts = squeeze(U_geostr(ii,jj,:)).*SWOT_radial_unit_x(ii,jj) + ...
                                 squeeze(V_geostr(ii,jj,:)).*SWOT_radial_unit_y(ii,jj);

            [Distance_from_SWOTcell,~,~] = m_idist(NORCAL.SWOT.lon{1}(ii,jj),NORCAL.SWOT.lat{1}(ii,jj), ...
                                                   Radial.LON,Radial.LAT);
            IND = [Distance_from_SWOTcell/1000 < KM_search];

            T_HFR_at_point = unique(Radial.TIME(IND));
            Ur_HFR_at_point = nan(length(T_HFR_at_point),1);
            % UrSTD_HFR_at_point = nan(length(T_HFR_at_point),1);
            tic
            for kk = 1:length(T_HFR_at_point)
                Ur_HFR_at_point(kk)    = mean(Radial.VEL(IND & Radial.TIME==T_HFR_at_point(kk))/100,"omitmissing");
                % UrSTD_HFR_at_point(ii) = std( Radial.VEL(IND & Radial.TIME==T_HFR_at_point(ii))/100,"omitmissing");
            end

            
            Ur_HFR_at_point_atSWOTtimes = interp1(T_HFR_at_point,Ur_HFR_at_point,NORCAL.SWOT.mean_time,'linear');;
            Radial_vel_rmsd(ii,jj) = rms(Ur_HFR_at_point_atSWOTtimes - SWOT_radial_vel_ts,'omitnan');

            [Radial_vel_corr_R,Radial_vel_corr_P] = ...
                corrcoef(Ur_HFR_at_point_atSWOTtimes(isfinite(Ur_HFR_at_point_atSWOTtimes) & isfinite(SWOT_radial_vel_ts)), ...
                SWOT_radial_vel_ts(isfinite(Ur_HFR_at_point_atSWOTtimes) & isfinite(SWOT_radial_vel_ts)));
            if isscalar(Radial_vel_corr_R)
                Radial_vel_corrR(ii,jj) = NaN;
                Radial_vel_corrP(ii,jj) = NaN;
            else
                Radial_vel_corrR(ii,jj) = Radial_vel_corr_R(1,2);
                Radial_vel_corrP(ii,jj) = Radial_vel_corr_P(1,2);
            end



            Ur_lowpass = ConvWithWindow(Ur_HFR_at_point);
            Ur_HFR_at_point_atSWOTtimes_lowpass = interp1(T_HFR_at_point,Ur_lowpass,NORCAL.SWOT.mean_time,'linear');;
            Radial_vel_lowpass_rmsd(ii,jj) = rms(Ur_HFR_at_point_atSWOTtimes_lowpass - SWOT_radial_vel_ts,'omitnan');

            [Radial_vel_lowpass_corr_R,Radial_vel_lowpass_corr_P] = ...
                corrcoef(Ur_HFR_at_point_atSWOTtimes_lowpass(isfinite(Ur_HFR_at_point_atSWOTtimes_lowpass) & isfinite(SWOT_radial_vel_ts)), ...
                         SWOT_radial_vel_ts(                 isfinite(Ur_HFR_at_point_atSWOTtimes_lowpass) & isfinite(SWOT_radial_vel_ts)));
            if isscalar(Radial_vel_lowpass_corr_R)
                Radial_vel_lowpass_corrR(ii,jj) = NaN;
                Radial_vel_lowpass_corrP(ii,jj) = NaN;
            else
                Radial_vel_lowpass_corrR(ii,jj) = Radial_vel_lowpass_corr_R(1,2);
                Radial_vel_lowpass_corrP(ii,jj) = Radial_vel_lowpass_corr_P(1,2);
            end
            % error % for testing
        else % Radial_vel_diff(ii,jj) = NaN
        end
    end
    disp(num2str(100*ii/size(NORCAL.SWOT.lat{1},1)))
end

mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
CurrentDir = dir(mfilePath);
cd(CurrentDir.folder)
cd ../data
save(['./Radial_vel_comparison_' Antenna_ID '.mat'], ...
      'Radial_vel_lowpass_rmsd','Radial_vel_lowpass_corrR','Radial_vel_lowpass_corrP',...
      'Radial_vel_rmsd','Radial_vel_corrR','Radial_vel_corrP','Antenna_ID')
disp('done')

%% Map correlation

close all
figure('Color','w')
tiledlayout(1,3,'TileSpacing','tight')
% tiledlayout(1,4,'TileSpacing','tight')

P_threshold = 0.05;

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  Radial_vel_rmsd);
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([0 1]*0.5)
CB.Label.String = 'RMSD (m/s)'; CB.Label.FontSize = 16;

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  Radial_vel_corrR .* [Radial_vel_corrP<P_threshold]./[Radial_vel_corrP<P_threshold]);
                  ...Radial_vel_corrR);
% m_scatter(real(unique(Radial.LON + 1i*Radial.LAT)),...
%           imag(unique(Radial.LON + 1i*Radial.LAT)),1,'k','filled')
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
CB.Label.String = ['U_{full radial}: corr. coef. (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

AX3 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  Radial_vel_lowpass_corrR .* [Radial_vel_lowpass_corrP<P_threshold]./[Radial_vel_lowpass_corrP<P_threshold]);
                  ...Radial_vel_lowpass_corrR);
% m_scatter(real(unique(Radial.LON + 1i*Radial.LAT)),...
%           imag(unique(Radial.LON + 1i*Radial.LAT)),1,'k','filled')
COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
CB.Label.String = ['U_{low pass radial}: corr. coef. (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

% AX4 = nexttile;
% M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
%              'longitudes',[-126 -122.] + [-.1 .1],...
%              'latitudes',[37 42] + [-.1 .1]);
% set(gcf,'color','w'); hold on
% m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
%                   [Radial_vel_lowpass_corrR .* [Radial_vel_lowpass_corrP<P_threshold]./[Radial_vel_lowpass_corrP<P_threshold]] - ...
%                   [Radial_vel_corrR .* [        Radial_vel_corrP<        P_threshold]./[Radial_vel_corrP<        P_threshold]]);
% COAST = m_gshhs_l('patch',0.5*[1 1 1],'FaceAlpha',0.5);
% m_grid('box','fancy', 'backgroundcolor','none');
% set(gca,'FontSize',16)
% CB = colorbar; clim([-1 1]*0.5)
% CB.Label.String = ['CorrCoef U_{low pass radial} - CorrCoef U_{full radial} (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;

colormap(turbo)
% colormap(bwr)

linkaxes([AX1 AX2 AX3],'xy')
% linkaxes([AX1 AX2 AX3 AX4],'xy')





figure
histogram_nan = @(IN,BINS) histogram(IN(isfinite(IN)),BINS,'Normalization','pdf');
tiledlayout(1,3,'TileSpacing','tight')
nexttile
histogram_nan(Radial_vel_rmsd,0:0.02:1); hold on
histogram_nan(Radial_vel_lowpass_rmsd,0:0.02:1)
xlabel('RMSD')
ylabel('PDF')
legend('SWOT vs. HFR (total)','SWOT vs. HFR (filtered)')
nexttile
histogram_nan(Radial_vel_corrR .* [Radial_vel_corrP<P_threshold]./[Radial_vel_corrP<P_threshold],...
              -1:0.05:1); hold on
histogram_nan(Radial_vel_lowpass_corrR .* [Radial_vel_lowpass_corrP<P_threshold]./[Radial_vel_lowpass_corrP<P_threshold],...
              -1:0.05:1)
xlabel('CorrCoef')
ylabel('PDF')
legend('SWOT vs. HFR (total)','SWOT vs. HFR (filtered)')
nexttile
histogram_nan([Radial_vel_lowpass_corrR .* [Radial_vel_lowpass_corrP<P_threshold]./[Radial_vel_lowpass_corrP<P_threshold]] - ...
              [Radial_vel_corrR .* [Radial_vel_corrP<P_threshold]./[Radial_vel_corrP<P_threshold]], ...
              -1:0.01:1)
xlabel('\DeltaCorrCoef (lowpass - full)')
ylabel('PDF')
legend('SWOT vs. HFR (total)','SWOT vs. HFR (filtered)')


%% Map all antennas' correlations together

% Custom Color Map inspired by cmocean('oxy')
% CMAP = nan(100,3);
% CMAP( 1:25, 1) = linspace(0.8,0.1,25)';
% CMAP( 1:25, 2) = linspace(0.8,0.1,25)';
% CMAP( 1:25, 3) = ones(            25,1);
% CMAP(26:75, 1) = linspace(0.8,0.1,50)';
% CMAP(26:75, 2) = linspace(0.8,0.1,50)';
% CMAP(26:75, 3) = linspace(0.8,0.1,50)';
% CMAP(76:100,1) = ones(            25,1);
% CMAP(76:100,2) = linspace(0.9,0.0,25)';
% CMAP(76:100,3) = linspace(0.5,0.0,25)';

% CMAP = nan(100,3);
% CMAP( 1:50, 1) = linspace(0.1,0.8,50)';
% CMAP( 1:50, 2) = linspace(0.1,0.8,50)';
% CMAP( 1:50, 3) = ones(            50,1);
% CMAP(51:75, 1) = linspace(1,0,    25)';
% CMAP(51:75, 2) = linspace(1,0,    25)';
% CMAP(51:75, 3) = linspace(1,0,    25)';
% CMAP(76:100,1) = ones(            25,1);
% CMAP(76:100,2) = linspace(0.9,0.0,25)';
% CMAP(76:100,3) = linspace(0.5,0.0,25)';

CMAP = bwr(100);
% CMAP = CMAP.*[abs([[.02:.02:2]-1]') ones(100,1) ones(100,1)];
% CMAP = CMAP.*[ones(100,1) ones(100,1) abs([[.01:.01:1]]')];

% %%
if exist('antenna_IDs','var')
    % error
else

    % antenna_IDs = {'SHEL','BRAG','PAFS'};
    antenna_IDs = {'TRIN','SHEL','BRAG','PAFS','GCVE'};

    % load(['../Dist_to_nearest_radial_' Antenna_ID '.mat'])
    % load(['../Radial_vel_comparison_' Antenna_ID '.mat'])

    for ii = 1:length(antenna_IDs)
        antenna_ID = antenna_IDs{ii};
        eval([antenna_ID ' = load(''./data/HFR_radials/Radial_vel_comparison_' antenna_ID '.mat'');'])
        eval([antenna_ID '.Origin_antenna = load(''./data/HFR_radials/Radial_data_allcalvaltimes_' antenna_ID '.mat'',''Origin_antenna'');'])
        eval([antenna_ID '.Origin_antenna = ' antenna_ID '.Origin_antenna.Origin_antenna;'])
    end

    % % % Full HFR
    maxCorrCoef = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(antenna_IDs));
    maxCorrCoef(:,:,1) = TRIN.Radial_vel_corrR .* [TRIN.Radial_vel_corrP<P_threshold]./[TRIN.Radial_vel_corrP<P_threshold];
    maxCorrCoef(:,:,2) = GCVE.Radial_vel_corrR .* [GCVE.Radial_vel_corrP<P_threshold]./[GCVE.Radial_vel_corrP<P_threshold];
    maxCorrCoef(:,:,3) = SHEL.Radial_vel_corrR .* [SHEL.Radial_vel_corrP<P_threshold]./[SHEL.Radial_vel_corrP<P_threshold];
    maxCorrCoef(:,:,4) = BRAG.Radial_vel_corrR .* [BRAG.Radial_vel_corrP<P_threshold]./[BRAG.Radial_vel_corrP<P_threshold];
    maxCorrCoef(:,:,5) = PAFS.Radial_vel_corrR .* [PAFS.Radial_vel_corrP<P_threshold]./[PAFS.Radial_vel_corrP<P_threshold];
    maxCorrCoef = max(maxCorrCoef,[],3,'omitnan');

    % % % Filtered HFR
    maxCorrCoef_lowpass = nan(size(NORCAL.SWOT.lon{1},1),size(NORCAL.SWOT.lon{1},2),length(antenna_IDs));
    maxCorrCoef_lowpass(:,:,1) = TRIN.Radial_vel_lowpass_corrR .* [TRIN.Radial_vel_lowpass_corrP<P_threshold]./[TRIN.Radial_vel_lowpass_corrP<P_threshold];
    maxCorrCoef_lowpass(:,:,2) = GCVE.Radial_vel_lowpass_corrR .* [GCVE.Radial_vel_lowpass_corrP<P_threshold]./[GCVE.Radial_vel_lowpass_corrP<P_threshold];
    maxCorrCoef_lowpass(:,:,3) = SHEL.Radial_vel_lowpass_corrR .* [SHEL.Radial_vel_lowpass_corrP<P_threshold]./[SHEL.Radial_vel_lowpass_corrP<P_threshold];
    maxCorrCoef_lowpass(:,:,4) = BRAG.Radial_vel_lowpass_corrR .* [BRAG.Radial_vel_lowpass_corrP<P_threshold]./[BRAG.Radial_vel_lowpass_corrP<P_threshold];
    maxCorrCoef_lowpass(:,:,5) = PAFS.Radial_vel_lowpass_corrR .* [PAFS.Radial_vel_lowpass_corrP<P_threshold]./[PAFS.Radial_vel_lowpass_corrP<P_threshold];
    maxCorrCoef_lowpass = max(maxCorrCoef_lowpass,[],3,'omitnan');
end

close all
figure('Color','w')
tiledlayout(1,2,'TileSpacing','tight')
P_threshold = 0.05;

AX1 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  maxCorrCoef);
m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
          maxCorrCoef,[0.5 0.5],'k','LineWidth',1);
m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
          maxCorrCoef,[0.75 0.75],'k');
% m_scatter(real(unique(Radial.LON + 1i*Radial.LAT)),...
%           imag(unique(Radial.LON + 1i*Radial.LAT)),1,'k','filled')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
for ii = 1:length(antenna_IDs) % unfiltered
    antenna_ID = antenna_IDs{ii};
    m_plot(real(eval([antenna_ID '.Origin_antenna'])), ...
           imag(eval([antenna_ID '.Origin_antenna'])), 'o', ...
           'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
    m_text(real(eval([antenna_ID '.Origin_antenna'])) + 0.1, ...
           imag(eval([antenna_ID '.Origin_antenna'])), antenna_ID)
end
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
% CB.Label.String = ['$u_\mathrm{radial}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT,HFR),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
CB.Label.Interpreter = 'latex';
colormap('turbo')

AX2 = nexttile;
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
             'longitudes',[-126 -122.] + [-.1 .1],...
             'latitudes',[37 42] + [-.1 .1]);
set(gcf,'color','w'); hold on
m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                  maxCorrCoef_lowpass);
m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
          maxCorrCoef_lowpass,[0.5 0.5],'k','LineWidth',1);
m_contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
          maxCorrCoef_lowpass,[0.75 0.75],'k');
% m_scatter(real(unique(Radial.LON + 1i*Radial.LAT)),...
%           imag(unique(Radial.LON + 1i*Radial.LAT)),1,'k','filled')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
for ii = 1:length(antenna_IDs) % unfiltered
    antenna_ID = antenna_IDs{ii};
    m_plot(real(eval([antenna_ID '.Origin_antenna'])), ...
           imag(eval([antenna_ID '.Origin_antenna'])), 'o', ...
           'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k')
end
set(gca,'FontSize',16)
CB = colorbar; clim([-1 1])
% CB.Label.String = ['$u_\mathrm{radial, lowpass}$ corr. coef. (where P$<$' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
CB.Label.String = ['CorrCoef$_\mathrm{radial}$(SWOT,HFR$_\mathrm{LP}$),  P$<$' num2str(P_threshold) '']; CB.Label.FontSize = 16;
CB.Label.Interpreter = 'latex';
% colormap('turbo')
colormap(CMAP)

set(gcf,'Position',[-919   238   882   540])
linkaxes([AX1 AX2],'xy')

% %%
% figure(1)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/radial_SWOT_corrcoef_2panel.pdf',...
% 'BackgroundColor','none','ContentType','vector')

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
                  maxCorrCoef_lowpass - maxCorrCoef);
% m_scatter(real(unique(Radial.LON + 1i*Radial.LAT)),...
%           imag(unique(Radial.LON + 1i*Radial.LAT)),1,'k','filled')
COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
m_grid('box','fancy', 'backgroundcolor','none');
for ii = 1:length(antenna_IDs) % unfiltered
    antenna_ID = antenna_IDs{ii};
    m_plot(real(eval([antenna_ID '.Origin_antenna'])), ...
           imag(eval([antenna_ID '.Origin_antenna'])), 'o', ...
           'MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    m_text(real(eval([antenna_ID '.Origin_antenna'])) + 0.1, ...
           imag(eval([antenna_ID '.Origin_antenna'])), antenna_ID)
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
histogram(maxCorrCoef_lowpass - maxCorrCoef,[-1:0.01:1],'Normalization','pdf'); hold on
text(-1,1,['Mean = ' num2str(mean(maxCorrCoef_lowpass(:) - maxCorrCoef(:),'omitnan'))])
text(-1,0.9,['Median = ' num2str(median(maxCorrCoef_lowpass(:) - maxCorrCoef(:),'omitnan'))])
text(-1,0.8,['STDev = ' num2str(std(maxCorrCoef_lowpass(:) - maxCorrCoef(:),'omitnan'))])


% %%
% figure(1)
% exportgraphics(gcf,...
% '/Users/kachelein/Documents/JPL/papers/my_work/SWOT_HFR_Analysis/figures/draft/radial_SWOT_corrcoef_diff.pdf',...
% 'BackgroundColor','none','ContentType','vector')

%% All antennas

figure('Color','w')
tiledlayout(2,length(antenna_IDs),'TileSpacing','tight')

for ii = 1:length(antenna_IDs) % unfiltered
    antenna_ID = antenna_IDs{ii};
    eval(['AX' num2str(ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[-126 -122.] + [-.1 .1],...
        'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval([antenna_ID '.Radial_vel_corrR .* [' ...
                            antenna_ID '.Radial_vel_corrP<P_threshold]./[' ...
                            antenna_ID '.Radial_vel_corrP<P_threshold]']));
    m_plot(real(eval([antenna_ID '.Origin_antenna'])), imag(eval([antenna_ID '.Origin_antenna'])), 'k.','MarkerSize',25)
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1])
    title(['BML: ' antenna_ID])
end
CB.Label.String = ['Vel_{radial}: corr. coef. (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
for ii = 1:length(antenna_IDs) % filtered
    antenna_ID = antenna_IDs{ii};
    eval(['AX' num2str(length(antenna_IDs) + ii) ' = nexttile;'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[-126 -122.] + [-.1 .1],...
        'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval([antenna_ID '.Radial_vel_lowpass_corrR .* [' ...
                            antenna_ID '.Radial_vel_lowpass_corrP<P_threshold]./[' ...
                            antenna_ID '.Radial_vel_lowpass_corrP<P_threshold]']));
    m_plot(real(eval([antenna_ID '.Origin_antenna'])), imag(eval([antenna_ID '.Origin_antenna'])), 'k.','MarkerSize',25)
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1])
    title(['BML: ' antenna_ID])
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

for ii = 1:length(antenna_IDs)
    antenna_ID = antenna_IDs{ii};
    eval(['AX' num2str(ii) ' = nexttile([2,1]);'])
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',[-126 -122.] + [-.1 .1],...
        'latitudes',[37 42] + [-.1 .1]);
    set(gcf,'color','w'); hold on
    m_pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},...
                      eval([antenna_ID '.Radial_vel_lowpass_corrR' ...
                    ' .* [' antenna_ID '.Radial_vel_lowpass_corrP<P_threshold]./[' ... dot comment for low-P areas
                            antenna_ID '.Radial_vel_lowpass_corrP<P_threshold]' ... dot comment for low-P areas
                            ' - ' ...
                            antenna_ID '.Radial_vel_corrR' ...
                    ' .* [' antenna_ID '.Radial_vel_corrP<P_threshold]./[' ... dot comment for low-P areas
                            antenna_ID '.Radial_vel_corrP<P_threshold]' ... dot comment for low-P areas
                            ]));
    m_plot(real(eval([antenna_ID '.Origin_antenna'])), imag(eval([antenna_ID '.Origin_antenna'])), 'k.','MarkerSize',25)
    COAST = m_gshhs_i('patch',0.5*[1 1 1],'FaceAlpha',0.5);
    m_grid('box','fancy', 'backgroundcolor','none');
    set(gca,'FontSize',16)
    CB = colorbar; clim([-1 1]*0.5)
    CB.Label.String = ['U_{radial, lowpass}: corr. coef. (where P<' num2str(P_threshold) ')']; CB.Label.FontSize = 16;
    title(['BML: ' antenna_ID])
end
colormap(bwr)
linkaxes([AX1 AX2 AX3],'xy')

nexttile;
histogram(maxCorrCoef_lowpass - maxCorrCoef,[-1:0.01:1],'Normalization','pdf')
xlabel('CorrCoef_{lowpass} - CorrCoef_{full}')
xlabel('PDF')
title({['Mean = ' num2str(mean(maxCorrCoef_lowpass(:) - maxCorrCoef(:),'omitnan'))] ; ...
       ['Median = ' num2str(median(maxCorrCoef_lowpass(:) - maxCorrCoef(:),'omitnan'))] ; ...
       ['\sigma = ' num2str(std(maxCorrCoef_lowpass(:) - maxCorrCoef(:),'omitnan'))]})

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
%% Along-radial spectra

% This comparison is to see if there is a noise floor in HFR radial
% velocities along (almost) straight lines.

close all

dTHETA = 0.65; % width of angle interval
THETA = -170 + 3; % should be 3 + 5*N

IND = Radial.ANGLE > THETA - dTHETA/2 & Radial.ANGLE < THETA + dTHETA/2;
R_dist_ = Radial.R_dist(IND);
ANGLE_ = Radial.ANGLE(IND);
TIME_ = Radial.TIME(IND);
LON_ = Radial.LON(IND);
LAT_ = Radial.LAT(IND);
VEL_ = Radial.VEL(IND);

R_dist_unique = unique(R_dist_);
ANGLE_unique = unique(ANGLE_);
TIME_unique = unique(TIME_);
LON_unique = real(unique(LON_ + 1i*LAT_));
LAT_unique = imag(unique(LON_ + 1i*LAT_));
VEL_unique = unique(VEL_);

% Filter out the points from neighboring lines that leak in:
% Sort by radial distance:
[~, i_R_unique_sorted] = sort(R_dist_unique);
LON_unique_sorted = LON_unique(i_R_unique_sorted);
LAT_unique_sorted = LAT_unique(i_R_unique_sorted);

% Get rid of very close data:
% Start at the beginning and check the next point: if it's too far off in
% angle, discard them.
i_R_unique_sorted_ = [i_R_unique_sorted(1)];
LastGoodIndex      = [i_R_unique_sorted(1)];
dTHETA_ii = angle([LON_unique_sorted(2) + 1i*LAT_unique_sorted(2)] - ...
                  [LON_unique_sorted(1) + 1i*LAT_unique_sorted(1)])*180/pi;
ii = 0; delta_i = 0;
for jj = 1:[length(i_R_unique_sorted) - 1]
    ii = i_R_unique_sorted(jj);
    if ii == LastGoodIndex
        delta_i = 1;
        if [ii + delta_i] <= [length(i_R_unique_sorted)]
            dTHETA_ii = angle([LON_unique_sorted(ii+delta_i) + 1i*LAT_unique_sorted(ii+delta_i)] - ...
                              [LON_unique_sorted(ii) + 1i*LAT_unique_sorted(ii)])*180/pi;
        else
        end
        while [dTHETA_ii > THETA + 1.2*dTHETA] || [dTHETA_ii < THETA - 1.2*dTHETA]
            delta_i = delta_i + 1;
            if [ii + delta_i] <= [length(i_R_unique_sorted)]
                dTHETA_ii = angle([LON_unique_sorted(ii+delta_i) + 1i*LAT_unique_sorted(ii+delta_i)] - ...
                                  [LON_unique_sorted(ii) + 1i*LAT_unique_sorted(ii)])*180/pi;
            else
                dTHETA_ii = THETA;
            end
        end
        if [ii + delta_i] <= [length(i_R_unique_sorted)]
            dTHETA_ii = angle([LON_unique_sorted(ii+delta_i) + 1i*LAT_unique_sorted(ii+delta_i)] - ...
                              [LON_unique_sorted(ii) + 1i*LAT_unique_sorted(ii)])*180/pi;
        else
            dTHETA_ii = angle([LON_unique_sorted(end) + 1i*LAT_unique_sorted(end)] - ...
                              [LON_unique_sorted(ii) + 1i*LAT_unique_sorted(ii)])*180/pi;
        end
        i_R_unique_sorted_ = [i_R_unique_sorted_ ; i_R_unique_sorted(ii)];
        if [ii + delta_i] <= [length(i_R_unique_sorted)]
            LastGoodIndex = i_R_unique_sorted(ii+delta_i);
        else
            LastGoodIndex = i_R_unique_sorted(end);
        end
    else
        % disp('skipped')
    end
    % dTHETA_ii
end

LON_unique_sorted = LON_unique(i_R_unique_sorted_);
LAT_unique_sorted = LAT_unique(i_R_unique_sorted_);


for tt = 1:length(TIME_unique)
    IND_time = TIME_ == TIME_unique(tt);
    Sum_Instantaneous_Alongtrack_Data(tt) = sum(IND_time);
end


figure
plot(Sum_Instantaneous_Alongtrack_Data,'.-')
ylim([0 length(unique(LON_ + 1i*LAT_))])
%%
close all

figure
plot(Sum_Instantaneous_Alongtrack_Data,'.-')
ylim([0 length(unique(LON_ + 1i*LAT_))])


figure
plot(real(unique(LON_ + 1i*LAT_)),imag(unique(LON_ + 1i*LAT_)),'r*-'); hold on
plot(LON_unique_sorted,LAT_unique_sorted,'b.-');
HH = [ ones(size(LON_unique_sorted)), ...
       LON_unique_sorted - mean(LON_unique_sorted,'omitnan'), ...
      [LON_unique_sorted - mean(LON_unique_sorted,'omitnan')].^2];
plot(LON_unique_sorted, HH*[HH\LAT_unique_sorted],'o--')
plot(real(Radial.Origin_antenna),imag(Radial.Origin_antenna),'sb')
% axis equal

% figure
% plot(real(unique(LON_ + 1i*LAT_)),imag(unique(LON_ + 1i*LAT_)),'.-'); hold on
% HH = [ones(size(real(unique(LON_ + 1i*LAT_)))), ...
%       real(unique(LON_ + 1i*LAT_)) - mean(real(unique(LON_ + 1i*LAT_)),'omitnan'), ...
%       [real(unique(LON_ + 1i*LAT_)) - mean(real(unique(LON_ + 1i*LAT_)),'omitnan')].^2];
% plot(real(unique(LON_ + 1i*LAT_)), HH*[HH\imag(unique(LON_ + 1i*LAT_))],'.-')

% figure
% plot(diff(unique(R_dist_)),'.-')
% 
% figure
% plot(unique(ANGLE_),'.-')

% Plot a single along-view radial velocity data series:
figure
DATENUM = datenum('09-May-2023 12:00:00');

tt = dsearchn(TIME_unique,DATENUM);
IND_time = TIME_ == TIME_unique(tt);

plot(R_dist_(IND_time), VEL_(IND_time), '.-' ); hold on



% Estimate the power spectrum:
Min_Points = 20;
NUFFT = nan(Min_Points,1);
II = 1;
for tt = 1:length(TIME_unique)
    IND_time = TIME_ == TIME_unique(tt);
    if sum(IND_time) >= Min_Points
        DIST = m_lldist(LON_(IND_time),LAT_(IND_time));
        DIST = [0 ; cumsum(DIST)];
        NUFFT(:,II) = nufft(VEL_(IND_time),DIST,(0:(Min_Points-1))/Min_Points);
        II = II + 1;
    else
    end
end
figure
k_Ny = 1/[0.5*mode(diff(DIST))];
semilogy(k_Ny*2*(0:(Min_Points-1))/Min_Points,sum(abs(NUFFT).^2,2),'.-')

cascade_figures(4)

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

function [AA,dd_fit] = fit_6term_2d(xx,yy,dd,nn)
xx.*yy.*dd; % size check
Collimate = @(IN) IN(:);
dd_fit = nan(size(dd));
AA = nan(size(dd,1), size(dd,2), 6);
II = size(xx,2);
JJ = size(xx,1);
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
AA = AA.*repmat(isfinite(dd),1,1,size(AA,3))./repmat(isfinite(dd),1,1,size(AA,3));
dd_fit = dd_fit.*isfinite(dd)./isfinite(dd);
end