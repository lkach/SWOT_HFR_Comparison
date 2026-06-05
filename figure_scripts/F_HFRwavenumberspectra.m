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

% load('./SWOT_and_HFR_velocities_8pix.mat')
load('./SWOT_and_HFR_velocities_8pix_weighted.mat')
% Name                              Size                  Bytes  Class              Attributes
% 
% CG_Iteration                      1x1                       8  double                       
% ConvWithWindow                    1x1                      32  function_handle              
% PIXELS                            1x1                       8  double                       
% SWOT_time                        89x1                     712  double                       
% U_cyclogeostr_Nit               285x69x89            14001480  double                       
% U_geostr                        285x69x89            14001480  double                       
% U_rot_filtered_SWOTtimes         94x56x89             3747968  double                       
% U_rot_unfiltered_SWOTtimes       94x56x89             3747968  double                       
% Ucg_SWOT_HFRgrid_all             94x56x89             3747968  double                       
% Ug_SWOT_HFRgrid_all              94x56x89             3747968  double                       
% V_cyclogeostr_Nit               285x69x89            14001480  double                       
% V_geostr                        285x69x89            14001480  double                       
% V_rot_filtered_SWOTtimes         94x56x89             3747968  double                       
% V_rot_unfiltered_SWOTtimes       94x56x89             3747968  double                       
% Vcg_SWOT_HFRgrid_all             94x56x89             3747968  double                       
% Vg_SWOT_HFRgrid_all              94x56x89             3747968  double                       
% WINDOW                           36x1                     288  double 

load('./UV_rot_filtered.mat')

close all

%% Animate 6km HFR snapshots to get a feel for what they are

% close all
% 
% figure
% set(gcf,'Position',[-1437         100         495         870])
% for ii = 1:length(NORCAL.HFR.time)
%     imagesc(NORCAL.HFR.u(:,:,ii))
%     title(num2str(ii))
%     pause(0.1)
%     clim([-1 1])
% end

%% Examine how complete the snapshots are
% and limit later spectral analysis to these times

close all

N_slice = numel(NORCAL.HFR.u(:,:,1));
N_time = length(NORCAL.HFR.time);

figure
subplot(211)
% plot(NORCAL.HFR.time/24, ...
plot(sum(reshape(isfinite(NORCAL.HFR.u),N_slice,N_time),1,'omitnan')'/N_slice,'.-')
subplot(212)
histogram(sum(reshape(isfinite(NORCAL.HFR.u),N_slice,N_time),1,'omitnan')'/N_slice,0:0.01:0.51)

%% Examine how long vertical and horizontal strips are

close all

VerticalStripLengths = [];
HorizontalStripLengths = [];
Collimate = @(IN) IN(:);

for ii = 1:N_time
    % if sum(Collimate(isfinite(NORCAL.HFR.u(:,:,ii))))/N_slice > 0.4
        VerticalStripLengths   = [VerticalStripLengths;   sum(isfinite(NORCAL.HFR.u(:,:,ii)),1)'];
        HorizontalStripLengths = [HorizontalStripLengths; sum(isfinite(NORCAL.HFR.u(:,:,ii)),2) ];
    % else
    % end
end

figure
histogram(VerticalStripLengths,  0:100); hold on
histogram(HorizontalStripLengths,0:100);
legend('Vertical Strip Lengths','Horizontal Strip Lengths')

figure
imagesc(sum(isfinite(NORCAL.HFR.u),3)./size(NORCAL.HFR.u,3))
clim([0 1])
colormap("turbo")
colorbar
set(gcf,'Position',[-862   271   383   540])

% space steps (km)
dy = abs(median(diff(NORCAL.HFR.lat))*111);
dx = abs(median(diff(NORCAL.HFR.lon))*111*cosd(40));


%% Test applying NUFFT to get spectral estimates in the x and y directions

close all

% space steps (km)
dy = abs(median(diff(NORCAL.HFR.lat))*111);
dx = abs(median(diff(NORCAL.HFR.lon))*111*cosd(40));

% Establish the set of wavenumbers to fit to:
k_vec = [0:[1/[20*dx]]:[1/[2*dx]]]';
l_vec = [0:[1/[30*dy]]:[1/[2*dy]]]';

% Use my custom function "nunanspectrum", provided with this software/data
% bundle:
% [SPEC, Freq, Err] = nunanspectrum(TS, T, TIME_UNITS, ...
% 'Segments',SEGS, 'Window','hann', 'Freq',l_vec);

% Example:
figure
imagesc(NORCAL.HFR.u(:,:,1070))
set(gcf,'Position',[-1437         100         495         870])

figure
plot(NORCAL.HFR.u(:,12:15,1070))

figure
[SPEC, ~, ~] = nunanspectrum(NORCAL.HFR.u(:,12,1070), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
hold on
nunanspectrum(NORCAL.HFR.u(:,13,1070), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
              'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
nunanspectrum(NORCAL.HFR.u(:,14,1070), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
              'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
nunanspectrum(NORCAL.HFR.u(:,15,1070), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
              'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);


% White Noise Test
figure
[SPEC, ~, ~] = nunanspectrum(randn(size(NORCAL.HFR.u(:,12,1070))), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
hold on
SPEC_mat_v = [];
for ii = 1:100
[SPEC, ~, ~] = nunanspectrum(randn(size(NORCAL.HFR.u(:,13,1070))), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
SPEC_mat_v = [SPEC_mat_v , SPEC];
end
loglog(l_vec,mean(SPEC_mat_v,2),'k','LineWidth',2)
title('White noise test - nunanspectrum, randn')

% White Noise Test (1 fewer datum)
figure
[SPEC, ~, ~] = nunanspectrum(randn(size(NORCAL.HFR.u(2:end,12,1070))), 0:dx:[dx*[size(NORCAL.HFR.u,1) - 1]], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
hold on
SPEC_mat_v = [];
for ii = 1:100
[SPEC, ~, ~] = nunanspectrum(randn(size(NORCAL.HFR.u(2:end,12,1070))), 0:dx:[dx*[size(NORCAL.HFR.u,1) - 1]], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
SPEC_mat_v = [SPEC_mat_v , SPEC];
end
loglog(l_vec,mean(SPEC_mat_v,2),'k','LineWidth',2)
title('White noise test (1 fewer datum) - nunanspectrum, randn')


% Red Noise Test
figure
[SPEC, ~, ~] = nunanspectrum(AR_make(0.9,length(NORCAL.HFR.u(:,12,1070))), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
hold on
SPEC_mat_v = [];
for ii = 1:100
[SPEC, ~, ~] = nunanspectrum(AR_make(0.9,length(NORCAL.HFR.u(:,13,1070))), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
SPEC_mat_v = [SPEC_mat_v , SPEC];
end
loglog(l_vec,mean(SPEC_mat_v,2),'k','LineWidth',2)
title('Red noise test - nunanspectrum, AR\_make')

%% More tests

close all

% Red Noise to white noise Test via differentiation
figure
[SPEC, ~, ~] = nunanspectrum(diff(AR_make(0.9,length(NORCAL.HFR.u(:,12,1070)))), 0:dx:[dx*[-1+size(NORCAL.HFR.u,1)]], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
hold on
SPEC_mat_v = [];
for ii = 1:100
[SPEC, ~, ~] = nunanspectrum(diff(AR_make(0.9,length(NORCAL.HFR.u(:,13,1070)))), 0:dx:[dx*[-1+size(NORCAL.HFR.u,1)]], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
SPEC_mat_v = [SPEC_mat_v , SPEC];
end
loglog(l_vec,mean(SPEC_mat_v,2),'k','LineWidth',2)
title('Red-to-white noise test - nunanspectrum, diff(AR\_make)')


% Red Noise with floor Test
WHITE_FLOOR = 0.5;
figure
[SPEC, ~, ~] = nunanspectrum(WHITE_FLOOR*randn(size(NORCAL.HFR.u(:,12,1070))) + AR_make(0.99,length(NORCAL.HFR.u(:,12,1070))), 0:dx:[dx*[size(NORCAL.HFR.u,1)]], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
hold on
SPEC_mat_v = [];
for ii = 1:100
[SPEC, ~, ~] = nunanspectrum(WHITE_FLOOR*randn(size(NORCAL.HFR.u(:,13,1070))) + AR_make(0.99,length(NORCAL.HFR.u(:,13,1070))), 0:dx:[dx*[size(NORCAL.HFR.u,1)]], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
SPEC_mat_v = [SPEC_mat_v , SPEC];
end
loglog(l_vec,mean(SPEC_mat_v,2),'k','LineWidth',2)
title('White noise + Red noise test - nunanspectrum, randn + AR\_make')

% Red Noise with floor to white noise with rise Test via differentiation
figure
[SPEC, ~, ~] = nunanspectrum(diff(WHITE_FLOOR*randn(size(NORCAL.HFR.u(:,12,1070))) + AR_make(0.99,length(NORCAL.HFR.u(:,12,1070)))), 0:dx:[dx*[-1+size(NORCAL.HFR.u,1)]], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
hold on
SPEC_mat_v = [];
for ii = 1:100
[SPEC, ~, ~] = nunanspectrum(diff(WHITE_FLOOR*randn(size(NORCAL.HFR.u(:,13,1070))) + AR_make(0.99,length(NORCAL.HFR.u(:,13,1070)))), 0:dx:[dx*[-1+size(NORCAL.HFR.u,1)]], 'km', ...
               'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
SPEC_mat_v = [SPEC_mat_v , SPEC];
end
loglog(l_vec,mean(SPEC_mat_v,2),'k','LineWidth',2)
title('[White noise + Red noise] to [blue + white] test - nunanspectrum, diff(randn + AR\_make)')

%% Along vertical

% This variable is true if we wish to also evaluate the rotational,
% filtered HFR velocities in the wavenumber domain (they are used elsewhere
% to compare in the time domain):
test_rotational_filtered_hfr = true;

% % % Establish the set of wavenumbers to fit to:
% k_vec = [0:[1/[30*dx]]:[1/[2*dx]]]';
% l_vec = [0:[1/[30*dy]]:[1/[2*dy]]]';

% % % % % % % Analyze the first 17 columns, which tend to be good,
% % % % % % % using only >=81 long columns:

SPEC_mat_u_l = [];
SPEC_mat_v_l = [];
SPEC_mat_u_rf_l = [];
SPEC_mat_v_rf_l = [];
Length_Data = [];
XRANGE = 1:17;
for ii = 1:N_time
    for jj = XRANGE
        HFR_Data_u     = NORCAL.HFR.u(:,jj,ii);
        HFR_Data_v     = NORCAL.HFR.v(:,jj,ii);

        % HFR_Data_v     = 0.5*(NORCAL.HFR.u(:,jj,ii).^2 + NORCAL.HFR.v(:,jj,ii).^2);
        % HFR_Data_v     = cumsum(randn(length(NORCAL.HFR.u(:,jj,ii)),1));%.*[isfinite(NORCAL.HFR.u(:,jj,ii))]./[isfinite(NORCAL.HFR.u(:,jj,ii))]; % test case: red spectrum with no white noise
        % HFR_Data_v     = cumsum(randn(91,1)) + 0.5*randn(91,1);%$
        % HFR_Data_v     = randn(length(NORCAL.HFR.u(:,jj,ii)),1);%.*[isfinite(NORCAL.HFR.u(:,jj,ii))]./[isfinite(NORCAL.HFR.u(:,jj,ii))]; % test case: white noise
        % HFR_Data_v     = randn(90,1);%$

        HFR_Data_u     = trim_nans(HFR_Data_u);
            Length_Data = [Length_Data ; length(HFR_Data_u)];
        HFR_Location = 0:dy:[dy*length(HFR_Data_u)];

        HFR_Data_v     = trim_nans(HFR_Data_v);
            Length_Data = [Length_Data ; length(HFR_Data_v)];
        HFR_Location = 0:dy:[dy*length(HFR_Data_v)];

        % if length(HFR_Data) == 90 && sum(~isfinite(HFR_Data)) == 0 % NONUNIFORM FFT
            % [SPEC_u, ~, ~] = nunanspectrum(HFR_Data_u,HFR_Location, 'km', ...
            %                'Segments',3, 'Window','hanning', 'Freq',l_vec, 'Plot',false, 'PlotSegments',false);
            % [SPEC_v, ~, ~] = nunanspectrum(HFR_Data_v,HFR_Location, 'km', ...
            %                'Segments',3, 'Window','hanning', 'Freq',l_vec, 'Plot',false, 'PlotSegments',false);

        % if [length(HFR_Data_v) >= 81] && ~mod(length(HFR_Data_v),3) && sum(~isfinite(HFR_Data_v)) == 0 % UNIFORM FFT
        if [length(HFR_Data_v) >= 81] && sum(~isfinite(HFR_Data_v)) == 0
            HFR_Data_u = HFR_Data_u(1:81);
            HFR_Data_v = HFR_Data_v(1:81);
            [SPEC_u, l_vec, ~] = nanspectrum(HFR_Data_u, dy, 'km', 3, '.-',false,0,'hanning');
            [SPEC_v, ~,     ~] = nanspectrum(HFR_Data_v, dy, 'km', 3, '.-',false,0,'hanning');

            SPEC_mat_u_l = [SPEC_mat_u_l , SPEC_u];
            SPEC_mat_v_l = [SPEC_mat_v_l , SPEC_v];
        else
        end

        % Rotational filtered HFR data, must be evaluated in a different
        % script first to save space here, or load('./UV_rot_filtered.mat')
        if test_rotational_filtered_hfr
            HFR_Data_u_rf  = U_rot_filtered(:,jj,ii); % ^ see above
            HFR_Data_v_rf  = V_rot_filtered(:,jj,ii);

            HFR_Data_u_rf = trim_nans(HFR_Data_u_rf);
            HFR_Data_v_rf = trim_nans(HFR_Data_v_rf);

            if [length(HFR_Data_u_rf) >= 81] && [length(HFR_Data_v_rf) >= 81] && sum(~isfinite(HFR_Data_u_rf)) == 0
                HFR_Data_u_rf = HFR_Data_u_rf(1:81);
                HFR_Data_v_rf = HFR_Data_v_rf(1:81);
                [SPEC_u_rf, ~, ~] = nanspectrum(HFR_Data_u_rf, dy, 'km', 3, '.-',false,0,'hanning');
                [SPEC_v_rf, ~, ~] = nanspectrum(HFR_Data_v_rf, dy, 'km', 3, '.-',false,0,'hanning');

                SPEC_mat_u_rf_l = [SPEC_mat_u_rf_l , SPEC_u_rf];
                SPEC_mat_v_rf_l = [SPEC_mat_v_rf_l , SPEC_v_rf];
            else
            end
        else
        end




    end

    if ~mod(ii,22)
        disp(100*ii/N_time)
    else
    end
end
% SPEC_mat_u_l = SPEC_mat_u_l/size(SPEC_mat_u_l,2);
% SPEC_mat_v_l = SPEC_mat_v_l/size(SPEC_mat_v_l,2);

close all

figure
subplot(121)
loglog(l_vec,SPEC_mat_v_l(:,1:10:end),'-');hold on
loglog(l_vec,median(SPEC_mat_v_l,2),'k*-','LineWidth',2)
loglog(l_vec,[median(SPEC_mat_v_l(1,:))/l_vec(2)^-2]*[l_vec.^-2],'g--','LineWidth',2)
xlabel('Meridional wavenumber l (km^-1)')
ylabel('PSD ([m/2]^4 [km^-1]^-1)')
subplot(122)
loglog(l_vec,prctile(SPEC_mat_v_l',[50]),'ko-','LineWidth',2); hold on
loglog(l_vec,prctile(SPEC_mat_v_l',[10 90]),'k--','LineWidth',1)

Suu_l_10_50_90 = prctile(SPEC_mat_u_l',[10 50 90]);
Svv_l_10_50_90 = prctile(SPEC_mat_v_l',[10 50 90]);
Suu_l_mean = mean(SPEC_mat_u_l,2);
Svv_l_mean = mean(SPEC_mat_v_l,2);
if test_rotational_filtered_hfr
    Suu_rf_l_10_50_90 = prctile(SPEC_mat_u_rf_l',[10 50 90]);
    Svv_rf_l_10_50_90 = prctile(SPEC_mat_v_rf_l',[10 50 90]);
    Suu_rf_l_mean = mean(SPEC_mat_u_rf_l,2);
    Svv_rf_l_mean = mean(SPEC_mat_v_rf_l,2);
else
end

% Scatter of indivual spectra around their average (ratios in log space)
figure;histogram(log10(SPEC_mat_v_l ./ mean(SPEC_mat_v_l,2)),'Normalization','pdf')
hold on
% Compare the whole distribution to that of a single frequency:
mean_SPEC_mat = mean(SPEC_mat_v_l,2);
histogram(log10(SPEC_mat_v_l(5,:) ./ ...
          mean_SPEC_mat (5) ),'Normalization','pdf')

% Histogram of how long the vertical strips are:
figure
histogram(Length_Data,0:100)
legend('length of vertical strips')

% Map of where you are looking:
figure
Analysis_Block = zeros(size(NORCAL.HFR.v(:,:,1)));
Analysis_Block(:,XRANGE) = 1;
imagesc(sum(isfinite(NORCAL.HFR.v),3)./size(NORCAL.HFR.v,3) + Analysis_Block)
set(gcf,'Position',[-1390 288 344 540])


% % % % % % % Analyze the first 15 columns, which tend to be good:
% figure
% SPEC_mat = [];
% for ii = 1:N_time
%     % if sum(Collimate(isfinite(NORCAL.HFR.u(:,:,ii))))/N_slice > 0.49
%     for jj = 1:15
%         if sum(isfinite(NORCAL.HFR.v(:,jj,ii)))/length(NORCAL.HFR.v(:,jj,ii)) > 0.95
%             [SPEC, ~, ~] = nunanspectrum(NORCAL.HFR.v(:,jj,ii), 0:dy:[dy*size(NORCAL.HFR.u,1)], 'km', ...
%                            'Segments',3, 'Window','rectwin', 'Freq',l_vec, 'Plot',true, 'PlotSegments',false);
%             hold on
%             SPEC_mat = [SPEC_mat , SPEC];
%             % error
%         else
%         end
%     end
% 
%     if ~mod(ii,22)
%         disp(100*ii/N_time)
%     else
%     end
% end
% loglog(l_vec,mean(SPEC_mat,2),'k','LineWidth',2)

% % % % % % % Analyze columns 12 through 15, which tend to be good:
% figure
% SPEC_mat = [];
% for ii = 1:N_time
%     % if sum(Collimate(isfinite(NORCAL.HFR.u(:,:,ii))))/N_slice > 0.49
%     if sum(isfinite(NORCAL.HFR.v(:,15,ii)))/length(NORCAL.HFR.v(:,15,ii)) > 0.95
%         [SPEC1, ~, ~] = nunanspectrum(NORCAL.HFR.v(:,12,ii), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
%                        'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
%         hold on
%         [SPEC2, ~, ~] = nunanspectrum(NORCAL.HFR.v(:,13,ii), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
%                       'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
%         [SPEC3, ~, ~] = nunanspectrum(NORCAL.HFR.v(:,14,ii), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
%                       'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
%         [SPEC4, ~, ~] = nunanspectrum(NORCAL.HFR.v(:,15,ii), 0:dx:[dx*size(NORCAL.HFR.u,1)], 'km', ...
%                       'Segments',2, 'Window','hanning', 'Freq',l_vec, 'Plot',true);
%         SPEC_mat = [SPEC_mat , SPEC1 , SPEC2 , SPEC3 , SPEC4];
%     else
%     end
%     if ~mod(ii,22)
%         disp(100*ii/N_time)
%     else
%     end
% end
% loglog(l_vec,mean(SPEC_mat,2),'k','LineWidth',2)

%% Along horizontal

close all

% % % Establish the set of wavenumbers to fit to:
% k_vec = [0:[1/[20*dx]]:[1/[2*dx]]]';
% l_vec = [0:[1/[30*dy]]:[1/[2*dy]]]';

figure
SPEC_mat_u_k = [];
SPEC_mat_v_k = [];
SPEC_mat_u_rf_k = [];
SPEC_mat_v_rf_k = [];
Length_Data = [];
YRANGE = 50:94;
for ii = 1:N_time
    % if sum(Collimate(isfinite(NORCAL.HFR.u(:,:,ii))))/N_slice > 0.49
    for jj = YRANGE
        HFR_Data_u     = NORCAL.HFR.u(jj,:,ii);
        HFR_Data_v     = NORCAL.HFR.v(jj,:,ii);

        % HFR_Data_v     = diff(NORCAL.HFR.v(jj,:,ii));
        % HFR_Data_v     = diff(cumsum(randn(41,1)));
        % HFR_Data_v     = 0.5*(NORCAL.HFR.u(jj,:,ii).^2 + NORCAL.HFR.v(jj,:,ii).^2);
        
        HFR_Data_u     = trim_nans(HFR_Data_u);
            Length_Data = [Length_Data ; length(HFR_Data_u)];
        HFR_Location = 0:dx:[dx*length(HFR_Data_u)];

        HFR_Data_v     = trim_nans(HFR_Data_v);
            Length_Data = [Length_Data ; length(HFR_Data_v)];
        HFR_Location = 0:dx:[dx*length(HFR_Data_v)];
        % if length(HFR_Data_v) >= 38 % && length(HFR_Data) <= 42
        %     [SPEC_u, ~, ~] = nunanspectrum(HFR_Data_u,HFR_Location, 'km', ...
        %                    'Segments',2, 'Window','rectwin', 'Freq',k_vec, 'Plot',false, 'PlotSegments',false);
        %     [SPEC_v, ~, ~] = nunanspectrum(HFR_Data_v,HFR_Location, 'km', ...
        %                    'Segments',2, 'Window','rectwin', 'Freq',k_vec, 'Plot',false, 'PlotSegments',false);
        
        if [length(HFR_Data_v) >= 38] && sum(~isfinite(HFR_Data_v)) == 0
            HFR_Data_u = HFR_Data_u(1:38);
            HFR_Data_v = HFR_Data_v(1:38);
            [SPEC_u, k_vec, ~]    = nanspectrum(HFR_Data_u, dx, 'km', 2, '.-',false,0,'hanning');
            [SPEC_v, ~,     ~]    = nanspectrum(HFR_Data_v, dx, 'km', 2, '.-',false,0,'hanning');

            SPEC_mat_u_k = [SPEC_mat_u_k , SPEC_u];
            SPEC_mat_v_k = [SPEC_mat_v_k , SPEC_v];
        else
        end

        % Rotational filtered HFR data, must be evaluated in a different
        % script first to save space here:
        if test_rotational_filtered_hfr
            HFR_Data_u_rf = U_rot_filtered(jj,:,ii);
            HFR_Data_v_rf = V_rot_filtered(jj,:,ii);

            HFR_Data_u_rf = trim_nans(HFR_Data_u_rf);
            HFR_Data_v_rf = trim_nans(HFR_Data_v_rf);

            if [length(HFR_Data_u_rf) >= 38] && [length(HFR_Data_v_rf) >= 38] && sum(~isfinite(HFR_Data_v)) == 0
                HFR_Data_u_rf = HFR_Data_u_rf(1:38);
                HFR_Data_v_rf = HFR_Data_v_rf(1:38);
                [SPEC_u_rf, ~, ~] = nanspectrum(HFR_Data_u_rf, dx, 'km', 2, '.-',false,0,'hanning');
                [SPEC_v_rf, ~, ~] = nanspectrum(HFR_Data_v_rf, dx, 'km', 2, '.-',false,0,'hanning');

                SPEC_mat_u_rf_k = [SPEC_mat_u_rf_k , SPEC_u_rf];
                SPEC_mat_v_rf_k = [SPEC_mat_v_rf_k , SPEC_v_rf];
            else
            end
        else
        end

    end

    if ~mod(ii,22)
        disp(100*ii/N_time)
    else
    end
end
% SPEC_mat_u_k = SPEC_mat_u_k/size(SPEC_mat_u_k,2);
% SPEC_mat_v_k = SPEC_mat_v_k/size(SPEC_mat_v_k,2);
% %%
figure
subplot(121)
loglog(k_vec,SPEC_mat_v_k(:,1:10:end),'-');hold on
loglog(k_vec,mean(SPEC_mat_v_k,2),'k*-','LineWidth',2)
loglog(k_vec,[mean(SPEC_mat_v_k(1,:))/k_vec(2)^-2]*[k_vec.^-2],'g--','LineWidth',2)
xlabel('Meridional wavenumber k (km^-1)')
ylabel('PSD ([m/2]^4 [km^-1]^-1)')
subplot(122)
loglog(k_vec,prctile(SPEC_mat_v_k',[50]),'ko-','LineWidth',2); hold on
loglog(k_vec,prctile(SPEC_mat_v_k',[10 90]),'k--','LineWidth',1)

Suu_k_10_50_90 = prctile(SPEC_mat_u_k',[10 50 90]);
Svv_k_10_50_90 = prctile(SPEC_mat_v_k',[10 50 90]);
Suu_k_mean = mean(SPEC_mat_u_k,2);
Svv_k_mean = mean(SPEC_mat_v_k,2);
if test_rotational_filtered_hfr
    Suu_rf_k_10_50_90 = prctile(SPEC_mat_u_rf_k',[10 50 90]);
    Svv_rf_k_10_50_90 = prctile(SPEC_mat_v_rf_k',[10 50 90]);
    Suu_rf_k_mean = mean(SPEC_mat_u_rf_k,2);
    Svv_rf_k_mean = mean(SPEC_mat_v_rf_k,2);
else
end

% figure;histogram(log10(SPEC_mat ./ mean(SPEC_mat,2)),'Normalization','pdf')
% hold on
% histogram(log10(SPEC_mat(5,:) ./ mean_SPEC_mat(5)),'Normalization','pdf')

% %% Histogram of how long the horizontal strips are:
figure
histogram(Length_Data,0:100)
legend('length of horizontal strips')

% Map of where you are looking:
figure
Analysis_Block = zeros(size(NORCAL.HFR.v(:,:,1)));
Analysis_Block(YRANGE,:) = 1;
imagesc(sum(isfinite(NORCAL.HFR.v),3)./size(NORCAL.HFR.v,3) + Analysis_Block)
set(gcf,'Position',[-1390 288 344 540])

%% Plot the average spectra together

close all
CO = colororder;

figure('Color','w')
loglog(k_vec,Suu_k_10_50_90(2,:),    'o-','LineWidth',2,'Color',CO(1,:)); hold on
loglog(k_vec,Suu_k_10_50_90([1 3],:),'--','LineWidth',1,'Color',CO(1,:))
loglog(k_vec,Svv_k_10_50_90(2,:),    'o-','LineWidth',2,'Color',CO(2,:)); hold on
loglog(k_vec,Svv_k_10_50_90([1 3],:),'--','LineWidth',1,'Color',CO(2,:))

loglog(l_vec,Suu_l_10_50_90(2,:),    'o-','LineWidth',2,'Color',CO(1,:)); hold on
loglog(l_vec,Suu_l_10_50_90([1 3],:),'--','LineWidth',1,'Color',CO(1,:))
loglog(l_vec,Svv_l_10_50_90(2,:),    'o-','LineWidth',2,'Color',CO(2,:)); hold on
loglog(l_vec,Svv_l_10_50_90([1 3],:),'--','LineWidth',1,'Color',CO(2,:))


figure('Color','w')
subplot(121)
% loglog(k_vec,Svv_k_10_50_90(2,:),    'o-','LineWidth',2,'Color',CO(1,:)); hold on
loglog(k_vec,Svv_k_mean,    'o-','LineWidth',2,'Color',CO(1,:)); hold on
loglog(k_vec,Svv_k_10_50_90([1 3],:),'--','LineWidth',1,'Color',CO(1,:),'HandleVisibility','off')
% loglog(l_vec,Svv_l_10_50_90(2,:),    'o-','LineWidth',2,'Color',CO(2,:))
loglog(l_vec,Svv_l_mean,    'o-','LineWidth',2,'Color',CO(2,:))
loglog(l_vec,Svv_l_10_50_90([1 3],:),'--','LineWidth',1,'Color',CO(2,:),'HandleVisibility','off')
LEG = legend('$S_{vv}(k)$','$S_{vv}(l)$');
LEG.Interpreter = 'latex';
LEG.FontSize = 16;
set(gca,'FontSize',16)
xlabel('wavenumber (cpkm)','Interpreter','latex')
ylabel('PSD (m s$^{-1}$ cpkm$^{-1}$)','Interpreter','latex')

subplot(122)
% loglog(1./k_vec,Svv_k_10_50_90(2,:),    'o-','LineWidth',2,'Color',CO(1,:)); hold on
loglog(1./k_vec,Svv_k_mean,    'o-','LineWidth',2,'Color',CO(1,:)); hold on
loglog(1./k_vec,Svv_k_10_50_90([1 3],:),'--','LineWidth',1,'Color',CO(1,:),'HandleVisibility','off')
% loglog(1./l_vec,Svv_l_10_50_90(2,:),    'o-','LineWidth',2,'Color',CO(2,:))
loglog(1./l_vec,Svv_l_mean,    'o-','LineWidth',2,'Color',CO(2,:))
loglog(1./l_vec,Svv_l_10_50_90([1 3],:),'--','LineWidth',1,'Color',CO(2,:),'HandleVisibility','off')
LEG = legend('$S_{vv}(k^{-1})$','$S_{vv}(l^{-1})$');
LEG.Interpreter = 'latex';
LEG.FontSize = 16;
set(gca,'FontSize',16)
xlabel('wavelength (km)','Interpreter','latex')
xticks([10:10:50 100 150])
ylabel('PSD (m s$^{-1}$ cpkm$^{-1}$)','Interpreter','latex')

set(gcf,'Position',[1         257        1251         540])

%%
% Plot time-domain covariance to estimate decorrelation scales for the
% purposes of estimating degrees of freedom for spectral error later:
figure
reshaped_u = reshape(NORCAL.HFR.u,[size(NORCAL.HFR.u,1)*size(NORCAL.HFR.u,2) size(NORCAL.HFR.u,3)]);
reshaped_v = reshape(NORCAL.HFR.v,[size(NORCAL.HFR.v,1)*size(NORCAL.HFR.v,2) size(NORCAL.HFR.v,3)]);
HFR_fracgooddata = sum(isfinite(reshaped_v),2)/size(reshaped_v,2);
II = 0;
CORRu = 0;
CORRv = 0;
for ii = 1:length(HFR_fracgooddata)
    if HFR_fracgooddata(ii) > 0.9
        CORRu_ii = cov_gaps(reshaped_u(ii,:),reshaped_u(ii,:),100);
        CORRv_ii = cov_gaps(reshaped_v(ii,:),reshaped_v(ii,:),100);
        CORRu = CORRu + CORRu_ii;
        CORRv = CORRv + CORRv_ii;
        II = II + 1;
    else
    end
end
CORRu = CORRu/II;
CORRv = CORRv/II;
plot(0:[length(CORRu)-1],CORRu,'.-'); hold on
plot(0:[length(CORRv)-1],CORRv,'.-')
% 8 Hours seems reasonable (the rebound is clearly tidal)

%% Analogy for SWOT velocities


close all

% Limit analysis to these indices: (137:277,[12:27,42:58])

SPEC_mat_ug_l = [];
SPEC_mat_vg_l = [];
Length_Data_SWOTvel = [];
XRANGE = [12:27,42:55];
YRANGE = 137:277;
dy_swot = 2;
for ii = 1:size(U_geostr,3)
    for jj = XRANGE
        SWOTvel_Data_u     = U_geostr(YRANGE,jj,ii);
        SWOTvel_Data_v     = V_geostr(YRANGE,jj,ii);

        Length_Data_SWOTvel = [Length_Data_SWOTvel ; length(SWOTvel_Data_u)];
        SWOT_Location = 0:dy_swot:[dy_swot*length(SWOTvel_Data_u)];

        % if sum(~isfinite(HFR_Data)) == 0 % NONUNIFORM FFT
        %    [SPEC_u, ~, ~] = nunanspectrum(SWOTvel_Data_u,SWOT_Location, 'km', ...
        %                   'Segments',3, 'Window','hanning', 'Freq',l_vec, 'Plot',false, 'PlotSegments',false);
        %    [SPEC_v, ~, ~] = nunanspectrum(SWOTvel_Data_v,SWOT_Location, 'km', ...
        %                   'Segments',3, 'Window','hanning', 'Freq',l_vec, 'Plot',false, 'PlotSegments',false);

        if sum(~isfinite(SWOTvel_Data_v)) <= 2 % == 0
            % SWOTvel_Data_u = SWOTvel_Data_u;
            % SWOTvel_Data_v = SWOTvel_Data_v;
            [SPEC_u, l_vec_swot, ~] = nanspectrum(SWOTvel_Data_u, dy_swot, 'km', 3, '.-',false,0,'hanning');
            [SPEC_v, ~,          ~] = nanspectrum(SWOTvel_Data_v, dy_swot, 'km', 3, '.-',false,0,'hanning');

            SPEC_mat_ug_l = [SPEC_mat_ug_l , SPEC_u];
            SPEC_mat_vg_l = [SPEC_mat_vg_l , SPEC_v];
        else
        end
    end


    % disp(100*ii/size(U_geostr,3))

end
% SPEC_mat_ug_l = SPEC_mat_ug_l/size(SPEC_mat_ug_l,2);
% SPEC_mat_vg_l = SPEC_mat_vg_l/size(SPEC_mat_vg_l,2);

figure
subplot(221)
loglog(l_vec_swot,SPEC_mat_ug_l(:,1:10:end),'-');hold on
loglog(l_vec_swot,median(SPEC_mat_ug_l,2),'k*-','LineWidth',2)
loglog(l_vec_swot,[median(SPEC_mat_ug_l(1,:))/l_vec_swot(2)^-2]*[l_vec_swot.^-2],'g--','LineWidth',2)
loglog(l_vec_swot,[median(SPEC_mat_ug_l(1,:))/l_vec_swot(2)^-3]*[l_vec_swot.^-3],'g:','LineWidth',2)
xlabel('Meridional wavenumber l (km^-1)')
ylabel('PSD ([m/2]^4 [km^-1]^-1)')
ylim(10.^[-7 1])
title('Weighted') %$#
subplot(222)
loglog(l_vec_swot,prctile(SPEC_mat_ug_l',[50]),'ko-','LineWidth',2); hold on
loglog(l_vec_swot,prctile(SPEC_mat_ug_l',[10 90]),'k--','LineWidth',1)
ylim(10.^[-7 1])
title('S_{u_gu_g}')

subplot(223)
loglog(l_vec_swot,SPEC_mat_vg_l(:,1:10:end),'-');hold on
loglog(l_vec_swot,median(SPEC_mat_vg_l,2),'k*-','LineWidth',2)
loglog(l_vec_swot,[median(SPEC_mat_vg_l(1,:))/l_vec_swot(2)^-2]*[l_vec_swot.^-2],'g--','LineWidth',2)
loglog(l_vec_swot,[median(SPEC_mat_vg_l(1,:))/l_vec_swot(2)^-3]*[l_vec_swot.^-3],'g:','LineWidth',2)
xlabel('Meridional wavenumber l (km^-1)')
ylabel('PSD ([m/2]^4 [km^-1]^-1)')
ylim(10.^[-7 1])
subplot(224)
loglog(l_vec_swot,prctile(SPEC_mat_vg_l',[50]),'ko-','LineWidth',2); hold on
loglog(l_vec_swot,prctile(SPEC_mat_vg_l',[10 90]),'k--','LineWidth',1)
ylim(10.^[-7 1])
title('S_{v_gv_g}')

Sugug_l_10_50_90 = prctile(SPEC_mat_ug_l',[10 50 90]); Sugug_l_mean = mean(SPEC_mat_ug_l,2);
Svgvg_l_10_50_90 = prctile(SPEC_mat_vg_l',[10 50 90]); Svgvg_l_mean = mean(SPEC_mat_vg_l,2);

% Scatter of indivual spectra around their average (ratios in log space)
figure;histogram(log10(SPEC_mat_vg_l ./ mean(SPEC_mat_vg_l,2)),'Normalization','pdf')
hold on
% Compare the whole distribution to that of a single frequency:
mean_SPEC_mat = mean(SPEC_mat_vg_l,2);
histogram(log10(SPEC_mat_vg_l(5,:) ./ ...
          mean_SPEC_mat (5) ),'Normalization','pdf')

% Histogram of how long the vertical strips are:
figure
histogram(Length_Data_SWOTvel,0:300)
legend('length of vertical strips')

% Map of where you are looking:
figure
Analysis_Block = zeros(size(U_geostr(:,:,1)));
Analysis_Block(YRANGE,XRANGE) = 1;
imagesc(sum(isfinite(U_geostr),3)./size(U_geostr,3) + Analysis_Block)
set(gcf,'Position',[-1390 288 344 540])

% Show the angle of the swath in this region (~73 degrees is typical)
figure
subplot(1,2,1)
plot(angle(   [NORCAL.SWOT.lon{1}(1:[end-1],:) - NORCAL.SWOT.lon{1}(2:end,:)] + ...
           1i*[NORCAL.SWOT.lat{1}(1:[end-1],:) - NORCAL.SWOT.lat{1}(2:end,:)])*180/pi)
title('angles')
subplot(1,2,2)
scatter(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},'.');axis equal
hold on; plot([-124, -123],[37, 37 + (1*sind(73)/sind(90-73))],'k*-')
title('tracks and a 73 degree angled line for comparison')

%% Analogy for SWOT velocities that are gridded to the HFR grid

figure
SPEC_mat_ug_gridded_l = [];
SPEC_mat_vg_gridded_l = [];
for i_t = 1:size(Ug_SWOT_HFRgrid_all,3)
    for i_x = 1:size(Ug_SWOT_HFRgrid_all,2)
        SWOT_Data_Ug_gridded = trim_nans(Ug_SWOT_HFRgrid_all(:,i_x,i_t));
        SWOT_Data_Vg_gridded = trim_nans(Vg_SWOT_HFRgrid_all(:,i_x,i_t));
        if [length(SWOT_Data_Ug_gridded) >= 42] && sum(~isfinite(SWOT_Data_Ug_gridded)) == 0
            plot(SWOT_Data_Ug_gridded,'-');hold on
            [SPEC_u,l_vec_swot_gridded] = nanspectrum(SWOT_Data_Ug_gridded(1:42), dy, 'km', 2, '.-',false,0,'hanning');
            [SPEC_v,                 ~] = nanspectrum(SWOT_Data_Vg_gridded(1:42), dy, 'km', 2, '.-',false,0,'hanning');
            SPEC_mat_ug_gridded_l = [SPEC_mat_ug_gridded_l , SPEC_u];
            SPEC_mat_vg_gridded_l = [SPEC_mat_vg_gridded_l , SPEC_v];
        else
        end
    end
end
% SPEC_mat_ug_gridded_l = SPEC_mat_ug_gridded_l/size(SPEC_mat_ug_gridded_l,2);
% SPEC_mat_vg_gridded_l = SPEC_mat_vg_gridded_l/size(SPEC_mat_vg_gridded_l,2);


Sugug_gridded_l_10_50_90 = prctile(SPEC_mat_ug_gridded_l',[10 50 90]);
Svgvg_gridded_l_10_50_90 = prctile(SPEC_mat_vg_gridded_l',[10 50 90]);

Sugug_gridded_l_mean = mean(SPEC_mat_ug_gridded_l,2);
Svgvg_gridded_l_mean = mean(SPEC_mat_vg_gridded_l,2);

% SWOT_Data_Ug_gridded = trim_nans(Ug_SWOT_HFRgrid_all(:,17,1));
% figure
% plot(SWOT_Data_Ug_gridded,'.-')

%% Along vertical spectra of filtered, rotational HFR

% % % Establish the set of wavenumbers to fit to:
% k_vec = [0:[1/[30*dx]]:[1/[2*dx]]]';
% l_vec_rotfilt = [0:[1/[30*dy]]:[1/[2*dy]]]';

% % % % % % % Analyze the first 15 columns, which tend to be good,
% % % % % % % using only >=90 long columns:

SPEC_mat_urotfilt_l = [];
SPEC_mat_vrotfilt_l = [];
Length_Data = [];
XRANGE = 1:17;
for ii = 1:size(U_rot_filtered_SWOTtimes,3)
    for jj = XRANGE
        HFR_Data_u     = U_rot_filtered_SWOTtimes(:,jj,ii);
        HFR_Data_v     = V_rot_filtered_SWOTtimes(:,jj,ii);

        HFR_Data_u     = trim_nans(HFR_Data_u);
            Length_Data = [Length_Data ; length(HFR_Data_u)];
        HFR_Location = 0:dy:[dy*length(HFR_Data_u)];

        HFR_Data_v     = trim_nans(HFR_Data_v);
            Length_Data = [Length_Data ; length(HFR_Data_v)];
        HFR_Location = 0:dy:[dy*length(HFR_Data_v)];

        % if length(HFR_Data) == 90 && sum(~isfinite(HFR_Data)) == 0 % NONUNIFORM FFT
            % [SPEC_u, ~, ~] = nunanspectrum(HFR_Data_u,HFR_Location, 'km', ...
            %                'Segments',3, 'Window','hanning', 'Freq',l_vec_rotfilt, 'Plot',false, 'PlotSegments',false);
            % [SPEC_v, ~, ~] = nunanspectrum(HFR_Data_v,HFR_Location, 'km', ...
            %                'Segments',3, 'Window','hanning', 'Freq',l_vec_rotfilt, 'Plot',false, 'PlotSegments',false);

        % if [length(HFR_Data_v) >= 81] && ~mod(length(HFR_Data_v),3) && sum(~isfinite(HFR_Data_v)) == 0 % UNIFORM FFT
        if [length(HFR_Data_u) >= 81] && [length(HFR_Data_v) >= 81] && sum(~isfinite(HFR_Data_v)) == 0
            HFR_Data_u = HFR_Data_u(1:81);
            HFR_Data_v = HFR_Data_v(1:81);
            [SPEC_u, l_vec_rotfilt, ~] = nanspectrum(HFR_Data_u, dy, 'km', 3, '.-',false,0,'hanning');
            [SPEC_v, ~,     ~] = nanspectrum(HFR_Data_v, dy, 'km', 3, '.-',false,0,'hanning');

            SPEC_mat_urotfilt_l = [SPEC_mat_urotfilt_l , SPEC_u];
            SPEC_mat_vrotfilt_l = [SPEC_mat_vrotfilt_l , SPEC_v];
        else
        end
    end

    if ~mod(ii,22)
        disp(100*ii/size(U_rot_filtered_SWOTtimes,3))
    else
    end
end

close all

figure
subplot(121)
loglog(l_vec_rotfilt,SPEC_mat_vrotfilt_l(:,1:10:end),'-');hold on
loglog(l_vec_rotfilt,median(SPEC_mat_vrotfilt_l,2),'k*-','LineWidth',2)
loglog(l_vec_rotfilt,[median(SPEC_mat_vrotfilt_l(1,:))/l_vec_rotfilt(2)^-2]*[l_vec_rotfilt.^-2],'g--','LineWidth',2)
xlabel('Meridional wavenumber l (km^-1)')
ylabel('PSD ([m/2]^4 [km^-1]^-1)')
subplot(122)
loglog(l_vec_rotfilt,prctile(SPEC_mat_vrotfilt_l',[50]),'ko-','LineWidth',2); hold on
loglog(l_vec_rotfilt,prctile(SPEC_mat_vrotfilt_l',[10 90]),'k--','LineWidth',1)

Suu_rotfilt_l_10_50_90 = prctile(SPEC_mat_urotfilt_l',[10 50 90]); Suu_rotfilt_l_mean = mean(SPEC_mat_urotfilt_l,2);
Svv_rotfilt_l_10_50_90 = prctile(SPEC_mat_vrotfilt_l',[10 50 90]); Svv_rotfilt_l_mean = mean(SPEC_mat_vrotfilt_l,2);

% Scatter of indivual spectra around their average (ratios in log space)
figure;histogram(log10(SPEC_mat_vrotfilt_l ./ mean(SPEC_mat_vrotfilt_l,2)),'Normalization','pdf')
hold on
% Compare the whole distribution to that of a single frequency:
mean_SPEC_mat = mean(SPEC_mat_vrotfilt_l,2);
histogram(log10(SPEC_mat_vrotfilt_l(5,:) ./ ...
          mean_SPEC_mat (5) ),'Normalization','pdf')

% Histogram of how long the vertical strips are:
figure
histogram(Length_Data,0:100)
legend('length of vertical strips')

% Map of where you are looking:
figure
Analysis_Block = zeros(size(NORCAL.HFR.v(:,:,1)));
Analysis_Block(:,XRANGE) = 1;
imagesc(sum(isfinite(NORCAL.HFR.v),3)./size(NORCAL.HFR.v,3) + Analysis_Block)
set(gcf,'Position',[-1390 288 344 540])

%% Along horizontal spectra of filtered, rotational HFR

close all

% % % Establish the set of wavenumbers to fit to:
% k_vec = [0:[1/[20*dx]]:[1/[2*dx]]]';
% l_vec = [0:[1/[30*dy]]:[1/[2*dy]]]';

figure
SPEC_mat_urotfilt_k = [];
SPEC_mat_vrotfilt_k = [];
Length_Data = [];
YRANGE = 50:94;
for ii = 1:size(U_rot_filtered_SWOTtimes,3)
    % if sum(Collimate(isfinite(NORCAL.HFR.u(:,:,ii))))/N_slice > 0.49
    for jj = YRANGE
        HFR_Data_u     = U_rot_filtered_SWOTtimes(jj,:,ii);
        HFR_Data_v     = V_rot_filtered_SWOTtimes(jj,:,ii);
        
        HFR_Data_u     = trim_nans(HFR_Data_u);
            Length_Data = [Length_Data ; length(HFR_Data_u)];
        HFR_Location = 0:dx:[dx*length(HFR_Data_u)];

        HFR_Data_v     = trim_nans(HFR_Data_v);
            Length_Data = [Length_Data ; length(HFR_Data_v)];
        HFR_Location = 0:dx:[dx*length(HFR_Data_v)];
        % if length(HFR_Data_v) >= 38 % && length(HFR_Data) <= 42
        %     [SPEC_u, ~, ~] = nunanspectrum(HFR_Data_u,HFR_Location, 'km', ...
        %                    'Segments',2, 'Window','rectwin', 'Freq',k_vec, 'Plot',false, 'PlotSegments',false);
        %     [SPEC_v, ~, ~] = nunanspectrum(HFR_Data_v,HFR_Location, 'km', ...
        %                    'Segments',2, 'Window','rectwin', 'Freq',k_vec, 'Plot',false, 'PlotSegments',false);
        
        if [length(HFR_Data_u) >= 38] && [length(HFR_Data_v) >= 38] && sum(~isfinite(HFR_Data_v)) == 0
            HFR_Data_u = HFR_Data_u(1:38);
            HFR_Data_v = HFR_Data_v(1:38);
            [SPEC_u, k_vec_rotfilt,  ~]    = nanspectrum(HFR_Data_u, dx, 'km', 2, '.-',false,0,'hanning');
            [SPEC_v,             ~,  ~]    = nanspectrum(HFR_Data_v, dx, 'km', 2, '.-',false,0,'hanning');

            SPEC_mat_urotfilt_k = [SPEC_mat_urotfilt_k , SPEC_u];
            SPEC_mat_vrotfilt_k = [SPEC_mat_vrotfilt_k , SPEC_v];
        else
        end
    end

    disp(100*ii/size(U_rot_filtered_SWOTtimes,3))
end
% %%
figure
subplot(121)
loglog(k_vec_rotfilt,SPEC_mat_vrotfilt_k(:,1:10:end),'-');hold on
loglog(k_vec_rotfilt,mean(SPEC_mat_vrotfilt_k,2),'k*-','LineWidth',2)
loglog(k_vec_rotfilt,[mean(SPEC_mat_vrotfilt_k(1,:))/k_vec_rotfilt(2)^-2]*[k_vec_rotfilt.^-2],'g--','LineWidth',2)
xlabel('Meridional wavenumber k (km^-1)')
ylabel('PSD ([m/2]^4 [km^-1]^-1)')
subplot(122)
loglog(k_vec_rotfilt,prctile(SPEC_mat_vrotfilt_k',[50]),'ko-','LineWidth',2); hold on
loglog(k_vec_rotfilt,prctile(SPEC_mat_vrotfilt_k',[10 90]),'k--','LineWidth',1)

Suu_rotfilt_k_10_50_90 = prctile(SPEC_mat_urotfilt_k',[10 50 90]); Suu_rotfilt_k_mean = mean(SPEC_mat_urotfilt_k,2);
Svv_rotfilt_k_10_50_90 = prctile(SPEC_mat_vrotfilt_k',[10 50 90]); Svv_rotfilt_k_mean = mean(SPEC_mat_vrotfilt_k,2);

% %% Histogram of how long the horizontal strips are:
figure
histogram(Length_Data,0:100)
legend('length of horizontal strips')

% Map of where you are looking:
figure
Analysis_Block = zeros(size(NORCAL.HFR.v(:,:,1)));
Analysis_Block(YRANGE,:) = 1;
imagesc(sum(isfinite(NORCAL.HFR.v),3)./size(NORCAL.HFR.v,3) + Analysis_Block)
set(gcf,'Position',[-1390 288 344 540])

%% ERROR ESTIMATION:
% % In short, do this calculation for 95% confidence ratio:
% EffectiveSegments = 2*SEGMENTS - 1;% degrees of freedom; for overlapping
%                      % Hann windowed segments, this is the total number of
%                      % segments, i.e.:
%                      % 2*SEGMENTS - 1
%                      % This may only work perfectly with the Hann window.
%                      % With other windows, it's not entirely clear, but my
%                      % guess is, if the window is "Hann-like", it will not
%                      % make a big difference.
% 
% err_high = 2*EffectiveSegments/chi2inv(.05/2,2*EffectiveSegments);
% err_low = 2*EffectiveSegments/chi2inv(1-.05/2,2*EffectiveSegments);
% err = [err_low err_high];
EffectiveSegments_HFR_k  =         2*size(SPEC_mat_u_k,         2)*[[dx     /(20)].^2]*[1/8];
EffectiveSegments_HFR_l  =         3*size(SPEC_mat_u_l,         2)*[[dy     /(20)].^2]*[1/8];
EffectiveSegments_SWOT_l =         3*size(SPEC_mat_ug_l,        2)*[[dy_swot/(20)].^2]*[1/1];
EffectiveSegments_SWOT_gridded_l = 2*size(SPEC_mat_ug_gridded_l,2)*[[dy     /(20)].^2]*[1/1];
% % % ^ Assume decorrelation scale of 20 km so as to be conservative with
% % % uncertainty. Assume decorrelation time of 24 hours in order for HFR
% % % to be comparable to SWOT's daily sampling.
% 
err_high_HFR_k = 2*EffectiveSegments_HFR_k/chi2inv(.05/2,2*EffectiveSegments_HFR_k);
err_low_HFR_k = 2*EffectiveSegments_HFR_k/chi2inv(1-.05/2,2*EffectiveSegments_HFR_k);
err_HFR_k = [err_low_HFR_k err_high_HFR_k]/err_high_HFR_k;
%
err_high_HFR_l = 2*EffectiveSegments_HFR_l/chi2inv(.05/2,2*EffectiveSegments_HFR_l);
err_low_HFR_l = 2*EffectiveSegments_HFR_l/chi2inv(1-.05/2,2*EffectiveSegments_HFR_l);
err_HFR_l = [err_low_HFR_l err_high_HFR_l]/err_high_HFR_l;
%
err_high_SWOT_l = 2*EffectiveSegments_SWOT_l/chi2inv(.05/2,2*EffectiveSegments_SWOT_l);
err_low_SWOT_l = 2*EffectiveSegments_SWOT_l/chi2inv(1-.05/2,2*EffectiveSegments_SWOT_l);
err_SWOT_l = [err_low_SWOT_l err_high_SWOT_l]/err_high_SWOT_l;
%
err_high_SWOT_gridded_l = 2*EffectiveSegments_SWOT_gridded_l/chi2inv(.05/2,2*EffectiveSegments_SWOT_gridded_l);
err_low_SWOT_gridded_l = 2*EffectiveSegments_SWOT_gridded_l/chi2inv(1-.05/2,2*EffectiveSegments_SWOT_gridded_l);
err_SWOT_gridded_l = [err_low_SWOT_gridded_l err_high_SWOT_gridded_l]/err_high_SWOT_gridded_l;
%
% % Add to the plot:
% loglog(1.00*[1 1]*10^-1,err_HFR_k *10^-3,'k','LineWidth',1,'HandleVisibility','off')
% loglog(1.20*[1 1]*10^-1,err_HFR_l *10^-3,'k','LineWidth',1,'HandleVisibility','off')
% loglog(1.44*[1 1]*10^-1,err_SWOT_l*10^-3,'k','LineWidth',1,'HandleVisibility','off')

error('Forced stop by user.')

%% PLOT: Combine and make publication-worthy plot

k_cutoff = 1/18; % wavenumber of vertical dashed line
XLIM = 1./[100 10];
YLIM = [2*10^-3 10.^0];
    YLIM = [10^-6.5 10.^0];

close all
% % % % % % % % % % % % % % % % % % % % % % % TWO-PANEL:

figure('Color','w')
% tiledlayout(1,2,"TileSpacing","tight")

% AX1 = nexttile;
AX1 = subplot(1,2,1);
loglog(k_vec,Suu_k_mean,'.-',         'LineWidth',1,'MarkerSize',10,'Color',[0 114 178  ]/255); hold on
loglog(l_vec,Suu_l_mean,'.-',         'LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0 ]/255);
    % loglog(k_vec_rotfilt,Suu_rotfilt_k_mean,'.:', 'LineWidth',2,'MarkerSize',10,'Color',[0 114 178  ]/255,'HandleVisibility','off'); hold on
    % loglog(l_vec_rotfilt,Suu_rotfilt_l_mean,'.:', 'LineWidth',2,'MarkerSize',10,'Color',[213 94/2 0 ]/255,'HandleVisibility','off'); hold on
    loglog(k_vec_rotfilt,Suu_rf_k_mean,'.:', 'LineWidth',2,'MarkerSize',10,'Color',[0 114 178  ]/255,'HandleVisibility','off'); hold on
    loglog(l_vec_rotfilt,Suu_rf_l_mean,'.:', 'LineWidth',2,'MarkerSize',10,'Color',[213 94/2 0 ]/255,'HandleVisibility','off'); hold on
loglog(l_vec_swot,Sugug_l_mean,'.-',  'LineWidth',1,'MarkerSize',10,'Color',[0 0 0      ]/255);
    loglog(l_vec_swot_gridded,Sugug_gridded_l_mean,'.-',  'LineWidth',1,'MarkerSize',10,'Color',[150 150 150]/255);
loglog([1 1]*k_cutoff,YLIM,'k--')
loglog(1.000*[1 1]*10^-1,err_HFR_k         *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.200*[1 1]*10^-1,err_HFR_l         *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.440*[1 1]*10^-1,err_SWOT_l        *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.728*[1 1]*10^-1,err_SWOT_gridded_l*10^-3,'k','LineWidth',1,'HandleVisibility','off')
text(1.05*k_cutoff,10^-1.5,['1/(' num2str(1/k_cutoff) ' km)'],'Interpreter','latex','FontSize',12)
LEG1 = legend('$\langle S_{uu}(k)\rangle$',...
              '$\langle S_{uu}(l)\rangle$',...
              '$\langle S_{u_\mathrm{g}u_\mathrm{g}}(l^\star)\rangle$',...
              '$\langle S_{u''_\mathrm{g}u''_\mathrm{g}}(l)\rangle$');
LEG1.Interpreter = 'latex';
LEG1.FontSize = 16;
LEG1.Position = [0.1419    0.1874    0.1686    0.2871]; % [0.1433    0.1885    0.1686    0.2159]; % [0.1494    0.1850    0.1852    0.2159];
xlabel('$k,l,l^\star$ (cpkm)','Interpreter','latex')
ylabel('PSD (m$^2$ s$^{-2}$ cpkm$^{-1}$)','Interpreter','latex')
AX1.FontSize = 12;
AX1.XLim = [0.9*min([k_vec(1) l_vec(1) l_vec_swot(1)]), 1.1*max([k_vec(end) l_vec(end) l_vec_swot(end)])];
AX1.YLim = YLIM;
grid on

% AX2 = nexttile;
AX2 = subplot(1,2,2);
loglog(k_vec,Svv_k_mean,'.-',         'LineWidth',1,'MarkerSize',10,'Color',[0 114 178  ]/255); hold on
loglog(l_vec,Svv_l_mean,'.-',         'LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0 ]/255);
    % loglog(k_vec_rotfilt,Svv_rotfilt_k_mean,'.:', 'LineWidth',2,'MarkerSize',10,'Color',[0 114 178  ]/255,'HandleVisibility','off'); hold on
    % loglog(l_vec_rotfilt,Svv_rotfilt_l_mean,'.:', 'LineWidth',2,'MarkerSize',10,'Color',[213 94/2 0 ]/255,'HandleVisibility','off'); hold on
    loglog(k_vec_rotfilt,Svv_rf_k_mean,'.:', 'LineWidth',2,'MarkerSize',10,'Color',[0 114 178  ]/255,'HandleVisibility','off'); hold on
    loglog(l_vec_rotfilt,Svv_rf_l_mean,'.:', 'LineWidth',2,'MarkerSize',10,'Color',[213 94/2 0 ]/255,'HandleVisibility','off'); hold on
loglog(l_vec_swot,Svgvg_l_mean,'.-',  'LineWidth',1,'MarkerSize',10,'Color',[0 0 0      ]/255);
    loglog(l_vec_swot_gridded,Svgvg_gridded_l_mean,'.-',  'LineWidth',1,'MarkerSize',10,'Color',[150 150 150]/255);
loglog([1 1]*k_cutoff,YLIM,'k--')
loglog(1.000*[1 1]*10^-1,err_HFR_k         *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.200*[1 1]*10^-1,err_HFR_l         *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.440*[1 1]*10^-1,err_SWOT_l        *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.728*[1 1]*10^-1,err_SWOT_gridded_l*10^-3,'k','LineWidth',1,'HandleVisibility','off')
text(1.05*k_cutoff,10^-1.5,['1/(' num2str(1/k_cutoff) ' km)'],'Interpreter','latex','FontSize',12)
LEG2 = legend('$\langle S_{vv}(k)\rangle$',...
             '$\langle S_{vv}(l)\rangle$',...
             '$\langle S_{v_\mathrm{g}v_\mathrm{g}}(l^\star)\rangle$',...
             '$\langle S_{v''_\mathrm{g}v''_\mathrm{g}}(l)\rangle$');
LEG2.Interpreter = 'latex';
LEG2.FontSize = 16;
LEG2.Position = [0.4833    0.1874    0.1686    0.2871]; % [0.4833    0.1885    0.1686    0.2159]; % [0.5374    0.1850    0.1852    0.2159];
xlabel('$k,l,l^\star$ (cpkm)','Interpreter','latex')
% ylabel('PSD (m$^2$ s$^{-2}$ cpkm$^{-1}$)','Interpreter','latex')
AX2.FontSize = 12;
AX2.XLim = [0.9*min([k_vec(1) l_vec(1) l_vec_swot(1)]), 1.1*max([k_vec(end) l_vec(end) l_vec_swot(end)])];
AX2.YLim = YLIM;
AX2.YTickLabel = {};
grid on

set(gcf,'Position',[-1181         521         702         290])

AX1.Position = AX1.Position - [ 0.001 -0.05 0 0 ];
AX2.Position = AX2.Position - [ 0.1   -0.05 0 0 ];

% % % % % % % % % % % % % % % % % % % % % % % ONE PANEL:

% figure('Color','w')
% loglog(k_vec,Suu_k_mean,'.-',         'LineWidth',1,'MarkerSize',10,'Color',[0 114 178  ]/255); hold on
% loglog(k_vec,Svv_k_mean,'.-',         'LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0 ]/255);
% loglog(l_vec,Suu_l_mean,'.-',         'LineWidth',1,'MarkerSize',10,'Color',[86 180 233 ]/255); hold on
% loglog(l_vec,Svv_l_mean,'.-',         'LineWidth',1,'MarkerSize',10,'Color',[204 121 167]/255);
% loglog(l_vec_swot,Sugug_l_mean,'.--k','LineWidth',1,'MarkerSize',10,'Color',[0 0 0      ]/255);
% loglog(l_vec_swot,Svgvg_l_mean,'.--k','LineWidth',1,'MarkerSize',10,'Color',[230 159 0  ]/255);
% AX1 = gca;
%     AX1.FontSize = 12;
%     % AX1.XLim = [0.9*min([k_vec(2) l_vec(2)]), 1.1*max([k_vec(end) l_vec(end)])];
%     AX1.XLim = [0.9*min([k_vec(1) l_vec(1) l_vec_swot(1)]), 1.1*max([k_vec(end) l_vec(end) l_vec_swot(end)])];
%     % AX1.XLim = XLIM;
%     AX1.YLim = YLIM;
% loglog([1 1]*k_cutoff,YLIM,'k--')
% text(1.05*k_cutoff,10^-1.5,['1/(' num2str(1/k_cutoff) ' km)'],'Interpreter','latex','FontSize',12)
% LEG = legend('$\langle S_{uu}(k)\rangle$','$\langle S_{vv}(k)\rangle$',...
%              '$\langle S_{uu}(l)\rangle$','$\langle S_{vv}(l)\rangle$',...
%              '$\langle S_{u_gu_g}(l^\star)\rangle$','$\langle S_{v_gv_g}(l^\star)\rangle$');
% LEG.Interpreter = 'latex';
% LEG.FontSize = 16;
% LEG.Position = [0.1836    0.1584    0.3164    0.2723];
% xlabel('$k,l,l^\star$ (cpkm)','Interpreter','latex')
% ylabel('PSD (m$^2$ s$^{-2}$ cpkm$^{-1}$)','Interpreter','latex')
% 
% set(gcf,'Position',[-896   364   417   447])

% % % % % % % % % % % % % % % % % % % % % % % 

% figure('Color','w')
% tiledlayout(2,1,"TileSpacing","tight")
% 
% AX1 = nexttile;
% loglog(k_vec,Suu_k_mean,'.-','LineWidth',1,'MarkerSize',10,'Color',[0 114 178]/255); hold on
% loglog(k_vec,Svv_k_mean,'.-','LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0]/255);
% % loglog(k_vec,Svv_k_10_50_90([1 3],:),'-r','LineWidth',1);
% loglog([1 1]*k_cutoff,YLIM,'k--')
%     AX1.FontSize = 12;
%     % AX1.XLim = [0.9*min([k_vec(2) l_vec(2)]), 1.1*max([k_vec(end) l_vec(end)])];
%     AX1.XLim = [0.9*min([k_vec(1) l_vec(1) l_vec_swot(1)]), 1.1*max([k_vec(end) l_vec(end) l_vec_swot(end)])];
%     % AX1.XLim = XLIM;
%     AX1.YLim = YLIM;
% LEG = legend('$\langle S_{uu}(k)\rangle$','$\langle S_{vv}(k)\rangle$','');
% LEG.Interpreter = 'latex';
% xlabel('$k$ (km$^{-1}$)','Interpreter','latex')
% ylabel('PSD (m$^2$ s$^{-2}$ cpkm$^{-1}$)','Interpreter','latex')
% 
% AX2 = nexttile;
% loglog(l_vec,Suu_l_mean,'.-','LineWidth',         1,'MarkerSize',10,'Color',[0 114 178]/255); hold on
% loglog(l_vec,Svv_l_mean,'.-','LineWidth',         1,'MarkerSize',10,'Color',[213 94/2 0]/255);
% loglog(l_vec_swot,Sugug_l_mean,'.-k','LineWidth', 1,'MarkerSize',10,'Color',[86 180 233]/255);
% loglog(l_vec_swot,Svgvg_l_mean,'.-k','LineWidth', 1,'MarkerSize',10,'Color',[204 121 167]/255);
%     AX2.FontSize = 12;
%     % AX2.XLim = [0.9*min([k_vec(2) l_vec(2)]), 1.1*max([k_vec(end) l_vec(end)])];
%     AX2.XLim = [0.9*min([k_vec(1) l_vec(1) l_vec_swot(1)]), 1.1*max([k_vec(end) l_vec(end) l_vec_swot(end)])];
%     % AX2.XLim = XLIM;
%     AX2.YLim = YLIM;
% loglog([1 1]*k_cutoff,YLIM,'k--')
% text(1.05*k_cutoff,10^-1.5,['1/(' num2str(1/k_cutoff) ' km)'],'Interpreter','latex','FontSize',12)
% LEG = legend('$\langle S_{uu}(l)\rangle$','$\langle S_{vv}(l)\rangle$','');
% LEG.Interpreter = 'latex';
% xlabel('$l$ (km$^{-1}$)','Interpreter','latex')
% ylabel('PSD (m$^2$ s$^{-2}$ cpkm$^{-1}$)','Interpreter','latex')
% 
% set(gcf,'Position',[757   173   344   503])

% %%
save_boolean = input('Do you want to save this figure? 1 (yes) or 0 (no)     ');
if save_boolean
    figure(1)
    exportgraphics(gcf,...
    '../figures/draft/F_vel_wavenumberspectra_HFR_SWOT.pdf',...
    'BackgroundColor','none','ContentType','vector')
    disp('Image saved')

else
end

%% PLOT: Publication-worthy plot using median spectra, not mean

k_cutoff = 1/18; % wavenumber of vertical dashed line
XLIM = 1./[100 10];
YLIM = [2*10^-3 10.^0];
    YLIM = [10^-6 10.^0];

% close all
% % % % % % % % % % % % % % % % % % % % % % % TWO-PANEL:

figure('Color','w')
% tiledlayout(1,2,"TileSpacing","tight")

% AX1 = nexttile;
AX1 = subplot(1,2,1);
loglog(k_vec,Suu_k_10_50_90(2,:),'.-',         'LineWidth',1,'MarkerSize',10,'Color',[0 114 178  ]/255); hold on
loglog(l_vec,Suu_l_10_50_90(2,:),'.-',         'LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0 ]/255);
    loglog(k_vec_rotfilt,Suu_rotfilt_k_10_50_90(2,:),'.--', 'LineWidth',1,'MarkerSize',10,'Color',[0 114 178  ]/255,'HandleVisibility','off'); hold on
    loglog(l_vec_rotfilt,Suu_rotfilt_l_10_50_90(2,:),'.--', 'LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0 ]/255,'HandleVisibility','off'); hold on
    loglog(k_vec_rotfilt,Suu_rf_k_10_50_90(2,:),'.:', 'LineWidth',1,'MarkerSize',10,'Color',[0 114 178  ]/255,'HandleVisibility','off'); hold on
    loglog(l_vec_rotfilt,Suu_rf_l_10_50_90(2,:),'.:', 'LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0 ]/255,'HandleVisibility','off'); hold on
loglog(l_vec_swot,Sugug_l_10_50_90(2,:),'.-',  'LineWidth',1,'MarkerSize',10,'Color',[0 0 0      ]/255);
    loglog(l_vec_swot_gridded,Sugug_gridded_l_10_50_90(2,:),'.-',  'LineWidth',1,'MarkerSize',10,'Color',[150 150 150]/255);
loglog([1 1]*k_cutoff,YLIM,'k--')
loglog(1.00*[1 1]*10^-1,err_HFR_k *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.20*[1 1]*10^-1,err_HFR_l *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.44*[1 1]*10^-1,err_SWOT_l*10^-3,'k','LineWidth',1,'HandleVisibility','off')
text(1.05*k_cutoff,10^-1.5,['1/(' num2str(1/k_cutoff) ' km)'],'Interpreter','latex','FontSize',12)
LEG1 = legend('median $S_{uu}(k)$',...
             'median $S_{uu}(l)$',...
             'median $S_{u_\mathrm{g}u_\mathrm{g}}(l^\star)$');
LEG1.Interpreter = 'latex';
LEG1.FontSize = 16;
LEG1.Position = [0.1433    0.1885    0.1686    0.2159]; % [0.1494    0.1850    0.1852    0.2159];
xlabel('$k,l,l^\star$ (cpkm)','Interpreter','latex')
ylabel('PSD (m$^2$ s$^{-2}$ cpkm$^{-1}$)','Interpreter','latex')
AX1.FontSize = 12;
AX1.XLim = [0.9*min([k_vec(1) l_vec(1) l_vec_swot(1)]), 1.1*max([k_vec(end) l_vec(end) l_vec_swot(end)])];
AX1.YLim = YLIM;
grid on

% AX2 = nexttile;
AX2 = subplot(1,2,2);
loglog(k_vec,Svv_k_10_50_90(2,:),'.-',         'LineWidth',1,'MarkerSize',10,'Color',[0 114 178  ]/255); hold on
loglog(l_vec,Svv_l_10_50_90(2,:),'.-',         'LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0 ]/255);
    loglog(k_vec_rotfilt,Svv_rotfilt_k_10_50_90(2,:),'.--', 'LineWidth',1,'MarkerSize',10,'Color',[0 114 178  ]/255,'HandleVisibility','off'); hold on
    loglog(l_vec_rotfilt,Svv_rotfilt_l_10_50_90(2,:),'.--', 'LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0 ]/255,'HandleVisibility','off'); hold on
    loglog(k_vec_rotfilt,Svv_rf_k_10_50_90(2,:),'.:', 'LineWidth',1,'MarkerSize',10,'Color',[0 114 178  ]/255,'HandleVisibility','off'); hold on
    loglog(l_vec_rotfilt,Svv_rf_l_10_50_90(2,:),'.:', 'LineWidth',1,'MarkerSize',10,'Color',[213 94/2 0 ]/255,'HandleVisibility','off'); hold on
loglog(l_vec_swot,Svgvg_l_10_50_90(2,:),'.-',  'LineWidth',1,'MarkerSize',10,'Color',[0 0 0      ]/255);
    loglog(l_vec_swot_gridded,Svgvg_gridded_l_10_50_90(2,:),'.:',  'LineWidth',1,'MarkerSize',10,'Color',[150 150 150]/255);
loglog([1 1]*k_cutoff,YLIM,'k--')
loglog(1.00*[1 1]*10^-1,err_HFR_k *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.20*[1 1]*10^-1,err_HFR_l *10^-3,'k','LineWidth',1,'HandleVisibility','off')
loglog(1.44*[1 1]*10^-1,err_SWOT_l*10^-3,'k','LineWidth',1,'HandleVisibility','off')
text(1.05*k_cutoff,10^-1.5,['1/(' num2str(1/k_cutoff) ' km)'],'Interpreter','latex','FontSize',12)
LEG2 = legend('median $S_{vv}(k)$',...
             'median $S_{vv}(l)$',...
             'median $S_{v_\mathrm{g}v_\mathrm{g}}(l^\star)$');
LEG2.Interpreter = 'latex';
LEG2.FontSize = 16;
LEG2.Position = [0.4833    0.1885    0.1686    0.2159]; % [0.5374    0.1850    0.1852    0.2159];
xlabel('$k,l,l^\star$ (cpkm)','Interpreter','latex')
% ylabel('PSD (m$^2$ s$^{-2}$ cpkm$^{-1}$)','Interpreter','latex')
AX2.FontSize = 12;
AX2.XLim = [0.9*min([k_vec(1) l_vec(1) l_vec_swot(1)]), 1.1*max([k_vec(end) l_vec(end) l_vec_swot(end)])];
AX2.YLim = YLIM;
AX2.YTickLabel = {};
grid on

set(gcf,'Position',[-1181         521         702         290])

AX1.Position = AX1.Position - [ 0.001 -0.05 0 0 ];
AX2.Position = AX2.Position - [ 0.1   -0.05 0 0 ];

%%
%%
%% Comparison of wavenumber spectra across different SWOT fitting window sizes:

% CG_Iteration                      1x1                       8  double                       
% ConvWithWindow                    1x1                      32  function_handle              
% PIXELS                            1x1                       8  double                       
% SWOT_time                        89x1                     712  double                       
% U_cyclogeostr_Nit               285x69x89            14001480  double                       
% U_geostr                        285x69x89            14001480  double                       
% U_rot_filtered_SWOTtimes         94x56x89             3747968  double                       
% U_rot_unfiltered_SWOTtimes       94x56x89             3747968  double                       
% Ucg_SWOT_HFRgrid_all             94x56x89             3747968  double                       
% Ug_SWOT_HFRgrid_all              94x56x89             3747968  double                       
% V_cyclogeostr_Nit               285x69x89            14001480  double                       
% V_geostr                        285x69x89            14001480  double                       
% V_rot_filtered_SWOTtimes         94x56x89             3747968  double                       
% V_rot_unfiltered_SWOTtimes       94x56x89             3747968  double                       
% Vcg_SWOT_HFRgrid_all             94x56x89             3747968  double                       
% Vg_SWOT_HFRgrid_all              94x56x89             3747968  double                       
% WINDOW                           36x1                     288  double

% SWOT_VEL_FILES = dir('./SWOT_and_HFR_velocities_*pix.mat');
SWOT_VEL_FILES = dir('./SWOT_and_HFR_velocities_*pix.mat');
SWOT_cell = cell(5,length(SWOT_VEL_FILES));

SWOTw_VEL_FILES = dir('./SWOT_and_HFR_velocities_*pix_weighted.mat');
SWOTw_cell = cell(5,length(SWOT_VEL_FILES));

for hh = 1:length(SWOT_VEL_FILES)
    SWOT_struct = load(SWOT_VEL_FILES(hh).name,'U_geostr','V_geostr');
    SWOTw_struct = load(SWOTw_VEL_FILES(hh).name,'U_geostr','V_geostr');
    if hh == 2
        MAX_BAD = 4;
    else
        MAX_BAD = 2;
    end

    SPEC_mat_ug_l_ = [];
    SPEC_mat_vg_l_ = [];
    wSPEC_mat_ug_l_ = [];
    wSPEC_mat_vg_l_ = [];
    XRANGE_ = [12:27,42:55];
    YRANGE_ = 137:277;
    dy_swot = 2;
    
    for ii = 1:size(SWOT_struct.U_geostr,3)
        for jj = XRANGE_
            SWOTvel_Data_u_     = SWOT_struct.U_geostr(YRANGE_,jj,ii);
            SWOTvel_Data_v_     = SWOT_struct.V_geostr(YRANGE_,jj,ii);
            if sum(~isfinite(SWOTvel_Data_v_)) <= MAX_BAD %  >= 0.9
                % SWOTvel_Data_u_ = SWOTvel_Data_u_;
                % SWOTvel_Data_v_ = SWOTvel_Data_v_;
                [SPEC_u_, l_vec_swot, ~] = nanspectrum(SWOTvel_Data_u_(1:end), dy_swot, 'km', 3, '.-',false,0,'hanning');
                [SPEC_v_, ~,          ~] = nanspectrum(SWOTvel_Data_v_(1:end), dy_swot, 'km', 3, '.-',false,0,'hanning');

                SPEC_mat_ug_l_ = [SPEC_mat_ug_l_ , SPEC_u_];
                SPEC_mat_vg_l_ = [SPEC_mat_vg_l_ , SPEC_v_];
            else
            end
        end

        if ~mod(ii,11)
            disp([hh round(100*ii/size(SWOT_struct.U_geostr,3))])
        else
        end

        for jj = XRANGE_
            SWOTwvel_Data_u_     = SWOTw_struct.U_geostr(YRANGE_,jj,ii);
            SWOTwvel_Data_v_     = SWOTw_struct.V_geostr(YRANGE_,jj,ii);
            if sum(~isfinite(SWOTwvel_Data_v_)) <= MAX_BAD %  >= 0.9
                SWOTwvel_Data_u_ = SWOTwvel_Data_u_;
                SWOTwvel_Data_v_ = SWOTwvel_Data_v_;
                [SPEC_u_, l_vec_swot, ~] = nanspectrum(SWOTwvel_Data_u_(1:end), dy_swot, 'km', 3, '.-',false,0,'hanning');
                [SPEC_v_, ~,          ~] = nanspectrum(SWOTwvel_Data_v_(1:end), dy_swot, 'km', 3, '.-',false,0,'hanning');

                wSPEC_mat_ug_l_ = [wSPEC_mat_ug_l_ , SPEC_u_];
                wSPEC_mat_vg_l_ = [wSPEC_mat_vg_l_ , SPEC_v_];
            else
            end
        end

        if ~mod(ii,11)
            disp([hh round(100*ii/size(SWOTw_struct.U_geostr,3))])
        else
        end
    end

    Sugug_l_10_50_90_ = prctile(SPEC_mat_ug_l_',[10 50 90]);
    Svgvg_l_10_50_90_ = prctile(SPEC_mat_vg_l_',[10 50 90]);
    Sugug_l_mean_ = mean(SPEC_mat_ug_l_,2,'omitnan');
    Svgvg_l_mean_ = mean(SPEC_mat_vg_l_,2,'omitnan');

    SWOT_cell{1,hh} = SWOT_VEL_FILES(hh).name;
    SWOT_cell{2,hh} = Sugug_l_10_50_90_;
    SWOT_cell{3,hh} = Svgvg_l_10_50_90_;
    SWOT_cell{4,hh} = Sugug_l_mean_;
    SWOT_cell{5,hh} = Svgvg_l_mean_;
    
    SWOT_cell{6,hh} = SWOT_struct.U_geostr(:,:,33);
    SWOT_cell{7,hh} = SWOT_struct.V_geostr(:,:,33);

    % % % % % % % 

    wSugug_l_10_50_90 = prctile(wSPEC_mat_ug_l_',[10 50 90]);
    wSvgvg_l_10_50_90 = prctile(wSPEC_mat_vg_l_',[10 50 90]);
    wSugug_l_mean = mean(wSPEC_mat_ug_l_,2,'omitnan');
    wSvgvg_l_mean = mean(wSPEC_mat_vg_l_,2,'omitnan');

    SWOTw_cell{1,hh} = SWOTw_VEL_FILES(hh).name;
    SWOTw_cell{2,hh} = wSugug_l_10_50_90;
    SWOTw_cell{3,hh} = wSvgvg_l_10_50_90;
    SWOTw_cell{4,hh} = wSugug_l_mean;
    SWOTw_cell{5,hh} = wSvgvg_l_mean;
    
    SWOTw_cell{6,hh} = SWOTw_struct.U_geostr(:,:,33);
    SWOTw_cell{7,hh} = SWOTw_struct.V_geostr(:,:,33);
end

%% Plot different SWOT wavenumber spectra:

close all

k_cutoff = 1/18;
YLIM = 10.^[-7 0];

figure
subplot(1,2,1)
loglog(l_vec_swot,SWOT_cell{5,1},'.-','LineWidth',2,'Markersize',20); hold on
loglog(l_vec_swot,SWOT_cell{5,2},'.-','LineWidth',2,'Markersize',20);
loglog(l_vec_swot,SWOT_cell{5,3},'.-','LineWidth',2,'Markersize',20);
loglog(l_vec_swot,SWOT_cell{5,4},'.-','LineWidth',2,'Markersize',20);
loglog(l_vec_swot,SWOT_cell{5,5},'.-','LineWidth',2,'Markersize',20);
loglog(l_vec_swot,SWOT_cell{5,6},'.-','LineWidth',2,'Markersize',20);
loglog([1 1]*k_cutoff,YLIM,'k--')
LEG_in = {};
for ii = 1:size(SWOT_cell,2)
    LEG_in{ii} = replace(SWOT_cell{1,ii},{'SWOT_and_HFR_velocities_','.mat','_w'},{'','','\_w'});
end
legend(LEG_in)
ylim(10.^[-6 0])
% % % % % % % 
subplot(1,2,2)
loglog(l_vec_swot,SWOTw_cell{5,1},'.-','LineWidth',2,'Markersize',20); hold on
loglog(l_vec_swot,SWOTw_cell{5,2},'.-','LineWidth',2,'Markersize',20);
loglog(l_vec_swot,SWOTw_cell{5,3},'.-','LineWidth',2,'Markersize',20);
loglog(l_vec_swot,SWOTw_cell{5,4},'.-','LineWidth',2,'Markersize',20);
loglog(l_vec_swot,SWOTw_cell{5,5},'.-','LineWidth',2,'Markersize',20);
loglog(l_vec_swot,SWOTw_cell{5,6},'.-','LineWidth',2,'Markersize',20);
loglog([1 1]*k_cutoff,YLIM,'k--')
LEG_in = {};
for ii = 1:size(SWOTw_cell,2)
    LEG_in{ii} = replace(SWOTw_cell{1,ii},{'SWOT_and_HFR_velocities_','.mat','_w'},{'','','\_w'});
end
legend(LEG_in)
ylim(10.^[-6 0])


figure
tiledlayout(2,3,"TileSpacing","tight")
for ii = 1:size(SWOT_cell,2)
    II = dsearchn([2:2:12]',...
                  str2num(replace(SWOT_cell{1,ii},...
                                  {'SWOT_and_HFR_velocities_','pix_weighted.mat','pix.mat'},...
                                  {'','',''}))...
                  );
    nexttile(II)
    loglog(l_vec_swot,SWOT_cell{     4,ii},'.-','LineWidth',2,'Markersize',20); hold on
    LEG_in_{1} = replace(SWOT_cell{  1,ii},{'SWOT_and_HFR_velocities_','.mat','_w'},{'','','\_w'});
    loglog(l_vec_swot,SWOTw_cell{    4,ii},'.-','LineWidth',2,'Markersize',20); hold on
    LEG_in_{2} = replace(SWOTw_cell{ 1,ii},{'SWOT_and_HFR_velocities_','.mat','_w'},{'','','\_w'});

    % loglog(k_vec_rotfilt,Suu_rf_k_mean,'.-', 'LineWidth',1,'MarkerSize',10,'HandleVisibility','on'); hold on
    loglog(l_vec_rotfilt,Suu_rf_l_mean,'.-', 'LineWidth',1,'MarkerSize',10,'HandleVisibility','on'); hold on
    % loglog(k_vec_rotfilt,Svv_rf_k_mean,'.-', 'LineWidth',1,'MarkerSize',10,'HandleVisibility','on'); hold on
    loglog(l_vec_rotfilt,Svv_rf_l_mean,'.-', 'LineWidth',1,'MarkerSize',10,'HandleVisibility','on'); hold on

    LEG_in_{3} = '<S(l)> of rot. l.p. HFR u';
    LEG_in_{4} = '<S(l)> of rot. l.p. HFR v';

    loglog([1 1]*k_cutoff,YLIM,'k--')
    LEG = legend(LEG_in_,'Location','southwest');
end

% %%

figure
tiledlayout(1,size(SWOT_cell,2),"TileSpacing","tight")
Change_unanalyzed_pixels = 0*ones(size(SWOT_cell{7,1}));
Change_unanalyzed_pixels(YRANGE_,XRANGE_) = 1;
for ii = 1:size(SWOT_cell,2)
    nexttile
    imagesc(SWOT_cell{7,ii} ./ Change_unanalyzed_pixels)
    axis equal
    title(LEG_in{ii})
    clim([-1 1])
end
xlabel('Smoothed with no weighting (Original)')

figure
tiledlayout(1,size(SWOTw_cell,2),"TileSpacing","tight")
Change_unanalyzed_pixels = 0*ones(size(SWOTw_cell{7,1}));
Change_unanalyzed_pixels(YRANGE_,XRANGE_) = 1;
for ii = 1:size(SWOT_cell,2)
    nexttile
    imagesc(SWOTw_cell{7,ii} ./ Change_unanalyzed_pixels)
    axis equal
    title(LEG_in{ii})
    clim([-1 1])
end
xlabel('Smoothed with edge weighting')


figure
Change_unanalyzed_pixels = nan(size(SWOT_cell{7,1}));
Change_unanalyzed_pixels(YRANGE_,XRANGE_) = 1;
for ii = 1:size(SWOT_cell,2)
    histogram(SWOT_cell{7,ii} .* Change_unanalyzed_pixels,'Normalization','pdf'); hold on
end
legend(LEG_in)

figure
for ii = 1:size(SWOT_cell,2)
    plot(SWOT_cell{7,ii}(YRANGE_,20),'.-'); hold on
end
legend(LEG_in)


%% Supplemental plot to show different spectral levels of HFR vs. SWOT

% SWOTw_cell{1,hh} = SWOTw_VEL_FILES(hh).name;
% SWOTw_cell{2,hh} = wSugug_l_10_50_90;
% SWOTw_cell{3,hh} = wSvgvg_l_10_50_90;
% SWOTw_cell{4,hh} = wSugug_l_mean;
% SWOTw_cell{5,hh} = wSvgvg_l_mean;

close all

clear LEG_in_
YLIM = 10.^[-8 0.1];

figure('Color','w')
tiledlayout(2,3,"TileSpacing","tight")
for ii = 1:size(SWOT_cell,2)
    II = dsearchn([2:2:12]',...
                  str2num(replace(SWOT_cell{1,ii},...
                                  {'SWOT_and_HFR_velocities_','pix_weighted.mat','pix.mat'},...
                                  {'','',''}))...
                  );
    nexttile(II)
    loglog(l_vec_swot,SWOT_cell{     4,ii},'.-','LineWidth',1,'Markersize',10); hold on
    LEG_in_{1} = [num2str(1 + 2*str2double(replace(SWOT_cell{  1,ii},{'SWOT_and_HFR_velocities_','pix.mat'},{'',''}))) ' pixels'];
    loglog(l_vec_swot,SWOTw_cell{    4,ii},'.-','LineWidth',1,'Markersize',10); hold on
    LEG_in_{2} = [num2str(1 + 2*str2double(replace(SWOTw_cell{ 1,ii},{'SWOT_and_HFR_velocities_','pix_weighted.mat'},{'',''}))) ' pixels, weighted'];

    loglog(l_vec_rotfilt,Suu_rf_l_mean,'.-k', 'LineWidth',1,'MarkerSize',10,'HandleVisibility','on'); hold on

    LEG_in_{3} = 'Rotational Low-pass HFR';

    loglog([1 1]*k_cutoff,YLIM,'k--','HandleVisibility','off')
    LEG = legend(LEG_in_,'Location','southwest','BackgroundAlpha',0.5,'Interpreter','latex');

    if II == 1
        ylabel({'Power Spectral Density';'($m^{2} s^{-2}$ cpkm$^{-1}$)'},'Interpreter','latex')
        title('$S_{uu}(l)$, $S_{uu}(l^\star)$','Interpreter','latex')
    elseif II == 4
        xlabel({'Wavenumber (cpkm)'},'Interpreter','latex')
    else
    end

    ylim(YLIM)
    set(gca,'FontSize',14)
end
% set(gcf,'Position',[154   119   812   574])
% set(gcf,'Position',[154   119   857   574])
set(gcf,'Position',[154   119   903   574])

figure('Color','w')
tiledlayout(2,3,"TileSpacing","tight")
for ii = 1:size(SWOT_cell,2)
    II = dsearchn([2:2:12]',...
                  str2num(replace(SWOT_cell{1,ii},...
                                  {'SWOT_and_HFR_velocities_','pix_weighted.mat','pix.mat'},...
                                  {'','',''}))...
                  );
    nexttile(II)
    loglog(l_vec_swot,SWOT_cell{     5,ii},'.-','LineWidth',1,'Markersize',10); hold on
    LEG_in_{1} = [num2str(1 + 2*str2double(replace(SWOT_cell{  1,ii},{'SWOT_and_HFR_velocities_','pix.mat'},{'',''}))) ' pixels'];
    loglog(l_vec_swot,SWOTw_cell{    5,ii},'.-','LineWidth',1,'Markersize',10); hold on
    LEG_in_{2} = [num2str(1 + 2*str2double(replace(SWOTw_cell{ 1,ii},{'SWOT_and_HFR_velocities_','pix_weighted.mat'},{'',''}))) ' pixels, weighted'];

    loglog(l_vec_rotfilt,Svv_rf_l_mean,'.-k', 'LineWidth',1,'MarkerSize',10,'HandleVisibility','on'); hold on

    LEG_in_{3} = 'Rotational Low-pass HFR';

    loglog([1 1]*k_cutoff,YLIM,'k--','HandleVisibility','off')
    LEG = legend(LEG_in_,'Location','southwest','BackgroundAlpha',0.5,'Interpreter','latex');

    if II == 1
        ylabel({'Power Spectral Density';'($m^{2} s^{-2}$ cpkm$^{-1}$)'},'Interpreter','latex')
        title('$S_{vv}(l)$, $S_{vv}(l^\star)$','Interpreter','latex')
    elseif II == 4
        xlabel({'Wavenumber (cpkm)'},'Interpreter','latex')
    else
    end

    ylim(YLIM)
    set(gca,'FontSize',14)
end
% set(gcf,'Position',[154   119   812   574])
% set(gcf,'Position',[154   119   857   574])
set(gcf,'Position',[200   119   903   574])

%%
figure(2)
exportgraphics(gcf,...
'../figures/draft/F_vel_wavenumberspectra_HFR_SWOT_6panel_v.pdf',...
'BackgroundColor','none','ContentType','vector')
disp('Image saved')


%%
%%
%%
%% Auxiliary functions

function Input_trimmed = trim_nans(Input)
if sum(isfinite(Input)) == 0
    Input_trimmed = [];
else
Indices        = 1:length(Input);
WhereIsFinite  = isfinite(Input);
IndicesFinite  = Indices(WhereIsFinite);
i_first_finite = IndicesFinite(1);
i_last_finite  = IndicesFinite(end);
Input_trimmed = Input(i_first_finite:i_last_finite);
end
end

%% cov_gaps
% [Cov, Lags, Cov_std, n_Cov] = cov_gaps(A,B,maxlags)
% 
% Gives the true cross covariance between two equally sized vectors. This
% does not use any tricks to make the calculation fast, but rather it uses
% the definition of covariance:
% 
% C_ab(dt) = <A(t)B(t + dt)>
% 
% Where both A and B are de-meaned vectors.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% IN:   A =         Nx1 vector (may be 1xN, but outputs will always be columns).
%                   This may have NaN or Inf for undefined values, but what is key
%                   is that the data (whether finite or not) be even spaced.
% IN:   B =         Vector of the same size as A. If A = B, then this calculates
%                   autocovariance.
% IN:   maxlags =   (Optional, scalar) Maximum lag considered for calculating Cov.
%                   The default value is length(A) - 1.
%                   Note that Lags(end) = maxlags.
% 
% OUT:  Cov =       maxlags x 1 vector, covariance between A and B. See the
%                   equation for C_ab above.
% OUT:  Lags =      (Optional, maxlags x 1) The corresponding lags of the
%                   entries in Cov.
% OUT:  Cov_std =   (Optional, maxlags x 1) The STD of the values that went
%                   into calculating each entry in Cov.
% OUT:  n_Cov =     (Optional, maxlags x 2) The first column is the number of
%                   pairs A(t)B(t + dt) for each dt which were finite (i.e.
%                   neither A(t) nor B(t + dt) were NaN or Inf). The second
%                   column is the number of pairs that would have been used if
%                   there were no NaN's or Inf's.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [Cov, varargout] = cov_gaps(A,B,varargin)

if size(A) ~= size(B)
    error('A and B must be the same size.')
else
end

if isrow(A)
    A = A'; B = B';
else
end

N = length(A);

A = A - nanmean(A);
B = B - nanmean(B);

if nargin == 2
    maxlags = N - 1;
elseif nargin == 3
    maxlags = varargin{1};
else
    error('Incorrect number of inputs; please see documentation.')
end

Cov = zeros(1 + maxlags,1); % initialize the crosscovariance
Lags = zeros(1 + maxlags,1); % initialize the corresponding lags
Cov_std = zeros(1 + maxlags,1); % the std of the points that went into each element in pre-shaped "Cov" (starting with zero-lag)
n_Cov = zeros(1 + maxlags,2); % how many pairs were finite and how many pairs could have been finite if perfect
for i=1:(1 + maxlags) % lag = (i-1)*dt
    Crosscov_ABi = zeros(N - i + 1,1);
    Crosscov_BAi = zeros(N - i + 1,1);
    for j=1:(N - i + 1)
        Crosscov_ABi(j) = A(j)*B(j+i-1);
        Crosscov_BAi(j) = B(j)*A(j+i-1);
    end
    Crosscov_i = [Crosscov_ABi ; Crosscov_BAi];
    Lags(i) = i - 1;
    Cov(i) = nanmean(Crosscov_i); % autocov
    Cov_std(i) = nanstd(Crosscov_i); % std of points that made each element of autocov
    n_Cov(i,1) = sum(isfinite(Crosscov_i)); % how many pairs were finite
    n_Cov(i,2) = 2*(N - i + 1); % how many pairs could have been finite if perfect
end

% % Use "Cov_unfolded" (others are for completeness/clarity) for
% % estimating (cross)spectra by taking the fft.
% Cov_unfolded = [Cov;flip(Cov(2:(end-1)))];
% Cov_std_unfolded = [Cov_std;flip(Cov_std(2:(end-1)))];
% Lags_unfolded = [Lags;flip(Lags(2:(end-1)))];

if nargout == 1
elseif nargout == 2
    varargout{1} = Lags;
elseif nargout == 3
    varargout{1} = Lags;
    varargout{2} = Cov_std;
elseif nargout == 4
    varargout{1} = Lags;
    varargout{2} = Cov_std;
    varargout{3} = n_Cov;
elseif nargout == 0
else
    error('The number of outputs must be 1, 2, 3, or 4.')
end

% % This appears to taper just the right amount, at least for white noise:
% foo = randn(1000,1);
% [Cov, Lags, Cov_std, n_Cov] = cov_gaps(foo,foo);
% plot(Lags,Cov.*(n_Cov(:,1)/n_Cov(1,1)).^.5,'.:')

end


%% Discarded
% close all
% 
% % % % Establish the set of wavenumbers to fit to:
% % k_vec = [0:[1/[20*dx]]:[1/[2*dx]]]';
% % l_vec = [0:[1/[30*dy]]:[1/[2*dy]]]';
% 
% figure
% SPEC_mat_urotfilt_l = [];
% SPEC_mat_vrotfilt_l = [];
% XRANGE = 1:17;
% for ii = 1:size(U_rot_filtered_SWOTtimes,3)
%     for jj = XRANGE
%         HFR_Data_urotfilt = U_rot_filtered_SWOTtimes(:,jj,ii);
%         HFR_Data_vrotfilt = V_rot_filtered_SWOTtimes(:,jj,ii);
%         HFR_Data_urotfilt = trim_nans(HFR_Data_urotfilt);
%         HFR_Data_vrotfilt = trim_nans(HFR_Data_vrotfilt);
% 
%         if [length(HFR_Data_urotfilt) >= 81] && [length(HFR_Data_vrotfilt) >= 81] && sum(~isfinite(HFR_Data_urotfilt)) == 0
%             [SPEC_urotfilt, l_vec_rotfilt, ~] = nanspectrum(HFR_Data_urotfilt(1:81), dy, 'km', 3, '.-',false,0,'hanning');
%             [SPEC_vrotfilt,              ~, ~] = nanspectrum(HFR_Data_vrotfilt(1:81), dy, 'km', 3, '.-',false,0,'hanning');
%             SPEC_mat_urotfilt_l = [SPEC_mat_urotfilt_l , SPEC_urotfilt];
%             SPEC_mat_vrotfilt_l = [SPEC_mat_vrotfilt_l , SPEC_vrotfilt];
%         else
%         end
%     end
% 
%     disp(100*ii/size(U_rot_filtered_SWOTtimes,3))
% 
% end
% % %%
% figure
% subplot(121)
% loglog(l_vec_rotfilt,      SPEC_mat_urotfilt_l(:,1:10:end),'-');hold on
% loglog(l_vec_rotfilt, mean(SPEC_mat_urotfilt_l,2),'k*-','LineWidth',2)
% loglog(l_vec_rotfilt,[mean(SPEC_mat_urotfilt_l(1,:))/l_vec_rotfilt(2)^-2]*[l_vec_rotfilt.^-2],'g--','LineWidth',2)
% xlabel('Zonal wavenumber l (km^-1)')
% ylabel('PSD ([m/2]^4 [km^-1]^-1)')
% subplot(122)
% loglog(l_vec_rotfilt,prctile(SPEC_mat_urotfilt_l',[50]),'ko-','LineWidth',2); hold on
% loglog(l_vec_rotfilt,prctile(SPEC_mat_urotfilt_l',[10 90]),'k--','LineWidth',1)
% 
% Suu_rotfilt_l_10_50_90 = prctile(SPEC_mat_urotfilt_l,[10 50 90]); Suu_rotfilt_l_mean = mean(SPEC_mat_urotfilt_l,2);
% Svv_rotfilt_l_10_50_90 = prctile(SPEC_mat_vrotfilt_l,[10 50 90]); Svv_rotfilt_l_mean = mean(SPEC_mat_vrotfilt_l,2);
% 
% % Map of where you are looking:
% figure
% Analysis_Block = zeros(size(U_rot_filtered_SWOTtimes(:,:,1)));
% Analysis_Block(:,XRANGE) = 1;
% imagesc(sum(isfinite(U_rot_filtered_SWOTtimes),3)./size(U_rot_filtered_SWOTtimes,3) + Analysis_Block)
% set(gcf,'Position',[-1390 288 344 540])