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

% Load previously-calculated geostrophic velocity estimates:
VELS = load('./SWOT_and_HFR_velocities_8pix.mat');

close all

%% Visualize time coverage of SST

for ii = 1:NORCAL.SST.mean_datenum
    SST_COVERAGE = [];
end

figure
plot(NORCAL.SST.mean_datenum,SST_COVERAGE,'.-')

%%

% Manually verified:
SST_is_actually_good =[...
    20; ... %
    23; ...
    32; ...
    33; ... %
    35; ...
    38; ...
    40; ... %
    42; ...
    80; ...
    81; ... %
    82; ...
    86; ...
   110; ... %
   116; ...
   117; ...
   130; ... %
   131; ...
   132; ...
   134; ... %
   135; ...
   143];

% DATE                     =   INDEX
% 
% '26-Apr-2023 20:30:25'   =   80
% '26-Apr-2023 22:10:26'   =   81
% '27-Apr-2023 10:29:10'   =   82
% '27-Apr-2023 21:51:27'   =   86

% '09-May-2023 21:26:01'   =   20
% '11-May-2023 20:48:44'   =   23    } Not actually very good
% '28-May-2023 22:09:57'   =   110
% '30-May-2023 21:32:03'   =   116   } same data
% '30-May-2023 21:32:04'   =   32    }
% '31-May-2023 09:51:05'   =   117
% '31-May-2023 21:13:12'   =   33

% '01-Jun-2023 20:54:44'   =   35
% '02-Jun-2023 20:35:15'   =   38
% '19-Jun-2023 20:16:35'   =   130
% '19-Jun-2023 21:56:55'   =   131   } same data
% '19-Jun-2023 21:56:56'   =   40    }
% '20-Jun-2023 10:15:52'   =   132
% '20-Jun-2023 21:37:58'   =   134   } same data
% '20-Jun-2023 21:37:59'   =   42    }
% '21-Jun-2023 09:57:00'   =   135

% '01-Jul-2023 21:31:45'   =   143

% Copied from a CSV sheet from the IOOS, so it includes East Coast stations
% as well.
% <https://hfradar.ioos.us/hfrnet/sitediag/stationList.php>
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

%% % % % Select from SST indices
% Dates_List = NORCAL.SST.mean_datenum([80 81 82 86]); % April
% Dates_List = NORCAL.SST.mean_datenum([20 23]); % Early May
% Dates_List = NORCAL.SST.mean_datenum([110 116 117 33]); % Late May
% Dates_List = NORCAL.SST.mean_datenum([35 38]); % Early June
% Dates_List = NORCAL.SST.mean_datenum([130 131 132 134 135]); % Late June
% Dates_List = NORCAL.SST.mean_datenum([143]); % July (one day)
% 
Dates_List = NORCAL.SST.mean_datenum([20 116 35]); % 3x3 for publication
% Dates_List = NORCAL.SST.mean_datenum([20 35]); % 2x3 for publication if 3x3 is too much
% Dates_List = NORCAL.SST.mean_datenum([20]); % 1x3 for publication if 2x3 is too much
% 
[Dates_List,SortInd] = sort(Dates_List);
Dates_CharArray = datestr(Dates_List);
USER_PICKED_DATES = {};
for ii = 1:size(Dates_List,1)
    USER_PICKED_DATES{ii} = Dates_CharArray(ii,:);
end
% error

LatRange = [38 42];
LonRange = [min(NORCAL.HFR.lon(:)) -122.5];
Antenna_Locations = Antenna_Locations(real(Antenna_Locations) > min(LonRange) & ...
                                      real(Antenna_Locations) < max(LonRange) & ...
                                      imag(Antenna_Locations) > min(LatRange) & ...
                                      imag(Antenna_Locations) < max(LatRange));

FONTSIZE = 24;

% % % Manually chosen:
% USER_PICKED_DATES = {'2023-05-10 02:00:00' , ...
%                      '2023-05-19 02:00:00' , ...
%                      '2023-05-29 02:00:00'};
% USER_PICKED_DATES = {'2023-05-10 02:00:00'}; % One row to test

close all
figure('Color','w')
tiledlayout(length(USER_PICKED_DATES),3,'TileSpacing','tight')
Alphabet = 'abcdefghijklmnopqrstuvwxyz';
LetterLabels = Alphabet(1:[3*length(USER_PICKED_DATES)]);

for jj = 1:length(USER_PICKED_DATES)

    USER_PICKED_DATE = USER_PICKED_DATES{jj};

    ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum(USER_PICKED_DATE));
    SWOT_time = NORCAL.SWOT.mean_time(ti_swot);
    ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, SWOT_time);
    ti_SST = dsearchn(NORCAL.SST.mean_datenum,datenum(USER_PICKED_DATE));
    SST_t = NORCAL.SST.mean_datenum(ti_SST);
    
    eval(['AX' num2str(1 + 3*[jj-1]) ' = nexttile;']) % subplot(131) % SWOT SSHA
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',LonRange + [-.1 .1],...
        'latitudes',LatRange + [-.1 .1]);
    set(gcf,'color','w')
    COAST = m_gshhs_i('patch',0.8*[1 1 1]);
    m_grid('box','fancy', 'backgroundcolor','none'); hold on
    set(gcf,'Position',[262    69   962   728])
    m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, 100* ...
              [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
              [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
              [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]);
    CB = colorbar; CB.Label.String = '$\eta''$ (cm)'; CB.Label.Interpreter = 'latex';
    CB.FontSize = 12; CB.Label.FontSize = 16;
    clim([-1 1]*15)
    m_text(-123.4, 41.7, ['(' LetterLabels(1 + 3*[jj-1]) ')'], 'fontsize', FONTSIZE, 'Interpreter', 'latex')
    title(['SWOT time = ' datestr(SWOT_time)])
    
    eval(['AX' num2str(2 + 3*[jj-1]) ' = nexttile;']) % subplot(132) % HFR VORTICITY AND ARROWS
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',LonRange + [-.1 .1],...
        'latitudes',LatRange + [-.1 .1]);
    set(gcf,'color','w')
    COAST = m_gshhs_i('patch',0.8*[1 1 1]);
    m_grid('box','fancy', 'backgroundcolor','none'); hold on
    set(gcf,'Position',[262    69   962   728])
    % warning('Requires "vorticity2D" function')
    % m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, vorticity2D(NORCAL.HFR.u(:,:,ti_hfr), NORCAL.HFR.v(:,:,ti_hfr), 6000, 6000));
    m_pcolor_centered(NORCAL.HFR.LON, NORCAL.HFR.LAT, ...
                          [curl(NORCAL.HFR.LON.*cosd(NORCAL.HFR.LAT), NORCAL.HFR.LAT , NORCAL.HFR.u(:,:,ti_hfr), NORCAL.HFR.v(:,:,ti_hfr))*[1/111000]]./...
                          fcor_degrees_radps(NORCAL.HFR.LAT));
    M_PLOT = m_plot(real(Antenna_Locations),imag(Antenna_Locations),'k.','MarkerSize',25,'MarkerFaceColor','k');
    CB = colorbar; CB.Label.String = '$\zeta/f$'; CB.Label.Interpreter = 'latex';
    CB.FontSize = 12; CB.Label.FontSize = 16;
    clim([-1 1.0001]*0.3)
    m_text(-123.4, 41.7, ['(' LetterLabels(2 + 3*[jj-1]) ')'], 'fontsize', FONTSIZE, 'Interpreter', 'latex')
    m_quiver(NORCAL.HFR.LON,NORCAL.HFR.LAT,NORCAL.HFR.u(:,:,ti_hfr),NORCAL.HFR.v(:,:,ti_hfr),'Color','k');
    title(['HFR time = '  datestr(T0 + NORCAL.HFR.time(ti_hfr)/24)])
    
    eval(['AX' num2str(3 + 3*[jj-1]) ' = nexttile;']) % subplot(133) % VIIRS SST
    M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
        'longitudes',LonRange + [-.1 .1],...
        'latitudes',LatRange + [-.1 .1]);
    set(gcf,'color','w')
    COAST = m_gshhs_i('patch',0.8*[1 1 1]);
    m_grid('box','fancy', 'backgroundcolor','none'); hold on
    set(gcf,'Position',[262    69   962   728])
    m_pcolor_centered(NORCAL.SST.longitude{ti_SST}, NORCAL.SST.latitude{ti_SST}, ... flagged up to some number
                      NORCAL.SST.sst{ti_SST}.*[[NORCAL.SST.qual_sst{ti_SST}<3]./ ...
                                               [NORCAL.SST.qual_sst{ti_SST}<3]]);
    m_contour(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, 100* ...
              [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
              [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
              [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}],[-15:5:15],'k');
    CB = colorbar; CB.Label.String = 'SST ($^\circ$C)'; CB.Label.Interpreter = 'latex';
    CB.FontSize = 12; CB.Label.FontSize = 16;
    clim([9 14])
    m_text(-123.4, 41.7, ['(' LetterLabels(3 + 3*[jj-1]) ')'], 'fontsize', FONTSIZE, 'Interpreter', 'latex')
    title(['SST time = '  datestr(SST_t)])

    % eval(['colormap(AX' num2str(1 + 3*[jj-1]) ',cmocean(''balance''));']);
    eval(['colormap(AX' num2str(1 + 3*[jj-1]) ',''bwr'');']);
    eval(['colormap(AX' num2str(2 + 3*[jj-1]) ',flip(cmocean(''balance'')));']);
    eval(['colormap(AX' num2str(3 + 3*[jj-1]) ',''turbo'');']);

end

% linkaxes([AX1 AX2 AX3 AX4 AX5 AX6 AX7 AX8 AX9],'xy')
% linkaxes([AX1 AX2 AX3 AX4 AX5 AX6],'xy')
linkaxes([AX1 AX2 AX3],'xy')

% % % For LatRange = [38 41]:
% set(gcf,'Position',[-1512 79 926 898]) % 3x3
% set(gcf,'Position',[-1512 83 926 615]) % 2x3 and 1x3

% % % For LatRange = [38 42]:
set(gcf,'Position',[-1371 83 785 894]) % 3x3
% set(gcf,'Position',[-1609 83 1023 894]) % 2x3 and 1x3

%% Save image

% figure(1)
% exportgraphics(gcf,...
% ['../figures/draft/F1_' num2str(jj*3) 'panels.pdf'],...
% 'BackgroundColor','none','ContentType','vector')

%% Close-up of southern eddy with different vectors, SST background

LonRange = [-124.75 -123.5];
LatRange = [38.5 39.25];
% LonRange = [-130 -120];
% LatRange = [35 45];

close all
figure('Color','w')
tiledlayout(1,2,'TileSpacing','tight')

USER_PICKED_DATE = USER_PICKED_DATES{1};
ti_swot = dsearchn(NORCAL.SWOT.mean_time,datenum(USER_PICKED_DATE));
SWOT_time = NORCAL.SWOT.mean_time(ti_swot);
ti_hfr = dsearchn(T0 + NORCAL.HFR.time/24, SWOT_time);
ti_SST = dsearchn(NORCAL.SST.mean_datenum,datenum(USER_PICKED_DATE));
SST_t = NORCAL.SST.mean_datenum(ti_SST);

max_swot_speed = max(abs(Collimate(VELS.U_geostr(:,:,ti_swot) + 1i*VELS.V_geostr(:,:,ti_swot))));
max_hfr_speed  = max(abs(Collimate(VELS.U_rot_filtered_SWOTtimes(:,:,ti_swot) + 1i*VELS.V_rot_filtered_SWOTtimes(:,:,ti_swot))));

nexttile
M_P = m_proj('Miller',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
    'longitudes',LonRange + [-.1 .1],...
    'latitudes',LatRange + [-.1 .1]);
set(gcf,'color','w')
COAST = m_gshhs_i('patch',0.8*[1 1 1]); hold on
set(gcf,'Position',[262    69   962   728])
m_pcolor_centered(NORCAL.SST.longitude{ti_SST}, NORCAL.SST.latitude{ti_SST}, ... flagged up to some number
                  NORCAL.SST.sst{ti_SST}.*[[NORCAL.SST.qual_sst{ti_SST}<3]./ ...
                  [NORCAL.SST.qual_sst{ti_SST}<3]]);
m_quiver(NORCAL.SWOT.lon{ti_swot},NORCAL.SWOT.lat{ti_swot},...
         ...VELS.U_geostr(:,:,ti_swot),VELS.V_geostr(:,:,ti_swot),...
         VELS.U_cyclogeostr_Nit(:,:,ti_swot),VELS.V_cyclogeostr_Nit(:,:,ti_swot),...
         2,'Color','k')
m_grid('box','fancy', 'backgroundcolor','none');
CB = colorbar; CB.Label.String = 'SST ($^\circ$C)'; CB.Label.Interpreter = 'latex';
CB.FontSize = 12; CB.Label.FontSize = 16;
clim([11 13.5])
m_text(-123.4, 41.7, ['(' LetterLabels(3 + 3*[jj-1]) ')'], 'fontsize', FONTSIZE, 'Interpreter', 'latex')
title({['SST time = '  datestr(SST_t)] , ...
       ['SWOT time = ' datestr(SWOT_time)]})

nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
    'longitudes',LonRange + [-.1 .1],...
    'latitudes',LatRange + [-.1 .1]);
set(gcf,'color','w')
COAST = m_gshhs_i('patch',0.8*[1 1 1]); hold on
set(gcf,'Position',[262    69   962   728])
m_pcolor_centered(NORCAL.SST.longitude{ti_SST}, NORCAL.SST.latitude{ti_SST}, ... flagged up to some number
                  NORCAL.SST.sst{ti_SST}.*[[NORCAL.SST.qual_sst{ti_SST}<3]./ ...
                  [NORCAL.SST.qual_sst{ti_SST}<3]]);
m_quiver(NORCAL.HFR.LON,NORCAL.HFR.LAT,...
         VELS.U_rot_filtered_SWOTtimes(:,:,ti_swot), VELS.V_rot_filtered_SWOTtimes(:,:,ti_swot),...
         [2]*max_hfr_speed/max_swot_speed,'Color','k')
% ^ Note that [2] would need to be [2/3] to scale to be the same as the
% first panel, but that also leads to hard-to-see vectors here. Make sure
% to note the factor of ~3 in the manuscript.
m_grid('box','fancy', 'backgroundcolor','none');
CB = colorbar; CB.Label.String = 'SST ($^\circ$C)'; CB.Label.Interpreter = 'latex';
CB.FontSize = 12; CB.Label.FontSize = 16;
clim([11 13.5])
m_text(-123.4, 41.7, ['(' LetterLabels(3 + 3*[jj-1]) ')'], 'fontsize', FONTSIZE, 'Interpreter', 'latex')
title(['HFR time = '  datestr(T0 + NORCAL.HFR.time(ti_hfr)/24)])

colormap(turbo)

% %%
% figure(1)
% exportgraphics(gcf,...
% ['../figures/draft/F_vectors_over_SST.pdf'],...
% 'BackgroundColor','none','ContentType','vector')


%%
%%
%%
%% Auxiliary functions
%%
% This gradient function was written to account for rotated bases, like
% with SWOT, which are not accounted fo r in the default MATLAB gradient
% function:
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
%%
% Coriolis frequency in rad s^-1, given input of degrees latitude:
function FCOR = fcor_degrees_radps(LAT)
    % sideral day = 23.9344696 hours = 86164.1 seconds
    FCOR = 2*[2*pi/(86164.1)]*sind(LAT);
end
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
% % Remove a planar fit from each square before smoothing, then put it back:
% function BB = smoothdata2_gaussian_deplaned(AA)
% end

% 1D interpolate over nan, used in the following function:
function TSinterp = interp_over_nan(T,TS,varargin)
    if nargin == 2
        method = 'linear';
    elseif nargin == 3
        method = varargin{1};
    else
        error('Only 2 or 3 arguments are expected.')
    end
    TSnonan = TS(isfinite(TS));
    % get rid of leading and trailing NaN's:
    if length(TSnonan)/length(TS) < 0.5
        TSinterp = TS;
    else
        % TS(1) = TSnonan(1);
        % TS(end) = TSnonan(end);
        TSnonan = TS(isfinite(TS));
        Tnonan = T(isfinite(TS));
        TSinterp = interp1(Tnonan,TSnonan,T,method);
    end    
end
%%
% Interpolate over the nadir gap before smoothing (only works for SWOT
% passes where N-S is the vertical direction):
function BB = smoothdata2_interpnadirgap(XX,YY,AA,nn)
    % Error will occur if any of the inputs are different sizes:
    XX.*YY.*AA;
    for ii = 1:size(AA,1)
        AA(ii,:) = interp_over_nan(XX(ii,:),AA(ii,:),'linear');
    end
    BB = smoothdata2(AA,'gaussian',nn,'omitnan');
end
%%
% Needed to plot streamlines:
function OUT = nan2zero(IN)
    OUT = IN;
    OUT(~isfinite(OUT)) = 0;
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
