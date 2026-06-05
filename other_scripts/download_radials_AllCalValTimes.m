% Examine radial data directly

% Determine antena curator and name at:
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

% % % Example antenna:
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

%% CODAR QA/QC manual:
% https://cdn.ioos.noaa.gov/media/2017/12/cos_qaqc_document.pdf

%% 

% Desired download directory:
DOWNLOAD_DIR = '/mnt/flow/swot/Analysis_Luke/HFR_NORCAL/';
cd(DOWNLOAD_DIR)

% % Load SWOT flyover times to match (only start and finish at this point):
% NORCAL.SWOT = load('../NORCAL_SWOTdata_CCS.mat','time');
% NORCAL.SWOT.mean_time = nan(length(NORCAL.SWOT.time),1);
% for ii = 1:length(NORCAL.SWOT.mean_time)
%     NORCAL.SWOT.mean_time(ii) = mean(NORCAL.SWOT.time{ii}(:),'omitnan');
% end
% time_start = floor(NORCAL.SWOT.mean_time(1));
% time_end   = floor(NORCAL.SWOT.mean_time(end)) + 1;
time_start = 738984; % '08-Apr-2023'
time_end = 739077; % '10-Jul-2023'
time_vector = time_start:[1/24]:time_end;
list_of_times = datestr(time_vector,'yyyy-mm-dd HH');

% Template for download link (North to South):

% BML_antenna = 'SHEL';
% BML_antenna = 'BRAG';
% BML_antenna = 'PAFS';

BML_antennas = {'SHEL','BRAG','PAFS'};

% Full path to wget on my machine; yours may vary:
% WGET = '/opt/homebrew/bin/wget';
WGET = '/usr/bin/wget';

%% Download only what you need:
% 1) wget the file from the server
% 2) save only the part you need
% 3) overwrite the remaining text file with the next one to save drive
%    space


for ant_i = 1:length(BML_antennas)
    BML_antenna = BML_antennas{ant_i};
    URL_template = ['https://www.ncei.noaa.gov/data/oceans/ndbc/hfradar/radial/' ...
        'YYYY/YYYYMM/BML/' BML_antenna '/RDL_i_BML_' BML_antenna '_YYYY_MM_DD_HH00.ruv'];
    TIME = [];
    LON = [];
    LAT = [];
    VEL = [];

    for i_url = 1:length(list_of_times)
        URL_i = replace(URL_template,{'YYYY','MM','DD','HH'},...
            {list_of_times(i_url,1:4),...
             list_of_times(i_url,6:7),...
             list_of_times(i_url,9:10),...
             list_of_times(i_url,12:13)});
        STATUS = system([WGET ' "' URL_i '" -O ' DOWNLOAD_DIR 'temp.ruv'])

        if STATUS % i.e. wget failed
            % Nothing ever happens
        else

            DATA_struct = importdata([DOWNLOAD_DIR 'temp.ruv'],' ',55);
            if i_url == 1
                Variable_names = [DATA_struct.textdata{end-1} '   Time(UTC)'];
                Origin_antenna = str2num(replace(replace(DATA_struct.textdata{10},{'%Origin:  '},{'1i*'}),{' '},{' + '}));
            else
            end
            TimeStamp = datenum( ...
                replace( ...
                replace(DATA_struct.textdata{7},{'%TimeStamp: ','  '},{'','$'}), ...
                {' ','$'},{'-',' '}) , 'yyyy-mm-dd HH');
            DATA = [DATA_struct.data, TimeStamp*ones(size(DATA_struct.data,1),1)];

            % % % Decide what data to throw away:
            % rows 6 and 7 are "Spatial Quality" and "Temporal Quality", which
            % are standard deviations (see manual). They appear mostly
            % uncorrelated, so I filter out data with abs(SQ + i*TQ) < 95th percentile,
            % which is about 50 cm/s. This cutoff is arbitrary.
            SQTQ_Cutoff = 60;
            % row 10 minus row 11 is the max minus min velocity that contribute
            % to the average; again filter by using only diff_max_min < 95th percentile,
            % which is about 92 cm/s. This cutoff is arbitrary.
            DiffMaxMin_Cutoff = 92;
            IND = [abs(DATA(:,6) + 1i*DATA(:,7)) < SQTQ_Cutoff & ...
                [DATA(:,10) - DATA(:,11)] < DiffMaxMin_Cutoff];

            TIME = [TIME ; DATA(IND,21)];
            LON = [LON ; DATA(IND,1)];
            LAT = [LAT ; DATA(IND,2)];
            VEL = [VEL ; DATA(IND,18)];

            disp(datestr(TimeStamp))
        end
    end
    save([DOWNLOAD_DIR 'Radial_data_allcalvaltimes_' BML_antenna '.mat'],...
         'TIME','LON','LAT','VEL','Origin_antenna','Variable_names','SQTQ_Cutoff','DiffMaxMin_Cutoff')
end


