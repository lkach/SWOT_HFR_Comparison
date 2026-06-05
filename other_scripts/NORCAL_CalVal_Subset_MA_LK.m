% Used on JPL Eddy server to subset data from where all available data are
% stored.
% 
% Originally written by Matthew Archer (329B) for the same purpose,
% modified by Luke Kachelein.

% ls -alht /mnt/flow/swot/KaRIn/SWOT_L2_LR_SSH_2.0_Expert/*Expert_4*.nc | head
% ^ Yields Apr 7, 2023 as the earliest day in the cal/val record

% vim /mnt/flow/swot/Analysis_Luke/HFR_NORCAL/NORCAL_CalVal_Subset_MA_LK.m

%% Load and Subset SWOT data
clear;clc
% Get file list
files_4xx = dir('/mnt/flow/swot/KaRIn/SWOT_L2_LR_SSH_2.0_AllToDate/*Expert_4*.nc');
files_5xx = dir('/mnt/flow/swot/KaRIn/SWOT_L2_LR_SSH_2.0_AllToDate/*Expert_5*.nc');
% ^ Use Expert to get SWH
% For some reason, the first of these files is corrupted, so eliminate it manually:
files_4xx = files_4xx(2:end);

files = files_4xx;
% Combine 400 and 500 SWOT records from the cal/val orbit into one
% structure:
for ii = 1:length(files_5xx)
    files(end+1) = files_5xx(ii);
end

lenfiles = length(files); % get # of files
% Define California Current System (CCS) area
% ln = 50; ls = 30; le = -116; lw = -133; % lat/lon limits for CCS
% ln = 38.5; ls = 33; le = -116; lw = -133; % lat/lon limits for CCS
ln = 42; ls = 37; le = -122; lw = -125.4; % lat/lon limits for CCS
count=0;
for ii = 1:lenfiles % foreach file (i.e. for each SWOT pass)
    disp([num2str(ii) ' out of ' num2str(lenfiles)]) % Display progress
    fname = [files(ii).folder '/' files(ii).name]; % Construct filename
    % Load in Lon/Lat data
    lon1 = ncread(fname,'longitude');lon1(lon1>180) = lon1(lon1>180)-360;
    lat1 = ncread(fname,'latitude');
    % Conditional indexing - find only CCS region
    cswath = 35; % get center track of swath
    lonc = lon1(cswath,:); latc = lat1(cswath,:); % reduce lon/lat to a 1D track for indexing purposes
    indc = find(lonc >= lw & lonc <= le & latc >= ls & latc <= ln); % index to retain only swot data in CCS
    if length(indc) > 1 % if there's data in region
        count = count + 1;
        % Subset to region (creates a cell array)
        name{count} = files(ii).name;
        time{count} = ncread(fname,'time',indc(1),length(indc))/60/60/24 + datenum(2000,1,1);
        lon{count} = lon1(:,indc);
        lat{count} = lat1(:,indc);
        % Load in and save only key variables from SWOT files
        ssha_karin_2{count} = ncread(fname,'ssha_karin_2',[1 indc(1)],[size(lon1,1) length(indc)]);
        ssha_karin_2_qual{count} = ncread(fname,'ssha_karin_2_qual',[1 indc(1)],[size(lon1,1) length(indc)]);
        mean_sea_surface_cnescls{count} = ncread(fname,'mean_sea_surface_cnescls',[1 indc(1)],[size(lon1,1) length(indc)]);
        internal_tide_hret{count} = ncread(fname,'internal_tide_hret',[1 indc(1)],[size(lon1,1) length(indc)]);
        ocean_tide_fes{count} = ncread(fname,'ocean_tide_fes',[1 indc(1)],[size(lon1,1) length(indc)]);
        distance_to_coast{count} = ncread(fname,'distance_to_coast',[1 indc(1)],[size(lon1,1) length(indc)]);

        % % Added by LK
        % lon_full{count} = lon1;
        % lat_full{count} = lat1;
        % time_full{count} = ncread(fname,'time')/60/60/24 + datenum(2000,1,1);
        % ssha_full{count} = ncread(fname,'ssha_karin_2');
        % flag_full{count} = ncread(fname,'ssha_karin_2_qual');
        height_cor_xover{count} = ncread(fname,'height_cor_xover',[1 indc(1)],[size(lon1,1) length(indc)]);
        height_cor_xover_qual{count} = ncread(fname,'height_cor_xover_qual',[1 indc(1)],[size(lon1,1) length(indc)]);

        % % Added by LK later than the above, in order to add SWH information from the Expert files:
        swh_karin{count} =      ncread(fname,'swh_karin',[1 indc(1)],[size(lon1,1) length(indc)]);
        swh_karin_qual{count} = ncread(fname,'swh_karin_qual',[1 indc(1)],[size(lon1,1) length(indc)]);
        cross_track_distance{count} = ncread(fname,'cross_track_distance',[1 indc(1)],[size(lon1,1) length(indc)]);
    end
end
save('/mnt/flow/swot/Analysis_Luke/SWOT_HFR_Overlap/NORCAL_SWOTdata_CCS.mat', ...
    'name', 'time', 'lon', 'lat', 'ssha_karin_2', 'ssha_karin_2_qual', 'mean_sea_surface_cnescls', 'internal_tide_hret', 'ocean_tide_fes', 'distance_to_coast', ...
    ... %'ssha_full', 'lat_full', 'lon_full', 'flag_full', 'time_full', ...
    'height_cor_xover', 'height_cor_xover_qual', 'swh_karin', 'swh_karin_qual', 'cross_track_distance' ...
    );
%%
error('Forced error to allow only file saving to be run, and not perform plotting.')
%% Load in subsetted data, perform ad-hoc quality control, and plot it up
clear;clc
warning off
% Load subsetted SWOT data
load SWOTdata_CCS
% Load coast data
load coastlines
% Define California Current System (CCS) area
ln = 45; ls = 25; le = -110; lw = -140; % lat/lon limits for CCS
figure % initialize figure window
%
for ii = 1:length(time) % for each pass
    disp([num2str(ii) ' out of ' num2str(length(time))]) % Display progress
    % Grab variables at time i
    sshai = ssha_karin_2{ii};
    d2ca = distance_to_coast{ii};
    flaga = ssha_karin_2_qual{ii};
    lona = lon{ii};
    lata = lat{ii};
    % Ad-hoc Quality Control
    sshai(d2ca > 10000000 | d2ca <= 1000 | isnan(d2ca) ) = NaN; % Remove data near coast
    sshai(sshai > 3 | sshai < -3) = NaN; % Remove data with unrealistic amplitudes
    sshai(flaga > 200) = NaN; % ssha(flaga~=0) = NaN; % Use flag variable to remove bad data
    sshf = detrend(sshai','omitnan')'; % detrend to remove long wavelength error
    % Plot coastlines
    plot(coastlon,coastlat),hold on
    % Plot SWOT swath pass
    pcolor(lona,lata,sshf), shading flat, hold on
    caxis([-.35 .35])
    % Clean figure
    set(gcf,'color','w'),set(gca,'Layer','Top','fontsize',20), box on, grid on
    xlim([lw le]); ylim([ls ln])
    Y=get(gca,'ylim');set(gca,'dataaspectratio',[1 cos((Y(2)+Y(1))/2/180*pi) 1])
end
