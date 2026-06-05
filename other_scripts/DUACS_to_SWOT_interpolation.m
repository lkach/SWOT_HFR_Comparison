%% DUACS Product
% "Global Ocean Gridded L4 Sea Surface Heights And Derived Variables
% Reprocessed 1993 Ongoing"

% Downloaded from:
% https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_MY_008_047/download?dataset=cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D_202411

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

% % % % % % % % % % %
% SWOT:
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
SSHA = [];
for ti_swot = 1:length(NORCAL.SWOT.ssha_karin_2)
    SSHA(:,:,ti_swot) = [NORCAL.SWOT.ssha_karin_2{ti_swot} + NORCAL.SWOT.height_cor_xover{ti_swot}] .* ...
                        [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}]./ ...
                        [~NORCAL.SWOT.ssha_karin_2_qual{ti_swot} .* ~NORCAL.SWOT.height_cor_xover_qual{ti_swot}];
end
NORCAL.SWOT.mean_SSHA = mean(SSHA,3,'omitnan'); % clear SSHA

NORCAL.GEOSTR = load('SWOT_and_HFR_velocities_8pix_weighted.mat');

% % % % % % % % % % %
% HFR
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

% % % % % % % % % % %
% 1/4 degree
DUACS_4 = ncreadall('./DUACS/c3s_obs-sl_glo_phy-ssh_my_twosat-l4-duacs-0.25deg_P1D_1773452684889.nc');
% Global Attributes:
%            Conventions      = "CF-1.11"
%            title            = "DT merged two satellites Global Ocean Gridded SSALTO/DUACS Sea Surface Height L4 product and derived variables"
%            institution      = "CLS, CNES"
%            source           = "Altimetry measurements"
%            history          = "2024/09/23 20:26:01 pva_axp@trex014.sis.cnes.fr Import depuis MSLA.nc"
%            contact          = "http://climate.copernicus.eu/c3s-user-service-desk"
%            references       = "http://climate.copernicus.eu"
%            comment          = "Sea Surface Height measured by Altimetry and derived variables"
%            subset:source    = "ARCO data downloaded from the Marine Data Store using the MyOcean Data Portal"
%            subset:productId = "SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057"
%            subset:datasetId = "c3s_obs-sl_glo_phy-ssh_my_twosat-l4-duacs-0.25deg_P1D_202411"
%            subset:date      = "2026-03-14T01:44:44.890Z"

% 1/8 degree
DUACS_8 = ncreadall('./DUACS/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D_1773450091510.nc');
% Global Attributes:
%            Conventions      = "CF-1.11"
%            title            = "DT merged all satellites Global Ocean Gridded SSALTO/DUACS Sea Surface Height L4 product and derived variables"
%            institution      = "CLS, CNES"
%            source           = "Altimetry measurements"
%            history          = "2024-10-23 12:55:06Z: Creation"
%            contact          = "servicedesk.cmems@mercator-ocean.eu"
%            references       = "http://marine.copernicus.eu"
%            comment          = "Sea Surface Height measured by Altimetry and derived variables"
%            subset:source    = "ARCO data downloaded from the Marine Data Store using the MyOcean Data Portal"
%            subset:productId = "SEALEVEL_GLO_PHY_L4_MY_008_047"
%            subset:datasetId = "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D_202411"
%            subset:date      = "2026-03-14T01:01:31.510Z"

% % % % % % % % % % %
% Flip the x- and y-axes because the NC files arrange them inconveniently:
VARS = {'adt','err_sla','err_ugosa','err_vgosa','flag_ice','sla',...
        'ugos','ugosa','vgos','vgosa'};
STRUCTS = {'DUACS_4','DUACS_8'};
for ii = 1:length(STRUCTS)
    for jj = 1:length(VARS)
        eval([STRUCTS{ii} '.' VARS{jj} ' = permute(' STRUCTS{ii} '.' VARS{jj} ',[2 1 3]);']);
    end
end

[DUACS_4.LON, DUACS_4.LAT] = meshgrid(DUACS_4.longitude, DUACS_4.latitude);
[DUACS_8.LON, DUACS_8.LAT] = meshgrid(DUACS_8.longitude, DUACS_8.latitude);

% Both have the "time" variable with units:
% units         = "seconds since 1970-01-01 00:00:00"

time0 = datenum('1970-01-01 00:00:00');
DUACS_4.time = time0 + DUACS_4.time/[24*60*60];
DUACS_8.time = time0 + DUACS_8.time/[24*60*60];

%% Examine visually

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
ti_swot   = dsearchn(NORCAL.SWOT.mean_time, datenum(USER_PICKED_DATE));
ti_duacs4 = dsearchn(DUACS_4.time, NORCAL.SWOT.mean_time(ti_swot));
ti_duacs8 = dsearchn(DUACS_8.time, NORCAL.SWOT.mean_time(ti_swot));

close all

% figure('Color','w')
% pcolor_centered(DUACS_4.LON,DUACS_4.LAT,DUACS_4.sla(:,:,ti_duacs4));
% axis tight; axis equal

figure('Color','w')
tiledlayout(1,3,"TileSpacing","tight")

nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
    'longitudes',[-126.5 -121.5],...
    'latitudes',[36.5 42.5]);
COAST = m_gshhs_l('patch',0.8*[1 1 1]);
m_grid('box','fancy', 'backgroundcolor','none'); hold on
set(gcf,'Position',[262    69   962   728])
m_pcolor_centered(DUACS_4.LON,DUACS_4.LAT,DUACS_4.sla(:,:,ti_duacs4));
CB = colorbar; CB.Label.String = 'SLA (m)'; CB.Label.Interpreter = 'latex';
CB.FontSize = 12; CB.Label.FontSize = 16;
clim([-1 1]*0.2)
colormap(turbo)
m_text(-123.4, 41.7, '1/4$^\circ$', 'fontsize', 18, 'Interpreter', 'latex')
title(['DUACS SLA 1/4^o at ' datestr(DUACS_4.time(ti_duacs4))])

nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
    'longitudes',[-126.5 -121.5],...
    'latitudes',[36.5 42.5]);
COAST = m_gshhs_l('patch',0.8*[1 1 1]);
m_grid('box','fancy', 'backgroundcolor','none'); hold on
set(gcf,'Position',[262    69   962   728])
m_pcolor_centered(DUACS_8.LON,DUACS_8.LAT,DUACS_8.sla(:,:,ti_duacs8));
CB = colorbar; CB.Label.String = 'SLA (m)'; CB.Label.Interpreter = 'latex';
CB.FontSize = 12; CB.Label.FontSize = 16;
clim([-1 1]*0.2)
colormap(turbo)
m_text(-123.4, 41.7, '1/8$^\circ$', 'fontsize', 18, 'Interpreter', 'latex')
title(['DUACS 1/8^o at ' datestr(DUACS_8.time(ti_duacs8))])

nexttile
M_P = m_proj('Lambert',... <--- ALTERNATIVE, MORE LIMITED BOUNDS
    'longitudes',[-126.5 -121.5],...
    'latitudes',[36.5 42.5]);
COAST = m_gshhs_l('patch',0.8*[1 1 1]);
m_grid('box','fancy', 'backgroundcolor','none'); hold on
set(gcf,'Position',[262    69   962   728])
m_pcolor_centered(NORCAL.SWOT.lon{ti_swot}, NORCAL.SWOT.lat{ti_swot}, ...
                  SSHA(:,:,30));
CB = colorbar; CB.Label.String = 'SSHA (m)'; CB.Label.Interpreter = 'latex';
CB.FontSize = 12; CB.Label.FontSize = 16;
clim([-1 1]*0.2)
colormap(turbo)
m_text(-123.4, 41.7, '2 km', 'fontsize', 18, 'Interpreter', 'latex')
title(['SWOT SSHA\_Karin\_2 at ' datestr(NORCAL.SWOT.mean_time(ti_swot))])

%% Interpolate to the HFR grid

% Shorter variable name for SWOT grid
LON_SWOT = NORCAL.SWOT.lon{1}; % They're the same across time steps
LAT_SWOT = NORCAL.SWOT.lat{1};

% Shorter variable name for HFR grid
LON_HFR = NORCAL.HFR.LON;
LAT_HFR = NORCAL.HFR.LAT;

% NORCAL.SWOT
%         cross_track_distance: {1×89 cell}
%            distance_to_coast: {1×89 cell}
%             height_cor_xover: {1×89 cell}
%        height_cor_xover_qual: {1×89 cell}
%           internal_tide_hret: {1×89 cell}
%                          lat: {1×89 cell}
%                          lon: {1×89 cell}
%     mean_sea_surface_cnescls: {1×89 cell}
%                         name: {1×89 cell}
%               ocean_tide_fes: {1×89 cell}
%                 ssha_karin_2: {1×89 cell}
%            ssha_karin_2_qual: {1×89 cell}
%                    swh_karin: {1×89 cell}
%               swh_karin_qual: {1×89 cell}
%                         time: {1×89 cell}
%                    mean_time: [89×1 double]
%                    mean_SSHA: [285×69 double]

% NORCAL.GEOSTR
%                   CG_Iteration: 1
%                 ConvWithWindow: @(IN)nanconv(squeeze(IN),WINDOW,'same','edge')
%                         PIXELS: 8
%                      SWOT_time: [89×1 double]
%              U_cyclogeostr_Nit: [285×69×89 double]
%                       U_geostr: [285×69×89 double]
%       U_rot_filtered_SWOTtimes: [94×56×89 double]
%     U_rot_unfiltered_SWOTtimes: [94×56×89 double]
%           Ucg_SWOT_HFRgrid_all: [94×56×89 double]
%            Ug_SWOT_HFRgrid_all: [94×56×89 double]
%              V_cyclogeostr_Nit: [285×69×89 double]
%                       V_geostr: [285×69×89 double]
%       V_rot_filtered_SWOTtimes: [94×56×89 double]
%     V_rot_unfiltered_SWOTtimes: [94×56×89 double]
%           Vcg_SWOT_HFRgrid_all: [94×56×89 double]
%            Vg_SWOT_HFRgrid_all: [94×56×89 double]
%                         WINDOW: [36×1 double]

% DUACS_4 = 
%               time: [106×1 double]
%           latitude: [20×1 single]
%          longitude: [16×1 single]
%                adt: [20×16×106 double]
%            err_sla: [20×16×106 double]
%          err_ugosa: [20×16×106 double]
%          err_vgosa: [20×16×106 double]
%           flag_ice: [20×16×106 double]
%                sla: [20×16×106 double]
%     tpa_correction: [106×1 double]
%               ugos: [20×16×106 double]
%              ugosa: [20×16×106 double]
%               vgos: [20×16×106 double]
%              vgosa: [20×16×106 double]
%                LON: [20×16 single]
%                LAT: [20×16 single]

% NORCAL.HFR
%            hdop: [94×56×2211 single]
%        time_run: [2211×1 double]
%            time: [2211×1 double]
%             lat: [94×1 single]
%             lon: [56×1 single]
%           wgs84: -127
%               u: [94×56×2211 single]
%               v: [94×56×2211 single]
%     u_mean_2023: [94×56 single]
%     v_mean_2023: [94×56 single]
%             LON: [94×56 single]
%             LAT: [94×56 single]



% warning('Actually, look to the files that contain "U_geostr" and related variables, not the plain SSHA.')

DUACS_4_struct = struct;
% DUACS_4_struct.readme = ['Note that the variable names are not necessarily accurate,' ...
%     ' but rather are named for the closest analogue to SWOT variables so that they can be dropped' ...
%     ' into existing scripts written to compare HFR and SWOT. For example, DUACS has "sla", which I' ...
%     ' assign to the SWOT variable "ssha_karin_2".'];
DUACS_4_struct.LON_SWOT = LON_SWOT;
DUACS_4_struct.LAT_SWOT = LAT_SWOT;
DUACS_4_struct.LON_HFR = LON_HFR;
DUACS_4_struct.LAT_HFR = LAT_HFR;
DUACS_4_struct.time = DUACS_4.time;
DUACS_4_struct.U_geostr = [];
DUACS_4_struct.V_geostr = [];
for ii = 1:length(DUACS_4_struct.time)
    % vgos =  "Absolute geostrophic velocity: meridian component"
    % vgosa = "Geostrophic velocity anomalies: meridian component"
    DUACS_4_struct.U_geostr(:,:,ii) = interp2(DUACS_4.LON,DUACS_4.LAT,DUACS_4.ugos(:,:,ii),DUACS_4_struct.LON_SWOT,DUACS_4_struct.LAT_SWOT,'linear');
    DUACS_4_struct.V_geostr(:,:,ii) = interp2(DUACS_4.LON,DUACS_4.LAT,DUACS_4.vgos(:,:,ii),DUACS_4_struct.LON_SWOT,DUACS_4_struct.LAT_SWOT,'linear');

    DUACS_4_struct.U_geostr_HFRgrid(:,:,ii) = interp2(DUACS_4.LON,DUACS_4.LAT,DUACS_4.ugos(:,:,ii),DUACS_4_struct.LON_HFR,DUACS_4_struct.LAT_HFR,'linear');
    DUACS_4_struct.V_geostr_HFRgrid(:,:,ii) = interp2(DUACS_4.LON,DUACS_4.LAT,DUACS_4.vgos(:,:,ii),DUACS_4_struct.LON_HFR,DUACS_4_struct.LAT_HFR,'linear');
end

DUACS_8_struct = struct;
% DUACS_8_struct.readme = ['Note that the variable names are not necessarily accurate,' ...
%     ' but rather are named for the closest analogue to SWOT variables so that they can be dropped' ...
%     ' into existing scripts written to compare HFR and SWOT. For example, DUACS has "sla", which I' ...
%     ' assign to the SWOT variable "ssha_karin_2".'];
DUACS_8_struct.LON_SWOT = LON_SWOT;
DUACS_8_struct.LAT_SWOT = LAT_SWOT;
DUACS_8_struct.LON_HFR = LON_HFR;
DUACS_8_struct.LAT_HFR = LAT_HFR;
DUACS_8_struct.time = DUACS_8.time;
DUACS_8_struct.U_geostr = [];
DUACS_8_struct.V_geostr = [];
for ii = 1:length(DUACS_8_struct.time)
    % vgos =  "Absolute geostrophic velocity: meridian component"
    % vgosa = "Geostrophic velocity anomalies: meridian component"
    DUACS_8_struct.U_geostr(:,:,ii) = interp2(DUACS_8.LON,DUACS_8.LAT,DUACS_8.ugos(:,:,ii),DUACS_8_struct.LON_SWOT,DUACS_8_struct.LAT_SWOT,'linear');
    DUACS_8_struct.V_geostr(:,:,ii) = interp2(DUACS_8.LON,DUACS_8.LAT,DUACS_8.vgos(:,:,ii),DUACS_8_struct.LON_SWOT,DUACS_8_struct.LAT_SWOT,'linear');

    DUACS_8_struct.U_geostr_HFRgrid(:,:,ii) = interp2(DUACS_8.LON,DUACS_8.LAT,DUACS_8.ugos(:,:,ii),DUACS_8_struct.LON_HFR,DUACS_8_struct.LAT_HFR,'linear');
    DUACS_8_struct.V_geostr_HFRgrid(:,:,ii) = interp2(DUACS_8.LON,DUACS_8.LAT,DUACS_8.vgos(:,:,ii),DUACS_8_struct.LON_HFR,DUACS_8_struct.LAT_HFR,'linear');
end



%%
close all

ti_duacs = 70;
ti_swot = dsearchn(NORCAL.SWOT.mean_time, DUACS_4.time(ti_duacs));
clickBool = false; % Set to false if you want to predefine the examined points' coordinates

% figure
% subplot(2,2,1);
% pcolor_centered(DUACS_4_struct.LON_SWOT,DUACS_4_struct.LAT_SWOT,DUACS_4_struct.U_geostr(:,:,ii));
% shading flat; axis tight; axis equal
% subplot(2,2,2);
% pcolor_centered(DUACS_4_struct.LON_HFR,DUACS_4_struct.LAT_HFR,DUACS_4_struct.U_geostr_HFRgrid(:,:,ii));
% shading flat; axis tight; axis equal
% subplot(2,2,3);
% pcolor_centered(DUACS_8_struct.LON_SWOT,DUACS_8_struct.LAT_SWOT,DUACS_8_struct.U_geostr(:,:,ii));
% shading flat; axis tight; axis equal
% subplot(2,2,4);
% pcolor_centered(DUACS_8_struct.LON_HFR,DUACS_8_struct.LAT_HFR,DUACS_8_struct.U_geostr_HFRgrid(:,:,ii));
% shading flat; axis tight; axis equal

figure
tiledlayout(2,3,"TileSpacing","tight")

% subplot(3,4,[1 5 9])
% imagesc(DUACS_8_struct.U_geostr_HFRgrid(:,:,ii));
% [x_click,y_click] = ginput(1); x_click = round(x_click); y_click = round(y_click);

nexttile([2 1]) % subplot(1,4,[1])
pcolor_centered(DUACS_8_struct.LON_HFR,DUACS_8_struct.LAT_HFR,DUACS_8_struct.U_geostr_HFRgrid(:,:,ti_duacs)); hold on
pcolor_centered(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},         NORCAL.GEOSTR.U_geostr(:,:,ti_swot));
% contour(NORCAL.SWOT.lon{1},NORCAL.SWOT.lat{1},                 NORCAL.GEOSTR.U_geostr(:,:,ti_swot),'k');
axis equal; axis tight
clim([-1 1]*0.5)
% 
MARKER = {'*k','+k'}
if clickBool
    [lon_click ,lat_click ] = ginput(2);
else
    lon_click = [-124.55072, -124.17593];
    lat_click = [40.0131416, 38.8804016];
end
y_click = []; x_click = [];
y_click_swot = []; x_click_swot = [];
for kk = 1:length(lon_click)
    DIST_MAT = abs([DUACS_8_struct.LON_HFR + 1i*DUACS_8_struct.LAT_HFR] - [lon_click(kk) + 1i*lat_click(kk)]);
    DIST_MAT = DIST_MAT == min(DIST_MAT(:));
    % [~,y_click] = min(DIST_MAT ); y_click = unique(y_click);
    % [~,x_click] = min(DIST_MAT'); x_click = unique(x_click);
    ROWS_COLUMNS = [1:size(DUACS_8_struct.LON_HFR,1)]' + 1i*[1:size(DUACS_8_struct.LON_HFR,2)];
    y_click = [y_click, real(ROWS_COLUMNS(DIST_MAT))];
    x_click = [x_click, imag(ROWS_COLUMNS(DIST_MAT))];
    DIST_MAT_SWOT = abs([NORCAL.SWOT.lon{1} + 1i*NORCAL.SWOT.lat{1}] - [lon_click(kk) + 1i*lat_click(kk)]);
    DIST_MAT_SWOT = DIST_MAT_SWOT == min(DIST_MAT_SWOT(:));
    ROWS_COLUMNS = [1:size(NORCAL.SWOT.lon{1},1)]' + 1i*[1:size(NORCAL.SWOT.lon{1},2)];
    y_click_swot = [y_click_swot, real(ROWS_COLUMNS(DIST_MAT_SWOT))];
    x_click_swot = [x_click_swot, imag(ROWS_COLUMNS(DIST_MAT_SWOT))];
    plot(lon_click(kk),lat_click(kk),MARKER{kk})
end



% nexttile([1 2]) % subplot(1,4,[2:4])
% plot(NORCAL.HFR.time/24 + T0,squeeze(  NORCAL.HFR.u(y_click,x_click,:)),'.-','Color',[1 1 1]*0.9,'LineWidth',1.5); hold on
% 
% plot(NORCAL.SWOT.mean_time, squeeze(NORCAL.GEOSTR.U_rot_filtered_SWOTtimes(y_click,x_click,:)), 'k-','LineWidth',1.5)
% plot(NORCAL.SWOT.mean_time, squeeze(NORCAL.GEOSTR.U_rot_unfiltered_SWOTtimes(y_click,x_click,:)), 'k--','LineWidth',1.5,'HandleVisibility','off')
% 
% plot(DUACS_4_struct.time, squeeze( DUACS_4_struct.U_geostr_HFRgrid(y_click,x_click,:)), 'b.-','LineWidth',1.5);
% plot(DUACS_8_struct.time, squeeze( DUACS_8_struct.U_geostr_HFRgrid(y_click,x_click,:)), 'r.-','LineWidth',1.5)
% plot(NORCAL.SWOT.mean_time, squeeze(NORCAL.GEOSTR.U_geostr(y_click_swot,x_click_swot,:)), 'g.-','LineWidth',1.5)
% datetick
% legend('HFR','HFR, rot., filtered, @ SWOT times','DUACS 1/4^o','DUACS 1/8^o','SWOT')
% set(gcf,'Position',[-1919         473        1920         499])

for kk = 1:2
    y_clickN = y_click(kk);
    x_clickN = x_click(kk);
    y_clickN_swot = y_click_swot(kk);
    x_clickN_swot = x_click_swot(kk);
    nexttile([1 2]) % Second point's time series
    plot(NORCAL.HFR.time/24 + T0,squeeze(  NORCAL.HFR.v(y_clickN,x_clickN,:)),'.-','Color',[1 1 1]*0.9,'LineWidth',1.5); hold on

    plot(NORCAL.SWOT.mean_time, squeeze(NORCAL.GEOSTR.V_rot_filtered_SWOTtimes(y_clickN,x_clickN,:)), 'k-','LineWidth',1.5)
    % plot(NORCAL.SWOT.mean_time, squeeze(NORCAL.GEOSTR.V_rot_unfiltered_SWOTtimes(y_clickN,x_clickN,:)), 'k--','LineWidth',1.5,'HandleVisibility','off')

    plot(DUACS_4_struct.time, squeeze( DUACS_4_struct.V_geostr_HFRgrid(y_clickN,x_clickN,:)), 'b.-','LineWidth',1.5);
    plot(DUACS_8_struct.time, squeeze( DUACS_8_struct.V_geostr_HFRgrid(y_clickN,x_clickN,:)), 'r.-','LineWidth',1.5)
    plot(NORCAL.SWOT.mean_time, squeeze(NORCAL.GEOSTR.V_geostr(y_clickN_swot,x_clickN_swot,:)), 'g.-','LineWidth',1.5)
    xlim([datenum('2023-04-01') datenum('2023-08-01')])
    datetick
    legend('HFR','HFR, rot., filtered, @ SWOT times','DUACS 1/4^o','DUACS 1/8^o','SWOT')
end

set(gcf,'Position',[-1919         473        1920         499])



%% Wavenumber spectra of DUACS (just 1/8 degree)

warning('Likely not necessary, this will be a steep spectrum anyway.')

%% Save the gridded output for later convenient use:

save(['./DUACS_4_gridded_to_HFR_and_SWOT_grids.mat'],'-struct','DUACS_4_struct')
save(['./DUACS_8_gridded_to_HFR_and_SWOT_grids.mat'],'-struct','DUACS_8_struct')

%% Auxiliary functions
%% ncreadall(NC_file) loads all of the variables in a .nc and names them in
% MATLAB the same as they are named in the file. WARNING: There is no way
% to tell if the .nc file is very large (i.e. if its contents would fill of
% the memory of the MATLAB workspace). Therefore, be aware of the size of
% the file before using this function.
% 
% IN:   NC_file  = string of the .nc file name, e.g. 'data.nc'
%       EXCLUDED = (optional) cell of strings of names of variables to be
%                  excluded
% 
% OUT:  var_struct = OPTIONAL, this is a structure with the variables
%                    stored in it. If you choose not to have an output
%                    for ncreadall, then the variables in NC_file simply
%                    go to your workspace separately, rather than being
%                    bundled into a structure.
% Luke Kachelein - February 21, 2018
% 
% For use to meet course requirements and perform academic research at
% degree-granting institutions only. It is not available for government,
% commercial, or other organizational use, unless granted permission by the
% original author.
function var_struct = ncreadall(NC_file,varargin)
if nargin == 1
    EXCLUDED = {''};
elseif nargin == 2
    EXCLUDED = varargin{1};
else
    error('One or two inputs expected.')
end
if sum(NC_file((end-2):end) == '.nc') == 3
else
    help ncreadall
    error('NCREADALL only accepts .nc files. You must terminate the string argument "NC_file" with ".nc". Read the instructions.')
end
INFO = ncinfo(NC_file);
variables_struct = INFO.Variables; % this will be a struct
N = length(variables_struct); % number of variables
if length(EXCLUDED) == 1 && strcmp(EXCLUDED{1},'')
    variables_cell = cell(1,N);
else
    variables_cell = {};
end
n_var = 1;
for n = 1:N
    % Exclude variables if request:
    if ~sum(strcmp(variables_struct(n).Name,EXCLUDED)) % 0 if the n'th variable is in EXCLUDED
        variables_cell{n_var} = variables_struct(n).Name;
        n_var = n_var + 1;
    else
    end
end
% Redefine N after eliminating possible undesired variables:
N = length(variables_cell);
% If no output is expected, just define variables in workspace
if nargout == 0
    for n = 1:N
        eval([variables_cell{n},' = ncread(''',NC_file,''',''',variables_cell{n},''');']);
        eval(['assignin(''base'',''',variables_cell{n},''',',variables_cell{n},');']);
    end
% If one is expected, define a structure
elseif nargout == 1
    var_struct = struct;
    for n = 1:N
        eval(['var_struct.',variables_cell{n},' = ncread(''',NC_file,''',''',variables_cell{n},''');']);
    end
else
    help ncreadall
    error('NCREADALL allows only one or zero output arguments. Read the instructions.')
end
end
%% PCOLOR_CENTERED(first 3 basic pcolor arguments)
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