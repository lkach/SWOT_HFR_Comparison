%% The following needs to be run on the server because it will take a very long time to run
addpath ~; lkaddpath
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
CurrentDir = dir(mfilePath);
cd(CurrentDir.folder)
Antenna_ID = 'GCVE';
load(['./workspace_' Antenna_ID '.mat'],...
     'NORCAL',...
     'SWOT_radial_unit_x',...
     'SWOT_radial_unit_y',...
     'Radial'...
     )
clear U_geostr V_geostr
load(['./Dist_to_nearest_radial_' Antenna_ID '.mat'])

SWOTQuadraticFitDir = dir(['./SWOT_and_HFR_velocities_*pix_weighted.mat']);

%% 
%% Calculate the difference in velocity at each radar point over time,
%  and the correlation coefficient between the two.
%% First figure out how far the closest radial is from each SWOT point,
% which will be used in the next step.
% This should take ~30 minutes to complete:
if exist('Dist_to_nearest_radial','var')
    % error('You already calculated these variables, which take a long time to calculate.')
    disp('done part 1 (pre-loaded Dist_to_nearest_radial)')
else
disp('Calculating distance to nearest radials on SWOT grid.')
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
% cd ../data
save(['./Dist_to_nearest_radial_' Antenna_ID '.mat'],'Dist_to_nearest_radial','Antenna_ID')
disp('done part 1')
end
%% Now compare the SWOT and HFR radial velocities head-to-head
% Make this on the SWOT grid, with NaN where there are not enough nearby
% radar points to give useful information. I have to do it as a loop; maybe
% there is a more efficient way, but I only need an answer.
% if exist('Radial_vel_rmsd','var')
%     error('You already calculated these variables, which take a long time to calculate.')
% else
% end
KM_search = 3.0; % search radius for averaging radial velocities to compare to a SWOT grid point
SWOT_mean_time = NORCAL.SWOT.mean_time; % For naming simplicity


for pixel_i = 1:length(SWOTQuadraticFitDir)
    load([SWOTQuadraticFitDir(pixel_i).folder '/' ...
          SWOTQuadraticFitDir(pixel_i).name]);
    
    Radial_g_rmsd = nan(size(NORCAL.SWOT.lat{1}));
    Radial_g_corrR = nan(size(NORCAL.SWOT.lat{1}));
    Radial_g_corrP = nan(size(NORCAL.SWOT.lat{1}));
    Radial_g_lowpass_rmsd = nan(size(NORCAL.SWOT.lat{1}));
    Radial_g_lowpass_corrR = nan(size(NORCAL.SWOT.lat{1}));
    Radial_g_lowpass_corrP = nan(size(NORCAL.SWOT.lat{1}));
    Radial_cg_rmsd = nan(size(NORCAL.SWOT.lat{1}));
    Radial_cg_corrR = nan(size(NORCAL.SWOT.lat{1}));
    Radial_cg_corrP = nan(size(NORCAL.SWOT.lat{1}));
    Radial_cg_lowpass_rmsd = nan(size(NORCAL.SWOT.lat{1}));
    Radial_cg_lowpass_corrR = nan(size(NORCAL.SWOT.lat{1}));
    Radial_cg_lowpass_corrP = nan(size(NORCAL.SWOT.lat{1}));
    
    WINDOW = hanning(36); WINDOW = WINDOW/sum(WINDOW); % low pass
    ConvWithWindow = @(IN) conv(squeeze(IN),WINDOW,'same');
    
    for ii = 1:size(NORCAL.SWOT.lat{1},1)
        for jj = 1:size(NORCAL.SWOT.lat{1},2)
            % ii = 124; jj = 50; %$
            if sum(isfinite(U_geostr(ii,jj,:)),3)/size(U_geostr,3) >= 0.9 & Dist_to_nearest_radial(ii,jj) < 5*KM_search*1000
    
                SWOTg_radial_vel_ts = squeeze(U_geostr(ii,jj,:)).*SWOT_radial_unit_x(ii,jj) + ...
                                      squeeze(V_geostr(ii,jj,:)).*SWOT_radial_unit_y(ii,jj);
                SWOTcg_radial_vel_ts = squeeze(U_cyclogeostr_Nit(ii,jj,:)).*SWOT_radial_unit_x(ii,jj) + ...
                                       squeeze(V_cyclogeostr_Nit(ii,jj,:)).*SWOT_radial_unit_y(ii,jj);
    
                % [Distance_from_SWOTcell,~,~] = m_idist(NORCAL.SWOT.lon{1}(ii,jj),NORCAL.SWOT.lat{1}(ii,jj), ...
                %                                        Radial.LON,Radial.LAT);
                %                                        % Old, slow way
    
                % New, fast way that quickly eliminates most of the obviously
                % unneeded m_idist calculations:
                DistInDeg = abs([NORCAL.SWOT.lon{1}(ii,jj) + 1i*NORCAL.SWOT.lat{1}(ii,jj)] - [Radial.LON + 1i*Radial.LAT]);
                Distance_from_SWOTcell = nan(size(Radial.LON));
                [Distance_from_SWOTcell_,~,~] = m_idist(NORCAL.SWOT.lon{1}(ii,jj),NORCAL.SWOT.lat{1}(ii,jj), ...
                                                        Radial.LON(DistInDeg<0.2),Radial.LAT(DistInDeg<0.2));
                Distance_from_SWOTcell(DistInDeg<0.2) = Distance_from_SWOTcell_;
    
                IND = [Distance_from_SWOTcell/1000 < KM_search];
    
                T_HFR_at_point = unique(Radial.TIME(IND));
                if length(T_HFR_at_point) > 10 % this just makes sure that T_HFR_at_point isn't empty or nearly so
                Ur_HFR_at_point = nan(length(T_HFR_at_point),1);
                % UrSTD_HFR_at_point = nan(length(T_HFR_at_point),1);
                
                Radial_VEL_IND = Radial.VEL(IND);
                TIME_at_IND = Radial.TIME(IND);
                for kk = 1:length(T_HFR_at_point)
                    Ur_HFR_at_point(kk) = mean(Radial_VEL_IND(TIME_at_IND==T_HFR_at_point(kk))/100,'omitnan');
                end

                Ur_HFR_at_point_atSWOTtimes = interp1(T_HFR_at_point,Ur_HFR_at_point,SWOT_mean_time,'linear');;
                
                % Full Radials vs. G vel.
                Radial_g_rmsd(ii,jj) = rms(Ur_HFR_at_point_atSWOTtimes - SWOTg_radial_vel_ts,'omitnan');
                [corr_R,corr_P] = ...
                    corrcoef(Ur_HFR_at_point_atSWOTtimes(isfinite(Ur_HFR_at_point_atSWOTtimes) & isfinite(SWOTg_radial_vel_ts)), ...
                             SWOTg_radial_vel_ts(        isfinite(Ur_HFR_at_point_atSWOTtimes) & isfinite(SWOTg_radial_vel_ts)));
                if isscalar(corr_R)
                    Radial_g_corrR(ii,jj) = NaN;
                    Radial_g_corrP(ii,jj) = NaN;
                else
                    Radial_g_corrR(ii,jj) = corr_R(1,2);
                    Radial_g_corrP(ii,jj) = corr_P(1,2);
                end
    
                % Full Radials vs. CG vel.
                Radial_cg_rmsd(ii,jj) = rms(Ur_HFR_at_point_atSWOTtimes - SWOTcg_radial_vel_ts,'omitnan');
                [corr_R,corr_P] = ...
                    corrcoef(Ur_HFR_at_point_atSWOTtimes(isfinite(Ur_HFR_at_point_atSWOTtimes) & isfinite(SWOTcg_radial_vel_ts)), ...
                             SWOTcg_radial_vel_ts(       isfinite(Ur_HFR_at_point_atSWOTtimes) & isfinite(SWOTcg_radial_vel_ts)));
                if isscalar(corr_R)
                    Radial_cg_corrR(ii,jj) = NaN;
                    Radial_cg_corrP(ii,jj) = NaN;
                else
                    Radial_cg_corrR(ii,jj) = corr_R(1,2);
                    Radial_cg_corrP(ii,jj) = corr_P(1,2);
                end
    
    
                % Define Low Pass Radials
                Ur_lowpass = ConvWithWindow(Ur_HFR_at_point);
                Ur_HFR_at_point_atSWOTtimes_lowpass = interp1(T_HFR_at_point,Ur_lowpass,SWOT_mean_time,'linear');;

                % Low Pass Radials vs. G vel.
                Radial_g_lowpass_rmsd(ii,jj) = rms(Ur_HFR_at_point_atSWOTtimes_lowpass - SWOTg_radial_vel_ts,'omitnan');
                [lowpass_corr_R,lowpass_corr_P] = ...
                    corrcoef(Ur_HFR_at_point_atSWOTtimes_lowpass(isfinite(Ur_HFR_at_point_atSWOTtimes_lowpass) & isfinite(SWOTg_radial_vel_ts)), ...
                             SWOTg_radial_vel_ts(                isfinite(Ur_HFR_at_point_atSWOTtimes_lowpass) & isfinite(SWOTg_radial_vel_ts)));
                if isscalar(lowpass_corr_R)
                    Radial_g_lowpass_corrR(ii,jj) = NaN;
                    Radial_g_lowpass_corrP(ii,jj) = NaN;
                else
                    Radial_g_lowpass_corrR(ii,jj) = lowpass_corr_R(1,2);
                    Radial_g_lowpass_corrP(ii,jj) = lowpass_corr_P(1,2);
                end

                % Low Pass Radials vs. CG vel.
                Radial_cg_lowpass_rmsd(ii,jj) = rms(Ur_HFR_at_point_atSWOTtimes_lowpass - SWOTcg_radial_vel_ts,'omitnan');
                [lowpass_corr_R,lowpass_corr_P] = ...
                    corrcoef(Ur_HFR_at_point_atSWOTtimes_lowpass(isfinite(Ur_HFR_at_point_atSWOTtimes_lowpass) & isfinite(SWOTcg_radial_vel_ts)), ...
                             SWOTcg_radial_vel_ts(               isfinite(Ur_HFR_at_point_atSWOTtimes_lowpass) & isfinite(SWOTcg_radial_vel_ts)));
                if isscalar(lowpass_corr_R)
                    Radial_cg_lowpass_corrR(ii,jj) = NaN;
                    Radial_cg_lowpass_corrP(ii,jj) = NaN;
                else
                    Radial_cg_lowpass_corrR(ii,jj) = lowpass_corr_R(1,2);
                    Radial_cg_lowpass_corrP(ii,jj) = lowpass_corr_P(1,2);
                end
    
                else
    
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
    % cd ../data
    save(['./Radial_vel_comparison_' Antenna_ID '_quadfit' num2str(PIXELS) 'pix_weighted.mat'], ...
          'Radial_g_lowpass_rmsd','Radial_g_lowpass_corrR','Radial_g_lowpass_corrP',...
          'Radial_g_rmsd','Radial_g_corrR','Radial_g_corrP',...
          'Radial_cg_lowpass_rmsd','Radial_cg_lowpass_corrR','Radial_cg_lowpass_corrP',...
          'Radial_cg_rmsd','Radial_cg_corrR','Radial_cg_corrP',...
          'Antenna_ID','KM_search','PIXELS')
    disp('done')
end

%%
