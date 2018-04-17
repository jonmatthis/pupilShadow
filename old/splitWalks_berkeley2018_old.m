
close all
clear all
restoredefaultpath


if ispc
    cd E:\Dropbox\UTexas\pupilShadow
    addpath(genpath(cd))
    
    dataPath = 'E:\Dropbox\UTexas\OpticFlowProject';
    cd(dataPath);
    
elseif ismac
    
    cd /Users/matthis/Dropbox/UTexas/pupilShadow
    addpath(genpath(cd))
    
    dataPath = '/Users/matthis/Dropbox/UTexas/KateBerkeley2018';

    cd(dataPath);
end


%%Add the folders relevant to the experiment to the path
addpath(genpath(cd))


%%Add the folders relevant to the experiment to the path
addpath(genpath(cd))
path = cd;


%% Load the session specific details

spotcheck = false;

sessionID = '2018-04-05-S03';


% condID = 'Woodchips';
condID = 'Binocular';


dataPath = strcat(dataPath,filesep,sessionID);
shadowPath = strcat(dataPath,filesep,condID,'/Shadow/');
pupilDataPath = strcat(dataPath,filesep, condID,'/Pupil/');


splitWalks_date = datetime;

cd(dataPath)
cd('OutputFiles')

load(strcat(condID,'.mat'));
%%

saveOutData = true;
debug = true;
spotcheck = false;

%%


%%
%%% clip out data relevant to each walk, zeroing and rotating as necessary

clear allWalks
walks = sesh.walks;
comXYZ = squeeze(shadow_fr_mar_dim(:,1,:));
for ii = 1:length(walks)
    
    
    ii
    
    thisWalk = [];
    
    %% load various bits of data into struct
    thisWalk.avg_fps = round(mean(diff(syncedUnixTime).^-1));
    thisWalk.shadowMarkerNames = shadowMarkerNames;
    thisWalk.calibDist = calibDist;
    thisWalk.legLength = sesh.legLength;
    thisWalk.vorCalibErr = vorCalibErr;
    
    thisWalk.splitWalks_date = datetime;
    thisWalk.processData_date = processData_date ;
    
    thisWalk.sessionID = sessionID;
    
    %% comXYZ
    thisWalk.comXYZ = comXYZ(walks(ii,1): walks(ii,2),:);
    
    zCom = thisWalk.comXYZ(1, :); %this is the original comXYZ start point, used to zero other data.
    
    thisWalk.comXYZ(:, 1) = thisWalk.comXYZ(:, 1) - zCom(1); %zero X data
    thisWalk.comXYZ(:, 3) = thisWalk.comXYZ(:, 3) - zCom(3); %zero Z data
    
    %these won't change for the other data downstream
    pt0 = thisWalk.comXYZ(end, [1 3]); %original endpoint
    pt1 = [1000 0]; %positive-X vector
    origin = thisWalk.comXYZ(1, [1 3]); %startpoint
    
    
    %the data to be rotated
    X = thisWalk.comXYZ(:,1); %original X
    Z = thisWalk.comXYZ(:,3); %original Z (Y)
    
    disp('rotating COM')
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    thisWalk.comXYZ = [x_r' thisWalk.comXYZ(:,2) z_r'];
    
    %% ground fixations (R)
    thisWalk.rGazeGroundIntersection = rGazeGroundIntersection(walks(ii,1): walks(ii,2),:);
    
    thisWalk.rGazeGroundIntersection(:, 1) = thisWalk.rGazeGroundIntersection(:, 1) - zCom(1); %zero X data
    thisWalk.rGazeGroundIntersection(:, 3) = thisWalk.rGazeGroundIntersection(:, 3) - zCom(3); %zero Z data
    
    X = thisWalk.rGazeGroundIntersection(:,1); %original X
    Z = thisWalk.rGazeGroundIntersection(:,3); %original Z (Y)
    
    disp('rotating rGazeGroundIntersection')
    
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    thisWalk.rGazeGroundIntersection = [x_r' thisWalk.rGazeGroundIntersection(:,2) z_r'];
    
    %% ground fixations (L)
    thisWalk.lGazeGroundIntersection = lGazeGroundIntersection(walks(ii,1): walks(ii,2),:);
    
    thisWalk.lGazeGroundIntersection(:, 1) = thisWalk.lGazeGroundIntersection(:, 1) - zCom(1); %zero X data
    thisWalk.lGazeGroundIntersection(:, 3) = thisWalk.lGazeGroundIntersection(:, 3) - zCom(3); %zero Z data
    
    X = thisWalk.lGazeGroundIntersection(:,1); %original X
    Z = thisWalk.lGazeGroundIntersection(:,3); %original Z (Y)
    
    disp('rotating lGazeGroundIntersection')
    
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    thisWalk.lGazeGroundIntersection = [x_r' thisWalk.lGazeGroundIntersection(:,2) z_r'];
    
    %% rEyeballCenterXYZ
    thisWalk.rEyeballCenterXYZ = rEyeballCenterXYZ(walks(ii,1): walks(ii,2),:);
    
    thisWalk.rEyeballCenterXYZ(:, 1) = thisWalk.rEyeballCenterXYZ(:, 1) - zCom(1); %zero X data
    thisWalk.rEyeballCenterXYZ(:, 3) = thisWalk.rEyeballCenterXYZ(:, 3) - zCom(3); %zero Z data
    
    X = thisWalk.rEyeballCenterXYZ(:,1); %original X
    Z = thisWalk.rEyeballCenterXYZ(:,3); %original Z (Y)
    
    disp('rotating rEyeballCenterXYZ')
    
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    
    thisWalk.rEyeballCenterXYZ = [x_r' thisWalk.rEyeballCenterXYZ(:,2) z_r'];
    
    %% lEyeballCenterXYZ
    thisWalk.lEyeballCenterXYZ = lEyeballCenterXYZ(walks(ii,1): walks(ii,2),:);
    
    thisWalk.lEyeballCenterXYZ(:, 1) = thisWalk.lEyeballCenterXYZ(:, 1) - zCom(1); %zero X data
    thisWalk.lEyeballCenterXYZ(:, 3) = thisWalk.lEyeballCenterXYZ(:, 3) - zCom(3); %zero Z data
    
    X = thisWalk.lEyeballCenterXYZ(:,1); %original X
    Z = thisWalk.lEyeballCenterXYZ(:,3); %original Z (Y)
    
    disp('rotating lEyeballCenterXYZ')
    
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    
    thisWalk.lEyeballCenterXYZ = [x_r' thisWalk.lEyeballCenterXYZ(:,2) z_r'];
    
    
    
    %% shadow data
    thisWalk.shadow_fr_mar_dim = shadow_fr_mar_dim(walks(ii,1): walks(ii,2),:,:);
    
    s = thisWalk.shadow_fr_mar_dim;
    
    disp('rotating Marker Data')
    for mm = 1:length(s(1,:,1))
        thisM = squeeze(s(:,mm,:));
        
        thisM(:, 1) = thisM(:, 1) - zCom(1); %zero X data
        thisM(:, 3) = thisM(:, 3) - zCom(3); %zero X data
        
        X = thisM(:,1); %original X
        Z = thisM(:,3); %original Z (Y)
        
        
        [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
        hold on
        s(:,mm,:) = [x_r' thisM(:,2) z_r'];
    end
    
    if spotcheck; dbstack; keyboard; end; hold off
    
    thisWalk.shadow_fr_mar_dim = s;
    
    %% rGazeXYZ
    thisWalk.rGazeXYZ= rGazeXYZ(walks(ii,1): walks(ii,2),:);
    
    
    X = thisWalk.rGazeXYZ(:,1); %original X
    Z = thisWalk.rGazeXYZ(:,3); %original Z (Y)
    
    disp('rotating rGazeXYZ')
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    thisWalk.rGazeXYZ = [x_r' thisWalk.rGazeXYZ(:,2) z_r'];
    
    
    %% lGazeXYZ
    thisWalk.lGazeXYZ= lGazeXYZ(walks(ii,1): walks(ii,2),:);
    
    
    X = thisWalk.lGazeXYZ(:,1); %original X
    Z = thisWalk.lGazeXYZ(:,3); %original Z (Y)
    
    disp('rotating lGazeXYZ')
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    thisWalk.lGazeXYZ = [x_r' thisWalk.lGazeXYZ(:,2) z_r'];
    
    %%%this stuff isn't correctly rotated, I think it's in eye camera
    %%%coordinates or soemthing??
% % %     %% rEye_sphCenCam
% % %     thisWalk.rEye_sphCenCamXYZ(:,1)= rEye_sphCenCam_x(walks(ii,1): walks(ii,2),:);
% % %     thisWalk.rEye_sphCenCamXYZ(:,2)= rEye_sphCenCam_y(walks(ii,1): walks(ii,2),:);
% % %     thisWalk.rEye_sphCenCamXYZ(:,3)= rEye_sphCenCam_z(walks(ii,1): walks(ii,2),:);
% % %     
% % %     X = thisWalk.rEye_sphCenCamXYZ(:,1); %original X
% % %     Z = thisWalk.rEye_sphCenCamXYZ(:,3); %original Z (Y)
% % %     
% % %     disp('rotating rEye_sphCenCamXYZ')
% % %     [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
% % %     
% % %     if spotcheck; dbstack; keyboard; end
% % %     
% % %     thisWalk.rEye_sphCenCamXYZ = [x_r' thisWalk.rEye_sphCenCamXYZ(:,2) z_r'];
% % %         
% % %     %% lEye_sphCenCam
% % %     thisWalk.lEye_sphCenCamXYZ(:,1)= lEye_sphCenCam_x(walks(ii,1): walks(ii,2),:);
% % %     thisWalk.lEye_sphCenCamXYZ(:,2)= lEye_sphCenCam_y(walks(ii,1): walks(ii,2),:);
% % %     thisWalk.lEye_sphCenCamXYZ(:,3)= lEye_sphCenCam_z(walks(ii,1): walks(ii,2),:);
% % %     
% % %     X = thisWalk.lEye_sphCenCamXYZ(:,1); %original X
% % %     Z = thisWalk.lEye_sphCenCamXYZ(:,3); %original Z (Y)
% % %     
% % %     disp('rotating lEye_sphCenCamXYZ')
% % %     [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
% % %     
% % %     if spotcheck; dbstack; keyboard; end
% % %     
% % %     thisWalk.lEye_sphCenCamXYZ = [x_r' thisWalk.lEye_sphCenCamXYZ(:,2) z_r'];
% % %     
% % %         %% rEye_pupCirCen
% % %     thisWalk.rEye_pupCircCenXYZ(:,1)= rEye_pupCircCen_x(walks(ii,1): walks(ii,2),:);
% % %     thisWalk.rEye_pupCircCenXYZ(:,2)= rEye_pupCircCen_y(walks(ii,1): walks(ii,2),:);
% % %     thisWalk.rEye_pupCircCenXYZ(:,3)= rEye_pupCircCen_z(walks(ii,1): walks(ii,2),:);
% % %     
% % %     X = thisWalk.rEye_pupCircCenXYZ(:,1); %original X
% % %     Z = thisWalk.rEye_pupCircCenXYZ(:,3); %original Z (Y)
% % %     
% % %     disp('rotating rEye_pupCircCenXYZ')
% % %     [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
% % %     
% % %     if spotcheck; dbstack; keyboard; end
% % %     
% % %     thisWalk.rEye_pupCircCenXYZ = [x_r' thisWalk.rEye_pupCircCenXYZ(:,2) z_r'];
% % %     
% % %     
% % %     %% lEye_pupCirCen
% % %     thisWalk.lEye_pupCircCenXYZ(:,1)= lEye_pupCircCen_x(walks(ii,1): walks(ii,2),:);
% % %     thisWalk.lEye_pupCircCenXYZ(:,2)= lEye_pupCircCen_y(walks(ii,1): walks(ii,2),:);
% % %     thisWalk.lEye_pupCircCenXYZ(:,3)= lEye_pupCircCen_z(walks(ii,1): walks(ii,2),:);
% % %     
% % %     X = thisWalk.lEye_pupCircCenXYZ(:,1); %original X
% % %     Z = thisWalk.lEye_pupCircCenXYZ(:,3); %original Z (Y)
% % %     
% % %     disp('rotating lEye_pupCircCenXYZ')
% % %     [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
% % %     
% % %     if spotcheck; dbstack; keyboard; end
% % %     
% % %     thisWalk.lEye_pupCircCenXYZ = [x_r' thisWalk.lEye_pupCircCenXYZ(:,2) z_r'];
% % %     
    
    %% headVecX_fr_xyz
    thisWalk.headVecX_fr_xyz= headVecX_fr_xyz(walks(ii,1): walks(ii,2),:);
    
    X = thisWalk.headVecX_fr_xyz(:,1); %original X
    Z = thisWalk.headVecX_fr_xyz(:,3); %original Z (Y)
    
    disp('rotating headVecX_fr_xyz')
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    thisWalk.headVecX_fr_xyz = [x_r' thisWalk.headVecX_fr_xyz(:,2) z_r'];
    
     %% headVecY_fr_xyz
    thisWalk.headVecY_fr_xyz= headVecY_fr_xyz(walks(ii,1): walks(ii,2),:);
    
    X = thisWalk.headVecY_fr_xyz(:,1); %original X
    Z = thisWalk.headVecY_fr_xyz(:,3); %original Z (Y)
    
    disp('rotating headVecY_fr_xyz')
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    thisWalk.headVecY_fr_xyz = [x_r' thisWalk.headVecY_fr_xyz(:,2) z_r'];
    
     %% headVecZ_fr_xyz
    thisWalk.headVecZ_fr_xyz= headVecZ_fr_xyz(walks(ii,1): walks(ii,2),:);
    
    X = thisWalk.headVecZ_fr_xyz(:,1); %original X
    Z = thisWalk.headVecZ_fr_xyz(:,3); %original Z (Y)
    
    disp('rotating headVecZ_fr_xyz')
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    thisWalk.headVecZ_fr_xyz = [x_r' thisWalk.headVecZ_fr_xyz(:,2) z_r'];
    %% step data
    theseStepIDs = steps_HS_TO_StanceLeg_XYZ(:,1)>=walks(ii,1) & steps_HS_TO_StanceLeg_XYZ(:,1)<=walks(ii,2);
    
    thisWalk.steps_HS_TO_StanceLeg_XYZ = steps_HS_TO_StanceLeg_XYZ(theseStepIDs,:);
    
    
    thisWalk.steps_HS_TO_StanceLeg_XYZ(:, 4) = thisWalk.steps_HS_TO_StanceLeg_XYZ(:, 4) - zCom(1); %zero X data
    thisWalk.steps_HS_TO_StanceLeg_XYZ(:, 6) = thisWalk.steps_HS_TO_StanceLeg_XYZ(:, 6) - zCom(3); %zero Z data
    
    X = thisWalk.steps_HS_TO_StanceLeg_XYZ(:,4); %original X
    Z = thisWalk.steps_HS_TO_StanceLeg_XYZ(:,6); %original Z (Y)
    
    disp('rotating Steps')
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin, debug );
    
    if spotcheck; dbstack; keyboard; end
    
    thisWalk.steps_HS_TO_StanceLeg_XYZ(:,[4:6]) = [x_r' thisWalk.steps_HS_TO_StanceLeg_XYZ(:,5) z_r'];
    
    thisWalk.steps_HS_TO_StanceLeg_XYZ(:,1) = thisWalk.steps_HS_TO_StanceLeg_XYZ(:,1) - walks(ii,1);
    thisWalk.steps_HS_TO_StanceLeg_XYZ(:,2) = thisWalk.steps_HS_TO_StanceLeg_XYZ(:,2) - walks(ii,1);
    
    
    %% add non-rotated data
    
    
    thisWalk.rEye_norm_pos_x = rEye.norm_pos_x(walks(ii,1): walks(ii,2),:);
    thisWalk.rEye_norm_pos_y = rEye.norm_pos_y(walks(ii,1): walks(ii,2),:);
    
    thisWalk.lEye_norm_pos_x = lEye.norm_pos_x(walks(ii,1): walks(ii,2),:);
    thisWalk.lEye_norm_pos_y = lEye.norm_pos_y(walks(ii,1): walks(ii,2),:);
    
    thisWalk.rEye_norm_pos_x = rEye_pupRadius(walks(ii,1): walks(ii,2),:);
    thisWalk.lEye_norm_pos_x = lEye_pupRadius(walks(ii,1): walks(ii,2),:);

    
    thisWalk.frames = walks(ii,1): walks(ii,2);
    thisWalk.syncedUnixTime = syncedUnixTime(walks(ii,1): walks(ii,2));
    
    
    
% % % %     %% correct for alignment error
% % % %     [gazeTheta rho] = cart2pol(thisWalk.rGazeXYZ(:,1), thisWalk.rGazeXYZ(:,3));
% % % %     
% % % %     figure(483)
% % % %     histogram(gazeTheta,'Normalization','probability')
% % % %     xlim([-.5 .5])
% % % %     hold on
% % % %     
% % % %     [thisWalk] = correctAlignmentError_opt(thisWalk);
% % % %     
% % % %     [gazeTheta rho] = cart2pol(thisWalk.rGazeXYZ(:,1), thisWalk.rGazeXYZ(:,3));
% % % %     
% % % %     histogram(gazeTheta,'Normalization','probability')
% % % %     xlim([-.5 .5])
% % % %     hold off

% [gazeTheta gazeRho] = cart2pol(thisWalk.rGazeXYZ(:,1), thisWalk.rGazeXYZ(:,3));
% figure(3453)
% histogram(gazeTheta,'Normalization','probability')
% xlim([-.5 .5])
% hold on
%
% [thisWalk.rGazeXYZ(:,1), thisWalk.rGazeXYZ(:,3)] = pol2cart(gazeTheta-nanmean(gazeTheta), gazeRho);
%
% [gazeThetaCorr gazeRho] = cart2pol(thisWalk.rGazeXYZ(:,1), thisWalk.rGazeXYZ(:,3));
% histogram(gazeThetaCorr,'Normalization','probability')

    
    %% load a buncha individual marker datums into the struct
    
    rHeelID = find(strcmp('RightHeel', shadowMarkerNames));
    thisWalk.rHeelXYZ = squeeze(s(:,rHeelID,:)); % pull out rHeelID marker
    
    rToeID = find(strcmp('RightToe', shadowMarkerNames));
    thisWalk.rToeXYZ = squeeze(s(:,rToeID,:)); % pull out rHeelID marker
    
    rFootID = find(strcmp('RightFoot', shadowMarkerNames));
    thisWalk.rFootXYZ = squeeze(s(:,rFootID,:)); % pull out RightFoot marker
    
    
    lHeelID = find(strcmp('LeftHeel', shadowMarkerNames));
    thisWalk.lHeelXYZ = squeeze(s(:,lHeelID,:)); % pull out rHeelID marker
    
    lToeID = find(strcmp('LeftToe', shadowMarkerNames));
    thisWalk.lToeXYZ = squeeze(s(:,lToeID,:)); % pull out rHeelID marker000
    
    lFootID = find(strcmp('LeftFoot', shadowMarkerNames));
    thisWalk.lFootXYZ = squeeze(s(:,lFootID,:)); % pull out LeftFoot marker
    
    
    hTopID = find(strcmp('HeadTop', shadowMarkerNames));
    thisWalk.hTopXYZ= squeeze(s(:,hTopID,:)); % pull out LeftFoot marker
    
    hC1ID = find(strcmp('HeadTop', shadowMarkerNames));
    thisWalk.hC1XYZ= squeeze(s(:,hC1ID,:)); % pull out LeftFoot marker
    
    thisWalk.hCenXYZ = (thisWalk.hTopXYZ + thisWalk.hC1XYZ)/2;
    
    
    
    
    
    
    %%
    allWalks{ii} = thisWalk;
    
    figure(123213)
    
    if mod(ii,2); %% plot even/odd walks on differnet subplots (because even plots are rotated 180 from odd ones)
        subplot(311);
    else
        subplot(312)
    end
    
    plot(thisWalk.comXYZ(:,1), thisWalk.comXYZ(:,3))
    
    hold on
    plot(thisWalk.rGazeGroundIntersection(:,1), thisWalk.rGazeGroundIntersection(:,3),'r.')
    plot(thisWalk.lGazeGroundIntersection(:,1), thisWalk.lGazeGroundIntersection(:,3),'b.')
    
    plot(thisWalk.comXYZ(1,1), thisWalk.comXYZ(1,3),'gp')
    plot(thisWalk.comXYZ(end,1), thisWalk.comXYZ(end,3),'rp')
    plot(thisWalk.steps_HS_TO_StanceLeg_XYZ(:,4), thisWalk.steps_HS_TO_StanceLeg_XYZ(:,6),'o')
    
    plot(thisWalk.steps_HS_TO_StanceLeg_XYZ(1,4), thisWalk.steps_HS_TO_StanceLeg_XYZ(1,6),'go')
    plot(thisWalk.steps_HS_TO_StanceLeg_XYZ(end,4), thisWalk.steps_HS_TO_StanceLeg_XYZ(end,6),'ro')
    
    axis equal
    
    if mod(ii,2);
        title(strcat(condID,' - odd walks'))
    else
        title(strcat(condID,' - even walks'))
    end
    
    title(condID)
    
    subplot(313)
    plot(comXYZ)
    hold on
    plot(walks(ii,1):walks(ii,2), comXYZ(walks(ii,1):walks(ii,2),3),'o')
    drawnow
    beep
    
    %%
end


%%
if saveOutData
    
    walks = allWalks;
    walkIDs = sesh.walkID;
    save('allWalks.mat','walks')
    
end

