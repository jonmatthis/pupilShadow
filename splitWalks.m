function[walks] = splitWalks(sessionPath,takeID,walks,walk_names, varargin)
% function[walks] = splitWalksFun(sessionPath,takeID,walks,walk_names, varargin)
%
% Loads a mat file from data already processed from an experiment involving
% pupil labs eyetracker and a shadow walking suit.  It takes the processed
% data arrays and separates them into individual walks.
%
% Inputs:
% sessionPath - where the data is located
% sessionID - name of the session
% takeID - name of the take
% useEye - which eyes were used in this experiment both = [1,1], right =
%    [1,0], left = [0,1]
% sesh - this is a struct which contains useful information about the
%    session being split. (TODO: write down exactly what info it must have)
% walks - nx2 array; n is the number of walks. The first column is where
%    each walk starts and the second is where each walk ends.
% walk_names - nx1 cell array containing the name/label for each walk.
%
% optional Name/Variable input pairs:
%
% 'Transfer' - data to transfer to individual walks
% 'Split' - data that needs to be split up and put into individual walks
% 'Rotate' - data that needs to be split up, then rotated so the starting
%     point is at the origin, then put into individual walks
% 'CorrectAlignmentFcn' - A function that accepts the current walks struct
%     and does some small corrections on the alignment
%% define some path locations
outputPath = [sessionPath, filesep, 'OutputFiles'];
takePath = [sessionPath filesep takeID];
shadowDataPath = [takePath filesep 'Shadow' filesep];
pupilDataPath = [takePath filesep 'Pupil' filesep];

% load data
fn = [sessionPath filesep 'OutputFiles' filesep takeID '.mat'];
disp(['loading data: ' fn '...']);
load(fn);
%% process optional function parameters
defaultTransfer = { 'sessionID','takeID', 'subID', 'shadowMarkerNames', 'processData_date',...
    'framerate', 'legLength','calibDist', 'px2mmScale'};
lEyeData = {'lVorCalibErr'};
rEyeData = {'rVorCalibErr'};
if useEye(1), defaultTransfer = [defaultTransfer, rEyeData]; end
if useEye(2), defaultTransfer = [defaultTransfer, lEyeData]; end


defaultSplit = {'syncedUnixTime','gaze.norm_pos_x','gaze.norm_pos_y'...
    'headAccXYZ','chestAccXYZ','hipsAccXYZ','worldFrameIndex'};
lEyeData = {'lEye.theta','lEye.phi','lEye.norm_pos_x',...
    'lEye.norm_pos_y','lEye.circle_3d_radius','lEye.blinks'};
rEyeData = {'rEye.theta','rEye.phi','rEye.norm_pos_x',...
    'rEye.norm_pos_y','rEye.circle_3d_radius','rEye.blinks'};
if useEye(1), defaultSplit = [defaultSplit, rEyeData]; end
if useEye(2), defaultSplit = [defaultSplit, lEyeData]; end


defaultRotate = {'comXYZ','headVecX_fr_xyz','headVecY_fr_xyz','headVecZ_fr_xyz'};
rEyeData = {'rGazeGroundIntersection','rEyeballCenterXYZ','rGazeXYZ'};
lEyeData = {'lGazeGroundIntersection', 'lEyeballCenterXYZ','lGazeXYZ'};
if useEye(1), defaultRotate = [defaultRotate, rEyeData]; end
if useEye(2), defaultRotate = [defaultRotate, lEyeData]; end


p = inputParser();
p.addParameter('Transfer',defaultTransfer);
p.addParameter('Split',defaultSplit);
p.addParameter('Rotate',defaultRotate);
% p.addParameter('CorrectAlignmentFcn',@correctAlignmentError_opt);
p.addParameter('CorrectAlignmentFcn',@correctAlignmentError_lookUpTable);
parse(p,varargin{:});

data_to_transfer = p.Results.Transfer;
data_to_split = p.Results.Split;
data_to_split_rotate = p.Results.Rotate;
corrAlign = p.Results.CorrectAlignmentFcn;

%%

comXYZ = squeeze(shadow_fr_mar_dim(:,1,:));

for ww = 1:size(walks,1)
    
    thisWalk = [];
    thisWalk.ww = ww;
    thisWalk.splitWalks_date = datetime;
    thisWalk.name = walk_names{ww};
    thisWalk.eyes = useEye;
    thisWalk.frames = walks(ww,1): walks(ww,2);
    
    % transfer data
    for dd = 1:length(data_to_transfer)
        varName = data_to_transfer{dd};
        disp(['Adding ' varName '...']);
        varValue = eval(varName);
        idx = find(varName=='.');
        varName(idx) = '_';
        thisWalk.(varName) = varValue;
    end
    
    
    % prepare to rotate a bunch of data based on center of mass
    thisWalk.comXYZ=comXYZ(walks(ww,1):walks(ww,2),:);
    origin = [thisWalk.comXYZ(1, 1) 0 thisWalk.comXYZ(1, 3)]; %startpoint
    
    if strcmp(thisWalk.name,'vor') % don't rotate VOR data
        pt0 = calibPoint([1 3]); %for the VOR frames, use the Calib point to define the "end point"
        pt1 = pt0; %don't rotate VOR data
        thisWalk.isThisVORCalibrationData = true;
    else
        pt0 = thisWalk.comXYZ(end, [1 3]); %original endpoint
        pt1 = [1000 0]; %positive-X vector
        thisWalk.isThisVORCalibrationData = false;
    end
    
    % rotate data
    for dd = 1:length(data_to_split_rotate)
        varName = data_to_split_rotate{dd};
        disp(['Rotating ' varName '...']);
        varValue = eval(varName);
        
        
        X = varValue(walks(ww,1): walks(ww,2), 1) - origin(1); %original X (translated to origin)
        Y = varValue(walks(ww,1): walks(ww,2), 2);
        Z = varValue(walks(ww,1): walks(ww,2), 3) - origin(3); %original Z (translated to origin)
        
        
        debug = true;
        [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin([1 3]), debug);
        thisWalk.(varName) = [x_r', Y, z_r'];
    end
    
    % rotate shadowmarker data (it's special)
    disp('Rotating Marker Data...')
    thisWalk.shadow_fr_mar_dim = shadow_fr_mar_dim(walks(ww,1): walks(ww,2),:,:);
    s = thisWalk.shadow_fr_mar_dim;
    for mm=1:size(s,2)
        thismarker = squeeze(s(:,mm,:));
        thismarker(:,[1,3]) = thismarker(:,[1,3])-origin([1 3]);
        X = thismarker(:,1);
        Z = thismarker(:,3);
        [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin([1 3]), 0 );
        s(:,mm,:) = [x_r' thismarker(:,2) z_r'];
    end
    thisWalk.shadow_fr_mar_dim = s;
    
    % load some indvidual data markers in separately (e.g., RightHeel ->
    % rHeelXYZ)
    markers = {'RightHeel','RightToe','RightFoot','LeftHeel','LeftToe',...
        'LeftFoot','HeadTop'};
    for mm=1:length(markers)
        id = (strcmp(markers{mm},shadowMarkerNames));
        upperCase = find(markers{mm} >= 'A' & markers{mm} <= 'Z');
        newName = [lower(markers{mm}(upperCase(1))),markers{mm}(upperCase(2):end),'XYZ'];
        thisWalk.(newName) = squeeze(s(:,id,:));
    end
    
    %rotate step data
    disp('Rotating Steps...')
    theseStepIDs = steps_HS_TO_StanceLeg_XYZ(:,1)>=walks(ww,1) & steps_HS_TO_StanceLeg_XYZ(:,1)<=walks(ww,2);
    thisWalk.steps_HS_TO_StanceLeg_XYZ = steps_HS_TO_StanceLeg_XYZ(theseStepIDs,:);
    thisWalk.steps_HS_TO_StanceLeg_XYZ(:, [4,6]) = thisWalk.steps_HS_TO_StanceLeg_XYZ(:, [4,6]) - origin([1 3]); %zero data
    X = thisWalk.steps_HS_TO_StanceLeg_XYZ(:,4); %original X
    Z = thisWalk.steps_HS_TO_StanceLeg_XYZ(:,6); %original Z (Y)
    [x_r, z_r] = rotateFromV0toV1(X, Z, pt0, pt1, origin([1 3]),0 );
    thisWalk.steps_HS_TO_StanceLeg_XYZ(:,4:6) = [x_r' thisWalk.steps_HS_TO_StanceLeg_XYZ(:,5) z_r'];
    thisWalk.steps_HS_TO_StanceLeg_XYZ(:,1:2) = thisWalk.steps_HS_TO_StanceLeg_XYZ(:,1:2) - walks(ww,1:2);
    
    % add unrotated split data
    for dd = 1:length(data_to_split)
        varName = data_to_split{dd};
        disp(['Adding ' varName '...']);
        varValue = eval(varName);
        idx = find(varName=='.');
        varName(idx) = '_';
        thisWalk.(varName) = varValue(walks(ww,1): walks(ww,2),:);
    end
    
    
    thisWalk.ww = ww; %need this info for corrAlignLookUpTable
    thisWalk.sesh = sesh;
    thisWalk.trialType = sesh.trialType{ww};
    
    % correct for alignment error should this be somewhere else?
    thisWalk_orig = thisWalk;
    thisWalk = corrAlign(thisWalk_orig);
    
    
    %%% final debug plot. Gaze data should be well aligned with the
    %%% COM path
    figure(sum(double(sessionID))+sum(double(takeID))) %dumb way to get a unique number for each session/take combo
    %     if mod(ii,2); %% plot even/odd walks on differnet subplots (because even plots are rotated 180 from odd ones)
    %         subplot(311);
    %     else
    %         subplot(312)
    %     end
    subplot(round(length(walks(2:end,:))/2),2,ww)
    
    
    hold on
    plot(thisWalk.rGazeGroundIntersection(:,1), thisWalk.rGazeGroundIntersection(:,3),'r.')
    plot(thisWalk.lGazeGroundIntersection(:,1), thisWalk.lGazeGroundIntersection(:,3),'b.')
    
    plot(thisWalk.comXYZ(1,1), thisWalk.comXYZ(1,3),'gp')
    plot(thisWalk.comXYZ(end,1), thisWalk.comXYZ(end,3),'rp')
    
    plot(thisWalk.steps_HS_TO_StanceLeg_XYZ(1,4), thisWalk.steps_HS_TO_StanceLeg_XYZ(1,6),'go')
    plot(thisWalk.steps_HS_TO_StanceLeg_XYZ(end,4), thisWalk.steps_HS_TO_StanceLeg_XYZ(end,6),'ro')
    
        plot(thisWalk.comXYZ(:,1), thisWalk.comXYZ(:,3),'g','LineWidth',2)

    axis equal
    
    title(strcat(sessionID,{' walk# '},num2str(ww),'-',thisWalk.trialType))
    
    drawnow
    
    
    
    
    allWalks{ww} = thisWalk; %push thisWalk into the allWalks cell array
    
end

disp(['Saving: ' outputPath filesep takeID '_allWalks.mat'])
save([outputPath filesep takeID '_allWalks.mat'],'allWalks','sesh');