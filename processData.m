function[] = processData(sessionID,takeID,sessionPath,useEye,sesh)
% function[] = processDataFun(sessionID,takeID,sessionPath,useEye,sesh)
%
% Loads data from pupil labs eye tracker and shadow mobile motion capture
% system.  Aligns data from both systems in order to calculate gaze ground
% intersections.
%
% Inputs:
% sessionID - name of the session
% takeID - name of the take
% sessionPath - where the data is located
% useEye - which eyes were used in this experiment both = [1,1], right =
%    [1,0], left = [0,1]
% sesh - this is a struct which contains useful information about the
%    session being split. (TODO: write down exactly what info it must have)

processData_date = datetime;
outputPath = [sessionPath, filesep, 'OutputFiles'];
takePath = [sessionPath filesep takeID];
shadowTakeName = sesh.shadowTakeName;

subID = sesh.subID;
vorFrames = sesh.vorFrames;
calibFrame = vorFrames(1);
legLength = sesh.legLength;
shadowDataPath = [takePath filesep 'Shadow' filesep];
pupilDataPath = [takePath filesep 'Pupil' filesep];


%% get pupil start time date & world time stamps

% PUPIL

worldTimestamps = readNPY([pupilDataPath 'world_timestamps.npy']);

delimiter = ','; formatSpec = '%q%q%[^\n\r]';

fileID = fopen([pupilDataPath 'info.csv'],'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
key = dataArray{:, 1};
value = dataArray{:, 2};
res = str2double(strsplit(value{strcmp(key, 'World Camera Resolution' )},'x'));
resWidth = res(1);
resHeight = res(2);


%% resample pupil data to enforce constant framerate (also grab porX porY)
% PUPIL

pupilExportPath = [pupilDataPath filesep 'exports' filesep sesh.pupilExportDir];
%%%%% WE NEED TO CHOOSE THE CONVENTION FOR THIS %%%%%%%%%%%

pupilUnixStartTime = str2double(value{strcmp(key, 'Start Time (System)')});
prefRate = 120;

pupTable = readtable(strcat(pupilExportPath,filesep,'pupil_positions.csv'));
gazeTable = readtable(strcat(pupilExportPath,filesep,'gaze_positions.csv'));

startTime = max([pupTable.timestamp(1) gazeTable.timestamp(1)]);
endTime = min([pupTable.timestamp(end) gazeTable.timestamp(end)]);

skipThese = {'id','index' 'method','model_id', 'timestamp'};
if useEye(1)
    rEye = fixSamplingRate(pupTable(pupTable.id==0,:),startTime,endTime,...
        pupilUnixStartTime, prefRate, skipThese);
    % give a couple things more meaningful names
    rEye.pupRadius = rEye.circle_3d_radius;
    
end
if useEye(2)
    lEye = fixSamplingRate(pupTable(pupTable.id==1,:),startTime,endTime,...
        pupilUnixStartTime, prefRate, skipThese);
    % give a couple things more meaningful names
    lEye.pupRadius = lEye.circle_3d_radius;
end

skipThese = {'base_data','index','timestamp','gaze_normal1_z'};
gaze = fixSamplingRate(gazeTable,startTime,endTime,pupilUnixStartTime,...
    prefRate,skipThese);
pupUnixTime = gaze.unixTimestamp;

% porX and porY
[index, norm_pos_x, norm_pos_y] = pull_index_normposx_normposy(strcat(pupilExportPath,filesep,'gaze_positions.csv'));
[porX porY] =  downsamplegaze(norm_pos_x,norm_pos_y,index,resHeight,resWidth);


%% Find frame of world video that most closely matches each rEye, lEye and gaze data frame
worldUnixTimestamp = (worldTimestamps-worldTimestamps(1)) + pupilUnixStartTime;

disp('Finding World Frame Indices')
if useEye(1)
    rEye.index = knnsearch(worldUnixTimestamp,rEye.unixTimestamp);
end
if useEye(2)
    lEye.index = knnsearch(worldUnixTimestamp,lEye.unixTimestamp);
end
gaze.index = knnsearch(worldUnixTimestamp,gaze.unixTimestamp);

%% Load Shadow Data - SHADOW ONLY
disp('Loading Shadow Data ...')

streamFilename = [shadowDataPath, filesep, shadowTakeName, '_stream.csv'];
streamData = readtable(streamFilename);  % this is way better than importdata because it auto fills some stuff for us

% if isfield(streamData,'Head_time')==0
%     if isempty(streamData.time)==1
%         streamData.timestamp = streamData.Body_time;
%         
%     else
%         streamData.timestamp = streamData.time;     % shadow V3
%     end
% else
%     streamData.timestamp = streamData.Head_time;  % shadow V2
% end

if any(strcmp('Body_time', streamData.Properties.VariableNames))
    streamData.time = streamData.Body_time;
elseif any(strcmp('Head_time', streamData.Properties.VariableNames))
    streamData.time = streamData.Head_time;
elseif any(strcmp('timestamp', streamData.Properties.VariableNames))
   streamData.time = streamData.timestamp; 
else
    streamData.time = streamData.time;
end

streamData.timestamp = streamData.time;
%% Find ShadoUnixStartTime - SHADOW ONLY

takeStruct = xml2struct([shadowDataPath, filesep,'take.mTake']); % turns out '.mTake' files are just secret XML's!! :O
shadowStartDateTime = takeStruct.take.start;
t = strsplit(shadowStartDateTime.Text,{'-','T',':','Z'});
shadowStartDateTime = datetime(str2num(t{1}),...
    str2num(t{2}),...
    str2num(t{3}),...
    str2num(t{4}),...
    str2num(t{5}),...
    str2num(t{6}));
shadowUnixStartTime = posixtime(shadowStartDateTime);

%% Shadow Resampling - SHADOW ONLY
skipThese = {'markerData','og_time'};
shadowDataResamp = fixSamplingRate(streamData,streamData.time(1),...
    streamData.time(end), shadowUnixStartTime, prefRate,skipThese);

shadowUnixTime = shadowDataResamp.timestamp + shadowUnixStartTime;

%% Build a skeleton (linear translation data for all markers) - SHADOW ONLY

xx = find(endsWith(streamData.Properties.VariableNames,'ltx'));
yy = find(endsWith(streamData.Properties.VariableNames,'lty'));
zz = find(endsWith(streamData.Properties.VariableNames,'ltz'));

shadowMarkerNames = streamData.Properties.VariableNames(xx);
for i=1:length(shadowMarkerNames)
    shadowMarkerNames{i} = shadowMarkerNames{i}(1:end-4);
    
    % shadow V3
    if strcmp(shadowMarkerNames{i},'HeadEnd')==1
        shadowMarkerNames{i}='HeadTop';
    end
end

streamShadow_fr_mar_dim = nan(size(shadowDataResamp.unixTimestamp,1), length(xx), 3);

for mm=1:length(xx)
    streamShadow_fr_mar_dim(:,mm,1) = shadowDataResamp.(streamData.Properties.VariableNames{xx(mm)})*10;     %%% multiply by 10 to convert from cm to mm
    streamShadow_fr_mar_dim(:,mm,2) = shadowDataResamp.(streamData.Properties.VariableNames{yy(mm)})*10;
    streamShadow_fr_mar_dim(:,mm,3) = shadowDataResamp.(streamData.Properties.VariableNames{zz(mm)})*10;
end

shadowDataResamp.markerData = streamShadow_fr_mar_dim;

%% Trim pupil and Shadow data

disp('Trimming pupil & shadow data to be the same chunks of time...');
firstFrameTime = max([pupUnixTime(1) shadowUnixTime(1)]);
lastFrameTime = min([pupUnixTime(end) shadowUnixTime(end)]);

if lastFrameTime - firstFrameTime < 0
    warning('Something went wrong with the timestamps...')
    keyboard
end

trimPupFrames = pupUnixTime<firstFrameTime |  pupUnixTime>lastFrameTime;
trimShadowFrames = shadowUnixTime<firstFrameTime | shadowUnixTime>lastFrameTime;

% the trimming might be off by a frame if things aren't perfectly aligned
% We'll remove the extra frame.
if sum(~trimPupFrames)==sum(~trimShadowFrames)+1
    idx = find(~trimPupFrames,1,'last');
    trimPupFrames(idx) = true;
elseif sum(~trimPupFrames)+1==sum(~trimShadowFrames)
    idx = find(~trimShadowFrames,1,'last');
    trimShadowFrames(idx) = true;
end

if useEye(1)
    rEye = trimFrames(rEye,trimPupFrames,{});
end

if useEye(2)
    lEye = trimFrames(lEye,trimPupFrames,{});
end
gaze = trimFrames(gaze,trimPupFrames,{});

skipThese = {'timestamp','markerData'};
shadowDataTrimmed = trimFrames(shadowDataResamp,trimShadowFrames,skipThese);
shadowDataTrimmed.markerData = shadowDataTrimmed.markerData(~trimShadowFrames,:,:);

% UM why do we have rt and st????
rt = pupUnixTime;
rt(trimPupFrames) = [];

st = shadowUnixTime;
st(trimShadowFrames) = [];
syncedUnixTime = mean([st rt]')';
framerate = round(mean(diff(syncedUnixTime).^-1));

shadowRAW_fr_mar_dim = shadowDataTrimmed.markerData;

%% SH Pull out head/chest/hips acceleration -- SHADOW ONLY

headAccXYZ = [shadowDataTrimmed.Head_lax, shadowDataTrimmed.Head_lay, shadowDataTrimmed.Head_laz];
chestAccXYZ = [shadowDataTrimmed.Chest_lax, shadowDataTrimmed.Chest_lay, shadowDataTrimmed.Chest_laz];
hipsAccXYZ = [shadowDataTrimmed.Hips_lax, shadowDataTrimmed.Hips_lay, shadowDataTrimmed.Hips_laz];


%% Find Steps -- SHADOW ONLY
disp('Finding steps')
wRaw.shadow_fr_mar_dim = shadowRAW_fr_mar_dim;
wRaw.shadowMarkerNames = shadowMarkerNames;
wRaw.avg_fps = mean(diff(syncedUnixTime).^-1);
% wRaw.walks = walks;

[allSteps_HS_TO_StanceLeg] = ZeniStepFinder(wRaw);

%% Try to fix 'skateboarding' problem by pinning the feet to the ground during each step
%SHADOW ONLY
disp('Fixing Skateboards')
[shadow_fr_mar_dim] = fixSkateboarding_kb(wRaw, allSteps_HS_TO_StanceLeg);
comXYZ = squeeze(shadow_fr_mar_dim(:,1,:));

%% build step_TO_HS_ft_XYZ variable
% SHADOW ONLY

% get l/r heel information
rHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('RightHeel', shadowMarkerNames),:)); % pull out lHeel marker
lHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('LeftHeel', shadowMarkerNames),:)); % pull out lHeel marker


steps_HS_TO_StanceLeg_XYZ = nan(length(allSteps_HS_TO_StanceLeg), 6);
idx = allSteps_HS_TO_StanceLeg(:,3) == 1; % right foot on ground
steps_HS_TO_StanceLeg_XYZ(idx,:) = [allSteps_HS_TO_StanceLeg(idx,1:3),rHeelXYZ(allSteps_HS_TO_StanceLeg(idx,1),1:3)];
steps_HS_TO_StanceLeg_XYZ(~idx,:) = [allSteps_HS_TO_StanceLeg(~idx,1:3),lHeelXYZ(allSteps_HS_TO_StanceLeg(~idx,1),1:3)];



if true %%debug plot of sorts
    figure(6484); clf;
    
    plot(steps_HS_TO_StanceLeg_XYZ(idx,4), steps_HS_TO_StanceLeg_XYZ(idx,6), 'ro','MarkerFaceColor','r')
    hold on
    plot(steps_HS_TO_StanceLeg_XYZ(~idx,4), steps_HS_TO_StanceLeg_XYZ(~idx,6), 'bo','MarkerFaceColor','b')
    
end
axis equal
hold off

%% Make Head rotation matrices
% SHADOW ONLY

HeadGqw = shadowDataTrimmed.Head_Gqw;
HeadGqx = shadowDataTrimmed.Head_Gqx;
HeadGqy = shadowDataTrimmed.Head_Gqy;
HeadGqz = shadowDataTrimmed.Head_Gqz;


headGlobalQuat_wxyz = normalize(quaternion(HeadGqw, HeadGqx, HeadGqy, HeadGqz));% + quaternion(sqrt(.5), sqrt(.5), 0, 0)); %rotate by 90 degs about X axis to make Z point up

headRotMat_row_col_fr = headGlobalQuat_wxyz.RotationMatrix;


%% find eye positions in Shadow reference frame - BOTH
[ rEyeballCenterXYZ, lEyeballCenterXYZ,worldCamCenterXYZ ] = findEyePositions(headGlobalQuat_wxyz, shadow_fr_mar_dim, shadowMarkerNames,  calibFrame);

%% calc calib mat points - BOTH
[ calibPoint ] = calcCalibPoint( shadow_fr_mar_dim, shadowMarkerNames, calibFrame);


%% find head vecotrs
[ headVecX_fr_xyz, headVecY_fr_xyz, headVecZ_fr_xyz] = findHeadVecs(headGlobalQuat_wxyz, shadow_fr_mar_dim, shadowMarkerNames,  calibFrame, calibPoint, vorFrames);

for ii = 1:length(headGlobalQuat_wxyz)
   
    basisX(ii,:) = quat2rotm(headGlobalQuat_wxyz(ii).e')*[1;0;0];
    basisY(ii,:) = quat2rotm(headGlobalQuat_wxyz(ii).e')*[0;1;0];
    basisZ(ii,:) = quat2rotm(headGlobalQuat_wxyz(ii).e')*[0;0;1];

end

%% calibrate yr eyeballs (and yr WorldCamGaze guy)!
% VOR FRAME METHOD - find camera alignment (i.e. the rotations needed for to make gaze vector align with calibration points during vorFrames)
if useEye(1)
    %right eye first
    vData.dataType                =       1; %1 = right eye, 2 = left eye, 3 = world cam
    vData.calibPoint            = calibPoint;
    vData.eyeballCenterXYZ      = rEyeballCenterXYZ(vorFrames,:);
    
    vData.confidence            = rEye.confidence(vorFrames);
    vData.eye_pupCircCen_x      = rEye.circle_3d_center_x(vorFrames,:);
    vData.eye_pupCircCen_y      = rEye.circle_3d_center_y(vorFrames,:);
    vData.eye_pupCircCen_z      = rEye.circle_3d_center_z(vorFrames,:);
    
    vData.eye_sphCenCam_x     = rEye.sphere_center_x(vorFrames,:);
    vData.eye_sphCenCam_y     = rEye.sphere_center_y(vorFrames,:);
    vData.eye_sphCenCam_z     = rEye.sphere_center_z(vorFrames,:);
    
    vData.headRotMat_row_col_fr        = headRotMat_row_col_fr(:,:,vorFrames);
    vData.shadow_fr_mar_dim     = squeeze(shadow_fr_mar_dim(vorFrames,:,:));
    vData.rHeelXYZ              = rHeelXYZ(vorFrames,:);
    vData.lHeelXYZ              = lHeelXYZ(vorFrames,:);
    vData.shadowMarkerNames     = shadowMarkerNames;
    vData.plotDebug             = true;
    
    vorAlignLossFun = @(camAlignEulerGuess) vorPupilAlignErrFun_eyeCam(vData, camAlignEulerGuess);
    
    initialCamEulerGuess = [0 0 0]; %starting guess for camAlignRotMat
    opts = optimset('Display', 'iter', 'MaxFunEvals',5000, 'PlotFcns',{@optimplotx, @optimplotfval,@optimplotfirstorderopt});
    [camAlignEuler, rVorCalibErr] = fminunc(vorAlignLossFun, initialCamEulerGuess, opts);
    camAlignQuat= quaternion.eulerangles('123',camAlignEuler(1),camAlignEuler(2),camAlignEuler(3));
    rEyeAlignRotMat = camAlignQuat.RotationMatrix;
end

if useEye(2)
    % left eye
    clear vData
    vData.dataType                =       2; %1 = right eye, 2 = left eye, 3 = world cam
    vData.calibPoint            = calibPoint;
    vData.eyeballCenterXYZ      = lEyeballCenterXYZ(vorFrames,:);
    
    vData.confidence            = lEye.confidence(vorFrames);
    vData.eye_pupCircCen_x      = lEye.circle_3d_center_x(vorFrames,:);
    vData.eye_pupCircCen_y      = lEye.circle_3d_center_y(vorFrames,:);
    vData.eye_pupCircCen_z      = lEye.circle_3d_center_z(vorFrames,:);
    
    vData.eye_sphCenCam_x     = lEye.sphere_center_x(vorFrames,:);
    vData.eye_sphCenCam_y     = lEye.sphere_center_y(vorFrames,:);
    vData.eye_sphCenCam_z     = lEye.sphere_center_z(vorFrames,:);
    
    vData.headRotMat_row_col_fr        = headRotMat_row_col_fr(:,:,vorFrames);
    vData.shadow_fr_mar_dim     = squeeze(shadow_fr_mar_dim(vorFrames,:,:));
    vData.rHeelXYZ              = rHeelXYZ(vorFrames,:);
    vData.lHeelXYZ              = lHeelXYZ(vorFrames,:);
    vData.shadowMarkerNames     = shadowMarkerNames;
    vData.plotDebug             = true;
    
    vorAlignLossFun = @(camAlignEulerGuess) vorPupilAlignErrFun_eyeCam(vData, camAlignEulerGuess);
    
    initialCamEulerGuess = [0 0 0]; %starting guess for camAlignRotMat
    opts = optimset('Display', 'iter', 'MaxFunEvals',5000, 'PlotFcns',{@optimplotx, @optimplotfval,@optimplotfirstorderopt});
    [camAlignEuler, lVorCalibErr] = fminunc(vorAlignLossFun, initialCamEulerGuess, opts);
    camAlignQuat= quaternion.eulerangles('123',camAlignEuler(1),camAlignEuler(2),camAlignEuler(3));
    lEyeAlignRotMat = camAlignQuat.RotationMatrix;
end

%%  calc px2mmScale (from 'fixPupilVid.m) - 2018-04-19
% % % % JSM Woodchips
% tapePxX = [849,1063;835,1060;843,1064;855,1067;856,1068;857,1071;862,1073;844,1059;849,1063;566,798;552,782;821,1040;831,1046]
% tapePxY = [815,820;1011,1012;937,941;729.000000000000,735.000000000000;734.000000000000,736;558.000000000000,562.000000000000;555.000000000000,560.000000000000;779,781;768.000000000000,770;793,781;792,779;774,773;766.000000000000,767.000000000000]

% % % JSM Rocks
% tapePxX =[885,1094;874,1084;878,1088;879,1087;882,1092;881,1091;882,1093;884,1092;881,1091;880,1091;882,1093;882,1094;880,1096;876,1096;875,1098];
% tapePxY =[819,803;814,800;813,799;812,797;808,792;810,795;811,797;811,796;808,795;811,799;823,808;881,868;933,919;984,970;1053,1037];
%
%         %[x,y] pixel coordinates of either end of .5m tape at ~vorFrames(1)
%         calibTape1 = [tapePxX(:,1) tapePxY(:,1)];
%         calibTape2 = [tapePxX(:,2) tapePxY(:,2)];
%         scaleDist = sqrt( (calibTape1(:,1)-calibTape2(:,1)).^2 + (calibTape1(:,2)-calibTape2(:,2)).^2); %euclidean distance between pt1 & pt2 (size of 0.5m tape in pixel units)
%
%         scaleSize = 500;
%         px2mmScale = scaleSize/mean(scaleDist);

if ~(resWidth == 1920 && resHeight == 1080)
    disp('Be Aware - The px2mmScale value here assumes 1920x1080 resolution, and this video is some other resolution.')
    keyboard
end

px2mmScale =2.3232;

%% Calculate Gaze Vectors

if useEye(1)
    [calibDist,rGazeXYZ] = calculateGazeVectors(rEye,rEyeAlignRotMat,rEyeballCenterXYZ,...
        calibFrame,calibPoint,headRotMat_row_col_fr);
end

if useEye(2)
    [calibDist,lGazeXYZ] = calculateGazeVectors(lEye,lEyeAlignRotMat,lEyeballCenterXYZ,...
        calibFrame,calibPoint,headRotMat_row_col_fr);
end


%% Calculate Gaze/Ground intersections
if useEye(1)
    disp('rGazeGroundIntersections')
    [ rGazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, rGazeXYZ, rEyeballCenterXYZ );
end

if useEye(2)
    disp('lGazeGroundIntersections')
    [ lGazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, lGazeXYZ, lEyeballCenterXYZ );
end


%% Get saccades (Jon)
if useEye(1)
    [saccFrames] = findSaccades(rad2deg(rEye.phi), rad2deg(rEye.theta));
end

if useEye(2)
    [saccFrames] = findSaccades(rad2deg(lEye.phi), rad2deg(lEye.theta));
end

%% Find blinks - corralignmenterror is dependent on this...

% use eye traces to find blinks
if useEye(1)
    rEye.theta_blink = removeBlinks(rEye.theta);
    rEye.phi_blink = removeBlinks(rEye.phi);
    rEye.blinks = rEye.theta_blink | rEye.phi_blink;
end

if useEye(2)
    lEye.theta_blink = removeBlinks(lEye.theta);
    lEye.phi_blink = removeBlinks(lEye.phi);
    lEye.blinks = lEye.theta_blink | lEye.phi_blink;
end
% correct for alignment error

%%

if useEye(1)
    worldFrameIndex = rEye.index;
else
    worldFrameIndex = lEye.index;
end

%% Save out all the variables
if ~exist(outputPath,'file')
    mkdir(outputPath)
end
disp('Saving out mat file')
save([outputPath filesep takeID '.mat'])
