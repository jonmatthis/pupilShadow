function [ rEye, lEye, gaze ] = fixPupilSamplingRate( pupilExportPath, prefRate, pupilUnixStartTime )
%FIXPUPILSAMPLINGRATE Summary of this function goes here
%  resample PupilLabs data at a constant rate
% path = The datapath of the pupil data
% prefRate = desired sampling rate in Hz


pupTable = readtable(strcat(pupilExportPath,'/pupil_positions.csv'));
gazeTable = readtable(strcat(pupilExportPath,'/gaze_positions.csv'));

startTime = max([pupTable.timestamp(1) gazeTable.timestamp(1)]);
endTime = min([pupTable.timestamp(end) gazeTable.timestamp(end)]);

desTimestamp = (startTime-startTime):1/prefRate:(endTime-startTime);

numFrames = length(desTimestamp)-prefRate; %the desired number of frames after resampling (minus 1 second, to avoid problems arising from minor mismatches in the various streams)

%% Do the thing for the Pupil Positions Data (info from the left (1) and right (0) eye camera screens) 
pupVarNames = pupTable.Properties.VariableNames;

id = pupTable.id; %Which eyeball does this row reference? 0 = Right eye, 1 = Left eye


origTimestamps_all = pupTable.timestamp;

origTimestamps_eye0 = origTimestamps_all(id==0);
origTimestamps_eye1 = origTimestamps_all(id==1);


rEye.index_orig = pupTable.index(id==0);
lEye.index_orig = pupTable.index(id==1);

skipThese = {'id','index' 'method','model_id', 'timestamp'};
for vv = 1:length(pupVarNames)
    thisVarName = pupVarNames{vv};
    thisVar = pupTable.(thisVarName);
    
    if isempty(strmatch(thisVarName, skipThese)) %resample data
        [resampThisVarEye0, r_timestamp] = resample(thisVar(id==0), origTimestamps_eye0, prefRate,'pchip');
        [resampThisVarEye1, l_timestamp] = resample(thisVar(id==1), origTimestamps_eye1, prefRate,'pchip');
        
        rEye.timestamp = r_timestamp(1:numFrames);
        lEye.timestamp = l_timestamp(1:numFrames);
        
        resampThisVarEye0 =   resampThisVarEye0(1:numFrames);
        resampThisVarEye1 =   resampThisVarEye1(1:numFrames);
        
        rEye.(thisVarName) = resampThisVarEye0;
        lEye.(thisVarName) = resampThisVarEye1;
    
       
    elseif ~strcmp(thisVarName,'timestamp')
        
        rEye.(thisVarName) = thisVar(id ==0);
        lEye.(thisVarName) = thisVar(id ==1);
        
    end
end

rEye.unixTimestamp = ((rEye.timestamp - rEye.timestamp(1)) + pupilUnixStartTime);
lEye.unixTimestamp = ((lEye.timestamp - lEye.timestamp(1)) + pupilUnixStartTime);


%% Do the thing for the Gaze Data (Crosshairs on the world camera screen)

gazeVarNames = gazeTable.Properties.VariableNames;


origTimestamps_gaze = gazeTable.timestamp;


gaze.index_orig = gazeTable.index;

skipThese = {'base_data','index','timestamp'};
for vv = 1:length(gazeVarNames)
    thisVarName = gazeVarNames{vv};
    thisVar = gazeTable.(thisVarName);
    
    if isempty(strmatch(thisVarName, skipThese)) %resample data
        [resampThisVarGaze, g_timestamp] = resample(thisVar, origTimestamps_gaze, prefRate,'pchip');
        
        gaze.timestamp = g_timestamp(1:numFrames);
        
        resampThisVarGaze =   resampThisVarGaze(1:numFrames);
        
        gaze.(thisVarName) = resampThisVarGaze;    
       
    elseif ~strcmp(thisVarName,'timestamp')
        
        gaze.(thisVarName) = thisVar;
        
    end
end

gaze.unixTimestamp = ((gaze.timestamp - gaze.timestamp(1)) + pupilUnixStartTime);

