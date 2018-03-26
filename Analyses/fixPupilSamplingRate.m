function [ rEye, lEye ] = fixPupilSamplingRate( pupilExportPath, prefRate, pupilUnixStartTime )
%FIXPUPILSAMPLINGRATE Summary of this function goes here
%  resample PupilLabs data at a constant rate
% path = The datapath of the pupil data
% prefRate = desired sampling rate in Hz

beep

pupTable = readtable(strcat(pupilExportPath,'/pupil_positions.csv'));

varNames = pupTable.Properties.VariableNames;

id = pupTable.id; %Which eyeball does this row reference? 0 = Right eye, 1 = Left eye


origTimestamps_all = pupTable.timestamp;

origTimestamps_eye0 = origTimestamps_all(id==0);
origTimestamps_eye1 = origTimestamps_all(id==1);



skipThese = {'id', 'index', 'method','model_id', 'timestamp'};
for vv = 1:length(varNames)
    thisVarName = varNames{vv};
    thisVar = pupTable.(thisVarName);
    
    if isempty(strmatch(thisVarName, skipThese)) %resample data
        [resampThisVarEye0, r_timestamp] = resample(thisVar(id==0), origTimestamps_eye0, prefRate,'pchip');
        [resampThisVarEye1, l_timestamp] = resample(thisVar(id==1), origTimestamps_eye1, prefRate,'pchip');
        
        rEye.timestamp = r_timestamp(1:min([length(r_timestamp) length(l_timestamp)]));
        lEye.timestamp = l_timestamp(1:min([length(r_timestamp) length(l_timestamp)]));
        
        resampThisVarEye0 =   resampThisVarEye0(1:min([length(resampThisVarEye0) length(resampThisVarEye1)]));
        resampThisVarEye1 =   resampThisVarEye1(1:min([length(resampThisVarEye0) length(resampThisVarEye1)]));
        
        rEye.(thisVarName) = resampThisVarEye0;
        lEye.(thisVarName) = resampThisVarEye1;
        
    elseif ~strcmp(thisVarName,'timestamp')
        
        rEye.(thisVarName) = thisVar(id ==0);
        lEye.(thisVarName) = thisVar(id ==1);
        
    end
end

rEye.unixTimestamp = ((rEye.timestamp - rEye.timestamp(1)) + pupilUnixStartTime);
lEye.unixTimestamp = ((lEye.timestamp - lEye.timestamp(1)) + pupilUnixStartTime);

beep

% %% debug plots
% x = pupTable.norm_pos_x;
% y = pupTable.norm_pos_y;
% 
% figure(22344566)
% subplot(211)
% plot(origTimestamps_eye0, x(id==0),'.-');
% hold on
% plot(rEye.timestamp, rEye.norm_pos_x,'.-')
% 
% subplot(212)
% plot(origTimestamps_eye0, y(id==0),'.-');
% hold on
% plot(rEye.timestamp, rEye.norm_pos_y,'.-')
