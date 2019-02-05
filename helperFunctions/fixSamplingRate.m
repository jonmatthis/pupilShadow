function [ resampledTable ] = fixSamplingRate_( pupTable, startTime,endTime,unixStartTime,prefRate,skipThese )
%FIXPUPILSAMPLINGRATE Summary of this function goes here
%  resample PupilLabs data at a constant rate
% path = The datapath of the pupil data
% prefRate = desired sampling rate in Hz

desTimestamp = (startTime-startTime):1/prefRate:(endTime-startTime);

numFrames = length(desTimestamp)-prefRate; %the desired number of frames after resampling (minus 1 second, to avoid problems arising from minor mismatches in the various streams)

VarNames = pupTable.Properties.VariableNames;

origTimestamps = pupTable.timestamp;

for vv = 1:length(VarNames)
    thisVarName = VarNames{vv};
    thisVar = pupTable.(thisVarName);
    
    if isempty(strmatch(thisVarName, skipThese)) %resample data
        if ~any(isnan(thisVar))

            [resampThisVar, timestamp] = resample(thisVar, origTimestamps, prefRate,'pchip');

            if numFrames > length(timestamp)
                'This is a weird difference between eyes in this recording!!!! hopefully nothing fucks up'
                numFrames = length(timestamp);
            end

                    
            resampThisVar =   resampThisVar(1:numFrames);
            resampledTable.(thisVarName) = resampThisVar;
        end
       
    elseif ~strcmp(thisVarName,'timestamp')
        
        resampledTable.(thisVarName) = thisVar;
        
    end
    
    
    
end

resampledTable.timestamp = timestamp(1:numFrames);
resampledTable.unixTimestamp = ((resampledTable.timestamp - resampledTable.timestamp(1)) + unixStartTime);

