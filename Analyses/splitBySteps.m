function [dataByStepNorm,dataByStepRaw] = splitBySteps(data,steps)
%SPLITBYSTEPS split the input data up according to step events. Send it
%back both as a cell array, and as a matrix of data nomalized to be the
%same length

dataByStepRaw= [];

hs = steps(:,1); %heelstrike frames

for ss = 5:length(hs)-5
    
    thisStepData= data(hs(ss):(hs(ss+1)-1));
    thisStepData = thisStepData(~isnan(thisStepData));
    dataByStepRaw{end+1} = thisStepData;
    
    stepDur(ss) = length(thisStepData);
       
end

meanStepDur = nanmean(stepDur);


newIndex = linspace(1,meanStepDur,round(meanStepDur)); %resample data to be 100 data points across the mean step duration

for nn = 1:length(dataByStepRaw)
   
    thisStepRaw = dataByStepRaw{nn};
    
    origIndex = 1:length(thisStepRaw);
    
    
    thisStepNorm = interp1(origIndex, thisStepRaw, newIndex);
    
    dataByStepNorm(nn,:) = thisStepNorm;
    
end

beep
