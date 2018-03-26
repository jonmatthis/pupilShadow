function [ shadowDataResamp] = fixShadowSamplingRate( shadowDataRaw, prefRate )
%FIXPUPILSAMPLINGRATE Summary of this function goes here
%  resample PupilLabs data at a constant rate
% path = The datapath of the pupil data
% prefRate = desired sampling rate in Hz

og_time = shadowDataRaw.og_time;

newTimestamp = [0:1/prefRate:og_time(end)]';

fNames = fieldnames(shadowDataRaw);

skipThese = {'markerData','og_time'}; %resample all the non marker data first
for vv = 1:length(fNames)
    thisVarName = fNames{vv};
    thisVar = shadowDataRaw.(thisVarName);
    
    if isempty(strmatch(thisVarName, skipThese))
        [resampThisVar] = interp1(og_time, thisVar, newTimestamp);

        shadowDataResamp.(thisVarName) = resampThisVar;
    end
end

shadowDataResamp.timestamp = newTimestamp;%get new timestamps


% %%%debug
% figure(43958)
% plot(shadowDataRaw.og_time, shadowDataRaw.HeadGqw,'.-');hold on 
% plot(shadowDataRaw.og_time, shadowDataRaw.HeadGqx,'.-')
% plot(shadowDataRaw.og_time, shadowDataRaw.HeadGqy,'.-')
% plot(shadowDataRaw.og_time, shadowDataRaw.HeadGqz,'.-')
% 
% 
% plot(shadowDataResamp.timestamp, shadowDataResamp.HeadGqw,'o-','MarkerSize',3);hold on 
% plot(shadowDataResamp.timestamp, shadowDataResamp.HeadGqx,'o-','MarkerSize',3)
% plot(shadowDataResamp.timestamp, shadowDataResamp.HeadGqy,'o-','MarkerSize',3)
% plot(shadowDataResamp.timestamp, shadowDataResamp.HeadGqz,'o-','MarkerSize',3)
% 

c3d_fr_mar_dimRaw = shadowDataRaw.markerData;
c3d_fr_mar_dimResamp = [];

for mm = 1:length(squeeze(c3d_fr_mar_dimRaw(1,:,1))) % resample the marker data now
    thisMar_xyz = squeeze(c3d_fr_mar_dimRaw(:,mm,:));
    
    for dd = 1:length(squeeze(c3d_fr_mar_dimRaw(1,mm,:))) % loop through each dimension
%         [resampThisMar_xyz(:,dd), newTimestamp] = resample(thisMar_xyz(:,dd), og_time, prefRate);
        [resampThisMar_xyz(:,dd)] = interp1(og_time, thisMar_xyz(:,dd), newTimestamp);

    end
    c3d_fr_mar_dimResamp(:,mm,:) = resampThisMar_xyz;
end

shadowDataResamp.markerData = c3d_fr_mar_dimResamp;
% %%%%%debug
% figure(4395)
% plot(shadowDataRaw.og_time, squeeze(c3d_fr_mar_dimRaw(:,1,:)),'.-');hold on 
% plot(shadowDataResamp.timestamp, squeeze(c3d_fr_mar_dimResamp(:,1,:)),'o-','MarkerSize',3);
% hold off



