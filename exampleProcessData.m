%% example 1
clear;

sessionID = '2018-04-05-S03';
takeID = 'Binocular';
sessionPath = '/Volumes/Walking4TB/Data/2018-04-05-S03';
useEye = [1;0];
sessionFunction = @loadSessionInfo;

sesh = sessionFunction(sessionID,takeID);
% process data
processData(sessionID,takeID,sessionPath,useEye,sesh);

% split walks
load([sessionPath filesep 'OutputFiles' filesep takeID],'gaze');
[walks, walk_names] = loadWalks_Berkeley(...
    [sessionPath filesep takeID filesep 'Pupil' filesep 'exports' filesep sesh.pupilExportDir],...
    gaze);

splitWalks(sessionPath,takeID,walks,walk_names)


%% example 2

sessionID = '2018-01-31_JAW';
takeID = 'Rocks';
sessionPath = ['/Volumes/Walking4TB/OpticFlowProject/Data/' sessionID];
useEye = [1;1];
sessionFunction = @loadSessionInfo_opticflow;

% process data
sesh = sessionFunction(sessionID,takeID);
processData(sessionID,takeID,sessionPath,useEye,sesh);


% walks and walk_names are necessary for splitWalks
walks = sesh.walks;
walk_names = repmat({takeID},size(walks,1),1);
walk_names{1} = 'vor';
splitWalks(sessionPath,takeID,walks,walk_names)
