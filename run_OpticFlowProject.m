close all
clear all
restoredefaultpath

[ret, name] = system('hostname');
name = name(1:end-1);
name = string(name);

switch name
    case 'MATHISPCWIN10'
        repoPath = 'D:\Dropbox\ResearchProjects\pupilShadow';
        basePath = 'D:\Dropbox\ResearchProjects\OpticFlowProject\Data';
    case 'DESKTOP-L29LOMC'
        repoPath ='C:\Users\jon\Dropbox\ResearchProjects\pupilShadow';
        basePath = 'C:\Users\jon\Dropbox\ResearchProjects\OpticFlowProject\Data';
    case 'karl-G551JW'
        repoPath = '/home/karl/pupilShadow/';
        basePath = '/home/karl/Dropbox/OpticFlowProject/Data';
end

addpath(genpath(repoPath))

numSubs = 3;
numConds = 2;
for subNum = 4:4
    
    switch subNum
        case 1
            sessionID = '2018-01-23_JSM';
        case 2
            sessionID = '2018-01-26_JAC';
        case 3
            sessionID = '2018-01-31_JAW';
        case 4
            sessionID = '2019-02-27_JSM';
    end
    
    
    
    for condNum = 2:numConds
        
        switch condNum
            case 1
                takeID = 'Woodchips';
            case 2
                takeID = 'Rocks';
        end
        
        sessionPath = [basePath filesep sessionID];
        useEye = [1;1];
        sessionFunction = @loadSessionInfo_opticflow;
        
        
        
        % process data
        sesh = sessionFunction(sessionID,takeID);
        
      
        
        processData(sessionID,takeID,sessionPath,useEye,sesh);
        
        
        
        % walks and walk_names are necessary for splitWalks
        walks = sesh.walks;
        walk_names = repmat({takeID},size(walks,1),1);
        % walk_names{1} = 'vor'; %this needs to be done differently
        splitWalks(sessionPath,takeID,walks,walk_names);
        
    end%numConds
end%numSubs
