close all
clear all
restoredefaultpath

switch getenv('computername')
    case 'MATHISPCWIN10' %Jon's desktop PC
        repoPath = 'D:\Dropbox\ResearchProjects\pupilShadow';
        basePath = 'D:\Dropbox\ResearchProjects\OpticFlowProject\Data';
                addpath(genpath('D:\Dropbox\ResearchProjects\toolboxes')); %add necessary toolboxes to path

    case 'DESKTOP-L29LOMC' %Jon's windows laptop
        repoPath ='C:\Users\jon\Dropbox\ResearchProjects\pupilShadow';
        basePath = 'C:\Users\jon\Dropbox\ResearchProjects\OpticFlowProject\Data';
        addpath(genpath('C:\Users\jon\Dropbox\ResearchProjects\toolboxes')); %add necessary toolboxes to path
end

addpath(genpath(repoPath))
%% 
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
