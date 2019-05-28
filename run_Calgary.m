close all
clear all
restoredefaultpath

if isunix
    [ret, name] = system('hostname');
    name = name(1:end-1);
    name = string(name);
elseif ispc || ismac
    name = getenv('computername');
end

switch name
    case 'MATHISPCWIN10' %Jon's desktop PC
        repoPath = 'D:\Dropbox\ResearchProjects\pupilShadow';
        basePath = 'D:\Dropbox\ResearchProjects\OpticFlowProject\Data';
        addpath(genpath('D:\Dropbox\ResearchProjects\toolboxes')); %add necessary toolboxes to path
    case 'DESKTOP-L29LOMC' %Jon's windows laptop
        repoPath ='C:\Users\jon\Dropbox\ResearchProjects\pupilShadow';
        basePath = 'C:\Users\jon\Dropbox\ResearchProjects\OpticFlowProject\Data';
        addpath(genpath('C:\Users\jon\Dropbox\ResearchProjects\toolboxes')); %add necessary toolboxes to path
    case 'karl-G551JW'
        assert(exist('/home/karl/mexopencv', 'dir')==7, 'Laser skeletons require MexOpenCV to function')
        addpath('/home/karl/mexopencv/')
        addpath('/home/karl/mexopencv/opencv_contrib/')
        addpath(genpath('/home/karl/toolboxes/'));
        repoPath = '/home/karl/pupilShadow/';
        basePath = '/home/karl/Dropbox/OpticFlowProject/Data';
        %          basePath = '/media/karl/44CD-7F85/OpticFlowProject/Data';
end

addpath(genpath(repoPath))
%%

%%

numSubs = 5;
numConds = 3;
for subNum = 5;%1:4
    
    switch subNum
        case 1
            sessionID = '2018-01-23_JSM';
        case 2
            sessionID = '2018-01-26_JAC';
        case 3
            sessionID = '2018-01-31_JAW';
        case 4
            sessionID = '2019-03-07_JPL';
        case 5
            sessionID = 'PreCalgaryTest';
            
    end
    
    
    
    for condNum = 3;%1:numConds
        
        switch condNum
            case 1
                takeID = 'Rocks';
            case 2
                takeID = 'Woodchips';
            case 3
                takeID = 'test';
        end
        
        sessionPath = [basePath filesep sessionID];
        useEye = [1;1]; %use both eyes [true; true;] or just one? [false;  true] or [true; false]
        sessionFunction = @loadSessionInfo_Calgary; %name of the loadSessionInfo file we're going to use
        
        
        
        % process data
        sesh = sessionFunction(sessionID,takeID); %'sesh' is a struct with a bunch of information needed for the rest of the code (i.e. what were the calibration/walk frames, etc)
        
        
        
        processData(sessionID,takeID,sessionPath,useEye,sesh);
        
        
        % walks and walk_names are necessary for splitWalks
        walks = sesh.walks;
        walk_names = repmat({takeID},size(walks,1),1);
        walk_names{end+1} = 'vor';
        splitWalks(sessionPath,takeID,walks,walk_names);
        
    end%numConds
end%numSubs
