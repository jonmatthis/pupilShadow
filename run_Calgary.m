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

numSubs = 1;
numConds = 1;

for subNum = 1;%1:4
    
    switch subNum
        case 1
            sessionID = '2019-05-29_CalgaryWorkshopData';
    end
    
    
    
    for condNum = 1;%1:numConds
        
        switch condNum
            case 1
                takeID = 'Frisbee';
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
        %         walk_names{end+1} = 'vor';
        splitWalks(sessionPath,takeID,walks,walk_names);
        
        %% LASER SKELETON!!!!
        
        %%add paths to mexopencv
        assert(exist('C:/dev/mexopencv', 'dir')==7, 'Laser skeletons require MexOpenCV to function :`(')
        addpath('C:/dev/mexopencv')
        addpath('C:/dev/mexopencv/opencv_contrib/')
        
        %load camera parameters
        camParamPath = [repoPath filesep 'laserSkeletonCode' filesep 'camParam1080_20180626.mat'];
        load(camParamPath);
        camParams = cameraParams1080_2;
        
        %load allWalks cell array (from 'splitWalks.m')
        allWalksPath =[sessionPath filesep 'OutputFiles' filesep takeID '_allWalks.mat'];
        load(allWalksPath)
        
        walkNum = 2; %which walk to play
        w = allWalks{walkNum}; %'w' is a struct that contains all data relevant to walk#walkNum
        
        playLaserSkeleton(w,camParams,basePath)
    end%numConds
end%numSubs
