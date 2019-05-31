close all
clear all
restoredefaultpath


fullPath = mfilename('fullpath');

repoPath =  fullPath(1:end-12);
dataPath = 'C:\Users\jon\Dropbox\ResearchProjects\OpticFlowProject\Data';
toolboxPath = 'C:\Users\jon\Dropbox\ResearchProjects\toolboxes';

addpath(genpath(toolboxPath)); %add necessary toolboxes to path
addpath(genpath(repoPath)) %add this repository to the path
%%

%%

numSubs = 1;
numConds = 1;

for subNum = 1;%1:4
    
    switch subNum
        case 1
            sessionID = '2019-05-29_CalgaryWorkshopData';
        case 2
            sessionID = '2019-05-30_CalgaryExtraSkeleton';
    end
    
    
    
    for takeNum = 1;%1:numConds
        
        switch takeNum
            case 1
                takeID = 'Frisbee';
            case 2
                takeID = 'LeafStomping';
            case 3
                takeID = 'RandomPlaytime';
        end
        
        sessionPath = [dataPath filesep sessionID];
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
        
        playLaserSkeleton(w,camParams,dataPath)
    end%numConds
end%numSubs
