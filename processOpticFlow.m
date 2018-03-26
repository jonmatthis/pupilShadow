
clearvars
close all

cd E:\Dropbox\UTexas\OpticFlowProject

addpath(genpath(cd))

sub = 'JSM';
cond = 'Woodchips';
framesPath = strcat('E:\OpticFlowFrames',filesep,sub,filesep,cond);

opticFlowObj = opticalFlowFarneback;

load('paramStruct.mat');
camParams = cameraParameters(paramStruct);

%%


recordVid = false;
plotVid = false;
t = [];
flowCell = {};
frameRef = 0;

cd(framesPath)
startFrame = 1;
endFrame = numel(dir('world')) - 2; %play til the end of the video

for ff = startFrame:endFrame
    tic
    
    disp(strcat('PROGRESS:', num2str(ff-startFrame),'-of-',num2str(endFrame-startFrame),...
        '- estimated time remaining ~',num2str((nanmean(t) * (endFrame-ff))/60 ),'mins - Mean Frame Dur ~',num2str(nanmean(t)),...
        '- RecordVid is set to:',num2str(recordVid)))
    
    clf
    f = gcf;
    f.Units = 'normalized';
    f.Position = [1 0.1 1 0.83];
    
    axes('Position',[0 0 1 1])
    
    thisWorldFrameFID = strcat(framesPath,filesep,'world',filesep,'world-',sprintf('%09d',ff),'.png');
    
    thisFrameRaw = imread(thisWorldFrameFID);
    thisFrameRGB = undistortImage(thisFrameRaw, camParams);
    
    thisFrameGray = rgb2gray(thisFrameRGB);
    
    flowRaw = estimateFlow(opticFlowObj, thisFrameGray);
    
    if plotVid
        imshow(thisFrameRGB)
        hold on
        plot(flowRaw, 'DecimationFactor', [10 10],'ScaleFactor',10)
        
        drawnow
    end
    
    cd(strcat(framesPath,filesep,'matFiles'))
    
    
    
    frameName= strcat('flowFrame_',sprintf('%08d',ff));
    save(strcat(frameName,'.mat'),'flowRaw')
    
    t(end+1) = toc;
    
end



