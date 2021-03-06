% load 'D:\Dropbox\ResearchProjects\OpticFlowProject\Data\PreCalgaryTest\OutputFiles\test_allWalks.mat'; %path to the 'allWalks' file
% load 'D:\Dropbox\ResearchProjects\pupilShadow\laserSkeleton\camParam1080_20180626.mat'; %path to the camera intrinsics file

% repoPath = 'D:\Dropbox\ResearchProjects\pupilShadow';
% basePath = 'D:\Dropbox\ResearchProjects\OpticFlowProject\Data'; %data folder that contains the session data folders
% 
% addpath(genpath('D:\Dropbox\ResearchProjects\toolboxes')); %add necessary toolboxes to path

function runLaserSkeleton(repoPath, basePath, camParams)

repoPath ='C:\Users\jon\Dropbox\ResearchProjects\pupilShadow';
basePath = 'C:\Users\jon\Dropbox\ResearchProjects\OpticFlowProject\Data';
addpath(genpath('C:\Users\jon\Dropbox\ResearchProjects\toolboxes')); %add necessary toolboxes to path
addpath(genpath(repoPath))

%%add paths to mexopencv
assert(exist('C:/dev/mexopencv', 'dir')==7, 'Laser skeletons require MexOpenCV to function')
addpath('C:/dev/mexopencv')
addpath('C:/dev/mexopencv/opencv_contrib/')

walkNum = 1; %which walk to play
w = allWalks{walkNum};
camParams = cameraParams1080_2;
playLaserSkeleton(w,camParams,basePath)
end