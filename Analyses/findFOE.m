
function [FOExy, XY, maxKDE] = findFOE(flow,frameRGB,top,left,bottom,right, maxKDE)

debugPlot = false;


numRows = 1;
numCols = 3;

if debugPlot
    figure(99984)
    clf
    subplot(numRows, numCols,1)
    imagesc(frameRGB);hold on
    plot(flow,'DecimationFactor',[40 40],'ScaleFactor',8)
    
    h1 = findobj(gca,'Type','Quiver');
    h1.MaxHeadSize = .001;
    h1.Marker = '.';
    h1.MarkerEdgeColor = 'k';
    
    h1.Color = 'r';
    h1.LineWidth = 1;
    
    axis equal
    axis([0 size(frameRGB,2) 0 size(frameRGB,1)])
end



%% Wind demo
% close all
% load wind
% [sx,sy] = meshgrid(80:130,20:50);
% streamline(stream2(x(:,:,5),y(:,:,5),u(:,:,5),v(:,:,5),sx,sy));


%% Invert flow field
nvx = -flow.Vx;
nvy = -flow.Vy;



meshSize = 100;
[starty,startx] = meshgrid(top+1:meshSize:bottom, left+1:meshSize:right);%  full Image

% numStreams = 8; %the number of streamlines you'll make on each row and column
% [starty,startx] = meshgrid(linspace(top+1, bottom, numStreams), linspace(left+1,right,numStreams));%  full Image


if debugPlot
    subplot(numRows, numCols, 2)
    m1 = plot(startx, starty,'k.');
    hold on
    axis ij
    axis equal
    axis([0 size(frameRGB,2) 0 size(frameRGB,1)])
end
%%

XY = stream2(nvx, nvy,startx,starty, [1 1000]);

% z = zeros(size(starty));

% for ll = 1:length(XY)
%
%     ll
%     if isempty(XY{ll})
%         continue
%     end
%
%     Aw=hist3(XY{ll}(end,:), 'edges', {1:meshSize:size(nvx,1), 1:meshSize:size(nvy,2)})';
%     z = z+Aw;
%     %     imagesc(z)
%     %     colorbar
%     %     drawnow
%
%
%
% end


% [value, location] = max(z(:));
% [R,C] = ind2sub(size(A),location);


if debugPlot
    s1 = streamline(XY);
    
    for ll = 1:length(s1)
        s1(ll).LineWidth = .5;
        %         s1(ll).Color = 'k';
    end
    
end


flowEnds =  [];

lowest_empt=0;
lowest_ind=length(XY)+1;

for bb = 1:length(XY)
    if isempty(XY{bb})&&lowest_empt==0
        lowest_empt=1;
        lowest_ind=bb;
    end
end


for ll = 1:lowest_ind-1;
    
    flowEnds(end+1,:) = XY{ll}(end,:);
    
    if debugPlot
        plot(XY{ll}(end,1)+rand , XY{ll}(end,2)+rand , 'rp')
    end
    
end


%%%%Animate particles
%  str1 = streamparticles(XY,1000,'Animate',50);
%%
% figure(3939)
% strAx = gca; 
% axis ij 
% axis equal
% xlim([0 size(frameRGB,2)])
% ylim([0 size(frameRGB,1)])
% str1 = streamparticlesVideo(imAx, XY,5000, 'Animate',5 )
% beep



%%

%                 XY=[dy{cc,kLag}, dx{cc,kLag}];% NOTE - swapping X for Y to make the plot easier

nflowEnds=size(flowEnds,1);


% our smoothing kernel
xax=[-500:5:500];
[xx,yy]=meshgrid(xax);
sigma=50; % std dev of smoothing kernel (variance of gaussian)
smkern=exp(-.5 * (xx.^2 + yy.^2)/sigma^2); %smoothed with an isotropic gaussian kernel with std dev of sigma
smkern=smkern/sum(smkern(:));

%surf(smkern,'EdgeColor','none'); %check out yer kernerl, if you're into that


n=size(flowEnds,1);


% step 1) 2D histogram of all fixations (eyepositions) in our
% xy grid
Aw=hist3(fliplr(flowEnds), 'edges', {1:size(frameRGB,1), 1:size(frameRGB,2)});
% step 2) smooth it a bit with a kernel that sums to 1 (so we're not changing the vales)
KDEraw=conv2(Aw, smkern, 'same');

% KDE=KDEraw/numel(startx); %normalize by the number of starting particles, if you're into that
KDE=KDEraw;

%find location of peak of prob dist
maxVal = max(KDE(:));
[FOEy, FOEx] = find(KDE == maxVal);
FOEy=FOEy(1);
FOEx=FOEx(1);

    FOExy = [FOEx, FOEy];


if debugPlot
    subplot(numRows, numCols, 3)
    %
    % s1 = surf(KDE);
    % s1.EdgeColor = 'none';
    % colorbar
    imagesc(KDE);
    hold on
    
    caxis([ 0 .5]) %set colorbar to show range from 0 to 1 (rather than the default range of min(KDE(:)) to max(KDE(:)))
    
    colorbar
    axis equal
       axis([0 size(frameRGB,2) 0 size(frameRGB,1)])
    
    subplot(numRows, numCols, 1) %%plot FOE on Subplot 1
    plot3(FOEx, FOEy, maxVal,'kp','MarkerFaceColor', 'm', 'MarkerSize' ,24)

    drawnow
    
end


maxKDE(end+1) = max(KDE(:));

% figure(3434)
% subplot(211)
% imagesc(KDE);
% hold on

% caxis([ 0 numel(startx)/2]) %set colorbar to show range from 0 to 1 (rather than the default range of min(KDE(:)) to max(KDE(:)))

% subplot(212)
% plot(maxKDE,'-o')
% hold on
% drawnow


%%
end


%     %%
%     figure
%     xlim([0 640])
%     ylim([0 480])
% %     imagesc(frameRGB);
%     hold on
%
%     o1 = imagesc(thisFlow.Orientation);
%     hold on
% %     o1.AlphaData = .3;
%     colormap hsv
%
%     plot(thisFlow,'DecimationFactor',[10 10],'ScaleFactor',4)%, 'MarkerEdgeColor','r','ShowArrowHead','off','LineWidth',2);
%
%     h1 = findobj(gca,'Type','Quiver');
%     h1.MaxHeadSize = .3;
%     h1.Marker = '.';
%     h1.MarkerEdgeColor = 'k';
%
%
%
%

