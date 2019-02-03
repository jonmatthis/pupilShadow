function [saccFrames] = findSaccades(porX, porY)
% takes in porX and porY (in degrees visual angle) and spits out a vector
% of bools where 'true' values denote frames where subject is saccading.



%% stupid/simple velocity/acceleration based saccade detector 

gVraw = zeros (length(porX),1);
for fr = 2:length(porX)
    gVraw(fr) = sqrt( (porX(fr-1) - porX(fr)).^2 + (porY(fr-1) - porY(fr)).^2);
end

gV = smooth(fillmissing(gVraw,'linear'));
gV  = gV-median(gV)

medAcc = median(abs(diff(gV))); %median acceleration value of gaze
step = medAcc*5; %any step bigger than this is probably a saccade

saccFrames = logical(zeros(size(gV)));
for fr = 2:length(gV)-1
    %%%if the distance between the current value (gaze velocity) and the previous value (and the next value) is
    %%%greater than STEP, label this as a saccade frame
    %%i.e., this is an acceleration threshold
    if abs(gV(fr)-gV(fr-1))>step && abs(gV(fr)-gV(fr+1))>step
        saccFrames(fr) = true;%this is a saccade frame
    else
        saccFrames(fr) = false;%this is not a saccadeFrame
    end
end

saccFramesRaw = saccFrames;
saccFrames =  logical(round(smooth(saccFrames,5))); % my slop-city way to get rid of low-duration fixations and saccades


%%
debug = true;
if debug
-mean(porX)
figure(53298);clf
subplot(2,1,1)
t = 1:length(porX);
plot(t,porX-mean(porX),'.-')
hold on 
plot(t,porY-mean(porY),'.-')
plot(t,porY-mean(porY),'.-')

plot(t(saccFrames), porX(saccFrames)-mean(porX),'ro','MarkerSize',4,'MarkerFaceColor','r')
plot(t(saccFrames), porY(saccFrames)-mean(porY),'ro','MarkerSize',4,'MarkerFaceColor','r')

subplot(2,1,2)
plot(t,gV, '.-')
hold on 
plot(t(saccFrames), gV(saccFrames),'ro','MarkerSize',4,'MarkerFaceColor','r')
end
