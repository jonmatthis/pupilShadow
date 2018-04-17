    function [rEyeOffset, lEyeOffset] = calcTemporalOffset(headVecX_fr_xyz, rEye, lEye, wiggleFrames)
%% align yr wiggly bits
%%% Find a temporal offset to align shadow and pupil streams using the

rEyeX = rEye.norm_pos_x(wiggleFrames);
rEyeY = rEye.norm_pos_y(wiggleFrames);

lEyeX = lEye.norm_pos_x(wiggleFrames);
lEyeY = lEye.norm_pos_y(wiggleFrames);

headX = headVecX_fr_xyz(wiggleFrames,1);
headY = headVecX_fr_xyz(wiggleFrames,2);
headZ = headVecX_fr_xyz(wiggleFrames,3);


figure(3244)
 plot(headX,'.-');
 hold on 
 plot(headY,'.-');
 plot(headZ,'.-');
 
 plot(rEyeX(108:end),'ro-')
 plot(rEyeY(108:end),'r.-')
 
 plot(lEyeX(108:end),'bo-')
 plot(lEyeY(108:end),'b.-')
 drawnow
 
 
 minPeakDistance = 20;
 
 [headXpks, headXlocs] = findpeaks(headX);
 [headYpks, headYlocs] = findpeaks(headY);
 [headZpks, headZlocs] = findpeaks(headZ);

 
 [rEyeXpks, rEyeXlocs] = findpeaks(rEyeX);
 [rEyeYpks, rEyeYlocs] = findpeaks(EyeY);
 

