function [corrAlignErr] = correctAlignmentErrorFcn(gazeXYZ, comXYZ, camXYZ, rHeelXYZ, lHeelXYZ, frames, plotDebug, thetaGuess)

% gaze originating from camXYZ (inertial reference frame)
gaze_guess(:,1) = gazeXYZ(:,1)+camXYZ(:,1);
gaze_guess(:,2) = gazeXYZ(:,2)+camXYZ(:,2);
gaze_guess(:,3) = gazeXYZ(:,3)+camXYZ(:,3);

% center g on comXYZ (comXYZ reference frame)
gaze_guess(:,1) = gaze_guess(:,1)-comXYZ(:,1);
gaze_guess(:,2) = gaze_guess(:,2)-comXYZ(:,2);
gaze_guess(:,3) = gaze_guess(:,3)-comXYZ(:,3);

%center camXYZ on comXYZ (comXYZ reference frame
cam(:,1) = camXYZ(:,1)-comXYZ(:,1);
cam(:,2) = camXYZ(:,2)-comXYZ(:,2);
cam(:,3) = camXYZ(:,3)-comXYZ(:,3);





for rr = 1:length(gaze_guess) %rotate gaze and cam by thetaGuess around their new origin (i.e. the COM)
    
    gaze_guess(rr,1) = ...
        gaze_guess(rr,1) * cos(thetaGuess)+... %x*cos(theta)
        gaze_guess(rr,3) * sin(thetaGuess);    %y*sin(theta)
    
    gaze_guess(rr,3) = ...
        -gaze_guess(rr,1) * sin(thetaGuess)+... %x*cos(theta)
        gaze_guess(rr,3) * cos(thetaGuess);    %y*sin(theta)
    
    cam(rr,1) = ...
        cam(rr,1) * cos(thetaGuess)+... %x*cos(theta)
        cam(rr,3) * sin(thetaGuess);    %y*sin(theta)
    
    cam(rr,3) = ...
        -cam(rr,1) * sin(thetaGuess)+... %x*cos(theta)
        cam(rr,3) * cos(thetaGuess);    %y*sin(theta)
    
end

%%put 'em back into an inertial reference frame
gaze_guess(:,1) = gaze_guess(:,1)+comXYZ(:,1);
gaze_guess(:,2) = gaze_guess(:,2)+comXYZ(:,2);
gaze_guess(:,3) = gaze_guess(:,3)+comXYZ(:,3);

cam(:,1) = cam(:,1)+comXYZ(:,1);
cam(:,2) = cam(:,2)+comXYZ(:,2);
cam(:,3) = cam(:,3)+comXYZ(:,3);

[ gazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gaze_guess, cam );

if plotDebug
    figure(342332); clf
    
    
    plot(comXYZ(:,1), comXYZ(:,3))
    hold on
    
    plot(cam(:,1), cam(:,3))
    
    
    plot(gazeGroundIntersection(:,1),gazeGroundIntersection(:,3),'.')
    
    axis([min(comXYZ(:,1))-1000 max(comXYZ(:,1))+1000 min(comXYZ(:,3))-1000 max(comXYZ(:,3))+1000  ])

   axis equal 
end



[gazeTheta rho] = cart2pol(gazeXYZ(:,1), gazeXYZ(:,3));

figure(483)
histogram(gazeTheta,'Normalization','probability')
xlim([-.5 .5])

corrAlignErr = nanmean(gazeTheta);

% 
% parfor ee = frames
%     
% % for ee = frames
%     thisDist = nan(length(gazeGroundIntersection),1);
%     for cc = ee:ee+300
%         
%         if cc > length(comXYZ)
%             break
%         end
%         
%         thisDist(end+1) = pdist([gazeGroundIntersection(ee,[1 3]); comXYZ(cc,[1 3])]);
%         
%     end
%     
%     err(ee) = min(thisDist);
%     
% end

% figure(324324)
% plot(err,'-o')
% drawnow
% corrAlignErr = nansum(err)/length(comXYZ);







