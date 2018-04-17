function  [thisWalk_fixed] = correctAlignmentError_opt(thisWalk_orig)

thisWalk_fixed = thisWalk_orig;
for iter = 1:2 % 1 = Right EYE, 2 = left eye
    %%
    if iter == 1
        gazeXYZ = thisWalk_orig.rGazeXYZ;
        eyeCenterXYZ = thisWalk_orig.rEyeballCenterXYZ;
        groundFix = thisWalk_orig.rGazeGroundIntersection;
    elseif iter == 2
        gazeXYZ = thisWalk_orig.lGazeXYZ;
        eyeCenterXYZ = thisWalk_orig.lEyeballCenterXYZ;
        groundFix = thisWalk_orig.lGazeGroundIntersection;
    end
    
    
    shadow_fr_mar_dim = thisWalk_orig.shadow_fr_mar_dim;
    shadowMarkerNames = thisWalk_orig.shadowMarkerNames;
    
    comXYZ = squeeze(thisWalk_orig.shadow_fr_mar_dim(:,1,:));
    
    
    %%
    
    groundFix_z = groundFix - comXYZ;
    
    
    [comTheta, comRho] = cart2pol(comXYZ(:,1), comXYZ(:,3));
    [groundFixTheta, groundFixRho] = cart2pol(groundFix_z(:,1), groundFix_z(:,3));
    
    
    
    thetaErr = comTheta-groundFixTheta; %trying to get gaze to be more or less equally distributed around COM path
    
    medianErr =  nanmedian(thetaErr); %median angular difference between gaze and COM <<--rotate everythign by -ThisTheta

    meanErr =  nanmean(thetaErr); %"Switching your 'means' to 'medians' is a sign of desperation" <<-- My own advice, apt in this case (meanErr was underrotating in some cases, I guess due to outliers?)

    figure(483);
    subplot(1,2,iter)
    histogram(thetaErr,'Normalization','probability')
    xlim([-pi pi])
    hold on
    histogram(thetaErr-medianErr,'Normalization','probability')
    hold off

    corrAlignTheta = -medianErr;
    %% get to rotatin'
    
    
    % center g on comXYZ (comXYZ reference frame)
    g_z(:,1) = gazeXYZ(:,1)-comXYZ(:,1);
    g_z(:,2) = gazeXYZ(:,3)-comXYZ(:,3);
    
    %center camXYZ on comXYZ (comXYZ reference frame
    c_z(:,1) = eyeCenterXYZ(:,1)-comXYZ(:,1);
    c_z(:,2) = eyeCenterXYZ(:,3)-comXYZ(:,3);
    
    [gTheta, gRho] = cart2pol(g_z(:,1), g_z(:,2));
    [g_z(:,1), g_z(:,2)] = pol2cart(gTheta-corrAlignTheta, gRho); %rotate by -CorrAlignTheta
    
    [cTheta, cRho] = cart2pol(c_z(:,1), c_z(:,2));
    [c_z(:,1), c_z(:,2)] = pol2cart(cTheta-corrAlignTheta, cRho); %rotate by -CorrAlignTheta
    %
    %
    %
    %
    %
    %
    %         for rr = 1:length(g) %rotate gaze and cam by thetaGuess around their new origin (i.e. the COM)
    %
    %             g(rr,1) = ...
    %                 g(rr,1) * cos(corrAlignTheta)... %x*cos(theta)
    %                 -g(rr,2) * sin(corrAlignTheta);    %-y*sin(theta)
    %
    %             g(rr,2) = ...
    %                 g(rr,1) * sin(corrAlignTheta)+... %x*sin(theta)
    %                 g(rr,2) * cos(corrAlignTheta);    %y*cos(theta)
    %
    %
    %
    %             c(rr,1) = ...
    %                 c(rr,1) * cos(corrAlignTheta)... %x*cos(theta)
    %                 -c(rr,2) * sin(corrAlignTheta);    %-y*sin(theta)
    %
    %             c(rr,2) = ...
    %                 c(rr,1) * sin(corrAlignTheta)... %x*sin(theta)
    %                 +c(rr,2) * cos(corrAlignTheta);    %y*cos(theta)
    %
    %
    %
    %         end
    %
    %
    %
    % revert gaze to right coord system
    gazeXYZ(:,1) = g_z(:,1)+comXYZ(:,1);
    gazeXYZ(:,3) = g_z(:,2)+comXYZ(:,3);
    
    % revert camXYZ to right coord system
    eyeCenterXYZ(:,1) = c_z(:,1)+comXYZ(:,1);
    eyeCenterXYZ(:,3) = c_z(:,2)+comXYZ(:,3);
    
    % recalculate ground fixations
    
    rHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('RightHeel', shadowMarkerNames),:)); % pull out rHeelID marker
    
    rToeXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('RightToe', shadowMarkerNames),:)); % pull out rHeelID marker
    
    lHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('LeftHeel', shadowMarkerNames),:)); % pull out rHeelID marker
    
    lToeXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('LeftToe', shadowMarkerNames),:)); % pull out rHeelID marker000
    
    
    
    if iter == 1
        disp('calckin up some rGazeGroundIntersections')
        [ thisWalk_fixed.rGazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeCenterXYZ );
        thisWalk_fixed.rCorrAlignTheta = corrAlignTheta;
        
        thisWalk_fixed.rGazeXYZ = gazeXYZ;
        thisWalk_fixed.rEyeballCenterXYZ = eyeCenterXYZ;
    elseif iter == 2
        disp('calckin up some lGazeGroundIntersections')
        [ thisWalk_fixed.lGazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeCenterXYZ );
        thisWalk_fixed.lCorrAlignTheta = corrAlignTheta;
        
        thisWalk_fixed.lGazeXYZ = gazeXYZ;
        thisWalk_fixed.lEyeballCenterXYZ = eyeCenterXYZ;
    end
    
    
end

% rotate shadow data
%
% s = nan(size(shadow_fr_mar_dim));
%
% % zero everything (i.e. set origin to comXYZ)
% for ff = 1:length(comXYZ)
%      if mod(ff,1000) == 0; disp(strcat(num2str(ff),' - inCorrAlignErrorWhatver')); end
%
%     for mm = 1:numel(shadow_fr_mar_dim(1,:,1)) %m = Marker
%
%         s(ff,mm,1) = shadow_fr_mar_dim(ff,mm,1) - comXYZ(ff,1);
%         s(ff,mm,2) = shadow_fr_mar_dim(ff,mm,2) - comXYZ(ff,2);
%         s(ff,mm,3) = shadow_fr_mar_dim(ff,mm,3) - comXYZ(ff,3);
%     end
%
% end
%
%
%
%
% %rotate shadow by thetaGuess around their new origin (i.e. the COM)
% for rr = 1:length(comXYZ)
%
%     for mm = 1:numel(shadow_fr_mar_dim(1,:,1)) %m = Marker
%
%         s(rr,mm,1) = ...
%             s(rr,mm,1) * cos(corrAlignTheta)+... %x*cos(theta)
%             s(rr,mm,3) * sin(corrAlignTheta);    %y*sin(theta)
%
%         s(rr,mm,3) = ...
%             -s(rr,mm,1) * sin(corrAlignTheta)+... %x*cos(theta)
%             s(rr,mm,3) * cos(corrAlignTheta);    %y*sin(theta)
%     end
% end
%
% %do, like, the opposite of zero-ing (i.e. put te shadow dat back into inertial frame)
% for ff = 1:length(comXYZ)
%
%
%     for mm = 1:numel(shadow_fr_mar_dim(1,:,1)) %m = Marker
%
%         s(ff,mm,1) = s(ff,mm,1) + comXYZ(ff,1);
%         s(ff,mm,2) = s(ff,mm,2) + comXYZ(ff,2);
%         s(ff,mm,3) = s(ff,mm,3) + comXYZ(ff,3);
%     end
%
% end
%
% thisWalk_orig.shadow_ds_fr_mar_dim = s;


thisWalk_fixed.corrAlignTheta = corrAlignTheta;


%%
% figure(1566); 
% plot(thisWalk_fixed.comXYZ(:,1),thisWalk_fixed.comXYZ(:,3),'k.-')
% hold on 
% plot(thisWalk_fixed.rGazeGroundIntersection(:,1),thisWalk_fixed.rGazeGroundIntersection(:,3),'m.-')
% plot(thisWalk_fixed.lGazeGroundIntersection(:,1),thisWalk_fixed.lGazeGroundIntersection(:,3),'g.-')
% drawnow



