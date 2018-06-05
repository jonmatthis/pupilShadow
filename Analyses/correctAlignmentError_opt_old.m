function  [thisWalk] = correctAlignmentError_opt(thisWalk_orig)


for iter = 1:2 % 1 = Right EYE, 2 = left eye
    %%
    if iter == 1
        gazeXYZ = thisWalk_orig.rGazeXYZ;
        eyeCenterXYZ = thisWalk_orig.rEyeballCenterXYZ;
        
    elseif iter == 2
        gazeXYZ = thisWalk_orig.lGazeXYZ;
        eyeCenterXYZ = thisWalk_orig.lEyeballCenterXYZ;
        
    end
    
    walks = thisWalk_orig.walks;
    
    shadow_fr_mar_dim = thisWalk_orig.shadow_fr_mar_dim;
    shadowMarkerNames = thisWalk_orig.shadowMarkerNames;
    
    comXYZ = squeeze(thisWalk_orig.shadow_fr_mar_dim(:,1,:));
    
    
    %%
    
    frames = 1:length(comXYZ);%:(walks(1,1):walks(1,2),:));
    
    rHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('RightHeel', shadowMarkerNames),:)); % pull out lHeel marker
    lHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('LeftHeel', shadowMarkerNames),:)); % pull out lHeel marker
    
    6jyr
    %%
    
    if walks(end) > length(gazeXYZ) %this function was called by "split walks", so use the whole datastream
        com = comXYZ;
        
        [comTheta, comRho] = cart2pol(com(:,1), com(:,3));
        [gazeTheta gazeRho] = cart2pol(gazeXYZ(:,1), gazeXYZ(:,3));
        
    elseif walks(end) <= length(gazeXYZ) %this function was called by "processData," so only use a subset of the data (will be refined later in SplitWalks)
        com = comXYZ(walks(1,1):walks(1,2),:) - comXYZ(walks(1,1),:) ;
        
        [comTheta, comRho] = cart2pol(com(:,1), com(:,3));
        [gazeTheta gazeRho] = cart2pol(gazeXYZ(walks(1,1):walks(1,2),1), gazeXYZ(walks(1,1):walks(1,2),3));
    end

thetaErr = comTheta-gazeTheta; %trying to get gaze to be more or less equally distributed around COM path

modeErr =  mode(round(thetaErr,3));

figure(483);clf
histogram(thetaErr,'Normalization','probability')
xlim([-pi pi])
hold on
histogram(thetaErr-modeErr,'Normalization','probability')

% corrAlignTheta = -modeErr;
corrAlignTheta = 0;
    
    %% get to rotatin'
    
    
    % gaze originating from camXYZ (inerial reference frame)
    g(:,1) = gazeXYZ(:,1)+eyeCenterXYZ(:,1);
    g(:,2) = gazeXYZ(:,3)+eyeCenterXYZ(:,3);
    
    
    % center g on comXYZ (comXYZ reference frame)
    g(:,1) = g(:,1)-comXYZ(:,1);
    g(:,2) = g(:,2)-comXYZ(:,3);
    
    %center camXYZ on comXYZ (comXYZ reference frame
    c(:,1) = eyeCenterXYZ(:,1)-comXYZ(:,1);
    c(:,2) = eyeCenterXYZ(:,3)-comXYZ(:,3);
    
    
    
%     
%     for rr = 1:length(g) %rotate gaze and cam by thetaGuess around their new origin (i.e. the COM)
%         
%         g(rr,1) = ...
%             g(rr,1) * cos(corrAlignTheta)... %x*cos(theta)
%             -g(rr,2) * sin(corrAlignTheta);    %-y*sin(theta)
%         
%         g(rr,2) = ...
%             g(rr,1) * sin(corrAlignTheta)+... %x*sin(theta)
%             g(rr,2) * cos(corrAlignTheta);    %y*cos(theta)
%         
%         
%         
%         c(rr,1) = ...
%             c(rr,1) * cos(corrAlignTheta)... %x*cos(theta)
%             -c(rr,2) * sin(corrAlignTheta);    %-y*sin(theta)
%         
%         c(rr,2) = ...
%             c(rr,1) * sin(corrAlignTheta)... %x*sin(theta)
%             +c(rr,2) * cos(corrAlignTheta);    %y*cos(theta)
%         
%         
%         
%     end
%     
    
    
    % revert gaze to right coord system
    g(:,1) = g(:,1)+comXYZ(:,1);
    g(:,2) = g(:,2)+comXYZ(:,3);
    
    g(:,1) = g(:,1)-eyeCenterXYZ(:,1);
    g(:,2) = g(:,2)-eyeCenterXYZ(:,3);
    
    gazeXYZ(:,1) = g(:,1);
    gazeXYZ(:,3) = g(:,2);
    
    thisWalk_orig.gazeXYZ = gazeXYZ;
    
    % revert camXYZ to right coord system
    c(:,1) = c(:,1)+comXYZ(:,1);
    c(:,2) = c(:,2)+comXYZ(:,3);
    
    eyeCenterXYZ(:,1) = c(:,1);
    eyeCenterXYZ(:,3) = c(:,2);
    
    
    if iter == 1
        thisWalk_orig.rGazeXYZ = gazeXYZ;
        thisWalk_orig.rEyeballCenterXYZ = eyeCenterXYZ;
        
    elseif iter == 2
        
        thisWalk_orig.lGazeXYZ = gazeXYZ;
        thisWalk_orig.lEyeballCenterXYZ = eyeCenterXYZ;
    end
    
    
    
    % recalculate ground fixations
    
    rHeelID = find(strcmp('RightHeel', shadowMarkerNames));
    rHeelXYZ = squeeze(shadow_fr_mar_dim(:,rHeelID,:)); % pull out rHeelID marker
    
    rToeID = find(strcmp('RightToe', shadowMarkerNames));
    rToeXYZ = squeeze(shadow_fr_mar_dim(:,rToeID,:)); % pull out rHeelID marker
    
    lHeelID = find(strcmp('LeftHeel', shadowMarkerNames));
    lHeelXYZ = squeeze(shadow_fr_mar_dim(:,lHeelID,:)); % pull out rHeelID marker
    
    lToeID = find(strcmp('LeftToe', shadowMarkerNames));
    lToeXYZ = squeeze(shadow_fr_mar_dim(:,lToeID,:)); % pull out rHeelID marker000
    
    
    
    if iter == 1
        disp('calckin up some rGazeGroundIntersections')
        [ thisWalk_orig.rGazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeCenterXYZ );
        thisWalk_orig.rCorrAlignTheta = corrAlignTheta;
    elseif iter == 2
        disp('calckin up some lGazeGroundIntersections')
        [ thisWalk_orig.lGazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeCenterXYZ );
        thisWalk_orig.lCorrAlignTheta = corrAlignTheta;
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


thisWalk_orig.corrAlignTheta = corrAlignTheta;
thisWalk = thisWalk_orig;




