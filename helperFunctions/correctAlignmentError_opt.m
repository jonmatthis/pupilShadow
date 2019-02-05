function  [thisWalk_fixed] = correctAlignmentError_opt_(thisWalk_orig)

thisWalk_fixed = thisWalk_orig;
figure(8543);clf

for optIter = find(thisWalk_orig.eyes(:))' % 1 = Right EYE, 2 = left eye
    corrAlignLossFun = @(corrAlignGuess) correctAlignmentError_errFun(thisWalk_orig, optIter, corrAlignGuess);
    opts = optimset('Display', 'iter', 'MaxFunEvals',5000);%, 'PlotFcns',{@optimplotx, @optimplotfval});
    %     [corrAlignTheta] = fminunc(corrAlignLossFun, initialCorrAlignGuess, opts);
    guessBound = pi/2;
    if optIter == 1
        [rCorrAlignTheta, rErr] = fminbnd(corrAlignLossFun, -guessBound, guessBound, opts);
    elseif optIter == 2
        [lCorrAlignTheta, lErr] = fminbnd(corrAlignLossFun, -guessBound, guessBound, opts);
    end
end


if (thisWalk_orig.eyes(1) && ~thisWalk_orig.eyes(2))
    corrAlignTheta = rCorrAlignTheta;
    'left'
elseif thisWalk_orig.eyes(2) && ~thisWalk_orig.eyes(1)
    corrAlignTheta = lCorrAlignTheta;
    'right'
elseif sum(thisWalk_orig.lEye_blinks) >= sum(thisWalk_orig.rEye_blinks)
    corrAlignTheta = rCorrAlignTheta;
    'left'
elseif sum(thisWalk_orig.rEye_blinks) >= sum(thisWalk_orig.lEye_blinks)
    corrAlignTheta = lCorrAlignTheta;
    'right'
else
    'problem'
end



%%

for iter = find(thisWalk_orig.eyes(:))'
    
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
        thisWalk_fixed.rCorrAlignTheta = rCorrAlignTheta;
        
        thisWalk_fixed.rGazeXYZ = gazeXYZ;
        thisWalk_fixed.rEyeballCenterXYZ = eyeCenterXYZ;
        
        newGroundFix = thisWalk_fixed.rGazeGroundIntersection;
        origGroundFix = thisWalk_orig.rGazeGroundIntersection;
        
        
    elseif iter == 2
        disp('calckin up some lGazeGroundIntersections')
        [ thisWalk_fixed.lGazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeCenterXYZ );
        thisWalk_fixed.lCorrAlignTheta = lCorrAlignTheta;
        
        thisWalk_fixed.lGazeXYZ = gazeXYZ;
        thisWalk_fixed.lEyeballCenterXYZ = eyeCenterXYZ;
        
        
        newGroundFix = thisWalk_fixed.lGazeGroundIntersection;
        origGroundFix = thisWalk_orig.lGazeGroundIntersection;
        
    end
    
    
    
    
    
    newGroundFix_z  = newGroundFix - comXYZ;
    origGroundFix_z = origGroundFix - comXYZ;
    
    
    [comTheta, comRho] = cart2pol(comXYZ(:,1), comXYZ(:,3));
    [newGroundFixTheta, newGroundFixRho] = cart2pol(newGroundFix_z(:,1), newGroundFix_z(:,3));
    [origGroundFixTheta, origGroundFixRho] = cart2pol(origGroundFix_z(:,1), origGroundFix_z(:,3));
    
    newThetaErr = comTheta-newGroundFixTheta; %trying to get gaze to be more or less equally distributed around COM path
    origThetaErr = comTheta-origGroundFixTheta; %trying to get gaze to be more or less equally distributed around COM path
    
    
    subplot(1,2,iter)
    if sum(isnan(origThetaErr))< numel(origThetaErr)
        polarhistogram(origThetaErr,40,'Normalization','probability','DisplayName','comTheta-groundFixTheta (original)')
    end
    hold on
    if sum(isnan(newThetaErr))< numel(newThetaErr)
        polarhistogram(newThetaErr,40,'Normalization','probability','DisplayName','comTheta-groundFixTheta ("fixed")')
    end
    l  =   legend;
    %     l.Location = 'south';
    
    if iter == 1
        title('Right Eye gaze alignment correction')
    elseif iter == 2
        title('Left Eye gaze alignment correction')
    end
    hold off
    
    drawnow
    
    
    
end
thisWalk_fixed.corrAlignTheta = corrAlignTheta;
end






%% %%%%%%%%%%%%%
%%%%%% error function for the alignment error correction estimation - Tries
%%%%%% to minimize the distance between ("Good," i.e. not implicated by the badDots bools) gaze ponts and the com path
%%%%%%

function [err] = correctAlignmentError_errFun(thisWalk_orig, iter, corrAlignGuess)

corrAlignTheta = corrAlignGuess;

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

% revert gaze to right coord system
gazeXYZ(:,1) = g_z(:,1)+comXYZ(:,1);
gazeXYZ(:,3) = g_z(:,2)+comXYZ(:,3);

% revert camXYZ to right coord system
eyeCenterXYZ(:,1) = c_z(:,1)+comXYZ(:,1);
eyeCenterXYZ(:,3) = c_z(:,2)+comXYZ(:,3);


%%%%%%%%%%%%%%%%%%%%%%%%
%%% evaluate how good this alignment is by measuring mean distance of
%%% gaze from COM path
%%%%%%%%%%%%%%%%%%%%%%%%

% recalculate ground fixations

rHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('RightHeel', shadowMarkerNames),:)); % pull out rHeelID marker

rToeXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('RightToe', shadowMarkerNames),:)); % pull out rHeelID marker

lHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('LeftHeel', shadowMarkerNames),:)); % pull out rHeelID marker

lToeXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('LeftToe', shadowMarkerNames),:)); % pull out rHeelID marker000



if iter == 1
    [ gazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeCenterXYZ );
    
elseif iter == 2
    [gazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeCenterXYZ );
end



groundFix_z = gazeGroundIntersection - comXYZ;


[groundFixTheta, groundFixRho] = cart2pol(groundFix_z(:,1), groundFix_z(:,3));

badDots = abs(groundFixTheta)>pi/1.5; %remove frames where gaze/ground intersection is more than +/- (pi/8) from straight ahead (assumes the misalignment is small, which is usually true?)
badDots = badDots | groundFixRho > 5e3; %remove frames whre gaze/ground intersection is more than 5e4 mm away from gaze point
badDots = badDots | isnan(groundFixTheta);

%     %%%^%%%debug dem bad ol' dots
%     polarplot(groundFixTheta, groundFixRho,'b.','DisplayName','all gaze/ground dots')
%     rlim([0 2e4])
%     hold on
%     polarplot(groundFixTheta(badDots), groundFixRho(badDots),'r.','DisplayName','bad gaze/ground dots :(')
%     l = legend;
%     l.Location = 'south';
%




%%%%remove bad gaze/ground dots!
gazeGroundIntersection(badDots,:) = [];
comXYZ(badDots,:) = [];

[~,dists] = knnsearch(gazeGroundIntersection, comXYZ);

err = nanmean(dists);



end

