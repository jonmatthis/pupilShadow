function  [thisWalk_fixed] = correctAlignmentError_lookUpTable(thisWalk_orig)

thisWalk_fixed = thisWalk_orig;
figure(8543);clf

corrAlignTheta = corrAlignLookup(thisWalk_orig.sessionID, thisWalk_orig.takeID, thisWalk_orig.ww); %this'll return a NaN if you haven't set the proper value for this walk yet



stillLooking = true; %have you found your favorite correction value yet?
firstRun = true;
%%
while stillLooking
    
    if ~isnan(corrAlignTheta) && firstRun %if it's the first run of the loop and corrAlignTheta isn't a Nan - Congrats! You're good to go! Hopefully!
        stillLooking = false;
    end
    
    firstRun = false;
    
    if isnan(corrAlignTheta) %if the lookup table returned a NaN, set corrAlignTheta to 0 and start the manual correction loop
        corrAlignTheta = 0;
    end
    
    for iter = 1:2
        
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
            %             disp('calckin up some rGazeGroundIntersections')
            [ thisWalk_fixed.rGazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeCenterXYZ );
            %             thisWalk_fixed.rCorrAlignTheta = rCorrAlignTheta;
            
            thisWalk_fixed.rGazeXYZ = gazeXYZ;
            thisWalk_fixed.rEyeballCenterXYZ = eyeCenterXYZ;
            
            if strcmp(thisWalk_orig.trialType, 'Fix') % For 'Distant Fixation' trials, use GazeXYZ instead of GazeGroundIntersections, becuase sub's gaze vectors rarely intersect groundplane in that condition
                newGroundFix = thisWalk_fixed.rGazeXYZ;
                origGroundFix = thisWalk_orig.rGazeXYZ;
            else
                newGroundFix = thisWalk_fixed.rGazeGroundIntersection;
                origGroundFix = thisWalk_orig.rGazeGroundIntersection;
            end
            
        elseif iter == 2
            %             disp('calckin up some lGazeGroundIntersections')
            [ thisWalk_fixed.lGazeGroundIntersection] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeCenterXYZ );
            %             thisWalk_fixed.lCorrAlignTheta = lCorrAlignTheta;
            
            thisWalk_fixed.lGazeXYZ = gazeXYZ;
            thisWalk_fixed.lEyeballCenterXYZ = eyeCenterXYZ;
            
            if strcmp(thisWalk_orig.trialType, 'Fix') % For 'Distant Fixation' trials, use GazeXYZ instead of GazeGroundIntersections, becuase sub's gaze vectors rarely intersect groundplane in that condition
                newGroundFix = thisWalk_fixed.lGazeXYZ;
                origGroundFix = thisWalk_orig.lGazeXYZ;
            else
                newGroundFix = thisWalk_fixed.lGazeGroundIntersection;
                origGroundFix = thisWalk_orig.lGazeGroundIntersection;
            end
            
        end
        
        
        
        
        
        newGroundFix_z  = newGroundFix - comXYZ;
        origGroundFix_z = origGroundFix - comXYZ;
        
        
        %         [comTheta, comRho] = cart2pol(comXYZ(:,1), comXYZ(:,3));
        %         [newGroundFixTheta, newGroundFixRho] = cart2pol(newGroundFix_z(:,1), newGroundFix_z(:,3));
        %         [origGroundFixTheta, origGroundFixRho] = cart2pol(origGroundFix_z(:,1), origGroundFix_z(:,3));
        %
        %
        %
        %         newThetaErr = comTheta-newGroundFixTheta; %trying to get gaze to be more or less equally distributed around COM path
        %         origThetaErr = comTheta-origGroundFixTheta; %trying to get gaze to be more or less equally distributed around COM path
        %
        %
        %         subplot(1,2,iter)
        %         polarhistogram(origThetaErr,'Normalization','probability','DisplayName','comTheta-groundFixTheta (original)')
        %         hold on
        %         polarhistogram(newThetaErr,'Normalization','probability','DisplayName','comTheta-groundFixTheta ("fixed")')
        %         l  =   legend;
        %         l.Location = 'south';
        %
        %         if iter == 1
        %             title('Right Eye gaze alignment correction')
        %         elseif iter == 2
        %             title('Left Eye gaze alignment correction')
        %         end
        %         hold off
        if iter == 2
            clf
            subplot(2,1,1)
            
            plot(thisWalk_orig.comXYZ(:,1), thisWalk_orig.comXYZ(:,3),'m')
            
            hold on
            plot(thisWalk_orig.rGazeGroundIntersection(:,1), thisWalk_orig.rGazeGroundIntersection(:,3),'r.')
            plot(thisWalk_orig.lGazeGroundIntersection(:,1), thisWalk_orig.lGazeGroundIntersection(:,3),'b.')
            title('Original Data - NOTE: Y axis scale magnifies offset error ')
            
            subplot(2,1,2)
            
            plot(thisWalk_fixed.comXYZ(:,1), thisWalk_fixed.comXYZ(:,3),'m')
            
            hold on
            plot(thisWalk_fixed.rGazeGroundIntersection(:,1), thisWalk_fixed.rGazeGroundIntersection(:,3),'r.')
            plot(thisWalk_fixed.lGazeGroundIntersection(:,1), thisWalk_fixed.lGazeGroundIntersection(:,3),'b.')
            title('rotated Data - Positive corrAlignTheta shifts move gaze points down. Try to get red/blue dots to align with COM path (magenta line)')
            
            drawnow
        end
        
        
    end
    
    if stillLooking
        
        disp('......')
        dbstack
        disp([thisWalk_orig.sessionID, ' ', thisWalk_orig.takeID, 'Walk# ',num2str(thisWalk_orig.ww) ' trial type: ' thisWalk_orig.trialType ])
        
        disp(['Current corrAlignTheta value is: ' num2str(corrAlignTheta) ', hows it look?'])
        newVal = input('Input new corrAlignTheta value (enter 99 to break loop):');
        disp('......')
        if newVal == 99
            stillLooking = false;
        else
            corrAlignTheta = newVal;
        end
    end
    
end
thisWalk_fixed.corrAlignTheta = corrAlignTheta;
end





%%% this is where we'll store the corrAlignTheta values that we're happy
%%% with

function [corrAlignTheta] = corrAlignLookup(sessionID, takeID, ww)
corrAlignTheta = nan;

switch sessionID
    case '2018-01-23_JSM'
        switch takeID
            case 'Woodchips'
                corrTable_ww = [ 0 .36 .07 .25 .15 .15]; %should be 7 values for Woodchips data (6 walks, plus the VOR at the beginning)
            case'Rocks'
                corrTable_ww = [-.35 -.3 -.35 -.3 -.35 -.32]; %should be 7 values for walks Data(6 walks, plus the VOR at the beginning)
        end
    case '2018-01-31_JAW'
        switch takeID
            case 'Woodchips'
                corrTable_ww = [ -.34 -.18 -.3 -.18 -.38 -.18 ]; %should be 7 values for Woodchips data (6 walks, plus the VOR at the beginning)
            case'Rocks'
                corrTable_ww = [ -.55 -.4 -.5 -.35 -.45 -.35]; %should be 7 values for walks Data(6 walks, plus the VOR at the beginning)
        end
        
    case '2018-01-26_JAC'
        switch takeID
            case 'Woodchips'
                corrTable_ww = [ -.3 -.12 -.3 .05 -.32 -.1 ]; %should be 7 values for Woodchips data (6 walks, plus the VOR at the beginning)
            case'Rocks'
                corrTable_ww = [ -.45 -.45 -.35 -.4 -.35 -.48 ]; %should be 7 values for walks Data(6 walks, plus the VOR at the beginning)
        end
    case '2019-02-27_JSM' %test
        switch takeID
              case'Rocks'
                corrTable_ww = [0.3 .24 .31 .23 0 0]; %test
        end
end

corrAlignTheta = corrTable_ww(ww);


end





