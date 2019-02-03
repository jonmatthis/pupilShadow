function [shadowFixed_fr_mar_dim] = fixSkateboarding(wRaw, allSteps_HS_TO_StanceLeg)
%Try to fix the thing where shadow data appears to "skateboard" because it
%fails to recognize a step. In this function, I'ma try to pin each foot to
%the ground at heelstrike and subtract the 'false' motion from the rest of
%the skeleton. Hopefully, that'll fix that up. We'll see :/

wFixed = wRaw;

shadowRaw_fr_mar_dim = wRaw.shadow_fr_mar_dim;

shadowMarkerNames = wRaw.shadowMarkerNames;

hs = allSteps_HS_TO_StanceLeg(:,1);
stanceLeg = allSteps_HS_TO_StanceLeg(:,3);


shadowFixed_fr_mar_dim = shadowRaw_fr_mar_dim; %we'll be fixing this one

%%marker ID's to assist with skeleton plotting
lLeg = [2 3 4 5 6 7 5];
rLeg = [2 8 9 10 11 12 10];
tors = [2 13 14 15 26 27 28];
lArm = [15 16 17 26 17 18 19 20];
rArm = [15 21 22 26 22 23 24 25];


rheel = strcmp('RightHeel', shadowMarkerNames);
lheel = strcmp('LeftHeel', shadowMarkerNames);
figure(146789)

length(hs)

for ss = 1:length(hs)-1
%     disp(strcat({'Anti-Skateboarding Step# '},num2str(ss),{' of '}, num2str(length(hs))))
        if stanceLeg(ss) == 1
            thisFootholdXYZ = squeeze(shadowFixed_fr_mar_dim(hs(ss),rheel,:)); % pull out lHeel marker;
        elseif stanceLeg(ss) == 2
            thisFootholdXYZ = squeeze(shadowFixed_fr_mar_dim(hs(ss),lheel,:)); % pull out lHeel marker;
        end

    
    fr = hs(ss):hs(ss+1)-1;
  
    % correct for slippage during a step

    if stanceLeg(ss) == 1
        currHeelXYZ =   squeeze(shadowFixed_fr_mar_dim(fr,rheel,:)); % pull out lHeel marker
    elseif stanceLeg(ss) == 2
        currHeelXYZ =   squeeze(shadowFixed_fr_mar_dim(fr,lheel,:)); % pull out lHeel marker
    end


        
        %%find error between current heel position and heel position at
        %%heel strike, and add that error to all the subsequent skeleton
        %%data. Theoretically, this should pin the heel to it's position at
        %%HS
        
        slippageXYZ = thisFootholdXYZ' - currHeelXYZ;
        
                
%         for mm = 1:length(shadowMarkerNames)
           shadowFixed_fr_mar_dim(fr,:, 1) = shadowFixed_fr_mar_dim(fr,:, 1) + slippageXYZ(:,1);
           shadowFixed_fr_mar_dim(fr,:, 2) = shadowFixed_fr_mar_dim(fr,:, 2) + slippageXYZ(:,2);
           shadowFixed_fr_mar_dim(fr,:, 3) = shadowFixed_fr_mar_dim(fr,:, 3) + slippageXYZ(:,3);

           shadowFixed_fr_mar_dim(fr(end)+1:end,:, 1) = shadowFixed_fr_mar_dim(fr(end)+1:end,:, 1) + slippageXYZ(end,1);
           shadowFixed_fr_mar_dim(fr(end)+1:end,:, 2) = shadowFixed_fr_mar_dim(fr(end)+1:end,:, 2) + slippageXYZ(end,2);
           shadowFixed_fr_mar_dim(fr(end)+1:end,:, 3) = shadowFixed_fr_mar_dim(fr(end)+1:end,:, 3) + slippageXYZ(end,3);

           
%         end
        
        plotSkel = false;
        
        if plotSkel
        %%Plot raw skeleton data
        clf
        plot3(shadowFixed_fr_mar_dim(fr,1:28,1),shadowFixed_fr_mar_dim(fr,1:28,2),shadowFixed_fr_mar_dim(fr,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
        hold on
        
        plot3(shadowFixed_fr_mar_dim(fr,lLeg,1),shadowFixed_fr_mar_dim(fr,lLeg,2),shadowFixed_fr_mar_dim(fr,lLeg,3),'b','LineWidth',2)
        plot3(shadowFixed_fr_mar_dim(fr,rLeg,1),shadowFixed_fr_mar_dim(fr,rLeg,2),shadowFixed_fr_mar_dim(fr,rLeg,3),'r','LineWidth',2)
        plot3(shadowFixed_fr_mar_dim(fr,tors,1),shadowFixed_fr_mar_dim(fr,tors,2),shadowFixed_fr_mar_dim(fr,tors,3),'g','LineWidth',2)
        plot3(shadowFixed_fr_mar_dim(fr,lArm,1),shadowFixed_fr_mar_dim(fr,lArm,2),shadowFixed_fr_mar_dim(fr,lArm,3),'b','LineWidth',2)
        plot3(shadowFixed_fr_mar_dim(fr,rArm,1),shadowFixed_fr_mar_dim(fr,rArm,2),shadowFixed_fr_mar_dim(fr,rArm,3),'r','LineWidth',2)
        
        
        
        rHeelXYZ = squeeze(shadowFixed_fr_mar_dim(:,strcmp('RightHeel', shadowMarkerNames),:)); % pull out lHeel marker;
        lHeelXYZ = squeeze(shadowFixed_fr_mar_dim(:,strcmp('LeftHeel', shadowMarkerNames),:)); % pull out lHeel marker;

        if stanceLeg(ss) == 1
            plot3(rHeelXYZ(hs(ss),1), rHeelXYZ(hs(ss),2), rHeelXYZ(hs(ss),3),'rp','MarkerSize',12,'MarkerFaceColor','r')
            
        elseif stanceLeg(ss) == 2
            plot3(lHeelXYZ(hs(ss),1), lHeelXYZ(hs(ss),2), lHeelXYZ(hs(ss),3),'bp','MarkerSize',12,'MarkerFaceColor','b')
        end
            
        axis equal
        title(num2str(fr))
        xlabel('x');ylabel('y'); zlabel('z');
        
        a = gca;
        comXYZ = squeeze(shadowFixed_fr_mar_dim(:,1,:));

        a.CameraTarget = [comXYZ(fr,1), comXYZ(fr,2), comXYZ(fr,3)]; %point figure 'camera' at COM
        %         a.CameraPosition = a.CameraTarget + [-1800 1800 2000]; %set camera position
        %         a.CameraViewAngle = 80;
        a.CameraUpVector = [ 0 1 0];
        a.Position = [0 0 1 1];
        
        
        drawnow
    end
        

       
       if mod(ss, 100) == 0
%            disp(fr)
       clf
       plot(squeeze(shadowRaw_fr_mar_dim(:,1,:)), '-','DisplayName','Raw Skel Data')
       hold on
       plot(squeeze(shadowFixed_fr_mar_dim(:,1,:)), '.-','DisplayName','Fixed Skel Data')
       legend('show')
       title('Skateboarding *is* a crime')
       drawnow
       end

end

figure(33423)
plot(squeeze(shadowRaw_fr_mar_dim(:,1,:)), '-','DisplayName','Raw Skel Data')
hold on
plot(squeeze(shadowFixed_fr_mar_dim(:,1,:)), '.-','DisplayName','Fixed Skel Data')
legend('show')
title('Skateboarding *is* a crime')
drawnow