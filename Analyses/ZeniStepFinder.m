
% function [rightTO_HS, leftTO_HS] = ZeniStepFinder(shadow_fr_mar_dim,shadowMarkerNames, avg_fps)

% function [allStancePhases] = ZeniStepFinder(w)
function [allSteps_HS_TO_StanceLeg] = ZeniStepFinder(w)

shadow_fr_mar_dim = w.shadow_fr_mar_dim;
shadowMarkerNames = w.shadowMarkerNames;
avg_fps = mean(w.avg_fps);
walks = w.walks;


% The basic idea is to first situate the data by subtracting the x & y
% coordinates of the Root marker from each marker position at each frame,
% essentially setting the origin to the subject's root.
%
% The velocity of each foot in the X direction is then calculated. Zero
% crossings (within a band from Zero) are either toe-offs or heel strikes.
% Positive to Negative denotes heel strike, Negative to Positive denotes
% Toe off.


%working for threshold and working for new data and then check signs in directional walking

%define frames of interes as only frames when the person is moving --> for
%loop
%x=body velocity in 2d not zeroed because need movement and not relative
%movement
%frInt = x>___
%X(frInt) = new X
%interesct function
% framesOfInterest = 11000:13500;
c3dData = shadow_fr_mar_dim;

% Pull out right and left ankle data before they are zeroed and downsampled
% want all locations/velocities to be sum of z and x axes


rFootID = find(strcmp('RightFoot', shadowMarkerNames));
rAnkX_noZ = squeeze(c3dData(:,rFootID,1)); % pull out rFootID marker as X vector
rAnkY_noZ = squeeze(c3dData(:,rFootID,2)); % pull out rFootID marker as Y vector
rAnkZ_noZ = squeeze(c3dData(:,rFootID,3)); % pull out rFootID marker as Z vector

lFootID = find(strcmp('LeftFoot', shadowMarkerNames));
lAnkX_noZ = squeeze(c3dData(:,lFootID,1)); % pull out lFootID marker as X vector
lAnkY_noZ = squeeze(c3dData(:,lFootID,2)); % pull out lFootID marker as Y vector
lAnkZ_noZ = squeeze(c3dData(:,lFootID,3)); % pull out lFootID marker as Z vector

%define root as 'Body'
BodyID = find(strcmp('Hips',shadowMarkerNames)); %find body marker
root = squeeze(c3dData(:,BodyID,:));





% zero the data first (it's not already zeroed since not downsampled)
for f = 1:numel(c3dData(:,1,1)) %f = Frame
    for m = 1:numel(c3dData(1,:,1)) %m = Marker
        c3dData(f,m,1) = c3dData(f,m,1) - root(f,1);
        c3dData(f,m,2) = c3dData(f,m,2) - root(f,2);
        c3dData(f,m,3) = c3dData(f,m,3) - root(f,3);
    end
end


rThighXZ = squeeze(c3dData(:,strcmp('RightThigh', shadowMarkerNames), [1 3])); % pull out lFootID marker as X vector

[rThighTheta, rThighRho] = cart2pol(rThighXZ(:,1), rThighXZ(:,2));

for f = 1:numel(c3dData(:,1,1)) %f = Frame
    for m = 1:numel(c3dData(1,:,1)) %m = Marker
        
        X = c3dData(f,m,1);
        Z = c3dData(f,m,3);
        
        %         x_r = c3dData(f,m,1)* cos(-rThighTheta(f))- ...
        %             c3dData(f,m,3)* sin(-rThighTheta(f));
        %
        %         z_r = c3dData(f,m,1)* sin(-rThighTheta(f))+...
        %             c3dData(f,m,3)* cos(-rThighTheta(f));
        %
        
        x_r = ...
            X * cos(rThighTheta(f))+... %x*cos(theta)
            Z * sin(rThighTheta(f));    %y*sin(theta)
        
        
        z_r = ...
            -X * sin(rThighTheta(f))+... %x*cos(theta)
            Z * cos(rThighTheta(f));    %y*sin(theta)
        
        
        c3dData(f,m,1) = x_r;
        c3dData(f,m,3) = z_r;
        
    end
end


rThighXZrot = squeeze(c3dData(:,strcmp('RightThigh', shadowMarkerNames), [1 3])); % pull out lFootID marker as X vector




rAnkX = squeeze(c3dData(:,rFootID,1)); % pull out rFootID marker as X vector
rAnkY = squeeze(c3dData(:,rFootID,2)); % pull out rFootID marker as Y vector
rAnkZ = squeeze(c3dData(:,rFootID,3)); % pull out rFootID marker as Z vector

lAnkX = squeeze(c3dData(:,lFootID,1)); % pull out lFootID marker as X vector
lAnkY = squeeze(c3dData(:,lFootID,2)); % pull out lFootID marker as Y vector
lAnkZ = squeeze(c3dData(:,lFootID,3)); % pull out lFootID marker as Z vector


walkDir = zeros(length(root),1);
for ww = 1:size(walks,1)
    %
    %     if mean(sum(diff(root(walks(ww,1):walks(ww,2),[1 3])),2)) > 0 %if mean Root Vel in this walk is positive, walkDir is positive
    %         dir = 1;
    %     elseif mean(sum(diff(root(walks(ww,1):walks(ww,2),[1 3])),2)) < 0 %if mean Root Vel in this walk is positive, walkDir is positive
    %         dir = -1;
    %     end
    %
    dir = -1; %I think this is all you need now that I'm doing the thigh rotation thing, but based on my experience, there is an equal or greater chance that I am just fucking everything up :D
    
    walkDir(walks(ww,1):walks(ww,2)) = dir;
end

%
% %find frames wherein subject is moving
% rootVel = diff(root(:,1))*mean(avg_fps)/1000 + diff(root(:,3))*mean(avg_fps)/1000;
%
% order = 4;
% cutoff = .3;
%
% posRootVel = butterLowZero(order,cutoff,mean(avg_fps),rootVel);
% posRootVel(posRootVel<0) = nan;
%
% negRootVel = butterLowZero(order,cutoff,mean(avg_fps),rootVel);
% negRootVel(negRootVel>0) = nan;
%
%
%
% thresh = .3;
% moving = posRootVel>thresh | -negRootVel>thresh;
%
% walkDir = +(posRootVel>=0 | -negRootVel>0); %convert logical to double
% walkDir(rootVel<0) = -1; %tag walking direction with  a 1 or a -1
%
% walkDir = walkDir .* moving;

%% Right and Left Ankle markers
% add Xvel to Zvel to get 2d vel

order = 4;
cutoff = 3;

rAnkVelX = [0; diff(rAnkX)];
rAnkVelZ = [0; diff(rAnkZ)];
rAnkVel =  rAnkVelZ*mean(avg_fps)/1000; %now that I'm rotating the data, all the action should be happening in the Z direction. <<-Hey, FutureJon. Nice to see you. The problem is probably here <<--
% rAnkVel = plus(rAnkVelX, rAnkVelZ)*mean(avg_fps)/1000; <<-- this is what it used to look like ^^
rAnkVelRaw = rAnkVel.* walkDir;
rAnkVel = butterLowZero(order, cutoff, mean(avg_fps), rAnkVel);
rAnkVel = rAnkVel .* walkDir; %flip ankVel when sub is walking in "negative" direction. Also zero out non-moving frames

lAnkVelX = [0; diff(lAnkX)];
lAnkVelZ = [0; diff(lAnkZ)];
lAnkVel = lAnkVelZ*mean(avg_fps)/1000;
lAnkVelRaw = lAnkVel.* walkDir;
lAnkVel = butterLowZero(order, cutoff, mean(avg_fps), lAnkVel);
lAnkVel = lAnkVel .* walkDir; %flip ankVel when sub is walking in "negative" direction. Also zero out non-moving frames




lAnkVelX_noZ = [0; diff(lAnkX_noZ)];
lAnkVelZ_noZ = [0; diff(lAnkZ_noZ)];
lAnkVel_noZ = plus(lAnkVelX_noZ, lAnkVelZ_noZ)*mean(avg_fps)/1000;
lAnkVel_noZ = lAnkVel_noZ .* walkDir; %flip ankVel when sub is walking in "negative" direction. Also zero out non-moving frames



% % %
% % % %%correct for weird patches where "correctGroundSlip did its thing by linearly interpolating sectoins where Ankel Acceleration == 0. That's
% % % %%right, I'm correcting for the correction, cuz that's how I roll :-/
% % %
% % % rAnkAcc = [0; diff(rAnkVel)];
% % % lAnkAcc = [0; diff(lAnkVel)];
% % %
% % % rAnkAcc(abs(rAnkAcc) < 1e-10) = 0;
% % % lAnkAcc(abs(lAnkAcc) < 1e-10) = 0;
% % %
% % % rAnkAcc(~moving) = nan;
% % % lAnkAcc(~moving) = nan;
% % %
% % % wonkyFrames = lAnkAcc==0;
% % %
% % % wonkyFrameStarts  = find( diff([0; wonkyFrames]) == 1)-2;
% % % wonkyFrameEnds  = find( diff([0; wonkyFrames]) == -1)+2;
% % %
% % % % %%debug plot
% % % % plot(rAnkVel,'ro-')
% % % % hold on
% % % % plot(lAnkVel,'bo-')
% % %
% % % for ff = 1:length(wonkyFrameStarts) % replace wonky frames with linearly interpolated velocity
% % %
% % %     rAnkVel(wonkyFrameStarts(ff):wonkyFrameEnds(ff)) = linspace(rAnkVel(wonkyFrameStarts(ff)), rAnkVel(wonkyFrameEnds(ff)), length(wonkyFrameStarts(ff):wonkyFrameEnds(ff)));
% % %
% % %     lAnkVel(wonkyFrameStarts(ff):wonkyFrameEnds(ff)) = linspace(lAnkVel(wonkyFrameStarts(ff)), lAnkVel(wonkyFrameEnds(ff)), length(wonkyFrameStarts(ff):wonkyFrameEnds(ff)));
% % % end
% % %
% % % % plot(rAnkVel,'r.-')
% % % % hold on
% % % % plot(lAnkVel,'b.-')
% % % % hold off


%Find frames where sub is walking


% thresh = .1;
% startFrame = 0;
% endFrame = 0;
%
% for ii = 1:length(posRootVel)
%
%     if(posRootVel(ii))>thresh & startFrame ==0
%         startFrame = ii;
%         ii = ii + 200;
%     end
%
%     if(posRootVel(ii))<thresh & startFrame > 0 & endFrame == 0
%         endFrame = ii;
%         break
%     end
%
% end
%
% if endFrame == 0
%     endFrame = length(posRootVel);
% end
%
% % %%% Debug Plot
% figure;clf;
% plot(posRootVel)
% hold on
% plot([startFrame endFrame], posRootVel([startFrame endFrame]), 'o')
%
startFrame = 1;
endFrame = length(root);

rTO = [];
rHS = [];

lTO = [];
lHS = [];


for i = startFrame+1:endFrame-1
    
    %Right foot toe offs
    if rAnkVel(i-1) < 0 && rAnkVel(i) > 0
        rTO(end+1) = i;
        
    end
    
    %Left foot Toe offs
    if (lAnkVel(i-1) < 0 && lAnkVel(i) > 0)
        lTO(end+1) = i;
        
    end
    
    
    
    
    if (rAnkVel(i-1) > 0 && rAnkVel(i) < 0)
        rHS(end+1) = i;
    end
    
    %Left foot heel strikes
    if (lAnkVel(i-1) > 0 && lAnkVel(i) < 0)
        lHS(end+1) = i;
        
    end
    
    
    
end


rTO = rTO';
lTO = lTO';
rHS = rHS';
lHS = lHS';


% Uncomment for debugging Plots
figure(43242)
clf
title('Raw steps raised, cleaned steps on 0 line')
subplot(2,1,1)
refline(0,0); hold on;
plot(rAnkVelRaw,'m.-','MarkerSize',2)
plot(rAnkVel,'r-o','MarkerSize',2)
plot(rTO,.01,'ko','MarkerFaceColor','m','MarkerSize',8)
plot(rHS,.01,'kp','MarkerFaceColor','m','MarkerSize',8)
grid on

subplot(212)
refline(0,0); hold on;
% plot(lAnkVelRaw,'c.-','MarkerSize',2)
plot(lAnkVel, 'b-o','MarkerSize',2)
plot(lTO(:,1),.01,'ko','MarkerFaceColor','g','MarkerSize',8)
plot(lHS(:,1),.01,'kp','MarkerFaceColor','g','MarkerSize',8)
grid on


figure(432421)
clf
refline(0,0); hold on;
plot(rAnkVelRaw,'m.-')
plot(rAnkVel,'r.-')
plot(rTO(:,1),.02,'ko','MarkerFaceColor','m','MarkerSize',8)
plot(rHS(:,1),.02,'kp','MarkerFaceColor','m','MarkerSize',8)
grid on

refline(0,0); hold on;
plot(lAnkVelRaw,'c.-')
plot(lAnkVel, 'b.-')
plot(lTO(:,1),.02,'ko','MarkerFaceColor','g','MarkerSize',8)
plot(lHS(:,1),.02,'kp','MarkerFaceColor','g','MarkerSize',8)
grid on





%% Clean up

%%% start with (temporally) linear sequence of all gait events, tagged with
%%% foot (R = 1, L = 2) and type (Toe off = 3, Heel strike = 4)


rTO(:,2) = 1; %'1' means this is something the right leg is doing
rHS(:,2) = 1;

rTO(:,3) = 3; % '3' means it's a Toe Off (p.s. using 3 and 4 b/c 1 and 2 are used to mean "right" and "left" elsewhere in this fxn)
rHS(:,3) = 4; % '4' means it's a Heel Strike

lTO(:,2) = 2; %'2' means this is something the right leg is doing
lHS(:,2) = 2;

lTO(:,3) = 3; % '3' means it's a Toe Off
lHS(:,3) = 4; % '4' means it's a Heel Strike


allGaitEvents = sortrows([ rTO; rHS; lTO; lHS]);

%remove Same-Leg Gait events that are below a threshold
minEventDur = round(.200/(1/avg_fps)); %200ms (0.2sec) = 20frames @~100fps

for gg = 1:length(allGaitEvents)-1
    
    if allGaitEvents(gg,2) == allGaitEvents(gg+1,2) % if these two subsequent events occur on the same leg
        if (allGaitEvents(gg+1,1)- allGaitEvents(gg,1)) < minEventDur %and the events are separated by less time than threshold..
            allGaitEvents(gg:gg+1,:) = [nan nan nan; nan nan nan]; %NaN 'em out!!
        end
    end
end

allGaitEvents( isnan(allGaitEvents(:,1)),:) = [];


%tag suspicious steps that don't conform to normal gait pattern (i.e. TO-HS with one leg, folloed by TO-HS with the other leg)
while allGaitEvents(1,3) == 4 %it's easier if we start with a toe off
    allGaitEvents(1,:) = [];
end

allGaitEventsRaw = allGaitEvents;

iter = 0;
suspect = zeros(length(allGaitEvents),1);

while sum(suspect)>0 || iter == 0
    
    suspect = zeros(length(allGaitEvents),1);
    iter = iter+1;
    
    for ee = 1:length(allGaitEvents)-1
        
        %% Sequence should be ...TO(L)-HS(L)-TO(R)-HS(R)-TO(R)-HS(L)-TO(L)....
        
        %If this event is  TO, AND the next event is EITHER not a HS, OR not on the same leg (P & (Q or R)).... SUSPECT!!
        if allGaitEvents(ee,3) ==3 &&  (allGaitEvents(ee+1,3) ~= 4 || allGaitEvents(ee,2) ~= allGaitEvents(ee+1,2))
            suspect(ee+1) = 1;
        end
        
        
        if allGaitEvents(ee,3) ==4 && (allGaitEvents(ee+1,3) ~= 3 || allGaitEvents(ee,2) == allGaitEvents(ee+1,2))
            suspect(ee+1) = 1;
        end
    end
    allGaitEvents(logical(suspect),:) = [];
    
end

%% Debug plot
rTO_0 = allGaitEvents( [allGaitEvents(:,2) == 1 & allGaitEvents(:,3) == 3],:);
rHS_0 = allGaitEvents( [allGaitEvents(:,2) == 1 & allGaitEvents(:,3) == 4],:);
lTO_0 = allGaitEvents( [allGaitEvents(:,2) == 2 & allGaitEvents(:,3) == 3],:);
lHS_0 = allGaitEvents( [allGaitEvents(:,2) == 2 & allGaitEvents(:,3) == 4],:);


figure(432421)

plot(rTO_0(:,1),0,'mo','MarkerFaceColor','k','MarkerSize',12)
plot(rHS_0(:,1),0,'mp','MarkerFaceColor','k','MarkerSize',12)
grid on

plot(lTO_0(:,1),0,'go','MarkerFaceColor','k','MarkerSize',12)
plot(lHS_0(:,1),0,'gp','MarkerFaceColor','k','MarkerSize',12)
grid on


%% build output variable(s)
%
% %% build "allStancePhases"
% rAllEvents = allGaitEvents(allGaitEvents(:,2)==1, :);
% lAllEvents = allGaitEvents(allGaitEvents(:,2)==2, :);
%
% if ~( sum(diff(rAllEvents(:,3))==0) == 0) || ~( sum(diff(lAllEvents(:,3))==0) == 0)  %if there are any TO-TO or HS-HS's in either variable, stop
%     disp('Problem in Zeni')
%     dbstack
%     keyboard
% end
%
% while rAllEvents(1,3) ~=4; %make sure the first entry is a HS
%     rAllEvents(1,:) = [];
% end
%
%
% while rAllEvents(end,3) ~=3; %make sure the first entry is a TO
%     rAllEvents(end,:) = [];
% end
%
%
% while lAllEvents(1,3) ~=4; %make sure the first entry is a HS
%     lAllEvents(1,:) = [];
% end
%
%
% while lAllEvents(end,3) ~=3; %make sure the first entry is a TO
%     lAllEvents(end,:) = [];
% end
%
% rHSTO(:,1) = rAllEvents(rAllEvents(:,3)==4); %r HS's
% rHSTO(:,2) = rAllEvents(rAllEvents(:,3)==3); %r TO's
% rHSTO(:,3) = ones(length(rHSTO),1)*1; %tag right leg events with a '1'
%
%
% lHSTO(:,1) = lAllEvents(rAllEvents(:,3)==4); %lHS's
% lHSTO(:,2) = lAllEvents(rAllEvents(:,3)==3); %lTO's
% lHSTO(:,3) = ones(length(lHSTO),1)*2; %tag left leg events with a '3'
%
%
% allStancePhases = sortrows([rHSTO; lHSTO]);
%
% %delete steps that are longer than 1.5 seconds
%
% stanceDur = allStancePhases(:,2) - allStancePhases(:,1);
% maxStanceDur = round(1.5/(1/avg_fps)); %1.5s = 150frames @~100fps
%
% allStancePhases(stanceDur > maxStanceDur, :) = [];
%

%% build allSteps
ss= 0;
% if allGaitEvents(1,3) == 3;
%     ss= ss+1; %skip the first entry if it is a TO
% end

allSteps_HS_TO_StanceLeg = [];

while ss < length(allGaitEvents)-1
    ss = ss+1;
    
    
    if allGaitEvents(ss:ss+1,3) == [4; 3]
        
        allSteps_HS_TO_StanceLeg(end+1,1) = allGaitEvents(ss   , 1); %%this Heelstrike (beginning of double support)
        allSteps_HS_TO_StanceLeg(end,3)   = allGaitEvents(ss   , 2); %this stance leg (which leg is doing the heel striking, 1 = r, 2 = L)
        allSteps_HS_TO_StanceLeg(end,2)   = allGaitEvents(ss+1 , 1); %%this Toe Off (Beginning of single support)
        
        if allGaitEvents(ss+1 , 2) ==  allSteps_HS_TO_StanceLeg(end,3) %HS and TO should happen with differnt feet
            disp('Problem in Zeni')
            dbstack
            keyboard
        end
        
        
        ss = ss+1;
    end
    
end




%% Debug plot
rHS_1 =  allSteps_HS_TO_StanceLeg( [ allSteps_HS_TO_StanceLeg(:,3) == 1],1);
lHS_1 =  allSteps_HS_TO_StanceLeg( [ allSteps_HS_TO_StanceLeg(:,3) == 2],1);

rTO_1 =  allSteps_HS_TO_StanceLeg( [ allSteps_HS_TO_StanceLeg(:,3) == 2],2);%rTOs come second on rows starting with a lHS
lTO_1 =  allSteps_HS_TO_StanceLeg( [ allSteps_HS_TO_StanceLeg(:,3) == 1],2);%and vice versa

figure(432421)


plot(rTO_1(:,1),0,'ro','MarkerFaceColor','r','MarkerSize',4)
plot(rHS_1(:,1),0,'ro','MarkerFaceColor','r','MarkerSize',4)
grid on

plot(lTO_1(:,1),0,'go','MarkerFaceColor','g','MarkerSize',4)
plot(lHS_1(:,1),0,'go','MarkerFaceColor','g','MarkerSize',4)
grid on

for i = 1:length(allSteps_HS_TO_StanceLeg)-1
    ssp(i) = allSteps_HS_TO_StanceLeg(i+1,1)- allSteps_HS_TO_StanceLeg(i,2);
end

figure(3432432)
title('Step Durations')
plot(diff(allSteps_HS_TO_StanceLeg(:,1)),'-p');
hold on
plot(ssp,'-v')
plot(diff(allSteps_HS_TO_StanceLeg(:,[1 2]),1,2),'-o')
plot(allSteps_HS_TO_StanceLeg(:,3),'-s')

legend('Full Step (HS-HS)',  'SingleSupport (TO-HS)', 'DoubleSupport( HS-TO)','StanceLeg (1=R, 2=L)')
ylim([0, 100])
hold off







% %% Debug plot
% rTO_1 =  allStancePhases( [ allStancePhases(:,3) == 1],2);
% rHS_1 =  allStancePhases( [ allStancePhases(:,3) == 1],1);
%
% lTO_1 =  allStancePhases( [ allStancePhases(:,3) == 2],2);
% lHS_1 =  allStancePhases( [ allStancePhases(:,3) == 2],1);
%
%
% figure(432421)
%
% plot(rTO_1(:,1),0,'ro','MarkerFaceColor','r','MarkerSize',4)
% plot(rHS_1(:,1),0,'ro','MarkerFaceColor','r','MarkerSize',4)
% grid on
%
% plot(lTO_1(:,1),0,'ro','MarkerFaceColor','r','MarkerSize',4)
% plot(lHS_1(:,1),0,'ro','MarkerFaceColor','r','MarkerSize',4)
% grid on
%

% %%
%
%
% %Uncomment for Debug Plots
% figure(43242)
% subplot(2,1,1)
% plot(rAnkVel,'r-o','MarkerSize',2)
% hold on
% plot(allStancePhases(allStancePhases(:,3) == 1, 2),0,'ro','MarkerFaceColor','k','MarkerSize',8)
% plot(allStancePhases(allStancePhases(:,3) == 1, 1),0,'rp','MarkerFaceColor','k','MarkerSize',8)
% grid on
%
% subplot(212)
% plot(lAnkVel, 'b-o','MarkerSize',2)
% hold on
% plot(allStancePhases(allStancePhases(:,3) == 2, 2),0,'bo','MarkerFaceColor','k','MarkerSize',8)
% plot(allStancePhases(allStancePhases(:,3) == 2, 1),0,'bp','MarkerFaceColor','k','MarkerSize',8)
% grid on
% hold off

