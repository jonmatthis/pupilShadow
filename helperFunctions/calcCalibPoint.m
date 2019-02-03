 function [ calibPoint ] = calcCalibPoint( shadow_fr_mar_dim, shadowMarkerNames, calibFrame)
%CALCCALIBMATPOINTS calculate location of calibration points on calibration
%carpet thing. The subject stands at/near/above? (0,0) on ground plane.
%Calibrates to point at ([0,1000] (1m in front of them). Other points are
%at [0,1400], [-450,1400].[450,1400], and [(+/-?)450,1000]

thisFrame = calibFrame;

rHeelID = find(strcmp('RightHeel', shadowMarkerNames));
rHeelXYZ = squeeze(shadow_fr_mar_dim(:,rHeelID,:)); % pull out rHeelID marker

rToeID = find(strcmp('RightToe', shadowMarkerNames));
rToeXYZ = squeeze(shadow_fr_mar_dim(:,rToeID,:)); % pull out rHeelID marker

rAnkID = find(strcmp('RightFoot', shadowMarkerNames));
rAnkXYZ = squeeze(shadow_fr_mar_dim(:,rAnkID,:)); 

lHeelID = find(strcmp('LeftHeel', shadowMarkerNames));
lHeelXYZ = squeeze(shadow_fr_mar_dim(:,lHeelID,:)); % pull out rHeelID marker

lToeID = find(strcmp('LeftToe', shadowMarkerNames));
lToeXYZ = squeeze(shadow_fr_mar_dim(:,lToeID,:)); % pull out rHeelID marker

lAnkID = find(strcmp('LeftFoot', shadowMarkerNames));
lAnkXYZ = squeeze(shadow_fr_mar_dim(:,lAnkID,:)); 

headID = find(strcmp('HeadTop', shadowMarkerNames));
headXYZ = squeeze(shadow_fr_mar_dim(:,headID,:)); % This will be all Zeros if zero-ing is turned on


%% grab your mean heel and toe
grHei = min([ rHeelXYZ(thisFrame,2) lHeelXYZ(thisFrame,2)]); %the Z values denote vertical position

toeMean = mean([rToeXYZ(thisFrame,:); lToeXYZ(thisFrame,:)]);
toeMean(2) = grHei;

heelMean = mean([rHeelXYZ(thisFrame,:); lHeelXYZ(thisFrame,:)]);
heelMean(2) = grHei;

ankMean = mean([rAnkXYZ(thisFrame,:); lAnkXYZ(thisFrame,:)]);
ankMean(2) = grHei;

%% %%%debug - Plot skeleton and mean heel/toe/ank postion

ii = calibFrame;
figure(9384);clf
plot3(rHeelXYZ(ii,1), rHeelXYZ(ii,2), rHeelXYZ(ii,3),'r^')
hold on
plot3(rToeXYZ(ii,1), rToeXYZ(ii,2), rToeXYZ(ii,3),'ro')

plot3(lHeelXYZ(ii,1), lHeelXYZ(ii,2), lHeelXYZ(ii,3),'b^')
plot3(lToeXYZ(ii,1), lToeXYZ(ii,2), lToeXYZ(ii,3),'bo')

plot3(heelMean(1), heelMean(2), heelMean(3), 'ko')
plot3(toeMean(1), toeMean(2), toeMean(3), 'ko')
plot3(ankMean(1), ankMean(2), ankMean(3), 'kp')

lLeg = [2 3 4 5 6 7 5];
rLeg = [2 8 9 10 11 12 10];
tors = [2 13 14 15 26 27 28];
lArm = [15 16 17 26 17 18 19 20];
rArm = [15 21 22 26 22 23 24 25];

plot3(shadow_fr_mar_dim(thisFrame,1:28,1),shadow_fr_mar_dim(thisFrame,1:28,2),shadow_fr_mar_dim(thisFrame,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
hold on


plot3(shadow_fr_mar_dim(thisFrame,lLeg,1),shadow_fr_mar_dim(thisFrame,lLeg,2),shadow_fr_mar_dim(thisFrame,lLeg,3),'b','LineWidth',2)
plot3(shadow_fr_mar_dim(thisFrame,rLeg,1),shadow_fr_mar_dim(thisFrame,rLeg,2),shadow_fr_mar_dim(thisFrame,rLeg,3),'r','LineWidth',2)
plot3(shadow_fr_mar_dim(thisFrame,tors,1),shadow_fr_mar_dim(thisFrame,tors,2),shadow_fr_mar_dim(thisFrame,tors,3),'g','LineWidth',2)
plot3(shadow_fr_mar_dim(thisFrame,lArm,1),shadow_fr_mar_dim(thisFrame,lArm,2),shadow_fr_mar_dim(thisFrame,lArm,3),'b','LineWidth',2)
plot3(shadow_fr_mar_dim(thisFrame,rArm,1),shadow_fr_mar_dim(thisFrame,rArm,2),shadow_fr_mar_dim(thisFrame,rArm,3),'r','LineWidth',2)
axis equal




g_x = meshgrid(rHeelXYZ(ii,1)-2000:500:rHeelXYZ(ii,1)+2000);
g_z = meshgrid(rHeelXYZ(ii,2)-2000:500:rHeelXYZ(ii,2)+2000)';
g_y = ones(size(g_x)) * min([rHeelXYZ(ii,2) lHeelXYZ(ii,2) ]);

surface(g_x, g_y, g_z,'FaceColor','w','EdgeColor','k')

%% make ankle-toe vectors
X = [ankMean(1) toeMean(1)];
Y = [ankMean(2) toeMean(2)];
Z = [ankMean(3) toeMean(3)];
plot3(X,Y,Z,'r-','LineWidth',2)

%%%put Toe vector in ankle-centered frame, and normalize it
normToeInAnkRef  = (toeMean-ankMean)./norm(toeMean-ankMean);
toeInAnkRefScaled = normToeInAnkRef * 1400;

calibPoint = toeInAnkRefScaled + ankMean;

plot3([ankMean(1) calibPoint(1)],[ankMean(2) calibPoint(2)],[ankMean(3) calibPoint(3)],'rp-')

