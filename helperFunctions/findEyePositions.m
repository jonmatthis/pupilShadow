function [ rEyeballCenterXYZ, lEyeballCenterXYZ, worldCamCenterXYZ ] = findEyePositions( headGlobalQuat_wxyz, shadow_fr_mar_dim, shadowMarkerNames, calibFrame,shadowVersion)
%FINDSCENECAMERAPOSITION calculate the location oof the eyeballs relative
%to subject's head markers. This is the point from which the gaze vectors should
%originate.




%%
hTop =  squeeze(shadow_fr_mar_dim(:,strcmp('HeadTop', shadowMarkerNames),1:3)); %marker at Head Top
hC1 =   squeeze(shadow_fr_mar_dim(:,strcmp('Head', shadowMarkerNames),1:3)); %marker at C1 vertebra

lToe = squeeze(shadow_fr_mar_dim(:,strcmp('LeftToe', shadowMarkerNames),1:3));
rToe = squeeze(shadow_fr_mar_dim(:,strcmp('RightToe', shadowMarkerNames),1:3));



hCen = (hTop+hC1)/2; 

hCenCalib = hCen(calibFrame,:);


rEyeToeCalib_XYZ = [rToe(calibFrame,1) hCenCalib(2) rToe(calibFrame,3) ]; % As a starting guess - assume that the eyes are direction above the toes at the height of the head center... yeah... 
lEyeToeCalib_XYZ = [lToe(calibFrame,1) hCenCalib(2) lToe(calibFrame,3) ]; % As a starting guess - assume that the eyes are direction above the toes at the height of the head center... yeah... 



%%% scale eyeToe guess until it has an interocular separation of 70mm

rEyeCalib_XYZ_z = rEyeToeCalib_XYZ - hCenCalib;
lEyeCalib_XYZ_z = lEyeToeCalib_XYZ - hCenCalib;

iod = pdist([rEyeCalib_XYZ_z; lEyeCalib_XYZ_z]); %interocular distance

while pdist([rEyeCalib_XYZ_z; lEyeCalib_XYZ_z]) > 70
    rEyeCalib_XYZ_z = rEyeCalib_XYZ_z * .99;
    lEyeCalib_XYZ_z = lEyeCalib_XYZ_z * .99;
    iod(end+1) = pdist([rEyeCalib_XYZ_z; lEyeCalib_XYZ_z]); %interocular distance
end
% 
% % Debug Plot
% 
% figure(2)
% sphRes = 20;
% r = ones(sphRes, sphRes);
% [th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
% [x1,y1,z1] = sph2cart(th, phi, r);
% 
% 
% lLeg = [2 3 4 5 6 7 5];
% rLeg = [2 8 9 10 11 12 10];
% tors = [2 13 14 15 26 27 28];
% lArm = [15 16 17 26 17 18 19 20];
% rArm = [15 21 22 26 22 23 24 25];
% 
% i = calibFrame;
% 
% hx = hTop(i,1);
% hy = hTop(i,2);
% hz = hTop(i,3);
% 
% plot3(shadow_fr_mar_dim(i,1:28,1),shadow_fr_mar_dim(i,1:28,2),shadow_fr_mar_dim(i,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
% hold on
% 
% 
% plot3(shadow_fr_mar_dim(i,lLeg,1),shadow_fr_mar_dim(i,lLeg,2),shadow_fr_mar_dim(i,lLeg,3),'b','LineWidth',2)
% plot3(shadow_fr_mar_dim(i,rLeg,1),shadow_fr_mar_dim(i,rLeg,2),shadow_fr_mar_dim(i,rLeg,3),'r','LineWidth',2)
% plot3(shadow_fr_mar_dim(i,tors,1),shadow_fr_mar_dim(i,tors,2),shadow_fr_mar_dim(i,tors,3),'g','LineWidth',2)
% plot3(shadow_fr_mar_dim(i,lArm,1),shadow_fr_mar_dim(i,lArm,2),shadow_fr_mar_dim(i,lArm,3),'b','LineWidth',2)
% plot3(shadow_fr_mar_dim(i,rArm,1),shadow_fr_mar_dim(i,rArm,2),shadow_fr_mar_dim(i,rArm,3),'r','LineWidth',2)
% 
% 
% plot3(hx,hy, hz,'mh', 'MarkerSize', 20)
% 
% plot3(hTop(i,1), hTop(i,2), hTop(i,3),'gp', 'MarkerSize', 15)
% plot3(hC1(i,1), hC1(i,2), hC1(i,3),'kp', 'MarkerSize', 15)
% 
% plot3(rEyeToeCalib_XYZ(1),rEyeToeCalib_XYZ(2), rEyeToeCalib_XYZ(3),'rp')
% plot3(lEyeToeCalib_XYZ(1),lEyeToeCalib_XYZ(2), lEyeToeCalib_XYZ(3),'bp')
% 
% plot3(rEyeCalib_XYZ_z(1)+hCen(calibFrame, 1),rEyeCalib_XYZ_z(2)+hCen(calibFrame, 2), rEyeCalib_XYZ_z(3)+hCen(calibFrame, 3),'ro')
% plot3(lEyeCalib_XYZ_z(1)+hCen(calibFrame, 1),lEyeCalib_XYZ_z(2)+hCen(calibFrame, 2), lEyeCalib_XYZ_z(3)+hCen(calibFrame, 3),'bo')
% 
% 
% 
% %     surface(x1+hx,y1+hy,z1+hz,'FaceColor', 'none','EdgeColor',[.9 .9 .9])
% bx =   shadow_fr_mar_dim(i,1,1);
% by =   shadow_fr_mar_dim(i,1,2);
% bz =   shadow_fr_mar_dim(i,1,3);
% 
% axis([-2000+bx 2000+bx -2000+by 2000+by -2500+bz 2000+bz])
% 
% 
% axis equal
% title(num2str(i))
% set(gca,'CameraUpVector',[0 1 0])
% xlabel('x');ylabel('y'); zlabel('z');
% 
% hold off
% drawnow








% find unit vectors f head reference frame
headRotMat_row_col_fr = headGlobalQuat_wxyz.RotationMatrix;

calibRotMat = headRotMat_row_col_fr(:,:,calibFrame);

hux = calibRotMat * [1; 0; 0];
huy = calibRotMat * [0; 1; 0];
huz = calibRotMat * [0; 0; 1];

rEye_ux = dot(rEyeCalib_XYZ_z, hux);
rEye_uy = dot(rEyeCalib_XYZ_z, huy);
rEye_uz = dot(rEyeCalib_XYZ_z, huz);

lEye_ux = dot(lEyeCalib_XYZ_z, hux);
lEye_uy = dot(lEyeCalib_XYZ_z, huy);
lEye_uz = dot(lEyeCalib_XYZ_z, huz);


% % Debug Plot
% 
% figure(2)
% sphRes = 20;
% r = ones(sphRes, sphRes);
% [th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
% [x1,y1,z1] = sph2cart(th, phi, r);
% 
% 
% lLeg = [2 3 4 5 6 7 5];
% rLeg = [2 8 9 10 11 12 10];
% tors = [2 13 14 15 26 27 28];
% lArm = [15 16 17 26 17 18 19 20];
% rArm = [15 21 22 26 22 23 24 25];
% 
%  i = calibFrame;
%     
%     hx = hTop(i,1);
%     hy = hTop(i,2);
%     hz = hTop(i,3);
%     
%     plot3(shadow_fr_mar_dim(i,1:28,1),shadow_fr_mar_dim(i,1:28,2),shadow_fr_mar_dim(i,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
%     hold on
%     
%     
%     plot3(shadow_fr_mar_dim(i,lLeg,1),shadow_fr_mar_dim(i,lLeg,2),shadow_fr_mar_dim(i,lLeg,3),'b','LineWidth',2)
%     plot3(shadow_fr_mar_dim(i,rLeg,1),shadow_fr_mar_dim(i,rLeg,2),shadow_fr_mar_dim(i,rLeg,3),'r','LineWidth',2)
%     plot3(shadow_fr_mar_dim(i,tors,1),shadow_fr_mar_dim(i,tors,2),shadow_fr_mar_dim(i,tors,3),'g','LineWidth',2)
%     plot3(shadow_fr_mar_dim(i,lArm,1),shadow_fr_mar_dim(i,lArm,2),shadow_fr_mar_dim(i,lArm,3),'b','LineWidth',2)
%     plot3(shadow_fr_mar_dim(i,rArm,1),shadow_fr_mar_dim(i,rArm,2),shadow_fr_mar_dim(i,rArm,3),'r','LineWidth',2)
%     
%     
%         plot3(hx,hy, hz,'mh', 'MarkerSize', 20)
%         
%     plot3(hTop(i,1), hTop(i,2), hTop(i,3),'gp', 'MarkerSize', 15)
%     plot3(hC1(i,1), hC1(i,2), hC1(i,3),'kp', 'MarkerSize', 15)
% 
%         
%     plot3(headCalibVec(:,1),headCalibVec(:,2), headCalibVec(:,3),'rp-', 'MarkerSize', 15)
%     plot3(headCalibVecScaled(:,1),headCalibVecScaled(:,2), headCalibVecScaled(:,3),'bp-')
%     
%         plot3(calibPoints(:,1),calibPoints(:,2), calibPoints(:,3),'h')
% 
%     %     surface(x1+hx,y1+hy,z1+hz,'FaceColor', 'none','EdgeColor',[.9 .9 .9])
%       bx =   shadow_fr_mar_dim(i,1,1);
%       by =   shadow_fr_mar_dim(i,1,2);
%       bz =   shadow_fr_mar_dim(i,1,3);
%    
%     axis([-2000+bx 2000+bx -2000+by 2000+by -2500+bz 2000+bz])
%     
%        
%     axis equal
%     title(num2str(i))
%     set(gca,'CameraUpVector',[0 1 0])
%     xlabel('x');ylabel('y'); zlabel('z');
%     
%     hold off
%     drawnow
%     
% keyboard


%%



%make a head vector from IMU rotation matrix
rEyeballCenterXYZ = nan(length(headRotMat_row_col_fr),3);
lEyeballCenterXYZ = nan(length(headRotMat_row_col_fr),3);
for ii = 2:length(headRotMat_row_col_fr)
    if mod(ii,1000) == 0 ; disp(strcat({'Finding Eyeballs: '},num2str(ii),{' of '}, num2str(length(headRotMat_row_col_fr)))); end
    
    rEyeballCenterXYZ(ii,:) = headRotMat_row_col_fr(:,:,ii) * [rEye_ux; rEye_uy; rEye_uz]; %camera is location in head coordinate frame (where the origin is (headTopXYZ + HeadC1XYZ)/2)
    lEyeballCenterXYZ(ii,:) = headRotMat_row_col_fr(:,:,ii) * [lEye_ux; lEye_uy; lEye_uz]; %camera is location in head coordinate frame (where the origin is (headTopXYZ + HeadC1XYZ)/2)
end

rEyeballCenterXYZ = rEyeballCenterXYZ+hCen;
lEyeballCenterXYZ = lEyeballCenterXYZ+hCen;

worldCamCenterXYZ(:,1) = rEyeballCenterXYZ(:,1)*(1/3) + lEyeballCenterXYZ(:,1)*(1/3) + hTop(:,1)*(1/3);
worldCamCenterXYZ(:,2) = rEyeballCenterXYZ(:,2)*(1/3) + lEyeballCenterXYZ(:,2)*(1/3) + hTop(:,2)*(1/3);
worldCamCenterXYZ(:,3) = rEyeballCenterXYZ(:,3)*(1/3) + lEyeballCenterXYZ(:,3)*(1/3) + hTop(:,3)*(1/3);

% %% Debug Plot - play video of head orientation axes, cam position is black star
% 
% figure(2)
% sphRes = 20;
% r = ones(sphRes, sphRes);
% [th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
% [x1,y1,z1] = sph2cart(th, phi, r);
% 
% 
% lLeg = [2 3 4 5 6 7 5];
% rLeg = [2 8 9 10 11 12 10];
% tors = [2 13 14 15 26 27 28];
% lArm = [15 16 17 26 17 18 19 20];
% rArm = [15 21 22 26 22 23 24 25];
% 
% lToeID = 6;
% lToeXYZ = squeeze(shadow_fr_mar_dim(:,lToeID,:)); % pull out rHeelID marker000
% 
% for i = calibFrame:3:length(headZ)
%     
%     hx = hTop(i,1);
%     hy = hTop(i,2);
%     hz = hTop(i,3);
%     
%     plot3(shadow_fr_mar_dim(i,1:28,1),shadow_fr_mar_dim(i,1:28,2),shadow_fr_mar_dim(i,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
%     hold on
%     
%     
%     plot3(shadow_fr_mar_dim(i,lLeg,1),shadow_fr_mar_dim(i,lLeg,2),shadow_fr_mar_dim(i,lLeg,3),'b','LineWidth',2)
%     plot3(shadow_fr_mar_dim(i,rLeg,1),shadow_fr_mar_dim(i,rLeg,2),shadow_fr_mar_dim(i,rLeg,3),'r','LineWidth',2)
%     plot3(shadow_fr_mar_dim(i,tors,1),shadow_fr_mar_dim(i,tors,2),shadow_fr_mar_dim(i,tors,3),'g','LineWidth',2)
%     plot3(shadow_fr_mar_dim(i,lArm,1),shadow_fr_mar_dim(i,lArm,2),shadow_fr_mar_dim(i,lArm,3),'b','LineWidth',2)
%     plot3(shadow_fr_mar_dim(i,rArm,1),shadow_fr_mar_dim(i,rArm,2),shadow_fr_mar_dim(i,rArm,3),'r','LineWidth',2)
%     
%     
%     plot3([hx camXYZ(i,1)],[hy camXYZ(i,2)], [hz camXYZ(i,3)],'k-p');
%     plot3(camXYZ([i-100:i],1), camXYZ([i-100:i],2), camXYZ([i-100:i],3),'k-');
%     
%     plot3([hx headX(i,1)+hx],[hy headX(i,2)+hy], [hz headX(i,3)+hz],'r-');
%     plot3([hx headY(i,1)+hx],[hy headY(i,2)+hy], [hz headY(i,3)+hz],'g-');
%     plot3([hx headZ(i,1)+hx],[hy headZ(i,2)+hy], [hz headZ(i,3)+hz],'b-');
%         
%     plot3(calibPoints(:,1),calibPoints(:,2), calibPoints(:,3),'h')
%     
%     %     surface(x1+hx,y1+hy,z1+hz,'FaceColor', 'none','EdgeColor',[.9 .9 .9])
%       bx =   shadow_fr_mar_dim(i,1,1);
%       by =   shadow_fr_mar_dim(i,1,2);
%       bz =   shadow_fr_mar_dim(i,1,3);
%    
%     axis([-2000+bx 2000+bx -2000+by 2000+by -2500+bz 2000+bz])
%         
%     g_x = meshgrid(-1000+bx:100:1000+bx);
%     g_z = meshgrid(-1000+bz:100:1000+bz)';
%     grHeight(i) = lToeXYZ(i,2);  %the Y values denote vertical position
%     g_y = ones(length(g_x),length(g_x)) * grHeight(i);
%     surface(g_x, g_y, g_z,'FaceColor','w','EdgeColor',[.3 .3 .3])
%     
%     axis equal
%     title(num2str(i))
%     set(gca,'CameraUpVector',[0 1 0])
%     xlabel('x');ylabel('y'); zlabel('z');
%     view (-181, -52)
%     
%     hold off
%     drawnow
%     
% end













% Debug Plot

figure(9999393)
sphRes = 20;
r = ones(sphRes, sphRes);
[th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
[x1,y1,z1] = sph2cart(th, phi, r);


if shadowVersion == 2
    lLeg = [2 3 4 5 6 7 5];
    rLeg = [2 8 9 10 11 12 10];
    tors = [2 13 14 15 26 27 28];
    lArm = [15 16 17 26 17 18 19 20];
    rArm = [15 21 22 26 22 23 24 25];
numMarkers = 28;
elseif shadowVersion == 3    
    lLeg = [2 9 10 11 14 12 13 12 11];
    rLeg = [2 3 4  5  8  6  7  6  5];
    tors = [2 15 16 17 30 31 32];
    lArm = [17 24 25 30 25 26 27 28 29];
    rArm = [17 18 19 30 19 20 21 22 23];
    numMarkers = 32;
end

 ii = calibFrame;
    
    hx = hTop(ii,1);
    hy = hTop(ii,2);
    hz = hTop(ii,3);
    
    plot3(shadow_fr_mar_dim(ii,1:numMarkers,1),shadow_fr_mar_dim(ii,1:numMarkers,2),shadow_fr_mar_dim(ii,1:numMarkers,3),'ko','MarkerFaceColor','k','MarkerSize',4)
    hold on
    
    
    plot3(shadow_fr_mar_dim(ii,lLeg,1),shadow_fr_mar_dim(ii,lLeg,2),shadow_fr_mar_dim(ii,lLeg,3),'b','LineWidth',2)
    plot3(shadow_fr_mar_dim(ii,rLeg,1),shadow_fr_mar_dim(ii,rLeg,2),shadow_fr_mar_dim(ii,rLeg,3),'r','LineWidth',2)
    plot3(shadow_fr_mar_dim(ii,tors,1),shadow_fr_mar_dim(ii,tors,2),shadow_fr_mar_dim(ii,tors,3),'g','LineWidth',2)
    plot3(shadow_fr_mar_dim(ii,lArm,1),shadow_fr_mar_dim(ii,lArm,2),shadow_fr_mar_dim(ii,lArm,3),'b','LineWidth',2)
    plot3(shadow_fr_mar_dim(ii,rArm,1),shadow_fr_mar_dim(ii,rArm,2),shadow_fr_mar_dim(ii,rArm,3),'r','LineWidth',2)
    
    
        plot3(hx,hy, hz,'mh', 'MarkerSize', 20)
        
    plot3(hTop(ii,1), hTop(ii,2), hTop(ii,3),'gp', 'MarkerSize', 15)
    plot3(hC1(ii,1), hC1(ii,2), hC1(ii,3),'kp', 'MarkerSize', 15)

        

    
        plot3([hCen(ii,1) rEyeballCenterXYZ(ii,1)], [hCen(ii,2) rEyeballCenterXYZ(ii,2)], [hCen(ii,3) rEyeballCenterXYZ(ii,3)],'rp-')
        plot3([hCen(ii,1) lEyeballCenterXYZ(ii,1)], [hCen(ii,2) lEyeballCenterXYZ(ii,2)], [hCen(ii,3) lEyeballCenterXYZ(ii,3)],'bp-')
        plot3([hCen(ii,1) worldCamCenterXYZ(ii,1)], [hCen(ii,2) worldCamCenterXYZ(ii,2)], [hCen(ii,3) worldCamCenterXYZ(ii,3)],'kp-')
        

    %     surface(x1+hx,y1+hy,z1+hz,'FaceColor', 'none','EdgeColor',[.9 .9 .9])
      bx =   shadow_fr_mar_dim(ii,1,1);
      by =   shadow_fr_mar_dim(ii,1,2);
      bz =   shadow_fr_mar_dim(ii,1,3);
   
    axis([-2000+bx 2000+bx -2000+by 2000+by -2500+bz 2000+bz])
    
       view(-173, -43);
    axis equal
    title(num2str(ii))
    set(gca,'CameraUpVector',[0 1 0])
    xlabel('x');ylabel('y'); zlabel('z');
    
    hold off
    drawnow
    


