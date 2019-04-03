function [ headVecX_fr_xyz, headVecY_fr_xyz, headVecZ_fr_xyz]  = findHeadVecs( headGlobalQuat_wxyz, shadow_fr_mar_dim, shadowMarkerNames, calibFrame, calibPoint, vorFrames)
%FIND HEAD VECS
%
% Find head orientation vectors
% X - pointing out the front of the head
% Y - pointing out the top of the head
% Z - pointing out the right of the head




%%
hTop =  squeeze(shadow_fr_mar_dim(:,strcmp('HeadTop', shadowMarkerNames),1:3)); %marker at Head Top
hC1 =   squeeze(shadow_fr_mar_dim(:,strcmp('Head', shadowMarkerNames),1:3)); %marker at C1 vertebra

lToe = squeeze(shadow_fr_mar_dim(:,strcmp('LeftToe', shadowMarkerNames),1:3));
rToe = squeeze(shadow_fr_mar_dim(:,strcmp('RightToe', shadowMarkerNames),1:3));



hCen = (hTop+hC1)/2;

hCenCalib = hCen(calibFrame,:);


headCalib_vec = normr(calibPoint- hCenCalib); %  use headCen - CalibPoint vec to tell which part of head is front
headYCalib_XYZ_z = normr(hTop(calibFrame,:)- hCenCalib); %  Head Y is just the top of the head
headZCalib_XYZ_z = cross(headCalib_vec,headYCalib_XYZ_z); % Head Z is coming out of the right of the head
headXCalib_XYZ_z = cross(headYCalib_XYZ_z,headZCalib_XYZ_z); %define head X based on head Y and Z

headVecScale = 1;
headXCalib_XYZ_z = headXCalib_XYZ_z *headVecScale;
headYCalib_XYZ_z = headYCalib_XYZ_z *headVecScale;
headZCalib_XYZ_z = headZCalib_XYZ_z *headVecScale;

%%%make sure head vectors are orthogonal
assert(dot(headXCalib_XYZ_z,headYCalib_XYZ_z) + dot(headXCalib_XYZ_z,headZCalib_XYZ_z) + dot(headZCalib_XYZ_z,headYCalib_XYZ_z) < eps*10, 'Your head vectors arent othogonal')
%% debug plot - show head vectors only
figure(543)
clf
debugPlot = true;
[x y z] = sphere;
mesh(x,y,z,'FaceColor','none')
hold on
plot3(0,0,0,'kp')
plot3([0 headXCalib_XYZ_z(1)], [0 headXCalib_XYZ_z(2)], [0 headXCalib_XYZ_z(3)],'r-')
plot3([0 headYCalib_XYZ_z(1)], [0 headYCalib_XYZ_z(2)], [0 headYCalib_XYZ_z(3)],'g-')
plot3([0 headZCalib_XYZ_z(1)], [0 headZCalib_XYZ_z(2)], [0 headZCalib_XYZ_z(3)],'b-')
axis equal
set(gca,'CameraUpVector',[0 1 0])
xlabel('x');ylabel('y'); zlabel('z');

%%

% find unit vectors f head reference frame
headRotMat_row_col_fr = headGlobalQuat_wxyz.RotationMatrix;

calibRotMat = headRotMat_row_col_fr(:,:,calibFrame);


% find unit vectors in head reference frame
headRotMat_row_col_fr = headGlobalQuat_wxyz.RotationMatrix;

calibRotMat = headRotMat_row_col_fr(:,:,calibFrame);

ux = calibRotMat * [1; 0; 0];
uy = calibRotMat * [0; 1; 0];
uz = calibRotMat * [0; 0; 1];

%defining these vectors relative to unit vectors allows you to rotate them
%using the headRotMat.... or something. I really need to learn how math
%works... :-/
headX_ux = dot(headXCalib_XYZ_z, ux);
headX_uy = dot(headXCalib_XYZ_z, uy);
headX_uz = dot(headXCalib_XYZ_z, uz);

headY_ux = dot(headYCalib_XYZ_z, ux);
headY_uy = dot(headYCalib_XYZ_z, uy);
headY_uz = dot(headYCalib_XYZ_z, uz);

headZ_ux = dot(headZCalib_XYZ_z, ux);
headZ_uy = dot(headZCalib_XYZ_z, uy);
headZ_uz = dot(headZCalib_XYZ_z, uz);


debug = true;
if debug
    % Debug Plot
    
    figure(5432)
    sphRes = 20;
    r = ones(sphRes, sphRes);
    [th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
    [x1,y1,z1] = sph2cart(th, phi, r);
    
    
    lLeg = [2 3 4 5 6 7 5];
    rLeg = [2 8 9 10 11 12 10];
    tors = [2 13 14 15 26 27 28];
    lArm = [15 16 17 26 17 18 19 20];
    rArm = [15 21 22 26 22 23 24 25];
    
    i = calibFrame;
    
    hx = hCen(i,1);
    hy = hCen(i,2);
    hz = hCen(i,3);
    
    plot3(shadow_fr_mar_dim(i,1:28,1),shadow_fr_mar_dim(i,1:28,2),shadow_fr_mar_dim(i,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
    hold on
    
    
    plot3(shadow_fr_mar_dim(i,lLeg,1),shadow_fr_mar_dim(i,lLeg,2),shadow_fr_mar_dim(i,lLeg,3),'b','LineWidth',2)
    plot3(shadow_fr_mar_dim(i,rLeg,1),shadow_fr_mar_dim(i,rLeg,2),shadow_fr_mar_dim(i,rLeg,3),'r','LineWidth',2)
    plot3(shadow_fr_mar_dim(i,tors,1),shadow_fr_mar_dim(i,tors,2),shadow_fr_mar_dim(i,tors,3),'g','LineWidth',2)
    plot3(shadow_fr_mar_dim(i,lArm,1),shadow_fr_mar_dim(i,lArm,2),shadow_fr_mar_dim(i,lArm,3),'b','LineWidth',2)
    plot3(shadow_fr_mar_dim(i,rArm,1),shadow_fr_mar_dim(i,rArm,2),shadow_fr_mar_dim(i,rArm,3),'r','LineWidth',2)
    
    
    plot3(hx,hy, hz,'mh', 'MarkerSize', 20)
    
    plot3(hTop(i,1), hTop(i,2), hTop(i,3),'gp', 'MarkerSize', 15)
    plot3(hC1(i,1), hC1(i,2), hC1(i,3),'kp', 'MarkerSize', 15)
    
    
    %     plot3([hx hux(1)+hx], [hy hux(2)+hy], [hz hux(3)+hz],'rp-', 'MarkerSize', 15)
    %     plot3([hx huy(1)+hx], [hy huy(2)+hy], [hz huy(3)+hz],'gp-', 'MarkerSize', 15)
    %     plot3([hx huz(1)+hx], [hy huz(2)+hy], [hz huz(3)+hz],'bp-', 'MarkerSize', 15)
    
    
    plot3([hx headXCalib_XYZ_z(1)+hx], [hy headXCalib_XYZ_z(2)+hy], [hz headXCalib_XYZ_z(3)+hz],'r-o','LineWidth',2)
    plot3([hx headYCalib_XYZ_z(1)+hx], [hy headYCalib_XYZ_z(2)+hy], [hz headYCalib_XYZ_z(3)+hz],'g-o','LineWidth',2)
    plot3([hx headZCalib_XYZ_z(1)+hx], [hy headZCalib_XYZ_z(2)+hy], [hz headZCalib_XYZ_z(3)+hz],'b-o','LineWidth',2)
    
    plot3(calibPoint(1),calibPoint(2), calibPoint(3),'h')
    
    %     surface(x1+hx,y1+hy,z1+hz,'FaceColor', 'none','EdgeColor',[.9 .9 .9])
    bx =   shadow_fr_mar_dim(i,1,1);
    by =   shadow_fr_mar_dim(i,1,2);
    bz =   shadow_fr_mar_dim(i,1,3);
    
    [xx,zz] = meshgrid(-2000+bx:100:2000+bx, -2500+bz:100:2000+bz);
    surface(xx,zeros(size(xx)),zz,'FaceColor','w')
    axis([-2000+bx 2000+bx -2000+by 2000+by -2500+bz 2000+bz])
    
    
    axis equal
    title(['Head Vectors on CalibFrame ' num2str(i)])
    set(gca,'CameraUpVector',[0 1 0])
    xlabel('x');ylabel('y'); zlabel('z');
    
    hold off
    drawnow
    
end


%%



%make a head vector from IMU rotation matrix
headVecX_fr_xyz = nan(length(headRotMat_row_col_fr),3);
headVecY_fr_xyz = nan(length(headRotMat_row_col_fr),3);
headVecZ_fr_xyz = nan(length(headRotMat_row_col_fr),3);

for fr = 2:length(headRotMat_row_col_fr)
    
    headVecX_fr_xyz(fr,:) = headRotMat_row_col_fr(:,:,fr) * [headX_ux; headX_uy; headX_uz];
    headVecY_fr_xyz(fr,:) = headRotMat_row_col_fr(:,:,fr) * [headY_ux; headY_uy; headY_uz];
    headVecZ_fr_xyz(fr,:) = headRotMat_row_col_fr(:,:,fr) * [headZ_ux; headZ_uy; headZ_uz];
    
end

% %put headVecs in a world-centered reference frame.
% headVecX_fr_xyz = headVecX_fr_xyz+hCen;
% headVecY_fr_xyz = headVecY_fr_xyz+hCen;
% headVecZ_fr_xyz = headVecZ_fr_xyz+hCen;

% Debug Plot

figure(9999393)
clf
sphRes = 20;
r = ones(sphRes, sphRes);
[th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
[x1,y1,z1] = sph2cart(th, phi, r);


lLeg = [2 3 4 5 6 7 5];
rLeg = [2 8 9 10 11 12 10];
tors = [2 13 14 15 26 27 28];
lArm = [15 16 17 26 17 18 19 20];
rArm = [15 21 22 26 22 23 24 25];

fr = calibFrame;

hx = hCen(fr,1);
hy = hCen(fr,2);
hz = hCen(fr,3);

plot3(shadow_fr_mar_dim(fr,1:28,1),shadow_fr_mar_dim(fr,1:28,2),shadow_fr_mar_dim(fr,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
hold on


plot3(shadow_fr_mar_dim(fr,lLeg,1),shadow_fr_mar_dim(fr,lLeg,2),shadow_fr_mar_dim(fr,lLeg,3),'b','LineWidth',2)
plot3(shadow_fr_mar_dim(fr,rLeg,1),shadow_fr_mar_dim(fr,rLeg,2),shadow_fr_mar_dim(fr,rLeg,3),'r','LineWidth',2)
plot3(shadow_fr_mar_dim(fr,tors,1),shadow_fr_mar_dim(fr,tors,2),shadow_fr_mar_dim(fr,tors,3),'g','LineWidth',2)
plot3(shadow_fr_mar_dim(fr,lArm,1),shadow_fr_mar_dim(fr,lArm,2),shadow_fr_mar_dim(fr,lArm,3),'b','LineWidth',2)
plot3(shadow_fr_mar_dim(fr,rArm,1),shadow_fr_mar_dim(fr,rArm,2),shadow_fr_mar_dim(fr,rArm,3),'r','LineWidth',2)


plot3(hx,hy, hz,'mh', 'MarkerSize', 20)

plot3(hTop(fr,1), hTop(fr,2), hTop(fr,3),'gp', 'MarkerSize', 15)
plot3(hC1(fr,1), hC1(fr,2), hC1(fr,3),'kp', 'MarkerSize', 15)



headVecScale = 1000;

plot3([0 headVecX_fr_xyz(fr,1)*headVecScale]+hx, [0 headVecX_fr_xyz(fr,2)*headVecScale]+hy, [0 headVecX_fr_xyz(fr,3)*headVecScale]+hz,'r-o','LineWidth',2)
plot3([0 headVecY_fr_xyz(fr,1)*headVecScale]+hx, [0 headVecY_fr_xyz(fr,2)*headVecScale]+hy, [0 headVecY_fr_xyz(fr,3)*headVecScale]+hz,'g-o','LineWidth',2)
plot3([0 headVecZ_fr_xyz(fr,1)*headVecScale]+hx, [0 headVecZ_fr_xyz(fr,2)*headVecScale]+hy, [0 headVecZ_fr_xyz(fr,3)*headVecScale]+hz,'b-o','LineWidth',2)


plot3(calibPoint(1), calibPoint(2), calibPoint(3), 'kp')
%     surface(x1+hx,y1+hy,z1+hz,'FaceColor', 'none','EdgeColor',[.9 .9 .9])
bx =   shadow_fr_mar_dim(fr,1,1);
by =   shadow_fr_mar_dim(fr,1,2);
bz =   shadow_fr_mar_dim(fr,1,3);

[xx,zz] = meshgrid(-2000+bx:100:2000+bx, -2500+bz:100:2000+bz);
surface(xx,zeros(size(xx)),zz,'FaceColor','w')

axis([-2000+bx 2000+bx -2000+by 2000+by -2500+bz 2000+bz])

view(-173, -43);
axis equal
title(num2str(fr))
set(gca,'CameraUpVector',[0 1 0])
xlabel('x');ylabel('y'); zlabel('z');

hold off
drawnow

%% play a video of the VOR frames, if it comes to that
debugVid = false;

if debugVid
    
    figure(9999393)
    clf
    sphRes = 20;
    r = ones(sphRes, sphRes);
    [th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
    [x1,y1,z1] = sph2cart(th, phi, r);
    
    
    lLeg = [2 3 4 5 6 7 5];
    rLeg = [2 8 9 10 11 12 10];
    tors = [2 13 14 15 26 27 28];
    lArm = [15 16 17 26 17 18 19 20];
    rArm = [15 21 22 26 22 23 24 25];
    
    for fr = vorFrames(1):10:vorFrames(end)
        
        hx = hCen(fr,1);
        hy = hCen(fr,2);
        hz = hCen(fr,3);
        
        plot3(shadow_fr_mar_dim(fr,1:28,1),shadow_fr_mar_dim(fr,1:28,2),shadow_fr_mar_dim(fr,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
        hold on
        
        
        plot3(shadow_fr_mar_dim(fr,lLeg,1),shadow_fr_mar_dim(fr,lLeg,2),shadow_fr_mar_dim(fr,lLeg,3),'b','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,rLeg,1),shadow_fr_mar_dim(fr,rLeg,2),shadow_fr_mar_dim(fr,rLeg,3),'r','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,tors,1),shadow_fr_mar_dim(fr,tors,2),shadow_fr_mar_dim(fr,tors,3),'g','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,lArm,1),shadow_fr_mar_dim(fr,lArm,2),shadow_fr_mar_dim(fr,lArm,3),'b','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,rArm,1),shadow_fr_mar_dim(fr,rArm,2),shadow_fr_mar_dim(fr,rArm,3),'r','LineWidth',2)
        
        
        plot3(hx,hy, hz,'mh', 'MarkerSize', 20)
        
        plot3(hTop(fr,1), hTop(fr,2), hTop(fr,3),'gp', 'MarkerSize', 15)
        plot3(hC1(fr,1), hC1(fr,2), hC1(fr,3),'kp', 'MarkerSize', 15)
        
        
        
        headVecScale = 1000;
        
        plot3([0 headVecX_fr_xyz(fr,1)*headVecScale]+hx, [0 headVecX_fr_xyz(fr,2)*headVecScale]+hy, [0 headVecX_fr_xyz(fr,3)*headVecScale]+hz,'r-o','LineWidth',2)
        plot3([0 headVecY_fr_xyz(fr,1)*headVecScale]+hx, [0 headVecY_fr_xyz(fr,2)*headVecScale]+hy, [0 headVecY_fr_xyz(fr,3)*headVecScale]+hz,'g-o','LineWidth',2)
        plot3([0 headVecZ_fr_xyz(fr,1)*headVecScale]+hx, [0 headVecZ_fr_xyz(fr,2)*headVecScale]+hy, [0 headVecZ_fr_xyz(fr,3)*headVecScale]+hz,'b-o','LineWidth',2)
        
        plot3(calibPoint(1), calibPoint(2), calibPoint(3), 'kp')
        %     surface(x1+hx,y1+hy,z1+hz,'FaceColor', 'none','EdgeColor',[.9 .9 .9])
        bx =   shadow_fr_mar_dim(fr,1,1);
        by =   shadow_fr_mar_dim(fr,1,2);
        bz =   shadow_fr_mar_dim(fr,1,3);
        
        [xx,zz] = meshgrid(-2000+bx:100:2000+bx, -2500+bz:100:2000+bz);
        surface(xx,zeros(size(xx)),zz,'FaceColor','w')
        
        axis([-2000+bx 2000+bx -2000+by 2000+by -2500+bz 2000+bz])
        
        view(-173, -43);
        axis equal
        title(num2str(fr))
        set(gca,'CameraUpVector',[0 1 0])
        xlabel('x');ylabel('y'); zlabel('z');
        
        hold off
        drawnow
    end
end


