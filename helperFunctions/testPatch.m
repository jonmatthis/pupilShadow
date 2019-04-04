%% play a video of the VOR frames, if it comes to that

% first load a Rocks.mat / Woodchips.mat


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

headXCalib_XYZ_z = normr(headXCalib_XYZ_z);
headYCalib_XYZ_z = normr(headYCalib_XYZ_z);
headZCalib_XYZ_z = normr(headZCalib_XYZ_z);


debugVid = true;

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
        
        line([0 patchTopLeft(fr,1)]+hx, [0 patchTopLeft(fr,2)]+hy, [0 patchTopLeft(fr,3)]+hz,'color','k','linewidth',2);
        line([0 patchTopRight(fr,1)]+hx, [0 patchTopRight(fr,2)]+hy, [0 patchTopRight(fr,3)]+hz,'color','k','linewidth',2);
        line([0 patchBottomLeft(fr,1)]+hx, [0 patchBottomLeft(fr,2)]+hy, [0 patchBottomLeft(fr,3)]+hz,'color','k','linewidth',2);
        line([0 patchBottomRight(fr,1)]+hx, [0 patchBottomRight(fr,2)]+hy, [0 patchBottomRight(fr,3)]+hz,'color','k','linewidth',2);
 
        
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