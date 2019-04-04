function [patchTopLeft patchBottomLeft patchBottomRight patchTopRight] = getCamFrustum(headVecX_fr_xyz,headVecY_fr_xyz,headVecZ_fr_xyz,gaze_norm_pos_x,gaze_norm_pos_y,...
    px2mmScale,calibDist,rGazeXYZ,lGazeXYZ,rEyeballCenterXYZ,lEyeballCenterXYZ,resHeight,resWidth,shadow_fr_mar_dim, calibFrame, calibPoint, vorFrames,shadowMarkerNames)

porX = gaze_norm_pos_x*resWidth;
porY = gaze_norm_pos_y*resHeight;

%% head bases

bx = normr(headVecX_fr_xyz);
by = normr(headVecY_fr_xyz);
bz = normr(headVecZ_fr_xyz);


 %% find cam orientation
        
        p2m = px2mmScale;
        calibDist = calibDist;
        
        % unaligned 3D gaze vector based on norm_pos X/Y
        unalign_vec = [p2m*(porX-resWidth/2)...
            p2m*(porY-resHeight/2)...
            -repmat(calibDist,[length(porX) 1])];
        
        
        % unit vectors
        unalign_vec = normr(unalign_vec);
        
        % vectors that we want to align thoe to (calculated from 3D eye
        % balls)
        align_vec_R = rGazeXYZ-rEyeballCenterXYZ;
        align_vec_L = lGazeXYZ-lEyeballCenterXYZ;
        
        % average these for a cylcopean gaze vector
        align_vec = 0.5*(align_vec_R+align_vec_L);
        
        % normalize
        align_vec = normr(align_vec);
        
        % randexer
        dex = randperm(size(unalign_vec,1),1000);
        
        % set up sets of points we want to align
        p1 = unalign_vec;
        p2 = align_vec;
        
        % for each head orientation, inverse rotate the 3D eyeball vectors
        % so that they are now in head coordinates
        for ii = 1:length(bx)
            
            R = [bx(ii,:)' by(ii,:)' bz(ii,:)'];
            
            p2(ii,:) = (inv(R)*p2(ii,:)')';
            
        end
        
        p1 = p1(dex,:);
        p2 = p2(dex,:);
        
        % optimize to search for rotation angles
        save('gpp.mat','p1','p2');
        x0 = zeros(3,1);
        x = fminsearch(@errFunPoints,x0);
        
        % camera to head rotation matrix
        R2 = eul2rotm(x');
        %%
        clear frust frust_pre
        
        % frustum (previously patch)
        frust_pre = [-resWidth/2*p2m -resHeight/2*p2m -calibDist;...
            -resWidth/2*p2m resHeight/2*p2m -calibDist;...
            resWidth/2*p2m resHeight/2*p2m -calibDist;...
            resWidth/2*p2m -resHeight/2*p2m -calibDist];
        
        % figure out frustum for each frame by first rotation by the camera
        % to head rotation matrix, then the particular frames head rotation
        % matrix (frustum in world now)
        for ii = 1:length(bx)
            R1 = [bx(ii,:)' by(ii,:)' bz(ii,:)'];
            
            frust(:,:,ii) = (R1*R2*frust_pre')';
        end

        patchTopLeft = squeeze(frust(1,:,:))';
        patchBottomLeft = squeeze(frust(2,:,:))';
        patchBottomRight = squeeze(frust(3,:,:))';
        patchTopRight = squeeze(frust(4,:,:))';
        
    %% play a video of the VOR frames, if it comes to that

% first load a Rocks.mat / Woodchips.mat


hTop =  squeeze(shadow_fr_mar_dim(:,strcmp('HeadTop', shadowMarkerNames),1:3)); %marker at Head Top
hC1 =   squeeze(shadow_fr_mar_dim(:,strcmp('Head', shadowMarkerNames),1:3)); %marker at C1 vertebra

lToe = squeeze(shadow_fr_mar_dim(:,strcmp('LeftToe', shadowMarkerNames),1:3));
rToe = squeeze(shadow_fr_mar_dim(:,strcmp('RightToe', shadowMarkerNames),1:3));



hCen = (hTop+hC1)/2;

hCenCalib = hCen(calibFrame,:);


headCalib_vec = normr(calibPoint- hCenCalib); %  use headCen - CalibPoint vec to tell which part of head is front



debugVid = true;

if debugVid
    
    figure(91836)
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
    
    for fr = vorFrames(1)

%     for fr = vorFrames(1):10:vorFrames(end)
        
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



end

