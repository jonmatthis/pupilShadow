function playLaserFaceSkeletonMonster(w)

shadow_fr_mar_dim = w.shadow_fr_mar_dim;
calibDist = w.calibDist;
comXYZ = w.comXYZ;
rGazeGroundIntersection = w.rGazeGroundIntersection;
lGazeGroundIntersection = w.lGazeGroundIntersection;

rEyeballCenterXYZ = w.rEyeballCenterXYZ;
lEyeballCenterXYZ = w.lEyeballCenterXYZ;


rHeelXYZ = w.rHeelXYZ;
lHeelXYZ = w.lHeelXYZ;


rGazeXYZ = w.rGazeXYZ;
lGazeXYZ = w.lGazeXYZ;

headVecX_fr_xyz = w.headVecX_fr_xyz;
headVecY_fr_xyz = w.headVecY_fr_xyz;
headVecZ_fr_xyz = w.headVecZ_fr_xyz;

steps_HS_TO_StanceLeg_XYZ = w.steps_HS_TO_StanceLeg_XYZ;
%% %%% make sphere thingy fr eyeball guys
sphRes = 20;
r = 35;%mean(rEye.sphere_radius); %p.s. it's 12mm, but let's blow 'em up a bit for ... visibilitiy... 8D
[th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
[x1,y1,z1] = sph2cart(th, phi, r);

normScale = calibDist;
plotSkel = true;

lLeg = [2 3 4 5 6 7 5];
rLeg = [2 8 9 10 11 12 10];
tors = [2 13 14 15 26 27 28];
lArm = [15 16 17 26 17 18 19 20];
rArm = [15 21 22 26 22 23 24 25];

comXYZ = squeeze(shadow_fr_mar_dim(:,1,:));


frames = 1:10:length(comXYZ);
% plot (hypothetical) groundplane

xSpan = [min(rGazeGroundIntersection(frames,1))-5000, max(rGazeGroundIntersection(frames,1))+5000];
zSpan = [min(rGazeGroundIntersection(frames,3))-5000, max(rGazeGroundIntersection(frames,3))+5000];


res      = 100; % resultion for the meshgrid
[groundPlane_x, groundPlane_z] = meshgrid(xSpan(1):res:xSpan(2), zSpan(1):res:zSpan(2));


groundPlane_y = ones(size(groundPlane_x));
groundPlane_color = ones(size(groundPlane_x));




figure(1254);clf
% set(gcf,'Position',[1921 121 1920 979])



for ii = frames
    ii
    cla
    clf
    %%eyeball centers in shadow coordinats(not to be confused with "rEye_sphCen_x", which are in pupil camera coords)
    rEx = rEyeballCenterXYZ(ii,1);
    rEy = rEyeballCenterXYZ(ii,2);
    rEz = rEyeballCenterXYZ(ii,3);
    
    lEx = lEyeballCenterXYZ(ii,1);
    lEy = lEyeballCenterXYZ(ii,2);
    lEz = lEyeballCenterXYZ(ii,3);
    
    %     %%pull out the l and r eye sphere centers for this frame
    %     rCx = rEye_sphCenCam_x(ii);
    %     rCy = rEye_sphCenCam_y(ii);
    %     rCz = rEye_sphCenCam_z(ii);
    %
    %     lCx = lEye_sphCenCam_x(ii);
    %     lCy = lEye_sphCenCam_y(ii);
    %     lCz = lEye_sphCenCam_z(ii);
    %
    
    grHeight(ii) = min([rHeelXYZ(ii,2) lHeelXYZ(ii,2) ]);
    
    % right eye
    r1 =  mesh(x1+rEx, y1+rEy, z1+rEz);
    r1.FaceColor = [1 .9 .9];
    r1.EdgeColor = 'k';
    r1.EdgeAlpha = 0.1;
    hold on
    
    
    %     %%% Plot circular patch for pupil - centered on pupilNorm (code jacked from - https://www.mathworks.com/matlabcentral/fileexchange/26588-plot-circle-in-3d)
    %     thisRPupCenter = [rEye_pupCircCenXYZ(ii,1)-rCx rEye_pupCircCenXYZ(ii,2)-rCy rEye_pupCircCenXYZ(ii,3)-rCz] ;
    %     thisRPupNormal = thisRPupCenter*normScale;
    %     thisRPupRadius = rEye_pupRadius(ii);
    %
    %     if ~isnan(thisRPupNormal)
    %         theta=0:.1:2*pi;
    %         v=null(thisRPupNormal);
    %         points=repmat(thisRPupCenter',1,size(theta,2))+thisRPupRadius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    %         patch(points(1,:)+rEx, points(2,:)+rEy, points(3,:)+rEz ,'r');
    %     end
    %%%%
    %
    %     plot3([0+rEx thisRPupCenter(1)+rEx],...
    %         [0+rEy thisRPupCenter(2)+rEy],...
    %         [0+rEz thisRPupCenter(3)+rEz],'k-','LineWidth',2)
    %
    %
    %     plot3([rEx thisRPupNormal(1)+rEx],...
    %         [rEy thisRPupNormal(2)+rEy],...
    %         [rEz thisRPupNormal(3)+rEz],'m-')
    %
    %     plot3([ thisRPupNormal(1)*normScale+rEx],...
    %         [ thisRPupNormal(2)*normScale+rEy],...
    %         [r thisRPupNormal(3)*normScale+rEz],'kp')
    
    plot3([rEx rGazeXYZ(ii,1)],...
        [rEy rGazeXYZ(ii,2)],...
        [rEz rGazeXYZ(ii,3)], 'm-','LineWidth',2)
    
    plot3(rGazeGroundIntersection(ii,1),...
        rGazeGroundIntersection(ii,2),...
        rGazeGroundIntersection(ii,3),'kp','MarkerFaceColor','r','MarkerSize',12)
    
    %     plot3(rGazeGroundIntersection(frames,1),...
    %         rGazeGroundIntersection(frames,2),...
    %         rGazeGroundIntersection(frames,3),'-r')
    
    
    % left eye
    l1 =  mesh(x1+lEx, y1+lEy, z1+lEz);
    l1.FaceColor = [.9 .9 1];
    l1.EdgeColor = 'none';
    
    hold on
    
    
    %     %%% Plot circular patch for pupil - centered on pupilNorm (code jacked from - https://www.mathworks.com/matlabcentral/fileexchange/26588-plot-circle-in-3d)
    %     thisLPupCenter = [lEye_pupCircCenXYZ(ii,1)-lCx lEye_pupCircCenXYZ(ii,2)-lCy lEye_pupCircCenXYZ(ii,3)-lCz] ;
    %     thisLPupNormal = thisLPupCenter*1.3;
    %     thisLPupRadius = lEye_pupRadius(ii);
    %
    %     if ~isnan(thisLPupNormal)
    %         theta=0:.1:2*pi;
    %         v=null(thisLPupNormal);
    %         points=repmat(thisLPupCenter',1,size(theta,2))+thisLPupRadius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    %         patch(points(1,:)+lEx, points(2,:)+lEy, points(3,:)+lEz ,'b');
    %     end
    %     %%%%
    
    %     plot3([0+lEx thisLPupCenter(1)+lEx],...
    %         [0+lEy thisLPupCenter(2)+lEy],...
    %         [0+lEz thisLPupCenter(3)+lEz],'k-','LineWidth',2)
    %
    %
    %     plot3([lEx thisLPupNormal(1)*normScale+lEx],...
    %         [lEy thisLPupNormal(2)*normScale+lEy],...
    %         [lEz thisLPupNormal(3)*normScale+lEz],'c-','LineWidth',2)
    %
    %     plot3([thisLPupNormal(1)*normScale+lEx],...
    %         [thisLPupNormal(2)*normScale+lEy],...
    %         [thisLPupNormal(3)*normScale+lEz],'kp')
    
    plot3([lEx lGazeXYZ(ii,1)],...
        [lEy lGazeXYZ(ii,2)],...
        [lEz lGazeXYZ(ii,3)], 'c-','LineWidth',2)
    
    plot3(lGazeGroundIntersection(ii,1),...
        lGazeGroundIntersection(ii,2),...
        lGazeGroundIntersection(ii,3),'kp','MarkerFaceColor','b','MarkerSize',12)
    
    %     plot3(lGazeGroundIntersection(frames,1),...
    %         lGazeGroundIntersection(frames,2),...
    %         lGazeGroundIntersection(frames,3),'-b')
    
    
    if plotSkel
        %%%Plotcherself a skeleetoon
        plot3(shadow_fr_mar_dim(ii,1:28,1),shadow_fr_mar_dim(ii,1:28,2),shadow_fr_mar_dim(ii,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
        hold on
        
        
        plot3(shadow_fr_mar_dim(ii,lLeg,1),shadow_fr_mar_dim(ii,lLeg,2),shadow_fr_mar_dim(ii,lLeg,3),'b','LineWidth',2)
        plot3(shadow_fr_mar_dim(ii,rLeg,1),shadow_fr_mar_dim(ii,rLeg,2),shadow_fr_mar_dim(ii,rLeg,3),'r','LineWidth',2)
        plot3(shadow_fr_mar_dim(ii,tors,1),shadow_fr_mar_dim(ii,tors,2),shadow_fr_mar_dim(ii,tors,3),'g','LineWidth',2)
        plot3(shadow_fr_mar_dim(ii,lArm,1),shadow_fr_mar_dim(ii,lArm,2),shadow_fr_mar_dim(ii,lArm,3),'b','LineWidth',2)
        plot3(shadow_fr_mar_dim(ii,rArm,1),shadow_fr_mar_dim(ii,rArm,2),shadow_fr_mar_dim(ii,rArm,3),'r','LineWidth',2)
        
        %plot head axes
        hx = shadow_fr_mar_dim(ii,28,1);
        hy = shadow_fr_mar_dim(ii,28,2);
        hz = shadow_fr_mar_dim(ii,28,3);
        
        plot3([ hx headVecX_fr_xyz(ii,1)*1000+hx], [hy headVecX_fr_xyz(ii,2)*1000+hy],[hz headVecX_fr_xyz(ii,3)*1000+hz],'r-','LineWidth',3)
        plot3([ hx headVecY_fr_xyz(ii,1)*1000+hx], [hy headVecY_fr_xyz(ii,2)*1000+hy],[hz headVecY_fr_xyz(ii,3)*1000+hz],'g-','LineWidth',3)
        plot3([ hx headVecZ_fr_xyz(ii,1)*1000+hx], [hy headVecZ_fr_xyz(ii,2)*1000+hy],[hz headVecZ_fr_xyz(ii,3)*1000+hz],'b-','LineWidth',3)
        
        bx =   shadow_fr_mar_dim(ii,1,1);
        by =   shadow_fr_mar_dim(ii,1,2);
        bz =   shadow_fr_mar_dim(ii,1,3);
        
        %%% plot foothold locations
        rFootholds = steps_HS_TO_StanceLeg_XYZ(steps_HS_TO_StanceLeg_XYZ(:,3) == 1 ,:);
        lFootholds = steps_HS_TO_StanceLeg_XYZ(steps_HS_TO_StanceLeg_XYZ(:,3) == 2 ,:);
        
        rFootholds(rFootholds(:,1)<ii-2000 | rFootholds(:,1)>ii+2000,:) = [];
        lFootholds(lFootholds(:,1)<ii-2000 | lFootholds(:,1)>ii+2000,:) = [];
        
        %   plot vertical projection of foothold locations onto groundplane
        
        plot3(rFootholds(:,4), ones(length(rFootholds(:,1)))*grHeight(ii), rFootholds(:,6),'ko','MarkerSize', 9, 'MarkerFaceColor','r')
        plot3(lFootholds(:,4), ones(length(lFootholds(:,1)))*grHeight(ii), lFootholds(:,6),'ko','MarkerSize', 9, 'MarkerFaceColor','c')
        
        
        %plot gaussianly burnt groundplane
        sigma = 7500;
        meanGazeGround = mean([lGazeGroundIntersection(ii,:); rGazeGroundIntersection(ii,:)]);
        gaussian        = 1./sqrt(2*pi*sigma).*exp(-1./(2*sigma).*( (groundPlane_z-meanGazeGround(3)).^2 + (groundPlane_x-meanGazeGround(1)).^2));
        gaussianNorm    = gaussian ./ max(max(gaussian));
        
        if ~isnan(gaussianNorm)
            groundPlane_color = groundPlane_color + gaussianNorm; %add 2d gaussian for this frame's gaze/ground intersection ground plane
        end
        
        %         g_x = meshgrid(-10e4:500:10e4) + comXYZ(ii,1);
        %         g_y = ones(size(g_x)) * min([rHeelXYZ(ii,2) lHeelXYZ(ii,2) ]);
        %         g_z = meshgrid(-10e4:500:10e4)' + comXYZ(ii,3);
        
        s1 = surface(groundPlane_x , groundPlane_y*grHeight(ii), groundPlane_z, groundPlane_color  );
        s1.LineStyle = 'none';
        s1.FaceColor = 'interp';
        
        %         CT = cbrewer('div', 'Spectral', 64);
        %         colormap(flipud(CT));
        colormap jet
        caxis([0 5])
        
    end
    %     view(-173, -43);
    axis equal
    title(num2str(ii))
    %     set(gca,'CameraUpVector',[0 1 0])
    xlabel('x');ylabel('y'); zlabel('z');
    %     axis([-5000+bx 5000+bx -5000+by 5000+by -5000+bz 5000+bz])
    
    a = gca;
    a.CameraTarget = [comXYZ(ii,1), comXYZ(ii,2), comXYZ(ii,3)]; %point figure 'camera' at COM
    a.CameraPosition = a.CameraTarget + [-1800 1800 2000]; %set camera position
    a.CameraViewAngle = 80;
    a.CameraUpVector = [ 0 1 0];
    a.Position = [0 0 1 1];
    
    hold off
    drawnow
    
end

