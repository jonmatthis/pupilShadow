function [camAlignErr] = vorPupilAlignErrFun_worldCam(vData, camAlignGuess_xyz)
%VORPUPILALIGNERRFUN Summary of this function goes here
%   Detailed explanation goes here

camQuatGuess = quaternion.eulerangles('123',camAlignGuess_xyz(1),camAlignGuess_xyz(2),camAlignGuess_xyz(3));
camRotMatGuess = camQuatGuess.RotationMatrix;

%%unpack vData
headRotMat_row_col_fr = vData.headRotMat_row_col_fr;

worldCamCenterXYZ = vData.worldCamCenterXYZ;

gazeXYZ = vData.gazeXYZ;

rHeelXYZ = vData.rHeelXYZ;
lHeelXYZ = vData.lHeelXYZ;

patchXYZ = vData.patchXYZ;

patchX =  vData.patchXYZ(1,:);
patchY =  vData.patchXYZ(2,:);
patchZ =  vData.patchXYZ(3,:);

patchNormX =   mean(patchX(1:4));
patchNormY =   mean(patchY(1:4));
patchNormZ =   mean(patchZ(1:4));


calibPoint = vData.calibPoint;

confidence = vData.confidence;

%normalize its length
for ll = 1:length(gazeXYZ)
    gazeXYZ(ll,:) = gazeXYZ(ll,:)/norm(gazeXYZ(ll,:));
end

%multiply it by your desired length ;)
calibDist = pdist([worldCamCenterXYZ(1,:); calibPoint]); %myboy pythag
gazeXYZ = gazeXYZ*calibDist;
%%%%
%%%%%%% This part's important - Rotate gaze vector by [this guess at the proper alignment matrix], prior to resituating  the origin on on the eyeball
%%%%


for rr = 1:length(gazeXYZ)
    
    thisET_frame_unrot = camRotMatGuess * [gazeXYZ(rr,1); gazeXYZ(rr,2); gazeXYZ(rr,3)];
    thisETframe = headRotMat_row_col_fr(:,:,rr) * thisET_frame_unrot;
    
    headOrVec(rr,:) =     headRotMat_row_col_fr(:,:,rr) * [2e3; 0;0];
    
    gazeXYZ(rr,:) = thisETframe;
    
    
    %rotate the vidFramePatch
    for j = 1:5 %it takes 5 vectors to define a square patch
        thisPatchUnrot(:,j) =  (camRotMatGuess * patchXYZ(:,j));
        thisPatch(:,j) = headRotMat_row_col_fr(:,:,rr)  * thisPatchUnrot(:,j);
    end
    
    patchX(rr,:) = thisPatch(1,:);
    patchY(rr,:) = thisPatch(2,:);
    patchZ(rr,:) = thisPatch(3,:);
end




% add the eyeball center (in shadow/world coordiates) to translate origin of gaze vector onto the shadow eyeball
gazeXYZ(:,1) = gazeXYZ(:,1)+ worldCamCenterXYZ(:,1);
gazeXYZ(:,2) = gazeXYZ(:,2)+ worldCamCenterXYZ(:,2);
gazeXYZ(:,3) = gazeXYZ(:,3)+ worldCamCenterXYZ(:,3);

gazeXYZ(abs(sum(diff(gazeXYZ)'))>40,:) = nan;

confThresh = .8;
gazeXYZ(confidence < confThresh,:) = nan;

% [ gazeGroundIntersectionXYZ] = calcGroundFixations( rHeelXYZ, lHeelXYZ, gazeXYZ, eyeballCenterXYZ );

% gazeXYZ = gazeXYZ + worldCamCenterXYZ;

for cc = 1:length(gazeXYZ)
    % error for each frame == distance between gaze/ground intersection and the calibration point
    err(cc) = sqrt( (calibPoint(1) - gazeXYZ(cc,1))^2 + (calibPoint(2) - gazeXYZ(cc,2))^2 +(calibPoint(3) - gazeXYZ(cc,3))^2 );
%     err(cc) = sqrt( (calibPoint(1) - gazeGroundIntersectionXYZ(cc,1))^2 + (calibPoint(2) - gazeGroundIntersectionXYZ(cc,2))^2 +(calibPoint(3) - gazeGroundIntersectionXYZ(cc,3))^2 );

end


camAlignErr = nansum(err)/length(err);

if vData.plotDebug
    
    %%unpack vData
    shadow_fr_mar_dim = vData.shadow_fr_mar_dim;
    
    
    %%%%% make sphere thingy fr eyeball guys
    sphRes = 20;
    r =12; %your eyeballs have a radius of 12 mm
    [th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
    [x1,y1,z1] = sph2cart(th, phi, r);
    
    normScale = 1;
    plotSkel = true;
    
    lLeg = [2 3 4 5 6 7 5];
    rLeg = [2 8 9 10 11 12 10];
    tors = [2 13 14 15 26 27 28];
    lArm = [15 16 17 26 17 18 19 20];
    rArm = [15 21 22 26 22 23 24 25];
    
    
    figure(987)
    
    if camAlignErr < 0;
        frames = 1:10:length(gazeXYZ);
    else
        frames = 1;
    end
    
    for ii = frames%:length(eye_sphCenCam_x);
        subplot(1,2,1)
        cla
        
        %%eyeball centers in shadow coordinats(not to be confused with "rEye_sphCen_x", which are in pupil camera coords)
        ex = worldCamCenterXYZ(ii,1);
        ey = worldCamCenterXYZ(ii,2);
        ez = worldCamCenterXYZ(ii,3);
        
        
      
        
        
        
        % right eye
        e1 =  mesh(x1+ex, y1+ey, z1+ez);
        e1.EdgeColor = 'k';
        e1.EdgeAlpha = 0.1;
        hold on
        

        
        %%%%
        switch vData.dataType
            case 1
                gzCol = 'r';
            case 2
                gzCol = 'b';
            case 3
                gzCol = 'k';
        end
        
        plot3([0+ex nanmean(gazeXYZ(:,1))],...
            [0+ey nanmean(gazeXYZ(:,2))],...
            [0+ez nanmean(gazeXYZ(:,3))],'-s','LineWidth',2,'Color',gzCol)

        
        plot3(gazeXYZ(:,1),...
            gazeXYZ(:,2),...
            gazeXYZ(:,3),'-','Color',gzCol)

        
        
%         plot3(patchX(ii,:)+ex, patchY(ii,:)+ey, patchZ(ii,:)+ez,'k','LineWidth',3)

%         plot3(gazeGroundIntersectionXYZ(:,1),...
%             gazeGroundIntersectionXYZ(:,2),...
%             gazeGroundIntersectionXYZ(:,3),'m.')
%         
%         plot3(gazeGroundIntersectionXYZ(ii,1),...
%             gazeGroundIntersectionXYZ(ii,2),...
%             gazeGroundIntersectionXYZ(ii,3),'kp','MarkerFaceColor','r','MarkerSize',12)
%         
        %%%plot yer calibration point
        plot3(calibPoint(1), calibPoint(2), calibPoint(3), 'hr', 'MarkerFaceColor','r','MarkerSize',12)
        
        %%%Plotcherself a skeleetoon
        plot3(shadow_fr_mar_dim(ii,1:28,1),shadow_fr_mar_dim(ii,1:28,2),shadow_fr_mar_dim(ii,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
        hold on
        
        
        plot3(shadow_fr_mar_dim(ii,lLeg,1),shadow_fr_mar_dim(ii,lLeg,2),shadow_fr_mar_dim(ii,lLeg,3),'b','LineWidth',2)
        plot3(shadow_fr_mar_dim(ii,rLeg,1),shadow_fr_mar_dim(ii,rLeg,2),shadow_fr_mar_dim(ii,rLeg,3),'r','LineWidth',2)
        plot3(shadow_fr_mar_dim(ii,tors,1),shadow_fr_mar_dim(ii,tors,2),shadow_fr_mar_dim(ii,tors,3),'g','LineWidth',2)
        plot3(shadow_fr_mar_dim(ii,lArm,1),shadow_fr_mar_dim(ii,lArm,2),shadow_fr_mar_dim(ii,lArm,3),'b','LineWidth',2)
        plot3(shadow_fr_mar_dim(ii,rArm,1),shadow_fr_mar_dim(ii,rArm,2),shadow_fr_mar_dim(ii,rArm,3),'r','LineWidth',2)
        
        
        
        bx =   shadow_fr_mar_dim(ii,1,1);
        by =   shadow_fr_mar_dim(ii,1,2);
        bz =   shadow_fr_mar_dim(ii,1,3);
        
        grndSize = 3e3;
        g_x = meshgrid(rHeelXYZ(ii,1)-grndSize:100:rHeelXYZ(ii,1)+grndSize);
        g_z = meshgrid(rHeelXYZ(ii,2)-grndSize:100:rHeelXYZ(ii,2)+grndSize)';
        g_y = ones(size(g_x)) * min([rHeelXYZ(ii,2) lHeelXYZ(ii,2) ]);
        
        surface(g_x, g_y, g_z,'FaceColor','none','EdgeColor','k')
        
        
        
        axis equal
        title(num2str(ii))
        %     set(gca,'CameraUpVector',[0 1 0])
        xlabel('x');ylabel('y'); zlabel('z');
        axis([-grndSize+bx grndSize+bx -grndSize+by grndSize+by -grndSize+bz grndSize+bz])
        
        a = gca;
        a.CameraTarget = [ex ey ez]; %point figure 'camera' at COM
        a.CameraPosition = a.CameraTarget + [-1800 1800 2000]; %set camera position
        a.CameraViewAngle = 80;
        a.CameraUpVector = [ 0 1 0];
%         a.Position = [0 0 1 1];
        
        hold off
        
        subplot(122)
        plot(gazeXYZ);
        rx = refline(0,calibPoint(1)); rx.Color = 'r';
        ry = refline(0,calibPoint(2)); ry.Color = 'r';
        rz = refline(0,calibPoint(3)); rz.Color = 'r';
        drawnow
    end
end