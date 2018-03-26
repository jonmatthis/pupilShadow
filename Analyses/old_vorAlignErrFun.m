function [camAlignErr] = vorAlignErrFun(gazeXYZ, headGlobalRotMat, rHeelXYZ , lHeelXYZ , eyeballCenterXYZ, shadow_fr_mar_dim,  shadowMarkerNames, calibPoint, plotDebug, writerObj, camEulerGuess)

camQuatGuess = quaternion.eulerangles('123',camEulerGuess(1),camEulerGuess(2),camEulerGuess(3));
camRotMatGuess = camQuatGuess.RotationMatrix;

%PIn Gze to origin before rotating
gazeXYZ_zeroed(:,1) = gazeXYZ(:,1) - eyeballCenterXYZ(:,1);
gazeXYZ_zeroed(:,2) = gazeXYZ(:,2) - eyeballCenterXYZ(:,2);
gazeXYZ_zeroed(:,3) = gazeXYZ(:,3) - eyeballCenterXYZ(:,3);






for ii = 1:length(gazeXYZ)-10

    thisET_frame_unrot = camRotMatGuess * [gazeXYZ_zeroed(ii,1); gazeXYZ_zeroed(ii,2); gazeXYZ_zeroed(ii,3)];
    thisETframe = headGlobalRotMat(:,:,ii+10) * thisET_frame_unrot;
    
    headOrVec(ii,:) =     headGlobalRotMat(:,:,ii) * [2e3; 0;0];

    gazeXYZ_zeroed(ii,:) = thisETframe;
   
end

%Move gaze back to eyeball centered reference frame
gazeXYZ(:,1) = gazeXYZ_zeroed(:,1) + eyeballCenterXYZ(:,1);
gazeXYZ(:,2) = gazeXYZ_zeroed(:,2) + eyeballCenterXYZ(:,2);
gazeXYZ(:,3) = gazeXYZ_zeroed(:,3) + eyeballCenterXYZ(:,3);

[ groundFixation] = calcGroundFixations(  rHeelXYZ , lHeelXYZ ,  gazeXYZ , eyeballCenterXYZ  );



% for cc = 1:length(gazeXYZ)
%     err(cc) = sqrt(pdist([gazeXYZ(cc,:)+camXYZ(cc,:); calibPoints(1,:)])^2);
% end
for cc = 1:length(gazeXYZ)
    err(cc) = sqrt( (calibPoint(1) - gazeXYZ(cc,1))^2 + (calibPoint(2) - gazeXYZ(cc,2))^2 +(calibPoint(3) - gazeXYZ(cc,3))^2 );
end


camAlignErr = nansum(err)/length(err);

if plotDebug
    
    figure(2223)
    clf
    

    
    lLeg = [2 3 4 5 6 7 5];
    rLeg = [2 8 9 10 11 12 10];
    tors = [2 13 14 15 26 27 28];
    lArm = [15 16 17 26 17 18 19 20];
    rArm = [15 21 22 26 22 23 24 25];
   
    
    
    ex = eyeballCenterXYZ(end,1);
    ey = eyeballCenterXYZ(end,2);
    ez = eyeballCenterXYZ(end,3);
    
    subplot(1,2,1)

    
    plot3(gazeXYZ(end,1), gazeXYZ(end,2), gazeXYZ(end,3),'bp');     hold on
    
    
    
    %Plot skeleton
    plot3(shadow_fr_mar_dim(end,1:28,1),shadow_fr_mar_dim(end,1:28,2),shadow_fr_mar_dim(end,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
    
    
    plot3(shadow_fr_mar_dim(end, lLeg,1),shadow_fr_mar_dim(end, lLeg,2),shadow_fr_mar_dim(end, lLeg,3),'b','LineWidth',2)
    plot3(shadow_fr_mar_dim(end, rLeg,1),shadow_fr_mar_dim(end, rLeg,2),shadow_fr_mar_dim(end, rLeg,3),'r','LineWidth',2)
    plot3(shadow_fr_mar_dim(end, tors,1),shadow_fr_mar_dim(end, tors,2),shadow_fr_mar_dim(end, tors,3),'g','LineWidth',2)
    plot3(shadow_fr_mar_dim(end, lArm,1),shadow_fr_mar_dim(end, lArm,2),shadow_fr_mar_dim(end, lArm,3),'b','LineWidth',2)
    plot3(shadow_fr_mar_dim(end, rArm,1),shadow_fr_mar_dim(end, rArm,2),shadow_fr_mar_dim(end, rArm,3),'r','LineWidth',2)
    
    %plot presumed eyeball location
    plot3(eyeballCenterXYZ(end, 1),eyeballCenterXYZ(end, 2),eyeballCenterXYZ(end, 3),'kp')
    
    
    
    
    plot3( [ex gazeXYZ(end,  1)],  [ey gazeXYZ(end,  2)], [ez gazeXYZ(end,  3)],'-mo','MarkerSize',5,'MarkerFaceColor','m','LineWidth',4)
    
   
    
    plot3(groundFixation(:,1), groundFixation(:, 2), groundFixation(:, 3),'-')
    

    
    g_x = meshgrid(rHeelXYZ(ii,1)-4000:500:rHeelXYZ(ii,1)+4000);
    g_y = meshgrid(rHeelXYZ(ii,2)-4000:500:rHeelXYZ(ii,2)+4000)';
    g_z = ones(size(g_x)) * min([rHeelXYZ(end,3) lHeelXYZ(end,3) ]);
    
    surface(g_x, g_y, g_z,'FaceColor','w','EdgeColor','k')
    

    
    
    plot3(0, 0, min([rHeelXYZ(end,3) lHeelXYZ(end,3) ]),'rx','MarkerSize',10)
    
    %plot calibration points
    plot3(calibPoint(1),calibPoint(2),calibPoint(3),'ko','MarkerSize',8, 'MarkerFaceColor','k');

    
    axis([-3000+rHeelXYZ(ii,1) 3000+rHeelXYZ(ii,1) -1500+rHeelXYZ(ii,2) 2000+rHeelXYZ(ii,2) -3000+rHeelXYZ(ii,3) 3000+rHeelXYZ(ii,3) ])
    
    daspect([1,1,1])
    xlabel('x');ylabel('z'); zlabel('y'); %label Z axis as Y and vice versa, because that's how it's displayed
    
    title(strcat('Camera Rotation (Euler Angles): ', num2str(camEulerGuess,3)))
    
    hold off
 

    
    
    subplot(122)
    plot(calibPoint(1), calibPoint(2),'ko','MarkerSize',8, 'MarkerFaceColor','k');
    hold on

    plot(groundFixation(:,1), groundFixation(:, 2),'-')
        plot(calibPoint(1), calibPoint(2),'kp','MarkerSize',8, 'MarkerFaceColor','m');

        title(strcat('Gaze on Ground, Mean Error: ',num2str(camAlignErr,3),'mm'))
        
    axis equal
    grid on 
    grid minor
    
    axis([calibPoint(1)-5000 calibPoint(1)+5000 calibPoint(2)-5000 calibPoint(2)+5000])
    
    hold off
    f = gcf;
%     f.Position = [1 41 1920 1084];
    drawnow
    pause(.1)
%     if ~isempty(writerObj)
%         frame = getframe(gcf);
%         writeVideo(writerObj,frame);
%     end
end















