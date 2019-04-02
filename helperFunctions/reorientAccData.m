function  [hipsAccXYZ_rot, chestAccXYZ_rot, headAccXYZ_rot] =  reorientAccData(inData ,debug)

%%unpack inData struct
varNames = fieldnames(inData);
for i=1:length(varNames)
    eval([varNames{i} '=inData.' varNames{i} ';']);
end

%% define head center (the origin of the head vectors)
hTopXYZ =  squeeze(shadow_fr_mar_dim(:,strcmp('HeadTop', shadowMarkerNames),1:3)); %marker at Head Top
hC1XYZ =   squeeze(shadow_fr_mar_dim(:,strcmp('Head', shadowMarkerNames),1:3)); %marker at C1 vertebra
hCenXYZ = (hTopXYZ+hC1XYZ)/2;% head center

hipsXYZ =   squeeze(shadow_fr_mar_dim(:,strcmp('Hips', shadowMarkerNames),1:3)); %marker at C1 vertebra

chestXYZ =   squeeze(shadow_fr_mar_dim(:,strcmp('Chest', shadowMarkerNames),1:3)); %marker at C1 vertebra


% acc vectors in the calibration frame
headAccCalibXYZ = headAccXYZ(calibFrame,:);
chestAccCalibXYZ = chestAccXYZ(calibFrame,:);
hipsAccCalibXYZ = hipsAccXYZ(calibFrame,:);

%vertical projection of relevant marker tells you which way the gravity vector should point
headDownVectorXYZ   = [hCenXYZ(calibFrame,1)  0 hCenXYZ(calibFrame,3)]  - hCenXYZ(calibFrame,:);
chestDownVectorXYZ  = [chestXYZ(calibFrame,1) 0 chestXYZ(calibFrame,3)] - chestXYZ(calibFrame,:);
hipsDownVectorXYZ   = [hipsXYZ(calibFrame,1)  0 hipsXYZ(calibFrame,3)]  - hipsXYZ(calibFrame,:);

%rotate acceleration vector to align with gravity/down vector
%%%HEAD
[headAccCalibTh,headAccCalibPhi, headAccCalibRho] = cart2sph(headAccCalibXYZ(1),headAccCalibXYZ(2), headAccCalibXYZ(3));
[headDownTh,headDownPhi, headDownRho] = cart2sph(headDownVectorXYZ(1),headDownVectorXYZ(2), headDownVectorXYZ(3));

headDownRotTh = headDownTh-headAccCalibTh; %difference between Acc theta and down-vector on calib frame
headDownRotPhi = headDownPhi-headAccCalibPhi; %difference between Acc phi and down-vector on calib frame

[headAccTh,headAccPhi, headAccRho] = cart2sph(headAccXYZ(:,1),headAccCalibXYZ(:,2), headAccCalibXYZ(:,3));

[headAccXYZ_rot(:,1),headAccXYZ_rot(:,2),headAccXYZ_rot(:,3)]=sph2cart((headAccTh+headDownRotTh), headAccPhi+headDownRotPhi,headAccRho); %rotate all acc vectors to align them with gravity vertical, and rotate them so that forward is head-forward

%%%CHEST
[chestAccCalibTh,chestAccCalibPhi, chestAccCalibRho] = cart2sph(chestAccCalibXYZ(1),chestAccCalibXYZ(2), chestAccCalibXYZ(3));
[chestDownTh,chestDownPhi, chestDownRho] = cart2sph(chestDownVectorXYZ(1),chestDownVectorXYZ(2), chestDownVectorXYZ(3));

chestDownRotTh = chestDownTh-chestAccCalibTh; %difference between Acc theta and down-vector on calib frame
chestDownRotPhi = chestDownPhi-chestAccCalibPhi; %difference between Acc phi and down-vector on calib frame

[chestAccTh,chestAccPhi, chestAccRho] = cart2sph(chestAccXYZ(:,1),chestAccCalibXYZ(:,2), chestAccCalibXYZ(:,3));

[chestAccXYZ_rot(:,1),chestAccXYZ_rot(:,2),chestAccXYZ_rot(:,3)]=sph2cart((chestAccTh+chestDownRotTh), chestAccPhi+chestDownRotPhi,chestAccRho); %rotate all acc vectors to align them with gravity vertical

%%%HIPS
[hipsAccCalibTh,hipsAccCalibPhi, hipsAccCalibRho] = cart2sph(hipsAccCalibXYZ(1),hipsAccCalibXYZ(2), hipsAccCalibXYZ(3));
[hipsDownTh,hipsDownPhi, hipsDownRho] = cart2sph(hipsDownVectorXYZ(1),hipsDownVectorXYZ(2), hipsDownVectorXYZ(3));

hipsDownRotTh = hipsDownTh-hipsAccCalibTh; %difference between Acc theta and down-vector on calib frame
hipsDownRotPhi = hipsDownPhi-hipsAccCalibPhi; %difference between Acc phi and down-vector on calib frame

[hipsAccTh,hipsAccPhi, hipsAccRho] = cart2sph(hipsAccXYZ(:,1),hipsAccCalibXYZ(:,2), hipsAccCalibXYZ(:,3));

[hipsAccXYZ_rot(:,1),hipsAccXYZ_rot(:,2),hipsAccXYZ_rot(:,3)]=sph2cart((hipsAccTh+hipsDownRotTh), hipsAccPhi+hipsDownRotPhi,hipsAccRho); %rotate all acc vectors to align them with gravity vertical


%%
if debug
    %%
    figure(8544)
    clf
    lLeg = [2 3 4 5 6 7 5];
    rLeg = [2 8 9 10 11 12 10];
    tors = [2 13 14 15 26 27 28];
    lArm = [15 16 17 26 17 18 19 20];
    rArm = [15 21 22 26 22 23 24 25];
    
    for fr = calibFrame%-1000:10:calibFrame+30000
        
        
        
        plot3(shadow_fr_mar_dim(fr,1:28,1),shadow_fr_mar_dim(fr,1:28,2),shadow_fr_mar_dim(fr,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
        hold on
        
        
        plot3(shadow_fr_mar_dim(fr,lLeg,1),shadow_fr_mar_dim(fr,lLeg,2),shadow_fr_mar_dim(fr,lLeg,3),'b','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,rLeg,1),shadow_fr_mar_dim(fr,rLeg,2),shadow_fr_mar_dim(fr,rLeg,3),'r','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,tors,1),shadow_fr_mar_dim(fr,tors,2),shadow_fr_mar_dim(fr,tors,3),'g','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,lArm,1),shadow_fr_mar_dim(fr,lArm,2),shadow_fr_mar_dim(fr,lArm,3),'b','LineWidth',2)
        plot3(shadow_fr_mar_dim(fr,rArm,1),shadow_fr_mar_dim(fr,rArm,2),shadow_fr_mar_dim(fr,rArm,3),'r','LineWidth',2)
        
        
        %plot head axes
        hx = hCenXYZ(fr,1);
        hy = hCenXYZ(fr,2);
        hz = hCenXYZ(fr,3);
        
        headVecScale = 500;
        
        plot3([0 headVecX_fr_xyz(fr,1)*headVecScale]+hx, [0 headVecX_fr_xyz(fr,2)*headVecScale]+hy, [0 headVecX_fr_xyz(fr,3)*headVecScale]+hz,'r-o','MarkerSize',3,'LineWidth',2)
        plot3([0 headVecY_fr_xyz(fr,1)*headVecScale]+hx, [0 headVecY_fr_xyz(fr,2)*headVecScale]+hy, [0 headVecY_fr_xyz(fr,3)*headVecScale]+hz,'g-o','MarkerSize',3,'LineWidth',2)
        plot3([0 headVecZ_fr_xyz(fr,1)*headVecScale]+hx, [0 headVecZ_fr_xyz(fr,2)*headVecScale]+hy, [0 headVecZ_fr_xyz(fr,3)*headVecScale]+hz,'b-o','MarkerSize',3,'LineWidth',2)
        
        accScale =1000;
%         quiver3(hx, hy, hz, headAccXYZ(fr,1)*accScale,  headAccXYZ(fr,2)*accScale,  headAccXYZ(fr,3)*accScale,'AutoScale','off')
        quiver3(hx, hy, hz, headAccXYZ_rot(fr,1)*accScale,  headAccXYZ_rot(fr,2)*accScale,  headAccXYZ_rot(fr,3)*accScale,'LineWidth',3,'AutoScale','off')

%         quiver3(chestXYZ(fr,1), chestXYZ(fr,2), chestXYZ(fr,3), chestAccXYZ(fr,1)*accScale,  chestAccXYZ(fr,2)*accScale,  chestAccXYZ(fr,3)*accScale,'AutoScale','off')
        quiver3(chestXYZ(fr,1), chestXYZ(fr,2), chestXYZ(fr,3), chestAccXYZ_rot(fr,1)*accScale,  chestAccXYZ_rot(fr,2)*accScale,  chestAccXYZ_rot(fr,3)*accScale,'LineWidth',3,'AutoScale','off')

%         quiver3(hipsXYZ(fr,1), hipsXYZ(fr,2), hipsXYZ(fr,3), hipsAccXYZ(fr,1)*accScale,  hipsAccXYZ(fr,2)*accScale,  hipsAccXYZ(fr,3)*accScale,'AutoScale','off')
        quiver3(hipsXYZ(fr,1), hipsXYZ(fr,2), hipsXYZ(fr,3), hipsAccXYZ_rot(fr,1)*accScale,  hipsAccXYZ_rot(fr,2)*accScale,  hipsAccXYZ_rot(fr,3)*accScale,'LineWidth',3,'AutoScale','off')

        %plot calib point
        plot3(calibPoint(1), calibPoint(2), calibPoint(3), 'kp', 'MarkerFaceColor','y', 'MarkerSize',16)
        
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
        view(-182, -67)
        drawnow
    end%%
end




