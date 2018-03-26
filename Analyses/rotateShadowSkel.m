function [ shadow_fr_mar_dim_raw ] = rotateShadowSkel( shadow_fr_mar_dim_raw, shadowMarkerNames, headVec_fr_xyz, calibFrame )
%ROTATESHADOWSKEL Suiiary of this function goes here
%   Detailed explanation goes here
beep
hTop =  squeeze(shadow_fr_mar_dim_raw(:,strcmp('HeadTop', shadowMarkerNames),1:3)); %marker at Head Top
hC1 =   squeeze(shadow_fr_mar_dim_raw(:,strcmp('Head', shadowMarkerNames),1:3)); %marker at C1 vertebra

rHeelXYZ =   squeeze(shadow_fr_mar_dim_raw(:,strcmp('RightHeel', shadowMarkerNames),1:3));
lHeelXYZ =   squeeze(shadow_fr_mar_dim_raw(:,strcmp('LeftHeel', shadowMarkerNames),1:3));



hCen = (hTop+hC1)/2;

for ii = calibFrame:10:calibFrame+1e4
    ii
    lLeg = [2 3 4 5 6 7 5];
    rLeg = [2 8 9 10 11 12 10];
    tors = [2 13 14 15 26 27 28];
    lArm = [15 16 17 26 17 18 19 20];
    rArm = [15 21 22 26 22 23 24 25];
    
    figure(334);clf
    
    plot3(shadow_fr_mar_dim_raw(ii,1:28,1),shadow_fr_mar_dim_raw(ii,1:28,2),shadow_fr_mar_dim_raw(ii,1:28,3),'ko','MarkerFaceColor','k','MarkerSize',4)
    hold on
    
    plot3(shadow_fr_mar_dim_raw(ii,lLeg,1),shadow_fr_mar_dim_raw(ii,lLeg,2),shadow_fr_mar_dim_raw(ii,lLeg,3),'b','LineWidth',2)
    plot3(shadow_fr_mar_dim_raw(ii,rLeg,1),shadow_fr_mar_dim_raw(ii,rLeg,2),shadow_fr_mar_dim_raw(ii,rLeg,3),'r','LineWidth',2)
    plot3(shadow_fr_mar_dim_raw(ii,tors,1),shadow_fr_mar_dim_raw(ii,tors,2),shadow_fr_mar_dim_raw(ii,tors,3),'Color',[.5 .5 .5],'LineWidth',2)
    plot3(shadow_fr_mar_dim_raw(ii,lArm,1),shadow_fr_mar_dim_raw(ii,lArm,2),shadow_fr_mar_dim_raw(ii,lArm,3),'b','LineWidth',2)
    plot3(shadow_fr_mar_dim_raw(ii,rArm,1),shadow_fr_mar_dim_raw(ii,rArm,2),shadow_fr_mar_dim_raw(ii,rArm,3),'r','LineWidth',2)
    
    %%plot head vector
    hold on
    plot3([hCen(ii,1) hCen(ii,1)+100], [hCen(ii,2) hCen(ii,2)],       [hCen(ii,3) hCen(ii,3)],     '-rp') %x-vector, aka magnetic north
    plot3([hCen(ii,1) hCen(ii,1)],     [hCen(ii,2) hCen(ii,2)+100],   [hCen(ii,3) hCen(ii,3)],     '-gp') %y-vector, aka West
    plot3([hCen(ii,1) hCen(ii,1)],     [hCen(ii,2) hCen(ii,2)],       [hCen(ii,3) hCen(ii,3)+100], '-bp') %z-vector, aka vertical
    
    %plot head orientation
    plot3([hCen(ii,1) hCen(ii,1)+headVec_fr_xyz(ii,1)],[hCen(ii,2) hCen(ii,2)+headVec_fr_xyz(ii,2)],[hCen(ii,3) hCen(ii,3)+headVec_fr_xyz(ii,3)],'-mo','LineWidth',2)
    plot3([hCen(ii,1)+headVec_fr_xyz(ii-100:ii+100,1)],[hCen(ii,2)+headVec_fr_xyz(ii-100:ii+100,2)],[hCen(ii,3)+headVec_fr_xyz(ii-100:ii+100,3)],'-mo','LineWidth',2)
    
    g_x = meshgrid(rHeelXYZ(ii,1)-500:100:rHeelXYZ(ii,1)+500);
    g_y = meshgrid(rHeelXYZ(ii,2)-500:100:rHeelXYZ(ii,2)+500)';
    g_z = ones(size(g_x)) * min([rHeelXYZ(3) lHeelXYZ(3) ]);
    
    surface(g_x, g_y, g_z,'FaceColor','w','EdgeColor','k')
    
    
    
    axis equal
    drawnow
end
beep
