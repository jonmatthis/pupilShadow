function[calibDist,gazeXYZ] = calculateGazeVectors(theEye,theEyeAlignRotMat,theEyeballCenterXYZ,calibFrame,calibPoint,headRotMat_row_col_fr)

% right eye
gazeXYZ = [theEye.circle_3d_center_x theEye.circle_3d_center_y theEye.circle_3d_center_z] ...  Take your "PupilCircleCenter" (in 3D Eye camera coordinate system, units are mm)
    -[theEye.sphere_center_x  theEye.sphere_center_y  theEye.sphere_center_z];        %Subtract EyeSphereCenter (in eye camera coordiates) >> Origin is now at the center of the EyeSphere in camera coords

%normalize its length
gazeXYZ = gazeXYZ./(vecnorm(gazeXYZ')');

%multiply it by your desired length ;)
calibDist = pdist([theEyeballCenterXYZ(calibFrame,:); calibPoint]); %myboy pythag
gazeXYZ = gazeXYZ*calibDist*10;

%%%%
%%%%%%% Rotate gaze vector by the camera alignment matrix from the 
%%%%%%% "vorAlignLossFun' (aka theEyeAlignRotMat) & head orientation 
%%%%%%% (in that order), prior to resituating  the origin on on the eyeball

for rr = 1:length(gazeXYZ)    
    thisET_frame_unrot = theEyeAlignRotMat * gazeXYZ(rr,1:3)';
    thisETframe = headRotMat_row_col_fr(:,:,rr) * thisET_frame_unrot;
    gazeXYZ(rr,:) = thisETframe;
end

% add the eyeball center (in shadow/world coordinates) to translate origin of gaze vector onto the shadow eyeball
gazeXYZ = gazeXYZ+ theEyeballCenterXYZ;
