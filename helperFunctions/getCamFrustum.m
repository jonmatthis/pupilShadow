function [frust] = getCamFrustum(headVecX_fr_xyz,headVecY_fr_xyz,headVecZ_fr_xyz,gaze_norm_pos_x,gaze_norm_pos_y,...
    px2mmScale,calibDist,rGazeXYZ,lGazeXYZ,rEyeballCenterXYZ,lEyeballCenterXYZ,resHeight,resWidth);

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





