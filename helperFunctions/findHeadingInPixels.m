function [shortTermHeading_normPosX,shortTermHeading_normPosY, longTermHeading_normPosX,longTermHeading_normPosY] = findHeadingInPixels(hData)


%% blow apart hData into its component variables
varNames = fieldnames(hData);
for i=1:length(varNames)
    eval([varNames{i} '=hData.' varNames{i} ';']);
end
%% preallocate some bad boyz

longTermHeading_normPosX  = nan(length(shadow_fr_mar_dim),1); %this will be based on where the sub will be at the end of the walk.
longTermHeading_normPosY  = longTermHeading_normPosX;

shortTermHeading_normPosX  = nan(length(shadow_fr_mar_dim),1); %this will be based on the subjects location 3 steps ahead of the current frame
shortTermHeading_normPosY = shortTermHeading_normPosX;

longTermHeadingXYZ  = nan(length(shadow_fr_mar_dim), 3); %this will be based on where the sub will be at the end of the walk.
shortTermHeadingXYZ  = nan(length(shadow_fr_mar_dim),3); %this will be based on the subjects location 3 steps ahead of the current frame

%% find head center
hTop =  squeeze(shadow_fr_mar_dim(:,strcmp('HeadTop', shadowMarkerNames),1:3)); %marker at Head Top
hC1 =   squeeze(shadow_fr_mar_dim(:,strcmp('Head', shadowMarkerNames),1:3)); %marker at C1 vertebra
hCen = (hTop+hC1)/2;

%% find short/long term heading in each frame in XYZ world coordinates
for fr = 1:length(hCen)
    
    if mod(fr,1000)==0
        disp(['Finding heading in Pixel Coords, frame: ', num2str(fr), ' of ', num2str(length(hCen))])
    end
    
    %%%%find head positoin 3 seconds into the future
    shortTermHeadingXYZ(fr,:) = hCen(min([length(hCen) fr+(framerate*3)]),:); %short term heading is where the head will be 3 seconds from now
    
    %%% find where the head will be at the end of the current walk
    walkEnds = walks(:,2);
    walkEnds(walkEnds<fr) = []; %snip off walkEnds that are less than this frame, so walkEnds(1) is the frame we need
    
    if ~isempty(walkEnds) %if there is a walkEnd left, walkEnds(1) is our guy
        longTermHeadingXYZ(fr,:) = hCen(walkEnds(1),:);
    else %if there are no walkEnds left, use the end of the headTrajectory instead
        longTermHeadingXYZ(fr,:) = hCen(end,:);
    end
    
    
    %% find pixel location of short and long term heading on each frame in loop
    
    
    hx = hCen(fr,1);
    hy = hCen(fr,2);
    hz = hCen(fr,3);
    
    thisTopLeft = patchTopLeft(fr,:) + hCen(fr,:);
    thisBottomLeft = patchBottomLeft(fr,:) + hCen(fr,:);
    thisTopRight = patchTopRight(fr,:) + hCen(fr,:);
    thisBottomRight = patchBottomRight(fr,:) + hCen(fr,:);
    
    %create the image plane (using geom3d toolbox, here and below)
    imagePlane = createPlane(thisTopLeft, thisTopRight, thisBottomLeft);
    
    %find intersection between shortTermHeading and image plane
    shortTermHeadingEdge = [hCen(fr,:) shortTermHeadingXYZ(fr,:)];
    shortTermHeadingIntersectXYZ = intersectLinePlane(shortTermHeadingEdge, imagePlane);
    
    %ditto for long term heading point
    longTermHeadingEdge = [hCen(fr,:) longTermHeadingXYZ(fr,:)];
    longTermHeadingIntersectXYZ = intersectLinePlane(longTermHeadingEdge, imagePlane);
    
    % project intersection point onto top and left edge of image plane to find screen X and Y position, respectively
    
    imagePlaneTopEdge_z = thisTopRight-thisTopLeft; %top edge of image plane in TopLeftCorner-centered reference frame
    imagePlaneLeftEdge_z = thisBottomLeft-thisTopLeft; %top edge of image plane in TopLeftCorner-centered reference frame
    
    
    
    
    shortTermHeadingIntersectXYZ_z = shortTermHeadingIntersectXYZ-thisTopLeft;
    longTermHeadingIntersectXYZ_z = longTermHeadingIntersectXYZ-thisTopLeft;
    
    %find normalized projection of short/long term heading onto Top/Left edge of screen to find normalized X and Y screen coordinates
    shortTermHeading_normPosX(fr) = dot(imagePlaneTopEdge_z/norm(imagePlaneTopEdge_z), shortTermHeadingIntersectXYZ_z/norm(imagePlaneTopEdge_z));
    shortTermHeading_normPosY(fr) = dot(imagePlaneLeftEdge_z/norm(imagePlaneLeftEdge_z), shortTermHeadingIntersectXYZ_z/norm(imagePlaneLeftEdge_z));
    
    longTermHeading_normPosX(fr) = dot(imagePlaneTopEdge_z/norm(imagePlaneTopEdge_z), longTermHeadingIntersectXYZ_z/norm(imagePlaneTopEdge_z));
    longTermHeading_normPosY(fr) = dot(imagePlaneLeftEdge_z/norm(imagePlaneLeftEdge_z), longTermHeadingIntersectXYZ_z/norm(imagePlaneLeftEdge_z));
    
    showDebug = false;%I can't think of a good summary debug plot, but the booled-out videos show a nice enough frame by frame
    if showDebug && mod(fr,10)==0
        figure(743);clf
        subplot(122)
        plot3(0,0,0,'rp')
        hold on
        plot3([0 imagePlaneTopEdge_z(1)], [0 imagePlaneTopEdge_z(2)], [0 imagePlaneTopEdge_z(3)],'c-o')%top edge of image plane
        plot3([0 imagePlaneTopEdge_z(1)*shortTermHeading_normPosX(fr)], [0 imagePlaneTopEdge_z(2)*shortTermHeading_normPosX(fr)], [0 imagePlaneTopEdge_z(3)*shortTermHeading_normPosX(fr)],'b-o','LineWidth',2)%top edge of image plane
        plot3([0 imagePlaneTopEdge_z(1)*longTermHeading_normPosX(fr)], [0 imagePlaneTopEdge_z(2)*longTermHeading_normPosX(fr)], [0 imagePlaneTopEdge_z(3)*longTermHeading_normPosX(fr)],'g--o','LineWidth',2)%top edge of image plane
        
        
        plot3([0 imagePlaneLeftEdge_z(1)], [0 imagePlaneLeftEdge_z(2)], [0 imagePlaneLeftEdge_z(3)],'m-o')
        plot3([0 imagePlaneLeftEdge_z(1)*shortTermHeading_normPosY(fr)], [0 imagePlaneLeftEdge_z(2)*shortTermHeading_normPosY], [0 imagePlaneLeftEdge_z(3)*shortTermHeading_normPosY(fr)],'r-o','LineWidth',2)%top edge of image plane
        plot3([0 imagePlaneLeftEdge_z(1)*longTermHeading_normPosY(fr)], [0 imagePlaneLeftEdge_z(2)*longTermHeading_normPosY(fr)], [0 imagePlaneLeftEdge_z(3)*longTermHeading_normPosY(fr)],'k--o','LineWidth',2)%top edge of image plane
        
        plot3(shortTermHeadingIntersectXYZ_z(1), shortTermHeadingIntersectXYZ_z(2), shortTermHeadingIntersectXYZ_z(3),'mh')
        plot3(longTermHeadingIntersectXYZ_z(1), longTermHeadingIntersectXYZ_z(2), longTermHeadingIntersectXYZ_z(3),'kp')
        
        axis equal
        
        
        subplot(121)
        
        lLeg = [2 3 4 5 6 7 5];
        rLeg = [2 8 9 10 11 12 10];
        tors = [2 13 14 15 26 27 28];
        lArm = [15 16 17 26 17 18 19 20];
        rArm = [15 21 22 26 22 23 24 25];
        
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
        
        
        
        
        line([0 patchTopLeft(fr,1)]+hx, [0 patchTopLeft(fr,2)]+hy, [0 patchTopLeft(fr,3)]+hz,'color','k','linewidth',2);
        line([0 patchTopRight(fr,1)]+hx, [0 patchTopRight(fr,2)]+hy, [0 patchTopRight(fr,3)]+hz,'color','k','linewidth',2);
        line([0 patchBottomLeft(fr,1)]+hx, [0 patchBottomLeft(fr,2)]+hy, [0 patchBottomLeft(fr,3)]+hz,'color','k','linewidth',2);
        line([0 patchBottomRight(fr,1)]+hx, [0 patchBottomRight(fr,2)]+hy, [0 patchBottomRight(fr,3)]+hz,'color','k','linewidth',2);
        
        thisFrustum = [patchTopLeft(fr,:); patchTopRight(fr,:); patchBottomRight(fr,:);patchBottomLeft(fr,:); patchTopLeft(fr,:)];
        plot3(thisFrustum(:,1)+hx, thisFrustum(:,2)+hy, thisFrustum(:,3)+hz,'r-o','LineWidth',3)
        
        plot3(thisTopLeft(1), thisTopLeft(2), thisTopLeft(3),'gp', 'MarkerSize',20)
        plot3(thisBottomLeft(1), thisBottomLeft(2), thisBottomLeft(3),'mp', 'MarkerSize',20)
        
        plot3(thisTopRight(1), thisTopRight(2), thisTopRight(3),'bp', 'MarkerSize',20)
        plot3(thisBottomRight(1), thisBottomRight(2), thisBottomRight(3),'kp', 'MarkerSize',20)
        
        
        plot3(shortTermHeadingIntersectXYZ(1), shortTermHeadingIntersectXYZ(2), shortTermHeadingIntersectXYZ(3),'mh')
        
        plot3(longTermHeadingIntersectXYZ(1), longTermHeadingIntersectXYZ(2), longTermHeadingIntersectXYZ(3),'kp')
        
        bx =   shadow_fr_mar_dim(fr,1,1);
        by =   shadow_fr_mar_dim(fr,1,2);
        bz =   shadow_fr_mar_dim(fr,1,3);
        
        [xx,zz] = meshgrid(-2000+bx:100:2000+bx, -2500+bz:100:2000+bz);
        surface(xx,zeros(size(xx)),zz,'FaceColor','w')
        
        axis([-4000+bx 4000+bx -4000+by 4000+by -4500+bz 4000+bz])
        
        view(-173, -43);
        axis equal
        title(num2str(fr))
        set(gca,'CameraUpVector',[0 1 0])
        xlabel('x');ylabel('y'); zlabel('z');
        
        drawnow
        
    end
end

