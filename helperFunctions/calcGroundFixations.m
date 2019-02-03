function [ gazeGroundIntersection ] = calcGroundFixations(  rHeelXYZ, lHeelXYZ,  gazeXYZ, eyeballCenterXYZ )


%CALCGROUNDFIXATIONS calculate intersections between gaze vector (GazeXYZ,
%which orignates from CamXYZ and the ground plane (defined as flat plane at
%height of lowest foot marker)



%check if camXYZ is a single [3,1] vector or a list of multiple vectors
if min(size(eyeballCenterXYZ)) == 1
    frames = 1;
else
    frames = 1:length(eyeballCenterXYZ);
end

 gazeGroundIntersection = nan(length(frames),3);
 check = nan(length(frames),1);
 
 
for ii = frames
    if min(size(eyeballCenterXYZ)) == 1 %if only only one frame is input (this if/else statement probably isn't necessary, but whatevs)
        
        grHeight = min([rHeelXYZ(2) lHeelXYZ(2) ]); %the Y values denote vertical position
        p0 = eyeballCenterXYZ; %endpoint #1 (origin of Gaze Vector = eyeballcenter)
        p1(1) = gazeXYZ(1); %endpoint #2 Endpoint of Gaze Vector
        p1(2) = gazeXYZ(2);
        p1(3) = gazeXYZ(3); %(doing it this way because the vectors come in different orientations (3x1 or 1x3)sometimes, and this is easier than, ya know, being smart)
    else
        grHeight = min([ rHeelXYZ(ii,2) lHeelXYZ(ii,2)]); %the Y values denote vertical position
        p0 = eyeballCenterXYZ(ii,:); %endpoint #1 (origin of Gaze Vector = eyeballcenter)
        p1 = gazeXYZ(ii,:); %endpoint #2  Endpoint of Gaze Vector
        
    end
    
    V0 = [0 grHeight 0 ]; % V) - any point on the ground plane
    n = [0 grHeight+10 0 ]; %normal Vector of the ground plane (I think?)
    
    
    %plane_line_intersect computes the intersection of a plane and a segment(or
    %a straight line)
    % Inputs:
    %       n: normal vector of the Plane
    %       V0: any point that belongs to the Plane
    %       p0: end point 1 of the segment p0p1
    %       p1:  end point 2 of the segment p0p1
    %Outputs:
    %      I    is the point of interection
    %     Check is an indicator:
    %      0 => disjoint (no intersection)
    %      1 => the plane intersects P0P1 in the unique point I
    %      2 => the segment lies in the plane
    %      3=>the intersection lies outside the segment P0P1
    
    [gazeGroundIntersection(ii,:),check(ii)] = plane_line_intersect(n,V0,p0,p1);
    
    
    if check(ii) ~= 1
        gazeGroundIntersection(ii,:) = [nan nan nan]; %this guy doesn't interesect with ol' Mr. Ground :(
    end
    
    debugPlot = false;
    if debugPlot
        % % %
        % % % debug plot
        % % %
        figure(3)
        cla
        plot3(rHeelXYZ(ii,1),rHeelXYZ(ii,2),rHeelXYZ(ii,3),'ro','MarkerFaceColor','r')
        hold on
        
        plot3(lHeelXYZ(ii,1),lHeelXYZ(ii,2),lHeelXYZ(ii,3),'bo','MarkerFaceColor','b')
        
        plot3(p0(1), p0(2), p0(3),'mo')
        plot3(p1(1), p1(2), p1(3),'mp')
        
        plot3([p0(1) p1(1)], [p0(2) p1(2)], [p0(3) p1(3)],'m-')
        
        plot3(gazeGroundIntersection(ii,1),gazeGroundIntersection(ii,2),gazeGroundIntersection(ii,3),'p','MarkerSize',8,'MarkerFaceColor','r')
        
        plot3(eyeballCenterXYZ(ii,1),eyeballCenterXYZ(ii,2), eyeballCenterXYZ(ii,3),'kp')
        
        g_x = meshgrid(rHeelXYZ(ii,1)-1000:100:rHeelXYZ(ii,1)+1000);
        
        g_y = ones(size(g_x)) * grHeight;
        
        g_z = meshgrid(rHeelXYZ(ii,2)-1000:100:rHeelXYZ(ii,2)+1000)';
        
        surface(g_x, g_y, g_z,'FaceColor','w','EdgeColor','k')
        
%         axis([rHeelXYZ(ii,1)-2000 rHeelXYZ(ii,1)+2000 rHeelXYZ(ii,2)-2000 rHeelXYZ(ii,2)+2000 rHeelXYZ(ii,3)-2000 rHeelXYZ(ii,3)+2000])
        drawnow
        pause(.1)
        hold off
    end
    
end

