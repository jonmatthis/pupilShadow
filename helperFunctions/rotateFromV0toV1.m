function [x_r, y_r, theta] = rotateFromV0toV1(X, Y, pt0, pt1, origin, debug)

%%%%%[x_r, y_r] = rotateFromV0toV1(X, Y, pt0, pt1, origin, rotDir)
% % % rotates X and Y data by amount defined by angle between 
% % % vec0 (defined as [origin pt0]) and vec1(defined as [origin pt1])



if length(X) ~= length(Y)
    disp('X and Y must be the same length')
    return
end
% 

%% pt 1 & pt 2 define reference line of points to be rotated... I think ? 


vec_orig = pt0 - origin;
vec_goal = pt1 ;

[theta_orig, ~] = cart2pol(vec_orig(1), vec_orig(2));
[theta_goal, ~] = cart2pol(vec_goal(1), vec_goal(2));


theta = theta_goal - theta_orig;


theta = -theta; %Don't ask me, I just work here. Check the bottom playground cell if you want to check on this


%% rotate Data

for m = 1:length(X)
    
        x_r(m) = ...
            X(m) * cos(theta)+... %x*cos(theta)
            Y(m) * sin(theta);    %y*sin(theta)
        
        
        y_r(m) = ...
            -X(m) * sin(theta)+... %x*cos(theta)
             Y(m) * cos(theta);    %y*sin(theta)
    
end

if debug
    figure(8989)
    plot(X,Y,'-ob')
    hold on
%     plot(X(1),Y(1),'gp','MarkerSize',12,'MarkerFaceColor','g')
%     plot([origin(1) pt0(1)], [origin(2) pt0(2)],'g-p')

    plot(x_r, y_r,'-or')

    plot(x_r(1),y_r(1),'mp','MarkerSize',12,'MarkerFaceColor','m')
%     plot([origin(1) pt1(1)], [origin(2) pt1(2)],'m-p')
%     plot([origin(1) pt1(1)]*10, [origin(2) pt1(2)]*10,'m:p')
    
    axis equal
    hold off

end