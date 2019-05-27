function [error_out] = errFunPoints(p1,p2,xguess)


% load('gpp.mat');

R = eul2rotm(xguess');





p_test = (R*p1')';

% p_test = p_test + x(1:3)';

[x,y,z] = sphere;

figure(1)
clf
plot3(p_test(:,1),p_test(:,2),p_test(:,3),'go');
hold on
plot3(p2(:,1),p2(:,2),p2(:,3),'ro');
mesh(x,y,z);
axis equal
lim = 1.2;
axis([-lim lim -lim lim -lim lim])
view(37,14)
error_out = sum(vecnorm([p2 - p_test]'));





end