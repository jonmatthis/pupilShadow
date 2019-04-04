function [error_out] = errFunPoints(x)


load('gpp.mat');

R = eul2rotm(x');





p_test = (R*p1')';

% p_test = p_test + x(1:3)';


figure(1)
clf
plot3(p_test(:,1),p_test(:,2),p_test(:,3),'go');
hold on
plot3(p2(:,1),p2(:,2),p2(:,3),'ro');
axis equal

error_out = sum(vecnorm([p2 - p_test]'));





end