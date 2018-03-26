%
% Convert a quaternion to a set of Euler angles. Requires
% an input quaternion of unit length. 
%
% Returns a 3-vector [x y z] in x-y-z rotation order.
%
% Quaternion is defined by a 4-vector, q = [w x*i y*j z*k].
%   euler = quaternion_to_euler(q);
%
% @file    tools/matlab/quaternion_to_euler.m
% @author  Luke Tokheim, luke@motionnode.com
% @version 2.0
%

function [result] = quaternion_to_euler(q)
  if ~isvector(q) || (4 ~= length(q)),
    error('Input quaternion is not a vector of length 4.');
  end

  if abs(1 - (sqrt(dot(q, q)))) > 1e-6,
    error('Input quaternion is not of unit length.');
  end

  q1 = q(1);
  q2 = q(2);
  q3 = q(3);
  q4 = q(4);

  q1_2 = q1*q1;
  q2_2 = q2*q2;
  q3_2 = q3*q3;
  q4_2 = q4*q4;

  result = zeros(3,1);
  result(1) = atan2(2 * (q3*q4 + q1*q2), q1_2 - q2_2 - q3_2 + q4_2);
  result(2) = asin(-2 * (q2*q4 - q1*q3));
  result(3) = atan2(2 * (q2*q3 + q1*q4), q1_2 + q2_2 - q3_2 - q4_2);
end
% function [result] = quaternion_to_euler(q)
