%
% Convert a quaternion to a 3-by-3 rotation matrix. This
% does not require a unit-length quaternion, but it does
% require a non-zero quaternion.
%
% Returns a 3-by-3 matrix in row major order.
%
% Quaternion is defined by a 4-vector, q = [w x*i y*j z*k].
%   rotate = quaternion_to_R3_rotation(q);
%
% Ported from the Boost.Quaternion library at:
%   http://www.boost.org/libs/math/quaternion/HSO3.hpp
%
% @file    tools/matlab/quaternion_to_R3_rotation.m
% @author  Luke Tokheim, luke@motionnode.com
% @version 2.0
%

function [result] = quaternion_to_R3_rotation(q)
  if ~isvector(q) | (4 ~= length(q)),
    error('Input quaternion is not a vector of length 4.');
  end

  a = q(1);
  b = q(2);
  c = q(3);
  d = q(4);

  aa = a*a;
  ab = a*b;
  ac = a*c;
  ad = a*d;
  bb = b*b;
  bc = b*c;
  bd = b*d;
  cc = c*c;
  cd = c*d;
  dd = d*d;

  norme_carre = aa+bb+cc+dd;

  result = eye(3,3);
  if norme_carre > 1e1*eps,
    result(1,1) = (aa+bb-cc-dd)/norme_carre;
    result(1,2) = 2*(-ad+bc)/norme_carre;
    result(1,3) = 2*(ac+bd)/norme_carre;
    result(2,1) = 2*(ad+bc)/norme_carre;
    result(2,2) = (aa-bb+cc-dd)/norme_carre;
    result(2,3) = 2*(-ab+cd)/norme_carre;
    result(3,1) = 2*(-ac+bd)/norme_carre;
    result(3,2) = 2*(ab+cd)/norme_carre;
    result(3,3) = (aa-bb-cc+dd)/norme_carre;
  end

end
% function [result] = quaternion_to_R3_rotation(q)
