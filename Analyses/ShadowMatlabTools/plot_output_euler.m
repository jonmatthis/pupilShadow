%
% Plot an "output" data file. Packed 32-bit IEEE floating
% point numbers, little-endian, 4 per sample.
%
% Returns the Mx3 matrix of Euler angle data. Each row is a
% single set of Euler angles, specified in degrees in x-y-z
% rotation order. e = [x y z].
%
% Returns the Mx4 matrix of sample data. Each row is a single 
% sample of "output" data, global quaternion orientation data.
% q = [w x*i y*j z*k], |q| = 1
%
% @file    tools/matlab/plot_output_euler.m
% @author  Luke Tokheim, luke@motionnode.com
% @version 2.0
%

function [data,output] = plot_output_euler(filename)
  output = import_binary(filename, 'float32', 4);

  data = [];
  for i=[1:size(output,1)],
    % Normalize the quaternion.
    magnitude = sqrt(dot(output(i,:), output(i,:)));
    output(i,:) = output(i,:) ./ magnitude;
      
    data(i,:) = quaternion_to_euler(output(i,:))';
  end
  
  % Convert radians to degrees.
  data = data * (180/pi);
  
  clf;
  plot(data);
end
% function [data,output] = plot_output_euler(filename)
