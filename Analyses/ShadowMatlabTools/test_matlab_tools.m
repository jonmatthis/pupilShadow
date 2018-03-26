%
% Test script for all Matlab tool scripts. Also, show
% basic usage of all of the scripts. The scripts have
% been simplified to maximize compatiblity with older
% versions of Matlab and GNU Octave.
%
% @file    tools/matlab/test_matlab_tools.m
% @author  Luke Tokheim, luke@motionnode.com
% @version 2.0
%

% function [] = test_matlab_tools()
  % Read a "raw" binary file. Plot the raw sensor data and
  % return it.
  path = pwd;
  filename_raw = strcat(path,'\test_data\raw.bin');
  raw = plot_raw(filename_raw);

  % Read a "sensor" binary file. Plot the calibrated sensor
  % data and return it.
  filename_sensor = strcat(path,'\test_data\sensor.bin');
  sensor = plot_sensor(filename_sensor);
  
  % Read a "output" binary file. Plot the output orientation
  % data as a Euler angles and return it.
  filename_output = strcat(path,'\test_data\output.bin');
  [euler,output] = plot_output_euler(filename_output);

  
  % Convert a single quaternion to a set of Euler angles.
  euler_1 = quaternion_to_euler(output(1,:));
  % Convert a single quaternion to a 3-by-3 rotation matrix.
  matrix_1 = quaternion_to_R3_rotation(output(1,:));
% end
% function [] = test_matlab_tools()
