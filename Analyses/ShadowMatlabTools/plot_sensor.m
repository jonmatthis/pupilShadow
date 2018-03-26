%
% Plot a "sensor" data file. Packed 32-bit IEEE floating
% point numbers, little-endian, 3, 4, 9, or 10 per sample.
%
% Returns the MxN matrix of sample data. Each row is a single
% sample of "sensor" data, calibrated sensor data.
% [accelerometer [magnetometer gyroscope] [temperature]] =
%   [ax ay az [mx my mz gx gy gz] [t]]
%
% Note that temperature output is not available on all
% versions of the data.
%
% @file    tools/matlab/plot_sensor.m
% @author  Luke Tokheim, luke@motionnode.com
% @version 2.0
%

function [data] = plot_sensor(filename)
  data = plot_stream(filename, 'float32');
end
% function [data] = plot_sensor(filename, has_temperature)
