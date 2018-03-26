%
% Plot a "raw" data file. Packed 16-bit signed integers,
% little-endian, 3, 4, 9, or 10 per sample.
%
% Returns the MxN matrix of sample data. Each row is a single 
% sample of "raw" data, un-processed sensor data.
% [accelerometer [magnetometer gyroscope] [temperature]] =
%   [ax ay az [mx my mz gx gy gz] [t]]
%
% Note that temperature output is not available on all
% versions of the data.
%
% @file    tools/matlab/plot_raw.m
% @author  Luke Tokheim, luke@motionnode.com
% @version 2.0
%

function [data] = plot_raw(filename)
  data = plot_stream(filename, 'int16');
end
% function [data] = plot_raw(filename)
