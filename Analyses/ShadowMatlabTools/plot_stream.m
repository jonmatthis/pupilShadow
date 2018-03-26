%
% Plot a Motion Service stream data file. This will either be a
% "raw" data file or a "sensor" data file. The data_type
% input parameter should either be 'float32' or 'int16'.
%
% Returns the MxN matrix of sample data. Each row is a single
% sample of stream data.
% [accelerometer [magnetometer gyroscope] [temperature]] =
%   [ax ay az [mx my mz gx gy gz] [t]]
%
% Detect accelerometer only streams at run time and remove
% zero data in the magnetometer and gyroscope channels.
%
% Note that temperature output is not available on all
% versions of the data.
%
% @file    tools/matlab/plot_stream.m
% @author  Luke Tokheim, luke@motionnode.com
% @version 2.4
%

function [data] = plot_stream(filename, data_type)
  if nargin < 2,
    data_type = 'float32';
  end
  nchannel = [9,10];

  [data,header] = import_binary(filename, data_type, nchannel);

  if header.version > 0,
    if nargin < 2 || ~strcmp('int16', data_type),
      buffer = zeros(size(data, 1), 9);
      for i=1:3,
        index = header_find_channel_index(header, 1, bitshift(1, 7 + i));
        
        buffer(:, (1:3) + (i - 1) * 3) = data(:, (0:2) + index);
      end
      
      data = buffer;
    else
      buffer = zeros(size(data, 1), 9);
      for i=1:3,
        index = header_find_channel_index(header, 1, bitshift(1, 11 + i));
        
        buffer(:, (1:3) + (i - 1) * 3) = data(:, (0:2) + index);
      end
      
      data = buffer;
    end
  else
    % Detect data files with temperature output if the
    % number of elements in ambiguous.
    n_element = prod(size(data));
    if (mod(n_element, 10) == 0),
      if (mod(n_element, 9) == 0),
        if (size(data,2) == 10) && (data(1,10) ~= 0),
          data = reshape(data', 9, n_element / 9)';
        elseif (size(data,2) == 9) && (data(2,1) == 0),      
          data = reshape(data', 10, n_element / 10)';
        end
      end

      % Duplicate the second temperature sample in the logged
      % data. The first sample is zero to allow for file
      % detection.
      if (size(data,2) == 10) && (data(1,10) == 0),
        data(1,10) = data(2,10);
      end
    end
  end
  
  clf;
  nplot = 3;
  
  % Detect accelerometer only data stream. Remove
  % zero data in the magnetometer and gyroscope
  % channels.
  if (size(data,1) > 0) && (sum(data(1,4:9)) == 0),
    data(:,4:9) = [];
    nplot = 1;
  end
  
  x = linspace(0, size(data,1) - 1, size(data,1));
  if header.version > 0,
    x = x * header.h;
  end
  
  % Accelerometer.
  subplot(nplot,1,1);
  plot(x, data(:,1:3));
  
  if size(data,2) >= 6,
    % Magnetometer.
    subplot(nplot,1,2);
    plot(x, data(:,4:6));
  end
  
  if size(data,2) >= 9,
    % Gyroscope.
    subplot(nplot,1,3);
    plot(x, data(:,7:9));
  end
end
% function [data] = plot_stream(filename, data_type)
