%
% Read a binary data file of packed samples into an MxN
% matrix. N is defined by the nchannel argument. M is
% defined by N, the length of the data file, and the type
% of each data sample.
%
% Assumes that all data samples are little-endian.
%
% For example, a "raw" data file is 16-bit signed integers,
% little-endian, 9 or 10 values per sample.
%   data = import_binary(filename, 'int16', [9,10]); 
%
% @file    tools/matlab/import_binary.m
% @version 2.4
%

function [data, header] = import_binary(filename, type, nchannel)
  if nargin < 2,
    type = 'float32';
  end
  if nargin < 3,
    nchannel = 9;
  end

  fid = fopen(filename, 'rb', 'ieee-le');
  
  header = fread_take_stream_header(fid);
  if header.version > 0,
    nchannel = header.frame_stride / 4;
    type = 'float32';
  end
  
  buffer = fread(fid, inf, type);
  fclose(fid);

  n = nchannel(1);
  for i=(1:length(nchannel)),
    if mod(size(buffer,1), nchannel(i)) == 0,
      n = nchannel(i);
    end
  end

  if (n > 0) && mod(size(buffer,1), n) ~= 0,
    error('invalid number of samples in input buffer');
  end

  nsamples = size(buffer,1) / n;
  data = reshape(buffer, n, nsamples)';
end
% function [data] = import_binary(filename, type, nchannel)
