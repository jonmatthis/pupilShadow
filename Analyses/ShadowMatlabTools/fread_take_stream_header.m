%
% @file    tools/matlab/fread_take_stream_header.m
% @version 2.4
%

function [s, offset] = fread_take_stream_header(fid)
  s = struct();
  s.version = 0;
  
  offset = 0;

  % Detect the header by our take file magic bytes.
  magic = int32(fread(fid, 2, 'int32'));
  if ~isequal(magic, [-8882056; 87652969]),
    % No header present, rewind.
    frewind(fid);
    return;
  end

  s.magic = magic;  
  s.version = fread(fid, 1, 'int32');
  s.uuid = int32(fread(fid, 4, 'int32'));
  s.num_node = fread(fid, 1, 'uint32');
  s.frame_stride = fread(fid, 1, 'uint32');
  s.num_frame = fread(fid, 1, 'uint32');
  s.channel_mask = uint32(fread(fid, 1, 'uint32'));
  s.h = fread(fid, 1, 'float32');
  s.location = fread(fid, 3, 'float32');
  s.geomagnetic = fread(fid, 3, 'float32');
  s.tv_sec = fread(fid, 2, 'uint32');
  s.tv_usec = fread(fid, 1, 'uint32');

  % Padding. Reserved to 128 bytes.
  fread(fid, 11, 'int32');

  % Two integers per node. The node key and its channel mask.
  s.node_header = uint32(fread(fid, 2 * s.num_node, 'uint32'));
  s.node_header = reshape(s.node_header, 2, s.num_node)';

  % Optionally report back the offset in bytes to get past the header
  if nargout > 1,
    offset = ftell(fid); 
  end
end
% function [s, offset] = fread_take_stream_header(fid)
