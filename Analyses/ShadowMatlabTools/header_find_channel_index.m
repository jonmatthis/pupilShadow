%
% @file    tools/matlab/header_find_channel_index.m
% @version 2.4
%
function [index] = header_find_channel_index(header, node_index, node_channel)
  channel_stride = [
    4, 4, 4, 3, 3, 3, 3, 4, 3, 3, 3, 1, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1
  ];

  index = 0;
  itr = 1;
  for i=1:header.num_node,
    mask = header.node_header(i, 2);
    if 0 == mask,
      continue;
    end
    
    for j=1:length(channel_stride),
      channel = bitshift(1, j - 1);
      if bitand(mask, channel),
        if (i == node_index + 1) && (node_channel == channel),
          index = itr;
          break;
        end
        itr = itr + channel_stride(j);
      end
    end
    
    if index > 0,
      break;
    end
  end
end
