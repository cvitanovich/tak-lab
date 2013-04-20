function mtlcorhe(filename)


%function mtlcorhe(filename)


[filetype,info_blocks,...
 n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,...
 comment1, comment2] = mtlrh(filename)

%added Dec 4, 1995
%get direction_matrix
if n_directions > 1
direction_matrix=mtlrdir(filename,n_directions);
end;
%end addition

sampling_rate=24000
last_line=sampling_rate/(2*n_lines)*n_lines-1

mtlwh(filename,filetype,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);


%added Dec 4, 1995
%rewrite direction_matrix
if n_directions > 1
status=mtlwdir(filename,direction_matrix);
end;
%end addition

return;
