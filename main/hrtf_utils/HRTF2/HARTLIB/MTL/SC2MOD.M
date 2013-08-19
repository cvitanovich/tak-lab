function sc2mod(list)
file_list=readlist(list);
[n_files,str_len] = size(file_list);
for n=1:n_files
filename=deblank(file_list(n,:));
[file_type_ref,info_blocks,n_channels,n_lines_ref, sampling_rate_ref,...
 first_line_ref,last_line_ref,n_directions,comment1, comment2] =mtlrh(filename);

n_lines_new=2*n_lines_ref;
mtlwh(filename,file_type_ref,info_blocks,n_channels,n_lines_new, sampling_rate_ref,...
       first_line_ref,last_line_ref,n_directions,...
       comment1, comment2);
end;
disp('ready');
