function repair(list_name);
file_list=readlist(list_name);

[n_files,str_len] = size(file_list)
for k=1:n_files
filename=deblank(file_list(k,:));
filename
[filetype,info_blocks,...
           n_channels,n_lines, sampling_rate,...
           first_line,last_line,n_directions,...
           comment1, comment2] = mtlrh(filename);
n_lines=2*n_lines

mtlwh(filename,filetype,info_blocks,...
           n_channels,n_lines, sampling_rate,...
           first_line,last_line,n_directions,...
           comment1, comment2);

end
