function lsplot(listname)

% LSPLOT plot loudspeaker transfer functions
%
%function lsplot(listname)

file_list=readlist(listname);
[n_files, str_len]=size(file_list);
i=1;
filename=deblank(file_list(i,:));
[filetype,info_blocks,...
           n_channels,n_lines, sampling_rate,...
           first_line,last_line,n_directions,...
           comment1, comment2] = mtlrh(filename);

H=zeros(n_files,n_lines);
for i=1:n_files
filename=deblank(file_list(i,:));

H(i,:)=mtlrch(filename,1);
end;

magplot(H',sampling_rate);grcntrl;
end;
