function [Hl,Hr,azi,ele]=readsc2(filename);
% READSC2 reads SC2 file
%         [Hl,Hr,azi,ele]=READSC2(filename);


[filetype,info_blocks,n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,comment1, comment2] = readh(filename);


direction_matrix=read_dir(filename,n_directions);
ele=direction_matrix(1,:);
azi=direction_matrix(2,:);

[fid,msg]=fopen(filename,'r','ieee-le');
if fid<0 error(msg); end;

N=2*n_lines;

fseek(fid,info_blocks*256,'bof');
H=fread(fid,inf,'float32');
Hl=H(1:N); Hr=H(N+1:2*N);
Hl=Hl(1:2:N)+i*Hl(2:2:N);
Hr=Hr(1:2:N)+i*Hr(2:2:N);

fclose(fid);
