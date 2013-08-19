function mtliaerb(file1, outfile, frq)

% function mtliaerb(infile, outfile, frq)

% read head file1
[file_type1,info_blocks1,n_channels1,n_lines1, sampling_rate1,...
first_line1,last_line1,n_directions1,...
           comment1_1, comment2_1] = mtlrh(file1);

% compare header

index=0;
%if direction_matrix
if ((info_blocks1-2) > 0)
   direction_matrix=mtlrdir(file1,n_directions1);
   if (index ~= 0) direction_matrix(:) = direction_matrix(:,index); end;
   % if time_delay exists read time_delay
end;


% calculate output_header

n_byte=n_directions1*4;
n_blocks=ceil(n_byte/256);

filetype=11;
info_blocks=info_blocks1+n_blocks;
n_channels=n_directions1
n_lines=length(frq);
sampling_rate=sampling_rate1
first_line=min(frq);
last_line=last_line1
n_directions=n_directions1;
comment1=comment1_1;
comment2=comment2_1;
filename=outfile;


 
% write output_header

mtlwh(filename,filetype,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);


status=mtlwdir(filename,direction_matrix);
 
mtlwf(filename, frq);

for n=1:n_directions
 clc
  disp(['MTLIA: Calculating direction ', num2str(n), ' of ', num2str(n_channels)] );

 channel1=mtlrch(file1,(n-1)*2+1); 
 channel2=mtlrch(file1,n*2);

 channel_out = channel1./channel2;

 itd=ia_tdiff(channel1, channel2, frq, sampling_rate);
 ild=ia_ldiff(channel1, channel2, frq, sampling_rate);

 mtlwch(outfile,ild, itd ,n);
end;
return;

