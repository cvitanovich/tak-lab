mtlinterpol(infile,outfile,new_dir)


%mtlinterpol(infile,outfile,new_dir)


[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
first_line,last_line,n_directions,...
           comment1, comment2] = mtlrh(infile)


% compare header

index=0;
%if direction_matrix
if ((info_blocks-2) > 0)
   direction_matrix=mtlrdir(infile,n_directions);
   if (index ~= 0) direction_matrix(:) = direction_matrix(:,index); end;
   % if time_delay exists read time_delay
   if (( (info_blocks-2)*256 - n_directions*8) == (n_channels*4))
      [delay_l,delay_r]=mtlrdel(infile,n_directions,binaural);
     if (index ~= 0) 
        delay_l = delay_l(index); 
        if (binaural == 1) delay_r = delay_r(index); end;
     end;
   end;
end;


out_file_type=file_type;
out_info_blocks=info_blocks;
out_n_channels=n_channels;
out_n_lines=n_lines;
out_sampling_rate=sampling_rate;
out_first_line=first_line;
out_last_line=last_line;
out_n_directions=n_directions;
out_comment1=comment1;
out_comment2=comment2;

mtlwh(outfile,out_filetype,out_info_blocks,out_n_channels,out_n_lines,... out_sampling_rate, out_first_line,out_last_line,out_n_directions,...
out_comment1, out_comment2);



impulse=interpolate (ele,azi,direction_matrix,filename)

%readh
%read direction matrix d
p=[ele;azi]
dist=sphdist(p,d)
%if in dist 0 find that position directions allready exits
[dist_sort,i]=sort(dist)

dist_sum=sum(dist_sort(1:3))

ir=ones(1,n_lines)*(1+i);
for n=1:3
 actual_dist=dist_sum-dist_sort(n)
 old_ir=mtlrch(filename,i(n))  
 %verschiebung um Anfangszeitpunkt 
 ir=ir+old_ir*actual_dist/dist_sum; 
 %Mittelung der Anfangszeitpunkte
 %verschiebung
end

write impulse response


