function mtlinterpol(infile,outfile,new_dir)


%mtlinterpol(infile,outfile,new_dir)


[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
first_line,last_line,n_directions,...
 comment1, comment2] = mtlrh(infile);



binaural=0;
if ((n_channels/n_directions) == 2) binaural=1; end;


index=0;
%if direction_matrix
if ((info_blocks-2) > 0)
   direction_matrix=mtlrdir(infile,n_directions);
   if (index ~= 0) direction_matrix(:) = direction_matrix(:,index); end;
   % if time_delay exists read time_delay
   if (( (info_blocks-2)*256 - n_directions*8) >= (n_channels*4))
           [delay_l,delay_r]=mtlrdel(infile,n_directions);
           
          if (index ~= 0) 
      %  delay_l = delay_l(index); 
    %    if (binaural == 1) delay_r = delay_r(index); end;
     end;
   end;
end;
out_ft=file_type;
out_info_blocks=info_blocks;

out_n_lines=n_lines;
out_sampling_rate=sampling_rate;
out_first_line=first_line;
out_last_line=last_line;
new_ndir=size(new_dir);
n_dir=new_ndir(2);
out_n_directions=n_dir;
out_comment1='test';
out_comment2='test';
out_n_ch=out_n_directions*(binaural+1);
mtlwh(outfile,out_ft,out_info_blocks,out_n_ch,out_n_lines,... 
out_sampling_rate, out_first_line,out_last_line,...
out_n_directions,'test', 'test');



if (n_dir > 0)
  n_byte=n_dir*2*4;
  sn_del=length(delay_l);
  if (sn_del > 0) n_byte=n_byte+n_dir*4*(binaural+1); end; 
  n_blocks=ceil(n_byte/256);
end;

out_info_blocks=2+n_blocks;

mtlwh(outfile,out_ft,out_info_blocks,out_n_ch,out_n_lines,... 
out_sampling_rate, out_first_line,out_last_line,...
out_n_directions,out_comment1, out_comment2);


%if direction_matrix
[fid,message] = fopen(outfile,'a','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

if ((n_blocks) > 0)
 fseek(fid,512,'bof');
 fwrite(fid,zeros(n_blocks*256,1),'float32');
 fclose(fid);
 
 n_delay_l=zeros(1,out_n_directions);
 n_delay_r=zeros(1,out_n_directions);

 status=mtlwdir(outfile,new_dir);
 if(length(delay_l) > 0) mtlwdel(outfile, n_delay_l,n_delay_r); end;
end


for n=1:out_n_directions
clc
  disp(['MTLINTERPOL: Calculating direction ', num2str(n), ' of ', num2str(out_n_directions)] );

p=new_dir(:,n)

dist=sphdist(p,direction_matrix);

 if any(dist==0) 
    index=find(dist==0);
 
index
    channel1=mtlrch(infile,(index-1)*(binaural+1)+1);
    mtlwch(outfile,channel1,(n-1)*(binaural+1)+1);

    n_del_l(n)=delay_l(index);
    if (binaural==1) 
    channel2=mtlrch(infile,index*(1+binaural));
    mtlwch(outfile,channel2,n*(binaural+1));
  
    n_del_r(n)=delay_r(index);
 
   end
else

[dist_sort,i]=sort(dist);
dist_sort(1:4)
direction_matrix(1:2,i(1:4))

dist_sum=sum(dist_sort(1:4));


mean_del_r=0;
mean_del_l=0;
sum1=zeros(n_lines,1);
sum2=zeros(n_lines,1);

for j=1:4
 actual_dist=dist_sum-dist_sort(j);
 del_l=delay_l(i(j));
 mean_del_l=fix(mean_del_l+del_l*actual_dist/(3*dist_sum));
 n_del_l(n)=mean_del_l;
 old_channel1=mtlrch(infile,(i(j)-1)*(binaural+1)+1);
 sum1(1:n_lines-del_l+1)=sum1(1:n_lines-del_l+1)+old_channel1(del_l:n_lines)*actual_dist/(3*dist_sum);
 if (binaural==1) 
  del_r=delay_r(i(j));
  mean_del_r=fix(mean_del_r+del_r*actual_dist/(3*dist_sum));
  n_del_r(n)=mean_del_r;
  old_channel2=mtlrch(infile,i(j)*(1+binaural));
    sum2(1:n_lines-del_r+1)=sum2(1:n_lines-del_r+1)+old_channel2(del_r:n_lines)*actual_dist/(3*dist_sum);
 end
end
sum1=sum1';
sum1=[zeros(1,mean_del_l), sum1(1:n_lines-mean_del_l)];
mtlwch(outfile,sum1,(n-1)*(binaural+1)+1);

 if (binaural==1) 
  sum2=sum2';
  sum2=[zeros(1,mean_del_r), sum2(1:n_lines-mean_del_r)];
  mtlwch(outfile,sum2,n*(1+binaural));
 end
end
end

[fid,message] = fopen(outfile,'a','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

if(length(n_del_l) > 0) mtlwdel(outfile, n_del_r,n_del_r); end;
end



