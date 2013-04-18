function fresample(infile,outfile,old_sf,new_sf,stereo)

%function fresample(infile,outfile,old_sf,new_sf,stereo)

ratio=old_sf/new_sf;
block=old_sf;
i=zeros(1,100);
[o,b]=resample(i,new_sf,old_sf);
a=1;
[ifid,message]=fopen(infile,'r');
if ifid==-1          % file can't be opened
  disp(message);
  return;
end;


[ofid,message]=fopen(outfile,'w');
if ofid==-1          % file can't be opened
  disp(message);
  return;
end;


[input,count]=fread(ifid,block,'int16');

zi1=filtic(b,a,count/2);
zi2=filtic(b,a,count/2);
zi=filtic(b,a,count/2);

while(count>0) 
if (stereo)
  i1=input(1:2:length(input));
  [o1,zf1]=filter (b,a,i1,zi1);
  o1=o1(1:ratio:length(o1));
  zi1=zf1;
  i2=input(2:2:length(input));
  [o2,zf2]=filter (b,a,i2,zi2);
  o2=o2(1:ratio:length(o2));
  zi2=zf2;

  output=zeros(1,count/ratio);
  
  output(1:2:length(o1)*2)= o1;
  output(2:2:length(o1)*2)= o2;
 else   
   [output,zf]=filter(b,a,input,zi);
   zi=zf;
  output=output(1:2:length(output));
 end
 
 fwrite(ofid,output,'int16');

 [input,count]=fread(ifid,block,'int16');
end;



fclose(ifid);
fclose(ofid);
end;
