function [data,main_chunk,length1,chunk_type,sub_chunk,length2,format,modus,sample_fq,byte_p_sec,...
          byte_p_spl,bit_p_spl,data_chunk,data_length]=readwav(fname,dispmode);


if nargin<2
  dispmode=0;
end;

[fid,msg]=fopen(fname,'r');
if fid<0 error(msg); end;

main_chunk=fread(fid,4,'char');
length1=fread(fid,1,'ulong');
chunk_type=fread(fid,4,'char');
sub_chunk=fread(fid,4,'char');
length2=fread(fid,1,'ulong');
format=fread(fid,1,'ushort');
modus=fread(fid,1,'ushort');
sample_fq=fread(fid,1,'ulong');
byte_p_sec=fread(fid,1,'ulong');
byte_p_spl=fread(fid,1,'ushort');
bit_p_spl=fread(fid,1,'ushort');
data_chunk=fread(fid,4,'char');
data_length=fread(fid,1,'ulong');

if byte_p_spl==1
  [data,cnt]=fread(fid,inf,'uchar');
elseif byte_p_spl==2
  [data,cnt]=fread(fid,inf,'ushort');
else
  error('invalid byte_p_spl value');
end;
if ~iseven(cnt) data=data(1:length(data)-1); end;

if modus==2
  L=length(data);
  data1=data(1:2:L);
  data2=data(2:2:L);
  data=[data1(:),data2(:)];
end;

if dispmode==1
  disp(setstr(main_chunk'));
  disp(length1);
  disp(setstr(chunk_type'));
  disp(setstr(sub_chunk'));
  disp(length2);
  disp(format);
  disp(modus);
  disp(sample_fq);
  disp(byte_p_sec);
  disp(byte_p_spl);
  disp(bit_p_spl);
  disp(setstr(data_chunk'));
  disp(data_length);
end;




