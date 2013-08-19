function status=mtlseek(fid, offset, origin)

if (offset > 32767)
 n_seeks=fix(offset/32767);
 remainder=rem(offset,32767);
 status=fseek(fid,0,'bof');
 if (status<0) error('mtlseek not succesful -1'); end;
 for n=1:n_seeks
  status= fseek(fid,32767,'cof');
  if (status<0) error('mtlseek not succesful -2' ); end;
 end
 status=fseek(fid,remainder,'cof');
 if (status<0) error('mtlseek not succesful -3' ); end;

else
 status=fseek(fid,offset,origin);
 if (status<0) error('mtlseek not succesful -4' ); end;

end
return;
