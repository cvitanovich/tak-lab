function  s2m(filename)

%function  mtlorbo(filename)

SIGM_BLOCK_LEN = 256;
DOS_EOF = 26;

SO_EOF = 79;
SO_TYPE = 80;
SO_N = 128;
SO_CHA = 84;

INT_TYPE = 1;
FLOAT_TYPE = 2;
COMPLEX_TYPE = 3;
ADDA_TYPE = 4;


BIFI_MAGIC_NUMBER = hex2dec('7caffee1');

% status mask

BIFI_MEASURED     = 1;
BIFI_X_TESTSIGNAL = 2;
BIFI_SPEAKER_EQ   = 4;
BIFI_TOTAL_EQ     = 8;

BIFI_RES_NUMBER = 5;

BIFI_RES_AS_MEASURED         = 1;
BIFI_RES_CONSTANT_AZIMUTH    = 2;
BIFI_RES_CONSTANT_DENSITY    = 3;
BIFI_RES_OLD_FULL_SPHERE     = 4;
BIFI_RES_WITH_DIRECTION_INFO = 5;


%open file for reading and appending new data

[fid,message] = fopen(filename,'r+','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;






stat=fseek(fid,220,'bof');
if (stat < 0) error('ERROR: mtlwh, magic number, fseek'); return; end;
%count = fwrite(fid,BIFI_MAGIC_NUMBER,'ulong');
count = fwrite(fid,122,'ulong');

%if (count < length(BIFI_MAGIC_NUMBER)) error('ERROR: mtlwh, magic number, fwrite'); return; end;

end;



status=1;
stat=fseek(fid,232,'bof');
if (stat < 0) error('ERROR: mtlwh, status, fseek'); return; end;
count = fwrite(fid,status,'ushort');
if (count < length(status)) error('ERROR: mtlwh, status , fwrite'); return; end;

end;


res_type=5
stat=fseek(fid,234,'bof');
if (stat < 0) error('ERROR: mtlwh, res_type, fseek'); return; end;
count = fwrite(fid,res_type,'ushort');
if (count < length(res_type)) error('ERROR: mtlwh, res_type , fwrite'); return; end;

end;

channels_per_dir=2;
stat=fseek(fid,238,'bof');
if (stat < 0) error('ERROR: mtlwh, channels per direction, fseek'); return; end;
count = fwrite(fid,channels_per_dir,'ushort');
if (count < length(channels_per_dir)) error('ERROR: mtlwh, channels per direction , fwrite'); return; end;

end;





%close file

stat=fclose(fid);
if (stat < 0) error('ERROR: mtlwh, close file'); return; end;

return;



