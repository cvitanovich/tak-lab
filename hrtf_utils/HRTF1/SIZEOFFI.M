function s=sizeoffi(filename)

% SIZEOFFI returns size of file in byte 
%
%	s=SIZEOFFILE(filename)
%       filename: name of file
%	s: 	  size of file in bytes


% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date: 25.09.1995
% Modified: 19.1.1996 Klaus Hartung
% Modified: 7.7.1996 Klaus Hartung
% 	   help improved
 
% calls		toolbox
% no calls to other routines

if (nargin < 1) error('Insufficient number of arguments'); end;

fid=fopen(filename);     
if (fid<0) 
  error('D: fopen'); s=-1; 
  return; 
end;

status=fseek(fid,0,'eof');    
       
if (status<0) error('D: fseek'); s=-1; fclose(fid); return; end;

s=ftell(fid);             
if (s<0) 
  error('D: ftell'); 
  s=-1; 
  message=ferror(fid);
  error(message);
end;

fclose(fid);

% end of sizeoff.m 
    
