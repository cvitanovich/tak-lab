function [flag,fname]=chkfile(fname,extension)
% CHKFILE checks for the existence of the specified file
%         [flag,fname]=chkfile(fname,extension)
%         A warning message is diplayed if the file does not exist.
%         parameters:
%           fname     --> file main name
%           extension --> file name extension, i.e '.mat'
%           fname     <-- becomes a '?' if file does not exist, remains
%                         unchanged if file exists
%           flag      --> returns 1 if file exists, 0 else

if exist(sprintf('%s%s',fname,extension))~=2
  flag=0;
  disp(sprintf('file %s not found, press any key ...',fname));
  fname='?';
  pause;
else
  flag=1;
end;
