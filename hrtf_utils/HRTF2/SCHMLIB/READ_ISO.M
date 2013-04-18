function [isoloaded,criterium,filetype,sampling_line,value,delta,...
          iso_points,iso_array] = read_iso(filename,isoloaded,criterium,...
          filetype,sampling_line,value,delta,iso_points,iso_array)

% function [isoloaded,criterium,filetype,sampling_line,value,delta,...
%          iso_points,iso_array] = read_iso(filename,isoloaded,criterium,...
%          filetype,sampling_line,value,delta,iso_points,iso_array)
%
% READ_ISO reads the datas of an iso-file and returns the following variables:
%
%          - isoloaded: number of loaded isolines
%          - criterium: number of identification for the kind of isoline
%                       (magnitude, phase, phase delay, group delay)
%          - filetype: number of identification between impulse reponses
%                      and transfer functions 
%          - sampling_line: the time- or frequency value in [s] or [Hz],
%                           where the signal is observed
%          - value: vector including the choosen values of the isolines
%          - delta: vector including the tolerances of value vector
%          - iso_points: vector including the number of isopoints belonging
%                        to each isoline. The size of the vector is the
%                        variable 'isoloaded'.
%          - iso_array: matrix including the isopoints. Each column
%                       represents an isopoint. The rows include the
%                       following informations:
%                        1) first row: x-coordinate of first isoline
%                        2) second row: y-coordinate of the first isoline
%                        3) third row: z-coordinate of the first isoline
%                        4) fourth row: x-coordinate of the second isoline
%                        5) ...
%
%          If several isolines are loaded, all of then must have the same
%          criterium, filetype, and sampling line. If not the program aborts
%          with an error message. The informations about an isoline are
%          displayed. The attributes of a new isoline are appended to the
%          already existing isolines. For this reason input and output
%          arguments are the same. The files must have the ISOSIG - format.
%
%          For further information see the ISOSIG manual.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 26.08.1994 Joerg Schmatter
%
% For evaluating and displaying spatial lines of the same level another
% method is used, than the with C++ provided methods. The ISOSIG format is
% not used at the moment !


[fid,message] = fopen(filename,'r','ieee-le');
disp(['file: ',filename]);
if fid==-1
  disp(message);
  return;
end;

eof = fread(fid,1,'char');
eof_int = abs(eof);
if (eof_int ~= 26)
  fclose(fid);
  error('choosen file isn''t readable');
end;

if isoloaded==0
 isoloaded = isoloaded + 1;
 criterium = fread(fid,1,'short');
 filetype = fread(fid,1,'short');
 sampling_line = fread(fid,1,'float32');
 value(isoloaded) = fread(fid,1,'float32');
 delta(isoloaded) = fread(fid,1,'float32');
 iso_points(isoloaded) = fread(fid,1,'short');
 isofield = fread(fid,3*iso_points(isoloaded),'float32');
 isofield = reshape(isofield,3,iso_points(isoloaded));
 fclose(fid);
else
  crit = fread(fid,1,'short');
  if criterium~=crit
    error('same criterium neccessary for multiple isolines');
  end;
  f_type = fread(fid,1,'short');
  if filetype~=f_type
    error('same filetype neccessary for multiple isolines');
  end;
  s_line = fread(fid,1,'float32');
  if sampling_line~=s_line
    error('same sampling line neccessary for multiple isolines');
  end;

  isoloaded = isoloaded + 1;

  value(isoloaded) = fread(fid,1,'float32');
  delta(isoloaded) = fread(fid,1,'float32');
  iso_points(isoloaded) = fread(fid,1,'short');
  isofield = fread(fid,3*iso_points(isoloaded),'float32');
  isofield = reshape(isofield,3,iso_points(isoloaded));
  fclose(fid);
end;

if (isoloaded==1)
  iso_array = isofield;
else
   max_n_o_points = max(iso_points);               % matrices must have the same dimensions
   row = size(iso_array,1);
   col = size(iso_array,2);
   actual_n_o_points = iso_points(isoloaded);
   if actual_n_o_points>col                        % new isoline has more points
      append = zeros(row,actual_n_o_points-col);
      iso_array = [iso_array,append];
   elseif actual_n_o_points<col                    % new isoline has less points
      append = zeros(row,col-actual_n_o_points);
      isofield = [isofield,append];
   end;
   iso_array = [iso_array;isofield];
end;

if criterium==0
  disp('isolines of magnitude');
elseif criterium==1
  disp('isolines of phase');
elseif criterium==2
  disp('isolines of phase delay');
elseif criterium==3
  disp('isolines of group delay');
end;

if filetype==2
  disp('file type: impulse response');
  disp(['sampling line: ' num2str(sampling_line) ' sec']);
  disp(['absolute value: ' num2str(value(isoloaded))]);
  disp(['delta: ' num2str(delta(isoloaded))]);
else 
  disp('file type: transfer function');
  disp(['sampling line: ' num2str(sampling_line) ' Hz']);
  disp(['value: ' num2str(value(isoloaded)) ' dB']);
  disp(['delta: ' num2str(delta(isoloaded)) ' dB']);
end;
