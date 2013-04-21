function status = exist1(S,scr)

%function status = exist1('S',scr)
%
% extends MATLAB's exist.m for use with structures
% see exist.m for 'help'
% also supported for individual fields of a structure
% full of kluges!!!

warning('off')
X = [''' '''];
SS = [X(1) S X(3)];
pt = findstr(S,'.');
len = length(S);
if isempty(pt)		pt = len +1;	end

if nargin ==1						% test to see if it is there at all
   try
      temp = evalin('caller',[ 'exist(' SS ')' ]);
   catch
      status = 0;
      return
   end
elseif nargin ==2				% test to see if there as a particular type
   if strcmp(scr,'var')
      try
         temp = evalin('caller',[ 'exist(' SS ' ,' [X(1) 'var' X(3)] ')' ]);
      catch
         status = 0;
         return
      end
   elseif strcmp(scr,'builtin')
      status = evalin('caller',[ 'exist(' SS ',' [X(1) 'builtin' X(3)] ')' ]);
	   return
   elseif strcmp(scr,'file')
      status = evalin('caller',[ 'exist(' SS ',' [X(1) 'file' X(3)] ')' ]);
  		 return
   elseif strcmp(scr,'dir')
     status = evalin('caller',[ 'exist(' SS ',' [X(1) 'dir' X(3)] ')' ]);
   	return
   else
      help exist;		disp('exist2 2nd argment error');
   return
   end
else										% needs 1 or 2 argins
   disp('exist2 requires either one or two argins')
   help exist2
   return
end


if temp == 0								% claims it doesn't exist; is it a structure?
   S_name = S(1:pt-1);			% structure name
   F_name = S(pt+1:len);		% field name
   if ~isempty(F_name)			% if a field name was supplied in argin1
      try
      if isfield(evalin('caller',S_name),F_name)
      F = getfield(evalin('caller',S_name),F_name);
      if ~isempty(F);			% field with data
         status = 1;
      else							% empty field
         status = 0;
      end
     else							% not a field
      status = 0;
   	end
   	catch
      status = 0;
    end
  else								% no field name is it a structure?
     try
        status = evalin('caller',['isstruct(' S_name ')']);
     catch
        status = 0;
     end
  end
else									% exists and not a structure
   status = temp;
end

warning('on')