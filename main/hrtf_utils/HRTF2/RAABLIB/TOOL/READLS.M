function [h,attr]=readls(catname,number,mode,display);
% READLS    Read measured loudspeaker impulse response from catalogue.
%           [h,attr]=READLS(catname,number,mode,display);
%           catname   --> catalogue name
%           number    --> number of loudspeaker (1..11)
%           mode      --> 'raw' returns measured data directly
%                         'norm' or [] returns data normalized
%                         (peak value=1) and with delay time removed.
%           display   --> 'on' shows some informations about the data,
%                       'off' or [] does nothing.
%           h         <-- vector with measured data
%           attr      <-- vector with informations
%                         attr(1)=delay time of data
%                         attr(2)=normalization (peak) value

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 05.09.1994 A. Raab

if nargin<4              % set default values
  display='off';
  if nargin<3
    mode='n';
  end;
end;

if isempty(mode) mode='n'; end;
if isempty(display) display='off'; end;

mode=find(mode(1)==['nr']);
if isempty(mode) error('invalid mode'); end;

modename=['norm';'raw '];


% load catalogue
%
if ~exist(sprintf('%s.mat',catname)) error('catalogue not found'); end;
eval(['load ' catname]);

% load data file and validate specified number
%
maxnum=size(flt,1);
if isempty(find(number==[1:maxnum])) error('number not in catalogue'); end;
filename=deblank(setstr(flt(number,:)));
s=sprintf('h=readsr1(''%s%s'');',datadir,filename);
eval(s);

% prepare output data
%
attr=[delay(number);peak(number)];
if mode==1
  h=h-offset(number);
  h=h./peak(number);
  h=[h(attr(1):length(h)) ; zeros(attr(1)-1,1)];
end;

% print information if required
%
if strcmp(display,'on')
  disp(sprintf('catalogue:%s, speaker no. %i, %s mode',...
       catname,number, modename(mode,:)));
  disp(sprintf('initial delay=%g samples',attr(1)));
  disp(sprintf('first value (%s)=%g, peak value (raw)=%g',...
       modename(mode,:),h(1),attr(2)));
end;
