function [h,attr]=readhp(catname,channel,mode,display);
% READHP    Read measured headphone impulse response 
%           for given catalogue.
%           [h,attr]=READHP(catname,channel,mode,display);
%           catname   --> catalogue name
%           channel   --> 'left' ([]) or 'right' (default 'left')
%           mode      --> 'raw' returns measured data directly
%                         'norm' or [] returns data normalized 
%                         (peak value=1) and with delay time and offset 
%                         removed.
%           display   --> 'on' shows some informations about the data,
%                         'off' or [] does nothing.
%           h         <-- vector with measured data
%           attr      <-- vector with informations
%                         attr(1)=delay time of data
%                         attr(2)=normalization (peak) value

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 23.08.1994 A. Raab

if nargin<5             % set default values
  mode='n';
  if nargin<4
    channel='l';
  end;
end;

if isempty(channel) channel='l'; end;
if isempty(mode) mode='n'; end;

channel=find(channel(1)==['lr']);
if isempty(channel) error('invalid channel'); end;

mode=find(mode(1)==['nr']);
if isempty(mode) error('invalid mode'); end;

modename=['norm';'raw '];
chname=['left ';'right'];

% load catalogue
%
if ~exist(sprintf('%s.mat',catname)) error('catalogue not found'); end;
eval(['load ' catname]);                

% load impulse response
%
s=sprintf('%s%s',datadir,hpfile);
if ~exist(s) error(sprintf('file %s not found',s)); end;
s=sprintf('[hl,hr]=readsr2(''%s'');',s);
eval(s);

% Check file format
%
if exist('hpoffset')==0 
  disp('warning: old catalogue file format, missing headphone information');
  hpoffset=[0,0];
end;
if exist('hpdelay')==0
  disp('warning: old catalogue file format, missing headphone information');
  hpdelay=[1,1];
end;
if exist('hppeak')==0
  disp('warning: old catalogue file format, missing headphone information');
  hppeak=[max(abs(hl)) max(abs(hr))];
end;

% prepare output data
%
if channel==2 h=hr(:); else h=hl(:); end;
attr=[hpdelay(channel) hppeak(channel)];
if mode==1
  h=h-hpoffset(channel);
  h=h./hppeak(channel);
  h=[h(attr(1):length(h)) ; zeros(attr(1)-1,1)];
end;

% print information if required
%
if strcmp(display,'on')
  disp(sprintf('catalogue:%s, hp %s channel, %s mode',...
       catname,chname(channel,:),modename(mode,:)));
  disp(sprintf('first value (%s)=%g, peak value (raw)=%g',...
       modename(mode,:),h(1),peak));
end;
