function [h,attr]=readdata(catname,azimuth,elevation,channel,mode,display);
% READDATA  Read measured impulse response from catalogue.
%           [h,attr]=READDATA(catname,azimuth,elevation,channel,mode,display);
%           catname   --> catalogue name
%           azimuth   --> angle in degrees (0,15,30,...,345)
%           elevation --> ignored in this version
%           channel   --> 'left' ([]) or 'right' (default 'left')
%           mode      --> 'raw' returns measured data directly
%                         'norm' or [] returns data normalized 
%                         (peak value=1) and with delay time and reflections removed.
%           display   --> 'on' shows some informations about the data,
%                         'off' or [] does nothing.
%           h         <-- vector with measured data
%           attr      <-- vector with informations
%                         attr(1)=delay time of data
%                         attr(2)=normalization (peak) value
%                         attr(3)=interaural delay (not implemented)
%                         attr(4)=start index for fading (only for newest catalogues)
%                         attr(5)=stop index for fading (only for newest catalogues)

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

% validate specified direction and get corresponding filename
%
p=find(flt(:,1)==azimuth & flt(:,2)==elevation);
if isempty(p) error('direction not in catalogue'); end;
p=p(1);
filename=deblank(setstr(flt(p,3:size(flt,2))));
s=sprintf('[hl,hr]=read%s(''%s%s.%s'');',filetype,datadir,filename,filetype);
eval(s);

% prepare output data
%
if channel==2 h=hr(:); else h=hl(:); end;

attr=[delay(channel,p);peak(channel,p);idelay(p)];
if mode==1
  h=h-offset(channel,p);
  h=h./peak(channel,p);
  h=[h(attr(1):length(h)) ; zeros(attr(1)-1,1)];
  if exist('fade1')==1 & exist('fade2')==1
    h=fade(h,0,0,fade1,fade2);
    attr=[attr;fade1;fade2];
  end;
end;

% print information if required
%
if strcmp(display,'on')
  disp(sprintf('catalogue:%s, %s channel, %s mode, azimuth=%g, elevation=%g',...
       catname,chname(channel,:),modename(mode,:),azimuth,elevation));
  disp(sprintf('initial delay=%g samples, max ccf=%g samples',...
       attr(1),attr(3)));
  disp(sprintf('first value (%s)=%g, peak value (raw)=%g',...
       modename(mode,:),h(1),attr(2)));
end;ÿ
