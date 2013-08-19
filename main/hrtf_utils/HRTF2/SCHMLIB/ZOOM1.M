function zoom1(xaxis)

% function zoom1(xaxis)
%
% ZOOM1 is an user menu for selecting ranges in a figure including only
%       one plot. If MATLAB runs under WINDOWS the following menu is
%       attached to a plot window:
%	Set X-Range	: use mouse to select displayed x-range
%	Set Y-Range	: use mouse to select displayed y-range
%	Set XY-Range    : use mouse to select displayed x- and y-range
%	Autoscale       : autoscale plot
%	Grid on/off     : set grid on/off
%
%       The input argument 'xaxis' is a vector including the points of
%       sampling.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

%if computer == 'PCWIN'
  item = uimenu('Label','Axis-Control');
         uimenu(item,'Label','Set X-Range ','Callback','zoom1_(1,xaxis)');
         uimenu(item,'Label','Set Y-Range ','Callback','zoom1_(2,xaxis)');
         uimenu(item,'Label','Set XY-Range','Callback','zoom1_(3,xaxis)');
         uimenu(item,'Label','Autoscale','Callback','zoom1_(4,xaxis)');
         uimenu(item,'Label','Grid on/off','Callback','zoom1_(5,xaxis)');
%end;
