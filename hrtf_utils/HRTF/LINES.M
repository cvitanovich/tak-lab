function lines(x,y,linetype)
%LINES  Plot discrete sequence data.
%       LINES(Y) plots the data sequence Y as stems from the x-axis
%       terminated with dots for the data value.
%       LINES(X,Y) plots the data sequence Y at the values specfied
%       in X.
%       There is an optional final string argument to specify a line-type
%       for the stems of the data sequence.  E.g. LINES(X,Y,'-.') or
%       LINES(Y,':').
%
%       See also PLOT, BAR, STAIRS, STEM.

%       Copyright (c) 1984-93 by The MathWorks, Inc.

%       Changed 23.05.1994 by A. Raab
%       Original file was stem.m
%       I have removed the awful circles!

n = length(x);
if nargin == 1
	y = x(:)';
	x = 1:n;
	linetype = '-';
elseif nargin == 2
	if isstr(y)
		linetype = y;
		y = x(:)';
		x = 1:n;
	else
		x = x(:)';
		y = y(:)';
		linetype = '-';
	end
elseif nargin == 3
	x = x(:)';
	y = y(:)';
end
xx = [x;x;nan*ones(size(x))];
yy = [zeros(1,n);y;nan*ones(size(y))];
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = strcmp(next,'add') & strcmp('add',lower(get(gcf,'NextPlot')));
plot(x,y,'y.',xx(:),yy(:),['y' linetype])
q = axis;hold on;plot([q(1) q(2)],[0 0],'w-');
if ~hold_state, set(cax,'NextPlot',next); end

