function l_contr(cs)

% function l_contr(cs)
%
% L_CONTR adds height labels to the current contour plot using the
%         contour structure 'cs', the returning output from the CONTOUR
%         routine. The function is based on the M-file 'clabel()', but uses
%         less points of labeling. The position of the labels is choosen
%         randomly.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 02.08.1994 Joerg Schmatter

if min(size(cs)) > 2
	error('First input must be a valid contour description matrix.')
end
cax = gca;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

hold on;
[mcs, ncs] = size(cs);

k = 0; n = 0; flip = 0;

while (1)
  % Select a labeling point randomly.
   
  k = k + n + 1; if k > ncs, break, end
  c = cs(1, k); n = cs(2, k);
  r = rand; while (r > 1 | r < 0) r = rand; end
  j = fix(r.* (n - 1)) + 1;
  if flip, j = n - j; end
  flip = ~flip;
  x1 = cs(1, j+k); y1 = cs(2, j+k);
  x2 = cs(1, j+k+1); y2 = cs(2, j+k+1);
  xx = (x1 + x2) ./ 2; yy = (y1 + y2) ./ 2;   % Test was here; removed.

  % Label the point.
  if r>0.8
    s = [' ' sprintf('%0.3g', c)];
    plot(xx, yy, '+');
    text(xx, yy, s, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');
  end;
end

if ~hold_state, set(cax,'NextPlot',next); end

