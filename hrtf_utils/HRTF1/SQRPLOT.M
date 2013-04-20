function sqrplot()
% SQRPLOT produces equal sized axis for a graph.

ax=axis;
minxy=min([ax(1) ax(3)]);
maxxy=max([ax(2) ax(4)]);
axis([minxy maxxy minxy maxxy]);
