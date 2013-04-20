function okay=isstable(p);
% ISSTABLE      Check polynom with coefficients in p for roots 
%               outside or on the unit circle. 
%               Return 1 if none, 0 else.

okay=~any(abs(lattice(p))>=1);ÿ
