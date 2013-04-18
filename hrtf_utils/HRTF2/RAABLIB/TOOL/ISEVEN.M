function bool=iseven(n);
% ISEVEN return 1 if the input number is even, 0 else
%        bool=ISEVEN(n)
bool=2*fix(n/2)==n;
