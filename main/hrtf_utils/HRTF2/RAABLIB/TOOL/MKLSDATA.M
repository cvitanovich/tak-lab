function [H1, h1, H2, ofs] = mklsdata(N, M, L, h)
% MKLSDATA make data arrays for least-square approximation
%          [H1, h1, H2, ofs] = MKLSDATA(N, M, L, h)

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 27.05.1994 A. Raab

h=h(:);
ofs = 1;
if N-M+1 < 0                   % insert zeros if necessary
  ofs  = -(N-M)+1;             % store start of original data
  h = [zeros((ofs-1),1); h];   % insert N-M zeros
end;
%
% form all vectors and matrices
%
% h1=[h(N+1) ... h(N+L)]'
%
h1 = h((ofs+N+1):(ofs+N+L));
%
% form matrix H1
%        | h(N)     h(N-1)   ... h(N-M+1) |
%        | h(N+1)   h(N)     ... h(N-M+2) |
% H1 =   | h(N+2)   h(N+1)   ... h(N-M+3) |
%        |   :       :             :      |
%        | h(N+L-1) h(N+L-2) ... h(N+L-M) |
% note: N-M+1 might be negative, but we have inserted some zeros,
% see above
%
for col=0:M-1
  H1 = [H1 h((ofs+N-col):(ofs+N+L-1-col))];
end;
%
% form matrix H2
%        | h(0) 0      0      ... 0      |
%        | h(1) h(0)   0      ... 0      |
% H2 =   | h(2) h(1)   h(0)   ... 0      |
%        |  :    :      :         :      |
%        | h(N) h(N-1) h(N-2) ... h(N-M) |
%

K=min([N M]);
for col=0:K
  t = [zeros(col,1) ; h(ofs:(ofs+N-col))];    % temporary column vector
  H2 = [H2  t];
end;                                          % H2 is now K x K
clear t;
if (M-N) > 0                     % fill the remaining
  H2 = [H2 zeros(N+1,M-N)];       % columns with zero
end;
