function [H, h] = mklsrdat(N, M, L, r, zlag)
% MKLSRDATA make data arrays for least-square approximation with
%                autocorrelation data
% [H, h] = MKLSRDAT(N, M, L, r, zlag)

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 27.05.1994 A. Raab

% form all vectors and matrices
%
% h=[r(N+1) ... r(N+L)]'
%
h = r((zlag+N+1):(zlag+N+L))';  
%
% form matrix H1
%        | r(N)     r(N-1)   ... r(N-M+1) |
%        | r(N+1)   r(N)     ... r(N-M+2) |
% H =    | r(N+2)   r(N+1)   ... r(N-M+3) |
%        |   :        :            :      |
%        | r(N+L-1) r(N+L-2) ... r(N+L-M) |
%
for col=0:M-1
  H = [H r((zlag+N-col):(zlag+N+L-1-col))'];
end;

