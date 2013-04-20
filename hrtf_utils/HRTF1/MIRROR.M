function out = mirror(n_lines,in)

% function out = mirror(n_lines,in)
%
% MIRROR conjugate complex mirroring of the matrix 'in'.
%        'n_lines' is the number of columns of the matrix. If n_lines is even the 
%        resulting output matrix has 2 * n_lines columns, if n_lines is odd
%        the resulting output matrix, has 2 * (n_lines-1) columns.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter

out = [];
if rem(n_lines,2)~=0
  out(1:n_lines) = in;
  out(n_lines+1:2*n_lines-1) = conj(fliplr(in(1:n_lines-1)));
else
  out(1:n_lines) = in;
  out(n_lines+1:2*n_lines) = conj(fliplr(in(1:n_lines)));
end;
