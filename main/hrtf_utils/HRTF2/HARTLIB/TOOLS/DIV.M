function c = div(a,b)

% DIV calculates c=a/b for each value of each column of a. b has one column.

% written by Klaus Hartung
% Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 
% Date 27.09.1995

r=1./b;
c = diag(sparse(r)) * a;

