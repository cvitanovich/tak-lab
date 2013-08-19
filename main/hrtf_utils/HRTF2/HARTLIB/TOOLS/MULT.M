function c = mult(a,b)

% MULT calculates c=a*b for each value of each column of a. b has one column

% written by Klaus Hartung
% Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 
% Date 27.09.1995

c = diag(sparse(b)) * a;

