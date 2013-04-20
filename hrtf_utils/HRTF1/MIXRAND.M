function x=mixrand(x)
% MIXRAND resorts the elements in vector x in a random way.
%         x=MIXRAND(x)

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 13.06.1994 A. Raab

L=length(x);
ri=round((L-1)*rand(1,L+1)+1);  % random index to vector x
for i=1:L
  t=x(ri(i));                   % exchange two elements
  x(ri(i))=x(ri(i+1));          % at random position
  x(ri(i+1))=t;
end;
