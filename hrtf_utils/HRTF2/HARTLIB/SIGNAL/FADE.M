function [h] = fade(h, b1, e1, b2, e2);
%FADE   Fading function. Syntax h = FADE(h, b1, e1, b2, e2);
%       h(1:b1)         --> set zero;
%       h(b1:e1)        --> fade on with squared cosine function f1
%       h(e1:b2)        --> unchanged
%       h(b2:e2)        --> fade off with squared cosine function f2
%       h(e2:length(h)) --> zero
%
%       f1 = cos(pi/2*(e1-i)/(e1-b1))^2
%       f2 = cos(pi/2*(i-b2)/(e2-b2))^2

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 20.06.1994 A. Raab

t=size(h);
h=h(:)';
h(1:b1) = zeros(size(h(1:b1)));
r=(b1+1):(e1-1);
h(r)=h(r).*cos(pi/2*(e1-r)/(e1-b1)).^2;
r=(b2+1):(e2-1);
h(r) = h(r).*cos(pi/2*(r-b2)/(e2-b2)).^2;
h(e2:length(h)) = zeros(size(h(e2:length(h))));
if size(h)~=t h=h'; end;
