function en=energy2(x,w)

% ENERGY2 calculate energy of array, integrate power in window
%
%        ENERGY2(X) returns energy of x
%        ENERGY2(X,W) returns energy in window w shifted over x


% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektortechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%
% Date: 29.5.96
% Modified: 7.7.96
%        improved help
%  

% Calls		toolbox
%
% power         hartung/tools

if (nargin < 1) error('No time signal defined'); end;



p=power(x);
[m,n] = size(x);

if (nargin<2)
   en = zeros(1,n);
   en = sum(p);
  elseif ((w>1) & (w<m))
    en = zeros(m-w+1,n);
    disp('energy_: Computing sum in window for each position. This will take a while...'); 
    for l = 1:m-w+1
%    clc
%    disp(['energy_: Computing line ' num2str(l) ' of ' num2str(m-w+1)]);
      en(l,:) = sum(p(l:l+w-1,:));
    end
  elseif (w==1) en=p;

  elseif (w==m) en=sum(p);

 else
   error(' 0 < W <= m; [m,n] = size(x).'); 
end;


% end of file energy2.m