function y=optwin(x,w,m_max)

% OPTWIN find the most powerful part 
%
%        y=optwin(x,w,m)
%        y=optwin(x);
%        y=optwin(x,w);
%        optwin returns the most powerful part 
%        of x in a rectangular window of length w in a 
%        new vector with length m

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date 25.09.1995
% Modified: 7.7.96
%        improved help

% Calls 	Toolbox
% energy2	

if (nargin < 1) error('Insufficient number of arguments'); end;
if (nargin < 2) w=1; end;
[m,n] = size(x);
if (nargin < 3) m_max=m; end;


e = energy2(x,w);

[e_max,k_e_max] = max(e);

k_e_max_min = min(k_e_max);

y=zeros(m_max,n);



if (~any(k_e_max==1))
  new_k1 = k_e_max - k_e_max_min + 1;
  new_k2 = new_k1+w-1;

  [max_new_k2, maxpos]= max(new_k2);

  if ( max_new_k2 > m_max) error(['m_max ' int2str(max_new_k2) ' at ' int2str(maxpos) ' is smaller than ' int2str(m_max) '.' ]);  end;
  
  for n_a=1:n
    % clc
    % disp(['optwin: Computing channel ' num2str(n_a) '/' num2str(n)]);
      y(new_k1(n_a):new_k2(n_a),n_a) = x(k_e_max(n_a):k_e_max(n_a)+w-1,n_a);
  end

 else
%   disp('No touch');
   y=x(1:m_max,:);
end;


% end of file optwin.m 