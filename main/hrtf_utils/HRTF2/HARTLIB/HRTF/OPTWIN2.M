function y=optwin2(x,w,m_max)

% OPTWIN2 find the most powerful part 
%
%        y=optwin(x,w,m)
%        y=optwin(x);
%        y=optwin(x,w);
%        optwin returns the most powerful part 
%        of x in a rectangular window of length w in a 
%        new vector with length m
%        w is ignored only the real length is considered

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

%wl=128;
wl=w;
FIR_ORDER=zeros(n,1);
N1=zeros(n,1);
N2=zeros(n,1);
e=zeros(m+wl,n);
for n_count=1:n
   disp(['Calculating channel ' int2str(n_count) ' of ' int2str(n) '.'])
   [FIR_ORDER(n_count),N1(n_count),N2(n_count),e(:,n_count)]=ordfir(x(:,n_count),wl,0.02,1.0);
   disp([int2str(n_count)  ':' int2str(FIR_ORDER(n_count))]);
   if (FIR_ORDER(n_count)>wl)
      FIR_ORDER(n_count)=N1(n_count);
      disp(['Corrected: ' int2str(n_count)  ':' int2str(FIR_ORDER(n_count))]);
   end
   if (FIR_ORDER(n_count)>wl)
      FIR_ORDER(n_count)=N2(n_count);
      disp(['Corrected: ' int2str(n_count)  ':' int2str(FIR_ORDER(n_count))]);
   end     
   if (FIR_ORDER(n_count)>wl)
      disp(['Something wrong with Impulse Response  ' int2str(n_count)  '.']);
   end 
end

   irplot(e(:,n_count),44100);

   N1
N2
end
FIR_ORDER_MAX=max(FIR_ORDER)
if (FIR_ORDER_MAX>wl) error('Length of impulse response greater than wl. Contact Klaus Hartung.'); end;
e2=zeros(n+FIR_ORDER_MAX,n); 


for n_count=1:n
	ed= energy2(x(:,n_count),fix(FIR_ORDER(n_count)));
	e2(1:length(ed),n_count)=ed;
end;


[e_max,k_e_max] = max(e2);

k_e_max=k_e_max-2;
k_e_max_min = min(k_e_max);

y=zeros(m_max,n);



if (~any(k_e_max==1))
  new_k1 = k_e_max - k_e_max_min + 1;
%  new_k2 = new_k1+w-1;
  new_k2 = new_k1+FIR_ORDER'-1;
  [max_new_k2, maxpos]= max(new_k2);

  if ( max_new_k2 > m_max) error(['m_max ' int2str(max_new_k2) ' at ' int2str(maxpos) ' is smaller than ' int2str(m_max) '.' ]);  end;
  
  for n_a=1:n
    % clc
    % disp(['optwin: Computing channel ' num2str(n_a) '/' num2str(n)]);
      y(new_k1(n_a):new_k2(n_a),n_a) = x(k_e_max(n_a):k_e_max(n_a)+FIR_ORDER(n_a)-1,n_a);
  end

 else
%   disp('No touch');
   y=x(1:m_max,:);
end;


% end of file optwin.m 
