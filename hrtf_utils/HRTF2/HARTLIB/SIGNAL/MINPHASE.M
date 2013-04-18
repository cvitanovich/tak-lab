function H_min = minphase(H_mix, t, p)

% MINPHASE calculates the minimum phase impulse response or transfer function
% function H_min = minphase(H_mix, t, p)
% initial delays are pre
% H_mix: mixed phase impulse response or transfer function
% t: domain of input and output data t=1 time domain, t=0 frequency domain
% p: preserve initial delay(1), pure minimum phase signale might be not causal


% (p) Klaus Hartung 
% Lehrstuhl für allgemeine Elektrotechnik und Akustik
% Ruhr-Universität Bochum
% 28.07.1995


if (nargin<=1) t=0; end;
if (nargin<=2) p=0; end;

if (t==0) h_mix=rifft(H_mix); end;
  else
     h_mix=H_mix;
     H_mix=rfft(h_mix);
  end;

%irplot(h_mix',44100);

len=length(h_mix);
h_min=zeros(1,len);


if (p==1)  [me,pe]=max(etc(h_mix)); end;

A=abs(H_mix);

B_ht=imag(hilbert(A));

B_ht=B_ht-B_ht(1);

H_min=A.*(cos(B_ht)+j*sin(B_ht));

%phasplot(H_min', 44100);

h_min=rifft(H_min);

%irplot(h_min',44100);


%irplot(h_min',44100);

if (p==1) 

  x=[0:1:len/2-1 -len/2:1:-1];
  x=x+len/2+1;

   h_min=h_min(x);

  [me,pe_m]=max(etc(h_min));
  
%  irplot(h_min',44100);

  x=[pe_m-pe+1:1:len 1:1:pe_m-pe];
  h_min=h_min(x);
%  irplot(h_min',44100);
end;
  [me,pe_m]=max(etc(h_min));

if (t==0) H_min=rfft(h_min); end;
else H_min=h_min; end;


end;

