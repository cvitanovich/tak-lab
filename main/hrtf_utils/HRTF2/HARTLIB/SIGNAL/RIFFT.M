function x = rifft(X,nf)

% RIFFT inverse fft of symmetric spectrum 
% 	rifft(X,nf)
%       rifft(X)
% 	Performs a IFFT of a symmetric spectrum of real valued time signals
% 	rifft(X) assumes that Nyquist frequency is included
% 	rifft(X,nf) nf=0 if not included
% 


% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%

% Date 25.09.1995
% Modified: 26.09.1995 by Klaus Hartung
% Modified: 9.7.96 Klaus Hartung 

% Calls 	Toolbox
% even		hartung/tools
% 

if (nargin < 1) error('No Spectrum defined'); end;
if (nargin < 2) nf=1; end;
	
[m,n] = size(X);

if nf   
	Y = zeros(2*m-2,1);
	x=zeros(2*m-2,n); 
else 
	Y = zeros(2*m-1,1); 
	x=zeros(2*m-1,n);
end



%n_x   = ceil((m-1)*2)

if (m==1)  X=X(:); [m,n] = size(X); end;

   
if(even(m))     
%   disp('Set Nyquist to 0.')
%    Y  = zeros(m+1,n); 
%    Y(1:m,:)  = X ;
%disp('RIFFT: Warning assuming nyquist frequency is 0')
    nf = 1;   
end

%disp('mirroring');
for nc=1:n
   if nf                         % with Nyquist frequency      
%    disp('WIth Nyquist');
%    X   = X/2;
%    X(1) = X(1)*2;
%    X(m) = X(m)*2;
%     Y    = zeros(2*m-2,1);
     if (even(m))      Y    = [X(1:m,nc) ; 0; conj(X(m:-1:2,nc))];
	x=zeros(2*m,n);

       else            Y    = [X(1:m,nc) ; conj(X(m-1:-1:2,nc))]; end;
   else                        % without Nyquist frequency    
%    disp('Without Nyquist');
%    X    = X/2;
%    X(1) = X(1)*2;
%    Y   = zeros(2*m-1,1);
     Y    = [X ; conj(X(m:-1:2,nc))];
   end
%disp('fft');

   x(:,nc) = real(ifft(Y));
end

%x = x*n_x;


% end of file rifft.m