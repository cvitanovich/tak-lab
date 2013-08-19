function X = rfft(x,len)

% RFFT calculates the spectrum of a real valued time series
%
% 	X = rfft(x)
% 	X = rfft(x,len)
%
% 	Performs a FFT on real valued x vectors or matrices. The symmetric part of the
% 	transform is cut and only one half is preserved.
% 	For len see FFT.

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date 25.09.1995
% Modified: 26.09.1995  Klaus Hartung
% Modified: 9.7.96  Klaus Hartung

% Calls 	Toolbox
% no calls

if (nargin < 1) error('No time series defined'); end;

[m,n] = size(x);

if (m==1) x = x'; [m,n] = size(x); end;

if (nargin<2) len = m; end;

X = fft(x,len);
NumUniquePts = ceil((len+1)/2);

X = X(1:NumUniquePts,:);

% correct amplitude

%X=X*2;

%X(1,:)=X(1,:)/2;

%if even(len)
% X(NumUniquePts,:)=X(NumUniquePts,:)/2;
%end

%X=X/len;


% end of file rfft.m 