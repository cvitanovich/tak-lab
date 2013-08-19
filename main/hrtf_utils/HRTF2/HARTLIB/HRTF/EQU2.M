function EQ_IR = equ2(IR, IR2, Wn, Fs)

% EQU2  Equalization of  HRTF with inverse TF (e.g. headphone).
%	eq_ir=EQU2(ir, ir2, Wn, Fs) 
%	eq_ir=EQU2(ir, ir2, Wn)    
%	convolves the impulse responses in ir with 
%	the inverse impulse responses of ir2. ir2 is 
%	a matrix with a 2-channel IR. Map is
%	a matrix with the corresponding directions
%	of ir. Wn is a vector with the corner frequencies
%	of the target bandpass. Fs is the sampling 
%	frequency (default 44100).

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektortechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%
% Date: 28.06.96  
% Modified: 6.7.96 Klaus Hartung

% Calls		toolbox
% lsqinv3	hartung/tools


if (nargin < 4)
	Fs = 44100;
	if (nargin < 3)	
		error ('Insufficient number of arguments');
	end;
end;

disp('Calculating inverse Filter');

INV_IR2 = lsqinv6(IR2/max(max(IR2)),Wn,Fs,0);
disp('ready');

[m,n_channels]= size(IR);
[m2,n2_channels]=size(INV_IR2);


for nc=1:1:n_channels/2
	clc
	disp('Calculating equalized IRs');
	disp(['channel: ' num2str(nc) '/' num2str(n_channels/2)]);
	disp('   ');
	EQ_IR(:,nc*2-1) = conv(IR(:,nc*2-1),INV_IR2(:,1));
	EQ_IR(:,nc*2)   = conv(IR(:,nc*2),  INV_IR2(:,2));
end
%irplot(EQ_IR,44100);

disp('EQU ready.')
return

% end of equ.m 
