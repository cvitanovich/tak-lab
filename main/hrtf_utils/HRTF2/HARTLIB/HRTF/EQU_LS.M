function EQ_IR = equ_ls(IR, LSIR, MAP, Wn, Fs)

% EQU_LS equalize HRTF with inverse loudspeaker TF.
%
%	eq_ir=EQU_LS(ir, lsir, map, Wn, Fs) 
%	eq_ir=EQU_LS(ir, lsir, map, Wn)    
%	convolves the impulse responses in ir with 
%	the inverse impulse response of lsir. map is
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

if (nargin < 5)
	Fs = 4400;
	if (nargin < 4)	
		error ('Insufficient number of arguments');
	end;
end;
 
disp('Calculating inverse Filter ...');


LSIR2=optwin(LSIR,88,128);
LSIR=LSIR2/max(max(LSIR));

INV_LSIR = lsqinv3(LSIR,256,Wn,Fs);

disp('...ready');

[m,n_channels]= size(IR);
[m2,n2_channels]=size(INV_LSIR);


for nc=1:1:n_channels/2
	clc
	disp('Calculating equalized IRs');
	disp(['channel: ' num2str(nc) '/' num2str(n_channels/2)]);
	disp(' ');
        elevation = MAP(1,nc);
	ls_number = (elevation+10)/10+1;
	EQ_IR(:,nc*2-1) = conv(IR(:,nc*2-1),INV_LSIR(:,ls_number));
	EQ_IR(:,nc*2)   = conv(IR(:,nc*2),  INV_LSIR(:,ls_number));
	
end

clc

disp('EQU_LS ready.')

% end of equ_ls.m