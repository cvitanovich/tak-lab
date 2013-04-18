function EQ = equ_ls(IR, LSIR, MAP, Wn,Fs)

% EQUALIZE equalize HRTF with inverse loudspeaker TF.
%
% function EQ = equalize(IR, LSIR, MAP,Wn,Fs)
%
%	eq_ir=EQUALIZE(ir, lsir, map, Wn, Fs) 
%	convolves the impulse responses in ir with 
%	the inverse impulse response of lsir. map is
%	a matrix with the corresponding directions
%	of ir. Wn is a vector with the corner frequencies
%	of the target bandpass. Fs is the sampling 
%	frequency (default 44100).

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%

disp('Calculating inverse Filter');

INV_LSIR = lsqinv3(LSIR,512,Wn,Fs);

[m,n_channels]= size(IR);
[m2,n2_channels]=size(INV_LSIR);

disp('Calculating equalized IRs');

for nc=1:1:n_channels/2

        	elevation = MAP(1,nc);
		ls_number = (elevation+10)/10+1;
		EQ(:,nc*2-1) = conv(IR(:,nc*2-1),INV_LSIR(:,ls_number));
		EQ(:,nc*2)   = conv(IR(:,nc*2),  INV_LSIR(:,ls_number));
	
end


return
