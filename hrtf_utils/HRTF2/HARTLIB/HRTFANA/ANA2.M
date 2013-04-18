function ana2(file)


% ANA2 analysis of HRTF
% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 


% Klaus Hartung
% Dec. 95

filename=[file 'mic'];
eval([' load ' filename]);
 
[m,n] = size(HRTF);
% [mi,ni] = size(HRIR);
mi=(m-1)*2;

b_array = [4:22];
disp(['ANA']);

disp(['Allocation of memory'])

CCSP=zeros(m,n/2)+i*zeros(m,n/2);
NBCCSP=zeros(size(CCSP));

ILD=zeros(length(b_array),n/2);
ITD=zeros(length(b_array),n/2);
LEVEL=zeros(length(b_array),n);

disp(['Calculation']);

disp(['Multiplication CCSP'])
CCSP= HRTF(:,1:2:n) .* conj(HRTF(:,2:2:n));
disp(['End of ']);


disp(['h_lp'])
h_lp=firls(256,[0 1000 1500 Fs]/Fs,[1 1 0 0]);
H_lp=rfft(h_lp,(m-1)*2);
disp(['End h_lp']);


for b=1:length(b_array)

	disp(['Calculating band ' num2str(b) ' ' num2str(b_array(b))]);
         
	f=bark2hz(b_array(b))

	h=auditory_filter(f,Fs);
        H=rfft(h,mi);
        
        disp(['time']);


        NBCCSP = mult(CCSP,abs(H).^2);

        if (f< 1600)
                 CCF=rifft(NBCCSP);
        else
                 CCF=real(ifft([4*NBCCSP; zeros(m-2,n/2)]));
        end;
        
        len=length(CCF);
        x=[0:1:len/2-1   -len/2:1:-1];

        [maximum,maxpos]=max(abs(CCF));

        ITD(b,:)=-x(maxpos)/Fs;   


        disp('level');

        NB = mult(HRTF,H);

    	I_l=sum(abs(NB(:,1:2:n)).^2);
	I_r=sum(abs(NB(:,2:2:n)).^2);

	L_l=10*log10(I_l);
	L_r=10*log10(I_r);

        LEVEL(b,1:2:244)=L_l;
	LEVEL(b,2:2:244)=L_r; 
	ILD(b,:)=L_l-L_r;
    

end;

clear I_L I_R L_l L_r NB maximum maxpos x len CCF NBCCSP CCSP H h f b h_lp H_lp  m n
outname=[file 'ana'];

eval(['save ' outname ' LEVEL ILD ITD MAP'])
end
