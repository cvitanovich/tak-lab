% MON_ANA ?

[m,n] = size(HRTF);
[mi,ni] = size(HRIR);

b_array = [4:22];
Level=zeros(length(b_array),n);

for b=1:length(b_array)

         
	f=bark2hz(b_array(b));

	disp(['Calculating band ' num2str(b) ', ' num2str(b_array(b)) ' Bark,  f = ' num2str(f) ' Hz' ]);

	h=auditory_filter(f);
        H=rfft(h,mi);
        
        NB = mult(HRTF,H);
    	Intensity=sum(abs(NB).^2);
	Level(b,:)=10*log10(Intensity);
    
end;

