function itd=ia_tdiff(H_l, H_r, frq, sampling_frq)

% itd=ia_tdiff(H_l, H_r, selected_frq, sampling_frq)

% Calculates the interaural time delay of a critical band centered a frq
% time delay is positiv if left ear is leading, negative for right ear leading
% frq might be a vector with several frequencies

sz_l=size(H_l);
number_of_lines=length(H_l);

if (sz_l(1) ~= 1) H_l=reshape(H_l,1,number_of_lines); end;
sz_r=size(H_r);
if (sz_r(1) ~= 1) H_r=reshape(H_r,1,number_of_lines); end;

itd=[];

%filter_order=(number_of_lines-1)/4;
filter_order=(number_of_lines-1);
len_frq=length(frq);

ccsp = (H_r.*conj(H_l));


for i=1:len_frq

	roex=earfilt2(frq(i), number_of_lines, filter_order,sampling_frq);

	ccsp = ccsp .* abs(roex).^2;

	if (frq < 1600)
		 ccf=rifft(ccsp);
	else
		 ccf=abs(ifft([4*ccsp, zeros(1,number_of_lines-2)]));
	end;

	ccf=real(ccf);

	s_ccf=[];

	len=length(ccf);
        x=[0:1:len/2-1 -len/2:1:-1];

%	s_ccf(1:len/2)=ccf(len/2+1:len);
%	s_ccf(len/2+1:len)=ccf(1:len/2);

	[maximum,maxpos]=max(abs(ccf));

%	[maximum,maxpos]=max(abs(s_ccf));

%	maxpos=maxpos-length(s_ccf)/2-1;

%	itd(i)=maxpos/sampling_frq;
        itd(i)=x(maxpos)/sampling_frq;   
end;

sz=size(itd);
if(sz(1) ~= 1) itd=reshape(itd,1,length(itd)); end;

end;