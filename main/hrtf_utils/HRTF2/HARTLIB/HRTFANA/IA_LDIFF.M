function ild=ia_ldiff(H_l, H_r, frq, sampling_frq)


% ild=ia_ldiff(H_l, H_r, selected_frq, sampling_frq)

% Calculates the interaural level difference  of a critical band centered a frq
% level difference is positiv if left ear is leading, 
% negative for right ear leading
% frq might be a vector with several frequencies


sz_l=size(H_l);
number_of_lines=length(H_l);

if (sz_l(1) ~= 1) H_l=reshape(H_l,1,number_of_lines); end;
sz_r=size(H_r);
if (sz_r(1) ~= 1) H_r=reshape(H_r,1,number_of_lines); end;

ild=[];
%filter_order=(number_of_lines-1)/4;
filter_order=(number_of_lines-1);


len_frq=length(frq);

for i=1:len_frq
roex=earfilt2(frq(i), number_of_lines, filter_order,sampling_frq);

H_l=roex.*H_l;
H_r=roex.*H_r;

I_l=sum(abs(H_l).^2);
I_r=sum(abs(H_r).^2);

L_l=10*log10(I_l);
L_r=10*log10(I_r);

ild(i)=L_l-L_r;
end;

sz=size(ild);
if(sz(1) ~= 1) ild=reshape(ild,1,length(ild)); end;


end;
