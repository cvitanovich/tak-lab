function preplow()

% PREPLOW correct directions for  upside-down measurements

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 


%  28.6.96
load HEADL

MAP_L(1,:)=-MAP_L(1,:);

index = find_ele(MAP_L,10);
index=[index find_ele(MAP_L,0);

for i=1:length(index)
	[m,n]=size(MAP);
	[m2,n2]=size(HRIR_L);
	if (n2 ~= n) = error('Size of MAP and HRIR_L do not match'); return; end;
        new_index=find(1:n~=index);
	HRIR_L=HRIR_L(:,new_index);
	MAP_L=MAP_L(:,new_index);
end;
