function addlow(file)

% ADDLOWW merge the lower hemisphere to an upper hemisphere catalogue
%         
%         addlow(file) adds the lower hemisphere catalogue ../hksu/hksu.cat
%         to <file>cat.mat. The result is stored in <file>Xcat.mat and
%         <file>Xcat.srx.
% 	

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%
% Date: 28.06.96  
% Modified: 9.7.96 Klaus Hartung

% Calls 	Toolbox
% mat2mtl	hartung/hrtf
% fexist	hartung/tools

if (nargin < 1) error('No filename defined.'); end;
if ~fexist('../hksu/hksucat.mat') error (' ../hksu/hksucat.mat does not not exist'); end;
 
eval(['load ../hksu/hksucat.mat']);

MAP_L(1,:)=-MAP(1,:);
MAP_L(2,:)=MAP(2,:);

index = find_ele(MAP_L,-10)
index=[index find_ele(MAP_L,0)]
index=[index find_ele(MAP_L,10)]



[m,n]=size(MAP)
[m2,n2]=size(HRIR)
if (n2/2 ~= n) error('Size of MAP and HRIR_L do not match'); return; end;


new_index=1:n;
for i=1:length(index)
        new_index=new_index(find(new_index~=index(i)));	
end;

for k=1:length(new_index)
	HRIR_L(:,2*k-1)=HRIR(:,new_index(k)*2-1);
	HRIR_L(:,2*k)  =HRIR(:,new_index(k)*2);
end

MAP_L=MAP_L(:,new_index);

eval(['load ' file 'cat'])
        
[m,N_CHANNELS_L]=size(HRIR_L);

MAP=[MAP MAP_L];
HRIR=[HRIR HRIR_L];
N_CHANNELS=N_CHANNELS+N_CHANNELS_L;
N_DIR = N_CHANNELS/2;

newfile=[file(1:4) 'Xcat'];
eval(['save ' newfile ' HRIR MAP N_DIR N_CHANNELS Fs N_LINES']);
mat2mtl(newfile);

end

% end of file addlow.m