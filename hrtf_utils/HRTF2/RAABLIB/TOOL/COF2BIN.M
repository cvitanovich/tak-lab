function cof2bin(coffile,az,el,binfileprefix)
% COF2BIN writes the coefficients of two FIR
%         filters to two binary files.
%         COF2BIN(coffile,az,el,binfileprefix)

eval(['load ' coffile]);

dirpos=find(dirs(:,1)==az & dirs(:,2)==el);
cofpos=2*(dirpos-1)+1;
cofleft=B(cofpos+0,:);
cofright=B(cofpos+1,:);

writebin([binfileprefix '.l'],cofleft,'float32','ieee-be');
writebin([binfileprefix '.r'],cofright,'float32','ieee-be');

disp([cofleft(:),cofright(:)]);
