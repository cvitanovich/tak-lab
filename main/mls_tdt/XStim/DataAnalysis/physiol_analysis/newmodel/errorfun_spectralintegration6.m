function err = errorfun_spectralintegration6(V,...
   IA_fullspec,IA_piece1,IA_piece2,IA_piece3,IA_piece4,index1,index2,index3,index4)
%This function takes only ILDAlone mean arrays
%errorfun_spectralintegration6(V,...
%   IA_fullspec,IA_piece1,IA_piece2,IA_piece3,IA_piece4,index1,index2,index3,index4)
%V:               variables to fit
%IA_fullspec:     measured ILDAlone RS with broadband noise
%IA_piece1:       measured ILDAlone RS, lo freq, left side
%IA_piece2:       measured ILDAlone RS, lo freq, right side
%IA_piece3:       measured ILDAlone RS, hi freq, left side
%IA_piece4:       measured ILDAlone RS, hi freq, right side

lofreq = zeros(size(IA_fullspec,1),size(IA_fullspec,2));
hifreq = zeros(size(IA_fullspec,1),size(IA_fullspec,2));

lofreq(index1) = V(1) * IA_piece1;
lofreq(index2) = V(2) * IA_piece2;
hifreq(index3) = V(3) * IA_piece3;
hifreq(index4) = V(4) * IA_piece4;

result = (lofreq + hifreq);
result = V(5) .* logsig(result);

err = sum( (result - IA_fullspec).^2);

return