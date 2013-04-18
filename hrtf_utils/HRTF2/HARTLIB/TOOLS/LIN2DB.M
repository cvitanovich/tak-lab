function y=lin2dB(x)

%Klaus Hartung
%14.12.96

y=zeros(size(x));

index=find(x>0);
y(index)=20.*log10(x(index));

return
