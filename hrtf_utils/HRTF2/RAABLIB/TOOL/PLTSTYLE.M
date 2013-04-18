function [so,nc]=pltstyle(si);
% PLTSTYLE Make correct plot style string
%          so=pltstyle(si); 
%          See MAGPLOT.M for usage.

n=length(si);
so='';
nc=0;
i=1;
q=0;
for i=1:n
  if isletter(si(i))
        if q==0 q=1; so=[so si(i)];      
    elseif q==1 q=1; so=[so '  ' si(i)];
    elseif q==2 q=1; so=[so ' ' si(i)]; 
    elseif q==3 q=1; so=[so si(i)];      end;
    nc=nc+1;
  else
        if q==0 error('invalid linetype');
    elseif q==1 q=2; so=[so si(i)];
    elseif q==2 q=3; so=[so si(i)];
    elseif q==3 error('invalid linetype'); end;
  end;
end;
if q==1 so=[so '  ']; elseif q==2 so=[so ' ']; end;  




  
  