function [a2] = dezero(a1)
%Function dzero
%To remove zero values from an array, as in an array of REACTION TIMES
%[a2] = dezero(a1)

   count = 0;
	for n=1:length(a1)
   	if(a1(n) ~=0)
      	count = count+1;
      	temp(count) = a1(n);
   	end
	end
   
   if (~count)
      a2 = 0;
   else
      a2 = temp';
   end
   
return
