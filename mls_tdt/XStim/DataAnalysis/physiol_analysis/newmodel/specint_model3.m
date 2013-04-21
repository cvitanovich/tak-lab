clear; close all

load 'd:\mlspezio\matlab\save\Neuron'
%specify model
%func = 'result = (1/sum(V)) * (V(1) * IA_bandwidths(1,:)) + (V(2) * IA_bandwidths(2,:));';
[err,func] = errorfun_spectralintegration3;

count = 0;
for cell_num = 19:20
   count = count + 1;
   IA_fullspec = Neuron{cell_num}.ia_meansurf{1};
   for band = 1:2
      IA_bandwidths(band,:) = Neuron{cell_num}.ia_meansurf{2*band};
   end
   %minimize the model error
   Options = optimset('Display','iter','MaxIter',500,'TolFun',0.002,'MaxFunEvals',1000);
   V = [1 1]; %coefficients for model
   [V,fval,exitflag,output] = fminsearch('errorfun_spectralintegration3',V,Options,...
      IA_fullspec,IA_bandwidths);
   eval(func);
   temp = corrcoef(result,IA_fullspec);
   cc(count) = temp(1,2);
   vars(count,:) = V;
end

