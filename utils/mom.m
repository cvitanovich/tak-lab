function [moments] = mom(input,n)

%function [moments] = mom(input,n)
% argin 'n' specifies which single moment to return
% default is to return all four
% calculates moments to fourth order:
%		all normalized by appropriate RMS level
%		mean (DC POWER)
%			for normal curve == 0
%		std (AC POWER)
%			for normal curve == 1
%		AC skewness
%			negative value -> skewness to left of mean
%			positive value -> skewness to right of mean
%		AC kurtosis (kurtosis = 4th moment - 3)
%			negative kurtosis -> platykurtosis (more intermediate values)
%			positive kurtosis -> leptokurtosis (more values near mean & tails)
%		Note: Pumplin uses Z = sqrt(mom(4) - 1) to compare
%				modulations


A = mean(input);
B = sqrt(mean((input - A).^2));				% RMS
if B > 0
   C = mean((input - A).^3) / (B .^3);
   D = mean((input - A).^4) / (B .^4);
else
   C = 0;
   D = 0;
end

moments = [A B C D];
if exist('n') & (n==1 | n==2 | n==3 | n==4)
   moments = moments(n);
end
