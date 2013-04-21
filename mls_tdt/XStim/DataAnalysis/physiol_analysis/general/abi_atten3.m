function end_fr = abi_atten3(initialrs, dbatten);
% ABI_ATTEN3  attenutate signal based on average loudness
%             this one written for pure tone responses
%            takes an input response surface (e.g. space test)
%            and a set of attenutations for each point (same size
%            matrix as response surface.  For each point in the resp
%            surface, estimates the degree of response change due to a 
%            abi change of given magnitude and returns a new response
%            surface which reflects these estimates.  dbatten is attenuation
%            relative to db used during initialrs test, which is taken as zero.
%            Derived from abi_atten.  
%            This version is based on pure-tone ILD tuning curves
%            Tuning curves are sigmoidal in ABI
%            and linear in terms of initial response strength.
%            NOTE: you may have to add or subtract values in order to
%                  get your ABL values to hit the ABL response curve at the
%                  correct location.  Do a histogram of all ABL values in your
%                  HRTFs and compare with the ABL tuning curve.
% 						 use the following command to view the underlying ABL curve:
%                  plot(-30:30, abi_atten3(ones(61,1),-40:30))

%abax = -50:.1:10;
%abi_curve = sigmoid(abax,.30,-9);
%abi_curve(abi_curve>1) = ones(size(abi_curve(abi_curve>1)));

% unwrap the input matrices so we can opperate on them
[n m] = size(initialrs);
initialrs = initialrs(:);
dbatten = dbatten(:);

% initially the abi curve intercepts zero at about .6 and then we
% normalize it so that it crosses zero at the value 1 so that when we 
% multiply by the value from this curve, it will do nothing when the atten
% is zero and will decrease output proportionately as lower abi values are 
% specified

% I think the values slope: .18 and intercept: -2 best fit the data
% however, I have played with these values to affect training results
if 1
	abi_curve_vals = sigmoid(dbatten, .14, -15);
else
	abi_curve_vals = sigmoid(dbatten, .18, -15);
end;

end_fr = initialrs.*abi_curve_vals;
end_fr_temp = end_fr;

% rewrap the input output argument
end_fr = ones(n,m);
end_fr(:) = end_fr_temp;
