function end_fr = abi_atten2(initialrs, dbatten, spont_rate);
% ABI_ATTEN2  attenutate signal based on average loudness
%            takes an input response surface (e.g. space test)
%            and a set of attenutations for each point (same size
%            matrix as response surface.  For each point in the resp
%            surface, estimates the degree of response change due to a 
%            abi change of given magnitude and returns a new response
%            surface which reflects these estimates.  dbatten is attenuation
%            relative to db used during initialrs test, which is taken as zero.
%            Derived from abi_atten.  This one uses a similar algorithm,
%            but the parameters were derived by comparing space tests done
%            with and without abi equalization.  Tuning curves are sigmoidal in ABI
%            and linear in terms of initial response strength.

if (nargin<3) spont_rate = 0; end;

%abax = -50:.1:10;
%abi_curve = sigmoid(abax,.30,-9);
%abi_curve(abi_curve>1) = ones(size(abi_curve(abi_curve>1)));

% unwrap the input matrices so we can opperate on them
[n m] = size(initialrs);
initialrs = initialrs(:);
dbatten = dbatten(:);

end_fr = initialrs.*sigmoid(dbatten,.3,-9);
end_fr_temp = end_fr;

% rewrap the input output argument
end_fr = ones(n,m);
end_fr(:) = end_fr_temp;
