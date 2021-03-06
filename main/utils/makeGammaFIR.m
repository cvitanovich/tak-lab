function [coefs]=makeGammaFIR(Fs,cF,species)

% Change the EarQ and ERB if you wish to use a different ERB scale.
switch species
    case 'human'
        EarQ = 9.26449;				% Glasberg and Moore Parameters
        minBW = 24.7;
        order = 1;
        ERB = ((cF/EarQ).^order + minBW^order).^(1/order);
    case 'owl'                      % Koppl parameters
        EarQ = 1.9 * 0.074 * (cF .^ .504);
        ERB = cF ./ EarQ;
end

cF = cF(:);
T=1/Fs;
B=1.019*2*pi*ERB;

A0 = T;
A2 = 0;
B0 = 1;
B1 = -2*cos(2*cF*pi*T)./exp(B*T);
B2 = exp(-2*B*T);

A11 = -(2*T*cos(2*cF*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cF*pi*T)./ ...
    exp(B*T))/2;
A12 = -(2*T*cos(2*cF*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cF*pi*T)./ ...
    exp(B*T))/2;
A13 = -(2*T*cos(2*cF*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cF*pi*T)./ ...
    exp(B*T))/2;
A14 = -(2*T*cos(2*cF*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cF*pi*T)./ ...
    exp(B*T))/2;

gain = abs((-2*exp(4*i*cF*pi*T)*T + ...
    2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
    (cos(2*cF*pi*T) - sqrt(3 - 2^(3/2))* ...
    sin(2*cF*pi*T))) .* ...
    (-2*exp(4*i*cF*pi*T)*T + ...
    2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
    (cos(2*cF*pi*T) + sqrt(3 - 2^(3/2)) * ...
    sin(2*cF*pi*T))).* ...
    (-2*exp(4*i*cF*pi*T)*T + ...
    2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
    (cos(2*cF*pi*T) - ...
    sqrt(3 + 2^(3/2))*sin(2*cF*pi*T))) .* ...
    (-2*exp(4*i*cF*pi*T)*T + 2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
    (cos(2*cF*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cF*pi*T))) ./ ...
    (-2 ./ exp(2*B*T) - 2*exp(4*i*cF*pi*T) +  ...
    2*(1 + exp(4*i*cF*pi*T))./exp(B*T)).^4);

allfilts = ones(length(cF),1);
fcoefs = [A0*allfilts A11 A12 A13 A14 A2*allfilts B0*allfilts B1 B2 gain];

% 512 points:
impulse = [1 zeros(1, 511)];
y1 = filter([A0/gain A11/gain A2/gain], [B0 B1 B2], impulse);
y2 = filter([A0 A12 A2], [B0 B1 B2], y1);
y3 = filter([A0 A13 A2], [B0 B1 B2], y2);
response = filter([A0 A14 A2], [B0 B1 B2], y3);

% scale response for unity gain:
coefs = response./sum(abs(response));