function [rbest phi_best] = lonoisenoiseNorm(time,Fs,s)

matlabpool(8);

global minfreq
global maxfreq
global dur
global len
global range
global mag
global Fsamp

fftw('planner','exhaustive');
Fsamp = Fs;

prior = -inf;

minfreq = 2000;
maxfreq = 11000;

dur = time/1000;    %now dur is in sec
len = round(dur*Fsamp); %number of samples

minfreq = round(((minfreq + 1)/Fsamp) * len);
maxfreq = round(((maxfreq + 1)/Fsamp) * len);
range = maxfreq - minfreq + 1;
% mag spectrum = 1 between set frequencies:
mag = zeros(len,1);
mag(minfreq:maxfreq) = ones(range,1);

for k = 1:1
    phi = zeros(len,1);
    phi(minfreq:maxfreq) = (rand(s,range,1) - 0.5) * (2*pi);
    g0 = Grad2(phi);
    x = min1(phi,g0);
    phi = phi + x*g0;

    noise = real(ifft(mag .* ( (cos(phi)) + (sqrt(-1)* sin(phi)) )))';
    %noise = (noise/max(abs(noise)));
    %noise = (abs(hilbert(noise)));
    noise = noise - mean(noise);
    noise = noise./std(noise);
    %mu0 = mean(noise);
    %n0= length(noise);
    r = std(abs(hilbert(noise))); 
    %var(noise)/var(abs(hilbert(noise))); %(1/n0)*sum((noise - mu0).^4) / ((1/n0)*sum((noise - mu0).^2))^2; 

    rlast = r;
    g1 = Grad2(phi);
    flag = 1;
    
    while flag == 1
        tic
        noise = real(ifft(mag .* ( (cos(phi)) + (sqrt(-1)* sin(phi)) )))';
        %noise = (noise/max(abs(noise)));
        %noise = abs(hilbert(noise));
        %mu0 = mean(noise);
        %n0= length(noise);
        noise = noise - mean(noise);
        noise = noise./std(noise);
        old = std(abs(hilbert(noise))); 
        %var(noise)/var(abs(hilbert(noise))); %(1/n0)*sum((noise - mu0).^4) / ((1/n0)*sum((noise - mu0).^2))^2;
        
        [x y] = min2(phi,g0,g1);
        phi = phi + x*g0 + y*g1;
        g0 = g1;
        g1 = Grad2(phi);
        noise = real(ifft(mag .* ( (cos(phi)) + (sqrt(-1)* sin(phi)) )))';
        noise = noise - mean(noise);
        noise = noise./std(noise);

        %(noise/max(abs(noise))); 
        %         noise = abs(hilbert(noise));
        %mu0 = mean(noise);
        %n0 = length(noise);
        new = std(abs(hilbert(noise)))
	%var(noise)/var(abs(hilbert(noise))) %(1/n0)*sum((noise - mu0).^4) / ((1/n0)*sum((noise - mu0).^2))^2;
        toc
        if new <= (old - 0.00001)
        else
            flag = 0;
        end
    end

    if new > prior
        prior = new;
        phi_best = phi;
    end

end

rbest = prior;

    % Usage: g = Grad(fun, x0)
    %     fun: name of the mu0ltidimensional scalar function
    %          (string). This function takes a vector argument of
    %          length n and returns a scalar.
    %     x0: point of interest (vector of length n)
    %     g: column vector containing the gradient of fun at x0. The
    %        size(g) = size(x)

    function g = Grad2(phi)
    % |delta(i)| is relative to |x0(i)|
    
    len = length(phi);
    delta = 0.01;
    g = zeros(length(phi),1);

    
    parfor p=minfreq:maxfreq
        tmp = zeros(len,1);
        tmp(p) = delta;
        noi = zeros(length(phi),1);
        noi = real( ifft( mag .* ((cos(phi+tmp)) + (sqrt(-1)*sin(phi+tmp))) ));
        noi = noi - mean(noi);
        noi = noi./std(noi);
        %        noi = noi/(abs(max(noi)));
        %         noi = abs(hilbert(noi));
%         mu = mean(noi);
%         n = length(noi);
        f1 = std(abs(hilbert(noi))); %var(noi)/var(abs(hilbert(noi))); %f1 = (1/n)*sum((noi - mu).^4) / ((1/n)*sum((noi - mu).^2))^2;
        
        noi = real(ifft(mag .* ((cos(phi-tmp)) + (sqrt(-1)*sin(phi-tmp))) ));
        noi = noi - mean(noi);
        noi = noi./std(noi);
        %noi = noi/(abs(max(noi)));
%         noi = abs(hilbert(noi));
        %mu = mean(noi);
        %n = length(noi);
        f2 = std(abs(hilbert(noi))); %var(noi)/var(abs(hilbert(noi))); %f2 = (1/n)*sum((noi - mu).^4) / ((1/n)*sum((noi - mu).^2))^2;
        g(p) = (f1-f2)/(2*delta);
    end

    end

    function x = min1(phi, g0)

    % phi' = phi + x*gradient
    x=0;
    r = zeros(201,1);
    
    noi = zeros(length(phi),201);

    parfor i = 1:201
        noi = zeros(length(phi),1);
        noi = real(ifft(mag .* ( cos( phi + ((1/100)*i-(101/100)) .*g0 ) + (sqrt(-1) * sin( phi + ((1/100)*i-(101/100)) .*g0 ) ) ) ));
        %noi = noi/abs(max(noi));
        %noi = abs(hilbert(noi));
        %mu = mean(noi);
        %n = length(noi);
        %r(i) = var(noi)/var(abs(hilbert(noi)));
        noi = noi - mean(noi);
        noi = noi./std(noi);
        r(i) = std(abs(hilbert(noi)));
    end
        

    lowest = inf;
    for k2 = 1 : length(r)
        if r(k2) < lowest
            lowest = k2;
        end
    end
    
    x = (1/100)*lowest - (101/100);
    
    end
    
    function [xout yout] = min2(phi, g0, g1)

    %phi = phi + x*g0 + y*g1;
    
    %xy = randn(s,1,1000);

    %newvals = zeros(1000,1);
    lowest = inf;
    newvals = zeros(1000,1);
    x = zeros(1000,1);
    y = zeros(1000,1);
    parfor i2=1:1000
        x(i2) = randn(s,1,1);
        y(i2) = randn(s,1,1);
        noi = real(ifft(mag .* ( (cos(phi + x(i2)*g0 + y(i2)*g1)) + (sqrt(-1)* sin(phi + x(i2)*g0 + y(i2)*g1)) )));
        noi = noi - mean(noi);
        noi = noi./std(noi);
%        noi = noi/abs(max(noi));
%         noi = abs(hilbert(noi));
%        mu = mean(noi);
%        n = length(noi);
        newvals(i2) = std(abs(hilbert(noi))); %var(noi)/var(abs(hilbert(noi)));
%        newvals(i2,2) = x;
 %       newvals(i2,3) = y;
    end

    lowest = find(newvals == min(newvals));
    xout = x(lowest(1));
    yout = y(lowest(1));
    
    end
end 
