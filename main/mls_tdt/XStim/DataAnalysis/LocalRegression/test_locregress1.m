clear
numpoints = 100;
alpha = .1;
mode = 2;
num_iterations = 20;
conv_value = 0.001;
X = 0:(2*pi)/(numpoints-1):2*pi;
Y = hist(randn(1,1000),numpoints);


[result1,residual] = locregress1(X,Y,alpha,mode);
[result2,residual,K,V,weights] = locrbregress1(X,Y,alpha,mode,num_iterations,conv_value);
plot(X,Y,'r.')
hold on
plot(X,result1,'m');
plot(X,result2,'g');
