function sigmoid_curve = sigmoid(x_axis,base,x_shift)
%sigmoid_curve = sigmoid(x_axis,base,x_shift)
%Function to produce a sigmoid, squashing-like, function
%Set to make "x_shift" the threshhold of the function (y ~ 0.05)
%x_axis: the values of the independent variable
%base:	the value of the base used in the exponential form of the sigmoid
%x_shift:the distance that the midpoint should be shifted from 0 (+/- in x_axis)

thresh_value = 0.1;
sigmoid_curve = 1 - 1./((1+base.^(x_axis - x_shift)));

threshhold = min(find(sigmoid_curve > thresh_value));
halfpoint = max(find(sigmoid_curve <= 0.5));

shift_more = x_axis(halfpoint) - x_axis(threshhold);

sigmoid_curve = 1 - 1./((1+base.^(x_axis - x_shift - shift_more)));

%plot(x_axis,sigmoid_curve)

return