function [result,residual] = locregress1(X,Y,alpha,mode)
%locregress1(X,Y,alpha,mode)
%Reference: Clive Loader, _Local Regression and Likelihood_(New York, NY: Springer, 1999)
%alpha:				0 -> 1, determines bandwidth
%mode:				0: local, 1: linear, 2: quadratic

%Check to see that inputs are correct
if( (size(X,1) > 1 & size(X,2) > 1) | (size(Y,1) > 1 & size(Y,2) > 1) )
   error('Requires 2 vectors of the same length.')
end
if(length(X) ~= length(Y))
   error('Requires 2 vectors of the same length.')
end
if(nargin < 5)
   peak_pos = 0;
   if(nargin < 4)
      mode = 1;
      if(nargin < 3)
         alpha = 0.2;
      end
   end
end
if(size(X,2) > 1)
   X = X';
end
if(size(Y,2) > 1)
   Y = Y';
end

X_mat = repmat(X,1,size(X,1));
Y_mat = repmat(Y,1,size(Y,1));

%Get the kth distance to use as h (h is the "bandwidth"; it determines how smooth the curve is)
%This implements the "nearest neighbor" method for estimating h
kth_distance = fix(alpha * length(X));

%Calculate all distances, then sort each row from smallest to largest
dist_X = (abs(X_mat - X_mat'))';
diff_X = (X_mat - X_mat')';
[dist_X_sort,sort_index] = sort(dist_X,1);
X_mat_sort = X_mat(sort_index);
Y_mat_sort = Y_mat(sort_index);
diff_X_sort = diff_X(sort_index);

%Calculate weights, using h as specified by the kth distance from each X
h_mat = repmat(dist_X_sort(kth_distance,:),size(X,1),1);
U = dist_X_sort ./ h_mat';
weights = (1 - U.^3).^3; %tricubic method

neg_index = find(weights < 0);
weights(neg_index) = 0;
clear neg_index

%Calculate pieces of the local quadratic regression estimate

%1. Local Constant Regression Term (a0)
a0 = sum(Y_mat_sort .* weights,1) ./ sum(weights,1); %a0
nan_index = find(isnan(a0));
a0(nan_index) = 0;

%2. Local Linear Regression Term (a1 * (xi - x))
X_bar_w = sum(weights .* X_mat_sort,1) ./ sum(weights,1);
X_bar_w_mat = repmat(X_bar_w',1,size(X_mat_sort,2));
diff_w = X_mat_sort - X_bar_w_mat;
a1 = (sum(weights .* diff_w .* Y_mat_sort,1)./sum(weights .* (diff_w.^2),1));
nan_index = find(isnan(a1));
a1(nan_index) = 0;


%3. Local Quadratic Regression Term (1/2 * a2 * (xi - x)^2)
a2 = (sum(weights .* diff_w.^2 .* Y_mat_sort,1)./sum(weights .* (diff_w.^4),1));
nan_index = find(isnan(a2));
a2(nan_index) = 0;



%Put it all together to get the result
switch mode
case 0,
   result = a0;
case 1,
   result = a0 + (a1 .* (diff_w(:,1))');
case 2,
   result = a0 + (a1 .* (diff_w(:,1))') + (a2 .* ((diff_w(:,1))').^2);
end

residual = Y' - result;

return