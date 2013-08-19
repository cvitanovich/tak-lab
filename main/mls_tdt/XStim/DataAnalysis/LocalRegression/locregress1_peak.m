function [result,residual,peak_pos] = locregress1_peak(X,Y,alpha,mode,peak_pos)
%Reference: Clive Loader, _Local Regression and Likelihood_(New York, NY: Springer, 1999)
%alpha:				0 -> 1, determines bandwidth
%mode:				0: local, 1: linear, 2: quadratic
%peak_pos:			specify a position in X that is to be weighted very heavily

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

%Downweight selected values if they don't match some threshhold relative to the Y value at the peak_pos
if(peak_pos > 0)
   extra_weights = ones(size(weights));
   nXpos = fix(0.30 * length(X));
   peakYval = Y(peak_pos);
   threshhold = 0.7;
   for X_num = 1:size(X_mat_sort,2)
      Xval_index(X_num) = find(sort_index(:,X_num) == peak_pos);
      Yval_index = ...
         find(Y_mat_sort(max(1,Xval_index(X_num)-nXpos):...
         min(Xval_index(X_num)+nXpos,length(X)),X_num) < (threshhold * peakYval));
      Yval_index2 = ...
         find(Y_mat_sort(max(1,Xval_index(X_num)-nXpos):...
         min(Xval_index(X_num)+nXpos,length(X)),X_num) >= (threshhold * peakYval));
      for Y_num = 1:length(Yval_index)
         first_ind = Xval_index(X_num)-nXpos + Yval_index(Y_num) - 1;
         if(first_ind <= 0) first_ind = 1; end
         extra_weights(first_ind,X_num) = 2000;
      end
      for Y_num = 1:length(Yval_index2)
         first_ind = Xval_index(X_num)-nXpos + Yval_index2(Y_num) - 1;
         if(first_ind <= 0) first_ind = 1; end
         extra_weights(first_ind,X_num) = 1;
      end
      clear Yval_index Yval_index2
   end
   weights = weights .* extra_weights;
end
   
neg_index = find(weights < 0);
weights(neg_index) = 0;
clear neg_index

%Calculate pieces of the local quadratic regression estimate

%1. Local Constant Regression Term (a0)
a0 = sum(Y_mat_sort .* weights,1) ./ sum(weights,1); %a0

%2. Local Linear Regression Term (a1 * (xi - x))
X_bar_w = sum(weights .* X_mat_sort,1) ./ sum(weights,1);
X_bar_w_mat = repmat(X_bar_w',1,size(X_mat_sort,2));
diff_w = X_mat_sort - X_bar_w_mat;
a1 = (sum(weights .* diff_w .* Y_mat_sort,1)./sum(weights .* (diff_w.^2),1));

%3. Local Quadratic Regression Term (1/2 * a2 * (xi - x)^2)
a2 = (sum(weights .* diff_w.^2 .* Y_mat_sort,1)./sum(weights .* (diff_w.^4),1));


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