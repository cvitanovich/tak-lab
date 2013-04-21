function [beta, delta, tstats, stats] = regress_mls(X,y,inmodel,alpha);
% Local function for fitting stepwise regression model.
% Taken from FitModel in Stepwise
[n, p] = size(X);
colidx = 1:p;
notin  = colidx;
notin(inmodel) = [];

if nargin == 4
   sigprob = 1 - alpha/2;
else
   sigprob = (0.975).^(1/p);
end

df = n-2:-1:n-p-1; %degrees of freedom
crit = tinv(sigprob,df);


% Do initial statistical calculations.
mx = mean(X);
X = (X - mx(ones(n,1),:));
y = y - mean(y);
sx = std(X);
scaledx = X./sx(ones(n,1),:);

if ~isempty(inmodel)
   [Q,R]=qr(X(:,inmodel),0);
   b = R\(Q'*y);
   r = y - X(:,inmodel)*b;
   RI = R\eye(size(R));
   xtxid = (sum(RI'.*RI'))';
else
   r = y;
end
terms = length(inmodel);
if n-terms-1 == 0
   mse = 0;
else
   mse = r'*r./(n-terms-1);
end
rmse = sqrt(mse);
tss = y'*y;
rss = y'*y-r'*r;
if terms == 0
   prob = 1;
   r2   = 0;
   f    = 0;
elseif mse == 0
   f = Inf;
   prob = 0;
   r2 = 1;
   seb = 0;
   t = Inf;
else
   f   = (rss./terms)./mse;
   prob = 1-fcdf(f,terms,(n-terms-1));
   r2  = rss./tss;
   seb = sqrt(xtxid)*rmse;
   t = b./seb;
end

beta = zeros(p,1);
tstats = zeros(p,1);
delta = zeros(p,1);
if ~isempty(inmodel)
   delta(inmodel) = seb*crit(terms);
   beta(inmodel) = b;
   tstats(inmodel) = t;
end


for index = notin
   Xnew = [X(:,inmodel) X(:,index)];
   [Qadd, Radd] = qr(Xnew,0);
   bnew = Radd\(Qadd'*y);
   rnew = y - Xnew*bnew;
   RInew = Radd\eye(terms+1,terms+1);
   xtxidnew = (sum(RInew'.*RInew'))';
   if (n == terms+2)
     rmse1 = Inf;
   else
     rmse1 = sqrt(sum(rnew.*rnew)./(n-terms-2));
   end
   seb = sqrt(xtxidnew(terms+1))*rmse1;
   
   beta(index) = bnew(terms+1);
   tstats(index) = bnew(terms+1)./seb;
   delta(index) = seb*crit(terms+1);
end
delta(isnan(delta)) = 0;
stats = [rmse r2 f prob];
