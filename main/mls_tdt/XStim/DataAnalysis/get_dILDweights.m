function	[final_coef, pervar] = get_dILDweights(dILD, IAnSpikes, cF, dir, flags);

%function	[final_coef, pervar] = get_dILDweights(dILD, IAnSpikes, cF, dir, flags);
%
% requires argins:
%		dILD - ILD distance organized as [n_cF x nLocs]
%		IAnSpikes - mean number of spikes for each location using ILD-alone HRTFs
%		cF - center frequencies of bands
% optionally:	dir
%					flags = [pflag]
% returns the weights for each frequency band and the percent variance explained
%
% calculate the curves for Activity vs. d_ILD,
% find the principal frequency components of those curves
% transform the Activity/d_ILD curves by the principal components to get 1 estimate of the ILDAlone
%			response surface for each principal component
% perform multiple regression using those estimates as X's and the measured ILDAlone response surface as the Y
% to find the frequency weights, calculate a weighted sum of principal
%			components (each component weighted by its regression coefficient,
%			then use those frequency weights to calculate an estimate of the ILDAlone
%			response surface from the Activity vs. d_ILD curves
%
% also included: regress_mls, to perform regression using a specified subset of X's (from MATLAB's stepwise).

n_cF = length(cF);
use_stepwise = 0;

if nargin < 5
   pflag = 0;
else
   pflag = flags(1);
end

%Calculate Activity vs. ILD_distance curves
      alpha = 0.15;
      thresh = 0.3;
      for ib = 1:n_cF
			[sort_dILD, ind_dILD] = sort(dILD(ib,:));			% sort ILD distance
		 	[temp,reverse_ind_dILD] = sort(ind_dILD);     % get reverse sort index       
  
       [tempia,resia] = locregress1(sort_dILD,IAnSpikes(ind_dILD),alpha,0);		% local regression of dILD on Spikes
       tempia_allfreqs(ib,:) = tempia;
       resia_allfreqs(ib,:) = resia;
       bbif_meansurf(ib,:) = tempia(reverse_ind_dILD);
       %disp(['Finished bb simulation for ' num2str(cF(ib)) ' Hz'])
      end
      
%Perform PCR (Principal Components Regression) on bbif_meansurf
      [pc,score,latent,tsq] = princomp(bbif_meansurf');
      frac_latent = latent/sum(latent);
      temp = find(frac_latent >= 0.05);	 %Keep all principal components contributing >= threshhold fractional variance
      npc = length(temp);						 %number of PCs kept
      pc_totvar = sum(frac_latent(temp));
      disp([num2str(frac_latent(temp))])
      disp(['percent variance explained by ' num2str(npc) ' PCs: ' num2str(pc_totvar)])
      
%Transform bbif_meansurf by its PCs
      for pc_num = 1:npc
         pc_bbif_meansurf(:,pc_num) = (sum(repmat(pc(:,pc_num),1,size(bbif_meansurf,2)) .* bbif_meansurf,1))';
      end
      
%Perform multiple regression using the PCs
alpha = 0.05;
      if use_stepwise
         stepwise(pc_bbif_meansurf,IAnSpikes',1:size(pc_bbif_meansurf,2),alpha);
         keyboard
      else
         [b1, delta1, tstats1, stats1] = regress_mls(pc_bbif_meansurf,IAnSpikes',[1:size(pc_bbif_meansurf,2)],alpha);        
         inmodel = [];
         bint = [b1-delta1 b1+delta1];
         for b_num = 1:length(b1)
            if(sign(bint(b_num,1)) == sign(bint(b_num,2))) inmodel(b_num) = b_num; end
         end
         ind = find(inmodel ~= 0);
         inmodel = inmodel(ind);
      end
      [b2, delta, tstats, stats] = regress_mls(pc_bbif_meansurf,IAnSpikes',inmodel,alpha);
      beta = zeros(size(b2));
      beta(inmodel) = b2(inmodel);
      
      final_coef = pc(:,1:npc) * beta; 		%scale the PC's with the regression coefficients
      
      Test_meanarray = final_coef' * bbif_meansurf;
      temp = corrcoef(IAnSpikes,Test_meanarray);
      pervar = temp(1,2)^2;
      disp(['percent variance explained (retrospective): ' num2str(pervar)])
      
      temp = pc_bbif_meansurf(:,1) * frac_latent(1);
      for i = 2:npc
         temp = temp + pc_bbif_meansurf(:,i) * frac_latent(i);
      end
      
if pflag
   plot_diam(IAnSpikes,dir,1);
   title('ILD alone')
   plot_diam(temp,dir,1);
   title(['reconstructed ILD alone  : ' num2str(pervar) ' of variance explained'])
end
      
      
%%%%%%%%%%%%%%%%%%%%%%%%%      
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
