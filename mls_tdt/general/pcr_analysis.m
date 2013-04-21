function [PCR_meansurf,...
      PCR_diamond,...
      final_coef,...
      freq_bands,...
      pervar,...
      bbif_meansurf,...
      actd_curve_all,...
      new_bbif_meansurf,...
      new_ildaxis,...
      new_bbif_meansurf2,...
      new_ildaxis2,...
      pc,...
      latent,...
      pc_bbif_meansurf,...
      pc_totvar,...
      beta,...
   X,Y,IA_bestloc_ind,ILD_matrix_focus] = pcr_analysis(HRTFinfo,...
   ILD_matrix,...
   IA_meansurf,...
   IA_locs,...
   IA_bestloc_ind,...
   bandlimits,...
   bandwidth,...
   frac_points,...
   locreg_alpha,...
   locreg_mode,...
   mr_alpha,...
   stepwise_flag)   
%Function to carry out Principal Component Regression Analysis on ILD-Alone Data in order to
%predict the broadband ILD-freq response surface and to predict broadband ILD-Alone responses
%Variables on input:
%HRTFinfo: 		a structure with the HRTF header info, location_matrix and hrtf_freqs
%ILD_matrix:		the ILD matrix, with locations (columns) corresponding to HRTFinfo.location_matrix
%IA_meansurf: 	the mean vector of ILD-alone responses
%IA_locs: 		locations corresponding to IA_meansurf, in that order
%IA_bestloc_ind:	in case there are multiple top spots, you can specify the one of your choosing		
%bandlimits:		limits in Hz of broadband analysis
%bandwidth:		width (in octaves) of frequency band analysis for ILD-distance
%frac_points:		fraction of the total IA_meansurf points to include for parameterization
%locreg_alpha:	alpha for local regression (0 to 1)
%locreg_mode:		mode of local regression, 0 (default), 1 or 2 order (see locregress1)
%mr_alpha:		alpha for Multiple Regression (0.05 or smaller)
%stepwise_flag:	0 for automated MR, 1 for stepwise MR

new_ildaxis = [-30:1:30];

if(size(IA_locs,1) < 3)
   IA_locs = IA_locs';
end

%Get those locations in the ILDmatrix that match the measured ILDAlone locations
for num_loc = 1:size(IA_locs,1)
   ILDmat_index(num_loc) =...
      max(find(HRTFinfo.location_matrix(1,:) == IA_locs(num_loc,1) &...
      HRTFinfo.location_matrix(2,:) == IA_locs(num_loc,2)));
end
ILD_matrix_focus = ILD_matrix(:,ILDmat_index);


%Get the optimal ILD-alone location
if(isempty(IA_bestloc_ind))
   IA_bestloc_ind = max(find(IA_meansurf == max(IA_meansurf)));
end


%*****Begin Analysis for Octave bands

%Determine Frequency Bands to use in the MR analysis
clear freq_bands
freq_bands = get_freqbands(bandlimits,bandwidth);
bb_ILDf_freqaxis = freq_bands(:,1);


%Calculate distance_ILD matrix, unsigned and signed
disp('Calculating signed distance matrix...')
sign_flag = 1;
dist_ild_signed = get_ilddist(ILD_matrix_focus,HRTFinfo.hrtf_freqs,...
   IA_bestloc_ind,freq_bands,sign_flag);
dist_ild = dist_ild_signed;
[dist_ild_s_sort,dsort_s_index] = sort(dist_ild_signed,2);
%***Important: reverse the sorting for the comparison
[junk,reverse_sort_index] = sort(dsort_s_index,2);
clear junk

%Get the indices for inclusion in the Activity vs. d_ILD calculation
randtemp = randperm(size(ILD_matrix_focus,2));
keep_inds = randtemp(1:round(frac_points*size(ILD_matrix_focus,2)));
ind_num_in = 0;
ind_num_out = 0;
for loc = 1:size(IA_meansurf,2)
   if(~all(loc - keep_inds))
      ind_num_in = ind_num_in + 1;
      train_inds(ind_num_in) = loc;
   else
      ind_num_out = ind_num_out + 1;
      test_inds(ind_num_out) = loc;
   end
end
if(frac_points == 1) test_inds = train_inds; end

%Calculate Activity vs. ILD_distance curves
for freq_num = 1:size(freq_bands,1)
   %Determine how to sort the training indices, and how to reverse the sort
   [y,key_index] = sort(dist_ild_signed(freq_num,train_inds));
   [y,rev_key_index] = sort(key_index);
   %Find max and min d_ILDs, and use them to get the curve
   [y,min_ind] = min(dist_ild_signed(freq_num,:));
   [y,max_ind] = max(dist_ild_signed(freq_num,:));
   
   if(all(train_inds - min_ind) & all(train_inds - max_ind)) %max and min are NOT in the training set
      [actd_curve,resia] = ...
         locregress1([dist_ild_signed(freq_num,min_ind) ...
            sort(dist_ild_signed(freq_num,train_inds)) ...
            dist_ild_signed(freq_num,max_ind)],...
         IA_meansurf([min_ind train_inds(key_index) max_ind]),locreg_alpha);
      temp_curve = actd_curve(2:length(actd_curve)-1);
      d_ild_all{freq_num} = [dist_ild_signed(freq_num,min_ind) ...
            sort(dist_ild_signed(freq_num,train_inds)) ...
            dist_ild_signed(freq_num,max_ind)];
      actd_curve_all{freq_num} = actd_curve;
      actd_curve_use(freq_num,:) = temp_curve(rev_key_index);
   elseif(all(train_inds - min_ind)) %min is NOT in the training set
      [actd_curve,resia] = ...
         locregress1([dist_ild_signed(freq_num,min_ind) ...
            sort(dist_ild_signed(freq_num,train_inds))],...
         IA_meansurf([min_ind train_inds(key_index)]),locreg_alpha);
      temp_curve = actd_curve(2:length(actd_curve));
      d_ild_all{freq_num} = [dist_ild_signed(freq_num,min_ind) ...
            sort(dist_ild_signed(freq_num,train_inds))];
      actd_curve_all{freq_num} = actd_curve;
      actd_curve_use(freq_num,:) = temp_curve(rev_key_index);
   elseif(all(train_inds - max_ind)) %max is NOT in the training set
      [actd_curve,resia] = ...
         locregress1([sort(dist_ild_signed(freq_num,train_inds)) ...
            dist_ild_signed(freq_num,max_ind)],...
         IA_meansurf([train_inds(key_index) max_ind]),locreg_alpha);
      temp_curve = actd_curve(1:length(actd_curve)-1);
      d_ild_all{freq_num} = [sort(dist_ild_signed(freq_num,train_inds)) ...
            dist_ild_signed(freq_num,max_ind)];
      actd_curve_all{freq_num} = actd_curve;
      actd_curve_use(freq_num,:) = temp_curve(rev_key_index);
   else %both max and min ARE in the training set
      [actd_curve,resia] = ...
         locregress1(sort(dist_ild_signed(freq_num,train_inds)),...
         IA_meansurf(train_inds(key_index)),locreg_alpha);
      d_ild_all{freq_num} = sort(dist_ild_signed(freq_num,train_inds));
      actd_curve_all{freq_num} = actd_curve;
      actd_curve_use(freq_num,:) = actd_curve(rev_key_index);
   end
   
   disp(['Finished bb simulation for ' num2str(freq_bands(freq_num,1)) '-' num2str(freq_bands(freq_num,2)) ' Hz'])
end%end loop over frequencies

%Putative broadband dILD-freq surface
bbif_meansurf = actd_curve_use;

%Perform PCR (Principal Components Regression) on bbif_meansurf
[pc,score,latent,tsq] = princomp(bbif_meansurf');
frac_latent = latent/sum(latent);
temp = find(frac_latent >= 0.01) %Keep all principal components contributing >= threshhold variance
npc = length(temp) %number of PCs kept
%npc = 1;
pc_totvar = sum(frac_latent(temp));

%Transform bbif_meansurf by its PCs
for pc_num = 1:npc
   pc_bbif_meansurf(:,pc_num) = (sum(repmat(pc(:,pc_num),1,size(bbif_meansurf,2)) .* bbif_meansurf,1))';
end
      
%Center and scale the Y and X input in preparation for Multiple Regression
Y = IA_meansurf; X = pc_bbif_meansurf(:,1:npc);
syy = sum((Y - mean(Y)).^2);
Ycs = (Y - mean(Y))/syy;
for xnum = 1:size(X,2)
   sxx(xnum) = sum((X(:,xnum) - mean(X(:,xnum))).^2);
   X(:,xnum) = (X(:,xnum) - mean(X(:,xnum)))./sxx(xnum);
end

%Perform multiple regression using the PCs
if(stepwise_flag)
   stepwise(X,Y(train_inds)',1:npc,mr_alpha);
   keyboard
else
   [b1, delta1, tstats1, stats1] = ...
      regress_mls(X,Y(train_inds)',[1:npc],mr_alpha);        
   inmodel = [];
   bint = [b1-delta1 b1+delta1];
   for b_num = 1:length(b1)
      if(sign(bint(b_num,1)) == sign(bint(b_num,2))) inmodel(b_num) = b_num; end
   end
   inmodel = dezero(inmodel);
end
[b2, delta, tstats, stats] = regress_mls(X,Y(train_inds)',inmodel,mr_alpha);
beta = zeros(size(b2));
beta(inmodel) = b2(inmodel);

final_coef = pc(:,1:npc) * beta; %scale the PC's with the regression coefficients
final_coef = final_coef/max(final_coef);

%Test agreement on remaining locations
for freq_num = 1:size(freq_bands,1)
   activities(freq_num,:) = interp1(d_ild_all{freq_num},actd_curve_all{freq_num},dist_ild_signed(freq_num,:));
end

PCR_meansurf = final_coef' * activities;
temp = corrcoef(IA_meansurf(test_inds),PCR_meansurf(test_inds));
pervar = temp(1,2)^2;
%Make the diamond
[PCR_azi,PCR_ele,temp] = ...
   array2diamond(PCR_meansurf',IA_locs');

%Interpolate the missing values in the ILDAlone measurement
[AZ EL] = meshgrid(PCR_azi,PCR_ele);
% generate mask for missing points
missmask = NaN*ones(size(temp));
i = 1;
for az = -90:5:90;
   for el = -90+abs(az):5:90-abs(az)
      if (~(IA_locs(1,:)==el & IA_locs(2,:)==az))
         missmask(AZ==az & EL==el) = 1;
      else
         missmask(AZ==az & EL==el) = 0;
      end; 
   end;
end;
% replace all valid locations presently containing NaN's with zeros
ri = (missmask==1);
temp(ri) = zeros(size(temp(ri)));
% generate interpolation function - uses surrounding 4 squares to generate missing square
intfun = [[0 1 0]; [1 0 1]; [0 1 0]];
temp(find(isnan(temp))) = zeros(size(find(isnan(temp))));
intval_old = conv2(temp, intfun,'same')/4;
intval = intval_old.*missmask;
PCR_diamond = intval + temp;

      
%Generate the predicted Broadband ILD-Freq Response Surface
clear bbif_meansurf;
bbif_meansurf = activities;
%Interpolate to fit Activity/ILD curves onto newildaxis
if(size(bb_ILDf_freqaxis,2) > 1) bb_ILDf_freqaxis = bb_ILDf_freqaxis'; end
new_bbif_meansurf = griddata(repmat(bb_ILDf_freqaxis,1,size(bbif_meansurf,2)),...
   dist_ild_signed,...
   bbif_meansurf,...
   bb_ILDf_freqaxis,...
   new_ildaxis,'nearest');

new_bbif_meansurf = (new_bbif_meansurf)';

%Make sure we're not interpolating out of bounds
for freq_num = 1:length(bb_ILDf_freqaxis)
   tempind1 = find(new_ildaxis < min(dist_ild_signed(freq_num,:)));
   tempind2 = find(new_ildaxis > max(dist_ild_signed(freq_num,:)));
   new_bbif_meansurf(freq_num,[tempind1 tempind2]) = NaN;
end

%Apply frequency weights
for freq_num = 1:length(bb_ILDf_freqaxis)
   new_bbif_meansurf(freq_num,:) = ...
      new_bbif_meansurf(freq_num,:) * final_coef(freq_num);
end

%Transform by correcting for signed RMS ILD at optimal location
[new_bbif_meansurf2,new_ildaxis2] = ...
   iftransform_regress(new_bbif_meansurf,...
   bb_ILDf_freqaxis,...
   new_ildaxis,...
   ILD_matrix_focus(:,IA_bestloc_ind),...
   HRTFinfo);

return