%mr_analysis2
%Script to carry out a Multiple Regression modelling and analysis of ILDAlone activity
%Dependent variable: ILDAlone activity
%Independent variables:	Tonal activity at constant frequency, all locations

clear; close all;

%set parameters
bird_number = 899;
begin_cells = 1;
tot_num_cells = 24;
end_cells = begin_cells + tot_num_cells - 1;

%Graphics parameters
bkgndcolor = [0.85 0.85 0.85];
fig_ht = 10;
fig_wd = 8;

%Printing
print_flag = 0;

%Directories & files
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\HRTF_files\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];
get_HRTF = 0;
get_ITDmatrix = 0;
load 'd:\mlspezio\matlab\save\899HRTF' ILD_matrix ABL_matrix HRTFinfo

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'
clear ITD_matrix ABI_matrix test_numbers reps

load d:\mlspezio\matlab\save\tiaregressvar tiaregressvar
for cell_num = [12:20 22:24]
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
      
      disp(['Processing cell # ' num2str(cell_num)])
      
      %Choose the elements for analysis:     
      ILDf_freqaxis  = Neuron{cell_num}.tif_freqaxis;
      ILDf_ildaxis   = Neuron{cell_num}.tif_ildaxis;
      ILDf_meansurf  = Neuron{cell_num}.tif_meansurf;
      ILDf_meanarray = Neuron{cell_num}.tia_meanarray;
      
      IA_meansurf = Neuron{cell_num}.ia_meansurf{1};
      IA_stdsurf  = Neuron{cell_num}.ia_stdsurf{1};
      IA_diamond  = Neuron{cell_num}.ia_diamond{1};
      IA_locs		= Neuron{cell_num}.ia_locs{1};
      IA_azi		= Neuron{cell_num}.ia_azi{1};
      IA_ele		= Neuron{cell_num}.ia_ele{1};
      
   %Get those locations in the ILDmatrix that match the measured ILDAlone locations
   for num_loc = 1:size(IA_locs,1)
      ILDmat_index(num_loc) =...
         max(find(HRTFinfo.location_matrix(1,:) == IA_locs(num_loc,1) &...
         HRTFinfo.location_matrix(2,:) == IA_locs(num_loc,2)));
   end
   ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
   
   %Get the minimum and maximum frequency indices
   [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - ILDf_freqaxis(1)));
   [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - ...
      ILDf_freqaxis(length(ILDf_freqaxis))));
   
   IA_bestloc_ind(cell_num) = max(find(IA_meansurf == max(IA_meansurf)));
   if(cell_num == 1) IA_bestloc_ind(cell_num) = 129; end
   if(cell_num == 4) IA_bestloc_ind(cell_num) = 198; end
   if(cell_num == 7) IA_bestloc_ind(cell_num) = 163; end
   if(cell_num == 9) IA_bestloc_ind(cell_num) = 282; end
   if(cell_num == 11) IA_bestloc_ind(cell_num) = 127; end
   if(cell_num == 13) IA_bestloc_ind(cell_num) = 130; end
   if(cell_num == 16) IA_bestloc_ind(cell_num) = 159; end
   if(cell_num == 18) IA_bestloc_ind(cell_num) = 78; end
   if(cell_num == 24) IA_bestloc_ind(cell_num) = 128; end
   
   %Calculate the simulated ILDAlone RS for all frequency bands
   disp('Simulating ILDAlone RSs for all frequency bands...')
   for freq_num = 1:length(ILDf_freqaxis)
      temp_ILDf_meansurf = zeros(size(ILDf_meansurf));
      temp_ILDf_meansurf(freq_num,:) = ILDf_meansurf(freq_num,:);
      [temp,dirs,full_rs,new_f_ax,new_ild_ax] = ...
         if2space(ILD_matrix_focus',...
         HRTFinfo.hrtf_freqs,...
         ILDmat_index,...
         temp_ILDf_meansurf,...
         ILDf_freqaxis,...
         ILDf_ildaxis);
      %This step makes sure that we don't have completely 0 rows
      if(~any(temp))
         temp(1) = 0.01;
      end
      tia_meansurf(freq_num,:) = temp;
      disp(['Finished simulation for ' num2str(ILDf_freqaxis(freq_num)) ' Hz'])
   end
   
   
   %Perform multiple regression
   alpha = 0.05; %confidence interval criterion
   stepwise(tia_meansurf',IA_meansurf',1:length(ILDf_freqaxis),alpha);
   keyboard
   
      %Only take those coefficients that contribute significantly
      newbeta = beta;
      for freq_num = 1:length(ILDf_freqaxis)
         if( (sign(betaci(freq_num,1)) ~= sign(betaci(freq_num,2))) | ~all(out - freq_num)) %out is a vector of freq_bands not included in the model
            newbeta(freq_num) = 0;
         end
      end
      
      %Calculate the IA RS using the regression analysis
      Regress_tia_meanarray = newbeta' * tia_meansurf;
      
      IA_locs = IA_locs';
      [Regress_tia_azi,Regress_tia_ele,temp] = ...
         array2diamond(Regress_tia_meanarray,IA_locs);
      
      %Interpolate the missing values in the ILDAlone measurement
      [AZ EL] = meshgrid(Regress_tia_azi,Regress_tia_ele);
      
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
      Regress_tia_diamond = intval + temp;
      
      [y,minfreqind] = min(abs(HRTFinfo.hrtf_freqs - min(ILDf_freqaxis)));
      [y,maxfreqind] = min(abs(HRTFinfo.hrtf_freqs - max(ILDf_freqaxis)));
      
      cctemp = corrcoef(IA_meansurf,ILDf_meanarray(ILDmat_index));
      temp = corrcoef(Regress_tia_meanarray,Neuron{cell_num}.ia_meansurf{1});
      tiaregressvar(cell_num) = temp(1,2)^2;
   
   end %end if

end %end for loop
