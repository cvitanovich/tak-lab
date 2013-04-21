function [net,New_BP_surface] = ...
   nnmodel_1(bird_number,side_of_brain,test_numbers,hrtf_file,get_hrtf)
%Function to model the transformation of Tonal to 1/3 BP ILD/Freq surface 
%This is a one-layer network that has the following architecture:
%        Input --> Layer 1 (Linear) --> Output
%Input is the Tonal ILD/Freq RS, Output is the BP ILD/Freq RS (modeled)

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];

%1. Get information from bird's HRTF catalogue - for use in predicting ILD Alone Surfaces
   %Header information
   [filetype,...
         info_blocks,...
         n_channels,n_lines, sampling_rate,...
         first_line,last_line,num_locations,...
         comment1, comment2] = mtlrh([hrtf_directory hrtf_file]);
   %Location information
   temp = mtlrdir([hrtf_directory hrtf_file]);
   location_matrix = sph2dbl(temp);
   clear temp
   
   hrtf_freqs = [first_line:(last_line-first_line)/(n_lines-1):last_line]';
   
if(get_hrtf == 1)
   %Transfer Functions
   for loc = 1:num_locations
      temp = mtlrch([hrtf_directory hrtf_file],2*loc-1);
      left_raw(:,loc) = temp;
      left = 20*log10(abs(temp));
      clear temp
      temp = mtlrch([hrtf_directory hrtf_file],2*loc);
      right_raw(:,loc) = temp;
      right = 20*log10(abs(temp));
      clear temp	
      ILD_matrix(:,loc) = right - left;
      ABI_matrix(:,loc) = right + left ./2;
         
      if(mod(loc,10) == 0)
         disp(['Finished location ' num2str(loc)])
      end
   end
   eval(['save ' hrtf_directory num2str(bird_number) 'ildmatrix ILD_matrix left_raw right_raw']);
   eval(['save ' hrtf_directory num2str(bird_number) 'abimatrix ABI_matrix']);
else
   eval(['load ' hrtf_directory num2str(bird_number) 'ildmatrix.mat']);
   eval(['load ' hrtf_directory num2str(bird_number) 'abimatrix.mat']);
end

% eliminate location of AZ=0, EL=90 so that it won't give us problems
oklocs = location_matrix(1,:)~=90;
location_matrix = location_matrix(:,oklocs);
ILD_matrix = ILD_matrix(:,oklocs);

nhrtfs = size(ILD_matrix,2);
count = 0;
for cell_num = 1:size(test_numbers,1)

%2. Get the mean response surface of the Tonal ILD/Freq data
[Tonal_meansurf, Tonal_stdsurf, Tonal_dim2vals, Tonal_dim1vals, Tonal_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell_num,1), 1, 0);
%freq_tuning = max(Tonal_meansurf,[],2);
%freq_tuning = freq_tuning./max(freq_tuning);
nfreqs = length(Tonal_dim1vals);
%Generate the PREDICTED ILD Alone surface from the Tonal Data
[Tonal_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
   hrtf_freqs, location_matrix, Tonal_meansurf, Tonal_dim1vals, Tonal_dim2vals);
[Tonal_azi,Tonal_ele,Tonal_diamond] = array2diamond(Tonal_prediction,location_matrix);
Tonal_diamond = Tonal_diamond./max(max(Tonal_diamond));

%3. Get the mean response surface of the BP ILD/Freq data
[BP_meansurf, BP_stdsurf, BP_dim2vals, BP_dim1vals, BP_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell_num,2), 1, 0);

%Generate the PREDICTED ILD Alone surface from the BP Data
[BP_prediction, dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
   hrtf_freqs, location_matrix, BP_meansurf, BP_dim1vals, BP_dim2vals);
[BP_azi,BP_ele,BP_diamond] = array2diamond(BP_prediction,location_matrix);
BP_diamond = BP_diamond./max(max(BP_diamond));

%3. Get the mean response surface of the ILDAlone data and the ILDAlone diamond
[ILDAlone_meansurf, ILDAlone_stdsurf, ILDAlone_dim2vals, ILDAlone_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell_num,3), 1, 0);
max_ILDAlone_meansurf = max(ILDAlone_meansurf);
[ILDAlone_azi,ILDAlone_ele,temp] = array2diamond(ILDAlone_meansurf,locs');
   %Interpolate the missing values in the ILDAlone measurement
   [AZ EL] = meshgrid(ILDAlone_azi, ILDAlone_ele);

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(locs(:,1)==el & locs(:,2)==az))
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
	
   ILDAlone_diamond = intval + temp;
   ILDAlone_diamond = ILDAlone_diamond./max(max(ILDAlone_diamond));
   
   
   %4. Generate input pattern (Tonal_meansurf) and target (ILDAlone_diamond)
   if(length(Tonal_dim1vals) == 21 & length(BP_dim1vals) == 21 &...
         length(Tonal_dim2vals) == 16 & length(Tonal_dim2vals) == 16)
      count = count +1;
      Input_Pattern(:,count) = Tonal_meansurf(:);
      %Input_Pattern = repmat(Tonal_meansurf(:),1,nhrtfs);
      %Input_Pattern = Input_Pattern';
      %[locs temp] = diam2array(ILDAlone_azi, ILDAlone_ele, ILDAlone_diamond, location_matrix);
      Target_Pattern(:,count) = BP_meansurf(:);
      %Target_Pattern = repmat(temp,length(Tonal_meansurf(:)),1);
   end

end %end loop over cells

%5. Set up the neural network model
net = newff([min(min(Input_Pattern)) max(max(Input_Pattern))],...
   [size(Input_Pattern,1)],{'purelin'},'trainrp');
net.inputs{1}.size = length(Input_Pattern);
net.b{1} = zeros(length(Input_Pattern),1);
net.biasConnect = 0;
%net.trainFcn = 'trainwb';
%net.adaptFcn = 'adaptwb';
%net.inputWeights{:,:}.learnFcn = 'learngdm';
%net.layerWeights{:,:}.learnFcn = 'learngdm';


%6. Train the network
net.trainParam.show = 5;
net.TrainParam.epochs = 50;
net.TrainParam.goal = 1e-4;
[net tr] = train(net,Input_Pattern',Target_Pattern);

  
  %7. Make the relevant outputs
  for cell_num = 1:size(Input_Pattern,2)
[temp] = sim(net,Input_Pattern(:,cell_num));
temp2 = zeros(size(BP_meansurf,1),size(BP_meansurf,2));
temp2(:) = temp;
New_BP_meansurf(:,:,cell_num) = temp;
%Generate the PREDICTED ILD Alone surface from the Modeled BP Data
[New_BP_prediction(cell_num,:), dirs, full_rs, new_f_ax, new_ild_ax] = if2space(ILD_matrix',...
   hrtf_freqs, location_matrix, New_BP_meansurf(:,:,cell_num), Tonal_dim1vals, Tonal_dim2vals);
[New_BP_azi,New_BP_ele,New_BP_diamond(:,:,cell_num)] = array2diamond(New_BP_prediction(cell_num,:),location_matrix);
New_BP_diamond(:,:,cell_num) = New_BP_diamond(:,:,cell_num)./max(max(New_BP_diamond(:,:,cell_num)));
end
return
figure
subplot(1,2,1)
plotdiam(ILDAlone_azi,ILDAlone_ele,ILDAlone_diamond);
colorbar
title('Measured')

subplot(1,2,2)
plotdiam(New_BP_azi,New_BP_ele,New_BP_diamond);
colorbar
title('Modeled')

figure
subplot(3,1,1)
plotsurf(Tonal_dim1vals,Tonal_dim2vals,Tonal_meansurf');
colorbar
title('Tonal, Measured')

subplot(3,1,2)
plotsurf(BP_dim1vals,BP_dim2vals,BP_meansurf');
colorbar
title('BP, Measured')

subplot(3,1,3)
plotsurf(Tonal_dim1vals,Tonal_dim2vals,New_BP_meansurf');
colorbar
title('BP, Modeled')


return

