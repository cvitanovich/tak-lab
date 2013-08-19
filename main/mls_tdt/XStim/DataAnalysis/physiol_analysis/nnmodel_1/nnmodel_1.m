function [net,New_ILDAlone_diamond,New_ILDf_surface] = ...
   nnmodel_1(bird_number,side_of_brain,test_numbers,hrtf_file,get_hrtf)
%Function to model the transformation of a space-tuned cell's Tonal ILD/Freq 
%Response Surface in order to produce an ILD/Freq RS that will predict the ILDAlone RS
%This is a two-layer network that has the following architecture:
%        Input --> Layer 1 (Linear) --> Layer 2 (Linear) --> Output
%Input is the Tonal ILD/Freq RS, Output is the ILDAlone RS

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

%2. Get the mean response surface of the Tonal ILD/Freq data
[Tonal_meansurf, Tonal_stdsurf, Tonal_ildaxis, Tonal_freqaxis, Tonal_testpars, spont_spikes, spont_dur, nincl_reps] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(1), 1, 0);
%freq_tuning = max(Tonal_meansurf,[],2);
%freq_tuning = freq_tuning./max(freq_tuning);
nfreqs = length(Tonal_freqaxis);

%3. Get the mean response surface of the ILDAlone data and the ILDAlone diamond
[ILDAlone_meansurf, ILDAlone_stdsurf, ILDAlone_dim2vals, ILDAlone_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(2), 1, 0);
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
Input_Pattern = Tonal_meansurf(:)';
%Input_Pattern = repmat(Tonal_meansurf(:),1,nhrtfs);
%Input_Pattern = Input_Pattern';
[locs temp] = diam2array(ILDAlone_azi, ILDAlone_ele, ILDAlone_diamond, location_matrix);
Target_Pattern = temp';
%Target_Pattern = repmat(temp,length(Tonal_meansurf(:)),1);

%5. Set up the neural network model
net = newff([min(Input_Pattern) max(Input_Pattern)],...
   [length(Input_Pattern) length(Target_Pattern)],{'purelin' 'purelin'},'trainrp');
net.inputs{1}.size = length(Input_Pattern);
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
[temp] = sim(net,Input_Pattern');
New_ILDAlone = temp;
clear temp
temp = eval([net.layers{1}.transferFcn '(net.IW{1}*Input_Pattern'' + net.b{1})']);
New_ILDf_surface = zeros(length(Tonal_freqaxis),length(Tonal_ildaxis));
New_ILDf_surface(:) = temp;

figure
subplot(1,2,1)
plotdiam(ILDAlone_azi,ILDAlone_ele,ILDAlone_diamond);
colorbar
title('Measured')

subplot(1,2,2)
[New_ILDAlone_azi,New_ILDAlone_ele,New_ILDAlone_diamond] = array2diamond(New_ILDAlone,location_matrix);
plotdiam(New_ILDAlone_azi,New_ILDAlone_ele,New_ILDAlone_diamond);
colorbar
title('Modeled')

figure
subplot(2,1,1)
plotsurf(Tonal_freqaxis,Tonal_ildaxis,Tonal_meansurf');
colorbar
title('Measured')

subplot(2,1,2)
plotsurf(Tonal_freqaxis,Tonal_ildaxis,New_ILDf_surface');
colorbar
title('Modeled')


return

