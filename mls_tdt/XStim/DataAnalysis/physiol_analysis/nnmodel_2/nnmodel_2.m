%Script to set up and train a neural network model of ICCls or ICx cells, 
%using inputs from tonal ILD/freq-tuned cells and weights that are a function
%of the Euclidean distance between tonal cells' ILD/freq tuning and the optimal ILD/freq
%spectrum of the ICCls or ICx cell

clear

%Load network made by nnmodel_2_setup
load 'd:\mlspezio\matlab\save\ICspacenet_init'
net_temp = ICspacenet;
clear ICspacenet

%****************
%Specify Input (P)
%****************
%set parameters
bird_number = 899;
cmap = jet;
linecolor = 'white';
bkgndcolor = [0.85 0.85 0.85];
%Directories & files
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\HRTF_files\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];
get_HRTF = 0;
get_ITDmatrix = 0;

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'
clear ITD_matrix ABI_matrix

load d:\mlspezio\matlab\save\nnet_24b
net = nnet_24b;
for cell_num = 12:24
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')))
      
      ICspacenet{cell_num} = net_temp;
      numfreqs = length(Neuron{cell_num}.tif_freqaxis);
      numILDs  = length(Neuron{cell_num}.tif_ildaxis);
      
      ICspacenet{cell_num}.layers{1}.size = numfreqs*numILDs; %Layer of tonal ILD/freq-tuned units
      ICspacenet{cell_num}.layers{1}.dimensions = [numfreqs numILDs];
      
      %Get the contour for the SRF of the neuron
      if(any(strcmp(fieldnames(Neuron{cell_num}),'ts_meansurf')));
         [rf_center_az, rf_center_el, srflim, rf_area, rf_i, space_contour] = ...      
            srfstats(Neuron{cell_num}.ts_azi{1},...
            Neuron{cell_num}.ts_ele{1},...
            Neuron{cell_num}.ts_diamond{1}, 0.5, 0);
      end
      
      %Get those locations in the ILDmatrix that match the measured ILDAlone locations
      for num_loc = 1:size(Neuron{cell_num}.ia_locs{1},1)
         ILDmat_index(num_loc) =...
            max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(num_loc,1) &...
            HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(num_loc,2)));
      end
      ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
      
      %Find the best ILDAlone RS location
      IA_bestloc_ind(cell_num) = max(find(Neuron{cell_num}.ia_meansurf{1} == max(Neuron{cell_num}.ia_meansurf{1})));
      if(cell_num == 1) IA_bestloc_ind(cell_num) = 129; end
      if(cell_num == 4) IA_bestloc_ind(cell_num) = 198; end
      if(cell_num == 7) IA_bestloc_ind(cell_num) = 163; end
      if(cell_num == 9) IA_bestloc_ind(cell_num) = 282; end
      if(cell_num == 11) IA_bestloc_ind(cell_num) = 127; end
      if(cell_num == 13) IA_bestloc_ind(cell_num) = 130; end
      if(cell_num == 16) IA_bestloc_ind(cell_num) = 159; end
      if(cell_num == 18) IA_bestloc_ind(cell_num) = 78; end
      if(cell_num == 24) IA_bestloc_ind(cell_num) = 128; end
      
      %Select frequency range matching the tif_meansurf
      minfreq = min(Neuron{cell_num}.tif_freqaxis); maxfreq = max(Neuron{cell_num}.tif_freqaxis);
      [y,ind_minfreq] = min(abs(HRTFinfo.hrtf_freqs - minfreq));
      [y,ind_maxfreq] = min(abs(HRTFinfo.hrtf_freqs - maxfreq));
      
      ILD_matrix_new = interp1(HRTFinfo.hrtf_freqs(ind_minfreq:ind_maxfreq),...
         ILD_matrix_focus(ind_minfreq:ind_maxfreq,:),...
         Neuron{cell_num}.tif_freqaxis,'spline');
      
      %Set input ranges
      minILD = min(min(ILD_matrix_new)); maxILD = max(max(ILD_matrix_new));
      ICspacenet{cell_num}.inputs{1}.range = [minILD*ones(numfreqs*numILDs,1) maxILD*ones(numfreqs*numILDs,1)];
      
      %Set inputs and targets - each input/target pair is an ILD spectrum for a location
      %and the unit's response to that location
      [FREQs,ILDs] = meshgrid(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis);
      tempP = zeros(length(Neuron{cell_num}.tif_ildaxis),length(Neuron{cell_num}.tif_freqaxis),size(ILD_matrix_new,2));
      for loc = 1:size(ILD_matrix_new,2)
         for freq = 1:size(ILD_matrix_new,1)
            [y,ildindex]  = min(abs(Neuron{cell_num}.tif_ildaxis  - ILD_matrix_new(freq,loc)));
            tempP(ildindex,freq,loc) = 1;
         end
         %P is the input we use. Each column is a location containing numfreqs*numILDs values
         %1's are in the appropriate matrix cells for a location's ILD/freq spectrum
         P(:,loc) = reshape(tempP(:,:,loc),...
            length(Neuron{cell_num}.tif_ildaxis)*length(Neuron{cell_num}.tif_freqaxis),1);
      end
      TonalP = zeros(numfreqs*numILDs,numfreqs*numILDs);
      for inputnum = 1:size(TonalP,2)
         TonalP(inputnum,inputnum) = 1;
      end
      P = [P TonalP]; %Include Tonal inputs for training set
      ICspacenet{cell_num}.inputs{1}.range = minmax(P);
      
      T = Neuron{cell_num}.ia_meansurf{1}; %ILDAlone targets
      b = ones(1,numfreqs*numILDs);
      b(:) = Neuron{cell_num}.tif_meansurf';
      T = [T b]; clear b;
      
      %Set the input weights to 1 for all identical ILD/Freq values, 0 otherwise
      b = ones(numILDs*numfreqs,numILDs*numfreqs);
      ICspacenet{cell_num}.IW{1} = spdiags(b,0,numILDs*numfreqs,numILDs*numfreqs); %Set all appropriate input weights to 1
      clear b
      
      %Do the following if the input weights are to be set according to the cell's tonal response
      if(1)
         b = ones(numfreqs*numILDs,1);
         b(:) = Neuron{cell_num}.tif_meansurf'; c = ICspacenet{cell_num}.IW{1};
         ICspacenet{cell_num}.IW{1} = spdiags(b,0,c); %Replace main diag in c with b
         clear b c
      end
      
      %Initialize the weights from Layer 1 to Layer 1
      %rand('state',sum(100*clock));
      %tempmat = rand(numfreqs*numILDs,numfreqs*numILDs);
      %tempmat(find(tempmat < 0.5)) = -1; tempmat(find(tempmat >= 0.5)) = 1;
      %ICspacenet{cell_num}.LW{1,1} = 0.01*tempmat;
      %ILDs_vec = zeros(numILDs*numfreqs,1);
      %ILDs_vec(:) = ILDs;
      %ILDwidth = 10;
      %tempmat = zeros(numILDs*numfreqs,numILDs*numfreqs);
      %for ild1 = 1:length(ILDs_vec)
      %   dist = sqrt((ILDs_vec - ILDs_vec(ild1)).^2);
      %   tempmat(:,ild1) = 1 - (dist/ILDwidth);
      %end
      %ICspacenet{cell_num}.LW{1,1} = tempmat; clear tempmat;
      
      %Initialize the weights from Layer 1 to Layer 2,
      %based on distance measure from optimal ILDs
      ILDwidth = 10;
      for freq = 1:numfreqs
         for ild = 1:numILDs
            dist = sqrt((ILDs(ild,freq) -...
               ILD_matrix_new(freq,IA_bestloc_ind(cell_num))).^2);
            a(ild,freq) = 1 - (dist/ILDwidth);
         end
      end
      index = find(Neuron{cell_num}.tif_meansurf' ~= 0);
      b = zeros(numILDs,numfreqs);
      b(index) = 1;
      ICspacenet{cell_num}.LW{2,1}(:) = a;
      clear a;
      %ICspacenet{cell_num}.LW{2,1}(:) = a.*b;
      %ICspacenet{cell_num}.LW{2,1}(:) = Neuron{cell_num}.tif_meansurf';
      %ICspacenet{cell_num}.LW{2,1} = ones(1,numILDs*numfreqs);
      
      clear tempP ILD_matrix_focus
      
      %Train network
      figure(50);
      [net{cell_num},tr{cell_num}] = train(ICspacenet{cell_num},P,T);
      nnet_24b{cell_num} = net{cell_num};
      save d:\mlspezio\matlab\save\nnet_24b nnet_24b
      
      %Plot
      fig_ht = 9.5;
      fig_wd = 7.5;
      mainfig = figure;
      set(mainfig,'Units','inches');
      set(mainfig,'Position',[0 0 fig_wd fig_ht]);
      
      %Simulate network and calculate r^2
      IA_nnpredict = sim(net{cell_num},P(:,1:359));
      subplot(3,3,1)
      set(gca,'FontSize',8);
      plot(Neuron{cell_num}.ia_meansurf{1},IA_nnpredict,'.');
      cc = corrcoef(Neuron{cell_num}.ia_meansurf{1},IA_nnpredict);
      vari = cc(2,1)^2;
      xlabel('Measured ILDAlone activity for 359 locations','FontSize',8);
      ylabel('Modeled ILDAlone activity for 359 locations','FontSize',8);
      title(['ILDAlone activity comparison, r^{2} = ' num2str(vari)],'FontSize',8);
      
      %Measured ILDAlone Diamond
      subplot(3,3,2)
      set(gca,'FontSize',8);
      axis square
      plotdiam(Neuron{cell_num}.ia_azi{1},Neuron{cell_num}.ia_ele{1},Neuron{cell_num}.ia_diamond{1});
      colorbar
      colormap(cmap);
      xlabel('Azimuth (deg)','FontSize',8);
      ylabel('Elevation (deg)','FontSize',8);
      title(['Measured ILDAlone for Cell #' num2str(cell_num)],'FontSize',8);
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
         hold off;
      end

      %Modeled ILDAlone Diamond
      %make the diamond for the prediction
      [IA_nnpredict_azi,IA_nnpredict_ele,temp] = ...
         array2diamond(IA_nnpredict,Neuron{cell_num}.ia_locs{1}');
      %Interpolate the missing values in the ILDAlone measurement
      [AZ EL] = meshgrid(IA_nnpredict_azi, IA_nnpredict_ele);
      
      % generate mask for missing points
      missmask = NaN*ones(size(temp));
      i = 1;
      for az = -90:5:90;
         for el = -90+abs(az):5:90-abs(az)
            if (~(Neuron{cell_num}.ia_locs{1}(:,1)==el & Neuron{cell_num}.ia_locs{1}(:,2)==az))
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
      IA_nnpredict_diam = intval + temp;
      %PLOT
      subplot(3,3,3)
      axis square
      plotdiam(IA_nnpredict_azi,IA_nnpredict_ele,IA_nnpredict_diam);
      colorbar
      colormap(cmap);
      xlabel('Azimuth (deg)','FontSize',8);
      ylabel('Elevation (deg)','FontSize',8);
      title('Modeled ILDAlone','FontSize',8);
      if(exist('space_contour'))
         hold on;
         plot(space_contour(1,:),space_contour(2,:),linecolor,'LineWidth',1.5);
         hold off;
      end
      
      %Simulate Tonal response and calculate r^2
      z = sim(net{cell_num},TonalP);
      negind = find(z < 0); z(negind) = 0;
      subplot(3,3,4)
      set(gca,'FontSize',8);
      a = zeros(numfreqs*numILDs,1);
      a(:) = Neuron{cell_num}.tif_meansurf';
      plot(a,z,'.');
      cc = corrcoef(a,z);
      vari = cc(2,1)^2;
      xlabel('Measured Tonal activity','FontSize',8);
      ylabel('Modeled Tonal activity','FontSize',8);
      title(['Tonal activity comparison, r^{2} = ' num2str(vari)],'FontSize',8);
      
      %Measured Tonal ILD/Freq surface
      subplot(3,3,5)
      set(gca,'FontSize',8);
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,Neuron{cell_num}.tif_meansurf');
      hold on
      colorbar
      colormap(cmap);
      plot(Neuron{cell_num}.tif_freqaxis,ILD_matrix_new(:,IA_bestloc_ind(cell_num)),'yellow',...
         'LineWidth',2);
      axis([min(Neuron{cell_num}.tif_freqaxis) max(Neuron{cell_num}.tif_freqaxis)...
            min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
      xlabel('Freq (Hz)','FontSize',8);
      ylabel('ILD (dB)','FontSize',8);
      title('Measured Tonal ILD/Freq','FontSize',8);

      %Modeled Tonal ILD/Freq surface
      subplot(3,3,6)
      set(gca,'FontSize',8);
      tif_nnpredict = zeros(numILDs,numfreqs);
      tif_nnpredict(:) = z;
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,tif_nnpredict);
      hold on
      colorbar
      colormap(cmap);
      plot(Neuron{cell_num}.tif_freqaxis,ILD_matrix_new(:,IA_bestloc_ind(cell_num)),'yellow',...
         'LineWidth',2);
      axis([min(Neuron{cell_num}.tif_freqaxis) max(Neuron{cell_num}.tif_freqaxis)...
            min(Neuron{cell_num}.tif_ildaxis)  max(Neuron{cell_num}.tif_ildaxis)]);
      xlabel('Freq (Hz)','FontSize',8);
      ylabel('ILD (dB)','FontSize',8);
      title('Modeled Tonal ILD/Freq','FontSize',8);
      
      %Trained input weights
      subplot(3,3,8)
      set(gca,'FontSize',8);
      a = zeros(numILDs,numfreqs);
      a(:) = diag(net{cell_num}.IW{1});
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,a);
      shading flat
      colorbar
      colormap(cmap)
      xlabel('Freq (Hz)','FontSize',8); ylabel('ILD (dB)','FontSize',8);
      title('Input weights, after training','FontSize',8)
      
      %Trained Layer 1-->2 weights
      subplot(3,3,9)
      set(gca,'FontSize',8);
      a = zeros(numILDs,numfreqs);
      a(:) = net{cell_num}.LW{2,1};
      plotsurf(Neuron{cell_num}.tif_freqaxis,Neuron{cell_num}.tif_ildaxis,a);
      shading flat
      colorbar
      colormap(cmap)
      xlabel('Freq (Hz)','FontSize',8); ylabel('ILD (dB)','FontSize',8);
      title('Weights for layer 1-->2, after training','FontSize',8)
      hold on
      plot(Neuron{cell_num}.tif_freqaxis,ILD_matrix_new(:,IA_bestloc_ind(cell_num)),'k',...
         'LineWidth',2);
      
      clear P T space_contour IA_nnpredict b z a
   end %end if for tif_meansurf
end %end loop over cells
      


