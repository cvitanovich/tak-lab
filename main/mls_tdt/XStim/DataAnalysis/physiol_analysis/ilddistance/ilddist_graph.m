%Script to graphically display the ILD distance idea

clear
%set parameters
bird_number = 899;
cmap = 1-gray;
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
nntwarn off

%load the Neuron structures
load 'd:\mlspezio\matlab\save\Neuron_24'
clear ITD_matrix ABI_matrix

%Limit frequencies used
minfreq = 2000; maxfreq = 11000;
[y,minfreqind] = min(abs(HRTFinfo.hrtf_freqs - minfreq));
[y,maxfreqind] = min(abs(HRTFinfo.hrtf_freqs - maxfreq));


for cell_num = 6:6
   
   %Find the best ILDAlone RS location
   IA_bestloc_ind(cell_num) = max(find(Neuron{cell_num}.ia_meansurf{1} == ...
      max(Neuron{cell_num}.ia_meansurf{1})));
   if(cell_num == 1)  IA_bestloc_ind(cell_num) = 129; end
   if(cell_num == 4)  IA_bestloc_ind(cell_num) = 198; end
   if(cell_num == 7)  IA_bestloc_ind(cell_num) = 163; end
   if(cell_num == 9)  IA_bestloc_ind(cell_num) = 282; end
   if(cell_num == 11) IA_bestloc_ind(cell_num) = 127; end
   if(cell_num == 13) IA_bestloc_ind(cell_num) = 130; end
   if(cell_num == 16) IA_bestloc_ind(cell_num) = 159; end
   if(cell_num == 18) IA_bestloc_ind(cell_num) =  78; end
   if(cell_num == 24) IA_bestloc_ind(cell_num) = 128; end
   
	%Get those locations in the ILDmatrix that match the measured ILDAlone locations
	for num_loc = 1:size(Neuron{cell_num}.ia_locs{1},1)
   	ILDmat_index(num_loc) =...
      	max(find(HRTFinfo.location_matrix(1,:) == Neuron{cell_num}.ia_locs{1}(num_loc,1) &...
      	HRTFinfo.location_matrix(2,:) == Neuron{cell_num}.ia_locs{1}(num_loc,2)));
	end
	ILD_matrix_focus = ILD_matrix(:,ILDmat_index);
   
   
   fig_ht = 9.5; fig_wd = 7.5; mainfig = figure; set(mainfig,'Units','inches'); 
   set(mainfig,'Position',[0 0 fig_wd fig_ht]);
   
   subplot('Position',[0.05 0.65 0.25 0.25])
   set(gca,'FontSize',8)
   plotsurf(Neuron{cell_num}.ia_azi{1},Neuron{cell_num}.ia_ele{1},Neuron{cell_num}.ia_diamond{1});
   axis square
   colormap(cmap); colorbar;
   xlabel('Azimuth (\circ)'); ylabel('Elevation (\circ)');
   freqs(1) = Neuron{cell_num}.ia_header{1}(16);
   freqs(2) = Neuron{cell_num}.ia_header{1}(17);
   title(str2mat(['Measured ILDAlone RS, cell # ' num2str(cell_num)],...
      ['freqs = (' num2str(freqs(1)) ',' num2str(freqs(2)) ')']));
   
   subplot('Position',[0.38 0.65 0.25 0.25])
   set(gca,'FontSize',8)
   plotsurf(Neuron{cell_num}.ia_azi{2},Neuron{cell_num}.ia_ele{2},Neuron{cell_num}.ia_diamond{2});
   axis square
   colormap(cmap); colorbar;
   xlabel('Azimuth (\circ)'); ylabel('Elevation (\circ)');
   if(cell_num < 21)
      freqs(1,1) = Neuron{cell_num}.ia_header{2}(14) - (Neuron{cell_num}.ia_header{2}(15)/2);
      freqs(1,2) = Neuron{cell_num}.ia_header{2}(14) + (Neuron{cell_num}.ia_header{2}(15)/2);
   else
      freqs(1,1) = Neuron{cell_num}.ia_header{2}(16);
      freqs(1,2) = Neuron{cell_num}.ia_header{2}(17);
   end
   title(str2mat(['Measured ILDAlone RS, cell # ' num2str(cell_num)],...
      ['freqs = (' num2str(freqs(1,1)) ',' num2str(freqs(1,2)) ')']));
   
   subplot('Position',[0.71 0.65 0.25 0.25])
   set(gca,'FontSize',8)
   plotsurf(Neuron{cell_num}.ia_azi{4},Neuron{cell_num}.ia_ele{4},Neuron{cell_num}.ia_diamond{4});
   axis square
   colormap(cmap); colorbar;
   xlabel('Azimuth (\circ)'); ylabel('Elevation (\circ)');
   if(cell_num < 21)
      freqs(2,1) = Neuron{cell_num}.ia_header{4}(14) - (Neuron{cell_num}.ia_header{4}(15)/2);
      freqs(2,2) = Neuron{cell_num}.ia_header{4}(14) + (Neuron{cell_num}.ia_header{4}(15)/2);
   else
      freqs(2,1) = Neuron{cell_num}.ia_header{4}(16);
      freqs(2,2) = Neuron{cell_num}.ia_header{4}(17);
   end
   title(str2mat(['Measured ILDAlone RS, cell # ' num2str(cell_num)],...
      ['freqs = (' num2str(freqs(2,1)) ',' num2str(freqs(2,2)) ')']));
   
   switch cell_num
   case 19,
      az_limits(2,:) = [-90 -55]; az_limits(1,:) = [-40  0];
      el_limits(2,:) = [-20   0]; el_limits(1,:) = [ 30 80];
   case 20,
      az_limits(1,:) = [-90 -50]; az_limits(2,:) = [40  60];
      el_limits(1,:) = [-20   0]; el_limits(2,:) = [-15 -5];
   case 21,
      az_limits(1,:) = [ 35  40]; az_limits(2,:) = [-40  -5];
      el_limits(1,:) = [-45 -30]; el_limits(2,:) = [-65 -55];
   case 22,
      az_limits(1,:) = [-40 -10]; az_limits(2,:) = [20 50];
      el_limits(1,:) = [-65 -50]; el_limits(2,:) = [25 35];
   case 24,
      az_limits(1,:) = [-20 -5]; az_limits(2,:) = [50 60];
      el_limits(1,:) = [5 15]; el_limits(2,:) = [-30 -10];
   end %end switch
   
   %Calculate the distance matrix
   disp('Calculating ILD distances...')
   [dist_ild(:,:,cell_num)] = ...
      get_ilddist_all(ILD_matrix_focus(minfreqind:maxfreqind,:),IA_bestloc_ind(cell_num));
   
   disp('Calculating correlations...')
   ILD_freqs = HRTFinfo.hrtf_freqs(minfreqind:maxfreqind);
   for n = 1:length(ILD_freqs)
    	temp = corrcoef(Neuron{cell_num}.ia_meansurf{1},dist_ild(n,:,cell_num));
    	cc_actfreq(n,cell_num) = temp(1,2);
   end
   cc_actfreq(:,cell_num) = abs(cc_actfreq(:,cell_num));
   
   %Parse the ILD matrix & plot
   for n = 1:2
      if(n == 1)subplot('Position',[0.07 0.37 0.8 0.25]);end
      if(n == 2)subplot('Position',[0.07 0.05 0.8 0.25]);end
      hold on
      set(gca,'FontSize',8);
      [haxes,hline1,hline2] = ...
         plotyy(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
         ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind(cell_num)),...
         HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
         cc_actfreq(:,cell_num));
      axes(haxes(1));
      set(gca,'YColor','black','Ytick',[-40:10:40],'Box','off');
      xlim([2000 11000]); ylim([-45 45]);
      axes(haxes(2));
      xlim([2000 11000]); ylim([0 0.8]);
      set(gca,'FontSize',8,'YColor','magenta','Ytick',[0:0.2:0.8]);
      set(hline2,'Color','magenta','LineWidth',3);
      line(2000:1000:11000,0.4*ones(1,10),'Color','black','LineStyle','--');
      ylabel('Correlation (r)','Rotation',270);
      axes(haxes(1));
      for f = 1:2
         [freqindex,locindex,parsed_ILDmat] = ...
            parse_HRTF(ILD_matrix,...
            HRTFinfo.hrtf_freqs,...
            HRTFinfo.location_matrix,...
            freqs(f,:),...
            az_limits(n,:),...
            el_limits(n,:));
         if(f == 1)
            plot(HRTFinfo.hrtf_freqs(freqindex),parsed_ILDmat,'g--','LineWidth',1);
            hold on
         end
         if(f == 2)
            plot(HRTFinfo.hrtf_freqs(freqindex),parsed_ILDmat,'r-.','LineWidth',1);
            line(freqs(f,1)*ones(1,10),-45:90/9:45);
         end
      end
      xlabel('Frequency (Hz)','FontSize',8); ylabel('ILD (dB)','FontSize',8);
      title(['Location range: AZ = (' ...
            num2str(az_limits(n,1)) ',' ...
            num2str(az_limits(n,2)) '), EL = (' ...
            num2str(el_limits(n,1)) ',' ...
            num2str(el_limits(n,2)) ')']);
      axes(haxes(2));
   end
         
         
         
   
end
