clear

load d:\mlspezio\matlab\save\Neuron_24
clear ITD_matrix ILD_matrix ABI_matrix test_numbers

restrict = 0;
frac_act = [0.7 0.2];

bins = [0:.1:300];
fig_ht = 9.5; fig_wd = 7.5;      

for cell_num = 1:24
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   mainfig = figure; set(mainfig,'Units','inches'); 
   set(mainfig,'Position',[0 0 fig_wd fig_ht]);
   subplot(2,1,1)
   set(gca,'FontSize',8);
	%Tonal
	data_array = Neuron{cell_num}.tif_dataarray;
   if(restrict)
      
      tif_index1 = find(Neuron{cell_num}.tif_meansurf >= ...
         frac_act(1) * max(max(Neuron{cell_num}.tif_meansurf)));
      [r,c] = ind2sub(size(Neuron{cell_num}.tif_meansurf),tif_index1);
      
      tif_index1a = [];
      for rnum = 1:length(r)
         tif_index1a = [tif_index1a;find(data_array(:,4) == r(rnum) & data_array(:,5) == c(rnum))];
      end
      [i1,i2] = ind2sub(size(data_array),tif_index1a);
      
      %Calculate isi info
      [tif_isi1a{cell_num},tif_mean_isi1a(cell_num),tif_std_isi1a(cell_num),tif_Cv1a(cell_num)] = ...
         get_isi(data_array(i1,:),1,3);
      
      n = hist(data_array(i1,1),bins);
      n1 = n/length(tif_index1a);
      %n = n/max(n);
      plot(bins,n1,'r-');
      hold on
      
      minval = min(min(Neuron{cell_num}.tif_meansurf));
      if(minval == 0) minval = 0.01; end
      tif_index2 = find(Neuron{cell_num}.tif_meansurf <= ...
         minval + (frac_act(2) * minval));
      clear r c
      [r,c] = ind2sub(size(Neuron{cell_num}.tif_meansurf),tif_index2);
      tif_index2a = [];
      for rnum = 1:length(r)
         tif_index2a = [tif_index2a;find(data_array(:,4) == r(rnum) & data_array(:,5) == c(rnum))];
      end
      [i1,i2] = ind2sub(size(data_array),tif_index2a);
      
      %Calculate isi info
      [tif_isi2a{cell_num},tif_mean_isi2a(cell_num),tif_std_isi2a(cell_num),tif_Cv2a(cell_num)] = ...
         get_isi(data_array(i1,:),1,3);
      
      n = hist(data_array(i1,1),bins);
      %n = n/max(n);
      n2 = n/length(tif_index2a);
      plot(bins,n2,'b--');
      line(100*ones(1,10),0:max(max([n1 n2]))/9:max(max([n1 n2])),'Color','m');
      title(str2mat(['Cell # ' num2str(cell_num) ', Tonal Data'],...
         ['Red >= ' num2str(100*frac_act(1)) ...
            '% MaxAct, Blue <= MinAct + ' num2str(100*frac_act(2)) '%']));
      
   else
      n = hist(data_array(:,1),bins);
      %n = n/max(n);
      plot(bins,n);
      title(['Cell # ' num2str(cell_num) ', Tonal Data (all trials)'])
      
      %Calculate isi info
      [tif_isi{cell_num},tif_mean_isi(cell_num),tif_std_isi(cell_num),tif_Cv(cell_num)] = ...
         get_isi(data_array,1,3);
   end
   
   
   %ILDAlone
   clear data_array
   data_array = Neuron{cell_num}.ia_dataarray;
   set(mainfig,'Position',[0 0 fig_wd fig_ht]);
   subplot(2,1,2)
   set(gca,'FontSize',8);
   if(restrict)
      
      ia_index1 = find(Neuron{cell_num}.ia_meansurf{1} >= ...
         frac_act(1) * max(max(Neuron{cell_num}.ia_meansurf{1})));
      
      ia_index1a = [];
      for rnum = 1:length(ia_index1)
         ia_index1a = [ia_index1a;...
               find(data_array(:,4) == Neuron{cell_num}.ia_locs{1}(ia_index1(rnum),1) &...
               data_array(:,5) == Neuron{cell_num}.ia_locs{1}(ia_index1(rnum),2))];
      end
      [i1,i2] = ind2sub(size(data_array),ia_index1a);
      
      %Calculate isi info
      [ia_isi1a{cell_num},ia_mean_isi1a(cell_num),ia_std_isi1a(cell_num),ia_Cv1a(cell_num)] = ...
         get_isi(data_array(i1,:),1,3);
      
      n = hist(data_array(i1,1),bins);
      %n = n/max(n);
      n1 = n/length(ia_index1a);
      plot(bins,n1,'r-');
      hold on
      
      minval = min(min(Neuron{cell_num}.ia_meansurf{1}));
      if(minval == 0) minval = 0.01; end
      ia_index2 = find(Neuron{cell_num}.ia_meansurf{1} <= ...
         minval + (frac_act(2) * minval));
      ia_index2a = [];
      for rnum = 1:length(ia_index2)
         ia_index2a = [ia_index2a;...
               find(data_array(:,4) == Neuron{cell_num}.ia_locs{1}(ia_index2(rnum),1) &...
               data_array(:,5) == Neuron{cell_num}.ia_locs{1}(ia_index2(rnum),2))];
      end
      [i1,i2] = ind2sub(size(data_array),ia_index2a);
      
      %Calculate isi info
      [ia_isi2a{cell_num},ia_mean_isi2a(cell_num),ia_std_isi2a(cell_num),ia_Cv2a(cell_num)] = ...
         get_isi(data_array(i1,:),1,3);
      
      n = hist(data_array(i1,1),bins);
      %n = n/max(n);
      n2 = n/length(ia_index2a);
      plot(bins,n2,'b--');
      line(100*ones(1,10),0:max(max([n1 n2]))/9:max(max([n1 n2])),'Color','m');
      title(str2mat(['Cell # ' num2str(cell_num) ', IA Data'],...
         ['Red >= ' num2str(100*frac_act(1)) ...
            '% MaxAct, Blue <= MinAct + ' num2str(100*frac_act(2)) '%']));
      
   else
      n = hist(data_array(:,1),bins);
      %n = n/max(n);
      plot(bins,n);
      title(['Cell # ' num2str(cell_num) ', IA Data (all trials)'])
      
      %Calculate isi info
      [ia_isi{cell_num},ia_mean_isi(cell_num),ia_std_isi(cell_num),ia_Cv(cell_num)] = ...
         get_isi(data_array,1,3);
   end
   
   
   end %if for tif_meansurf
end %loop over cells