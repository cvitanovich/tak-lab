clear

test_numbers = [
%[Thr  SR ABI ITD ABI/F TIF BIF IA1 IA2 IA3 IA4 IA5 IA6 IA7 IA8 IA9 IA10 IA11 IA12 TS1 TS2 TS3 TS4 TS5 TS6 TS7 TS8 TS9 TS10 TS11 TS12 ]
  -57  -1  -1  65    -1  66  68  70  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  71  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -62  -1  74  75    -1  76  77  79  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  80  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1  -1  85    -1  87  88  89  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  90  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -72  -1  -1  91    -1  92  93  94  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  95  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -80  -1  96  97    -1  98 100 101  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -72  -1 103 104    -1 105 106 107  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 108 109    -1 110 111 112  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 113  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 121 122    -1 123 124 125  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 126  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1  -1 127    -1 129 130 131  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 132  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -73  -1 133 134    -1 135 136 137  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 138  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -73  -1 141 142    -1 143 144 145  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 146  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 152 153    -1 162 158 159  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 161  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 165 166    -1 167 168 169  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 170  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 171 172    -1 173 174 175  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 176  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 206 208    -1 209 210 211  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 207 212 213  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70 220 214 216   223 217 218 219  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 221  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 228 230   234 232 233 231  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 229  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 242 248    -1 249 250 247  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 246  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
];

restrict = 1;
frac_act = [0.8 1];

bins = [0:300];
cache_files_directory = 'd:\mlspezio\matlab\scripts\physiol_analysis\899analysis\899_cache';
for cell_num = 1:size(test_numbers,1)
%for cell_num = 1:1
   %Tonal
   data_file_wildcard = sprintf('899rSITE*T%03g.MAT', test_numbers(cell_num,6));
   cachedir = eval(['dir(''' cache_files_directory '\' data_file_wildcard ''')']);
   if (size(cachedir,1)>1) error('Bizarre Error: More than one cache file matches this test number.'); end;
   if (size(cachedir,1)~=0) 
      data_file = cachedir(1).name;
      data_file = [cache_files_directory '\' data_file];
      disp(['Loading cached data from file: ' data_file])
		eval(['load ' data_file ';']);
   else
      mean_resp_surf = -1;
      disp('Error: Cached file doesnt exist, run proc_raw first.');
      return;
   end;
   

   fig_ht = 9.5; fig_wd = 7.5;      
   mainfig = figure; set(mainfig,'Units','inches'); 
   set(mainfig,'Position',[0 0 fig_wd fig_ht]);
   subplot(3,1,1)
   set(gca,'FontSize',8);
   if(restrict)
      
      tif_index1 = find(Neuron{cell_num}.tif_meansurf >= ...
         frac_act(1) * max(max(Neuron{cell_num}.tif_meansurf)));
      [r,c] = ind2sub(size(Neuron{cell_num}.tif_meansurf),tif_index1);
      
      tif_index1a = [];
      for rnum = 1:length(r)
         tif_index1a = [tif_index1a;find(data_array(:,4) == r(rnum) & data_array(:,5) == c(rnum))];
      end
      [i1,i2] = ind2sub(size(data_array),tif_index1a);
      
      n = hist(data_array(i1,1),bins);
      n1 = n/length(tif_index1a);
      %n = n/max(n);
      plot(bins,n1,'r-');
      hold on
      
      minval = min(min(Neuron{cell_num}.tif_meansurf));
      if(minval == 0) minval = 0.01; end
      tif_index2 = find(Neuron{cell_num}.tif_meansurf <= ...
         minval + ...
         (frac_act(2) * minval));
      [r,c] = ind2sub(size(Neuron{cell_num}.tif_meansurf),tif_index2);
      tif_index2a = [];
      for rnum = 1:length(r)
         tif_index2a = [tif_index2a;find(data_array(:,4) == r(rnum) & data_array(:,5) == c(rnum))];
      end
      [i1,i2] = ind2sub(size(data_array),tif_index2a);
      
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
   end
   
   %BP
   clear data_array
   data_file_wildcard = sprintf('899rSITE*T%03g.MAT', test_numbers(cell_num,7));
   cachedir = eval(['dir(''' cache_files_directory '\' data_file_wildcard ''')']);
   if (size(cachedir,1)>1) error('Bizarre Error: More than one cache file matches this test number.'); end;
   if (size(cachedir,1)~=0) 
      data_file = cachedir(1).name;
      data_file = [cache_files_directory '\' data_file];
      disp(['Loading cached data from file: ' data_file])
      eval(['load ' data_file ';']);
   else
      mean_resp_surf = -1;
      disp('Error: Cached file doesnt exist, run proc_raw first.');
      return;
   end;
   
   subplot(3,1,2)
   set(gca,'FontSize',8);
   if(restrict)
      
      bif_index1 = find(Neuron{cell_num}.bif_meansurf >= ...
         frac_act(1) * max(max(Neuron{cell_num}.bif_meansurf)));
      [r,c] = ind2sub(size(Neuron{cell_num}.bif_meansurf),bif_index1);
      
      bif_index1a = [];
      for rnum = 1:length(r)
         bif_index1a = [bif_index1a;find(data_array(:,4) == r(rnum) & data_array(:,5) == c(rnum))];
      end
      [i1,i2] = ind2sub(size(data_array),bif_index1a);
      
      n = hist(data_array(i1,1),bins);
      n1 = n/length(bif_index1a);
      %n = n/max(n);
      plot(bins,n1,'r-');
      hold on
      
      minval = min(min(Neuron{cell_num}.bif_meansurf));
      if(minval == 0) minval = 0.01; end
      bif_index2 = find(Neuron{cell_num}.bif_meansurf <= ...
         minval + ...
         (frac_act(2) * minval));
      [r,c] = ind2sub(size(Neuron{cell_num}.bif_meansurf),bif_index2);
      bif_index2a = [];
      for rnum = 1:length(r)
         bif_index2a = [bif_index2a;find(data_array(:,4) == r(rnum) & data_array(:,5) == c(rnum))];
      end
      [i1,i2] = ind2sub(size(data_array),bif_index2a);
      
      n = hist(data_array(i1,1),bins);
      %n = n/max(n);
      n2 = n/length(bif_index2a);
      plot(bins,n2,'b--');
      line(100*ones(1,10),0:max(max([n1 n2]))/9:max(max([n1 n2])),'Color','m');
      title(str2mat(['Cell # ' num2str(cell_num) ', BP Data'],...
         ['Red >= ' num2str(100*frac_act(1)) ...
            '% MaxAct, Blue <= MinAct + ' num2str(100*frac_act(2)) '%']));
      
   else
      n = hist(data_array(:,1),bins);
      %n = n/max(n);
      plot(bins,n);
      title(['Cell # ' num2str(cell_num) ', BP Data (all trials)'])
   end
   
   
   %ILDAlone
   clear data_array
   data_file_wildcard = sprintf('899rSITE*T%03g.MAT', test_numbers(cell_num,8));
   cachedir = eval(['dir(''' cache_files_directory '\' data_file_wildcard ''')']);
   if (size(cachedir,1)>1) error('Bizarre Error: More than one cache file matches this test number.'); end;
   if (size(cachedir,1)~=0) 
      data_file = cachedir(1).name;
      data_file = [cache_files_directory '\' data_file];
      disp(['Loading cached data from file: ' data_file])
      eval(['load ' data_file ';']);
   else
      mean_resp_surf = -1;
      disp('Error: Cached file doesnt exist, run proc_raw first.');
      return;
   end;
   
   subplot(3,1,3)
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
      
      n = hist(data_array(i1,1),bins);
      %n = n/max(n);
      n1 = n/length(ia_index1a);
      plot(bins,n1,'r-');
      hold on
      
      minval = min(min(Neuron{cell_num}.ia_meansurf{1}));
      if(minval == 0) minval = 0.01; end
      ia_index2 = find(Neuron{cell_num}.ia_meansurf{1} <= ...
         minval + ...
         (frac_act(2) * minval));
      ia_index2a = [];
      for rnum = 1:length(ia_index2)
         ia_index2a = [ia_index2a;...
               find(data_array(:,4) == Neuron{cell_num}.ia_locs{1}(ia_index2(rnum),1) &...
               data_array(:,5) == Neuron{cell_num}.ia_locs{1}(ia_index2(rnum),2))];
      end
      [i1,i2] = ind2sub(size(data_array),ia_index2a);
      
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
   end
   xlabel('Time (ms)'); ylabel('# spikes');
end
   
   
