clear
pack

load d:\mlspezio\matlab\save\Neuron_24_psth
clear ITD_matrix ILD_matrix ABI_matrix test_numbers

%Paramters
bird_number = 899;
side_of_brain = 'r';
begin_cell = 1;
end_cell = 1;
bandwidth = 1/3;

%Directories & files
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\HRTF_files\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];
get_HRTF = 0;
get_ITDmatrix = 0;

[L_matrix,R_matrix,ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = get_HRTFinfo2(bird_number,...
   hrtf_file,...
   hrtf_file2,...
   get_HRTF,...
   get_ITDmatrix);


for cell_num = begin_cell:end_cell
   
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')))
      
      disp(['Processing cell # ' num2str(cell_num)])
      
      disp('Initializing neural network...')
      ILDf_freqaxis = Neuron{cell_num}.tif_freqaxis;
      ILDf_ildaxis  = Neuron{cell_num}.tif_ildaxis;
      ILDnet = dynsetup2(length(ILDf_freqaxis),length(ILDf_freqaxis)*length(ILDf_ildaxis));
      ILDnet = init(ILDnet);
		ILDnet.IW{1,1} = diag(ones(1,length(ILDf_freqaxis)),0);
		ILDnet.IW{2,2} = diag(ones(1,length(ILDf_freqaxis)),0);
      
      disp('Making inputs...')
      P = make_inputs(Neuron{cell_num},HRTFinfo,L_matrix,R_matrix,bandwidth);
      
      disp('Making targets...')
      T = make_targs(Neuron{cell_num});
      
      clear L_matrix R_matrix ABI_matrix
      disp('Training...')
      [ILDnet_new,tr] = train(ILDnet,P,T);
      
      O = sim(ILDnet_new,P);
      
      for n = 1:length(Neuron{cell_num}.psth_bins)
         T_mat(n,:) = T{n};
         O_mat(n,:) = O{n};
      end
      
      for n = 1:size(T_mat,2)
         temp = corrcoef(T_mat(:,n),O_mat(:,n));
         r(n) = temp(1,2);
      end
      
      figure
      plot(r,'.');
      
   end
end

