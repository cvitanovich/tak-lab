%Engage_IAMultRegress
new_ildaxis = -30:30;
vert_pos = [0.50 0.05];

%load e:\spezio\data\DATA894L5.mat
clear Neuron ILD_matrix ILD_matrix_focus

iamr_datafile = get(hdatafile,'String');
load(iamr_datafile);
if(~isstruct(DATA)) return; end

%Run the regression with ILD-Frequency or ILDAlone?
if(get(hregressia,'Value'))
   iamr_testnum = str2num(get(hiatestnum,'String'));
else
   iamr_testnum = str2num(get(hiftestnum,'String'));
end

if_testnum   = str2num(get(hiftestnum,'String'));
iamr_sitenum = XStimParams.recording_site;
if_sitenum   = XStimParams.recording_site;
%iamr_sitenum = 8;
%if_sitenum = 8;
%**************************************************
%First, perform data reduction into avg. spike rate
time_window = [105 205]; %ms
reps_to_include = [];
Neuron.ia_locs = DATA.Site{iamr_sitenum}.test{iamr_testnum}.params.locations;
Neuron.ia_datamatrix = DATA.Site{iamr_sitenum}.test{iamr_testnum}.datamatrix;
[Neuron.ia_mean, Neuron.ia_std] = ...
   proc_test_space(Neuron.ia_datamatrix,...
   Neuron.ia_locs,...
   time_window,...
   reps_to_include);

IA_meansurf = Neuron.ia_mean;
[Neuron.ia_azi,Neuron.ia_ele,temp] = ...
   array2diamond(Neuron.ia_mean,Neuron.ia_locs);
%Interpolate the missing values in the ILDAlone measurement
[AZ EL] = meshgrid(Neuron.ia_azi, Neuron.ia_ele);

% generate mask for missing points
missmask = NaN*ones(size(temp));

i = 1;
for az = -90:5:90;
   for el = -90+abs(az):5:90-abs(az)
      if (~(Neuron.ia_locs(1,:)==el & Neuron.ia_locs(2,:)==az))
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
Neuron.ia_diamond = intval + temp;

%****************************************************
%Now, get the ILD-Frequency response surface
if(if_testnum ~= 0)
   time_window = [105 205]; %ms
   reps_to_include = [];
   Neuron.if_freqaxis = ...
      DATA.Site{if_sitenum}.test{if_testnum}.params.lofreq:...
      (DATA.Site{if_sitenum}.test{if_testnum}.params.hifreq-DATA.Site{if_sitenum}.test{if_testnum}.params.lofreq)/...
      (DATA.Site{if_sitenum}.test{if_testnum}.params.numfreqs-1):DATA.Site{if_sitenum}.test{if_testnum}.params.hifreq;
   Neuron.if_ildaxis = ...
      DATA.Site{if_sitenum}.test{if_testnum}.params.loild:...
      (DATA.Site{if_sitenum}.test{if_testnum}.params.hiild-DATA.Site{if_sitenum}.test{if_testnum}.params.loild)/...
      (DATA.Site{if_sitenum}.test{if_testnum}.params.numilds-1):DATA.Site{if_sitenum}.test{if_testnum}.params.hiild;
   Neuron.if_datamatrix = DATA.Site{if_sitenum}.test{if_testnum}.datamatrix;
   [Neuron.if_mean, Neuron.if_std] = ...
      proc_test(Neuron.if_datamatrix,...
      Neuron.if_freqaxis,...
      Neuron.if_ildaxis,...
      time_window,...
      reps_to_include);
end   

%****************************************************
%Calculate d_ILD
iamr_hrtffile = get(hhrtffile,'String');
load(iamr_hrtffile);

IA_locs = Neuron.ia_locs;

[PCR_meansurf,...
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
   [2000 11000],...
   1/12,...
   1,...
   0.15,...
   0,...
   0.001,...
   0);

%Limit frequencies used
minfreq = 2000; maxfreq = 11000;
[y,minfreqind] = min(abs(HRTFinfo.hrtf_freqs - minfreq));
[y,maxfreqind] = min(abs(HRTFinfo.hrtf_freqs - maxfreq));
   
   for freq_num1 = 1:size(freq_bands,1)
      for freq_num2 = 1:size(freq_bands,1)
         temp = corrcoef(actd_curve_all{freq_num1},actd_curve_all{freq_num2});
         corrmat(freq_num1,freq_num2) = temp(1,2);
      end
   end
   
   figure
   pcolor(freq_bands(:,1),freq_bands(:,1),corrmat);
   colorbar
   
      
      [y,minfreqind] = min(abs(HRTFinfo.hrtf_freqs - min(freq_bands(:,1))));
      [y,maxfreqind] = min(abs(HRTFinfo.hrtf_freqs - max(freq_bands(:,1))));
      
      hf = figure;
      set(hf,'Units','inches',...
         'Position',[0 0 7 5]);
      
      %IF
      if(if_testnum ~= 0)
      hildf1 = subplot('Position',[0.05 vert_pos(1) 0.90 0.35]);
      set(hildf1,'FontSize',8);
      if(isfield(Neuron,'gif_freqaxis'))
         plotsurf(Neuron.gif_freqaxis,Neuron.gif_ildaxis,Neuron.gif_mean');
         title('ILD-Freq, Gammatones');
      else
         plotsurf(Neuron.if_freqaxis,Neuron.if_ildaxis,Neuron.if_mean');
         title('ILD-Freq, Tones');
      end
      if(isfield(Neuron,'ia_azi'))
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind),...
            'm',...
            'LineWidth',1.5);
      end
      if(isfield(Neuron,'ts_azi'))
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind),...
            'w',...
            'LineWidth',1.5);
      end
      xlim([2000 11000]);
      colorbar
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      end
      
      %BBIF
      hildf2 = subplot('Position',[0.05 vert_pos(2) 0.90 0.35]);
      set(hildf2,'FontSize',8);
      bb_ILDf_freqaxis = freq_bands(:,1);
      plotsurf(bb_ILDf_freqaxis,new_ildaxis2,(new_bbif_meansurf2)');
      if(isfield(Neuron,'ia_azi'))
         hold on
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,IA_bestloc_ind),...
            'm',...
            'LineWidth',1.5);
      end
      if(isfield(Neuron,'ts_azi'))
         plot(HRTFinfo.hrtf_freqs(minfreqind:maxfreqind),...
            ILD_matrix_focus(minfreqind:maxfreqind,TS_bestloc_ind),...
            'w',...
            'LineWidth',1.5);
      end
      xlim([2000 11000]);
      colorbar
      xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
      
      pervar
