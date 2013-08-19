%Engage_altIR

%*******************************************************************************
%	The altIR Test operation
%*******************************************************************************

if(XStimParams.space_flag == 1)
    if isempty(FN.space)        % FN not yet picked
  	    [FN.space,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File');
  	    if(FN.space_path ~= 0)
   	    set(H.spacefile,'String',[FN.space_path FN.space]);
  	    end
        set(H.spacefilepb,'Value',0);
    end
	hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space]));
elseif(XStimParams.ildalone_flag == 1)
    if isempty(FN.ildalone)        % FN not yet picked
        [FN.ildalone,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
	    if(FN.ILA_path ~= 0)
   	        set(H.ildalonefile,'String',[FN.ILA_path FN.ildalone]);
        end
        set(H.ildalonepb,'Value',0);
    end
	hrtfdirmat = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone]));
elseif(XStimParams.itdalone_flag == 1)
    if isempty(FN.itdalone)        % FN not yet picked
        [FN.itdalone,FN.ILA_path] = uigetfile('*.*','Select ITD-alone HRTF File');
	    if(FN.ITA_path ~= 0)
   	        set(H.itdalonefile,'String',[FN.ITA_path FN.itdalone]);
        end
        set(H.itdalonepb,'Value',0);
    end
	hrtfdirmat = sph2dbl(mtlrdir([FN.ITA_path FN.itdalone]));
    if isempty(FN.ildalone)        % FN not yet picked
        [FN.ildalone,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
	    if(FN.ILA_path ~= 0)
   	        set(H.ildalonefile,'String',[FN.ILA_path FN.ildalone]);
        end
        set(H.ildalonepb,'Value',0);
    end
end

%Put parameters into XStimParams
test_val = get(H.test_type,'Value');
test_type = get(H.test_type,'String');
XStimParams.test_type = deblank(test_type(test_val,:));
clear test_type
if get(H.stim_type,'Value') ~= 9			% not from file
	XStimParams.bandwidth = str2num(get(H.bandwidth,'String'));
	XStimParams.bandlimits(1) = str2num(get(H.bandlimits1,'String'));
	XStimParams.bandlimits(2) = str2num(get(H.bandlimits2,'String'));
else
    if isempty(FN.stim)         %Stimulus from file
     [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
    end
end
XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));

%%%%
clear locind

%Specify DAMA buffers
clear BUF
BUF.L1				= 1;
BUF.R1				= 2;
BUF.L2				= 3;
BUF.R2				= 4;
BUF.playseq_L1		= 5;
BUF.playseq_R1		= 6;
BUF.playspec1		= 7;
BUF.playseq_L2		= 8;
BUF.playseq_R2		= 9;
BUF.playspec2		= 10;

%Make play sequence buffers
S232('allot16',BUF.playseq_L1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.L1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_L1);
      
S232('allot16',BUF.playseq_R1,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.R1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_R1);

S232('allot16',BUF.playseq_L2,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.L2);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_L2);
      
S232('allot16',BUF.playseq_R2,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.R2);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_R2);

%Make play specification buffer
S232('allot16',BUF.playspec1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L1);
S232('make',1,BUF.playseq_R1);
S232('make',2,0);
S232('qpop16',BUF.playspec1);

S232('allot16',BUF.playspec2,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L2);
S232('make',1,BUF.playseq_R2);
S232('make',2,0);
S232('qpop16',BUF.playspec2);


%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
DUR_silence = 100; 					%ms
silence_len = (DUR_silence * round(TDT.Fs/1000));
%Add a piece of silence after stimulus 
DUR_silence2 = 50; 					%ms
silence_len2 = (DUR_silence2 * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = XStimParams.curr_stimdur;
if(XStimParams.space_flag == 0) 					%ILDAlone Test
   S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
   S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
   S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
   S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
elseif(XStimParams.space_flag == 1) 			%fully-cued Test
   S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
   S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
   S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
   S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
end   


S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
if(XStimParams.space_flag == 0) 					%ILDAlone Test
   PD1npts = (silence_len2 + silence_len + DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen  + TDT.ephonefiltlen + TDT.hrtffiltlen;
elseif(XStimParams.space_flag == 1) 			%fully-cued Test
    PD1npts = (silence_len2 + silence_len + DUR*(round(TDT.Fs/1000))) + TDT.ephonefiltlen + TDT.hrtffiltlen;
end
S232('PD1npts',1,PD1npts);

%Get Earphone filters
if isempty(FN.ephone)
    EphoneFileDiagBox;
end
ephonefname = [FN.ephone_path FN.ephone];
ephonefilt_left  = (mtlrch(ephonefname,1))';
ephonefilt_right = (mtlrch(ephonefname,2))';

%Load Earphone filters
dspid_left = 0; dspid_right = 1;
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('0xFFF'));
S232('dropall');
%Make connections for left ear
S232('PD1addsimp',1,S232('DSPout',dspid_left),S232('DAC',0)); %DSPout to DAC0
S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left)); %IB to DSPin
%Make connections for right ear
S232('PD1addsimp',1,S232('DSPout',dspid_right),S232('DAC',1));
S232('PD1specIB',1,S232('IB',1),S232('DSPin',dspid_right));
%Load left      
S232('pushf',ephonefilt_left,length(ephonefilt_left));
S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.ephonescale,1.0,1);
%Load right
S232('pushf',ephonefilt_right,length(ephonefilt_right));
S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.ephonescale,1.0,1);

%set(H.ephonefile,'Enable','off');
%set(H.ephoneuseit,'Enable','off');
set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

ITD = XStimParams.curr_ITD;
if(abs(ITD) > 250) return; end

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

ISI = XStimParams.test_ISI;
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); 				%correct for ITD filtlength
ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));				%correct for ephonefilt

%Get the HRTF spectra indices
if ~get(H.pickerflag,'value')
	locscriptfid = fopen([FN_script_path FN.script],'r');
	XStimParams.locations = fscanf(locscriptfid,'%i %i',[2 inf]);
	fclose(locscriptfid);
else
   XStimParams.locations = GUI.locations1';
   if ~length(XStimParams.locations)
   	set(H.pickerfig,'Color', [.1 .2 .8]);
%   	set(H.picker_error,'visible','on');
   	pause;
     	XStimParams.locations = GUI.locations1';
 		set(H.picker_error,'visible','off');
   	set(H.pickerfig,'Color', [.8 .8 .8]);
   end
end
for locnum = 1:size(XStimParams.locations,2)
	locind(locnum) = max1(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
	 	hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
end

%%%%% Loop to make the stimuli we'll use
remreps = 1;
set(H.buildplay,'String','Building Stimuli');
set(H.remreps,'String',num2str(remreps));
repnum = 1;
numtrials = size(XStimParams.locations,2);
finalspikematrix = zeros(1,numtrials);

while (exist1('H.altIRfig') & (repnum <= 1))
   %Randomize the stimuli
   randseq = randperm(numtrials);
   
   trialnum = 1;
   while (exist1('H.altIRfig') & (trialnum <= numtrials))
      set(H.buildplay,'BackgroundColor','red');
      %Check for pause by user
      while (exist1('H.altIRfig') & get(H.pausealtIR,'Value'))
         pause(0);
         if(~exist1('H.altIRfig')) return; end         
         set(H.exitaltIR,'Visible','on');
         set(H.resetaltIR,'Visible','on');
         if(exist1('H.altIRfig') & get(H.resetaltIR,'Value') == 1)
            set(H.resetaltIR,'Value',0);
            set(H.pausealtIR,'Value',0);
            Reset_altIR;
            return;
         end
         if isempty(XStimParams.locations)
             Reset_altIR;
             return;
         end
      end
      %Make the stimuli
      switch get(H.stim_type,'Value')      
      case 3,
         %Use Concurrent ToneStack if specified: 10 frequencies
         tonestackfreqs = ...
            XStimParams.bandwidth(1):...
            round((XStimParams.bandwidth(2) - XStimParams.bandwidth(1))/9):...
            XStimParams.bandwidth(2);
         [trial_left,trial_right] = MakeToneStack(TDT.Fs,tonestackfreqs,XStimParams.curr_stimdur);
      case 6,
         %Narrowband Noise
         [trial_left] = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
         trial_right = trial_left;
      case 7,
         %Narrowband Noise Stack
         [trial_left] = MakeNBNoiseStack(TDT.Fs,XStimParams.curr_stimdur,XStimParams.spec_NBNoiseStackBands);
         trial_right = trial_left;
      case 8,
         %Broadband Noise
         [trial_left] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
         trial_right = trial_left;
      case 9,
         %Stimulus from file
         if(~exist1('stim_from_file'))
            fid = fopen([FN.stim_path FN.stim],'r');
            stim_from_file = fread(fid,inf,'float');
            fclose(fid);
         end
         [trial_left] = stim_from_file;
         if(size(trial_left,1) > 1)
            trial_left = trial_left';
         end
         trial_right = trial_left;
      otherwise,
         set(H.stim_type,'Value',8);
         [trial_left] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
         trial_right = trial_left;
         disp('Stim type not supported for altIR.  Reset to BROADBAND');
      end
      
      if get(H.stim_type,'Value') ~= 9
      	%Ramp the stimuli
      	ramp_time = 5; %ms
      	[trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
      	[trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
      end
      
      %%%%%%% removed zero padding from here to lower down  8/29/02
      
      %Apply ITD filtering if conducting ILDAlone Space Test
      if(XStimParams.space_flag == 0)
         itdleft = 0; itdright = 0;
         ITD = round(str2num(get(H.ITD,'String')));
         if(ITD < 0)
            itdleft = 0;
            itdright = abs(ITD);
         elseif(ITD > 0)
            itdleft = abs(ITD);
            itdright = 0;
         end
         if(trialnum == 1)
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdleft)]);
            eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
            eval(['load ' FN.ITD_path 'itdfilt' num2str(itdright)]);
            eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
         end
         trial_left = conv(trial_left,itd_filt_left);
         trial_right = conv(trial_right,itd_filt_right);
      end
      
    
      %scale stimuli
      ABAval = 0.5*(max(abs(trial_left)) + max(abs(trial_right)));
      scalefact = TDT.scalevalue/ABAval;
      temp_left = round(scalefact*trial_left);
      temp_right = round(scalefact*trial_right);
      
      %Apply HRTF filtering
      if(XStimParams.space_flag == 1)
         hrtf_left = mtlrch([FN.space_path FN.space],(2*locind(randseq(trialnum)))-1);
         hrtf_right = mtlrch([FN.space_path FN.space],2*locind(randseq(trialnum)));
         % alter mag of IRs
         if XStimParams.alter_left
             if XStimParams.alterMag
                 hrtf_leftX = altmagIR3 (hrtf_left, XStimParams.alterMag_Minbandlimits, XStimParams.alterMag_Maxbandlimits, XStimParams.alterMag_dB);
             end
             if XStimParams.alterPhase
                 hrtf_leftX = altphaseIR3 (hrtf_left, XStimParams.alterPhase_Minbandlimits, XStimParams.alterPhase_Maxbandlimits, XStimParams.alterPhase_rad);
             end
         end   
         % alter phase of IRs     
         if XStimParams.alter_right
             if XStimParams.alterMag
                 hrtf_rightX = altmagIR3 (hrtf_right, XStimParams.alterMag_Minbandlimits, XStimParams.alterMag_Maxbandlimits, XStimParams.alterMag_dB);
             end
             if XStimParams.alterPhase
                 hrtf_rightX = altphaseIR3 (hrtf_right, XStimParams.alterPhase_Minbandlimits, XStimParams.alterPhase_Maxbandlimits, XStimParams.alterPhase_rad);
             end
         end   

         trial_left = conv(temp_left,hrtf_left);
         trial_right = conv(temp_right,hrtf_right);
        if XStimParams.alter_left
         trial_leftX = conv(temp_left,hrtf_leftX);              %%%% added 8/30/02
         ScaleFactor = mom(trial_left,2) / mom(trial_leftX,2);
         trial_left = trial_leftX * ScaleFactor;
        end             
        if XStimParams.alter_right   
         trial_rightX = conv(temp_right,hrtf_rightX);              %%%% added 8/30/02
         ScaleFactor = mom(trial_right,2) / mom(trial_rightX,2);
         trial_right = trial_rightX * ScaleFactor;
        end
      elseif(XStimParams.ildalone_flag == 1)
         hrtf_left = mtlrch([FN.ILA_path FN.ildalone],(2*locind(randseq(trialnum)))-1);
         hrtf_right = mtlrch([FN.ILA_path FN.ildalone],2*locind(randseq(trialnum)));
         
         % alter mag of IRs
         if XStimParams.alter_left
             if XStimParams.alterMag
                 hrtf_leftX = altmagIR3 (hrtf_left, XStimParams.alterMag_Minbandlimits, XStimParams.alterMag_Maxbandlimits, XStimParams.alterMag_dB);
             end
             if XStimParams.alterPhase
                 hrtf_leftX = altphaseIR3 (hrtf_left, XStimParams.alterPhase_Minbandlimits, XStimParams.alterPhase_Maxbandlimits, XStimParams.alterPhase_rad);
             end
         end   
         % alter phase of IRs    
         if XStimParams.alter_right
             if XStimParams.alterMag
                 hrtf_rightX = altmagIR3 (hrtf_right, XStimParams.alterMag_Minbandlimits, XStimParams.alterMag_Maxbandlimits, XStimParams.alterMag_dB);
             end
             if XStimParams.alterPhase
                 hrtf_rightX = altphaseIR3 (hrtf_right, XStimParams.alterPhase_Minbandlimits, XStimParams.alterPhase_Maxbandlimits, XStimParams.alterPhase_rad);
             end
         end   
            
         trial_left = conv(temp_left,hrtf_left);
         trial_right = conv(temp_right,hrtf_right);
         if XStimParams.alter_left
         trial_leftX = conv(temp_left,hrtf_leftX);              %%%% added 8/30/02
         ScaleFactor = mom(trial_left,2) / mom(trial_leftX,2);
         trial_left = trial_leftX * ScaleFactor;
        end             
        if XStimParams.alter_right   
         trial_rightX = conv(temp_right,hrtf_rightX);              %%%% added 8/30/02
         ScaleFactor = mom(trial_right,2) / mom(trial_rightX,2);
         trial_right = trial_rightX * ScaleFactor;
        end
     end
      
      % remove DC offset
      trial_left = trial_left - round(mean(trial_left));
      trial_right = trial_right - round(mean(trial_right));

%%%%% added padding here from above 8/29/02      
       %Add in the leading silent period
      trial_left =  [zeros(1,silence_len) trial_left];
      trial_right = [zeros(1,silence_len) trial_right];
     
      %Add in the trailing silent period
      trial_left =  [trial_left zeros(1,silence_len2)];
      trial_right = [trial_right zeros(1,silence_len2)];
%%%%%      

     %pad with zeros for ephonefilters
      filttrial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
      filttrial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
      
      if length(filttrial_left) ~= length(filttrial_right)
          disp('They do not match ?!!!')
      end
      
      % save stimuli to disk
      if(exist1('H.altIRfig'));
         S232('push16',filttrial_left,length(filttrial_left));
         S232('qpop16',BUF.L1);
         fname = ['stimbuf_left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
               '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
         evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
      	temp_left = dama2pc(BUF.L1);
         S232('push16',filttrial_right,length(filttrial_right));
         S232('qpop16',BUF.R1);
         fname = ['stimbuf_right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
               '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
         evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
         temp_right = dama2pc(BUF.R1);
      end
      
      
      %Plot PSD
      if(exist1('H.altIRfig') & ~exist1('H.psdaxes'))
         figure(H.altIRfig);
         H.psdaxes = axes;
         set(H.psdaxes,'Visible','off');
      end
      if(exist1('H.altIRfig') & get(H.plotpsd,'Value') == 1)
         axes(H.psdaxes);
         set(H.psdaxes,'Visible','on',...
            'Position',[0.1 0.04 0.8 0.18]);
         axis square
         [pl,f] = psd(temp_left,2048,TDT.Fs); [pr,f] = psd(temp_right,2048,TDT.Fs);
         H.pl = plot(f,10*log10(abs(pl)),'g');
         hold on
         H.pr = plot(f,10*log10(abs(pr)),'r');
         xlabel('Frequency (Hz)');
         ILDval = (20*log10(max(abs(temp_left)))) - (20*log10(max(abs(temp_right))));
         ABLval = 0.5*((20*log10(max(abs(temp_left)))) + (20*log10(max(abs(temp_right)))));
         title(['ABL = ' num2str(ABLval) ', ILD = ' num2str(ILDval)]);
         set(H.psdaxes,'YLim',[0 130]);
         hold off
         grid on
         pause(1);
      elseif(exist1('H.altIRfig') & exist1('H.pl') & get(H.plotpsd,'Value') == 0)
         set(H.psdaxes,'Visible','off');
         set(H.pl,'Visible','off');
         set(H.pr,'Visible','off');
      end
      
      %Debug
      debug = 0;
      if(~exist1('H.debugaxes') & debug)
         figure
         H.debugaxes = axes;
      end
      if(debug)
         axes(H.debugaxes)
         plot(temp_left,'g'); hold on; plot(temp_right,'r');
         hold off
      end
      
      remtrials = numtrials - trialnum;
      set(H.remtrials,'String',num2str(remtrials));
      trialnum = trialnum + 1;
      set(H.buildplay,'BackgroundColor','yellow');
      pause(0);
   end 										%end loop over trials
      
   remreps = 1;
   set(H.remreps,'String',num2str(remreps));
   repnum = repnum + 1;
   pause(0);
end 											%end loop over reps
%%% finished making stimuli

%Begin playing sounds   
set(H.buildplay,'String','Playing Stimuli');
set(H.buildplay,'BackgroundColor','yellow');
set(H.remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
buffcycle = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.altIRfig') & get(H.recorddata,'Value'))
   update_dataFN;
end

%loop for reps
while (exist1('H.altIRfig') & (repnum <= XStimParams.numreps))
   %Randomize the stimuli
   randseq = randperm(numtrials);   
   trialnum = 1;
   spikes_trial = [];
   EL_trial = [];
   AZ_trial = [];
   repnum_trial = [];
   Nspikes = [];
   
   % loop for trials
   while (exist1('H.altIRfig') & (trialnum <= numtrials))

      %Check for pause by user
      while (exist1('H.altIRfig') & get(H.pausealtIR,'Value'))
         pause(0);
         if(~exist1('H.altIRfig')) return; end         
         set(H.exitaltIR,'Visible','on');
         set(H.resetaltIR,'Visible','on');
         if(exist1('H.altIRfig') & get(H.resetaltIR,'Value') == 1)
            set(H.resetaltIR,'Value',0);
            set(H.pausealtIR,'Value',0);
            Reset_altIR;
            return;
         end
         if isempty(XStimParams.locations)
             Reset_altIR;
             return;
         end
    end
      
      %Wait till PD1 is finished
      while S232('PD1status',1) usec_delay(1000); end

  if(trialnum <= numtrials)
         if(buffcycle == 1);				% load first buffers
         fname = ['stimbuf_left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
               '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
         evalstr = ['S232(''disk2dama16'',BUF.L1,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
         fname = ['stimbuf_right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
               '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
         evalstr = ['S232(''disk2dama16'',BUF.R1,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
      elseif(buffcycle == 2);				% load second buffers
         fname = ['stimbuf_left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
               '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
         evalstr = ['S232(''disk2dama16'',BUF.L2,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
         fname = ['stimbuf_right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
               '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
         evalstr = ['S232(''disk2dama16'',BUF.R2,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
      end
  end
   
      %Wait till PD1 is finished
      while S232('PD1status',1) usec_delay(1000); end

      S232('PD1stop',1);
      %keyboard
      %Check for pause by user
      while (exist1('H.altIRfig') & get(H.pausealtIR,'Value'))
         pause(0);
         if(~exist1('H.altIRfig')) return; end         
         set(H.exitaltIR,'Visible','on');
         set(H.resetaltIR,'Visible','on');
         if(exist1('H.altIRfig') & get(H.resetaltIR,'Value') == 1)
            set(H.resetaltIR,'Value',0);
            set(H.pausealtIR,'Value',0);
            return;
         end
         if isempty(XStimParams.locations)
             Reset_altIR;
             return;
         end
      end
      
      %Stop the m110 and get spikes
      if(trialnum > 1)					% first trial just for loading sounds (sounds played after this)
         m110dx( C_.STOP);
         spikes = m110dx( C_.DATA, 1000); 			% Take 100 spikes max
         ind = find(spikes ~= 0); 						% Get clock events that are spikes
         spikes = spikes(ind);
         ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
         if(exist1('H.altIRfig') & ~isempty(spikes)) 
            spikes = spikes(ind);
            spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
            EL_trial = [EL_trial;hrtfdirmat(1,locind(randseq(trialnum-1)))* ones(size(spikes))];
            AZ_trial = [AZ_trial;hrtfdirmat(2,locind(randseq(trialnum-1)))* ones(size(spikes))];
            repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
            Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
         end
      end
      
      if(exist1('H.altIRfig'))
         
         if(exist1('H.altIRfig') & buffcycle == 1)
            S232('seqplay',BUF.playspec1);
         elseif(exist1('H.altIRfig') & buffcycle == 2)
            S232('seqplay',BUF.playspec2);
         end
         S232('PD1arm',1);
         
         %Send trigger
         %Set up MII
         m100x( C_.INIT );
         m110dx( C_.INIT );
         m110dx( C_.CLOCK, mii_us_per_sample);
         m110dx( C_.MODE, M_.PST );
         
         pause(ISI/1000);
         %%%%%
         if (trialnum <= numtrials)
            %Start clock
            m110dx( C_.START);
            %Send pulse: PD1 GO!
            m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
         end
         %Check for pause by user
         while (exist1('H.altIRfig') & get(H.pausealtIR,'Value'))
            pause(0);
            if(~exist1('H.altIRfig')) return; end         
            set(H.exitaltIR,'Visible','on');
            set(H.resetaltIR,'Visible','on');
            if(exist1('H.altIRfig') & get(H.resetaltIR,'Value') == 1)
               set(H.resetaltIR,'Value',0);
               set(H.pausealtIR,'Value',0);
               Reset_altIR;
               return;
            end
         end
         if isempty(XStimParams.locations)
             Reset_altIR;
             return;
         end
      end
      
      if(buffcycle == 1)
         buffcycle = 2;
      else
        buffcycle = 1;
      end
      
      if(trialnum > 1)
         finalspikematrix(randseq(trialnum-1)) = ...
            finalspikematrix(randseq(trialnum-1)) + ...
            length(spikes);
      end
      
      remtrials = numtrials - trialnum;
      set(H.remtrials,'String',num2str(remtrials));
      trialnum = trialnum + 1;
      pause(0);
   end %end loop over trials
   
   %Plot Spike Rate Data
   interimspikerate = finalspikematrix/repnum;
   if(exist1('H.altIRfig') & ~exist1('H.finalspikeratefig'))
      H.finalspikeratefig = figure('Position',[700 20 550 500],...
         'Name','altIR Test Spike Rate Plot',...
         'NumberTitle','off');
      H.spikeaxes = axes;
   end
   figure(H.finalspikeratefig)
%   [azi,ele,interimdiam] = array2diamond(interimspikerate,XStimParams.locations);
%   plotdiam(azi,ele,interimdiam);
   plotdiam1(XStimParams.locations, interimspikerate);
   set(H.spikeaxes,'Color','black');
   xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
   colorbar
   
   %Record Data
    if(exist1('H.altIRfig') & get(H.recorddata,'Value'))
       datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
       record_data3(XStimParams,datamatrix);
    end
    
   remreps = XStimParams.numreps - repnum;
   set(H.remreps,'String',num2str(remreps));
   repnum = repnum + 1;
   pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.finalspikeratefig)
set(H.finalspikeratefig,'Name','Final Plot for altIR Test');
%[azi,ele,finaldiam] = array2diamond(finalspikematrix,XStimParams.locations);
%plotdiam(azi,ele,finaldiam);
plotdiam1(XStimParams.locations, interimspikerate);
set(H.spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
      ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))],...
   'FontSize',8);
colorbar

set(H.buildplay,'String','Build/Play status');
set(H.exitaltIR,'Visible','on');
set(H.resetaltIR,'Visible','on');

% increment test number
if(exist1('H.altIRfig') & get(H.recorddata,'Value'))
	XStimParams.testnum = XStimParams.testnum +1;
	set(H.testnum, 'String', num2str(XStimParams.testnum));
	update_dataFN;
end
