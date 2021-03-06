%Engage_ILDFreq

%*******************************************************************************
%	The ILDfreq operation
%*******************************************************************************

set(H.exitildfreq,'Visible','off');
set(H.resetildfreq,'Visible','off');


%Put parameters into XStimParams
test_val = get(H.test_type,'Value');
test_type = get(H.test_type,'String');
XStimParams.test_type = deblank(test_type(test_val,:));
clear test_type
XStimParams.loild = str2num(get(H.lowild,'String'));
XStimParams.hiild = str2num(get(H.highild,'String'));
XStimParams.numilds = str2num(get(H.numilds,'String'));
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.loc_flag = get(H.locuseit,'Value');
XStimParams.ephone_flag = get(H.ephoneuseit,'Value');

lowfreq = str2num(get(H.lowfreq,'String'));
highfreq = str2num(get(H.highfreq,'String'));
numfreqs = str2num(get(H.numfreqs,'String'));
lowild = str2num(get(H.lowild,'String'));
highild = str2num(get(H.highild,'String'));
numilds = str2num(get(H.numilds,'String'));

%Setup for freqs and ilds used
freqss = lowfreq:(highfreq-lowfreq)/(numfreqs-1):highfreq;
ilds  = lowild:(highild-lowild)/(numilds-1):highild;
[dim1mat,dim2mat] = meshgrid(freqss,ilds);
FREQSS = dim1mat(:); ILDS = dim2mat(:);


numreps = str2num(get(H.numreps,'String'));
numtrials = length(freqss) * length(ilds);

clear BUF

%Specify DAMA buffers
BUF.left_1			= 1;
BUF.right_1		= 2;
%isibuffer 				= 3;
BUF.playseq_left_1			= 4;
BUF.playseq_right_1		= 5;
BUF.playspecbuf_1			= 6;
BUF.left_2			= 7;
BUF.right_2		= 8;
BUF.playseq_left_2			= 9;
BUF.playseq_right_2		= 10;
BUF.playspecbuf_2			= 11;

%Make play sequence buffers
S232('allot16',BUF.playseq_left_1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.left_1);
S232('make',1,1);
S232('make',2,0);
%S232('make',2,isibuffer);
%S232('make',2,1);
%S232('make',4,0);
S232('qpop16',BUF.playseq_left_1);
      
S232('allot16',BUF.playseq_right_1,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.right_1);
S232('make',1,1);
S232('make',2,0);
%S232('make',2,isibuffer);
%S232('make',2,1);
%S232('make',4,0);
S232('qpop16',BUF.playseq_right_1);

S232('allot16',BUF.playseq_left_2,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.left_2);
S232('make',1,1);
S232('make',2,0);
%S232('make',2,isibuffer);
%S232('make',2,1);
%S232('make',4,0);
S232('qpop16',BUF.playseq_left_2);
      
S232('allot16',BUF.playseq_right_2,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.right_2);
S232('make',1,1);
S232('make',2,0);
%S232('make',2,isibuffer);
%S232('make',2,1);
%S232('make',4,0);
S232('qpop16',BUF.playseq_right_2);


%Make play specification buffer
S232('allot16',BUF.playspecbuf_1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_left_1);
S232('make',1,BUF.playseq_right_1);
S232('make',2,0);
S232('qpop16',BUF.playspecbuf_1);

S232('allot16',BUF.playspecbuf_2,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_left_2);
S232('make',1,BUF.playseq_right_2);
S232('make',2,0);
S232('qpop16',BUF.playspecbuf_2);


%Check to see if earphone filters are being used
ephoneflag = 0;
if(exist1('H.fig') & XStimParams.ephone_flag == 1)
   ephoneflag = 1;
   ephonefname = [FN.ephone_path FN.ephone];
   ephonefilt_left  = (mtlrch(ephonefname,1))';
   ephonefilt_right = (mtlrch(ephonefname,2))';
end

%Check to see if location filters are being used
locflag = 0;
if(exist1('H.fig') & XStimParams.loc_flag == 1)
   locflag = 1;
   locfname = [FN.loc_path FN.loc];
   dirmat = sph2dbl(mtlrdir(locfname));
   locnum = find(dirmat(1,:) == XStimParams.loc_azel(2) & dirmat(2,:) == XStimParams.loc_azel(1));
   locfilt_left  = (mtlrch(locfname,2*locnum-1))';
   locfilt_right = (mtlrch(locfname,2*locnum))';
end

%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
DUR_silence = XStimParams.silence_lead; %ms
silence_len = (DUR_silence * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = str2num(get(H.DUR,'String'));
if(ephoneflag == 0)
   S232('allot16',BUF.left_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
   S232('allot16',BUF.right_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
   S232('allot16',BUF.left_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
   S232('allot16',BUF.right_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen);
elseif(ephoneflag == 1 & locflag == 0)
   S232('allot16',BUF.left_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
   S232('allot16',BUF.right_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
   S232('allot16',BUF.left_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
   S232('allot16',BUF.right_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen);
elseif(ephoneflag == 1 & locflag == 1)
   S232('allot16',BUF.left_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
   S232('allot16',BUF.right_1,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
   S232('allot16',BUF.left_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
   S232('allot16',BUF.right_2,(silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
end

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
if(exist1('H.fig') & ephoneflag == 1 & locflag == 1)
   S232('PD1npts',1,(silence_len + DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen  + ...
   TDT.ephonefiltlen + TDT.hrtffiltlen);
elseif(exist1('H.fig') & ephoneflag == 1 & locflag == 0)
   S232('PD1npts',1,(silence_len + DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen  + TDT.ephonefiltlen);
elseif(exist1('H.fig') & ephoneflag == 0)
   S232('PD1npts',1,(silence_len + DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen);
end

%Load Earphone filters
if(exist1('H.ildfreqfig') & XStimParams.loc_flag == 1 & XStimParams.ephone_flag == 1) %Use earphone and location filters
   dspid_left_loc = 0; dspid_right_loc = 1;
   dspid_left_ephone = 2; dspid_right_ephone = 3;
   S232('PD1clrsched',1);
   S232('PD1nstrms',1,2,0);
   S232('PD1resetDSP',1,hex2dec('0xFFF'));
   S232('dropall');
   %Make connections for left ear
   S232('PD1addsimp',1,S232('DSPout',dspid_left_ephone),S232('DAC',0)); %ephone to DAC1
   S232('PD1addsimp',1,S232('DSPout',dspid_left_loc),S232('DSPin',dspid_left_ephone)); %loc to ephone
   S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left_loc)); %IB to loc
   %Make connections for right ear
   S232('PD1addsimp',1,S232('DSPout',dspid_right_ephone),S232('DAC',1)); %ephone to DAC1
   S232('PD1addsimp',1,S232('DSPout',dspid_right_loc),S232('DSPin',dspid_right_ephone)); %loc to ephone
   S232('PD1specIB',1,S232('IB',1),S232('DSPin',dspid_right_loc)); %IB to loc
   
   %Load left
   S232('pushf',ephonefilt_left,length(ephonefilt_left));
   S232('PreLoadRaw',1,S232('DSPid',dspid_left_ephone),'MONO','STACK','','',TDT.ephonescale,1.0,1);
   S232('dropall');
   S232('pushf',locfilt_left,length(locfilt_left));
   S232('PreLoadRaw',1,S232('DSPid',dspid_left_loc),'MONO','STACK','','',1.0,1.0,1);
   S232('dropall');
   %Load right
   S232('pushf',ephonefilt_right,length(ephonefilt_right));
   S232('PreLoadRaw',1,S232('DSPid',dspid_right_ephone),'MONO','STACK','','',TDT.ephonescale,1.0,1);
   S232('dropall');
   S232('pushf',locfilt_right,length(locfilt_right));
   S232('PreLoadRaw',1,S232('DSPid',dspid_right_loc),'MONO','STACK','','',1.0,1.0,1);
   S232('dropall');
elseif(exist1('H.ildfreqfig') & XStimParams.loc_flag == 0 & XStimParams.ephone_flag == 1) %Use earphone filters
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
   while (exist1('H.ildfreqfig') & S232('APactive'))  pause(0); end
   S232('pushf',ephonefilt_left,length(ephonefilt_left));
   while (exist1('H.ildfreqfig') & S232('APactive'))  pause(0); end
   S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.ephonescale,1.0,1);
   %Load right
   S232('pushf',ephonefilt_right,length(ephonefilt_right));
   while (exist1('H.ildfreqfig') & S232('APactive'))  pause(0); end
   S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.ephonescale,1.0,1);
   while (exist1('H.ildfreqfig') & S232('APactive'))  pause(0); end
elseif(exist1('H.ildfreqfig') & get(H.ephoneuseit,'Value') == 0)
   S232('PD1clrsched',1);
   S232('PD1nstrms',1,2,0);
   S232('PD1specIB',1,S232('IB',0),S232('DAC',0));
   S232('PD1specIB',1,S232('IB',1),S232('DAC',1));
end


set(H.ephonefile,'Enable','off');
set(H.ephoneuseit,'Enable','off');
set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; %microsecond per sample
mii_separation = 100; %only take events separated by 100 samples (i.e., 1 ms)


ITD = str2num(get(H.ITD,'String'));
if(abs(ITD) > 250) return; end

ABL = str2num(get(H.ABL,'String'));
if(ABL < -110) return; end
S232('PA4atten',1,abs(ABL)-20);
S232('PA4atten',2,abs(ABL)-20);

ISI = str2num(get(H.ISI,'String'));
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); %correct for ITD filtlength
if(exist1('H.fig') & get(H.ephoneuseit,'Value')) %correct for ephone filtlength
   ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));
end

%Loop to make the stimuli we'll use
remreps = 1;
set(H.buildplay,'String','Building Stimuli');
set(H.remreps,'String',num2str(remreps));
repnum = 1;
buffcycle = 1;
finalspikematrix = zeros(size(dim1mat));
while (exist1('H.ildfreqfig') & (repnum <= 1))
   %Randomize the stimuli
   randseq = randperm(numtrials);
   
   %Loop to make stimuli
   trialnum = 1;
   while (exist1('H.ildfreqfig') & (trialnum <= numtrials))
      set(H.buildplay,'BackgroundColor','red');
      %Check for pause by user
      while (exist1('H.ildfreqfig') & get(H.pauseildfreq,'Value'))
         pause(0);
         if(~exist1('H.ildfreqfig')) return; end         
         set(H.exitildfreq,'Visible','on');
         set(H.resetildfreq,'Visible','on');
         if(exist1('H.ildfreqfig') & get(H.resetildfreq,'Value') == 1)
            set(H.resetildfreq,'Value',0);
            set(H.pauseildfreq,'Value',0);
            Reset_ILDfreq;
            return;
         end
      end
      %Make the stimuli
      clear tone tonestack rawtrial_left rawtrial_right trial_left trial_right
      clear itdtrial_left itdtrial_right
      switch get(H.stim_type,'Value')      
      case 1,
         %Pure Tones Get freq & ild for trial
         freq(repnum,trialnum) = FREQSS(randseq(trialnum));
         ild(repnum,trialnum) = ILDS(randseq(trialnum));
         %Make tone at specified frequency and ild - NOT GETTING ILD by ATTENUATORS
         %Lright - Lleft = 20*log10(Ampright/Ampleft), Hartmann p. 29
         mult_fact = 10^(ild(repnum,trialnum)/20);
         tone = MakeTone(TDT.Fs,freq(repnum,trialnum),DUR);
         tone = tone/max(abs(tone));
         trial_left  = tone/sqrt(mult_fact);
         trial_right = sqrt(mult_fact)*tone;
      case 2,
         %GammaTones Get freq & ild for trial
         freq(repnum,trialnum) = FREQSS(randseq(trialnum));
         ild(repnum,trialnum) = ILDS(randseq(trialnum));
         %Make tone at specified frequency and ild - NOT GETTING ILD by ATTENUATORS
         %Lright - Lleft = 20*log10(Ampright/Ampleft), Hartmann p. 29
         mult_fact = 10^(ild(repnum,trialnum)/20);
         gammatone = MakeGammaTone(TDT.Fs,freq(repnum,trialnum),DUR);
         gammatone = gammatone/max(abs(gammatone));
         trial_left  = gammatone/sqrt(mult_fact);
         trial_right = sqrt(mult_fact)*gammatone;
      case 3,
         %Use Concurrent ToneStack if specified
         freq(repnum,trialnum) = FREQSS(randseq(trialnum));
         ild(repnum,trialnum) = ILDS(randseq(trialnum));
         tonestackfreqs = freqss(find(freqss ~= freq(repnum,trialnum)));
         [trial_left,trial_right] = MakeToneStack(TDT.Fs,...
            tonestackfreqs,...
            DUR,...
            freq(repnum,trialnum),...
            ild(repnum,trialnum));
      case 4,
         %Use Concurrent GammaToneStack if specified
         freq(repnum,trialnum) = FREQSS(randseq(trialnum));
         ild(repnum,trialnum) = ILDS(randseq(trialnum));
         tonestackfreqs = freqss(find(freqss ~= freq(repnum,trialnum)));
         [trial_left,trial_right] = MakeGammaToneStack(TDT.Fs,...
            tonestackfreqs,...
            DUR,...
            freq(repnum,trialnum),...
            ild(repnum,trialnum));
      case 8,
         %Use Broadband Noise and shift ILDs of narrow frequency bands
         freq(repnum,trialnum) = FREQSS(randseq(trialnum));
         ild(repnum,trialnum) = ILDS(randseq(trialnum));
         [trial_left,trial_right] = MakeBBNoiseLR(TDT.Fs,DUR,freq(repnum,trialnum),ild(repnum,trialnum));
      otherwise,
         set(H.stim_type,'Value',2);
         freq(repnum,trialnum) = FREQSS(randseq(trialnum));
         ild(repnum,trialnum) = ILDS(randseq(trialnum));
         [trial_left,trial_right] = MakeBBNoiseLR(TDT.Fs,DUR,freq(repnum,trialnum),ild(repnum,trialnum));
        disp('Stim type not supported for ILDfreq.  Reset to BROADBAND.');    
     end
      
      %Ramp the stimuli
      ramp_time = 5; %ms
      [trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
      [trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
      
      %Add in the silent period
      trial_left =  [zeros(1,silence_len) trial_left];
      trial_right = [zeros(1,silence_len) trial_right];
      
      %Apply ITD filtering
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
         eval(['load ' FN.itd_path 'itdfilt' num2str(itdleft)]);
         eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
         eval(['load ' FN.itd_path 'itdfilt' num2str(itdright)]);
         eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
      end
      itdtrial_left = conv(trial_left,itd_filt_left);
      itdtrial_right = conv(trial_right,itd_filt_right);
      
      %Correct for filters
      if(locflag == 1 & ephoneflag == 1)
         itdtrial_left = [itdtrial_left zeros(1,TDT.ephonefiltlen) zeros(1,TDT.hrtffiltlen)];
         itdtrial_right = [itdtrial_right zeros(1,TDT.ephonefiltlen) zeros(1,TDT.hrtffiltlen)];
      elseif(locflag == 0 & ephoneflag == 1)
         itdtrial_left = [itdtrial_left zeros(1,TDT.ephonefiltlen)];
         itdtrial_right = [itdtrial_right zeros(1,TDT.ephonefiltlen)];
      end
      
      %Fill stimulus buffers
      itdtrial_left = itdtrial_left - round(mean(itdtrial_left));
      itdtrial_right = itdtrial_right - round(mean(itdtrial_right));
      ABAval = 0.5*(max(abs(itdtrial_left)) + max(abs(itdtrial_right)));
      scalefact = TDT.scalevalue/ABAval;
      itdtrial_left = round(scalefact*itdtrial_left);
      itdtrial_right = round(scalefact*itdtrial_right);
      
      if(exist1('H.ildfreqfig') & buffcycle == 1);
         S232('push16',itdtrial_left,length(itdtrial_left));
         S232('qpop16',BUF.left_1);
         fname = ['BUF.left_' num2str(FREQSS(randseq(trialnum))) '_' num2str(ILDS(randseq(trialnum)))];
         evalstr = ['S232(''dama2disk16'',BUF.left_1,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
      	 temp_left = dama2pc(BUF.left_1);
         S232('push16',itdtrial_right,length(itdtrial_right));
         S232('qpop16',BUF.right_1);
         fname = ['BUF.right_' num2str(FREQSS(randseq(trialnum))) '_' num2str(ILDS(randseq(trialnum)))];
         evalstr = ['S232(''dama2disk16'',BUF.right_1,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
         temp_right = dama2pc(BUF.right_1);
      end
      
      
      %Plot PSD
      if(exist1('H.ildfreqfig') & ~exist1('H.psdaxes'))
         figure(H.ildfreqfig);
         H.psdaxes = axes;
         set(H.psdaxes,'Visible','off');
      end
      if(exist1('H.ildfreqfig') & get(H.plotpsd,'Value') == 1)
         axes(H.psdaxes);
         set(H.psdaxes,'Visible','on',...
            'Position',[0.1 0.07 0.8 0.2]);
         axis square
         [pl,f] = psd(temp_left,2048,TDT.Fs); [pr,f] = psd(temp_right,2048,TDT.Fs);
         H.pl = plot(f,10*log10(abs(pl)),'g');
         hold on
         H.pr = plot(f,10*log10(abs(pr)),'r');
         xlabel('Frequency (Hz)');
         ILDval = (20*log10(max(abs(temp_right)))) - (20*log10(max(abs(temp_left))));
         ABLval = 0.5*((20*log10(max(abs(temp_left)))) + (20*log10(max(abs(temp_right)))));
         title(['ABL = ' num2str(ABLval) ', ILD = ' num2str(ILDS(randseq(trialnum)))]);
         set(H.psdaxes,'YLim',[0 130]);
         hold off
         grid on
         pause(1)
      elseif(exist1('H.ildfreqfig') & exist1('H.pl') & get(H.plotpsd,'Value') == 0)
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
         %[pl,f] = psd(temp_left,2048,TDT.Fs); [pr,f] = psd(temp_right,2048,TDT.Fs);
         %plot(f,10*log10(abs(pl)),'g');
         %hold on
         %plot(f,10*log10(abs(pr)),'r');
         hold off
         %keyboard
      end
      
      remtrials = numtrials - trialnum;
      set(H.remtrials,'String',num2str(remtrials));
      trialnum = trialnum + 1;
      set(H.buildplay,'BackgroundColor','yellow');
      pause(0);
   end %end loop over trials
   
   
      
   remreps = 1;
   set(H.remreps,'String',num2str(remreps));
   repnum = repnum + 1;
   pause(0);
end %end loop over reps

%Put up raster plot
position = [300 0 950 950];
xvals = freqss;
yvals = ilds;
xtext = 'Freq (Hz)';
ytext = 'ILD (dB)';
xoffset = 5;
yoffset = 5;
tot_dur = DUR_silence + DUR + 50;
stim_dur = DUR_silence + DUR;
startcount = DUR_silence;
endcount = DUR_silence + DUR;
row_space = 5;
col_space = .2*tot_dur;

[H.rasterfig,H.rastaxes] = makerasterfig(position,...
   xvals,...
   yvals,...
   xtext,...
   ytext,...
   xoffset,...
   yoffset,...
   numreps,...
   tot_dur,...
   stim_dur,...
   startcount,...
   endcount,...
   row_space,...
   col_space);

%Begin loop over experimental blocks (i.e., reps) and trials   
set(H.buildplay,'String','Playing Stimuli');
set(H.buildplay,'BackgroundColor','yellow');
set(H.remreps,'String',num2str(numreps));
repnum = 1;
buffcycle = 1;
datamatrix = [];

if(exist1('H.ildfreqfig') & get(H.recorddata,'Value'))
   XStimParams.testnum = XStimParams.testnum + 1;
end
while (exist1('H.ildfreqfig') & (repnum <= numreps))
   %Randomize the stimuli
   randseq = randperm(numtrials);   
   trialnum = 1;
   spikes_trial = [];
   freqs_trial = [];
   ilds_trial = [];
   repnum_trial = [];
   while (exist1('H.ildfreqfig') & (trialnum <= numtrials+1))
      
      %Check for pause by user
      while (exist1('H.ildfreqfig') & get(H.pauseildfreq,'Value'))
         pause(0);
         if(~exist1('H.ildfreqfig')) return; end         
         set(H.exitildfreq,'Visible','on');
         set(H.resetildfreq,'Visible','on');
         if(exist1('H.ildfreqfig') & get(H.resetildfreq,'Value') == 1)
            set(H.resetildfreq,'Value',0);
            set(H.pauseildfreq,'Value',0);
            Reset_ILDfreq;
            return;
         end
      end

      
      if(trialnum <= numtrials)
      if(buffcycle == 1);
         fname = ['BUF.left_' num2str(FREQSS(randseq(trialnum))) '_' num2str(ILDS(randseq(trialnum)))];
         evalstr = ['S232(''disk2dama16'',BUF.left_1,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
         fname = ['BUF.right_' num2str(FREQSS(randseq(trialnum))) '_' num2str(ILDS(randseq(trialnum)))];
         evalstr = ['S232(''disk2dama16'',BUF.right_1,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
      elseif(buffcycle == 2);
         fname = ['BUF.left_' num2str(FREQSS(randseq(trialnum))) '_' num2str(ILDS(randseq(trialnum)))];
         evalstr = ['S232(''disk2dama16'',BUF.left_2,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
         fname = ['BUF.right_' num2str(FREQSS(randseq(trialnum))) '_' num2str(ILDS(randseq(trialnum)))];
         evalstr = ['S232(''disk2dama16'',BUF.right_2,' ...
               [' ''' stimuli_dir fname ''' '] ...
               ',0);'];
         eval(evalstr);
      end
      end
      
      %Wait till PD1 is finished
      while S232('PD1status',1) usec_delay(1000); end
      S232('PD1stop',1);
      %Check for pause by user
      while (exist1('H.ildfreqfig') & get(H.pauseildfreq,'Value'))
         pause(0);
         if(~exist1('H.ildfreqfig')) return; end         
         set(H.exitildfreq,'Visible','on');
         set(H.resetildfreq,'Visible','on');
         if(exist1('H.ildfreqfig') & get(H.resetildfreq,'Value') == 1)
            set(H.resetildfreq,'Value',0);
            set(H.pauseildfreq,'Value',0);
            return;
         end
      end
      pause(0.05); %wait an additional 50 ms to collect spikes
      
      %Stop the m110 and get spikes
      if(trialnum > 1)
         m110dx( C_.STOP);
         spikes = m110dx( C_.DATA, 1000); %Take 100 spikes max
         ind = find(spikes ~= 0); %Get clock events that are spikes
         spikes = spikes(ind);
         ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
         if(exist1('H.ildfreqfig') & ~isempty(spikes)) 
            spikes = spikes(ind);
            spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
            freqs_trial = [freqs_trial;find(freqss == FREQSS(randseq(trialnum-1))) * ones(size(spikes))];
            ilds_trial = [ilds_trial;find(ilds == ILDS(randseq(trialnum-1))) * ones(size(spikes))];
            repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
         end 
      end
      
      if(exist1('H.ildfreqfig'))
         
         if(exist1('H.ildfreqfig') & buffcycle == 1)
            S232('seqplay',BUF.playspecbuf_1);
         elseif(exist1('H.ildfreqfig') & buffcycle == 2)
            S232('seqplay',BUF.playspecbuf_2);
         end
         S232('PD1arm',1);
         
         %Send trigger
         %Set up MII
         m100x( C_.INIT );
         m110dx( C_.INIT );
         m110dx( C_.CLOCK, mii_us_per_sample);
         m110dx( C_.MODE, M_.PST );
         
         pause(ISI/1000);
         if (trialnum <= numtrials)
            %Start clock
            m110dx( C_.START);
            %Send pulse: PD1 GO!
            m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
         end
      end
      
      if(buffcycle == 1)
         buffcycle = 2;
      else
         buffcycle = 1;
      end
      
      if(trialnum > 1)
         finalspikematrix(find(ilds == ILDS(randseq(trialnum-1))),...
            find(freqss == FREQSS(randseq(trialnum-1)))) = ...
            finalspikematrix(find(ilds == ILDS(randseq(trialnum-1))),...
            find(freqss == FREQSS(randseq(trialnum-1)))) + ...
            length(spikes);
      end
      
      if(trialnum <= numtrials)
         remtrials = numtrials - trialnum;
         set(H.remtrials,'String',num2str(remtrials));
      end
      trialnum = trialnum + 1;
      pause(0);
   end %end loop over trials
   
   %Plot all the spikes from trials
   plotraster(H.rasterfig,...
      spikes_trial,...
      freqs_trial,...
      ilds_trial,...
      repnum,...
      tot_dur,...
      stim_dur,...
      numreps,...
      xoffset,...
      yoffset,...
      row_space,...
      col_space);
   
    %Record Data
    if(exist1('H.ildfreqfig') & get(H.recorddata,'Value'))
       datamatrix = [datamatrix;[spikes_trial repnum_trial freqs_trial ilds_trial]];
       record_data(XStimParams,datamatrix);
    end
    
   remreps = numreps - repnum;
   set(H.remreps,'String',num2str(remreps));
   repnum = repnum + 1;
   pause(0);
end %end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/numreps;

H.finalspikeratefig = figure('Position',[700 20 550 500],...
   'Name','Final Test Results',...
   'NumberTitle','off');
H.temp = subplot('Position',[0.1 0.60 0.30 0.30]);
plot(freqss,max(finalspikematrix,[],1));
set(H.temp,'FontSize',8,'XLim',[min(freqss) max(freqss)]);
xlabel('Frequency (Hz)'); ylabel('Spikes/stim');
title('Frequency Tuning');
H.temp = subplot('Position',[0.65 0.60 0.30 0.30]);
plot(ilds,max(finalspikematrix,[],2));
set(H.temp,'FontSize',8,'XLim',[min(ilds) max(ilds)]);
xlabel('ILD (dB)'); ylabel('Spikes/stim');
title('ILD Tuning');
H.temp = subplot('Position',[0.1 0.1 0.80 0.35]);
set(H.temp,'FontSize',8);
plotsurf(freqss,ilds,finalspikematrix);
xlabel('Frequency (Hz)'); ylabel('ILD (dB)');
title('ILD-FREQ');
colorbar

set(H.buildplay,'String','Build/Play status');
set(H.exitildfreq,'Visible','on');
set(H.resetildfreq,'Visible','on');
