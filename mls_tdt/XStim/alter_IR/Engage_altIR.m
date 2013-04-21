function [] = Engage_altIR()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_altIR
% VERSION to use 1/12th octave bands instead of 1/3rd
% debugging pflag
pflag = 0;
%*******************************************************************************
%	The altIR Test operation
%*******************************************************************************

%%% re. altmag and altphase, will need to get rid of min and max bands
%%% and implement find(XStimParams.Mag_octlistBIN ~=0)
%%% and find(XStimParams.Phase_octlistBIN ~=0)

stimuli_dir = FN.temp_stim_path;
fclose all;
if ~XStimParams.buildOnly     % buildOnly keeps adding files to folder (clear these from Xstim menu)
    eval(['delete ' stimuli_dir '*.*;']);
end

if(XStimParams.space_flag == 1)
    if isempty(FN.space_std) | ~isempty(strfind(FN.space_std,'eq'))        % FN not yet picked
        [FN.space_std,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.std preferred)');
        if(FN.space_path ~= 0)
            set(H.altIR_spacefile,'String',[FN.space_path FN.space_std]);
        end
        set(H.altIR_spacefilepb,'Value',0);
        FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
    end
    %%%%%%%
    if FN.HRTFfiletype(1,2) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space_std]));
    elseif FN.HRTFfiletype(1,2) == 2
        dir = 0;
        eval(['load -mat ' FN.space_path FN.space_std]);
        TF1_space = TF1;    TF2_space = TF2;
        hrtfdirmat = dir;
        clear dir TF1 TF2
    else
        disp(['space HRTFfiletype incorrect'])
        return
    end
    str1 = 'fc';
elseif(XStimParams.ildalone_flag == 1)    
    if isempty(FN.ildalone_std)        % FN not yet picked
        [FN.ildalone_std,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.altIR_ildalonefile,'String',[FN.ILA_path FN.ildalone_std]);
        end
        set(H.altIR_ildalonepb,'Value',0);
        FN.HRTFfiletype(2,2) = testHRTFfiletype(FN.ILA_path, FN.ildalone_std);
    end
    if FN.HRTFfiletype(2,2) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone_std]));
    elseif FN.HRTFfiletype(2,2) == 2
        dir = 0;
        eval(['load -mat ' FN.ILA_path FN.ildalone_std]);
        TF1_ila = TF1;    TF2_ila = TF2;
        hrtfdirmat = dir;
        clear dir TF1 TF2
    else
        disp(['ildalone HRTFfiletype incorrect'])
        return
    end
    str1 = 'ila';
elseif(XStimParams.itdalone_flag == 1)
    if isempty(FN.itdalone_std)        % FN not yet picked
        [FN.itdalone_std,FN.ITA_path] = uigetfile('*.*','Select ITD-alone HRTF File');
        if(FN.ITA_path ~= 0)
            set(H.altIR_itdalonefile,'String',[FN.ITA_path FN.itdalone_std]);
        end
        set(H.altIR_itdalonepb,'Value',0);
        FN.HRTFfiletype(3,2) = testHRTFfiletype(FN.ITA_path, FN.itdalone_std);
    end
    if FN.HRTFfiletype(3,2) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ITA_path FN.itdalone_std]));
    elseif FN.HRTFfiletype(3,2) == 2
        dir = 0;
        eval(['load -mat ' FN.ITA_path FN.itdalone_std]);
        TF1_ita = TF1;    TF2_ita = TF2;
        hrtfdirmat = dir;
        clear dir TF1 TF2
    else
        disp(['itdalone HRTFfiletype incorrect'])
        return
    end
    
    % also need FN.ildalone for this test
    if isempty(FN.ildalone_std)        % FN not yet picked
        [FN.ildalone_std,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.altIR_ildalonefile,'String',[FN.ILA_path FN.ildalone_std]);
        end
        set(H.altIR_ildalonepb,'Value',0);
        FN.HRTFfiletype(2,2) = testHRTFfiletype(FN.ILA_path, FN.ildalone_std);
    end
    %%%%%%%
    if FN.HRTFfiletype(2,2) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone_std]));
    elseif FN.HRTFfiletype(2,2) == 2
        dir = 0;
        eval(['load -mat ' FN.ILA_path FN.ildalone_std]);
        TF1_ila = TF1;    TF2_ila = TF2;
        dir_ILA = dir;
        clear dir TF1 TF2
    else
        disp(['ildalone HRTFfiletype incorrect'])
        return
    end
    str1 = 'ita';
end

if ~XStimParams.alter_setbands
    if XStimParams.alterMag
        if isempty(find(XStimParams.Mag_octlistBIN ~=0))
            for i = 1:28
                set(H.altIR_Mag_octlist(i),'backgroundcolor','m');
            end
            disp('Please check frequencies');
            pause
            for i = 1:28
                set(H.altIR_Mag_octlist(i),'backgroundcolor',[.8 .8 .8]);
            end
            return
        end
    end
    if XStimParams.alterPhase
        if isempty(find(XStimParams.Phase_octlistBIN ~=0))
            for i = 1:28
                set(H.altIR_Phase_octlist(i),'backgroundcolor','m');
            end
            disp('Please check frequencies');
            pause
            for i = 1:28
                set(H.altIR_Phase_octlist(i),'backgroundcolor',[.8 .8 .8]);
            end
            return
        end
    end
end

%Put parameters into XStimParams
XStimParams.test_type = 'AltIR';
clear test_type
if get(H.stim_type,'Value') ~= 9			% not from file
    XStimParams.bandwidth = str2num(get(H.altIR_bandwidth,'String'));
    XStimParams.bandlimits(1) = str2num(get(H.altIR_bandlimits1,'String'));
    XStimParams.bandlimits(2) = str2num(get(H.altIR_bandlimits2,'String'));
else
    if isempty(FN.stim)         %Stimulus from file
        [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
    end
end
XStimParams.curr_ITD = str2num(get(H.altIR_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.altIR_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.altIR_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.altIR_ISI,'String'));
XStimParams.numreps = str2num(get(H.altIR_numreps,'String'));
XStimParams.reset_flag = 0;

% write diary
if get(H.altIR_recorddata,'Value')
    str2 = [];
    if XStimParams.alter_left;        str2 = [str2 '  left: '];    end
    if XStimParams.alter_right        str2 = [str2 '  right: '];    end
    if XStimParams.alterPhase        str2 = [str2 '  ' num2str(XStimParams.alterPhase_rad) ' phase'];    end
    if XStimParams.alterMag        str2 = [str2 '  ' num2str(XStimParams.alterMag_dB) ' mag'];    end
    tempstr = ['    ' str1 ' altIR-test: ' num2str(size(XStimParams.locations,2)) ' locations;     ' str2 '  atten: ' num2str(abs(XStimParams.curr_ABL))];
    update_diary
end

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
DUR_silence = XStimParams.silence_lead; 					%ms
silence_len = (DUR_silence * round(TDT.Fs/1000));
%Add a piece of silence after stimulus 
DUR_silence2 = XStimParams.silence_trail; 					%ms
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
if FN.HRTFfiletype(5) == 1
    ephonefname = [FN.ephone_path FN.ephone];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
else
    eval(['load -mat ' FN.ephone_path FN.ephone]);
    ephonefilt_left  = TF1;
    ephonefilt_right = TF2;
    clear TF1 TF2 dir
end

%Load Earphone filters
dspid_left = 0; dspid_right = 1;
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));
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
XStimParams.locations = GUI.locations1';
if ~length(XStimParams.locations)
    set(H.pickerfig,'Color', [.1 .2 .8]);
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end
for locnum = 1:size(XStimParams.locations,2)
    locind(locnum) = max1(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
end

% save XStimParams for this test
Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_AltIR_' str1 ' Temp_params;'])
clear Temp_params str1

%%%%% Loop to make the stimuli we'll use
% saving (or using saved) altIRs or hrtfs asif needed
if get(H.altIR_saveHRTFflag,'Value')
    if exist(['' char(FN.HRTF_path) char(FN.saveHRTF) '.mat' ''],'file')
        delete(['' char(FN.HRTF_path) char(FN.saveHRTF) '.mat' ''])
    end
    alter_left = XStimParams.alter_left;    alter_right = XStimParams.alter_right;
    alterPhase = XStimParams.alterPhase;    alterPhase_rad = XStimParams.alterPhase_rad;
    alterMag = XStimParams.alterMag;    alterMag_dB = XStimParams.alterMag_dB;
    Mag_octlistBIN = XStimParams.Mag_octlistBIN;    Phase_octlistBIN = XStimParams.Phase_octlistBIN;
    
    locations = XStimParams.locations;
    eval(['save ' char(FN.HRTF_path) char(FN.saveHRTF) ' hrtfdirmat locind alter* '...
            'locations Phase_octlistBIN Mag_octlistBIN;']);
end        
% use saved altered HRTFs
if get(H.altIR_use_saveHRTFflag,'Value') |get(H.altIR_use_saveALTflag,'Value')
    eval(['load ' char(FN.HRTF_path) char(FN.saveHRTF) ' hrtfdirmat locind alter* '...
            'locations Phase_octlistBIN Mag_octlistBIN;']);
    XStimParams.alter_left = alter_left;    XStimParams.alter_right = alter_right;
    XStimParams.alterPhase =  alterPhase;    XStimParams.alterPhase_rad = alterPhase_rad;
    XStimParams.alterMag = alterMag;    XStimParams.alterMag_dB = alterMag_dB;
    XStimParams.Mag_octlistBIN = Mag_octlistBIN;    XStimParams.Phase_octlistBIN = Phase_octlistBIN;
    
    for i = 1:length(Mag_octlistBIN)-1
        set(H.altIR_Mag_octlist(i),'value', Mag_octlistBIN(i));
        set(H.altIR_Phase_octlist(i),'value', Phase_octlistBIN(i));
    end
    
    XStimParams.locations = locations;
    GUI.locations1 = locations';
    clear alter* locations Mag_octlistBIN Phase_octlistBIN
    refresh_spaceaxis
end        

remreps = 1;
set(H.altIR_buildplay,'String','Building Stimuli');
set(H.altIR_remreps,'String',num2str(remreps));
repnum = 1;
numtrials = size(XStimParams.locations,2);
finalspikematrix = zeros(1,numtrials);

while (exist1('H.altIRfig') & (repnum <= 1))
    %Randomize the stimuli
    if ~get(H.altIR_use_saveHRTFflag,'Value')
        randseq = randperm(numtrials);
    end
    
    trialnum = 1;
    while (exist1('H.altIRfig') & (trialnum <= numtrials))
        set(H.altIR_buildplay,'BackgroundColor','red');
        %Check for pause by user
        if pause_check; return; end
        clear ALTout* hrtf_*
        % get az and el for this trial for naming purposes below:
        el = hrtfdirmat(1,locind(randseq(trialnum)));
        az = hrtfdirmat(2,locind(randseq(trialnum)));
        if el<0
            el = ['n' num2str(abs(el))];
        else
            el = ['p' num2str(el)];
        end
        if az<0
            az = ['n' num2str(abs(az))];
        else
            az = ['p' num2str(az)];
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
        
        % remove DC offset (added here on 3/17/03)
        trial_left = trial_left - mom(trial_left,1);
        trial_right = trial_right - mom(trial_right,1);
        
        % set ACpower==1
        trial_left = trial_left /mom(trial_left,2);
        trial_right = trial_right /mom(trial_right,2);
        
        if get(H.stim_type,'Value') ~= 9        %Ramp the stimuli
            ramp_time = 5; %ms
            [trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
            [trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
        end
        
        % remove DCoffset
        trial_left = trial_left - mom(trial_left,1);
        trial_right = trial_right - mom(trial_right,1);
        
        %Apply ITD filtering if conducting ILDAlone Space Test
        if(XStimParams.space_flag == 0)
            itdleft = 0; itdright = 0;
            ITD = round(str2num(get(H.altIR_ITD,'String')));
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
        
        %Apply HRTF filtering   
        if(XStimParams.space_flag == 1)                 % FC
            if FN.HRTFfiletype(1,2) == 1
                hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind(randseq(trialnum)))-1);
                hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind(randseq(trialnum)));
            else
                hrtf_left = TF1_space(locind(randseq(trialnum)),:)';
                hrtf_right = TF2_space(locind(randseq(trialnum)),:)';
            end
            str1 = 'FC';           
            % alter mag of IRs
            if XStimParams.alter_left
                % use saved altered HRTFs
                if get(H.altIR_use_saveHRTFflag,'Value') | get(H.altIR_use_saveALTflag,'Value')
                    fname = ['hrtf_left_' el '_' az];
                    fname0 = ['ALTout_' el '_' az];
                    eval(['load ' char(FN.HRTF_path) char(FN.saveHRTF) ' ' fname ' ' fname0]);
                    eval(['hrtf_leftX = ' fname ';']);
                    eval(['ALTin = ' fname0 ';']);
                else            
                    if XStimParams.alter_setbands
                        if XStimParams.alterMag
                            hrtf_leftX = altmagIR3 (hrtf_left, XStimParams.alterMag_Minbandlimits, XStimParams.alterMag_Maxbandlimits, XStimParams.alterMag_dB);
                        end
                        if XStimParams.alterPhase
                            hrtf_leftX = altphaseIR3 (hrtf_left, XStimParams.alterPhase_Minbandlimits, XStimParams.alterPhase_Maxbandlimits, XStimParams.alterPhase_rad);
                        end
                    else            % use 1/12 octave bands
                        if XStimParams.alterMag
                            [hrtf_leftX, temp, ALTout{1}] = altmagIR5 (hrtf_left, find(XStimParams.Mag_octlistBIN ~=0), XStimParams.alterMag_dB,[],XStimParams.Twelths_to_combine);
                        end
                        if XStimParams.alterPhase
                            [hrtf_leftX, temp, ALTout{2}] = altphaseIR5 (hrtf_left, find(XStimParams.Phase_octlistBIN ~=0), XStimParams.alterPhase_rad);        %removed last two argins: [],XStimParams.Twelths_to_combine
                        end
                    end    
                end
                if get(H.altIR_use_saveALTflag,'Value')
                    if XStimParams.alterMag
                        [hrtf_leftX, temp, ALTout{1}] = altmagIR5 (hrtf_left, find(XStimParams.Mag_octlistBIN ~=0), XStimParams.alterMag_dB, ALTin{1},XStimParams.Twelths_to_combine);
                    end
                    if XStimParams.alterPhase
                        [hrtf_leftX, temp, ALTout{2}] = altphaseIR5 (hrtf_left, find(XStimParams.Phase_octlistBIN ~=0), XStimParams.alterPhase_rad, ALTin{2});  % removed last argin: XStimParams.Twelths_to_combine
                    end
                end
            end   
            if XStimParams.alter_right
                if get(H.altIR_use_saveHRTFflag,'Value') | get(H.altIR_use_saveALTflag,'Value')
                    fname = ['hrtf_right_' el '_' az];
                    fname0 = ['ALTout_' el '_' az];
                    eval(['load ' char(FN.HRTF_path) char(FN.saveHRTF) ' ' fname ' ' fname0]);
                    eval(['hrtf_rightX = ' fname ';']);
                    eval(['ALTin = ' fname0 ';']);
                else            
                    if XStimParams.alter_setbands
                        if XStimParams.alterMag
                            hrtf_rightX = altmagIR3 (hrtf_right, XStimParams.alterMag_Minbandlimits, XStimParams.alterMag_Maxbandlimits, XStimParams.alterMag_dB);
                        end
                        if XStimParams.alterPhase
                            hrtf_rightX = altphaseIR3 (hrtf_right, XStimParams.alterPhase_Minbandlimits, XStimParams.alterPhase_Maxbandlimits, XStimParams.alterPhase_rad);
                        end
                    else                    % use 1/12 octave bands
                        if XStimParams.alterMag
                            [hrtf_rightX, temp, ALTout{3}] = altmagIR5 (hrtf_right, find(XStimParams.Mag_octlistBIN ~=0), XStimParams.alterMag_dB,[],XStimParams.Twelths_to_combine);
                        end
                        if XStimParams.alterPhase
                            [hrtf_rightX, temp, ALTout{4}] = altphaseIR5 (hrtf_right, find(XStimParams.Phase_octlistBIN ~=0), XStimParams.alterPhase_rad);      % removed last two argins: [],XStimParams.Twelths_to_combine
                        end
                    end
                end
                if get(H.altIR_use_saveALTflag,'Value')
                    if XStimParams.alterMag
                        [hrtf_rightX, temp, ALTout{3}] = altmagIR5 (hrtf_right, find(XStimParams.Mag_octlistBIN ~=0), XStimParams.alterMag_dB, ALTin{3},XStimParams.Twelths_to_combine);
                    end
                    if XStimParams.alterPhase
                        [hrtf_rightX, temp, ALTout{4}] = altphaseIR5 (hrtf_right, find(XStimParams.Phase_octlistBIN ~=0), XStimParams.alterPhase_rad, ALTin{4});        % removed last argin: XStimParams.Twelths_to_combine
                    end
                end
            end   
            
            temp_left = trial_left;
            temp_right = trial_right;
            trial_left = conv(trial_left,hrtf_left);
            trial_right = conv(trial_right,hrtf_right);
            
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif(XStimParams.ildalone_flag == 1)                  % ILA
            if FN.HRTFfiletype(2,2) == 1      
                hrtf_left = mtlrch([FN.ILA_path FN.ildalone_std],(2*locind(randseq(trialnum)))-1);
                hrtf_right = mtlrch([FN.ILA_path FN.ildalone_std],2*locind(randseq(trialnum)));
            else
                hrtf_left = TF1_ila(locind(randseq(trialnum)),:)';
                hrtf_right = TF2_ila(locind(randseq(trialnum)),:)';
            end
            str1 = 'ILA';
            if XStimParams.alter_left
                % use saved altered HRTFs or ALTerations
                if get(H.altIR_use_saveHRTFflag,'Value') | get(H.altIR_use_saveALTflag,'Value')
                    fname = ['hrtf_left_' el '_' az];
                    fname0 = ['ALTout_' el '_' az];
                    eval(['load ' char(FN.HRTF_path) char(FN.saveHRTF) ' ' fname ' ' fname0]);
                    eval(['hrtf_leftX = ' fname ';']);
                    eval(['ALTin = ' fname0 ';']);
                else            
                    if XStimParams.alter_setbands
                        if XStimParams.alterMag
                            hrtf_leftX = altmagIR3 (hrtf_left, XStimParams.alterMag_Minbandlimits, XStimParams.alterMag_Maxbandlimits, XStimParams.alterMag_dB);
                        end
                        if XStimParams.alterPhase
                            hrtf_leftX = altphaseIR3 (hrtf_left, XStimParams.alterPhase_Minbandlimits, XStimParams.alterPhase_Maxbandlimits, XStimParams.alterPhase_rad);
                        end
                    else        % use 1/12 octave bands
                        if XStimParams.alterMag
                            [hrtf_leftX, temp, ALTout{1}] = altmagIR5 (hrtf_left, find(XStimParams.Mag_octlistBIN ~=0), XStimParams.alterMag_dB,[],XStimParams.Twelths_to_combine);
                        end
                        if XStimParams.alterPhase
                            [hrtf_leftX, temp, ALTout{2}] = altphaseIR5 (hrtf_left, find(XStimParams.Phase_octlistBIN ~=0), XStimParams.alterPhase_rad);        % removed last two argins: v[],XStimParams.Twelths_to_combine
                        end
                    end  
                end
                if get(H.altIR_use_saveALTflag,'Value')
                    if XStimParams.alterMag
                        [hrtf_leftX, temp, ALTout{1}] = altmagIR5 (hrtf_left, find(XStimParams.Mag_octlistBIN ~=0), XStimParams.alterMag_dB, ALTin{1},XStimParams.Twelths_to_combine);
                    end
                    if XStimParams.alterPhase
                        [hrtf_leftX, temp, ALTout{2}] = altphaseIR5 (hrtf_left, find(XStimParams.Phase_octlistBIN ~=0), XStimParams.alterPhase_rad, ALTin{2});      % removed last argin: ,XStimParams.Twelths_to_combine
                    end
                end
            end   
            
            if XStimParams.alter_right
                % use saved altered HRTFs
                if get(H.altIR_use_saveHRTFflag,'Value') | get(H.altIR_use_saveALTflag,'Value')
                    fname = ['hrtf_right_' el '_' az];
                    fname0 = ['ALTout_' el '_' az];
                    eval(['load ' char(FN.HRTF_path) char(FN.saveHRTF) ' ' fname ' ' fname0]);
                    eval(['hrtf_rightX = ' fname ';']);
                    eval(['ALTin = ' fname0 ';']);
                else            
                    if XStimParams.alter_setbands
                        if XStimParams.alterMag
                            hrtf_rightX = altmagIR3 (hrtf_right, XStimParams.alterMag_Minbandlimits, XStimParams.alterMag_Maxbandlimits, XStimParams.alterMag_dB);
                        end
                        if XStimParams.alterPhase
                            hrtf_rightX = altphaseIR3 (hrtf_right, XStimParams.alterPhase_Minbandlimits, XStimParams.alterPhase_Maxbandlimits, XStimParams.alterPhase_rad);
                        end
                    else       %use 1/12 octave bands
                        if XStimParams.alterMag
                            [hrtf_rightX, temp, ALTout{3}] = altmagIR5 (hrtf_right, find(XStimParams.Mag_octlistBIN ~=0), XStimParams.alterMag_dB,[],XStimParams.Twelths_to_combine);
                        end
                        if XStimParams.alterPhase
                            [hrtf_rightX, temp, ALTout{4}] = altphaseIR5 (hrtf_right, find(XStimParams.Phase_octlistBIN ~=0), XStimParams.alterPhase_rad);      % removed last two argins: [],XStimParams.Twelths_to_combine
                        end
                    end
                end
                if get(H.altIR_use_saveALTflag,'Value')
                    if XStimParams.alterMag
                        [hrtf_rightX, temp, ALTout{3}] = altmagIR5 (hrtf_right, find(XStimParams.Mag_octlistBIN ~=0), XStimParams.alterMag_dB, ALTin{3},XStimParams.Twelths_to_combine);
                    end
                    if XStimParams.alterPhase
                        [hrtf_rightX, temp, ALTout{4}] = altphaseIR5 (hrtf_right, find(XStimParams.Phase_octlistBIN ~=0), XStimParams.alterPhase_rad, ALTin{4});        % removed last argin: XStimParams.Twelths_to_combine
                    end
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
        
        % for debugging
        if pflag
            figure;
            subplot(1,2,1); hold on
            plot(hrtf_left,'b')
            if exist('hrtf_leftX')      plot(hrtf_leftX,'r');   end
            subplot(1,2,2); hold on
            plot(hrtf_right,'b')
            if exist('hrtf_rightX')     plot(hrtf_rightX,'r');  end
            pause
            close
        end
        
        % remove DC offset
        trial_left = trial_left - round(mean(trial_left));
        trial_right = trial_right - round(mean(trial_right));
        
        % adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
        % presented through DSPs (*.eq) at 0,0 (broadband)
        trial_left = trial_left * TDT.hrtf_Lfactor;
        trial_right = trial_right * TDT.hrtf_Rfactor;
        
        %%%%%   scale 3/7/07
        trial_left = trial_left * TDT.scaleFactor;
        trial_right = trial_right * TDT.scaleFactor;
        
        %%%%% added padding here from above on 8/29/02      
        % Add in the leading silent period
        trial_left =  [zeros(1,silence_len) trial_left];
        trial_right = [zeros(1,silence_len) trial_right];
        
        % Add in the trailing silent period
        trial_left =  [trial_left zeros(1,silence_len2)];
        trial_right = [trial_right zeros(1,silence_len2)];
        %%%%%      
        
        % pad with zeros for ephonefilters
        filttrial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
        filttrial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
        
        if length(filttrial_left) ~= length(filttrial_right)
            disp('They do not match ?!!!')
        end
        
        % save altered HRTFs to disk
        if get(H.altIR_saveHRTFflag,'Value')
            fname0 = ['ALTout_' el '_' az];
            eval([fname0 '= ALTout;']);
            if XStimParams.alter_left
                fname = ['hrtf_left_' el '_' az];
                eval([fname ' = hrtf_leftX;']);
                eval(['save ' char(FN.HRTF_path) char(FN.saveHRTF) ' ' fname ' ' fname0 ' -append;']);
            end             
            if XStimParams.alter_right   
                fname = ['hrtf_right_' el '_' az];
                eval([fname ' = hrtf_rightX;']);
                eval(['save ' char(FN.HRTF_path) char(FN.saveHRTF) ' ' fname ' ' fname0 ' -append;']);
            end
        end        
        
        % save stimuli to disk
        if(exist1('H.altIRfig'));
            S232('push16',filttrial_left,length(filttrial_left));
            S232('qpop16',BUF.L1);
            fname = ['altIR_' str1 '.left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
                    '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
            evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
                    [' ''' stimuli_dir fname ''' '] ...
                    ',0);'];
            eval(evalstr);
            temp_left = dama2pc(BUF.L1);
            S232('push16',filttrial_right,length(filttrial_right));
            S232('qpop16',BUF.R1);
            fname = ['altIR_' str1 '.right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
                    '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
            evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
                    [' ''' stimuli_dir fname ''' '] ...
                    ',0);'];
            eval(evalstr);
            temp_right = dama2pc(BUF.R1);
        end
        
        
        %Plot PSD
        if(exist1('H.altIRfig') & ~exist1('H.altIR_psdaxes'))
            figure(H.altIRfig);
            H.altIR_psdaxes = axes;
            set(H.altIR_psdaxes,'Visible','off');
        end
        if(exist1('H.altIRfig') & get(H.altIR_plotpsd,'Value') == 1)
            axes(H.altIR_psdaxes);
            set(H.altIR_psdaxes,'Visible','on',...
                'Position',[0.1 0.04 0.8 0.18]);
            axis square
            [pl,f] = psd(trial_left,2048,TDT.Fs); [pr,f] = psd(trial_right,2048,TDT.Fs);
            H.altIR_pl = plot(f,10*log10(abs(pl)),'g');
            hold on
            H.altIR_pr = plot(f,10*log10(abs(pr)),'r');
            xlabel('Frequency (Hz)');
            ILDval = (20*log10(max(abs(trial_left)))) - (20*log10(max(abs(trial_right))));
            ABLval = 0.5*((20*log10(max(abs(trial_left)))) + (20*log10(max(abs(trial_right)))));
            title(['ABL = ' num2str(ABLval) ', ILD = ' num2str(ILDval)]);
            set(H.altIR_psdaxes,'YLim',[0 130]);
            hold off
            grid on
            pause(1);
        elseif(exist1('H.altIRfig') & exist1('H.altIR_pl') & get(H.altIR_plotpsd,'Value') == 0)
            set(H.altIR_psdaxes,'Visible','off');
            set(H.altIR_pl,'Visible','off');
            set(H.altIR_pr,'Visible','off');
        end
        
        remtrials = numtrials - trialnum;
        set(H.altIR_remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        set(H.altIR_buildplay,'BackgroundColor','yellow');
        pause(0);
    end 										%end loop over trials
    
    remreps = 1;
    set(H.altIR_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
end 											%end loop over reps
%%% finished making stimuli

if ~XStimParams.buildOnly
    
    %Begin playing sounds   
    set(H.altIR_buildplay,'String','Playing Stimuli');
    set(H.altIR_buildplay,'BackgroundColor','yellow');
    set(H.altIR_remreps,'String',num2str(XStimParams.numreps));
    repnum = 1;
    buffcycle = 1;
    datamatrix = [];
    
    % increment testnumber
    if(exist1('H.altIRfig') & get(H.altIR_recorddata,'Value'))
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
        tic
        while (exist1('H.altIRfig') & (trialnum <= numtrials))
            
            %Wait till PD1 is finished
            while S232('PD1status',1) usec_delay(1000); end
            
            if(trialnum <= numtrials)
                if(buffcycle == 1);				% load first buffers
                    fname = ['altIR_' str1 '.left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
                            '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
                    evalstr = ['S232(''disk2dama16'',BUF.L1,' ...
                            [' ''' stimuli_dir fname ''' '] ...
                            ',0);'];
                    eval(evalstr);
                    fname = ['altIR_' str1 '.right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
                            '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
                    evalstr = ['S232(''disk2dama16'',BUF.R1,' ...
                            [' ''' stimuli_dir fname ''' '] ...
                            ',0);'];
                    eval(evalstr);
                elseif(buffcycle == 2);				% load second buffers
                    fname = ['altIR_' str1 '.left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
                            '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
                    evalstr = ['S232(''disk2dama16'',BUF.L2,' ...
                            [' ''' stimuli_dir fname ''' '] ...
                            ',0);'];
                    eval(evalstr);
                    fname = ['altIR_' str1 '.right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
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
            
            while toc < ISI/1000 
                if pause_check    return; end
            end
            
            %Check for pause by user
            if pause_check  return; end
            
            %%%%%
            if (trialnum <= numtrials)
                %Start clock
                m110dx( C_.START);
                %Send pulse: PD1 GO!
                m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
                tic
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
            set(H.altIR_remtrials,'String',num2str(remtrials));
            trialnum = trialnum + 1;
            pause(0);
        end %end loop over trials
        
        %Plot Spike Rate Data
        interimspikerate = finalspikematrix/repnum;
        if(exist1('H.altIRfig') & ~exist1('H.altIR_finalspikeratefig'))
            H.altIR_finalspikeratefig = figure('Position',[700 20 550 500],...
                'Name','altIR Test Spike Rate Plot',...
                'NumberTitle','off');
            H.altIR_spikeaxes = axes;
        end
        figure(H.altIR_finalspikeratefig)
        plotdiam1(XStimParams.locations, interimspikerate);
        set(H.altIR_spikeaxes,'Color','black');
        xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
        colorbar
        
        %Record Data
        if(exist1('H.altIRfig') & get(H.altIR_recorddata,'Value'))
            datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
            record_data3(XStimParams,datamatrix);
        end
        
        remreps = XStimParams.numreps - repnum;
        set(H.altIR_remreps,'String',num2str(remreps));
        repnum = repnum + 1;
        pause(0);
    end 		%end loop over reps
    
    %Plot final spike rate figure
    finalspikematrix = finalspikematrix/XStimParams.numreps;
    figure(H.altIR_finalspikeratefig)
    set(H.altIR_finalspikeratefig,'Name','Final Plot for altIR Test');
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.altIR_spikeaxes,'Color','black');
    locmaxspikes = find(finalspikematrix == max(finalspikematrix));
    xlabel('Azimuth'); ylabel('Elevation');
    title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
            ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))],...
        'FontSize',8);
    colorbar
end     %if ~XStimParams.buildOnly

set(H.altIR_buildplay,'String','Build/Play status');
set(H.exitaltIR,'Visible','on');
set(H.resetaltIR,'Visible','on');

% increment test number
if(exist1('H.altIRfig') & get(H.altIR_recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end

%%%%%%%%%
function [flag] = pause_check

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

flag = 0;
%Check for pause by user
while (exist1('H.altIRfig') & get(H.pausealtIR,'Value'))
    pause(0);
    if(~exist1('H.altIRfig')) return; end         
    set(H.exitaltIR,'Visible','on');
    set(H.resetaltIR,'Visible','on');
    if(exist1('H.altIRfig') & get(H.resetaltIR,'Value') == 1)
        set(H.resetaltIR,'Value',0);
        set(H.pausealtIR,'Value',0);
        Reset_altIR;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_altIR;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end