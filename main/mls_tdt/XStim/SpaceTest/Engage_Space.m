function [] = Engage_space()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_space

%*******************************************************************************
%	The space Test operation
% altered version to eliminate double buffering
% should use HRIRs without earphone equalization - usually called *.std
%*******************************************************************************

stimuli_dir = FN.temp_stim_path;
fclose all;
if ~XStimParams.buildOnly     % buildOnly keeps adding files to folder (clear these from Xstim menu)
    eval(['delete ' stimuli_dir '*.*;']);
end

if XStimParams.HiDynamicRange
    scaleFactor = 60;
else
    scaleFactor = TDT.scaleFactor;
end

% reset non-allowed stim_type to BroadBand
if isempty(find([1 2 6 8 9] == get(H.stim_type,'Value')))
    set(H.stim_type,'Value',8);
    disp('Unsupported stimulus type for space test; reset to BROADBAND');
end

XStimParams.curr_ITD = str2num(get(H.space_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.space_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.space_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space_ISI,'String'));
XStimParams.numreps = str2num(get(H.space_numreps,'String'));
XStimParams.reset_flag = 0;

% load filter coeffs & other params for ABLalone test
if XStimParams.ABLalone_flag
    Fs = 30000;
    cF = round(1000*exp(([40:168]/48)*log(2)))'; 
    n_cF = length(cF);
    fcoefs = Make_ERBFiltA(Fs,cF);
    maxFactor = .00003764*cF(n_cF)+.6236;
    Factor = maxFactor ./ (.00003764*cF+.6236);
    Factormat = repmat(Factor,1,255);
    fftpts = 2048;
    freq = 0:15000/((fftpts/2)-1):15000;
    freq_ind = nearest_index(freq,cF);
    endpt = max1(freq_ind);
    startpt = min1(freq_ind);
    Xpart1 = startpt:endpt;
    Xpart2 = (fftpts + 2 - endpt):(fftpts + 2 - startpt);
end

%Specify DAMA buffers
clear BUF
BUF.L1				= 1;
BUF.R1				= 2;
BUF.playseq_L1		= 5;
BUF.playseq_R1		= 6;
BUF.playspec1		= 7;

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

%Make play specification buffer
S232('allot16',BUF.playspec1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L1);
S232('make',1,BUF.playseq_R1);
S232('make',2,0);
S232('qpop16',BUF.playspec1);

% Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
DUR_silence = XStimParams.silence_lead; 					%ms
silence_len = (DUR_silence * round(TDT.Fs/1000));
% Add a piece of silence after stimulus 
DUR_silence2 = XStimParams.silence_trail; 					%ms
silence_len2 = (DUR_silence2 * round(TDT.Fs/1000));

if get(H.stim_type,'Value') == 9
    if isempty(FN.stim)         %Stimulus from file
        [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
    end
    D = dir([FN.stim_path FN.stim]);
    XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
end

%Make Stimulus buffers
DUR = XStimParams.curr_stimdur;
if XStimParams.ildalone_flag 				%ILDAlone Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.ephonefiltlen + TDT.itdfiltlen;
    XStimParams.test_type = 'Space ILA';     str1 = 'ILA';
    disp('This is an ILD-ALONE space test')
elseif XStimParams.itdalone_flag 			%ITDAlone Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + (2* TDT.hrtffiltlen) + TDT.ephonefiltlen;
    XStimParams.test_type = 'Space ITA';        str1 = 'ITA';
    disp('This is an ITD-ALONE space test')
elseif XStimParams.space_flag 			%fully-cued Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.ephonefiltlen;
    XStimParams.test_type = 'Space FC';     str1 = 'FC';
    disp('This is a FULLY CUED space test')
elseif XStimParams.ABLalone_flag			%ABL-alone Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.ephonefiltlen;
    XStimParams.test_type = 'Space ABA';        str1 = 'ABA';
    disp('This is an ABL-alone space test')
end   

S232('allot16',BUF.L1,Npts_totalplay);
S232('allot16',BUF.R1,Npts_totalplay);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,Npts_totalplay);

%Get Earphone filters
if FN.HRTFfiletype(6) == 1
    ephonefname = [FN.ephone_path FN.ephone2];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
else
    eval(['load -mat ' FN.ephone_path FN.ephone2]);
    ephonefilt_left  = TF1;
    ephonefilt_right = TF2;
    clear TF1 TF2 dir
end

%Load Earphone filters
dspid_left = 0; dspid_right = 1;
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));           % used to be '0xFFF'
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
ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));			%correct for ephonefilt

%Get the HRTF indices
XStimParams.locations = GUI.locations1';
if ~length(XStimParams.locations)
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    pause;
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end

% save XStimParams for this test
if get(H.space_recorddata,'Value')
    tempstr = ['    ' str1 'space     atten: ' num2str(abs(ABL)) 'dB'];
    update_diary
end

Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_Space_' str1 ' Temp_params;'])
clear Temp_params

if XStimParams.space_flag | XStimParams.ABLalone_flag
    if isempty(FN.space_std) | ~isempty(strfind(FN.space_std,'eq'))        % FN not yet picked
        [FN.space_std,FN.space_path] = uigetfile('*.*','Select Fully-cued HRTF File (*.std preferred)');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space_std]);
        end
        set(H.spacefilepb,'Value',0);
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
elseif XStimParams.ildalone_flag
    if isempty(FN.ildalone_std)        % FN not yet picked
        [FN.ildalone_std,FN.ILA_path] = uigetfile('*.*','Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.space_ildalonefile,'String',[FN.ILA_path FN.ildalone_std]);
        end
        set(H.space_ildalonepb,'Value',0);
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

elseif XStimParams.itdalone_flag
    if isempty(FN.itdalone_std)        % FN not yet picked
        [FN.itdalone_std,FN.ITA_path] = uigetfile('*.*','Select ITD-alone HRTF File');
        if(FN.ITA_path ~= 0)
            set(H.space_itdalonefile,'String',[FN.ITA_path FN.itdalone_std]);
        end
        set(H.space_itdalonepb,'Value',0);
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
            set(H.space_ildalonefile,'String',[FN.ILA_path FN.ildalone_std]);
        end
        set(H.space_ildalonepb,'Value',0);
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
end
clear locind
for locnum = 1:size(XStimParams.locations,2)
    locind(locnum) = max(find(hrtfdirmat(1,:) == XStimParams.locations(1,locnum) &...
        hrtfdirmat(2,:) == XStimParams.locations(2,locnum)));
end

% update display
setInfo_space;

% get reference Lref, Rref and ABL for ABLalone test
if XStimParams.ABLalone_flag
    ind00 = max(find(hrtfdirmat(1,:) == XStimParams.el & hrtfdirmat(2,:) == XStimParams.az));
    if FN.HRTFfiletype(1,2) == 1
        Lref = mtlrch([FN.space_path FN.space_std],(2*ind00)-1);
        Rref = mtlrch([FN.space_path FN.space_std],2*ind00);
    else
        Lref = TF1_space(ind00,:)';
        Rref = TF2_space(ind00,:)';
    end
    tempL = ERBFilterBankB(Lref, fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
    tempR = ERBFilterBankB(Rref, fcoefs) .* Factormat;
    [temp ablref] = calclevel_time(tempL,tempR, cF);
    clear tempL tempR temp
end

%%%%% Loop to make the stimuli we'll use
remreps = 1;
set(H.space_status,'String','Status: building stimuli');
set(H.space_status,'BackgroundColor','blue');

set(H.space_remreps,'String',num2str(remreps));
repnum = 1;
numtrials = size(XStimParams.locations,2);
finalspikematrix = zeros(1,numtrials);
clear stim_from_file

%Randomize the stimuli
randseq = randperm(numtrials);

trialnum = 1;
while (exist1('H.spacefig') & (trialnum <= numtrials))
    set(H.space_status,'BackgroundColor','yellow');
    %Check for pause by user
    if pause_check    return; end
    
    %Make the stimuli
    switch get(H.stim_type,'Value')      
        case 1            % tone
            [trial_left] = MakeTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur);
            trial_right = trial_left;
        case 2         % GammaTones
            trial_left = MakeGammaTone(TDT.Fs,XStimParams.curr_freq,XStimParams.curr_stimdur);
            trial_right = trial_left;
        case 3
            %Use Concurrent ToneStack if specified: 10 frequencies
            tonestackfreqs = ...
                XStimParams.bandwidth(1):...
                round((XStimParams.bandwidth(2) - XStimParams.bandwidth(1))/9):...
                XStimParams.bandwidth(2);
            [trial_left,trial_right] = MakeToneStack(TDT.Fs,tonestackfreqs,XStimParams.curr_stimdur);
        case 6
            %Narrowband Noise
            [trial_left] = MakeNBNoise(TDT.Fs,XStimParams.curr_stimdur,XStimParams.bandlimits);
            trial_right = trial_left;
        case 7
            %Narrowband Noise Stack
            [trial_left] = MakeNBNoiseStack(TDT.Fs,XStimParams.curr_stimdur,XStimParams.spec_NBNoiseStackBands);
            trial_right = trial_left;
        case 8
            %Broadband Noise
            [trial_left] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
            trial_right = trial_left;
        case 9
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
            disp('Unsupported stimulus type for space test; reset to BROADBAND');
    end
    
    % remove any DCoffset
    trial_left = trial_left - mom(trial_left,1);
    trial_right = trial_right - mom(trial_right,1);
    
    % normalize to ACPower (added 3/5/07 as in 2-source)
    if ~XStimParams.HiDynamicRange
        trial_left = trial_left / mom(trial_left,2);
        trial_right = trial_right / mom(trial_right,2);
    end

    % modulate stim1
    if get(H.stim_type,'Value') ~= 9 & ~strcmp(XStimParams.mod_type,'None')
        Envelope = make_env(DUR, XStimParams.mod_type, XStimParams.mod_depth(1), XStimParams.mod_freq(1), XStimParams.mod_phase(1));
        trial_right = trial_right .* Envelope(:)';
        trial_left = trial_left .* Envelope(:)';
    end
    
    if get(H.stim_type,'Value') ~= 9
        %Ramp the stimuli
        ramp_time = 5; %ms
        [trial_left] = ramp_sound(trial_left,TDT.Fs,ramp_time);
        [trial_right] = ramp_sound(trial_right,TDT.Fs,ramp_time);
    end
    
     % remove any DCoffset
    trial_left = trial_left - mom(trial_left,1);
    trial_right = trial_right - mom(trial_right,1);

    %Apply ITD filtering if conducting ILDAlone space Test
    if(XStimParams.ildalone_flag == 1)
        itdleft = 0; itdright = 0;
        ITD = round(str2num(get(H.space_ITD,'String')));
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
    
    % Apply ILD filtering if conducting ITDalone Test
    if(XStimParams.itdalone_flag == 1)
        if(trialnum == 1)
            if FN.HRTFfiletype(2,2) == 1
                dir_ILA = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone_std]));
            end
            ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
            if isempty('ILAind')
                disp('Could not find ILA location in HRTF file');
                return
            end
            if FN.HRTFfiletype(2,2) == 1
                eval(['ila_filt_left = mtlrch(''' FN.ILA_path FN.ildalone_std ''', ' num2str(ILAind * 2-1) ');']);
                eval(['ila_filt_right = mtlrch(''' FN.ILA_path FN.ildalone_std ''', ' num2str(ILAind * 2) ');']);
            else
                ila_filt_left = TF1_ila(ILAind,:);
                ila_filt_right = TF2_ila(ILAind,:);
            end
        end
        trial_left = conv(trial_left,ila_filt_left);
        trial_right = conv(trial_right,ila_filt_right);
    end
    
    %Apply HRTF filtering
    if(XStimParams.space_flag == 1)
        if FN.HRTFfiletype(1,2) == 1
            hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind(randseq(trialnum)))-1);
            hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind(randseq(trialnum)));
        else
            hrtf_left = TF1_space(locind(randseq(trialnum)),:);
            hrtf_right = TF2_space(locind(randseq(trialnum)),:);
        end
        str1 = 'FC';
    elseif(XStimParams.ABLalone_flag == 1)
        % get ABLfactor for this location
        if FN.HRTFfiletype(1,2) == 1
            tempL = mtlrch([FN.space_path FN.space_std],(2*locind(randseq(trialnum)))-1);
            tempR = mtlrch([FN.space_path FN.space_std],2*locind(randseq(trialnum)));
        else
            tempL = TF1_space(locind(randseq(trialnum)),:);
            tempR = TF2_space(locind(randseq(trialnum)),:);
        end
        tempL = ERBFilterBankB(tempL, fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
        tempR = ERBFilterBankB(tempR, fcoefs) .* Factormat;
        [ildx ablx] = calclevel_time(tempL,tempR, cF);
        ABLfactor = (10 .^(ablx/20)) ./ (10 .^(ablref/20));
        ABLfactor_long = (interp1(freq(freq_ind),ABLfactor,freq(startpt:endpt)))';
        % apply ABLfactor to reference location
        FT_L = fft(Lref,fftpts);
        FT_L(Xpart1)= FT_L(Xpart1) .* (ABLfactor_long);       	% positive freqs
        FT_L(Xpart2)= FT_L(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
        hrtf_left = real(ifft(FT_L));
        hrtf_left = hrtf_left(1:255);
        FT_R = fft(Rref,fftpts);
        FT_R(Xpart1)= FT_R(Xpart1) .* (ABLfactor_long);       	% positive freqs
        FT_R(Xpart2)= FT_R(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
        hrtf_right = real(ifft(FT_R));
        hrtf_right = hrtf_right(1:255);
        clear FT_L FT_R ABLfactor* 
        
    elseif(XStimParams.ildalone_flag == 1)
        if FN.HRTFfiletype(2,2) == 1      
            hrtf_left = mtlrch([FN.ILA_path FN.ildalone_std],(2*locind(randseq(trialnum)))-1);
            hrtf_right = mtlrch([FN.ILA_path FN.ildalone_std],2*locind(randseq(trialnum)));
        else
            hrtf_left = TF1_ila(locind(randseq(trialnum)),:);
            hrtf_right = TF2_ila(locind(randseq(trialnum)),:);
        end
        str1 = 'ILA';
    elseif(XStimParams.itdalone_flag == 1)
        if FN.HRTFfiletype(3,2) == 1      
            hrtf_left = mtlrch([FN.ITA_path FN.itdalone_std],(2*locind(randseq(trialnum)))-1);
            hrtf_right = mtlrch([FN.ITA_path FN.itdalone_std],2*locind(randseq(trialnum)));
        else
            hrtf_left = TF1_ita(locind(randseq(trialnum)),:);
            hrtf_right = TF2_ita(locind(randseq(trialnum)),:);        
        end
            str1 = 'ITA';
    end
    
    trial_left = conv(trial_left, hrtf_left);
    trial_right = conv(trial_right, hrtf_right);

    % remove DC offset
    trial_left = trial_left - round(mean(trial_left));
    trial_right = trial_right - round(mean(trial_right));

    % adjust headphone stimuli pre-convolved with hrtfs (*.std) to match those
    % presented through DSPs (*.eq) at 0,0 (broadband)
    trial_left = trial_left * TDT.hrtf_Lfactor;
    trial_right = trial_right * TDT.hrtf_Rfactor;
        
    %scale stimuli 3/7/07
    % note scale factor changes with XStimParams.HiDynamicRange
    trial_left = trial_left * scaleFactor;
    trial_right = trial_right * scaleFactor;
    
    %Add in the leading silent period
    trial_left =  [zeros(1,silence_len) trial_left];
    trial_right = [zeros(1,silence_len) trial_right];
    
    %Add in the trailing silent period
    trial_left =  [trial_left zeros(1,silence_len2)];
    trial_right = [trial_right zeros(1,silence_len2)];
    
    %pad with zeros for ephonefilters
    trial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
    trial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
    
    % save stimuli to disk
    S232('push16',trial_left,length(trial_left));
    S232('qpop16',BUF.L1);
    fname = ['space1_' str1 '.left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
            '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
    evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
            [' ''' FN.temp_stim_path fname ''' '] ',0);'];
    eval(evalstr);
    %temp_left = dama2pc(BUF.L1);
    S232('push16',trial_right,length(trial_right));
    S232('qpop16',BUF.R1);
    fname = ['space1_' str1 '.right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
            '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
    evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
            [' ''' FN.temp_stim_path fname ''' '] ',0);'];
    eval(evalstr);
    %temp_right = dama2pc(BUF.R1);
    
    % set up for next trial                        
    remtrials = numtrials - trialnum;
    set(H.space_remtrials,'String',num2str(remtrials));
    trialnum = trialnum + 1;
    set(H.space_status,'BackgroundColor','blue');
    pause(0);
end 										%end loop over trials
%%%%%%%%%%%%%%%%%%% finished making stimuli
if ~XStimParams.buildOnly
%Begin playing sounds   
set(H.space_status,'String','Status: playing stimuli');
set(H.space_status,'BackgroundColor','green');
set(H.space_remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
buffcycle = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.spacefig') & get(H.space_recorddata,'Value'))
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
%loop for reps
while (exist1('H.spacefig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli
    randseq = randperm(numtrials);   
    trialnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    tic
    while (exist1('H.spacefig') & (trialnum <= numtrials) )
        %Check for pause by user
        if pause_check    return; end
        
        % load stimuli
        fname = ['space1_' str1 '.left_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
                '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
        eval(evalstr);
        evalstr = ['S232(''qpop16'''  ',BUF.L1);'];
        eval(evalstr);
        fname = ['space1_' str1 '.right_' num2str(hrtfdirmat(1,locind(randseq(trialnum)))) ...
                '_' num2str(hrtfdirmat(2,locind(randseq(trialnum))))];
        evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
        eval(evalstr);
        evalstr = ['S232(''qpop16'''  ',BUF.R1);'];
        eval(evalstr);
                
        S232('seqplay',BUF.playspec1);
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

        %Start clock
        m110dx( C_.START);
        %Send pulse: PD1 GO!
        m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
        tic
        
        pause(DUR/1000+.2);
        
        while(S232('PD1status',1)) usec_delay(1000);  end
        S232('PD1stop',1);
       
        %Stop the m110 and get spikes
        m110dx( C_.STOP);
        spikes = m110dx( C_.DATA, XStimParams.curr_stimdur*2); 			% Take 2*XStimParams.curr_stimdur spikes max
        ind = find(spikes ~= 0); 						% Get clock events that are spikes
        spikes = spikes(ind);
        ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
        if(exist1('H.spacefig') & ~isempty(spikes)) 
            spikes = spikes(ind);
            spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
            EL_trial = [EL_trial;hrtfdirmat(1,locind(randseq(trialnum)))* ones(size(spikes))];
            AZ_trial = [AZ_trial;hrtfdirmat(2,locind(randseq(trialnum)))* ones(size(spikes))];
            repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
            Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
        end
        finalspikematrix(randseq(trialnum)) = finalspikematrix(randseq(trialnum)) + length(spikes);
        if pause_check    return; end
    
        remtrials = numtrials - trialnum;
        set(H.space_remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end %end loop over trials
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.spacefig') & ~exist1('H.space_finalspikeratefig'))
        H.space_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','space Test Spike Rate Plot',...
            'NumberTitle','off');
        H.space_spikeaxes = axes;
    end
    figure(H.space_finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.space_spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    %Record Data
    if(exist1('H.spacefig') & get(H.space_recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.space_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.space_finalspikeratefig)
set(H.space_finalspikeratefig,'Name','Final Plot for space Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.space_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
        ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))], 'FontSize',8);
colorbar

% plot weights
if get(H.space_plotwts,'Value')
    set(H.space_status,'String','Status: calculating weights');
    set(H.space_status,'BackgroundColor','red');
    plot_MLSweights(finalspikematrix,locind);
    set(H.space_plotwts,'value',0);
end
end         % if ~XStimParams.buildOnly

set(H.space_status,'String','Status: results');
set(H.space_status,'BackgroundColor','blue');

set(H.exitspace,'Visible','on');
set(H.resetspace,'Visible','on');

% increment test number
if(exist1('H.spacefig') & get(H.space_recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end


%%%%%%%%%%%
function [Envelope] = make_env(DUR,mod_type,param1,param2,param3)

% param1 used as follows:
%    'tone' or 'LP noise' or 'File': mod_depth
% param2 used as follows:
%    'tone or 'LP noise' : mod_freq
% param3 used as follows
%    'tone' : addition to modulation starting phase (0 - pi)

global H
global TDT
global FN

if nargin < 3   param3 = 0; end

Npts = DUR*(round(TDT.Fs/1000));
Inc = 1/TDT.Fs;
switch mod_type
    case 'Sq wave'
        T = 0:Inc:(DUR/1000 - Inc);
        Tone = (param1 / 2)* square(2 * pi * param2 .* T + (.75 * 2 * pi + param3));
        Envelope = Tone + (1-param1/2);
    case 'Tone'
        T = 0:Inc:(DUR/1000 - Inc);
        Tone = (param1 / 2)* sin(2 * pi * param2 .* T + (.75 * 2 * pi + param3));
        Envelope = Tone + (1-param1/2);
    case 'LP Noise'
        LP_noise = m_noi(5, param2, Npts/(TDT.Fs), TDT.Fs/2);
        LP_noise = (param1 / 2)* (LP_noise / max1(LP_noise));
        Envelope = LP_noise + (1-param1/2);
    case 'File'			
        if(~exist1('mod_from_file'))
            fid = fopen(FN.mod,'r');
            mod_from_file = fread(fid,inf,'float');
            fclose(fid);
            while length(mod_from_file) ~= Npts
                [FN.mod,FN.mod_path] = uigetfile('*.*','Select Envelope File');
                if(FN.mod_path ~= 0)
                    set(H.space_modfile,'String',[FN.mod_path FN.mod]);
                end
                fid = fopen(FN.mod,'r');
                mod_from_file = fread(fid,inf,'float');
                fclose(fid);
            end
            mod_from_file = mod_from_file - mean(mod_from_file);
            mod_from_file = (param1 / 2)* (mod_from_file / max1(mod_from_file));
        end
        Envelope = mod_from_file + (1-param1/2);
    otherwise
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
while (exist1('H.spacefig') & get(H.pausespace,'Value'))
    pause(0);
    if(~exist1('H.spacefig')) return; end         
    set(H.exitspace,'Visible','on');
    set(H.resetspace,'Visible','on');
    if(exist1('H.spacefig') & get(H.resetspace,'Value') == 1)
        set(H.resetspace,'Value',0);
        set(H.pausespace,'Value',0);
        Reset_space;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_space;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end