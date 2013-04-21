function [] = Engage_space_FPR()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_space_FPR

%*******************************************************************************
%	The space_FPR Test operation
% no double buffering
% use HRIRs without earphone equalization - usually called *.std
%
% alterred Jan 11, 2011 to keep DATA from multiple reps of same sound(s)
%*******************************************************************************
 colors = [ ...
         179 199 255; ...
        200 150 255; ...
        0    0  255; ...
        216 41  0; ...
        255 199 179;...
        255 150 200;...
        255   0    0;...
        199 255 179;...
        200 255 150;...
        0   255 0]/255;

fclose all;

if XStimParams.HiDynamicRange
    scaleFactor = 60;
else
    scaleFactor = TDT.scaleFactor;
end

% reset stim_type to file (#9)
set(H.stim_type,'Value',9);

XStimParams.curr_ITD = str2num(get(H.space_FPR_ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.space_FPR_ABL,'String'));
XStimParams.curr_stimdur = str2num(get(H.space_FPR_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space_FPR_ISI,'String'));
XStimParams.numreps = str2num(get(H.space_FPR_numreps,'String'));
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

% check to be sure stimFNs chosen
if isempty(FN.FPR)         %Stimulus from file
    set(H.space_FPR_initFN,'value',1)
    setinfo_space_FPR
end
nFiles = length(FN.FPR);
for iFile=1:nFiles
D = dir([FN.stim_path2 FN.FPR{iFile}]);
B(iFile) = D.bytes/4;
end
DUR = max1(B) / round(TDT.Fs/1000);
XStimParams.curr_stimdur = DUR;

%Make Stimulus buffers
% Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
DUR_silence = XStimParams.silence_lead; 					%ms
silence_len = (DUR_silence * round(TDT.Fs/1000));
% Add a piece of silence after stimulus 
DUR_silence2 = XStimParams.silence_trail; 					%ms
silence_len2 = (DUR_silence2 * round(TDT.Fs/1000));

if XStimParams.ildalone_flag 				%ILDAlone Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.ephonefiltlen + TDT.itdfiltlen;
    XStimParams.test_type = 'space_FPR ILA';     str1 = 'ila';
    disp('This is an ILD-ALONE space_FPR test')
elseif XStimParams.itdalone_flag 			%ITDAlone Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + (2* TDT.hrtffiltlen) + TDT.ephonefiltlen;
    XStimParams.test_type = 'space_FPR ITA';        str1 = 'ita';
    disp('This is an ITD-ALONE space_FPR test')
elseif XStimParams.space_flag 			%fully-cued Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.ephonefiltlen;
    XStimParams.test_type = 'space_FPR FC';     str1 = 'fc';
    disp('This is a FULLY CUED space_FPR test')
elseif XStimParams.ABLalone_flag			%ABL-alone Test
    Npts_totalplay = (silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.hrtffiltlen + TDT.ephonefiltlen;
    XStimParams.test_type = 'space_FPR ABA';        str1 = 'aba';
    disp('This is an ABL-alone space_FPR test')
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
if get(H.space_FPR_recorddata,'Value')
    tempstr = ['    ' str1 ' space_FPR-test: ' num2str(size(XStimParams.locations,2)) ' locations;     atten: ' num2str(abs(ABL))];
    update_diary
end

Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_space_FPR_' str1 ' Temp_params;'])
clear Temp_params str1

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
setInfo_space_FPR;

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
% numtrials, trialnum reflects all iFiles and all iLocs
remreps = 1;
set(H.space_FPR_status,'String','Status: building stimuli');
set(H.space_FPR_status,'BackgroundColor','yellow');
set(H.space_FPR_status,'ForegroundColor','blue');
set(H.space_FPR_remreps,'String',num2str(remreps));
repnum = 1;
nLocs = size(XStimParams.locations,2);
nFiles = size(FN.FPR,2);
numtrials = nLocs * nFiles;
finalspikematrix = zeros(1,numtrials);

%Randomize stimFNs
FN.FPR = FN.FPR(randperm(nFiles));

% clear out temporary directory
delete( [FN.temp_stim_path '\*.*'])

trialnum = 1;
for iFile = 1:nFiles
    %Check for pause by user
    if pause_check    return; end
    % get stimulus from file
    if strcmp('mat',FN.FPR{iFile}(end-2:end))
        eval(['load ' FN.stim_path2 FN.FPR{iFile}])
        trial_left  = stim;
    else
        fid = fopen([FN.stim_path2 FN.FPR{iFile}],'r');
        trial_left = fread(fid,inf,'float');
        fclose(fid);
    end
    trial_left = trial_left(:)' - mom(trial_left,1);
    
    % trailing pad stim length of short stims
    if length(trial_left)/TDT.Fs *1000 < DUR
        npts = ceil(DUR*30 - length(trial_left));
        trial_left = [trial_left zeros(1,npts)];
    end
    
    % normalize to ACPower
    if ~XStimParams.HiDynamicRange
        trial_left = trial_left / mom(trial_left,2);
    end
    trial_right = trial_left;
    
    %Apply ITD filtering if conducting ILDAlone Test
    if(XStimParams.ildalone_flag == 1)
        itdleft = 0; itdright = 0;
        ITD = round(str2num(get(H.space_FPR_ITD,'String')));
        if(ITD < 0)
            itdright = abs(ITD);
        elseif(ITD > 0)
            itdleft = abs(ITD);
        end
        if(iFile == 1)
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
        if(iFile == 1)
            if FN.HRTFfiletype(2,2) == 1
                dir_ILA = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone_std]));
            end
            ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
            if isempty('ILAind')
                error('engage397: Could not find ILA location in HRTF file');
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
    trial_left0 = trial_left;
    trial_right0 = trial_right;
    
    for iLoc = 1:nLocs        %Apply HRTF filtering
        if pause_check    return; end        %Check for pause by user
        if(XStimParams.space_flag == 1)
            if FN.HRTFfiletype(1,2) == 1
                hrtf_left = mtlrch([FN.space_path FN.space_std],(2*locind(iLoc))-1);
                hrtf_right = mtlrch([FN.space_path FN.space_std],2*locind(iLoc));
            else
                hrtf_left = TF1_space(locind(iLoc),:);
                hrtf_right = TF2_space(locind(iLoc),:);
            end
            str1 = 'FC';
        elseif(XStimParams.ABLalone_flag == 1)
            % get ABLfactor for this location
            if FN.HRTFfiletype(1,2) == 1
                tempL = mtlrch([FN.space_path FN.space_std],(2*locind(iLoc))-1);
                tempR = mtlrch([FN.space_path FN.space_std],2*locind(iLoc));
            else
                tempL = TF1_space(locind(iLoc),:);
                tempR = TF2_space(locind(iLoc),:);
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
                hrtf_left = mtlrch([FN.ILA_path FN.ildalone_std],(2*locind(iLoc))-1);
                hrtf_right = mtlrch([FN.ILA_path FN.ildalone_std],2*locind(iLoc));
            else
                hrtf_left = TF1_ila(locind(iLoc),:);
                hrtf_right = TF2_ila(locind(iLoc),:);
            end
            str1 = 'ILA';
        elseif(XStimParams.itdalone_flag == 1)
            if FN.HRTFfiletype(3,2) == 1      
                hrtf_left = mtlrch([FN.ITA_path FN.itdalone_std],(2*locind(iLoc))-1);
                hrtf_right = mtlrch([FN.ITA_path FN.itdalone_std],2*locind(iLoc));
            else
                hrtf_left = TF1_ita(locind(iLoc),:);
                hrtf_right = TF2_ita(locind(iLoc),:);        
            end
            str1 = 'ITA';
        end
        
        trial_left = conv(trial_left0, hrtf_left);
        trial_right = conv(trial_right0, hrtf_right);
        
        % remove DC offset
        trial_left = trial_left - mean(trial_left);
        trial_right = trial_right - mean(trial_right);
        
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

        if length(trial_left) ~= Npts_totalplay | length(trial_right) ~= Npts_totalplay
            error('engage499: filesize incorrect')
        end
        
        % save stimuli to disk
        S232('push16',trial_left,length(trial_left));
        S232('qpop16',BUF.L1);
        fname = ['FPR' num2str(iFile) '_' str1 '.left_' ...
                num2str(hrtfdirmat(1,locind(iLoc))) ...
                '_' num2str(hrtfdirmat(2,locind(iLoc)))];
        evalstr = ['S232(''dama2disk16'',BUF.L1,' [' ''' FN.temp_stim_path fname ''' '] ',0);'];
        eval(evalstr);
        
        S232('push16',trial_right,length(trial_right));
        S232('qpop16',BUF.R1);
        fname = ['FPR' num2str(iFile) '_' str1 '.right_' ...
                num2str(hrtfdirmat(1,locind(iLoc))) ...
                '_' num2str(hrtfdirmat(2,locind(iLoc)))];
        evalstr = ['S232(''dama2disk16'',BUF.R1,' [' ''' FN.temp_stim_path fname ''' '] ',0);'];
        eval(evalstr);
        
        % set up for next trial                        
        remtrials = numtrials - trialnum;
        set(H.space_FPR_remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end     % iLoc
end 	% iFile

%%%%%%%%%%%%%%%%%%% finished making stimuli

%Begin playing sounds   
set(H.space_FPR_status,'String','Status: playing stimuli');
set(H.space_FPR_status,'BackgroundColor','green');
set(H.space_FPR_remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.space_FPRfig') & get(H.space_FPR_recorddata,'Value'))
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
%loop for reps
    trialnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    FN_trial = [];
    Nspikes = [];
while (exist1('H.space_FPRfig') & (repnum <= XStimParams.numreps))
    
    tic
    for iFile = 1:nFiles
        randseq = randperm(nLocs);   
        for iLoc = 1:nLocs
            %Check for pause by user
            if pause_check    return; end
            
            %%% load stimuli
            %% push left file onto stack & check filesize
            fname = ['FPR' num2str(iFile) '_' str1 '.left_' ...
                    num2str(hrtfdirmat(1,locind(randseq(iLoc)))) ...
                    '_' num2str(hrtfdirmat(2,locind(randseq(iLoc))))];
            evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
            eval(evalstr);
            if S232('topsize') ~= Npts_totalplay
                error('engage567: LEFT file to play is incorrect size')
            end
            %% pop to DAMA
            eval(['S232(''qpop16'''  ',BUF.L1);']);
            
            %% push right file onto stack & check filesize            
            fname = ['FPR' num2str(iFile) '_' str1 '.right_' ...
                    num2str(hrtfdirmat(1,locind(randseq(iLoc)))) ...
                    '_' num2str(hrtfdirmat(2,locind(randseq(iLoc))))];
            evalstr = ['S232(''pushdisk16''' ',' [' ''' FN.temp_stim_path fname ''' '] ');'];
            eval(evalstr);
            if S232('topsize') ~= Npts_totalplay
                error('engage579: RIGHT file to play is incorrect size')
            end
            %% pop to DAMA
            eval(['S232(''qpop16'''  ',BUF.R1);']);
            
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
            
            %Stop the m110 and get spiketimes
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, round(XStimParams.curr_stimdur*2)); 			% Take 2*XStimParams.curr_stimdur spikes max
            ind = find(spikes ~= 0); 						% Get clock events that are spikes
            spikes = spikes(ind);
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            if(exist1('H.space_FPRfig') & ~isempty(spikes)) 
                spikes = spikes(ind);
                spikes_trial = [spikes_trial; spikes/(1000/mii_us_per_sample)];
                EL_trial = [EL_trial; hrtfdirmat(1,locind(randseq(iLoc)))* ones(size(spikes))];
                AZ_trial = [AZ_trial; hrtfdirmat(2,locind(randseq(iLoc)))* ones(size(spikes))];
                repnum_trial = [repnum_trial; repnum * ones(size(spikes))];
                FN_trial = [FN_trial; iFile* ones(size(spikes))];
                Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
            end
            ind = (iLoc-1)*nFiles + iFile;
            finalspikematrix(ind) = finalspikematrix(ind) + length(spikes);
            if pause_check    return; end
            
            remtrials = numtrials - trialnum;
            set(H.space_FPR_remtrials,'String',num2str(remtrials));
            trialnum = trialnum + 1;
            pause(0);
        end             % iLocs
        
        %Record Data
        if(exist1('H.space_FPRfig') & get(H.space_FPR_recorddata,'Value'))
            datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial FN_trial]];
            record_data3(XStimParams,datamatrix);
        end
    end                 % iFiles
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.space_FPRfig') & ~exist1('H.space_FPR_finalspikeratefig'))
        H.space_FPR_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','space_FPR Test Spike Rate Plot',...
            'NumberTitle','off');
        H.space_FPR_spikeaxes = axes;
    end
    figure(H.space_FPR_finalspikeratefig); hold off
    %plotdiam1(XStimParams.locations, interimspikerate);
    for iLoc = 1:nLocs
       ind = [1:nFiles] + (iLoc-1)*nFiles;
       plot(1:nFiles,finalspikematrix(ind),'color',colors(rem(iLoc,10)+1,:))
       hold on
    end
    set(H.space_FPR_spikeaxes,'Color','black');
    xlabel('file#'); ylabel('spikerate'); title(['Rep # ' num2str(repnum)]);
            
    remreps = XStimParams.numreps - repnum;
    set(H.space_FPR_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.space_FPR_finalspikeratefig); hold off
set(H.space_FPR_finalspikeratefig,'Name','Final Plot for space_FPR Test');
    for iLoc = 1:nLocs
       ind = [1:nFiles] + (iLoc-1)*nFiles;
       plot(1:nFiles,finalspikematrix(ind),'color',colors(rem(iLoc,10)+1,:))
       hold on
    end
set(H.space_FPR_spikeaxes,'Color','black');
xlabel('file#'); ylabel('spikerate');
title('good stuff - finished!', 'FontSize',8);

set(H.space_FPR_status,'String','Status: results');
set(H.space_FPR_status,'BackgroundColor','blue');

set(H.exitspace_FPR,'Visible','on');
set(H.resetspace_FPR,'Visible','on');

% increment test number
if(exist1('H.space_FPRfig') & get(H.space_FPR_recorddata,'Value'))
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
while (exist1('H.space_FPRfig') & get(H.pausespace_FPR,'Value'))
    pause(0);
    if(~exist1('H.space_FPRfig')) return; end         
    set(H.exitspace_FPR,'Visible','on');
    set(H.resetspace_FPR,'Visible','on');
    if(exist1('H.space_FPRfig') & get(H.resetspace_FPR,'Value') == 1)
        set(H.resetspace_FPR,'Value',0);
        set(H.pausespace_FPR,'Value',0);
        Reset_space_FPR;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_space_FPR;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end