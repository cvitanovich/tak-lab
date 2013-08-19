function [] = Engage_MaskedSpace()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

% Initialize Application and get AP2 and XBUS locks
if(S232('S2init', 0, 'INIT_SECONDARY', 20000) == 0)
    disp('Cannot initialize a secondary process')
    return;
end

if(S232('APlock', 100, 0) == 0)
    disp('Cannot acquire lock on AP2 Card')
    s232('S2close');
    return;
end

if(S232('XBlock', 100, 0) == 0)
    disp('Cannot acquire lock on X Bus')
    s232('APunlock', 0);
    s232('S2close');
    return;
end

S232('trash');

%Engage_MaskedSpace
%*******************************************************************************
%	The MaskedSpace Test operation
% calls a DLL to run masker then probe from many locations
% should use HRIRs with earphone equalization - usually called *.eq
%*******************************************************************************

stimuli_dir = FN.temp_stim_path;
fclose all;

% Check for scaleFactor
scaleFactor = 60;
if ~XStimParams.HiDynamicRange
    BN = questdlg('Use HighDynamic Range?','Dynamic Range is set to low','High','Low (normal)','High');
    if ~isempty(findstr(BN,'Low'))    scaleFactor = TDT.scaleFactor;      end
end

% reset stim_type to BroadBand
set(H.stim_type,'Value',8);

XStimParams.curr_ABL = str2num(get(H.MaskedSpace_ABL,'String'));
XStimParams.curr_stimdur2 = str2num(get(H.MaskedSpace_probeDUR,'String'));
XStimParams.curr_stimdur = str2num(get(H.MaskedSpace_maskDUR,'String'));
XStimParams.numreps = str2num(get(H.MaskedSpace_numreps,'String'));
XStimParams.reset_flag = 0;
XStimParams.ramp = 1;
XStimParams.mask_add_dB = str2num(get(H.MaskedSpace_mask_add_dB,'string'));
XStimParams.probe_add_dB = str2num(get(H.MaskedSpace_probe_add_dB,'string'));

% values for assembling stimuli
probe_pts = XStimParams.curr_stimdur2 *30;
mask_pts = XStimParams.curr_stimdur *30;
nptsTotalPlay = probe_pts + mask_pts;
ramp_flag = XStimParams.ramp;
ramp_usecs = 1000 * 10;
record_spikes = get(H.MaskedSpace_recorddata,'value');
nSpikes = round(nptsTotalPlay/100);
SRATE = 10^6/TDT.Fs;

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

% specify DAMA buffers
clear BUF
BUF.PLAY_SPEC   = 1;
BUF.IB0_SEQ     = 2;
BUF.IB1_SEQ     = 3;
BUF.LHRTF       = 10;
BUF.RHRTF       = 11;
BUF.ENDZEROS    = 12;
BUF.MASK        = 15;
BUF.PROBE       = 16;
BUF.ZEROS       = 17;
BUF.RAMPUP      = 22;
BUF.RAMPDOWN    = 23;

% atten
ABL = abs(XStimParams.curr_ABL);
if(ABL > 110) return; end
atten = ABL-20;

%Get the HRTF index for the single probe location
XStimParams.locations = GUI.locations1';
while size(XStimParams.locations,2) ~=1
    set(H.pickerfig,'Color', [.1 .2 .8]);
    set(H.picker_error,'visible','on');
    input('Enter single probe location')
    XStimParams.locations = GUI.locations1';
    set(H.picker_error,'visible','off');
    set(H.pickerfig,'Color', [.8 .8 .8]);
end

% get HRTF indices for masker location(s)
Probeind = [];       % index of probe loc into array of masker locs
while isempty(Probeind)
    eval(['load ' FN.current_path 'Locations_current;'])
    GUI.locations1 = locations;
    setinfo_spacePicker
    set(H.pickerfig,'Name','Masker Space Picker')
    input('Enter masker locations')
    XStimParams.locations2 = GUI.locations1';
    if ~length(XStimParams.locations2)
        set(H.pickerfig,'Color', [.1 .2 .8]);
        set(H.picker_error,'visible','on');
        input('Enter masker locations')
        XStimParams.locations2 = GUI.locations1';
        set(H.picker_error,'visible','off');
        set(H.pickerfig,'Color', [.8 .8 .8]);
    end
    % check to see that probe location is also one of the masker locations
    Probeind = find(XStimParams.locations2(1,:)==XStimParams.locations(1,1) & XStimParams.locations2(2,:)==XStimParams.locations(2,1));
    if isempty(Probeind)
        input('masker locations must include probe location')
    end
    clc
end

% HRTF file picked?
if XStimParams.space_flag
    while exist1([FN.space_path FN.space_eq]) ~=2
        [FN.space_eq,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select Fully-cued HRTF File *.eq preferred');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space_eq]);
        end
        FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);
    end
    disp('This is a FULLY-CUED MaskedSpace test')
    str1='FC ';
    XStimParams.HRTF_FN = FN.space_eq;
    filetype = FN.HRTFfiletype(1,1);
else
    while exist1([FN.space_path FN.ablequal_eq]) ~=2
        [FN.ablequal_eq,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select ABLequal FC HRTF File  *.eq preferred');
        if(FN.space_path ~= 0)
            set(H.ABLequalfile,'String',[FN.space_path FN.ablequal_eq]);
        end
        FN.HRTFfiletype(7,1) = testHRTFfiletype(FN.space_path, FN.ablequal_eq);
    end
    disp('This is an ABLequal FC MaskedSpace test')
    str1='ABL= ';
    XStimParams.HRTF_FN = FN.ablequal_eq;
    filetype = FN.HRTFfiletype(7,1);
end   

% load HRTFs
if filetype == 1
    dir = sph2dbl(mtlrdir([FN.space_path XStimParams.HRTF_FN]));
    nLocs = size(dir,2);
    for iloc = 1:nLocs
        TF1(iloc,:) = mtlrch([FN.space_path XStimParams.HRTF_FN],iloc*2-1);
        TF2(iloc,:) = mtlrch([FN.space_path XStimParams.HRTF_FN],iloc*2);
    end
elseif filetype == 2
    dir = 0;
    eval(['load -mat ' FN.space_path XStimParams.HRTF_FN]);
else
    disp(['space HRTFfiletype incorrect'])
    return
end

% keep only those HRTFlocs for use
nLocs = size(XStimParams.locations2,2);
clear hrtfdirmat
for iloc = 1:nLocs
    ind = max(find(dir(1,:) == XStimParams.locations2(1,iloc) &...
        dir(2,:) == XStimParams.locations2(2,iloc)));
    TF1_MaskedSpace(iloc,:) = TF1(ind,:);
    TF2_MaskedSpace(iloc,:) = TF2(ind,:);
    hrtfdirmat(:,iloc) = dir(:,ind);
end
clear TF1 TF2 dir
% index to RF location
locnum = find(hrtfdirmat(1,:)== XStimParams.locations(1,1)  & hrtfdirmat(2,:) == XStimParams.locations(2,1));

% load HRTF coefs **************/
TF1 = TF1_MaskedSpace';
TF1 = TF1(:)';
TF2 = TF2_MaskedSpace';
TF2 = TF2(:)';
S232('dropall');
S232('pushf',TF1, 255*nLocs);
S232('dpush',255);
S232('value',0.0);
S232('cat');
S232('allotf', BUF.LHRTF, 255*(nLocs+1));
S232('qpopf', BUF.LHRTF );

S232('dropall');
S232('pushf',TF2, 255*nLocs);
S232('dpush',255);
S232('value',0.0);
S232('cat');
S232('allotf', BUF.RHRTF, 255*(nLocs+1));
S232('qpopf', BUF.RHRTF );
nLocs = nLocs +1;           % added on a zeros HRTF to do noG condition
hrtfdirmat(:,nLocs) = [ 1 1 ]';
clear TF*

% save XStimParams for this test
if get(H.MaskedSpace_recorddata,'Value')
    tempStruct.record_spikes = 1;
    tempstr = ['    ' str1 'MaskedSpace     atten: ' num2str(abs(ABL)) 'dB'];
    update_diary
else
    tempStruct.record_spikes = 0;
end

Temp_params = XStimParams;
eval(['save ' FN.temp_stim_path 'XStimParams_MaskedSpace_' str1 ' Temp_params;'])
clear Temp_params str1

% update display
setInfo_MaskedSpace;

%Begin playing sounds   
set(H.MaskedSpace_status,'String','Status: playing stimuli');
set(H.MaskedSpace_status,'BackgroundColor','green');
set(H.MaskedSpace_remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
buffcycle = 1;

% increment testnumber
if(exist1('H.MaskedSpacefig') & get(H.MaskedSpace_recorddata,'Value'))
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.MaskedSpace_FN,'String',FN.data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up buffers
S232('allot16',BUF.PLAY_SPEC,10);
S232('allot16',BUF.IB0_SEQ,10);
S232('allot16',BUF.IB1_SEQ,10);

% play specification list
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.IB0_SEQ);
S232('make',1,BUF.IB1_SEQ);
S232('make',2,0);
S232('qpop16',BUF.PLAY_SPEC);

% playsequence for IB0 (MASK)
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.MASK);
S232('make',1,1);
S232('make',2,BUF.ZEROS);
S232('make',3,probe_pts/30);
S232('make',4,BUF.ENDZEROS);
S232('make',5,1);
S232('make',6,0);
S232('qpop16',BUF.IB0_SEQ);

% playsequence for IB1 (PROBE)
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.ZEROS);
S232('make',1,mask_pts/30);
S232('make',2,BUF.PROBE);
S232('make',3,1);
S232('make',4,BUF.ENDZEROS);
S232('make',5,1);
S232('make',6,0);
S232('qpop16',BUF.IB1_SEQ);

% ZEROS	***************/
S232('allot16', BUF.ZEROS, 30);
S232('dpush',30);
S232('value',0);
S232('qpop16', BUF.ZEROS );

% ENDZEROS   *************/
S232('allot16', BUF.ENDZEROS, 254);
S232('dpush',254);
S232('value',0);
S232('qpop16', BUF.ENDZEROS );

% MASKER   ***************/
S232('dropall');
if (ramp_flag)
    if(mask_pts * SRATE >= ramp_usecs*2)
        ramppts = floor(ramp_usecs/SRATE);
        S232('dpush',ramppts);
        S232('fill',0,1/ramppts);
        S232('dpush',mask_pts-ramppts);
        S232('value',1);
        S232('cat');
        
        S232('dpush',mask_pts-ramppts);
        S232('value',1);
        S232('dpush',ramppts);
        S232('fill',1,-1/ramppts);
        S232('cat');
        S232('mult');
    else
        ramppts = mask_pts/2;
        S232('dpush',ramppts);
        S232('fill',0,1/ramppts);
        
        S232('dpush',ramppts);
        S232('fill',1,-1/mask_pts);
        S232('cat');
    end
end

S232('dpush',mask_pts);
S232('value',0);

S232('flat');		% flat noise
if (ramp_flag)
    S232('mult');		% multiply with ramp
end
S232('scale',scaleFactor * (10^(XStimParams.mask_add_dB /20)));
S232('allot16', BUF.MASK, mask_pts);
S232('qpop16',BUF.MASK);

% PROBE ****************/
S232('dropall');
if (ramp_flag)
    if(probe_pts * SRATE >= ramp_usecs*2)
        ramppts = floor(ramp_usecs/SRATE);
        S232('dpush',ramppts);
        S232('fill',0,1/ramppts);
        S232('dpush',probe_pts-ramppts);
        S232('value',1);
        S232('cat');
        
        S232('dpush',probe_pts-ramppts);
        S232('value',1);
        S232('dpush',ramppts);
        S232('fill',1,-1/ramppts);
        S232('cat');
        S232('mult');
    else
        ramppts = probe_pts/2;
        S232('dpush',ramppts);
        S232('fill',0,1/ramppts);
        
        S232('dpush',ramppts);
        S232('fill',1,-1/probe_pts);
        S232('cat');
    end
end

S232('dpush',probe_pts);
S232('value',0);
S232('flat');		    % flat noise
if (ramp_flag)          % multiply with ramp
    S232('mult');	
end

S232('scale',scaleFactor * 10^(XStimParams.probe_add_dB /20));
S232('allot16', BUF.PROBE, probe_pts );
S232('qpop16', BUF.PROBE );

% setup PD1	****************/
S232('PD1clear',1);
S232('PD1srate',1,SRATE);
S232('PD1npts',1,nptsTotalPlay+254);

S232('PD1resetDSP',1,hex2dec('FFF'));
S232('PD1clrsched',1);
S232('PD1nstrms',1, 2, 0);

sources(1) = S232('DSPout',0);          % DSP[0] & DSP[2] to DAC[0] LEFT
sources(2) = S232('DSPout',2);
scales(1) = 1.0;
scales(2) = 1.0;
S232('PD1addmult',1, sources, scales, 2, S232('DAC',0));

sources(1) = S232('DSPout',1);          % DSP[1] & DSP[3] to DAC[1] RIGHT
sources(2) = S232('DSPout',3);
scales(1) = 1.0;
scales(2) = 1.0;
S232('PD1addmult',1, sources, scales, 2, S232('DAC',1));

S232('PD1addsimp',1, S232('IREG',0), S232('DSPin',0));		% IREG[0] to DSP[0]	LEFT
S232('PD1addsimp',1, S232('IREG',0), S232('DSPin',1));		% IREG[0] to DSP[1]	RIGHT
S232('PD1specIB',1, S232('IB',0),   S232('IREG',0));		% MASKER IB[0] to IREG[0]

S232('PD1addsimp',1, S232('IREG',1), S232('DSPin',2));		% IREG[1] to DSP[2]	LEFT
S232('PD1addsimp',1, S232('IREG',1), S232('DSPin',3));		% IREG[1] to DSP[3]	RIGHT
S232('PD1specIB',1, S232('IB',1),   S232('IREG',1));		% PROBE IB[1] to IREG[1]

% preload DSPs with PROBE HRTFs
% locum = index to probe location
S232('dropall');
S232('qpushpartf',BUF.LHRTF,(locnum-1)*255+1,255);
S232('PreLoadRaw',1, S232('DSPid',2), 'MONO', 'STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);

S232('dropall');
S232('qpushpartf',BUF.RHRTF,(locnum-1)*255+1,255);
S232('PreLoadRaw',1, S232('DSPid',3), 'MONO', 'STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);

m100x( C_.INIT );
%if (record_spikes)
%    m110dx( C_.INIT );
%    m110dx( C_.CLOCK, 10);
%    m110dx( C_.MODE, M_.PST );
%end

% set attenuation
S232('PA4atten',1,atten);
S232('PA4atten',2,atten);

datamatrix = [];
finalspikematrix = zeros(1,nLocs);

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
%loop for reps
for iRep = 1:XStimParams.numreps
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    %Randomize the stimuli
    randseq = randperm(nLocs);
    tempdir = hrtfdirmat(:,randseq);
    
    %Check for pause by user
    if pause_check    return; end
    
    %initialize spikes array
    spikes = zeros(nLocs,nSpikes);
    tic;
    
    %%%%%%%%%%%%%%%%%%%%%%% loop for locations
    for iLoc = 1:nLocs+1
        if (iLoc <= nLocs)
            if randseq(iLoc) == nLocs
                S232('dpush',255);
                S232('value',0.0);
                S232('dpush',255);
                S232('value',0.0);
            else
                hrtfpos = (randseq(iLoc)-1)*255+1;
                S232('dropall');
                S232('qpushpartf',BUF.RHRTF,hrtfpos,255);
                S232('qpushpartf',BUF.LHRTF,hrtfpos,255);
            end
        end
        
        %%%%%%% delay
        while (toc < (.15+nptsTotalPlay/TDT.Fs))   end
        
        if (iLoc>1)
            S232('PD1stop',1);
            m110dx( C_.STOP);
            spikes(iLoc,:) = m110dx( C_.DATA, nSpikes)'; 
        end
        if pause_check    return; end
        
        if (iLoc <= nLocs)
            S232('PreLoadRaw',1, S232('DSPid',0), 'MONO', 'STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);
            S232('PreLoadRaw',1, S232('DSPid',1), 'MONO', 'STACK','','',TDT.hrtf_scale,TDT.hrtf_scale,1);
            tic
            S232('seqplay',BUF.PLAY_SPEC);
            m110dx( C_.INIT );
            m110dx( C_.CLOCK, 10);
            m110dx( C_.MODE, M_.PST );
            S232('PD1arm',1);
            m110dx( C_.START );
            m101x( C_.DATA,M_.BIT,M_.PULSE,0 ); % Use port 0 for the pulse
        end	
        
    end	
    %%%%%%%%%%%%%%%%%%%%%%
    
    % process spikes array
    for iLoc = 1:nLocs
        temp = spikes(iLoc,1:max1(find(spikes(iLoc,:)~=0)))';    % strip off excess
        nspikes =0;
        if ~isempty(temp)
            ind0 = [1; (find(abs(diff(temp)) > mii_separation))+1];
            temp = temp(ind0); nspikes = length(temp);
            spikes_trial = [spikes_trial;temp/(1000/mii_us_per_sample)];
            EL_trial = [EL_trial;hrtfdirmat(1,randseq(iLoc))* ones(size(temp))];
            AZ_trial = [AZ_trial;hrtfdirmat(2,randseq(iLoc))* ones(size(temp))];
            repnum_trial = [repnum_trial;iRep * ones(size(temp))];
            Nspikes = [Nspikes; length(temp) * ones(size(temp))];
        end
        finalspikematrix(randseq(iLoc)) = finalspikematrix(randseq(iLoc)) + nspikes;
    end
    
    if pause_check    return; end
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.MaskedSpacefig') & ~exist1('H.MaskedSpace_finalspikeratefig'))
        H.MaskedSpace_finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','MaskedSpace Test Spike Rate Plot',...
            'NumberTitle','off');
        H.MaskedSpace_spikeaxes = axes;
    end
    figure(H.MaskedSpace_finalspikeratefig)
    plotdiam1(hrtfdirmat, interimspikerate);
    %set(H.MaskedSpace_spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(iRep)]);
    colorbar
    
    %Record Data
    if(exist1('H.MaskedSpacefig') & get(H.MaskedSpace_recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    set(H.MaskedSpace_remreps,'String',num2str(XStimParams.numreps - iRep));
    pause(0);
end 									%end loop over reps

S232('PA4mute',1);
S232('PA4mute',2);		

S232('trash');
S232('dropall');

S232('APunlock',0);
S232('XBunlock',0);
S232('S2close');


%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.MaskedSpace_finalspikeratefig)
set(H.MaskedSpace_finalspikeratefig,'Name','Final Plot for MaskedSpace Test');
plotdiam1(hrtfdirmat, interimspikerate);
set(H.MaskedSpace_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(hrtfdirmat(1,locmaxspikes)) ...
        ', AZ = ' num2str(hrtfdirmat(2,locmaxspikes))], 'FontSize',8);
colorbar

set(H.MaskedSpace_status,'String','Status: results');
set(H.MaskedSpace_status,'BackgroundColor','blue');

set(H.exitMaskedSpace,'Visible','on');
set(H.resetMaskedSpace,'Visible','on');

% increment test number
if(exist1('H.MaskedSpacefig') & get(H.MaskedSpace_recorddata,'Value'))
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
while (exist1('H.MaskedSpacefig') & get(H.pauseMaskedSpace,'Value'))
    pause(0);
    if(~exist1('H.MaskedSpacefig')) return; end         
    set(H.exitMaskedSpace,'Visible','on');
    set(H.resetMaskedSpace,'Visible','on');
    if(exist1('H.MaskedSpacefig') & get(H.resetMaskedSpace,'Value') == 1)
        set(H.resetMaskedSpace,'Value',0);
        set(H.pauseMaskedSpace,'Value',0);
        Reset_MaskedSpace;    flag = 1;
        return;
    end
    if isempty(XStimParams.locations)
        Reset_MaskedSpace;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end