function [] = Engage_MaskedSpace()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

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

% make structure for calling DLL
clear tempStruct
tempStruct.samplingPeriod = 33.3;
tempStruct.probe_pts = XStimParams.curr_stimdur2 *30;
tempStruct.mask_pts = XStimParams.curr_stimdur *30;
tempStruct.noiScale = scaleFactor;
tempStruct.nptsTotalPlay = tempStruct.probe_pts + tempStruct.mask_pts;
tempStruct.ramp_flag = XStimParams.ramp;
tempStruct.mask_add_dB = XStimParams.mask_add_dB;
tempStruct.probe_add_dB = XStimParams.probe_add_dB;

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

% atten
ABL = abs(XStimParams.curr_ABL);
if(ABL > 110) return; end
tempStruct.latten = ABL-20;
tempStruct.ratten = ABL-20;

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
tempStruct.nLocs = size(XStimParams.locations2,2);
clear hrtfdirmat
for iloc = 1:tempStruct.nLocs
    ind = max(find(dir(1,:) == XStimParams.locations2(1,iloc) &...
        dir(2,:) == XStimParams.locations2(2,iloc)));
    TF1_MaskedSpace(iloc,:) = TF1(ind,:);
    TF2_MaskedSpace(iloc,:) = TF2(ind,:);
    hrtfdirmat(:,iloc) = dir(:,ind);
end
clear TF1 TF2 dir

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
datamatrix = [];
finalspikematrix = zeros(1,tempStruct.nLocs);

% increment testnumber
if(exist1('H.MaskedSpacefig') & get(H.MaskedSpace_recorddata,'Value'))
    update_dataFN;
end

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
%loop for reps
while (exist1('H.MaskedSpacefig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli
    randseq = randperm(tempStruct.nLocs);
    TF1 = TF1_MaskedSpace(randseq,:)';
    TF1 = single(TF1(:)');
    TF2 = TF2_MaskedSpace(randseq,:)';
    TF2 = single(TF2(:)');
    tempdir = hrtfdirmat(:,randseq);
    tempStruct.locnum = find(tempdir(1,:)== XStimParams.locations(1,1)  & tempdir(2,:) == XStimParams.locations(2,1));
        
    %Check for pause by user
    if pause_check    return; end
    play2_record2b_rove3(tempStruct,TF1,TF2);
    
    %Stop the m110 and get spikes
    m110dx_old( C_.STOP);
    spikes = m110dx_old( C_.DATA, XStimParams.curr_stimdur*tempStruct.nLocs); 			% Take 2*XStimParams.curr_stimdur spikes max
    % stim onsets denoted by -1 -1 -1 -1
    spikes = spikes(1:max1(find(spikes~=0))+1);    % strip off excess
    badspikes=[];
    i = 1;
    while i <= length(spikes) && spikes(i)~=0
        test = 0;
        while spikes(i) == -1
            test=test+1;
            i=i+1;
        end
        switch mod(test,4)
            case 0
            case 1
                badspikes = [badspikes; i-1];
            case 2
                badspikes = [badspikes; i-1; i-2];
            case 3
                badspikes = [badspikes; i-1; i-2; i-3];
        end
        while spikes(i)>0
            i=i+1;
        end
    end
   
    spikes = spikes(setdiff(1:length(spikes)-1,badspikes));
    ind = find(spikes==-1);
    ind= ind(4:4:length(ind));
    
    temp = tempStruct.nLocs+1-length(ind);
    if temp>0
        ind = [ind; length(spikes)];
        spikes = [spikes; -ones(temp,1)];
    elseif (temp<0)
        error('too many -1 in spikes');
    end
    
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    ispike = 1;

    for iloc = 1:tempStruct.nLocs
        temp = spikes(ind(iloc)+1:ind(iloc+1)-1);
        if ~isempty(temp)
            ind0 = [1; (find(abs(diff(temp)) > mii_separation))+1];
            nspikes = length(ind0);
            spikes_trial(ispike:ispike+nspikes-1,1) = temp(ind0) /(1000/mii_us_per_sample);
            EL_trial(ispike:ispike+nspikes-1,1) = ones(nspikes,1) * tempdir(1,iloc);
            AZ_trial(ispike:ispike+nspikes-1,1) = ones(nspikes,1) * tempdir(2,iloc);
            repnum_trial(ispike:ispike+nspikes-1,1) = ones(nspikes,1) * repnum;
            Nspikes(ispike:ispike+nspikes-1,1) = ones(nspikes,1) * nspikes;
            ispike = ispike+nspikes;
        else
            nspikes = 0;
        end
        finalspikematrix(randseq(iloc)) = finalspikematrix(randseq(iloc)) + nspikes;
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
    plotdiam1(XStimParams.locations2, interimspikerate);
    %set(H.MaskedSpace_spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    %Record Data
    if(exist1('H.MaskedSpacefig') & get(H.MaskedSpace_recorddata,'Value'))
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
        record_data3(XStimParams,datamatrix);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.MaskedSpace_remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.MaskedSpace_finalspikeratefig)
set(H.MaskedSpace_finalspikeratefig,'Name','Final Plot for MaskedSpace Test');
plotdiam1(XStimParams.locations2, interimspikerate);
set(H.MaskedSpace_spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(XStimParams.locations2(1,locmaxspikes)) ...
        ', AZ = ' num2str(XStimParams.locations2(2,locmaxspikes))], 'FontSize',8);
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