% SetInfo_space_FPR
% combines all the other 'set...' scripts into one

% uses .STD HRIRs, therefore make sure earphone filters are also selected
   XStimParams.ephone_flag = 1;
   set(H.ephoneuseit,'Value',1);

 
if get(H.space_FPR_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.space_FPR_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.space_FPR_DynamicRangepb,'BackgroundColor','red');
else
    XStimParams.HiDynamicRange = 0;
    set(H.space_FPR_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.space_FPR_DynamicRangepb,'BackgroundColor','yellow');
end
   
XStimParams.buildOnly = 0;
   
if get(H.spaceflag,'Value') & ~XStimParams.space_flag
    XStimParams.space_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
elseif get(H.space_ABLaloneflag,'Value') & ~XStimParams.ABLalone_flag
    XStimParams.space_flag = 0;
    XStimParams.ABLalone_flag = 1;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 0;
elseif get(H.space_ildaloneflag,'Value') & ~XStimParams.ildalone_flag
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 1;
    XStimParams.itdalone_flag = 0;
    XStimParams.ABLalone_flag = 0;
elseif get(H.space_itdaloneflag,'Value') & ~XStimParams.itdalone_flag
    XStimParams.space_flag = 0;
    XStimParams.ildalone_flag = 0;
    XStimParams.itdalone_flag = 1;
    XStimParams.ABLalone_flag = 0;
end
set(H.spaceflag,'Value',XStimParams.space_flag);
set(H.space_ildaloneflag,'value',XStimParams.ildalone_flag);
set(H.space_itdaloneflag,'value',XStimParams.itdalone_flag);
set(H.space_ABLaloneflag,'value',XStimParams.ABLalone_flag);

% control visibility of file picking for ILA, ITA or fully-cued
if get(H.spaceflag,'Value')
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
    set(H.space_ildalonefilehdr,'Enable','off');
    set(H.space_ildalonefile,'Enable','off');
    set(H.space_ildalonepb,'Enable','off');  
    set(H.space_itdalonefilehdr,'Enable','off');
    set(H.space_itdalonefile,'Enable','off');
    set(H.space_itdalonepb,'Enable','off');  
    set(H.space_FPR_ITDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.space_FPR_ILDtxt,'BackgroundColor',[.8 .8 .8]);
elseif get(H.space_ABLaloneflag,'Value')
    set(H.spacefilehdr,'Enable','on');
    set(H.spacefile,'Enable','on');
    set(H.spacefilepb,'Enable','on');   
    set(H.space_ildalonefilehdr,'Enable','off');
    set(H.space_ildalonefile,'Enable','off');
    set(H.space_ildalonepb,'Enable','off');  
    set(H.space_itdalonefilehdr,'Enable','off');
    set(H.space_itdalonefile,'Enable','off');
    set(H.space_itdalonepb,'Enable','off');  
    set(H.space_FPR_ITDtxt,'BackgroundColor','green');
    set(H.space_FPR_ILDtxt,'BackgroundColor','green');
elseif get(H.space_ildaloneflag,'Value')
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.space_ildalonefilehdr,'Enable','on');
    set(H.space_ildalonefile,'Enable','on');
    set(H.space_ildalonepb,'Enable','on');
    set(H.space_itdalonefilehdr,'Enable','off');
    set(H.space_itdalonefile,'Enable','off');
    set(H.space_itdalonepb,'Enable','off');  
    set(H.space_FPR_ITDtxt,'BackgroundColor','green');
    set(H.space_FPR_ILDtxt,'BackgroundColor',[.8 .8 .8]);
elseif get(H.space_itdaloneflag,'Value')
    set(H.spacefilehdr,'Enable','off');
    set(H.spacefile,'Enable','off');
    set(H.spacefilepb,'Enable','off');
    set(H.space_ildalonefilehdr,'Enable','off');
    set(H.space_ildalonefile,'Enable','off');
    set(H.space_ildalonepb,'Enable','off');
    set(H.space_itdalonefilehdr,'Enable','on');
    set(H.space_itdalonefile,'Enable','on');
    set(H.space_itdalonepb,'Enable','on');  
    set(H.space_FPR_ITDtxt,'BackgroundColor',[.8 .8 .8]);
    set(H.space_FPR_ILDtxt,'BackgroundColor','green');
end

% check for ILA pushbutton (file-picker)
if get(H.space_ildalonepb,'Value')
    [FN.ildalone_std,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
    if(FN.ILA_path ~= 0)
        set(H.space_ildalonefile,'String',FN.ildalone_std);
    end
    set(H.space_ildalonepb,'Value',0);
    FN.HRTFfiletype(2,2) = testHRTFfiletype(FN.ILA_path, FN.ildalone_std);
end

% check for ITA pushbutton (file-picker)
if get(H.space_itdalonepb,'Value')
    [FN.itdalone_std,FN.ITA_path] = uigetfile([FN.ITA_path '*.*'],'Select ITD-alone HRTF File');
    if(FN.ITA_path ~= 0)
        set(H.space_itdalonefile,'String',FN.itdalone_std);
    end
    set(H.space_itdalonepb,'Value',0);
    FN.HRTFfiletype(3,2) = testHRTFfiletype(FN.ITA_path, FN.itdalone_std);
end

% check for fully-cued pushbutton (file-picker)
if get(H.spacefilepb,'Value')
    [FN.space_std,FN.space_path] = uigetfile([FN.HRTF_path '*.*'],'Select Fully-cued HRTF File (*.std preferred)');
    if(FN.space_path ~= 0)
        set(H.spacefile,'String',FN.space_std);
    end
    set(H.spacefilepb,'Value',0);
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end

% set reference locations for ILA
if get(H.spaceflag,'Value') | get(H.space_ABLaloneflag,'Value') | get(H.space_itdaloneflag,'Value')
    XStimParams.az = str2num(get(H.space_FPR_az,'String'));
    XStimParams.el = str2num(get(H.space_FPR_el,'String'));
end

% reset curr_ILD and curr_ITD
XStimParams.curr_ITD = str2num(get(H.space_FPR_ITD,'string'));
XStimParams.curr_ABL = str2num(get(H.space_FPR_ABL,'string'));

% set numreps
XStimParams.numreps = str2num(get(H.space_FPR_numreps,'String'));

%Set space_FPR PickerFlag
XStimParams.picker_flag = 1;
if exist1('H.pickerfig') 	
    if isempty(H.pickerfig)
        XStimParams.locations = [];
        H.pickerfig = space_Picker;
    end
else
    XStimParams.locations = [];
    H.pickerfig = spacePicker;  
end

% increment test number
if get(H.space_FPR_inc_testnum,'Value')
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.space_FPR_inc_testnum,'Value',0);
end
set(H.space_FPRdata_FN,'String', FN.data);


%%%%%%%%%%%%%%%%%%
% initialize FNs
if get(H.space_FPR_initFN,'value')
    FN.FPR = [];
    for i = 1:4    
        set(H.space_FPR_stimFN(i),'string','')        
    end

    [temp,FN.stim_path2] = uigetfile([FN.stim_path2 '*.*'],'Select stimulus PATH');
    
    % get directory of FN.stim_path2
    D = dir(FN.stim_path2);
    nFiles = length(D)-2;
    
    for iFile = 1:nFiles
        tempFN{iFile} = D(iFile+2).name;
    end
    
    % choose which stimFNs to include
    if exist('tempFN')
        ind = listdlg('PromptString','Select DATA Files',...
            'ListString', tempFN,...
            'InitialValue',[1:nFiles]);
        clear stimFN
        for i = 1:length(ind)
            FN.FPR{i} = tempFN{ind(i)};  % re-organized in engage for each rep
        end
    end
    nFiles = length(FN.FPR);
    % update display
    set(H.space_FPR_stimFN(1),'string',strvcat(FN.FPR(1:min1([nFiles,6]))))
    if nFiles>6
        set(H.space_FPR_stimFN(2),'string',strvcat(FN.FPR(7:min1([nFiles,12]))))
    end
    if nFiles>12
        set(H.space_FPR_stimFN(3),'string',strvcat(FN.FPR(13:min1([nFiles,18]))))
    end
    if nFiles>18
        set(H.space_FPR_stimFN(4),'string',strvcat(FN.FPR(19:min1([nFiles,24]))))
    end
    
    for iFile=1:nFiles
        D = dir([FN.stim_path2 FN.FPR{iFile}]);
        B(iFile) = D.bytes/4;
    end
    DUR = max1(B) / round(TDT.Fs/1000);
    XStimParams.curr_stimdur = DUR;

    clear ind tempFN D nFiles temp
    % update stimulus duration
    set(H.space_FPR_DUR, 'String',num2str(XStimParams.curr_stimdur));
end


eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim