% SetInfo_Composite
stimuli_dir = FN.temp_stim_path;

% ClearStim Directory
if get(H.composite_ClearStimspb,'value')
    eval(['delete ' stimuli_dir '*.*;']);    
end

% initialize stimulus FNs
if get(H.composite_initFNpb,'value')
    set(H.composite_status,'string','Initializing FNs');
    set(H.composite_initFNpb,'value',0);
    % get directory of MRGpath
    D = dir(stimuli_dir);
    nFiles = length(D);
    if nFiles <=2       % only directories, no stim files
        disp('no files available')
    else
        % find left ear files (assume right ear files also exist)
        clear FNlist tempFN
        ind0=[]; n =0;
        for iFile = 3:nFiles
            if findstr('left',D(iFile).name)
                n = n+1;    ind0 = [ind0 iFile];
                tempFN{n} = D(iFile).name;
            end
        end
        
        % choose which datafiles to include
        if exist('tempFN')
            ind = listdlg('PromptString','Select DATA Files',...
                'ListString', tempFN,...
                'InitialValue',[1:n]);
            for n=1:length(ind)
                FNlist(n).name = deal(D(ind0(ind(n))).name);
            end
            nFiles = length(FNlist);
        end
        
        % reset text to indicate # stims chosen
        set(H.composite_num_stims_txt,'string',[num2str(nFiles) ' stimuli chosen']);
        
        % determine longest stim duration
        max_dur = max1([D(ind0(ind)).bytes])/2;
        min_dur = min1([D(ind0(ind)).bytes])/2;
        
        % identify in FNlist, which type of file each one is
        for iFile = 1:nFiles
            if findstr(FNlist(iFile).name,'ABL')  FNlist(iFile).type = 1;   end 
            if findstr(FNlist(iFile).name,'ITD')  FNlist(iFile).type = 2;   end 
            if findstr(FNlist(iFile).name,'ILD')  FNlist(iFile).type = 3;   end 
            if findstr(FNlist(iFile).name,'space1_')        % type 4-6
                if findstr(FNlist(iFile).name,'FC')
                    FNlist(iFile).type = 4;
                elseif findstr(FNlist(iFile).name,'ILA')
                    FNlist(iFile).type = 5;
                elseif findstr(FNlist(iFile).name,'ITA')
                    FNlist(iFile).type = 6;
                end
            end 
            if findstr(FNlist(iFile).name,'space2')         % type 7-9
                if findstr(FNlist(iFile).name,'FC')
                    FNlist(iFile).type = 7;
                elseif findstr(FNlist(iFile).name,'ILA')
                    FNlist(iFile).type = 8;
                elseif findstr(FNlist(iFile).name,'ITA')
                    FNlist(iFile).type = 9;
                end
            end 
            if findstr(FNlist(iFile).name,'altIR')          % type 10-12  for now.
                if findstr(FNlist(iFile).name,'FC')
                    FNlist(iFile).type = 10;
                elseif findstr(FNlist(iFile).name,'ILA')
                    FNlist(iFile).type = 11;
                elseif findstr(FNlist(iFile).name,'ITA')    % not currently available
                    FNlist(iFile).type = 12;
                end
            end 
        end
        clear ind ind0 n D nFiles tempFN
    end
    if exist('max_dur')
        set(H.composite_DUR, 'String', num2str(max_dur/round(TDT.Fs/1000)))
    end
end

XStimParams.buildOnly = 1;

% uses .STD HRIRs, therefore make sure earphone filters are also required
XStimParams.ephone_flag = 1;
set(H.ephoneuseit,'Value',1);

% set numreps
XStimParams.numreps = str2num(get(H.composite_numreps,'String'));

% increment test number
if(exist1('H.composite_fig') & get(H.composite_inc_testnum,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.composite_inc_testnum,'Value',0);
end
set(H.xstim_recorddata_FN, 'String', FN.data);

% update stimulus duration
XStimParams.curr_stimdur = str2num(get(H.composite_DUR, 'String'));

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])

% choose test types
if get(H.composite_test_type,'value')~=1;
    % location on GUI
    F.str_loc = F.str_loc + F.ind_ht;
    % number of this test
    Test_num = Test_num+1;
    
    temp = get(H.composite_test_type,'String');
    Temp_type = deblank(temp(get(H.composite_test_type,'Value'),:));
    % run appropriate GUI
    switch Temp_type
        case 'None'
        case 'ABL'
            abltest;
        case 'ITD'
            itdtest;
        case 'ILD'
            ildtest;
        case 'Space'
            space;
        case 'Space2'
            space2;
        case 'AltIR'
            alt_IR_XXX;  
    end
        
   set(H.composite_test_type,'value',1)
end

if XStimParams.reset_flag
        % determine test type
    XStimParams.test_type = Temp_type;
    str1 = [Temp_type ' '];
    switch Temp_type
        case {'Space', 'Space2', 'AltIR'}
            if XStimParams.ildalone_flag
                str1 = [str1 'ILA'];
            elseif XStimParams.itdalone_flag
                str1 = [str1 'ITA'];
            else
                str1 = [str1 'FC'];
            end
    end
    % save dataFN, params...
    switch deblank(str1)
        case {'ABL', 'abl'}
            Test_array(Test_num).type = 1;
            eval(['load ' FN.temp_stim_path 'XStimParams_ABL'])
        case {'ITD', 'itd'}
            Test_array(Test_num).type = 2;
            eval(['load ' FN.temp_stim_path 'XStimParams_ITD'])
        case {'ILD', 'ild'}
            Test_array(Test_num).type = 3;
            eval(['load ' FN.temp_stim_path 'XStimParams_ILD'])
        case {'Space FC', 'space FC'}
            Test_array(Test_num).type = 4;
            eval(['load ' FN.temp_stim_path 'XStimParams_Space_fc'])
        case {'Space ILA', 'space ILA'}
            Test_array(Test_num).type = 5;
            eval(['load ' FN.temp_stim_path 'XStimParams_Space_ila'])
        case {'Space ITA', 'space ITA'}
            Test_array(Test_num).type = 6;
            eval(['load ' FN.temp_stim_path 'XStimParams_Space_ita'])
        case {'Space2 FC', 'space2 FC'}
            Test_array(Test_num).type = 7;
            eval(['load ' FN.temp_stim_path 'XStimParams_Space2_fc'])
        case {'Space2 ILA', 'space2 ILA'}
            Test_array(Test_num).type = 8;
            eval(['load ' FN.temp_stim_path 'XStimParams_Space2_ila'])
        case {'Space2 ITA', 'space2 ITA'}
            Test_array(Test_num).type = 9;
            eval(['load ' FN.temp_stim_path 'XStimParams_Space2_ita'])
        case {'AltIR FC', 'altIR FC'}
            Test_array(Test_num).type = 10;
            eval(['load ' FN.temp_stim_path 'XStimParams_AltIR_fc'])
        case {'AltIR ILA', 'altIR ILA'}
            Test_array(Test_num).type = 11;
            eval(['load ' FN.temp_stim_path 'XStimParams_AltIR_ila'])
        case {'AltIR ITA', 'altIR ITA'}
            Test_array(Test_num).type = 12;
            eval(['load ' FN.temp_stim_path 'XStimParams_AltIR_ita'])
    end
    Test_array(Test_num).params = Temp_params;
    %Test_array(Test_num).params.locations = GUI.locations1';
    Test_array(Test_num).repnum = [];
    Test_array(Test_num).stimnum = [];
    Test_array(Test_num).spikes = [];
    Test_array(Test_num).Nspikes = [];
    Test_array(Test_num).param1 = [];
    Test_array(Test_num).param2 = [];
    Test_array(Test_num).FN = FN.data;
    clear Temp_params
    
% display test type and dataFN
    H.composite_test_type_txt(Test_num) = uicontrol('Parent',H.composite_fig,...
        'Style','text',...
        'Units','pixels',...
        'ForegroundColor','blue',...
       'Position',[30 F.composite_ht-F.str_loc 80 F.ind_ht],...
        'String', str1);
    H.composite_dataFN(Test_num) = uicontrol('Parent',H.composite_fig,...
        'Style','text',...
        'Units','pixels',...
        'ForegroundColor','blue',...
       'Position',[120 F.composite_ht-F.str_loc 100 F.ind_ht],...
        'String', FN.data);
    % increment testnum
    if(exist1('H.composite_fig'))
        XStimParams.testnum = XStimParams.testnum +1;
        set(H.testnum, 'String', num2str(XStimParams.testnum));
        update_dataFN;
        set(H.xstim_recorddata_FN, 'String', FN.data);
    end
    clear Temp_type
    XStimParams.reset_flag = 0;
end