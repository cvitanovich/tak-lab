% SetInfo_space_RIR
% combines all the other 'set...' scripts into one

% uses .STD HRIRs, therefore make sure earphone filters are also selected
XStimParams.ephone_flag = 1;
set(H.ephoneuseit,'Value',1);


if get(H.space_RIR_DynamicRangepb,'value')
    XStimParams.HiDynamicRange = 1;
    set(H.space_RIR_DynamicRangepb,'string','push for Normal Dynamic Range');
    set(H.space_RIR_DynamicRangepb,'BackgroundColor','red');
else
    XStimParams.HiDynamicRange = 0;
    set(H.space_RIR_DynamicRangepb,'string','push for High Dynamic Range');
    set(H.space_RIR_DynamicRangepb,'BackgroundColor','yellow');
end

XStimParams.buildOnly = 0;


XStimParams.space_flag = 1;
XStimParams.ildalone_flag = 0;
XStimParams.itdalone_flag = 0;
XStimParams.ABLalone_flag = 0;

% set numreps
XStimParams.numreps = str2num(get(H.space_RIR_numreps,'String'));

% increment test number
if get(H.space_RIR_inc_testnum,'Value')
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.space_RIR_inc_testnum,'Value',0);
end
set(H.space_RIRdata_FN,'String', FN.data);


%%%%%%%%%%%%%%%%%%
% use old stimFNs
if get(H.space_RIR_useSavedFN,'value')
    str1 = 'rir';
    eval(['load ' FN.temp_stim_path 'XStimParams_space_RIR_' str1 ' Temp_FN_RIR'])
    FN.RIR = Temp_FN_RIR;
end
% initialize FNs
if get(H.space_RIR_initFN,'value')
    FN.RIR = [];
    set(H.space_RIR_stimFN,'string','')        
    
    [temp,FN.stim_path2] = uigetfile([FN.stim_path2 '*.*'],'Select stimulus PATH');
    
    % get directory of FN.stim_path2
    D = dir(FN.stim_path2);
    nFiles = length(D)-2;
    
    for iFile = 1:nFiles
        tempFN{iFile} = D(iFile+2).name;
    end
    
    % choose which stimFNs to include
    if exist('tempFN')
        ind = listdlg('PromptString','Select DATA Files', ...
            'ListSize', [400 500], ...
            'ListString', tempFN, ...
            'InitialValue',[1:nFiles]);
        clear stimFN
        for i = 1:length(ind)
            FN.RIR{i} = tempFN{ind(i)};  % re-organized in engage f
        end
    end
end

if get(H.space_RIR_useSavedFN,'value') | get(H.space_RIR_initFN,'value')
    
    nFiles = length(FN.RIR);
    clear repNum
    XStimParams.locations1 = [];
    for ifile = 1:nFiles
        fname = char(FN.RIR(ifile));
        ind = strfind(fname,'rep')+3;
        repNum(ifile) = str2num(fname(ind:ind+1));
        
        % update location matrix
        ind = findstr(fname,'el')-4;
        ele = str2num(fname(ind+1:ind+3));
        if strcmp(fname(ind),'n')   ele=-ele;   end
        ind = findstr(fname,'az')-4;
        azi = str2num(fname(ind+1:ind+3));
        if strcmp(fname(ind),'n')   azi = -azi; end
        
        if isempty(XStimParams.locations1)
            XStimParams.locations1 = [ele; azi];
        else
            if isempty(find(XStimParams.locations1(1,:)==ele & XStimParams.locations1(2,:)==azi))
                XStimParams.locations1 = [XStimParams.locations1 [ele; azi]];
            end
        end
    end
    [repNum, indRIRsort] = sort(repNum);
    FN.RIR = FN.RIR(indRIRsort);
    XStimParams.numreps = length(unique(repNum));
    set(H.space_RIR_numreps,'string',num2str(XStimParams.numreps))
    
    % display 1st picked file
    set(H.space_RIR_stimFN(1),'string',FN.RIR(1))
    
    % get locations and nptsTotalPlay
    eval(['load ' FN.stim_path2 char(FN.RIR(1)) ' nptsTotalPlay STIMdir*'])
    XStimParams.RIRpts2use = nptsTotalPlay;
    XStimParams.locations2 = STIMdir;           % all possible locations
    XStimParams.curr_stimdur = round(nptsTotalPlay/30);
    
    % update stimulus duration
    set(H.space_RIR_DUR, 'String',num2str(XStimParams.curr_stimdur));
    % update #files picked
    set(H.space_RIR_nFilestxt, 'String', [num2str(length(FN.RIR)) ' files picked ' ]);
    clear ind tempFN D nFiles temp
end

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim