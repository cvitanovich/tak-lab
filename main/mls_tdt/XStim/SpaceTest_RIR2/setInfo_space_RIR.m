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

%%%%%%%%%%%%%%%%%%
% initialize RIR FNs
if get(H.space_RIR_RIRfilepb,'value')
    FN.RIR = [];
    set(H.space_RIR_RIRFN,'string','')
    
    [temp,FN.stim_path2] = uigetfile([FN.stim_path2 '*.*'],'Select RIR PATH');
    
    % get directory of FN.stim_path2
    D = dir(FN.stim_path2);
    nFiles = length(D)-2;
    
    for ifile = 1:nFiles
        tempFN{ifile} = D(ifile+2).name;
    end
    
    % choose which RIR FNs to include
    if exist('tempFN')
        ind = listdlg('PromptString','Select RIR Files',...
            'ListString', tempFN,...
            'InitialValue',[1:nFiles]);
        clear stimFN
        for i = 1:length(ind)
            FN.RIR{i} = tempFN{ind(i)};  % re-organized in engage for each rep
        end
    end
    nFiles = length(FN.RIR);
    % update display
    set(H.space_RIR_RIRFN,'string',strvcat(FN.RIR))
    
    for ifile=1:nFiles
        eval(['S = load(''' FN.stim_path2 FN.RIR{ifile} ''');'])
        [nLocs(ifile) RIRpts(ifile)] = size(S.TF1);
        if ISFIELD(S,'TF1')
            RIR1{ifile} = S.TF1;
        else
            error(['no TF1 field in ' FN.RIR{ifile}])
        end
        if ISFIELD(S,'TF2')
            RIR2{ifile} = S.TF2;
        else
            error(['no TF2 field in ' FN.RIR{ifile}])
        end
        if ISFIELD(S,'dir')
            dirMat{ifile} = S.dir;
        else
            error(['no dir field in ' FN.RIR{ifile}])
        end
    end
    set(H.space_RIR_RIRpts, 'String',strvcat(num2str(RIRpts')));
    
    prompt = strvcat('nPts for the current RIRs:','',num2str(RIRpts'),'','How many RIR pts to use?');
    A = inputdlg(prompt,'choose RIR #pts to use');
    XStimParams.RIRpts2use = str2num(A{1});
    set(H.space_RIR_RIRpts2use, 'String',A{1});
    
    ifile = min1(find(nLocs==max1(nLocs)));
    [R C] = size(dirMat{ifile});
    if C==2
        GUI.locations1 = dirMat{ifile};
    else
        GUI.locations1 = dirMat{ifile}';
    end
    refresh_spaceaxis
    % locations first loaded into GUI
    % then copied to XStimParams.locations2 for storage
    % then picking allowed
    % picked are then placed in XStimParams.locations
    XStimParams.locations2 = GUI.locations1';        
    
    clear ind tempFN D nFiles temp R C S A
end

if ~strcmp(get(H.space_RIR_RIRpts2use, 'String'),num2str(XStimParams.RIRpts2use))
    XStimParams.RIRpts2use = str2num(get(H.space_RIR_RIRpts2use, 'String'));
end

% set numreps
XStimParams.numreps = str2num(get(H.space_RIR_numreps,'String'));

% set ramppts
XStimParams.ramppts = str2num(get(H.space_RIR_ramppts,'String'))*30;

% set AM params
XStimParams.mod_freq(1) = str2num(get(H.space_RIR_mod_freq,'String'));
XStimParams.mod_depth(1) = str2num(get(H.space_RIR_mod_depth,'String'));
XStimParams.mod_phase(1) = str2num(get(H.space_RIR_mod_phase,'String'));

%Set space_RIR PickerFlag
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
if get(H.space_RIR_inc_testnum,'Value')
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum))
    update_dataFN;
    set(H.space_RIR_inc_testnum,'Value',0);
end
set(H.space_RIRdata_FN,'String', FN.data);

eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim