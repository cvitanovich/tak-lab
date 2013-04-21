% update dataFile name

if XStimParams.testnum < 10
    FN.data = [num2str(XStimParams.bird_number)...
            XStimParams.recording_site...
            '00'...
            num2str(XStimParams.testnum)];
elseif XStimParams.testnum < 100
    FN.data = [num2str(XStimParams.bird_number)...
            XStimParams.recording_site...
            '0'...
            num2str(XStimParams.testnum)];
elseif XStimParams.testnum < 1000
    FN.data = [num2str(XStimParams.bird_number)...
            XStimParams.recording_site...
            num2str(XStimParams.testnum)];
else
    FN.data = [];
    disp('error making data file name');
    return
end
set(H.xstim_recorddata_FN,'String',FN.data);

% check to see if file exists

if exist([FN.data_path FN.data '.mat'],'file') == 2
    set(H.fig,'Color', [.1 .2 .8]);
    drawnow
    R = input([FN.data ' already exists. \n Would you like to over write it (0 or 1)? \n']);
    set(H.fig,'Color', [.8 .8 .8]);
    if R
        eval(['delete ' FN.data_path FN.data '.mat']);
    else
        
        XStimParams.testnum = XStimParams.testnum +1;
        set(H.testnum, 'String', num2str(XStimParams.testnum))
        update_dataFN;
    end
end
