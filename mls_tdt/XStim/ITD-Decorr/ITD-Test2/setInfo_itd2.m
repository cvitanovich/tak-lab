%setInfo_itd

% increment test number
if(exist1('H.itdfig') & get(H.inc_testnum,'Value'))
	XStimParams.testnum = XStimParams.testnum +1;
	set(H.testnum, 'String', num2str(XStimParams.testnum))
   update_dataFN;
   set(H.inc_testnum,'Value',0);
end


set(H.stim_type,'Value',9);
set(H.stim_type2,'Value',9);
% FNstim pushbuttons
if get(H.itd_FNpb(1),'value')
    set(H.itd_FNpb(1),'value',0)
    FN.stim = 'nosuchfile';
    while exist1([FN.stim_path FN.stim]) ~=2
        [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.*'],'Select stim#1 File');
        if(FN.stim_path ~= 0)
            set(H.stim_filename,'String',FN.stim);
            set(H.itd_FNtxt(1),'String',FN.stim);
        end
    end
end
    
    
if get(H.itd_FNpb(2),'value')
    set(H.itd_FNpb(2),'value',0)
    FN.stim2 = 'nosuchfile';
    while exist1([FN.stim_path2 FN.stim2]) ~=2
        [FN.stim2,FN.stim_path2] = uigetfile([FN.stim_path2 '*.*'],'Select stim#1 File');
        if(FN.stim_path2 ~= 0)
            set(H.stim_filename2,'String',FN.stim2);
            set(H.itd_FNtxt(2),'String',FN.stim2);
        end
    end
end


% set up noises from file
set(H.stim_filename,'String',FN.stim);
set(H.stim_filename2,'String',FN.stim2);
D = dir([FN.stim_path FN.stim]);
XStimParams.curr_stimdur = round(1000 * (D.bytes/4) / TDT.Fs);
set(H.DUR, 'String', num2str(XStimParams.curr_stimdur));
    
%setITDstep
low_itd  = str2num(get(H.lowitd,'String'));
high_itd = str2num(get(H.highitd,'String'));
if(high_itd > 250)
   high_itd = 250;
   set(H.highitd,'String',num2str(high_itd));
end
if(low_itd < -250)
   low_itd = -250;
   set(H.lowitd,'String',num2str(low_itd));
end
XStimParams.hiitd = high_itd;
XStimParams.loitd = low_itd;

numitds = str2num(get(H.numitds,'String'));

%set(H.stepitd,'String',num2str(round((high_itd-low_itd)/(numitds-1))));

% check out modulation parameters for stim
XStimParams.mod_depth(1) = str2num(get(H.ITDtest2mod_depth,'String'));
mod_type = get(H.ITDtest2mod_type,'String');
mod_num = get(H.ITDtest2mod_type,'Value');
XStimParams.mod_type = deblank(mod_type(mod_num,:));

switch XStimParams.mod_type
    case 'Tone'
        set(H.ITDtest2mod_txt,'String', 'Freq:');
        set(H.ITDtest2mod_freq,'Visible','on');
        set(H.ITDtest2mod_phase,'Visible','on');
        set(H.ITDtest2mod_txtA,'Visible','on');
        set(H.ITDtest2mod_depth,'Visible','on');
        set(H.ITDtest2mod_txtB,'Visible','on');
        set(H.ITDtest2mod_pb,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.ITDtest2mod_freq,'String'));
        XStimParams.mod_phase(1) = str2num(get(H.ITDtest2mod_phase,'String'));
    case 'File'
        set(H.ITDtest2mod_pb,'Visible','on');
        set(H.ITDtest2mod_freq,'Visible','off');
        set(H.ITDtest2mod_txtA,'Visible','on');
        set(H.ITDtest2mod_depth,'Visible','on');
        set(H.ITDtest2mod_phase,'Visible','off');
        set(H.ITDtest2mod_txtB,'Visible','off');
        if get(H.ITDtest2mod_pb,'Value')
            [mod_from_file, FN.mod, FN.mod_path] = Get_Stimulus_File(XStimParams, FN.mod_path);
            set(H.ITDtest2mod_pb,'Value',0);
        end
        set(H.ITDtest2mod_txt,'String', FN.mod);
        
    case 'LP Noise'
        set(H.ITDtest2mod_txt,'String', 'CutOff Freq:');
        set(H.ITDtest2mod_freq,'Visible','on');
        set(H.ITDtest2mod_txtA,'Visible','on');
        set(H.ITDtest2mod_depth,'Visible','on');
        set(H.ITDtest2mod_phase,'Visible','off');
        set(H.ITDtest2mod_pb,'Visible','off');
        set(H.ITDtest2mod_txtB,'Visible','off');
        XStimParams.mod_freq(1) = str2num(get(H.ITDtest2mod_freq,'String'));
        
    case 'None'
        set(H.ITDtest2mod_txt,'String', 'no mod  ');
        set(H.ITDtest2mod_freq,'Visible','off');
        set(H.ITDtest2mod_phase,'Visible','off');
        set(H.ITDtest2mod_pb,'Visible','off');
        set(H.ITDtest2mod_txtB,'Visible','off');
        set(H.ITDtest2mod_txtA,'Visible','off');
        set(H.ITDtest2mod_depth,'Visible','off');
        
    otherwise
end



% save current parameter values
eval(['save ' FN.current_path 'XStimParams_current XStimParams;'])
eval(['save ', FN.current_path  'FN_current FN;'])
update_XStim