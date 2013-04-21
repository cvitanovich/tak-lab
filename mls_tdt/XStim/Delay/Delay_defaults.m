function Delay_defaults(DefaultStr)
% function Delay_defaults(DefaultStr)
% DefaultStr default settings name
% Load default variables and update controls

global XStimParams;
global H;

XStimParams.Delay_default = get(H.delay_defaults,'Value');
DefaultStr = get(H.delay_defaults,'String');

menuStr = {'Onset_PE','OnRamps','LagXs','NoLagSeg',' ','NaNs','OffRamps'};

switch XStimParams.Delay_default
    case 1 % short prelim
        disp([char(menuStr(1)) ' (60 trials)']);
        
        XStimParams.SilentDelay =1;
        XStimParams.uncorrel =1;
        XStimParams.SilentLag =0;
        XStimParams.randOnsetPerms =1;
        XStimParams.DelayMods = [1 2 3]; % 1=normal, 2=uncorr, 3=corr no lead segment, 4=corr no lag segment
        XStimParams.DelayTimes =    [2 5 10 15 20 30 50 200 ...
                                    -2 -5 -10 -15 -20 -30 -50 -200];
        XStimParams.ramp_timeS_on = [2.5 NaN NaN NaN NaN NaN];
        XStimParams.ramp_timeS_off = [2.5 NaN NaN NaN NaN NaN];
        XStimParams.lag_seg_ext = [0 NaN NaN NaN NaN NaN NaN NaN];
        XStimParams.OnOff_mode = 1;  % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
        %XStimParams.Delay_default = 1;
    case 2 
        disp([char(menuStr(2))  ' (96 trials)']);
        
        XStimParams.SilentDelay =0;
        XStimParams.uncorrel =0;
        XStimParams.SilentLag =0;
        XStimParams.randOnsetPerms =1;
        XStimParams.DelayMods = [1];  % 1=normal, 2=uncorr, 3=corr no lead segment, 4=corr no lag segment
        XStimParams.DelayTimes =    [2 5 10 15 20 30 50 200 ...
                                    -2 -5 -10 -15 -20 -30 -50 -200];
        XStimParams.ramp_timeS_on = [0 2.5 10 20 30 50];
        XStimParams.ramp_timeS_off = [2.5 NaN NaN NaN NaN NaN];
        XStimParams.lag_seg_ext = [0 NaN NaN NaN NaN NaN NaN NaN];
        XStimParams.OnOff_mode = 1; % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
        %XStimParams.Delay_default = 1;
    case 3
        disp([char(menuStr(3))  ' (26 trials)']);
        
        XStimParams.SilentDelay =0;
        XStimParams.uncorrel =0;
        XStimParams.SilentLag =0;
        XStimParams.randOnsetPerms =1;
        XStimParams.DelayMods = [1]; % 1=normal, 2=uncorr, 3=corr no lead segment, 4=corr no lag segment
        XStimParams.DelayTimes =    [2 5 10 15 20 30 50 200 ...
                                    -2 -5 -10 -15 -20 -30 -50 -200];
        XStimParams.ramp_timeS_on = [2.5 NaN NaN NaN NaN NaN];
        XStimParams.ramp_timeS_off = [2.5 NaN NaN NaN NaN NaN];
        XStimParams.lag_seg_ext = [-10 -15 -20 -25 10 15 20 25];
        XStimParams.OnOff_mode = 1; % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
        %XStimParams.Delay_default = 1;
        
    case 4
        disp([char(menuStr(4))  ' (56 trials)']);
        
        XStimParams.SilentDelay =0;
        XStimParams.uncorrel =0;
        XStimParams.SilentLag =1;
        XStimParams.randOnsetPerms =1;
        XStimParams.DelayMods = [1 4]; % 1=normal, 2=uncorr, 3=corr no lead segment, 4=corr no lag segment
        XStimParams.DelayTimes =    [2 5 10 15 20 30 50 200 ...
                                    -2 -5 -10 -15 -20 -30 -50 -200];
        XStimParams.ramp_timeS_on = [2.5 50 NaN NaN NaN NaN];
        XStimParams.ramp_timeS_off = [2.5 NaN NaN NaN NaN NaN];
        XStimParams.lag_seg_ext = [0 NaN NaN NaN NaN NaN NaN NaN];
        XStimParams.OnOff_mode = 1; % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
        %XStimParams.Delay_default = 1;
        
        
    case 5
        return;
        
    case 6
        disp([char(menuStr(6))  ' (1 trial)']);
        
        XStimParams.SilentDelay =0;
        XStimParams.uncorrel =0;
        XStimParams.SilentLag =0;
        XStimParams.randOnsetPerms =1;
        XStimParams.DelayMods = [1]; % 1=normal, 2=uncorr, 3=corr no lead segment, 4=corr no lag segment
        XStimParams.DelayTimes =    [2 NaN NaN NaN NaN NaN NaN NaN ...
                                    NaN NaN NaN NaN NaN NaN NaN NaN];
        XStimParams.ramp_timeS_on = [2.5 NaN NaN NaN NaN NaN];
        XStimParams.ramp_timeS_off = [2.5 NaN NaN NaN NaN NaN];
        XStimParams.lag_seg_ext = [0 NaN NaN NaN NaN NaN NaN NaN];
        XStimParams.OnOff_mode = 1; % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
        %XStimParams.Delay_default = 1;
        
    case 7 
        disp([char(menuStr(7))  ' (96 trials)']);
        
        XStimParams.SilentDelay =0;
        XStimParams.uncorrel =0;
        XStimParams.SilentLag =0;
        XStimParams.randOnsetPerms =1;
        XStimParams.DelayMods = [1];  % 1=normal, 2=uncorr, 3=corr no lead segment, 4=corr no lag segment
        XStimParams.DelayTimes =    [2 5 10 15 20 30 50 200 ...
                                    -2 -5 -10 -15 -20 -30 -50 -200];
        XStimParams.ramp_timeS_on = [2.5 NaN NaN NaN NaN NaN];
        XStimParams.ramp_timeS_off = [0 2.5 10 20 30 50];
        XStimParams.lag_seg_ext = [0 NaN NaN NaN NaN NaN NaN NaN];
        XStimParams.OnOff_mode = 3; % 1=scr1 is rand, 2=both src's rand, 3=scr2 is rand
        %XStimParams.Delay_default = 1;
        
        
    otherwise
        warning('no defaults were set');
end

% update settings
set(H.uncorrel, 'Value', XStimParams.uncorrel);
set(H.silentdelay, 'Value', XStimParams.SilentDelay);
set(H.silentlag, 'Value', XStimParams.SilentLag);

set(H.randOnsetPerms, 'Value', XStimParams.randOnsetPerms);

set(H.delay01, 'String', num2str(XStimParams.DelayTimes(1)));
set(H.delay02, 'String', num2str(XStimParams.DelayTimes(2)));
set(H.delay03, 'String', num2str(XStimParams.DelayTimes(3)));
set(H.delay04, 'String', num2str(XStimParams.DelayTimes(4)));
set(H.delay05, 'String', num2str(XStimParams.DelayTimes(5)));
set(H.delay06, 'String', num2str(XStimParams.DelayTimes(6)));
set(H.delay07, 'String', num2str(XStimParams.DelayTimes(7)));
set(H.delay08, 'String', num2str(XStimParams.DelayTimes(8)));
set(H.delay09, 'String', num2str(XStimParams.DelayTimes(9)));
set(H.delay10, 'String', num2str(XStimParams.DelayTimes(10)));
set(H.delay11, 'String', num2str(XStimParams.DelayTimes(11)));
set(H.delay12, 'String', num2str(XStimParams.DelayTimes(12)));
set(H.delay13, 'String', num2str(XStimParams.DelayTimes(13)));
set(H.delay14, 'String', num2str(XStimParams.DelayTimes(14)));
set(H.delay15, 'String', num2str(XStimParams.DelayTimes(15)));
set(H.delay16, 'String', num2str(XStimParams.DelayTimes(16)));

set(H.onramp01, 'String', num2str(XStimParams.ramp_timeS_on(1)));
set(H.onramp02, 'String', num2str(XStimParams.ramp_timeS_on(2)));
set(H.onramp03, 'String', num2str(XStimParams.ramp_timeS_on(3)));
set(H.onramp04, 'String', num2str(XStimParams.ramp_timeS_on(4)));
set(H.onramp05, 'String', num2str(XStimParams.ramp_timeS_on(5)));
set(H.onramp06, 'String', num2str(XStimParams.ramp_timeS_on(6)));

set(H.DelayOnOff_mode, 'Value', XStimParams.OnOff_mode);

set(H.offramp01, 'String', num2str(XStimParams.ramp_timeS_off(1)));
set(H.offramp02, 'String', num2str(XStimParams.ramp_timeS_off(2)));
set(H.offramp03, 'String', num2str(XStimParams.ramp_timeS_off(3)));
set(H.offramp04, 'String', num2str(XStimParams.ramp_timeS_off(4)));
set(H.offramp05, 'String', num2str(XStimParams.ramp_timeS_off(5)));
set(H.offramp06, 'String', num2str(XStimParams.ramp_timeS_off(6)));

set(H.lagext01, 'String', num2str(XStimParams.lag_seg_ext(1)));
set(H.lagext02, 'String', num2str(XStimParams.lag_seg_ext(2)));
set(H.lagext03, 'String', num2str(XStimParams.lag_seg_ext(3)));
set(H.lagext04, 'String', num2str(XStimParams.lag_seg_ext(4)));
set(H.lagext05, 'String', num2str(XStimParams.lag_seg_ext(5)));
set(H.lagext06, 'String', num2str(XStimParams.lag_seg_ext(6)));
set(H.lagext07, 'String', num2str(XStimParams.lag_seg_ext(7)));
set(H.lagext08, 'String', num2str(XStimParams.lag_seg_ext(8)));

set(H.delay_defaults, 'Value', XStimParams.Delay_default);

%SetInfo_Delay;