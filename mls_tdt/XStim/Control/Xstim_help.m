function [] = Xstim_help

% Xstim_help
% function to provide help from main Xstim menu
% 

global H

val = get(H.Xstim_help,'Value');

switch val
   
case 1      % Timing Hangups
    DLGNAME = 'Timing Hangups';
    HELPSTRING{1} = 'When playing out stimuli, the rhythm sometimes skips.';
    HELPSTRING{2} = 'This appears to be due to too much traffic on the internet.';
    HELPSTRING{3} = '';
    HELPSTRING{4} = 'RT CLICK on MY COMPUTER. GO TO PROPERTIES, DEVICE MANAGER. ';
    HELPSTRING{5} = 'GO TO NETWORK ADAPTERS - FAST ETHERNET ADAPTER - PROPERTIES';
    HELPSTRING{6} = 'CHECK DISABLE IN THIS HARDWARE PROFILE';
    HELPSTRING{7} = 'NOW, UNPLUG THE INTERNET CABLE and REBOOT';
    HELPSTRING{8} = '';
    HELPSTRING{9} = 'IF this does not work, remove the network card from the computer';
    HELPSTRING{10} = '';
    HELPSTRING{11} = 'TO RETURN TO NORMAL FUNCTION, REVERSE ALL THIS.';
case 2      % HRTF files
    DLGNAME = 'HRTF files';
    HELPSTRING{1} = 'TRYING TO FIGURE OUT WHICH HRTF FILES TO USE?';
    HELPSTRING{2} = '';
    HELPSTRING{3} = 'First decide whether to use separate compensation for the earphones';
    HELPSTRING{4} = ' and ear canal, with ET_ECcs.inv,';
    HELPSTRING{5} = ' or whether you want these built in to the HRTFs.';
    HELPSTRING{6} = '';
    HELPSTRING{7} = 'For SEPARATE COMPENSATION, ';
    HELPSTRING{8} = 'load the ET_ECcs.inv file for your bird to the DSPs.';
    HELPSTRING{9} = 'Then use: ';
    HELPSTRING{10} = 'out8a_ (or *.std) for fully-cued tests';
    HELPSTRING{11} = 'out16a (or *.ILA) for ILD-alone tests';
    HELPSTRING{12} = 'out18a (or *.ITA) for ITD-alone tests';
    HELPSTRING{13} = '';
    HELPSTRING{14} = 'For BUILT-IN COMPENSATION, use:';
    HELPSTRING{15} = 'out12 (or *.eq) for fully-cued tests';
    HELPSTRING{16} = 'out16d (or *.ILA) for ILD-alone tests';
    HELPSTRING{17} = 'out18 (or *.ITA) for ITD-alone tests';
otherwise
    DLGNAME = 'SORRY';
    HELPSTRING{1} = 'MICROSOFT IS STILL DEVELOPING ANOTHER WAY TO RUIN YOUR WHOLE DAY.';
   set(H.datatools,'Value',1);
end

HELPDLG(HELPSTRING,DLGNAME);