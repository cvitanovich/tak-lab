%globals_var
%
% sets up global variables
global FN
global H

% added HRTFfiletype[1-6] on Sept 02, 2004
%
% 1-6, in order:
%   space, ila, ita, loc, ephone, ephone2
%
% filetypes:    0) untested
%               1) traditional binary
%               2) new *.mat with TRF1 and TF2
%               999) non-readable
%
% SHOULD probably run Globals_FN_update after this
% to go with current bird_number

home_drive = 'e:\';

FN.current_path = [home_drive 'kip\matlab\mls_tdt\XStim\control\' ];
FN.itd_path = [home_drive 'kip\matlab\mls_tdt\XStim\ITDfilters\' ];
FN.HRTF_path = [home_drive 'kip\HRTFdata\' ];
FN.data_path = [home_drive 'kip\datastor\' ];
FN.script_path = [home_drive 'kip\matlab\mls_tdt\XStim\space-test\stimuli\' ];
FN.temp_stim_path = [home_drive 'kip\noisetc\temp\' ];
FN.temp_stim_path2 = 'h:\kip_temp\';
FN.stim_path = [home_drive 'kip\noisetc\' ];
FN.stim_path2 = [home_drive 'kip\noisetc\' ];
FN.stim_path3 = [home_drive 'kip\noisetc\' ];
FN.ephone_path = [home_drive 'kip\HRTFdata\' ];
FN.ILA_path = [home_drive 'kip\HRTFdata\' ];
FN.ITA_path = [home_drive 'kip\HRTFdata\' ];
FN.ITD_path = [home_drive 'kip\matlab\mls_tdt\XStim\ITDFilters\' ];
FN.space_path = [home_drive 'kip\HRTFdata\' ];
FN.loc_path = [home_drive 'kip\HRTFdata\' ];
FN.home = [home_drive 'kip\matlab\mls_tdt\Xstim\' ];
FN.stim = 'rnd1_30.noi';
FN.data = [];
FN.stim2 = 'l_75_c01.noi';
FN.stim3 = 'l_75_b01.noi';

if 0
FN.stim4L = 'l_75_b01.noi';
FN.stim5L = 'l_75_b01.noi';
FN.stim6L = 'l_75_b01.noi';
FN.stim7L = 'l_75_b01.noi';
FN.stim8L = 'l_75_b01.noi';
FN.stim9L = 'l_75_b01.noi';
FN.stim10L = 'l_75_b01.noi';
FN.stim4R = 'l_75_b01.noi';
FN.stim5R = 'l_75_b01.noi';
FN.stim6R = 'l_75_b01.noi';
FN.stim7R = 'l_75_b01.noi';
FN.stim8R = 'l_75_b01.noi';
FN.stim9R = 'l_75_b01.noi';
FN.stim10R = 'l_75_b01.noi';
end

FN.mod_path = [home_drive 'kip\noisetc\' ];
FN.mod_path2 = [home_drive 'kip\noisetc\' ];
FN.mod_path3 = [home_drive 'kip\noisetc\' ];
FN.mod = [];
FN.mod2 = [];
FN.mod3 = [];
FN.saveHRTF = [];

%FN.script = 'checkn90.scr';
FN.space_eq = '938AD.eq.mat';
FN.ildalone_eq = '938AD.ila.mat';
FN.itdalone_eq = '938AD.ita.mat';
FN.space_std = '938AD.std.mat';
FN.ildalone_std = '938AD.ila.std.mat';
FN.itdalone_std = '938AD.ita.std.mat';
FN.loc = '938AD.eq.mat';
FN.ephone2 = '938ET_EC_inv.mat';
FN.ablequal_eq = '';
FN.ablequal_std = '';

FN.HRTFfiletype = zeros(6,2);
FN.ephone = 'onesX.imp';
%FN.ephone2 = 'onesX.imp';

% list of stim filenames organized by each rep
FN.FRP = [];

% get rid of all old handles
clear global H

% reset readme file
readme = [];

% save to disk
%eval(['save ' FN.current_path 'FN_current FN;'])
