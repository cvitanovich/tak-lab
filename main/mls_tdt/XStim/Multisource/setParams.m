global GUI
global H
global XStimParams
global FN


sp_pos_incr = 20;
sp_tot_ht = 0;
sp_ind_ht = 20;
sp_Two_source_wd = 1200; 
sp_Two_source_ht = size(XStimParams.sub_locations, 1)*30;

% define figure
H.setparamsfig = figure('Units', 'pixels',...
	'Position', [10 420 1215 sp_Two_source_ht],...
	'Name', 'Set/View Stimulus Params', ...
	'NumberTitle', 'off');
 
sp_tot_ht = sp_tot_ht + sp_ind_ht;

if size(XStimParams.sub_locations, 1) > 0
	for count1 = 1:size(XStimParams.sub_locations, 1)
		sp_tot_ht = sp_tot_ht + sp_ind_ht;
		uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[10 sp_Two_source_ht-sp_tot_ht 50 sp_ind_ht],...
			'FontWeight','bold',...	
			'Visible','on',...   
			'String', ['stim ' num2str(count1)]);
		H.Twosrc_stim_type(count1) = uicontrol('Parent', H.setparamsfig,...
			'Style','popup',...
			'Units','pixels',...
			'Position',[60 sp_Two_source_ht-sp_tot_ht 120 sp_ind_ht],...
			'BackgroundColor','white',...
			'String',...
			'Tone|Gammatone|Tone Stack|GammaTone Stack|BandPassed Noise|NarrowBand Noise|NarrowBand Noise Stack|BroadBand Noise|File|None',...
			'Value', H.stim_type(count1),...
			'Callback','setInfo_setParams');
		H.Twosrc_stim_filenamepb(count1) = uicontrol('Parent', H.setparamsfig,...
			'Style','pushbutton',...
			'Units','pixels',...
			'Position',[180 sp_Two_source_ht-sp_tot_ht 25 sp_ind_ht],...
			'Enable','off',...    
			'BackgroundColor','red',...
			'ForegroundColor','white',...
			'String', 'file',...
			'Callback','setInfo_setParams');
		H.Twosrc_stim_filename(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[205 sp_Two_source_ht-sp_tot_ht 75 sp_ind_ht],...
			'Visible','on',...  
			'ForegroundColor','blue',...
			'String', '');
% 		H.nostim_txt(count1) = uicontrol('Parent',H.setparamsfig,...
% 			'Style','text',...
% 			'Units','pixels',...
% 			'Position',[180 sp_Two_source_ht-sp_tot_ht 90 sp_ind_ht],...
% 			'Visible','off',...   
% 			'String', 'no stim');
		H.freq_txt(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[280 sp_Two_source_ht-sp_tot_ht 30 sp_ind_ht],...
			'Visible','on',...   
			'String', 'freq');
		H.Twosrc_freq(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','edit',...
			'Units','pixels',...
			'Position',[310 sp_Two_source_ht-sp_tot_ht 40 sp_ind_ht],...
			'Visible','on',...    
			'String', num2str(XStimParams.freq(count1)),...
			'Callback','setInfo_setParams');
		H.factor_txt(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[350 sp_Two_source_ht-sp_tot_ht 90 sp_ind_ht],...
			'String', 'weighting factor');
		H.factor(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','edit',...
			'Units','pixels',...
			'Position',[440 sp_Two_source_ht-sp_tot_ht 40 sp_ind_ht],...
			'String', num2str(XStimParams.factor(count1)),...
			'Callback','setInfo_setParams');
		H.bandwidth_txt(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[480 sp_Two_source_ht-sp_tot_ht 90 sp_ind_ht],...
			'String', 'bandwidth (oct)');
		H.bandwidth(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','edit',...
			'Units','pixels',...
			'Position',[570 sp_Two_source_ht-sp_tot_ht 40 sp_ind_ht],...
			'String', num2str(XStimParams.bandwidth(count1)),...
			'Callback','setInfo_setParams');
		H.limit_txt(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[610 sp_Two_source_ht-sp_tot_ht 70 sp_ind_ht],...
			'String', 'limits low/high');
		H.bandlimits1(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','edit',...
			'Units','pixels',...
			'Position',[680 sp_Two_source_ht-sp_tot_ht 40 sp_ind_ht],...
			'String', num2str(XStimParams.bandlimits1(count1)),...
			'Callback','setInfo_setParams');
		H.bandlimits2(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','edit',...
			'Units','pixels',...
			'Position',[720 sp_Two_source_ht-sp_tot_ht 40 sp_ind_ht],...
			'String', num2str(XStimParams.bandlimits2(count1)),...
			'Callback','setInfo_setParams');
%%%%%%%%%%%%%%%%%%%
% 		switch XStimParams.mod_type(count1)
% 			case 1
% 				H.mod_num(count1) = 1;
% 			case 2
% 				H.mod_num(count1) = 2;
% 			case 3
% 				H.mod_num(count1) = 3;
% 			case 4
% 				H.mod_num(count1) = 4;
% 			otherwise
% 				H.mod_num(count1) = 4;
% 		end
		uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[760 sp_Two_source_ht-sp_tot_ht 60 sp_ind_ht],...
			'FontWeight','bold',...	
			'Visible','on',...   
			'String', ['mod ' num2str(count1)]);
		H.Twosrc_mod_type(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','popup',...
			'Units','pixels',...
			'Position',[820 sp_Two_source_ht-sp_tot_ht 50 sp_ind_ht],...
			'BackgroundColor','White',...
			'String',...
			'Tone|LP Noise|File|None',...
			'Value',  H.mod_type(count1),...
			'Callback','setInfo_setParams');
		H.Twosrc_mod_pb(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','pushbutton',...
			'Units','pixels',...
			'Position',[870 sp_Two_source_ht-sp_tot_ht 25 sp_ind_ht],...
			'BackgroundColor','red',...
			'ForegroundColor','white',...
			'String', 'file',...
			'Enable','off',...    	
			'Visible','on',...
			'Callback','setInfo_setParams');
		H.mod_file_txt(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'ForegroundColor','blue',...
			'Position',[895 sp_Two_source_ht-sp_tot_ht 75 sp_ind_ht],...
			'Visible','on',...
			'String', '');
% 		H.nomod_txt(count1) = uicontrol('Parent',H.setparamsfig,...
% 			'Style','text',...
% 			'Units','pixels',...
% 			'ForegroundColor','blue',...
% 			'Position',[1095 sp_Two_source_ht-sp_tot_ht 100 sp_ind_ht],...
% 			'Visible','off',...
% 			'String', '');
		H.mod_freq_txt(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[970 sp_Two_source_ht-sp_tot_ht 30 sp_ind_ht],...
			'String', 'freq');
		H.Twosrc_mod_freq(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','edit',...
			'Units','pixels',...
			'Enable','off',...
			'Position',[1000 sp_Two_source_ht-sp_tot_ht 50 sp_ind_ht],...
			'String', num2str(XStimParams.mod_freq(count1)),...
			'Callback','setInfo_setParams');
		H.mod_depth_txt = uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[1050 sp_Two_source_ht-sp_tot_ht 30 sp_ind_ht],...
			'String', 'depth');
		H.Twosrc_mod_depth(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','edit',...
			'Units','pixels',...
			'Enable','off',...
			'Position',[1080 sp_Two_source_ht-sp_tot_ht 40 sp_ind_ht],...
			'String', num2str(XStimParams.mod_depth(count1)),...
			'Callback','setInfo_setParams');
		H.mod_phase_txt = uicontrol('Parent',H.setparamsfig,...
			'Style','text',...
			'Units','pixels',...
			'Position',[1120 sp_Two_source_ht-sp_tot_ht 40 sp_ind_ht],...
			'String', 'phase');
		H.Twosrc_mod_phase(count1) = uicontrol('Parent',H.setparamsfig,...
			'Style','edit',...
			'Units','pixels',...
			'Enable','off',...
			'Position',[1160 sp_Two_source_ht-sp_tot_ht 30 sp_ind_ht],...
			'String', num2str(XStimParams.mod_phase(count1)),...
			'Callback','setInfo_setParams');
	end
end

sp_tot_ht = sp_tot_ht + sp_ind_ht;
sp_tot_ht = sp_tot_ht + sp_ind_ht;

uicontrol('Style','text',...
    'Units','pixels',...
    'Position',[10 sp_Two_source_ht-sp_tot_ht 80 sp_ind_ht],...
    'String', 'stim dur (ms)');
H.DUR = uicontrol('Style','edit',...
    'Units','pixels',...
    'Position',[90 sp_Two_source_ht-sp_tot_ht 40 sp_ind_ht],...
    'String', num2str(XStimParams.curr_stimdur),...
    'Callback','setInfo_setParams');

setInfo_setParams;

 