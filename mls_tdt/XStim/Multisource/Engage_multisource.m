function [] = Engage_multisource()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_multisource

%*******************************************************************************
%	The Two_source Test operation
%*******************************************************************************
% NOTE: (Mar 4, 2003) problem with writing stim files discovered:
%   locind has multiple rows and this was not taken into account when writing files
%   so actual file written (and played out) was ceil(locind(randseq(trialnum))/2)
%   
% NOTE: (Jun 1, 2003) problem with make_env that always used mod_freq and mod_depth
%   from stimulus number1 (so all stims had same mod depth and freq
%   changed to use input param1,param2 to carry the correct depth and freq
%   also: eliminated 2nd buffercycle and associated buffers...
%
% NOTE: (July 29, 2004) added in ABL-alone and removed filecache, itd2,3 and az2,3 and
% el2,3 assuming these will always be the same as #1
% removed location script ability (forcing spacepicker)

rand('state',sum(100*clock));

stimuli_dir = FN.temp_stim_path;
fclose all;
eval(['delete ' stimuli_dir '*.*;']);

% check if filt files assigned
if exist([FN.ephone_path FN.ephone2]) ~= 2
    ephonefilediagbox;
end
if XStimParams.space_flag |  XStimParams.ABLalone_flag
    while exist1([FN.space_path FN.space]) ~=2 | isempty(strfind(FN.space,'std'))
        [FN.space,FN.space_path] = uigetfile([FN.space_path '*.*'],'Select *.STD HRTF File');
        if(FN.space_path ~= 0)
            set(H.spacefile,'String',[FN.space_path FN.space]);
        end
        FN.HRTFfiletype(1) = testHRTFfiletype(FN.space_path, FN.space);
    end
    if XStimParams.space_flag
        disp('This is a Multi-source SPACE test')
    else
        disp('This is an ABL-alone SPACE test')
    end
elseif XStimParams.ildalone_flag
    while exist1([FN.ILA_path FN.ildalone]) ~=2
        [FN.ildalone,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.ildalonefile,'String',[FN.ILA_path FN.ildalone]);
        end      
        FN.HRTFfiletype(2) = testHRTFfiletype(FN.ILA_path, FN.ildalone);
    end
    disp('This is an ILD-ALONE test')
elseif XStimParams.itdalone_flag
    while exist1([FN.ILA_path FN.ildalone]) ~=2
        [FN.ildalone,FN.ILA_path] = uigetfile([FN.ILA_path '*.*'],'Select ILD-alone HRTF File');
        if(FN.ILA_path ~= 0)
            set(H.ildalonefile,'String',[FN.ILA_path FN.ildalone]);
        end
        FN.HRTFfiletype(2) = testHRTFfiletype(FN.ILA_path, FN.ildalone);
    end
    while exist1([FN.ITA_path FN.itdalone]) ~=2
        [FN.itdalone,FN.ITA_path] = uigetfile([FN.ITA_path '*.*'],'Select ITD-alone HRTF File');
        if(FN.ITA_path ~= 0)
            set(H.itdalonefile,'String',[FN.ITA_path FN.itdalone]);
        end
        FN.HRTFfiletype(3) = testHRTFfiletype(FN.ITA_path, FN.itdalone);
    end
    disp('This is an ITD-ALONE test')
end

%Put parameters into XStimParams
XStimParams.curr_ITD = str2num(get(H.ITD,'String'));
XStimParams.curr_ABL = str2num(get(H.ABL,'String'));
%XStimParams.curr_stimdur = str2num(get(H.DUR,'String'));
XStimParams.test_ISI = str2num(get(H.ISI,'String'));
XStimParams.numreps = str2num(get(H.numreps,'String'));
XStimParams.reset_flag = 0;

if XStimParams.ABLalone_flag
    % load filter coeffs & other params for ABLalone test
    Fs = 30000;
    cF = round(1000*exp(([40:168]/48)*log(2)))'; 
    n_cF = length(cF);
    fcoefs = Make_ERBFiltA(Fs,cF);
    maxFactor = .00003764*cF(n_cF)+.6236;
    Factor = maxFactor ./ (.00003764*cF+.6236);
    Factormat = repmat(Factor,1,255);
    fftpts = 2048;
    freq = 0:15000/((fftpts/2)-1):15000;
    freq_ind = nearest_index(freq,cF);
    endpt = max1(freq_ind);
    startpt = min1(freq_ind);
    Xpart1 = startpt:endpt; 
    Xpart2 = (fftpts + 2 - endpt):(fftpts + 2 - startpt);
end

% if (get(H.stim_type,'Value') == 9) & isempty(FN.stim)
%     [FN.stim,FN.stim_path] = uigetfile([FN.stim_path '*.noi'],'Select stimulus File');
% end

%Specify DAMA buffers
clear BUF
BUF.L1				= 1;
BUF.R1				= 2;
BUF.L2				= 3;
BUF.R2				= 4;
BUF.playseq_L1		= 5;
BUF.playseq_R1		= 6;
BUF.playspec1		= 7;

%Make play sequence buffers
S232('allot16',BUF.playseq_L1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.L1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_L1);

S232('allot16',BUF.playseq_R1,10);
S232('dpush', 10);
S232('value',0);
S232('make',0,BUF.R1);
S232('make',1,1);
S232('make',2,0);
S232('qpop16',BUF.playseq_R1);

%Make play specification buffer
S232('allot16',BUF.playspec1,10);
S232('dpush',10);
S232('value',0);
S232('make',0,BUF.playseq_L1);
S232('make',1,BUF.playseq_R1);
S232('make',2,0);
S232('qpop16',BUF.playspec1);

%Add a piece of silence prior to stimulus to calculate spontaneous rate, 3/22/01
silence_len = (XStimParams.silence_lead * round(TDT.Fs/1000));
%Add a piece of silence after stimulus 
silence_len2 = (XStimParams.silence_trail * round(TDT.Fs/1000));

%Make Stimulus buffers
DUR = XStimParams.curr_stimdur;
if XStimParams.ildalone_flag 					%ILDAlone Test
    S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.itdfiltlen + TDT.ephonefiltlen + TDT.hrtffiltlen);
elseif XStimParams.itdalone_flag 				%ITDAlone Test
    S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2);
    S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2);
    S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2);
    S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen*2);
elseif XStimParams.space_flag | XStimParams.ABLalone_flag	    %fully-cued Test
    S232('allot16',BUF.L1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.R1,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.L2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
    S232('allot16',BUF.R2,(silence_len2 + silence_len + DUR*round(TDT.Fs/1000)) + TDT.ephonefiltlen + TDT.hrtffiltlen);
end   

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
if XStimParams.ildalone_flag 					%ILDAlone Test
    S232('PD1npts',1,(silence_len2 + silence_len + DUR*(round(TDT.Fs/1000))) + TDT.itdfiltlen  + TDT.ephonefiltlen + TDT.hrtffiltlen);
elseif XStimParams.itdalone_flag 				%ITDAlone Test
    S232('PD1npts',1,(silence_len2 + silence_len + DUR*(round(TDT.Fs/1000)))  + TDT.ephonefiltlen + (2 * TDT.hrtffiltlen));
elseif XStimParams.space_flag |  XStimParams.ABLalone_flag			% fully-cued or ABL-alone Test
    S232('PD1npts',1,(silence_len2 + silence_len + DUR*(round(TDT.Fs/1000))) + TDT.ephonefiltlen + TDT.hrtffiltlen);
end

%Get Earphone filters
if FN.HRTFfiletype(6) == 1
    ephonefname = [FN.ephone_path FN.ephone2];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
elseif FN.HRTFfiletype(6) == 2
    eval(['load ' FN.ephone_path FN.ephone2]);
    ephonefilt_left  = TF1;
    ephonefilt_right = TF2;
else
    disp('ephone filters of incorrect type');
end

%Load Earphone filters
dspid_left = 0; dspid_right = 1;
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));           % was '0xFFF'
S232('dropall');
%Make connections for left ear
S232('PD1addsimp',1,S232('DSPout',dspid_left),S232('DAC',0)); %DSPout to DAC0
S232('PD1specIB',1,S232('IB',0),S232('DSPin',dspid_left)); %IB to DSPin
%Make connections for right ear
S232('PD1addsimp',1,S232('DSPout',dspid_right),S232('DAC',1));
S232('PD1specIB',1,S232('IB',1),S232('DSPin',dspid_right));
%Load left      
S232('pushf',ephonefilt_left,length(ephonefilt_left));
S232('PreLoadRaw',1,S232('DSPid',dspid_left),'MONO','STACK','','',TDT.ephonescale,1.0,1);
%Load right
S232('pushf',ephonefilt_right,length(ephonefilt_right));
S232('PreLoadRaw',1,S232('DSPid',dspid_right),'MONO','STACK','','',TDT.ephonescale,1.0,1);

set(H.locfile,'Enable','off');
set(H.locAZ,'Enable','off');
set(H.locEL,'Enable','off');
set(H.locuseit,'Enable','off');

%Set MII parameters
mii_us_per_sample = 10; 							%microsecond per sample
mii_separation = 100; 								%only take events separated by 100 samples (i.e., 1 ms)

ITD = XStimParams.curr_ITD;
if(abs(ITD) > 250) return; end

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

ISI = XStimParams.test_ISI;

%Get all HRTF spectra indices

% 6/28/05 START

% XStimParams.locations = GUI.locations1';
% if ~size(XStimParams.locations,2)
%     set(H.pickerfig,'Color', [.1 .2 .8]);
%     set(H.picker_error,'visible','on');
%     pause;
%     XStimParams.locations = GUI.locations1';
%     set(H.picker_error,'visible','off');
%     set(H.pickerfig,'Color', [.8 .8 .8]);
% end


clear locations;
eval(['load ' FN.current_path 'Locations_current1;'])
XStimParams.locations1 = locations';
clear locations;
eval(['load ' FN.current_path 'Locations_current2;'])
XStimParams.locations2 = locations';		
											
% 6/28/05 END

%%%%%%%%%% load all HRTFs with HRTFfiletype == 2
dir = 0;
if FN.HRTFfiletype(1) == 2
    eval(['load -mat ' FN.space_path FN.space]);
    TF1_space = TF1; TF2_space = TF2;
    dir_space = dir;
end
if FN.HRTFfiletype(2) == 2
    eval(['load -mat ' FN.ILA_path FN.ildalone]);
    TF1_ila = TF1; TF2_ila = TF2;
    dir_ILA = dir;
end
if FN.HRTFfiletype(3) == 2
    eval(['load -mat ' FN.ITA_path FN.itdalone]);
    TF1_ita = TF1; TF2_ita = TF2;
    dir_ITA = dir;
end
clear dir TF1 TF2

if XStimParams.space_flag | XStimParams.ABLalone_flag
    if FN.HRTFfiletype(1) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.space_path FN.space]));
    elseif FN.HRTFfiletype(1) == 2
        hrtfdirmat = dir_space;
    end
elseif XStimParams.ildalone_flag
    if FN.HRTFfiletype(2) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone]));
    elseif FN.HRTFfiletype(2) == 2
        hrtfdirmat = dir_ILA;
    end
elseif XStimParams.itdalone_flag
    if FN.HRTFfiletype(3) == 1
        hrtfdirmat = sph2dbl(mtlrdir([FN.ITA_path FN.itdalone]));
    elseif FN.HRTFfiletype(3) == 2
        hrtfdirmat = dir_ITA;
    end
end
%%%%%%%%%%%

% find locations for 2 or 3 sources
clear locind locations
numlocs = 0;

% 6/28/05 START

% offset_el = XStimParams.offset_el(2)/2;
% offset_az = XStimParams.offset_az(2)/2;
% offset_el3 = XStimParams.offset_el(3)/2;
% offset_az3 = XStimParams.offset_az(3)/2;
% if ~XStimParams.ThreeStims
%     for locnum = 1:size(XStimParams.locations,2)
%         loc1 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum) + offset_el) &...
%             hrtfdirmat(2,:) == (XStimParams.locations(2,locnum)) + offset_az));
%         loc2 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum) - offset_el) &...
%             hrtfdirmat(2,:) == (XStimParams.locations(2,locnum)) - offset_az));
%         if ~isempty(loc1) & ~isempty(loc2)
%             numlocs = numlocs+1;
%             locind(:,numlocs) = [loc1; loc2];
%             locations(:,numlocs) = XStimParams.locations(:,locnum);
%         end
%     end
% else
%     for locnum = 1:size(XStimParams.locations,2)
%         loc1 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum) + offset_el) &...
%             hrtfdirmat(2,:) == (XStimParams.locations(2,locnum)) + offset_az));
%         loc2 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum) - offset_el) &...
%             hrtfdirmat(2,:) == (XStimParams.locations(2,locnum)) - offset_az));
%         loc3 = max(find(hrtfdirmat(1,:) == (XStimParams.locations(1,locnum) + offset_el3) &...
%             hrtfdirmat(2,:) == (XStimParams.locations(2,locnum)) + offset_az3));
%         if ~isempty(loc1) & ~isempty(loc2) & ~isempty(loc3)
%             numlocs = numlocs+1;
%             locind(:,numlocs) = [loc1; loc2; loc3];
%             locations(:,numlocs) = XStimParams.locations(:,locnum);
%         end
%     end
%     
% end

% first express XStimParams.locations1 (the shape of the stimulus) as differences from the left highest point
maxy = max(XStimParams.locations1(1, :));
minx = min(XStimParams.locations1(2, find(XStimParams.locations1(1, :) == maxy)));
for pointnum = 1:size(XStimParams.locations1, 2)
	XStimParams.locations1(:, pointnum) = XStimParams.locations1(:, pointnum) - [maxy; minx];
end

% then assuming elevation is in the 1st row of hrtfdirmat
numlocs = size(XStimParams.locations1, 2);
numtrials = 0;

for locnum = 1:size(XStimParams.locations2, 2)
	count1 = 0;
	for pointnum = 1:size(XStimParams.locations1, 2)
		count2 = abs(XStimParams.locations2(1, locnum) + XStimParams.locations1(1, pointnum));
		count3 = abs(XStimParams.locations2(2, locnum) + XStimParams.locations1(2, pointnum));
		if (count2 + count3) <= 90
			locarray(pointnum) = max(find(hrtfdirmat(1, :) == (XStimParams.locations2(1, locnum)) +... 
			XStimParams.locations1(1, pointnum) & hrtfdirmat(2, :) == (XStimParams.locations2(2, locnum)) +... 
			XStimParams.locations1(2, pointnum)));
			count1 = count1 + 1;
		end
	end
	if count1 == size(XStimParams.locations1, 2)
		numtrials = numtrials + 1;
		locind(:, numtrials) = locarray';
		locations(:, numtrials) = XStimParams.locations2(:, locnum);
	end
end

if exist('locations') & numtrials
    XStimParams.locations = locations;
else
    error('no matching locations - reset and try again')
end

% 6/28/05 END

% get reference Lref, Rref and ABL for ABLalone test (always use just first ref location)
if XStimParams.ABLalone_flag 
    ind00 = max(find(hrtfdirmat(1,:) == XStimParams.el & hrtfdirmat(2,:) == XStimParams.az));
    if FN.HRTFfiletype(1) == 1
        Lref = mtlrch([FN.space_path FN.space],(2*ind00)-1);
        Rref = mtlrch([FN.space_path FN.space],2*ind00);
    else
        Lref = TF1_space(ind00,:)';
        Rref = TF2_space(ind00,:)';
    end
    tempL = ERBFilterBankB(Lref, fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
    tempR = ERBFilterBankB(Rref, fcoefs) .* Factormat;
    [ildref ablref] = calclevel_time(tempL,tempR, cF);
    %ITDref = calcitd(tempL,tempR, cF, Fs, ones(n_cF,1));
    clear tempL tempR temp
end

%%%%%%%%%%%%%%%%%%%% make the stimuli we'll use.  in this section are major modifications from 6/28/05.
remreps = 1;
set(H.status,'String','Status: Building Stimuli');
set(H.remreps,'String',num2str(remreps));
repnum = 1;
numtrials = size(XStimParams.locations,2);
finalspikematrix = zeros(1,numtrials);

%Randomize the stimuli
randseq = randperm(numtrials);
trialnum = 1;

while (exist1('H.multi_sourcefig') & (trialnum <= numtrials))
    set(H.status,'BackgroundColor','red');
    %Check for pause by user
    if pause_check
		return;
	end
 
	for current_loc = 1:numlocs
		%Make the stimulus
		
        %a few variables are used in non-XStimParam form and are not set by default elsewhere
        H.stim_type(current_loc) = XStimParams.stim_type_ted(current_loc);
        if H.stim_type(current_loc) == 9
            FN.stim_ted(current_loc) = XStimParams.stimfile(current_loc);
        end
        H.mod_type(current_loc) = XStimParams.mod_type_ted(current_loc);	
        if H.mod_type == 3
            FN.mod_ted(current_loc) = XStimParams.modfile(current_loc);
        end
                   
        if H.stim_type(current_loc) ~= 9
			source1_L = get_stim(XStimParams.freq(current_loc), current_loc, current_loc);
		else
			source1_L = get_stim(FN.stim_path, char(FN.stim_ted(current_loc)), current_loc);
		end
		source1_R = source1_L;

		% modulate stim1
		if XStimParams.mod_type_ted(current_loc) ~= 4
			Envelope = make_env(DUR, XStimParams.mod_type_ted(current_loc), XStimParams.mod_depth(current_loc), XStimParams.mod_freq(current_loc), XStimParams.mod_phase(current_loc), current_loc);
			source1_R = source1_R .* Envelope(:)';
			source1_L = source1_L .* Envelope(:)';
		end

		if H.stim_type(current_loc) ~= 9      %Ramp the stimuli
			ramp_time = 5; %ms
			[source1_L] = ramp_sound(source1_L,TDT.Fs,ramp_time);
			[source1_R] = ramp_sound(source1_R,TDT.Fs,ramp_time);
		end

		% remove any DCoffset
		source1_L = source1_L - mom(source1_L,1);
		source1_R = source1_R - mom(source1_R,1);

		%Apply ITD filtering if conducting ILDAlone Two_source Test
		if(XStimParams.ildalone_flag == 1)
			itdleft = 0; itdright = 0;
			ITD = round(str2num(get(H.ITD,'String')));
			if(ITD < 0)
				itdleft = 0;
				itdright = abs(ITD);
			elseif(ITD > 0)
				itdleft = abs(ITD);
				itdright = 0;
			end
			if(trialnum == 1)
				eval(['load ' FN.ITD_path 'itdfilt' num2str(itdleft)]);
				eval(['itd_filt_left = itd_filt' num2str(itdleft) ';']);
				eval(['load ' FN.ITD_path 'itdfilt' num2str(itdright)]);
				eval(['itd_filt_right = itd_filt' num2str(itdright) ';']);
			end
			source1_L = conv(source1_L,itd_filt_left);
			source1_R = conv(source1_R,itd_filt_right);
		end

		% Apply ILD filtering if conducting ITDalone Test
		if(XStimParams.itdalone_flag == 1)
			if(current_loc == 1)
                if(trialnum == 1)
					%%%%%%%%%%%%%%
					if FN.HRTFfiletype(2) == 1
						dir_ILA = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone]));
						ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
						if isempty(ILAind)
							disp('Could not find ILA location in HRTF file');
							return
						end
						eval(['ila_filt_left(1, :, :) = mtlrch(''' FN.ILA_path FN.ildalone ''', ' num2str(ILAind * 2-1) ');']);
						eval(['ila_filt_right(1, :, :) = mtlrch(''' FN.ILA_path FN.ildalone ''', ' num2str(ILAind * 2) ');']);
					else
						ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
						if isempty(ILAind)
							disp('Could not find ILA location in HRTF file');
							return
						end
                        ila_filt_left(1, :, :) = TF1_ila(ILAind,:);
						ila_filt_right(1, :, :) = TF2_ila(ILAind,:);
					end
					%%%%%%%%%%%%
				end
        		source1_L = conv(source1_L,ila_filt_left(1, :, :));
				source1_R = conv(source1_R,ila_filt_right(1, :, :));
            else
                if(trialnum == 1)
                    ILAind = max(find(dir_ILA(1,:) == XStimParams.el & dir_ILA(2,:) == XStimParams.az));
                    if isempty(ILAind)
                        disp('Could not find ILA location in HRTF file');
                        return
                    end
                    if FN.HRTFfiletype(2) == 1
                        eval(['ila_filt_left(current_loc, :, :) = mtlrch(''' FN.ILA_path FN.ildalone ''', ' num2str(ILAind * 2-1) ');']);
                        eval(['ila_filt_right(current_loc, :, :) = mtlrch(''' FN.ILA_path FN.ildalone ''', ' num2str(ILAind * 2) ');']);
                    else
                        ila_filt_left(current_loc, :, :) = TF1_ila(ILAind,:);
                        ila_filt_right(current_loc, :, :) = TF2_ila(ILAind,:);
                    end   
                end
                source1_L = conv(source1_L,ila_filt_left(current_loc, :, :));
                source1_R = conv(source1_R,ila_filt_right(current_loc, :, :));
            end
        end

		%Add in the leading silent period
		source1_L =  [zeros(1,silence_len) source1_L];
		source1_R = [zeros(1,silence_len) source1_R];

		%Add in the trailing silent period
		source1_L =  [source1_L zeros(1,silence_len2)];
		source1_R = [source1_R zeros(1,silence_len2)];
       
		%Apply HRTF filtering
		if(XStimParams.space_flag == 1)
			if FN.HRTFfiletype(1) == 1
				hrtf_left = mtlrch([FN.space_path FN.space],(2*locind(current_loc,randseq(trialnum)))-1);
				hrtf_right = mtlrch([FN.space_path FN.space],2*locind(current_loc,randseq(trialnum)));
			else
				hrtf_left = TF1_space(locind(current_loc,randseq(trialnum)),:);
				hrtf_right = TF2_space(locind(current_loc,randseq(trialnum)),:);
			end
		elseif(XStimParams.ABLalone_flag == 1)
			% get ABLfactor for this location
			if FN.HRTFfiletype(1) == 1
				tempL = mtlrch([FN.space_path FN.space],(2*locind(current_loc,randseq(trialnum)))-1);
				tempR = mtlrch([FN.space_path FN.space],2*locind(current_loc,randseq(trialnum)));
			else
				tempL = TF1_space(locind(current_loc,randseq(trialnum)),:);
				tempR = TF2_space(locind(current_loc,randseq(trialnum)),:);
			end
			tempL = ERBFilterBankB(tempL, fcoefs) .* Factormat;		% has dimensions n_cF x length(noi)
			tempR = ERBFilterBankB(tempR, fcoefs) .* Factormat;
			[ildx ablx] = calclevel_time(tempL,tempR, cF);
			%ITDx = calcitd(tempL,tempR, cF, Fs, ones(n_cF,1));
			ABLfactor = (10 .^(ablx/20)) ./ (10 .^(ablref/20));
			%%%%%%%%%%%%%%%%%% next line
			ABLfactor_long = (interp1(freq(freq_ind),ABLfactor,freq(startpt:endpt)))';
			% apply ABLfactor to reference location
			FT_L = fft(Lref,fftpts);
			FT_L(Xpart1) = FT_L(Xpart1) .* (ABLfactor_long);       	% positive freqs
			FT_L(Xpart2) = FT_L(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
			hrtf_left = real(ifft(FT_L));
			hrtf_left = hrtf_left(1:255);
			FT_R = fft(Rref,fftpts);
			FT_R(Xpart1)= FT_R(Xpart1) .* (ABLfactor_long);       	% positive freqs
			FT_R(Xpart2)= FT_R(Xpart2) .* flipud(ABLfactor_long);   % negative freqs
			hrtf_right = real(ifft(FT_R));
			hrtf_right = hrtf_right(1:255);
			clear FT_L FT_R ABLfactor*         
		elseif(XStimParams.ildalone_flag == 1)
			if FN.HRTFfiletype(2) == 1
				hrtf_left = mtlrch([FN.ILA_path FN.ildalone],(2*locind(current_loc,randseq(trialnum)))-1);
				hrtf_right = mtlrch([FN.ILA_path FN.ildalone],2*locind(current_loc,randseq(trialnum)));
			else
				hrtf_left = TF1_ila(locind(current_loc,randseq(trialnum)),:);
				hrtf_right = TF2_ila(locind(current_loc,randseq(trialnum)),:);
			end
		elseif(XStimParams.itdalone_flag == 1)
			if FN.HRTFfiletype(3) == 1
				hrtf_left = mtlrch([FN.ITA_path FN.itdalone],(2*locind(current_loc,randseq(trialnum)))-1);
				hrtf_right = mtlrch([FN.ITA_path FN.itdalone],2*locind(current_loc,randseq(trialnum)));
			else
				hrtf_left = TF1_ita(locind(current_loc,randseq(trialnum)),:);
				hrtf_right = TF2_ita(locind(current_loc,randseq(trialnum)),:);
			end
		end
		source1_L = conv(source1_L,hrtf_left);
		source1_R = conv(source1_R,hrtf_right);

		% scale
		source1_L = source1_L - round(mean(source1_L));
		source1_R = source1_R - round(mean(source1_R));
		if XStimParams.ildalone_flag | XStimParams.itdalone_flag
			ABAval = 0.5*(mom(source1_L,2) + mom(source1_R,2));
			scalefact = TDT.scalevalue/ABAval;
			source1_L = round(scalefact*source1_L);
			source1_R = round(scalefact*source1_R);
		end

		if current_loc == 1
			trial_left = source1_L * XStimParams.factor(current_loc);
			trial_right = source1_R * XStimParams.factor(current_loc);
		else
			trial_left = trial_left + (source1_L * XStimParams.factor(current_loc));
			trial_right = trial_right + (source1_R * XStimParams.factor(current_loc));
		end
	end
		
	%%%%%%%%%%%%%%%%% add the weighted sources
    trial_left = trial_left / numlocs;
    trial_right = trial_right / numlocs;
       

    % re-scale
    ABAval = 0.5*(mom(trial_left,2) + mom(trial_right,2));
    scalefact = TDT.scalevalue/ABAval;
    trial_left = round(scalefact*trial_left);
    trial_right = round(scalefact*trial_right);
    
    %pad with zeros
    filttrial_left = [trial_left zeros(1,TDT.ephonefiltlen)];
    filttrial_right = [trial_right zeros(1,TDT.ephonefiltlen)];
    
    % save stims to disk with name of loc1
    if(exist1('H.multi_sourcefig'));
        S232('push16',filttrial_left,length(filttrial_left));
        S232('qpop16',BUF.L1);
        fname = ['stimbuf_left_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum)))) ...
                '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum))))];
        evalstr = ['S232(''dama2disk16'',BUF.L1,' ...
                [' ''' stimuli_dir fname ''' ']   ',0);'];
        eval(evalstr);
        temp_left = dama2pc(BUF.L1);
        S232('push16',filttrial_right,length(filttrial_right));
        S232('qpop16',BUF.R1);
        fname = ['stimbuf_right_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum)))) ...
                '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum))))];
        evalstr = ['S232(''dama2disk16'',BUF.R1,' ...
                [' ''' stimuli_dir fname ''' ']   ',0);'];
        eval(evalstr);
        temp_right = dama2pc(BUF.R1);
    end
    
    remtrials = numtrials - trialnum;
    set(H.remtrials,'String',num2str(remtrials));
    trialnum = trialnum + 1;
    set(H.status,'BackgroundColor','blue');
    pause(0);
end 										%end loop over trials

%%%%%%%%%%%%%%%%%%%%%%%%%%% finished making sounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin playing sounds   
set(H.status,'String','Status: Playing Stimuli');
set(H.status,'BackgroundColor','green');
set(H.remreps,'String',num2str(XStimParams.numreps));
repnum = 1;
datamatrix = [];

% increment testnumber
if(exist1('H.multi_sourcefig') & get(H.recorddata,'Value'))
    update_dataFN;
end

%loop for reps
while (exist1('H.multi_sourcefig') & (repnum <= XStimParams.numreps))
    %Randomize the stimuli
    randseq = randperm(numtrials);   
    trialnum = 1;
    spikes_trial = [];
    EL_trial = [];
    AZ_trial = [];
    repnum_trial = [];
    Nspikes = [];
    
    % loop for trials
    tic
    while (exist1('H.multi_sourcefig') & (trialnum <= numtrials+1))
        
        %Check for pause by user
        if pause_check  return; end
        
        %Wait till PD1 is finished
        while S232('PD1status',1) usec_delay(1000); end
        
        if(trialnum <= numtrials)
            fname = ['stimbuf_left_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum)))) ...
                    '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum))))];
            evalstr = ['S232(''disk2dama16'',BUF.L1,'  [' ''' stimuli_dir fname ''' '] ',0);'];
            eval(evalstr);
            fname = ['stimbuf_right_' num2str(hrtfdirmat(1,locind(1,randseq(trialnum)))) ...
                    '_' num2str(hrtfdirmat(2,locind(1,randseq(trialnum))))];
            evalstr = ['S232(''disk2dama16'',BUF.R1,' [' ''' stimuli_dir fname ''' '] ',0);'];
            eval(evalstr);
        end
        
        %Wait till PD1 is finished
        while S232('PD1status',1) usec_delay(1000); end
        S232('PD1stop',1);
        
        %Stop the m110 and get spikes
        if(trialnum > 1)							% first trial just for loading
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, 1000); 			% Take 100 spikes max
            ind = find(spikes ~= 0); 						% Get clock events that are spikes
            spikes = spikes(ind);
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            if(exist1('H.multi_sourcefig') & ~isempty(spikes)) 
                spikes = spikes(ind);
                spikes_trial = [spikes_trial;spikes/(1000/mii_us_per_sample)];
                EL_trial = [EL_trial;locations(1,randseq(trialnum-1))* ones(size(spikes))];
                AZ_trial = [AZ_trial;locations(2,randseq(trialnum-1))* ones(size(spikes))];
                repnum_trial = [repnum_trial;repnum * ones(size(spikes))];
                Nspikes = [Nspikes; length(spikes) * ones(size(spikes))];
            end
        end
        
         %Check for pause by user
        if pause_check  return; end

        if(exist1('H.multi_sourcefig') & (trialnum <= numtrials))
            
            S232('seqplay',BUF.playspec1);
            S232('PD1arm',1);
            
            %Send trigger
            %Set up MII
            m100x( C_.INIT );
            m110dx( C_.INIT );
            m110dx( C_.CLOCK, mii_us_per_sample);
            m110dx( C_.MODE, M_.PST );
            
            if (trialnum <= numtrials)
                while toc < ISI/1000     end
                %Start clock
                m110dx( C_.START);
                %Send pulse: PD1 GO!
                m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
                tic
            end
            
        end
        
        if(trialnum > 1)
            finalspikematrix(randseq(trialnum-1)) = ...
                finalspikematrix(randseq(trialnum-1)) + ...
                length(spikes);
        end
        
        remtrials = numtrials - trialnum +1;
        set(H.remtrials,'String',num2str(remtrials));
        trialnum = trialnum + 1;
        pause(0);
    end %end loop over trials
    
    %Plot Spike Rate Data
    interimspikerate = finalspikematrix/repnum;
    if(exist1('H.multi_sourcefig') & ~exist1('H.finalspikeratefig'))
        H.finalspikeratefig = figure('Position',[700 20 550 500],...
            'Name','Multisource Test Spike Rate Plot',...
            'NumberTitle','off');
        H.spikeaxes = axes;
    end
    figure(H.finalspikeratefig)
    plotdiam1(XStimParams.locations, interimspikerate);
    set(H.spikeaxes,'Color','black');
    xlabel('Azimuth'); ylabel('Elevation'); title(['Rep # ' num2str(repnum)]);
    colorbar
    
    %Record Data
    if(exist1('H.multi_sourcefig') & get(H.recorddata,'Value'))
        tempseq{repnum} = randseq;
        datamatrix = [datamatrix;[Nspikes spikes_trial repnum_trial EL_trial AZ_trial]];
        record_data3(XStimParams,datamatrix, tempseq);
    end
    
    remreps = XStimParams.numreps - repnum;
    set(H.remreps,'String',num2str(remreps));
    repnum = repnum + 1;
    pause(0);
end 									%end loop over reps

%Plot final spike rate figure
finalspikematrix = finalspikematrix/XStimParams.numreps;
figure(H.finalspikeratefig)
set(H.finalspikeratefig,'Name','Final Plot for Multisource Test');
plotdiam1(XStimParams.locations, interimspikerate);
set(H.spikeaxes,'Color','black');
locmaxspikes = find(finalspikematrix == max(finalspikematrix));
xlabel('Azimuth'); ylabel('Elevation');
title(['Maximum Activity at EL = ' num2str(XStimParams.locations(1,locmaxspikes)) ...
        ', AZ = ' num2str(XStimParams.locations(2,locmaxspikes))],...
    'FontSize',8);
colorbar

set(H.status,'BackgroundColor','blue');
set(H.status,'String','Status: Results');
set(H.exitTwo_source,'Visible','on');
set(H.resetTwo_source,'Visible','on');

% increment test number
if(exist1('H.multi_sourcefig') & get(H.recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
end

%%%%%%%%%
function [stim] = get_stim(param1, param2, cloc)
% param1: XStimParams.freq or FN.stim_path
% param2: FN.stim

global H
global XStimParams
global TDT

switch H.stim_type(cloc)     
    case 1    %tone at specified frequency 
        stim = MakeTone(TDT.Fs,param1,XStimParams.curr_stimdur);
    case 2    %GammaTones
        stim = MakeGammaTone(TDT.Fs,param1,XStimParams.curr_stimdur);
    case 8    %Broadband Noise
        [stim] = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
    case 9    %Stimulus from file
        if(~exist1('stim_from_file'))
            fid = fopen([param1 param2],'r');
            stim = fread(fid,inf,'float');
            fclose(fid);
        end
        if(size(stim,1) > 1)
            stim = stim';
        end
    otherwise
        H.stim_type(cloc) = 8;
        stim = MakeBBNoise(TDT.Fs,XStimParams.curr_stimdur);
        disp('Stimulus type not supported for multisource tests.  Reset to BROADBAND');
        return
end

%%%%%%
function [Envelope] = make_env(DUR,mod_type,param1,param2,param3,param4)

% param1 used as follows:
%    'tone' or 'LP noise' or 'File': mod_depth
% param2 used as follows:
%    'tone or 'LP noise' : mod_freq
% param3 used as follows
%    'tone' : addition to modulation starting phase (0 - pi)

global H
global TDT
global FN

if nargin < 3   param3 = 0; end

Npts = DUR*(round(TDT.Fs/1000));
Inc = 1/TDT.Fs;
switch mod_type
    case 1
        T = 0:Inc:(DUR/1000 - Inc);
        Tone = (param1 / 2)* sin(2 * pi * param2 .* T + (.75 * 2 * pi + param3));
        Envelope = Tone + (1-param1/2);
    case 2
        LP_noise = m_noi(5, param2, Npts/(TDT.Fs), TDT.Fs/2);
        LP_noise = (param1 / 2)* (LP_noise / max1(LP_noise));
        Envelope = LP_noise + (1-param1/2);
    case 3			
        fid = fopen([FN.mod_path char(FN.mod_ted(param4))],'r');
        mod_from_file = fread(fid,inf,'float');
        fclose(fid);
        while length(mod_from_file) ~= Npts
            [tempstr1, FN.mod_path] = uigetfile('*.*','Select Envelope File');
            FN.mod_ted(param4) = cellstr(tempstr1);
            if(FN.mod_path ~= 0)
                set(H.modfile,'String',[FN.mod_path char(FN.mod_ted(param4))]);
            end
            fid = fopen(char(FN.mod_ted(param4)), 'r');
            mod_from_file = fread(fid,inf,'float');
            fclose(fid);
        end
        mod_from_file = mod_from_file - mean(mod_from_file);
        mod_from_file = (param1 / 2)* (mod_from_file / max1(mod_from_file));
        Envelope = mod_from_file + (1-param1/2);
    otherwise
end


%%%%%%%%%
function [flag] = pause_check

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

flag = 0;
%Check for pause by user
while (exist1('H.multi_sourcefig') & get(H.pauseTwo_source,'Value'))
    pause(0);
    if(~exist1('H.multi_sourcefig')) return; end         
    set(H.exitTwo_source,'Visible','on');
    set(H.resetTwo_source,'Visible','on');
    if(exist1('H.multi_sourcefig') & get(H.resetTwo_source,'Value') == 1)
        set(H.resetTwo_source,'Value',0);
        set(H.pauseTwo_source,'Value',0);
        Reset_multisource;   flag = 1;
        return;
    end
    if isempty(XStimParams.locations1)
        Reset_multisource;   flag=1;
        return;
    end
end

if XStimParams.reset_flag
    flag = 1;
    XStimParams.reset_flag = 0;
end

