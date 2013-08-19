function [] = Engage_space_RIR()

global H
global XStimParams
global TDT
global FN
global C_
global M_
global GUI

%Engage_space_RIR
%*******************************************************************************
%	The space_RIR Test operation
% use RIRs without earphone equalization - usually called *.std
% in saved DATA, the (ASCII of the) first letter of RoomType is saved as
% the first column in DATA. Also param3 is saved as a cellarray of
% roomtypes fpr each sound played
%
%*******************************************************************************
colors = [ ...
        179 199 255; ...
        200 150 255; ...
        0    0  255; ...
        216 41  0; ...
        255 199 179;...
        255 150 200;...
        255   0    0;...
        199 255 179;...
        200 255 150;...
        0   255 0]/255;

fclose all;

if XStimParams.HiDynamicRange
    scaleFactor = 60;
else
    scaleFactor = TDT.scaleFactor;
end

% reset stim_type to file (#9)
set(H.stim_type,'Value',9);

XStimParams.curr_stimdur = str2num(get(H.space_RIR_DUR,'String'));
XStimParams.test_ISI = str2num(get(H.space_RIR_ISI,'String'));
XStimParams.numreps = str2num(get(H.space_RIR_numreps,'String'));
XStimParams.reset_flag = 0;

%Specify DAMA buffers
clear BUF
BUF.L1				= 1;
BUF.R1				= 2;
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

% check to be sure stimFNs chosen
if isempty(FN.RIR)         %Stimulus from file
    set(H.space_RIR_initFN,'value',1)
    setinfo_space_RIR
end
nFiles = length(FN.RIR);
for ifile = 1:nFiles
    temp = char(FN.RIR(ifile));
    ind = strfind(temp,'rep')+3;
    RIRrepNum(ifile) = str2num(temp(ind:ind+1));
end
[RIRrepNum, indRIRsort] = sort(RIRrepNum);
FN.RIR = FN.RIR(indRIRsort);
XStimParams.numreps = length(unique(RIRrepNum));
set(H.space_RIR_numreps,'string',num2str(XStimParams.numreps))

%Make Stimulus buffers
% set Npts to length of file
Npts_totalplay = XStimParams.RIRpts2use;

XStimParams.test_type = 'space_RIR FC';     str1 = 'rir';
disp('This is a FULLY CUED space_RIR test')

S232('allot16',BUF.L1,Npts_totalplay);
S232('allot16',BUF.R1,Npts_totalplay);

S232('PD1clear',1);
S232('PD1srate', 1, 1e6/TDT.Fs);
S232('PD1npts',1,Npts_totalplay);

%Get Earphone filters
if FN.HRTFfiletype(6) == 1
    ephonefname = [FN.ephone_path FN.ephone2];
    ephonefilt_left  = (mtlrch(ephonefname,1))';
    ephonefilt_right = (mtlrch(ephonefname,2))';
else
    eval(['load -mat ' FN.ephone_path FN.ephone2]);
    ephonefilt_left  = TF1;
    ephonefilt_right = TF2;
    clear TF1 TF2 dir
end

%Load Earphone filters
dspid_left = 0; dspid_right = 1;
S232('PD1clrsched',1);
S232('PD1nstrms',1,2,0);
S232('PD1resetDSP',1,hex2dec('FFF'));           % used to be '0xFFF'
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

ABL = XStimParams.curr_ABL;
if(ABL < -110) return; end
S232('PA4atten',1,abs(ABL)-20);					% HB1 adds 20 dB attenuation
S232('PA4atten',2,abs(ABL)-20);

ISI = XStimParams.test_ISI;
ISI = ISI - (TDT.itdfiltlen/(TDT.Fs/1000)); 				%correct for ITD filtlength
ISI = ISI - (TDT.ephonefiltlen/(TDT.Fs/1000));			%correct for ephonefilt

% save XStimParams for this test
if get(H.space_RIR_recorddata,'Value')
    tempstr = ['    ' str1 ' space_RIR-test     atten: ' num2str(abs(ABL))];
    update_diary
end

Temp_params = XStimParams;
Temp_FN_RIR = FN.RIR;
eval(['save ' FN.temp_stim_path 'XStimParams_space_RIR_' str1 ' Temp_params Temp_FN_RIR'])
clear Temp_params Temp_FN_RIR str1

% update display
setInfo_space_RIR;

nLocs = size(XStimParams.locations1,2);
nFiles = size(FN.RIR,2);
RoomTemp = [];
for ifile = 1:nFiles
    fname = char(FN.RIR(ifile));
    ind = findstr(fname,'_');
    RoomTemp = [RoomTemp fname(ind(3)+1)];
end
RoomList = unique(RoomTemp);
nRooms = length(RoomList);

%Begin playing sounds   
set(H.space_RIR_status,'String','Status: playing stimuli');
set(H.space_RIR_status,'BackgroundColor','green');
set(H.space_RIR_remreps,'String',num2str(XStimParams.numreps));

% increment testnumber
if(exist1('H.space_RIRfig') & get(H.space_RIR_recorddata,'Value'))
    update_dataFN;
end

% initalize dataRaster figure
%Plot Spike Raster
if(exist1('H.finalspikerasterfig')) close(H.finalspikerasterfig);  H.finalspikerasterfig =[]; end
H.finalspikerasterfig = figure('Position',[700 20 550 500],...
    'Name','RIR Test Spike Raster',...
    'NumberTitle','off');
H.spikeaxes = axes;
set(H.spikeaxes,'Color','black');
hold on
axis([0 XStimParams.curr_stimdur + XStimParams.silence_lead + XStimParams.silence_trail 0 1])
set(H.spikeaxes,'Ytick',0:1/nLocs/nRooms:1);
ipt=1;
str = repmat(' ',1,round(40/(nRooms*nLocs-1)));
templabel2 = str;
for iLoc=1:nLocs
    for iRoom = 1:nRooms
        templabel{ipt} = ['loc' num2str(iLoc) ' ' RoomList(iRoom)];
        templabel2 = [templabel2 templabel{ipt}  str];
        ipt=ipt+1;    
    end
end
templabel{ipt} = ' ';
set(H.spikeaxes,'YtickLabel','');

ylabel(templabel2)
xlabel('time (ms)')
title('RIR spike raster plot')

plot([1 1]*XStimParams.silence_lead,[0 nFiles],'g')
plot([1 1]*XStimParams.silence_lead+1000,[0 nFiles],'g')
plot([1 1]*(XStimParams.silence_lead+XStimParams.curr_stimdur),[0 nFiles],'r')

%%%%%%%%%%%%%%%%%%%%%%%% main data collection
datamatrix = [];
spikes_trial = [];
EL_trial = [];
AZ_trial = [];
Room_trial = [];
FN_trial = [];
Nspikes = [];
seed_trial = [];
iLoc = zeros(nFiles,1);
spike_raster_matrix = cell(nFiles);

%%%% loop for reps (max of 50 reps- that's all the files available)
irep=0; jrep=0; jfile = 0; last_jFile = 0;
while (exist1('H.space_RIRfig') & (irep <= max1(RIRrepNum)))
    irep = irep+1;                              % for counting up through filenames
    indReps = find(RIRrepNum == irep);
    if ~isempty(indReps)
        jrep = jrep+1;                          % which rep is actually being played
        indReps = indReps(randperm(length(indReps)));
        tic
        for ifile = 1:length(indReps)
            jfile= jfile+1;
            %Check for pause by user
            if pause_check    return; end
            
            set(H.space_RIR_status,'String',['Loading file']);
            set(H.space_RIR_status,'BackgroundColor','blue');
            set(H.space_RIR_status,'ForegroundColor','yellow');
            
            %%% load stimuli
            fname = char(FN.RIR(indReps(ifile)));
            eval(['load ' FN.stim_path2 fname ' t* n* seed STIMparams'])
            
            ind = findstr(fname,'el')-4;
            ele = str2num(fname(ind+1:ind+3));
            if strcmp(fname(ind),'n')   ele=-ele;   end
            ind = findstr(fname,'az')-4;
            azi = str2num(fname(ind+1:ind+3));
            if strcmp(fname(ind),'n')   azi = -azi; end
            iLoc(jfile) = find(XStimParams.locations1(1,:)==ele & XStimParams.locations1(2,:)==azi);
            
            ind = findstr(fname,'_');
            RoomType{jfile} = fname(ind(3)+1:ind(4)-1);
            
            if nptsTotalPlay ~= Npts_totalplay
                error('number of points in stim file incorrect')
            end
            % load to buffers
            S232('push16',trial_left,Npts_totalplay);
            S232('qpop16',BUF.L1);
            S232('push16',trial_right,Npts_totalplay);
            S232('qpop16',BUF.R1);
            
            set(H.space_RIR_status,'String',['Playing: ' fname]);
            set(H.space_RIR_status,'BackgroundColor','blue');
            set(H.space_RIR_status,'ForegroundColor','green');
            
            S232('seqplay',BUF.playspec1);
            S232('PD1arm',1);
            
            %Send trigger
            %Set up MII
            m100x( C_.INIT );
            m110dx( C_.INIT );
            m110dx( C_.CLOCK, mii_us_per_sample);
            m110dx( C_.MODE, M_.PST );
            
            while toc < ISI/1000 
                if pause_check    return; end
            end
            
            %Start clock
            m110dx( C_.START);
            %Send pulse: PD1 GO!
            m101x( C_.DATA,M_.BIT,M_.PULSE,0); %Use port 0 for the pulse
            tic
            
            pause(Npts_totalplay/TDT.Fs + .2);
            
            while(S232('PD1status',1)) usec_delay(1000);  end
            S232('PD1stop',1);
            
            %Stop the m110 and get spiketimes
            m110dx( C_.STOP);
            spikes = m110dx( C_.DATA, round(XStimParams.curr_stimdur*2)); 			% Take 2*XStimParams.curr_stimdur spikes max
            ind = find(spikes ~= 0); 						% Get clock events that are spikes
            spikes = spikes(ind);
            ind = [1; (find(diff(spikes) > mii_separation))+1]; %Only take those events separated by >=1 ms
            if ~isempty(spikes) 
                spikes = spikes(ind);
                spike_raster_matrix{jfile} = spikes/(1000/mii_us_per_sample);
                
                Ones1 = ones(size(spikes));
                if exist1('H.space_RIRfig')
                    spikes_trial = [spikes_trial; spikes/(1000/mii_us_per_sample)];
                    EL_trial = [EL_trial; ele * Ones1];
                    AZ_trial = [AZ_trial; azi * Ones1];
                    Room_trial = [Room_trial; double(uint8(RoomType{jfile}(1)))* Ones1];
                    FN_trial = [FN_trial; indReps(ifile)* Ones1];
                    Nspikes = [Nspikes; length(spikes) * Ones1];
                    seed_trial = [seed_trial; seed * Ones1];
                end
            end
            
            if pause_check    return; end
            
            set(H.space_RIR_remtrials,'String',num2str(nFiles-jfile));
            pause(0);
            
            %Record Data
            if(exist1('H.space_RIRfig') & get(H.space_RIR_recorddata,'Value'))
                datamatrix = [Nspikes spikes_trial Room_trial EL_trial AZ_trial FN_trial seed_trial];
                record_data3(XStimParams,datamatrix,RoomType,STIMparams);
            end
        end                 % ifiles
    end                     % ~isempty
    
    figure(H.finalspikerasterfig);
    for ifile = last_jFile+1:jfile
        tempType = RoomType{ifile}(1);
        switch tempType
            case 'N'
                iRoom = 10;
            case 'S'
                iRoom = 2;
            case 'M'
                iRoom = 5;
            case 'H'
                iRoom = 4;
            otherwise
                iRoom = 6;
        end
        ypt = (iLoc(ifile)-1)/nLocs + (strfind(RoomList,tempType)-1)/nRooms/nLocs + (jrep-1)/XStimParams.numreps/nRooms/nLocs;
        
        temp = spike_raster_matrix{ifile}(spike_raster_matrix{ifile}<XStimParams.silence_lead | spike_raster_matrix{ifile}>XStimParams.silence_lead +XStimParams.curr_stimdur);
        plot(temp,ones(size(temp))*ypt,'b.')
        temp = spike_raster_matrix{ifile}(spike_raster_matrix{ifile}>XStimParams.silence_lead & spike_raster_matrix{ifile}<XStimParams.silence_lead +XStimParams.curr_stimdur);
        plot(temp,ones(size(temp))*ypt,'.','color', colors(iRoom,:))
    end
    last_jFile = jfile;
    drawnow
    
    set(H.space_RIR_remreps,'String',num2str(XStimParams.numreps - jrep));
    pause(0);
end 						% reps

set(H.space_RIR_status,'String','Status: results');
set(H.space_RIR_status,'BackgroundColor','blue');

set(H.exitspace_RIR,'Visible','on');
set(H.resetspace_RIR,'Visible','on');

% increment test number
if(exist1('H.space_RIRfig') & get(H.space_RIR_recorddata,'Value'))
    XStimParams.testnum = XStimParams.testnum +1;
    set(H.testnum, 'String', num2str(XStimParams.testnum));
    update_dataFN;
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
while (exist1('H.space_RIRfig') & get(H.pausespace_RIR,'Value'))
    pause(0);
    if(~exist1('H.space_RIRfig')) return; end         
    set(H.exitspace_RIR,'Visible','on');
    set(H.resetspace_RIR,'Visible','on');
    if(exist1('H.space_RIRfig') & get(H.resetspace_RIR,'Value') == 1)
        set(H.resetspace_RIR,'Value',0);
        set(H.pausespace_RIR,'Value',0);
        Reset_space_RIR;    flag = 1;
        return;
    end
end
if XStimParams.reset_flag ==1
    flag = 1;
    XStimParams.reset_flag = 0;
end