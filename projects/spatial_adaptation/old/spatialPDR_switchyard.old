function [varargout] = spatialPDR_switchyard(fcn,varargin)
% Call as: spatialPDR_switchyard('function name');
global PDR session

switch fcn
    case 'MTLreadCH'
        % PROTOTYPE: channel = MTLreadCH(index)
        index = varargin{1};
        channel = feval(fcn,index);
        varargout{1} = channel;
    case 'sphere2double'
        % PROTOTYPE: [dir1, dir2] = sphere2double(p1, p2)
        dirIn = varargin{1};
        dirOut = feval(fcn,dirIn);
        varargout{1} = dirOut;
    case 'setupCalib'
        flag_adapt=varargin{1};
        feval(fcn,flag_adapt);
    case 'setupTrialSeq'
        quit = feval(fcn);
        varargout{1} = quit;
    case 'makeTest'
        % PROTOTYPE: stim = makeTest(seedval,dur,minband,maxband,samprate,ramp,base_rms)
        seedval = varargin{1};
        dur = varargin{2};
        minband = varargin{3};
        maxband = varargin{4};
        samprate = varargin{5};
        ramp = varargin{6};
        base_rms = varargin{7};
        stim = feval(fcn);
        varargout{1} = stim;
    case 'run_calibration'
        quit_flag = feval(fcn);
        varargout{1} = quit_flag;
    case 'calc_dB'
        bird = varargin{1};
        levels = varargin{2};
        coords = varargin{3};
        scales = feval(fcn,bird,levels,coords);
        varargout{1} = scales;
    otherwise
        feval(fcn); % no special inputs or outputs
end

function setDefaults()
global PDR

PDR = struct(...                            % MAIN PARAMETERS:
    'calib',0,...                           % flag for calibrations
    'virtual',1, ...                        % flag for virtual sound presentation
    'record',1, ...                        % flag for recording pupillometer output
    'flag_adapt',1,...                     % flag for adapting stimulus
    'ntrials',50,...                       % # of trials in session
    'npretrials',5,...                     % # of trials before first test trial
    'n_test_trials',[],...                  % # of test trials 
    'buf_pts',16384,...                     % number of pts in each buffer
    'buf_dur',[],...                        % buffer duration in ms
    'isi_buf',3,...                         % no. buffers between trials
    'isi_time',[],...                       % ISI (seconds) between trials
    'decimationfactor',4, ...               % decimation factor for data collection
    'stim_Fs', 30000, ...                   % sampling rate (Hz) for stimuli
    'npts_totalplay',[],...                % total no. pts in the session
    'len_session',[], ...                   % length of session (in minutes)
    'starttime',[], ...                     % session start time
    'stoptime', [], ...                     % session stop time
    'code_path', 'C:\andrew\pdr\code\',...  % path to code
    'data_path', 'C:\andrew\pdr\data\', ...     % data path (normally: c:\andrew\pdr\data normally, calib: c:\andrew\pdr\calib_data)
    'base_atten',0,...                      % base attenuation for PA4
    'filename',[],...                       % filename for session data
    'bird_id',929,...
    ...
    ...                                     % ADAPTING STIMULUS PARAMETERS:
    'ADAPT_loc',[0 0],...                   % location of adapting stimulus
    'ADAPT_dB',[],...                       % TBD
    'ADAPT_species','owl',...               % species
    'ADAPT_type','gammatone',...            % type of sound for adapting stimulus
    'ADAPT_cF',6000,...                     % center frequency (Hz) (if relevant)
    'ADAPT_coefs',[],...                    % FIR coefficients to generate adapting stimulus
    'ADAPT_scale',30000,...                  % test scale for adaptor
    'ADAPT_ramp',[],...                     % ramp for trial segments
    'ADAPT_seed',47,...                     % seed for generating adaptor (using C code)
    ...
    ...                                     % TEST STIMULUS PARAMETERS:
    'TEST_soundtype','octave',...           % type of test sound ('gammatone', 'octave', 'broadband', etc.)
    'TEST_seed',47,...                      % seed value to for random number generator
    'TEST_base_rms',0.0350,...              % rms level for test sound
    'TEST_bandwidth',[4000 8000],...            % frequency range (min max)
    'TEST_dur',100,...                      % duration of test sounds (ms)
    'TEST_ramp',5,...                       % ramp duration for test sound (should be 5ms)
    'TEST_on_delay_pts',[],...              % delay before test sound onset (TBD)
    'TEST_trial_freq',3, ...                % test sound every x trials
    'TEST_trial_jitter',0, ...              % maximum jitter in isi_buf
    'TEST_scales',[1000 5000 15000 32760],...                 %[10 20 60 80 500],... % test scales (try 10, 20 ,40, 60, 80, 200, 300, 500, 1000)
    'TEST_sound',[],...                     % test sound stored here
    'TEST_scale_sequence',[],...            % sequence of test sound scales
    'TEST_nlocs',5,...                     % no. of test locations (default is 3)
    'TEST_locs',[0 0; 0 20; 0 40; 0 60; 0 80],...      % locations for test stimuli (el,az)
    'TEST_loc_sequence',[],...              % sequence of test sound locations
    ...
    ...                                     % HRTF PARAMETERS:
    'HRTF_directory','C:\andrew\pdr\HRTF_lib\',...  % directory of HRTF coefficient files
    'HRTF_fname','929AD_ABLequal.eq',...            % HRTF filename
    'HRTF_nchannels', 2,...                 % number of channels
    'HRTF_nlines', 255,...                  % number of lines (HRIR taps)
    'HRTF_sampling_rate', 30000);           % sampling rate (same as stimulus Fs)

% buffer duration (ms):
PDR.buf_dur = 1000*PDR.buf_pts/PDR.stim_Fs;

% Time (in seconds) between test stimuli:
PDR.isi_time = ( PDR.buf_dur*PDR.isi_buf+(PDR.buf_dur-PDR.TEST_dur)/2 )/1000;

q=clock;
y=num2str(q(1));y=y(3:4);
m=num2str(q(2));if size(m,2)<2;m=['0' m];end
d=num2str(q(3));if size(d,2)<2;d=['0' d];end
switch PDR.bird_id
    case 929
        LETTER = 'l';
    case 930
        LETTER = 'm';
    case 0
        LETTER = 't';
end

PDR.filename = [y m d LETTER 'a'];  %930=m 929=l 882=d 883=e  (a is index to experiment number for that day)
PDR.npts_totalplay = PDR.ntrials*(PDR.isi_buf+1)*PDR.buf_pts; % Calculate length of session!
PDR.len_session = (1/60)*(PDR.npts_totalplay/PDR.stim_Fs); % length of session in minutes
h=msgbox(['Session will last approximately ' num2str(PDR.len_session) ' minutes']);
uiwait(h)
cd(PDR.code_path);
tmp = pwd;
if tmp(end-3:end) == 'code'
    PDR.data_path = [tmp(1:end-4) 'data\'];
else
    warndlg('Something could be wrong with the path setup!')
end

% make test sound:
stim = makeTest(PDR.TEST_seed,PDR.TEST_dur,PDR.TEST_bandwidth(1),PDR.TEST_bandwidth(2),PDR.stim_Fs,PDR.TEST_ramp,PDR.TEST_base_rms);
PDR.TEST_sound = zeros(1,PDR.buf_pts);
on_delay_pts = floor((PDR.buf_pts - length(stim))/2);
PDR.TEST_on_delay_pts = on_delay_pts;
PDR.TEST_sound(on_delay_pts+1:on_delay_pts+length(stim)) = stim; % place stimulus in buffer mid-section
PDR.TEST_sound = (0.999)*PDR.TEST_sound ./ (max(abs(PDR.TEST_sound)));
PDR.TEST_sound = PDR.TEST_sound;
clear stim;

function quit = setupTrialSeq()
global PDR session

setup=1;
while setup
	% Setup randomized scale and loc sequences:
	scales = PDR.TEST_scales;
	
	mu = mean(scales);
	scales = sort(scales);
	
	% lowest test scales:
	lows = scales(1:round(length(scales)*.333)); % lower third
	scale_seq = zeros(1,PDR.ntrials+1);
	start = PDR.npretrials+1;
	stop = PDR.ntrials+1 - mod(PDR.ntrials-PDR.npretrials+1,PDR.TEST_trial_freq);
	nm = floor((stop - start)/PDR.TEST_trial_freq)+1;
	tmp=ceil(length(lows)*rand(1,nm));
	scale_seq(start:PDR.TEST_trial_freq:stop) = lows(tmp);
	
	% middle scales:
	mids = scales(round(length(scales)*.333)+1:round(length(scales)*.667)); % middle third
	% highest scales:
	highs = scales(round(length(scales)*.667)+1:end);
	% combine mids and highs
	mid_hi = [mids highs];
	
	% randomized sequence with this pattern:
	% (hi or med) low (hi or med) low (hi or med) ...
	FRQ = 2*PDR.TEST_trial_freq;
	start = PDR.npretrials+1;
	stop = PDR.ntrials+1 - mod(PDR.ntrials-PDR.npretrials+1,FRQ);
	nm = floor((stop - start)/FRQ)+1;
	tmp = ceil(length(mid_hi)*rand(1,nm));
	scale_seq(start:FRQ:stop) = mid_hi(tmp);
	
	PDR.TEST_scale_sequence = scale_seq;
	
	% create a sequence of locations that is randomized:
	start = PDR.npretrials+1;
	stop = PDR.ntrials+1 - mod(PDR.ntrials-PDR.npretrials+1,PDR.TEST_trial_freq);
	nm = floor((stop - start)/PDR.TEST_trial_freq)+1;
	
	% randomized sequence of locations:
	loc_seq = zeros(1,PDR.ntrials+1);
	for i = start:PDR.TEST_trial_freq:stop
        loc_seq(i) = round(1 + (PDR.TEST_nlocs - 1) * rand(1)); % gives a random integer indicating trial loc IDs
	end
	
	PDR.TEST_loc_sequence = loc_seq(1:PDR.ntrials+1);


    screen_size = get(0, 'ScreenSize');
	session.hTrialPlot = figure;
	set(session.hTrialPlot, 'Position', [0.1*screen_size(3) 0.1*screen_size(4) 0.7*screen_size(3) 0.7*screen_size(4)] );
	hold on;
	x=0;
	tmp = find(PDR.TEST_loc_sequence~=0);
    if length(tmp)<1
        warndlg('No trials with this config!','No Trials');
        quit = 1;
        close(session.hTrialPlot);
        return;
    end
    cnt=0;
    scale_list=sort(PDR.TEST_scales);
	for i=1:length(PDR.TEST_loc_sequence)
        x=x+1;
        if PDR.TEST_loc_sequence(i) ~=0
            cnt=cnt+1;
            loc = mod(PDR.TEST_loc_sequence(i),5);
            scale = PDR.TEST_scale_sequence(i);
            %r= scale;
            sz(i)=ceil(log10(scale)+1)*2+4;
            colr = find(scale_list==PDR.TEST_scale_sequence(i))/length(scale_list);
            plot(x,PDR.TEST_loc_sequence(i),...
                    'Marker','*',...
                    'MarkerFaceColor','none',...
                    'MarkerEdgeColor',[colr 0 1-colr],...
                    'MarkerSize',sz(i));
        end
     
	end
    
    PDR.n_test_trials = cnt;
	
	xlabel('Trial #');
	ylabel('Location (el,az)');
	y_tick_labels{1,PDR.TEST_nlocs} = [];
	for i=1:PDR.TEST_nlocs
        y_tick_labels{i} = ['(' num2str(PDR.TEST_locs(i,1)) ',' num2str(PDR.TEST_locs(i,2)) ')'];
	end
	set(gca,'YTick',1:PDR.TEST_nlocs,'YTickLabel',y_tick_labels);
	axis([0 PDR.ntrials+1 0 PDR.TEST_nlocs+1]);
	whitebg(gcf,'k');
	warning off MATLAB:QUESTDLG:stringMismatch;
	test=questdlg(['# Tests = ' num2str(PDR.n_test_trials) ' ... Acceptable trial sequence?'],'Trial Sequence Confirmation','YES','NO','QUIT','NO');
    
    if strcmp(test,'YES')
        setup = 0;
    end
    if strcmp(test,'QUIT')
        quit = 1;
        close;
        return;
    end
    if strcmp(test,'NO')
        close;
    end
end
quit = 0;
return;

function MTLreadHDR()
global PDR
% Adapted from MTLRH
% MTLreadHDR reads the header of a SIGMAP and MTLSIG file.
%
% filetype: 3 -> frq-domain
%           2 -> time-domain
% info_blocks: number of info blocks
% n_channels: number of channels in file
% n_lines: number of lines per channel
% sampling_rate: sampling rate in kHz
% first_line: value of first line ( in ms for time domain
%                                   in kHz for frq domain )
% last_line: value of last line
% n_directions: number of directions in file
% comment1: string with 80 characters
% comment2: string with 200 characters
%
% See also: MTLWH, MTLRCH, MTLWCH

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995

%open file for reading
%disp(HRTFfile);
HRTFfile = [PDR.HRTF_directory PDR.HRTF_fname];
[fid,message] = fopen(HRTFfile,'r','ieee-le');
if (fid==-1)         				% file can't be opened
    disp(message);
    filetype=-1;
    return;
end;


fseek(fid,0,'bof');

%test if data are stored in PC-format (EOF=26) or Workstation format(EOF=XX)

fseek(fid,79,'bof');
eof = fread(fid,1,'char');
eof_int = abs(eof);
if (eof_int ~= 26);
    fclose(fid);
    error('chosen file isn''t readable');
    filetype=-1;
    return;
end;

%check filetype

fseek(fid,80,'bof');
filetype = fread(fid,1,'short');
if ( (filetype<2)|(filetype>11) )
    fclose(fid);
    error('choosen file includes unknown filetype');
end;

status=fseek(fid,82,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek info_blocks'); filetype=-1;return, end;
PDR.HRTF_info_blocks = fread(fid,1,'short');

status=fseek(fid,84,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek n_channels'); filetype=-1;return, end;
PDR.HRTF_nchannels = fread(fid,1,'short');

status=fseek(fid,128,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek n_lines'); filetype=-1;return, end;
PDR.HRTF_nlines = fread(fid,1,'short');


%if the filetype is a linear signal (hrir or hrtf)
if (filetype<7)
    status=fseek(fid,132,'bof');
    if (status<0) error('ERROR: MTLreadHDR, fseek sampling_rate'); filetype=-1;return, end;
    PDR.HRTF_sampling_rate = fread(fid,1,'float32');
    
    status= fseek(fid,136,'bof');
    if (status<0) error('ERROR: MTLreadHDR, fseek first_line'); filetype=-1;return, end;
    PDR.HRTF_first_line = fread(fid,1,'float32');
    
    status=fseek(fid,140,'bof');
    if (status<0) error('ERROR: MTLreadHDR, fseek last_lines'); filetype=-1;return, end;
    PDR.HRTF_last_line = fread(fid,1,'float32');
end;


% if the file is sigmap format it does not contain 'n_directions' field
% assume that it has one direction

status=fseek(fid,220,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek n_directions'); filetype=-1;return, end;
PDR.HRTF_ndirections = fread(fid,1,'short');
if (PDR.HRTF_ndirections==0) PDR.HRTF_ndirections=1; end;

status=fseek(fid,0,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek comment1'); filetype=-1;return, end;
comment1 = setstr(fread(fid,80,'uchar')) ;
PDR.HRTF_comment1=comment1';

status=fseek(fid,222,'bof');
if (status<0) error('ERROR: MTLreadHDR, fseek comment2'); filetype=-1;return, end;
comment2 = setstr(fread(fid,200,'uchar'));
PDR.HRTF_comment2=comment2';

PDR.HRTF_filetype = filetype;
fclose(fid);
return;

function MTLreadDIR()
global PDR
% name changed from MTLRDIR to MTLreadDIR

% direction_matrix=mtlrdir(filename,n_directions)
%
% MTLreadDIR reads direction matrix of a mtlsig file
%
% direction_matrix: matrix with n-directions columns
%                   first line: elevation
%                   second line: azimuth
%
% filename: name of the mtlsig file
%
% optional arguments:
%
% n_directions: number of directions
%
%
% See also: MTLWH, MTLRCH, MTLWCH

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995
% Date: 21.4.1995



HRTFfile = [PDR.HRTF_directory PDR.HRTF_fname];

[fid,message] = fopen(HRTFfile,'r','ieee-le');
if fid==-1          % file can't be opened
    disp(message);
    direction_matrix=-1;
    return;
end;

if (nargin < 2)
    status=fseek(fid,220,'bof');
    if (status<0) error('ERROR: mtlrdir, fseek n_directions'); direction_matrix=-1;return, end;
    PDR.HRTF_ndirections = fread(fid,1,'short');
    if (PDR.HRTF_ndirections<=0) direction_matrix=-1;return, end;
end;


status=fseek(fid,512,'bof');
if (status<0) error('ERROR: mtlrdir, fseek n_directions'); direction_matrix=-1;return, end;

direction_matrix=fread(fid,2*PDR.HRTF_ndirections,'float32');

status=fclose(fid);
if (status<0) error('ERROR: mtlrdir, fseek n_directions'); direction_matrix=-1;return, end;

direction_matrix=reshape(direction_matrix,2,PDR.HRTF_ndirections);

PDR.HRTF_dir_matrix = direction_matrix;


return;

function [dir1, dir2] = sphere2double(p1, p2)
% renamed sphere2double but functionality is the same.
%
% Original file: SPH2DBL.m	
% Description: 
% Converts spherical coordinates into double polar coordinates
%		[dir1, dir2] = sph2dbl(p1, p2):  (in degrees)
%		input arguments: saz, sel
%		output arguments: daz, del
%	or:	[dir1] = sph2dbl(p1):  (in degrees)
%		input argument: direction matrix (2 x n_dirs; el then az)
%		output argument: direction matrix (2 x n_dirs; del then daz)

if (nargin == 1)
    el = p1(1,:).*pi/180;
    az = p1(2,:).*pi/180;
    %  disp('SPH2DBL: calculating double polar direction matrix');
elseif (nargin ==2)
    el = p1*pi/180;
    az = p2*pi/180;
    %  disp('SPHDBL: calculating double polar azimuth and elevation');
else
    error('SPH2DBL: incorrect number of input arguments');
    return;
end;

% convert to cartesian coordinates:
x = cos(el) .* sin(az);
y = sin(el);
z = cos(el) .* cos(az);

%convert to double polar coordinates:
del = (el);
daz = (asin(x));

[value,index] = find(daz > pi/2);
daz(index) = pi - daz(index);

[value,index] = find(daz < -pi/2);
daz(index) = -pi - daz(index);

if (nargin == 1)
    dir1 = zeros(size(p1));
    dir1(1,:) = round(del*180/pi);
    dir1(2,:) = round(daz*180/pi);
else
    dir1 = round(del*180/pi);
    dir2 = round(daz*180/pi);
end;

%disp (['el : ' num2str(p1) ' del : ' num2str(dir1)]);
%disp (['az : ' num2str(p2) ' daz : ' num2str(dir2)]);

function [channel] = MTLreadCH(index)
global PDR
% adapted from MTLRCH"

% MTLRCH reads one channel from a file (type 2 or 3)
%
% [channel]=MTLreadCH(index);
%
% var_string = name of output variable
%
% See also: MTLWCH, MTLWH, MTLreadHDR

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995

% if no index available read first channel

if (nargin == 0) index=1; end

file_type = PDR.HRTF_filetype;
HRTFfile = [PDR.HRTF_directory PDR.HRTF_fname];

iscomplex=0;

if (file_type == 2) iscomplex=0; end
if (file_type == 3) iscomplex=1; end

[fid,message] = fopen(HRTFfile,'r','ieee-le');
if fid==-1          % file can't be opened
    disp(message);
    channel=-1;
    return
end

pos=PDR.HRTF_info_blocks*256+PDR.HRTF_nlines*4*(index-1)*(iscomplex+1);

% fseek to channel
status=fseek(fid,pos,'bof');
if (status<0) error('ERROR: MTLreadCH,  fseek to channel not succesful' );channel=-1; return; end

%fread of channel
ichannel=fread(fid,PDR.HRTF_nlines*(iscomplex+1),'float32');
if (ichannel<0) error('ERROR: MTLreadCH, reading of channel not succesful' );channel=-1; end

% creating a complex vector if necessary
if (iscomplex==1)
    channel(1:n_lines)=ichannel(1:2:PDR.HRTF_nlines*2)+i*ichannel(2:2:PDR.HRTF_nlines*2);
else
    channel=ichannel;
end

channel = channel(:);

status = fclose(fid);
if (status < 0) error('ERROR: MTLreadCH, close file');channel=-1; return; end
var_string = channel; % sets global variable with name = var_string to value of channel

function makeGammaFIR()
global PDR

Fs = PDR.stim_Fs;
cF = PDR.ADAPT_cF; % center frequency (Hz)
params = PDR.ADAPT_species; % species

% Change the EarQ and ERB if you wish to use a different ERB scale.
switch params
    case 'human'
        EarQ = 9.26449;				% Glasberg and Moore Parameters
        minBW = 24.7;
        order = 1;
        ERB = ((cF/EarQ).^order + minBW^order).^(1/order);
    case 'owl'                      % Koppl parameters
        EarQ = 1.9 * 0.074 * (cF .^ .504);
        ERB = cF ./ EarQ;
end

cF = cF(:);
T=1/Fs;
B=1.019*2*pi*ERB;

A0 = T;
A2 = 0;
B0 = 1;
B1 = -2*cos(2*cF*pi*T)./exp(B*T);
B2 = exp(-2*B*T);

A11 = -(2*T*cos(2*cF*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cF*pi*T)./ ...
    exp(B*T))/2;
A12 = -(2*T*cos(2*cF*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cF*pi*T)./ ...
    exp(B*T))/2;
A13 = -(2*T*cos(2*cF*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cF*pi*T)./ ...
    exp(B*T))/2;
A14 = -(2*T*cos(2*cF*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cF*pi*T)./ ...
    exp(B*T))/2;

gain = abs((-2*exp(4*i*cF*pi*T)*T + ...
    2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
    (cos(2*cF*pi*T) - sqrt(3 - 2^(3/2))* ...
    sin(2*cF*pi*T))) .* ...
    (-2*exp(4*i*cF*pi*T)*T + ...
    2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
    (cos(2*cF*pi*T) + sqrt(3 - 2^(3/2)) * ...
    sin(2*cF*pi*T))).* ...
    (-2*exp(4*i*cF*pi*T)*T + ...
    2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
    (cos(2*cF*pi*T) - ...
    sqrt(3 + 2^(3/2))*sin(2*cF*pi*T))) .* ...
    (-2*exp(4*i*cF*pi*T)*T + 2*exp(-(B*T) + 2*i*cF*pi*T).*T.* ...
    (cos(2*cF*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cF*pi*T))) ./ ...
    (-2 ./ exp(2*B*T) - 2*exp(4*i*cF*pi*T) +  ...
    2*(1 + exp(4*i*cF*pi*T))./exp(B*T)).^4);

allfilts = ones(length(cF),1);
fcoefs = [A0*allfilts A11 A12 A13 A14 A2*allfilts B0*allfilts B1 B2 gain];

% 512 points:
impulse = [1 zeros(1, 511)];
y1 = filter([A0/gain A11/gain A2/gain], [B0 B1 B2], impulse);
y2 = filter([A0 A12 A2], [B0 B1 B2], y1);
y3 = filter([A0 A13 A2], [B0 B1 B2], y2);
response = filter([A0 A14 A2], [B0 B1 B2], y3);

% scale response for unity gain:
response = response./sum(abs(response));
PDR.ADAPT_coefs = response;

function err = setupCalib(flag_adapt)
% setup for calibration session
global PDR
err = 0;

PDR = struct(...                        % MAIN PARAMETERS:
    'calib',1,...                       % flag for calibrations
    'flag_adapt',flag_adapt,...                  % adaptor or test probe calibration?
    'virtual',1,...                     % Freefield: virtual = 0, Headphones: virtual = 1
    'record',0,...                      % flag for recording pdr trace
    'ntrials',[],...                   % no. of trials in a session
    'buf_pts',16384*10,...              % 32768*10 for 6711 ms buffers (with Fs = 48828 Hz)... permits 5 second calibration sounds
    'buf_dur',[],...                    % buffer duration in ms
    'isi_buf',1,...                     % number of empty buffers between trials (6 buffers gives an ISI of 3.73s using 100ms sounds, for example)
    'decimationfactor',4,...            % decimation factor for data collection
    'stim_Fs',30000,...                 % sampling rate for sound production
    'npts_totalplay',[],...            % total no. pts in the session
    'len_session',[], ...               % length of session (in minutes)
    'starttime',[], ...                 % session start time
    'stoptime', [], ...                 % session stop time
    'code_path', 'c:\andrew\pdr\code\',...              % path to code
    'data_path', 'c:\andrew\pdr\calib_data\',...        % flag indicates that AD recording (pupillometer) will be needed
    'base_atten',30,...                  % attenuation value to send to TDT (combining scale_val and atten gives a certain SPL in dB -- This needs to be calibrated!!!)
    'filename',[],...                   % file name for storing session data
    ...
    ...
    'CALIB_ear',[],...                  % left or right earphone
    'CALIB_nreps',2,...                 % repeat each stimulus x times
    'CALIB_soundtype',[],...            % test or adaptor
    'CALIB_comments',[],...             % extra info about sound
    'CALIB_sound',[],...                % calibration sound (should be a 5 second version of sounds used in the experiment, with the same rms level)
    'CALIB_sound_dur',5000,...          % 5 second sounds for calibrations
    'CALIB_bandwidth',[],...
    'CALIB_scales',30000:-2000:6000,... % scales to use for calibration
    'CALIB_scale_sequence',[],...       % sequence of scales
    'CALIB_locations',[],...            % a 2xnlocs vector of location coordinates (Az,El) ... Azimuths in top row, Elevations in bottom row
    'CALIB_coord_sequence',[],...       % sequence of location coordinates (Az,El) ... Should be selected from CALIB_locations
    'CALIB_locID_sequence',[],...         % sequence of location IDs (for HRIR selection)
    'CALIB_nlocs',[],...                 % calibrate for x locations
    'CALIB_nscales',[],...              % number of scales used
    'CALIB_ramp',5);                  % use 5ms ramp
    
if PDR.virtual
    PDR.HRTF_directory = 'C:\andrew\pdr\HRTF_lib\';
    PDR.HRTF_fname = '930AD_ABLequal.eq';
    PDR.HRTF_nchannels = 2;
    PDR.HRTF_nlines = 255;
    PDR.HRTF_sampling_rate = 30000;
end


% GENERAL SETUP (Trial Sequence, Etc.):

ear = menu('Left or right ear?','LEFT','RIGHT');

if ear==1
    PDR.CALIB_ear = 'LEFT';
else
    PDR.CALIB_ear = 'RIGHT';
end

Hwarn=warndlg(['Make sure only the ' PDR.CALIB_ear ' headphone is connected!']);
uiwait(Hwarn);

PDR.CALIB_nscales = length(PDR.CALIB_scales);

% Reminder: the HRTF files are accessed using (el,az) coordinates with
% line one = elevation and line two = azimuth!

PDR.CALIB_locations(2,:) = -80:20:80; % azimuths (2nd row for azimuths!)
PDR.CALIB_nlocs = length(PDR.CALIB_locations(2,:));
PDR.CALIB_locations(1,:) = zeros(1,PDR.CALIB_nlocs); % zero elevation throughout (elevation is first row!)
PDR.CALIB_nscales = length(PDR.CALIB_scales);
PDR.ntrials = PDR.CALIB_nlocs*PDR.CALIB_nscales*PDR.CALIB_nreps;

% setup location/scale sequences:
loc_order = randperm(PDR.CALIB_nlocs);

PDR.CALIB_scale_sequence = [];
PDR.CALIB_coord_sequence = NaN*ones(2,PDR.ntrials);
PDR.CALIB_locID_sequence = NaN*ones(1,PDR.ntrials);
trialsperloc = PDR.CALIB_nscales*PDR.CALIB_nreps;
start = 1;
stop = trialsperloc;
for i = 1:PDR.CALIB_nlocs
    PDR.CALIB_locID_sequence(1,start:stop) = loc_order(i)*ones(1,trialsperloc);
    PDR.CALIB_coord_sequence(1,start:stop) = PDR.CALIB_locations(1,loc_order(i))*ones(1,trialsperloc);
    PDR.CALIB_coord_sequence(2,start:stop) = PDR.CALIB_locations(2,loc_order(i))*ones(1,trialsperloc);
    scale_order = PDR.CALIB_scales(randperm(length(PDR.CALIB_scales)));
    for j=1:PDR.CALIB_nscales
        PDR.CALIB_scale_sequence = [PDR.CALIB_scale_sequence scale_order(j)*ones(1,PDR.CALIB_nreps)];
    end
    start = stop + 1;
    stop = stop + trialsperloc;
end


if flag_adapt

    PDR.CALIB_soundtype = 'adaptor';
    PDR.CALIB_comments = 'gammatone (6 kHz center frequency)';

    PDR.CALIB_sound_dur = PDR.buf_pts/PDR.stim_Fs;
    % Make FIR coefficients for Adaptor:
    PDR.ADAPT_cF = 6000;
    PDR.ADAPT_species = 'owl';
    spatialPDR_switchyard('makeGammaFIR');

end

% setup for calibrating test sounds
if ~flag_adapt
    
    PDR.CALIB_soundtype = 'test';
    PDR.CALIB_comments = 'octave band test probes (4-8kHz)';
    PDR.CALIB_seed = 47;
    PDR.CALIB_base_rms = 0.0350;
    PDR.CALIB_bandwidth = [4000 8000];

    % Create a 5 second version of the test sound:
    stim = makeTest(PDR.CALIB_seed,PDR.CALIB_sound_dur,PDR.CALIB_bandwidth(1),PDR.CALIB_bandwidth(2),PDR.stim_Fs,PDR.CALIB_ramp,PDR.CALIB_base_rms);
    PDR.CALIB_sound = zeros(1,PDR.buf_pts);
    on_delay_pts = floor((PDR.buf_pts - length(stim))/2);
    PDR.CALIB_on_delay_pts = on_delay_pts;
    PDR.CALIB_sound(on_delay_pts+1:on_delay_pts+length(stim)) = stim; % place stimulus in buffer mid-section
    PDR.CALIB_sound = (0.999)*PDR.CALIB_sound ./ (max(abs(PDR.CALIB_sound)));
end

PDR = orderfields(PDR); % order fields by ASCII dictionary order

% buffer duration (ms):
PDR.buf_dur = 1000*PDR.buf_pts/PDR.stim_Fs;

q=clock;
y=num2str(q(1));y=y(3:4);
m=num2str(q(2));if size(m,2)<2;m=['0' m];end
d=num2str(q(3));if size(d,2)<2;d=['0' d];end

cd(PDR.code_path);
tmp = pwd;
if tmp(end-3:end) == 'code'
    PDR.data_path = [tmp(1:end-4) 'calib_data\'];
else
    warndlg('Something could be wrong with the path setup!')
end

cd(PDR.data_path);
PDR.filename = [y m d 'calib_' PDR.CALIB_soundtype '_' PDR.CALIB_ear '_a'];
count = double('a'+0);
while exist ([PDR.filename '.mat'],'file');
    count = count + 1;
    if count > 122
        disp('There are already several files with similar names!');
        PDR.filename = input('Enter a unique filename for this session: ', 's');
        break;
    else
        PDR.filename(end) = char(count);
    end
end
PDR.npts_totalplay = PDR.ntrials*(PDR.isi_buf+1)*PDR.buf_pts; % Calculate length of session!
PDR.len_session = (1/60)*(PDR.npts_totalplay/PDR.stim_Fs); % length of session in minutes

screen_size = get(0, 'ScreenSize');
session.hTrialPlot = figure;
set(session.hTrialPlot, 'Position', [0.02*screen_size(3) 0.25*screen_size(4) 0.95*screen_size(3) 0.6*screen_size(4)] );
hold on;
x=0;

tmp = sort(PDR.CALIB_locations(2,:));
axis([0 PDR.ntrials+1 tmp(1)-10 tmp(end)+10]);
for i=1:length(PDR.CALIB_locID_sequence)
    x=x+1;
    r = PDR.CALIB_scale_sequence(i)/32768;
    plot(x,PDR.CALIB_coord_sequence(2,i),...
                'Marker','*',...
                'MarkerFaceColor','none',...
                'MarkerEdgeColor',[r 0 1-r],...
                'MarkerSize',4+round(15*r));
end

title('Calibration Trial Sequence');
xlabel('Trial #');
ylabel('Location (el,az)');
y_tick_labels{1,PDR.CALIB_nlocs} = [];
for i=1:PDR.CALIB_nlocs
    y_tick_labels{i} = ['(' num2str(PDR.CALIB_locations(1,i)) ',' num2str(PDR.CALIB_locations(2,i)) ')'];
end
set(gca,'YTick',PDR.CALIB_locations(2,:),'YTickLabel',y_tick_labels,'XTick',1:10:PDR.ntrials);
whitebg(gcf,'k');
warning off MATLAB:QUESTDLG:stringMismatch;
    
function stim = makeTest(seedval,dur,minband,maxband,samprate,ramp,base_rms)
% makes test sounds using specific parameters
global PDR

rampPts = ceil((ramp/1000)*samprate);
[temp,stim] = seeded_whnoise(minband,maxband,samprate,dur/1000,seedval);
stim = stim';  %Stimuli are scaled between 1 and -1 in whnoise and not rescaled here; 
%MUST BE A ROW VECTOR!!!!
npts = length(stim);
% put envelope on it
env = [0:1/(rampPts-1):1 ones(1,npts - 2*rampPts) 1:-1/(rampPts-1):0];
stim = stim .* env;
% scale to desired rms:
rms = mean(stim.^2)/sqrt(2);
ratio = rms/base_rms;
stim = stim/sqrt(ratio);

function [X,Xtime] = seeded_whnoise(minfreq,maxfreq,Fs,duration,seedval,phi)
%A script to produce noise with a flat spectrum between specified frequencies.
%Usage: [freq domain, time domain] = whnoise(min frequency, max frequency,Fs,durationor)
%durationor is used to determine length; len = duration * Fs;
%Borrowed from Kip Keller w_rnd.m

if (nargin < 3) Fs = 30000; end
if (nargin < 4) duration = .1; end
if (nargin < 5) randn('seed',47);
else
    randn('seed',seedval);
end
if (nargin < 6) phase = (rand(duration*Fs,1) - 0.5) * (2*pi);
else
phase = phi;
end

len = duration*Fs;

minfreq = round(((minfreq+1)/Fs) * len);
maxfreq = round(((maxfreq+1)/Fs) * len);
range = maxfreq-minfreq+1;

% mag spectrum = 1 between set frequencies:
mag = zeros(len,1);
mag(minfreq:maxfreq) = ones(range,1);

% random phase spectrum between set frequencies:
			%-0.5 moves phases from 0to1 range to -0.5to0.5 range, then have to multiply by 2pi instead of pi
% combine phase and magnitude:
X = mag .* ( (cos(phase)) + (i .* sin(phase)) );

% convert to time domain:
Xtime = real(ifft(X));
Xtime=Xtime/max(abs(Xtime));						%have to normalize so that max and mins are between -1 and 1 otherwise will saturate = bad


while 0
valnoise = Xtime;
valnoise = abs(hilbert(valnoise));
valnoise = valnoise*10000;
valnoise = valnoise-mean(valnoise);
figure
% Use next highest power of 2 greater than or equal to length(x) to calculate FFT.
nfft= 2^(nextpow2(length(valnoise))); 
% Take fft, padding with zeros so that length(fftx) is equal to nfft 
fftx = fft(valnoise,nfft); 
% Calculate the numberof unique points
NumUniquePts = ceil((nfft+1)/2); 
% FFT is symmetric, throw away second half 
fftx = fftx(1:NumUniquePts); 
% Take the magnitude of fft of x and scale the fft so that it is not a function of the length of x
mx = abs(fftx)/length(valnoise); 
% Take the square of the magnitude of fft of x. 
mx = mx.^2; 
% Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
% The DC component and Nyquist component, if it exists, are unique and should not be multiplied by 2.
if rem(nfft, 2) % odd nfft excludes Nyquist point
  mx(2:end) = mx(2:end)*2;
else
  mx(2:end -1) = mx(2:end -1)*2;
end
% This is an evenly spaced frequency vector with NumUniquePts points. 
f = (0:NumUniquePts-1)*Fs/nfft; 
% normalize mx
mx = mx/max(mx);

mx = 20*log10(mx);
% Generate the plot, title and labels. 
plot(f,mx);
axis([0 20000 -100 100])
title('Power Spectrum'); 
xlabel('Frequency (Hz)');
ylabel('Power');
end

return;

function quit_flag = run_calibration
% plays 5 second sounds through the TDT for calibrations
global PDR HRTF

fig0 = gcf;
title('TRIAL SEQUENCE (CALIBRATIONS)');
quit_flag = 0; % quit flag

lengthOFtrials=(round(PDR.npts_totalplay/PDR.stim_Fs/60*10))/10;  %should be in approximate minutes
b=clock; PDR.starttime(1:2)=b(4:5);

calibInfo = struct(...
            'latten',PDR.base_atten, ...
            'ratten',PDR.base_atten, ...
            'scale',[], ...
            'buf_pts',16384*10,... % 5 seconds
            'flag_adapt',PDR.flag_adapt,...
            'nreps',PDR.CALIB_nreps,...
            'HRTF_nlines',PDR.HRTF_nlines);

%write header information to file... saving global variables
PDR = orderfields(PDR); % order fields by ASCII dictionary order
save ([PDR.data_path PDR.filename '.mat'],'calibInfo','PDR','HRTF');

% initialize all SystemII hardware

if(S232('S2init',0,'INIT_PRIMARY',1000)==0)
   warndlg('FAILED to initialize a primary process','warning');
   return;
else
   disp('Initialized primary process');
end

if(S232('APlock',100, 0)==0)
   warndlg('FAILED to get AP lock for initialization','warning');
   return;
else
   disp('AP lock in initiation');
end
if(S232('XBlock',100, 0)==0)
   warndlg('FAILED to get XB lock for initialization','warning');
   S232('APunlock',0);
   warndlg('AP unlock in initiation','warning');
   return;
else
   disp('XB lock in initiation');
end

err = S232('getS2err');
switch err
case 0
case 1: warndlg('APOS error in initiation','warning'); return;
case 2: warndlg('XBUS error in initiation','warning'); return;
end


SRATE =  (1/PDR.stim_Fs)*10^6; % sampling rate for TDT

% PREPARE PD1 FOR CONVERSION:
s232('PD1clear',1);
s232('PD1fixbug',1);
s232('PD1nstrms',1, 2, 0);
s232('PD1srate',1,33.3);
s232('PD1npts',1, PDR.buf_pts);
s232('PD1clrsched',1);
s232('PD1addsimp',1, s232('IREG',0), s232('DAC',0));
s232('PD1specIB',1, s232('IB',0), s232('IREG',0));
s232('PD1addsimp',1, s232('IREG',1), s232('DAC',1));
s232('PD1specIB',1, s232('IB',1), s232('IREG',1));
s232('PD1setIO',1,.01,9.99,.01,9.99);

% SET PA4 ATTENUATION:
s232('PA4atten',1,calibInfo.latten);
s232('PA4atten',2,calibInfo.ratten);

% COUNTDOWN TO START:

in = inputdlg('Enter a countdown (in seconds) for each set of trials to begin (45 seconds suggested): ');
if isempty(in)
    quit_flag = 1;
    return;
end
start_cnt = str2num(in{1});
clear in

in = inputdlg('Enter how many trials you want to do at a time (20 suggested): ');
if isempty(in)
    quit_flag = 1;
    return;
end
step = str2num(in{1});
clear in

% for data input:
PDR.CALIB_SPLs = nan*ones(1,PDR.ntrials);
prompt = cell(1,step);
def = cell(1,step);

% MAIN LOOP:
% check for what type of calibration this is:
if ~calibInfo.flag_adapt
    ADAPT_coefs = zeros(1,255); % empty set of coefficients if calibrating test sounds
else
    ADAPT_coefs = PDR.ADAPT_coefs;
    PDR.CALIB_sound = zeros(1,PDR.buf_pts); % adaptor sound created using C code (FIR filtering)
end


i=0;
while i < PDR.ntrials

    % give user some warning and time to prepare:
    if i+step <= PDR.ntrials
        txt = ['Click okay to begin calibration trials ' num2str(i+1) ' through ' num2str(i+step)];
    else
        step = PDR.ntrials-i; % use a curtailed step size
        txt = ['Click okay to begin calibration trials ' num2str(i+1) ' through ' num2str(PDR.ntrials)];
    end
    
    uiwait(msgbox(txt,'Calibrations:','modal'));
    fig1 = figure; axis off;
    countdown = start_cnt;
    for ii = countdown:-1:0
        if i+step <= PDR.ntrials
            txt = ['Next ' num2str(step) ' trials begin in ' num2str(countdown) ' seconds!'];
        end
        countdown = countdown - 1;
	    t=text(.5,.5,txt,'HorizontalAlignment','center','FontSize',18,'EdgeColor',[1 0 0]);
        pause(1);
        delete(t);
    end
    close(fig1);
    
    box_color = [rand rand rand];
    clear prompt trials def
    for j=1:step
        
        % trial counter:
        i=i+1;
              
        % prepare TDT and select parameters:
        s232('trash');
        s232('dropall');
        calibInfo.scale = PDR.CALIB_scale_sequence(i);
        hrtfL = HRTF.CALIB_Left(:,PDR.CALIB_locID_sequence(i))';
        hrtfR = HRTF.CALIB_Right(:,PDR.CALIB_locID_sequence(i))';
        
        % Plot trial sequence info:
        figure(fig0); hold on;
        plot(i,PDR.CALIB_coord_sequence(2,i),'s','MarkerSize',4,'MarkerEdgeColor',box_color);
        drawnow;
        
        % run c code:
        spatialPDR_calib(calibInfo, single(hrtfL), single(hrtfR), single(ADAPT_coefs), single(PDR.CALIB_sound));
        pause(2) % 2 second pause between trials
        
        % for data input:
        prompt{j} = ['Trial #' num2str(i)];
        trials(j) = i;
        def{j} = '';
        
    end
    
    
    ok = 0;
    in=questdlg('Repeat those trials again?','REPEAT','YES','NO','NO');
    if strcmp(in,'YES')
        ok = 1;
        i = i-step;
    end
    clear in
    
    % for data input:
    while ~ok
        in = inputdlg(prompt, ['Enter SPLs for trials ' num2str(trials(1)) ' through ' num2str(trials(step)) ' : '], 1, def);
        if isempty(in)
            quit_flag = 1;
            return;
        end
        
        list = [];
        for k=1:step     
            PDR.CALIB_SPLs(trials(k)) = str2num(in{k});
            list = [list ',' in{k}];
        end
        list = list(2:end);
        txt = ['Are you sure these are the corresponding SPLs? (SPLs = ' list ')'];
        check = questdlg(txt,'Check Input: ','No');
        if strcmp(check,'Yes')
            ok = 1;
            break;
        end
        
        clear in
    end
    
end

cd(PDR.data_path);
%write header information to file... saving global variables
save ([PDR.data_path PDR.filename '.mat'],'PDR','HRTF');

try
	% plot SPLs for each location and fit a line:
	PDR.CALIB_regression_fit = cell(1,PDR.CALIB_nlocs);
	PDR.CALIB_rsquare_fit = cell(1,PDR.CALIB_nlocs);
	close all;
    screen_size = get(0, 'ScreenSize');
    fig2=figure;
    set(fig2, 'Position', [0.02*screen_size(3) 0.05*screen_size(4) 0.95*screen_size(3) 0.8*screen_size(4)] );
	dim=ceil(sqrt(PDR.CALIB_nlocs));
    a = 0.05; %alpha
    h=waitbar(0,'Please wait for plots...');
	for ii = 1:PDR.CALIB_nlocs
        subplot(dim,dim,ii);
        set(gca,'Color',[.3 .3 .3]);
        hold on;
        xlabel('log10(scale)','FontSize',8);
        ylabel('Intensity (dB SPL)','FontSize',8);
        trials = find(PDR.CALIB_locID_sequence == ii);
        xes = log10(PDR.CALIB_scale_sequence(trials));
        yes = PDR.CALIB_SPLs(trials);
        xrange=0:0.01:5;
        
        % regression with stats
        [PDR.CALIB_rsquare_fit{ii}, PDR.CALIB_regression_fit{ii}] = regress_stats(xes,yes,a,xrange);
        hold on;
        title_text = ['Loc (el,az): (' num2str(PDR.CALIB_locations(1,ii)) ',' num2str(PDR.CALIB_locations(2,ii)) ')'];
        title_text = [title_text ' (R^{2} = ' num2str(PDR.CALIB_rsquare_fit{ii}) ')'];
        title(title_text,'FontSize',8);
        hold off;
        waitbar(ii/PDR.CALIB_nlocs,h);
    end
	close(h)
	set(gcf,'InvertHardcopy','off');
	saveas(fig2,PDR.filename,'fig');
catch err
    msgbox('Could not create/save regression fits!!!','ERROR','error');
end

% RESET TDT:
s232('PD1clear', 1);
s232('XBunlock', 0);
s232('S2close');


b=clock;
PDR.stoptime(1:2)=b(4:5);
%write header information to file... saving global variables (AGAIN)
save ([PDR.data_path PDR.filename '.mat'],'PDR','HRTF');

cd ..\code
return;

function [RSQ, COEFFS] = regress_stats(xes,yes,alfa,xrange)
% plot data
hold on;
n=length(yes);
df=length(yes)-2;
for j = 1:n
    plot(xes,yes,'m*');
end

% determine fit coefficients and CI
[p, S] = polyfit(xes,yes,1);
[Y,DELTA]=polyconf(p,xrange,S,alfa);
m=p(1);
b=p(2);
yfit =  m.* xrange + b;
ypred = (m.*xes + b);
yresid = yes - ypred;
SSresid = sum(yresid.^2);
SStotal = (length(yes)-1) * var(yes);
RSQ = 1 - SSresid/SStotal;

% plot fit with CI
plot(xrange,Y,'c-');
plot(xrange,Y-DELTA,'r--');
plot(xrange,Y+DELTA,'r--');

% standard error calculation
SSxx = sum(xes.^2)-n*mean(xes)^2;
SSyy = sum(yes.^2)-n*mean(yes)^2;
SSxy = sum(xes.*yes)-n*mean(xes)*mean(yes);
s=sqrt((SSyy-m*SSxy)/(n-2));
SE_m = s/sqrt(SSxx); % standard error for slope
SE_b = s*sqrt((1/n)+mean(xes)^2/SSxx); % standard error for intercept

% determine t statistic
step = 0.01;
t=step;
cum=0;
while cum < (1-alfa)
    tes=-t:step:t;
    tmp=tpdf(tes,df);
    cum=sum(tmp)*step;
    if cum > 0.95
        break;
    else
        t=t+step;
    end
end

% determine coefficient CIs
COEFFS = zeros(3,3);
COEFFS = [m m-SE_m*t m+SE_m*t; b b-SE_b*t b+SE_b*t];

% set axis values
minX=min(xrange);
maxX=max(xrange);
minY=min(Y-DELTA);
maxY=max(Y+DELTA);
axis([minX maxX minY maxY]);

% plot info
eqn = ['Fit: y = mx+b = ' num2str(m) 'x + ' num2str(b)];
text(minX+0.1*(maxX-minX),minY+0.9*(maxY-minY),eqn,'FontSize',8);
mcoeff=['m = ' num2str(m) ' [' num2str(COEFFS(1,2)) ',' num2str(COEFFS(1,3)) ']'];
text(minX+0.1*(maxX-minX),minY+0.75*(maxY-minY),mcoeff,'FontSize',8);
bcoeff=['b = ' num2str(b) ' [' num2str(COEFFS(2,2)) ',' num2str(COEFFS(2,3)) ']'];
text(minX+0.1*(maxX-minX),minY+0.6*(maxY-minY),bcoeff,'FontSize',8);
hold off;