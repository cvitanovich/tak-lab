% estimate distribution of amplitudes for adaptor

% 100 sec sound
P.ADAPT_seed=47; P.buf_pts=30000000; P.stim_Fs=30000; P.TEST_ramp=5;
HRTF.directory='/Users/cvitanovich/Documents/MATLAB/tak-lab/HRTFs/Matlab_V6/';
HRTF.fname='1073AC_eq_ABLequal_norm.mat';
P.HRTF_fname=HRTF.fname;
P.HRTF_directory=HRTF.directory;
% preallocate a sound buffer:
tmp=NaN*ones(1,P.buf_pts);

%% FIND HRTF COEFFIENCTS WITH LARGEST AMPLITUDE
load([P.HRTF_directory P.HRTF_fname]);
largest=-Inf;
for j=1:length(direc)
    el=direc(1,j); az=direc(2,j);
    [HRTF.left,HRTF.right] = readHRTFdotMAT(P.HRTF_directory,P.HRTF_fname,el,az);
    amp = max(sum(abs(HRTF.left)),sum(abs(HRTF.right)));
    if amp > largest
        largest=amp;
    end
end
largest

%% LOAD HRTF COEFFICIENTS
if(strcmp(P.HRTF_fname((end-3):end),'.mat'))
    [HRTF.left,HRTF.right] = readHRTFdotMAT(P.HRTF_directory,P.HRTF_fname,el,az);
else % using other format (e.g. for 930 or 929)
    %read HRTF coefficients files:
    
    HRTF = MTLreadHDR(HRTF);
    HRTF = MTLreadDIR(HRTF);
    % THIS CONVERSION IS NOT NECESSARY WITH ANY OF THE BIRDS WE USE NOW:
    % %convert coordinates to double polar, only if not using 929 or 930 or ones
    % if (~strcmp(P.HRTF_fname(1:3),'930') & ~strcmp(P.HRTF_fname(1:3),'929'))
    %     P.HRTF_dir_matrix = sphere2double(P.HRTF_dir_matrix);
    % end
    
    HRTF.left = nan*ones(HRTF.nlines,1);
    HRTF.right = HRTF.left;
    direc = HRTF.dir_matrix; % NOTE: first row = Elevation and 2nd row = Azimuth !!!
    [r c]=find(direc(1,:)==0);
    idx=find(direc(2,c)==0); % using frontal location for calibration
    HRTF.left = MTLreadCH(idx*2-1, HRTF);
    HRTF.right = MTLreadCH(idx*2, HRTF);
end


%% MAKE GAMMATONE (cF=6kHz)
P.GTONE_cF = 6000; P.GTONE_species = 'owl';
P.GTONE_coefs = makeGammaFIR(P.stim_Fs,P.GTONE_cF,P.GTONE_species);

% initial state
adapt_state = P.ADAPT_seed;
% circular buffer for continuous filtered adaptor:
CIRC_BUFS.adaptor=zeros(1,(length(P.GTONE_coefs)+P.buf_pts));
% circular buffers for HRTF filtering (left/right):
CIRC_BUFS.left=CIRC_BUFS.adaptor; CIRC_BUFS.right=CIRC_BUFS.adaptor;

% get new set of pseudorandomly generated numbers and run through a FIR filter
rand('state',adapt_state);
new_buffer=rand(1,P.buf_pts);
adapt_state=rand('state');

[filtered_buffer, CIRC_BUFS.adaptor] = circ_fir(CIRC_BUFS.adaptor,new_buffer,P.GTONE_coefs);
% filter adaptor with left/right HRTF coefficients
[adapt_left, CIRC_BUFS.left] = circ_fir(CIRC_BUFS.left,filtered_buffer,HRTF.left);
[adapt_right, CIRC_BUFS.right] = circ_fir(CIRC_BUFS.right,filtered_buffer,HRTF.right);

% make a histogram of amplitudes
nbins=100;
[f x] = hist([adapt_left adapt_right],nbins);
figure(1)
bar(x,f/sum(f));hold on