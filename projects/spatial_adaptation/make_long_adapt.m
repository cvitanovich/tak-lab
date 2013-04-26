
% make adaptor in a loop
dur_sec =60*2; % 200 minutes
P.ADAPT_seed=47; P.buf_pts=3000000; P.stim_Fs=30000; P.TEST_ramp=5;
HRTF.directory='C:\andrew\CORE\tak-lab\HRTFs\Matlab_V6\';
HRTF.fname='1073AC_eq_ABLequal_norm.mat';
P.HRTF_fname=HRTF.fname;
P.HRTF_directory=HRTF.directory;
% preallocate a sound buffer:
tmp=NaN*ones(1,P.buf_pts);

el=0; az=0;
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

%hL=figure; hR=figure;
%h=figure;
mx=-Inf;

clear tm
for j=0:dur_sec
    tic
    % get new set of pseudorandomly generated numbers and run through a FIR filter
    rand('state',adapt_state);
    new_buffer=rand(1,P.buf_pts);
    adapt_state=rand('state');
    
    [filtered_buffer, CIRC_BUFS.adaptor] = circ_fir(CIRC_BUFS.adaptor,new_buffer,P.GTONE_coefs);
    % filter adaptor with left/right HRTF coefficients
    [adapt_left, CIRC_BUFS.left] = circ_fir(CIRC_BUFS.left,filtered_buffer,HRTF.left);
    [adapt_right, CIRC_BUFS.right] = circ_fir(CIRC_BUFS.right,filtered_buffer,HRTF.right);
    %t=j:(1/P.stim_Fs):j+1;
    %t=t(1:P.buf_pts);
    %figure(hL); hold on; plot(t,adapt_left)
    %figure(hR); hold on; plot(t,adapt_right)
    
    tmp=max(max(abs(adapt_left)),max(abs(adapt_right)));
    last_mx=mx;
    mx=max(mx,tmp);
    tm(j+1)=toc;
    mean_tm=mean(tm);
    if mx>last_mx
        fprintf('new max = %d\n',mx);
    end
    fprintf('Est %d seconds remaining\n\n',((dur_sec-j)/100)*mean_tm);
%     subplot(2,1,1); plot(CIRC_BUFS.left);
%     subplot(2,1,2); plot(CIRC_BUFS.right,'r');
%     drawnow
end
disp(num2str(mx));