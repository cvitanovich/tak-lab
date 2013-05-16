function [PDR,A] = load_filter_pdr_data(PDR,A,invert_flg,rec_ch)
    A.numpts=PDR.buf_pts/(2^PDR.decimationfactor);
    A.numbufs=PDR.npts_totalplay/(2^PDR.decimationfactor)/A.numpts;  
    fid=fopen([A.pname PDR.filename '_REC' num2str(rec_ch) A.ext]); A.raw_data=fread(fid,'short'); fclose(fid);
    A.trialtype=A.raw_data(1:A.numpts:size(A.raw_data));
    A.idx2trials=find(A.trialtype==1);
    A.idx2trials = A.idx2trials(1:end); 
    PDR.ntrials = length(A.idx2trials);
    PDRdata.trace=cell(length(A.idx2trials),1);
    % remove last trial if there's a buffer mismatch
    while size(A.raw_data,1) < ((A.numpts*(A.idx2trials(end)-1)+1) + A.bufs_after*A.numpts)
        A.idx2trials=A.idx2trials(1:(end-1));
    end
    len=length(A.idx2trials);
    TRIALS = nan(len,A.bufs_total*(A.numpts-A.infopts)); % temporary matrix for parsing trials

    t=0;
    while t < len
        t=t+1;
        q=1 + A.numpts*(A.idx2trials(t)-1);
        q=q - A.numpts*A.bufs_before; % include A.bufs_before buffers before trial buffer (AC 6/27/2012)
        for k=1:A.bufs_total
            startPt = 1+((k-1)*(A.numpts-A.infopts));
            stopPt = (k*(A.numpts-A.infopts));
            if (q+1) < 1
                TRIALS(t,startPt:stopPt) = nan(A.numpts-A.infopts,1);
            else
                TRIALS(t,startPt:stopPt) = A.raw_data(q+A.infopts:q+A.numpts-1);
            end
            q=q+A.numpts;
        end
        
        % sets values in "prior" buffers equal to first value in 1st
        % buffer (only for first stimulus trace):
        tmp = isnan(PDRdata.trace{t,1});
        tmp = find(tmp == 1);
        if ~isempty(tmp)
            TRIALS(t,tmp) = TRIALS(t,tmp(end)+1).*ones(1,length(tmp));
        end
        PDRdata.numpts(t)=A.bufs_total*(A.numpts-A.infopts);
    end

    PDRdata.trace = mat2cell(TRIALS,ones(1,len));

    for n = 1:length(A.idx2trials)
        P_data (n,:) = PDRdata.trace{n}(1,:)';
    end

    P_data = P_data';
    A.ntrials = length(A.idx2trials);
    y = zeros(A.bufs_total*(A.dec_pts-A.infopts),A.ntrials)*NaN;
    A.meandata = zeros(A.bufs_total*(A.dec_pts-A.infopts),A.ntrials)*NaN;
    %filter AC noise out of trace; normalize all trial DC levels
    h0 = waitbar(0,'Hang on, filtering data...');
    % 20 Hz low pass filter:
    FilterSpec = fdesign.lowpass('N,Fc,Ap,Ast',300,20,1,60,A.dec_Fs);
    FilterObj = design(FilterSpec,'equiripple');
    b=FilterObj.Numerator;
    for trial = 1:A.ntrials
        %Fs = 1000;
        %t = 0:1/Fs:1;
        %b = ones (1,94)/94;
        y (:,trial)= filtfilt(b,1, P_data(:,trial));
        A.meandata (:, trial) = (y(:, trial))-((((y(round(A.bufs_before*(A.dec_pts-A.infopts)+A.sndStart*(A.dec_pts-A.infopts)), trial))))); % changed to zero at sound onset (A.sndStart)
        P_data (:, trial) = ( P_data(:, trial))-(((( P_data(A.bufs_before*(A.dec_pts-A.infopts)+round(A.sndStart*(A.dec_pts-A.infopts)), trial)))));
        waitbar(trial/A.ntrials,h0);
    end
    close(h0)
    if(invert_flg)
        A.meandata = -1*A.meandata; %INVERTING ONLY!
    end