function [peak_altIR,peakHeight_altIR,peakWidth_altIR,peak_std,peakHeight_std,peakWidth_std,peakHeight_altIR_std,peakWidth_altIR_std,dir] = compare_mls_stims2(FN_HRTF,flag,dir_tested)

%function [peak_altIR,peakHeight_altIR,peakWidth_altIR,peak_std,peakHeight_std,peakWidth_std,peakHeight_altIR_std,peakWidth_altIR_std,dir] = compare_mls_stims2(FN_HRTF,flag,dir_tested);
% built 8/1/03 to compare Hanna's altIR stims with unaltered HRIRs
% version to look at 2sources (will test dir_tested to see if files exist)
% returns dir with all locations used

fftpts = 1024;
freq = 0:15000/(fftpts/2 -1):15000;

if nargin < 2
    flag = 1;
end
stimuli_dir = ['e:\kip\noisetc\temp' num2str(flag) '\'];
saveFN = ['e:\kip\temp\XXX' num2str(flag)];

if nargin < 3
    load e:\kip\matlab\mls_tdt\dir_tested
end

dir_HRTF = sph2dbl(mtlrdir(FN_HRTF));
[temp,nLocs] = size(dir_tested);

Fs = 30000;
% calc cF (center frequencies for filterbank) on 1/12th octave scale
cF = round(1000*exp(([12:40]/12)*log(2)))';
n_cF = length(cF);
% calc'd from the envelopes of Gaussian-filtered long_b10:
maxFactor = ones(size(cF)) * .00003764*cF(n_cF)+.6236;
% scaling factor relative to highest freq
Factor = maxFactor./(.00003764*cF + .6236);
fcoefs = Make_ERBFiltA(Fs,cF);
Factormat1 = repmat(Factor,1,3254);
Factormat2 = repmat(Factor,1,255);
iLocInner = 0;

for iLoc = 1:nLocs
    disp(['Calculating stim # ' num2str(iLoc) ' of ' num2str(nLocs)]);
    % checking to see if files exist
    fname = ['stimbuf_left_' num2str(dir_tested(1,iLoc)) '_' num2str(dir_tested(2,iLoc))];
    if exist1([stimuli_dir fname])                % if files exist
        iLocInner = iLocInner+1;
        dir(:,iLocInner) = dir_tested(:,iLoc);
        % calc for unalterred IRs - 2source files are saved to disk with name of loc1 (no offset required in testname)
        if 1
            fname = ['stimbuf_left_' num2str(dir_tested(1,iLoc)) '_' num2str(dir_tested(2,iLoc))];
            fid = fopen(['e:\kip\noisetc\temp0\' fname],'r');
            [Lnoi, count] = fread(fid, inf,'int16');
            Lnoi = Lnoi(3001:6254);
            fclose(fid);
            fname = ['stimbuf_right_' num2str(dir_tested(1,iLoc)) '_' num2str(dir_tested(2,iLoc))];
            fid = fopen(['e:\kip\noisetc\temp0\' fname],'r');
            [Rnoi, count] = fread(fid, inf,'int16');
            Rnoi = Rnoi(3001:6254);
            fclose(fid);
            
            tempL = ERBFilterBankB(Lnoi, fcoefs).* Factormat1;		% has dimensions n_cF x length(noi)
            tempR = ERBFilterBankB(Rnoi, fcoefs).* Factormat1;
            [ITD_std, time, IPD_std] = calcitd(tempL,tempR, cF, Fs, ones(size(cF)));
            
        else    
            
            ind = find(dir_HRTF(1,:) == dir_tested(1,iLoc)+15 & dir_HRTF(2,:) == dir_tested(2,iLoc)-15);
            tempL = ERBFilterBankB(mtlrch(FN_HRTF,ind*2-1), fcoefs) .* Factormat2;		% has dimensions n_cF x length(noi)
            tempR = ERBFilterBankB(mtlrch(FN_HRTF,ind*2), fcoefs) .* Factormat2;
            [ITD_std, time, IPD_std] = calcitd(tempL,tempR, cF, Fs, ones(size(cF)));
        end 
        peak_std(iLocInner) = find(ITD_std == max1(ITD_std));
        peakHeight_std(iLocInner) = ITD_std(peak_std(iLocInner));
        Height = peakHeight_std(iLocInner);
        % find width at half-height
        i = peak_std(iLocInner);
        while Height >= (peakHeight_std(iLocInner)/2) & (i>1)      %test to the right
            Height = ITD_std(i);
            i = i-1;
        end
        R = i;  i = peak_std(iLocInner);
        while Height >= (peakHeight_std(iLocInner)/2) & (i>1)      %test to the left
            Height = ITD_std(i);
            i = i+1;
        end
        L = i;  peakWidth_std(iLocInner) = L-R;
        
        %% calc for altered IRs
        %%%%% load stimuli
        fname = ['stimbuf_left_' num2str(dir_tested(1,iLoc)) '_' num2str(dir_tested(2,iLoc))];
        fid = fopen([stimuli_dir fname],'r');
        [Lnoi, count] = fread(fid, inf,'int16');
        Lnoi = Lnoi(3001:6254);
        fclose(fid);
        fname = ['stimbuf_right_' num2str(dir_tested(1,iLoc)) '_' num2str(dir_tested(2,iLoc))];
        fid = fopen([stimuli_dir fname],'r');
        [Rnoi, count] = fread(fid, inf,'int16');
        Rnoi = Rnoi(3001:6254);
        fclose(fid);
        
        tempL = ERBFilterBankB(Lnoi, fcoefs).* Factormat1;		% has dimensions n_cF x length(noi)
        tempR = ERBFilterBankB(Rnoi, fcoefs).* Factormat1;
        [ITD_altIR, time, IPD_altIR] = calcitd(tempL,tempR, cF, Fs, ones(size(cF)));
        
        % find altIR peak data at normal HRTF peak
        peakHeight_altIR_std(iLocInner) = ITD_altIR(peak_std(iLocInner));
        Height = peakHeight_altIR_std(iLocInner);
        % find width at half-height
        i = peak_std(iLocInner);
        while Height >= (peakHeight_altIR_std(iLocInner)/2) & (i>1)      %test to the right
            Height = ITD_altIR(i);
            i = i-1;
        end
        R = i;  i = peak_std(iLocInner);
        while Height >= (peakHeight_altIR_std(iLocInner)/2) & (i>1)      %test to the left
            Height = ITD_altIR(i);
            i = i+1;
        end
        L = i;  peakWidth_altIR_std(iLocInner) = L-R;
        
        
        % find altIR peak data at it's own peak
        peak_altIR(iLocInner) = find(ITD_altIR == max1(ITD_altIR));
        peakHeight_altIR(iLocInner) = ITD_altIR(peak_altIR(iLocInner));
        Height = peakHeight_altIR(iLocInner);
        % find width at half-height
        i = peak_altIR(iLocInner);
        while Height >= (peakHeight_altIR(iLocInner)/2) & (i>1)      %test to the right
            Height = ITD_altIR(i);
            i = i-1;
        end
        R = i;  i = peak_altIR(iLocInner);
        while Height >= (peakHeight_altIR(iLocInner)/2) & (i>1)      %test to the left
            Height = ITD_altIR(i);
            i = i+1;
        end
        L = i;  peakWidth_altIR(iLocInner) = L-R;
        
        eval(['save ' saveFN ' peak*'])
    end
end


plot_diam(peakH_altIR_std-peakH_std,dir,1);
figure; hold on
plot(ones(1,iLocInner),peakH_altIR_std-peakH_std,'g.')