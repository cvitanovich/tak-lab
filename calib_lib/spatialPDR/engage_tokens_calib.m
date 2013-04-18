function engage_tokens_calib
% Engage Calibration of Sound Tokens
global PDR KNOWLES TDT

tokens = {'OCTAVE','BBN','GTONE'};
for cnt=1:length(tokens)
    eval(['left_snd = PDR.' tokens{cnt} '_left']);
    eval(['right_snd = PDR.' tokens{cnt} '_right']);
    [lo, hi] = det_cutoffs(PDR.scales_2_try_for_cutoffs, left_snd, right_snd);
    scales=lo:((hi-lo)/PDR.nscales):hi;
    avg_rms = test_scales(scales, left_snd, right_snd);
    both_dbs = (1/KNOWLES.coeffs(2))*log( (avg_rms - KNOWLES.coeffs(3)) ./ KNOWLES.coeffs(1) );
    
    % REGRESSION FIT
    h=figure;
    hold on;
    xrange = log10(1:1:32760);
    xes = log10(scales); % log10(scales) for lag speaker
    yes = both_dbs;
    colr = [1 0 0];
    [rsq, coefs] = regress_stats(xes,yes,0.05,xrange,colr,1);
    eval(['PDR.' tokens{cnt} '_RSQ = rsq']);
    eval(['PDR.' tokens{cnt} '_COEFS= coefs']);
    % title string
    title_string=[];
    title_string{1} = ['TEST SOUND IS: ' tokens{cnt}];
    title_string{2} = ['RSQ = ' num2str(rsq)];
    title(title_string,'FontSize',8);
    % axes labels
    xlabel('log10(scales)','FontSize',8);
    ylabel('SPL (dB, ABL)','FontSize',8);
    hold off;
    % save this figure
    set(h,'InvertHardcopy','off');
    fname = [PDR.data_path PDR.filename tokens{cnt} '_fit'];
    saveas(h,fname,'fig');
end
%write header information to file... saving global variables
CALIB_PDR=PDR;
save ([CALIB_PDR.data_path CALIB_PDR.filename '.mat'],'CALIB_PDR');
str{1} = 'Variables saved in: ';
str{2} = [CALIB_PDR.data_path CALIB_PDR.filename '.mat'];
hMsg=msgbox(str); uiwait(hMsg);

function [lo, hi] = det_cutoffs(scales2try,left_snd,right_snd)
global KNOWLES
avg_rms = test_scales(scales2try,left_snd,right_snd);
% select cutoffs
both_dbs = (1/KNOWLES.coeffs(2))*log( (avg_rms - KNOWLES.coeffs(3)) ./ KNOWLES.coeffs(1) );
[lo, hi] = select_cutoffs(log10(scales_lead),both_dbs);
if hi>log10(32760)
    hi=log10(32760);
end
if lo<0
    lo=0;
end


function [lo, hi] = select_cutoffs(xes,yes)
% Plots data in a scatter plot and asks user to input lo/hi cuttoffs
% for fitting linear portion of data range
screen_size = get(0, 'ScreenSize');
hFig=figure;
set(hFig, 'Position', [0.02*screen_size(3) 0.05*screen_size(4) 0.8*screen_size(3) 0.8*screen_size(4)] );

scatter(xes,yes,'.');
prompt = {'Enter Lower Cutoff','Enter Upper Cutoff'};
dlg_title = 'Input Cutoffs:';
num_lines = 1;
def = {'',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
lo=str2num(answer{1});
hi=str2num(answer{2});
close(hFig);
drawnow;

function avg_rms = test_scales(scales,left_snd,right_snd)
global TDT
for j=1:length(scales)
    S232('dropall'); % clear stack
    S232('qpushf',left_snd);
    S232('scale',scales(j));
    S232('qpop16',TDT.BUF_L1);
    S232('qpushf',right_snd);
    S232('scale',scales(j));
    S232('qpop16',TDT.BUF_R1);
    S232('seqplay',TDT.PLAY_SPEC);
    % recording voltage
    s232('seqrecord',TDT.REC_SPEC);
    s232('PD1arm',1);
    s232('PD1go',1);
    while(s232('PD1status',1))
        % wait for TDT to finish
    end
    % pop recorded voltage traces into matlab
    if exist('left_tmp')
        clear left_tmp;
    end
    if exist('right_tmp')
        clear right_tmp;
    end
    S232('qpush16', TDT.REC_L1);
    S232('qpop16', left_tmp);
    S232('qpush16', TDT.REC_R1);
    S232('qpop16', right_tmp);
    % calculate rms voltage
    left_rms = sqrt(sum(left_tmp.^2)./length(left_tmp));
    right_rms = sqrt(sum(right_tmp.^2).length(right_tmp));
    avg_rms(j) = (left_rms + right_rms)/2;
end

% REGRESSION STATS
function [RSQ, COEFFS] = regress_stats(xes,yes,alfa,xrange,col,info_flg)
% plot data
n=length(yes);
df=length(yes)-2;
for j = 1:n
    plot(xes,yes,'.','Color',col);
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
hold on;
plot(xrange,Y,'c-','Color',col);
plot(xrange,Y-DELTA,'r--','Color',col);
plot(xrange,Y+DELTA,'r--','Color',col);
hold off;

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

if info_flg
    % plot info
    eqn = ['Fit: y = mx+b = ' num2str(m) 'x + ' num2str(b)];
    text(minX+0.1*(maxX-minX),minY+0.9*(maxY-minY),eqn,'FontSize',8);
    mcoeff=['m = ' num2str(m) ' [' num2str(COEFFS(1,2)) ',' num2str(COEFFS(1,3)) ']'];
    text(minX+0.1*(maxX-minX),minY+0.75*(maxY-minY),mcoeff,'FontSize',8);
    bcoeff=['b = ' num2str(b) ' [' num2str(COEFFS(2,2)) ',' num2str(COEFFS(2,3)) ']'];
    text(minX+0.1*(maxX-minX),minY+0.6*(maxY-minY),bcoeff,'FontSize',8);
end