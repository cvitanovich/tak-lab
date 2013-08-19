function AMStim()
% Note: Noise1 is lead and Noise2 is lag
% AMStim (Adapted from Caitlin/Brian's head turn exp't)
% Generate stimuli for LDS session

global PDR
PDR.RPs=[];
%//////////////////////////////////////////
% low pass filter design:
fc=150;
fs=PDR.stim_Fs;
order=2;
[Bs,As]=filt_butter(fc,fs,order);
% Make stimulus or stimuli
SR = PDR.stim_Fs;
% Increase duration by delay
% so that the stimuli can be gated/windowed (see below)
StimDur = PDR.SOUNDS_length + abs(PDR.SOUNDS_carrier_delay); % should be in ms, add delay

%for sound generation (creates extended BBN that is then curtailed to
%the desired duration)
% Use PAD = 10 for sounds with less than 250ms total duration!
% Otherwise use PAD = 1
if StimDur <= 250
    PAD = 10;
else
    PAD = 1;
end

StimPnts = round((StimDur/1000)*SR);
DelayPnts = round((abs(PDR.SOUNDS_carrier_delay)/1000)*SR);
DelayPnts = max(1,DelayPnts);

numsnds=length(PDR.SOUNDS_rand_states);

% create a family of sounds the same envelope parameters, but with
% different carriers
hWait = waitbar(0,'Generating Sounds');
for i0=1:numsnds
    
    % Make Noise #1
    x = getBBN(PDR.SOUNDS_rand_states(i0), SR, StimDur, PDR.SOUNDS_carrier_bandwidth(1), PDR.SOUNDS_carrier_bandwidth(2));
    % scale to RMS = 0.23591 !!!!!!!!
    x=rms_scale(x,0.23591);
    Noise1 = x(1:StimPnts)';
    % Scale to desired amplitude
    Noise1 = Noise1 * 50; % scale by 50 so that speaker scaling factors aren't too big (well, small)!
    Noise1 = Noise1 * 10^(PDR.SOUNDS_amplitude/20)*0.00002; % amplitude re: 20 uPa   
    
    % copy Noise #1 to make Noise #2
    Noise2 = Noise1;
    
    % Make Envelope #1
    Env1 = getBBN(PDR.SOUNDS_states(1), SR, PAD*StimDur, PDR.SOUNDS_env_bandwidth(1), PDR.SOUNDS_env_bandwidth(2));
    Env1 = Env1(1:StimPnts)';
    Env1 = Env1 - min(Env1); Env1 =  Env1 / max(Env1);  % normalize
    Env1 = Env1 .* (PDR.SOUNDS_env_depth/100); % AM depth (as a percent)
    Env1 = Env1 + 1 - (PDR.SOUNDS_env_depth/100);
    
    % copy Envelope #1 to make Envelope #2
    Env2 = Env1;
    
    % get de-correlateing envelope #1
    DCEnv1 = getBBN(PDR.SOUNDS_states(2), SR, PAD*StimDur, PDR.SOUNDS_env_bandwidth(1), PDR.SOUNDS_env_bandwidth(2));
    DCEnv1 = DCEnv1(1:StimPnts)';
    DCEnv1 = DCEnv1 - min(DCEnv1); DCEnv1 =  DCEnv1 / max(DCEnv1);  % normalize
    DCEnv1 = DCEnv1 .* (PDR.SOUNDS_env_depth/100); % AM depth (as a percent)
    DCEnv1 = DCEnv1 + 1 - (PDR.SOUNDS_env_depth/100);
    
    % get de-correlateing envelope #2
    DCEnv2 = getBBN(PDR.SOUNDS_states(3), SR, PAD*StimDur, PDR.SOUNDS_env_bandwidth(1), PDR.SOUNDS_env_bandwidth(2));
    DCEnv2 = DCEnv2(1:StimPnts)';
    DCEnv2 = DCEnv2 - min(DCEnv2); DCEnv2 =  DCEnv2 / max(DCEnv2);  % normalize
    DCEnv2 = DCEnv2 .* (PDR.SOUNDS_env_depth/100); % AM depth (as a percent)
    DCEnv2 = DCEnv2 + 1 - (PDR.SOUNDS_env_depth/100);
    
    % ** mix de-correlateing envelopes in with original envelopes **
    Env1DC = ((1-PDR.SOUNDS_mix).*DCEnv1) + ((PDR.SOUNDS_mix).*Env1);
    Env2DC = ((1-PDR.SOUNDS_mix).*DCEnv2) + ((PDR.SOUNDS_mix).*Env2);
    
    % !!! re-normalize the envelopes
    Env1DC = Env1DC - min(Env1DC); Env1DC =  Env1DC / max(Env1DC);  % normalize
    Env1DC = Env1DC .* (PDR.SOUNDS_env_depth/100); % AM depth (as a percent)
    Env1DC = Env1DC + 1 - (PDR.SOUNDS_env_depth/100);
    Env2DC = Env2DC - min(Env2DC); Env2DC =  Env2DC / max(Env2DC);  % normalize
    Env2DC = Env2DC .* (PDR.SOUNDS_env_depth/100); % AM depth (as a percent)
    Env2DC = Env2DC + 1 - (PDR.SOUNDS_env_depth/100);
    
    % Delay to make "lead" and "lag"
    % DO NOT SHIFT ENVELOPES!!! (only carriers have a delay)
    Noise1 = [Noise1 zeros(1,DelayPnts)];
    Noise2 = [zeros(1,DelayPnts) Noise2];
    
    % window (i.e., gate) out onset/offset time disparity
    % DO NOT SHIFT ENVELOPES!!! (only carriers have a delay)
    if(PDR.SOUNDS_carrier_delay>=0)
        N1 = Noise1(DelayPnts:StimPnts);
        N2 = Noise2(DelayPnts:StimPnts);
    else
        PDR.CARRIERS_SWITCHED=1;
        N2 = Noise1(DelayPnts:StimPnts);
        N1 = Noise2(DelayPnts:StimPnts);
    end
    Noise1=N1;
    Noise2=N2;
    % Fix envelope points to match carriers!
    Env1DC = Env1DC(1:length(Noise1));
    Env2DC = Env2DC(1:length(Noise2));
    
    % ramp envelopes on/off
    Env1DC = rampSounds(Env1DC, SR, PDR.SOUNDS_ramp); % stim envelope
    Env2DC = rampSounds(Env2DC, SR, PDR.SOUNDS_ramp); % stim envelope
    
    % set rms amplitude for envelopes:
    Env1DC=rms_scale(Env1DC,PDR.SOUNDS_rms);
    Env2DC=rms_scale(Env2DC,PDR.SOUNDS_rms);
    
    
    % CALCULATE CORRELATION COEFFICIENTS (Final LEAD/LAG Envelopes)
    % get hilbert envelope (lead):           
    env1(i0,:)=abs(hilbert(Noise1 .* Env1DC));
    % butterworth lowpass filter (lead):
    for y0=1:size(Bs,1); env1(i0,:)=filter(Bs(y0,:),As(y0,:),env1(i0,:)); end
    
    % get hilbert envelope (lag);
    env2(i0,:)=abs(hilbert(Noise2 .* Env2DC));
    % butterworth lowpass filter (lag):
    for y1=1:size(Bs,1); env2(i0,:)=filter(Bs(y1,:),As(y1,:),env2(i0,:)); end
    
    % calculate correlation coefficient (lead/lag) for filtered hilbert
    % envelopes (this is the true envelope after convolving with
    % carrier)
    
    R=corrcoef(env1(i0,:),env2(i0,:));
    Res(i0)=R(2,1);
    PDR.SOUNDS_env_correlations{i0} = R(2,1); % correlation coefficients

    % CALCULATE RP:
    [leadP, lagP] = SpikeProbAM(env1(i0,:),env2(i0,:),PDR.stim_Fs);
    PDR.RPs(i0)= leadP/(leadP+lagP);
    
    % multiply with envelopes
    Noise1 = Noise1 .* Env1DC;
    Noise2 = Noise2 .* Env2DC;
    
    LEAD_sounds{i0} = Noise1;
    LAG_sounds{i0} = Noise2;
    
    waitbar(i0/numsnds,hWait);
end

close(hWait);
if(PDR.SOUNDS_carrier_delay<0)
    h=warndlg('Lead/Lag Carriers Switched!!!');
    uiwait(h);
end
% find the most highly correlated sounds (low pass filtered hilbert
% envelopes)
clims=[-1 1];

R0=corrcoef(env1');
for q0=1:numsnds
    R0(q0,q0)=-inf;
end
tmp=reshape(R0,1,size(R0,1)*size(R0,2));
tmp=sort(tmp);
tmp=fliplr(tmp);
chk=0; cnt=0;
while chk==0
    cnt=cnt+1;
    [i,j]=find(R0>=tmp(cnt));
    if length(unique(i)) >= PDR.SOUNDS_num_carriers
        snds=unique(i);
        snds=snds(1:PDR.SOUNDS_num_carriers);
        PDR.SOUNDS_carrier_states=PDR.SOUNDS_rand_states(snds);
        env_new1=env1(snds,:);
        env_new2=env2(snds,:);
        chk=1;
        break;
    end
end

for q1=1:length(snds)
    % zero out LEAD sounds if testing MAAs
    if PDR.maa_test_flag
        PDR.LEAD_sounds{q1}=zeros(size(LEAD_sounds{snds(q1)}));
    else
        PDR.LEAD_sounds{q1}=LEAD_sounds{snds(q1)};
    end
    PDR.LAG_sounds{q1}=LAG_sounds{snds(q1)};
end

scrn=get(0,'ScreenSize');
hTemp=figure('name','hist','Position',[0.1*scrn(3) 0.1*scrn(4) 0.8*scrn(3) 0.8*scrn(4)]); hold on;

if PDR.SOUNDS_num_carriers > 1
    % correlations of sound tokens
    subplot(3,2,1);
    Rnew1=corrcoef(env_new1');
    Rnew2=corrcoef(env_new2');
    Ravg=(Rnew1+Rnew2)./2;
    
    imagesc(Ravg,[min(min(Ravg)) max(max(Ravg))]); colormap bone; colorbar
    axis square
    title('Avg. Correlations Between Sound Tokens');
end
if PDR.SOUNDS_num_carriers > 1
    if ~PDR.maa_test_flag
        % plot histogram of calculated lead/lag correlations
        % (Just the ones that will be used)
        subplot(3,2,3);
        hold on;
        binz=(min(Res(snds))-0.05):0.01:(max(Res(snds))+0.05);
        out=hist(Res(snds),binz);
        bar(binz,out);
        xlabel('Correlation Coefficient')
        ylabel('Frequency')
        title('Histogram of Lead/Lag Correlations (low pass filtered envelopes)')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','w','EdgeColor','none')
        axis([-1 1 0 ceil(1.1*max(out))]);
        %axis([(min(Res(snds))-0.05) (max(Res(snds))+0.05) 0 ceil(1.1*max(out))])
        set(gca,'Color',[0 0 0]);
        % hist of RPs
        subplot(3,2,5);
        hold on;
        binz=0:0.01:1;
        out=hist(PDR.RPs,binz);
        bar(binz,out);
        xlabel('RP');
        ylabel('Frequency');
        title('Histogram of Response Probabilities (low pass filtered envelopes)');
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r','EdgeColor','none')
        axis([0 1 0 ceil(1.1*max(out))])
        set(gca,'Color',[0 0 0]);
    else
        subplot(3,2,3);
        axis off
        text(0.15,0.35,'MAA EXPERIMENT','FontSize',30,'Color',[1 0 0])
    end
end
% lead and lag envelopes
subplot(3,2,2); hold on;
for q2=1:length(snds)
    % do not plot leads if just testing MAAs
    if ~PDR.maa_test_flag
        A=PDR.LEAD_sounds{q2};
        A_env=abs(hilbert(A));
        for y0=1:size(Bs,1); A_env=filter(Bs(y0,:),As(y0,:),A_env); end
        plot(A_env,'Color',[q2/length(snds) 0 1-q2/length(snds)],'LineWidth',2);
    end
    % plot low pass filtered lag envelopes:
    B=PDR.LAG_sounds{q2};
    B_env=abs(hilbert(B));
    for y0=1:size(Bs,1); B_env=filter(Bs(y0,:),As(y0,:),B_env); end
    plot(B_env,'Color',[q2/length(snds) 0 1-q2/length(snds)],'LineWidth',1,'LineStyle',':');
end
title('Low Pass Filtered Envelopes (lead/lag)');
set(gca,'Color',[0 0 0]);


% wait for user
uicontrol('Style', 'pushbutton', 'String', 'Close and Continue',...
    'Units','Normalized','Position', [0.6 0.1 0.25 0.1],...
    'BackgroundColor','g','Callback', 'close(''hist'')');
uiwait(hTemp)
