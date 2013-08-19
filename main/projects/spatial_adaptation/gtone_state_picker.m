function [resulting_rms, state_list, success]=gtone_state_picker(gtone_coefs,Fs,dur,buf_pts,nStates,hrtf_lt, hrtf_rt,target_rms)
% picks states that avoid clipping given desired parameters

hWait=waitbar(0,'Grabbing states for adaptor...');

% a random list of state values to pick from:
S=round(1000*rand(1,nStates*10));

cnt=0; cyc=0; success=1;
resulting_rms=[];
state_list=[];
while(cnt<nStates)
    cyc=cyc+1;
    % make a gammatone with state=cyc
    rand('state',S(cyc));
    snd=rand(1,ceil(dur*Fs+buf_pts));
    snd=filtfilt(gtone_coefs,1,snd);
    snd=snd./max(abs(snd));
    % convolve with hrtf coefficients for desired location
    lt_snd=filter(hrtf_lt,1,snd);
    rt_snd=filter(hrtf_rt,1,snd);
    % set rms level:
    avg_rms=(sqrt(mean(lt_snd.^2))+sqrt(mean(rt_snd.^2)))/2;
    lt_snd=(target_rms/avg_rms).*lt_snd;
    rt_snd=(target_rms/avg_rms).*rt_snd;
    mx=max(abs([lt_snd rt_snd]));
    if(mx<.9999)
        cnt=cnt+1;
        state_list=[state_list S(cyc)];
        resulting_rms(cnt)=(sqrt(mean(lt_snd.^2))+sqrt(mean(rt_snd.^2)))/2;
    end
    if(cyc>10*nStates)
        warndlg('Having trouble finding enough states!!!')
        success=-1;
        close(hWait)
        break;
    end
    waitbar(cnt/nStates,hWait);
end
close(hWait)