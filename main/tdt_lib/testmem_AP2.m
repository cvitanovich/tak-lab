% init TDT
TDT_init;
pause(1);
% determine max words in DAMA space:
s232('dropall');
w=s232('freewords'); % num 32-bit words
bytes=(w*4);
kbytes=(bytes/1024);
mbytes=kbytes/1024;
disp(['AP2 card has ' num2str(mbytes) ' MB of memory storage space.'])
disp('Trying to fill DAMA buffers with zeros...')
buf_words=(128*1024)/4; % 128K buffers
nbuffers=w/buf_words;
dbn=[];
tmp=NaN*ones(1,buf_words*nbuffers);
figure; hold on; times2try=[1 ...
        20 ...
        30 45 60 75 90 120 240 300 600 3600 14400]; 
lp=0;
cnt=NaN*ones(2,length(times2try));
for tm=times2try
    s232('trash'); s232('dropall');
    lp=lp+1;
    p_len=tm;
    
    for j=1:(nbuffers-1)
        s232('dropall');
        s232('dpush',buf_words);
        s232('value',0);
        % allot DAMA space
        dbn(j)=s232('_allot16',buf_words);
        s232('qpop16',dbn(j));
    end
    % fill remaining "stack space" with zeros:
    s232('dropall');
    s232('dpush',buf_words);
    s232('value',0);
    hWait=waitbar(0,['Pausing for ~' num2str(p_len) ' seconds!']);
    for t=1:p_len
        pause(1)
        waitbar(t/p_len,hWait);
    end
    close(hWait);
    
    % grab stored zeros and check for errors
    tmp(((nbuffers-1)*buf_words+1):(nbuffers*buf_words))=s232('pop16');
    for j=(nbuffers-1):-1:1
        s232('dropall');
        s232('qpush16',dbn(j));
        tmp(((j-1)*buf_words+1):(j*buf_words))=s232('pop16');
    end
    tmp2=find(tmp~=0);
    cnt(1,lp)=p_len;
    cnt(2,lp)=length(tmp2);
    scatter(cnt(1,lp),cnt(2,lp));
    %plot(tmp,'Color',[lp/length(times2try) 1 1-lp/length(times2try)]);
end