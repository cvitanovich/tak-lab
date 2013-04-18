function lsequal(displaymode,rphfile)
% LSEQUAL Compute random phase signals with loadspeaker equalization
%         for the AUDIODROM loudspeakers.
%         LSEQUAL(displaymode,rphfile)
%         diplaymode=1 shows the inverse filters, rphfile is the file-
%         name of the original random-phase signal (SIGMAPS SI1 format),
%         i.e. 'e:\dsv\ran_p10.si1'

if nargin==0 displaymode=0; end;

N=1024;
catname='lscat';                     % catalogue for loudspeakers
lsnumbers=[1:11];                    % execute all speakers

j=sqrt(-1);
M=256+128;                            % length for inverse filters


[rps,header]=readsi1(rphfile);         % random phase signal
RPS=fft(rps,N);
rpse=zeros(length(lsnumbers),N);       % allocate memory for equalized signals

for i=lsnumbers                        % for all speakers
  [h,attr]=readls(catname,i,'norm','on');
  h=fade(h,0,0,350,450);
  h=h*attr(2);                       % denormalize data
  hi=lsqinv(h,M,0);                  % compute inverse filter
  Hm=fft(h,N);                       % measured loudspeaker response
  Hi=fft(hi,N);
  A=max(abs(Hm(1:N/2).*Hi(1:N/2)));  % scaling factor;
  Hi=Hi/A;
  if displaymode==1
    Heq=Hi.*Hm; Heq=Heq(1:N/2);
    gtitle=sprintf('impulse response of inverse filter for speaker %i',i);
    irplot(hi,samplefrq,gtitle);           % impulse response of inverse filter
    gtitle=sprintf('Total transfer function loudspeaker %i',i); 
    magplot(Heq(:),samplefrq,gtitle);      % display equalized transfer function
    disp('press any key...');
    pause;
    close; close;
  end;
  RPSe=RPS(:).*Hi(:);                      % 'equalized' random phase signal
  rpse(i,:)=real(ifft(RPSe,N))';
end;

maxrpse=max(max(abs(rpse)));              % maximum of all signals
rpse=(2^15-1)*rpse/maxrpse;               % scale for integer format

for i=lsnumbers                           % write output data
  filename=sprintf('rph_ls%i.si1',i);
  writesi1(rpse(i,:),filename,header);
end;
