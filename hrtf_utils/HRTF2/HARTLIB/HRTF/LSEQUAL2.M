function lsequal2(displaymode,listname,nl,M,Wn)
% LSEQUAL2 Compute random phase signals with loadspeaker equalization
%         for the loudspeakers for HRTF measurement.
%         LSEQUAL(displaymode,listname,nl,M,Wn)
%         diplaymode=1 shows the inverse filters, 
%         listname = name of a file with names of mtlsig-files
%                    name of output mtlsig-files (inverse TF) will be ls_iX.sc1
%                     (X = number of file)
%         nl = number of lines in impulse response to use for computation (default nl=128)
%         M  = length of impules response of inverse system (default M=256)
%         Wn = cutoff frequency of low-pass filter 0< Wn< 1 performs a
%              if Wn is a two element vector, Wn = [W1 W2] bandpass-filter

%modified by Klaus Hartung 28.6.95

% uses lsqinv2 of Klaus Hartung

if nargin==0 displaymode=0; end;
if (nargin < 4)
   M=512;
  if (nargin < 3) 
    nl=128; 
  end;
end;

test=0;

%nl=128;
%M=128*1.5;                            % length for inverse filters

file_list=readlist(listname);

[n_files, str_len]=size(file_list);

i=1;
filename=deblank(file_list(i,:));

           [filetype,info_blocks,...
           n_channels,n_lines, sampling_rate,...
           first_line,last_line,n_directions,...
           comment1, comment2] = mtlrh(filename);

for i=1:n_files 
                     % for all speakers

  filename=deblank(file_list(i,:));
  N=n_lines;
  samplefrq=sampling_rate;

  H=mtlrch(filename,1);

  h=rifft(H);
  H=rfft(h);

  [h,t(i)]=h_max_en(h,nl);
  h=h(1:nl);

  if (exist('Wn')==1)  hi=lsqinv2(h,M,Wn);                  % compute inverse filter
    else hi=lsqinv2(h,M);
  end;

  hi=hi(:);
  hi=[hi;zeros(2*n_lines-M,1)];
  
  if (displaymode==1)
    heq=conv(h,hi);

    Heq=rfft(heq);

    gtitle=sprintf('impulse response of inverse filter for speaker %i',i);
    irplot(hi(:),samplefrq,gtitle);   grcntrl        % impulse response of inverse filter
    gtitle=sprintf('Total transfer function loudspeaker %i',i); 
 
   magplot(Heq(:),samplefrq,gtitle); grcntrl     % display equalized transfer function
       
    disp('press any key...');
    pause;

    close; close;
  end;

lines=2*n_lines;
filename=sprintf('ls_i%i.sr1',i);
filetype=2;
mtlwh(filename,filetype,info_blocks,n_channels,lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);
size(hi)
mtlwch(filename,hi,1);

if (test)
  hi1=mtlrch(filename,1);
  magplot(rfft(conv(hi1,h)),samplefrq,gtitle); grcntrl     % display equalized transfer function

  disp('press any key...');
  pause;
  
  close;
end;

end;

