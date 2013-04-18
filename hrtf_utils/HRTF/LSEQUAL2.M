function lsequal2(displaymode,listname)
% LSEQUAL2 Compute random phase signals with loadspeaker
%         equalization
%         for the loudspeakers for HRTF measurement.
%         LSEQUAL(displaymode,listname)
%         diplaymode=1 shows the inverse filters, 
%         listname = name of a file with names of mtlsig-files
%                    name of output mtlsig-files 
%                     (inverse TF) will be ls_iX.sc1
%                     (X = number of file)

%modified by Klaus Hartung 28.6.95

% uses lsqinv of A. Raab

if nargin==0 displaymode=0; end;


j=sqrt(-1);
M=256+128;                            % length for inverse filters


file_list=readlist(listname);

[n_files, str_len]=size(file_list);

for i=1:n_files 
                     % for all speakers
filename=deblank(file_list(i,:));
[filetype,info_blocks,...
           n_channels,n_lines, sampling_rate,...
           first_line,last_line,n_directions,...
           comment1, comment2] = mtlrh(filename);

N=n_lines;
samplefrq=sampling_rate;

H=mtlrch(filename,1);

h=rifft(H);
dummy=h(280:535);
h=h*0;
h(1:256)=dummy;

  hi=lsqinv(h,M,0);                  % compute inverse filter
  [m,n]=size(hi);
  hi_l=zeros(1,1024);
  hi_l(1:m)=hi;
  Hm=H';

  Hi=rfft(hi_l);
  A=max(abs(Hm.*Hi));  % scaling factor;
  Hi=Hi/A;
  hi_l=hi_l/A;  
if displaymode==1
    Heq=Hi.*Hm; Heq=Heq(1:N/2);
    gtitle=sprintf('impulse response of inverse filter for speaker %i',i);
    irplot(hi,samplefrq,gtitle);   grcntrl        
% impulse response of inverse filter


    gtitle=sprintf('Total transfer function loudspeaker %i',i); 
    magplot(Heq(:),samplefrq,gtitle); grcntrl     
% display equalized transfer function


    disp('press any key...');
    pause;
    close; close;
  end;

filename=sprintf('ls_i%i.sc1',i);
mtlwh(filename,filetype,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);
channel=Hi;
mtlwch(filename,channel,1);

end;

end;