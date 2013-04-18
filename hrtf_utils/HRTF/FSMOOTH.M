function Hs=fsmooth(H,width,fftlen)
% FSMOOTH octave or third octave smoothing of transfer function H.
%         Hs=FSMOOTH(H,width,fftlen)
%         H     --> transfer function to be smoothed
%         width --> 'third' or 'octave'
%         fftlen--> number of computed frequency samples

if min(size(H))>1 error('H must be a vector'); end;
p=find(fftlen/length(H)==[1,2]);
if isempty(p) error('length of H and fftlen does not match'); end;
if p==1 H=H(1:fftlen/2); end;
w=find(width(1)==['to']);
if isempty(w) error('invalid width'); end;
if width==1 K=2^(1/3); else K=2; end;
K=sqrt(K);

Hs=H(:);
nmax=length(H)-1;
for n=2:nmax
  n0=round(n/K);
  n1=round(n*K);
  n1=min([n1,nmax]);
  Hs(n+1)=sum(H(n0+1:n1+1))/(n1-n0); 
end;

if p==1 Hs=[Hs;flipud(Hs)]; end;
