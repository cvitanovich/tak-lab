function Hs=fsmooth2(H,w,fftlen)
% FSMOOTH2 full, third or twelth octave smoothing of transfer function H.
%         Hs = FSMOOTH2(H,width,fftlen)
%         H     --> transfer function (mag) to be smoothed
%         width --> width 12)'twelth' or 3)'third' or 1)'octave'
%         fftlen--> number of computed frequency samples

if min(size(H))>1 error('H must be a vector'); end;

p=find(fftlen/length(H)==[1,2]);
if isempty(p) error('length of H and fftlen does not match'); end;
if (p==1) H=H(1:fftlen/2); end;

if isempty(w) error('invalid width');end
K = 2^(1/w);
K=sqrt(K);

Hs=H(:);
nmax=length(H)-1;
for n=2:nmax
  n0 = round(n/K);
  n1 = round(n*K);
  n1 = min([n1,nmax]);
  Hs(n+1) = sum(H(n0+1:n1+1))/(n1-n0+1); 
end;

if p==1 Hs=[Hs;flipud(Hs)]; end;
