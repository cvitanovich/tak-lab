function hpequal(catname,outname,displaymode);
% HPEQUAL computes the inverse filters to the
%         measured headphone responses specified
%         by the given catname. The coefficients
%         of the inverse filters are written to
%         the specified outfile.
%         HPEQUAL(catname,outfile,displaymode)
%         displaymode=1 shows the results.

if nargin<3 displaymode=0; end;
fftlen=2048;

[hl,al]=readhp(catname,'left','norm','on');
hl=hl*al(2);
[hr,ar]=readhp(catname,'right','norm','on');
hr=hr*ar(2);

N=128;
tshift=96;
hlinv=lsqinv(hl,N,tshift,1);
hrinv=lsqinv(hr,N,tshift,1);

if displaymode==1
  figure;
  subplot(2,1,1);
  plot(hl);
  title('headphone, left');
  subplot(2,1,2);
  plot(hr);
  title('headphone, right');

  Hl=fft(hl,fftlen);
  Hr=fft(hr,fftlen);

  Hlinv=fft(hlinv,fftlen);
  Hrinv=fft(hrinv,fftlen);

  Heql=Hlinv.*Hl; Heql=Heql(1:fftlen/2);
  Heqr=Hrinv.*Hr; Heqr=Heqr(1:fftlen/2);

  Hl=Hl(1:fftlen/2);
  Hr=Hr(1:fftlen/2);
  magplot([Hl(:) Hr(:)],samplfrq,'Headphone transfer function');
  grcntrl;

  magplot([Heql(:) Heqr(:)],samplfrq,'Equalized headphone transfer function');
  grcntrl;

  figure;
  subplot(2,1,1);
  plot(hlinv);
  title('Least square inverse for headphone, left');
  subplot(2,1,2);
  plot(hrinv);
  title('Least square inverse for headphone, right');
end;

info=['LSI-filters for headphone in ',catname];
s=sprintf('save %s hlinv hrinv info',outname);
disp(s);
eval(s);
