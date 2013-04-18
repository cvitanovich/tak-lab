function grpdelay=tgroup(H,samplerate);
% TGROUP Group delay for vector H.
%        grpdelay=TGROUP(H,samplerate);
%        note: H should contain samples from zero frequency
%              up to HALF of the sampling frequency.
%        If samplerate is omitted, 1/2 will be used (for
%        compatibility with grpdelay.m).

%       d phi   d         I   R dI/dw - I dR/dw
% -tgr= ----- = -- arctan - = -----------------
%        dw     dw        R       R^2 + I^2
%
% with R=real(H), I=imag(H), w=omega

% Reference:
%    A. Fettweis
%    Elemente nachrichtentechnischer Systeme
%    B.G. Teubner Verlag, Stuttgart, 1990
%    p. 212, 307

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 25.06.1994 A. Raab, revised 27.08.1994, 16.09.1994 A. Raab

if nargin==1 samplerate=1; end;
if min(size(H))>1 error('H must be a vector'); end;
domega=pi*samplerate/length(H);
R=real(H);
I=imag(H);
dI=derive(I,domega);
dR=derive(R,domega);
grpdelay=-(R.*dI-I.*dR)./(abs(H).^2);
