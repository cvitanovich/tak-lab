function anarad(respfname,ILD,ITD,LEVEL,T,MAP,b_array, direc)

% ANARAD calculates weighting factors for ILD,ITD and LEVEL input and T output using radial basis functions

% ANARAD(respfname,ILD,ITD,LEVEL,T,MAP,b_array, dir)
% 
% ILD interaural level difference
% ITD interaural time difference
% LEVEL monaural log mag spectrum
% T response of the neuron
% MAP coordinates of parameter
% b_array contains numbers frequency groups
% respfname file with neuronal response
% dir: directory for output files

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date May 1996


respfname=lower(respfname);
respfname=strrep(respfname,'ca0','c')


direc =[direc '/radbas']

cont=0;
w_off=1; 
%	ILD	ITD	LL	LR	LL_N	LR_N
%x=[	1	0	0	0	0	0	;...
%	0	1	0	0	0	0	;...
x=[	1	1	0	0	0	0	;...
%	0	0	1	0	0	0	;...
	1	0	1	0	0	0	;...
	0	1	1	0	0	0	;...
	1	1	1	0	0	0	;...
%	0	0	0	1	0	0	;...
	1	0	0	1	0	0	;...
	0	1	0	1	0	0	;...
	1	1	0	1	0	0	;...
	0	0	1	1	0	0	;...
	1	0	1	1	0	0	;...
	0	1	1	1	0	0	;...
	1	1	1	1	0	0	;...
%
	0	0	0	0	1	0	;...
	1	0	0	0	1	0	;...
	0	1	0	0	1	0	;...
	1	1	0	0	1	0	;...
	0	0	0	0	0	1	;...
	1	0	0	0	0	1	;...
	0	1	0	0	0	1	;...
	1	1	0	0	0	1	;...
	0	0	0	0	1	1	;...
	1	0	0	0	1	1	;...
	0	1	0	0	1	1	;...
	1	1	0	0	1	1	;...	
]


res = 15;
% load T
[vT,pT]=max(T); eleT=MAP(1,pT); aziT=MAP(2,pT);

%scale T : 
T=T-min(T); T=T/max(T);

% HRTF

if (exist('ILD')==1)
	disp('ILD exists')
 else
  	ana
end


[n_fgrp,ndir] = size(ILD);

fg=1;
lg=19;
n_fgrp=lg-fg+1;

% normalize ILD ITD
disp('normalize ILD and ITD');

N_ILD = div(ILD(fg:lg,:), max(abs(ILD(fg:lg,:).')));
%N_ILD=normc(ILD(fg:lg,:));
N_ITD = div(ITD(fg:lg,:), max(abs(ITD(fg:lg,:).')));
%N_ITD=normc(ITD(fg:lg,:));
% normalize LEVEL

disp('normalize LEVEL');
[m,n]= size(LEVEL);

N_LEVEL=LEVEL;

%maxLEVEL = max(LEVEL);
%meanLEVEL=mean(LEVEL);
%for dir=1:n
%	%N_LEVEL(:,dir)=LEVEL(:,dir)-maxLEVEL(dir)*ones(m,1);
%	N_LEVEL(:,dir)=LEVEL(:,dir)-meanLEVEL(dir)*ones(m,1);
%end

%absminN_LEVEL=min(min(N_LEVEL));
%N_LEVEL=N_LEVEL-absminN_LEVEL;

%maxN_LEVEL=max(max(N_LEVEL));
%N_LEVEL=N_LEVEL/maxN_LEVEL;


% separate LEVEL in LL (left) LR (right)
disp('separate LEVEL');

LL=N_LEVEL(fg:lg,1:2:n);
LR=N_LEVEL(fg:lg,2:2:n);

meanLL=mean(LL.').';
meanLR=mean(LR.').';

for dir=1:n/2
	LL(:,dir)=LL(:,dir)-meanLL;
	LR(:,dir)=LR(:,dir)-meanLR;
end

mLL=max(max(LL));
mLR=max(max(LR));

LL=LL/mLL;
LR=LR/mLR;

LL_N=normc(LL(:,:));
LR_N=normc(LR(:,:));

[n_cond,n_col]=size(x);


for cond=1:n_cond
	fname=[respfname '.c' int2str(cond)]
	pname=[respfname setstr(97+cond-1) '.eps'] 

 	disp(['calculating condition ' int2str(cond)]);

	P=[]; n_param=0;
	%P=[x(cond,1)*N_ILD(fg:lg,:); x(cond,2)*N_ITD(fg:lg,:); x(cond,3)*LL(fg:lg,:); x(cond,4)*LR(fg:lg,:)];
        if (x(cond,1)==1) P=[P;N_ILD(:,:)]; n_param=n_param+1; end;
        if (x(cond,2)==1) P=[P;N_ITD(:,:)]; n_param=n_param+1; end;
        if (x(cond,3)==1) P=[P;   LL(:,:)]; n_param=n_param+1; end;
        if (x(cond,4)==1) P=[P;   LR(:,:)]; n_param=n_param+1; end;
	if (x(cond,5)==1) P=[P;   LL_N(:,:)];n_param=n_param+1;end;	
	if (x(cond,6)==1) P=[P;   LR_N(:,:)];n_param=n_param+1;end;
%	n_param=n_param-1;

	niter=25;
	nmaxneuron=n_param;
	sse_goal=1;
	spread=NaN;
	dp=[niter nmaxneuron sse_goal spread]
	[W,b,W2,b2,nr,tr]=solverb(P,T,dp);
	a=simurb(P,W,b,W2,b2);
	
	[va,pa]=max(a); elea=MAP(1,pa); azia=MAP(2,pa);
	
	e=T-a;
	me=mean(e)
        msqe=mean(e.^2)
size(W)
figure
hintonw(W);
figure
hintonw(W2);
pause
if (w_off==0)
	w_r = reshape(W,n_fgrp,n_param);

	counter=1;
	if(x(cond,1)==1) w_ILD=w_r(:,counter); counter=counter+1; end
	if(x(cond,2)==1) w_ITD=w_r(:,counter); counter=counter+1; end
	if(x(cond,3)==1) w_LL =w_r(:,counter); counter=counter+1; end
	if(x(cond,4)==1) w_LR =w_r(:,counter); counter=counter+1; end
	if(x(cond,5)==1) w_LL_N =w_r(:,counter);counter=counter+1;end
	if(x(cond,6)==1) w_LR_N =w_r(:,counter);counter=counter+1;end

end
        cond_vec=x(cond,:)

	f=figure('NumberTitle','off','Name',fname);
	set(f,'Position',[50 50 700 500]);
	subplot(3,2,1); 
	axis off;
	text(0.1,1.0,['Name: ' respfname]);
	text(0.1,0.5,['MSE: ' num2str(msqe)]);
	
	drawnow
	subplot(3,2,2);
	axis off

	text(0,1.0, ['ILD:']), 		text(0.5,1.0,int2str(cond_vec(1)));
	text(0,0.8, ['ITD:']), 		text(0.5,0.8,int2str(cond_vec(2)));
	text(0,0.6, ['Left:']),		text(0.5,0.6,int2str(cond_vec(3)));
	text(0,0.4, ['Right:']),	text(0.5,0.4,int2str(cond_vec(4)));
	text(0,0.2, ['Left norm:']),	text(0.5,0.2,int2str(cond_vec(5)));
	text(0,0.0, ['Right norm:']),	text(0.5,0.0,int2str(cond_vec(6)));
  		
	drawnow	
	subplot(3,2,3);
	flatplot(T,MAP,cont,1,res,1,0,0,0,['Single Unit Response T, Max: ' num2str(eleT) ' Ele ' num2str(aziT) ' Azi'],1);
	% text(100,-37.5,['Max: ' num2str(eleT) ' Ele ' num2str(aziT) ' Azi']);
	subplot(3,2,4);
	flatplot(a,MAP,cont,1,res,1,0,0,0,['Simulation a, Max: ' num2str(elea) ' Ele ' num2str(azia) ' Azi'],1);
	%text(100,-37.5,['Max: ' num2str(elea) ' Ele ' num2str(azia) ' Azi']);

	subplot(3,2,5);
	flatplot(abs(e),MAP,cont,1,res,1,0,0,0,'Error T-a',1);
    
	drawnow;     
if (w_off==0)

	subplot(3,2,6);
	plot(b_array(fg:lg),w_r);
	label_str=zeros(1,1);
	for k=fg:lg label=int2str(bark2hz(b_array(k))); label_str(k-fg+1,1:length(label))=label; end;	
	label_str=label_str(1:4:n_fgrp,:);
	h=gca;
	set(h,'XTick',b_array(fg:4:lg),'XTickLabels',label_str);

	legstr=[];
 	legp=1;
	if (x(cond,1)==1) legstr(legp,1:3)='ILD'; legp=legp+1; end;
	if (x(cond,2)==1) legstr(legp,1:3)='ITD'; legp=legp+1; end;
	if (x(cond,3)==1) legstr(legp,1:4)='Left'; legp=legp+1;end;
	if (x(cond,4)==1) legstr(legp,1:5)='Right';legp=legp+1;end;
	if (x(cond,5)==1) legstr(legp,1:9)='Left norm';legp=legp+1;end;
	if (x(cond,6)==1) legstr(legp,1:10)='Right norm';legp=legp+1;end;

        legp=legp-1;
	
	if (legp==1) legend(legstr(1,:)); end
	if (legp==2) legend(legstr(1,:),legstr(2,:)); end
	if (legp==3) legend(legstr(1,:),legstr(2,:),legstr(3,:)); end
	if (legp==4) legend(legstr(1,:),legstr(2,:),legstr(3,:),legstr(4,:)); end
   
	ymax=max(max(w_r));     ymin=min(min(w_r));
	xmin=min(min(b_array(fg:lg))); xmax=max(max(b_array(fg:lg)));
	axis([xmin xmax ymin ymax]);
	xlabel('Frequency [Hz]');
	ylabel('Weighting Factor');
	drawnow;
end
	action=['save ' direc '/' fname ' a W b e P T cond_vec cond ']
	eval(action);
	action=['print -depsc '  direc '/' pname]
	action=strrep(action,'/','\')
	eval(action);
	close(f);
end
