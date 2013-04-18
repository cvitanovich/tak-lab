%function [w_ILD, w_ITD, w_LL, w_LR] = anares2 (ILD,ITD,LEVEL,T,MAP,b_array)
function [a_ILD, a_ITD, a_LL, a_LR] = anares2 (ILD,ITD,LEVEL,T,MAP,b_array)

% ANARES2 Calculates weigting factors for ILD ITD and LEVEL superposition for one neuron solving linear equations
%
% ANARES2 (ILD,ITD,LEVEL,T,MAP,b_array)
%
% ILD interaural level difference
% ITD interaural time difference
% LEVEL monaural log mag spectrum
% T response of the neuron
% MAP coordinates of parameter
% b_array contains numbers frequency groups


[n_fgrp,ndir] = size(ILD);

n_taps_l=5;
n_taps_t=5;

fg=14;
lg=22;

res = 15;
n_iter=5;


% load T

[vT,pT]=max(T); eleT=MAP(1,pT); aziT=MAP(2,pT);

%scale T : T=T+min(T); T=T/max(T);


n_fgrp=lg-fg+1;

% increase resolution of ILD ITD
%p=[];
%for k=1:n_fgrp
%   [ILD_h,MAP_h]=sphereinterp(ILD(k,:),MAP,5);
%   [m,n]=size(ILD_h);
%   temp=reshape(ILD_h,1,m*n);
%  p=[p; temp];   	
%end

%ILD=p;
%size(ILD)
%pause

% normalize ILD ITD
disp('normalize ILD and ITD');


fi=find(fg==b_array);
li=find(lg==b_array);

ITD_NORM=ITD(fi:li,:)/max(max(abs(ITD(fi:li,:))));
ILD_NORM=ILD(fi:li,:)/max(max(abs(ILD(fi:li,:))));

P_ITD=zeros(n_fgrp*n_taps_t,122);
P_ILD=zeros(n_fgrp*n_taps_t,122);

% place coding of time difference
for m=1:122
	for k=1:n_fgrp
		[N_ITD(:,k),new_delay]=value2place(ITD_NORM(k,m),n_taps_t);
		[N_ILD(:,k),new_delay]=value2place(ILD_NORM(k,m),n_taps_l);
	end
%surf(N_ILD)
%pause
P_ITD(:,m)=reshape(N_ITD,n_fgrp*n_taps_t,1);
P_ILD(:,m)=reshape(N_ILD,n_fgrp*n_taps_l,1);


end

% normalize LEVEL

disp('normalize LEVEL');
[m,n]= size(LEVEL);

N_LEVEL=LEVEL;


% separate LEVEL in LL (left) LR (right)

disp('separate LEVEL');

LL=N_LEVEL(fi:li,1:2:n);
LR=N_LEVEL(fi:li,2:2:n);

meanLL=mean(LL.').';
meanLR=mean(LR.').';

for dir=1:n/2
	LL(:,dir)=LL(:,dir)-meanLL;
	LR(:,dir)=LR(:,dir)-meanLR;
end

mLL=max(max(LL));
mLR=max(max(LR));

%LL=LL/mLL;
%LR=LR/mLR;

LL_N=normc(LL(:,:))/4;
LR_N=normc(LR(:,:));

%[n_cond,n_col]=size(x);

%ing condition ' int2str(cond)]);

	P=[]; n_param=0;
	P=[P;P_ILD(:,:)]; n_param=n_param+n_taps_l;
size(P) 
        P=[P;P_ITD(:,:)]; n_param=n_param+n_taps_t;
        P=[P;   LL_N(:,:)]; n_param=n_param+1;
        P=[P;   LR_N(:,:)]; n_param=n_param+1;
%	n_param=n_param-1;


	[W,b]=solvelin(P,T);
	a=simulin(P,W,b);
	
	[va,pa]=max(a); elea=MAP(1,pa); azia=MAP(2,pa);
	
	e=T-a;
	me=mean(e)
        msqe=mean(e.^2)
	w_r = reshape(W,n_fgrp,n_param);        


	counter=1;
	w_ILD=w_r(:,counter:counter+n_taps_l-1); counter=counter+n_taps_l;
  	w_ILD=reshape(w_ILD,n_fgrp,n_taps_l);  

	w_ITD=w_r(:,counter:counter+n_taps_t-1); counter=counter+n_taps_t;
        w_ITD=reshape(w_ITD,n_fgrp,n_taps_t);  

%        [max_val,pos_del]=max(w_ITD');

	w_LL =w_r(:,counter); counter=counter+1;
	w_LR =w_r(:,counter); counter=counter+1;


cont=0;
res=5;
figure
colormap('jet')
	a_ILD=simulin(P_ILD,reshape(w_ILD,1,n_fgrp*n_taps_l),b);
subplot(2,2,1)
	flatplot(a_ILD,MAP,cont,1,res,1,0,0,0,['Simulation a_ILD, '],1);

	a_ITD=simulin(P_ITD,reshape(w_ITD,1,n_fgrp*n_taps_t),b);
subplot(2,2,2)
	flatplot(a_ITD,MAP,cont,1,res,1,0,0,0,['Simulation a_ITD, '],1);

	a_LL=simulin(LL_N,reshape(w_LL,1,n_fgrp),b);
subplot(2,2,3)
	flatplot(a_LL,MAP,cont,1,res,1,0,0,0,['Simulation a_LR, '],1);

	a_LR=simulin(LR_N,reshape(w_LR,1,n_fgrp),b);
subplot(2,2,4)
	flatplot(a_LR,MAP,cont,1,res,1,0,0,0,['Simulation a_LR, '],1);


%	a_LL=simulin(P,W,b);
%	a_LR=simulin(P,W,b);



%	f=figure('NumberTitle','off','Name',fname);
	f=figure;
	set(f,'Position',[50 50 700 500]);
cont=0;
res=3;

%	subplot(3,2,1); 
%	axis off;
%	text(0.1,1.0,['Name: ' respfname]);
%	text(0.1,0.5,['MSE: ' num2str(msqe)]);
	
%	drawnow
%	subplot(3,2,2);
%	axis off

%	text(0,1.0, ['ILD:']), 		text(0.5,1.0,int2str(cond_vec(1)));
%	text(0,0.8, ['ITD:']), 		text(0.5,0.8,int2str(cond_vec(2)));
%	text(0,0.6, ['Left:']),		text(0.5,0.6,int2str(cond_vec(3)));
%	text(0,0.4, ['Right:']),	text(0.5,0.4,int2str(cond_vec(4)));
%	text(0,0.2, ['Left norm:']),	text(0.5,0.2,int2str(cond_vec(5)));
%	text(0,0.0, ['Right norm:']),	text(0.5,0.0,int2str(cond_vec(6)));
  		
	drawnow	

	subplot(2,2,1);
colormap('jet')
	flatplot(T,MAP,cont,1,res,1,0,0,0,['Single Unit Response T, Max: ' num2str(eleT) ' Ele ' num2str(aziT) ' Azi'],1,-10,20);
	% text(100,-37.5,['Max: ' num2str(eleT) ' Ele ' num2str(aziT) ' Azi']);
	
        subplot(2,2,2);
	flatplot(a,MAP,cont,1,res,1,0,0,0,['Simulation a, Max: ' num2str(elea) ' Ele ' num2str(azia) ' Azi'],1,-10,20);
	%text(100,-37.5,['Max: ' num2str(elea) ' Ele ' num2str(azia) ' Azi']);

	subplot(2,2,3);
	flatplot(abs(e),MAP,cont,1,res,1,0,0,0,'Error T-a',1,-10,20);
    	drawnow; 

	figure
	
	subplot(2,1,1)
	mesh(w_ITD);

	subplot(2,1,2)
	mesh(w_ILD);
figure
	subplot(2,1,1)
	plot(w_LL);

	subplot(2,1,2)
	plot(w_LR);


%	subplot(3,2,6);
%	plot(b_array(fg:lg),w_r);
%	label_str=zeros(1,1);
%	for k=fg:lg label=int2str(bark2hz(b_array(k))); label_str(k-fg+1,1:length(label))=label; end;	
%	label_str=label_str(1:4:n_fgrp,:);
%	h=gca;
%	set(h,'XTick',b_array(fg:4:lg),'XTickLabels',label_str);

%	legstr=[];
% 	legp=1;
%	legstr(legp,1:3)='ILD'; legp=legp+1;
%	legstr(legp,1:3)='ITD'; legp=legp+1;
%	legstr(legp,1:4)='Left'; legp=legp+1;
%	legstr(legp,1:5)='Right';legp=legp+1;
%	if (x(cond,5)==1) legstr(legp,1:9)='Left norm';legp=legp+1;end;
%	if (x(cond,6)==1) legstr(legp,1:10)='Right norm';legp=legp+1;end;

%        legp=legp-1;
%	
%	if (legp==1) legend(legstr(1,:)); end
%	if (legp==2) legend(legstr(1,:),legstr(2,:)); end
%	if (legp==3) legend(legstr(1,:),legstr(2,:),legstr(3,:)); end
%	if (legp==4) legend(legstr(1,:),legstr(2,:),legstr(3,:),legstr(4,:)); end
   
%	ymax=max(max(w_r));     ymin=min(min(w_r));
%	xmin=min(min(b_array(fg:lg))); xmax=max(max(b_array(fg:lg)));
%	axis([xmin xmax ymin ymax]);
%	xlabel('Frequency [Hz]');
%	ylabel('Weighting Factor');
%	drawnow;

%	action=['save ' direc '/' fname ' a W b e P T cond_vec cond ']
%	eval(action);
%	action=['print -depsc '  direc '/' pname]
%	c=computer;
%	if strcmp(c,'PCWIN') 
%		action=strrep(action,'/','\')
%	end
%	eval(action);

%	close(f);
end
