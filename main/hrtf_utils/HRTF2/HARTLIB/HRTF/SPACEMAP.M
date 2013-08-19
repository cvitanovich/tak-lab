function spacemap(ILD,ITD,LEVEL,T,MAP,b_array)

% SPACEMAP derives weighting factors for ILD, ITD and monaural spectrum (LEVEL)
%
% ILD interaural level difference
% ITD interaural time difference
% LEVEL monaural log mag spectrum
% T response of neuron
% MAP sperical coordinates for ILD ITD LEVEL and T
% b_array frequency groups


fg=13; % first frequency group
lg=23; % last  frequency group

n_taps_t = 5;
n_taps_l = 5;

delta_t= 100;
delta_l= 4;

[n_fgrp,ndir] = size(ILD);


n_fgrp=lg-fg+1;


% place coding of time difference

[N_ITD,new_delay]=mso(ITD,n_taps_t,0,delta_t);


% place coding of time difference

[N_ITD,new_delay]=mso(ILD,n_taps_l,0,delta_l);


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


% Calculation of P

	P=[]; n_param=0;
	P=[P;N_ILD(:,:)]; n_param=n_param+n_taps;
        P=[P;N_ITD(:,:)]; n_param=n_param+n_taps;
        P=[P;   LL(:,:)]; n_param=n_param+1;
        P=[P;   LR(:,:)]; n_param=n_param+1;


% Calculation of weighting factors W,b and test of modell a

	[W,b]=solvelin(P,T);
	a=simulin(P,W,b);


% Evaluation of fit

	[va,pa]=max(a); elea=MAP(1,pa); azia=MAP(2,pa);
	
	e=T-a;
	me=mean(e)
        msqe=mean(e.^2)

% Extraction of weighting factors for ILD, ITD, LL and LR

	w_r = reshape(W,n_fgrp,n_param);        
	counter=1;
	w_ILD=w_r(:,counter:counter+n_taps-1); counter=counter+n_taps_l;
        w_ITD=reshape(w_ILD,n_fgrp,n_taps);  
	w_ITD=w_r(:,counter:counter+n_taps-1); counter=counter+n_taps;
        w_ITD=reshape(w_ITD,n_fgrp,n_taps);  
%        [max_val,pos_del]=max(w_ITD');
	w_LL =w_r(:,counter); counter=counter+1;
	w_LR =w_r(:,counter); counter=counter+1;


% plotting results

	set(f,'Position',[50 50 700 500]);
	subplot(3,2,1); 
	axis off;
	text(0.1,1.0,['Name: ' respfname]);
	text(0.1,0.5,['MSE: ' num2str(msqe)]);
	
	drawnow
	subplot(3,2,2);
	axis off

%	text(0,1.0, ['ILD:']), 		text(0.5,1.0,int2str(cond_vec(1)));
%	text(0,0.8, ['ITD:']), 		text(0.5,0.8,int2str(cond_vec(2)));
%	text(0,0.6, ['Left:']),		text(0.5,0.6,int2str(cond_vec(3)));
%	text(0,0.4, ['Right:']),	text(0.5,0.4,int2str(cond_vec(4)));
%	text(0,0.2, ['Left norm:']),	text(0.5,0.2,int2str(cond_vec(5)));
%	text(0,0.0, ['Right norm:']),	text(0.5,0.0,int2str(cond_vec(6)));
  		
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

% plot weighting factors
	figure(2)

	subplot(2,2,1)
	surf(w_ITD);
	subplot(2,2,2)
	surf(w_ILD);
	subplot(2,2,3)
	plot(w_LL);
	subplot(2,2,4);
	plot(w_LR);
