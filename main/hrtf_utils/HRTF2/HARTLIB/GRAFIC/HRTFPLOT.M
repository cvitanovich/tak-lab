function hrtfplot(index)

global MAP HRTF HRIR filename filetype Fs channels_per_dir fmin fmax lmin lmax

size(HRTF)
size(HRIR)

if (length(index) ~=1) return; end;
 
if (filetype==2)
	for k=1:channels_per_dir
		HRIR(:,k)=mtlrch(filename,channels_per_dir*index-(channels_per_dir-k));
	end
	HRTF=rfft(HRIR,4096);
end;
if (filetype==3);
	for k=1:channels_per_dir
		HRTF(:,k)=mtlrch(filename,channels_per_dir*index-(channels_per_dir-k));	
	end
	[ns,nc]=size(HRTF);
        HRIR=rifft(HRTF);
	if (ns < 4096) HRTF=rfft(HRIR,4096);end
end

ele=round(MAP(1,index));
azi=round(MAP(2,index));

g=figure(2)
set(g,'Position',[0,0,650,400]);
subplot(1,1,1);
magplot(HRTF(:,1:2),Fs,['Transfer Function,  Ele.: ' int2str(ele) '  Azi.: ' int2str(azi)],'gr',fmin,fmax,1);
set(gca,'YLim',[lmin,lmax]);
legend('Left','Right');
%subplot(2,1,2);
%magplot(HRTF(:,2),Fs,'HRTF-right','r',100,15500,1);
%phasplot(HRTF(:,:index*2),Fs);
%phdeplot(HRTF(:,index*2-1:index*2),Fs);
%tgrpplot(HRTF(:,index*2-1:index*2),Fs);
%subplot(2,2,2);
%irplot(HRIR(:,1),Fs,'IR-left','b',0,5,1);
%subplot(2,2,4);
%irplot(HRIR(:,2),Fs,'IR-right','b',0,5,1);

end;