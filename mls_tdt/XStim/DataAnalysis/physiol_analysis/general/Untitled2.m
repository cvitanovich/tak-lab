freq_num = 26;
colors = 'bkgrcmybkgrcmybkgrcmybkgrcmy';
d_range = [5 20];


figure
hold
plot(dist_ild_s_sort(freq_num,:),bbif_meansurf(freq_num,dsort_s_index(freq_num,:)),'r')
plot(dist_ild_signed(freq_num,:),IA_meansurf,'.')
figure
hold
plot(dist_ild_s_sort(freq_num,:),bbif_meansurf(freq_num,dsort_s_index(freq_num,:)),'r')
plot(new_ildaxis,(new_bbif_meansurf(freq_num,:)),'b')
plot(new_ildaxis2,(new_bbif_meansurf2(freq_num,:))','g')
figure
hold
plot(dist_ild_s_sort(freq_num,:),resia_allfreqs(freq_num,:),'go')
[temp,res] = locregress1(dist_ild_s_sort(freq_num,:),...
   resia_allfreqs(freq_num,:),alpha);
plot(dist_ild_s_sort(freq_num,:),temp,'k')

index = find(dist_ild_signed(freq_num,:) >= d_range(1) &...
   dist_ild_signed(freq_num,:) <= d_range(2));
figure; hold
for ind_num = 1:length(index)
   plot(HRTFinfo.hrtf_freqs,ILD_matrix_focus(:,index(ind_num)),colors(ind_num));
   ht = text(HRTFinfo.hrtf_freqs(205),ILD_matrix_focus(205,index(ind_num)),num2str(IA_meansurf(index(ind_num))));
   set(ht,'Color',colors(ind_num),'FontWeight','bold');
   ht = text(HRTFinfo.hrtf_freqs(688),ILD_matrix_focus(688,index(ind_num)),num2str(IA_meansurf(index(ind_num))));
   set(ht,'Color',colors(ind_num),'FontWeight','bold');
end
hp = plot(HRTFinfo.hrtf_freqs,ILD_matrix_focus(:,IA_bestloc_ind),'k--');
set(hp,'LineWidth',2);
xlim([2000 11000])

