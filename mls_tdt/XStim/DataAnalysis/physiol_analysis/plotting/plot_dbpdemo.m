function [] = plot_dbpdemo(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r);
      %a = ildcurve;
      %b = Neuron{cell_num}.tif_freqaxis;
      %c = Neuron{cell_num}.tif_ildaxis;
      %d = Neuron{cell_num}.tif_meansurf;
      %e = bp_optimal_ildf_surface;
      %f = space_contour;
      %g = IA_dbpmodel_azi;
      %h = IA_dbpmodel_ele;
      %i = ILD_matrix_new(:,IA_bestloc_ind);
      %j = Neuron{cell_num}.ia_meansurf{1};
      %k = Neuron{cell_num}.ia_diamond{1};
      %l = IA_dbpmodel_meansurf;
      %m = IA_dbpmodel_diam;
      %n = cell_num;
      %o = Neuron{cell_num}.tia_meanarray;
      %p = Neuron{cell_num}.tia_diamond;
      %q = Neuron{cell_num}.tia_azi;
      %r = Neuron{cell_num}.tia_ele;
      
if(size(f,1) == 1 & size(f,2) == 1)
   f = [0
      0];
end
cmap = 1-gray;
color_iabest = 'yellow';
linew_iabest = 2;      
fig_ht = 9.5;
fig_wd = 7.5;
mainfig = figure;
set(mainfig,'Units','inches');
set(mainfig,'Position',[0 0 fig_wd fig_ht]);

h1 = subplot('Position',[0.07 0.7 0.4 0.25]);
hold on
set(h1,'Units','normalized','FontSize',8);
axis([min(b) max(b) min(c) max(c)]);
plotsurf(b,c,d');
colormap(cmap);
plot(b,i,color_iabest,'LineWidth',linew_iabest);
xlabel('Frequency (Hz)')
ylabel('ILD (dB)')
title(['Measured ILD/freq, cell # ' num2str(n)])
text(b(2),i(2),'\downarrow','VerticalAlignment','Bottom','Color','yellow');
text(b(2),i(2)+4,'Best loc in ILDAlone RS','VerticalAlignment','Bottom','Color','yellow','FontSize',14);

h2 = subplot('Position',[0.55 0.7 0.4 0.25]);
hold on
set(h2,'Units','normalized','FontSize',8);
axis([min(b) max(b) min(c) max(c)]);
plotsurf(b,c,e');
colormap(cmap);
plot(b,i,color_iabest,'LineWidth',linew_iabest);
colorbar
temp = corrcoef(d(:),e(:)); vari_explained = temp(1,2)^2;
xlabel('Frequency (Hz)')
ylabel('ILD (dB)')
title(['Modeled ILD/freq, r^{2} = ' num2str(vari_explained)])

d = d/max(max(d)); e = e/max(max(e));
err_ildf = d - e;
h3 = subplot('Position',[0.07 0.35 0.4 0.25]);
hold on
set(h3,'Units','normalized','FontSize',8);
axis([min(b) max(b) min(c) max(c)]);
plotsurf(b,c,err_ildf');
colormap(cmap);
colorbar
plot(b,i,color_iabest,'LineWidth',linew_iabest);
temp = corrcoef(d(:),e(:)); vari_explained = temp(1,2)^2;
xlabel('Frequency (Hz)')
ylabel('ILD (dB)')
title('Normalized Measured - Modeled, ILD/freq')

h4 = subplot('Position',[0.08 0.05 0.25 0.25]);
hold on
axis square
set(h4,'Units','normalized','FontSize',8);
plotsurf(g,h,k);
colormap(cmap);
xlabel('Azimuth (\circ)')
ylabel('Elevation (\circ)')
title('Measured ILDAlone RS')
plot(f(1,:),f(2,:),'white','LineWidth',2);

h5 = subplot('Position',[0.4 0.05 0.25 0.25]);
hold on
axis square
set(h5,'Units','normalized','FontSize',8);
plotsurf(q,r,p);
colormap(cmap);
temp = corrcoef(j,o); vari_explained = temp(1,2)^2;
xlabel('Azimuth (\circ)')
ylabel('Elevation (\circ)')
title(str2mat('Line-integral Modeled ILDAlone RS',['r^{2} = ' num2str(vari_explained)]));
plot(f(1,:),f(2,:),'white','LineWidth',2);

h6 = subplot('Position',[0.73 0.05 0.25 0.25]);
hold on
axis square
set(h6,'Units','normalized','FontSize',8);
plotsurf(g,h,m);
colormap(cmap);
temp = corrcoef(j,l); vari_explained = temp(1,2)^2;
xlabel('Azimuth (\circ)')
ylabel('Elevation (\circ)')
title(str2mat('Euston BkPrp Modeled ILDAlone RS',['r^{2} = ' num2str(vari_explained)]));
plot(f(1,:),f(2,:),'white','LineWidth',2);
