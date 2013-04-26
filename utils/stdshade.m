function stdshade(amatrix,alpha,acolor,Xes,smoothval)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23
if exist('acolor','var')==0 || isempty(acolor)
    acolor=['r' 'k']; 
end

if exist('F','var')==0 || isempty(Xes); 
    Xes=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smoothval); smoothval=1; end
else smoothval=1;
end

if ne(size(Xes,1),1)
    Xes=Xes';
end

amean=smooth(nanmean(amatrix),smoothval)';
astd=nanstd(amatrix); % to get std shading
% astd=nanstd(amatrix)/sqrt(size(amatrix,1)); % to get sem shading
if exist('alpha','var')==0 || isempty(alpha) 
    fill([Xes fliplr(Xes)],[amean+astd fliplr(amean-astd)],acolor(1),'linestyle','none');
    acolor='k';
else fill([Xes fliplr(Xes)],[amean+astd fliplr(amean-astd)],acolor(1), 'FaceAlpha', alpha,'linestyle','none');    
end

if ishold==0
    check=true; else check=false;
end

hold on; plot(Xes,amean,acolor(2),'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end



