function mtlcont(data,azi,ele,n)

% plots contour lines on sphere of mtlplots

azi=azi/180*pi;
ele=ele/180*pi;

dist=zeros(1,length(azi));

max_data=max(data);
min_data=min(data);

delta_data=max_data-min_data;

step_width=delta_data/(n+1);

for i=1:(n+1)
	index=find(data == min_data*(i-1)*step_width);
	azi_i=azi(index); 
	ele_i=ele(index);
	dist=sphdist([0,0],[ele_i,azi_i]);
	[min_dist,start]=min(dist);
	

end
