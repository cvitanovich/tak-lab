function m=medianp(MAP)

front=find_azi(MAP,0);
[se,index]=sort(front);
front=front(index);

rear=find_azi(MAP,180);
d90=mtlfd(MAP,180,90);
k=find(rear ~= d90);
rear=rear(k)
[se,index]=sort(rear);
rear=rear(index)
rear=fliplr(rear)

m=[front rear];

end
