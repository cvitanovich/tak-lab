function up2lo(listname)

list=readlist(listname);
[m,n]=size(list);

for k=1:m
	o=list(k,:);
	l=lower(o);
	unix(['mv ' o ' ' l]);
end
