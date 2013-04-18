function rename(listname,s2,s3)

list=readlist(listname);
[m,n]=size(list);

for k=1:m
	o=list(k,:);
        s=strrep(o,s2,s3);
	unix(['mv ' o ' ' s]);
end
