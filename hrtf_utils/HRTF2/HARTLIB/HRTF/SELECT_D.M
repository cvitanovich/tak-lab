function [e,a]=select_dir(select_flag)

global MAP map_stat min_ele max_ele min_azi max_azi min_res_azi min_res_ele n_ele n_azi

if (nargin <1)
   select_flag='exclusive';
end;

if (nargin==1)
 if ~isstr(select_flag) error('Choose exclusive, and'); end; 
end;

[l,c]=find(map_stat==30);

if select_flag=='exclusive'
    if (length(l) >0) map_stat(l,c)=15; end;
end



p=get(gca,'CurrentPoint');
a=round(p(1,1));
e=round(p(1,2));

l=round((e-min_ele)/min_res_ele)+1;
c=round((a-min_azi)/min_res_azi)+1;

e=(l-1)*min_res_ele+min_ele;
a=(c-1)*min_res_azi+min_azi;

if (map_stat(l,c)~=0)
	if (select_flag=='exclusive')   map_stat(l,c)=30; end;
	%if (select_flag=='compare')     map_stat(l,c)=45; end;
   	image([min_azi,max_azi],[max_ele,min_ele],map_stat);
else
	e=[]; a=[];
end;

return;

end;