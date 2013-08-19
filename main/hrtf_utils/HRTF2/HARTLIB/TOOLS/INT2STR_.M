function str = int2str_(a)


b=sprintf('%3i',a);
str=strrep(b,' ','_')
return;
end;
