function mtlproc(listname)

%function mtlproc(listname)

l_w = 0.003;
n_samples=256;
l=readlist(listname);

filename=strrep(listname,'.lst','.cx')
mtlmerge(l,filename);

mtlcorhe(filename);

fadename=strrep(filename,'.cx','.cxf')
mtlfade(filename,fadename);

irname=strrep(fadename,'.cxf','.rx')
mtlifft(fadename,irname);

cname=strrep(irname,'.rx','.rxc')
mtlcut(irname,cname,l_w, n_samples);

cleancatname=strrep(cname,'.rxc','.cxx')
mtlfft(cname,cleancatname);
