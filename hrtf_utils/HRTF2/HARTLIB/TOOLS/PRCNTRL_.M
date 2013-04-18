function prcntrl_(mode)

if mode==1
  [fname,fpath]=uiputfile('*.eps','Print to file');
  if fname~=0 
    s=sprintf('print -depsc %s%s',fpath,fname);
    eval(s);
  end;
elseif mode==2 
 [fname,fpath]=uiputfile('*.djp','Print to file');
  if fname~=0 
    s=sprintf('print -dcdeskjet %s%s',fpath,fname);
    eval(s);
  end;
elseif mode==3
 [fname,fpath]=uiputfile('*.hjp','Print to file');
  if fname~=0 
    s=sprintf('print -dcdjcolor %s%s',fpath,fname);
    eval(s);
  end;
end;

