function prcntrl_(mode)

if mode==1
  [fname,fpath]=uiputfile('*.ps','Print to file');
  if fname~=0 
    s=sprintf('print -dps %s%s',fpath,fname);
    eval(s);
  end;
elseif mode==2 | mode==3
  print -dwin -v
end;

ÿ