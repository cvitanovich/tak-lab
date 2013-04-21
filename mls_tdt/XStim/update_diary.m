% update_diary
% added June 2007
% called by each engage file
% create tempstr to write to diary identifying each test

FN.diary = [FN.data_path FN.data(1:4) 'diary.mat'];
if exist1(FN.diary)    
    load(FN.diary);
end

if exist1('STR')
    STRlen = length(STR);
else
    STRlen = 0;
end

temp = clock;
strmin = num2str(temp(5));
if length(strmin)==1    strmin = ['0' strmin];  end
tempstr0 = ['  '  num2str(temp(4)) ':' strmin '  '];
STR{STRlen+1} = [tempstr0 FN.data '      ' tempstr];

save(FN.diary,'STR')
clear tempstr* STRlen STR