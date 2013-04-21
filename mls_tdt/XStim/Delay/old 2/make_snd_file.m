function make_snd_file
    % make noise for preliminary experiments
    stim = MakeBBNoise(30000,150);
    [fileStr,pathStr,FilterIndex] = uiputfile('*.*','save stim');
	if(FilterIndex==0), return, end
	filepathStr=[pathStr fileStr];
    save(filepathStr, 'stim');