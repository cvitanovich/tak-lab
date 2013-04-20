if(0)

    p='C:\andrew\CORE\tak-lab\HRTFs\Matlab_V7';
    setenv('NEWPATH','C:\andrew\CORE\tak-lab\HRTFs\Matlab_V6\');
    
    s=what(p);
    
    nfiles=length(s.mat);
    
    for cnt=1:nfiles
        save tmp.mat
        clearvars -except s cnt
        setenv('FNAME',s.mat{cnt});
        load(getenv('FNAME'));
        clear s cnt
        save(['C:\andrew\CORE\tak-lab\HRTFs\Matlab_V6\' getenv('FNAME')],'-v6');
        load tmp.mat
    end
    
    delete('tmp.mat');

end