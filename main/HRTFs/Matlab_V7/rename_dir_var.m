if(0)
    s=what(pwd);
    nfiles=length(s.mat);
    
    for cnt=1:nfiles
        save tmp.mat
        clearvars -except s cnt
        setenv('FNAME',s.mat{cnt});
        load(getenv('FNAME'));
        clear s cnt
        if(exist('dir','var'))
            direc=dir;
            clear('dir')
        end
        save(getenv('FNAME'),'-v7');
        load tmp.mat
    end
    delete('tmp.mat')
end