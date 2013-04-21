function MakePEstim()
% MakePEstim: creates PEStim files in case rand function changes
    % /////////////////////////////////////
    % set SaveFile to 1 in function PEstim(...
    %/////////////////////////////////////
    %PEStimSeeds(); % make seeds once
    %return;
    
    RandStates = [...
    73 167 617 364 662 593 538 194 853 610 294 ...
    479 71 105 162 770 143 116 252 101 377 706 ...
    273 574 915 661 935 355 530 540 220 232 886 ...
    70 65 571 35 339 87 281 795 283 974 248 ...
    995 936 769 943 127 224];

    for s=1:size(RandStates,2)
        PEstim(RandStates(s),300,300,0,0,0,30000); % 300 ms noises
    end
    
end

    
function PEStimSeeds()
    unique=1; % unique numbers
    a = 1; b = 999;
    for s=1:50
        seeds(1,s)=round(a + (b-a) * rand(1));
    end
    for s=1:50
        for k=1:50
            if seeds(s)== seeds(k) & s~= k
               unique=0
               disp(seeds(k));
            end
        end
    end
    if(unique)
        disp(seeds);
    end
end