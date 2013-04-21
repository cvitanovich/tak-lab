%Test_type: Callback for Test Type popupmenu

val = get(H.test_type,'Value');

test_type = get(H.test_type,'String');
stim_val = get(H.test_type,'Value');
XStimParams.test_type = deblank(test_type(stim_val,:));

%   'String', ['None|Search for Cells|ABL|ITD|ILD|FREQ|McSpace|Mc2Source|Space|' ...
%'Space2|Space3|2Source|AltIR|AM|ITD2|ITD_decorr|Composite|BMLD|LRsounds|Adapt1|Adapt2|Space4'],...

disp('HUH????????')
switch val
    case 1,     % Search for cells
        searchcells;
    case 2,     % ABL Test
        abltest;
    case 3,     % ITD Test
        itdtest;
    case 4,     % ILD Test
        ildtest;
    case 5,     % FREQ Test
        freqtest;
    case 6,     % McSpace
        McSpace;
    case 7
        Mc2Source;
    case 8,     % Space or ILDAlone
        space;
    case 9,     % space2
        space2;
    case 10, 	 % Space3 or ILDAlone3
        space3;
    case 11		 % basic 2 source test
        Two_Source;
    case 12		 % alterred IR (Hanna's)
        alt_IR_XXX;  
    case 13
        AM;
    case 14
        ITD2;
    case 15
        ITD_decorr;
    case 16
        Composite; 
    case 17
        BMLD;
    case 18
        LRsounds;
    case 19
        Adapt1;
    case 20
        Adapt2;
    case 21
        Space4;
    case 24             % plays nReps of each file chosen
        Space_FilePerRep;
    case 25
        Long_SAMs;
    case 26
        MaskedSpace;
    otherwise,
        disp('Not implemented')
        set(H.test_type,'Value',1);
end
eval(['save ' FN.current_path 'XStimParams_current XStimParams']);