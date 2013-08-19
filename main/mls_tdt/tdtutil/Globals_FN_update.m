% Globals_FN_update
% updates FN to current bird and checks for HRTFfiletype
% requires loading XStimParams and FN first

FN.space_eq = '938AD.eq.mat';
FN.ildalone_eq = '938AD.ila.mat';
FN.itdalone_eq = '938AD.ita.mat';
FN.space_std = '938AD.std.mat';
FN.ildalone_std = '938AD.ila.std.mat';
FN.itdalone_std = '938AD.ita.std.mat';
FN.loc = '938AD.eq.mat';
FN.ephone2 = '938ET_EC_inv.mat';
FN.ablequal_eq = '938AD_ABLequal.eq.mat';
FN.ablequal_std = '938AD_ABLequal.std.mat';

switch XStimParams.bird_number
    case 933             % note: incomplete here for ILA and ITA
        FN.space_eq = '933AD.eq.mat';
        FN.ildalone_eq = '938AD.ila.mat';
        FN.itdalone_eq = '938AD.ita.mat';
        FN.space_std = '933AD.std.mat';
        FN.ildalone_std = '938AD.ila.std.mat';
        FN.itdalone_std = '938AD.ita.std.mat';
        FN.loc = '933AD.eq.mat';
        FN.ephone2 = '933ET_EC_inv.mat';
    case 938
        FN.space_eq = '938AD.eq.mat';
        FN.ildalone_eq = '938AD.ila.mat';
        FN.itdalone_eq = '938AD.ita.mat';
        FN.space_std = '938AD.std.mat';
        FN.ildalone_std = '938AD.ila.std.mat';
        FN.itdalone_std = '938AD.ita.std.mat';
        FN.loc = '938AD.eq.mat';
        FN.ephone2 = '938ET_EC_inv.mat';
        FN.ablequal_eq = '938AD_ABLequal.eq.mat';
        FN.ablequal_std = '938AD_ABLequal.std.mat';
    case 942
        FN.space_eq = '942AD.eq.mat';
        %FN.ildalone_eq = '938AD.ila.mat';
        %FN.itdalone_eq = '938AD.ita.mat';
        FN.space_std = '942AD.std.mat';
        %FN.ildalone_std = '938AD.ila.std.mat';
        %FN.itdalone_std = '938AD.ita.std.mat';
        FN.loc = '942AD.eq.mat';
        FN.ephone2 = '942_ET_EC__inv.mat';
        FN.ablequal_eq = '942AD_ABLequal.eq.mat';
        FN.ablequal_std = '942AD_ABLequal.std.mat';
    case 943
        FN.space_eq = '943AD.eq.mat';
        %FN.ildalone_eq = '938AD.ila.mat';
        %FN.itdalone_eq = '938AD.ita.mat';
        FN.space_std = '943AD.std.mat';
        %FN.ildalone_std = '938AD.ila.std.mat';
        %FN.itdalone_std = '938AD.ita.std.mat';
        FN.loc = '943AD.eq.mat';
        FN.ephone2 = '943_ET_EC__inv.mat';
        FN.ablequal_eq = '943AD_ABLequal.eq.mat';
        FN.ablequal_std = '943AD_ABLequal.std.mat';
    case 945
        FN.space_eq = '945AD.eq.mat';
        FN.ildalone_eq = '945AD.ila.mat';
        FN.itdalone_eq = '945AD.ita.mat';
        FN.space_std = '945AD.std.mat';
        FN.ildalone_std = '945AD.ila.std.mat';
        FN.itdalone_std = '945AD.ita.std.mat';
        FN.loc = '945AD.eq.mat';    
        FN.ephone2 = '945ET_EC_inv.mat';
    case 966
        FN.space_eq = '966AD.eq.mat';
        FN.ildalone_eq = '966AD.ila.eq.mat';
        FN.itdalone_eq = '966AD.ita.eq.mat';
        FN.space_std = '966AD.std.mat';
        FN.ildalone_std = '966AD.ila.std.mat';
        FN.itdalone_std = '966AD.ita.std.mat';
        FN.loc = '966AD.eq.mat';    
        FN.ephone2 = '938ET_EC_inv.mat';
    case 1007
        FN.space_eq = '1007BD.eq.mat';
        %FN.ildalone_eq = '1007BD.ila.mat';
        %FN.itdalone_eq = '1007BD.ita.mat';
        FN.space_std = '1007BD.std.mat';
        FN.ildalone_std = '1007BD.ila.std.mat';
        FN.itdalone_std = '1007BD.ita.std.mat';        
        FN.loc = '1007BD.eq.mat';
        FN.ephone2 = '1007_ET_EC__inv.mat';
        FN.ablequal_eq = '1007BD_ablequal.eq.mat';
        FN.ablequal_std = '1007BD_ablequal.std.mat';
    case 1027
        FN.space_eq = '1027AD.eq.mat';
        FN.ildalone_eq = '1027AD.ila.eq.mat';
        FN.itdalone_eq = '1027AD.ita.eq.mat';
        FN.space_std = '1027AD.std.mat';
        FN.ildalone_std = '1027AD.ila.std.mat';
        FN.itdalone_std = '1027AD.ita.std.mat';        
        FN.loc = '1027AD.eq.mat';
        FN.ephone2 = '1027AD_ET_EC__inv.mat';
        FN.ablequal_eq = '1027AD_ablequal.eq.mat';
        FN.ablequal_std = '1027AD_ablequal.std.mat';
    case 1029
        FN.space_eq = '1029AD.eq.mat';
        FN.ildalone_eq = '1029AD.ila.eq.mat';
        FN.itdalone_eq = '1029AD.ita.eq.mat';
        FN.space_std = '1029AD.std.mat';
        FN.ildalone_std = '1029AD.ila.std.mat';
        FN.itdalone_std = '1029AD.ita.std.mat';        
        FN.loc = '1029AD.eq.mat';
        FN.ephone2 = '1029AD_ET_EC__inv.mat';
        FN.ablequal_eq = '1029AD_ablequal.eq.mat';
        FN.ablequal_std = '1029AD_ablequal.std.mat';
    case 1046
        FN.space_eq = '1046AC.eq.mat';
        FN.ildalone_eq = '1046AC_ila.eq.mat';
        FN.itdalone_eq = '1046AC_ita.eq.mat';
        FN.space_std = '1046AC.std.mat';
        FN.ildalone_std = '1046AC_ila.std.mat';
        FN.itdalone_std = '1046AC_ita.std.mat';        
        FN.loc = '1046AC.eq.mat';
        FN.ephone2 = '1046AC_ET_EC__inv.mat';
        FN.ablequal_eq = '1046AC_ablequal.eq.mat';
        FN.ablequal_std = '1046AC_ablequal.std.mat';
    case 1047
        FN.space_eq = '1047AD.eq.mat';
        FN.ildalone_eq = '1047AD_ila.eq.mat';
        FN.itdalone_eq = '1047AD_ita.eq.mat';
        FN.space_std = '1047AD.std.mat';
        FN.ildalone_std = '1047AD_ila.std.mat';
        FN.itdalone_std = '1047AD_ita.std.mat';        
        FN.loc = '1047AD.eq.mat';
        FN.ephone2 = '1047AD_ET_EC__inv.mat';
        FN.ablequal_eq = '1047AD_ablequal.eq.mat';
        FN.ablequal_std = '1047AD_ablequal.std.mat';
    case 1049
        %FN.space_eq = '1049AE.eq.mat';
        FN.ildalone_eq = '1049AE_ila.eq.mat';
        FN.itdalone_eq = '1049AE_ita.eq.mat';
        %FN.space_std = '1049AE.std.mat';
        FN.ildalone_std = '1049AE_ila.std.mat';
        FN.itdalone_std = '1049AE_ita.std.mat';        
        %FN.loc = '1049AE.eq.mat';
        %FN.ephone2 = '1049AE_ET_EC__inv.mat';
        FN.ablequal_eq = '1049AE_ablequal.eq.mat';
        FN.ablequal_std = '1049AE_ablequal.std.mat';
        
        FN.space_eq = '1049CDF.eq.mat';
        FN.space_std = '1049CDF.std.mat';
        FN.loc = '1049CDF.eq.mat';
        FN.ephone2 = '1049CDF_ET_EC__inv.mat';
        
    otherwise
        FN.space_eq = '938AD.eq.mat';
        FN.ildalone_eq = '938AD.ila.mat';
        FN.itdalone_eq = '938AD.ita.mat';
        FN.space_std = '938AD.std.mat';
        FN.ildalone_std = '938AD.ila.std.mat';
        FN.itdalone_std = '938AD.ita.std.mat';
        FN.loc = '938AD.eq.mat';
        FN.ephone2 = '938ET_EC_inv.mat';
end


% check for HRTFfiletypes (added Sept 02, 2004)
% filetypes:    0) untested
%               1) traditional binary
%               2) new *.mat with TRF1 and TF2
%               999) non-readable
% see also Globals_var.m
if ~isempty(FN.space_eq)
    FN.HRTFfiletype(1,1) = testHRTFfiletype(FN.space_path, FN.space_eq);
end
if ~isempty(FN.ildalone_eq)
    FN.HRTFfiletype(2,1) = testHRTFfiletype(FN.ILA_path, FN.ildalone_eq);
end
if ~isempty(FN.itdalone_eq)
    FN.HRTFfiletype(3,1) = testHRTFfiletype(FN.ITA_path, FN.itdalone_eq);
end
if ~isempty(FN.space_std)
    FN.HRTFfiletype(1,2) = testHRTFfiletype(FN.space_path, FN.space_std);
end
if ~isempty(FN.ildalone_std)
    FN.HRTFfiletype(2,2) = testHRTFfiletype(FN.ILA_path, FN.ildalone_std);
end
if ~isempty(FN.itdalone_std)
    FN.HRTFfiletype(3,2) = testHRTFfiletype(FN.ITA_path, FN.itdalone_std);
end
if ~isempty(FN.loc)
    FN.HRTFfiletype(4) = testHRTFfiletype(FN.loc_path, FN.loc);
end
if ~isempty(FN.ephone)
    FN.HRTFfiletype(5) = testHRTFfiletype(FN.ephone_path, FN.ephone);
end
if ~isempty(FN.ephone2)
    FN.HRTFfiletype(6) = testHRTFfiletype(FN.ephone_path, FN.ephone2);
end
if ~isempty(FN.ablequal_eq)
    FN.HRTFfiletype(7,1) = testHRTFfiletype(FN.space_path, FN.ablequal_eq);
end
if ~isempty(FN.ablequal_std)
    FN.HRTFfiletype(7,2) = testHRTFfiletype(FN.space_path, FN.ablequal_std);
end
