% MAKE_DATA prepares a series of tests for later analysis.  Stores data in rapid 
%           retrieval format.
% all data format:
% bird side(0=L) site cluster itd ild-freq true_space itd_abi_space abi_space
% the variable whichtests specifies which tests in the all_data structure are spatially tuned
% cells based on their spatial tuning curves


all_data = [
%           itd	if		space		iaspace	aspace
[890 1 4  1 47		48		46			0			0]; % probably tectum, very deep cell
[890 1 15 1 119	122	118		0			0]; % too shallow to be IC, found before tonotopy started
[890 1 18 1 145	150	149		0			0];
[890 1 18 2 145	150	149		0			0];
[890 1 19 1 154	161	160		0			0];
[890 1 22 1 166	168	169		0			0]; % shallow, before tonotopy (suspect)
[890 1 24 1 175	179	181		0			0];
[890 0 3  1 6		14		15			0			0];
[890 0 3  2 28		20		22			0			0];
[890 0 6  1 45		59		58			0			0];
[890 0 8  1 65		70		69			0			0];
[890 0 9  1 81		84		89			0			0];
[890 0 12 1 97		101	98			0			0];
[890 0 17 1 118	132	120		0			0];
[890 0 18 1 133	136	137		0			0];
[890 0 20 1 139	147	146		0			0];
[890 0 21 1 149	152	150		0			0];
[891 1 2  1 8		12		13			0			0];
[891 1 3  1 14		18		15			0			0];
[891 1 5  1 26		32		31			0			0];
[891 1 6  1 33		40		34			0			0];
[891 1 7  2 43		51		47			0			0];
[891 1 7  3 48		49		47			0			0];
[891 1 9  1 53		56		54			0			0];
[891 1 13 1 75		80		77			0			0];
[891 1 13 2 75		80		77			0			0];
[891 1 14 1	86		93		91			0			0];
[891 1 15 1 101	98		100		0			0];
[891 1 16 1 116	121	117		0			0];
[891 1 17 1 129	131	130		0			0];
[891 1 18 1 132	134	133		0			0];
[891 1 21 1	142	148	149		145		143];
[891 1 27 1	165	171 	166		169		167]; 
[891 1 27 2	172	171 	174		0			173]; 
[891 1 28 1	175	178 	176		177		0]; 
[891 1 30 1 182	188	183		187		0];
[891 1 31 1	189	193	191		190		0];
[891 1 32 1 195	202	200		0			0];
[891 1 33 1 204	207	205		209		0];
[891 1 34 1 210	215	211		214		0];
[891 1 35 1 217	224	220		225		0];
[891 1 39 1 234	242	235		238		239];
[891 1 43 1 248	250	249		0			0];
[891 1 44 1 252	253	254		0			0];
[891 0 1  1	1		4		2			0			0];
[891 0 3  1	14		16		15			0			0]
[891 0 6  1 24		23		22			0			0];
[891 0 8  1 28		31		29			0			0];
[891 0 8  1 0     44		45			0			0];
[891 0 10 1 49		54		50			0			0];
[891 0 12 1 57		60		59			0			0];
[891 0 15 1 65		67		66			0			0];
[891 0 16 1 71		74		72			0			0];
[891 0 16 2 73		74		72			0			0];
[891 0 17 1	79		83		78			0			0]; % ITD tuning lame
[891 0 18 1 84		86		85			0			0]; % this cells is a high-freq extreme ILD responder
[891 0 19 1 88		90		94			0			0];
[891 0 20 1	97		98		95			0			0];
[891 0 21 1	103	105	104		0			0];
[891 0 23 1	116	118	114		0			0];
[891 0 26 1 136	141	137		0			0];
[891 0 26 1 132	135	133		0			0];
[891 0 30 1 165	166	164		0			0];
[891 0 32 1	171	173	170		172		0];  % interesting cell, demonstrates importance of loss of itd tuning for extreme ILD's
[891 0 33 1	177	182	179		183		0];
[891 0 35 1	189	192	190		0			0];
[891 0 36 1 194	201	196		200		0];
[891 0 41 1	209	215	211		212		0];
[891 0 45 1 229	232	231		230		0];
[891 0 46 1 235	238	236		237		0];
[891 0 47 1 239	244	241		242		0];
[891 0 48 1 259	257	258		260		0];
[891 0 49 1 262	264	261		263		0];
[897 1 1  1 1		7		6			9			0];
[897 1 3  1 13		17		15			18			0];
[897 1 7  1 26		32		27			28			29];
[897 1 8  1 33		37		34			36			0];
[897 1 10 1 44		50		42			43			48];
[897 1 10 2 44		50		42			43			48];
[897 1 11 1 53		57		54			55			0];
[897 1 12 1 58		62		59			60			0];
[897 1 15 1 77		81		88			0			0];
[897 1 16 1 93		97		96			0			0];
[897 1 18 1 109	113	110		111		0];
[897 1 19 1 117	120	121		0			0];
[897 1 22 1 126	130	127		129		0];
[897 1 22 2 132	134	131		133		135];
[897 1 23 1 141	145	142		143		149];
[897 1 24 1 154	159	157		158		0];
[897 1 26 1 163	166	168		165		167];
[897 1 27 1 175	180	177		178		0];
[897 1 35 1 198	201	199		200		202];
[897 1 37 1 208	212	210		214		0];
[897 1 38 1 215	218	216		217		0];
[897 1 41 1 225	228	226		227		0];
[897 1 48 1 242	249	252		246		251];
[897 1 49 1 253	258	254		257		0];
[897 1 51 1 263	268	264		266		0];
[897 1 51 2 263	268	264		266		0];
[897 1 52 1 271	273	277		272		278];
[897 1 53 1 285	288	280		287		0];
[897 1 54 1 292	294	290		293		0];
[897 0 1  1 1		3		4			0			0];
[897 0 3  1 16		21		20			0			0];
[897 0 6  1 26		28		27			0			0];
[897 0 8  1 31		34		32			35			0];
[897 0 10 1 38		52		49			40			43];
[897 0 12 1 58		63		60			62			67];
[897 0 13 1 72		73		71			74			0];
[897 0 14 1 77		79		76			78			0];
[897 0 15 1 82		85		86			84			87];
[897 0 16 1 88		91		90			92			0];
[897 0 18 1 95		98		96			97			100];
[897 0 19 1 101	104  102			103		0];
[897 0 20 1 110	109  107			108		0];
[897 0 21 1 112	115  113			114		116];
[897 0 23 1 123	130  125			126		131];
[897 0 25 1 136	140  138			139		141];
[897 0 26 1 150	154  152			153		155];
[897 0 27 1 163	166  164			165		0];
[897 0 27 2 163	166  164			165		0]
[897 0 28 1 171	173  169			172		0];
[897 0 30 1 177	181  179			180		185];
[897 0 31 1 189	192  190			191		193]
[897 0 36 1 210	218  211			212		215];
[897 0 38 1 229	226  223			225		227];
[897 0 39 1 240	242  238			241		0];
[897 0 40 1 257	251  256			250		253];
[897 0 43 1 270	271  269			272		0];
[897 0 44 1 274	277  273			276		0]; 
];

whichtests = [     
     1
     8
    20
    32
    34
    35
    48
    49
    53
    60
    63
    67
    68
    72
    73
    76
    77
    78
    79
    80
    81
    87
    88
    89
    90
    91
    92
    95
    96
   100
   104
   108
   110
   112
   114
   116
   117
   118
   119
   120
   121
   123
   125
   126
   127
   128
   130
];


for i = 1:size(all_data,1)
%for i = 76:size(all_data,1)
 
   bird = all_data(i,1);
   if (all_data(i,2)==0) side = 'L';
   else side = 'R';
   end;
   cl = all_data(i,4);
   itd_test = all_data(i,5);
   sptest = all_data(i,7:9);
   for j = 1:3
	   if (sptest(j)~=0)
	      [mean_space_surf std_resp_surf locs dim1 space_tpars{i,j} spont_count_space(i,j) spont_dur_space(i,j) reps_space(i,j) sr_space(i,j)] = proc_test(bird,side,sptest(j),cl,0);
			
			% interpolate
			[a, e, data] = array2diamond(mean_space_surf, locs');	
			[a, e, std_space{i,j}] = array2diamond(std_resp_surf, locs');	
				
			[AZ EL] = meshgrid(a, e);
			
			% interpolate missing values
			
			% generate mask for missing points
			missmask = NaN*ones(size(data));
			
			for az = -90:5:90;
			  for el = -90+abs(az):5:90-abs(az)
			    if (~(locs(:,1)==el & locs(:,2)==az))
			       missmask(AZ==az & EL==el) = 1;
			    else
			       missmask(AZ==az & EL==el) = 0;
			    end; 
			  end;
			end;
			
			% replace all valid locations presently containing NaN's with zeros
			ri = (missmask==1);
			data(ri) = zeros(size(data(ri)));
			
			% generate interpolation function
			intfun = [[0 1 0]; [1 0 1]; [0 1 0]];
			intval = conv2(data, intfun,'same')/4;
			intval = intval.*missmask;
			intval(find(isnan(intval))) = zeros(size(find(isnan(intval))));
			
			spacedata = intval + data;	      
	      
	      azi{i} = a;
	      ele{i} = e;
	      space_test{i,j} = spacedata;
	      
	   end;  % if sptest
   end;
   if (itd_test~=0)
      [mean_itd_curve{i} std_itd{i} null itd_axis{i} itd_tpars{i} spont_count_itd(i) spont_dur_itd(i) reps_itd(i) sr_itd(i)] = proc_test(bird,side,itd_test,cl,0);
   end;
   ildfreq_test = all_data(i,6);
   if (ildfreq_test~=0)
      [ildfreq_surf{i} std_ildf{i} ild_axis{i} freq_axis{i} ildfpars{i} spont_count_ildf(i) spont_dur_ildf(i) reps_ildf(i) sr_ildf(i)] = ...
                      proc_test(bird,side,ildfreq_test,cl,0);
   end;
end;

save space_data all_data space_test azi ele mean_itd_curve itd_axis ...
                ildfreq_surf ild_axis freq_axis ildfpars...
                 space_tpars itd_tpars ...
                spont_count_space spont_count_itd spont_count_ildf ...
                spont_dur_space   spont_dur_itd   spont_dur_ildf ...
					 std_space std_itd  std_ildf  reps_space reps_itd reps_ildf ...
                sr_itd sr_ildf sr_space
                
disp('DATA SAVED, OK TO QUIT');
