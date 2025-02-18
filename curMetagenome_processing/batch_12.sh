#!/bin/bash
NXF_VER=23.10.1 ./nextflow run ikmb/tofu-maapo -r 6887e5c -profile custom -c tofu.config --sra '[SRR1622928,SRR1622930,SRR1622931,SRR1622932,SRR1622933,SRR1622934,SRR1622935,SRR1622936,SRR1622937,SRR1622938,SRR1623045,SRR1623046,SRR1623047,SRR1623048,SRR1623049,SRR1623050,SRR1623051,SRR1623065,SRR1623074,SRR1623075,SRR1623158,SRR1623169,SRR1623175,SRR1623176,SRR1623185,SRR1623194,SRR1623208,SRR1623209,SRR1623215,SRR1623224,SRR1624794,SRR1624795,SRR1624804,SRR1624810,SRR1624819,SRR1624836,SRR1624837,SRR1624846,SRR1624855,SRR1624864,SRR1630934,SRR1630951,SRR1630952,SRR1630961,SRR1630970,SRR1630985,SRR1630988,SRR1631005,SRR1631006,SRR1631017,SRR1631030,SRR1631036,SRR1631055,SRR1631056,SRR1631070,SRR1631086,SRR1631094,SRR1631102,SRR1631103,SRR1631109,SRR1631244,SRR1631250,SRR1631256,SRR1631262,SRR1631276,SRR1631285,SRR1631293,SRR1631300,SRR1631307,SRR1631308,SRR1631336,SRR1631353,SRR1631354,SRR1631360,SRR1631366,SRR1631372,SRR1631378,SRR1631385,SRR1631391,SRR1631402,SRR1631434,SRR1631448,SRR1631461,SRR1631462,SRR1631468,SRR1631474,SRR1631483,SRR1631484,SRR1631488,SRR1631494,SRR1631500,SRR1631511,SRR1631517,SRR1631518,SRR1631525,SRR1631536,SRR1631537,SRR1631551,SRR1631557,SRR1631561,SRR1631666,SRR1631667,SRR1631673,SRR1631685,SRR1631699,SRR1631705,SRR1631706,SRR1631710,SRR1631716,SRR1631722,SRR1632296,SRR1632302,SRR1632313,SRR1632320,SRR1632321,SRR1632328,SRR1632332,SRR1632339,SRR1632346,SRR1632352,SRR1632358,SRR1632365,SRR1632379,SRR1632388,SRR1632394,SRR1632400,SRR1632401,SRR1632405,SRR1632412,SRR1632422,SRR1632430,SRR1632431,SRR1632437,SRR1632456,SRR1632457,SRR1632468,SRR1632474,SRR1632475,SRR1632481,SRR1632485,SRR1632785,SRR1632789,SRR1632795,SRR1632815,SRR1632816,SRR1632822,SRR1632834,SRR1632835,SRR1632841,SRR1632848,SRR1632977,SRR1632988,SRR1632994,SRR1632995,SRR1633006,SRR1633010,SRR1633011,SRR1633020,SRR1633021,SRR1633033,SRR1633146,SRR1633147,SRR1633158,SRR1633162,SRR1633170,SRR1633177,SRR3726381,SRR3726380,SRR3726379,SRR3726378,SRR3726377,SRR3726376,SRR3726375,SRR3726374,SRR3726373,SRR3726372,SRR3726371,SRR3726370,SRR3726369,SRR3726368,SRR3726367,SRR3726366,SRR3726365,SRR3726364,SRR3726363,SRR3726362,SRR3726361,SRR3726360,SRR3726359,SRR3726358,SRR3726357,SRR3726356,SRR3726355,SRR3726354,SRR3726353,SRR3726352,SRR3726351,SRR3726350,SRR3726349,SRR3726348,SRR3726347,SRR3726346,SRR3726345,SRR3726344,SRR3726343,SRR3726342,SRR3726339,SRR3726338,SRR3726337,SRR3726336,SRR3726335,SRR7658632,SRR7658646,SRR7658687,SRR7658611,SRR7658584,SRR7658644,SRR7658609,SRR7658610,SRR7658690,SRR7658588,SRR7658596,SRR7658615,SRR7658653,SRR7658672,SRR7658667,SRR7658620,SRR7658688,SRR7658676,SRR7658580,SRR7658678,SRR7658585,SRR7658627,SRR7658657,SRR7658663,SRR7658633,SRR7658617,SRR7658677,SRR7658635,SRR7658682,SRR7658624,SRR7658686,SRR7658680,SRR7658606,SRR7658618,SRR7658589,SRR7658619,SRR7658668,SRR7658603,SRR7658642,SRR7658665,SRR7658581,SRR7658616,SRR7658599,SRR7658626,SRR7658684,SRR7658628,SRR7658683,SRR7658604,SRR7658651,SRR7658637,SRR7658583,SRR7658631,SRR7658661,SRR7658636,SRR7658634,SRR7658602,SRR7658622,SRR7658689,SRR7658674,SRR7658592,SRR7658600,SRR7658650,SRR7658685,SRR7658594,SRR7658605,SRR7658598,SRR7658614,SRR7658625,SRR7658607,SRR7658586,SRR7658664,SRR7658591,SRR7658655,SRR7658579,SRR7658645,SRR7658679,SRR7658681,SRR7658593,SRR7658660,SRR7658649,SRR7658590,SRR7658675,SRR7658613,SRR7658662,SRR7658652,SRR7658658,SRR7658608,SRR7658587,SRR7658601,SRR7658629,SRR7658597,SRR7658639,SRR7658666,SRR7658647,SRR7658638,SRR7658621,SRR7658612,SRR7658643,SRR7658630,SRR7658670,SRR7658659,SRR7658671,SRR7658654,SRR7658669,SRR7658623,SRR7658595,SRR7658641,SRR7658640,SRR7658582,SRR7658673,SRR7658656,SRR7658648,SRR2938595,SRR2938594,SRR2938551,SRR2938550,SRR2938549,SRR2938548,SRR2938547,SRR2938546,SRR2938545,SRR2938544,SRR2938543,SRR2938542,SRR2938541,SRR2938540,SRR2938539,SRR2938538,SRR2938537,SRR2938536,SRR2938535,SRR2938534,SRR2938533,SRR2938532,SRR2938531,SRR2938530,SRR2938529,SRR2938528,SRR2938527,SRR2938526,SRR2938525,SRR2938524,SRR2938523,SRR2938521,SRR2938520,SRR2938519,SRR2938518,SRR2938516,SRR2938515,SRR2938514,SRR2938513,SRR2938512,SRR2938511,SRR2938510,SRR2938509,SRR2938508,SRR2938507,SRR2938506,SRR2938501,SRR2938500,SRR2938499,SRR2938498,SRR2938497,SRR2938496,SRR2938495,SRR2938494,SRR2938493,SRR2938492,SRR2938491,SRR2938490,SRR2938489,SRR2938488,SRR2938487,SRR2938486,SRR2938485,SRR2938484,SRR2938483,SRR2938482,SRR2938481,SRR2938480,SRR2938479,SRR2938478,SRR2938477,SRR2938476,SRR2938475,SRR2938474,SRR2938473,SRR2938472,SRR2938471,SRR2938470,SRR2938469,SRR2938468,SRR2938467,SRR2938466,SRR2938465,SRR2938464,SRR2938463,SRR2938462,SRR2938461,SRR2938460,SRR2938459,SRR2938458,SRR2938457,SRR2938456,SRR2938455,SRR2938454,SRR2938453,SRR2938452,SRR2938451,SRR2938450,SRR2938449,SRR2938448,SRR2938447,SRR2938446,SRR2938445,SRR2938444,SRR2938443,SRR2938442,SRR2938441,SRR2938440,SRR2938439,SRR2938438,SRR2938437,SRR2938436,SRR2938435,SRR2938434,SRR2938433,SRR2938432,SRR2938431,SRR2938430,SRR2938429,SRR2938428,SRR2938427,SRR2938426,SRR2938425,SRR2938424,SRR2938423,SRR2938422,SRR2938421,SRR2938420,SRR2938419,SRR2938418,SRR2938417,SRR2938416,SRR2938415,SRR2938414,SRR2938413,SRR2938412,SRR2938411,SRR2938402,SRR2938401,SRR2938400,SRR2938399,SRR2938398,SRR2938397,SRR2938396,SRR2938395,SRR2938394,SRR2938393,SRR2938392,SRR2938391,SRR2938388,SRR2938387,SRR2938386,SRR2938385,SRR2938384,SRR2938383,SRR2938382,SRR2938381,SRR2938380,SRR2938379,SRR2938378,SRR2938377,SRR2938376,SRR2938375,SRR2938374,SRR2938373,SRR2938372,SRR2938371,SRR2938370,SRR2938369,SRR2938368,SRR2938367,SRR2938366,SRR2938365,SRR2938364,SRR2938363,SRR2938362,SRR2938361,SRR2938360,SRR2938357,SRR2938356,SRR2938355,SRR2938354,SRR2938353,SRR2938352,SRR2938351,SRR2938350,SRR2938349,SRR2938348,SRR2938347,SRR2938338,SRR2938337,SRR9033724,SRR9033725,SRR9033722,SRR9033723,SRR9033715,SRR9033716,SRR9033718,SRR9033754,SRR9033753,SRR9033750,SRR9033749,SRR9033752,SRR9033751,SRR9033746,SRR9033745,SRR9033748,SRR9033747,SRR9033720,SRR9033721,SRR9033732,SRR9033727,SRR9033743,SRR9033742,SRR9033737,SRR9033760,SRR9033719,SRR9033738,SRR341581,SRR341582,SRR341583,SRR341584,SRR341585,SRR341586,SRR341587,SRR341588,SRR341654,SRR341589,SRR341655,SRR341656,SRR341657,SRR341658,SRR341590,SRR341591,SRR341592,SRR341593,SRR341594,SRR341595,SRR341596,SRR341597,SRR341659,SRR341598,SRR341660,SRR341661,SRR341662,SRR341663,SRR341664,SRR341665,SRR341666,SRR341667,SRR341668,SRR341669,SRR341670,SRR341671,SRR341672,SRR341673,SRR341674,SRR341599,SRR341675,SRR341676,SRR341677,SRR341600,SRR341678,SRR341601,SRR341602,SRR341603,SRR341604,SRR341605,SRR341606,SRR341679,SRR341680,SRR341607,SRR341681,SRR341682,SRR341608,SRR341683,SRR341609,SRR341684,SRR341685,SRR341610,SRR341686,SRR341687,SRR341611,SRR341612,SRR341688,SRR341613,SRR341614,SRR341615,SRR341689,SRR341616,SRR341617,SRR341618,SRR341619,SRR341620,SRR341690,SRR341621,SRR341622,SRR341623,SRR341691,SRR341692,SRR341624,SRR341693,SRR341625,SRR341694,SRR341626,SRR341695,SRR341627,SRR341696,SRR341697,SRR341628,SRR341629,SRR341698,SRR341630,SRR341631,SRR341632,SRR341699,SRR341700,SRR341633,SRR341701,SRR341634,SRR341702,SRR341703,SRR341704,SRR341705,SRR341706,SRR341707,SRR341635,SRR341636,SRR341637,SRR341708,SRR341638,SRR341639,SRR341709,SRR341640,SRR341710,SRR341641,SRR341711,SRR341642,SRR341643,SRR341712,SRR341713,SRR341644,SRR341714,SRR341645,SRR341646,SRR341715,SRR341716,SRR341647,SRR341648,SRR341717,SRR341718,SRR341719,SRR341720,SRR341649,SRR341721,SRR341722,SRR341650,SRR341651,SRR341723,SRR341652,SRR341724,SRR341725,SRR341653,SRR413737,SRR413761,SRR413770,SRR413773,SRR413680,SRR413686,SRR413694,SRR413700,SRR413708,SRR413714,SRR413721,SRR413728,SRR413733,SRR413734,SRR413735,SRR413736,SRR413738,SRR413739,SRR413740,SRR413741,SRR413742,SRR413743,SRR413744,SRR1778450,SRR413745,SRR413746,SRR413747,SRR413748,SRR413749,SRR413750,SRR413751,SRR413752,SRR413754,SRR413755,SRR413756,SRR413757,SRR413758,SRR413759,SRR413760,SRR413762,SRR413763,SRR413764,SRR413765,SRR413766,SRR413767,SRR413768,SRR413769,SRR413771,SRR413772,SRR413673,SRR413674,SRR413675,SRR413676,SRR413677,ERR528329,ERR527046,ERR527047,ERR527048,ERR527050,ERR527051,ERR527052,ERR527055,ERR527057,ERR527059,ERR527061,ERR527062,ERR527063,ERR527064,ERR527065,ERR527066,ERR527067,ERR527069,ERR527070,ERR527071,ERR527072,ERR527073,ERR527075,ERR527076,ERR527077,ERR527078,ERR527079,ERR527080,ERR527081,ERR527082,ERR527084,ERR527085,ERR527087,ERR527088,ERR527089,ERR527090,ERR527091,ERR527093,ERR527095,ERR527097,ERR527098,ERR527100,ERR527101,ERR527103,ERR527104,ERR527105,ERR527106,ERR527107,ERR527108,ERR527110,ERR527111,ERR527113,ERR527114,ERR527116,ERR527118,ERR527120,ERR527121,ERR527122,ERR527123,ERR527124,ERR527125,ERR527126,ERR527127,ERR527128,ERR527129,ERR527130,ERR528330,ERR527131,ERR527132,ERR527133,ERR527134,ERR527135,ERR527136,ERR527137,ERR527138,ERR527139,ERR528331,ERR528332,ERR528333,ERR532377,ERR528335,ERR528336,ERR528337,ERR527006,ERR527007,ERR527009,ERR527010,ERR527012,ERR527013,ERR527016,ERR527017,ERR527018,ERR527020,ERR527022,ERR527023,ERR527025,ERR527026,ERR527027,ERR527028,ERR527029,ERR527030,ERR527031,ERR527032,ERR527033,ERR527034,ERR527035,ERR527036,ERR527037,ERR527039,ERR527040,ERR527041,ERR527042,ERR527043,ERR527045,ERR527140,ERR527142,ERR527143,ERR527144,ERR527146,ERR527147,ERR527148,ERR527149,ERR527150,ERR527152,ERR527153,ERR527155,ERR527156,ERR527157,ERR527159,ERR527160,ERR527162,ERR527163,ERR527165,ERR527168,ERR527170,ERR527171,ERR527172,ERR527175,ERR527176,ERR527180,ERR527181,ERR527182,ERR527183,ERR527186,ERR527187,ERR527188,ERR527189,ERR527190,ERR527192,ERR527194,ERR527195,ERR527196,ERR527197,ERR527198,ERR527199,ERR527201,ERR527203,ERR527204,ERR527206,ERR527207,ERR527208,ERR527209,ERR527210,ERR527211,ERR527213,ERR527215,ERR527216,ERR527217,ERR527218,ERR528722,ERR528723,ERR528724,ERR528725,ERR528726,ERR528728,ERR528730,ERR528732,ERR528734,ERR528735,ERR528736,ERR528737,ERR528739,ERR528741,ERR528742,ERR528743,ERR528744,ERR528745,ERR528746,ERR528747,ERR528748,ERR528749,ERR532378,ERR532379,ERR532380,ERR532381,ERR532382,ERR532383,ERR532384,ERR532385,ERR532387,ERR532388,ERR532389,ERR532618,ERR532390,ERR532391,ERR532392,ERR532393,ERR532394,ERR532619,ERR532620,ERR532621,ERR532622,ERR528291,ERR528293,ERR528295,ERR528296,ERR528298,ERR528300,ERR528302,ERR528303,ERR528304,ERR528306,ERR528307,ERR528309,ERR528310,ERR528311,ERR528312,ERR528313,ERR528314,ERR528316,ERR528319,ERR528320,ERR528323,ERR528324,ERR528325,ERR528327,ERR528328,SRR1927149,SRR1929408,SRR1929484,SRR1929485,SRR1929563,SRR1930121,SRR1929574,SRR1930122,SRR1930123,SRR1930128,SRR1930132,SRR1930133,SRR1930134,SRR1930136,SRR1930138,SRR1930140,SRR1930141,SRR1930142,SRR1930143,SRR1930144,SRR1930145,SRR1930149,SRR1930176]' --exact-matches --outdir batch_12 --apikey "${NCBI_API_KEY}" --sylph_db /work_beegfs/sukmb465/projects/TOFUpaper/sylph_db/gtdb-r220-c200-dbv1.syldb --sylph_processing
