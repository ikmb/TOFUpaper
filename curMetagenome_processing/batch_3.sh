#!/bin/bash
NXF_VER=23.10.1 ./nextflow run ikmb/tofu-maapo -r 6887e5c -profile custom -c tofu.config --sra '[SRR5405977,SRR5405976,SRR5405975,SRR5405974,SRR5405973,SRR5405972,SRR5405971,SRR5405970,SRR5405969,SRR5405968,SRR5405967,SRR5405966,SRR5405965,SRR5405964,SRR5405963,SRR5405962,SRR5405961,SRR5405960,SRR5405959,SRR5405958,SRR5405957,SRR5405956,SRR5405955,SRR5405954,SRR5405953,SRR5405952,SRR5405951,SRR5405950,SRR5405949,SRR5405948,SRR5405947,SRR5405946,SRR5405945,SRR5405944,SRR5405943,SRR5405942,SRR5405941,SRR5405940,SRR5405939,SRR5405938,SRR5405937,SRR5405936,SRR5405935,SRR5405934,SRR5405933,SRR5405932,SRR5405931,SRR5405930,SRR5405929,SRR5405928,SRR5405927,SRR5405926,SRR5405925,SRR5405924,SRR5405923,SRR5405922,SRR5405921,SRR5405920,SRR5405919,SRR5405918,SRR5405917,SRR5405916,SRR5405915,SRR5405914,SRR5405913,SRR5405912,SRR5405911,SRR5405910,SRR5405909,SRR5405908,SRR5405907,SRR5405906,SRR5405905,SRR5405904,SRR5405903,SRR5405902,SRR5405901,SRR5405900,SRR5405899,SRR5405898,SRR5405897,SRR5405896,SRR5405895,SRR5405894,SRR5405893,SRR5405892,SRR5405891,SRR5405890,SRR5405889,SRR5405888,SRR5405887,SRR5405886,SRR5405885,SRR5405884,SRR5405883,SRR5405882,SRR5405881,SRR5405880,SRR5405879,SRR5405878,SRR5405877,SRR5405876,SRR5405875,SRR5405874,SRR5405873,SRR5405872,SRR5405871,SRR5405870,SRR5405869,SRR5405868,SRR5405867,SRR5405866,SRR5405865,SRR5405864,SRR5405863,SRR5405862,SRR5405861,SRR5405860,SRR5405859,SRR5405858,SRR5405857,SRR5405856,SRR5405855,SRR5405854,SRR5405853,SRR5405852,SRR5405851,SRR5405850,SRR5405849,SRR5405848,SRR5405847,SRR5405846,SRR5405845,SRR5405844,SRR5405843,SRR5405842,SRR5405841,SRR5405840,SRR5405839,SRR5405838,SRR5405837,SRR5405836,SRR5405835,SRR5405834,SRR5405833,SRR5405832,SRR5405831,SRR5405830,SRR5405829,SRR5405828,SRR5405827,SRR5405826,SRR5405825,SRR5405824,SRR5405823,SRR5405822,SRR5405821,SRR5405820,SRR5405819,SRR5405818,SRR5405817,SRR5405816,SRR5405815,SRR5405814,SRR5405813,SRR5405812,SRR5405811,SRR5405810,SRR5405809,SRR5405808,SRR5405807,SRR5405806,SRR5405805,SRR5405804,SRR5405803,SRR5405802,SRR5405801,SRR5405800,SRR5405799,SRR5405798,SRR5405797,SRR5405796,SRR5405795,SRR5405794,SRR5405793,SRR5405792,SRR5405791,SRR5405790,SRR5405789,SRR5405788,SRR5405787,SRR5405786,SRR5405785,SRR5405784,SRR5405783,SRR5405782,SRR5405781,SRR5405780,SRR5405779,SRR5405778,SRR5405777,SRR5405776,SRR5405775,SRR5405774,SRR5405773,SRR5405772,SRR5405771,SRR5405770,SRR5405769,SRR5405768,SRR5405767,SRR5405766,SRR5405765,SRR5405764,SRR5405763,SRR5405762,SRR5405761,SRR5405760,SRR5405759,SRR5405758,SRR5405757,SRR5405756,SRR5405755,SRR5405754,SRR5405753,SRR5405752,SRR5405751,SRR5405750,SRR5405749,SRR5405748,SRR5405747,SRR5405746,SRR5405745,SRR5405744,SRR5405743,SRR5405742,SRR5405741,SRR5405740,SRR5405739,SRR5405738,SRR5405737,SRR5405736,SRR5405735,SRR5405734,SRR5405733,SRR5405732,SRR5405731,SRR5405730,SRR5405729,SRR5405728,SRR5405727,SRR5405726,SRR5405725,SRR5405724,SRR5405723,SRR5405722,SRR5405721,SRR5405720,SRR5405719,SRR5405718,SRR5405717,SRR5405716,SRR5405715,SRR5405714,SRR5405713,SRR5405712,SRR5405711,SRR5405710,SRR5405709,SRR5405708,SRR5405707,SRR5405706,SRR5405705,SRR5405704,SRR5405703,SRR5405702,SRR5405701,SRR5405700,SRR5405699,SRR5405698,SRR5405697,SRR5405696,SRR5405695,SRR5405694,SRR5405693,SRR5405692,SRR5405691,SRR5405690,SRR5405689,SRR5405688,SRR5405687,SRR5405686,SRR5405685,SRR5405684,SRR5405683,SRR5405682,SRR5405681,SRR5405680,SRR5405679,SRR5405678,SRR5405677,SRR5405676,SRR5405675,SRR5405674,SRR5405673,SRR5405672,SRR5405671,SRR5405670,SRR5405669,SRR5405668,SRR5405667,SRR5405666,SRR5405665,SRR5405664,SRR5405663,SRR5405662,SRR5405661,SRR5405660,SRR5405659,SRR5405658,SRR5405657,SRR5405656,SRR5405655,SRR5405654,SRR5405653,SRR5405652,SRR5405651,SRR5405650,SRR5405649,SRR5405648,SRR5405647,SRR5405646,SRR5405645,SRR5405644,SRR5405643,SRR5405642,SRR5405641,SRR5405640,SRR5405639,SRR5405638,SRR5405637,SRR5405636,SRR5405635,SRR5405634,SRR5405633,SRR5405632,SRR5405631,SRR5405630,SRR5405629,SRR5405628,SRR5405627,SRR5405626,SRR5405625,SRR5405624,SRR5405623,SRR5405622,SRR5405621,SRR5405620,SRR5405619,SRR5405618,SRR5405617,SRR5405616,SRR5405615,SRR5405614,SRR5405613,SRR5405612,SRR5405611,SRR5405610,SRR5405609,SRR5405608,SRR5405607,SRR1519058,SRR1519059,SRR1519060,SRR1519061,SRR1519071,SRR1519072,SRR1519073,SRR1519075,SRR1519076,SRR1519077,SRR1518393,SRR1519078,SRR1518476,SRR1518536,SRR1518998,SRR1519031,SRR1519044,SRR1519045,SRR1519046,SRR1519079,SRR1519047,SRR1519048,SRR1519080,SRR1519082,SRR1519083,SRR1519084,SRR1519049,SRR1519050,SRR1519052,SRR1519055,SRR1519056,SRR1519057,SRR5279282,SRR5279276,SRR5279270,SRR5279263,SRR5279258,SRR5279308,SRR5279288,SRR5279287,SRR5279286,SRR5279285,SRR5279284,SRR5279283,SRR5279281,SRR5279280,SRR5279279,SRR5279278,SRR5279277,SRR5279275,SRR5279274,SRR5279273,SRR5279272,SRR5279271,SRR5279269,SRR5279268,SRR5279267,SRR5279266,SRR5279265,SRR5279264,SRR5279262,SRR5279261,SRR5279260,SRR5279259,SRR5279257,SRR5279256,SRR5279255,SRR5279254,SRR5279253,SRR5279252,SRR5279251,SRR5279250,SRR5279249,SRR5279248,SRR5279247,SRR5279246,SRR5279245,SRR5279244,SRR5279243,SRR5279242,SRR5279241,SRR5279313,SRR5279312,SRR5279311,SRR5279310,SRR5279309,SRR5279307,SRR5279306,SRR5279305,SRR5279304,SRR5279303,SRR5279302,SRR5279301,SRR5279300,SRR5279299,SRR5279298,SRR5279297,SRR5279296,SRR5279295,SRR5279294,SRR5279293,SRR5279292,SRR5279291,SRR5279290,SRR5279289,SRR5279240,SRR5279239,SRR5279238,SRR5279237,SRR5279236,SRR5279235,SRR5279234,SRR5279233,SRR5279232,SRR5279231,SRR5279230,SRR5279229,SRR5279228,SRR5279227,SRR5279226,SRR5279225,SRR5279224,SRR5279223,SRR5279222,SRR5279221,SRR5279220,SRR5279219,SRR5279218,SRR5279217,SRR1950743,SRR1950744,SRR1950713,SRR1950714,SRR1950715,SRR1950716,SRR1950745,SRR1950746,SRR1950747,SRR1950748,SRR1950749,SRR1950750,SRR1950751,SRR1950752,SRR1950717,SRR1950718,SRR1950719,SRR1950720,SRR1950753,SRR1950754,SRR1950755,SRR1950756,SRR1950757,SRR1950758,SRR1950759,SRR1950760,SRR1950761,SRR1950762,SRR1950763,SRR1950764,SRR1950765,SRR1950766,SRR1950767,SRR1950768,SRR1950769,SRR1950770,SRR1950771,SRR1950772,SRR1950790,SRR1950791,SRR1950773,SRR1950774,SRR1950775,SRR1950776,SRR1950777,SRR1950778,SRR1950782,SRR1950783,SRR1950784,SRR1950785,SRR1950721,SRR1950722,SRR1950723,SRR1950724,SRR1950786,SRR1950787,SRR1950788,SRR1950789,SRR1950725,SRR1950726,SRR1950727,SRR1950728,SRR1950779,SRR1950780,SRR1950729,SRR1950730,SRR1950731,SRR1950732,SRR1950733,SRR1950734,SRR1950735,SRR1950736,SRR1950737,SRR1950738,SRR1950739,SRR1950740,SRR1950741,SRR1950742,SRR8942345,SRR8942344,SRR8942347,SRR8942349,SRR8942350,SRR8942348,SRR8942351,SRR8942352,SRR8942354,SRR8942355,SRR8942353,SRR8942358,SRR8942357,SRR8942359,SRR8942361,SRR8942360,SRR8942362,SRR8942363,SRR8942364,SRR8942365,SRR8942366,SRR8942367,SRR8942369,SRR8942370,SRR8942375,SRR8942374,SRR8942371,SRR8942372,SRR8942373,SRR8942376,SRR8942377,SRR8942378,SRR8942379,SRR8942380,SRR8942381,SRR8942382,SRR8942383,SRR8942384,SRR8942385,SRR8942386,SRR8942388,SRR8942387,SRR8942389,SRR8942391,SRR8942390,SRR8942392,SRR8942393,SRR8942394,SRR8942395,SRR8942396,SRR8942397,SRR8942398,SRR8942399,SRR8942400,SRR8942402,SRR8942401,SRR8942403,SRR8942405,SRR8942406,SRR8942408,SRR8942407,SRR8942409,SRR8942410,SRR8942412,SRR8942411,SRR8942413,SRR8942414,SRR8942415,SRR8942416,SRR8942417,SRR8942420,SRR8942419,SRR8942418,SRR8942421,SRR8942422,SRR8942425,SRR8942423,SRR8942424,SRR8942426,SRR8942427,SRR8942429,SRR8942428,SRR8942430,SRR8942431,SRR8942432,SRR8942343,ERR1727744,ERR1727562,ERR1728278,ERR1727976,ERR1727736,ERR1727704,ERR1727752,ERR1727454,ERR1728212,ERR1727387,ERR1727654,ERR1727506,ERR1727586,ERR1728172,ERR1727475,ERR1727799,ERR1727832,ERR1728150,ERR1727335,ERR1728138,ERR1728338,ERR1728078,ERR1728274,ERR1727622,ERR1728377,ERR1727860,ERR1727630,ERR1728282,ERR1728385,ERR1727462,ERR1728000,ERR1728254,ERR1727700,ERR1728073,ERR1727816,ERR1727676,ERR1728146,ERR1727940,ERR1728114,ERR1728330,ERR1727650,ERR1727892,ERR1728224,ERR1727433,ERR1728322,ERR1728373,ERR1727534,ERR1728044,ERR1727415,ERR1728107,ERR1727582,ERR1727522,ERR1727466,ERR1727904,ERR1727445,ERR1727618,ERR1728159,ERR1727530,ERR1727411,ERR1727423,ERR1728167,ERR1727407,ERR1728102,ERR1728036,ERR1727606,ERR1727803,ERR1728310,ERR1728058,ERR1727784,ERR1728228,ERR1727792,ERR1727391,ERR1728354,ERR1727482,ERR1727913,ERR1727900,ERR1728216,ERR1728314,ERR1728040,ERR1727610,ERR1727780,ERR1727808,ERR1727291,ERR1727668,ERR1727471,ERR1727526,ERR1728110,ERR1727498,ERR1727287,ERR1727944,ERR1727594,ERR1728008,ERR1727479,ERR1727680,ERR1727856,ERR1728119,ERR1727419,ERR1727355,ERR1727371,ERR1728163,ERR1727343,ERR1727494,ERR1727988,ERR1727764,ERR1728342,ERR1727924,ERR1728298,ERR1727663,ERR1728326,ERR1727844,ERR1727598,ERR1728397,ERR1727510,ERR1728082,ERR1728061,ERR1728196,ERR1727692,ERR1728122,ERR1727347,ERR1727283,ERR1727840,ERR1728270,ERR1728054,ERR1727952,ERR1728090,ERR1728098,ERR1727972,ERR1727542,ERR1727696,ERR1727936,ERR1728188,ERR1727864,ERR1727359,ERR1727590,ERR1728142,ERR1727638,ERR1728401,ERR1727646,ERR1728032,ERR1728050,ERR1728020,ERR1727996,ERR1727518,ERR1727578,ERR1727716,ERR1727502,ERR1728016,ERR1727319,ERR1727968,ERR1728012,ERR1727724,ERR1727896,ERR1728246,ERR1727772,ERR1727812,ERR1727307,ERR1727339,ERR1727626,ERR1727732,ERR1728262,ERR1728290,ERR1727852,ERR1727672,ERR1727538,ERR1728155,ERR1728204,ERR1727303,ERR1727311,ERR1728180,ERR1727720,ERR1728126,ERR1727920,ERR1728306,ERR1728294,ERR1727828,ERR1727928,ERR1727880,ERR1727327,ERR1727546,ERR1727634,ERR1728065,ERR1727948,ERR1727375,ERR1728393,ERR1727667,ERR1727399,ERR1728028,ERR1728266,ERR1727614,ERR1727909,ERR1727437,ERR1727956,ERR1728357,ERR1728389,ERR1727876,ERR1727602,ERR1727836,ERR1728024,ERR1728069,ERR1727558,ERR1727554,ERR1727490,ERR1728250,ERR1727514,ERR1727379,ERR1727868,ERR1727980,ERR1728004,ERR1727658,ERR1727888,ERR1727642,ERR1728302,ERR1727776,ERR1728134,ERR1727299,ERR1728094,ERR1728238,ERR1728286,ERR1728346,ERR1727295,ERR1728086,ERR1727683,ERR1727795,ERR1727728,ERR1727848,ERR1728234,ERR1727429,ERR1727367,ERR1727748,ERR1727458,ERR1727486,ERR1727351,ERR1728176,ERR1727820,ERR1727574,ERR1727688,ERR1727331,ERR1727824,ERR1727550,ERR1728192,ERR1728334,ERR1727708,ERR1727403,ERR1727684,ERR1727964,ERR1727984,ERR1728200,ERR1727960,ERR1727383,ERR1727450,ERR1728369,ERR1727872,ERR1728350,ERR1728208,ERR1727428,ERR1727395,ERR1728381,ERR1728230,ERR1727992,ERR1728046,ERR1727441,ERR1727740,ERR1727323,ERR1727788,ERR1728130,ERR1728220,ERR1727662,ERR1727760,ERR1728258,ERR1727756,ERR1727566,ERR1727315,ERR1727570,ERR1728184,ERR1727712,ERR1728242,ERR1728361,ERR1727932,ERR1727884,ERR866583,ERR866584,ERR866585,ERR866591,ERR866592,ERR866596,ERR866601,ERR866608,ERR866609,ERR866561,ERR866562,ERR866563,ERR866564,ERR866565,ERR866566,ERR866567,ERR866568,ERR866569,ERR866570,ERR866571,ERR866572,ERR866573,ERR866574,ERR866575,ERR866576,ERR866577,ERR866578,ERR866579,ERR866580,ERR866581,ERR866582,ERR866586,ERR866587,ERR866588,ERR866590,ERR866593,ERR866594,ERR866595,ERR866597,ERR866598,ERR866599,ERR866600,ERR866602,ERR866604,ERR866605,ERR866606,ERR866607,SRR4074259,SRR8146952,SRR8146951,SRR8146937,SRR8146956,SRR4074260,SRR8146957,SRR4074287,SRR6367589,SRR4074351]' --exact-matches --outdir batch_3 --apikey "${NCBI_API_KEY}" --sylph_db /work_beegfs/sukmb465/projects/TOFUpaper/sylph_db/gtdb-r220-c200-dbv1.syldb --sylph_processing
