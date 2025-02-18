#!/bin/bash
NXF_VER=23.10.1 ./nextflow run ikmb/tofu-maapo -r 6887e5c -profile custom -c tofu.config --sra '[ERR3957791,ERR3957792,ERR3957793,ERR3957794,SRR5298222,SRR5298221,SRR5298220,SRR5298219,SRR5298218,SRR5298217,SRR5298216,SRR5298215,SRR5298214,SRR5298213,SRR5298212,SRR5298211,SRR5298210,SRR5298209,SRR5298208,SRR5298207,SRR5298206,SRR5298205,SRR5298204,SRR5298203,SRR5298202,SRR5298201,SRR5298200,SRR5298199,SRR5298198,SRR5298197,SRR5298196,SRR1980447,SRR1980437,SRR1980304,SRR1980459,SRR1980442,SRR1916033,SRR1910591,SRR1918822,SRR1915879,SRR1915993,SRR1916000,SRR1910580,SRR1980320,SRR1980291,SRR1980357,SRR1980342,SRR1980312,SRR1980225,SRR1980150,SRR1980144,SRR1980161,SRR1910609,SRR1910631,SRR1915992,SRR1915996,SRR1910718,SRR1916030,SRR1910574,SRR1910621,SRR1910612,SRR1910653,SRR1916024,SRR1980429,SRR1980448,SRR1980433,SRR1980443,SRR1980440,SRR1916016,SRR1980456,SRR1980451,SRR1980432,SRR1910668,SRR1980428,SRR1980163,SRR1980484,SRR1980149,SRR1980164,SRR1980439,SRR1980162,SRR1980319,SRR1980327,SRR1980004,SRR1980347,SRR1980318,SRR1980301,SRR1980290,SRR1980365,SRR1980288,SRR1980435,SRR1980427,SRR1980431,SRR1980452,SRR1980389,SRR1980309,SRR1980310,SRR1980298,SRR1980305,SRR1980353,SRR1980293,SRR1980302,SRR1980397,SRR1980311,SRR1980289,SRR1980333,SRR1980345,SRR1910643,SRR1910570,SRR1916005,SRR1918825,SRR1910596,SRR1980336,SRR1980294,SRR1980387,SRR1980316,SRR1980297,SRR1980287,SRR1980299,SRR1980286,SRR1980300,SRR1980295,SRR1910638,SRR1910697,SRR1916041,SRR1910678,SRR1910592,SRR1980724,SRR1910654,SRR1910605,SRR1918914,SRR1916017,SRR1916018,SRR1916014,SRR1918927,SRR1910832,SRR1918905,SRR1910864,SRR1980238,SRR1980434,SRR1980328,SRR1980450,SRR1980430,SRR1980160,SRR1980479,SRR1910634,SRR1915997,SRR1915989,SRR1980227,SRR1980476,SRR1980473,SRR1980477,ERR1474587,ERR1474586,ERR1474585,ERR1474584,ERR1474583,ERR1474582,ERR1474581,ERR1474580,ERR1474579,ERR1474578,ERR1474577,ERR1474576,ERR1474575,ERR1474574,ERR1474573,ERR1474572,ERR1474571,ERR1474570,ERR1474568,ERR1474567,ERR1474566,ERR1474565,ERR1474564,ERR209908,ERR209911,ERR209914,ERR209916,ERR209918,ERR321079,ERR321081,ERR321083,ERR321086,ERR321088,ERR321091,ERR321098,ERR321100,ERR321103,ERR321105,ERR321108,ERR321110,ERR321112,ERR321114,ERR321116,ERR321118,ERR321120,ERR321122,ERR321124,ERR321126,ERR321128,ERR321130,ERR321132,ERR321134,ERR321136,ERR321138,ERR321140,ERR321142,ERR321144,ERR321146,ERR321148,ERR321150,ERR321152,ERR321154,ERR321156,ERR321158,ERR321160,ERR321162,ERR321164,ERR321166,ERR321168,ERR321170,ERR321172,ERR321174,ERR321176,ERR321178,ERR321180,ERR321182,ERR321184,ERR321186,ERR321188,ERR321190,ERR321192,ERR321194,ERR321196,ERR321198,ERR321200,ERR321202,ERR321204,ERR321206,ERR321208,ERR321210,ERR321213,ERR321215,ERR321217,ERR321219,ERR321221,ERR321223,ERR321225,ERR321226,ERR321229,ERR321231,ERR321233,ERR321235,ERR321237,ERR321239,ERR321241,ERR321243,ERR321245,ERR321247,ERR321251,ERR321252,ERR321256,ERR321260,ERR321261,ERR321268,ERR321274,ERR321280,ERR321286,ERR321292,ERR321298,ERR321302,ERR321307,ERR321311,ERR321316,ERR321321,ERR321326,ERR321333,ERR321338,ERR321345,ERR321351,ERR321353,ERR321358,ERR321364,ERR321367,ERR321371,ERR321376,ERR321379,ERR321385,ERR321389,ERR321392,ERR321393,ERR321397,ERR321398,ERR321403,ERR321407,ERR321412,ERR321416,ERR321420,ERR321424,ERR321428,ERR321432,ERR321436,ERR321437,ERR321441,ERR321445,ERR321449,ERR321453,ERR321457,ERR321461,ERR321462,ERR321463,ERR321464,ERR321465,ERR321466,ERR321467,ERR321468,ERR321470,ERR321471,ERR321475,ERR321479,ERR321480,ERR321482,ERR321483,ERR321484,ERR321485,ERR321486,ERR321487,ERR321488,ERR321489,ERR321490,ERR321491,ERR321492,ERR321493,ERR321494,ERR321496,ERR321497,ERR321502,ERR321505,ERR321509,ERR321513,ERR321515,ERR321519,ERR321520,ERR321521,ERR321522,ERR321523,ERR321524,ERR321528,ERR321529,ERR321530,ERR321531,ERR321532,ERR321533,ERR321537,ERR321541,ERR321542,ERR321543,ERR321544,ERR321546,ERR321547,ERR321548,ERR321550,ERR321551,ERR321552,ERR321553,ERR321554,ERR321555,ERR321556,ERR321557,ERR321558,ERR321559,ERR321560,ERR321561,ERR321562,ERR321563,ERR321564,ERR321565,ERR321567,ERR321568,ERR321569,ERR321570,ERR321571,ERR321572,ERR321573,ERR321574,ERR321575,ERR321576,ERR321577,ERR321578,ERR321579,ERR321580,ERR321581,ERR321582,ERR321583,ERR321584,ERR321586,ERR321587,ERR321588,ERR321589,ERR321590,ERR321591,ERR321592,ERR321593,ERR321594,ERR321595,ERR321596,ERR321597,ERR321598,ERR321599,ERR321600,ERR321601,ERR321602,ERR321603,ERR321604,ERR321605,ERR321606,ERR321607,ERR321608,ERR321609,ERR321610,ERR321611,ERR321612,ERR321613,ERR321614,ERR321615,ERR321616,ERR321617,ERR321618,ERR321619,ERR321620,ERR321621,ERR321622,ERR321623,ERR321624,ERR321625,ERR321626,ERR321627,ERR321628,ERR321629,ERR321630,ERR321631,ERR321632,ERR321633,ERR321634,ERR321635,ERR321636,ERR321637,ERR321638,ERR321639,ERR321640,ERR321641,ERR321642,ERR321643,ERR321644,ERR321645,ERR321646,ERR321647,ERR321648,ERR321649,ERR321650,ERR321651,ERR321652,ERR321653,ERR321654,ERR321655,ERR321656,ERR6281997,ERR6229432,ERR6282298,ERR6281944,ERR6281939,ERR6282070,ERR6282137,ERR6282383,ERR6281677,ERR6281658,ERR6281586,ERR6281941,ERR6279618,ERR6275675,ERR6275672,ERR6275677,ERR6275670,ERR6275669,ERR6275673,ERR6275667,ERR6275664,ERR6275665,ERR6275659,ERR6275666,ERR6275663,ERR6275668,ERR6275674,ERR6275652,ERR6271720,ERR6275671,ERR6279623,ERR6279621,ERR6279620,ERR6279647,ERR6279654,ERR6279652,ERR6279657,ERR6279669,ERR6279625,ERR6279628,ERR6279661,ERR6279650,ERR6279629,ERR6279651,ERR6279655,ERR6279649,ERR6279646,ERR6279648,ERR6279645,ERR6279626,ERR6279644,ERR6279659,ERR6279627,ERR6279670,ERR6279656,ERR6231616,ERR6230278,ERR6230340,ERR6231475,ERR6230326,ERR6233695,ERR6230245,ERR6275682,ERR6275683,ERR6275678,ERR6275679,ERR6279615,ERR6275681,ERR6279616,ERR6279619,ERR6279622,ERR6275684,ERR6279614,ERR6279624,ERR6275661,ERR6254871,ERR6254870,ERR6243879,ERR6248618,ERR6239138,ERR6243877,ERR6243878,ERR6243880,ERR6235717,ERR6236131,ERR6235716,ERR6239139,ERR6235718,ERR6233693,ERR6232687,ERR6233692,ERR6232688,ERR6233694,ERR6232526,ERR6232529,ERR6232686,ERR6232528,ERR6232525,ERR6232523,ERR6232522,ERR6232521,ERR6232518,ERR6232527,ERR6232425,ERR6232382,ERR6231548,ERR6232458,ERR6232524,ERR6275680,ERR6275676,ERR6281751,ERR6281505,ERR6281942,ERR6281433,ERR6281431,ERR6281448,ERR6281500,ERR6281342,ERR6281340,ERR6281432,ERR6280948,ERR6281434,ERR6281341,ERR6281343,ERR6280910,ERR6279683,ERR6281000,ERR6279682,ERR6280907,ERR6279680,ERR6279679,ERR6279678,ERR6279673,ERR6279676,ERR6279677,ERR6279681,ERR6279668,ERR6279672,ERR6279664,ERR6279666,ERR6279663,ERR6279662,ERR6279658,ERR6279660,ERR6279653,ERR6271719,ERR6271721,ERR6257997,ERR6273965,ERR6281499,ERR6281504,ERR6281518,ERR6281344,ERR6279675,ERR6280899,ERR6279671,ERR6279667,ERR6264280,ERR6267712,ERR6279674,ERR6267711,ERR6254869,ERR6279665,ERR6264281,ERR504979,ERR504980,ERR504985,ERR504981,ERR414512,ERR414513,ERR414514,ERR414515,ERR414516,ERR414517,ERR414519,ERR414520,ERR414521,ERR414522,ERR414523,ERR414524,ERR414525,ERR414526,ERR414527,ERR414528,ERR414529,ERR414531,ERR414532,ERR414533,ERR414535,ERR414536,ERR414537,ERR414539,ERR414541,ERR414542,ERR414543,ERR414544,ERR414545,ERR414546,ERR414547,ERR414548,ERR414549,ERR414551,ERR414553,ERR414554,ERR414555,ERR414556,ERR414557,ERR414558,ERR414559,ERR414560,ERR414561,ERR414562,ERR414563,ERR414564,ERR414565,ERR414567,ERR414568,ERR414569,ERR414570,ERR414571,ERR414572,ERR414573,ERR414574,ERR414575,ERR414576,ERR414577,ERR414578,ERR414579,ERR414581,ERR414582,ERR414583,ERR414584,ERR414585,ERR414586,ERR414587,ERR414588,ERR414589,ERR414590,ERR414591,ERR414592,ERR414593,ERR414594,ERR414595,ERR414596,ERR414597,ERR414598,ERR414599,ERR414600,ERR414601,ERR414602,ERR414603,ERR414605,ERR414606,ERR414607,ERR414608,ERR414610,ERR414611,ERR414612,ERR414613,ERR414614,ERR414615,ERR414616,ERR414617,ERR414618,ERR414619,ERR414620,ERR414621,ERR414622,ERR414623,ERR414624,ERR414625,ERR414626,ERR414627,ERR414628,ERR414629,ERR414631,ERR414632,ERR504986,ERR504987,ERR504982,ERR504988,ERR504983,ERR504984,ERR504989,ERR414633,ERR414634,ERR414635,ERR414636,ERR414637,ERR414639,ERR414638,ERR414640,ERR414641,ERR414642,ERR414643,ERR414644,ERR414645,ERR414646,ERR414647,ERR414648,ERR414649,ERR414650,ERR414652,ERR414651,ERR414653,ERR414654,ERR414655,ERR414656,ERR414657,ERR414658,ERR414659,ERR414660,ERR414661,ERR414662,ERR414664,ERR414663,ERR414665,ERR414666,ERR414668,ERR414669,ERR414670,ERR414671,ERR414672,ERR414674,ERR414673,ERR414675,ERR414676,ERR414677,ERR414678,ERR414679,ERR414680,ERR414681,ERR414682,ERR414683,ERR414684,ERR414685,ERR414686,ERR414687,ERR414688,ERR414689,ERR414690,ERR414691,ERR414693,ERR414694,ERR414696,ERR414697,ERR414698,ERR414699,ERR414700,ERR414702,ERR414704,ERR414705,ERR414706,ERR414708,ERR414709,ERR414710,ERR414711,ERR414713,ERR414714,ERR414716,ERR414718,ERR414719,ERR414721,ERR414722,ERR414723,ERR414724,ERR414725,ERR414726,ERR414727,ERR414728,ERR414729,ERR414730,ERR414731,ERR414732,ERR414733,ERR414734,ERR414735,ERR414736,ERR414737,ERR414738,ERR414739,ERR414740,ERR414741,ERR414742,ERR414743,ERR414744,ERR414745,ERR414746,ERR414747,ERR414748,ERR414749,ERR414750,ERR414751,ERR414752,ERR414753,ERR414754,ERR414755,ERR414756,ERR414757,ERR414758,ERR414759,ERR414760,ERR414761,ERR414762,ERR414763,ERR414764,ERR414765,ERR414766,ERR414767,ERR414768,ERR414769,ERR414771,ERR414773,ERR414774,ERR414776,ERR414777,ERR414778,ERR414780,ERR414781,ERR414783,ERR414784,ERR414787,ERR414790,ERR414791,ERR1398129,ERR1398173,ERR1398206,ERR1398242,ERR1398075,ERR1398153,ERR1398078,ERR1398214,ERR1398192,ERR1398138,ERR1398083,ERR1398161,ERR1398223,ERR1398253,ERR1398218,ERR1398178,ERR1398113,ERR1398180,ERR1398148,ERR1398144,ERR1398217,ERR1398169,ERR1398244,ERR1398118,ERR1398164,ERR1398243,ERR1398198,ERR1398127,ERR1398089,ERR1398100,ERR1398263,ERR1398126,ERR1398248,ERR1398237,ERR1398224,ERR1398114,ERR1398152,ERR1398101,ERR1398109,ERR1398213,ERR1398205,ERR1398111,ERR1398190,ERR1398071,ERR1398209,ERR1398165,ERR1398195,ERR1398246,ERR1398193,ERR1398188,ERR1398249,ERR1398091,ERR1398081,ERR1398096,ERR1398167,ERR1398072,ERR1398241,ERR1398201,ERR1398155,ERR1398108,ERR1398084,ERR1398086,ERR1398106,ERR1398150,ERR1398141,ERR1398231,ERR1398158,ERR1398234,ERR1398092,ERR1398171,ERR1398074,ERR1398160,ERR1398098,ERR1398215,ERR1398133,ERR1398087,ERR1398156,ERR1398240,ERR1398259,ERR1398222,ERR1398175,ERR1398080,ERR1398210,ERR1398140,ERR1398070,ERR1398183,ERR1398121,ERR1398194,ERR1398255,ERR1398219,ERR1398094,ERR1398139,ERR1398187,ERR1398196,ERR1398257,ERR1398147,ERR1398102,ERR1398132,ERR1398131,ERR1398174,ERR1398212,ERR1398232,ERR1398251,ERR1398095,ERR1398189,ERR1398211,ERR1398076,ERR1398235,ERR1398110,ERR1398135]' --exact-matches --outdir batch_7 --apikey "${NCBI_API_KEY}" --sylph_db /work_beegfs/sukmb465/projects/TOFUpaper/sylph_db/gtdb-r220-c200-dbv1.syldb --sylph_processing
