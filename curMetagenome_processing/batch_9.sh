#!/bin/bash
NXF_VER=23.10.1 ./nextflow run ikmb/tofu-maapo -r 6887e5c -profile custom -c tofu.config --sra '[SRR6028190,SRR6028189,SRR6028192,SRR6028191,SRR6028186,SRR6028185,SRR6028188,SRR6028187,SRR6028225,SRR6028224,SRR6028382,SRR6028383,SRR6028384,SRR6028385,SRR6028386,SRR6028537,SRR6028538,SRR6028551,SRR6028539,SRR6028540,SRR6028373,SRR6028372,SRR6028371,SRR6028370,SRR6028369,SRR6028368,SRR6028448,SRR6028447,SRR6028223,SRR6028329,SRR6028491,SRR6028492,SRR5155843,SRR6028489,SRR6028490,SRR6028495,SRR6028496,SRR6028493,SRR6028494,SRR6028497,SRR6028498,SRR6028289,SRR6028288,SRR6028291,SRR6028290,SRR6028293,SRR6028292,SRR6028295,SRR6028294,SRR6028297,SRR6028296,SRR6028171,SRR6028172,SRR6028173,SRR6028174,SRR6038220,SRR6038221,SRR6038218,SRR6038219,SRR6028165,SRR6028166,SRR6028311,SRR6028310,SRR6028309,SRR6028308,SRR6028315,SRR6028314,SRR6028313,SRR6028312,SRR6054644,SRR6054643,SRR6054463,SRR6054587,SRR6038212,SRR6038213,SRR6038355,SRR6038356,SRR6038559,SRR6038560,SRR6038359,SRR6038360,SRR6038357,SRR6038358,SRR6038557,SRR6038558,SRR6038463,SRR6038462,SRR6038465,SRR6038464,SRR6054537,SRR6054538,SRR6054539,SRR6054540,SRR6038467,SRR6038469,SRR6038468,SRR6038461,SRR6038460,SRR6038520,SRR6038521,SRR6054504,SRR6054503,SRR6054502,SRR6054501,SRR6038522,SRR6038523,SRR6038524,SRR6038525,SRR6038526,SRR6038183,SRR6038184,SRR6038185,SRR6054563,SRR6054464,SRR6054469,SRR6054470,SRR6054432,SRR6038304,SRR6038303,SRR6038302,SRR6038301,SRR6054434,SRR6054437,SRR6054436,SRR6054439,SRR6054438,SRR6054441,SRR6054440,SRR6054443,SRR6054474,SRR6054692,SRR6054693,SRR6054694,SRR6054695,SRR6054688,SRR6054689,SRR6054690,SRR6054691,SRR6054697,SRR6054698,SRR6038300,SRR6038299,SRR6038298,SRR6038297,SRR6054658,SRR6054657,SRR6054656,SRR6054655,SRR6038296,SRR6038295,SRR6038351,SRR6038352,SRR6054668,SRR6054445,SRR6054700,SRR6054607,SRR6054604,SRR6054444,SRR6054621,SRR6054495,SRR6038349,SRR6038350,SRR6038347,SRR6038348,SRR6038345,SRR6038346,SRR6038353,SRR6038354,SRR6038418,SRR6038417,SRR6038420,SRR6038419,SRR6054571,SRR6054572,SRR6038414,SRR6038413,SRR6038416,SRR6038415,SRR6054569,SRR6054564,SRR6054565,SRR6054611,SRR6054610,SRR6054609,SRR6054608,SRR6054615,SRR6054614,SRR6054613,SRR6054612,SRR6054606,SRR6054605,SRR6038422,SRR6038421,SRR6038476,SRR6038477,SRR6054639,SRR6054640,SRR6054637,SRR5165875,SRR6054638,SRR6054646,SRR6054647,SRR6054677,SRR6054676,SRR6054679,SRR6054678,SRR6038478,SRR6038479,SRR6038472,SRR6038473,SRR6038474,SRR6038475,SRR6038470,SRR6038471,SRR6054424,SRR6054425,SRR6038247,SRR6038246,SRR6038245,SRR6038244,SRR6054429,SRR6054430,SRR6054431,SRR6054461,SRR6054460,SRR6054459,SRR6054458,SRR6054457,SRR6054456,SRR6054455,SRR6038251,SRR6038250,SRR6038249,SRR6038248,SRR6038243,SRR6038242,SRR6038436,SRR6038437,SRR6054486,SRR6054483,SRR6054484,SRR6054476,SRR6054477,SRR6054517,SRR6054516,SRR5165828,SRR6038434,SRR6038435,SRR6038440,SRR6054520,SRR6054523,SRR6054522,SRR6054509,SRR6054508,SRR6054551,SRR6054552,SRR6054553,SRR6054554,SRR6054547,SRR6054548,SRR6054549,SRR6054550,SRR6038441,SRR6038361,SRR6038362,SRR6038444,SRR6038364,SRR6038423,SRR6038277,SRR6038199,SRR6054595,SRR6054594,SRR6054593,SRR6054592,SRR6054599,SRR6054598,SRR6054667,SRR6054666,SRR6054665,SRR6054664,SRR6054663,SRR6054662,SRR6054661,SRR6054660,SRR6054671,SRR6054670,SRR6038198,SRR6038201,SRR6038200,SRR6038203,SRR6038202,SRR6038294,SRR6038293,SRR6038276,SRR6054625,SRR6054626,SRR6054447,SRR6054446,SRR6054449,SRR6054448,SRR5164261,SRR6038275,SRR6038274,SRR6038273,SRR6038279,SRR6038278,SRR6038272,SRR6038282,SRR6038281,SRR6038504,SRR6038505,SRR6038506,SRR6054706,SRR6054416,SRR6054419,SRR6054420,SRR6054527,SRR6054526,SRR6054525,SRR6054524,SRR6054531,SRR6054530,SRR6054529,SRR6054528,SRR6054533,SRR6054532,SRR6054491,SRR6054492,SRR6054493,SRR6054494,SRR6054487,SRR6054488,SRR6054489,SRR6054490,SRR6054496,SRR6054497,SRR6054601,SRR6054600,SRR6038507,SRR6038500,SRR6038501,SRR6038502,SRR6054628,SRR6054467,SRR6054616,SRR6054659,SRR6054561,SRR6054562,SRR6038503,SRR6038508,SRR6038509,SRR6038447,SRR6038446,SRR6038449,SRR6038448,SRR6038443,SRR6054687,SRR6054686,SRR5165715,SRR6054585,SRR6054535,SRR6054683,SRR6054682,SRR6038442,SRR6038445,SRR5165708,SRR6038363,SRR6038439,SRR6038438,SRR6038545,SRR6038546,SRR6451751,SRR6451750,SRR6451749,SRR6451748,SRR6038543,SRR6038544,SRR6451745,SRR6451744,SRR6451729,SRR6451728,SRR6451731,SRR6451730,SRR6451725,SRR6451724,SRR6451727,SRR6038542,SRR6038539,SRR6038540,SRR6038537,SRR6451831,SRR6451832,SRR6451833,SRR6451834,SRR6451835,SRR6451836,SRR6451837,SRR6451828,SRR6451829,SRR6451827,SRR6451826,SRR6451825,SRR6038538,SRR6038322,SRR6038321,SRR6038320,SRR6451820,SRR6451819,SRR6451818,SRR6451808,SRR6451809,SRR6451806,SRR6038319,SRR6038318,SRR6038317,SRR6038316,SRR6038315,SRR6038324,SRR6038323,SRR6038204,SRR6451796,SRR6451799,SRR6451798,SRR6451801,SRR6451800,SRR6451803,SRR6451802,SRR6451795,SRR6451794,SRR6038205,SRR6038206,SRR6038207,SRR6038208,SRR6038209,SRR6038210,SRR6038211,SRR6038196,SRR6451784,SRR6451785,SRR6451767,SRR6451766,SRR6451765,SRR6451764,SRR6038197,SRR6038481,SRR5165628,SRR6038480,SRR6451768,SRR6451773,SRR6451772,SRR6451742,SRR6451743,SRR6451740,SRR6451741,SRR6451738,SRR6038483,SRR6038482,SRR6038485,SRR6038484,SRR6451735,SRR6451721,SRR6451720,SRR6451723,SRR6451722,SRR6451717,SRR6451716,SRR5164116,SRR6451719,SRR6451718,SRR6451715,SRR6451714,SRR6451758,SRR6451759,SRR6038489,SRR6038488,SRR6451754,SRR6451755,SRR6038377,SRR6038378,SRR6451762,SRR6451763,SRR6038381,SRR6038382,SRR6451775,SRR6451774,SRR6451781,SRR6451780,SRR6451779,SRR6451778,SRR6451783,SRR6451782,SRR6451816,SRR6451817,SRR6451814,SRR6451815,ERR4086424,ERR4563248,ERR4560729,ERR4563253,ERR4560734,ERR4563258,ERR4560738,ERR4088391,ERR4088396,ERR4560743,ERR4560748,ERR4560753,ERR4563263,ERR4088401,ERR4563268,ERR4088408,ERR4088416,ERR4560759,ERR4560764,ERR4088420,ERR4563273,ERR4560769,ERR4086431,ERR4088425,ERR4560774,ERR4563278,ERR4560779,ERR4086438,ERR4563282,ERR4086442,ERR4088434,ERR4563287,ERR4560783,ERR4086447,ERR4563291,ERR4560788,ERR4560793,ERR4086451,ERR4086465,ERR4563299,ERR4563304,ERR4560799,ERR4563311,ERR4563315,ERR4086469,ERR4563319,ERR4088440,ERR4563323,ERR4563329,ERR4563334,ERR4560807,ERR4086478,ERR4088444,ERR4560811,ERR4563335,ERR4086482,ERR4086485,ERR4563341,ERR4088449,ERR4560815,ERR4563346,ERR4088453,ERR4086493,ERR4086497,ERR4563350,ERR4086506,ERR4088458,ERR4560820,ERR4086514,ERR4560823,ERR4088462,ERR4088469,ERR4560828,ERR4086520,ERR4086525,ERR4563355,ERR4563360,ERR4086529,ERR4088471,ERR4088472,ERR4086533,ERR4088474,ERR4563370,ERR4086537,ERR4088478,ERR4560833,ERR4086538,ERR4086546,ERR4560838,ERR4563374,ERR4088484,ERR4086556,ERR4088489,ERR4086561,ERR4088494,ERR4563377,ERR4560842,ERR4088502,ERR4563382,ERR4563386,ERR4563391,ERR4560847,ERR4560852,ERR4088506,ERR4563395,ERR4560858,ERR4563399,ERR4086567,ERR4086572,ERR4563403,ERR4086579,ERR4086587,ERR4086594,ERR4086601,ERR4088511,ERR4088517,ERR4563407,ERR4560862,ERR4086607,ERR4563412,ERR4086612,ERR4563416,ERR4560866,ERR4086618,ERR4563420,ERR4086625,ERR4086632,ERR4560870,ERR4088526,ERR4086637,ERR4086643,ERR4088531,ERR4086648,ERR4088538,ERR4088543,ERR4086653,ERR4560877,ERR4563427,ERR4560884,ERR4088549,ERR4086661,ERR4086669,ERR4563432,ERR4086674,ERR4086679,ERR4086687,ERR4086691,ERR4088557,ERR4563437,ERR4560889,ERR4560894,ERR4560898,ERR4086698,ERR4560903,ERR4088562,ERR4560908,ERR4560912,ERR4086705,ERR4563444,ERR4563448,ERR4086711,ERR4563449,ERR4560916,ERR4088568,ERR4088574,ERR4560923,ERR4563454,ERR4563458,ERR4560929,ERR4560934,ERR4086715,ERR4563463,ERR4560941,ERR4563467,ERR4086719,ERR4088581,ERR4563471,ERR4088585,ERR4563477,ERR4088591,ERR4086724,ERR4086730,ERR4563482,ERR4560946,ERR4086734,ERR4088595,ERR4560951,ERR4560956,ERR4088600,ERR4086739,ERR4560961,ERR4088605,ERR4560965,ERR4560970,ERR4086747,ERR4086755,ERR4086760,ERR4563487,ERR4563492,ERR4086765,ERR4563500,ERR4563509,ERR4563513,ERR4560975,ERR4563518,ERR4088609,ERR4560979,ERR4560986,ERR4560992,ERR4560997,ERR4088620,ERR4088622,ERR4086779,ERR4561001,ERR4086786,ERR4086794,ERR4086803,ERR4561009,ERR4086808,ERR4086814,ERR4088627,ERR4088631,ERR4561013,ERR4563523,ERR4561018,ERR4561025,ERR4563529,ERR4561029,ERR4088638,ERR4086822,ERR4086832,ERR4561037,ERR4086839,ERR4086849,ERR4088648,ERR4563535,ERR4561042,ERR4088653,ERR4086857,ERR4086864,ERR4088657,ERR4086869,ERR4086873,ERR4086879,ERR4086884,ERR4086889,ERR4561052,ERR4088662,ERR4561057,ERR4088665,ERR4086899,ERR4088669,ERR4088673,ERR4088682,ERR4563540,ERR4086903,ERR4563542,ERR4086907,ERR4563547,ERR4561062,ERR4086914,ERR4561067,ERR4086920,ERR4563552,ERR4086929,ERR4088685,ERR4086934,ERR4563557,ERR4086939,ERR4563561,ERR4088690,ERR4561074,ERR4561079,ERR4563570,ERR4561084,ERR4561089,ERR4561093,ERR4086946,ERR4086953,ERR4563575,ERR4563580,ERR4086958,ERR4561097,ERR4561101,ERR4563584,ERR4086966,ERR4561106,ERR4563593,ERR4086970,ERR4563597,ERR4086978,ERR4086988,ERR4563604,ERR4086990,ERR4561111,ERR4563609,ERR4086993,ERR4088695,ERR4563613,ERR4561116,ERR4087000,ERR4563617,ERR4087004,ERR4087009,ERR4087016,ERR4563625,ERR4561120,ERR4561125,ERR4561130,ERR4087019,ERR4561135,ERR4561140,ERR4087024,ERR4087029,ERR4563629,ERR4563633,ERR4561148,ERR4563638,ERR4561153,ERR4087034,ERR4563643,ERR4563647,ERR4088700,ERR4087039,ERR4087043,ERR4561157,ERR4563651,ERR4561162,ERR4087047,ERR4088702,ERR4563656,ERR4087057,ERR4561166,ERR4088711,ERR4088716,ERR4563660,ERR4087061,ERR4087065,ERR4561171,ERR4087069,ERR4561178,ERR4087077,ERR4563662,ERR4561184,ERR4561188,ERR4561192,ERR4561196,ERR4087083,ERR4561201,ERR4561206,ERR4087090,ERR4563663,ERR4561210,ERR4561214,ERR4087095,ERR4088718,ERR4088723,ERR4088729,ERR4561219,ERR4561223,ERR4087105,ERR4087110,ERR4563671,ERR4561227,ERR4088736,ERR4561231,ERR4563675,ERR4087115,ERR4087118,ERR4087124,ERR4087128,ERR4563679,ERR4087134,ERR4563683,ERR4561235,ERR4561242,ERR4563688,ERR4561246,ERR4563695,ERR4561251,ERR4561256,ERR4563700,ERR4561260,ERR4563705,ERR4563710,ERR4561264,ERR4561269,ERR4563714,ERR4561274,ERR4563718,ERR4561279,ERR4561285,ERR4563723,ERR4561290,ERR4087240,ERR4087249,ERR4561295,ERR4561301,ERR4563728,ERR4087253,ERR4563733,ERR4561305,ERR4561311,ERR4563739,ERR4561316,ERR4561321,ERR4561325,ERR4561330,ERR4561334,ERR4561339,ERR4563743,ERR4563746,ERR4563750,ERR4561343,ERR4561349,ERR4561353,ERR4563754,ERR4561361,ERR4563759,ERR4561367,ERR4561372,ERR4561377,ERR4561383,ERR4563763,ERR4561387,ERR4561391,ERR4561395,ERR4561399,ERR4561405,ERR4561409,ERR4563767,ERR4561415,ERR4563772,ERR4561421,ERR4561428,ERR4561432,ERR4561437,ERR4561442,ERR4561447,ERR4561449,ERR4563776,ERR4563781,ERR4561453,ERR4561461,ERR4563785,ERR4561465,ERR4563790,ERR4563794,ERR4561469,ERR4561474,ERR4561478,ERR4563798,ERR4563803,ERR4561483,ERR4561489,ERR4563806,ERR4561493,ERR4563810,ERR4561497,ERR4563815,ERR4561502,ERR4563819,ERR4561506,ERR4563823,ERR4561513,ERR4561518,ERR4561522,ERR4561527,ERR4563828,ERR4561532,ERR4561537,ERR4561541,ERR4563833,ERR4087258,ERR4561545,ERR4563837,ERR4561549,ERR4563843,ERR4561553,ERR4563849,ERR4563854,ERR4561558,ERR4561563,ERR4563858,ERR4563862,ERR4563867,ERR4563876,ERR4563886,ERR4563894,ERR4563899,ERR4563904,ERR4561571,ERR4563910,ERR4563919,ERR4563928,ERR4563932,ERR4563939,ERR4563943,ERR4563947,ERR4563953,ERR4563958,ERR4563961,ERR4563965,ERR4563967,ERR4563971]' --exact-matches --outdir batch_9 --apikey "${NCBI_API_KEY}" --sylph_db /work_beegfs/sukmb465/projects/TOFUpaper/sylph_db/gtdb-r220-c200-dbv1.syldb --sylph_processing
