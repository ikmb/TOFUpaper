#!/bin/bash
NXF_VER=23.10.1 ./nextflow run ikmb/tofu-maapo -r 6887e5c -profile custom -c tofu.config --sra '[ERR1137298,ERR1137299,ERR1137300,ERR1137310,ERR1137311,ERR1137312,ERR1137313,ERR1137314,ERR1137036,ERR1137315,ERR1137317,ERR1137324,ERR1137346,ERR1136697,ERR1136742,ERR1137348,ERR1137354,ERR1136743,ERR1137355,ERR1137062,ERR1137063,ERR1137064,ERR1137356,ERR1137357,ERR1137358,ERR1137359,ERR1137360,ERR1137361,ERR1136744,ERR1137362,ERR1137072,ERR1137668,ERR1137713,ERR1136634,ERR1136635,ERR1137580,ERR1137611,ERR1136745,ERR1136640,ERR1137581,ERR1136642,ERR1137582,ERR1137583,ERR1136645,ERR1137685,ERR1136647,ERR1137613,ERR1137584,ERR1136746,ERR1137585,ERR1136651,ERR1136652,ERR1137586,ERR1136654,ERR1137587,ERR1136656,ERR1136657,ERR1136658,ERR1136659,ERR1136747,ERR1136668,ERR1136669,ERR1136748,ERR1136670,ERR1136671,ERR1136672,ERR1136673,ERR1136674,ERR1137592,ERR1137593,ERR1136677,ERR1137594,ERR1137615,ERR1136749,ERR1136680,ERR1136681,ERR1136682,ERR1137595,ERR1137596,ERR1136685,ERR1136686,ERR1136687,ERR1136688,ERR1136689,ERR1136750,ERR1137637,ERR1137643,ERR1137499,ERR1137437,ERR1137438,ERR1137439,ERR1137644,ERR1137441,ERR1137442,ERR1137443,ERR1137444,ERR1136752,ERR1137445,ERR1137446,ERR1137447,ERR1137448,ERR1137714,ERR1137450,ERR1137451,ERR1137452,ERR1137453,ERR1137454,ERR1137700,ERR1137645,ERR1137619,ERR1137457,ERR1137458,ERR1137715,ERR1137460,ERR1137461,ERR1137646,ERR1137717,ERR1136754,ERR1137647,ERR1137374,ERR1137375,ERR1137376,ERR1137377,ERR1137378,ERR1137379,ERR1137380,ERR1137363,ERR1137364,ERR1136755,ERR1137431,ERR1137432,ERR1137433,ERR1137434,ERR1137435,ERR1137436,ERR1137365,ERR1137366,ERR1137648,ERR1136756,ERR1137466,ERR1137649,ERR1137650,ERR1137651,ERR1137652,ERR1137471,ERR1137472,ERR1137623,ERR1137474,ERR1137624,ERR1136757,ERR1137653,ERR1137485,ERR1136758,ERR1137654,ERR1137487,ERR1137488,ERR1137489,ERR1137490,ERR1137500,ERR1137501,ERR1137502,ERR1137503,ERR1137504,ERR1136759,ERR1137505,ERR1137382,ERR1137506,ERR1137507,ERR1137508,ERR1137509,ERR1137510,ERR1137511,ERR1137512,ERR1137513,ERR1136760,ERR1137383,ERR1137655,ERR1137385,ERR1137386,ERR1137387,ERR1137625,ERR1137389,ERR1137390,ERR1137391,ERR1137392,ERR1136761,ERR1137669,ERR1137535,ERR1137536,ERR1137537,ERR1137538,ERR1137539,ERR1137540,ERR1137541,ERR1137542,ERR1137543,ERR1136762,ERR1137544,ERR1137545,ERR1137546,ERR1137514,ERR1137515,ERR1137516,ERR1137517,ERR1137518,ERR1137519,ERR1137520,ERR1137521,ERR1137522,ERR1137523,ERR1137524,ERR1137525,ERR1137526,ERR1137527,ERR1137528,ERR1137547,ERR1136764,ERR1137559,ERR1137560,ERR1137561,ERR1137562,ERR1137563,ERR1137529,ERR1137530,ERR1136766,ERR1137566,ERR1137598,ERR1137567,ERR1137568,ERR1137599,ERR1137569,ERR1137600,ERR1137601,ERR1137626,ERR1137670,ERR1137656,ERR1137671,ERR1137672,ERR1137673,ERR1137628,ERR1137674,ERR1137630,ERR1137658,ERR1137659,ERR1137675,ERR1137676,ERR1137633,ERR1137634,ERR1137635,ERR1137677,ERR1137636,ERR1137660,ERR1137661,ERR1137662,ERR1137678,ERR1137679,ERR1137680,ERR1137681,ERR1137682,ERR1137683,ERR1137663,ERR1137684,ERR1137665,ERR1137686,ERR1137687,ERR1137688,ERR1137689,ERR1137690,ERR1137691,ERR1137692,ERR1137693,ERR1137694,ERR1137695,ERR1137696,ERR1137698,ERR1137726,ERR1137229,ERR1136772,ERR1136773,ERR1136774,ERR1136775,ERR1137235,ERR1136777,ERR1136778,ERR1136779,ERR1136780,ERR1136781,ERR1136693,ERR1136801,ERR1136802,ERR1136803,ERR1137702,ERR1137703,ERR1136832,ERR1136842,ERR1137327,ERR1136844,ERR1136845,ERR1136846,ERR1136847,ERR1136848,ERR1136849,ERR1136850,ERR1136875,ERR1136877,ERR1136878,ERR1136879,ERR1136880,ERR1137575,ERR1136882,ERR1136922,ERR1136923,ERR1136885,ERR1136886,ERR1136887,ERR1136888,ERR1136889,ERR1136890,ERR1136891,ERR1136892,ERR1136893,ERR1137576,ERR1136895,ERR1137577,ERR1136897,ERR1136898,ERR1136899,ERR1136900,ERR1137397,ERR1137398,ERR1137399,ERR1137400,ERR1137708,ERR1137401,ERR1137402,ERR1137404,ERR1137405,ERR1137211,ERR1137423,ERR1137275,ERR1137217,ERR1137276,ERR1137277,ERR1137278,ERR1137279,ERR1137280,ERR1137281,ERR1137237,ERR1136695,ERR1137285,ERR1137286,ERR1137121,ERR1137291,ERR1137294,ERR1137133,ERR1137134,ERR1137238,ERR1137136,ERR1137239,ERR1137138,ERR1137139,ERR1137140,ERR1137141,ERR1137142,ERR1137163,ERR1137164,ERR1137295,ERR1137168,ERR1137609,ERR1137170,ERR1137171,ERR1137172,ERR1137173,ERR1137296,ERR1137244,ERR1137176,ERR1137177,ERR1137178,ERR1137179,ERR1136618,ERR1136619,ERR1136620,ERR1136629,ERR1136732,ERR1137711,ERR1136733,ERR1136734,ERR1136971,ERR1136972,ERR1136973,ERR1136974,ERR1136975,ERR1136976,ERR1137245,ERR1137246,ERR1137247,ERR1137248,ERR1136735,ERR1137699,ERR1136737,ERR1137270,ERR1136738,ERR1136739,ERR1137301,ERR1137302,ERR1137303,ERR1137024,ERR1137304,ERR1137305,ERR1137306,ERR1137307,ERR1137308,ERR1137309,ERR1136740,ERR1137316,ERR1137318,ERR1136741,ERR1137319,ERR1137320,ERR1137043,ERR1137321,ERR1137045,ERR1137322,ERR1137323,ERR1137048,ERR1137347,ERR1137349,ERR1137325,ERR1137350,ERR1137351,ERR1137352,ERR1137353,ERR1137059,ERR1136637,ERR1136639,ERR1137588,ERR1137589,ERR1137590,ERR1137591,ERR1136664,ERR1136665,ERR1137614,ERR1136667,ERR1137642,ERR1137616,ERR1136692,ERR1137617,ERR1137492,ERR1137493,ERR1137494,ERR1137495,ERR1137496,ERR1137497,ERR1136698,ERR1137716,ERR1137718,ERR1137476,ERR1137478,ERR1137479,ERR1137480,ERR1137481,ERR1137482,ERR1137483,ERR1137484,ERR1136699,ERR1137701,ERR1137548,ERR1137549,ERR1137550,ERR1137551,ERR1137552,ERR1137553,ERR1137554,ERR1137555,ERR1137556,ERR1137597,ERR1137558,ERR1136765,ERR1137531,ERR1137532,ERR1137533,ERR1136767,ERR1137631,ERR1136768,ERR1136769,ERR1136770,ERR1136771,ERR1137697,ERR1137719,ERR1137720,ERR1137721,ERR1137722,ERR1137723,ERR1137724,ERR1137725,ERR1137727,ERR1136701,ERR1136782,ERR1136783,ERR1136784,ERR1136785,ERR1136786,ERR1136787,ERR1136788,ERR1136789,ERR1136790,ERR1136791,ERR1110474,ERR1110309,ERR1110337,ERR1110338,ERR1110310,ERR1110311,ERR1110312,ERR1110313,ERR1110339,ERR1110314,ERR1110315,ERR1110318,ERR1110438,ERR1110362,ERR1110443,ERR1110444,ERR1110445,ERR1110446,ERR1110447,ERR1110454,ERR1110473,ERR1110279,ERR1110302,ERR1110303,ERR1110334,ERR1110304,ERR1110305,ERR1110306,ERR1110466,ERR1110335,ERR1110336,ERR1110308,ERR1110331,ERR1110332,ERR1110316,ERR1110317,ERR1110319,ERR1110340,ERR1110341,ERR1110320,ERR1110321,ERR1110322,ERR1110323,ERR1110297,ERR1110324,ERR1110325,ERR1110326,ERR1110327,ERR1110328,ERR1110467,ERR1110330,ERR1110342,ERR1110352,ERR1110433,ERR1110298,ERR1110434,ERR1110468,ERR1110353,ERR1110343,ERR1110344,ERR1110436,ERR1110437,ERR1110355,ERR1110345,ERR1110346,ERR1110299,ERR1110347,ERR1110439,ERR1110440,ERR1110348,ERR1110349,ERR1110350,ERR1110358,ERR1110351,ERR1110359,ERR1110333,ERR1110441,ERR1110360,ERR1110361,ERR1110442,ERR1110300,ERR1110448,ERR1110449,ERR1110450,ERR1110363,ERR1110364,ERR1110451,ERR1110366,ERR1110452,ERR1110453,ERR1110301,ERR1110469,ERR1110470,ERR1110471,ERR1110458,ERR1110459,ERR1110472,ERR1110462,ERR1110463,ERR1110464,ERR480457,ERR480459,ERR480463,ERR480465,ERR480469,ERR480471,ERR480475,ERR480479,ERR480480,ERR480485,ERR480489,ERR480493,ERR480495,ERR480499,ERR480501,ERR480505,ERR480507,ERR480509,ERR480513,ERR480516,ERR480519,ERR480524,ERR480528,ERR480532,ERR480534,ERR480536,ERR480538,ERR480540,ERR480541,ERR480548,ERR480552,ERR480553,ERR480556,ERR480560,ERR480562,ERR480566,ERR480569,ERR480571,ERR480573,ERR480577,ERR480579,ERR480582,ERR480585,ERR480588,ERR480592,ERR480594,ERR480596,ERR480598,ERR480600,ERR480602,ERR480606,ERR480609,ERR480612,ERR480615,ERR480618,ERR480622,ERR480627,ERR480629,ERR480632,ERR480635,ERR480637,ERR480640,ERR480643,ERR480644,ERR480649,ERR480651,ERR480654,ERR480657,ERR480658,ERR480662,ERR480664,ERR480668,ERR480670,ERR480672,ERR480673,ERR480678,ERR480682,ERR480685,ERR480687,ERR480689,ERR480692,ERR480696,ERR480698,ERR480701,ERR480703,ERR480705,ERR480712,ERR480714,ERR480717,ERR480720,ERR480722,ERR480726,ERR480728,ERR480733,ERR480736,ERR480737,ERR480740,ERR480744,ERR480746,ERR480749,ERR480750,ERR480754,ERR480758,ERR480760,ERR480765,ERR480769,ERR480771,ERR480774,ERR480779,ERR480781,ERR480784,ERR480789,ERR480790,ERR480794,ERR480797,ERR480798,ERR480803,ERR480807,ERR480811,ERR480813,ERR480816,ERR480820,ERR480822,ERR480824,ERR480826,ERR480830,ERR480834,ERR480836,ERR480838,ERR480841,ERR480842,ERR480846,ERR480848,ERR480850,ERR480852,ERR480854,ERR480857,ERR480861,ERR480865,ERR480869,ERR480870,ERR480877,ERR480878,ERR480881,ERR480884,ERR480885,ERR480888,ERR480889,ERR480892,ERR480895,ERR480900,ERR480904,ERR480908,ERR480909,ERR480911,ERR480914,ERR2855866,ERR2855940,ERR2855878,ERR2855895,ERR2855822,ERR2855862,ERR2855881,ERR2855801,ERR2855918,ERR2855950,ERR2855927,ERR2855842,ERR2855879,ERR2855847,ERR2855807,ERR2855912,ERR2855872,ERR2855936,ERR2855946,ERR2855827,ERR2855805,ERR2855861,ERR2855860,ERR2855919,ERR2855865,ERR2855885,ERR2855906,ERR2855938,ERR2855835,ERR2855811,ERR2855806,ERR2855902,ERR2855853,ERR2855849,ERR2855848,ERR2855896,ERR2855921,ERR2855933,ERR2855914,ERR2855839,ERR2855929,ERR2855792,ERR2855804,ERR2855944,ERR2855953,ERR2855812,ERR2855901,ERR2855934,ERR2855892,ERR2855870,ERR2855851,ERR2855856,ERR2855907,ERR2855786,ERR2855803,ERR2855877,ERR2855844,ERR2855864,ERR2855794,ERR2855899,ERR2855880,ERR2855857,ERR2855826,ERR2855924,ERR2855814,ERR2855882,ERR2855858,ERR2855852,ERR2855863,ERR2855875,ERR2855926,ERR2855939,ERR2855789,ERR2855840,ERR2855949,ERR2855916,ERR2855836,ERR2855843,ERR2855824,ERR2855788,ERR2855954,ERR2855819,ERR2855890,ERR2855859,ERR2855813,ERR2855931,ERR2855797,ERR2855911,ERR2855952,ERR2855831,ERR2855855,ERR2855796,ERR2855893,ERR2855917,ERR2855928,ERR2855915,ERR2855790,ERR2855845,ERR2855888,ERR2855898,ERR2855910,ERR2855854,ERR2855920,ERR2855876,ERR2855820,ERR2855850,ERR2855947,ERR2855817,ERR2855791,ERR2855841,ERR2855923,ERR2855810,ERR2855809,ERR2855922,ERR2855884,ERR2855823,ERR2855828,ERR2855945,ERR2855904,ERR2855816,ERR2855838,ERR2855935,ERR2855909,ERR2855815,ERR2855834,ERR2855833,ERR2855937,ERR2855925,ERR2855874,ERR2855868,ERR2855951,ERR2855825,ERR2855873,ERR2855932,ERR2855942,ERR2855837,ERR2855913,ERR2855887,ERR2855897,ERR2855889,ERR2855818,ERR2855894,ERR2855799,ERR2855802,ERR2855905,ERR2855891,ERR2855871,ERR2855832,ERR2855867,ERR2855908,ERR2855941,ERR2855798,ERR2855787,ERR2855886,ERR2855955,ERR2855943,ERR2855829,ERR2855821,ERR2855948,ERR2855900,ERR2855808,ERR2855903,ERR2855930,ERR2855846,ERR2855956,ERR2855795,ERR2855793,ERR2855883,ERR2855800,ERR2855869,ERR2855830]' --exact-matches --outdir batch_18 --apikey "${NCBI_API_KEY}" --sylph_db /work_beegfs/sukmb465/projects/TOFUpaper/sylph_db/gtdb-r220-c200-dbv1.syldb --sylph_processing
