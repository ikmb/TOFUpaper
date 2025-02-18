#!/bin/bash
NXF_VER=23.10.1 ./nextflow run ikmb/tofu-maapo -r 6887e5c -profile custom -c tofu.config --sra '[SRR8146992,SRR8146994,SRR6367590,SRR4074352,SRR4074353,SRR8146959,SRR8146991,SRR6367591,SRR8146960,SRR6367592,SRR6367598,SRR4074355,SRR8146965,SRR4074356,SRR4074357,SRR6367597,SRR4074261,SRR4074262,SRR4074263,SRR8146961,SRR8146954,SRR4074264,SRR6367599,SRR4074265,SRR6367600,SRR6367587,SRR4074266,SRR4074267,SRR8146955,SRR4074269,SRR4074278,SRR4074296,SRR4074305,SRR8146993,SRR6367588,SRR4074323,SRR4074332,SRR6367593,SRR6367594,SRR4074341,SRR6367595,SRR6367596,SRR8146938,SRR8146935,SRR8146936,SRR8146962,SRR8146987,SRR8146958,SRR8146974,SRR8146973,SRR8146988,SRR8146978,SRR8146966,SRR8146964,SRR8146972,SRR8146971,SRR8146977,SRR8146947,SRR8146943,SRR8146953,SRR8146970,SRR8146950,SRR8146969,SRR8146980,SRR8146979,SRR8146945,SRR8146986,SRR8146946,SRR8146944,SRR8146963,SRR8146968,SRR8146967,SRR8146981,SRR8146941,SRR8146976,SRR8146942,SRR8146939,SRR8146982,SRR8146975,SRR8146983,SRR8146940,SRR8146984,SRR8146985,SRR8146948,SRR8146949,SRR8146989,SRR8146990,SRR5898908,SRR5898909,SRR5898910,SRR5898911,SRR5898912,SRR5898913,SRR5898914,SRR5898915,SRR5898916,SRR5898917,SRR5898918,SRR5898919,SRR5898920,SRR5898921,SRR5898922,SRR5898923,SRR5898924,SRR5898925,SRR5898926,SRR5898927,SRR5898928,SRR5898929,SRR5898930,SRR5898931,SRR5898932,SRR5898933,SRR5898934,SRR5898935,SRR5898936,SRR5898937,SRR5898938,SRR5898939,SRR5898940,SRR5898941,SRR5898942,SRR5898943,SRR5898944,SRR5898945,SRR5898946,SRR5898947,SRR5898948,SRR5898949,SRR5898950,SRR5898951,SRR5898952,SRR5898953,SRR5898954,SRR5898955,SRR5898956,SRR5898957,SRR5898958,SRR5898959,SRR5898960,SRR5898961,SRR5898962,SRR5898963,SRR5898964,SRR5898965,SRR5898966,SRR5898967,SRR5898968,SRR5898969,SRR5898970,SRR5898971,SRR5898972,SRR5898973,SRR5898974,SRR5898975,SRR5898976,SRR5898977,SRR5898978,SRR5898979,SRR5898980,SRR5898981,SRR5898982,SRR5898983,SRR5898984,SRR5898985,SRR5898986,SRR5898987,SRR5898988,SRR5898989,SRR5898990,SRR5898991,SRR5898992,SRR5898993,SRR5898994,SRR5898995,SRR5898996,SRR5898997,SRR5898998,SRR5898999,SRR5899000,SRR5899001,SRR5899002,SRR5899003,SRR5899004,SRR5899005,SRR5899006,SRR5899007,SRR5899008,SRR5899009,SRR5899010,SRR5899011,SRR5899012,SRR5899013,SRR5899014,SRR5899015,SRR5899016,SRR5899017,ERR688505,ERR688506,ERR688507,ERR688508,ERR688509,ERR688510,ERR688511,ERR688512,ERR688513,ERR688514,ERR688515,ERR688516,ERR688517,ERR688519,ERR688520,ERR688521,ERR688522,ERR688523,ERR688524,ERR688525,ERR688526,ERR688527,ERR688528,ERR688529,ERR688530,ERR688531,ERR688532,ERR688533,ERR688534,ERR688535,ERR688536,ERR688537,ERR688538,ERR688539,ERR688540,ERR688541,ERR688542,ERR688543,ERR688544,ERR688545,ERR688546,ERR688547,ERR688548,ERR688549,ERR688550,ERR688551,ERR688552,ERR688553,ERR688554,ERR688555,ERR688557,ERR688558,ERR688559,ERR688560,ERR688561,ERR688562,ERR688563,ERR688564,ERR688565,ERR688566,ERR688568,ERR688569,ERR688570,ERR688571,ERR688572,ERR688573,ERR688574,ERR688575,ERR688576,ERR688577,ERR688578,ERR688579,ERR688580,ERR688581,ERR688582,ERR688583,ERR688567,ERR688587,ERR688584,ERR688585,ERR688586,ERR688588,ERR688589,ERR688590,ERR688591,ERR688592,ERR688593,ERR688594,ERR688595,ERR688596,ERR688597,ERR688598,ERR688599,ERR688600,ERR688601,ERR688602,ERR688603,ERR688604,ERR688605,ERR688606,ERR688607,ERR688608,ERR688609,ERR688610,ERR688611,ERR688612,ERR688613,ERR688614,ERR688615,ERR688616,ERR688617,ERR688618,ERR688619,ERR688620,ERR688621,ERR688622,ERR688623,ERR688624,ERR688625,ERR688626,ERR688627,ERR688628,ERR688629,ERR688630,ERR688631,ERR688632,ERR688633,ERR688634,ERR688635,ERR688636,ERR688637,ERR688638,ERR688639,ERR688640,ERR688641,ERR688642,ERR688643,ERR688644,ERR688645,ERR688646,ERR688647,ERR688648,ERR688649,ERR710424,ERR688650,ERR688651,ERR710425,ERR710426,ERR710427,ERR710428,ERR710429,ERR710430,ERR710431,ERR710432,SRR5274093,SRR5274092,SRR5274091,SRR5274089,SRR5274088,SRR5274087,SRR5274086,SRR5274085,SRR5274084,SRR5274083,SRR5274082,SRR5274081,SRR5274080,SRR5274079,SRR5274078,SRR5274077,SRR5274076,SRR5274075,SRR5274074,SRR5274073,SRR5274072,SRR5274071,SRR5274070,SRR5274069,SRR5274067,SRR5274066,SRR5274065,SRR5274064,SRR5274063,SRR5274062,SRR5274061,SRR5274060,SRR5274059,SRR5274058,SRR5274057,SRR5274056,SRR5274055,SRR5274054,SRR5274053,SRR5274052,SRR5274051,SRR5274050,SRR5274049,SRR5274048,SRR5274047,SRR5274046,SRR5274045,SRR5274044,SRR5274043,SRR5274042,SRR5274041,SRR5274040,SRR5274039,SRR5274038,SRR5274037,SRR5274036,SRR5274035,SRR5274034,SRR5274033,SRR5274032,SRR5274031,SRR5274030,SRR5274029,SRR5274028,SRR5274027,SRR5274026,SRR5274025,SRR5274024,SRR5274023,SRR5274022,SRR5274021,SRR5274020,SRR5274019,SRR5274018,SRR5274017,SRR5274016,SRR5274015,SRR5274014,SRR5274013,SRR5274012,SRR5274011,SRR5274010,SRR5274009,SRR5274008,SRR5274007,SRR5274006,SRR5274005,SRR5274004,SRR5274003,SRR5274002,SRR5274001,SRR5274000,SRR5273999,SRR5273998,SRR5273997,SRR5273996,SRR5273995,SRR5273994,SRR5273993,SRR5273992,SRR5273991,SRR5273990,SRR5273989,SRR5273988,SRR5273987,SRR5273986,SRR5273985,SRR5273984,SRR5273983,SRR5273982,SRR5273981,SRR5273980,SRR5273979,SRR5273978,SRR5273977,SRR5273976,SRR5273975,SRR5273974,SRR5273973,SRR5273972,SRR5273971,SRR5273970,SRR5273969,SRR5273968,SRR5273967,SRR5273966,SRR5273965,SRR5273964,SRR5273963,SRR5273962,SRR5273961,SRR5273960,SRR5273959,SRR5273958,SRR5273957,SRR5273956,SRR5273955,SRR5273954,SRR5273953,SRR5273952,SRR5273951,SRR5273950,SRR5273949,SRR5273948,SRR5273947,SRR5273946,SRR5273945,SRR5273944,SRR5273943,SRR5273942,SRR5273941,SRR5273940,SRR5273939,SRR5273938,SRR5273937,SRR5273936,SRR5273935,SRR5273934,SRR5273933,SRR5273932,SRR5273931,SRR5273930,SRR5273929,SRR5273928,SRR5273927,SRR5273926,SRR5273925,SRR5273924,SRR5273923,SRR5273922,SRR5273921,SRR5273920,SRR5273919,SRR5273918,SRR5273917,SRR5273916,SRR5273915,SRR5273914,SRR5273913,SRR5273912,SRR5273911,SRR5273910,SRR5273909,SRR5273908,SRR5273907,SRR5273906,SRR5273905,SRR5273904,SRR5273903,SRR5273902,SRR5273901,SRR5273900,SRR5273899,SRR5273898,SRR5273897,SRR5273896,SRR5273895,SRR5273894,SRR5273893,SRR5273892,SRR5273891,SRR5273890,SRR5273889,SRR5273888,SRR5273887,SRR5273886,SRR5273885,SRR5273884,SRR5273883,SRR5273882,SRR5273881,SRR5273880,SRR5273879,SRR5273878,SRR5930495,SRR5930497,SRR5930498,SRR5930499,SRR5930501,SRR5930502,SRR5930493,SRR5930494,SRR5930526,SRR5930527,SRR5930528,SRR5930521,SRR5930522,SRR5930523,SRR5930533,SRR5930534,SRR5930511,SRR5930510,SRR5930513,SRR5930512,SRR5930515,SRR5930514,SRR5930517,SRR5930516,SRR5930509,SRR5930508,SRR5930536,SRR5930520,SRR5930503,SRR5930531,SRR5930532,SRR5930529,SRR5930530,SRR5930518,SRR5930519,SRR5930507,SRR5930506,SRR5930505,SRR5930504,SRR9217428,SRR9217406,SRR9217452,SRR9217391,SRR9217414,SRR9217393,SRR9217457,SRR9217463,SRR9217392,SRR9217451,SRR9217419,SRR9217422,SRR9217492,SRR9217400,SRR9217389,SRR9217415,SRR9217399,SRR9217408,SRR9217423,SRR9217471,SRR9217475,SRR9217449,SRR9217435,SRR9217454,SRR9217402,SRR9217473,SRR9217490,SRR9217401,SRR9217433,SRR9217404,SRR9217396,SRR9217487,SRR9217489,SRR9217479,SRR9217464,SRR9217462,SRR9217403,SRR9217470,SRR9217388,SRR9217394,SRR9217438,SRR9217478,SRR9217416,SRR9217461,SRR9217472,SRR9217430,SRR9217450,SRR9217459,SRR9217474,SRR9217384,SRR9217431,SRR9217484,SRR9217437,SRR9217456,SRR9217441,SRR9217434,SRR9217483,SRR9217446,SRR9217494,SRR9217460,SRR9217420,SRR9217398,SRR9217413,SRR9217476,SRR9217440,SRR9217468,SRR9217397,SRR9217444,SRR9217465,SRR9217425,SRR9217442,SRR9217417,SRR9217410,SRR9217485,SRR9217405,SRR9217385,SRR9217481,SRR9217436,SRR9217480,SRR9217491,SRR9217455,SRR9217467,SRR9217469,SRR9217448,SRR9217412,SRR9217458,SRR9217493,SRR9217387,SRR9217421,SRR9217409,SRR9217383,SRR9217407,SRR9217495,SRR9217395,SRR9217443,SRR9217429,SRR9217386,SRR9217418,SRR9217466,SRR9217486,SRR9217445,SRR9217427,SRR9217411,SRR9217488,SRR9217482,SRR9217390,SRR9217447,SRR9217424,SRR9217477,SRR9217426,SRR9217453,SRR9217432,SRR9217439,ERR2162224,ERR2162223,ERR2162222,ERR2162221,ERR2162220,ERR2162219,ERR2162218,ERR2162217,ERR2162216,ERR2162215,ERR2162214,ERR2162213,ERR2162212,ERR2162211,ERR2162210,ERR2162209,ERR2162208,ERR2162207,ERR2162206,ERR2162205,ERR2162204,ERR2162203,ERR2162202,ERR2162201,ERR2162200,SRR8865580,SRR8865577,SRR8865593,SRR8865575,SRR8865594,SRR8865573,SRR8865579,SRR8865586,SRR8865597,SRR8865600,SRR8865591,SRR8865574,SRR8865578,SRR8865581,SRR8865598,SRR8865583,SRR8865601,SRR8865572,SRR8865595,SRR8865587,SRR8865592,SRR8865599,SRR8865584,SRR8865596,SRR8865576,SRR8865582,SRR8865588,SRR8865585,SRR8865590,SRR8865589,SRR5651391,SRR5651392,SRR5651393,SRR5651394,SRR5651395,SRR5651396,SRR5651397,SRR5651398,SRR5651400,SRR5651401,SRR5651402,SRR5651403,SRR5651404,SRR5651405,SRR5651406,SRR5651407,SRR5651408,SRR5651409,SRR5651410,SRR5651411,SRR5651412,SRR5651413,SRR5651414,SRR5651415,SRR5651416,SRR5651417,SRR5651418,SRR5651419,SRR5651420,SRR5651421,SRR5651422,SRR5651423,SRR5651424,SRR5651425,SRR5651426,SRR5651427,SRR5651428,SRR5651429,SRR5651430,SRR5651431,SRR5651432,SRR5651433,SRR5651434,SRR5651435,SRR5651436,SRR5651437,SRR5651438,SRR5651439,SRR5651440,SRR5651441,SRR5651442,SRR5651443,SRR5651444,SRR5651445,SRR5651446,SRR5651447,SRR5651448,SRR5651449,SRR5651450,SRR5651451,SRR5651452,SRR5651453,SRR5651455,SRR5651456,SRR5651457,SRR5651458,SRR5651459,SRR5651460,SRR5651461,SRR5651462,SRR5651463,SRR5651464,SRR5651465,SRR5651466,SRR5651467,SRR5651468,SRR5651469,SRR5651470,SRR5651471,SRR5651472,SRR5651473,SRR5651474,SRR5651475,SRR5651476,SRR5651478,SRR5651479,SRR5651480,SRR5651481,SRR5651482,SRR5650021,SRR5650022,SRR5650023,SRR5650024,SRR5650025,SRR5650026,SRR5650027,SRR5650028,SRR5650029,SRR5650030,SRR5650031,SRR5650032,SRR5650033,SRR5650034,SRR5650035,SRR5650036,SRR5650037,SRR5650038,SRR5650039,SRR5650040,SRR5650041,SRR5650042,SRR5650043,SRR5650044,SRR5650045,SRR5650046,SRR5650047,SRR5650048,SRR5650049,SRR5650050,SRR5650051,SRR5650052,SRR5650053,SRR5650054,SRR5650055,SRR5650056,SRR5650057,SRR5650058,SRR5650059,SRR5650060,SRR5650061,SRR5650062,SRR5650063,SRR5650064,SRR5650065,SRR5650066,SRR5650067,SRR5650068,SRR5650069,SRR5650070,SRR5650071,SRR5650072,SRR5650073,SRR5650074,SRR5650075,SRR5650076,SRR5650077,SRR5650078,SRR5650079,SRR5650080,SRR5650081,SRR5650082,SRR5650083,SRR5650084,SRR5650085,SRR5650086,SRR5650087,SRR5650088,SRR5650089,SRR5650090,SRR5650091,SRR5650092,SRR5650093,SRR5650094,SRR5650095,SRR5650096,SRR5650097,SRR5650098,SRR5650099,SRR5650100,SRR5650101,SRR5650102,SRR5650103,SRR5650104,SRR5650105,SRR5650106,SRR5650107,SRR5650108,SRR5650109,SRR5650110,SRR5650111,SRR5650112,SRR5650113,SRR5650114,SRR5650115,SRR5650116,SRR5650117,SRR5650118,SRR5650119,SRR5650120,SRR5650121,SRR5650122,SRR5650123,SRR5650124,SRR5650125,SRR5650126,SRR5650127,SRR5650128,SRR5650129,SRR5650130,SRR5650131,SRR5650132,SRR5650133,SRR5650134,SRR5650135,SRR5650136,SRR5650137,SRR5650138,SRR5650139,SRR5650140,SRR5650141,SRR5650142,SRR5650143,SRR5650144,SRR5650145,SRR5650146,SRR5650147,SRR5650148,SRR5650149,SRR5650150,SRR5650151,SRR5650152,SRR5650153,SRR5650154,SRR5650155,SRR5650156,SRR5650157,SRR5650158,SRR5650159]' --exact-matches --outdir batch_4 --apikey "${NCBI_API_KEY}" --sylph_db /work_beegfs/sukmb465/projects/TOFUpaper/sylph_db/gtdb-r220-c200-dbv1.syldb --sylph_processing
