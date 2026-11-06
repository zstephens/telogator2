import copy
import numpy as np
import os
import random
import re
import shutil

LEXICO_2_IND = {'chr1':1, 'chr2':2, 'chr3':3, 'chr10':10, 'chr11':11, 'chr12':12, 'chr19':19, 'chr20':20,
                'chr4':4, 'chr5':5, 'chr6':6, 'chr13':13, 'chr14':14, 'chr15':15, 'chr21':21, 'chr22':22,
                'chr7':7, 'chr8':8, 'chr9':9, 'chr16':16, 'chr17':17, 'chr18':18,
                'chrX':23,'chrY':24, 'chrM':25, 'chrF':26, 'chrB':27, 'chrU':28}

SORTED_CHR_LIST = [n[1] for n in sorted([(LEXICO_2_IND[k],k) for k in LEXICO_2_IND])]

T2T_CHROMSIZE = {'chr1':248387328,
                 'chr2':242696752,
                 'chr3':201105948,
                 'chr4':193574945,
                 'chr5':182045439,
                 'chr6':172126628,
                 'chr7':160567428,
                 'chr8':146259331,
                 'chr9':150617247,
                 'chr10':134758134,
                 'chr11':135127769,
                 'chr12':133324548,
                 'chr13':113566686,
                 'chr14':101161492,
                 'chr15':99753195,
                 'chr16':96330374,
                 'chr17':84276897,
                 'chr18':80542538,
                 'chr19':61707364,
                 'chr20':66210255,
                 'chr21':45090682,
                 'chr22':51324926,
                 'chrX':154259566,
                 'chrY':62460029}

LARGE_NUMBER = int(1e9)

RC_DICT = {'A':'T','C':'G','G':'C','T':'A','N':'N'}

FILE_SUFFIX = {'.fa':       'fasta',
               '.fasta':    'fasta',
               '.fa.gz':    'fasta',
               '.fasta.gz': 'fasta',
               '.fq':       'fastq',
               '.fastq':    'fastq',
               '.fq.gz':    'fastq',
               '.fastq.gz': 'fastq',
               '.bam':      'bam',
               '.cram':     'cram',
               '.p':        'pickle'}

REF_CHAR  = 'MX=D'
READ_CHAR = 'MX=I'
CLIP_CHAR = 'SH'

BLANK_CHR   = 'chrBq'
UNCLUST_CHR = 'chrUq'
UNCLUST_POS = 0

SMALL_NUMBER = 1e-6

INTERSTITIAL_TEL_REGIONS = [('chm13_chr2q',         51395,  52589,  'p-1194'),
                            ('chm13_chr3q',         47928,  49762,  'q-1834'),
                            ('chm13_chr4q',         20760,  24132,  'q-3372'),
                            ('chm13_chr5p',         40612,  41197,  'p-585'),
                            ('chm13_chr5p',         82715,  83331,  'q-616'),
                            ('chm13_chr7p',         113853, 114611, 'p-758'),
                            ('chm13_chr9q',         132277, 132996, 'q-719'),
                            ('chm13_chr10q',        9721,   10735,  'q-1014'),
                            ('chm13_chr11p',        201399, 202732, 'p-1333'),
                            ('chm13_chr11p',        236575, 237735, 'q-1160'),
                            ('chm13_chr12p',        193038, 193985, 'p-947'),
                            ('chm13_chr13q',        46628,  47839,  'q-1211'),
                            ('chm13_chr16q',        13379,  14408,  'p-1029'),
                            ('chm13_chr16q',        129363, 130082, 'q-719'),
                            ('chm13_chr17p',        46080,  46795,  'q-715'),
                            ('chm13_chr17p',        239337, 239738, 'q-401'),
                            ('chm13_chr17q',        113842, 114989, 'q-1147'),
                            ('chm13_chr17q',        122492, 123392, 'q-900'),
                            ('chm13_chr18p',        55369,  56037,  'p-668'),
                            #('chm13_chr18p',        109389, 110031, 'p-642'),
                            #('chm13_chr18p',        163440, 164072, 'p-632'),
                            ('chm13_chr18p',        217347, 217895, 'p-548'),
                            ('chm13_chr18q',        242147, 242638, 'q-491'),
                            ('chm13_chr19p',        196736, 197191, 'p-455'),
                            ('chm13_chr20p',        70733,  71134,  'p-401'),
                            ('chm13_chr20q',        146949, 147785, 'q-836'),
                            ('chm13_chr22q',        121287, 123128, 'q-1841'),
                            ('hg002-mat_chr2q',     43000,  44207,  'p-1207'),
                            ('hg002-mat_chr3q',     51623,  53438,  'q-1815'),
                            ('hg002-mat_chr5p',     41100,  41696,  'p-596'),
                            ('hg002-mat_chr7p',     112489, 113241, 'p-752'),
                            ('hg002-mat_chr9q',     52787,  53530,  'q-743'),
                            ('hg002-mat_chr10q',    14566,  15644,  'q-1078'),
                            ('hg002-mat_chr11p',    201461, 202721, 'p-1260'),
                            ('hg002-mat_chr11p',    227606, 228926, 'q-1320'),
                            ('hg002-mat_chr12p',    94816,  95624,  'q-808'),
                            ('hg002-mat_chr12p',    184495, 185517, 'p-1022'),
                            ('hg002-mat_chr13q',    46372,  47827,  'q-1455'),
                            ('hg002-mat_chr16q',    14345,  14965,  'p-620'),
                            ('hg002-mat_chr16q',    129691, 130425, 'q-734'),
                            ('hg002-mat_chr17p',    41487,  42586,  'q-1099'),
                            ('hg002-mat_chr17q',    114813, 115430, 'q-617'),
                            ('hg002-mat_chr17q',    122952, 123942, 'q-990'),
                            ('hg002-mat_chr18p',    55734,  56259,  'p-525'),
                            ('hg002-mat_chr18p',    97321,  97779,  'p-458'),
                            ('hg002-mat_chr18q',    243086, 243577, 'q-491'),
                            ('hg002-mat_chr19p',    196005, 196448, 'p-443'),
                            ('hg002-mat_chr20p',    69721,  70153,  'p-432'),
                            ('hg002-mat_chr20q',    92361,  92927,  'p-566'),
                            ('hg002-mat_chr22q',    120441, 122600, 'q-2159'),
                            ('hg002-pat_chr2q',     50480,  51689,  'p-1209'),
                            ('hg002-pat_chr3q',     52029,  53852,  'q-1823'),
                            ('hg002-pat_chr4q',     43057,  46971,  'q-3914'),
                            ('hg002-pat_chr5p',     40807,  41383,  'p-576'),
                            ('hg002-pat_chr7p',     112511, 113268, 'p-757'),
                            ('hg002-pat_chr9q',     53466,  54204,  'q-738'),
                            ('hg002-pat_chr10q',    33047,  33511,  'q-464'),
                            ('hg002-pat_chr10q',    62960,  63969,  'q-1009'),
                            ('hg002-pat_chr11p',    202419, 203649, 'p-1230'),
                            ('hg002-pat_chr11p',    228679, 229970, 'q-1291'),
                            ('hg002-pat_chr12p',    189938, 190875, 'p-937'),
                            ('hg002-pat_chr13q',    44340,  46864,  'q-2524'),
                            ('hg002-pat_chr16p',    116727, 117460, 'p-733'),
                            ('hg002-pat_chr16q',    26926,  27679,  'p-753'),
                            ('hg002-pat_chr16q',    28650,  29107,  'p-457'),
                            ('hg002-pat_chr16q',    128566, 129027, 'p-461'),
                            ('hg002-pat_chr16q',    243861, 244645, 'q-784'),
                            ('hg002-pat_chr17p',    44295,  45387,  'q-1092'),
                            ('hg002-pat_chr17p',    239946, 240346, 'q-400'),
                            ('hg002-pat_chr17q',    113867, 115016, 'q-1149'),
                            ('hg002-pat_chr17q',    122505, 123491, 'q-986'),
                            ('hg002-pat_chr18p',    56033,  56554,  'p-521'),
                            ('hg002-pat_chr18p',    97568,  98019,  'p-451'),
                            ('hg002-pat_chr18q',    243433, 243928, 'q-495'),
                            ('hg002-pat_chr19p',    192733, 193194, 'p-461'),
                            ('hg002-pat_chr20p',    70137,  70587,  'p-450'),
                            ('hg002-pat_chr20q',    147296, 148128, 'q-832'),
                            ('hg002-pat_chr22q',    120123, 121891, 'q-1768'),
                            ('cn1-mat_chr2q',       46850,  48003,  'p-1153'),
                            ('cn1-mat_chr3q',       45778,  47602,  'q-1824'),
                            ('cn1-mat_chr4q',       81454,  83712,  'q-2258'),
                            ('cn1-mat_chr5p',       46621,  47215,  'p-594'),
                            ('cn1-mat_chr6p',       101529, 102058, 'p-529'),
                            ('cn1-mat_chr7p',       117925, 118675, 'p-750'),
                            ('cn1-mat_chr9q',       132915, 133686, 'q-771'),
                            ('cn1-mat_chr10q',      21642,  22154,  'q-512'),
                            ('cn1-mat_chr10q',      51593,  52569,  'q-976'),
                            ('cn1-mat_chr11p',      78385,  78810,  'p-425'),
                            ('cn1-mat_chr12p',      194550, 195504, 'p-954'),
                            ('cn1-mat_chr13q',      47799,  49233,  'q-1434'),
                            ('cn1-mat_chr16q',      17964,  18810,  'p-846'),
                            ('cn1-mat_chr16q',      132919, 133687, 'q-768'),
                            ('cn1-mat_chr17p',      118752, 119527, 'p-775'),
                            ('cn1-mat_chr17p',      156080, 156728, 'q-648'),
                            ('cn1-mat_chr17q',      113097, 114130, 'q-1033'),
                            ('cn1-mat_chr17q',      121347, 122947, 'q-1600'),
                            ('cn1-mat_chr18p',      62714,  63248,  'p-534'),
                            ('cn1-mat_chr18p',      104376, 104830, 'p-454'),
                            ('cn1-mat_chr18q',      239031, 239525, 'q-494'),
                            ('cn1-mat_chr19p',      175290, 175739, 'p-449'),
                            ('cn1-mat_chr20p',      75511,  75964,  'p-453'),
                            ('cn1-mat_chr20q',      154007, 154939, 'q-932'),
                            ('cn1-mat_chr22q',      117526, 119917, 'q-2391'),
                            ('cn1-pat_chr2q',       50685,  52275,  'p-1590'),
                            ('cn1-pat_chr3q',       49239,  51051,  'q-1812'),
                            ('cn1-pat_chr4q',       40859,  44538,  'q-3679'),
                            ('cn1-pat_chr5p',       42490,  43063,  'p-573'),
                            ('cn1-pat_chr6p',       153633, 154088, 'p-455'),
                            ('cn1-pat_chr7p',       119005, 119759, 'p-754'),
                            ('cn1-pat_chr9q',       128654, 129396, 'q-742'),
                            ('cn1-pat_chr10q',      19516,  20022,  'q-506'),
                            ('cn1-pat_chr10q',      49478,  50469,  'q-991'),
                            ('cn1-pat_chr11p',      32284,  32929,  'q-645'),
                            ('cn1-pat_chr12p',      190475, 191522, 'p-1047'),
                            ('cn1-pat_chr13q',      46711,  48179,  'q-1468'),
                            ('cn1-pat_chr15q',      165996, 166415, 'q-419'),
                            ('cn1-pat_chr16q',      14644,  15175,  'p-531'),
                            ('cn1-pat_chr16q',      131228, 131976, 'q-748'),
                            ('cn1-pat_chr17p',      87567,  88716,  'q-1149'),
                            ('cn1-pat_chr17q',      105971, 106594, 'q-623'),
                            ('cn1-pat_chr17q',      114107, 115098, 'q-991'),
                            ('cn1-pat_chr18p',      61272,  61818,  'p-546'),
                            ('cn1-pat_chr18p',      102871, 103357, 'p-486'),
                            ('cn1-pat_chr18q',      239591, 240084, 'q-493'),
                            ('cn1-pat_chr19p',      180757, 181211, 'p-454'),
                            ('cn1-pat_chr20p',      73844,  74284,  'p-440'),
                            ('cn1-pat_chr20q',      143894, 144729, 'q-835'),
                            ('cn1-pat_chr22q',      114951, 115662, 'q-711'),
                            ('yao-mat_chr1p',       222017, 224458, 'p-2441'),
                            ('yao-mat_chr2q',       53642,  55221,  'p-1579'),
                            ('yao-mat_chr3q',       44469,  46299,  'q-1830'),
                            ('yao-mat_chr5p',       43509,  44070,  'p-561'),
                            ('yao-mat_chr6p',       149626, 150069, 'p-443'),
                            ('yao-mat_chr7p',       115045, 115792, 'p-747'),
                            ('yao-mat_chr9q',       128374, 129110, 'q-736'),
                            ('yao-mat_chr10q',      40719,  41194,  'q-475'),
                            ('yao-mat_chr10q',      70645,  71657,  'q-1012'),
                            ('yao-mat_chr11p',      204597, 205853, 'p-1256'),
                            ('yao-mat_chr11p',      230895, 232173, 'q-1278'),
                            ('yao-mat_chr12p',      193318, 194141, 'p-823'),
                            ('yao-mat_chr13q',      46258,  47710,  'q-1452'),
                            ('yao-mat_chr16q',      28584,  29576,  'p-992'),
                            ('yao-mat_chr16q',      129944, 130579, 'p-635'),
                            ('yao-mat_chr17p',      91659,  92805,  'q-1146'),
                            ('yao-mat_chr17q',      112530, 113580, 'q-1050'),
                            ('yao-mat_chr17q',      120773, 122382, 'q-1609'),
                            ('yao-mat_chr18p',      57917,  58447,  'p-530'),
                            ('yao-mat_chr18p',      99445,  99902,  'p-457'),
                            ('yao-mat_chr18q',      236905, 237396, 'q-491'),
                            ('yao-mat_chr19p',      176027, 176480, 'p-453'),
                            ('yao-mat_chr20q',      51224,  51649,  'q-425'),
                            ('yao-mat_chr20q',      89173,  89764,  'p-591'),
                            ('yao-mat_chr22q',      116404, 118645, 'q-2241'),
                            ('yao-pat_chr2p',       208749, 210886, 'q-2137'),
                            ('yao-pat_chr2q',       49351,  50420,  'p-1069'),
                            ('yao-pat_chr3q',       50610,  52408,  'q-1798'),
                            ('yao-pat_chr5p',       43336,  43930,  'p-594'),
                            ('yao-pat_chr7p',       118387, 119154, 'p-767'),
                            ('yao-pat_chr9q',       126172, 126937, 'q-765'),
                            ('yao-pat_chr11p',      188589, 189098, 'p-509'),
                            ('yao-pat_chr11p',      215253, 215679, 'q-426'),
                            ('yao-pat_chr12p',      191278, 192307, 'p-1029'),
                            ('yao-pat_chr13q',      46137,  47399,  'q-1262'),
                            ('yao-pat_chr16q',      8210,   9274,   'p-1064'),
                            ('yao-pat_chr16q',      124197, 124920, 'q-723'),
                            ('yao-pat_chr17p',      44046,  45185,  'q-1139'),
                            ('yao-pat_chr17q',      116349, 116962, 'q-613'),
                            ('yao-pat_chr17q',      124418, 125312, 'q-894'),
                            ('yao-pat_chr18p',      58419,  58958,  'p-539'),
                            ('yao-pat_chr18p',      100616, 101087, 'p-471'),
                            ('yao-pat_chr18q',      241791, 242280, 'q-489'),
                            ('yao-pat_chr19p',      177492, 177941, 'p-449'),
                            ('yao-pat_chr20p',      70859,  71308,  'p-449'),
                            ('yao-pat_chr20q',      146446, 147284, 'q-838'),
                            ('yao-pat_chr22q',      116258, 118959, 'q-2701'),
                            ('i002c-mat_chr2p',     204704, 205245, 'q-541'),
                            ('i002c-mat_chr2q',     50889,  51930,  'p-1041'),
                            ('i002c-mat_chr3q',     47454,  49264,  'q-1810'),
                            ('i002c-mat_chr5p',     43123,  43699,  'p-576'),
                            ('i002c-mat_chr7p',     114673, 115423, 'p-750'),
                            ('i002c-mat_chr9q',     129684, 130422, 'q-738'),
                            ('i002c-mat_chr10q',    58411,  58880,  'q-469'),
                            ('i002c-mat_chr10q',    88322,  89352,  'q-1030'),
                            ('i002c-mat_chr11p',    159945, 160454, 'p-509'),
                            ('i002c-mat_chr11p',    188278, 188682, 'q-404'),
                            ('i002c-mat_chr12p',    191361, 192440, 'p-1079'),
                            ('i002c-mat_chr13q',    44796,  46209,  'q-1413'),
                            ('i002c-mat_chr16q',    13110,  14162,  'p-1052'),
                            ('i002c-mat_chr16q',    129180, 129913, 'q-733'),
                            ('i002c-mat_chr17p',    81696,  82835,  'q-1139'),
                            ('i002c-mat_chr17q',    117008, 118408, 'q-1400'),
                            ('i002c-mat_chr17q',    125661, 127433, 'q-1772'),
                            ('i002c-mat_chr18p',    58512,  59043,  'p-531'),
                            ('i002c-mat_chr18p',    100061, 100518, 'p-457'),
                            ('i002c-mat_chr18q',    240872, 241361, 'q-489'),
                            ('i002c-mat_chr19p',    177247, 177691, 'p-444'),
                            ('i002c-mat_chr20q',    146075, 146903, 'q-828'),
                            ('i002c-mat_chr22q',    120628, 122175, 'q-1547'),
                            ('i002c-pat_chr1p',     222799, 224198, 'p-1399'),
                            ('i002c-pat_chr3q',     49196,  51004,  'q-1808'),
                            ('i002c-pat_chr5p',     42169,  42753,  'p-584'),
                            ('i002c-pat_chr6p',     88125,  88641,  'p-516'),
                            ('i002c-pat_chr7p',     114114, 114865, 'p-751'),
                            ('i002c-pat_chr9q',     131843, 132584, 'q-741'),
                            ('i002c-pat_chr10q',    45685,  46209,  'q-524'),
                            ('i002c-pat_chr10q',    75650,  76731,  'q-1081'),
                            ('i002c-pat_chr11p',    159291, 159798, 'p-507'),
                            ('i002c-pat_chr12p',    190265, 191016, 'p-751'),
                            ('i002c-pat_chr13q',    46230,  47687,  'q-1457'),
                            ('i002c-pat_chr16p',    112376, 113158, 'p-782'),
                            ('i002c-pat_chr16q',    10968,  12013,  'p-1045'),
                            ('i002c-pat_chr16q',    126917, 127657, 'q-740'),
                            ('i002c-pat_chr17p',    38576,  39313,  'q-737'),
                            ('i002c-pat_chr17q',    112889, 113896, 'q-1007'),
                            ('i002c-pat_chr17q',    121095, 122690, 'q-1595'),
                            ('i002c-pat_chr18p',    58698,  59231,  'p-533'),
                            ('i002c-pat_chr18p',    100278, 100765, 'p-487'),
                            ('i002c-pat_chr18q',    237796, 238281, 'q-485'),
                            ('i002c-pat_chr19p',    176369, 176819, 'p-450'),
                            ('i002c-pat_chr22q',    117877, 118894, 'q-1017'),
                            ('HG06807-mat_chr2q',   50940,  53204,  'p-2264'),
                            ('HG06807-mat_chr3q',   47529,  49353,  'q-1824'),
                            ('HG06807-mat_chr4q',   22630,  23153,  'q-523'),
                            ('HG06807-mat_chr5p',   47361,  47965,  'p-604'),
                            ('HG06807-mat_chr5p',   89647,  90260,  'q-613'),
                            ('HG06807-mat_chr7p',   119346, 120113, 'p-767'),
                            ('HG06807-mat_chr8p',   148015, 148499, 'p-484'),
                            ('HG06807-mat_chr9q',   43880,  44652,  'q-772'),
                            ('HG06807-mat_chr10q',  13968,  15052,  'q-1084'),
                            ('HG06807-mat_chr11p',  161198, 161707, 'p-509'),
                            ('HG06807-mat_chr11p',  189487, 189904, 'q-417'),
                            ('HG06807-mat_chr12p',  193520, 194486, 'p-966'),
                            ('HG06807-mat_chr13q',  38125,  39511,  'q-1386'),
                            ('HG06807-mat_chr16p',  200175, 200944, 'p-769'),
                            ('HG06807-mat_chr16q',  71353,  72115,  'q-762'),
                            ('HG06807-mat_chr17p',  55984,  56724,  'q-740'),
                            ('HG06807-mat_chr17q',  109363, 110700, 'q-1337'),
                            ('HG06807-mat_chr17q',  117980, 119107, 'q-1127'),
                            ('HG06807-mat_chr18p',  62628,  63145,  'p-517'),
                            ('HG06807-mat_chr18p',  103978, 104434, 'p-456'),
                            ('HG06807-mat_chr18q',  238947, 239430, 'q-483'),
                            ('HG06807-mat_chr19p',  178661, 179118, 'p-457'),
                            ('HG06807-mat_chr20p',  73318,  73722,  'p-404'),
                            ('HG06807-mat_chr20q',  144493, 145393, 'q-900'),
                            ('HG06807-mat_chr22q',  104036, 106413, 'q-2377'),
                            ('HG06807-mat_chrXp',   10970,  11561,  'p-591'),
                            ('HG06807-mat_chrXp',   142805, 143472, 'p-667'),
                            ('HG06807-pat_chr2p',   205268, 205711, 'q-443'),
                            ('HG06807-pat_chr2q',   43712,  44904,  'p-1192'),
                            ('HG06807-pat_chr3q',   45797,  47597,  'q-1800'),
                            ('HG06807-pat_chr4q',   21169,  22487,  'q-1318'),
                            ('HG06807-pat_chr5p',   42797,  43400,  'p-603'),
                            ('HG06807-pat_chr6p',   163613, 164115, 'p-502'),
                            ('HG06807-pat_chr7p',   123771, 124535, 'p-764'),
                            ('HG06807-pat_chr9q',   50737,  51475,  'q-738'),
                            ('HG06807-pat_chr10q',  18140,  18665,  'q-525'),
                            ('HG06807-pat_chr10q',  48115,  49196,  'q-1081'),
                            ('HG06807-pat_chr11p',  205345, 206582, 'p-1237'),
                            ('HG06807-pat_chr11p',  230871, 232077, 'q-1206'),
                            ('HG06807-pat_chr12p',  98418,  99005,  'q-587'),
                            ('HG06807-pat_chr12p',  193023, 193973, 'p-950'),
                            ('HG06807-pat_chr13q',  38626,  40481,  'q-1855'),
                            ('HG06807-pat_chr16q',  6284,   7331,   'p-1047'),
                            ('HG06807-pat_chr16q',  123957, 124684, 'q-727'),
                            ('HG06807-pat_chr17p',  52483,  53408,  'q-925'),
                            ('HG06807-pat_chr17q',  111026, 112047, 'q-1021'),
                            ('HG06807-pat_chr17q',  119288, 121253, 'q-1965'),
                            ('HG06807-pat_chr18p',  59762,  60301,  'p-539'),
                            ('HG06807-pat_chr18p',  101352, 101847, 'p-495'),
                            ('HG06807-pat_chr18q',  232745, 233228, 'q-483'),
                            ('HG06807-pat_chr19p',  176899, 177365, 'p-466'),
                            ('HG06807-pat_chr20q',  138167, 139219, 'q-1052'),
                            ('HG06807-pat_chr22q',  115323, 115894, 'q-571'),
                            #
                            ('t2t-mouse_chr1q',       354063, 355001, 'q-938'),
                            ('t2t-mouse_chr2q',       352076, 352693, 'q-617'),
                            ('t2t-mouse_chr4q',       59298,  59703,  'q-405'),
                            ('t2t-mouse_chr9p',       350444, 350956, 'q-512'),
                            ('t2t-mouse_chr10p',      484951, 485471, 'q-520'),
                            ('t2t-mouse_chr10q',      425417, 425725, 'q-308'),
                            ('t2t-mouse_chr11q',      370707, 371083, 'q-376'),
                            ('t2t-mouse_chr11q',      403503, 404152, 'q-649'),
                            ('t2t-mouse_chr12p',      150378, 151692, 'q-1314'),
                            ('t2t-mouse_chr12p',      152574, 153720, 'q-1146'),
                            ('t2t-mouse_chr12p',      435951, 436426, 'q-475'),
                            ('t2t-mouse_chr16p',      78840,  79184,  'q-344'),
                            ('t2t-mouse_chrXq',       188515, 188911, 'q-396')]


def check_aligner_exe(fn):
    file_exists = os.path.isfile(fn)
    exe_exists = shutil.which(fn)
    if file_exists or exe_exists:
        return True
    return False


def exists_and_is_nonzero(fn):
    if os.path.isfile(fn) and os.path.getsize(fn) > 0:
        return True
    return False


def makedir(d):
    if not os.path.isdir(d):
        os.mkdir(d)


def dir_exists(d):
    return os.path.isdir(d)


def rm(fn):
    if os.path.isdir(fn):
        os.rmdir(fn)
    elif os.path.isfile(fn):
        os.remove(fn)


def mv(source, dest):
    os.rename(source, dest)


def strip_paths_from_string(s):
    if os.path.sep in s:
        s_out = s.split(os.path.sep)[-1]
        if len(s_out) > 0:
            return s_out
        else:
            return '.' + os.path.sep
    else:
        return s


def RC(s):
    return ''.join([RC_DICT[n] for n in s[::-1]])


def get_file_type(fn):
    fnl = fn.lower()
    file_type = None
    for k in FILE_SUFFIX:
        if fnl[-len(k):] == k:
            file_type = FILE_SUFFIX[k]
            break
    if file_type is None:
        print(f'Error: unknown file suffix, cannot determine input type: {fn}')
        exit(1)
    is_gzipped = (fnl[-3:] == '.gz')
    return (file_type, is_gzipped)


def shuffle_seq(s, rng=None):
    if rng is None:
        rng = np.random.default_rng(random.getrandbits(64))
    arr = np.frombuffer(s.encode('ascii'), dtype=np.uint8).copy()
    rng.shuffle(arr)
    return arr.tobytes().decode('ascii')


def shuffle_seq_old(s):
    return ''.join(random.sample(s,len(s)))


def get_downsample_inds(total_items, downsample_amount):
    inds = list(range(total_items))
    random.shuffle(inds)
    return sorted(inds[:downsample_amount], reverse=True)


def parse_cigar(cigar):
    letters = re.split(r"\d+",cigar)[1:]
    numbers = [int(n) for n in re.findall(r"\d+",cigar)]
    startPos = 0
    if letters[0] in CLIP_CHAR:
        startPos = numbers[0]
    endClip = 0
    if len(letters) > 1 and letters[-1] in CLIP_CHAR:
        endClip = numbers[-1]
    adj  = 0
    radj = 0
    for i in range(len(letters)):
        if letters[i] in REF_CHAR:
            adj += numbers[i]
        if letters[i] in READ_CHAR:
            radj += numbers[i]
    return (startPos, adj, radj, endClip)


def parse_read(splt):
    #
    cigar = splt[5]
    flag  = int(splt[1])
    ref   = splt[2]
    pos   = int(splt[3])
    rdat  = splt[9]
    rnm   = splt[0]
    mapq  = int(splt[4])
    #
    is_unmapped = False
    if ref == '*' or flag & 4:
        is_unmapped = True
    #
    orientation = 'FWD'
    if flag & 16:
        orientation = 'REV'
    #
    (read_pos_1, read_pos_2) = (None, None)
    (pos1, pos2) = (None, None)
    if is_unmapped is False:
        cig_dat    = parse_cigar(cigar)
        read_pos_1 = cig_dat[0]
        read_pos_2 = cig_dat[0] + cig_dat[2]
        read_len   = cig_dat[0] + cig_dat[2] + cig_dat[3]
        #
        if orientation == 'REV':
            [read_pos_1, read_pos_2] = [read_len - read_pos_2, read_len - read_pos_1]
        if orientation == 'FWD':
            pos1, pos2 = pos, pos + cig_dat[1]
        elif orientation == 'REV':
            pos1, pos2 = pos + cig_dat[1], pos
    else:
        ref = '*'
        pos = 0
        pos1, pos2 = 0, len(rdat)
        read_pos_1, read_pos_2 = 0, len(rdat)
        mapq = 0
    #
    return [rnm, pos, read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat]


def repeated_matches_trimming(alns, min_read_span_after_trimming=200, strategy='mapq', print_debug=False):
    #
    # trim repeated matches in the same manner as pbmm2 (only affects read_pos coords)
    #
    if print_debug:
        print('- matches trimming')
        for n in alns:
            print(n[:7], 'rdat len:', len(n[7]))
        print()
    r_coords = [[alns[n][0], alns[n][1], n] for n in range(len(alns)) if (alns[n][0] is not None and alns[n][1] is not None)]
    # don't try to trim unmapped reads
    if len(r_coords) == 0:
        return alns
    clust    = cluster_ranges(r_coords)
    any_lap  = any([len(n) > 1 for n in clust])
    # we have nothing to trim
    if any_lap is False:
        return alns
    #
    while any_lap:
        flat_clust = []
        #
        for i in range(len(clust)):
            if len(clust[i]) > 1:
                max_span = (clust[i][0][0], clust[i][0][1], alns[clust[i][0][2]][6])
                for j in range(1,len(clust[i])):
                    if strategy == 'largest' and clust[i][j][1] - clust[i][j][0] > max_span[1] - max_span[0]:
                        max_span = (clust[i][j][0], clust[i][j][1], alns[clust[i][j][2]][6])
                    elif strategy == 'mapq' and alns[clust[i][j][2]][6] > max_span[2]:
                        max_span = (clust[i][j][0], clust[i][j][1], alns[clust[i][j][2]][6])
                if print_debug:
                    print('MAX_SPAN', max_span)
                #
                max_found = False
                del_list  = []
                for j in range(len(clust[i])):
                    (x, y, q) = (clust[i][j][0], clust[i][j][1], alns[clust[i][j][2]][6])
                    # we are max (use first, if multiple)
                    if x == max_span[0] and y == max_span[1]:
                        if strategy == 'largest':
                            if max_found:
                                del_list.append(j)
                            else:
                                max_found = True
                        elif strategy == 'mapq':
                            if q == max_span[2] and max_found is False:
                                max_found = True
                            else:
                                del_list.append(j)
                    # are we completely consumed by max?
                    elif x >= max_span[0] and y <= max_span[1]:
                        del_list.append(j)
                    # do we overhang on the left?
                    elif x < max_span[0] and y >= max_span[0]:
                        clust[i][j][1] = max_span[0]-1
                    # do we overhang on the right?
                    elif x <= max_span[1] and y > max_span[1]:
                        clust[i][j][0] = max_span[1]+1
                    # did the trimming make us too short?
                    if clust[i][j][1] - clust[i][j][0] < min_read_span_after_trimming:
                        del_list.append(j)
                #
                del_list = sorted(list(set(del_list)), reverse=True)
                for di in del_list:
                    del clust[i][di]
            flat_clust.extend(clust[i])
        #
        if len(flat_clust) == 0:
            return []
        #
        flat_clust = sorted(flat_clust)
        clust      = cluster_ranges(flat_clust)
        any_lap    = any([len(n) > 1 for n in clust])
    #
    alns_out = copy.deepcopy(alns)
    del_list = []
    keep_ind = {n[2]:(n[0],n[1]) for n in flat_clust}
    for i in range(len(alns_out)):
        if i in keep_ind:
            alns_out[i][0] = keep_ind[i][0]
            alns_out[i][1] = keep_ind[i][1]
        else:
            del_list.append(i)
    del_list = sorted(list(set(del_list)), reverse=True)
    for di in del_list:
        del alns_out[di]
    #
    if print_debug:
        for n in alns_out:
            print(n[:7])
        print()
    return alns_out


def annotate_interstitial_tel(mychr, mypos, buffer=2000):
    for tr in INTERSTITIAL_TEL_REGIONS:
        if mychr == tr[0] and mypos >= tr[1]-buffer and mypos <= tr[2]+buffer:
            return True
    return False


def cluster_list(input_list, delta, which_val=None):
    #
    # cluster a sorted list
    #
    # "which_val = None" - assume input is list of values to directly sort
    # "which_val = 1"    - assume input is list of tuples, use index 1 to sort
    #
    if which_val is None:
        prev_val = input_list[0]
    else:
        prev_val = input_list[0][which_val]
    out_list    = [[input_list[0]]]
    current_ind = 0
    for n in input_list[1:]:
        if which_val is None:
            my_dist = n - prev_val
        else:
            my_dist = n[which_val] - prev_val
        if my_dist <= delta:
            out_list[current_ind].append(n)
        else:
            current_ind += 1
            out_list.append([])
            out_list[current_ind].append(n)
        if which_val is None:
            prev_val = n
        else:
            prev_val = n[which_val]
    return out_list


def cluster_ranges(input_list):
    c = [[input_list[0]]]
    for i in range(1,len(input_list)):
        found_a_home = False
        for j in range(len(c)):
            for k in range(len(c[j])):
                if input_list[i][0] <= c[j][k][1] and input_list[i][1] >= c[j][k][0]:
                    c[j].append(input_list[i])
                    found_a_home = True
                    break
            if found_a_home:
                break
        if found_a_home is False:
            c.append([input_list[i]])
    return c


def posmax(seq):
    m = seq[0]
    index = 0
    for i,x in enumerate(seq):
        if x > m:
            m = x
            index = i
    return index


def compute_n50(readlens):
    sorted_lengths = sorted(readlens, reverse=True)
    total_length = sum(sorted_lengths)
    threshold = total_length * 0.5
    cumulative_sum = 0
    for length in sorted_lengths:
        cumulative_sum += length
        if cumulative_sum >= threshold:
            return length
    return None
