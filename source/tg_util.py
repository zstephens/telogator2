import copy
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

INTERSTITIAL_TEL_REGIONS = [('t2t-chm13_chr1p',       304961, 305976, 'p-1015'),
                            ('t2t-chm13_chr2q',       203415, 207846, 'q-4431'),
                            ('t2t-chm13_chr3q',       297928, 299750, 'q-1822'),
                            ('t2t-chm13_chr4q',       270743, 274194, 'q-3451'),
                            ('t2t-chm13_chr4q',       475056, 475384, 'q-328'),
                            ('t2t-chm13_chr5p',       40627,  41173,  'p-546'),
                            ('t2t-chm13_chr7p',       113850, 114599, 'p-749'),
                            ('t2t-chm13_chr7q',       171481, 171851, 'q-370'),
                            ('t2t-chm13_chr9q',       382275, 382996, 'q-721'),
                            ('t2t-chm13_chr10q',      259721, 260692, 'q-971'),
                            ('t2t-chm13_chr10q',      474791, 475117, 'q-326'),
                            ('t2t-chm13_chr11p',      201398, 202729, 'p-1331'),
                            ('t2t-chm13_chr11p',      236825, 237626, 'q-801'),
                            ('t2t-chm13_chr13q',      239330, 240643, 'p-1313'),
                            ('t2t-chm13_chr15q',      45318,  47645,  'p-2327'),
                            ('t2t-chm13_chr16p',      12413,  12769,  'p-356'),
                            ('t2t-chm13_chr16p',      393502, 393934, 'q-432'),
                            ('t2t-chm13_chr16q',      379371, 380083, 'q-712'),
                            ('t2t-chm13_chr17q',      372489, 373350, 'q-861'),
                            ('t2t-chm13_chr18p',      55375,  56037,  'p-662'),
                            #('t2t-chm13_chr18p',      109394, 110038, 'p-644'), # I believe these two are actual telomere boundaries and that the T2T assembly for this arm is weird
                            #('t2t-chm13_chr18p',      163440, 164072, 'p-632'),
                            ('t2t-chm13_chr18p',      217355, 217862, 'p-507'),
                            ('t2t-chm13_chr18p',      258935, 259414, 'p-479'),
                            ('t2t-chm13_chr18q',      185935, 186369, 'q-434'),
                            ('t2t-chm13_chr18q',      492150, 492625, 'q-475'),
                            ('t2t-chm13_chr19p',      196741, 197186, 'p-445'),
                            ('t2t-chm13_chr20p',      70728,  71130,  'p-402'),
                            ('t2t-chm13_chr20q',      396970, 397762, 'q-792'),
                            ('t2t-chm13_chrXq',       488237, 488573, 'q-336'),
                            ('t2t-chm13_chrYq',       484759, 485099, 'q-340'),
                            ('t2t-yao-pat_chr1p',     311484, 314168, 'p-2684'),
                            ('t2t-yao-pat_chr2q',     203776, 206012, 'q-2236'),
                            ('t2t-yao-pat_chr3q',     300608, 302408, 'q-1800'),
                            ('t2t-yao-pat_chr4q',     471741, 472064, 'q-323'),
                            ('t2t-yao-pat_chr5p',     43349,  43888,  'p-539'),
                            ('t2t-yao-pat_chr6p',     103777, 104090, 'p-313'),
                            ('t2t-yao-pat_chr7p',     118391, 119148, 'p-757'),
                            ('t2t-yao-pat_chr7q',     127394, 128033, 'p-639'),
                            ('t2t-yao-pat_chr7q',     168168, 168613, 'q-445'),
                            ('t2t-yao-pat_chr9q',     376177, 376941, 'q-764'),
                            ('t2t-yao-pat_chr10q',    232679, 233643, 'q-964'),
                            ('t2t-yao-pat_chr10q',    473867, 474225, 'q-358'),
                            ('t2t-yao-pat_chr11p',    188586, 189089, 'p-503'),
                            ('t2t-yao-pat_chr13q',    239739, 240809, 'p-1070'),
                            ('t2t-yao-pat_chr15q',    43923,  48791,  'p-4868'),
                            ('t2t-yao-pat_chr16p',    17291,  17649,  'p-358'),
                            ('t2t-yao-pat_chr16p',    397634, 398003, 'q-369'),
                            ('t2t-yao-pat_chr16q',    374204, 374934, 'q-730'),
                            ('t2t-yao-pat_chr17p',    44532,  44936,  'q-404'),
                            ('t2t-yao-pat_chr17q',    374462, 375317, 'q-855'),
                            ('t2t-yao-pat_chr18p',    58423,  58930,  'p-507'),
                            ('t2t-yao-pat_chr18p',    100618, 101088, 'p-470'),
                            ('t2t-yao-pat_chr18q',    491800, 492275, 'q-475'),
                            ('t2t-yao-pat_chr19p',    177493, 177932, 'p-439'),
                            ('t2t-yao-pat_chr20p',    70855,  71301,  'p-446'),
                            ('t2t-yao-pat_chr20q',    396464, 397261, 'q-797'),
                            ('t2t-yao-pat_chrYq',     484565, 484957, 'q-392'),
                            ('t2t-yao-mat_chr1p',     222021, 224472, 'p-2451'),
                            ('t2t-yao-mat_chr2q',     197376, 197908, 'q-532'),
                            ('t2t-yao-mat_chr2q',     198511, 198855, 'q-344'),
                            ('t2t-yao-mat_chr2q',     207244, 210395, 'q-3151'),
                            ('t2t-yao-mat_chr3q',     294469, 296289, 'q-1820'),
                            ('t2t-yao-mat_chr4q',     474305, 474637, 'q-332'),
                            ('t2t-yao-mat_chr5p',     43516,  44036,  'p-520'),
                            ('t2t-yao-mat_chr6p',     149630, 150073, 'p-443'),
                            ('t2t-yao-mat_chr7p',     115046, 115782, 'p-736'),
                            ('t2t-yao-mat_chr7q',     122517, 123220, 'p-703'),
                            ('t2t-yao-mat_chr7q',     163615, 163917, 'q-302'),
                            ('t2t-yao-mat_chr8p',     148553, 148919, 'p-366'),
                            ('t2t-yao-mat_chr9q',     378375, 379108, 'q-733'),
                            ('t2t-yao-mat_chr10q',    320646, 321604, 'q-958'),
                            ('t2t-yao-mat_chr10q',    476033, 476371, 'q-338'),
                            ('t2t-yao-mat_chr11p',    204599, 205844, 'p-1245'),
                            ('t2t-yao-mat_chr11p',    231059, 231847, 'q-788'),
                            ('t2t-yao-mat_chr13q',    239822, 241274, 'p-1452'),
                            ('t2t-yao-mat_chr15q',    55067,  56944,  'p-1877'),
                            ('t2t-yao-mat_chr16p',    16101,  16449,  'p-348'),
                            ('t2t-yao-mat_chr16p',    396436, 396795, 'q-359'),
                            ('t2t-yao-mat_chr17p',    92137,  92705,  'q-568'),
                            ('t2t-yao-mat_chr17q',    370800, 372338, 'q-1538'),
                            ('t2t-yao-mat_chr18p',    57920,  58414,  'p-494'),
                            ('t2t-yao-mat_chr18p',    99447,  99902,  'p-455'),
                            ('t2t-yao-mat_chr18q',    486900, 487380, 'q-480'),
                            ('t2t-yao-mat_chr19p',    176038, 176477, 'p-439'),
                            ('t2t-yao-mat_chr20p',    112045, 112425, 'p-380'),
                            ('t2t-yao-mat_chrXq',     484032, 484358, 'q-326'),
                            ('t2t-cn1-pat_chr1p',     315714, 317904, 'p-2190'),
                            ('t2t-cn1-pat_chr2q',     194199, 194989, 'q-790'),
                            ('t2t-cn1-pat_chr2q',     204257, 207405, 'q-3148'),
                            ('t2t-cn1-pat_chr3q',     299248, 301042, 'q-1794'),
                            ('t2t-cn1-pat_chr4q',     290851, 294737, 'q-3886'),
                            ('t2t-cn1-pat_chr4q',     467535, 467852, 'q-317'),
                            ('t2t-cn1-pat_chr5p',     42498,  43034,  'p-536'),
                            ('t2t-cn1-pat_chr6p',     153647, 154086, 'p-439'),
                            ('t2t-cn1-pat_chr7p',     119005, 119749, 'p-744'),
                            ('t2t-cn1-pat_chr7q',     122142, 123058, 'p-916'),
                            ('t2t-cn1-pat_chr7q',     163199, 163637, 'q-438'),
                            ('t2t-cn1-pat_chr8p',     146801, 147169, 'p-368'),
                            ('t2t-cn1-pat_chr9q',     378674, 379384, 'q-710'),
                            ('t2t-cn1-pat_chr10q',    299489, 300414, 'q-925'),
                            ('t2t-cn1-pat_chr10q',    468273, 468607, 'q-334'),
                            ('t2t-cn1-pat_chr13q',    240569, 241348, 'p-779'),
                            ('t2t-cn1-pat_chr15q',    415991, 416429, 'q-438'),
                            ('t2t-cn1-pat_chr16p',    14941,  15296,  'p-355'),
                            ('t2t-cn1-pat_chr16p',    402409, 402768, 'q-359'),
                            ('t2t-cn1-pat_chr16q',    381237, 381975, 'q-738'),
                            ('t2t-cn1-pat_chr17p',    88046,  88598,  'q-552'),
                            ('t2t-cn1-pat_chr17q',    364132, 364962, 'q-830'),
                            ('t2t-cn1-pat_chr18p',    61278,  61792,  'p-514'),
                            ('t2t-cn1-pat_chr18p',    102875, 103357, 'p-482'),
                            ('t2t-cn1-pat_chr18q',    489595, 490074, 'q-479'),
                            ('t2t-cn1-pat_chr19p',    180769, 181210, 'p-441'),
                            ('t2t-cn1-pat_chr20p',    73845,  74284,  'p-439'),
                            ('t2t-cn1-pat_chr20q',    393908, 394707, 'q-799'),
                            ('t2t-cn1-pat_chrYq',     485584, 485979, 'q-395'),
                            ('t2t-cn1-mat_chr1p',     308036, 310145, 'p-2109'),
                            ('t2t-cn1-mat_chr2q',     191100, 191598, 'q-498'),
                            ('t2t-cn1-mat_chr2q',     199993, 200580, 'q-587'),
                            ('t2t-cn1-mat_chr2q',     200624, 202924, 'q-2300'),
                            ('t2t-cn1-mat_chr3q',     295793, 297608, 'q-1815'),
                            ('t2t-cn1-mat_chr4q',     331277, 333968, 'q-2691'),
                            ('t2t-cn1-mat_chr5p',     46642,  47187,  'p-545'),
                            ('t2t-cn1-mat_chr6p',     101542, 102058, 'p-516'),
                            ('t2t-cn1-mat_chr7p',     117923, 118660, 'p-737'),
                            ('t2t-cn1-mat_chr7q',     124160, 124845, 'p-685'),
                            ('t2t-cn1-mat_chr7q',     164990, 165436, 'q-446'),
                            ('t2t-cn1-mat_chr8p',     163205, 163576, 'p-371'),
                            ('t2t-cn1-mat_chr9q',     382912, 383679, 'q-767'),
                            ('t2t-cn1-mat_chr10q',    301597, 302540, 'q-943'),
                            ('t2t-cn1-mat_chr10q',    470276, 470617, 'q-341'),
                            ('t2t-cn1-mat_chr11p',    78385,  78809,  'p-424'),
                            ('t2t-cn1-mat_chr11p',    271533, 272772, 'p-1239'),
                            ('t2t-cn1-mat_chr11p',    299210, 299873, 'q-663'),
                            ('t2t-cn1-mat_chr13q',    241548, 242260, 'p-712'),
                            ('t2t-cn1-mat_chr16p',    15301,  15664,  'p-363'),
                            ('t2t-cn1-mat_chr16p',    396132, 397382, 'q-1250'),
                            ('t2t-cn1-mat_chr16q',    382914, 383681, 'q-767'),
                            ('t2t-cn1-mat_chr17p',    118753, 119518, 'p-765'),
                            ('t2t-cn1-mat_chr17q',    371405, 372953, 'q-1548'),
                            ('t2t-cn1-mat_chr18p',    62719,  63214,  'p-495'),
                            ('t2t-cn1-mat_chr18p',    104377, 104831, 'p-454'),
                            ('t2t-cn1-mat_chr18q',    489042, 489516, 'q-474'),
                            ('t2t-cn1-mat_chr19p',    175292, 175736, 'p-444'),
                            ('t2t-cn1-mat_chr20p',    75513,  75962,  'p-449'),
                            ('t2t-cn1-mat_chr20q',    404004, 404919, 'q-915'),
                            ('t2t-cn1-mat_chrXq',     482159, 482542, 'q-383'),
                            ('t2t-hg002-pat_chr1p',   308031, 310436, 'p-2405'),
                            ('t2t-hg002-pat_chr2q',   195719, 196494, 'q-775'),
                            ('t2t-hg002-pat_chr2q',   205135, 207222, 'q-2087'),
                            ('t2t-hg002-pat_chr3q',   302033, 303840, 'q-1807'),
                            ('t2t-hg002-pat_chr4q',   292874, 296993, 'q-4119'),
                            ('t2t-hg002-pat_chr5p',   40816,  41346,  'p-530'),
                            ('t2t-hg002-pat_chr6p',   100988, 101312, 'p-324'),
                            ('t2t-hg002-pat_chr7p',   112511, 113264, 'p-753'),
                            ('t2t-hg002-pat_chr7q',   170245, 170609, 'q-364'),
                            ('t2t-hg002-pat_chr8p',   142609, 142978, 'p-369'),
                            ('t2t-hg002-pat_chr9q',   303473, 304204, 'q-731'),
                            ('t2t-hg002-pat_chr10q',  312962, 313929, 'q-967'),
                            ('t2t-hg002-pat_chr10q',  475070, 475396, 'q-326'),
                            ('t2t-hg002-pat_chr11p',  202420, 203644, 'p-1224'),
                            ('t2t-hg002-pat_chr11p',  228898, 229596, 'q-698'),
                            ('t2t-hg002-pat_chr13q',  237780, 240045, 'p-2265'),
                            ('t2t-hg002-pat_chr15q',  30073,  33964,  'p-3891'),
                            ('t2t-hg002-pat_chr16p',  116729, 117459, 'p-730'),
                            ('t2t-hg002-pat_chr16q',  493874, 494635, 'q-761'),
                            ('t2t-hg002-pat_chr17p',  44675,  45088,  'q-413'),
                            ('t2t-hg002-pat_chr17q',  372539, 373376, 'q-837'),
                            ('t2t-hg002-pat_chr18p',  56040,  56539,  'p-499'),
                            ('t2t-hg002-pat_chr18p',  97566,  98024,  'p-458'),
                            ('t2t-hg002-pat_chr18q',  493435, 493924, 'q-489'),
                            ('t2t-hg002-pat_chr19p',  192742, 193188, 'p-446'),
                            ('t2t-hg002-pat_chr20p',  70139,  70584,  'p-445'),
                            ('t2t-hg002-pat_chr20q',  397286, 398118, 'q-832'),
                            ('t2t-hg002-pat_chrYq',   488814, 489157, 'q-343'),
                            ('t2t-hg002-mat_chr1p',   307383, 309493, 'p-2110'),
                            ('t2t-hg002-mat_chr2q',   196939, 198709, 'q-1770'),
                            ('t2t-hg002-mat_chr3q',   301623, 303441, 'q-1818'),
                            ('t2t-hg002-mat_chr4q',   244162, 244478, 'q-316'),
                            ('t2t-hg002-mat_chr4q',   475263, 475599, 'q-336'),
                            ('t2t-hg002-mat_chr5p',   41110,  41659,  'p-549'),
                            ('t2t-hg002-mat_chr6p',   100277, 100584, 'p-307'),
                            ('t2t-hg002-mat_chr7p',   112492, 113234, 'p-742'),
                            ('t2t-hg002-mat_chr7q',   168658, 169015, 'q-357'),
                            ('t2t-hg002-mat_chr8p',   142897, 143258, 'p-361'),
                            ('t2t-hg002-mat_chr9q',   302785, 303526, 'q-741'),
                            ('t2t-hg002-mat_chr10q',  264566, 265613, 'q-1047'),
                            ('t2t-hg002-mat_chr10q',  476134, 476468, 'q-334'),
                            ('t2t-hg002-mat_chr11p',  201469, 202713, 'p-1244'),
                            ('t2t-hg002-mat_chr11p',  227919, 228541, 'q-622'),
                            ('t2t-hg002-mat_chr13q',  240915, 241506, 'p-591'),
                            ('t2t-hg002-mat_chr15q',  39412,  41921,  'p-2509'),
                            ('t2t-hg002-mat_chr16p',  12608,  12975,  'p-367'),
                            ('t2t-hg002-mat_chr16p',  394358, 394931, 'q-573'),
                            ('t2t-hg002-mat_chr16q',  379701, 380418, 'q-717'),
                            ('t2t-hg002-mat_chr17p',  41859,  42515,  'q-656'),
                            ('t2t-hg002-mat_chr17q',  372981, 373811, 'q-830'),
                            ('t2t-hg002-mat_chr18p',  55736,  56238,  'p-502'),
                            ('t2t-hg002-mat_chr18p',  97322,  97782,  'p-460'),
                            ('t2t-hg002-mat_chr18q',  493094, 493561, 'q-467'),
                            ('t2t-hg002-mat_chr19p',  196013, 196444, 'p-431'),
                            ('t2t-hg002-mat_chr20p',  69721,  70153,  'p-432'),
                            ('t2t-hg002-mat_chrXq',   489026, 489347, 'q-321'),
                            ('t2t-ksa001-pat_chr1p',  311101, 314227, 'p-3126'),
                            ('t2t-ksa001-pat_chr2q',  20691,  22954,  'q-2263'),
                            ('t2t-ksa001-pat_chr3q',  293651, 295472, 'q-1821'),
                            ('t2t-ksa001-pat_chr4q',  468503, 468843, 'q-340'),
                            ('t2t-ksa001-pat_chr5p',  44523,  45039,  'p-516'),
                            ('t2t-ksa001-pat_chr6p',  102663, 103180, 'p-517'),
                            ('t2t-ksa001-pat_chr7p',  119226, 119971, 'p-745'),
                            ('t2t-ksa001-pat_chr7q',  166741, 167111, 'q-370'),
                            ('t2t-ksa001-pat_chr8p',  148463, 148826, 'p-363'),
                            ('t2t-ksa001-pat_chr9q',  285843, 286577, 'q-734'),
                            ('t2t-ksa001-pat_chr10q', 272953, 273933, 'q-980'),
                            ('t2t-ksa001-pat_chr10q', 471196, 471526, 'q-330'),
                            ('t2t-ksa001-pat_chr11p', 145332, 145831, 'p-499'),
                            ('t2t-ksa001-pat_chr13q', 236279, 237155, 'p-876'),
                            ('t2t-ksa001-pat_chr15q', 45137,  46932,  'p-1795'),
                            ('t2t-ksa001-pat_chr16p', 120509, 121280, 'p-771'),
                            ('t2t-ksa001-pat_chr16q', 379720, 380446, 'q-726'),
                            ('t2t-ksa001-pat_chr17q', 368005, 368874, 'q-869'),
                            ('t2t-ksa001-pat_chr18p', 59456,  59962,  'p-506'),
                            ('t2t-ksa001-pat_chr18p', 101048, 101530, 'p-482'),
                            ('t2t-ksa001-pat_chr18q', 185118, 185525, 'q-407'),
                            ('t2t-ksa001-pat_chr19p', 201102, 201552, 'p-450'),
                            ('t2t-ksa001-pat_chr20p', 75012,  75457,  'p-445'),
                            ('t2t-ksa001-pat_chr20q', 392807, 393710, 'q-903'),
                            ('t2t-ksa001-pat_chrXq',  482873, 483224, 'q-351'),
                            ('t2t-ksa001-mat_chr1p',  315080, 316100, 'p-1020'),
                            ('t2t-ksa001-mat_chr3q',  296903, 298711, 'q-1808'),
                            ('t2t-ksa001-mat_chr5p',  48560,  49102,  'p-542'),
                            ('t2t-ksa001-mat_chr6p',  105008, 105325, 'p-317'),
                            ('t2t-ksa001-mat_chr7p',  179959, 180715, 'p-756'),
                            ('t2t-ksa001-mat_chr7q',  121159, 123062, 'p-1903'),
                            ('t2t-ksa001-mat_chr7q',  163669, 164038, 'q-369'),
                            ('t2t-ksa001-mat_chr9q',  379743, 380461, 'q-718'),
                            ('t2t-ksa001-mat_chr10q', 261377, 262342, 'q-965'),
                            ('t2t-ksa001-mat_chr10q', 463093, 463429, 'q-336'),
                            ('t2t-ksa001-mat_chr11p', 204638, 205983, 'p-1345'),
                            ('t2t-ksa001-mat_chr11p', 237114, 237808, 'q-694'),
                            ('t2t-ksa001-mat_chr12q', 183536, 184786, 'q-1250'),
                            ('t2t-ksa001-mat_chr13q', 235599, 236669, 'p-1070'),
                            ('t2t-ksa001-mat_chr15q', 55699,  56544,  'p-845'),
                            ('t2t-ksa001-mat_chr16p', 15275,  15615,  'p-340'),
                            ('t2t-ksa001-mat_chr16p', 397685, 398163, 'q-478'),
                            ('t2t-ksa001-mat_chr16q', 321480, 322217, 'q-737'),
                            ('t2t-ksa001-mat_chr17p', 51317,  51701,  'q-384'),
                            ('t2t-ksa001-mat_chr17q', 372401, 373954, 'q-1553'),
                            ('t2t-ksa001-mat_chr18p', 58495,  58999,  'p-504'),
                            ('t2t-ksa001-mat_chr18p', 100083, 100580, 'p-497'),
                            ('t2t-ksa001-mat_chr18q', 486846, 487330, 'q-484'),
                            ('t2t-ksa001-mat_chr19p', 175997, 176445, 'p-448'),
                            ('t2t-ksa001-mat_chr20p', 70708,  71148,  'p-440'),
                            ('t2t-ksa001-mat_chrXq',  482160, 482502, 'q-342'),
                            ('t2t-i002c-pat_chr1p',   222224, 223450, 'p-1226'),
                            ('t2t-i002c-pat_chr2q',   26291,  28267,  'q-1976'),
                            ('t2t-i002c-pat_chr3q',   300227, 302018, 'q-1791'),
                            ('t2t-i002c-pat_chr4q',   217634, 218664, 'q-1030'),
                            ('t2t-i002c-pat_chr4q',   472159, 472496, 'q-337'),
                            ('t2t-i002c-pat_chr5p',   41718,  42251,  'p-533'),
                            ('t2t-i002c-pat_chr6p',   87286,  87796,  'p-510'),
                            ('t2t-i002c-pat_chr7p',   112959, 113703, 'p-744'),
                            ('t2t-i002c-pat_chr7q',   129418, 130124, 'p-706'),
                            ('t2t-i002c-pat_chr7q',   170377, 170826, 'q-449'),
                            ('t2t-i002c-pat_chr8p',   147208, 147580, 'p-372'),
                            ('t2t-i002c-pat_chr9q',   381869, 382584, 'q-715'),
                            ('t2t-i002c-pat_chr10q',  327026, 328073, 'q-1047'),
                            ('t2t-i002c-pat_chr10q',  472442, 472782, 'q-340'),
                            ('t2t-i002c-pat_chr11p',  159024, 159522, 'p-498'),
                            ('t2t-i002c-pat_chr13q',  234168, 237698, 'p-3530'),
                            ('t2t-i002c-pat_chr15q',  52178,  53146,  'p-968'),
                            ('t2t-i002c-pat_chr16p',  112017, 112792, 'p-775'),
                            ('t2t-i002c-pat_chr16q',  377899, 378615, 'q-716'),
                            ('t2t-i002c-pat_chr17q',  372221, 373748, 'q-1527'),
                            ('t2t-i002c-pat_chr18p',  58706,  59214,  'p-508'),
                            ('t2t-i002c-pat_chr18p',  100280, 100770, 'p-490'),
                            ('t2t-i002c-pat_chr18q',  488885, 489367, 'q-482'),
                            ('t2t-i002c-pat_chr19p',  172656, 173092, 'p-436'),
                            ('t2t-i002c-pat_chrYq',   485562, 485896, 'q-334'),
                            ('t2t-i002c-mat_chr1p',   307144, 308818, 'p-1674'),
                            ('t2t-i002c-mat_chr2q',   203409, 207517, 'q-4108'),
                            ('t2t-i002c-mat_chr3q',   297456, 299258, 'q-1802'),
                            ('t2t-i002c-mat_chr4q',   226642, 229386, 'q-2744'),
                            ('t2t-i002c-mat_chr4q',   473025, 473358, 'q-333'),
                            ('t2t-i002c-mat_chr5p',   43130,  43664,  'p-534'),
                            ('t2t-i002c-mat_chr6p',   103119, 103429, 'p-310'),
                            ('t2t-i002c-mat_chr7p',   113785, 114538, 'p-753'),
                            ('t2t-i002c-mat_chr7q',   126348, 127049, 'p-701'),
                            ('t2t-i002c-mat_chr7q',   167302, 167743, 'q-441'),
                            ('t2t-i002c-mat_chr8p',   144971, 145339, 'p-368'),
                            ('t2t-i002c-mat_chr9q',   379703, 380416, 'q-713'),
                            ('t2t-i002c-mat_chr10q',  340473, 341454, 'q-981'),
                            ('t2t-i002c-mat_chr10q',  472756, 473080, 'q-324'),
                            ('t2t-i002c-mat_chr11p',  159138, 159639, 'p-501'),
                            ('t2t-i002c-mat_chr13q',  238350, 239413, 'p-1063'),
                            ('t2t-i002c-mat_chr15q',  47298,  60292,  'p-12994'),
                            ('t2t-i002c-mat_chr16p',  12459,  12802,  'p-343'),
                            ('t2t-i002c-mat_chr16p',  393293, 393655, 'q-362'),
                            ('t2t-i002c-mat_chr16q',  379190, 379905, 'q-715'),
                            ('t2t-i002c-mat_chr17q',  376153, 377711, 'q-1558'),
                            ('t2t-i002c-mat_chr18p',  57579,  58077,  'p-498'),
                            ('t2t-i002c-mat_chr18p',  99113,  99562,  'p-449'),
                            ('t2t-i002c-mat_chr18q',  491239, 491718, 'q-479'),
                            ('t2t-i002c-mat_chr19p',  176942, 177377, 'p-435'),
                            ('t2t-i002c-mat_chr20p',  72676,  73065,  'p-389'),
                            ('t2t-i002c-mat_chr20q',  396034, 396863, 'q-829'),
                            ('t2t-i002c-mat_chrXq',   485129, 485479, 'q-350'),
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


def shuffle_seq(s):
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
