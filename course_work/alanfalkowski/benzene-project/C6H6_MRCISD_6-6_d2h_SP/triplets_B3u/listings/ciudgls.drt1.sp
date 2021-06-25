1
 program ciudg      
 multireference single and double excitation configuration
 interaction based on the graphical unitary group approach.


 references:  h. lischka, r. shepard, f. b. brown, and i. shavitt,
                  int. j. quantum chem. s 15, 91 (1981).
              r. shepard, r. a. bair, r. a. eades, a. f. wagner,
                  m. j. davis, l. b. harding, and t. h. dunning,
                  int j. quantum chem. s 17, 613 (1983).
              r. ahlrichs, h.-j. boehm, c. ehrhardt, p. scharf,
                  h. schiffer, h. lischka, and m. schindler,
                  j. comp. chem. 6, 200 (1985).
              r. shepard, i. shavitt, r. m. pitzer, d. c. comeau, m. pepper
                  h. lischka, p. g. szalay, r. ahlrichs, f. b. brown, and
                  j.-g. zhao, int. j. quantum chem. symp. 22, 149 (1988).

 This Version of Program CIUDG is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              CIUDG       **
     **    PROGRAM VERSION:      2009-03.    **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************


================ Computing sorting integral file structure ================

                    -----z----- -----y----- -----x----- -----w----- ---total---

                CSFs      5628      629626    11139461    12299891    24074606
      internal walks     23025       49686       15792       17871      106374
valid internal walks      5628       41898        6259        6729       60514
 getinfoarray: info=                     6 :                     1
                  8192                  6552                  8192
                  5460                     0
 icd(3)=                156520 ci%nnlev=                 17391  l2rec=
                  8192  n2max=                  5460
 lcore1,lcore2=            2621230597            2621240873
 lencor,maxblo            2621440000                 60000
========================================
 current settings:
 minbl3        6615
 minbl4        6615
 locmaxbl3   259474
 locmaxbuf   129737
 maxbl3       60000
 maxbl3       60000
 maxbl4       60000
 maxbuf       30006
========================================

 sorted 4-external integrals:   537 records of integral w-combinations 
                                530 records of integral x-combinations of length= 30006
                                        have been written (wstat,xstat=   81.01   79.81)

 sorted 3-external integrals:   241 records of integral w-combinations 
                                233 records of integral x-combinations of length= 30006
                                        have been written (wstat,xstat=   76.35   78.11)
 Orig.  diagonal integrals:  1electron:       186
                             0ext.    :       342
                             2ext.    :      6048
                             4ext.    :     28392


 Orig. off-diag. integrals:  4ext.    :  13413453
                             3ext.    :   5586948
                             2ext.    :    941415
                             1ext.    :     74598
                             0ext.    :      2427
                             2ext. SO :         0
                             1ext. SO :         0
                             0ext. SO :         0
                             1electron:      2206


 Sorted integrals            3ext.  w :   5513468 x :   5439988
                             4ext.  w :  13045051 x :  12688639


Cycle #  1 sortfile size=  31457280(      60 records of   524288) #buckets=   2
distributed memory consumption per node=         0 available core2621240873
Cycle #  2 sortfile size=  31457280(      60 records of   524288) #buckets=   2
distributed memory consumption per node=         0 available core2621240873
 minimum size of srtscr:  30408704 WP (    58 records)
 maximum size of srtscr:  31457280 WP (    60 records)
diagi   file:      8 records  of   6144 WP each=>      49152 WP total
ofdgi   file:    169 records  of   6144 WP each=>    1038336 WP total
fil3w   file:    241 records  of  30006 WP each=>    7231446 WP total
fil3x   file:    233 records  of  30006 WP each=>    6991398 WP total
fil4w   file:    537 records  of  30006 WP each=>   16113222 WP total
fil4x   file:    530 records  of  30006 WP each=>   15903180 WP total
 compressed index vector length=                 35500
 echo of the input for program ciudg:
 ------------------------------------------------------------------------
  &input
  NTYPE = 0,
  GSET = 0,
   DAVCOR =10,
  NCOREL = 30
  NROOT = 2
  IVMODE = 3
  NBKITR = 1
  NVBKMN = 2
  RTOLBK = 1e-3,1e-3,
  NITER = 40
  NVCIMN = 4
  RTOLCI = 1e-3,1e-3,
  NVCIMX = 7
  NVRFMX = 7
  NVBKMX = 7
  IDEN  = 1
  CSFPRN = 10,
 /&end
 ------------------------------------------------------------------------
lodens (list->root)=  2
invlodens (root->list)= -1  1
 bummer (warning):resetting fileloc for seriel operation0
 USING SEGMENTS OF EQUAL SIZE

****************  list of control variables  ****************
 lvlprt =    0      nroot  =    2      noldv  =   0      noldhv =   0
 nunitv =    4      nbkitr =    1      niter  =  40      davcor =  10
 csfprn =   10      ivmode =    3      istrt  =   0      vout   =   0
 iortls =    0      nvbkmx =    7      ibktv  =  -1      ibkthv =  -1
 nvcimx =    7      icitv  =   -1      icithv =  -1      frcsub =   0
 nvbkmn =    2      nvcimn =    4      maxseg = 300      nrfitr =  30
 ncorel =   30      nvrfmx =    7      nvrfmn =   4      iden   =   1
 itran  =    0      froot  =    0      rtmode =   0      ncouple=   1
 skipso =    F      dalton2=    0      molcas =   0      finalv =   0
 finalw =    0      cosmocalc=   0    with_tsklst=   0
 nsegwx =    1     1     1     1
 nseg0x =    1     1     1     1
 nseg1x =    1     1     1     1
 nseg2x =    1     1     1     1
 nseg3x =    1     1     1     1
 nseg4x =    1     1     1     1
 no0ex  =      0    no1ex  =      0    no2ex  =     0    no3ex  =     0
 no4ex  =      0    nodiag =      0
 cdg4ex =    1      c3ex1ex=    1      c2ex0ex=   1
 fileloc=    0     0     0     0     0     0     0     1     1     1
 directhd=   1      noaqccshift_zyxw=      0
 critical_crit=-1.00000    critical_delta= 0.05000

 ctol   = 0.010000    lrtshift=1.000000    smalld =0.001000


 convergence tolerances of bk and full diagonalization steps
 root #       rtolbk        rtol
 ------      --------      ------
    1        1.000E-03    1.000E-03
    2        1.000E-03    1.000E-03
 Computing density:                    .drt1.state2
 Main memory management:
 global                1 DP per process
 vdisk                 0 DP per process
 stack                 0 DP per process
 core         2621439999 DP per process
 gapointer%node_offset(*)=                     0                     0
 gapointer%node_width(*)=                     0                     0

********** Integral sort section *************


 workspace allocation information: lencor=2621439999

 echo of the input for program cisrt:
 ------------------------------------------------------------------------
  &input
  maxbl3=60000
  maxbl4=60000
  &end
 ------------------------------------------------------------------------
 
 ( 6) listing file:                    ciudgls             
 ( 5) input file:                      cisrtin   
 (17) cidrt file:                      cidrtfl             
 (11) transformed integrals file:      moints    
 (12) diagonal integral file:          diagint             
 (13) off-diagonal integral file:      ofdgint             
 (31) 4-external w integrals file:     fil4w               
 (32) 4-external x integrals file:     fil4x               
 (33) 3-external w integrals file:     fil3w               
 (34) 3-external x integrals file:     fil3x               
 (21) scratch da sorting file:         srtscr              
 (12) 2-e integral file [fsplit=2]:    moints2   

 input integral file header information:
 Hermit Integral Program : SIFS version  compute-0-12      14:54:08.910 24-Jun-21
  cidrt_title                                                                    
 MO-coefficients from mcscf.x                                                    
  with dummy occupation 1.0 for active orbitals                                  
  total ao core energy =  201.595839628                                          
 MCSCF energy =    -230.596267333                                                
 SIFS file created by program tran.      compute-0-12      14:54:25.927 24-Jun-21

 input energy(*) values:
 energy( 1)=  2.015958396276E+02, ietype=   -1,    core energy of type: Nuc.Rep.

 total core energy =   2.015958396276E+02

 nsym = 8 nmot= 186

 symmetry  =    1    2    3    4    5    6    7    8
 slabel(*) =  Ag   B3u  B2u  B1g  B1u  B2g  B3g  Au 
 nmpsy(*)  =   37   37   29   29   16   16   11   11

 info(*) =          1      8192      6552      8192      5460         0

 orbital labels, i:molab(i)=
   1:tout:001   2:tout:002   3:tout:003   4:tout:004   5:tout:005   6:tout:006   7:tout:007   8:tout:008   9:tout:009  10:tout:010
  11:tout:011  12:tout:012  13:tout:013  14:tout:014  15:tout:015  16:tout:016  17:tout:017  18:tout:018  19:tout:019  20:tout:020
  21:tout:021  22:tout:022  23:tout:023  24:tout:024  25:tout:025  26:tout:026  27:tout:027  28:tout:028  29:tout:029  30:tout:030
  31:tout:031  32:tout:032  33:tout:033  34:tout:034  35:tout:035  36:tout:036  37:tout:037  38:tout:038  39:tout:039  40:tout:040
  41:tout:041  42:tout:042  43:tout:043  44:tout:044  45:tout:045  46:tout:046  47:tout:047  48:tout:048  49:tout:049  50:tout:050
  51:tout:051  52:tout:052  53:tout:053  54:tout:054  55:tout:055  56:tout:056  57:tout:057  58:tout:058  59:tout:059  60:tout:060
  61:tout:061  62:tout:062  63:tout:063  64:tout:064  65:tout:065  66:tout:066  67:tout:067  68:tout:068  69:tout:069  70:tout:070
  71:tout:071  72:tout:072  73:tout:073  74:tout:074  75:tout:075  76:tout:076  77:tout:077  78:tout:078  79:tout:079  80:tout:080
  81:tout:081  82:tout:082  83:tout:083  84:tout:084  85:tout:085  86:tout:086  87:tout:087  88:tout:088  89:tout:089  90:tout:090
  91:tout:091  92:tout:092  93:tout:093  94:tout:094  95:tout:095  96:tout:096  97:tout:097  98:tout:098  99:tout:099 100:tout:100
 101:tout:101 102:tout:102 103:tout:103 104:tout:104 105:tout:105 106:tout:106 107:tout:107 108:tout:108 109:tout:109 110:tout:110
 111:tout:111 112:tout:112 113:tout:113 114:tout:114 115:tout:115 116:tout:116 117:tout:117 118:tout:118 119:tout:119 120:tout:120
 121:tout:121 122:tout:122 123:tout:123 124:tout:124 125:tout:125 126:tout:126 127:tout:127 128:tout:128 129:tout:129 130:tout:130
 131:tout:131 132:tout:132 133:tout:133 134:tout:134 135:tout:135 136:tout:136 137:tout:137 138:tout:138 139:tout:139 140:tout:140
 141:tout:141 142:tout:142 143:tout:143 144:tout:144 145:tout:145 146:tout:146 147:tout:147 148:tout:148 149:tout:149 150:tout:150
 151:tout:151 152:tout:152 153:tout:153 154:tout:154 155:tout:155 156:tout:156 157:tout:157 158:tout:158 159:tout:159 160:tout:160
 161:tout:161 162:tout:162 163:tout:163 164:tout:164 165:tout:165 166:tout:166 167:tout:167 168:tout:168 169:tout:169 170:tout:170
 171:tout:171 172:tout:172 173:tout:173 174:tout:174 175:tout:175 176:tout:176 177:tout:177 178:tout:178 179:tout:179 180:tout:180
 181:tout:181 182:tout:182 183:tout:183 184:tout:184 185:tout:185 186:tout:186

 input parameters:
 prnopt=  0
 ldamin=   32767 ldamax=  524288 ldainc=    4096
 maxbuf=   30006 maxbl3=   60000 maxbl4=   60000 intmxo=    6144
  Using 32 bit compression 

 drt information:
  cidrt_title                                                                    
 nmotd = 192 nfctd =   6 nfvtc =   0 nmot  = 186
 nlevel = 186 niot  =  18 lowinl= 169
 orbital-to-level map(*)
   -1  -1 169 170 171 172   1   2   3   4   5   6   7   8   9  10  11  12  13  14
   15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  -1
   -1 173 174 175  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49
   50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  -1 176
  177 178  68  69  70  71  72  73  74  75  76  77  78  79  80  81  82  83  84  85
   86  87  88  89  90  91  92  93  -1 179 180  94  95  96  97  98  99 100 101 102
  103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 181 182
  183 121 122 123 124 125 126 127 128 129 130 131 132 133 184 134 135 136 137 138
  139 140 141 142 143 144 145 146 147 148 185 149 150 151 152 153 154 155 156 157
  158 186 159 160 161 162 163 164 165 166 167 168
 compressed map(*)
  169 170 171 172   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
   17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33 173 174 175
   34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53
   54  55  56  57  58  59  60  61  62  63  64  65  66  67 176 177 178  68  69  70
   71  72  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90
   91  92  93 179 180  94  95  96  97  98  99 100 101 102 103 104 105 106 107 108
  109 110 111 112 113 114 115 116 117 118 119 120 181 182 183 121 122 123 124 125
  126 127 128 129 130 131 132 133 184 134 135 136 137 138 139 140 141 142 143 144
  145 146 147 148 185 149 150 151 152 153 154 155 156 157 158 186 159 160 161 162
  163 164 165 166 167 168
 levsym(*)
    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
    1   1   1   1   1   1   1   1   1   1   1   1   1   2   2   2   2   2   2   2
    2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2
    2   2   2   2   2   2   2   3   3   3   3   3   3   3   3   3   3   3   3   3
    3   3   3   3   3   3   3   3   3   3   3   3   3   4   4   4   4   4   4   4
    4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4   4
    5   5   5   5   5   5   5   5   5   5   5   5   5   6   6   6   6   6   6   6
    6   6   6   6   6   6   6   6   7   7   7   7   7   7   7   7   7   7   8   8
    8   8   8   8   8   8   8   8   1   1   1   1   2   2   2   3   3   3   4   4
    5   5   5   6   7   8
 repartitioning mu(*)=
   2.  2.  2.  2.  2.  2.  2.  2.  2.  2.  2.  2.  0.  0.  0.  0.  0.  0.

 new core energy added to the energy(*) list.
 from the integral file: h1_core= -2.935437811656E+02

 indxdg: diagonal integral statistics.
 total number of integrals contributing to diagonal matrix elements:     34782
 number with all external indices:     28392
 number with half external - half internal indices:      6048
 number with all internal indices:       342

 indxof: off-diagonal integral statistics.
    4-external integrals: num=   13413453 strt=          1
    3-external integrals: num=    5586948 strt=   13413454
    2-external integrals: num=     941415 strt=   19000402
    1-external integrals: num=      74598 strt=   19941817
    0-external integrals: num=       2427 strt=   20016415

 total number of off-diagonal integrals:    20018841


 indxof(2nd)  ittp=   3 numx(ittp)=      941415
 indxof(2nd)  ittp=   4 numx(ittp)=       74598
 indxof(2nd)  ittp=   5 numx(ittp)=        2427

 intermediate da file sorting parameters:
 nbuk=   2 lendar=  524288 nipbk=  349524 nipsg=2620428783
 pro2e        1   17392   34783   52174   69565   69736   69907   87298  786346 1485394
  2009682 2017874 2023334 2045173

 pro2e:  19540430 integrals read in  3579 records.

 pro2e:         0 integrals 34-ext integrals skipped.
 pro1e        1   17392   34783   52174   69565   69736   69907   87298  786346 1485394
  2009682 2017874 2023334 2045173
 pro1e: eref =   -1.325027193261321E+02
 total size of srtscr:                    60  records of                 524288 
 WP =             251658240 Bytes
 !timer: first half-sort required        cpu_time=     1.078 walltime=     3.534

 new core energy added to the energy(*) list.
 from the hamiltonian repartitioning, eref= -1.325027193261E+02
 putdg        1   17392   34783   52174   58318  582606  932131   87298  786346 1485394
  2009682 2017874 2023334 2045173

 putf:       8 buffers of length    6144 written to file 12
 diagonal integral file completed.

 putd34:   537 records of integral w-combinations and
           530 records of integral x-combinations of length= 30006 have been written.
 wstat,xstat=   81.01   79.81
 prep4e:   547 blocks of linear combinations of 4-external integrals processed.
 number of sorted 4-external integrals   25733690
 number of original 4-external integrals 13413453


 putf34: external integral file complete. nfilw=    31 nfilx=    32 nrecw=   537 nrecx=   530 lbufp= 30006

 putd34:   241 records of integral w-combinations and
           233 records of integral x-combinations of length= 30006 have been written.
 wstat,xstat=   76.35   78.11
 prep3e:   249 blocks of linear combinations of 3-external integrals processed.
 number of sorted 3-external integrals   10953456
 number of original 3-external integrals  5586948


 putf34: external integral file complete. nfilw=    33 nfilx=    34 nrecw=   241 nrecx=   233 lbufp= 30006
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd    941415         3  19000402  19000402    941415  20018841
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd     74598         4  19941817  19941817     74598  20018841
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd      2427         5  20016415  20016415      2427  20018841

 putf:     169 buffers of length    6144 written to file 13
 off-diagonal files sort completed.
 !timer: second half-sort required       cpu_time=     1.156 walltime=     7.246
 !timer: cisrt complete                  cpu_time=     2.239 walltime=    10.786
 executing brd_struct for cisrtinfo
cisrtinfo:
bufszi  6144
 diagfile 4ext:   28392 2ext:    6048 0ext:     342
 fil4w,fil4x  :13413453 fil3w,fil3x : 5586948
 ofdgint  2ext:  941415 1ext:   74598 0ext:    2427so0ext:       0so1ext:       0so2ext:       0
buffer minbl4    6615 minbl3    6615 maxbl2    6618nbas:  33  34  26  27  13  15  10  10 maxbuf 30006
 CIUDG version 5.9.7 ( 5-Oct-2004)

 workspace allocation information: lcore=2621439999

 core energy values from the integral file:
 energy( 1)=  2.015958396276E+02, ietype=   -1,    core energy of type: Nuc.Rep.
 energy( 2)= -2.935437811656E+02, ietype=    6,   fcore energy of type: H1(*)   
 energy( 3)= -1.325027193261E+02, ietype=    5,   fcore energy of type: Vref(*) 

 total core repulsion energy = -2.244506608641E+02
 nmot  =   192 niot  =    18 nfct  =     6 nfvt  =     0
 nrow  =   161 nsym  =     8 ssym  =     2 lenbuf=  1600
 nwalk,xbar:     106374    23025    49686    15792    17871
 nvalwt,nvalw:    60514     5628    41898     6259     6729
 ncsft:        24074606
 total number of valid internal walks:   60514
 nvalz,nvaly,nvalx,nvalw =     5628   41898    6259    6729

 cisrt info file parameters:
 file number  12 blocksize   6144
 mxbld  30720
 nd4ext,nd2ext,nd0ext 28392  6048   342
 n4ext,n3ext,n2ext,n1ext,n0ext,n2int,n1int,n0int 13413453  5586948   941415    74598     2427        0        0        0
 minbl4,minbl3,maxbl2  6615  6615  6618
 maxbuf 30006
 number of external orbitals per symmetry block:  33  34  26  27  13  15  10  10
 nmsym   8 number of internal orbitals  18
 executing brd_struct for drt
 executing brd_struct for orbinf
 executing brd_struct for momap
 calcthrxt: niot,maxw1=                    18                 27524
 block size     0
 pthz,pthy,pthx,pthw: 23025 49686 15792 17871 total internal walks:  106374
 maxlp3,n2lp,n1lp,n0lp 27524     0     0     0
 orbsym(*)= 1 1 1 1 2 2 2 3 3 3 4 4 5 5 5 6 7 8
setref: retained number of references =    48
 setref: total/valid number of walks=                 23025
                  5628
 nmb.of records onel     1
 nmb.of records 2-ext   154
 nmb.of records 1-ext    13
 nmb.of records 0-ext     1
 nmb.of records 2-int     0
 nmb.of records 1-int     0
 nmb.of records 0-int     0
 ---------memory usage in DP -----------------
 < n-ex core usage >
     routines:
    fourex            69518
    threx            225157
    twoex             36190
    onex               8593
    allin              6144
    diagon            32097
               =======
   maximum           225157
 
  __ static summary __ 
   reflst              5628
   hrfspc              5628
               -------
   static->            5628
 
  __ core required  __ 
   totstc              5628
   max n-ex          225157
               -------
   totnec->          230785
 
  __ core available __ 
   totspc        2621439999
   totnec -          230785
               -------
   totvec->      2621209214

 number of external paths / symmetry
 vertex x    2038    2119    2056    2055    1469    1467    1413    1411
 vertex w    2206    2119    2056    2055    1469    1467    1413    1411
segment: free space=  2621209214
 reducing frespc by                282055 to             2620927159 
  for index/conft/indsym storage .
 resegmenting ...



                   segmentation summary for type all-internal
 -------------------------------------------------------------------------------
 seg.      no. of|    no. of|  starting|  internal|  starting|  starting|
  no.    internal|        ci|       csf|     walks|      walk|       DRT|
            paths|  elements|    number|     /seg.|    number|    record|
 -------------------------------------------------------------------------------
  Z 1       23024|      5628|         0|      5628|         0|         1|
 -------------------------------------------------------------------------------
  Y 2       49686|    629626|      5628|     41898|      5628|         2|
 -------------------------------------------------------------------------------
  X 3       15747|  11139461|    635254|      6259|     47526|         5|
 -------------------------------------------------------------------------------
  W 4       17529|  12299891|  11774715|      6729|     53785|         6|
 -------------------------------------------------------------------------------
max. additional memory requirements:index=       19972DP  conft+indsym=      167592DP  drtbuffer=       94491 DP

dimension of the ci-matrix ->>>  24074606

 executing brd_struct for civct
 gentasklist: ntask=                    20
                    TASKLIST
----------------------------------------------------------------------------------------------------
TASK# BRA# KET#  T-TYPE    DESCR.   SEGMENTTYPE    SEGEL              SEGCI          VWALKS   
----------------------------------------------------------------------------------------------------
     1  3   1    24      two-ext xz   2X  3 1   15747   23024   11139461       5628    6259    5628
     2  4   1    25      two-ext wz   2X  4 1   17529   23024   12299891       5628    6729    5628
     3  4   3    26      two-ext wx*  WX  4 3   17529   15747   12299891   11139461    6729    6259
     4  4   3    27      two-ext wx+  WX  4 3   17529   15747   12299891   11139461    6729    6259
     5  2   1    11      one-ext yz   1X  2 1   49686   23024     629626       5628   41898    5628
     6  3   2    15      1ex3ex yx    3X  3 2   15747   49686   11139461     629626    6259   41898
     7  4   2    16      1ex3ex yw    3X  4 2   17529   49686   12299891     629626    6729   41898
     8  1   1     1      allint zz    OX  1 1   23024   23024       5628       5628    5628    5628
     9  2   2     5      0ex2ex yy    OX  2 2   49686   49686     629626     629626   41898   41898
    10  3   3     6      0ex2ex xx*   OX  3 3   15747   15747   11139461   11139461    6259    6259
    11  3   3    18      0ex2ex xx+   OX  3 3   15747   15747   11139461   11139461    6259    6259
    12  4   4     7      0ex2ex ww*   OX  4 4   17529   17529   12299891   12299891    6729    6729
    13  4   4    19      0ex2ex ww+   OX  4 4   17529   17529   12299891   12299891    6729    6729
    14  2   2    42      four-ext y   4X  2 2   49686   49686     629626     629626   41898   41898
    15  3   3    43      four-ext x   4X  3 3   15747   15747   11139461   11139461    6259    6259
    16  4   4    44      four-ext w   4X  4 4   17529   17529   12299891   12299891    6729    6729
    17  1   1    75      dg-024ext z  OX  1 1   23024   23024       5628       5628    5628    5628
    18  2   2    76      dg-024ext y  OX  2 2   49686   49686     629626     629626   41898   41898
    19  3   3    77      dg-024ext x  OX  3 3   15747   15747   11139461   11139461    6259    6259
    20  4   4    78      dg-024ext w  OX  4 4   17529   17529   12299891   12299891    6729    6729
----------------------------------------------------------------------------------------------------
REDTASK #   1 TIME=  19.000 N=  1 (task/type/sgbra)=(   1/24/0) (
REDTASK #   2 TIME=  18.000 N=  1 (task/type/sgbra)=(   2/25/0) (
REDTASK #   3 TIME=  17.000 N=  1 (task/type/sgbra)=(   3/26/1) (
REDTASK #   4 TIME=  16.000 N=  1 (task/type/sgbra)=(   4/27/2) (
REDTASK #   5 TIME=  15.000 N=  1 (task/type/sgbra)=(   5/11/0) (
REDTASK #   6 TIME=  14.000 N=  1 (task/type/sgbra)=(   6/15/0) (
REDTASK #   7 TIME=  13.000 N=  1 (task/type/sgbra)=(   7/16/0) (
REDTASK #   8 TIME=  12.000 N=  1 (task/type/sgbra)=(   8/ 1/0) (
REDTASK #   9 TIME=  11.000 N=  1 (task/type/sgbra)=(   9/ 5/0) (
REDTASK #  10 TIME=  10.000 N=  1 (task/type/sgbra)=(  10/ 6/1) (
REDTASK #  11 TIME=   9.000 N=  1 (task/type/sgbra)=(  11/18/2) (
REDTASK #  12 TIME=   8.000 N=  1 (task/type/sgbra)=(  12/ 7/1) (
REDTASK #  13 TIME=   7.000 N=  1 (task/type/sgbra)=(  13/19/2) (
REDTASK #  14 TIME=   6.000 N=  1 (task/type/sgbra)=(  14/42/1) (
REDTASK #  15 TIME=   5.000 N=  1 (task/type/sgbra)=(  15/43/1) (
REDTASK #  16 TIME=   4.000 N=  1 (task/type/sgbra)=(  16/44/1) (
REDTASK #  17 TIME=   3.000 N=  1 (task/type/sgbra)=(  17/75/1) (
REDTASK #  18 TIME=   2.000 N=  1 (task/type/sgbra)=(  18/76/1) (
REDTASK #  19 TIME=   1.000 N=  1 (task/type/sgbra)=(  19/77/1) (
REDTASK #  20 TIME=   0.000 N=  1 (task/type/sgbra)=(  20/78/1) (
 initializing v-file: 1:              24074606

    ---------trial vector generation----------

    trial vectors will be created by: 

    (ivmode= 3) diagonalizing h in the reference space.                     

      4 vectors will be written to unit 11 beginning with logical record   1

            4 vectors will be created
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:      170946 2x:           0 4x:           0
All internal counts: zz :     1034788 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:           0    task #     2:           0    task #     3:           0    task #     4:           0
task #     5:           0    task #     6:           0    task #     7:           0    task #     8:      877649
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:      113879    task #    18:           0    task #    19:           0    task #    20:           0
 reference space has dimension      48
 dsyevx: computed roots 1 to    8(converged:   8)

    root           eigenvalues
    ----           ------------
       1        -230.6036188789
       2        -230.5761126991
       3        -230.4791346227
       4        -230.2539226247
       5        -230.2407052199
       6        -230.2295003924
       7        -230.2149063471
       8        -230.2056666469

 strefv generated    4 initial ci vector(s).
    ---------end of vector generation---------

 ufvoutnew: ... writing  recamt=                  5628

         vector  1 from unit 11 written to unit 49 filename cirefv              
 ufvoutnew: ... writing  recamt=                  5628

         vector  2 from unit 11 written to unit 49 filename cirefv              

 ************************************************************************
 beginning the bk-type iterative procedure (nzcsf=  5628)...
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:          24074606
 number of initial trial vectors:                         4
 number of initial matrix-vector products:                0
 maximum dimension of the subspace vectors:               7
 number of roots to converge:                             2
 number of iterations:                                    1
 residual norm convergence criteria:               0.001000  0.001000

          starting bk iteration   1

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:           0 xx:           0 ww:           0
One-external counts: yz :     3649668 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       20775 wz:       17915 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:           0    task #     4:           0
task #     5:     3088196    task #     6:           0    task #     7:           0    task #     8:      877649
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:           0 xx:           0 ww:           0
One-external counts: yz :     3649668 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       20775 wz:       17915 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:           0    task #     4:           0
task #     5:     3088196    task #     6:           0    task #     7:           0    task #     8:      877649
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:           0 xx:           0 ww:           0
One-external counts: yz :     3649668 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       20775 wz:       17915 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:           0    task #     4:           0
task #     5:     3088196    task #     6:           0    task #     7:           0    task #     8:      877649
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:           0 xx:           0 ww:           0
One-external counts: yz :     3649668 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       20775 wz:       17915 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:           0    task #     4:           0
task #     5:     3088196    task #     6:           0    task #     7:           0    task #     8:      877649
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1    -6.15295801
   ht   2    -0.00000000    -6.12545184
   ht   3     0.00000000    -0.00000000    -6.02847376
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.80326176

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4
 ref    1    1.00000      -4.867218E-13  -3.177458E-13  -1.285217E-14
 ref    2   4.867218E-13    1.00000       4.468370E-13  -8.239504E-14

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4
 ref    1    1.00000        1.00000       3.006257E-25   6.954121E-27

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1     1.00000000    -0.00000000    -0.00000000    -0.00000000
 ref:   2     0.00000000     1.00000000     0.00000000    -0.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.6036188789  2.3981E-14  8.9049E-01  1.5662E+00  1.0000E-03   
 mr-sdci #  1  2   -230.5761126991 -2.1316E-14  0.0000E+00  1.5699E+00  1.0000E-03   
 mr-sdci #  1  3   -230.4791346227  3.0198E-14  0.0000E+00  1.5472E+00  1.0000E-04   
 mr-sdci #  1  4   -230.2539226247 -6.9278E-14  0.0000E+00  1.5651E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.351509
time for cinew                         3.102188
time for eigenvalue solver             0.000000
time for vector access                 0.000000

 mr-sdci  convergence not reached after  1 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.6036188789  2.3981E-14  8.9049E-01  1.5662E+00  1.0000E-03   
 mr-sdci #  1  2   -230.5761126991 -2.1316E-14  0.0000E+00  1.5699E+00  1.0000E-03   
 mr-sdci #  1  3   -230.4791346227  3.0198E-14  0.0000E+00  1.5472E+00  1.0000E-04   
 mr-sdci #  1  4   -230.2539226247 -6.9278E-14  0.0000E+00  1.5651E+00  1.0000E-04   
 
    2 of the   5 expansion vectors are transformed.
    2 of the   4 matrix-vector products are transformed.

    2 expansion eigenvectors written to unit nvfile (= 11)
    2 matrix-vector products written to unit nhvfil (= 10)

 ************************************************************************
 beginning the ci iterative diagonalization procedure... 
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:          24074606
 number of initial trial vectors:                         2
 number of initial matrix-vector products:                2
 maximum dimension of the subspace vectors:               7
 number of roots to converge:                             2
 number of iterations:                                   40
 residual norm convergence criteria:               0.001000  0.001000

          starting ci iteration   1

 Final subspace hamiltonian 

                ht   1         ht   2
   ht   1    -6.15295801
   ht   2    -0.00000000    -6.12545184

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2
 ref    1   -1.00000       5.350006E-13
 ref    2  -5.351570E-13   -1.00000    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2
 ref    1    1.00000        1.00000    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2

          reference overlap matrix  block   1

                ci   1         ci   2
 ref:   1    -1.00000000     0.00000000
 ref:   2    -0.00000000    -1.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.6036188789  8.8818E-16  8.9049E-01  1.5662E+00  1.0000E-03   
 mr-sdci #  1  2   -230.5761126991 -8.8818E-16  0.0000E+00  1.5699E+00  1.0000E-03   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.112585
time for cinew                         2.762751
time for eigenvalue solver             0.000097
time for vector access                 0.000000

          starting ci iteration   2

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3
   ht   1    -6.15295801
   ht   2    -0.00000000    -6.12545184
   ht   3     0.89049446     0.00365723    -1.63521384

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3
 ref    1  -0.906497      -3.974108E-03   0.422193    
 ref    2  -3.568788E-03   0.999992       1.750336E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3
 ref    1   0.821750        1.00000       0.178250    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3
 ref:   1    -0.90649707    -0.00397411     0.42219340
 ref:   2    -0.00356879     0.99999210     0.00175034

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  2  1   -231.2403606511  6.3674E-01  4.4533E-02  3.4698E-01  1.0000E-03   
 mr-sdci #  2  2   -230.5761131480  4.4894E-07  0.0000E+00  1.5699E+00  1.0000E-03   
 mr-sdci #  2  3   -227.6681809233 -2.8110E+00  0.0000E+00  1.4060E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.103104
time for cinew                         2.731232
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   3

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1    -6.15295801
   ht   2    -0.00000000    -6.12545184
   ht   3     0.89049446     0.00365723    -1.63521384
   ht   4     0.05488038     0.00224163     0.12923542    -0.09090701

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4
 ref    1   0.907724       3.589515E-03   0.106819       0.405727    
 ref    2   3.184964E-03  -0.999993       1.575160E-03   1.306703E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4
 ref    1   0.823972       0.999999       1.141277E-02   0.164616    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1     0.90772368     0.00358951     0.10681895     0.40572719
 ref:   2     0.00318496    -0.99999283     0.00157516     0.00130670

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  3  1   -231.2750745779  3.4714E-02  3.2328E-03  9.5657E-02  1.0000E-03   
 mr-sdci #  3  2   -230.5761158697  2.7217E-06  0.0000E+00  1.5699E+00  1.0000E-03   
 mr-sdci #  3  3   -228.6409202931  9.7274E-01  0.0000E+00  1.2585E+00  1.0000E-04   
 mr-sdci #  3  4   -227.3787591337 -2.8752E+00  0.0000E+00  1.3478E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.150269
time for cinew                         2.649414
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   4

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.15295801
   ht   2    -0.00000000    -6.12545184
   ht   3     0.89049446     0.00365723    -1.63521384
   ht   4     0.05488038     0.00224163     0.12923542    -0.09090701
   ht   5     0.01763143     0.00073457    -0.01006727    -0.00266085    -0.00526557

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.908128       3.257147E-03   0.102567       2.237537E-02  -0.405306    
 ref    2  -2.779108E-03  -0.999984       3.783764E-03   3.183685E-03  -6.759959E-04

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.824703       0.999978       1.053423E-02   5.107932E-04   0.164274    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.90812750     0.00325715     0.10256664     0.02237537    -0.40530637
 ref:   2    -0.00277911    -0.99998368     0.00378376     0.00318369    -0.00067600

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  4  1   -231.2774593078  2.3847E-03  2.9550E-04  2.8962E-02  1.0000E-03   
 mr-sdci #  4  2   -230.5761607463  4.4877E-05  0.0000E+00  1.5698E+00  1.0000E-03   
 mr-sdci #  4  3   -228.7686878204  1.2777E-01  0.0000E+00  1.1600E+00  1.0000E-04   
 mr-sdci #  4  4   -228.0001302800  6.2137E-01  0.0000E+00  1.6333E+00  1.0000E-04   
 mr-sdci #  4  5   -227.3462039730  2.8955E+00  0.0000E+00  1.3293E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.161285
time for cinew                         2.800171
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   5

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.15295801
   ht   2    -0.00000000    -6.12545184
   ht   3     0.89049446     0.00365723    -1.63521384
   ht   4     0.05488038     0.00224163     0.12923542    -0.09090701
   ht   5     0.01763143     0.00073457    -0.01006727    -0.00266085    -0.00526557
   ht   6     0.00490269    -0.00055410    -0.00201893    -0.00077734    -0.00021220    -0.00052070

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.907966       2.449959E-03  -4.739617E-02  -0.129180      -9.197980E-02   0.384964    
 ref    2  -2.314794E-03  -0.999822      -1.509480E-02   1.064803E-02   8.292918E-04   2.816172E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.824408       0.999651       2.474250E-03   1.680086E-02   8.460972E-03   0.148205    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.90796621     0.00244996    -0.04739617    -0.12918003    -0.09197980     0.38496390
 ref:   2    -0.00231479    -0.99982238    -0.01509480     0.01064803     0.00082929     0.00281617

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  5  1   -231.2777049232  2.4562E-04  2.7234E-05  8.7153E-03  1.0000E-03   
 mr-sdci #  5  2   -230.5767562225  5.9548E-04  0.0000E+00  1.5690E+00  1.0000E-03   
 mr-sdci #  5  3   -228.9622797599  1.9359E-01  0.0000E+00  1.2095E+00  1.0000E-04   
 mr-sdci #  5  4   -228.3587014624  3.5857E-01  0.0000E+00  1.1626E+00  1.0000E-04   
 mr-sdci #  5  5   -227.9396071510  5.9340E-01  0.0000E+00  1.8325E+00  1.0000E-04   
 mr-sdci #  5  6   -227.2835106785  2.8328E+00  0.0000E+00  1.3281E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.191833
time for cinew                         2.950378
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   6

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.15295801
   ht   2    -0.00000000    -6.12545184
   ht   3     0.89049446     0.00365723    -1.63521384
   ht   4     0.05488038     0.00224163     0.12923542    -0.09090701
   ht   5     0.01763143     0.00073457    -0.01006727    -0.00266085    -0.00526557
   ht   6     0.00490269    -0.00055410    -0.00201893    -0.00077734    -0.00021220    -0.00052070
   ht   7    -0.00075169    -0.00004615     0.00024164    -0.00044819     0.00006097    -0.00001184    -0.00006062

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.907892      -4.162533E-04   1.592267E-02   0.141224       2.002951E-02   8.180584E-02  -0.385281    
 ref    2  -1.693190E-03  -0.993889      -0.107672      -2.098333E-02   9.976988E-03   4.025805E-03  -5.704049E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.824270       0.987816       1.184672E-02   2.038453E-02   5.007215E-04   6.708403E-03   0.148474    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.90789180    -0.00041625     0.01592267     0.14122403     0.02002951     0.08180584    -0.38528063
 ref:   2    -0.00169319    -0.99388901    -0.10767165    -0.02098333     0.00997699     0.00402580    -0.00570405

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  6  1   -231.2777342797  2.9357E-05  6.0606E-06  4.2355E-03  1.0000E-03   
 mr-sdci #  6  2   -230.5878662248  1.1110E-02  0.0000E+00  1.5489E+00  1.0000E-03   
 mr-sdci #  6  3   -229.6853299690  7.2305E-01  0.0000E+00  1.2254E+00  1.0000E-04   
 mr-sdci #  6  4   -228.5445557348  1.8585E-01  0.0000E+00  1.2016E+00  1.0000E-04   
 mr-sdci #  6  5   -228.1641926315  2.2459E-01  0.0000E+00  1.1502E+00  1.0000E-04   
 mr-sdci #  6  6   -227.9201882755  6.3668E-01  0.0000E+00  1.9624E+00  1.0000E-04   
 mr-sdci #  6  7   -227.2661653449  2.8155E+00  0.0000E+00  1.3361E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.208435
time for cinew                         4.152527
time for eigenvalue solver             0.000183
time for vector access                 0.000000

          starting ci iteration   7

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.82707342
   ht   2    -0.00000000    -6.13720536
   ht   3     0.00000000    -0.00000000    -5.23466910
   ht   4    -0.00000000     0.00000000    -0.00000000    -4.09389487
   ht   5    -0.00013423     0.00075949    -0.00201632    -0.00110633    -0.00001309

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.907817      -1.526546E-02  -3.359803E-02   9.052600E-02   0.103841    
 ref    2   8.222652E-04  -0.886050       0.459789      -5.027518E-02   2.886915E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.824132       0.785318       0.212535       1.072255E-02   1.161637E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.90781702    -0.01526546    -0.03359803     0.09052600     0.10384094
 ref:   2     0.00082227    -0.88605006     0.45978917    -0.05027518     0.02886915

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  7  1   -231.2777426272  8.3475E-06  4.5140E-06  3.2593E-03  1.0000E-03   
 mr-sdci #  7  2   -230.6895939033  1.0173E-01  0.0000E+00  1.3256E+00  1.0000E-03   
 mr-sdci #  7  3   -230.1881285634  5.0280E-01  0.0000E+00  1.3721E+00  1.0000E-04   
 mr-sdci #  7  4   -228.9010031565  3.5645E-01  0.0000E+00  1.3452E+00  1.0000E-04   
 mr-sdci #  7  5   -228.0796992535 -8.4493E-02  0.0000E+00  1.3216E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.157043
time for cinew                         2.702087
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   8

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.82707342
   ht   2    -0.00000000    -6.13720536
   ht   3     0.00000000    -0.00000000    -5.23466910
   ht   4    -0.00000000     0.00000000    -0.00000000    -4.09389487
   ht   5    -0.00013423     0.00075949    -0.00201632    -0.00110633    -0.00001309
   ht   6     0.00043956    -0.00020637     0.00217210     0.00119450     0.00000396    -0.00001335

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.907829      -1.128579E-02  -1.455047E-02   1.021841E-02   0.149933      -1.806253E-02
 ref    2  -5.699511E-04  -0.599670       0.797042      -6.660957E-02   2.312418E-02  -4.897000E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.824154       0.359731       0.635487       4.541251E-03   2.301449E-02   3.502355E-04

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.90782914    -0.01128579    -0.01455047     0.01021841     0.14993251    -0.01806253
 ref:   2    -0.00056995    -0.59966975     0.79704160    -0.06660957     0.02312418    -0.00489700

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  8  1   -231.2777492290  6.6018E-06  2.3530E-06  2.6826E-03  1.0000E-03   
 mr-sdci #  8  2   -230.9572185228  2.6762E-01  0.0000E+00  6.8474E-01  1.0000E-03   
 mr-sdci #  8  3   -230.3694218031  1.8129E-01  0.0000E+00  1.5970E+00  1.0000E-04   
 mr-sdci #  8  4   -229.3303186444  4.2932E-01  0.0000E+00  1.1464E+00  1.0000E-04   
 mr-sdci #  8  5   -228.3583634290  2.7866E-01  0.0000E+00  1.1716E+00  1.0000E-04   
 mr-sdci #  8  6   -227.8475237676 -7.2665E-02  0.0000E+00  1.7318E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.197510
time for cinew                         2.914673
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   9

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.82707342
   ht   2    -0.00000000    -6.13720536
   ht   3     0.00000000    -0.00000000    -5.23466910
   ht   4    -0.00000000     0.00000000    -0.00000000    -4.09389487
   ht   5    -0.00013423     0.00075949    -0.00201632    -0.00110633    -0.00001309
   ht   6     0.00043956    -0.00020637     0.00217210     0.00119450     0.00000396    -0.00001335
   ht   7    -0.00030072    -0.00038405    -0.00045267     0.00014173     0.00000064    -0.00000045    -0.00000367

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.907816       1.124511E-02  -1.510030E-02   2.148658E-03  -0.125293      -7.977679E-02  -2.978109E-02
 ref    2  -1.662790E-03   0.599579       0.779122       0.158060      -6.599661E-02   3.811455E-02   5.055615E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.824132       0.359621       0.607258       2.498771E-02   2.005396E-02   7.817055E-03   3.442838E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.90781589     0.01124511    -0.01510030     0.00214866    -0.12529329    -0.07977679    -0.02978109
 ref:   2    -0.00166279     0.59957871     0.77912159     0.15806042    -0.06599661     0.03811455     0.05055615

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  9  1   -231.2777520148  2.7858E-06  8.5996E-07  1.4699E-03  1.0000E-03   
 mr-sdci #  9  2   -231.0481169216  9.0898E-02  0.0000E+00  4.5108E-01  1.0000E-03   
 mr-sdci #  9  3   -230.3717099354  2.2881E-03  0.0000E+00  1.5754E+00  1.0000E-04   
 mr-sdci #  9  4   -229.5223640797  1.9205E-01  0.0000E+00  1.0665E+00  1.0000E-04   
 mr-sdci #  9  5   -228.5156804652  1.5732E-01  0.0000E+00  1.1940E+00  1.0000E-04   
 mr-sdci #  9  6   -227.9891829332  1.4166E-01  0.0000E+00  1.2112E+00  1.0000E-04   
 mr-sdci #  9  7   -227.6992487563  4.3308E-01  0.0000E+00  1.9628E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.224731
time for cinew                         4.157715
time for eigenvalue solver             0.000183
time for vector access                 0.000000

          starting ci iteration  10

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.82709115
   ht   2     0.00000000    -6.59745606
   ht   3    -0.00000000    -0.00000000    -5.92104907
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.07170322
   ht   5     0.00018857     0.00055308     0.00003338    -0.00073300    -0.00000247

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.907792       1.400575E-02   1.056759E-02   8.910922E-03  -1.866794E-02
 ref    2   2.598710E-03   0.655744      -0.633262      -0.374560       0.143132    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.824094       0.430196       0.401133       0.140375       2.083523E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.90779230     0.01400575     0.01056759     0.00891092    -0.01866794
 ref:   2     0.00259871     0.65574353    -0.63326249    -0.37456035     0.14313189

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 10  1   -231.2777530822  1.0674E-06  6.9849E-07  1.4424E-03  1.0000E-03   
 mr-sdci # 10  2   -231.1034827110  5.5366E-02  0.0000E+00  4.2597E-01  1.0000E-03   
 mr-sdci # 10  3   -230.4838768057  1.1217E-01  0.0000E+00  1.3038E+00  1.0000E-04   
 mr-sdci # 10  4   -229.7043926143  1.8203E-01  0.0000E+00  1.4403E+00  1.0000E-04   
 mr-sdci # 10  5   -228.8984534047  3.8277E-01  0.0000E+00  1.6684E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.166870
time for cinew                         2.818359
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  11

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.82709115
   ht   2     0.00000000    -6.59745606
   ht   3    -0.00000000    -0.00000000    -5.92104907
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.07170322
   ht   5     0.00018857     0.00055308     0.00003338    -0.00073300    -0.00000247
   ht   6     0.00015378     0.00060936    -0.00027900    -0.00019522     0.00000008    -0.00000121

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.907812      -5.816875E-03  -2.346572E-02   7.032287E-03  -2.133480E-03  -4.372937E-02
 ref    2   4.534802E-03  -0.739689       0.320764      -0.468155      -0.320981      -0.152229    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.824142       0.547174       0.103440       0.219218       0.103034       2.508587E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.90781159    -0.00581687    -0.02346572     0.00703229    -0.00213348    -0.04372937
 ref:   2     0.00453480    -0.73968946     0.32076380    -0.46815481    -0.32098146    -0.15222881

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 11  1   -231.2777545131  1.4309E-06  7.0714E-07  1.4150E-03  1.0000E-03   
 mr-sdci # 11  2   -231.1743977316  7.0915E-02  0.0000E+00  3.3465E-01  1.0000E-03   
 mr-sdci # 11  3   -230.7563023016  2.7243E-01  0.0000E+00  6.9737E-01  1.0000E-04   
 mr-sdci # 11  4   -229.7191739910  1.4781E-02  0.0000E+00  1.5993E+00  1.0000E-04   
 mr-sdci # 11  5   -229.2844404124  3.8599E-01  0.0000E+00  1.5007E+00  1.0000E-04   
 mr-sdci # 11  6   -227.7192081431 -2.6997E-01  0.0000E+00  1.6116E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.203735
time for cinew                         2.993408
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  12

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.82709115
   ht   2     0.00000000    -6.59745606
   ht   3    -0.00000000    -0.00000000    -5.92104907
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.07170322
   ht   5     0.00018857     0.00055308     0.00003338    -0.00073300    -0.00000247
   ht   6     0.00015378     0.00060936    -0.00027900    -0.00019522     0.00000008    -0.00000121
   ht   7    -0.00061820    -0.00054646    -0.00012660    -0.00022263     0.00000071    -0.00000003    -0.00000163

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.907760       1.311074E-02  -1.076867E-02   6.205015E-05  -2.770818E-02   6.184132E-02  -4.784890E-02
 ref    2   6.283671E-03   0.764462       0.190953       0.404594       0.433723       2.096013E-02  -0.152637    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.824068       0.584574       3.657899E-02   0.163696       0.188883       4.263676E-03   2.558753E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.90775999     0.01311074    -0.01076867     0.00006205    -0.02770818     0.06184132    -0.04784890
 ref:   2     0.00628367     0.76446227     0.19095295     0.40459397     0.43372272     0.02096013    -0.15263685

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 12  1   -231.2777555141  1.0011E-06  4.4908E-07  1.1048E-03  1.0000E-03   
 mr-sdci # 12  2   -231.2085338089  3.4136E-02  0.0000E+00  2.0847E-01  1.0000E-03   
 mr-sdci # 12  3   -230.8447253891  8.8423E-02  0.0000E+00  4.4297E-01  1.0000E-04   
 mr-sdci # 12  4   -229.7269461806  7.7722E-03  0.0000E+00  1.4388E+00  1.0000E-04   
 mr-sdci # 12  5   -229.5002373290  2.1580E-01  0.0000E+00  1.6724E+00  1.0000E-04   
 mr-sdci # 12  6   -228.2046578292  4.8545E-01  0.0000E+00  1.5622E+00  1.0000E-04   
 mr-sdci # 12  7   -227.7172137314  1.7965E-02  0.0000E+00  1.6426E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.227539
time for cinew                         4.293823
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  13

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.82709465
   ht   2     0.00000000    -6.75787294
   ht   3    -0.00000000    -0.00000000    -6.39406453
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.27628532
   ht   5    -0.00076193     0.00014303    -0.00053118    -0.00041543    -0.00000099

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.907791       7.016195E-03   2.204808E-02   2.636583E-02  -8.794815E-02
 ref    2  -7.324539E-03   0.785830      -0.171312      -0.365219      -0.104632    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.824138       0.617577       2.983395E-02   0.134080       1.868266E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.90779072     0.00701620     0.02204808     0.02636583    -0.08794815
 ref:   2    -0.00732454     0.78582963    -0.17131209    -0.36521939    -0.10463167

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 13  1   -231.2777559704  4.5624E-07  0.0000E+00  7.3057E-04  1.0000E-03   
 mr-sdci # 13  2   -231.2207629158  1.2229E-02  7.1530E-03  1.4434E-01  1.0000E-03   
 mr-sdci # 13  3   -230.8769860405  3.2261E-02  0.0000E+00  3.7324E-01  1.0000E-04   
 mr-sdci # 13  4   -229.8679723867  1.4103E-01  0.0000E+00  1.3986E+00  1.0000E-04   
 mr-sdci # 13  5   -228.2467060505 -1.2535E+00  0.0000E+00  1.5795E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.165894
time for cinew                         5.460693
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  14

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.82709465
   ht   2     0.00000000    -6.75787294
   ht   3    -0.00000000    -0.00000000    -6.39406453
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.27628532
   ht   5    -0.00076193     0.00014303    -0.00053118    -0.00041543    -0.00000099
   ht   6     0.08377128     0.00313777     0.02204406    -0.04401593     0.00000439    -0.01392652

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.907763       1.126283E-02   1.390635E-02   1.445164E-02   2.866839E-03  -0.129879    
 ref    2  -8.153540E-03   0.804552      -0.162211       0.217302       0.254366       7.662194E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.824101       0.647431       2.650573E-02   4.742920E-02   6.471008E-02   2.273946E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.90776349     0.01126283     0.01390635     0.01445164     0.00286684    -0.12987892
 ref:   2    -0.00815354     0.80455206    -0.16221080     0.21730244     0.25436560     0.07662194

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 14  1   -231.2777562291  2.5872E-07  0.0000E+00  4.6691E-04  1.0000E-03   
 mr-sdci # 14  2   -231.2298631877  9.1003E-03  4.2143E-03  8.3528E-02  1.0000E-03   
 mr-sdci # 14  3   -230.9086506951  3.1665E-02  0.0000E+00  3.4285E-01  1.0000E-04   
 mr-sdci # 14  4   -230.3366012940  4.6863E-01  0.0000E+00  8.8584E-01  1.0000E-04   
 mr-sdci # 14  5   -228.5642476106  3.1754E-01  0.0000E+00  1.5098E+00  1.0000E-04   
 mr-sdci # 14  6   -227.8185801907 -3.8608E-01  0.0000E+00  1.7379E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.206299
time for cinew                         2.992554
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  15

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.82709465
   ht   2     0.00000000    -6.75787294
   ht   3    -0.00000000    -0.00000000    -6.39406453
   ht   4    -0.00000000    -0.00000000     0.00000000    -5.27628532
   ht   5    -0.00076193     0.00014303    -0.00053118    -0.00041543    -0.00000099
   ht   6     0.08377128     0.00313777     0.02204406    -0.04401593     0.00000439    -0.01392652
   ht   7    -0.04622798    -0.01188446    -0.05909262    -0.14020581    -0.00003839    -0.00176373    -0.02391053

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.907764      -9.611091E-03   1.545166E-02  -8.508218E-03  -1.682053E-02   5.389292E-03   0.133323    
 ref    2  -9.378575E-03  -0.861951      -0.221005      -0.107272       4.727493E-03   0.277645      -0.119702    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.824124       0.743053       4.908193E-02   1.157957E-02   3.052795E-04   7.711558E-02   3.210359E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.90776441    -0.00961109     0.01545166    -0.00850822    -0.01682053     0.00538929     0.13332266
 ref:   2    -0.00937858    -0.86195147    -0.22100493    -0.10727154     0.00472749     0.27764462    -0.11970238

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 15  1   -231.2777564819  2.5278E-07  0.0000E+00  5.1589E-04  1.0000E-03   
 mr-sdci # 15  2   -231.2393665702  9.5034E-03  5.0527E-03  1.2052E-01  1.0000E-03   
 mr-sdci # 15  3   -231.0219916811  1.1334E-01  0.0000E+00  4.4457E-01  1.0000E-04   
 mr-sdci # 15  4   -230.7741507467  4.3755E-01  0.0000E+00  5.2850E-01  1.0000E-04   
 mr-sdci # 15  5   -229.6064269361  1.0422E+00  0.0000E+00  1.0667E+00  1.0000E-04   
 mr-sdci # 15  6   -228.5300274548  7.1145E-01  0.0000E+00  1.5323E+00  1.0000E-04   
 mr-sdci # 15  7   -227.7765286800  5.9315E-02  0.0000E+00  1.7510E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.222412
time for cinew                         4.118286
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  16

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.82709562
   ht   2    -0.00000000    -6.78870571
   ht   3    -0.00000000     0.00000000    -6.57133082
   ht   4     0.00000000    -0.00000000     0.00000000    -6.32348988
   ht   5     0.00691527    -0.04901450     0.01021195    -0.00928982    -0.00868925

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.907758       9.772034E-03  -1.220352E-02  -1.280506E-02  -1.025168E-02
 ref    2   1.019524E-02   0.889250       0.152945      -3.854961E-02  -3.566236E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.824129       0.790860       2.354100E-02   1.650042E-03   1.376901E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.90775817     0.00977203    -0.01220352    -0.01280506    -0.01025168
 ref:   2     0.01019524     0.88924969     0.15294468    -0.03854961    -0.03566236

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 16  1   -231.2777566059  1.2404E-07  0.0000E+00  3.5591E-04  1.0000E-03   
 mr-sdci # 16  2   -231.2458933250  6.5268E-03  2.7714E-03  8.1807E-02  1.0000E-03   
 mr-sdci # 16  3   -231.0857339673  6.3742E-02  0.0000E+00  2.8703E-01  1.0000E-04   
 mr-sdci # 16  4   -230.8242976797  5.0147E-02  0.0000E+00  3.5632E-01  1.0000E-04   
 mr-sdci # 16  5   -228.3047287220 -1.3017E+00  0.0000E+00  1.6628E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.164917
time for cinew                         2.734985
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  17

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.82709562
   ht   2    -0.00000000    -6.78870571
   ht   3    -0.00000000     0.00000000    -6.57133082
   ht   4     0.00000000    -0.00000000     0.00000000    -6.32348988
   ht   5     0.00691527    -0.04901450     0.01021195    -0.00928982    -0.00868925
   ht   6    -0.00919076    -0.02589388     0.08721388     0.04795282     0.00171114    -0.00908115

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.907745      -1.169912E-02  -5.499588E-03   1.838699E-02   2.917953E-02  -6.656211E-03
 ref    2  -1.083109E-02  -0.902962       7.011799E-02   7.337692E-03  -6.154772E-03   1.446622E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.824119       0.815478       4.946778E-03   3.919232E-04   8.893263E-04   2.535768E-04

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.90774531    -0.01169912    -0.00549959     0.01838699     0.02917953    -0.00665621
 ref:   2    -0.01083109    -0.90296243     0.07011799     0.00733769    -0.00615477     0.01446622

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 17  1   -231.2777566840  7.8091E-08  0.0000E+00  2.8431E-04  1.0000E-03   
 mr-sdci # 17  2   -231.2506571139  4.7638E-03  1.1971E-03  5.9138E-02  1.0000E-03   
 mr-sdci # 17  3   -231.1279812120  4.2247E-02  0.0000E+00  1.6719E-01  1.0000E-04   
 mr-sdci # 17  4   -230.8464422685  2.2145E-02  0.0000E+00  3.2723E-01  1.0000E-04   
 mr-sdci # 17  5   -229.9865074005  1.6818E+00  0.0000E+00  1.0275E+00  1.0000E-04   
 mr-sdci # 17  6   -227.6542395249 -8.7579E-01  0.0000E+00  1.5434E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.215576
time for cinew                         2.892700
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  18

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.82709562
   ht   2    -0.00000000    -6.78870571
   ht   3    -0.00000000     0.00000000    -6.57133082
   ht   4     0.00000000    -0.00000000     0.00000000    -6.32348988
   ht   5     0.00691527    -0.04901450     0.01021195    -0.00928982    -0.00868925
   ht   6    -0.00919076    -0.02589388     0.08721388     0.04795282     0.00171114    -0.00908115
   ht   7    -0.02619020    -0.01988509     0.01201959    -0.00106926    -0.00054600     0.00068391    -0.00228517

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.907747      -1.093531E-02  -7.449206E-03   1.628128E-02  -4.682364E-03   5.521214E-02   3.226436E-02
 ref    2   1.102401E-02  -0.906298       5.009119E-02   2.523297E-03  -5.189156E-03  -5.755817E-02  -4.484057E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.824126       0.821497       2.564618E-03   2.714472E-04   4.885187E-05   6.361323E-03   3.051666E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.90774680    -0.01093531    -0.00744921     0.01628128    -0.00468236     0.05521214     0.03226436
 ref:   2     0.01102401    -0.90629849     0.05009119     0.00252330    -0.00518916    -0.05755817    -0.04484057

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 18  1   -231.2777567054  2.1412E-08  0.0000E+00  2.1465E-04  1.0000E-03   
 mr-sdci # 18  2   -231.2519584722  1.3014E-03  4.3613E-04  3.3903E-02  1.0000E-03   
 mr-sdci # 18  3   -231.1359314809  7.9503E-03  0.0000E+00  1.2473E-01  1.0000E-04   
 mr-sdci # 18  4   -230.8554100959  8.9678E-03  0.0000E+00  3.0306E-01  1.0000E-04   
 mr-sdci # 18  5   -230.2950383162  3.0853E-01  0.0000E+00  7.6132E-01  1.0000E-04   
 mr-sdci # 18  6   -228.2545407971  6.0030E-01  0.0000E+00  1.4614E+00  1.0000E-04   
 mr-sdci # 18  7   -227.5226071294 -2.5392E-01  0.0000E+00  1.6747E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.220215
time for cinew                         4.315186
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  19

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.82709584
   ht   2     0.00000000    -6.80129761
   ht   3    -0.00000000     0.00000000    -6.68527062
   ht   4     0.00000000    -0.00000000     0.00000000    -6.40474923
   ht   5    -0.00343173     0.00562186     0.00987318    -0.00093329    -0.00094557

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.907744      -1.158480E-02   5.241041E-03  -1.835124E-02  -4.687429E-02
 ref    2  -1.108277E-02  -0.906411      -3.829223E-02   2.164830E-03   3.747163E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.824121       0.821716       1.493764E-03   3.414544E-04   3.601322E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.90774360    -0.01158480     0.00524104    -0.01835124    -0.04687429
 ref:   2    -0.01108277    -0.90641150    -0.03829223     0.00216483     0.03747163

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 19  1   -231.2777567123  6.8564E-09  0.0000E+00  1.7216E-04  1.0000E-03   
 mr-sdci # 19  2   -231.2523587626  4.0029E-04  1.4311E-04  1.9763E-02  1.0000E-03   
 mr-sdci # 19  3   -231.1402591662  4.3277E-03  0.0000E+00  9.1053E-02  1.0000E-04   
 mr-sdci # 19  4   -230.8590028426  3.5927E-03  0.0000E+00  2.9721E-01  1.0000E-04   
 mr-sdci # 19  5   -228.8685883174 -1.4264E+00  0.0000E+00  1.7531E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.165039
time for cinew                         2.703125
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  20

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.82709584
   ht   2     0.00000000    -6.80129761
   ht   3    -0.00000000     0.00000000    -6.68527062
   ht   4     0.00000000    -0.00000000     0.00000000    -6.40474923
   ht   5    -0.00343173     0.00562186     0.00987318    -0.00093329    -0.00094557
   ht   6     0.00931383    -0.00612538    -0.00224285     0.00210276    -0.00009838    -0.00030440

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.907745      -1.114758E-02   6.738347E-03  -1.641104E-02  -1.653894E-02   8.354242E-02
 ref    2  -1.110414E-02  -0.907072      -3.518982E-02   2.005832E-03   1.468193E-02  -7.572494E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.824124       0.822904       1.283729E-03   2.733456E-04   4.890956E-04   1.271360E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.90774480    -0.01114758     0.00673835    -0.01641104    -0.01653894     0.08354242
 ref:   2    -0.01110414    -0.90707207    -0.03518982     0.00200583     0.01468193    -0.07572494

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 20  1   -231.2777567140  1.7353E-09  0.0000E+00  1.6409E-04  1.0000E-03   
 mr-sdci # 20  2   -231.2525272521  1.6849E-04  2.9953E-05  9.3977E-03  1.0000E-03   
 mr-sdci # 20  3   -231.1422657776  2.0066E-03  0.0000E+00  6.8590E-02  1.0000E-04   
 mr-sdci # 20  4   -230.8618099901  2.8071E-03  0.0000E+00  2.8753E-01  1.0000E-04   
 mr-sdci # 20  5   -230.0197066304  1.1511E+00  0.0000E+00  9.7252E-01  1.0000E-04   
 mr-sdci # 20  6   -227.6571258388 -5.9741E-01  0.0000E+00  1.7830E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.180664
time for cinew                         3.104004
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  21

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.82709584
   ht   2     0.00000000    -6.80129761
   ht   3    -0.00000000     0.00000000    -6.68527062
   ht   4     0.00000000    -0.00000000     0.00000000    -6.40474923
   ht   5    -0.00343173     0.00562186     0.00987318    -0.00093329    -0.00094557
   ht   6     0.00931383    -0.00612538    -0.00224285     0.00210276    -0.00009838    -0.00030440
   ht   7    -0.00520059     0.00303686     0.00290314    -0.00023694    -0.00004580    -0.00001267    -0.00006401

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.907743      -1.145103E-02   5.556571E-03  -1.717568E-02   1.840129E-02  -5.800315E-02  -0.109224    
 ref    2  -1.111004E-02  -0.907052      -3.390096E-02   2.341704E-03   1.984246E-03   7.966487E-03   8.228892E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.824121       0.822874       1.180151E-03   3.004875E-04   3.425449E-04   3.427831E-03   1.870145E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.90774325    -0.01145103     0.00555657    -0.01717568     0.01840129    -0.05800315    -0.10922445
 ref:   2    -0.01111004    -0.90705160    -0.03390096     0.00234170     0.00198425     0.00796649     0.08228892

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 21  1   -231.2777567149  8.3544E-10  0.0000E+00  1.5542E-04  1.0000E-03   
 mr-sdci # 21  2   -231.2525610778  3.3826E-05  1.3751E-05  5.8085E-03  1.0000E-03   
 mr-sdci # 21  3   -231.1427566695  4.9089E-04  0.0000E+00  6.2353E-02  1.0000E-04   
 mr-sdci # 21  4   -230.8619895660  1.7958E-04  0.0000E+00  2.8795E-01  1.0000E-04   
 mr-sdci # 21  5   -230.3362917546  3.1659E-01  0.0000E+00  8.5361E-01  1.0000E-04   
 mr-sdci # 21  6   -228.5782875707  9.2116E-01  0.0000E+00  1.4810E+00  1.0000E-04   
 mr-sdci # 21  7   -227.5413865331  1.8779E-02  0.0000E+00  1.8413E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.225342
time for cinew                         4.204590
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  22

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.82709585
   ht   2    -0.00000000    -6.80190021
   ht   3    -0.00000000     0.00000000    -6.69209581
   ht   4    -0.00000000    -0.00000000     0.00000000    -6.41132870
   ht   5    -0.00384284    -0.00031924    -0.00004189    -0.00004714    -0.00003384

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.907743      -1.124986E-02  -6.318047E-03   1.532499E-02   9.532213E-02
 ref    2   1.111010E-02  -0.907113       3.374995E-02  -2.142727E-03  -2.221293E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.824121       0.822981       1.178977E-03   2.394465E-04   9.579723E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.90774327    -0.01124986    -0.00631805     0.01532499     0.09532213
 ref:   2     0.01111010    -0.90711334     0.03374995    -0.00214273    -0.02221293

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 22  1   -231.2777567149  1.4388E-13  0.0000E+00  1.5547E-04  1.0000E-03   
 mr-sdci # 22  2   -231.2525720137  1.0936E-05  4.0725E-06  3.3001E-03  1.0000E-03   
 mr-sdci # 22  3   -231.1429099143  1.5324E-04  0.0000E+00  6.0772E-02  1.0000E-04   
 mr-sdci # 22  4   -230.8627563300  7.6676E-04  0.0000E+00  2.8566E-01  1.0000E-04   
 mr-sdci # 22  5   -228.8445535728 -1.4917E+00  0.0000E+00  1.8726E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.155273
time for cinew                         2.677002
time for eigenvalue solver             0.000244
time for vector access                 0.000000

          starting ci iteration  23

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.82709585
   ht   2    -0.00000000    -6.80190021
   ht   3    -0.00000000     0.00000000    -6.69209581
   ht   4    -0.00000000    -0.00000000     0.00000000    -6.41132870
   ht   5    -0.00384284    -0.00031924    -0.00004189    -0.00004714    -0.00003384
   ht   6     0.00310096     0.00037771     0.00008038     0.00013621    -0.00000210    -0.00001061

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.907742      -1.143499E-02  -5.337002E-03   1.692907E-02  -4.440061E-02  -0.168039    
 ref    2   1.111166E-02  -0.907149       3.366778E-02  -1.970281E-03  -1.580674E-02   1.426688E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.824119       0.823050       1.162003E-03   2.904754E-04   2.221267E-03   2.844074E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.90774201    -0.01143499    -0.00533700     0.01692907    -0.04440061    -0.16803926
 ref:   2     0.01111166    -0.90714889     0.03366778    -0.00197028    -0.01580674     0.01426688

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 23  1   -231.2777567151  2.5225E-10  0.0000E+00  1.5280E-04  1.0000E-03   
 mr-sdci # 23  2   -231.2525774604  5.4467E-06  1.7417E-06  2.1251E-03  1.0000E-03   
 mr-sdci # 23  3   -231.1430549353  1.4502E-04  0.0000E+00  5.8403E-02  1.0000E-04   
 mr-sdci # 23  4   -230.8630931281  3.3680E-04  0.0000E+00  2.8585E-01  1.0000E-04   
 mr-sdci # 23  5   -230.2130030148  1.3684E+00  0.0000E+00  1.0430E+00  1.0000E-04   
 mr-sdci # 23  6   -227.7445060251 -8.3378E-01  0.0000E+00  1.7898E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.178955
time for cinew                         2.838379
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  24

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.82709585
   ht   2    -0.00000000    -6.80190021
   ht   3    -0.00000000     0.00000000    -6.69209581
   ht   4    -0.00000000    -0.00000000     0.00000000    -6.41132870
   ht   5    -0.00384284    -0.00031924    -0.00004189    -0.00004714    -0.00003384
   ht   6     0.00310096     0.00037771     0.00008038     0.00013621    -0.00000210    -0.00001061
   ht   7    -0.00222047     0.00033313    -0.00094546    -0.00021356    -0.00000438     0.00000003    -0.00000451

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.907741      -1.133577E-02  -5.532097E-03   1.507250E-02   1.397776E-02   0.110878      -0.198987    
 ref    2   1.111118E-02  -0.907151       3.364259E-02  -1.869438E-03  -9.153928E-03   1.959314E-02   1.021576E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.824118       0.823051       1.162428E-03   2.306750E-04   2.791721E-04   1.267790E-02   3.970008E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.90774129    -0.01133577    -0.00553210     0.01507250     0.01397776     0.11087837    -0.19898672
 ref:   2     0.01111118    -0.90715067     0.03364259    -0.00186944    -0.00915393     0.01959314     0.01021576

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 24  1   -231.2777567152  9.6357E-11  0.0000E+00  1.5101E-04  1.0000E-03   
 mr-sdci # 24  2   -231.2525792594  1.7990E-06  6.7722E-07  1.3138E-03  1.0000E-03   
 mr-sdci # 24  3   -231.1430619170  6.9818E-06  0.0000E+00  5.8349E-02  1.0000E-04   
 mr-sdci # 24  4   -230.8636911992  5.9807E-04  0.0000E+00  2.8343E-01  1.0000E-04   
 mr-sdci # 24  5   -230.5710077806  3.5800E-01  0.0000E+00  7.1848E-01  1.0000E-04   
 mr-sdci # 24  6   -228.1827712136  4.3827E-01  0.0000E+00  1.3192E+00  1.0000E-04   
 mr-sdci # 24  7   -227.7194857510  1.7810E-01  0.0000E+00  1.8407E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.191895
time for cinew                         4.094727
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  25

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1523481 2x:      354974 4x:       54886
All internal counts: zz :     1034788 yy:     8467233 xx:      769004 ww:      818396
One-external counts: yz :     3649668 yx:     5083171 yw:     5121536
Two-external counts: yy :     3234074 ww:      389580 xx:      417504 xz:       20775 wz:       17915 wx:      695182
Three-ext.   counts: yx :      972741 yw:     1009554

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:      109673
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.82709585
   ht   2    -0.00000000    -6.80191840
   ht   3    -0.00000000     0.00000000    -6.69240105
   ht   4     0.00000000    -0.00000000     0.00000000    -6.41303034
   ht   5    -0.00061241    -0.00007751    -0.00061805    -0.00008985    -0.00000172

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.907741      -1.139169E-02   5.225833E-03  -1.669692E-02  -9.979526E-02
 ref    2  -1.111154E-02  -0.907158      -3.365063E-02   1.721771E-03  -1.193882E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.824117       0.823066       1.159674E-03   2.817515E-04   1.010163E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.90774063    -0.01139169     0.00522583    -0.01669692    -0.09979526
 ref:   2    -0.01111154    -0.90715837    -0.03365063     0.00172177    -0.01193882

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 25  1   -231.2777567153  9.1140E-11  0.0000E+00  1.4989E-04  1.0000E-03   
 mr-sdci # 25  2   -231.2525798997  6.4025E-07  2.6041E-07  7.9436E-04  1.0000E-03   
 mr-sdci # 25  3   -231.1430800376  1.8121E-05  0.0000E+00  5.8464E-02  1.0000E-04   
 mr-sdci # 25  4   -230.8641436054  4.5241E-04  0.0000E+00  2.8351E-01  1.0000E-04   
 mr-sdci # 25  5   -229.1487203623 -1.4223E+00  0.0000E+00  1.6164E+00  1.0000E-04   
 
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.156250
time for cinew                         2.855957
time for eigenvalue solver             0.000000
time for vector access                 0.000000

 mr-sdci  convergence criteria satisfied after 25 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 25  1   -231.2777567153  9.1140E-11  0.0000E+00  1.4989E-04  1.0000E-03   
 mr-sdci # 25  2   -231.2525798997  6.4025E-07  2.6041E-07  7.9436E-04  1.0000E-03   
 mr-sdci # 25  3   -231.1430800376  1.8121E-05  0.0000E+00  5.8464E-02  1.0000E-04   
 mr-sdci # 25  4   -230.8641436054  4.5241E-04  0.0000E+00  2.8351E-01  1.0000E-04   
 mr-sdci # 25  5   -229.1487203623 -1.4223E+00  0.0000E+00  1.6164E+00  1.0000E-04   

####################CIUDGINFO####################

   ci vector at position   1 energy= -231.277756715303
   ci vector at position   2 energy= -231.252579899652

################END OF CIUDGINFO################

 
    2 of the   6 expansion vectors are transformed.
    2 of the   5 matrix-vector products are transformed.

    2 expansion eigenvectors written to unit nvfile (= 11)
    2 matrix-vector products written to unit nhvfil (= 10)


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 1) =      -231.2777567153

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
                                               18

                                          orbital     3    4    5    6   42   43   44   80   81   82  110  111  139  140  141  155  171
                                              182

                                         symmetry   ag   ag   ag   ag   b3u  b3u  b3u  b2u  b2u  b2u  b1g  b1g  b1u  b1u  b1u  b2g  b3g
                                              au 

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  3  1       2  0.087088                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +  
                                             +  
 z*  3  1       4  0.042994                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -        +  
                                             +  
 z*  3  1       7  0.588016                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +    +- 
                                                
 z*  3  1       8 -0.089483                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +       
                                             +- 
 z*  3  1       9  0.013003                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-        +  
                                             +  
 z*  3  1      10  0.210593                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +    +    +- 
                                                
 z*  3  1      11 -0.037202                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +    +       
                                             +- 
 z*  3  1      12 -0.634712                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +  
                                             +  
 z*  3  1      17 -0.049742                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-        +    +- 
                                                
 z*  3  1      22 -0.034971                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -        +-   +  
                                             +  
 z*  3  1      25 -0.016157                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +     -   +    +- 
                                                
 z*  3  1      29 -0.015087                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +         +-    - 
                                             +  
 z*  3  1      30  0.052606                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +         +-   +  
                                              - 
 z*  3  1      33 -0.012995                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -   +-   +  
                                             +  
 z*  3  1      34 -0.024348                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +    +-    - 
                                             +  
 z*  3  1      35  0.028139                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +    +-   +  
                                              - 
 z*  3  1      36 -0.050680                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +    +- 
                                             +- 
 z*  3  1      40  0.024055                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-        +-   +  
                                             +  
 z*  3  1      43  0.012861                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +     -   +-   +  
                                             +  
 z*  3  1      46 -0.026158                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +         +    +- 
                                             +- 
 y   3  1    5826  0.019986              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +-      
                                                
 y   3  1    5881 -0.021034              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +- 
                                                
 y   3  1    6134  0.017224              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-              -   +  
                                             +  
 y   3  1    6162  0.027037              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +     - 
                                             +  
 y   3  1    6164 -0.011061              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +     - 
                                             +  
 y   3  1    6488 -0.010950              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -   +          -   +  
                                             +  
 y   3  1    6892  0.017704              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-             +- 
                                                
 y   3  1    7145  0.036854              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -         -   +  
                                             +  
 y   3  1    7147 -0.016654              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -         -   +  
                                             +  
 y   3  1    7150 -0.010875              8( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -         -   +  
                                             +  
 y   3  1    7489  0.024116              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +         +     - 
                                              - 
 y   3  1    7491 -0.010603              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +         +     - 
                                              - 
 y   3  1    7655  0.014369              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -    -   +  
                                             +  
 y   3  1    7856  0.017985              2( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-   +- 
                                                
 y   3  1    7857 -0.059429              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-   +- 
                                                
 y   3  1    7859  0.023616              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-   +- 
                                                
 y   3  1    7862  0.013648              8( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-   +- 
                                                
 y   3  1    7898  0.019223              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-      
                                             +- 
 y   3  1  591437 -0.010382             13( b3u)   +    +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +- 
                                                

 ci coefficient statistics:
           rq > 0.1                3
      0.1> rq > 0.01              36
     0.01> rq > 0.001          27992
    0.001> rq > 0.0001        753659
   0.0001> rq > 0.00001      3934470
  0.00001> rq > 0.000001     8630901
 0.000001> rq               10727545
           all              24074606
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:      170946 2x:           0 4x:           0
All internal counts: zz :     1034788 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:       15706
  iref  icsf         v(icsf)             hv(icsf)
     1     1      0.003727462670     -0.860202851379
     2     2      0.087087589975    -20.088273299157
     3     3      0.004598270041     -1.060966123166
     4     4      0.042993836305     -9.917170097633
     5     5      0.002753831197     -0.635432152057
     6     6     -0.004736252422      1.092964637586
     7     7      0.588016270178   -135.596697564186
     8     8     -0.089483209584     20.640727848943
     9     9      0.013002587466     -2.999197552545
    10    10      0.210592553160    -48.562082523918
    11    11     -0.037202457993      8.581769513322
    12    12     -0.634711797804    146.364939161912
    13    13     -0.000208007582      0.047965245290
    14    14     -0.000964494408      0.223299911029
    15    15     -0.003349200152      0.773211880971
    16    16      0.001967529308     -0.454284348762
    17    17     -0.049742110322     11.475691004323
    18    18      0.001430669673     -0.330499176467
    19    19     -0.000333838771      0.077364499803
    20    20     -0.003878981214      0.895722857237
    21    21     -0.001632003032      0.376332376436
    22    22     -0.034970611525      8.066497784862
    23    23     -0.001272202805      0.293738836381
    24    24      0.000545336745     -0.125976446267
    25    25     -0.016156890237      3.724607545396
    26    26     -0.001841051907      0.425489710831
    27    27      0.009968048483     -2.300214513547
    28    28     -0.002089082909      0.482242219074
    29    29     -0.015087129346      3.479901810791
    30    30      0.052606062416    -12.135300132501
    31    31     -0.000467923846      0.106749001569
    32    32     -0.001862255772      0.430009154830
    33    33     -0.012994852238      3.000560472815
    34    34     -0.024348277121      5.618042354739
    35    35      0.028138670824     -6.492107843624
    36    36     -0.050679698994     11.691098313781
    37    37     -0.000757139432      0.174830801137
    38    38     -0.001759658841      0.406081734347
    39    39      0.000106578106     -0.024582254101
    40    40      0.024055279583     -5.548154376012
    41    41      0.003129295346     -0.722013047338
    42    42     -0.000801671259      0.185127634098
    43    43      0.012861469984     -2.966536792002
    44    44     -0.000910431155      0.210197522651
    45    45      0.001776870488     -0.410186058370
    46    46     -0.026157930073      6.033272203605
    47    47      0.007517451809     -1.734323616767
    48    48     -0.007740264173      1.784934300158

 number of reference csfs (nref) is    48.  root number (iroot) is  1.
 c0**2 =   0.82481622  c**2 (all zwalks) =   0.82536333

 pople ci energy extrapolation is computed with 30 correlated electrons.

 eref      =   -230.603174579159   "relaxed" cnot**2         =   0.824816223649
 eci       =   -231.277756715303   deltae = eci - eref       =  -0.674582136144
 eci+dv1   =   -231.395932561372   dv1 = (1-cnot**2)*deltae  =  -0.118175846069
 eci+dv2   =   -231.421032080991   dv2 = dv1 / cnot**2       =  -0.143275365688
 eci+dv3   =   -231.459668574171   dv3 = dv1 / (2*cnot**2-1) =  -0.181911858868
 eci+pople =   -231.440963093020   ( 30e- scaled deltae )    =  -0.837788513861


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 2) =      -231.2525798997

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
                                               18

                                          orbital     3    4    5    6   42   43   44   80   81   82  110  111  139  140  141  155  171
                                              182

                                         symmetry   ag   ag   ag   ag   b3u  b3u  b3u  b2u  b2u  b2u  b1g  b1g  b1u  b1u  b1u  b2g  b3g
                                              au 

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  3  1       2 -0.039329                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +  
                                             +  
 z*  3  1       5  0.012288                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +          - 
                                             +  
 z*  3  1       7  0.608492                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +    +- 
                                                
 z*  3  1       8 -0.038345                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +       
                                             +- 
 z*  3  1      10  0.169176                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +    +    +- 
                                                
 z*  3  1      12  0.619591                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +  
                                             +  
 z*  3  1      17 -0.074887                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-        +    +- 
                                                
 z*  3  1      18  0.012476                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-        +       
                                             +- 
 z*  3  1      20 -0.011001                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +    +    +- 
                                                
 z*  3  1      25 -0.025842                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +     -   +    +- 
                                                
 z*  3  1      29  0.155079                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +         +-    - 
                                             +  
 z*  3  1      34  0.070493                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +    +-    - 
                                             +  
 z*  3  1      35 -0.015328                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +    +-   +  
                                              - 
 z*  3  1      40 -0.030239                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-        +-   +  
                                             +  
 z*  3  1      43 -0.019547                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +     -   +-   +  
                                             +  
 z*  3  1      46 -0.035862                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +         +    +- 
                                             +- 
 z*  3  1      47 -0.010271                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +-   +  
                                             +  
 z*  3  1      48 -0.011208                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +    +    +- 
                                             +- 
 y   3  1    5841  0.013640              3( b3g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -   +  
                                                
 y   3  1    5881 -0.029816              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +- 
                                                
 y   3  1    5883  0.012102              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +- 
                                                
 y   3  1    6162 -0.027097              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +     - 
                                             +  
 y   3  1    6164  0.012016              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +     - 
                                             +  
 y   3  1    6176 -0.016141              2( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +    +  
                                              - 
 y   3  1    6177  0.052598              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +    +  
                                              - 
 y   3  1    6179 -0.020453              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +    +  
                                              - 
 y   3  1    6182 -0.011493              8( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +    +  
                                              - 
 y   3  1    6205 -0.010634              3( b3g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +- 
                                             +  
 y   3  1    6488  0.010096              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -   +          -   +  
                                             +  
 y   3  1    6892  0.011425              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-             +- 
                                                
 y   3  1    7145 -0.017717              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -         -   +  
                                             +  
 y   3  1    7173 -0.015954              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -        +     - 
                                             +  
 y   3  1    7489  0.014240              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +         +     - 
                                              - 
 y   3  1    9880  0.010915              9( b2u)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +          -   +  
                                             +  
 y   3  1   88487 -0.010317              8( b1g)   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-   +          -   +  
                                             +  
 w   3  111778816 -0.011244    3( b2g)   3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +       
                                                
 w   3  111795561 -0.011784    3( b2g)   3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +  
                                             +  

 ci coefficient statistics:
           rq > 0.1                4
      0.1> rq > 0.01              33
     0.01> rq > 0.001          28359
    0.001> rq > 0.0001        757788
   0.0001> rq > 0.00001      3625900
  0.00001> rq > 0.000001     8570004
 0.000001> rq               11092518
           all              24074606
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:      170946 2x:           0 4x:           0
All internal counts: zz :     1034788 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       34785    task #     2:       27571    task #     3:      346917    task #     4:      173458
task #     5:     3088196    task #     6:     4047341    task #     7:     4057693    task #     8:      877649
task #     9:     7773377    task #    10:      216687    task #    11:      216687    task #    12:      181831
task #    13:      181831    task #    14:       90915    task #    15:       90917    task #    16:       90917
task #    17:      113879    task #    18:      767180    task #    19:      100479    task #    20:       15706
  iref  icsf         v(icsf)             hv(icsf)
     1     1      0.006410024620     -1.478825081507
     2     2     -0.039329287934      9.070687267705
     3     3      0.006259168233     -1.443928273639
     4     4     -0.004715633866      1.085680399710
     5     5      0.012288435206     -2.835769785839
     6     6     -0.008381756353      1.934449339971
     7     7      0.608491652128   -140.302079520949
     8     8     -0.038345347320      8.843252728938
     9     9     -0.000007191070      0.000634134524
    10    10      0.169175907941    -39.004311860499
    11    11     -0.001933389298      0.444583342055
    12    12      0.619590513344   -142.860483943694
    13    13     -0.000310559070      0.071610553028
    14    14      0.000836074583     -0.193219664857
    15    15      0.001436251136     -0.331445214474
    16    16     -0.000779620735      0.179961974344
    17    17     -0.074887276297     17.271336158959
    18    18      0.012476148571     -2.877780008040
    19    19      0.001083976237     -0.250459860361
    20    20     -0.011000986139      2.537775254230
    21    21      0.003155616142     -0.728078559330
    22    22     -0.001263546192      0.292853727904
    23    23     -0.000813631561      0.187844212938
    24    24      0.000041965757     -0.009730295592
    25    25     -0.025841912321      5.957079304557
    26    26      0.004310681590     -0.993807685541
    27    27      0.006070299001     -1.400569394744
    28    28     -0.000913718803      0.210772054448
    29    29      0.155079173175    -35.766248978399
    30    30     -0.009045775899      2.083693529189
    31    31     -0.002170839603      0.499533109740
    32    32      0.000796367458     -0.183423386330
    33    33      0.002854875212     -0.660448335917
    34    34      0.070492610608    -16.259040644059
    35    35     -0.015327562081      3.534971714592
    36    36     -0.008969296067      2.065736493884
    37    37      0.000475407230     -0.109743034167
    38    38     -0.000266344105      0.061363242053
    39    39     -0.000608639808      0.140555086833
    40    40     -0.030239378010      6.975099725505
    41    41      0.003828409424     -0.883253654290
    42    42     -0.000792185822      0.182905998809
    43    43     -0.019547029922      4.509546951794
    44    44     -0.000212330870      0.049256565587
    45    45      0.001503045060     -0.347175365709
    46    46     -0.035862269790      8.272085125765
    47    47     -0.010271115400      2.369749506282
    48    48     -0.011207518566      2.585348198025

 number of reference csfs (nref) is    48.  root number (iroot) is  2.
 c0**2 =   0.82500678  c**2 (all zwalks) =   0.82566651

 pople ci energy extrapolation is computed with 30 correlated electrons.

 eref      =   -230.575490967075   "relaxed" cnot**2         =   0.825006783506
 eci       =   -231.252579899652   deltae = eci - eref       =  -0.677088932576
 eci+dv1   =   -231.371065869815   dv1 = (1-cnot**2)*deltae  =  -0.118485970164
 eci+dv2   =   -231.396198076535   dv2 = dv1 / cnot**2       =  -0.143618176884
 eci+dv3   =   -231.434862202940   dv3 = dv1 / (2*cnot**2-1) =  -0.182282303289
 eci+pople =   -231.416130327558   ( 30e- scaled deltae )    =  -0.840639360482
maximum overlap with reference    1(overlap= 0.90774)
weight of reference states=  0.8241

 information on vector: 1 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    2(overlap= 0.90716)
weight of reference states=  0.8231

 information on vector: 2 from unit 11 written to unit 48 filename civout              
 passed aftci ... 
 readint2: molcas,dalton2=                     0                     0
 files%faoints=aoints              
lodens (list->root)=  1  2
sifcfg setup: record length 4096 DP
# d1 elements per record  3272
# d2 elements per record  2730
  The MR-CISD density will be calculated.
 item #                     1 suffix=:.drt1.state1:
 read_civout: repnuc=  -224.450660864067     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 91% root-following 0
 MR-CISD energy:  -231.27775672    -6.82709585
 residuum:     0.00014989
 deltae:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.90774063    -0.01139169     0.00522583    -0.01669692    -0.09979526     0.11087837    -0.19898672     0.00000000
 ref:   2    -0.01111154    -0.90715837    -0.03365063     0.00172177    -0.01193882     0.01959314     0.01021576     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.90774063    -0.01139169     0.00522583    -0.01669692    -0.09979526     0.00000000     0.00000000     0.00000000
 ref:   2    -0.01111154    -0.90715837    -0.03365063     0.00172177    -0.01193882     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 computing MRCISD density
 densi: densityinfo%a4den=   1.00000000000000     
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    1
--------------------------------------------------------------------------------
================================================================================
   DYZ=   80766  DYX=  101959  DYW=  106443
   D0Z=   19955  D0Y=  162182  D0X=   19509  D0W=   21336
  DDZI=   40986 DDYI=  280518 DDXI=   35816 DDWI=   38640
  DDZE=       0 DDYE=   41898 DDXE=    6259 DDWE=    6729
================================================================================
Trace of MO density:    30.000000
   30  correlated and    12  frozen core electrons

          modens reordered block   1

               ag    1        ag    2        ag    3        ag    4        ag    5        ag    6        ag    7        ag    8
  ag    1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    3    0.00000        0.00000        1.98712      -3.035417E-05   1.113611E-03  -1.969814E-05   1.475387E-03  -7.786033E-05
  ag    4    0.00000        0.00000      -3.035417E-05    1.98075      -1.616371E-06  -2.182971E-03   1.233624E-04  -8.589133E-04
  ag    5    0.00000        0.00000       1.113611E-03  -1.616371E-06    1.97951       5.964333E-05   1.427275E-03  -1.403212E-04
  ag    6    0.00000        0.00000      -1.969814E-05  -2.182971E-03   5.964333E-05    1.97452      -1.123004E-04   1.735205E-03
  ag    7    0.00000        0.00000       1.475387E-03   1.233624E-04   1.427275E-03  -1.123004E-04   6.092263E-04  -1.979121E-06
  ag    8    0.00000        0.00000      -7.786033E-05  -8.589133E-04  -1.403212E-04   1.735205E-03  -1.979121E-06   2.206530E-04
  ag    9    0.00000        0.00000      -1.278903E-04   3.201460E-04   1.274714E-03  -2.425936E-04   6.055384E-04   3.479113E-06
  ag   10    0.00000        0.00000       2.796954E-04   2.253244E-03   1.287969E-04  -1.034638E-03   2.444594E-05  -3.969893E-04
  ag   11    0.00000        0.00000      -3.463828E-03  -2.656688E-04  -4.480291E-04   2.401285E-05  -9.755251E-04  -2.862220E-06
  ag   12    0.00000        0.00000       3.518598E-05  -2.493735E-04   3.581646E-04  -3.762960E-03  -2.061029E-06  -3.615660E-04
  ag   13    0.00000        0.00000       2.539784E-04   2.932914E-03  -1.470330E-04   1.691429E-03  -1.646496E-06  -1.800138E-04
  ag   14    0.00000        0.00000      -3.164824E-03  -2.627457E-04   3.461512E-03  -1.106323E-04  -1.234419E-03   3.748435E-06
  ag   15    0.00000        0.00000       4.620946E-04   1.598689E-03  -2.160537E-04   6.509083E-04  -3.300332E-06  -5.986942E-04
  ag   16    0.00000        0.00000      -1.063906E-04  -1.786430E-03   8.781936E-04  -6.072385E-03   6.399522E-07  -5.881680E-04
  ag   17    0.00000        0.00000      -3.640542E-03  -5.824378E-04   1.229545E-03  -2.653486E-04  -7.346438E-04   1.463427E-06
  ag   18    0.00000        0.00000      -5.563768E-04   1.875581E-03   8.064257E-04  -1.268102E-03   9.762518E-07   3.642264E-04
  ag   19    0.00000        0.00000      -3.540657E-04   1.498588E-04   1.397445E-03  -3.644389E-04   3.837203E-04   1.653364E-06
  ag   20    0.00000        0.00000       3.482954E-03  -1.651999E-05  -7.515446E-03   6.866236E-04   4.232980E-04   4.707330E-07
  ag   21    0.00000        0.00000       3.799397E-04   3.080224E-03  -1.563451E-04   1.645744E-03   8.981069E-07  -1.532611E-04
  ag   22    0.00000        0.00000       1.317822E-04   2.297848E-03  -5.302853E-04   1.464441E-03  -1.998761E-06   3.828244E-04
  ag   23    0.00000        0.00000      -1.368268E-04  -2.670671E-05  -4.379747E-04   1.231822E-03  -7.799104E-06   5.019305E-04
  ag   24    0.00000        0.00000      -4.576807E-04   2.535913E-03   4.797273E-04   8.102771E-04   1.028018E-05   9.812161E-05
  ag   25    0.00000        0.00000      -1.852067E-03   5.227734E-05   2.257715E-03  -5.562958E-04   7.148687E-04   4.214263E-06
  ag   26    0.00000        0.00000       3.441679E-05   1.046850E-03  -1.735711E-04  -8.022409E-05  -6.218819E-06   1.260753E-04
  ag   27    0.00000        0.00000      -1.469437E-03  -2.540804E-04   1.803888E-03  -2.279191E-04  -4.892020E-04  -5.910689E-07
  ag   28    0.00000        0.00000       3.758435E-04   2.269239E-03  -3.242666E-05  -1.072626E-03   1.090702E-04   3.149214E-06
  ag   29    0.00000        0.00000      -9.419924E-04   2.639653E-04  -2.961578E-04  -2.841594E-04  -3.955437E-04   2.266768E-06
  ag   30    0.00000        0.00000      -1.652474E-05   2.693363E-04   5.299307E-05   2.320405E-03  -6.137078E-07  -1.221777E-04
  ag   31    0.00000        0.00000       1.795721E-03   1.753118E-06   3.216997E-04   1.819490E-04  -1.540477E-04  -1.761403E-06
  ag   32    0.00000        0.00000      -1.349393E-03   4.192465E-06   4.474882E-03  -3.294595E-04   1.638784E-04  -7.446332E-07
  ag   33    0.00000        0.00000       1.191560E-05  -2.947196E-03   6.012689E-05  -3.792307E-03  -1.307519E-06   1.456148E-04
  ag   34    0.00000        0.00000      -2.778295E-04  -1.095643E-03   2.986079E-05   4.992333E-04  -7.125313E-07   1.486762E-04
  ag   35    0.00000        0.00000      -1.202200E-04   2.332033E-03  -3.307879E-05  -1.141506E-03   2.227461E-07  -1.276505E-06
  ag   36    0.00000        0.00000       5.266312E-05   9.276438E-04   2.684971E-04  -3.875803E-03  -1.381742E-05   6.648717E-06
  ag   37    0.00000        0.00000      -9.512298E-04   8.750006E-05  -2.620074E-03  -5.226583E-04   1.493877E-04   8.481511E-07
  ag   38    0.00000        0.00000      -9.678930E-06  -2.831383E-03   1.309500E-04  -2.907710E-03  -2.756764E-07  -1.023942E-04
  ag   39    0.00000        0.00000      -4.368706E-05   5.406194E-04   1.969046E-05   5.801510E-04   2.004608E-07   1.110096E-05

               ag    9        ag   10        ag   11        ag   12        ag   13        ag   14        ag   15        ag   16
  ag    3  -1.278903E-04   2.796954E-04  -3.463828E-03   3.518598E-05   2.539784E-04  -3.164824E-03   4.620946E-04  -1.063906E-04
  ag    4   3.201460E-04   2.253244E-03  -2.656688E-04  -2.493735E-04   2.932914E-03  -2.627457E-04   1.598689E-03  -1.786430E-03
  ag    5   1.274714E-03   1.287969E-04  -4.480291E-04   3.581646E-04  -1.470330E-04   3.461512E-03  -2.160537E-04   8.781936E-04
  ag    6  -2.425936E-04  -1.034638E-03   2.401285E-05  -3.762960E-03   1.691429E-03  -1.106323E-04   6.509083E-04  -6.072385E-03
  ag    7   6.055384E-04   2.444594E-05  -9.755251E-04  -2.061029E-06  -1.646496E-06  -1.234419E-03  -3.300332E-06   6.399522E-07
  ag    8   3.479113E-06  -3.969893E-04  -2.862220E-06  -3.615660E-04  -1.800138E-04   3.748435E-06  -5.986942E-04  -5.881680E-04
  ag    9   2.099123E-03   2.405633E-05  -4.915913E-04  -4.325282E-06  -1.094411E-05  -2.933356E-03  -2.472342E-05  -8.327518E-06
  ag   10   2.405633E-05   8.539443E-04  -2.326407E-05   5.295643E-04   5.563392E-04  -6.730190E-05   1.579604E-03   6.235217E-04
  ag   11  -4.915913E-04  -2.326407E-05   1.788892E-03   1.274041E-05   1.284731E-05   1.728806E-03   3.098599E-05   7.741708E-06
  ag   12  -4.325282E-06   5.295643E-04   1.274041E-05   7.532859E-04   1.695473E-05   1.489153E-06   6.480677E-04   1.551676E-03
  ag   13  -1.094411E-05   5.563392E-04   1.284731E-05   1.695473E-05   7.311041E-04   6.702489E-06   1.084832E-03  -5.336772E-04
  ag   14  -2.933356E-03  -6.730190E-05   1.728806E-03   1.489153E-06   6.702489E-06   5.801096E-03   1.215477E-05   1.027372E-05
  ag   15  -2.472342E-05   1.579604E-03   3.098599E-05   6.480677E-04   1.084832E-03   1.215477E-05   4.013303E-03   2.043077E-04
  ag   16  -8.327518E-06   6.235217E-04   7.741708E-06   1.551676E-03  -5.336772E-04   1.027372E-05   2.043077E-04   3.941926E-03
  ag   17   6.383406E-04  -2.086815E-05   1.787421E-03   9.997845E-06   3.725748E-06   6.171129E-04   4.529893E-06   7.071683E-06
  ag   18   6.708822E-06  -1.325975E-03  -2.426461E-05  -2.188721E-04  -9.919956E-04   1.910479E-05  -4.774794E-03   1.023508E-03
  ag   19   9.075543E-04   9.677176E-06  -3.063044E-04  -2.085253E-06  -4.671906E-06  -6.575865E-04  -1.919736E-05   3.673331E-06
  ag   20   1.447936E-03   2.695855E-05  -8.230819E-04  -3.087302E-06  -9.921627E-06  -4.165890E-03  -2.896400E-05  -1.224198E-05
  ag   21  -4.859358E-06   5.083137E-04   8.666298E-06  -1.127071E-04   9.453478E-04  -8.209938E-06   7.551762E-04  -1.002583E-03
  ag   22   1.009398E-08  -4.527995E-04  -4.223716E-06  -1.055075E-03   4.638225E-04   9.374396E-07  -4.924918E-04  -2.810145E-03
  ag   23  -6.897737E-06  -8.586141E-04  -5.634515E-06  -1.060329E-03  -1.531044E-04   1.819198E-05  -1.142299E-03  -2.590989E-03
  ag   24   2.037737E-05  -4.152802E-04  -1.878039E-05  -9.405876E-05  -3.170193E-04  -1.453768E-05  -2.145988E-03   6.239891E-04
  ag   25   1.698804E-03   2.603110E-05  -6.291186E-04  -1.324334E-05   2.271271E-06  -1.648351E-03   1.179403E-05  -3.975996E-05
  ag   26  -5.142827E-07  -1.898596E-04   9.451187E-06  -2.823577E-04  -4.852117E-05   1.494394E-05   4.252302E-05  -8.993702E-04
  ag   27  -4.391723E-04  -1.993939E-05   9.991242E-04   7.036607E-06   5.781590E-06   1.808348E-03   1.240721E-06   1.411050E-05
  ag   28  -1.486121E-04   8.647382E-05  -2.505708E-04  -2.722164E-04   5.604732E-04   1.424676E-04  -1.709181E-04  -1.105452E-03
  ag   29   5.510031E-04   1.346386E-05   9.237807E-04  -6.813483E-05   1.506652E-04  -5.565956E-04  -5.145243E-05  -2.955722E-04
  ag   30  -3.641274E-07   1.226062E-04   2.137555E-06   2.334249E-04   4.811415E-05  -7.181445E-06  -1.136034E-04   3.104235E-04
  ag   31  -7.019968E-04  -5.036940E-06   7.128535E-05   2.894834E-06   2.290262E-06   9.094569E-04   1.366675E-05   8.546584E-06
  ag   32   1.257295E-04   3.613519E-06  -1.870843E-04   1.737347E-06  -2.490484E-06   9.168503E-05  -3.798597E-06   8.212887E-06
  ag   33   3.181141E-06  -4.263006E-04  -5.523827E-06  -8.007415E-05  -4.211346E-04   4.662632E-06  -7.553511E-04   1.493852E-05
  ag   34  -1.841750E-08  -3.005793E-04  -4.669282E-06  -2.940424E-04  -2.372630E-04   6.212039E-06  -5.694703E-04  -5.242571E-04
  ag   35   1.962865E-07  -3.555582E-05  -1.702300E-06   9.384117E-06  -4.562672E-05  -2.392986E-07  -1.699796E-04   4.728977E-05
  ag   36  -3.112119E-05  -8.293502E-05   1.837179E-05   1.322915E-05  -9.017732E-05   6.577303E-05  -3.151497E-04   1.034043E-05
  ag   37   3.432271E-04  -1.153168E-06  -2.136062E-04  -1.200801E-07  -9.966140E-06  -6.852487E-04  -3.352705E-05  -2.558153E-06
  ag   38  -2.088663E-06   1.797058E-04   3.627676E-06   2.156176E-04   2.638214E-05   4.662079E-07   3.220520E-04   4.710287E-04
  ag   39   9.440282E-07  -4.775355E-05  -1.502457E-06   4.150977E-05  -8.574466E-05  -1.259298E-06  -1.246128E-04   2.401585E-04

               ag   17        ag   18        ag   19        ag   20        ag   21        ag   22        ag   23        ag   24
  ag    3  -3.640542E-03  -5.563768E-04  -3.540657E-04   3.482954E-03   3.799397E-04   1.317822E-04  -1.368268E-04  -4.576807E-04
  ag    4  -5.824378E-04   1.875581E-03   1.498588E-04  -1.651999E-05   3.080224E-03   2.297848E-03  -2.670671E-05   2.535913E-03
  ag    5   1.229545E-03   8.064257E-04   1.397445E-03  -7.515446E-03  -1.563451E-04  -5.302853E-04  -4.379747E-04   4.797273E-04
  ag    6  -2.653486E-04  -1.268102E-03  -3.644389E-04   6.866236E-04   1.645744E-03   1.464441E-03   1.231822E-03   8.102771E-04
  ag    7  -7.346438E-04   9.762518E-07   3.837203E-04   4.232980E-04   8.981069E-07  -1.998761E-06  -7.799104E-06   1.028018E-05
  ag    8   1.463427E-06   3.642264E-04   1.653364E-06   4.707330E-07  -1.532611E-04   3.828244E-04   5.019305E-04   9.812161E-05
  ag    9   6.383406E-04   6.708822E-06   9.075543E-04   1.447936E-03  -4.859358E-06   1.009398E-08  -6.897737E-06   2.037737E-05
  ag   10  -2.086815E-05  -1.325975E-03   9.677176E-06   2.695855E-05   5.083137E-04  -4.527995E-04  -8.586141E-04  -4.152802E-04
  ag   11   1.787421E-03  -2.426461E-05  -3.063044E-04  -8.230819E-04   8.666298E-06  -4.223716E-06  -5.634515E-06  -1.878039E-05
  ag   12   9.997845E-06  -2.188721E-04  -2.085253E-06  -3.087302E-06  -1.127071E-04  -1.055075E-03  -1.060329E-03  -9.405876E-05
  ag   13   3.725748E-06  -9.919956E-04  -4.671906E-06  -9.921627E-06   9.453478E-04   4.638225E-04  -1.531044E-04  -3.170193E-04
  ag   14   6.171129E-04   1.910479E-05  -6.575865E-04  -4.165890E-03  -8.209938E-06   9.374396E-07   1.819198E-05  -1.453768E-05
  ag   15   4.529893E-06  -4.774794E-03  -1.919736E-05  -2.896400E-05   7.551762E-04  -4.924918E-04  -1.142299E-03  -2.145988E-03
  ag   16   7.071683E-06   1.023508E-03   3.673331E-06  -1.224198E-05  -1.002583E-03  -2.810145E-03  -2.590989E-03   6.239891E-04
  ag   17   2.625268E-03  -3.946906E-06   4.381209E-04  -7.293089E-04   2.134177E-06  -5.654802E-06  -1.009121E-05  -9.873480E-07
  ag   18  -3.946906E-06   7.516371E-03   2.655179E-05   7.346985E-06  -5.484554E-04  -1.745962E-04   9.671899E-05   4.551067E-03
  ag   19   4.381209E-04   2.655179E-05   1.071829E-03  -5.230594E-04  -3.529639E-06  -2.648304E-06  -9.280609E-06   2.831470E-05
  ag   20  -7.293089E-04   7.346985E-06  -5.230594E-04   4.744895E-03   5.415091E-06   7.452478E-06   7.712100E-06   2.919450E-06
  ag   21   2.134177E-06  -5.484554E-04  -3.529639E-06   5.415091E-06   1.516204E-03   1.045224E-03   2.326606E-04  -2.497589E-04
  ag   22  -5.654802E-06  -1.745962E-04  -2.648304E-06   7.452478E-06   1.045224E-03   2.280103E-03   2.023150E-03  -3.215433E-04
  ag   23  -1.009121E-05   9.671899E-05  -9.280609E-06   7.712100E-06   2.326606E-04   2.023150E-03   2.675842E-03  -4.723693E-04
  ag   24  -9.873480E-07   4.551067E-03   2.831470E-05   2.919450E-06  -2.497589E-04  -3.215433E-04  -4.723693E-04   3.932492E-03
  ag   25   5.703334E-04  -5.652858E-05   1.179854E-03  -7.701185E-04   6.764683E-06   1.958131E-05   6.710001E-06  -1.504122E-05
  ag   26   1.255467E-05  -8.374400E-04  -1.947086E-06  -1.400814E-05   3.068875E-05   6.699650E-04   9.908653E-04  -9.590657E-04
  ag   27   1.226170E-03   2.032252E-05   3.346771E-04  -1.747160E-03   3.284810E-06  -4.555901E-06  -6.693623E-06   8.781699E-06
  ag   28  -3.693730E-04   5.128568E-04   4.147783E-05  -2.204749E-04   1.180470E-03   1.257988E-03   7.434332E-04   3.587008E-04
  ag   29   1.355710E-03   1.404108E-04  -1.544482E-04   8.410004E-04   3.194789E-04   3.350483E-04   1.923044E-04   1.000525E-04
  ag   30  -3.003758E-07   3.265867E-04  -4.135257E-06   1.534145E-05   2.618492E-04  -3.066414E-05   2.192007E-04   9.041890E-06
  ag   31  -4.928977E-04  -8.951993E-06  -7.149733E-04  -3.768658E-04  -2.557471E-06  -8.466120E-06  -6.520800E-06  -6.494791E-06
  ag   32   6.698547E-05   8.583491E-06   3.173612E-04  -1.177642E-03  -4.388624E-06  -7.642683E-06  -1.032348E-05   1.352960E-05
  ag   33  -4.241287E-07   5.703110E-05  -1.242454E-06  -1.399912E-07  -4.173911E-04   8.488647E-05   6.656992E-04  -9.706785E-04
  ag   34  -4.761215E-06   7.971962E-04   1.220994E-06  -1.653096E-06  -3.659827E-04   2.698085E-04   7.214970E-04   9.673343E-04
  ag   35  -2.786105E-06   4.018148E-04   7.936173E-07   3.253706E-06  -2.302324E-06   9.326647E-05   1.524719E-04   5.410032E-04
  ag   36   3.841633E-06   3.960620E-04  -1.110746E-05  -2.580230E-05  -2.026751E-05   1.402532E-04   2.699834E-04   5.135443E-05
  ag   37  -4.464052E-05   3.649057E-05   1.387483E-04   2.738702E-04  -1.536192E-06   1.431441E-05   2.252562E-05   1.207966E-05
  ag   38   2.383272E-06  -2.739053E-04  -1.897266E-06  -7.845243E-07  -5.640920E-05  -3.547074E-04  -4.072402E-04  -1.522887E-04
  ag   39  -5.960955E-07   2.140338E-04   8.099095E-07   2.594335E-06  -1.316114E-04  -1.800547E-04  -2.604513E-04   1.744315E-04

               ag   25        ag   26        ag   27        ag   28        ag   29        ag   30        ag   31        ag   32
  ag    3  -1.852067E-03   3.441679E-05  -1.469437E-03   3.758435E-04  -9.419924E-04  -1.652474E-05   1.795721E-03  -1.349393E-03
  ag    4   5.227734E-05   1.046850E-03  -2.540804E-04   2.269239E-03   2.639653E-04   2.693363E-04   1.753118E-06   4.192465E-06
  ag    5   2.257715E-03  -1.735711E-04   1.803888E-03  -3.242666E-05  -2.961578E-04   5.299307E-05   3.216997E-04   4.474882E-03
  ag    6  -5.562958E-04  -8.022409E-05  -2.279191E-04  -1.072626E-03  -2.841594E-04   2.320405E-03   1.819490E-04  -3.294595E-04
  ag    7   7.148687E-04  -6.218819E-06  -4.892020E-04   1.090702E-04  -3.955437E-04  -6.137078E-07  -1.540477E-04   1.638784E-04
  ag    8   4.214263E-06   1.260753E-04  -5.910689E-07   3.149214E-06   2.266768E-06  -1.221777E-04  -1.761403E-06  -7.446332E-07
  ag    9   1.698804E-03  -5.142827E-07  -4.391723E-04  -1.486121E-04   5.510031E-04  -3.641274E-07  -7.019968E-04   1.257295E-04
  ag   10   2.603110E-05  -1.898596E-04  -1.993939E-05   8.647382E-05   1.346386E-05   1.226062E-04  -5.036940E-06   3.613519E-06
  ag   11  -6.291186E-04   9.451187E-06   9.991242E-04  -2.505708E-04   9.237807E-04   2.137555E-06   7.128535E-05  -1.870843E-04
  ag   12  -1.324334E-05  -2.823577E-04   7.036607E-06  -2.722164E-04  -6.813483E-05   2.334249E-04   2.894834E-06   1.737347E-06
  ag   13   2.271271E-06  -4.852117E-05   5.781590E-06   5.604732E-04   1.506652E-04   4.811415E-05   2.290262E-06  -2.490484E-06
  ag   14  -1.648351E-03   1.494394E-05   1.808348E-03   1.424676E-04  -5.565956E-04  -7.181445E-06   9.094569E-04   9.168503E-05
  ag   15   1.179403E-05   4.252302E-05   1.240721E-06  -1.709181E-04  -5.145243E-05  -1.136034E-04   1.366675E-05  -3.798597E-06
  ag   16  -3.975996E-05  -8.993702E-04   1.411050E-05  -1.105452E-03  -2.955722E-04   3.104235E-04   8.546584E-06   8.212887E-06
  ag   17   5.703334E-04   1.255467E-05   1.226170E-03  -3.693730E-04   1.355710E-03  -3.003758E-07  -4.928977E-04   6.698547E-05
  ag   18  -5.652858E-05  -8.374400E-04   2.032252E-05   5.128568E-04   1.404108E-04   3.265867E-04  -8.951993E-06   8.583491E-06
  ag   19   1.179854E-03  -1.947086E-06   3.346771E-04   4.147783E-05  -1.544482E-04  -4.135257E-06  -7.149733E-04   3.173612E-04
  ag   20  -7.701185E-04  -1.400814E-05  -1.747160E-03  -2.204749E-04   8.410004E-04   1.534145E-05  -3.768658E-04  -1.177642E-03
  ag   21   6.764683E-06   3.068875E-05   3.284810E-06   1.180470E-03   3.194789E-04   2.618492E-04  -2.557471E-06  -4.388624E-06
  ag   22   1.958131E-05   6.699650E-04  -4.555901E-06   1.257988E-03   3.350483E-04  -3.066414E-05  -8.466120E-06  -7.642683E-06
  ag   23   6.710001E-06   9.908653E-04  -6.693623E-06   7.434332E-04   1.923044E-04   2.192007E-04  -6.520800E-06  -1.032348E-05
  ag   24  -1.504122E-05  -9.590657E-04   8.781699E-06   3.587008E-04   1.000525E-04   9.041890E-06  -6.494791E-06   1.352960E-05
  ag   25   3.610864E-03   2.184041E-05  -2.493847E-04  -2.034310E-05   1.176605E-04  -6.218913E-06  -2.715556E-04   1.243384E-03
  ag   26   2.184041E-05   8.500641E-04   3.571526E-06   9.538507E-05   2.745636E-05   5.897280E-05  -3.690333E-06   1.615150E-06
  ag   27  -2.493847E-04   3.571526E-06   1.524377E-03  -4.805879E-05   1.782578E-04  -2.327401E-06  -3.981252E-04   3.733436E-04
  ag   28  -2.034310E-05   9.538507E-05  -4.805879E-05   1.560710E-03   3.551276E-05   3.100570E-04   5.069147E-05   8.397728E-05
  ag   29   1.176605E-04   2.745636E-05   1.782578E-04   3.551276E-05   1.430242E-03   8.555695E-05  -2.002022E-04  -3.127551E-04
  ag   30  -6.218913E-06   5.897280E-05  -2.327401E-06   3.100570E-04   8.555695E-05   7.664983E-04   1.041366E-07  -6.421404E-06
  ag   31  -2.715556E-04  -3.690333E-06  -3.981252E-04   5.069147E-05  -2.002022E-04   1.041366E-07   1.155487E-03  -3.006184E-04
  ag   32   1.243384E-03   1.615150E-06   3.733436E-04   8.397728E-05  -3.127551E-04  -6.421404E-06  -3.006184E-04   1.737319E-03
  ag   33   1.591096E-05   3.188150E-04  -5.151025E-06  -3.374721E-04  -9.181987E-05   1.002327E-04   9.753535E-07   8.915185E-07
  ag   34  -3.524466E-06  -2.414481E-05  -1.583775E-06  -6.221221E-05  -1.973143E-05   2.519613E-05   4.760776E-07   1.357723E-06
  ag   35  -8.040425E-06  -2.347433E-04   2.945988E-06   3.929779E-04   1.037369E-04   2.646377E-04  -9.368643E-07  -2.686740E-06
  ag   36  -7.308945E-05   2.390447E-04   2.842067E-05   1.811998E-04   3.590246E-05   1.445606E-04   1.501252E-05  -4.270839E-05
  ag   37   7.890243E-04   2.161552E-05  -3.017501E-04  -7.681971E-06   1.119542E-04   1.333670E-05  -1.864549E-04   4.583950E-04
  ag   38  -4.557404E-06  -2.276252E-04   2.498374E-06  -9.809884E-05  -2.473835E-05   5.806291E-05   2.908305E-06  -1.577919E-06
  ag   39  -6.880147E-06  -1.327140E-04   1.339141E-06  -2.369875E-05  -5.995095E-06   4.709362E-06  -2.494954E-07  -5.081059E-07

               ag   33        ag   34        ag   35        ag   36        ag   37        ag   38        ag   39
  ag    3   1.191560E-05  -2.778295E-04  -1.202200E-04   5.266312E-05  -9.512298E-04  -9.678930E-06  -4.368706E-05
  ag    4  -2.947196E-03  -1.095643E-03   2.332033E-03   9.276438E-04   8.750006E-05  -2.831383E-03   5.406194E-04
  ag    5   6.012689E-05   2.986079E-05  -3.307879E-05   2.684971E-04  -2.620074E-03   1.309500E-04   1.969046E-05
  ag    6  -3.792307E-03   4.992333E-04  -1.141506E-03  -3.875803E-03  -5.226583E-04  -2.907710E-03   5.801510E-04
  ag    7  -1.307519E-06  -7.125313E-07   2.227461E-07  -1.381742E-05   1.493877E-04  -2.756764E-07   2.004608E-07
  ag    8   1.456148E-04   1.486762E-04  -1.276505E-06   6.648717E-06   8.481511E-07  -1.023942E-04   1.110096E-05
  ag    9   3.181141E-06  -1.841750E-08   1.962865E-07  -3.112119E-05   3.432271E-04  -2.088663E-06   9.440282E-07
  ag   10  -4.263006E-04  -3.005793E-04  -3.555582E-05  -8.293502E-05  -1.153168E-06   1.797058E-04  -4.775355E-05
  ag   11  -5.523827E-06  -4.669282E-06  -1.702300E-06   1.837179E-05  -2.136062E-04   3.627676E-06  -1.502457E-06
  ag   12  -8.007415E-05  -2.940424E-04   9.384117E-06   1.322915E-05  -1.200801E-07   2.156176E-04   4.150977E-05
  ag   13  -4.211346E-04  -2.372630E-04  -4.562672E-05  -9.017732E-05  -9.966140E-06   2.638214E-05  -8.574466E-05
  ag   14   4.662632E-06   6.212039E-06  -2.392986E-07   6.577303E-05  -6.852487E-04   4.662079E-07  -1.259298E-06
  ag   15  -7.553511E-04  -5.694703E-04  -1.699796E-04  -3.151497E-04  -3.352705E-05   3.220520E-04  -1.246128E-04
  ag   16   1.493852E-05  -5.242571E-04   4.728977E-05   1.034043E-05  -2.558153E-06   4.710287E-04   2.401585E-04
  ag   17  -4.241287E-07  -4.761215E-06  -2.786105E-06   3.841633E-06  -4.464052E-05   2.383272E-06  -5.960955E-07
  ag   18   5.703110E-05   7.971962E-04   4.018148E-04   3.960620E-04   3.649057E-05  -2.739053E-04   2.140338E-04
  ag   19  -1.242454E-06   1.220994E-06   7.936173E-07  -1.110746E-05   1.387483E-04  -1.897266E-06   8.099095E-07
  ag   20  -1.399912E-07  -1.653096E-06   3.253706E-06  -2.580230E-05   2.738702E-04  -7.845243E-07   2.594335E-06
  ag   21  -4.173911E-04  -3.659827E-04  -2.302324E-06  -2.026751E-05  -1.536192E-06  -5.640920E-05  -1.316114E-04
  ag   22   8.488647E-05   2.698085E-04   9.326647E-05   1.402532E-04   1.431441E-05  -3.547074E-04  -1.800547E-04
  ag   23   6.656992E-04   7.214970E-04   1.524719E-04   2.699834E-04   2.252562E-05  -4.072402E-04  -2.604513E-04
  ag   24  -9.706785E-04   9.673343E-04   5.410032E-04   5.135443E-05   1.207966E-05  -1.522887E-04   1.744315E-04
  ag   25   1.591096E-05  -3.524466E-06  -8.040425E-06  -7.308945E-05   7.890243E-04  -4.557404E-06  -6.880147E-06
  ag   26   3.188150E-04  -2.414481E-05  -2.347433E-04   2.390447E-04   2.161552E-05  -2.276252E-04  -1.327140E-04
  ag   27  -5.151025E-06  -1.583775E-06   2.945988E-06   2.842067E-05  -3.017501E-04   2.498374E-06   1.339141E-06
  ag   28  -3.374721E-04  -6.221221E-05   3.929779E-04   1.811998E-04  -7.681971E-06  -9.809884E-05  -2.369875E-05
  ag   29  -9.181987E-05  -1.973143E-05   1.037369E-04   3.590246E-05   1.119542E-04  -2.473835E-05  -5.995095E-06
  ag   30   1.002327E-04   2.519613E-05   2.646377E-04   1.445606E-04   1.333670E-05   5.806291E-05   4.709362E-06
  ag   31   9.753535E-07   4.760776E-07  -9.368643E-07   1.501252E-05  -1.864549E-04   2.908305E-06  -2.494954E-07
  ag   32   8.915185E-07   1.357723E-06  -2.686740E-06  -4.270839E-05   4.583950E-04  -1.577919E-06  -5.081059E-07
  ag   33   1.472520E-03  -9.979384E-05  -3.024557E-04   2.959747E-04   2.896656E-05   2.003676E-04  -1.940780E-04
  ag   34  -9.979384E-05   1.040466E-03   2.999652E-04   1.602133E-05   1.950769E-06  -6.406125E-05  -1.223579E-04
  ag   35  -3.024557E-04   2.999652E-04   8.601259E-04   2.367983E-05   2.637971E-06   1.264326E-04   2.291500E-04
  ag   36   2.959747E-04   1.602133E-05   2.367983E-05   4.079820E-04  -1.325534E-05   9.058779E-05  -2.077571E-05
  ag   37   2.896656E-05   1.950769E-06   2.637971E-06  -1.325534E-05   5.585454E-04   8.103228E-06  -2.050649E-06
  ag   38   2.003676E-04  -6.406125E-05   1.264326E-04   9.058779E-05   8.103228E-06   4.994466E-04   1.166777E-04
  ag   39  -1.940780E-04  -1.223579E-04   2.291500E-04  -2.077571E-05  -2.050649E-06   1.166777E-04   3.243597E-04

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98730243     1.98147408     1.97941805     1.97389910     0.01446501     0.01265214
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.01033253     0.00666770     0.00532211     0.00463655     0.00232900     0.00177666     0.00148538     0.00127368
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00092170     0.00077516     0.00046800     0.00046137     0.00029974     0.00025221     0.00013828     0.00011880
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00010063     0.00009414     0.00008577     0.00005038     0.00002131     0.00002029     0.00001484     0.00000808
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000639     0.00000550     0.00000272     0.00000128     0.00000021     0.00000010     0.00000003

          modens reordered block   1

               b3u   1        b3u   2        b3u   3        b3u   4        b3u   5        b3u   6        b3u   7        b3u   8
  b3u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   3    0.00000        0.00000        1.98451      -3.088796E-05   1.745385E-03  -3.992931E-04   9.846881E-05   1.031385E-03
  b3u   4    0.00000        0.00000      -3.088796E-05    1.97879       6.327188E-05   3.305085E-04  -1.031970E-03  -4.244584E-04
  b3u   5    0.00000        0.00000       1.745385E-03   6.327188E-05    1.97702      -2.582297E-03   7.361762E-05   2.680922E-03
  b3u   6    0.00000        0.00000      -3.992931E-04   3.305085E-04  -2.582297E-03   3.830055E-04  -1.754730E-06  -4.916407E-04
  b3u   7    0.00000        0.00000       9.846881E-05  -1.031970E-03   7.361762E-05  -1.754730E-06   1.027336E-04   3.418100E-06
  b3u   8    0.00000        0.00000       1.031385E-03  -4.244584E-04   2.680922E-03  -4.916407E-04   3.418100E-06   6.980272E-04
  b3u   9    0.00000        0.00000      -3.619658E-05  -4.961821E-04   2.893909E-03  -5.957967E-04   4.880613E-06   6.697976E-04
  b3u  10    0.00000        0.00000      -3.039879E-04   2.312506E-03  -1.548599E-04  -2.761608E-06  -3.231508E-04  -2.491641E-07
  b3u  11    0.00000        0.00000      -1.712902E-03   1.289991E-04  -2.211957E-05   2.465781E-04  -1.439204E-06  -4.982773E-04
  b3u  12    0.00000        0.00000       4.126757E-04  -1.753150E-03   1.710419E-04  -3.481847E-05   4.359178E-04   3.987876E-05
  b3u  13    0.00000        0.00000       1.606377E-04   6.295438E-05   6.435447E-04   9.255934E-04   8.063290E-06  -9.055996E-04
  b3u  14    0.00000        0.00000      -4.485474E-04   1.785085E-04  -8.588444E-05   3.282251E-07  -4.316704E-04  -8.324471E-06
  b3u  15    0.00000        0.00000       3.907250E-04  -8.375453E-04   4.478841E-03  -5.730995E-04  -4.434396E-06   7.940512E-04
  b3u  16    0.00000        0.00000      -2.254848E-03  -4.203826E-04   4.500484E-03   2.344133E-04  -2.785320E-07  -3.862534E-04
  b3u  17    0.00000        0.00000       9.943448E-04  -4.418403E-04   1.535875E-03  -1.908961E-04  -5.266944E-07   4.733170E-04
  b3u  18    0.00000        0.00000      -7.185500E-04  -4.001724E-03   5.756421E-04   6.766159E-06  -4.080986E-04  -1.265538E-05
  b3u  19    0.00000        0.00000      -5.574052E-04   6.717587E-05  -2.273755E-03  -3.274646E-04  -3.554315E-06   1.934449E-04
  b3u  20    0.00000        0.00000      -1.138973E-04   1.491094E-03  -3.952121E-04   3.228597E-06  -2.614757E-04  -6.541608E-06
  b3u  21    0.00000        0.00000      -6.332686E-04  -6.157362E-04   1.322368E-03  -7.062691E-04   4.093350E-08   9.166193E-04
  b3u  22    0.00000        0.00000      -4.327841E-04  -3.628507E-03   9.502308E-04  -2.235015E-07   2.161633E-04   2.627017E-06
  b3u  23    0.00000        0.00000      -4.989595E-04  -1.386714E-03   5.121586E-04  -2.150200E-06  -1.158516E-04   1.860634E-06
  b3u  24    0.00000        0.00000       1.669795E-03   2.513377E-04  -2.384903E-03   1.668975E-05   5.117047E-07   3.910201E-05
  b3u  25    0.00000        0.00000      -1.372226E-04  -2.102014E-04  -1.016656E-04  -3.442241E-04   1.100576E-06   4.686073E-04
  b3u  26    0.00000        0.00000       3.422661E-04   3.004462E-05  -9.920528E-04  -1.757441E-04   1.616469E-06   1.432488E-04
  b3u  27    0.00000        0.00000      -5.919283E-05   2.265599E-03   3.130948E-04   1.390737E-06   2.031764E-04  -1.648578E-06
  b3u  28    0.00000        0.00000      -8.125113E-04   1.211194E-04  -2.384279E-03   4.504246E-06  -4.096590E-06  -1.717849E-04
  b3u  29    0.00000        0.00000      -1.855628E-03   7.872779E-05   1.523730E-03   2.495424E-04  -2.252798E-06  -4.565220E-04
  b3u  30    0.00000        0.00000      -9.292524E-05  -1.351159E-03  -2.849468E-05  -4.209406E-06  -7.208103E-05   6.096897E-06
  b3u  31    0.00000        0.00000       2.299247E-04   1.309539E-04  -4.125564E-03   9.256903E-05  -1.031177E-06  -5.861270E-05
  b3u  32    0.00000        0.00000       1.139578E-03  -7.116561E-05   2.392050E-04  -6.593844E-05   2.164720E-07   9.704114E-05
  b3u  33    0.00000        0.00000       1.873997E-04  -2.404308E-03  -1.770965E-04  -4.030058E-08   9.831230E-05   7.859786E-07
  b3u  34    0.00000        0.00000      -1.807575E-03  -1.200168E-04  -9.667974E-04  -1.460900E-04   3.164994E-07   1.970640E-04
  b3u  35    0.00000        0.00000       2.603154E-04  -3.041564E-03   4.670420E-04  -6.170446E-07  -1.169225E-05   9.954938E-06
  b3u  36    0.00000        0.00000      -8.215743E-04  -1.147934E-03  -1.419615E-03  -1.278852E-06  -4.659070E-06  -2.191931E-05
  b3u  37    0.00000        0.00000      -2.043536E-03  -6.647869E-05   2.593740E-03  -4.845054E-05   8.557658E-08   6.157815E-05
  b3u  38    0.00000        0.00000       3.549565E-05   1.111245E-07   2.503925E-05   2.497427E-07   6.047514E-05   3.637695E-07
  b3u  39    0.00000        0.00000       2.660460E-05   1.056629E-03   1.854419E-05   4.195962E-07   1.097665E-05  -4.608264E-07

               b3u   9        b3u  10        b3u  11        b3u  12        b3u  13        b3u  14        b3u  15        b3u  16
  b3u   3  -3.619658E-05  -3.039879E-04  -1.712902E-03   4.126757E-04   1.606377E-04  -4.485474E-04   3.907250E-04  -2.254848E-03
  b3u   4  -4.961821E-04   2.312506E-03   1.289991E-04  -1.753150E-03   6.295438E-05   1.785085E-04  -8.375453E-04  -4.203826E-04
  b3u   5   2.893909E-03  -1.548599E-04  -2.211957E-05   1.710419E-04   6.435447E-04  -8.588444E-05   4.478841E-03   4.500484E-03
  b3u   6  -5.957967E-04  -2.761608E-06   2.465781E-04  -3.481847E-05   9.255934E-04   3.282251E-07  -5.730995E-04   2.344133E-04
  b3u   7   4.880613E-06  -3.231508E-04  -1.439204E-06   4.359178E-04   8.063290E-06  -4.316704E-04  -4.434396E-06  -2.785320E-07
  b3u   8   6.697976E-04  -2.491641E-07  -4.982773E-04   3.987876E-05  -9.055996E-04  -8.324471E-06   7.940512E-04  -3.862534E-04
  b3u   9   1.136184E-03  -2.233503E-06  -7.190726E-05   8.670714E-05  -2.085348E-03  -6.672788E-06   9.119749E-04  -2.461378E-04
  b3u  10  -2.233503E-06   1.035769E-03  -3.758269E-07  -1.421318E-03  -4.745120E-05   1.478149E-03   2.857502E-05  -4.392272E-06
  b3u  11  -7.190726E-05  -3.758269E-07   9.463490E-04  -8.879065E-06   3.192028E-04  -6.894573E-06   2.513845E-04   1.541302E-03
  b3u  12   8.670714E-05  -1.421318E-03  -8.879065E-06   1.999729E-03  -1.540987E-04  -2.140368E-03  -3.317911E-05  -9.505615E-05
  b3u  13  -2.085348E-03  -4.745120E-05   3.192028E-04  -1.540987E-04   6.141398E-03  -1.120459E-04   4.721818E-04   3.470365E-03
  b3u  14  -6.672788E-06   1.478149E-03  -6.894573E-06  -2.140368E-03  -1.120459E-04   2.628523E-03  -7.202167E-06  -7.973500E-05
  b3u  15   9.119749E-04   2.857502E-05   2.513845E-04  -3.317911E-05   4.721818E-04  -7.202167E-06   3.144606E-03   3.009095E-03
  b3u  16  -2.461378E-04  -4.392272E-06   1.541302E-03  -9.505615E-05   3.470365E-03  -7.973500E-05   3.009095E-03   6.243427E-03
  b3u  17  -1.356208E-04   6.661886E-06  -6.908910E-04  -5.374051E-05   1.352507E-03  -3.790804E-06   7.894237E-04   4.372406E-04
  b3u  18  -2.741947E-05   1.550578E-03   2.645144E-06  -2.446970E-03  -1.386265E-05   3.609416E-03   5.413167E-05   1.086362E-05
  b3u  19   9.795512E-04   2.254572E-05   2.489198E-04   7.507746E-05  -3.107699E-03   7.586211E-05  -3.940564E-04  -1.490246E-03
  b3u  20  -1.187918E-05   7.797411E-04  -2.800819E-06  -1.022161E-03  -2.602909E-05   6.929084E-04  -4.923124E-06  -1.835046E-05
  b3u  21   1.450113E-03   1.817463E-05   2.203828E-04   5.642213E-05  -1.667554E-03  -5.929755E-06   3.015413E-03   2.126777E-03
  b3u  22   8.849240E-07  -5.970211E-04  -3.970137E-06   6.733717E-04   2.959595E-05  -6.968164E-05  -8.961318E-06  -1.249353E-05
  b3u  23  -1.118171E-06   4.528786E-04  -4.395122E-06  -7.720644E-04  -2.444938E-05   1.159087E-03   1.773992E-05  -1.798135E-05
  b3u  24  -8.622374E-05  -2.955238E-06  -5.834989E-04   3.688989E-05  -1.193061E-03   2.087481E-05  -1.420765E-03  -3.106007E-03
  b3u  25   4.809140E-04   3.717670E-06  -5.236877E-04   3.782858E-05  -1.192573E-03   1.265385E-05  -1.372644E-04  -1.438737E-03
  b3u  26   2.951078E-04  -2.328876E-06   9.425090E-05   1.753261E-05  -3.039276E-04  -7.816463E-06   7.746525E-05   2.186137E-04
  b3u  27   9.217600E-06  -7.135615E-04   9.109885E-06   1.065986E-03   6.249013E-06  -1.228902E-03  -2.101665E-05   1.704380E-06
  b3u  28   4.026020E-04   1.514301E-05   5.964009E-04   3.041504E-05  -1.465740E-03   3.532425E-05   1.537080E-05   1.182811E-04
  b3u  29  -2.827317E-04   1.067125E-06   5.501672E-04  -2.058337E-05   4.368988E-04   4.848447E-06  -5.067472E-04   4.122039E-04
  b3u  30   4.385558E-06   1.680474E-04  -6.497092E-06  -1.863431E-04  -1.385160E-05  -9.762305E-05   7.227666E-06  -3.365851E-06
  b3u  31  -2.130676E-04   8.387061E-07  -3.054661E-05  -1.802686E-05   4.648197E-04  -3.700933E-06   8.972279E-06  -9.075957E-05
  b3u  32   1.141418E-04   1.002105E-06  -2.277361E-05   5.817368E-06  -1.698136E-04  -2.859244E-07   2.024691E-04  -3.348894E-04
  b3u  33   2.820353E-06  -3.342354E-04   1.407955E-06   5.163864E-04   1.326114E-05  -5.063310E-04  -8.542531E-06   2.274990E-06
  b3u  34   2.986885E-04   2.454853E-06  -5.222919E-05   1.944923E-05  -5.416866E-04  -5.880563E-07   4.351737E-04   4.591753E-05
  b3u  35  -2.826726E-05   8.351704E-05  -2.657456E-05  -1.972718E-04   1.285000E-04   3.273312E-04   3.411260E-05   3.325173E-05
  b3u  36   7.459971E-05   3.331491E-05   6.547398E-05  -6.465893E-05  -3.433554E-04   1.332220E-04  -6.229969E-05  -8.452454E-05
  b3u  37   7.552988E-05   8.504603E-07   3.813686E-05  -2.627663E-06   1.047232E-04  -3.803087E-06   2.909510E-04   4.418524E-04
  b3u  38   9.571098E-07  -2.052039E-04  -4.843116E-08   2.956760E-04   8.344638E-06  -3.402965E-04  -6.415089E-06  -6.582789E-08
  b3u  39  -4.403723E-07  -3.160563E-05   6.400202E-08   3.910149E-05   1.409857E-06  -5.526535E-05  -2.294901E-06  -5.230964E-07

               b3u  17        b3u  18        b3u  19        b3u  20        b3u  21        b3u  22        b3u  23        b3u  24
  b3u   3   9.943448E-04  -7.185500E-04  -5.574052E-04  -1.138973E-04  -6.332686E-04  -4.327841E-04  -4.989595E-04   1.669795E-03
  b3u   4  -4.418403E-04  -4.001724E-03   6.717587E-05   1.491094E-03  -6.157362E-04  -3.628507E-03  -1.386714E-03   2.513377E-04
  b3u   5   1.535875E-03   5.756421E-04  -2.273755E-03  -3.952121E-04   1.322368E-03   9.502308E-04   5.121586E-04  -2.384903E-03
  b3u   6  -1.908961E-04   6.766159E-06  -3.274646E-04   3.228597E-06  -7.062691E-04  -2.235015E-07  -2.150200E-06   1.668975E-05
  b3u   7  -5.266944E-07  -4.080986E-04  -3.554315E-06  -2.614757E-04   4.093350E-08   2.161633E-04  -1.158516E-04   5.117047E-07
  b3u   8   4.733170E-04  -1.265538E-05   1.934449E-04  -6.541608E-06   9.166193E-04   2.627017E-06   1.860634E-06   3.910201E-05
  b3u   9  -1.356208E-04  -2.741947E-05   9.795512E-04  -1.187918E-05   1.450113E-03   8.849240E-07  -1.118171E-06  -8.622374E-05
  b3u  10   6.661886E-06   1.550578E-03   2.254572E-05   7.797411E-04   1.817463E-05  -5.970211E-04   4.528786E-04  -2.955238E-06
  b3u  11  -6.908910E-04   2.645144E-06   2.489198E-04  -2.800819E-06   2.203828E-04  -3.970137E-06  -4.395122E-06  -5.834989E-04
  b3u  12  -5.374051E-05  -2.446970E-03   7.507746E-05  -1.022161E-03   5.642213E-05   6.733717E-04  -7.720644E-04   3.688989E-05
  b3u  13   1.352507E-03  -1.386265E-05  -3.107699E-03  -2.602909E-05  -1.667554E-03   2.959595E-05  -2.444938E-05  -1.193061E-03
  b3u  14  -3.790804E-06   3.609416E-03   7.586211E-05   6.929084E-04  -5.929755E-06  -6.968164E-05   1.159087E-03   2.087481E-05
  b3u  15   7.894237E-04   5.413167E-05  -3.940564E-04  -4.923124E-06   3.015413E-03  -8.961318E-06   1.773992E-05  -1.420765E-03
  b3u  16   4.372406E-04   1.086362E-05  -1.490246E-03  -1.835046E-05   2.126777E-03  -1.249353E-05  -1.798135E-05  -3.106007E-03
  b3u  17   1.266767E-03   4.410987E-05  -1.112548E-03  -3.155976E-06   1.177659E-04   2.018661E-05   1.713561E-05  -3.817602E-04
  b3u  18   4.410987E-05   6.637895E-03   6.281746E-05  -1.592429E-04   2.139520E-05   2.587615E-03   2.895004E-03  -3.961539E-05
  b3u  19  -1.112548E-03   6.281746E-05   1.955493E-03  -5.337961E-06   5.061471E-04   3.521655E-05   4.109228E-05   7.350556E-04
  b3u  20  -3.155976E-06  -1.592429E-04  -5.337961E-06   1.144344E-03  -2.475529E-05  -1.751477E-03  -3.148954E-04   1.881533E-05
  b3u  21   1.177659E-04   2.139520E-05   5.061471E-04  -2.475529E-05   4.936542E-03  -2.995532E-06   1.748542E-05  -1.549581E-03
  b3u  22   2.018661E-05   2.587615E-03   3.521655E-05  -1.751477E-03  -2.995532E-06   5.346272E-03   2.641753E-03  -2.606330E-05
  b3u  23   1.713561E-05   2.895004E-03   4.109228E-05  -3.148954E-04   1.748542E-05   2.641753E-03   2.115707E-03  -1.214119E-05
  b3u  24  -3.817602E-04  -3.961539E-05   7.350556E-04   1.881533E-05  -1.549581E-03  -2.606330E-05  -1.214119E-05   2.487912E-03
  b3u  25   2.549502E-04  -1.316167E-05   5.760441E-04   8.581414E-06  -5.094180E-04  -3.878122E-06   7.242442E-07   7.268863E-04
  b3u  26   6.437206E-05  -1.516558E-05   2.779147E-04   2.019476E-06  -6.220662E-04   3.371230E-06  -5.405697E-06  -2.225291E-04
  b3u  27  -2.247135E-05  -1.061474E-03   5.826910E-06  -7.511507E-04   1.420904E-05   1.927976E-03   4.399689E-04  -6.280922E-06
  b3u  28  -1.027654E-03   1.048633E-05   1.149930E-03   7.024591E-06   1.090239E-03  -3.850449E-05  -7.000010E-06  -9.161314E-05
  b3u  29  -5.037587E-04   1.004118E-05   2.651300E-04   7.434211E-06  -1.300510E-03   1.736004E-05   8.809300E-06   1.075532E-04
  b3u  30   4.788672E-06  -3.620100E-04  -8.339535E-06   3.508365E-04   1.339853E-05   2.284094E-04   2.967738E-04  -5.153173E-06
  b3u  31   2.525637E-05  -1.878174E-06  -4.720377E-05   3.598884E-06  -5.282798E-05   4.889740E-07   3.705134E-06   4.867165E-04
  b3u  32  -5.726083E-05  -1.812439E-06   1.544504E-04   3.208884E-06  -2.000635E-04  -1.536510E-06   3.620473E-06   8.965045E-04
  b3u  33  -5.436846E-06  -7.894060E-04  -7.123812E-06  -2.946190E-04  -8.973486E-06  -3.044402E-04  -7.437362E-04   8.332333E-06
  b3u  34   4.625567E-05  -6.077453E-06   2.587762E-04  -2.022103E-06   7.349246E-04  -4.612062E-06  -4.988754E-07  -1.607401E-04
  b3u  35   8.029721E-05   7.567264E-04  -5.373672E-05  -3.499863E-06  -1.030601E-04   1.765466E-04   3.279841E-04   2.484030E-05
  b3u  36  -1.878954E-04   2.949092E-04   1.658427E-04  -2.810661E-06   2.890222E-04   6.582929E-05   1.276625E-04  -7.843231E-05
  b3u  37   1.094727E-04   4.640282E-06  -8.538886E-05  -3.153107E-06   3.116549E-04   4.420759E-06   1.848372E-06  -4.754681E-04
  b3u  38  -1.812671E-06  -3.468202E-04  -4.110701E-06  -1.548598E-04  -3.766411E-06   2.820680E-04  -3.944223E-06   4.198518E-07
  b3u  39  -1.344912E-06  -1.774830E-04  -2.707301E-06   4.329790E-05  -3.048545E-06  -2.973415E-04  -1.443576E-04   2.531054E-06

               b3u  25        b3u  26        b3u  27        b3u  28        b3u  29        b3u  30        b3u  31        b3u  32
  b3u   3  -1.372226E-04   3.422661E-04  -5.919283E-05  -8.125113E-04  -1.855628E-03  -9.292524E-05   2.299247E-04   1.139578E-03
  b3u   4  -2.102014E-04   3.004462E-05   2.265599E-03   1.211194E-04   7.872779E-05  -1.351159E-03   1.309539E-04  -7.116561E-05
  b3u   5  -1.016656E-04  -9.920528E-04   3.130948E-04  -2.384279E-03   1.523730E-03  -2.849468E-05  -4.125564E-03   2.392050E-04
  b3u   6  -3.442241E-04  -1.757441E-04   1.390737E-06   4.504246E-06   2.495424E-04  -4.209406E-06   9.256903E-05  -6.593844E-05
  b3u   7   1.100576E-06   1.616469E-06   2.031764E-04  -4.096590E-06  -2.252798E-06  -7.208103E-05  -1.031177E-06   2.164720E-07
  b3u   8   4.686073E-04   1.432488E-04  -1.648578E-06  -1.717849E-04  -4.565220E-04   6.096897E-06  -5.861270E-05   9.704114E-05
  b3u   9   4.809140E-04   2.951078E-04   9.217600E-06   4.026020E-04  -2.827317E-04   4.385558E-06  -2.130676E-04   1.141418E-04
  b3u  10   3.717670E-06  -2.328876E-06  -7.135615E-04   1.514301E-05   1.067125E-06   1.680474E-04   8.387061E-07   1.002105E-06
  b3u  11  -5.236877E-04   9.425090E-05   9.109885E-06   5.964009E-04   5.501672E-04  -6.497092E-06  -3.054661E-05  -2.277361E-05
  b3u  12   3.782858E-05   1.753261E-05   1.065986E-03   3.041504E-05  -2.058337E-05  -1.863431E-04  -1.802686E-05   5.817368E-06
  b3u  13  -1.192573E-03  -3.039276E-04   6.249013E-06  -1.465740E-03   4.368988E-04  -1.385160E-05   4.648197E-04  -1.698136E-04
  b3u  14   1.265385E-05  -7.816463E-06  -1.228902E-03   3.532425E-05   4.848447E-06  -9.762305E-05  -3.700933E-06  -2.859244E-07
  b3u  15  -1.372644E-04   7.746525E-05  -2.101665E-05   1.537080E-05  -5.067472E-04   7.227666E-06   8.972279E-06   2.024691E-04
  b3u  16  -1.438737E-03   2.186137E-04   1.704380E-06   1.182811E-04   4.122039E-04  -3.365851E-06  -9.075957E-05  -3.348894E-04
  b3u  17   2.549502E-04   6.437206E-05  -2.247135E-05  -1.027654E-03  -5.037587E-04   4.788672E-06   2.525637E-05  -5.726083E-05
  b3u  18  -1.316167E-05  -1.516558E-05  -1.061474E-03   1.048633E-05   1.004118E-05  -3.620100E-04  -1.878174E-06  -1.812439E-06
  b3u  19   5.760441E-04   2.779147E-04   5.826910E-06   1.149930E-03   2.651300E-04  -8.339535E-06  -4.720377E-05   1.544504E-04
  b3u  20   8.581414E-06   2.019476E-06  -7.511507E-04   7.024591E-06   7.434211E-06   3.508365E-04   3.598884E-06   3.208884E-06
  b3u  21  -5.094180E-04  -6.220662E-04   1.420904E-05   1.090239E-03  -1.300510E-03   1.339853E-05  -5.282798E-05  -2.000635E-04
  b3u  22  -3.878122E-06   3.371230E-06   1.927976E-03  -3.850449E-05   1.736004E-05   2.284094E-04   4.889740E-07  -1.536510E-06
  b3u  23   7.242442E-07  -5.405697E-06   4.399689E-04  -7.000010E-06   8.809300E-06   2.967738E-04   3.705134E-06   3.620473E-06
  b3u  24   7.268863E-04  -2.225291E-04  -6.280922E-06  -9.161314E-05   1.075532E-04  -5.153173E-06   4.867165E-04   8.965045E-04
  b3u  25   1.363807E-03   6.086590E-04  -1.482305E-05  -3.957472E-04   1.365705E-04   3.370146E-06  -4.569258E-04   4.636605E-04
  b3u  26   6.086590E-04   1.276579E-03  -7.176494E-06  -3.093357E-04   3.513453E-04   4.166959E-06  -6.022421E-04   2.818450E-04
  b3u  27  -1.482305E-05  -7.176494E-06   1.672272E-03  -2.723045E-06   8.527380E-06   4.150476E-04   1.286679E-05  -9.389440E-06
  b3u  28  -3.957472E-04  -3.093357E-04  -2.723045E-06   1.473270E-03   1.911788E-04  -1.383864E-05   3.266264E-04  -2.823122E-04
  b3u  29   1.365705E-04   3.513453E-04   8.527380E-06   1.911788E-04   1.259427E-03  -5.588031E-06  -5.416550E-05   1.278523E-04
  b3u  30   3.370146E-06   4.166959E-06   4.150476E-04  -1.383864E-05  -5.588031E-06   6.469127E-04  -3.424830E-06   9.306975E-07
  b3u  31  -4.569258E-04  -6.022421E-04   1.286679E-05   3.266264E-04  -5.416550E-05  -3.424830E-06   1.211700E-03  -1.639490E-04
  b3u  32   4.636605E-04   2.818450E-04  -9.389440E-06  -2.823122E-04   1.278523E-04   9.306975E-07  -1.639490E-04   1.118823E-03
  b3u  33   1.688418E-06   1.858530E-06  -5.917499E-05   1.854785E-06  -8.372522E-07  -2.381812E-04  -9.107592E-07   1.335377E-06
  b3u  34   4.460919E-04   9.127236E-05  -1.027446E-06   9.592681E-05  -6.230686E-05   3.574643E-06  -3.850928E-04   2.083301E-04
  b3u  35   3.110537E-05  -2.883212E-05  -3.801972E-04  -7.206021E-05   5.825473E-05  -1.535854E-04   8.748913E-05   1.853099E-05
  b3u  36  -7.749693E-05   6.215405E-05  -1.476043E-04   2.057744E-04  -1.540243E-04  -5.791627E-05  -2.253703E-04  -4.842997E-05
  b3u  37  -6.395765E-05  -3.909929E-05   2.461581E-06   2.245796E-05   1.248800E-04  -5.422243E-07  -1.004858E-04  -3.692574E-04
  b3u  38  -9.849037E-07   7.080169E-07   2.379913E-04  -4.440356E-06   5.949967E-07  -1.136541E-06   8.278290E-07  -3.289848E-07
  b3u  39   4.436265E-07   6.613617E-07  -7.377637E-05   6.698142E-07  -9.513630E-07  -5.605492E-05  -4.281399E-07   8.553254E-07

               b3u  33        b3u  34        b3u  35        b3u  36        b3u  37        b3u  38        b3u  39
  b3u   3   1.873997E-04  -1.807575E-03   2.603154E-04  -8.215743E-04  -2.043536E-03   3.549565E-05   2.660460E-05
  b3u   4  -2.404308E-03  -1.200168E-04  -3.041564E-03  -1.147934E-03  -6.647869E-05   1.111245E-07   1.056629E-03
  b3u   5  -1.770965E-04  -9.667974E-04   4.670420E-04  -1.419615E-03   2.593740E-03   2.503925E-05   1.854419E-05
  b3u   6  -4.030058E-08  -1.460900E-04  -6.170446E-07  -1.278852E-06  -4.845054E-05   2.497427E-07   4.195962E-07
  b3u   7   9.831230E-05   3.164994E-07  -1.169225E-05  -4.659070E-06   8.557658E-08   6.047514E-05   1.097665E-05
  b3u   8   7.859786E-07   1.970640E-04   9.954938E-06  -2.191931E-05   6.157815E-05   3.637695E-07  -4.608264E-07
  b3u   9   2.820353E-06   2.986885E-04  -2.826726E-05   7.459971E-05   7.552988E-05   9.571098E-07  -4.403723E-07
  b3u  10  -3.342354E-04   2.454853E-06   8.351704E-05   3.331491E-05   8.504603E-07  -2.052039E-04  -3.160563E-05
  b3u  11   1.407955E-06  -5.222919E-05  -2.657456E-05   6.547398E-05   3.813686E-05  -4.843116E-08   6.400202E-08
  b3u  12   5.163864E-04   1.944923E-05  -1.972718E-04  -6.465893E-05  -2.627663E-06   2.956760E-04   3.910149E-05
  b3u  13   1.326114E-05  -5.416866E-04   1.285000E-04  -3.433554E-04   1.047232E-04   8.344638E-06   1.409857E-06
  b3u  14  -5.063310E-04  -5.880563E-07   3.273312E-04   1.332220E-04  -3.803087E-06  -3.402965E-04  -5.526535E-05
  b3u  15  -8.542531E-06   4.351737E-04   3.411260E-05  -6.229969E-05   2.909510E-04  -6.415089E-06  -2.294901E-06
  b3u  16   2.274990E-06   4.591753E-05   3.325173E-05  -8.452454E-05   4.418524E-04  -6.582789E-08  -5.230964E-07
  b3u  17  -5.436846E-06   4.625567E-05   8.029721E-05  -1.878954E-04   1.094727E-04  -1.812671E-06  -1.344912E-06
  b3u  18  -7.894060E-04  -6.077453E-06   7.567264E-04   2.949092E-04   4.640282E-06  -3.468202E-04  -1.774830E-04
  b3u  19  -7.123812E-06   2.587762E-04  -5.373672E-05   1.658427E-04  -8.538886E-05  -4.110701E-06  -2.707301E-06
  b3u  20  -2.946190E-04  -2.022103E-06  -3.499863E-06  -2.810661E-06  -3.153107E-06  -1.548598E-04   4.329790E-05
  b3u  21  -8.973486E-06   7.349246E-04  -1.030601E-04   2.890222E-04   3.116549E-04  -3.766411E-06  -3.048545E-06
  b3u  22  -3.044402E-04  -4.612062E-06   1.765466E-04   6.582929E-05   4.420759E-06   2.820680E-04  -2.973415E-04
  b3u  23  -7.437362E-04  -4.988754E-07   3.279841E-04   1.276625E-04   1.848372E-06  -3.944223E-06  -1.443576E-04
  b3u  24   8.332333E-06  -1.607401E-04   2.484030E-05  -7.843231E-05  -4.754681E-04   4.198518E-07   2.531054E-06
  b3u  25   1.688418E-06   4.460919E-04   3.110537E-05  -7.749693E-05  -6.395765E-05  -9.849037E-07   4.436265E-07
  b3u  26   1.858530E-06   9.127236E-05  -2.883212E-05   6.215405E-05  -3.909929E-05   7.080169E-07   6.613617E-07
  b3u  27  -5.917499E-05  -1.027446E-06  -3.801972E-04  -1.476043E-04   2.461581E-06   2.379913E-04  -7.377637E-05
  b3u  28   1.854785E-06   9.592681E-05  -7.206021E-05   2.057744E-04   2.245796E-05  -4.440356E-06   6.698142E-07
  b3u  29  -8.372522E-07  -6.230686E-05   5.825473E-05  -1.540243E-04   1.248800E-04   5.949967E-07  -9.513630E-07
  b3u  30  -2.381812E-04   3.574643E-06  -1.535854E-04  -5.791627E-05  -5.422243E-07  -1.136541E-06  -5.605492E-05
  b3u  31  -9.107592E-07  -3.850928E-04   8.748913E-05  -2.253703E-04  -1.004858E-04   8.278290E-07  -4.281399E-07
  b3u  32   1.335377E-06   2.083301E-04   1.853099E-05  -4.842997E-05  -3.692574E-04  -3.289848E-07   8.553254E-07
  b3u  33   8.357011E-04   8.282484E-07  -1.356813E-04  -5.382528E-05  -1.185401E-06   3.902849E-05  -9.854314E-05
  b3u  34   8.282484E-07   5.853160E-04  -3.983917E-05   1.040129E-04   2.258116E-05  -1.371627E-06  -3.728470E-07
  b3u  35  -1.356813E-04  -3.983917E-05   4.588928E-04   4.386047E-05  -4.513202E-06   7.493300E-05  -2.704361E-05
  b3u  36  -5.382528E-05   1.040129E-04   4.386047E-05   3.669546E-04   1.353845E-05   2.922871E-05  -1.065041E-05
  b3u  37  -1.185401E-06   2.258116E-05  -4.513202E-06   1.353845E-05   4.483137E-04  -3.488886E-07  -6.324546E-07
  b3u  38   3.902849E-05  -1.371627E-06   7.493300E-05   2.922871E-05  -3.488886E-07   2.535807E-04  -7.885650E-05
  b3u  39  -9.854314E-05  -3.728470E-07  -2.704361E-05  -1.065041E-05  -6.324546E-07  -7.885650E-05   1.095126E-04

Natural orbital populations,block 2
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98490684     1.97882182     1.97669075     0.01406890     0.01290694     0.01084934
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00851044     0.00470531     0.00331181     0.00211261     0.00200117     0.00125339     0.00080572     0.00075418
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00044342     0.00044180     0.00032355     0.00021807     0.00017112     0.00013091     0.00009768     0.00006196
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00004065     0.00003130     0.00002936     0.00001832     0.00001178     0.00000672     0.00000653     0.00000391
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000381     0.00000141     0.00000058     0.00000018     0.00000012     0.00000006     0.00000000

          modens reordered block   1

               b2u   1        b2u   2        b2u   3        b2u   4        b2u   5        b2u   6        b2u   7        b2u   8
  b2u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b2u   2    0.00000        1.98450       8.434328E-05  -1.732922E-03   5.698120E-04   1.468877E-04   1.386972E-03   5.643425E-05
  b2u   3    0.00000       8.434328E-05    1.97470      -4.118111E-05   4.225650E-05  -5.218430E-03   1.071384E-04   2.477055E-05
  b2u   4    0.00000      -1.732922E-03  -4.118111E-05    1.97703      -2.952882E-03   3.416621E-06  -3.043426E-03  -3.548948E-03
  b2u   5    0.00000       5.698120E-04   4.225650E-05  -2.952882E-03   3.880913E-04   8.166873E-06   4.988952E-04   6.005054E-04
  b2u   6    0.00000       1.468877E-04  -5.218430E-03   3.416621E-06   8.166873E-06   4.604341E-04   4.946097E-06   1.092516E-05
  b2u   7    0.00000       1.386972E-03   1.071384E-04  -3.043426E-03   4.988952E-04   4.946097E-06   7.075518E-04   6.769787E-04
  b2u   8    0.00000       5.643425E-05   2.477055E-05  -3.548948E-03   6.005054E-04   1.092516E-05   6.769787E-04   1.141376E-03
  b2u   9    0.00000       2.257642E-03  -3.316984E-04  -9.305515E-05   2.483446E-04   9.978467E-06   5.000614E-04   7.106358E-05
  b2u  10    0.00000       3.248098E-04   6.949005E-04   5.001289E-04  -9.390803E-04  -7.589945E-06  -9.291285E-04  -2.105835E-03
  b2u  11    0.00000      -2.263782E-04   8.501252E-03  -1.512953E-04  -1.158664E-05  -9.632508E-04   1.502373E-06  -2.427493E-05
  b2u  12    0.00000      -8.615687E-04  -9.219085E-04   4.589668E-03  -5.686894E-04  -5.564605E-06  -7.906359E-04  -8.925665E-04
  b2u  13    0.00000       2.487871E-03  -1.723598E-03   3.985382E-03   2.401183E-04   1.241729E-05   3.957637E-04   2.560779E-04
  b2u  14    0.00000      -1.637189E-03  -1.246728E-04   1.529992E-03  -1.886295E-04  -9.942455E-06  -4.693053E-04   1.466664E-04
  b2u  15    0.00000       1.013442E-03   4.277377E-04  -1.740705E-03  -3.328522E-04  -4.085743E-06  -2.053205E-04  -9.870991E-04
  b2u  16    0.00000      -1.922632E-04   6.497042E-03  -3.114686E-04   8.104448E-07  -9.783450E-04   1.607383E-05  -3.997435E-07
  b2u  17    0.00000      -9.922889E-05  -8.072164E-04   1.127016E-03  -7.065983E-04  -1.184438E-05  -9.204524E-04  -1.441924E-03
  b2u  18    0.00000       1.475724E-03  -1.000234E-03   2.223197E-03  -1.687181E-05   2.900385E-06   3.930801E-05  -8.595492E-05
  b2u  19    0.00000      -2.575235E-04  -3.641139E-04  -4.677083E-04   3.461503E-04   1.017975E-05   4.713718E-04   4.803141E-04
  b2u  20    0.00000       5.872508E-05   8.596059E-04   1.638831E-04  -5.885013E-06   4.047654E-04  -1.223870E-05  -8.690011E-06
  b2u  21    0.00000      -1.381906E-04   1.789815E-05   3.974340E-04   1.753496E-04   1.164872E-05   1.430539E-04   2.924375E-04
  b2u  22    0.00000      -9.945677E-04   4.223772E-05   2.439188E-03  -5.510078E-06  -3.129663E-06  -1.712882E-04   4.043665E-04
  b2u  23    0.00000       8.594263E-05  -7.027993E-03   1.962098E-04   6.054921E-07   5.400898E-04  -6.672036E-06   8.458859E-07
  b2u  24    0.00000       2.335582E-03  -4.062817E-05   1.621486E-03   2.508596E-04   7.402431E-06   4.581287E-04   2.822168E-04
  b2u  25    0.00000      -4.231587E-04   1.991472E-04  -4.455110E-03   9.422400E-05   1.235794E-06   6.055746E-05   2.160469E-04
  b2u  26    0.00000       8.378761E-04  -1.285926E-04  -4.018398E-04   6.626620E-05   5.453192E-07   9.798910E-05   1.130335E-04
  b2u  27    0.00000       6.831947E-05  -3.611879E-03   7.044313E-05   1.377281E-06   2.396914E-05   1.982200E-06   2.276614E-06
  b2u  28    0.00000       1.769055E-03  -1.218544E-05  -8.369287E-04  -1.465821E-04  -2.740884E-06  -1.979146E-04  -2.984768E-04
  b2u  29    0.00000      -8.773882E-04   2.171758E-05   1.363608E-03  -3.270847E-07  -2.179013E-07  -2.535568E-05   7.796779E-05
  b2u  30    0.00000       1.903346E-03  -2.095191E-04   2.569294E-03  -4.837205E-05  -1.674659E-07  -6.088774E-05  -7.553100E-05

               b2u   9        b2u  10        b2u  11        b2u  12        b2u  13        b2u  14        b2u  15        b2u  16
  b2u   2   2.257642E-03   3.248098E-04  -2.263782E-04  -8.615687E-04   2.487871E-03  -1.637189E-03   1.013442E-03  -1.922632E-04
  b2u   3  -3.316984E-04   6.949005E-04   8.501252E-03  -9.219085E-04  -1.723598E-03  -1.246728E-04   4.277377E-04   6.497042E-03
  b2u   4  -9.305515E-05   5.001289E-04  -1.512953E-04   4.589668E-03   3.985382E-03   1.529992E-03  -1.740705E-03  -3.114686E-04
  b2u   5   2.483446E-04  -9.390803E-04  -1.158664E-05  -5.686894E-04   2.401183E-04  -1.886295E-04  -3.328522E-04   8.104448E-07
  b2u   6   9.978467E-06  -7.589945E-06  -9.632508E-04  -5.564605E-06   1.241729E-05  -9.942455E-06  -4.085743E-06  -9.783450E-04
  b2u   7   5.000614E-04  -9.291285E-04   1.502373E-06  -7.906359E-04   3.957637E-04  -4.693053E-04  -2.053205E-04   1.607383E-05
  b2u   8   7.106358E-05  -2.105835E-03  -2.427493E-05  -8.925665E-04   2.560779E-04   1.466664E-04  -9.870991E-04  -3.997435E-07
  b2u   9   9.426179E-04  -3.147609E-04  -4.907547E-06   2.449010E-04   1.534619E-03  -6.898374E-04   2.458898E-04  -2.933811E-07
  b2u  10  -3.147609E-04   6.148213E-03   4.081567E-05  -4.849233E-04  -3.457298E-03  -1.367279E-03   3.100257E-03  -1.232563E-05
  b2u  11  -4.907547E-06   4.081567E-05   2.110816E-03  -5.872634E-06  -3.076774E-05  -5.279854E-06   3.042150E-05   2.199052E-03
  b2u  12   2.449010E-04  -4.849233E-04  -5.872634E-06   3.124791E-03   3.006106E-03   7.980168E-04  -4.010149E-04  -1.836153E-05
  b2u  13   1.534619E-03  -3.457298E-03  -3.076774E-05   3.006106E-03   6.236215E-03   4.488526E-04  -1.482631E-03  -1.236413E-05
  b2u  14  -6.898374E-04  -1.367279E-03  -5.279854E-06   7.980168E-04   4.488526E-04   1.274401E-03  -1.114295E-03  -5.846233E-08
  b2u  15   2.458898E-04   3.100257E-03   3.042150E-05  -4.010149E-04  -1.482631E-03  -1.114295E-03   1.938881E-03   6.546857E-06
  b2u  16  -2.933811E-07  -1.232563E-05   2.199052E-03  -1.836153E-05  -1.236413E-05  -5.846233E-08   6.546857E-06   2.419893E-03
  b2u  17   2.218517E-04   1.654834E-03   2.401669E-05   3.005136E-03   2.148787E-03   1.222501E-04   4.987311E-04  -6.139766E-06
  b2u  18   5.812959E-04  -1.184384E-03  -8.968118E-06   1.425235E-03   3.109297E-03   3.890618E-04  -7.297371E-04  -5.485026E-06
  b2u  19   5.241124E-04  -1.197490E-03  -1.538337E-05   1.536707E-04   1.453261E-03  -2.474985E-04  -5.791932E-04  -4.716828E-06
  b2u  20  -2.992369E-06   2.486341E-05  -8.740634E-04   4.787300E-07  -1.153982E-05   3.976043E-07   1.082785E-05  -1.078415E-03
  b2u  21  -9.489954E-05  -2.983477E-04  -2.339966E-05  -7.354676E-05  -2.227625E-04  -6.371673E-05  -2.730505E-04  -2.474996E-05
  b2u  22  -5.980796E-04  -1.460194E-03  -1.421820E-05  -5.291074E-06  -1.234247E-04   1.032333E-03  -1.142063E-03  -4.028622E-06
  b2u  23   4.415902E-06   5.443327E-06  -1.426775E-03   1.278434E-05   1.550218E-05  -2.997969E-06  -4.370450E-06  -1.687641E-03
  b2u  24   5.467959E-04  -4.392794E-04  -5.615397E-06  -5.083057E-04   4.063578E-04  -5.015239E-04   2.589489E-04   3.573345E-06
  b2u  25  -3.100393E-05  -4.700211E-04  -5.495752E-06   1.461217E-05  -8.432546E-05   2.690409E-05  -5.071519E-05  -1.518546E-06
  b2u  26   2.415134E-05  -1.689910E-04   1.176372E-06  -2.008312E-04   3.369803E-04   5.772694E-05  -1.527757E-04   8.004863E-06
  b2u  27   2.264129E-06  -5.698490E-06  -9.552091E-05  -7.930099E-07   1.562659E-05   1.244572E-06  -5.342583E-06  -2.895489E-04
  b2u  28  -5.187824E-05   5.459711E-04   6.249342E-06   4.272770E-04   4.283157E-05   4.297026E-05   2.608429E-04  -1.023013E-06
  b2u  29  -7.112447E-05  -3.628834E-04  -3.850378E-06   7.732795E-05   9.414922E-05   2.061230E-04  -1.719446E-04  -1.639645E-06
  b2u  30   3.868554E-05  -9.867157E-05  -1.136652E-06   2.864971E-04   4.380694E-04   1.071091E-04  -8.179616E-05  -1.707890E-06

               b2u  17        b2u  18        b2u  19        b2u  20        b2u  21        b2u  22        b2u  23        b2u  24
  b2u   2  -9.922889E-05   1.475724E-03  -2.575235E-04   5.872508E-05  -1.381906E-04  -9.945677E-04   8.594263E-05   2.335582E-03
  b2u   3  -8.072164E-04  -1.000234E-03  -3.641139E-04   8.596059E-04   1.789815E-05   4.223772E-05  -7.027993E-03  -4.062817E-05
  b2u   4   1.127016E-03   2.223197E-03  -4.677083E-04   1.638831E-04   3.974340E-04   2.439188E-03   1.962098E-04   1.621486E-03
  b2u   5  -7.065983E-04  -1.687181E-05   3.461503E-04  -5.885013E-06   1.753496E-04  -5.510078E-06   6.054921E-07   2.508596E-04
  b2u   6  -1.184438E-05   2.900385E-06   1.017975E-05   4.047654E-04   1.164872E-05  -3.129663E-06   5.400898E-04   7.402431E-06
  b2u   7  -9.204524E-04   3.930801E-05   4.713718E-04  -1.223870E-05   1.430539E-04  -1.712882E-04  -6.672036E-06   4.581287E-04
  b2u   8  -1.441924E-03  -8.595492E-05   4.803141E-04  -8.690011E-06   2.924375E-04   4.043665E-04   8.458859E-07   2.822168E-04
  b2u   9   2.218517E-04   5.812959E-04   5.241124E-04  -2.992369E-06  -9.489954E-05  -5.980796E-04   4.415902E-06   5.467959E-04
  b2u  10   1.654834E-03  -1.184384E-03  -1.197490E-03   2.486341E-05  -2.983477E-04  -1.460194E-03   5.443327E-06  -4.392794E-04
  b2u  11   2.401669E-05  -8.968118E-06  -1.538337E-05  -8.740634E-04  -2.339966E-05  -1.421820E-05  -1.426775E-03  -5.615397E-06
  b2u  12   3.005136E-03   1.425235E-03   1.536707E-04   4.787300E-07  -7.354676E-05  -5.291074E-06   1.278434E-05  -5.083057E-04
  b2u  13   2.148787E-03   3.109297E-03   1.453261E-03  -1.153982E-05  -2.227625E-04  -1.234247E-04   1.550218E-05   4.063578E-04
  b2u  14   1.222501E-04   3.890618E-04  -2.474985E-04   3.976043E-07  -6.371673E-05   1.032333E-03  -2.997969E-06  -5.015239E-04
  b2u  15   4.987311E-04  -7.297371E-04  -5.791932E-04   1.082785E-05  -2.730505E-04  -1.142063E-03  -4.370450E-06   2.589489E-04
  b2u  16  -6.139766E-06  -5.485026E-06  -4.716828E-06  -1.078415E-03  -2.474996E-05  -4.028622E-06  -1.687641E-03   3.573345E-06
  b2u  17   4.938588E-03   1.572485E-03   5.345477E-04  -2.038955E-05   6.248218E-04  -1.082231E-03   4.275795E-06  -1.306364E-03
  b2u  18   1.572485E-03   2.497181E-03   7.374380E-04  -4.938961E-06  -2.230383E-04  -9.532729E-05   7.058060E-06  -1.156584E-04
  b2u  19   5.345477E-04   7.374380E-04   1.373361E-03  -1.873955E-05   6.091927E-04  -4.006902E-04   3.217084E-06  -1.437989E-04
  b2u  20  -2.038955E-05  -4.938961E-06  -1.873955E-05   9.542968E-04  -9.087889E-06   1.060646E-05   4.907860E-04   7.817814E-06
  b2u  21   6.248218E-04  -2.230383E-04   6.091927E-04  -9.087889E-06   1.273617E-03  -3.106073E-04   1.057700E-05  -3.535995E-04
  b2u  22  -1.082231E-03  -9.532729E-05  -4.006902E-04   1.060646E-05  -3.106073E-04   1.472428E-03  -4.177709E-07  -1.894851E-04
  b2u  23   4.275795E-06   7.058060E-06   3.217084E-06   4.907860E-04   1.057700E-05  -4.177709E-07   2.054145E-03  -1.282593E-06
  b2u  24  -1.306364E-03  -1.156584E-04  -1.437989E-04   7.817814E-06  -3.535995E-04  -1.894851E-04  -1.282593E-06   1.258228E-03
  b2u  25  -4.480579E-05  -4.809527E-04   4.623006E-04  -1.414568E-05   6.053553E-04  -3.290788E-04   1.307667E-06  -5.705943E-05
  b2u  26   2.041929E-04   9.008586E-04   4.656491E-04  -1.539652E-05   2.789422E-04  -2.830997E-04  -1.299577E-05  -1.301934E-04
  b2u  27   8.656045E-06   2.606703E-05   1.340114E-05   2.202233E-04   1.100645E-05  -6.772723E-06   5.136046E-04  -4.024411E-06
  b2u  28   7.302303E-04   1.630215E-04  -4.436197E-04   4.464265E-06  -9.083961E-05  -9.507788E-05   2.706184E-07  -6.268930E-05
  b2u  29  -2.979013E-04  -7.821562E-05  -8.600394E-05  -3.949363E-07   6.779208E-05   2.162581E-04   8.796389E-07   1.657226E-04
  b2u  30   3.123531E-04   4.756727E-04   6.808449E-05  -1.783422E-06   3.775100E-05  -2.492896E-05   1.330601E-06   1.218016E-04

               b2u  25        b2u  26        b2u  27        b2u  28        b2u  29        b2u  30
  b2u   2  -4.231587E-04   8.378761E-04   6.831947E-05   1.769055E-03  -8.773882E-04   1.903346E-03
  b2u   3   1.991472E-04  -1.285926E-04  -3.611879E-03  -1.218544E-05   2.171758E-05  -2.095191E-04
  b2u   4  -4.455110E-03  -4.018398E-04   7.044313E-05  -8.369287E-04   1.363608E-03   2.569294E-03
  b2u   5   9.422400E-05   6.626620E-05   1.377281E-06  -1.465821E-04  -3.270847E-07  -4.837205E-05
  b2u   6   1.235794E-06   5.453192E-07   2.396914E-05  -2.740884E-06  -2.179013E-07  -1.674659E-07
  b2u   7   6.055746E-05   9.798910E-05   1.982200E-06  -1.979146E-04  -2.535568E-05  -6.088774E-05
  b2u   8   2.160469E-04   1.130335E-04   2.276614E-06  -2.984768E-04   7.796779E-05  -7.553100E-05
  b2u   9  -3.100393E-05   2.415134E-05   2.264129E-06  -5.187824E-05  -7.112447E-05   3.868554E-05
  b2u  10  -4.700211E-04  -1.689910E-04  -5.698490E-06   5.459711E-04  -3.628834E-04  -9.867157E-05
  b2u  11  -5.495752E-06   1.176372E-06  -9.552091E-05   6.249342E-06  -3.850378E-06  -1.136652E-06
  b2u  12   1.461217E-05  -2.008312E-04  -7.930099E-07   4.272770E-04   7.732795E-05   2.864971E-04
  b2u  13  -8.432546E-05   3.369803E-04   1.562659E-05   4.283157E-05   9.414922E-05   4.380694E-04
  b2u  14   2.690409E-05   5.772694E-05   1.244572E-06   4.297026E-05   2.061230E-04   1.071091E-04
  b2u  15  -5.071519E-05  -1.527757E-04  -5.342583E-06   2.608429E-04  -1.719446E-04  -8.179616E-05
  b2u  16  -1.518546E-06   8.004863E-06  -2.895489E-04  -1.023013E-06  -1.639645E-06  -1.707890E-06
  b2u  17  -4.480579E-05   2.041929E-04   8.656045E-06   7.302303E-04  -2.979013E-04   3.123531E-04
  b2u  18  -4.809527E-04   9.008586E-04   2.606703E-05   1.630215E-04  -7.821562E-05   4.756727E-04
  b2u  19   4.623006E-04   4.656491E-04   1.340114E-05  -4.436197E-04  -8.600394E-05   6.808449E-05
  b2u  20  -1.414568E-05  -1.539652E-05   2.202233E-04   4.464265E-06  -3.949363E-07  -1.783422E-06
  b2u  21   6.053553E-04   2.789422E-04   1.100645E-05  -9.083961E-05   6.779208E-05   3.775100E-05
  b2u  22  -3.290788E-04  -2.830997E-04  -6.772723E-06  -9.507788E-05   2.162581E-04  -2.492896E-05
  b2u  23   1.307667E-06  -1.299577E-05   5.136046E-04   2.706184E-07   8.796389E-07   1.330601E-06
  b2u  24  -5.705943E-05  -1.301934E-04  -4.024411E-06  -6.268930E-05   1.657226E-04   1.218016E-04
  b2u  25   1.215172E-03   1.647450E-04   3.632241E-06  -3.876986E-04   2.387501E-04  -9.939498E-05
  b2u  26   1.647450E-04   1.118965E-03   1.172481E-05  -2.051547E-04  -5.154039E-05   3.707299E-04
  b2u  27   3.632241E-06   1.172481E-05   6.489917E-04  -4.461670E-06  -1.064371E-06   9.079263E-06
  b2u  28  -3.876986E-04  -2.051547E-04  -4.461670E-06   5.842000E-04  -1.095069E-04   2.222848E-05
  b2u  29   2.387501E-04  -5.154039E-05  -1.064371E-06  -1.095069E-04   3.476932E-04  -1.364818E-05
  b2u  30  -9.939498E-05   3.707299E-04   9.079263E-06   2.222848E-05  -1.364818E-05   4.470607E-04

Natural orbital populations,block 3
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98490368     1.97670402     1.97480706     0.01408601     0.01083768     0.00651437     0.00469515
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00331280     0.00211335     0.00125196     0.00111467     0.00080258     0.00068218     0.00044245     0.00032346
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00020834     0.00017084     0.00013068     0.00009750     0.00003124     0.00002940     0.00002033     0.00001173
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00000672     0.00000383     0.00000345     0.00000141     0.00000012     0.00000006

          modens reordered block   1

               b1g   1        b1g   2        b1g   3        b1g   4        b1g   5        b1g   6        b1g   7        b1g   8
  b1g   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b1g   2    0.00000        1.98080       2.119815E-03  -1.181509E-03  -2.821300E-03  -2.539990E-04  -8.875266E-06   3.207100E-03
  b1g   3    0.00000       2.119815E-03    1.97448      -1.903515E-03  -1.273778E-03  -3.990924E-03   1.891201E-05  -1.519768E-03
  b1g   4    0.00000      -1.181509E-03  -1.903515E-03   2.221039E-04   3.997304E-04   3.635808E-04   2.184771E-07  -1.811624E-04
  b1g   5    0.00000      -2.821300E-03  -1.273778E-03   3.997304E-04   8.568553E-04   5.344458E-04  -1.374685E-07  -5.553967E-04
  b1g   6    0.00000      -2.539990E-04  -3.990924E-03   3.635808E-04   5.344458E-04   7.565821E-04   9.116071E-07  -1.874111E-05
  b1g   7    0.00000      -8.875266E-06   1.891201E-05   2.184771E-07  -1.374685E-07   9.116071E-07   9.283387E-06   8.895691E-07
  b1g   8    0.00000       3.207100E-03  -1.519768E-03  -1.811624E-04  -5.553967E-04  -1.874111E-05   8.895691E-07   7.295550E-04
  b1g   9    0.00000       2.758619E-03  -9.191275E-04  -6.033415E-04  -1.586626E-03  -6.568529E-04   1.414714E-06   1.078329E-03
  b1g  10    0.00000       1.273117E-03  -6.582873E-03   5.873292E-04   6.222633E-04   1.553457E-03   2.837305E-06   5.363891E-04
  b1g  11    0.00000      -3.884914E-04  -2.697020E-03  -3.631464E-04  -1.323559E-03  -2.168401E-04   4.223643E-06   9.828446E-04
  b1g  12    0.00000      -1.651597E-04   4.647877E-04   5.714185E-06   2.183140E-05   2.295283E-06   8.925593E-05  -1.678464E-05
  b1g  13    0.00000      -3.327854E-03   1.346950E-03   1.559231E-04   5.115385E-04  -1.100909E-04  -1.202477E-06  -9.478437E-04
  b1g  14    0.00000       2.012043E-03  -1.667695E-03   3.854773E-04   4.591478E-04   1.062568E-03   1.591089E-06   4.648784E-04
  b1g  15    0.00000      -5.416290E-04   1.225808E-03  -2.069298E-06  -1.202299E-06  -6.846958E-06   1.898438E-04  -1.884608E-06
  b1g  16    0.00000      -3.888817E-04  -1.630706E-03   5.003957E-04   8.581234E-04   1.054133E-03   2.358893E-06  -1.586832E-04
  b1g  17    0.00000       1.488742E-03   5.932277E-04   9.655491E-05   4.123059E-04   9.240129E-05  -5.140016E-08  -3.117090E-04
  b1g  18    0.00000       4.220022E-04  -7.768676E-04  -3.050424E-07  -3.111716E-06   9.863219E-07  -6.272113E-05   2.985059E-06
  b1g  19    0.00000       1.137051E-03  -3.930285E-04   1.265972E-04   1.895982E-04   2.855044E-04   2.830302E-07  -4.479539E-05
  b1g  20    0.00000       2.223193E-03   1.397608E-03   8.966365E-07  -8.861980E-05   2.757162E-04   9.384793E-07   5.813456E-04
  b1g  21    0.00000       3.232066E-04  -2.249259E-03  -1.218543E-04  -1.228919E-04  -2.336883E-04  -4.042728E-07   4.574915E-05
  b1g  22    0.00000       1.923228E-04  -2.213217E-04   1.841221E-07   4.439772E-07   4.013725E-07  -6.227633E-07  -8.784473E-07
  b1g  23    0.00000       2.823003E-03  -3.275431E-03  -1.471564E-04  -4.288520E-04  -8.140049E-05   3.461213E-07   4.224393E-04
  b1g  24    0.00000      -1.337837E-03  -3.165167E-04   1.493990E-04   3.017239E-04   2.953586E-04   3.662643E-07  -2.375969E-04
  b1g  25    0.00000      -1.057886E-04   6.747700E-05   1.136812E-07   5.861275E-08   8.050153E-07  -3.965976E-05   2.839535E-07
  b1g  26    0.00000       2.237160E-03   1.385623E-03  -2.460279E-06   3.364156E-05  -1.165474E-05   1.304458E-08  -4.494861E-05
  b1g  27    0.00000       9.319489E-04   3.730759E-03   6.455613E-06   8.201954E-05  -1.350687E-05  -1.578650E-07  -9.004725E-05
  b1g  28    0.00000       2.638257E-03  -2.858927E-03   1.015990E-04   1.783048E-04   2.140366E-04   1.068813E-07  -2.592333E-05
  b1g  29    0.00000      -9.868272E-05   3.207798E-05   1.374037E-07   2.669682E-07   2.939152E-07  -5.297433E-06   9.120385E-09
  b1g  30    0.00000      -5.020570E-04   5.371872E-04  -1.058916E-05  -4.641577E-05   4.239689E-05   1.888170E-07   8.480075E-05

               b1g   9        b1g  10        b1g  11        b1g  12        b1g  13        b1g  14        b1g  15        b1g  16
  b1g   2   2.758619E-03   1.273117E-03  -3.884914E-04  -1.651597E-04  -3.327854E-03   2.012043E-03  -5.416290E-04  -3.888817E-04
  b1g   3  -9.191275E-04  -6.582873E-03  -2.697020E-03   4.647877E-04   1.346950E-03  -1.667695E-03   1.225808E-03  -1.630706E-03
  b1g   4  -6.033415E-04   5.873292E-04  -3.631464E-04   5.714185E-06   1.559231E-04   3.854773E-04  -2.069298E-06   5.003957E-04
  b1g   5  -1.586626E-03   6.222633E-04  -1.323559E-03   2.183140E-05   5.115385E-04   4.591478E-04  -1.202299E-06   8.581234E-04
  b1g   6  -6.568529E-04   1.553457E-03  -2.168401E-04   2.295283E-06  -1.100909E-04   1.062568E-03  -6.846958E-06   1.054133E-03
  b1g   7   1.414714E-06   2.837305E-06   4.223643E-06   8.925593E-05  -1.202477E-06   1.591089E-06   1.898438E-04   2.358893E-06
  b1g   8   1.078329E-03   5.363891E-04   9.828446E-04  -1.678464E-05  -9.478437E-04   4.648784E-04  -1.884608E-06  -1.586832E-04
  b1g   9   4.031835E-03  -1.945281E-04   4.778925E-03  -7.701984E-05  -7.629031E-04  -5.034113E-04   2.865294E-06  -1.144110E-03
  b1g  10  -1.945281E-04   3.952575E-03   1.060175E-03  -1.953809E-05  -1.007167E-03   2.829012E-03  -1.463857E-05   2.566040E-03
  b1g  11   4.778925E-03   1.060175E-03   7.509975E-03  -1.005872E-04  -5.656304E-04   1.862955E-04   5.263901E-05  -8.784696E-05
  b1g  12  -7.701984E-05  -1.953809E-05  -1.005872E-04   9.670242E-04   9.204901E-06  -7.930554E-06   2.311380E-03   1.413660E-05
  b1g  13  -7.629031E-04  -1.007167E-03  -5.656304E-04   9.204901E-06   1.519572E-03  -1.050028E-03  -4.463019E-07  -2.210556E-04
  b1g  14  -5.034113E-04   2.829012E-03   1.862955E-04  -7.930554E-06  -1.050028E-03   2.308185E-03  -1.865563E-05   2.018664E-03
  b1g  15   2.865294E-06  -1.463857E-05   5.263901E-05   2.311380E-03  -4.463019E-07  -1.865563E-05   6.217326E-03   2.677083E-05
  b1g  16  -1.144110E-03   2.566040E-03  -8.784696E-05   1.413660E-05  -2.210556E-04   2.018664E-03   2.677083E-05   2.640920E-03
  b1g  17  -2.138008E-03  -6.379697E-04  -4.530245E-03   8.247889E-05   2.601446E-04  -3.295508E-04   1.908560E-05  -4.713612E-04
  b1g  18   1.207228E-05   5.876179E-06  -1.397964E-06  -1.080276E-03  -5.925151E-07   8.340634E-06  -3.638265E-03  -1.909831E-05
  b1g  19   4.081670E-05   9.106187E-04   8.368988E-04  -1.615867E-05  -4.161184E-05   6.863974E-04  -8.188476E-06   9.898558E-04
  b1g  20  -1.852692E-04   1.131605E-03  -5.426552E-04   9.826790E-06  -1.222402E-03   1.300834E-03   1.243089E-06   7.459124E-04
  b1g  21  -1.058329E-04  -3.115604E-04  -3.103010E-04   5.151058E-06  -2.542284E-04  -3.265974E-05   2.308729E-06   2.195874E-04
  b1g  22  -1.137201E-06   1.132849E-06  -6.587906E-06  -2.045504E-04   1.761984E-06   2.167671E-06  -1.068279E-03  -2.688335E-06
  b1g  23   7.651796E-04   2.184872E-05   7.424054E-05  -3.633432E-06  -4.191599E-04  -8.555873E-05  -7.670077E-06  -6.632531E-04
  b1g  24  -5.726410E-04   5.223074E-04  -8.004301E-04   1.535370E-05   3.683786E-04   2.721872E-04   2.446293E-06   7.176875E-04
  b1g  25  -1.043461E-06   1.541401E-06  -7.738633E-06  -2.392326E-04  -2.730552E-07   8.581501E-07  -7.654371E-05  -8.837263E-09
  b1g  26  -1.663459E-04  -5.099873E-05  -3.956246E-04   7.547706E-06   2.854879E-06   9.250323E-05   9.946763E-07   1.507812E-04
  b1g  27  -3.174003E-04  -1.335704E-05  -3.993294E-04   7.610563E-06   2.143694E-05   1.402532E-04   4.001240E-06   2.672890E-04
  b1g  28  -3.174566E-04   4.668416E-04  -2.628006E-04   2.229104E-06  -5.476523E-05   3.536550E-04  -7.471620E-06   4.007649E-04
  b1g  29  -6.176748E-07   4.540465E-07  -4.362564E-07   2.022631E-05  -2.886515E-07  -6.578546E-08   2.360465E-04   1.536476E-06
  b1g  30   1.205485E-04   2.415359E-04   2.073929E-04  -3.562319E-06  -1.318858E-04   1.824105E-04  -1.154305E-06   2.588362E-04

               b1g  17        b1g  18        b1g  19        b1g  20        b1g  21        b1g  22        b1g  23        b1g  24
  b1g   2   1.488742E-03   4.220022E-04   1.137051E-03   2.223193E-03   3.232066E-04   1.923228E-04   2.823003E-03  -1.337837E-03
  b1g   3   5.932277E-04  -7.768676E-04  -3.930285E-04   1.397608E-03  -2.249259E-03  -2.213217E-04  -3.275431E-03  -3.165167E-04
  b1g   4   9.655491E-05  -3.050424E-07   1.265972E-04   8.966365E-07  -1.218543E-04   1.841221E-07  -1.471564E-04   1.493990E-04
  b1g   5   4.123059E-04  -3.111716E-06   1.895982E-04  -8.861980E-05  -1.228919E-04   4.439772E-07  -4.288520E-04   3.017239E-04
  b1g   6   9.240129E-05   9.863219E-07   2.855044E-04   2.757162E-04  -2.336883E-04   4.013725E-07  -8.140049E-05   2.953586E-04
  b1g   7  -5.140016E-08  -6.272113E-05   2.830302E-07   9.384793E-07  -4.042728E-07  -6.227633E-07   3.461213E-07   3.662643E-07
  b1g   8  -3.117090E-04   2.985059E-06  -4.479539E-05   5.813456E-04   4.574915E-05  -8.784473E-07   4.224393E-04  -2.375969E-04
  b1g   9  -2.138008E-03   1.207228E-05   4.081670E-05  -1.852692E-04  -1.058329E-04  -1.137201E-06   7.651796E-04  -5.726410E-04
  b1g  10  -6.379697E-04   5.876179E-06   9.106187E-04   1.131605E-03  -3.115604E-04   1.132849E-06   2.184872E-05   5.223074E-04
  b1g  11  -4.530245E-03  -1.397964E-06   8.368988E-04  -5.426552E-04  -3.103010E-04  -6.587906E-06   7.424054E-05  -8.004301E-04
  b1g  12   8.247889E-05  -1.080276E-03  -1.615867E-05   9.826790E-06   5.151058E-06  -2.045504E-04  -3.633432E-06   1.535370E-05
  b1g  13   2.601446E-04  -5.925151E-07  -4.161184E-05  -1.222402E-03  -2.542284E-04   1.761984E-06  -4.191599E-04   3.683786E-04
  b1g  14  -3.295508E-04   8.340634E-06   6.863974E-04   1.300834E-03  -3.265974E-05   2.167671E-06  -8.555873E-05   2.721872E-04
  b1g  15   1.908560E-05  -3.638265E-03  -8.188476E-06   1.243089E-06   2.308729E-06  -1.068279E-03  -7.670077E-06   2.446293E-06
  b1g  16  -4.713612E-04  -1.909831E-05   9.898558E-04   7.459124E-04   2.195874E-04  -2.688335E-06  -6.632531E-04   7.176875E-04
  b1g  17   3.907113E-03  -3.033915E-05  -9.545455E-04   3.819372E-04  -2.418387E-06  -1.332906E-06   9.551007E-04   9.724165E-04
  b1g  18  -3.033915E-05   2.887140E-03   7.816132E-06  -4.446683E-06  -4.077050E-06   1.175811E-03   3.402124E-07  -1.027760E-05
  b1g  19  -9.545455E-04   7.816132E-06   8.505406E-04   1.018401E-04   6.385797E-05   2.375520E-06  -3.120486E-04  -2.708604E-05
  b1g  20   3.819372E-04  -4.446683E-06   1.018401E-04   1.570943E-03   3.152821E-04  -1.278359E-06   3.516350E-04  -6.269667E-05
  b1g  21  -2.418387E-06  -4.077050E-06   6.385797E-05   3.152821E-04   7.617811E-04   5.594286E-07  -1.038011E-04   2.597835E-05
  b1g  22  -1.332906E-06   1.175811E-03   2.375520E-06  -1.278359E-06   5.594286E-07   6.511554E-04   1.111306E-06   7.648486E-07
  b1g  23   9.551007E-04   3.402124E-07  -3.120486E-04   3.516350E-04  -1.038011E-04   1.111306E-06   1.468095E-03   1.012648E-04
  b1g  24   9.724165E-04  -1.027760E-05  -2.708604E-05  -6.269667E-05   2.597835E-05   7.648486E-07   1.012648E-04   1.042743E-03
  b1g  25   3.396570E-06  -5.961595E-04  -8.328368E-07   1.086429E-06   1.756259E-06  -4.603851E-04   3.296961E-06   3.445447E-06
  b1g  26   5.337477E-04  -4.046566E-06  -2.302951E-04   4.100016E-04   2.624206E-04   8.518417E-07   2.976982E-04   3.005563E-04
  b1g  27   5.351375E-05  -4.134806E-06   2.409136E-04   1.816551E-04   1.432380E-04  -8.256061E-07  -2.966173E-04   1.527720E-05
  b1g  28   1.449452E-04   3.585579E-06   2.278225E-04   9.531869E-05  -5.766340E-05   2.091291E-06   2.008504E-04   6.303889E-05
  b1g  29   1.737414E-06  -4.009904E-04  -6.652733E-07   4.769792E-07   1.160483E-06  -3.080988E-04   7.255574E-07   1.323244E-06
  b1g  30  -1.679997E-04   8.319473E-07   1.315248E-04   2.415399E-05  -4.052491E-06   2.473359E-07  -1.920882E-04   1.218440E-04

               b1g  25        b1g  26        b1g  27        b1g  28        b1g  29        b1g  30
  b1g   2  -1.057886E-04   2.237160E-03   9.319489E-04   2.638257E-03  -9.868272E-05  -5.020570E-04
  b1g   3   6.747700E-05   1.385623E-03   3.730759E-03  -2.858927E-03   3.207798E-05   5.371872E-04
  b1g   4   1.136812E-07  -2.460279E-06   6.455613E-06   1.015990E-04   1.374037E-07  -1.058916E-05
  b1g   5   5.861275E-08   3.364156E-05   8.201954E-05   1.783048E-04   2.669682E-07  -4.641577E-05
  b1g   6   8.050153E-07  -1.165474E-05  -1.350687E-05   2.140366E-04   2.939152E-07   4.239689E-05
  b1g   7  -3.965976E-05   1.304458E-08  -1.578650E-07   1.068813E-07  -5.297433E-06   1.888170E-07
  b1g   8   2.839535E-07  -4.494861E-05  -9.004725E-05  -2.592333E-05   9.120385E-09   8.480075E-05
  b1g   9  -1.043461E-06  -1.663459E-04  -3.174003E-04  -3.174566E-04  -6.176748E-07   1.205485E-04
  b1g  10   1.541401E-06  -5.099873E-05  -1.335704E-05   4.668416E-04   4.540465E-07   2.415359E-04
  b1g  11  -7.738633E-06  -3.956246E-04  -3.993294E-04  -2.628006E-04  -4.362564E-07   2.073929E-04
  b1g  12  -2.392326E-04   7.547706E-06   7.610563E-06   2.229104E-06   2.022631E-05  -3.562319E-06
  b1g  13  -2.730552E-07   2.854879E-06   2.143694E-05  -5.476523E-05  -2.886515E-07  -1.318858E-04
  b1g  14   8.581501E-07   9.250323E-05   1.402532E-04   3.536550E-04  -6.578546E-08   1.824105E-04
  b1g  15  -7.654371E-05   9.946763E-07   4.001240E-06  -7.471620E-06   2.360465E-04  -1.154305E-06
  b1g  16  -8.837263E-09   1.507812E-04   2.672890E-04   4.007649E-04   1.536476E-06   2.588362E-04
  b1g  17   3.396570E-06   5.337477E-04   5.351375E-05   1.449452E-04   1.737414E-06  -1.679997E-04
  b1g  18  -5.961595E-04  -4.046566E-06  -4.134806E-06   3.585579E-06  -4.009904E-04   8.319473E-07
  b1g  19  -8.328368E-07  -2.302951E-04   2.409136E-04   2.278225E-04  -6.652733E-07   1.315248E-04
  b1g  20   1.086429E-06   4.100016E-04   1.816551E-04   9.531869E-05   4.769792E-07   2.415399E-05
  b1g  21   1.756259E-06   2.624206E-04   1.432380E-04  -5.766340E-05   1.160483E-06  -4.052491E-06
  b1g  22  -4.603851E-04   8.518417E-07  -8.256061E-07   2.091291E-06  -3.080988E-04   2.473359E-07
  b1g  23   3.296961E-06   2.976982E-04  -2.966173E-04   2.008504E-04   7.255574E-07  -1.920882E-04
  b1g  24   3.445447E-06   3.005563E-04   1.527720E-05   6.303889E-05   1.323244E-06   1.218440E-04
  b1g  25   6.194003E-04   1.722108E-06   2.517584E-07  -1.750560E-07   2.881265E-04  -2.775631E-07
  b1g  26   1.722108E-06   8.588360E-04   2.179647E-05  -1.284265E-04   6.229550E-07  -2.283264E-04
  b1g  27   2.517584E-07   2.179647E-05   4.048933E-04  -9.094424E-05   1.470481E-07   2.195224E-05
  b1g  28  -1.750560E-07  -1.284265E-04  -9.094424E-05   4.987028E-04  -2.337991E-07   1.167762E-04
  b1g  29   2.881265E-04   6.229550E-07   1.470481E-07  -2.337991E-07   3.354689E-04  -6.878645E-08
  b1g  30  -2.775631E-07  -2.283264E-04   2.195224E-05   1.167762E-04  -6.878645E-08   3.235079E-04

Natural orbital populations,block 4
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98148362     1.97390147     0.01445169     0.01032791     0.00961380     0.00463589     0.00233174
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00176698     0.00148824     0.00092085     0.00077590     0.00046809     0.00029930     0.00021695     0.00013814
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00010132     0.00008610     0.00006588     0.00005051     0.00002028     0.00001924     0.00000813     0.00000640
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00000288     0.00000128     0.00000109     0.00000021     0.00000003     0.00000002

          modens reordered block   1

               b1u   1        b1u   2        b1u   3        b1u   4        b1u   5        b1u   6        b1u   7        b1u   8
  b1u   1    1.94732       1.185206E-03   1.659029E-03   4.265751E-04   8.586184E-04   3.649100E-05  -1.701114E-03   1.537338E-04
  b1u   2   1.185206E-03   0.466098       0.167415      -2.470240E-04  -1.251816E-03   3.947486E-04   1.535304E-03   4.390177E-03
  b1u   3   1.659029E-03   0.167415       6.294349E-02  -8.057986E-05   4.998518E-05  -5.248360E-06   1.246965E-03   1.751576E-03
  b1u   4   4.265751E-04  -2.470240E-04  -8.057986E-05   7.285001E-04   7.694787E-05   6.832865E-04   3.895912E-05  -6.410155E-04
  b1u   5   8.586184E-04  -1.251816E-03   4.998518E-05   7.694787E-05   2.043729E-03   2.577132E-05  -1.619104E-03  -3.508750E-05
  b1u   6   3.649100E-05   3.947486E-04  -5.248360E-06   6.832865E-04   2.577132E-05   1.007081E-03   4.119089E-06  -4.840344E-04
  b1u   7  -1.701114E-03   1.535304E-03   1.246965E-03   3.895912E-05  -1.619104E-03   4.119089E-06   2.201956E-03   2.283576E-05
  b1u   8   1.537338E-04   4.390177E-03   1.751576E-03  -6.410155E-04  -3.508750E-05  -4.840344E-04   2.283576E-05   1.442021E-03
  b1u   9   5.368826E-03  -6.881187E-04  -2.509894E-03  -2.156352E-04  -5.764477E-04  -4.840332E-05  -7.715797E-04   2.047179E-05
  b1u  10  -2.042007E-06   1.546980E-03   4.189882E-04   6.477423E-04   1.197755E-05   1.179417E-03   1.482339E-05  -9.806703E-05
  b1u  11   3.707761E-05   4.981609E-03   1.916865E-03  -5.500553E-04  -3.109090E-05  -2.532301E-04   2.046380E-05   9.946970E-04
  b1u  12   3.579619E-03  -7.068066E-04   1.342405E-04   4.706262E-05   1.751388E-03   1.263602E-05  -1.423408E-03  -1.701480E-05
  b1u  13  -1.217297E-04  -2.858212E-03  -1.019110E-03  -1.330435E-04   3.111558E-06  -4.068444E-04  -7.386601E-06  -9.575617E-04
  b1u  14   3.036624E-05   1.255810E-04   1.156942E-04  -1.995210E-04  -6.829666E-06  -7.143222E-04   1.884472E-06   3.734766E-04
  b1u  15   2.434814E-03  -4.053297E-04  -2.070091E-04  -1.080095E-06   2.522029E-04  -1.000984E-06  -6.448753E-04  -2.959776E-06
  b1u  16   1.810636E-05   5.868433E-04   1.962764E-04   3.532213E-05   1.932690E-06   6.127324E-05  -3.202603E-06  -1.650840E-04

               b1u   9        b1u  10        b1u  11        b1u  12        b1u  13        b1u  14        b1u  15        b1u  16
  b1u   1   5.368826E-03  -2.042007E-06   3.707761E-05   3.579619E-03  -1.217297E-04   3.036624E-05   2.434814E-03   1.810636E-05
  b1u   2  -6.881187E-04   1.546980E-03   4.981609E-03  -7.068066E-04  -2.858212E-03   1.255810E-04  -4.053297E-04   5.868433E-04
  b1u   3  -2.509894E-03   4.189882E-04   1.916865E-03   1.342405E-04  -1.019110E-03   1.156942E-04  -2.070091E-04   1.962764E-04
  b1u   4  -2.156352E-04   6.477423E-04  -5.500553E-04   4.706262E-05  -1.330435E-04  -1.995210E-04  -1.080095E-06   3.532213E-05
  b1u   5  -5.764477E-04   1.197755E-05  -3.109090E-05   1.751388E-03   3.111558E-06  -6.829666E-06   2.522029E-04   1.932690E-06
  b1u   6  -4.840332E-05   1.179417E-03  -2.532301E-04   1.263602E-05  -4.068444E-04  -7.143222E-04  -1.000984E-06   6.127324E-05
  b1u   7  -7.715797E-04   1.482339E-05   2.046380E-05  -1.423408E-03  -7.386601E-06   1.884472E-06  -6.448753E-04  -3.202603E-06
  b1u   8   2.047179E-05  -9.806703E-05   9.946970E-04  -1.701480E-05  -9.575617E-04   3.734766E-04  -2.959776E-06  -1.650840E-04
  b1u   9   2.827186E-03  -4.220456E-05   1.722679E-05  -4.074436E-04  -4.393609E-06   8.644141E-06   1.135982E-04  -1.653915E-06
  b1u  10  -4.220456E-05   1.702629E-03   8.492358E-05   1.784097E-06  -1.125210E-03  -9.322684E-04  -3.526100E-06   1.807179E-04
  b1u  11   1.722679E-05   8.492358E-05   1.250435E-03  -1.385585E-05  -7.192728E-04  -1.540048E-04  -3.190294E-06   1.024609E-04
  b1u  12  -4.074436E-04   1.784097E-06  -1.385585E-05   2.035789E-03   2.795249E-06  -6.731926E-06   2.759591E-04   2.541047E-06
  b1u  13  -4.393609E-06  -1.125210E-03  -7.192728E-04   2.795249E-06   2.119522E-03  -3.695127E-05   7.237240E-07  -4.912114E-05
  b1u  14   8.644141E-06  -9.322684E-04  -1.540048E-04  -6.731926E-06  -3.695127E-05   1.258999E-03   1.020936E-06  -2.950878E-04
  b1u  15   1.135982E-04  -3.526100E-06  -3.190294E-06   2.759591E-04   7.237240E-07   1.020936E-06   6.199861E-04   1.419309E-06
  b1u  16  -1.653915E-06   1.807179E-04   1.024609E-04   2.541047E-06  -4.912114E-05  -2.950878E-04   1.419309E-06   4.714075E-04

Natural orbital populations,block 5
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.94734774     0.52668693     0.00546086     0.00512610     0.00411605     0.00354596     0.00136057     0.00070322
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00048168     0.00047312     0.00034745     0.00020877     0.00012300     0.00004432     0.00003606     0.00000773

          modens reordered block   1

               b2g   1        b2g   2        b2g   3        b2g   4        b2g   5        b2g   6        b2g   7        b2g   8
  b2g   1    1.46727       2.234555E-03  -6.680289E-04   1.900520E-03  -2.326018E-04  -5.943801E-04  -3.819425E-03  -8.334088E-03
  b2g   2   2.234555E-03   9.972515E-04   2.116435E-05  -9.276875E-05  -8.615974E-04   4.168670E-05  -1.089006E-03  -1.432109E-04
  b2g   3  -6.680289E-04   2.116435E-05   1.495064E-03  -5.325772E-03   1.437584E-05   2.174498E-03   6.804217E-06  -1.293801E-05
  b2g   4   1.900520E-03  -9.276875E-05  -5.325772E-03   1.944256E-02  -4.088179E-05  -8.357478E-03  -6.817022E-06   5.911060E-05
  b2g   5  -2.326018E-04  -8.615974E-04   1.437584E-05  -4.088179E-05   1.280729E-03   1.627966E-05   1.326561E-03  -5.872864E-04
  b2g   6  -5.943801E-04   4.168670E-05   2.174498E-03  -8.357478E-03   1.627966E-05   4.048575E-03   2.738879E-07  -3.199557E-05
  b2g   7  -3.819425E-03  -1.089006E-03   6.804217E-06  -6.817022E-06   1.326561E-03   2.738879E-07   1.693476E-03  -2.612264E-04
  b2g   8  -8.334088E-03  -1.432109E-04  -1.293801E-05   5.911060E-05  -5.872864E-04  -3.199557E-05  -2.612264E-04   1.372894E-03
  b2g   9  -5.276626E-04   1.068994E-06   1.311079E-03  -5.373557E-03   1.872150E-05   2.957093E-03   1.826614E-05   4.594052E-06
  b2g  10  -8.811596E-03  -9.590032E-04  -3.517824E-05   1.682890E-04   5.151086E-04  -9.741203E-05   8.452786E-04   7.449047E-04
  b2g  11  -1.300209E-04   1.267001E-05   5.406293E-04  -2.038559E-03   9.022546E-07   8.080165E-04  -4.227004E-06  -6.877178E-06
  b2g  12   3.424742E-03  -4.150807E-04   3.714478E-06  -1.166607E-05   1.113188E-03   7.195212E-06   1.080929E-03  -1.211580E-03
  b2g  13  -1.144048E-05   5.838258E-06  -4.815143E-06  -3.980204E-04  -4.653090E-06   7.663688E-04  -5.848767E-06  -2.670991E-06
  b2g  14  -4.231329E-03  -2.991911E-05  -2.798111E-06   1.490407E-05   8.318090E-05  -1.057203E-05  -3.785833E-05   2.715130E-04
  b2g  15  -3.364033E-03  -4.537056E-05  -3.431126E-06   1.579915E-05   9.053605E-05  -9.147977E-06   3.621549E-04   2.852331E-04
  b2g  16   8.758346E-05   6.072142E-07  -1.117313E-05   1.393801E-05  -1.177194E-06   4.672214E-05  -2.019888E-06   1.566046E-07

               b2g   9        b2g  10        b2g  11        b2g  12        b2g  13        b2g  14        b2g  15        b2g  16
  b2g   1  -5.276626E-04  -8.811596E-03  -1.300209E-04   3.424742E-03  -1.144048E-05  -4.231329E-03  -3.364033E-03   8.758346E-05
  b2g   2   1.068994E-06  -9.590032E-04   1.267001E-05  -4.150807E-04   5.838258E-06  -2.991911E-05  -4.537056E-05   6.072142E-07
  b2g   3   1.311079E-03  -3.517824E-05   5.406293E-04   3.714478E-06  -4.815143E-06  -2.798111E-06  -3.431126E-06  -1.117313E-05
  b2g   4  -5.373557E-03   1.682890E-04  -2.038559E-03  -1.166607E-05  -3.980204E-04   1.490407E-05   1.579915E-05   1.393801E-05
  b2g   5   1.872150E-05   5.151086E-04   9.022546E-07   1.113188E-03  -4.653090E-06   8.318090E-05   9.053605E-05  -1.177194E-06
  b2g   6   2.957093E-03  -9.741203E-05   8.080165E-04   7.195212E-06   7.663688E-04  -1.057203E-05  -9.147977E-06   4.672214E-05
  b2g   7   1.826614E-05   8.452786E-04  -4.227004E-06   1.080929E-03  -5.848767E-06  -3.785833E-05   3.621549E-04  -2.019888E-06
  b2g   8   4.594052E-06   7.449047E-04  -6.877178E-06  -1.211580E-03  -2.670991E-06   2.715130E-04   2.852331E-04   1.566046E-07
  b2g   9   2.457813E-03  -3.437113E-05   5.082398E-04  -5.521254E-06   1.029007E-03  -1.450058E-06  -9.375932E-06   1.516464E-04
  b2g  10  -3.437113E-05   1.571905E-03  -2.110788E-05  -1.015118E-04  -3.512186E-05   2.377072E-04  -7.770983E-05  -4.484572E-06
  b2g  11   5.082398E-04  -2.110788E-05   5.836103E-04  -1.018607E-06  -2.369346E-04  -1.367049E-06  -1.760979E-06  -6.308049E-05
  b2g  12  -5.521254E-06  -1.015118E-04  -1.018607E-06   1.695745E-03  -4.101676E-06  -2.506806E-04  -5.260201E-05  -1.248974E-06
  b2g  13   1.029007E-03  -3.512186E-05  -2.369346E-04  -4.101676E-06   1.142555E-03  -3.008671E-06  -1.058497E-06   2.895973E-04
  b2g  14  -1.450058E-06   2.377072E-04  -1.367049E-06  -2.506806E-04  -3.008671E-06   6.576156E-04   2.078846E-04  -2.450777E-07
  b2g  15  -9.375932E-06  -7.770983E-05  -1.760979E-06  -5.260201E-05  -1.058497E-06   2.078846E-04   9.378610E-04  -1.288636E-06
  b2g  16   1.516464E-04  -4.484572E-06  -6.308049E-05  -1.248974E-06   2.895973E-04  -2.450777E-07  -1.288636E-06   3.557915E-04

Natural orbital populations,block 6
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.46741065     0.02642284     0.00497927     0.00298596     0.00227680     0.00114412     0.00060619     0.00045520
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00030526     0.00022082     0.00006165     0.00005832     0.00005389     0.00001154     0.00000524     0.00000127

          modens reordered block   1

               b3g   1        b3g   2        b3g   3        b3g   4        b3g   5        b3g   6        b3g   7        b3g   8
  b3g   1    1.45223      -3.607565E-04  -2.448002E-04  -7.465151E-04  -3.828409E-03  -7.763400E-03  -8.162214E-03   8.045046E-04
  b3g   2  -3.607565E-04   9.690278E-04   8.519338E-04   5.371537E-06   1.071045E-03   1.323260E-04   9.219328E-04  -4.407829E-06
  b3g   3  -2.448002E-04   8.519338E-04   1.274638E-03   8.297279E-07   1.319835E-03  -5.825146E-04   5.135285E-04  -2.321512E-06
  b3g   4  -7.465151E-04   5.371537E-06   8.297279E-07   2.400337E-03  -4.186263E-06   6.556285E-06   1.169683E-05  -2.511028E-03
  b3g   5  -3.828409E-03   1.071045E-03   1.319835E-03  -4.186263E-06   1.685776E-03  -2.550319E-04   8.358625E-04   4.400191E-06
  b3g   6  -7.763400E-03   1.323260E-04  -5.825146E-04   6.556285E-06  -2.550319E-04   1.375975E-03   7.302677E-04  -5.422397E-06
  b3g   7  -8.162214E-03   9.219328E-04   5.135285E-04   1.169683E-05   8.358625E-04   7.302677E-04   1.523727E-03  -9.471102E-06
  b3g   8   8.045046E-04  -4.407829E-06  -2.321512E-06  -2.511028E-03   4.400191E-06  -5.422397E-06  -9.471102E-06   3.337599E-03
  b3g   9   3.065670E-03   4.213245E-04   1.111588E-03  -2.259546E-08   1.080743E-03  -1.204878E-03  -8.306380E-05  -8.932978E-07
  b3g  10  -3.853445E-03   2.025915E-05   8.044039E-05   5.103001E-07  -4.662089E-05   2.621116E-04   2.267936E-04  -7.948441E-07
  b3g  11  -3.497176E-03   3.878877E-05   8.854471E-05  -4.865671E-06   3.569570E-04   2.887053E-04  -8.118294E-05   5.464125E-06

               b3g   9        b3g  10        b3g  11
  b3g   1   3.065670E-03  -3.853445E-03  -3.497176E-03
  b3g   2   4.213245E-04   2.025915E-05   3.878877E-05
  b3g   3   1.111588E-03   8.044039E-05   8.854471E-05
  b3g   4  -2.259546E-08   5.103001E-07  -4.865671E-06
  b3g   5   1.080743E-03  -4.662089E-05   3.569570E-04
  b3g   6  -1.204878E-03   2.621116E-04   2.887053E-04
  b3g   7  -8.306380E-05   2.267936E-04  -8.118294E-05
  b3g   8  -8.932978E-07  -7.948441E-07   5.464125E-06
  b3g   9   1.689823E-03  -2.460761E-04  -5.815294E-05
  b3g  10  -2.460761E-04   6.468071E-04   2.032242E-04
  b3g  11  -5.815294E-05   2.032242E-04   9.457292E-04

Natural orbital populations,block 7
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.45235675     0.00542260     0.00494992     0.00294488     0.00114399     0.00060501     0.00031459     0.00021921
              MO     9       MO    10       MO    11
  occ(*)=     0.00006139     0.00005266     0.00001143

          modens reordered block   1

               au    1        au    2        au    3        au    4        au    5        au    6        au    7        au    8
  au    1   0.540895      -1.784019E-03  -3.291934E-04   6.301479E-03  -1.505130E-03   7.075769E-05  -6.315584E-03   3.795567E-03
  au    2  -1.784019E-03   7.992693E-04   7.271944E-04  -7.087564E-04  -6.750816E-04  -2.251723E-07   6.023557E-04   1.328008E-04
  au    3  -3.291934E-04   7.271944E-04   1.029430E-03  -5.162666E-04  -1.192589E-03   1.356037E-06   2.698295E-04   4.118962E-04
  au    4   6.301479E-03  -7.087564E-04  -5.162666E-04   1.506063E-03   1.157860E-04   2.707436E-06  -1.037739E-03   9.703125E-04
  au    5  -1.505130E-03  -6.750816E-04  -1.192589E-03   1.157860E-04   1.707836E-03  -2.988186E-06   8.224796E-05  -1.130744E-03
  au    6   7.075769E-05  -2.251723E-07   1.356037E-06   2.707436E-06  -2.988186E-06   8.835930E-05  -3.015359E-06   3.766483E-06
  au    7  -6.315584E-03   6.023557E-04   2.698295E-04  -1.037739E-03   8.224796E-05  -3.015359E-06   1.277965E-03  -7.342718E-04
  au    8   3.795567E-03   1.328008E-04   4.118962E-04   9.703125E-04  -1.130744E-03   3.766483E-06  -7.342718E-04   2.132600E-03
  au    9   4.038919E-04  -2.007846E-04  -7.140543E-04   3.840460E-04   9.284306E-04  -1.513011E-06   1.514869E-04   4.043330E-05
  au   10   6.697510E-04   3.347593E-05   5.993073E-05  -1.653172E-04  -1.789229E-04   4.681322E-07  -1.032503E-04   4.815969E-05
  au   11   1.346509E-04  -5.835974E-07   7.586182E-08   4.719748E-06  -1.870249E-06   1.677492E-04  -3.844012E-06   5.950261E-06

               au    9        au   10        au   11
  au    1   4.038919E-04   6.697510E-04   1.346509E-04
  au    2  -2.007846E-04   3.347593E-05  -5.835974E-07
  au    3  -7.140543E-04   5.993073E-05   7.586182E-08
  au    4   3.840460E-04  -1.653172E-04   4.719748E-06
  au    5   9.284306E-04  -1.789229E-04  -1.870249E-06
  au    6  -1.513011E-06   4.681322E-07   1.677492E-04
  au    7   1.514869E-04  -1.032503E-04  -3.844012E-06
  au    8   4.043330E-05   4.815969E-05   5.950261E-06
  au    9   1.257659E-03  -2.943441E-04   1.242043E-06
  au   10  -2.943441E-04   4.708359E-04   4.741120E-08
  au   11   1.242043E-06   4.741120E-08   5.688515E-04

Natural orbital populations,block 8
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     0.54108167     0.00416142     0.00358955     0.00138106     0.00062156     0.00047766     0.00021396     0.00012489
              MO     9       MO    10       MO    11
  occ(*)=     0.00003902     0.00003558     0.00000779


 total number of electrons =   42.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        ag  partial gross atomic populations
   ao class       1ag        2ag        3ag        4ag        5ag        6ag 
    C1_ s       0.024115   1.965539   0.527976   0.474216  -0.012031   0.057586
    C1_ p       0.003468  -0.000704   0.037687   0.126377   0.351221   0.266435
    C1_ d       0.000470   0.000179   0.004979  -0.112259   0.038968  -0.026618
    C2_ s       1.965551   0.024164   1.077390   0.238548  -0.023655   0.028793
    C2_ p       0.000801   0.007425   0.079441   0.118499   0.711336   1.287252
    C2_ d      -0.000546   0.001050   0.010404  -0.052960   0.079240   0.010099
    H1_ s       0.001179   0.000814   0.077161   0.771902   0.277637   0.278641
    H1_ p       0.000351   0.000587   0.005767   0.031403   0.000429  -0.036975
    H2_ s       0.003058   0.000734   0.154865   0.383862   0.555347   0.138448
    H2_ p       0.001552   0.000211   0.011631   0.001886   0.000927  -0.029761
 
   ao class       7ag        8ag        9ag       10ag       11ag       12ag 
    C1_ s       0.004356   0.000460   0.000478  -0.000052   0.001113   0.000314
    C1_ p       0.001553   0.000771   0.003883   0.001400   0.000120   0.000348
    C1_ d       0.000400   0.000920  -0.001162   0.000823   0.000876   0.000179
    C2_ s       0.002148   0.000953   0.000239  -0.000091   0.002220   0.000173
    C2_ p       0.004279   0.001488   0.002439   0.002824   0.000235   0.002227
    C2_ d       0.001459   0.001844  -0.000072   0.001665   0.001744   0.000788
    H1_ s       0.000185   0.002194   0.002838  -0.000017  -0.000363  -0.000206
    H1_ p      -0.000053  -0.000120   0.000176   0.000048   0.000032   0.000061
    H2_ s       0.000079   0.004387   0.001435  -0.000027  -0.000717  -0.000104
    H2_ p       0.000059  -0.000245   0.000078   0.000095   0.000063   0.000857
 
   ao class      13ag       14ag       15ag       16ag       17ag       18ag 
    C1_ s      -0.000338   0.000609  -0.000143   0.000007   0.000075  -0.000096
    C1_ p       0.001329  -0.000043  -0.000132   0.000384   0.000236   0.000012
    C1_ d       0.000690   0.000293   0.000748   0.000096   0.000200   0.000157
    C2_ s      -0.000168   0.001213  -0.000069   0.000017   0.000038  -0.000049
    C2_ p       0.000628  -0.000084   0.000310   0.000789   0.000084   0.000729
    C2_ d       0.001292   0.000581   0.000671   0.000196   0.000135  -0.000119
    H1_ s      -0.001076  -0.000268  -0.000213  -0.000099  -0.000032   0.000000
    H1_ p       0.000305   0.000001   0.000257   0.000028   0.000049   0.000027
    H2_ s      -0.000540  -0.000528  -0.000109  -0.000202  -0.000017   0.000000
    H2_ p       0.000208   0.000002   0.000165   0.000057   0.000154   0.000113
 
   ao class      19ag       20ag       21ag       22ag       23ag       24ag 
    C1_ s       0.000087   0.000016   0.000045  -0.000141  -0.000013  -0.000015
    C1_ p       0.000169  -0.000007   0.000087   0.000028   0.000002  -0.000017
    C1_ d       0.000019   0.000046   0.000060   0.000088   0.000065   0.000027
    C2_ s       0.000022   0.000078   0.000023  -0.000285  -0.000007  -0.000027
    C2_ p       0.000052   0.000030  -0.000025   0.000061  -0.000080  -0.000033
    C2_ d       0.000073   0.000128   0.000121   0.000182   0.000188   0.000053
    H1_ s      -0.000030   0.000007   0.000094   0.000120  -0.000001   0.000052
    H1_ p      -0.000035   0.000050  -0.000112  -0.000014   0.000014  -0.000008
    H2_ s      -0.000021   0.000027   0.000045   0.000247  -0.000001   0.000104
    H2_ p       0.000130   0.000087  -0.000039  -0.000031  -0.000028  -0.000017
 
   ao class      25ag       26ag       27ag       28ag       29ag       30ag 
    C1_ s      -0.000162  -0.000011   0.000035   0.000038   0.000001  -0.000001
    C1_ p       0.000081  -0.000006   0.000015  -0.000025  -0.000002  -0.000015
    C1_ d       0.000036   0.000008   0.000026   0.000005   0.000000   0.000001
    C2_ s      -0.000085  -0.000015   0.000018   0.000019   0.000002  -0.000000
    C2_ p       0.000152  -0.000014   0.000008  -0.000026  -0.000003   0.000005
    C2_ d       0.000084   0.000013  -0.000020   0.000062  -0.000000   0.000001
    H1_ s       0.000014   0.000043  -0.000020   0.000011  -0.000002   0.000029
    H1_ p      -0.000003  -0.000003   0.000018   0.000002   0.000010  -0.000012
    H2_ s       0.000008   0.000085  -0.000010   0.000006  -0.000004   0.000015
    H2_ p      -0.000024  -0.000005   0.000016  -0.000042   0.000020  -0.000003
 
   ao class      31ag       32ag       33ag       34ag       35ag       36ag 
    C1_ s      -0.000001   0.000002   0.000003   0.000002   0.000012  -0.000001
    C1_ p       0.000007  -0.000006  -0.000001   0.000000  -0.000008   0.000005
    C1_ d      -0.000000   0.000000   0.000002  -0.000000  -0.000000  -0.000000
    C2_ s      -0.000001  -0.000001   0.000001   0.000018   0.000004  -0.000000
    C2_ p       0.000016  -0.000003  -0.000007  -0.000000  -0.000007   0.000001
    C2_ d      -0.000001  -0.000000  -0.000010  -0.000000   0.000000  -0.000000
    H1_ s      -0.000001  -0.000002   0.000002  -0.000007   0.000002  -0.000002
    H1_ p      -0.000000   0.000011   0.000001   0.000001   0.000000   0.000000
    H2_ s      -0.000002  -0.000001   0.000001  -0.000011   0.000001  -0.000001
    H2_ p      -0.000001   0.000008   0.000014   0.000003  -0.000000   0.000000
 
   ao class      37ag       38ag       39ag 
    C1_ s      -0.000000  -0.000000   0.000000
    C1_ p       0.000000  -0.000000  -0.000000
    C1_ d      -0.000000  -0.000000   0.000000
    C2_ s      -0.000000  -0.000000   0.000000
    C2_ p       0.000001  -0.000000  -0.000000
    C2_ d       0.000000  -0.000000   0.000000
    H1_ s      -0.000000   0.000000   0.000000
    H1_ p      -0.000000  -0.000000   0.000000
    H2_ s      -0.000000   0.000000   0.000000
    H2_ p      -0.000000  -0.000000   0.000000

                        b3u partial gross atomic populations
   ao class       1b3u       2b3u       3b3u       4b3u       5b3u       6b3u
    C1_ s       0.022924   1.952003   0.718032   0.079517   0.078158  -0.000084
    C1_ p      -0.000455   0.001641  -0.038660   0.216137   0.725249   0.001106
    C1_ d      -0.000102   0.000126   0.042428  -0.185859  -0.018025   0.001063
    C2_ s       1.958708   0.025923   0.366089   0.170240   0.039595  -0.000053
    C2_ p       0.013886   0.014358   0.210075   0.427650   0.560611   0.008277
    C2_ d       0.003923   0.004386  -0.029779  -0.371442  -0.040218   0.002320
    H1_ s       0.000617   0.000533   0.446605   0.551676   0.434541   0.000922
    H1_ p       0.000315   0.000354   0.040363  -0.002420  -0.016504  -0.000040
    H2_ s       0.000067   0.000410   0.224253   1.097962   0.217588   0.000465
    H2_ p       0.000116   0.000266   0.005501  -0.004640  -0.004304   0.000092
 
   ao class       7b3u       8b3u       9b3u      10b3u      11b3u      12b3u
    C1_ s       0.003400   0.000835   0.000212   0.000050  -0.000538   0.000145
    C1_ p       0.000358   0.001859   0.002726   0.001814  -0.000136  -0.000093
    C1_ d      -0.000252   0.000137  -0.001068   0.000913   0.001305   0.000347
    C2_ s       0.006734   0.000430   0.000423   0.000022  -0.000270   0.000071
    C2_ p       0.000693   0.002954   0.005445   0.003961   0.000796  -0.000558
    C2_ d      -0.000484   0.000215  -0.002145   0.000871   0.002099   0.001658
    H1_ s       0.000641   0.002976   0.000859  -0.002221  -0.000069   0.000045
    H1_ p       0.000192  -0.000061   0.000109   0.000020   0.000116  -0.000011
    H2_ s       0.001246   0.001506   0.001734  -0.001111  -0.000033   0.000022
    H2_ p       0.000379  -0.000002   0.000216   0.000386   0.000042   0.000485
 
   ao class      13b3u      14b3u      15b3u      16b3u      17b3u      18b3u
    C1_ s       0.000302   0.000495  -0.000221  -0.000004   0.000173   0.000048
    C1_ p       0.001453   0.000606  -0.000010  -0.000133   0.000055   0.000127
    C1_ d      -0.000573   0.000297   0.000090   0.000363   0.000126  -0.000027
    C2_ s       0.000601   0.000253  -0.000115   0.000006   0.000075   0.000116
    C2_ p       0.002911  -0.000062   0.000733  -0.000262  -0.000187   0.000362
    C2_ d      -0.001143   0.000115   0.000260   0.000728   0.000131  -0.000017
    H1_ s      -0.000855  -0.000455   0.000017   0.000044   0.000044  -0.000061
    H1_ p       0.000340   0.000113   0.000005  -0.000025  -0.000024   0.000011
    H2_ s      -0.001711  -0.000230   0.000010   0.000085   0.000019  -0.000117
    H2_ p       0.000676   0.000121   0.000036  -0.000048   0.000032  -0.000000
 
   ao class      19b3u      20b3u      21b3u      22b3u      23b3u      24b3u
    C1_ s      -0.000195  -0.000055  -0.000298  -0.000037   0.000006  -0.000064
    C1_ p       0.000045   0.000093   0.000013   0.000010   0.000004   0.000029
    C1_ d       0.000016   0.000021   0.000155   0.000001   0.000062   0.000057
    C2_ s      -0.000098  -0.000094  -0.000155  -0.000018   0.000003  -0.000132
    C2_ p       0.000416   0.000181   0.000490   0.000311  -0.000019   0.000057
    C2_ d       0.000046   0.000043  -0.000158  -0.000247   0.000035   0.000114
    H1_ s      -0.000031   0.000057   0.000083   0.000036  -0.000011   0.000017
    H1_ p       0.000014  -0.000048  -0.000013   0.000004   0.000016  -0.000016
    H2_ s      -0.000015   0.000114   0.000042   0.000018  -0.000006   0.000033
    H2_ p       0.000127  -0.000095   0.000012   0.000053   0.000007  -0.000033
 
   ao class      25b3u      26b3u      27b3u      28b3u      29b3u      30b3u
    C1_ s       0.000011   0.000009   0.000016   0.000007  -0.000007  -0.000008
    C1_ p      -0.000004  -0.000011  -0.000036  -0.000002   0.000004   0.000006
    C1_ d       0.000006  -0.000010   0.000008   0.000000  -0.000001   0.000000
    C2_ s       0.000018   0.000004   0.000010   0.000014  -0.000004  -0.000005
    C2_ p      -0.000007  -0.000026  -0.000100  -0.000004   0.000015   0.000020
    C2_ d       0.000013   0.000076   0.000094   0.000000  -0.000006  -0.000006
    H1_ s      -0.000012   0.000009   0.000037  -0.000003  -0.000003  -0.000005
    H1_ p       0.000013   0.000001   0.000003   0.000004   0.000012  -0.000001
    H2_ s      -0.000025   0.000004   0.000019  -0.000005  -0.000002  -0.000004
    H2_ p       0.000027  -0.000024  -0.000022   0.000007   0.000004   0.000009
 
   ao class      31b3u      32b3u      33b3u      34b3u      35b3u      36b3u
    C1_ s       0.000001   0.000030  -0.000004  -0.000000   0.000002  -0.000000
    C1_ p      -0.000006  -0.000000   0.000001   0.000002  -0.000002   0.000001
    C1_ d       0.000004  -0.000001  -0.000001  -0.000000  -0.000000  -0.000000
    C2_ s       0.000004   0.000012   0.000001  -0.000000   0.000003  -0.000001
    C2_ p      -0.000016  -0.000039  -0.000000   0.000002  -0.000003   0.000001
    C2_ d       0.000006  -0.000001  -0.000001  -0.000001  -0.000000  -0.000000
    H1_ s       0.000005   0.000001   0.000000  -0.000001   0.000000  -0.000000
    H1_ p      -0.000001   0.000000   0.000002   0.000000   0.000000  -0.000000
    H2_ s       0.000012   0.000001   0.000001  -0.000001   0.000000  -0.000000
    H2_ p      -0.000002   0.000000   0.000004   0.000001   0.000000  -0.000000
 
   ao class      37b3u      38b3u      39b3u
    C1_ s      -0.000000  -0.000000   0.000000
    C1_ p       0.000000  -0.000000  -0.000000
    C1_ d      -0.000000  -0.000000   0.000000
    C2_ s      -0.000000  -0.000000   0.000000
    C2_ p       0.000000   0.000000  -0.000000
    C2_ d       0.000000  -0.000000   0.000000
    H1_ s       0.000000   0.000000   0.000000
    H1_ p       0.000000   0.000000  -0.000000
    H2_ s       0.000000   0.000000   0.000000
    H2_ p      -0.000000   0.000000  -0.000000

                        b2u partial gross atomic populations
   ao class       1b2u       2b2u       3b2u       4b2u       5b2u       6b2u
    C1_ p       0.016768   0.146466   0.126917   0.650202   0.005176   0.001338
    C1_ d       0.005055  -0.036067  -0.020659   0.018178   0.001190   0.000107
    C2_ s       1.965412   1.101473   0.117849  -0.000001  -0.000142   0.001266
    C2_ p       0.007796   0.012494   1.160378   1.314629   0.004213   0.003465
    C2_ d       0.001844   0.045952  -0.037507   0.035779   0.002201   0.000240
    H1_ p       0.000103  -0.010189   0.002501  -0.014628   0.000075   0.000018
    H2_ s       0.001851   0.669198   0.650416  -0.000000   0.001397   0.004485
    H2_ p       0.001171   0.055577  -0.023191  -0.029352  -0.000023  -0.000082
 
   ao class       7b2u       8b2u       9b2u      10b2u      11b2u      12b2u
    C1_ p       0.000819   0.002020   0.000591  -0.000340  -0.000248   0.000341
    C1_ d       0.000875   0.000272   0.000968   0.000991  -0.000021  -0.000025
    C2_ s       0.000000   0.000092  -0.000830   0.000213   0.000753   0.000000
    C2_ p       0.001636   0.003736   0.000082  -0.000308   0.000789   0.000694
    C2_ d       0.001760   0.001503   0.002444   0.001016   0.000432  -0.000054
    H1_ p       0.000475   0.000251  -0.000010   0.000328   0.000043   0.000052
    H2_ s       0.000000  -0.003333  -0.000099   0.000067  -0.000688  -0.000000
    H2_ p       0.000950   0.000155   0.000168   0.000146   0.000192   0.000106
 
   ao class      13b2u      14b2u      15b2u      16b2u      17b2u      18b2u
    C1_ p       0.000488   0.000020  -0.000108   0.000262  -0.000021   0.000318
    C1_ d       0.000145   0.000125   0.000058   0.000025   0.000148  -0.000162
    C2_ s      -0.000333   0.000000   0.000247  -0.000293  -0.000000  -0.000444
    C2_ p       0.000231   0.000041  -0.000028   0.000200  -0.000040   0.000180
    C2_ d       0.000205   0.000251   0.000201   0.000036   0.000295   0.000155
    H1_ p       0.000022   0.000082   0.000023   0.000080  -0.000058   0.000013
    H2_ s       0.000027  -0.000000   0.000065  -0.000046   0.000000   0.000124
    H2_ p       0.000018   0.000164  -0.000015   0.000060  -0.000115  -0.000013
 
   ao class      19b2u      20b2u      21b2u      22b2u      23b2u      24b2u
    C1_ p       0.000210  -0.000015  -0.000014  -0.000054  -0.000001   0.000008
    C1_ d      -0.000165   0.000003   0.000054   0.000060  -0.000007  -0.000003
    C2_ s      -0.000064   0.000010   0.000013   0.000026   0.000000  -0.000011
    C2_ p       0.000119  -0.000002  -0.000023  -0.000082  -0.000003   0.000010
    C2_ d      -0.000081   0.000095   0.000012   0.000042  -0.000013  -0.000004
    H1_ p       0.000034  -0.000001  -0.000016  -0.000015   0.000015  -0.000001
    H2_ s       0.000056  -0.000018   0.000013   0.000057   0.000000  -0.000005
    H2_ p       0.000022   0.000024  -0.000007  -0.000004   0.000030   0.000017
 
   ao class      25b2u      26b2u      27b2u      28b2u      29b2u      30b2u
    C1_ p       0.000010  -0.000026   0.000001   0.000001   0.000000   0.000000
    C1_ d      -0.000005  -0.000000  -0.000000  -0.000001   0.000000  -0.000000
    C2_ s      -0.000013   0.000042  -0.000000  -0.000001  -0.000000  -0.000000
    C2_ p       0.000015  -0.000013   0.000003   0.000003   0.000000  -0.000000
    C2_ d      -0.000002  -0.000002  -0.000000  -0.000000   0.000000  -0.000000
    H1_ p       0.000006   0.000000  -0.000000   0.000001  -0.000000   0.000000
    H2_ s      -0.000008   0.000002  -0.000000  -0.000002   0.000000   0.000000
    H2_ p       0.000002   0.000000  -0.000000   0.000000   0.000000   0.000000

                        b1g partial gross atomic populations
   ao class       1b1g       2b1g       3b1g       4b1g       5b1g       6b1g
    C1_ p       0.003407   0.035210   0.760229   0.002351   0.000324   0.002644
    C1_ d       0.000348   0.002085   0.015479   0.000824   0.000341   0.000502
    C2_ s       1.989749   0.720334   0.087424   0.006468   0.000712   0.000000
    C2_ p       0.004512   0.205105   0.792288   0.003507   0.005997   0.005266
    C2_ d       0.001538  -0.167352  -0.032364   0.001032  -0.001585   0.001020
    H1_ p       0.000044  -0.009301  -0.007647   0.000057  -0.000006   0.000061
    H2_ s       0.000202   1.152810   0.417458   0.000261   0.004290   0.000000
    H2_ p       0.000200   0.042593  -0.058965  -0.000048   0.000256   0.000122
 
   ao class       7b1g       8b1g       9b1g      10b1g      11b1g      12b1g
    C1_ p       0.001370  -0.000019  -0.000264   0.000246  -0.000025   0.000472
    C1_ d       0.000453   0.000640   0.000695   0.000195   0.000025  -0.000131
    C2_ s       0.000474  -0.000514   0.000000  -0.000205   0.000114  -0.000142
    C2_ p       0.001233   0.001972  -0.000528  -0.000069   0.000344   0.000268
    C2_ d       0.000496   0.001357   0.001385   0.001223   0.000310   0.000170
    H1_ p       0.000554   0.000037   0.000160   0.000024   0.000086   0.000066
    H2_ s      -0.000313  -0.001610  -0.000000  -0.000328  -0.000051   0.000000
    H2_ p       0.000369   0.000469   0.000319   0.000402   0.000118   0.000073
 
   ao class      13b1g      14b1g      15b1g      16b1g      17b1g      18b1g
    C1_ p      -0.000007  -0.000045   0.000155  -0.000056   0.000071   0.000003
    C1_ d       0.000055   0.000061  -0.000146   0.000104   0.000046  -0.000021
    C2_ s       0.000110   0.000068   0.000000  -0.000020  -0.000244   0.000049
    C2_ p       0.000230   0.000107   0.000308  -0.000023   0.000151   0.000029
    C2_ d       0.000037   0.000120  -0.000294   0.000149   0.000082   0.000022
    H1_ p       0.000095   0.000011   0.000065  -0.000023  -0.000015   0.000005
    H2_ s      -0.000052   0.000138   0.000000  -0.000002   0.000023  -0.000031
    H2_ p       0.000001  -0.000162   0.000130   0.000009  -0.000013   0.000030
 
   ao class      19b1g      20b1g      21b1g      22b1g      23b1g      24b1g
    C1_ p       0.000184  -0.000008   0.000010  -0.000068  -0.000001  -0.000004
    C1_ d      -0.000187   0.000040  -0.000001   0.000094  -0.000000  -0.000008
    C2_ s       0.000000   0.000057  -0.000001  -0.000000   0.000001   0.000004
    C2_ p       0.000359  -0.000043  -0.000020  -0.000130  -0.000008  -0.000004
    C2_ d      -0.000365   0.000028   0.000003   0.000185   0.000000  -0.000001
    H1_ p       0.000025  -0.000029   0.000003  -0.000021   0.000002   0.000009
    H2_ s      -0.000000   0.000017   0.000044   0.000000  -0.000003   0.000003
    H2_ p       0.000050  -0.000012  -0.000018  -0.000041   0.000017   0.000007
 
   ao class      25b1g      26b1g      27b1g      28b1g      29b1g      30b1g
    C1_ p      -0.000002  -0.000001   0.000001   0.000000   0.000000   0.000000
    C1_ d       0.000000  -0.000000  -0.000001   0.000000   0.000000   0.000000
    C2_ s       0.000017  -0.000002  -0.000000  -0.000001   0.000000   0.000000
    C2_ p      -0.000014   0.000007   0.000001   0.000001  -0.000000   0.000000
    C2_ d      -0.000000  -0.000000  -0.000003  -0.000000   0.000000   0.000000
    H1_ p      -0.000000  -0.000000   0.000001  -0.000000   0.000000  -0.000000
    H2_ s       0.000003  -0.000003  -0.000000  -0.000000   0.000000   0.000000
    H2_ p       0.000000   0.000000   0.000002  -0.000000   0.000000  -0.000000

                        b1u partial gross atomic populations
   ao class       1b1u       2b1u       3b1u       4b1u       5b1u       6b1u
    C1_ p       0.632639   0.326780   0.000160   0.000718   0.000721   0.001375
    C1_ d       0.010507  -0.001392   0.001393   0.000359   0.000911   0.000019
    C2_ p       1.266700   0.164970   0.000527   0.003395   0.000098   0.000448
    C2_ d       0.021634   0.023082   0.002901   0.000495   0.001716   0.001682
    H1_ p       0.005318   0.008769   0.000152   0.000040   0.000464   0.000017
    H2_ p       0.010550   0.004478   0.000328   0.000119   0.000207   0.000005
 
   ao class       7b1u       8b1u       9b1u      10b1u      11b1u      12b1u
    C1_ p       0.000228   0.000008  -0.000025   0.000005   0.000001   0.000089
    C1_ d       0.000338   0.000017   0.000039   0.000023   0.000141   0.000086
    C2_ p       0.000106   0.000071   0.000247   0.000225   0.000006   0.000008
    C2_ d       0.000490   0.000043   0.000049   0.000064   0.000306   0.000053
    H1_ p       0.000129   0.000188   0.000148   0.000054  -0.000034  -0.000019
    H2_ p       0.000070   0.000377   0.000023   0.000102  -0.000072  -0.000008
 
   ao class      13b1u      14b1u      15b1u      16b1u
    C1_ p       0.000017  -0.000001   0.000015   0.000002
    C1_ d       0.000023  -0.000009   0.000019  -0.000002
    C2_ p      -0.000004   0.000003  -0.000004  -0.000000
    C2_ d       0.000086   0.000002   0.000001  -0.000001
    H1_ p       0.000001   0.000021  -0.000003   0.000006
    H2_ p      -0.000000   0.000029   0.000008   0.000003

                        b2g partial gross atomic populations
   ao class       1b2g       2b2g       3b2g       4b2g       5b2g       6b2g
    C1_ p       0.951225   0.008953   0.000931   0.001406   0.000189   0.000023
    C1_ d       0.001359  -0.000394   0.001794   0.000540   0.000269   0.000042
    C2_ p       0.472373   0.017983   0.000457   0.000679   0.000402   0.000012
    C2_ d       0.023390  -0.000807   0.001086   0.000274   0.000521   0.001061
    H1_ p       0.012654   0.000224   0.000477   0.000056   0.000299   0.000005
    H2_ p       0.006411   0.000464   0.000234   0.000031   0.000597   0.000002
 
   ao class       7b2g       8b2g       9b2g      10b2g      11b2g      12b2g
    C1_ p      -0.000002   0.000247   0.000028   0.000024   0.000005  -0.000024
    C1_ d       0.000051  -0.000132   0.000024   0.000155   0.000006   0.000059
    C2_ p      -0.000001   0.000493   0.000054   0.000014   0.000003  -0.000041
    C2_ d       0.000068  -0.000266   0.000043   0.000079   0.000048   0.000116
    H1_ p       0.000327   0.000037   0.000051  -0.000034  -0.000000  -0.000018
    H2_ p       0.000163   0.000076   0.000105  -0.000018   0.000000  -0.000035
 
   ao class      13b2g      14b2g      15b2g      16b2g
    C1_ p       0.000021   0.000002  -0.000000   0.000000
    C1_ d       0.000002  -0.000002  -0.000001   0.000000
    C2_ p       0.000007   0.000001  -0.000000   0.000001
    C2_ d       0.000015  -0.000001  -0.000002   0.000000
    H1_ p       0.000007   0.000008   0.000003  -0.000000
    H2_ p       0.000003   0.000004   0.000006  -0.000000

                        b3g partial gross atomic populations
   ao class       1b3g       2b3g       3b3g       4b3g       5b3g       6b3g
    C1_ d       0.013882   0.001812   0.000125  -0.000005   0.000691   0.000028
    C2_ p       1.410679   0.000000   0.001326   0.002088   0.000035  -0.000003
    C2_ d       0.009231   0.003610   0.002784   0.000778   0.000411   0.000091
    H2_ p       0.018565   0.000000   0.000714   0.000084   0.000007   0.000489
 
   ao class       7b3g       8b3g       9b3g      10b3g      11b3g
    C1_ d       0.000106  -0.000002   0.000033   0.000006  -0.000000
    C2_ p      -0.000000   0.000037   0.000005   0.000029   0.000003
    C2_ d       0.000208   0.000236   0.000023   0.000007  -0.000003
    H2_ p      -0.000000  -0.000052   0.000000   0.000010   0.000012

                        au  partial gross atomic populations
   ao class       1au        2au        3au        4au        5au        6au 
    C1_ d       0.018059   0.000717   0.001217   0.000230   0.000207   0.000052
    C2_ p       0.502407   0.001014   0.001703   0.000315  -0.000000   0.000003
    C2_ d       0.006465   0.001749   0.000656   0.000630   0.000415   0.000151
    H2_ p       0.014150   0.000682   0.000013   0.000207   0.000000   0.000271
 
   ao class       7au        8au        9au       10au       11au 
    C1_ d       0.000014   0.000056   0.000000   0.000012  -0.000000
    C2_ p       0.000115   0.000012   0.000012   0.000000   0.000001
    C2_ d       0.000113   0.000055   0.000029   0.000024  -0.000003
    H2_ p      -0.000028   0.000002  -0.000002  -0.000000   0.000010


                        gross atomic populations
     ao           C1_        C2_        H1_        H2_
      s         5.898941  11.875306   2.846629   5.679222
      p         5.391507  10.853825   0.017945   0.038875
      d        -0.201958  -0.400292   0.000000   0.000000
    total      11.088491  22.328838   2.864573   5.718097
 

 Total number of electrons:   42.00000000

 item #                     2 suffix=:.drt1.state2:
 read_civout: repnuc=  -224.450660864067     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 91% root-following 0
 MR-CISD energy:  -231.27775672    -6.82709585
 residuum:     0.00014989
 deltae:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  1max overlap with ref# 91% root-following 0
 MR-CISD energy:  -231.25257990    -6.80191904
 residuum:     0.00079436
 deltae:     0.00000064
 apxde:     0.00000026

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.90774063    -0.01139169     0.00522583    -0.01669692    -0.09979526     0.11087837    -0.19898672     0.00000000
 ref:   2    -0.01111154    -0.90715837    -0.03365063     0.00172177    -0.01193882     0.01959314     0.01021576     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.90774063    -0.01139169     0.00522583    -0.01669692    -0.09979526     0.00000000     0.00000000     0.00000000
 ref:   2    -0.01111154    -0.90715837    -0.03365063     0.00172177    -0.01193882     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 computing MRCISD density
 densi: densityinfo%a4den=   1.00000000000000     
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    2
--------------------------------------------------------------------------------
================================================================================
   DYZ=   80766  DYX=  101959  DYW=  106443
   D0Z=   19955  D0Y=  162182  D0X=   19509  D0W=   21336
  DDZI=   40986 DDYI=  280518 DDXI=   35816 DDWI=   38640
  DDZE=       0 DDYE=   41898 DDXE=    6259 DDWE=    6729
================================================================================
Trace of MO density:    30.000000
   30  correlated and    12  frozen core electrons

          modens reordered block   1

               ag    1        ag    2        ag    3        ag    4        ag    5        ag    6        ag    7        ag    8
  ag    1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    3    0.00000        0.00000        1.98725      -5.333264E-05   1.163792E-03  -7.496124E-05   1.441526E-03  -1.129008E-04
  ag    4    0.00000        0.00000      -5.333264E-05    1.98063       2.120165E-05  -2.311001E-03   3.162614E-04  -7.505878E-04
  ag    5    0.00000        0.00000       1.163792E-03   2.120165E-05    1.97942       1.183475E-04   1.155298E-03  -1.420899E-04
  ag    6    0.00000        0.00000      -7.496124E-05  -2.311001E-03   1.183475E-04    1.97426      -2.276588E-05   1.565626E-03
  ag    7    0.00000        0.00000       1.441526E-03   3.162614E-04   1.155298E-03  -2.276588E-05   6.189923E-04  -5.876573E-06
  ag    8    0.00000        0.00000      -1.129008E-04  -7.505878E-04  -1.420899E-04   1.565626E-03  -5.876573E-06   2.246893E-04
  ag    9    0.00000        0.00000       9.478491E-05   5.158548E-05   7.310511E-04  -1.729895E-04   5.955775E-04   5.998183E-06
  ag   10    0.00000        0.00000       3.401681E-04   1.970539E-03   4.050427E-05  -8.231461E-04   3.072157E-05  -4.037970E-04
  ag   11    0.00000        0.00000      -3.269389E-03  -6.121810E-04  -1.297059E-04  -1.897497E-04  -9.925363E-04   3.448587E-06
  ag   12    0.00000        0.00000       4.275903E-05  -3.768964E-04   3.900995E-04  -3.432601E-03   3.143683E-07  -3.658298E-04
  ag   13    0.00000        0.00000       4.184870E-04   2.720683E-03  -2.283797E-04   1.647725E-03   1.359245E-05  -1.865926E-04
  ag   14    0.00000        0.00000      -3.288767E-03  -2.512143E-04   4.144167E-03  -4.426522E-04  -1.226736E-03   1.387070E-06
  ag   15    0.00000        0.00000       4.275922E-04   1.018545E-03  -5.270834E-04   8.355843E-04  -1.012257E-05  -6.034579E-04
  ag   16    0.00000        0.00000      -1.424704E-04  -1.852577E-03   1.038877E-03  -5.492925E-03   5.337518E-07  -5.919999E-04
  ag   17    0.00000        0.00000      -3.216082E-03  -8.422834E-04   1.275778E-03  -5.460851E-04  -7.485269E-04   6.756832E-06
  ag   18    0.00000        0.00000      -5.169914E-04   2.337363E-03   1.129268E-03  -1.260833E-03   1.602335E-05   3.633106E-04
  ag   19    0.00000        0.00000      -2.510483E-04   4.975285E-05   1.153698E-03  -4.803996E-04   3.895475E-04  -7.521477E-07
  ag   20    0.00000        0.00000       3.578469E-03   5.888521E-05  -7.695305E-03   1.054040E-03   4.136198E-04   4.828983E-06
  ag   21    0.00000        0.00000       5.594317E-04   2.847564E-03  -2.753182E-04   1.580507E-03   2.315575E-05  -1.611339E-04
  ag   22    0.00000        0.00000       2.404299E-04   2.361729E-03  -6.033114E-04   1.106646E-03   6.722261E-06   3.828263E-04
  ag   23    0.00000        0.00000      -1.965827E-04   1.922169E-04  -4.204702E-04   9.901999E-04  -1.595223E-05   5.060678E-04
  ag   24    0.00000        0.00000      -5.198951E-04   2.641841E-03   5.156656E-04   8.322894E-04   1.023917E-05   9.756191E-05
  ag   25    0.00000        0.00000      -1.750149E-03  -4.700418E-05   1.753324E-03  -5.867448E-04   7.156471E-04   3.270872E-06
  ag   26    0.00000        0.00000       4.973958E-05   1.221030E-03  -1.811472E-04   3.325092E-05  -7.728217E-06   1.283446E-04
  ag   27    0.00000        0.00000      -1.359521E-03  -3.773294E-04   1.898040E-03  -4.076131E-04  -4.877350E-04  -7.947046E-08
  ag   28    0.00000        0.00000       3.710399E-04   2.188452E-03  -1.383506E-04  -1.168156E-03   1.212430E-04  -8.099373E-07
  ag   29    0.00000        0.00000      -6.752385E-04   2.307025E-04  -2.358409E-04  -4.036209E-04  -3.987300E-04   4.229295E-06
  ag   30    0.00000        0.00000      -4.887150E-05   2.743987E-04  -1.427707E-05   2.529289E-03  -1.411584E-06  -1.232611E-04
  ag   31    0.00000        0.00000       1.744627E-03   1.304938E-04   3.351150E-04   3.572251E-04  -1.596107E-04  -1.120384E-06
  ag   32    0.00000        0.00000      -1.467284E-03  -9.887123E-05   4.502548E-03  -5.160468E-04   1.659653E-04  -1.409301E-06
  ag   33    0.00000        0.00000      -1.789839E-05  -2.982552E-03   2.614068E-04  -4.090752E-03  -3.414198E-06   1.467077E-04
  ag   34    0.00000        0.00000      -3.376652E-04  -1.101412E-03   5.384105E-05   3.965589E-04  -9.957534E-06   1.515289E-04
  ag   35    0.00000        0.00000      -1.135619E-04   2.334499E-03  -8.669806E-05  -1.132407E-03  -2.271715E-06  -1.400010E-06
  ag   36    0.00000        0.00000       2.919074E-05   9.409949E-04   3.221830E-04  -3.959408E-03  -1.426912E-05   7.412974E-06
  ag   37    0.00000        0.00000      -9.866949E-04   5.632677E-05  -2.713102E-03  -5.507022E-04   1.504299E-04   1.210812E-06
  ag   38    0.00000        0.00000      -6.169663E-05  -3.064375E-03   1.815423E-04  -3.110400E-03  -1.128000E-06  -1.038322E-04
  ag   39    0.00000        0.00000      -3.167306E-05   5.845539E-04  -1.081949E-05   6.342194E-04   1.548386E-07   1.021115E-05

               ag    9        ag   10        ag   11        ag   12        ag   13        ag   14        ag   15        ag   16
  ag    3   9.478491E-05   3.401681E-04  -3.269389E-03   4.275903E-05   4.184870E-04  -3.288767E-03   4.275922E-04  -1.424704E-04
  ag    4   5.158548E-05   1.970539E-03  -6.121810E-04  -3.768964E-04   2.720683E-03  -2.512143E-04   1.018545E-03  -1.852577E-03
  ag    5   7.310511E-04   4.050427E-05  -1.297059E-04   3.900995E-04  -2.283797E-04   4.144167E-03  -5.270834E-04   1.038877E-03
  ag    6  -1.729895E-04  -8.231461E-04  -1.897497E-04  -3.432601E-03   1.647725E-03  -4.426522E-04   8.355843E-04  -5.492925E-03
  ag    7   5.955775E-04   3.072157E-05  -9.925363E-04   3.143683E-07   1.359245E-05  -1.226736E-03  -1.012257E-05   5.337518E-07
  ag    8   5.998183E-06  -4.037970E-04   3.448587E-06  -3.658298E-04  -1.865926E-04   1.387070E-06  -6.034579E-04  -5.919999E-04
  ag    9   2.103328E-03   2.172401E-05  -4.754005E-04  -5.158380E-06  -2.864670E-05  -2.926443E-03  -1.448848E-06  -1.033011E-05
  ag   10   2.172401E-05   8.663475E-04  -3.327265E-05   5.359028E-04   5.685142E-04  -6.395494E-05   1.590752E-03   6.273128E-04
  ag   11  -4.754005E-04  -3.327265E-05   1.819546E-03   8.881372E-06  -1.454549E-05   1.723754E-03   4.524402E-05   7.521015E-06
  ag   12  -5.158380E-06   5.359028E-04   8.881372E-06   7.587657E-04   2.254809E-05   3.079331E-06   6.511650E-04   1.557606E-03
  ag   13  -2.864670E-05   5.685142E-04  -1.454549E-05   2.254809E-05   7.462101E-04   1.125560E-05   1.093236E-03  -5.299329E-04
  ag   14  -2.926443E-03  -6.395494E-05   1.723754E-03   3.079331E-06   1.125560E-05   5.811516E-03   5.297132E-06   1.656117E-05
  ag   15  -1.448848E-06   1.590752E-03   4.524402E-05   6.511650E-04   1.093236E-03   5.297132E-06   4.033825E-03   1.991624E-04
  ag   16  -1.033011E-05   6.273128E-04   7.521015E-06   1.557606E-03  -5.299329E-04   1.656117E-05   1.991624E-04   3.949435E-03
  ag   17   6.505686E-04  -2.993809E-05   1.813407E-03   6.071743E-06  -1.929168E-05   6.188807E-04   1.213956E-05   4.780972E-06
  ag   18  -3.080807E-05  -1.328812E-03  -5.200968E-05  -2.157428E-04  -9.941394E-04   4.122695E-05  -4.793854E-03   1.038276E-03
  ag   19   9.088172E-04   1.393626E-05  -3.124497E-04  -1.071544E-07  -3.940305E-06  -6.447647E-04  -1.442315E-05   5.171192E-06
  ag   20   1.454092E-03   1.990900E-05  -8.158043E-04  -5.428347E-06  -1.714128E-05  -4.213017E-03  -2.619150E-05  -1.627110E-05
  ag   21  -3.389744E-05   5.225318E-04  -3.232094E-05  -1.049136E-04   9.631191E-04  -3.734408E-06   7.655173E-04  -9.946540E-04
  ag   22  -1.169659E-05  -4.513024E-04  -2.052165E-05  -1.056111E-03   4.669026E-04  -1.235123E-06  -4.874750E-04  -2.811986E-03
  ag   23   9.155327E-06  -8.642376E-04   1.173084E-05  -1.064757E-03  -1.628541E-04   1.072454E-05  -1.138907E-03  -2.593501E-03
  ag   24   1.204196E-05  -4.165715E-04  -1.689571E-05  -9.298104E-05  -3.222753E-04   3.936858E-06  -2.154130E-03   6.317180E-04
  ag   25   1.686752E-03   2.915228E-05  -6.278677E-04  -1.515760E-05   6.054415E-06  -1.608098E-03   1.767655E-05  -4.697630E-05
  ag   26   4.559749E-06  -1.917608E-04   1.201290E-05  -2.848870E-04  -5.120038E-05   5.556821E-06   5.065665E-05  -9.036570E-04
  ag   27  -4.297555E-04  -2.207148E-05   9.982980E-04   7.586419E-06   6.257295E-07   1.793537E-03  -3.257843E-06   1.754612E-05
  ag   28  -1.633028E-04   9.380215E-05  -2.718474E-04  -2.680339E-04   5.686521E-04   1.440747E-04  -1.621724E-04  -1.101033E-03
  ag   29   5.496896E-04   9.949984E-06   9.273946E-04  -6.948950E-05   1.452211E-04  -5.621660E-04  -5.092174E-05  -2.971512E-04
  ag   30   2.807572E-06   1.254243E-04   5.646504E-06   2.345847E-04   4.779061E-05  -6.752852E-06  -1.063277E-04   3.138187E-04
  ag   31  -7.007645E-04  -4.838179E-06   8.231387E-05   1.712006E-06   2.077304E-06   9.132525E-04   1.911045E-05   7.038060E-06
  ag   32   1.126544E-04   3.528455E-06  -1.916225E-04   1.750112E-06  -7.204776E-07   1.189233E-04  -1.075474E-05   8.101009E-06
  ag   33   7.471365E-06  -4.284821E-04  -1.166284E-06  -8.165668E-05  -4.212273E-04   3.417806E-06  -7.606738E-04   1.164531E-05
  ag   34   1.313000E-05  -3.077031E-04   1.375320E-05  -2.969348E-04  -2.466694E-04   1.131054E-05  -5.858264E-04  -5.226105E-04
  ag   35   2.959585E-06  -3.525399E-05   3.900443E-06   8.849673E-06  -4.701471E-05   4.307491E-06  -1.698727E-04   4.734356E-05
  ag   36  -2.925615E-05  -8.371400E-05   2.020859E-05   1.192707E-05  -9.098438E-05   6.424141E-05  -3.131799E-04   7.697535E-06
  ag   37   3.395149E-04  -1.549628E-06  -2.153625E-04  -1.600367E-06  -8.537309E-06  -6.764382E-04  -3.505869E-05  -6.527135E-06
  ag   38  -3.968893E-06   1.811680E-04   5.300698E-06   2.170235E-04   2.940321E-05   7.185816E-06   3.170449E-04   4.726741E-04
  ag   39  -3.439979E-07  -4.544438E-05  -1.453092E-06   4.244514E-05  -8.455237E-05  -1.175799E-06  -1.177062E-04   2.399875E-04

               ag   17        ag   18        ag   19        ag   20        ag   21        ag   22        ag   23        ag   24
  ag    3  -3.216082E-03  -5.169914E-04  -2.510483E-04   3.578469E-03   5.594317E-04   2.404299E-04  -1.965827E-04  -5.198951E-04
  ag    4  -8.422834E-04   2.337363E-03   4.975285E-05   5.888521E-05   2.847564E-03   2.361729E-03   1.922169E-04   2.641841E-03
  ag    5   1.275778E-03   1.129268E-03   1.153698E-03  -7.695305E-03  -2.753182E-04  -6.033114E-04  -4.204702E-04   5.156656E-04
  ag    6  -5.460851E-04  -1.260833E-03  -4.803996E-04   1.054040E-03   1.580507E-03   1.106646E-03   9.901999E-04   8.322894E-04
  ag    7  -7.485269E-04   1.602335E-05   3.895475E-04   4.136198E-04   2.315575E-05   6.722261E-06  -1.595223E-05   1.023917E-05
  ag    8   6.756832E-06   3.633106E-04  -7.521477E-07   4.828983E-06  -1.611339E-04   3.828263E-04   5.060678E-04   9.756191E-05
  ag    9   6.505686E-04  -3.080807E-05   9.088172E-04   1.454092E-03  -3.389744E-05  -1.169659E-05   9.155327E-06   1.204196E-05
  ag   10  -2.993809E-05  -1.328812E-03   1.393626E-05   1.990900E-05   5.225318E-04  -4.513024E-04  -8.642376E-04  -4.165715E-04
  ag   11   1.813407E-03  -5.200968E-05  -3.124497E-04  -8.158043E-04  -3.232094E-05  -2.052165E-05   1.173084E-05  -1.689571E-05
  ag   12   6.071743E-06  -2.157428E-04  -1.071544E-07  -5.428347E-06  -1.049136E-04  -1.056111E-03  -1.064757E-03  -9.298104E-05
  ag   13  -1.929168E-05  -9.941394E-04  -3.940305E-06  -1.714128E-05   9.631191E-04   4.669026E-04  -1.628541E-04  -3.222753E-04
  ag   14   6.188807E-04   4.122695E-05  -6.447647E-04  -4.213017E-03  -3.734408E-06  -1.235123E-06   1.072454E-05   3.936858E-06
  ag   15   1.213956E-05  -4.793854E-03  -1.442315E-05  -2.619150E-05   7.655173E-04  -4.874750E-04  -1.138907E-03  -2.154130E-03
  ag   16   4.780972E-06   1.038276E-03   5.171192E-06  -1.627110E-05  -9.946540E-04  -2.811986E-03  -2.593501E-03   6.317180E-04
  ag   17   2.651305E-03  -2.010290E-05   4.389885E-04  -7.267030E-04  -3.171102E-05  -1.764612E-05   9.554195E-06   5.597855E-06
  ag   18  -2.010290E-05   7.546622E-03   2.777552E-05  -1.504207E-05  -5.542517E-04  -1.849364E-04   8.823342E-05   4.568973E-03
  ag   19   4.389885E-04   2.777552E-05   1.088788E-03  -5.451824E-04  -5.137999E-06  -4.694143E-06  -5.429246E-06   3.558017E-05
  ag   20  -7.267030E-04  -1.504207E-05  -5.451824E-04   4.831141E-03  -1.755718E-06   7.469400E-06   1.044730E-05  -1.855738E-05
  ag   21  -3.171102E-05  -5.542517E-04  -5.137999E-06  -1.755718E-06   1.536996E-03   1.047035E-03   2.163136E-04  -2.617384E-04
  ag   22  -1.764612E-05  -1.849364E-04  -4.694143E-06   7.469400E-06   1.047035E-03   2.282347E-03   2.018857E-03  -3.313627E-04
  ag   23   9.554195E-06   8.823342E-05  -5.429246E-06   1.044730E-05   2.163136E-04   2.018857E-03   2.679432E-03  -4.692992E-04
  ag   24   5.597855E-06   4.568973E-03   3.558017E-05  -1.855738E-05  -2.617384E-04  -3.313627E-04  -4.692992E-04   3.951042E-03
  ag   25   5.701119E-04  -6.125523E-05   1.190615E-03  -8.194869E-04   9.489786E-06   2.382887E-05   1.345073E-05  -1.401436E-05
  ag   26   1.490652E-05  -8.570411E-04  -5.048902E-06  -4.067048E-06   2.552807E-05   6.692779E-04   9.925890E-04  -9.717304E-04
  ag   27   1.231981E-03   2.770613E-05   3.404748E-04  -1.732162E-03  -3.444701E-06  -8.680336E-06  -3.944706E-06   1.651229E-05
  ag   28  -3.846253E-04   5.023078E-04   4.232654E-05  -2.281307E-04   1.187872E-03   1.258345E-03   7.337419E-04   3.460649E-04
  ag   29   1.360354E-03   1.342130E-04  -1.587932E-04   8.575716E-04   3.128583E-04   3.335605E-04   1.964432E-04   9.577160E-05
  ag   30   6.478551E-06   3.216715E-04  -4.936370E-07   1.267368E-05   2.571397E-04  -3.630540E-05   2.163241E-04   1.364351E-05
  ag   31  -4.855704E-04  -1.786533E-05  -7.217664E-04  -3.897745E-04  -4.899108E-06  -9.672773E-06  -8.307404E-06  -1.034576E-05
  ag   32   5.681528E-05   2.200851E-05   3.187139E-04  -1.204793E-03  -2.968369E-07  -4.733289E-06  -8.849006E-06   2.114477E-05
  ag   33   3.588337E-06   6.560730E-05  -5.165098E-07  -1.519894E-07  -4.141827E-04   8.976337E-05   6.656007E-04  -9.644439E-04
  ag   34   1.165933E-05   8.251291E-04   7.695955E-06  -9.004277E-06  -3.770839E-04   2.661796E-04   7.290201E-04   9.938183E-04
  ag   35   3.763696E-06   4.053861E-04   5.201325E-06  -3.360084E-06  -6.591526E-06   9.055678E-05   1.552164E-04   5.485956E-04
  ag   36   7.837243E-06   3.913495E-04  -9.860992E-06  -2.576182E-05  -2.217612E-05   1.412622E-04   2.710239E-04   4.949204E-05
  ag   37  -4.730964E-05   3.814982E-05   1.402270E-04   2.643806E-04   1.150577E-06   1.789081E-05   2.513884E-05   1.242480E-05
  ag   38   3.412073E-06  -2.612009E-04  -4.856461E-07  -7.800235E-06  -5.010085E-05  -3.526016E-04  -4.095851E-04  -1.448805E-04
  ag   39   3.775600E-07   2.043418E-04   8.064067E-07   3.498338E-06  -1.320532E-04  -1.808671E-04  -2.609644E-04   1.690800E-04

               ag   25        ag   26        ag   27        ag   28        ag   29        ag   30        ag   31        ag   32
  ag    3  -1.750149E-03   4.973958E-05  -1.359521E-03   3.710399E-04  -6.752385E-04  -4.887150E-05   1.744627E-03  -1.467284E-03
  ag    4  -4.700418E-05   1.221030E-03  -3.773294E-04   2.188452E-03   2.307025E-04   2.743987E-04   1.304938E-04  -9.887123E-05
  ag    5   1.753324E-03  -1.811472E-04   1.898040E-03  -1.383506E-04  -2.358409E-04  -1.427707E-05   3.351150E-04   4.502548E-03
  ag    6  -5.867448E-04   3.325092E-05  -4.076131E-04  -1.168156E-03  -4.036209E-04   2.529289E-03   3.572251E-04  -5.160468E-04
  ag    7   7.156471E-04  -7.728217E-06  -4.877350E-04   1.212430E-04  -3.987300E-04  -1.411584E-06  -1.596107E-04   1.659653E-04
  ag    8   3.270872E-06   1.283446E-04  -7.947046E-08  -8.099373E-07   4.229295E-06  -1.232611E-04  -1.120384E-06  -1.409301E-06
  ag    9   1.686752E-03   4.559749E-06  -4.297555E-04  -1.633028E-04   5.496896E-04   2.807572E-06  -7.007645E-04   1.126544E-04
  ag   10   2.915228E-05  -1.917608E-04  -2.207148E-05   9.380215E-05   9.949984E-06   1.254243E-04  -4.838179E-06   3.528455E-06
  ag   11  -6.278677E-04   1.201290E-05   9.982980E-04  -2.718474E-04   9.273946E-04   5.646504E-06   8.231387E-05  -1.916225E-04
  ag   12  -1.515760E-05  -2.848870E-04   7.586419E-06  -2.680339E-04  -6.948950E-05   2.345847E-04   1.712006E-06   1.750112E-06
  ag   13   6.054415E-06  -5.120038E-05   6.257295E-07   5.686521E-04   1.452211E-04   4.779061E-05   2.077304E-06  -7.204776E-07
  ag   14  -1.608098E-03   5.556821E-06   1.793537E-03   1.440747E-04  -5.621660E-04  -6.752852E-06   9.132525E-04   1.189233E-04
  ag   15   1.767655E-05   5.065665E-05  -3.257843E-06  -1.621724E-04  -5.092174E-05  -1.063277E-04   1.911045E-05  -1.075474E-05
  ag   16  -4.697630E-05  -9.036570E-04   1.754612E-05  -1.101033E-03  -2.971512E-04   3.138187E-04   7.038060E-06   8.101009E-06
  ag   17   5.701119E-04   1.490652E-05   1.231981E-03  -3.846253E-04   1.360354E-03   6.478551E-06  -4.855704E-04   5.681528E-05
  ag   18  -6.125523E-05  -8.570411E-04   2.770613E-05   5.023078E-04   1.342130E-04   3.216715E-04  -1.786533E-05   2.200851E-05
  ag   19   1.190615E-03  -5.048902E-06   3.404748E-04   4.232654E-05  -1.587932E-04  -4.936370E-07  -7.217664E-04   3.187139E-04
  ag   20  -8.194869E-04  -4.067048E-06  -1.732162E-03  -2.281307E-04   8.575716E-04   1.267368E-05  -3.897745E-04  -1.204793E-03
  ag   21   9.489786E-06   2.552807E-05  -3.444701E-06   1.187872E-03   3.128583E-04   2.571397E-04  -4.899108E-06  -2.968369E-07
  ag   22   2.382887E-05   6.692779E-04  -8.680336E-06   1.258345E-03   3.335605E-04  -3.630540E-05  -9.672773E-06  -4.733289E-06
  ag   23   1.345073E-05   9.925890E-04  -3.944706E-06   7.337419E-04   1.964432E-04   2.163241E-04  -8.307404E-06  -8.849006E-06
  ag   24  -1.401436E-05  -9.717304E-04   1.651229E-05   3.460649E-04   9.577160E-05   1.364351E-05  -1.034576E-05   2.114477E-05
  ag   25   3.628982E-03   2.162807E-05  -2.509737E-04  -1.768263E-05   1.087180E-04  -4.417395E-06  -2.650688E-04   1.250914E-03
  ag   26   2.162807E-05   8.547047E-04   1.844839E-06   8.867116E-05   2.856090E-05   5.609652E-05  -1.925966E-06  -1.009262E-06
  ag   27  -2.509737E-04   1.844839E-06   1.528453E-03  -5.133188E-05   1.847010E-04   4.043076E-07  -4.060552E-04   3.658352E-04
  ag   28  -1.768263E-05   8.867116E-05  -5.133188E-05   1.562324E-03   3.217463E-05   3.023311E-04   4.703622E-05   9.156860E-05
  ag   29   1.087180E-04   2.856090E-05   1.847010E-04   3.217463E-05   1.435407E-03   8.658674E-05  -2.006372E-04  -3.232544E-04
  ag   30  -4.417395E-06   5.609652E-05   4.043076E-07   3.023311E-04   8.658674E-05   7.659212E-04  -2.914950E-06  -5.385279E-06
  ag   31  -2.650688E-04  -1.925966E-06  -4.060552E-04   4.703622E-05  -2.006372E-04  -2.914950E-06   1.163632E-03  -2.863048E-04
  ag   32   1.250914E-03  -1.009262E-06   3.658352E-04   9.156860E-05  -3.232544E-04  -5.385279E-06  -2.863048E-04   1.735950E-03
  ag   33   1.793521E-05   3.217752E-04  -4.470824E-06  -3.292324E-04  -8.870929E-05   9.757452E-05   4.508014E-07   1.515023E-06
  ag   34   1.062555E-06  -2.400292E-05   2.477042E-06  -6.416142E-05  -1.915292E-05   3.218923E-05  -2.724926E-07   2.346301E-06
  ag   35  -5.490040E-06  -2.370596E-04   5.209778E-06   3.879014E-04   1.027527E-04   2.678847E-04  -2.874654E-06  -6.399365E-07
  ag   36  -7.211279E-05   2.385524E-04   3.002523E-05   1.801677E-04   3.785708E-05   1.419488E-04   1.272150E-05  -4.142655E-05
  ag   37   7.924105E-04   2.213988E-05  -3.025413E-04  -4.260674E-06   1.087761E-04   1.327777E-05  -1.831431E-04   4.583198E-04
  ag   38  -3.989165E-06  -2.279902E-04   4.226058E-06  -9.210986E-05  -2.389954E-05   5.877888E-05   2.202350E-06   4.379242E-09
  ag   39  -9.048527E-06  -1.348668E-04   2.691443E-06  -2.577450E-05  -5.494494E-06   3.815815E-06  -1.718381E-06  -5.127072E-07

               ag   33        ag   34        ag   35        ag   36        ag   37        ag   38        ag   39
  ag    3  -1.789839E-05  -3.376652E-04  -1.135619E-04   2.919074E-05  -9.866949E-04  -6.169663E-05  -3.167306E-05
  ag    4  -2.982552E-03  -1.101412E-03   2.334499E-03   9.409949E-04   5.632677E-05  -3.064375E-03   5.845539E-04
  ag    5   2.614068E-04   5.384105E-05  -8.669806E-05   3.221830E-04  -2.713102E-03   1.815423E-04  -1.081949E-05
  ag    6  -4.090752E-03   3.965589E-04  -1.132407E-03  -3.959408E-03  -5.507022E-04  -3.110400E-03   6.342194E-04
  ag    7  -3.414198E-06  -9.957534E-06  -2.271715E-06  -1.426912E-05   1.504299E-04  -1.128000E-06   1.548386E-07
  ag    8   1.467077E-04   1.515289E-04  -1.400010E-06   7.412974E-06   1.210812E-06  -1.038322E-04   1.021115E-05
  ag    9   7.471365E-06   1.313000E-05   2.959585E-06  -2.925615E-05   3.395149E-04  -3.968893E-06  -3.439979E-07
  ag   10  -4.284821E-04  -3.077031E-04  -3.525399E-05  -8.371400E-05  -1.549628E-06   1.811680E-04  -4.544438E-05
  ag   11  -1.166284E-06   1.375320E-05   3.900443E-06   2.020859E-05  -2.153625E-04   5.300698E-06  -1.453092E-06
  ag   12  -8.165668E-05  -2.969348E-04   8.849673E-06   1.192707E-05  -1.600367E-06   2.170235E-04   4.244514E-05
  ag   13  -4.212273E-04  -2.466694E-04  -4.701471E-05  -9.098438E-05  -8.537309E-06   2.940321E-05  -8.455237E-05
  ag   14   3.417806E-06   1.131054E-05   4.307491E-06   6.424141E-05  -6.764382E-04   7.185816E-06  -1.175799E-06
  ag   15  -7.606738E-04  -5.858264E-04  -1.698727E-04  -3.131799E-04  -3.505869E-05   3.170449E-04  -1.177062E-04
  ag   16   1.164531E-05  -5.226105E-04   4.734356E-05   7.697535E-06  -6.527135E-06   4.726741E-04   2.399875E-04
  ag   17   3.588337E-06   1.165933E-05   3.763696E-06   7.837243E-06  -4.730964E-05   3.412073E-06   3.775600E-07
  ag   18   6.560730E-05   8.251291E-04   4.053861E-04   3.913495E-04   3.814982E-05  -2.612009E-04   2.043418E-04
  ag   19  -5.165098E-07   7.695955E-06   5.201325E-06  -9.860992E-06   1.402270E-04  -4.856461E-07   8.064067E-07
  ag   20  -1.519894E-07  -9.004277E-06  -3.360084E-06  -2.576182E-05   2.643806E-04  -7.800235E-06   3.498338E-06
  ag   21  -4.141827E-04  -3.770839E-04  -6.591526E-06  -2.217612E-05   1.150577E-06  -5.010085E-05  -1.320532E-04
  ag   22   8.976337E-05   2.661796E-04   9.055678E-05   1.412622E-04   1.789081E-05  -3.526016E-04  -1.808671E-04
  ag   23   6.656007E-04   7.290201E-04   1.552164E-04   2.710239E-04   2.513884E-05  -4.095851E-04  -2.609644E-04
  ag   24  -9.644439E-04   9.938183E-04   5.485956E-04   4.949204E-05   1.242480E-05  -1.448805E-04   1.690800E-04
  ag   25   1.793521E-05   1.062555E-06  -5.490040E-06  -7.211279E-05   7.924105E-04  -3.989165E-06  -9.048527E-06
  ag   26   3.217752E-04  -2.400292E-05  -2.370596E-04   2.385524E-04   2.213988E-05  -2.279902E-04  -1.348668E-04
  ag   27  -4.470824E-06   2.477042E-06   5.209778E-06   3.002523E-05  -3.025413E-04   4.226058E-06   2.691443E-06
  ag   28  -3.292324E-04  -6.416142E-05   3.879014E-04   1.801677E-04  -4.260674E-06  -9.210986E-05  -2.577450E-05
  ag   29  -8.870929E-05  -1.915292E-05   1.027527E-04   3.785708E-05   1.087761E-04  -2.389954E-05  -5.494494E-06
  ag   30   9.757452E-05   3.218923E-05   2.678847E-04   1.419488E-04   1.327777E-05   5.877888E-05   3.815815E-06
  ag   31   4.508014E-07  -2.724926E-07  -2.874654E-06   1.272150E-05  -1.831431E-04   2.202350E-06  -1.718381E-06
  ag   32   1.515023E-06   2.346301E-06  -6.399365E-07  -4.142655E-05   4.583198E-04   4.379242E-09  -5.127072E-07
  ag   33   1.465318E-03  -1.063565E-04  -3.027208E-04   2.966270E-04   2.937322E-05   1.957902E-04  -1.915336E-04
  ag   34  -1.063565E-04   1.048820E-03   3.086557E-04   1.754688E-05   2.087097E-06  -6.793450E-05  -1.195814E-04
  ag   35  -3.027208E-04   3.086557E-04   8.634647E-04   2.313419E-05   3.038340E-06   1.273978E-04   2.282771E-04
  ag   36   2.966270E-04   1.754688E-05   2.313419E-05   4.073543E-04  -1.283199E-05   9.006487E-05  -2.122283E-05
  ag   37   2.937322E-05   2.087097E-06   3.038340E-06  -1.283199E-05   5.583644E-04   7.775620E-06  -2.556304E-06
  ag   38   1.957902E-04  -6.793450E-05   1.273978E-04   9.006487E-05   7.775620E-06   4.983512E-04   1.181870E-04
  ag   39  -1.915336E-04  -1.195814E-04   2.282771E-04  -2.122283E-05  -2.556304E-06   1.181870E-04   3.231697E-04

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98744591     1.98141275     1.97932543     1.97357643     0.01453694     0.01266887
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.01034229     0.00671804     0.00542727     0.00465656     0.00235113     0.00180410     0.00150184     0.00126916
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00092245     0.00076621     0.00046925     0.00046562     0.00030068     0.00026606     0.00014108     0.00012175
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00010534     0.00009542     0.00008854     0.00005147     0.00002200     0.00002045     0.00001620     0.00000839
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000650     0.00000646     0.00000280     0.00000130     0.00000022     0.00000015     0.00000003

          modens reordered block   1

               b3u   1        b3u   2        b3u   3        b3u   4        b3u   5        b3u   6        b3u   7        b3u   8
  b3u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   3    0.00000        0.00000        1.98460      -1.986167E-05   1.740045E-03  -3.292177E-04   1.186246E-04   8.466196E-04
  b3u   4    0.00000        0.00000      -1.986167E-05    1.97867       9.369249E-05   3.754759E-04  -9.060223E-04  -5.012881E-04
  b3u   5    0.00000        0.00000       1.740045E-03   9.369249E-05    1.97693      -2.333723E-03   7.189526E-05   2.454463E-03
  b3u   6    0.00000        0.00000      -3.292177E-04   3.754759E-04  -2.333723E-03   3.875796E-04  -2.697093E-06  -4.981506E-04
  b3u   7    0.00000        0.00000       1.186246E-04  -9.060223E-04   7.189526E-05  -2.697093E-06   1.038718E-04   4.965548E-06
  b3u   8    0.00000        0.00000       8.466196E-04  -5.012881E-04   2.454463E-03  -4.981506E-04   4.965548E-06   7.073232E-04
  b3u   9    0.00000        0.00000      -1.791742E-05  -4.901563E-04   2.467888E-03  -5.993285E-04   5.495041E-06   6.750897E-04
  b3u  10    0.00000        0.00000      -3.609189E-04   1.927423E-03  -1.344023E-04  -7.753704E-07  -3.259727E-04  -3.890650E-06
  b3u  11    0.00000        0.00000      -1.387515E-03   7.455316E-05   3.696146E-06   2.515653E-04  -2.761479E-06  -5.048956E-04
  b3u  12    0.00000        0.00000       4.880489E-04  -1.242499E-03   7.937810E-05  -3.693546E-05   4.391546E-04   4.435406E-05
  b3u  13    0.00000        0.00000       7.960442E-05  -2.470311E-04   1.458721E-03   9.249278E-04   9.325982E-06  -9.042196E-04
  b3u  14    0.00000        0.00000      -5.483539E-04  -2.472986E-04   2.046179E-05  -1.172925E-06  -4.323989E-04  -7.855436E-06
  b3u  15    0.00000        0.00000       2.905124E-04  -1.221254E-03   4.418522E-03  -5.803414E-04  -2.454469E-06   8.065050E-04
  b3u  16    0.00000        0.00000      -1.961730E-03  -8.510899E-04   4.626041E-03   2.337346E-04   8.360136E-07  -3.820587E-04
  b3u  17    0.00000        0.00000       7.007844E-04  -4.625432E-04   1.606113E-03  -1.968086E-04   1.392299E-06   4.816376E-04
  b3u  18    0.00000        0.00000      -9.485072E-04  -4.323465E-03   8.277250E-04   1.297542E-06  -4.070121E-04  -8.459289E-06
  b3u  19    0.00000        0.00000      -3.740696E-04   1.972231E-04  -2.628515E-03  -3.264851E-04  -4.707046E-06   1.907987E-04
  b3u  20    0.00000        0.00000      -4.569627E-05   1.132064E-03  -4.208934E-04   7.497774E-06  -2.650365E-04  -1.298201E-05
  b3u  21    0.00000        0.00000      -7.913257E-04  -8.568783E-04   1.169925E-03  -7.082356E-04   1.384956E-06   9.234874E-04
  b3u  22    0.00000        0.00000      -6.415709E-04  -3.424597E-03   1.149952E-03  -3.665151E-06   2.171751E-04   5.492541E-06
  b3u  23    0.00000        0.00000      -5.898213E-04  -1.549924E-03   6.621455E-04  -2.417812E-06  -1.172579E-04  -1.427162E-07
  b3u  24    0.00000        0.00000       1.591811E-03   2.628185E-04  -2.209855E-03   1.804208E-05  -5.769380E-07   3.437888E-05
  b3u  25    0.00000        0.00000      -2.214796E-04  -1.163967E-04  -2.317430E-04  -3.489355E-04   1.214162E-06   4.740903E-04
  b3u  26    0.00000        0.00000       4.789970E-04   1.756436E-04  -1.306091E-03  -1.775182E-04   1.376872E-06   1.438122E-04
  b3u  27    0.00000        0.00000      -9.320947E-05   2.436207E-03   3.281498E-04   2.145965E-06   2.030744E-04  -2.307246E-06
  b3u  28    0.00000        0.00000      -6.565458E-04   7.348038E-05  -2.474697E-03   6.504633E-06  -4.534526E-06  -1.742182E-04
  b3u  29    0.00000        0.00000      -1.658474E-03   1.592884E-05   1.571357E-03   2.496075E-04  -2.601459E-06  -4.564177E-04
  b3u  30    0.00000        0.00000      -7.498853E-05  -1.532880E-03   4.316325E-05  -2.948031E-06  -7.354714E-05   4.534880E-06
  b3u  31    0.00000        0.00000       1.930600E-04   7.656991E-05  -4.031702E-03   9.392888E-05  -7.313054E-07  -6.184858E-05
  b3u  32    0.00000        0.00000       1.166338E-03  -1.369345E-04   2.929499E-04  -6.700569E-05  -4.539107E-07   9.731320E-05
  b3u  33    0.00000        0.00000       1.207541E-04  -2.321955E-03  -1.912059E-04  -9.948372E-07   9.956222E-05   2.167979E-06
  b3u  34    0.00000        0.00000      -1.875441E-03  -1.026571E-04  -1.005352E-03  -1.491934E-04   8.089882E-08   2.021525E-04
  b3u  35    0.00000        0.00000       2.894684E-04  -3.063849E-03   4.727847E-04  -7.958283E-07  -1.218429E-05   9.745393E-06
  b3u  36    0.00000        0.00000      -7.788541E-04  -1.134393E-03  -1.515492E-03  -2.172096E-06  -5.061468E-06  -2.051763E-05
  b3u  37    0.00000        0.00000      -2.063742E-03  -1.486389E-05   2.558786E-03  -4.876644E-05   4.189619E-07   6.363698E-05
  b3u  38    0.00000        0.00000       6.146723E-05   1.037972E-04  -1.441331E-05   1.174856E-07   6.063022E-05   7.260197E-07
  b3u  39    0.00000        0.00000       2.807370E-05   1.067240E-03   3.369946E-05   2.015932E-07   1.100708E-05   3.500894E-07

               b3u   9        b3u  10        b3u  11        b3u  12        b3u  13        b3u  14        b3u  15        b3u  16
  b3u   3  -1.791742E-05  -3.609189E-04  -1.387515E-03   4.880489E-04   7.960442E-05  -5.483539E-04   2.905124E-04  -1.961730E-03
  b3u   4  -4.901563E-04   1.927423E-03   7.455316E-05  -1.242499E-03  -2.470311E-04  -2.472986E-04  -1.221254E-03  -8.510899E-04
  b3u   5   2.467888E-03  -1.344023E-04   3.696146E-06   7.937810E-05   1.458721E-03   2.046179E-05   4.418522E-03   4.626041E-03
  b3u   6  -5.993285E-04  -7.753704E-07   2.515653E-04  -3.693546E-05   9.249278E-04  -1.172925E-06  -5.803414E-04   2.337346E-04
  b3u   7   5.495041E-06  -3.259727E-04  -2.761479E-06   4.391546E-04   9.325982E-06  -4.323989E-04  -2.454469E-06   8.360136E-07
  b3u   8   6.750897E-04  -3.890650E-06  -5.048956E-04   4.435406E-05  -9.042196E-04  -7.855436E-06   8.065050E-04  -3.820587E-04
  b3u   9   1.137589E-03  -3.221580E-06  -7.723238E-05   8.695241E-05  -2.082244E-03  -4.686246E-06   9.145768E-04  -2.495275E-04
  b3u  10  -3.221580E-06   1.042808E-03   2.960922E-06  -1.429585E-03  -5.249710E-05   1.480116E-03   2.311373E-05  -8.467414E-06
  b3u  11  -7.723238E-05   2.960922E-06   9.549143E-04  -1.401015E-05   3.348959E-04  -7.662196E-06   2.539667E-04   1.557524E-03
  b3u  12   8.695241E-05  -1.429585E-03  -1.401015E-05   2.010160E-03  -1.472811E-04  -2.142909E-03  -2.755384E-05  -9.141851E-05
  b3u  13  -2.082244E-03  -5.249710E-05   3.348959E-04  -1.472811E-04   6.181037E-03  -1.262897E-04   5.240237E-04   3.548962E-03
  b3u  14  -4.686246E-06   1.480116E-03  -7.662196E-06  -2.142909E-03  -1.262897E-04   2.630574E-03  -1.504654E-05  -9.329252E-05
  b3u  15   9.145768E-04   2.311373E-05   2.539667E-04  -2.755384E-05   5.240237E-04  -1.504654E-05   3.201314E-03   3.082640E-03
  b3u  16  -2.495275E-04  -8.467414E-06   1.557524E-03  -9.141851E-05   3.548962E-03  -9.329252E-05   3.082640E-03   6.356231E-03
  b3u  17  -1.292681E-04   1.172065E-06  -6.919317E-04  -4.680437E-05   1.363314E-03  -9.894758E-06   8.148326E-04   4.629105E-04
  b3u  18  -2.174538E-05   1.549941E-03  -1.439617E-06  -2.449691E-03  -3.792743E-05   3.619692E-03   4.146641E-05  -1.509413E-05
  b3u  19   9.793598E-04   2.667136E-05   2.444332E-04   6.962081E-05  -3.135484E-03   8.489174E-05  -4.237308E-04  -1.534148E-03
  b3u  20  -1.351998E-05   7.883767E-04   4.891312E-06  -1.032630E-03  -2.408047E-05   6.922957E-04  -5.478098E-06  -1.174406E-05
  b3u  21   1.443823E-03   1.377548E-05   2.186199E-04   6.190673E-05  -1.614711E-03  -1.528686E-05   3.054214E-03   2.185650E-03
  b3u  22   2.508359E-06  -5.953385E-04  -9.560891E-06   6.662592E-04   1.929262E-05  -4.711647E-05  -2.238094E-05  -3.529764E-05
  b3u  23   5.936589E-07   4.594612E-04  -2.144771E-06  -7.847707E-04  -3.311604E-05   1.174876E-03   8.260877E-06  -2.976066E-05
  b3u  24  -8.363888E-05   6.135064E-07  -5.862085E-04   3.277702E-05  -1.234420E-03   2.866673E-05  -1.458928E-03  -3.159188E-03
  b3u  25   4.861976E-04   4.395300E-06  -5.318964E-04   3.627608E-05  -1.212459E-03   1.820216E-05  -1.501056E-04  -1.471106E-03
  b3u  26   3.007067E-04  -8.624194E-07   9.806077E-05   1.490640E-05  -3.162403E-04  -2.514538E-06   7.845834E-05   2.136844E-04
  b3u  27   7.010372E-06  -7.117414E-04   7.636388E-06   1.062184E-03   1.046844E-05  -1.218955E-03  -2.475167E-05  -1.111062E-06
  b3u  28   3.991487E-04   1.620271E-05   5.965155E-04   2.976631E-05  -1.466188E-03   3.556277E-05   8.920189E-06   1.141302E-04
  b3u  29  -2.837756E-04   2.352372E-06   5.444499E-04  -2.179949E-05   4.203175E-04   8.357981E-06  -5.237794E-04   3.825534E-04
  b3u  30   3.208158E-06   1.729583E-04  -5.781902E-06  -1.930382E-04  -1.129061E-05  -8.897899E-05   5.969624E-06  -3.542261E-06
  b3u  31  -2.119604E-04  -6.109914E-07  -2.578846E-05  -1.520686E-05   4.572253E-04  -7.948161E-06   6.650337E-06  -8.784756E-05
  b3u  32   1.164865E-04   3.663135E-06  -2.380175E-05   1.930837E-06  -1.833899E-04   6.597122E-06   1.919977E-04  -3.509934E-04
  b3u  33   2.843368E-06  -3.390265E-04  -7.399501E-07   5.246529E-04   1.038374E-05  -5.173420E-04  -9.721030E-06  -1.419571E-06
  b3u  34   3.000088E-04   3.578167E-06  -5.944688E-05   1.801149E-05  -5.455515E-04   3.013727E-06   4.339858E-04   3.563176E-05
  b3u  35  -2.708597E-05   8.520036E-05  -2.497446E-05  -2.001797E-04   1.256559E-04   3.297328E-04   3.495196E-05   3.396669E-05
  b3u  36   7.536628E-05   3.468589E-05   6.482836E-05  -6.694545E-05  -3.403244E-04   1.351442E-04  -5.860063E-05  -8.139320E-05
  b3u  37   7.360204E-05  -2.650513E-07   3.382633E-05  -8.938829E-07   1.095868E-04  -5.529932E-06   2.938079E-04   4.411253E-04
  b3u  38   6.520146E-07  -2.050808E-04  -9.308536E-07   2.950208E-04   8.964271E-06  -3.365140E-04  -6.928655E-06  -1.093722E-06
  b3u  39  -6.404197E-07  -3.155636E-05  -1.040781E-06   3.924245E-05   3.568413E-06  -5.334930E-05  -4.844426E-07   1.034265E-06

               b3u  17        b3u  18        b3u  19        b3u  20        b3u  21        b3u  22        b3u  23        b3u  24
  b3u   3   7.007844E-04  -9.485072E-04  -3.740696E-04  -4.569627E-05  -7.913257E-04  -6.415709E-04  -5.898213E-04   1.591811E-03
  b3u   4  -4.625432E-04  -4.323465E-03   1.972231E-04   1.132064E-03  -8.568783E-04  -3.424597E-03  -1.549924E-03   2.628185E-04
  b3u   5   1.606113E-03   8.277250E-04  -2.628515E-03  -4.208934E-04   1.169925E-03   1.149952E-03   6.621455E-04  -2.209855E-03
  b3u   6  -1.968086E-04   1.297542E-06  -3.264851E-04   7.497774E-06  -7.082356E-04  -3.665151E-06  -2.417812E-06   1.804208E-05
  b3u   7   1.392299E-06  -4.070121E-04  -4.707046E-06  -2.650365E-04   1.384956E-06   2.171751E-04  -1.172579E-04  -5.769380E-07
  b3u   8   4.816376E-04  -8.459289E-06   1.907987E-04  -1.298201E-05   9.234874E-04   5.492541E-06  -1.427162E-07   3.437888E-05
  b3u   9  -1.292681E-04  -2.174538E-05   9.793598E-04  -1.351998E-05   1.443823E-03   2.508359E-06   5.936589E-07  -8.363888E-05
  b3u  10   1.172065E-06   1.549941E-03   2.667136E-05   7.883767E-04   1.377548E-05  -5.953385E-04   4.594612E-04   6.135064E-07
  b3u  11  -6.919317E-04  -1.439617E-06   2.444332E-04   4.891312E-06   2.186199E-04  -9.560891E-06  -2.144771E-06  -5.862085E-04
  b3u  12  -4.680437E-05  -2.449691E-03   6.962081E-05  -1.032630E-03   6.190673E-05   6.662592E-04  -7.847707E-04   3.277702E-05
  b3u  13   1.363314E-03  -3.792743E-05  -3.135484E-03  -2.408047E-05  -1.614711E-03   1.929262E-05  -3.311604E-05  -1.234420E-03
  b3u  14  -9.894758E-06   3.619692E-03   8.489174E-05   6.922957E-04  -1.528686E-05  -4.711647E-05   1.174876E-03   2.866673E-05
  b3u  15   8.148326E-04   4.146641E-05  -4.237308E-04  -5.478098E-06   3.054214E-03  -2.238094E-05   8.260877E-06  -1.458928E-03
  b3u  16   4.629105E-04  -1.509413E-05  -1.534148E-03  -1.174406E-05   2.185650E-03  -3.529764E-05  -2.976066E-05  -3.159188E-03
  b3u  17   1.277166E-03   3.743733E-05  -1.123467E-03  -7.447203E-06   1.409252E-04   1.987162E-05   1.336191E-05  -3.987389E-04
  b3u  18   3.743733E-05   6.680135E-03   7.770211E-05  -1.694646E-04  -2.269542E-06   2.645899E-03   2.935675E-03  -2.419539E-05
  b3u  19  -1.123467E-03   7.770211E-05   1.977649E-03  -4.810573E-06   4.740685E-04   4.193323E-05   4.850949E-05   7.612243E-04
  b3u  20  -7.447203E-06  -1.694646E-04  -4.810573E-06   1.159749E-03  -2.192837E-05  -1.765258E-03  -3.143414E-04   1.753936E-05
  b3u  21   1.409252E-04  -2.269542E-06   4.740685E-04  -2.192837E-05   4.967047E-03  -3.126818E-05  -1.825339E-06  -1.584049E-03
  b3u  22   1.987162E-05   2.645899E-03   4.193323E-05  -1.765258E-03  -3.126818E-05   5.392250E-03   2.675817E-03  -1.214025E-05
  b3u  23   1.336191E-05   2.935675E-03   4.850949E-05  -3.143414E-04  -1.825339E-06   2.675817E-03   2.150428E-03  -2.691792E-06
  b3u  24  -3.987389E-04  -2.419539E-05   7.612243E-04   1.753936E-05  -1.584049E-03  -1.214025E-05  -2.691792E-06   2.513964E-03
  b3u  25   2.545517E-04   5.527373E-07   5.876246E-04   5.942133E-06  -5.289503E-04   7.371664E-06   9.123336E-06   7.455775E-04
  b3u  26   5.941858E-05  -1.943039E-06   2.891493E-04   1.468494E-07  -6.255184E-04   1.463346E-05   2.933655E-06  -2.153595E-04
  b3u  27  -2.035590E-05  -1.040857E-03   3.698298E-06  -7.560528E-04   6.923647E-06   1.937228E-03   4.469559E-04  -4.878065E-06
  b3u  28  -1.030706E-03   5.926274E-06   1.154351E-03   9.017407E-06   1.082399E-03  -4.756350E-05  -1.247774E-05  -8.755674E-05
  b3u  29  -5.099231E-04   1.634307E-05   2.780049E-04   6.064190E-06  -1.316925E-03   2.263882E-05   1.234358E-05   1.261746E-04
  b3u  30   5.825062E-06  -3.493967E-04  -9.654453E-06   3.517051E-04   1.130850E-05   2.285738E-04   2.993469E-04  -4.485067E-06
  b3u  31   1.878736E-05  -1.187159E-05  -4.156430E-05   5.054756E-06  -5.129755E-05  -6.869038E-06  -1.503470E-06   4.827898E-04
  b3u  32  -6.180492E-05   1.219230E-05   1.648478E-04   1.726676E-06  -2.129570E-04   8.282931E-06   1.072980E-05   9.040686E-04
  b3u  33  -6.850320E-06  -8.153248E-04  -6.759945E-06  -2.931733E-04  -9.002194E-06  -3.238757E-04  -7.617431E-04   1.006546E-05
  b3u  34   5.013052E-05  -9.415616E-07   2.601438E-04  -2.954604E-06   7.324592E-04  -4.342363E-06  -1.600509E-07  -1.569373E-04
  b3u  35   7.836685E-05   7.630904E-04  -5.087980E-05  -2.366789E-06  -1.019705E-04   1.839409E-04   3.355552E-04   2.515537E-05
  b3u  36  -1.848487E-04   2.989882E-04   1.626767E-04  -1.954483E-06   2.916034E-04   6.879502E-05   1.313818E-04  -7.977904E-05
  b3u  37   1.133766E-04   1.463275E-06  -8.925387E-05  -4.150990E-06   3.167393E-04   3.036724E-06  -2.415436E-07  -4.749228E-04
  b3u  38  -9.992866E-07  -3.376935E-04  -4.497714E-06  -1.585078E-04  -4.529409E-06   2.889635E-04  -1.123868E-06   7.659524E-07
  b3u  39   5.763577E-07  -1.742245E-04  -4.524739E-06   4.072042E-05  -3.756054E-07  -2.962714E-04  -1.459047E-04   9.167419E-07

               b3u  25        b3u  26        b3u  27        b3u  28        b3u  29        b3u  30        b3u  31        b3u  32
  b3u   3  -2.214796E-04   4.789970E-04  -9.320947E-05  -6.565458E-04  -1.658474E-03  -7.498853E-05   1.930600E-04   1.166338E-03
  b3u   4  -1.163967E-04   1.756436E-04   2.436207E-03   7.348038E-05   1.592884E-05  -1.532880E-03   7.656991E-05  -1.369345E-04
  b3u   5  -2.317430E-04  -1.306091E-03   3.281498E-04  -2.474697E-03   1.571357E-03   4.316325E-05  -4.031702E-03   2.929499E-04
  b3u   6  -3.489355E-04  -1.775182E-04   2.145965E-06   6.504633E-06   2.496075E-04  -2.948031E-06   9.392888E-05  -6.700569E-05
  b3u   7   1.214162E-06   1.376872E-06   2.030744E-04  -4.534526E-06  -2.601459E-06  -7.354714E-05  -7.313054E-07  -4.539107E-07
  b3u   8   4.740903E-04   1.438122E-04  -2.307246E-06  -1.742182E-04  -4.564177E-04   4.534880E-06  -6.184858E-05   9.731320E-05
  b3u   9   4.861976E-04   3.007067E-04   7.010372E-06   3.991487E-04  -2.837756E-04   3.208158E-06  -2.119604E-04   1.164865E-04
  b3u  10   4.395300E-06  -8.624194E-07  -7.117414E-04   1.620271E-05   2.352372E-06   1.729583E-04  -6.109914E-07   3.663135E-06
  b3u  11  -5.318964E-04   9.806077E-05   7.636388E-06   5.965155E-04   5.444499E-04  -5.781902E-06  -2.578846E-05  -2.380175E-05
  b3u  12   3.627608E-05   1.490640E-05   1.062184E-03   2.976631E-05  -2.179949E-05  -1.930382E-04  -1.520686E-05   1.930837E-06
  b3u  13  -1.212459E-03  -3.162403E-04   1.046844E-05  -1.466188E-03   4.203175E-04  -1.129061E-05   4.572253E-04  -1.833899E-04
  b3u  14   1.820216E-05  -2.514538E-06  -1.218955E-03   3.556277E-05   8.357981E-06  -8.897899E-05  -7.948161E-06   6.597122E-06
  b3u  15  -1.501056E-04   7.845834E-05  -2.475167E-05   8.920189E-06  -5.237794E-04   5.969624E-06   6.650337E-06   1.919977E-04
  b3u  16  -1.471106E-03   2.136844E-04  -1.111062E-06   1.141302E-04   3.825534E-04  -3.542261E-06  -8.784756E-05  -3.509934E-04
  b3u  17   2.545517E-04   5.941858E-05  -2.035590E-05  -1.030706E-03  -5.099231E-04   5.825062E-06   1.878736E-05  -6.180492E-05
  b3u  18   5.527373E-07  -1.943039E-06  -1.040857E-03   5.926274E-06   1.634307E-05  -3.493967E-04  -1.187159E-05   1.219230E-05
  b3u  19   5.876246E-04   2.891493E-04   3.698298E-06   1.154351E-03   2.780049E-04  -9.654453E-06  -4.156430E-05   1.648478E-04
  b3u  20   5.942133E-06   1.468494E-07  -7.560528E-04   9.017407E-06   6.064190E-06   3.517051E-04   5.054756E-06   1.726676E-06
  b3u  21  -5.289503E-04  -6.255184E-04   6.923647E-06   1.082399E-03  -1.316925E-03   1.130850E-05  -5.129755E-05  -2.129570E-04
  b3u  22   7.371664E-06   1.463346E-05   1.937228E-03  -4.756350E-05   2.263882E-05   2.285738E-04  -6.869038E-06   8.282931E-06
  b3u  23   9.123336E-06   2.933655E-06   4.469559E-04  -1.247774E-05   1.234358E-05   2.993469E-04  -1.503470E-06   1.072980E-05
  b3u  24   7.455775E-04  -2.153595E-04  -4.878065E-06  -8.755674E-05   1.261746E-04  -4.485067E-06   4.827898E-04   9.040686E-04
  b3u  25   1.377704E-03   6.165423E-04  -1.385184E-05  -3.997756E-04   1.374638E-04   3.983417E-06  -4.571937E-04   4.714562E-04
  b3u  26   6.165423E-04   1.281427E-03  -6.055469E-06  -3.037196E-04   3.564965E-04   4.755801E-06  -6.018887E-04   2.876450E-04
  b3u  27  -1.385184E-05  -6.055469E-06   1.670702E-03  -5.709412E-06   9.911392E-06   4.110156E-04   1.235079E-05  -9.376531E-06
  b3u  28  -3.997756E-04  -3.037196E-04  -5.709412E-06   1.475653E-03   1.952217E-04  -1.591480E-05   3.317902E-04  -2.797563E-04
  b3u  29   1.374638E-04   3.564965E-04   9.911392E-06   1.952217E-04   1.263655E-03  -6.889132E-06  -4.606263E-05   1.344437E-04
  b3u  30   3.983417E-06   4.755801E-06   4.110156E-04  -1.591480E-05  -6.889132E-06   6.431027E-04  -3.661709E-06   9.564726E-07
  b3u  31  -4.571937E-04  -6.018887E-04   1.235079E-05   3.317902E-04  -4.606263E-05  -3.661709E-06   1.208503E-03  -1.645853E-04
  b3u  32   4.714562E-04   2.876450E-04  -9.376531E-06  -2.797563E-04   1.344437E-04   9.564726E-07  -1.645853E-04   1.120131E-03
  b3u  33  -1.156481E-07   1.603750E-06  -6.027727E-05   2.715190E-06  -9.594590E-07  -2.371363E-04  -1.140270E-06   2.414301E-06
  b3u  34   4.478177E-04   9.387212E-05  -2.548075E-06   9.488950E-05  -6.403975E-05   3.053193E-06  -3.827040E-04   2.091938E-04
  b3u  35   3.155007E-05  -2.870765E-05  -3.779207E-04  -6.963120E-05   6.002076E-05  -1.515247E-04   8.810016E-05   1.855265E-05
  b3u  36  -7.707860E-05   6.203570E-05  -1.472934E-04   2.028394E-04  -1.588307E-04  -5.714691E-05  -2.259727E-04  -4.852389E-05
  b3u  37  -6.853469E-05  -4.238000E-05   3.465316E-06   2.211051E-05   1.203440E-04  -9.014754E-07  -9.587083E-05  -3.694292E-04
  b3u  38  -1.399073E-06   8.290523E-07   2.386197E-04  -4.302439E-06  -9.224873E-08  -1.780436E-06   1.786467E-06  -8.300466E-07
  b3u  39   1.596570E-07  -3.030895E-07  -7.472756E-05   5.341387E-07  -1.534056E-06  -5.740142E-05  -6.126691E-08  -4.864095E-08

               b3u  33        b3u  34        b3u  35        b3u  36        b3u  37        b3u  38        b3u  39
  b3u   3   1.207541E-04  -1.875441E-03   2.894684E-04  -7.788541E-04  -2.063742E-03   6.146723E-05   2.807370E-05
  b3u   4  -2.321955E-03  -1.026571E-04  -3.063849E-03  -1.134393E-03  -1.486389E-05   1.037972E-04   1.067240E-03
  b3u   5  -1.912059E-04  -1.005352E-03   4.727847E-04  -1.515492E-03   2.558786E-03  -1.441331E-05   3.369946E-05
  b3u   6  -9.948372E-07  -1.491934E-04  -7.958283E-07  -2.172096E-06  -4.876644E-05   1.174856E-07   2.015932E-07
  b3u   7   9.956222E-05   8.089882E-08  -1.218429E-05  -5.061468E-06   4.189619E-07   6.063022E-05   1.100708E-05
  b3u   8   2.167979E-06   2.021525E-04   9.745393E-06  -2.051763E-05   6.363698E-05   7.260197E-07   3.500894E-07
  b3u   9   2.843368E-06   3.000088E-04  -2.708597E-05   7.536628E-05   7.360204E-05   6.520146E-07  -6.404197E-07
  b3u  10  -3.390265E-04   3.578167E-06   8.520036E-05   3.468589E-05  -2.650513E-07  -2.050808E-04  -3.155636E-05
  b3u  11  -7.399501E-07  -5.944688E-05  -2.497446E-05   6.482836E-05   3.382633E-05  -9.308536E-07  -1.040781E-06
  b3u  12   5.246529E-04   1.801149E-05  -2.001797E-04  -6.694545E-05  -8.938829E-07   2.950208E-04   3.924245E-05
  b3u  13   1.038374E-05  -5.455515E-04   1.256559E-04  -3.403244E-04   1.095868E-04   8.964271E-06   3.568413E-06
  b3u  14  -5.173420E-04   3.013727E-06   3.297328E-04   1.351442E-04  -5.529932E-06  -3.365140E-04  -5.334930E-05
  b3u  15  -9.721030E-06   4.339858E-04   3.495196E-05  -5.860063E-05   2.938079E-04  -6.928655E-06  -4.844426E-07
  b3u  16  -1.419571E-06   3.563176E-05   3.396669E-05  -8.139320E-05   4.411253E-04  -1.093722E-06   1.034265E-06
  b3u  17  -6.850320E-06   5.013052E-05   7.836685E-05  -1.848487E-04   1.133766E-04  -9.992866E-07   5.763577E-07
  b3u  18  -8.153248E-04  -9.415616E-07   7.630904E-04   2.989882E-04   1.463275E-06  -3.376935E-04  -1.742245E-04
  b3u  19  -6.759945E-06   2.601438E-04  -5.087980E-05   1.626767E-04  -8.925387E-05  -4.497714E-06  -4.524739E-06
  b3u  20  -2.931733E-04  -2.954604E-06  -2.366789E-06  -1.954483E-06  -4.150990E-06  -1.585078E-04   4.072042E-05
  b3u  21  -9.002194E-06   7.324592E-04  -1.019705E-04   2.916034E-04   3.167393E-04  -4.529409E-06  -3.756054E-07
  b3u  22  -3.238757E-04  -4.342363E-06   1.839409E-04   6.879502E-05   3.036724E-06   2.889635E-04  -2.962714E-04
  b3u  23  -7.617431E-04  -1.600509E-07   3.355552E-04   1.313818E-04  -2.415436E-07  -1.123868E-06  -1.459047E-04
  b3u  24   1.006546E-05  -1.569373E-04   2.515537E-05  -7.977904E-05  -4.749228E-04   7.659524E-07   9.167419E-07
  b3u  25  -1.156481E-07   4.478177E-04   3.155007E-05  -7.707860E-05  -6.853469E-05  -1.399073E-06   1.596570E-07
  b3u  26   1.603750E-06   9.387212E-05  -2.870765E-05   6.203570E-05  -4.238000E-05   8.290523E-07  -3.030895E-07
  b3u  27  -6.027727E-05  -2.548075E-06  -3.779207E-04  -1.472934E-04   3.465316E-06   2.386197E-04  -7.472756E-05
  b3u  28   2.715190E-06   9.488950E-05  -6.963120E-05   2.028394E-04   2.211051E-05  -4.302439E-06   5.341387E-07
  b3u  29  -9.594590E-07  -6.403975E-05   6.002076E-05  -1.588307E-04   1.203440E-04  -9.224873E-08  -1.534056E-06
  b3u  30  -2.371363E-04   3.053193E-06  -1.515247E-04  -5.714691E-05  -9.014754E-07  -1.780436E-06  -5.740142E-05
  b3u  31  -1.140270E-06  -3.827040E-04   8.810016E-05  -2.259727E-04  -9.587083E-05   1.786467E-06  -6.126691E-08
  b3u  32   2.414301E-06   2.091938E-04   1.855265E-05  -4.852389E-05  -3.694292E-04  -8.300466E-07  -4.864095E-08
  b3u  33   8.388446E-04   8.146338E-07  -1.386470E-04  -5.514582E-05  -1.660188E-06   3.995118E-05  -9.581676E-05
  b3u  34   8.146338E-07   5.858920E-04  -3.933657E-05   1.033577E-04   2.225778E-05  -2.213369E-06   1.898438E-09
  b3u  35  -1.386470E-04  -3.933657E-05   4.589355E-04   4.380742E-05  -4.061617E-06   7.431097E-05  -2.748498E-05
  b3u  36  -5.514582E-05   1.033577E-04   4.380742E-05   3.680653E-04   1.185993E-05   2.893660E-05  -1.093503E-05
  b3u  37  -1.660188E-06   2.225778E-05  -4.061617E-06   1.185993E-05   4.473329E-04  -3.138089E-07  -1.995270E-07
  b3u  38   3.995118E-05  -2.213369E-06   7.431097E-05   2.893660E-05  -3.138089E-07   2.535169E-04  -7.887622E-05
  b3u  39  -9.581676E-05   1.898438E-09  -2.748498E-05  -1.093503E-05  -1.995270E-07  -7.887622E-05   1.093687E-04

Natural orbital populations,block 2
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98499284     1.97870658     1.97661110     0.01436632     0.01300965     0.01084906
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00854709     0.00473831     0.00332326     0.00211244     0.00200716     0.00125676     0.00080141     0.00075327
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00044247     0.00044222     0.00032351     0.00021738     0.00017843     0.00013466     0.00009824     0.00006600
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00004397     0.00003158     0.00002927     0.00001878     0.00001183     0.00000683     0.00000673     0.00000410
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000394     0.00000144     0.00000067     0.00000019     0.00000011     0.00000006     0.00000000

          modens reordered block   1

               b2u   1        b2u   2        b2u   3        b2u   4        b2u   5        b2u   6        b2u   7        b2u   8
  b2u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b2u   2    0.00000        1.98461       1.847856E-04  -1.770648E-03   5.989913E-04   1.925933E-04   1.400538E-03   2.045417E-05
  b2u   3    0.00000       1.847856E-04    1.97435       1.714798E-05  -2.081408E-05  -4.998952E-03   1.214336E-06   1.592462E-06
  b2u   4    0.00000      -1.770648E-03   1.714798E-05    1.97687      -2.742477E-03   1.332929E-04  -2.828520E-03  -3.110268E-03
  b2u   5    0.00000       5.989913E-04  -2.081408E-05  -2.742477E-03   3.977872E-04   1.096381E-05   5.147125E-04   6.051595E-04
  b2u   6    0.00000       1.925933E-04  -4.998952E-03   1.332929E-04   1.096381E-05   4.664972E-04   9.418974E-06   1.214161E-05
  b2u   7    0.00000       1.400538E-03   1.214336E-06  -2.828520E-03   5.147125E-04   9.418974E-06   7.347222E-04   6.828488E-04
  b2u   8    0.00000       2.045417E-05   1.592462E-06  -3.110268E-03   6.051595E-04   1.214161E-05   6.828488E-04   1.145153E-03
  b2u   9    0.00000       2.257321E-03  -5.590145E-04  -1.067374E-04   2.634898E-04   1.482349E-05   5.292196E-04   7.260765E-05
  b2u  10    0.00000       2.997289E-04   1.010525E-03  -1.825533E-04  -9.305296E-04  -6.856064E-06  -9.154513E-04  -2.097941E-03
  b2u  11    0.00000      -2.387014E-04   8.116206E-03  -4.579510E-04  -1.509025E-05  -9.727488E-04  -3.533901E-06  -2.670986E-05
  b2u  12    0.00000      -7.204306E-04  -1.152943E-03   4.253902E-03  -5.849382E-04  -8.742463E-06  -8.156529E-04  -9.022619E-04
  b2u  13    0.00000       2.489971E-03  -2.248089E-03   3.983288E-03   2.399999E-04   1.376993E-05   4.002829E-04   2.492053E-04
  b2u  14    0.00000      -1.487959E-03  -1.254241E-04   1.532839E-03  -2.064763E-04  -1.422871E-05  -5.000954E-04   1.388638E-04
  b2u  15    0.00000       9.357105E-04   4.747711E-04  -2.036918E-03  -3.250411E-04  -1.813941E-06  -1.897849E-04  -9.868133E-04
  b2u  16    0.00000      -2.184929E-04   6.134491E-03  -6.055057E-04  -3.034237E-06  -9.861251E-04   1.002851E-05  -2.517652E-06
  b2u  17    0.00000       7.748663E-05  -7.796317E-04   7.015029E-04  -7.112188E-04  -1.457828E-05  -9.319496E-04  -1.434417E-03
  b2u  18    0.00000       1.471629E-03  -1.186769E-03   2.209800E-03  -1.949471E-05   2.295885E-06   3.488693E-05  -8.647612E-05
  b2u  19    0.00000      -3.879357E-04  -1.923367E-04  -3.978398E-04   3.544330E-04   1.090023E-05   4.838983E-04   4.862247E-04
  b2u  20    0.00000       9.319088E-05   1.213486E-03   1.959966E-04  -2.745856E-06   4.040241E-04  -5.629785E-06  -9.809322E-06
  b2u  21    0.00000      -1.024930E-04   1.212598E-04   6.674555E-04   1.770893E-04   1.258606E-05   1.432605E-04   2.984835E-04
  b2u  22    0.00000      -8.604483E-04   6.048529E-06   2.491221E-03  -1.194245E-05  -4.244796E-06  -1.813959E-04   3.994374E-04
  b2u  23    0.00000       1.031549E-04  -7.243823E-03   4.839919E-04   8.316655E-07   5.424589E-04  -7.165089E-06   2.841414E-06
  b2u  24    0.00000       2.240601E-03  -3.100767E-04   1.827909E-03   2.540458E-04   1.111874E-05   4.654111E-04   2.799681E-04
  b2u  25    0.00000      -4.592190E-04   4.060676E-04  -4.530977E-03   9.550777E-05   1.097853E-06   6.360636E-05   2.144563E-04
  b2u  26    0.00000       8.151093E-04  -8.254249E-05  -3.134796E-04   6.665192E-05   5.570177E-07   9.787091E-05   1.143952E-04
  b2u  27    0.00000       1.073483E-04  -3.739060E-03   1.611774E-04   1.739175E-06   2.665746E-05   1.961022E-06   3.670047E-06
  b2u  28    0.00000       1.898801E-03  -1.344743E-04  -8.890708E-04  -1.499288E-04  -2.767767E-06  -2.043223E-04  -2.984844E-04
  b2u  29    0.00000      -8.672323E-04   3.723739E-05   1.378809E-03  -9.085320E-07   5.694322E-07  -2.590511E-05   7.652342E-05
  b2u  30    0.00000       1.977722E-03  -3.411699E-04   2.693170E-03  -4.989985E-05   1.094255E-06  -6.424787E-05  -7.495035E-05

               b2u   9        b2u  10        b2u  11        b2u  12        b2u  13        b2u  14        b2u  15        b2u  16
  b2u   2   2.257321E-03   2.997289E-04  -2.387014E-04  -7.204306E-04   2.489971E-03  -1.487959E-03   9.357105E-04  -2.184929E-04
  b2u   3  -5.590145E-04   1.010525E-03   8.116206E-03  -1.152943E-03  -2.248089E-03  -1.254241E-04   4.747711E-04   6.134491E-03
  b2u   4  -1.067374E-04  -1.825533E-04  -4.579510E-04   4.253902E-03   3.983288E-03   1.532839E-03  -2.036918E-03  -6.055057E-04
  b2u   5   2.634898E-04  -9.305296E-04  -1.509025E-05  -5.849382E-04   2.399999E-04  -2.064763E-04  -3.250411E-04  -3.034237E-06
  b2u   6   1.482349E-05  -6.856064E-06  -9.727488E-04  -8.742463E-06   1.376993E-05  -1.422871E-05  -1.813941E-06  -9.861251E-04
  b2u   7   5.292196E-04  -9.154513E-04  -3.533901E-06  -8.156529E-04   4.002829E-04  -5.000954E-04  -1.897849E-04   1.002851E-05
  b2u   8   7.260765E-05  -2.097941E-03  -2.670986E-05  -9.022619E-04   2.492053E-04   1.388638E-04  -9.868133E-04  -2.517652E-06
  b2u   9   9.827569E-04  -3.060539E-04  -9.365270E-06   2.287107E-04   1.555880E-03  -7.223687E-04   2.669167E-04  -6.406430E-06
  b2u  10  -3.060539E-04   6.150139E-03   4.057953E-05  -5.129129E-04  -3.471032E-03  -1.379560E-03   3.112964E-03  -1.270707E-05
  b2u  11  -9.365270E-06   4.057953E-05   2.124785E-03  -9.545209E-07  -3.069433E-05  -6.411892E-07   2.813292E-05   2.209760E-03
  b2u  12   2.287107E-04  -5.129129E-04  -9.545209E-07   3.166685E-03   3.029955E-03   8.308556E-04  -4.175666E-04  -1.323169E-05
  b2u  13   1.555880E-03  -3.471032E-03  -3.069433E-05   3.029955E-03   6.276300E-03   4.471048E-04  -1.480795E-03  -1.348605E-05
  b2u  14  -7.223687E-04  -1.379560E-03  -6.411892E-07   8.308556E-04   4.471048E-04   1.307821E-03  -1.135992E-03   5.616728E-06
  b2u  15   2.669167E-04   3.112964E-03   2.813292E-05  -4.175666E-04  -1.480795E-03  -1.135992E-03   1.962349E-03   3.907983E-06
  b2u  16  -6.406430E-06  -1.270707E-05   2.209760E-03  -1.323169E-05  -1.348605E-05   5.616728E-06   3.907983E-06   2.427884E-03
  b2u  17   2.021390E-04   1.621127E-03   2.806335E-05   3.013553E-03   2.149647E-03   1.499965E-04   4.756466E-04  -1.362689E-06
  b2u  18   5.788250E-04  -1.198153E-03  -7.821159E-06   1.438916E-03   3.120879E-03   3.986840E-04  -7.393741E-04  -4.359507E-06
  b2u  19   5.343063E-04  -1.195265E-03  -1.565090E-05   1.425897E-04   1.458631E-03  -2.592269E-04  -5.783615E-04  -5.914367E-06
  b2u  20   6.900844E-06   2.808094E-05  -8.716871E-04  -2.250970E-06  -6.842749E-06  -7.075635E-06   1.592614E-05  -1.076491E-03
  b2u  21  -1.041017E-04  -3.046353E-04  -2.524232E-05  -8.435387E-05  -2.372433E-04  -5.631745E-05  -2.840191E-04  -2.673088E-05
  b2u  22  -6.065632E-04  -1.461976E-03  -1.304267E-05   9.076234E-06  -1.199940E-04   1.041515E-03  -1.151818E-03  -3.180335E-06
  b2u  23   2.875362E-06   9.577172E-07  -1.428443E-03   1.332063E-05   1.811061E-05  -1.315884E-06  -7.389105E-06  -1.687208E-03
  b2u  24   5.568807E-04  -4.228250E-04  -1.129205E-05  -5.209812E-04   3.951868E-04  -5.152009E-04   2.794486E-04  -2.088385E-06
  b2u  25  -2.690154E-05  -4.578022E-04  -5.022532E-06   1.044471E-05  -8.540919E-05   1.874341E-05  -4.235689E-05  -1.389154E-06
  b2u  26   2.233482E-05  -1.719749E-04   9.843038E-07  -1.993055E-04   3.372248E-04   5.969346E-05  -1.579412E-04   7.987610E-06
  b2u  27   1.675963E-06  -9.400781E-06  -9.977320E-05   1.328040E-07   1.804847E-05   2.400066E-06  -7.974097E-06  -2.926526E-04
  b2u  28  -6.067344E-05   5.414050E-04   6.035574E-06   4.294891E-04   3.583362E-05   5.165302E-05   2.580517E-04  -8.134338E-07
  b2u  29  -7.082016E-05  -3.570196E-04  -5.468362E-06   7.707921E-05   9.209990E-05   2.033654E-04  -1.669982E-04  -3.072430E-06
  b2u  30   3.306790E-05  -1.008066E-04  -3.598135E-06   2.872882E-04   4.324468E-04   1.117414E-04  -8.412039E-05  -3.717327E-06

               b2u  17        b2u  18        b2u  19        b2u  20        b2u  21        b2u  22        b2u  23        b2u  24
  b2u   2   7.748663E-05   1.471629E-03  -3.879357E-04   9.319088E-05  -1.024930E-04  -8.604483E-04   1.031549E-04   2.240601E-03
  b2u   3  -7.796317E-04  -1.186769E-03  -1.923367E-04   1.213486E-03   1.212598E-04   6.048529E-06  -7.243823E-03  -3.100767E-04
  b2u   4   7.015029E-04   2.209800E-03  -3.978398E-04   1.959966E-04   6.674555E-04   2.491221E-03   4.839919E-04   1.827909E-03
  b2u   5  -7.112188E-04  -1.949471E-05   3.544330E-04  -2.745856E-06   1.770893E-04  -1.194245E-05   8.316655E-07   2.540458E-04
  b2u   6  -1.457828E-05   2.295885E-06   1.090023E-05   4.040241E-04   1.258606E-05  -4.244796E-06   5.424589E-04   1.111874E-05
  b2u   7  -9.319496E-04   3.488693E-05   4.838983E-04  -5.629785E-06   1.432605E-04  -1.813959E-04  -7.165089E-06   4.654111E-04
  b2u   8  -1.434417E-03  -8.647612E-05   4.862247E-04  -9.809322E-06   2.984835E-04   3.994374E-04   2.841414E-06   2.799681E-04
  b2u   9   2.021390E-04   5.788250E-04   5.343063E-04   6.900844E-06  -1.041017E-04  -6.065632E-04   2.875362E-06   5.568807E-04
  b2u  10   1.621127E-03  -1.198153E-03  -1.195265E-03   2.808094E-05  -3.046353E-04  -1.461976E-03   9.577172E-07  -4.228250E-04
  b2u  11   2.806335E-05  -7.821159E-06  -1.565090E-05  -8.716871E-04  -2.524232E-05  -1.304267E-05  -1.428443E-03  -1.129205E-05
  b2u  12   3.013553E-03   1.438916E-03   1.425897E-04  -2.250970E-06  -8.435387E-05   9.076234E-06   1.332063E-05  -5.209812E-04
  b2u  13   2.149647E-03   3.120879E-03   1.458631E-03  -6.842749E-06  -2.372433E-04  -1.199940E-04   1.811061E-05   3.951868E-04
  b2u  14   1.499965E-04   3.986840E-04  -2.592269E-04  -7.075635E-06  -5.631745E-05   1.041515E-03  -1.315884E-06  -5.152009E-04
  b2u  15   4.756466E-04  -7.393741E-04  -5.783615E-04   1.592614E-05  -2.840191E-04  -1.151818E-03  -7.389105E-06   2.794486E-04
  b2u  16  -1.362689E-06  -4.359507E-06  -5.914367E-06  -1.076491E-03  -2.673088E-05  -3.180335E-06  -1.687208E-03  -2.088385E-06
  b2u  17   4.930519E-03   1.579103E-03   5.397537E-04  -2.502197E-05   6.188594E-04  -1.064638E-03   5.803236E-06  -1.325419E-03
  b2u  18   1.579103E-03   2.501870E-03   7.439973E-04  -6.306293E-06  -2.247166E-04  -8.909169E-05   1.004720E-05  -1.303683E-04
  b2u  19   5.397537E-04   7.439973E-04   1.388967E-03  -1.882259E-05   6.139607E-04  -4.065290E-04   5.603774E-06  -1.468418E-04
  b2u  20  -2.502197E-05  -6.306293E-06  -1.882259E-05   9.614366E-04  -9.833672E-06   9.796546E-06   4.846761E-04   1.195310E-05
  b2u  21   6.188594E-04  -2.247166E-04   6.139607E-04  -9.833672E-06   1.275846E-03  -2.990517E-04   1.253698E-05  -3.665194E-04
  b2u  22  -1.064638E-03  -8.909169E-05  -4.065290E-04   9.796546E-06  -2.990517E-04   1.475218E-03  -3.761371E-07  -1.947722E-04
  b2u  23   5.803236E-06   1.004720E-05   5.603774E-06   4.846761E-04   1.253698E-05  -3.761371E-07   2.051504E-03  -4.413775E-07
  b2u  24  -1.325419E-03  -1.303683E-04  -1.468418E-04   1.195310E-05  -3.665194E-04  -1.947722E-04  -4.413775E-07   1.268962E-03
  b2u  25  -4.114265E-05  -4.786202E-04   4.648972E-04  -1.501125E-05   6.050002E-04  -3.361195E-04   1.806466E-06  -5.096719E-05
  b2u  26   2.095416E-04   9.014372E-04   4.717538E-04  -1.686727E-05   2.821271E-04  -2.801412E-04  -1.130522E-05  -1.381071E-04
  b2u  27   1.013539E-05   2.806808E-05   1.596174E-05   2.172903E-04   1.276448E-05  -6.269906E-06   5.154712E-04  -4.127882E-06
  b2u  28   7.255827E-04   1.592533E-04  -4.474437E-04   3.671932E-06  -9.428723E-05  -9.128299E-05   6.087137E-08  -6.383456E-05
  b2u  29  -2.963950E-04  -7.902400E-05  -8.693093E-05  -5.878780E-07   6.707356E-05   2.117104E-04   1.862107E-06   1.698254E-04
  b2u  30   3.118911E-04   4.725801E-04   6.904770E-05  -2.837800E-06   3.743198E-05  -2.298821E-05   3.326686E-06   1.182390E-04

               b2u  25        b2u  26        b2u  27        b2u  28        b2u  29        b2u  30
  b2u   2  -4.592190E-04   8.151093E-04   1.073483E-04   1.898801E-03  -8.672323E-04   1.977722E-03
  b2u   3   4.060676E-04  -8.254249E-05  -3.739060E-03  -1.344743E-04   3.723739E-05  -3.411699E-04
  b2u   4  -4.530977E-03  -3.134796E-04   1.611774E-04  -8.890708E-04   1.378809E-03   2.693170E-03
  b2u   5   9.550777E-05   6.665192E-05   1.739175E-06  -1.499288E-04  -9.085320E-07  -4.989985E-05
  b2u   6   1.097853E-06   5.570177E-07   2.665746E-05  -2.767767E-06   5.694322E-07   1.094255E-06
  b2u   7   6.360636E-05   9.787091E-05   1.961022E-06  -2.043223E-04  -2.590511E-05  -6.424787E-05
  b2u   8   2.144563E-04   1.143952E-04   3.670047E-06  -2.984844E-04   7.652342E-05  -7.495035E-05
  b2u   9  -2.690154E-05   2.233482E-05   1.675963E-06  -6.067344E-05  -7.082016E-05   3.306790E-05
  b2u  10  -4.578022E-04  -1.719749E-04  -9.400781E-06   5.414050E-04  -3.570196E-04  -1.008066E-04
  b2u  11  -5.022532E-06   9.843038E-07  -9.977320E-05   6.035574E-06  -5.468362E-06  -3.598135E-06
  b2u  12   1.044471E-05  -1.993055E-04   1.328040E-07   4.294891E-04   7.707921E-05   2.872882E-04
  b2u  13  -8.540919E-05   3.372248E-04   1.804847E-05   3.583362E-05   9.209990E-05   4.324468E-04
  b2u  14   1.874341E-05   5.969346E-05   2.400066E-06   5.165302E-05   2.033654E-04   1.117414E-04
  b2u  15  -4.235689E-05  -1.579412E-04  -7.974097E-06   2.580517E-04  -1.669982E-04  -8.412039E-05
  b2u  16  -1.389154E-06   7.987610E-06  -2.926526E-04  -8.134338E-07  -3.072430E-06  -3.717327E-06
  b2u  17  -4.114265E-05   2.095416E-04   1.013539E-05   7.255827E-04  -2.963950E-04   3.118911E-04
  b2u  18  -4.786202E-04   9.014372E-04   2.806808E-05   1.592533E-04  -7.902400E-05   4.725801E-04
  b2u  19   4.648972E-04   4.717538E-04   1.596174E-05  -4.474437E-04  -8.693093E-05   6.904770E-05
  b2u  20  -1.501125E-05  -1.686727E-05   2.172903E-04   3.671932E-06  -5.878780E-07  -2.837800E-06
  b2u  21   6.050002E-04   2.821271E-04   1.276448E-05  -9.428723E-05   6.707356E-05   3.743198E-05
  b2u  22  -3.361195E-04  -2.801412E-04  -6.269906E-06  -9.128299E-05   2.117104E-04  -2.298821E-05
  b2u  23   1.806466E-06  -1.130522E-05   5.154712E-04   6.087137E-08   1.862107E-06   3.326686E-06
  b2u  24  -5.096719E-05  -1.381071E-04  -4.127882E-06  -6.383456E-05   1.698254E-04   1.182390E-04
  b2u  25   1.214368E-03   1.654835E-04   4.330182E-06  -3.855106E-04   2.393912E-04  -9.658369E-05
  b2u  26   1.654835E-04   1.119105E-03   1.305227E-05  -2.068259E-04  -5.215139E-05   3.697338E-04
  b2u  27   4.330182E-06   1.305227E-05   6.500145E-04  -4.776603E-06  -5.693606E-07   1.024577E-05
  b2u  28  -3.855106E-04  -2.068259E-04  -4.776603E-06   5.846546E-04  -1.077046E-04   2.176792E-05
  b2u  29   2.393912E-04  -5.215139E-05  -5.693606E-07  -1.077046E-04   3.481305E-04  -1.248657E-05
  b2u  30  -9.658369E-05   3.697338E-04   1.024577E-05   2.176792E-05  -1.248657E-05   4.459449E-04

Natural orbital populations,block 3
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98501368     1.97654084     1.97444752     0.01416315     0.01080706     0.00653822     0.00483177
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00334541     0.00211814     0.00126179     0.00111704     0.00079783     0.00068571     0.00044325     0.00032603
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00021581     0.00018373     0.00013569     0.00009901     0.00003171     0.00002931     0.00002053     0.00001213
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00000684     0.00000441     0.00000339     0.00000145     0.00000015     0.00000007

          modens reordered block   1

               b1g   1        b1g   2        b1g   3        b1g   4        b1g   5        b1g   6        b1g   7        b1g   8
  b1g   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b1g   2    0.00000        1.98076       2.157639E-03  -1.143675E-03  -2.720614E-03  -1.899643E-04  -1.274952E-06   3.042965E-03
  b1g   3    0.00000       2.157639E-03    1.97425      -1.795952E-03  -1.209950E-03  -3.636839E-03   5.137424E-05  -1.310083E-03
  b1g   4    0.00000      -1.143675E-03  -1.795952E-03   2.273661E-04   4.081205E-04   3.702207E-04   2.006371E-07  -1.881267E-04
  b1g   5    0.00000      -2.720614E-03  -1.209950E-03   4.081205E-04   8.719413E-04   5.436026E-04  -3.132350E-07  -5.689474E-04
  b1g   6    0.00000      -1.899643E-04  -3.636839E-03   3.702207E-04   5.436026E-04   7.661203E-04   8.527439E-07  -2.649128E-05
  b1g   7    0.00000      -1.274952E-06   5.137424E-05   2.006371E-07  -3.132350E-07   8.527439E-07   9.366009E-06   9.419900E-07
  b1g   8    0.00000       3.042965E-03  -1.310083E-03  -1.881267E-04  -5.689474E-04  -2.649128E-05   9.419900E-07   7.414208E-04
  b1g   9    0.00000       2.777005E-03  -1.009061E-03  -6.100199E-04  -1.604035E-03  -6.596223E-04   2.518592E-06   1.097991E-03
  b1g  10    0.00000       1.461558E-03  -5.882069E-03   5.938509E-04   6.273196E-04   1.564616E-03   3.044490E-06   5.310268E-04
  b1g  11    0.00000      -7.723439E-05  -2.845913E-03  -3.643702E-04  -1.336642E-03  -2.112059E-04   6.270593E-06   1.002572E-03
  b1g  12    0.00000      -8.741264E-05   7.024197E-04   4.948156E-06   1.880591E-05   7.137012E-07   9.015717E-05  -1.558844E-05
  b1g  13    0.00000      -3.019290E-03   1.153362E-03   1.648094E-04   5.278271E-04  -9.837576E-05  -1.116592E-06  -9.604513E-04
  b1g  14    0.00000       2.079093E-03  -1.225023E-03   3.868750E-04   4.585178E-04   1.065025E-03   1.618164E-06   4.629423E-04
  b1g  15    0.00000      -3.745693E-04   1.608019E-03  -4.377090E-06  -1.027676E-05  -1.098631E-05   1.924989E-04   3.167136E-06
  b1g  16    0.00000      -4.256643E-05  -1.506744E-03   5.045870E-04   8.623018E-04   1.060325E-03   2.305312E-06  -1.657034E-04
  b1g  17    0.00000       1.334431E-03   8.714800E-04   9.907280E-05   4.231758E-04   9.369938E-05  -1.438922E-06  -3.263921E-04
  b1g  18    0.00000       3.741809E-04  -8.317724E-04   1.214914E-06   2.855712E-06   3.685282E-06  -6.467273E-05  -7.709771E-07
  b1g  19    0.00000       1.294701E-03  -4.223800E-04   1.275433E-04   1.889327E-04   2.863622E-04   7.176621E-07  -4.481922E-05
  b1g  20    0.00000       2.092961E-03   1.546045E-03  -2.102887E-06  -9.364546E-05   2.711192E-04   7.298779E-07   5.820548E-04
  b1g  21    0.00000       2.598671E-04  -2.605741E-03  -1.239553E-04  -1.273531E-04  -2.359560E-04  -6.925164E-07   4.694269E-05
  b1g  22    0.00000       2.230315E-04  -1.410483E-04   3.930271E-07   1.718774E-06   7.766001E-07  -1.386451E-06  -2.232613E-06
  b1g  23    0.00000       2.749230E-03  -3.217233E-03  -1.471779E-04  -4.277633E-04  -8.209812E-05   1.546525E-07   4.198386E-04
  b1g  24    0.00000      -1.121827E-03  -2.650949E-04   1.524002E-04   3.083059E-04   2.998085E-04   9.110782E-08  -2.434497E-04
  b1g  25    0.00000      -1.363656E-04   1.660945E-05  -2.231749E-07  -1.174078E-06   1.951706E-07  -3.909981E-05   1.136847E-06
  b1g  26    0.00000       2.308550E-03   1.319098E-03  -1.828019E-06   3.483847E-05  -9.289817E-06  -3.546460E-07  -4.599194E-05
  b1g  27    0.00000       9.861759E-04   3.676093E-03   6.914285E-06   8.207657E-05  -1.237756E-05  -5.073677E-08  -9.086323E-05
  b1g  28    0.00000       2.710102E-03  -2.758143E-03   1.026959E-04   1.795725E-04   2.143080E-04   1.967321E-07  -3.055349E-05
  b1g  29    0.00000      -1.325332E-04  -2.309140E-05  -1.600467E-07  -8.135558E-07  -1.330918E-07  -4.911819E-06   9.301784E-07
  b1g  30    0.00000      -5.132839E-04   6.224710E-04  -1.081791E-05  -4.680680E-05   4.207945E-05   3.253614E-07   8.619009E-05

               b1g   9        b1g  10        b1g  11        b1g  12        b1g  13        b1g  14        b1g  15        b1g  16
  b1g   2   2.777005E-03   1.461558E-03  -7.723439E-05  -8.741264E-05  -3.019290E-03   2.079093E-03  -3.745693E-04  -4.256643E-05
  b1g   3  -1.009061E-03  -5.882069E-03  -2.845913E-03   7.024197E-04   1.153362E-03  -1.225023E-03   1.608019E-03  -1.506744E-03
  b1g   4  -6.100199E-04   5.938509E-04  -3.643702E-04   4.948156E-06   1.648094E-04   3.868750E-04  -4.377090E-06   5.045870E-04
  b1g   5  -1.604035E-03   6.273196E-04  -1.336642E-03   1.880591E-05   5.278271E-04   4.585178E-04  -1.027676E-05   8.623018E-04
  b1g   6  -6.596223E-04   1.564616E-03  -2.112059E-04   7.137012E-07  -9.837576E-05   1.065025E-03  -1.098631E-05   1.060325E-03
  b1g   7   2.518592E-06   3.044490E-06   6.270593E-06   9.015717E-05  -1.116592E-06   1.618164E-06   1.924989E-04   2.305312E-06
  b1g   8   1.097991E-03   5.310268E-04   1.002572E-03  -1.558844E-05  -9.604513E-04   4.629423E-04   3.167136E-06  -1.657034E-04
  b1g   9   4.068519E-03  -1.781668E-04   4.829910E-03  -6.264284E-05  -7.859968E-04  -4.894580E-04   4.353588E-05  -1.137344E-03
  b1g  10  -1.781668E-04   3.968582E-03   1.094609E-03  -1.841754E-05  -9.928644E-04   2.831290E-03  -9.913675E-06   2.570273E-03
  b1g  11   4.829910E-03   1.094609E-03   7.599685E-03  -7.517858E-05  -5.883967E-04   2.088649E-04   1.245821E-04  -7.726040E-05
  b1g  12  -6.264284E-05  -1.841754E-05  -7.517858E-05   9.752421E-04   9.859901E-06  -8.630999E-06   2.336879E-03   1.279533E-05
  b1g  13  -7.859968E-04  -9.928644E-04  -5.883967E-04   9.859901E-06   1.531530E-03  -1.042296E-03  -1.295590E-06  -2.054428E-04
  b1g  14  -4.894580E-04   2.831290E-03   2.088649E-04  -8.630999E-06  -1.042296E-03   2.305739E-03  -1.921755E-05   2.014468E-03
  b1g  15   4.353588E-05  -9.913675E-06   1.245821E-04   2.336879E-03  -1.295590E-06  -1.921755E-05   6.294956E-03   2.408024E-05
  b1g  16  -1.137344E-03   2.570273E-03  -7.726040E-05   1.279533E-05  -2.054428E-04   2.014468E-03   2.408024E-05   2.639921E-03
  b1g  17  -2.175294E-03  -6.505539E-04  -4.595573E-03   6.535596E-05   2.783584E-04  -3.390982E-04  -3.005998E-05  -4.651704E-04
  b1g  18  -1.415037E-05   2.750742E-06  -4.805481E-05  -1.098347E-03   7.126756E-07   8.711398E-06  -3.691251E-03  -1.765725E-05
  b1g  19   5.287888E-05   9.121148E-04   8.613299E-04  -1.256807E-05  -3.891353E-05   6.843051E-04   7.286657E-07   9.871545E-04
  b1g  20  -1.759249E-04   1.123146E-03  -5.339845E-04   6.672225E-06  -1.218863E-03   1.295316E-03  -6.229121E-06   7.374853E-04
  b1g  21  -9.822521E-05  -3.170800E-04  -3.067976E-04   1.646094E-06  -2.506979E-04  -3.957547E-05  -6.538895E-06   2.147200E-04
  b1g  22  -7.021945E-06  -3.459917E-07  -1.740535E-05  -2.138161E-04   3.230232E-06   1.690624E-06  -1.093857E-03  -2.592059E-06
  b1g  23   7.627459E-04   1.907924E-05   7.437275E-05  -7.108281E-06  -4.147910E-04  -8.829581E-05  -1.864820E-05  -6.591850E-04
  b1g  24  -5.905027E-04   5.240963E-04  -8.310008E-04   1.361884E-05   3.748070E-04   2.711608E-04  -2.349626E-06   7.243276E-04
  b1g  25   4.100416E-06   1.863263E-06   1.316160E-06  -2.352991E-04  -6.011163E-07   6.551067E-07  -6.718187E-05  -4.777203E-07
  b1g  26  -1.723490E-04  -4.859241E-05  -4.094666E-04   3.448902E-06   5.312507E-06   9.243609E-05  -1.028009E-05   1.543089E-04
  b1g  27  -3.149368E-04  -1.236006E-05  -3.974419E-04   8.290972E-06   2.444553E-05   1.396418E-04   5.466221E-06   2.654053E-04
  b1g  28  -3.096260E-04   4.669010E-04  -2.445658E-04   1.722163E-06  -4.610220E-05   3.507872E-04  -1.008675E-05   4.037304E-04
  b1g  29   3.183734E-06   6.450306E-07   5.397291E-06   2.373204E-05  -9.673524E-07  -3.916589E-07   2.452761E-04   7.009149E-07
  b1g  30   1.194587E-04   2.410253E-04   2.055178E-04  -1.483992E-06  -1.345353E-04   1.823451E-04   4.588275E-06   2.578928E-04

               b1g  17        b1g  18        b1g  19        b1g  20        b1g  21        b1g  22        b1g  23        b1g  24
  b1g   2   1.334431E-03   3.741809E-04   1.294701E-03   2.092961E-03   2.598671E-04   2.230315E-04   2.749230E-03  -1.121827E-03
  b1g   3   8.714800E-04  -8.317724E-04  -4.223800E-04   1.546045E-03  -2.605741E-03  -1.410483E-04  -3.217233E-03  -2.650949E-04
  b1g   4   9.907280E-05   1.214914E-06   1.275433E-04  -2.102887E-06  -1.239553E-04   3.930271E-07  -1.471779E-04   1.524002E-04
  b1g   5   4.231758E-04   2.855712E-06   1.889327E-04  -9.364546E-05  -1.273531E-04   1.718774E-06  -4.277633E-04   3.083059E-04
  b1g   6   9.369938E-05   3.685282E-06   2.863622E-04   2.711192E-04  -2.359560E-04   7.766001E-07  -8.209812E-05   2.998085E-04
  b1g   7  -1.438922E-06  -6.467273E-05   7.176621E-07   7.298779E-07  -6.925164E-07  -1.386451E-06   1.546525E-07   9.110782E-08
  b1g   8  -3.263921E-04  -7.709771E-07  -4.481922E-05   5.820548E-04   4.694269E-05  -2.232613E-06   4.198386E-04  -2.434497E-04
  b1g   9  -2.175294E-03  -1.415037E-05   5.287888E-05  -1.759249E-04  -9.822521E-05  -7.021945E-06   7.627459E-04  -5.905027E-04
  b1g  10  -6.505539E-04   2.750742E-06   9.121148E-04   1.123146E-03  -3.170800E-04  -3.459917E-07   1.907924E-05   5.240963E-04
  b1g  11  -4.595573E-03  -4.805481E-05   8.613299E-04  -5.339845E-04  -3.067976E-04  -1.740535E-05   7.437275E-05  -8.310008E-04
  b1g  12   6.535596E-05  -1.098347E-03  -1.256807E-05   6.672225E-06   1.646094E-06  -2.138161E-04  -7.108281E-06   1.361884E-05
  b1g  13   2.783584E-04   7.126756E-07  -3.891353E-05  -1.218863E-03  -2.506979E-04   3.230232E-06  -4.147910E-04   3.748070E-04
  b1g  14  -3.390982E-04   8.711398E-06   6.843051E-04   1.295316E-03  -3.957547E-05   1.690624E-06  -8.829581E-05   2.711608E-04
  b1g  15  -3.005998E-05  -3.691251E-03   7.286657E-07  -6.229121E-06  -6.538895E-06  -1.093857E-03  -1.864820E-05  -2.349626E-06
  b1g  16  -4.651704E-04  -1.765725E-05   9.871545E-04   7.374853E-04   2.147200E-04  -2.592059E-06  -6.591850E-04   7.243276E-04
  b1g  17   3.959371E-03   1.540398E-06  -9.710688E-04   3.777289E-04   2.514500E-06   6.136961E-06   9.522166E-04   1.000204E-03
  b1g  18   1.540398E-06   2.920413E-03   3.607832E-06   2.524629E-07   1.355752E-06   1.191893E-03   8.449723E-06  -8.652148E-06
  b1g  19  -9.710688E-04   3.607832E-06   8.550143E-04   9.423613E-05   6.180002E-05   1.975895E-06  -3.128629E-04  -3.043508E-05
  b1g  20   3.777289E-04   2.524629E-07   9.423613E-05   1.568425E-03   3.067947E-04  -1.240485E-06   3.457100E-04  -6.148074E-05
  b1g  21   2.514500E-06   1.355752E-06   6.180002E-05   3.067947E-04   7.600274E-04   1.567580E-06  -1.010943E-04   3.215170E-05
  b1g  22   6.136961E-06   1.191893E-03   1.975895E-06  -1.240485E-06   1.567580E-06   6.561087E-04   3.192275E-06   1.245256E-06
  b1g  23   9.522166E-04   8.449723E-06  -3.128629E-04   3.457100E-04  -1.010943E-04   3.192275E-06   1.459546E-03   1.071973E-04
  b1g  24   1.000204E-03  -8.652148E-06  -3.043508E-05  -6.148074E-05   3.215170E-05   1.245256E-06   1.071973E-04   1.052014E-03
  b1g  25  -2.418851E-06  -5.994063E-04  -7.228120E-07   2.283225E-07   7.063926E-07  -4.621426E-04   2.010225E-06   3.835782E-06
  b1g  26   5.479846E-04   2.728670E-06  -2.326836E-04   4.081203E-04   2.647977E-04   2.398177E-06   2.995409E-04   3.076623E-04
  b1g  27   5.612183E-05  -4.645254E-06   2.389993E-04   1.804252E-04   1.395170E-04  -1.040924E-06  -2.947989E-04   1.597067E-05
  b1g  28   1.337246E-04   6.794980E-06   2.290010E-04   8.986126E-05  -5.875769E-05   3.866168E-06   1.947520E-04   6.583687E-05
  b1g  29  -1.462828E-06  -4.057249E-04  -1.115134E-06   2.785381E-07   7.046260E-07  -3.097479E-04   5.902745E-07   1.407690E-06
  b1g  30  -1.675194E-04  -3.123586E-06   1.327610E-04   2.459871E-05  -2.225053E-06  -6.433431E-07  -1.900570E-04   1.204926E-04

               b1g  25        b1g  26        b1g  27        b1g  28        b1g  29        b1g  30
  b1g   2  -1.363656E-04   2.308550E-03   9.861759E-04   2.710102E-03  -1.325332E-04  -5.132839E-04
  b1g   3   1.660945E-05   1.319098E-03   3.676093E-03  -2.758143E-03  -2.309140E-05   6.224710E-04
  b1g   4  -2.231749E-07  -1.828019E-06   6.914285E-06   1.026959E-04  -1.600467E-07  -1.081791E-05
  b1g   5  -1.174078E-06   3.483847E-05   8.207657E-05   1.795725E-04  -8.135558E-07  -4.680680E-05
  b1g   6   1.951706E-07  -9.289817E-06  -1.237756E-05   2.143080E-04  -1.330918E-07   4.207945E-05
  b1g   7  -3.909981E-05  -3.546460E-07  -5.073677E-08   1.967321E-07  -4.911819E-06   3.253614E-07
  b1g   8   1.136847E-06  -4.599194E-05  -9.086323E-05  -3.055349E-05   9.301784E-07   8.619009E-05
  b1g   9   4.100416E-06  -1.723490E-04  -3.149368E-04  -3.096260E-04   3.183734E-06   1.194587E-04
  b1g  10   1.863263E-06  -4.859241E-05  -1.236006E-05   4.669010E-04   6.450306E-07   2.410253E-04
  b1g  11   1.316160E-06  -4.094666E-04  -3.974419E-04  -2.445658E-04   5.397291E-06   2.055178E-04
  b1g  12  -2.352991E-04   3.448902E-06   8.290972E-06   1.722163E-06   2.373204E-05  -1.483992E-06
  b1g  13  -6.011163E-07   5.312507E-06   2.444553E-05  -4.610220E-05  -9.673524E-07  -1.345353E-04
  b1g  14   6.551067E-07   9.243609E-05   1.396418E-04   3.507872E-04  -3.916589E-07   1.823451E-04
  b1g  15  -6.718187E-05  -1.028009E-05   5.466221E-06  -1.008675E-05   2.452761E-04   4.588275E-06
  b1g  16  -4.777203E-07   1.543089E-04   2.654053E-04   4.037304E-04   7.009149E-07   2.578928E-04
  b1g  17  -2.418851E-06   5.479846E-04   5.612183E-05   1.337246E-04  -1.462828E-06  -1.675194E-04
  b1g  18  -5.994063E-04   2.728670E-06  -4.645254E-06   6.794980E-06  -4.057249E-04  -3.123586E-06
  b1g  19  -7.228120E-07  -2.326836E-04   2.389993E-04   2.290010E-04  -1.115134E-06   1.327610E-04
  b1g  20   2.283225E-07   4.081203E-04   1.804252E-04   8.986126E-05   2.785381E-07   2.459871E-05
  b1g  21   7.063926E-07   2.647977E-04   1.395170E-04  -5.875769E-05   7.046260E-07  -2.225053E-06
  b1g  22  -4.621426E-04   2.398177E-06  -1.040924E-06   3.866168E-06  -3.097479E-04  -6.433431E-07
  b1g  23   2.010225E-06   2.995409E-04  -2.947989E-04   1.947520E-04   5.902745E-07  -1.900570E-04
  b1g  24   3.835782E-06   3.076623E-04   1.597067E-05   6.583687E-05   1.407690E-06   1.204926E-04
  b1g  25   6.177671E-04   9.718934E-07   3.331140E-08  -9.369835E-07   2.877856E-04   4.783490E-07
  b1g  26   9.718934E-07   8.614552E-04   2.123012E-05  -1.282190E-04   4.342063E-07  -2.274911E-04
  b1g  27   3.331140E-08   2.123012E-05   4.032138E-04  -8.917714E-05  -1.784115E-07   2.153837E-05
  b1g  28  -9.369835E-07  -1.282190E-04  -8.917714E-05   4.955271E-04  -1.115322E-06   1.169990E-04
  b1g  29   2.877856E-04   4.342063E-07  -1.784115E-07  -1.115322E-06   3.352313E-04   2.370356E-07
  b1g  30   4.783490E-07  -2.274911E-04   2.153837E-05   1.169990E-04   2.370356E-07   3.233849E-04

Natural orbital populations,block 4
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98143872     1.97365893     0.01464936     0.01033212     0.00974188     0.00465037     0.00234683
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00176124     0.00149489     0.00092388     0.00077183     0.00046834     0.00030105     0.00021663     0.00014225
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00010512     0.00008986     0.00006453     0.00005214     0.00002055     0.00001889     0.00000848     0.00000652
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00000294     0.00000131     0.00000107     0.00000022     0.00000003     0.00000002

          modens reordered block   1

               b1u   1        b1u   2        b1u   3        b1u   4        b1u   5        b1u   6        b1u   7        b1u   8
  b1u   1    1.92184      -5.654278E-02  -1.427397E-02  -5.740627E-04   9.111604E-04  -6.270175E-04  -2.064448E-03   7.852898E-04
  b1u   2  -5.654278E-02   0.504623       0.145292       2.102762E-03   1.376213E-03   1.435897E-03  -3.184053E-03  -1.869028E-03
  b1u   3  -1.427397E-02   0.145292       4.569026E-02   2.637435E-04   7.906878E-04  -2.342026E-05  -9.174030E-06  -4.952914E-05
  b1u   4  -5.740627E-04   2.102762E-03   2.637435E-04   7.885752E-04   8.761797E-05   7.194570E-04   3.101709E-06  -7.275367E-04
  b1u   5   9.111604E-04   1.376213E-03   7.906878E-04   8.761797E-05   2.080911E-03   2.878684E-05  -1.666074E-03  -2.569094E-05
  b1u   6  -6.270175E-04   1.435897E-03  -2.342026E-05   7.194570E-04   2.878684E-05   1.030682E-03  -2.495141E-05  -5.370033E-04
  b1u   7  -2.064448E-03  -3.184053E-03  -9.174030E-06   3.101709E-06  -1.666074E-03  -2.495141E-05   2.291500E-03   2.676718E-05
  b1u   8   7.852898E-04  -1.869028E-03  -4.952914E-05  -7.275367E-04  -2.569094E-05  -5.370033E-04   2.676718E-05   1.553936E-03
  b1u   9   5.590906E-03   2.550982E-03  -1.666685E-03  -1.762133E-04  -5.520416E-04  -9.781758E-06  -8.325231E-04   1.176793E-06
  b1u  10  -4.120753E-04  -1.884564E-04  -4.321041E-04   6.499861E-04   1.113933E-05   1.184427E-03  -7.822277E-06  -1.072799E-04
  b1u  11   1.918375E-04   2.018079E-03   9.580461E-04  -5.867777E-04  -2.869881E-05  -2.734372E-04   1.231312E-05   1.032829E-03
  b1u  12   3.296916E-03   2.324214E-04   4.287765E-04   5.227872E-05   1.763362E-03   1.497800E-05  -1.433669E-03  -1.183319E-05
  b1u  13  -7.676292E-05   1.132085E-03   3.811573E-04  -1.058193E-04  -7.534225E-06  -3.955133E-04   1.757267E-05  -9.907196E-04
  b1u  14   1.379730E-04  -1.128394E-03  -2.147013E-04  -2.142240E-04   3.505719E-06  -7.245184E-04  -3.352794E-06   3.936063E-04
  b1u  15   2.418074E-03   5.194606E-04   5.302417E-05   4.962880E-06   2.614098E-04   4.045556E-06  -6.596955E-04   2.736323E-07
  b1u  16  -2.021641E-05   6.635832E-04   1.572700E-04   3.886052E-05  -2.740318E-06   6.572725E-05  -3.759551E-06  -1.770704E-04

               b1u   9        b1u  10        b1u  11        b1u  12        b1u  13        b1u  14        b1u  15        b1u  16
  b1u   1   5.590906E-03  -4.120753E-04   1.918375E-04   3.296916E-03  -7.676292E-05   1.379730E-04   2.418074E-03  -2.021641E-05
  b1u   2   2.550982E-03  -1.884564E-04   2.018079E-03   2.324214E-04   1.132085E-03  -1.128394E-03   5.194606E-04   6.635832E-04
  b1u   3  -1.666685E-03  -4.321041E-04   9.580461E-04   4.287765E-04   3.811573E-04  -2.147013E-04   5.302417E-05   1.572700E-04
  b1u   4  -1.762133E-04   6.499861E-04  -5.867777E-04   5.227872E-05  -1.058193E-04  -2.142240E-04   4.962880E-06   3.886052E-05
  b1u   5  -5.520416E-04   1.113933E-05  -2.869881E-05   1.763362E-03  -7.534225E-06   3.505719E-06   2.614098E-04  -2.740318E-06
  b1u   6  -9.781758E-06   1.184427E-03  -2.734372E-04   1.497800E-05  -3.955133E-04  -7.245184E-04   4.045556E-06   6.572725E-05
  b1u   7  -8.325231E-04  -7.822277E-06   1.231312E-05  -1.433669E-03   1.757267E-05  -3.352794E-06  -6.596955E-04  -3.759551E-06
  b1u   8   1.176793E-06  -1.072799E-04   1.032829E-03  -1.183319E-05  -9.907196E-04   3.936063E-04   2.736323E-07  -1.770704E-04
  b1u   9   2.870328E-03  -7.809108E-07   1.805899E-05  -4.124638E-04  -3.481449E-05  -6.998250E-08   1.179635E-04   1.033229E-05
  b1u  10  -7.809108E-07   1.710060E-03   7.889213E-05   3.332231E-06  -1.130813E-03  -9.366125E-04   2.294794E-07   1.828474E-04
  b1u  11   1.805899E-05   7.889213E-05   1.269296E-03  -1.046771E-05  -7.258348E-04  -1.515483E-04   1.921527E-06   9.971640E-05
  b1u  12  -4.124638E-04   3.332231E-06  -1.046771E-05   2.038906E-03  -4.371284E-06  -2.214831E-06   2.773223E-04   5.420986E-07
  b1u  13  -3.481449E-05  -1.130813E-03  -7.258348E-04  -4.371284E-06   2.133576E-03  -3.985613E-05  -6.404403E-06  -4.323384E-05
  b1u  14  -6.998250E-08  -9.366125E-04  -1.515483E-04  -2.214831E-06  -3.985613E-05   1.265484E-03   3.753398E-06  -3.016219E-04
  b1u  15   1.179635E-04   2.294794E-07   1.921527E-06   2.773223E-04  -6.404403E-06   3.753398E-06   6.215033E-04   9.235102E-09
  b1u  16   1.033229E-05   1.828474E-04   9.971640E-05   5.420986E-07  -4.323384E-05  -3.016219E-04   9.235102E-09   4.753870E-04

Natural orbital populations,block 5
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.92434045     0.54435183     0.00613129     0.00545894     0.00411645     0.00380432     0.00142090     0.00076802
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00063341     0.00047996     0.00035207     0.00021179     0.00012527     0.00004824     0.00003700     0.00000792

          modens reordered block   1

               b2g   1        b2g   2        b2g   3        b2g   4        b2g   5        b2g   6        b2g   7        b2g   8
  b2g   1    1.47671       2.676286E-04   3.377108E-03  -1.240524E-02   3.365858E-04   5.400953E-03  -2.281449E-03  -4.952640E-03
  b2g   2   2.676286E-04   1.098609E-03  -2.156166E-05   6.790324E-05  -8.970075E-04  -3.107048E-05  -1.132116E-03  -1.971818E-04
  b2g   3   3.377108E-03  -2.156166E-05   1.302042E-03  -4.649080E-03   3.220795E-05   1.883383E-03   2.178126E-05  -2.019260E-06
  b2g   4  -1.240524E-02   6.790324E-05  -4.649080E-03   1.709943E-02  -1.097274E-04  -7.363240E-03  -6.341491E-05   2.313103E-05
  b2g   5   3.365858E-04  -8.970075E-04   3.220795E-05  -1.097274E-04   1.296604E-03   4.892806E-05   1.342919E-03  -5.745967E-04
  b2g   6   5.400953E-03  -3.107048E-05   1.883383E-03  -7.363240E-03   4.892806E-05   3.633180E-03   2.623573E-05  -1.867816E-05
  b2g   7  -2.281449E-03  -1.132116E-03   2.178126E-05  -6.341491E-05   1.342919E-03   2.623573E-05   1.707287E-03  -2.508801E-04
  b2g   8  -4.952640E-03  -1.971818E-04  -2.019260E-06   2.313103E-05  -5.745967E-04  -1.867816E-05  -2.508801E-04   1.395323E-03
  b2g   9   3.084358E-03  -4.640527E-05   1.127756E-03  -4.755750E-03   4.021946E-05   2.702894E-03   3.606786E-05   1.112289E-05
  b2g  10  -7.286684E-03  -1.021050E-03  -1.213309E-05   8.624864E-05   5.326241E-04  -6.310967E-05   8.636285E-04   7.727105E-04
  b2g  11   4.352472E-04   9.873097E-06   5.406829E-04  -2.037064E-03   1.535511E-06   8.053205E-04  -3.211015E-06  -9.246129E-06
  b2g  12   2.271957E-03  -4.078961E-04   6.446689E-06  -2.350729E-05   1.113598E-03   1.340669E-05   1.083260E-03  -1.213049E-03
  b2g  13   6.417102E-04  -3.928613E-06  -5.747841E-05  -2.232564E-04  -2.123146E-07   6.962367E-04  -1.871710E-06  -1.380115E-06
  b2g  14  -2.521424E-03  -5.047479E-05  -2.236045E-06   1.272244E-05   8.762959E-05  -9.179749E-06  -3.573423E-05   2.759088E-04
  b2g  15  -1.103352E-03  -5.843256E-05  -6.554122E-06   2.542210E-05   9.487995E-05  -1.150970E-05   3.616337E-04   2.789594E-04
  b2g  16  -3.004633E-05   1.253304E-07  -3.638244E-05   9.492906E-05  -1.121497E-06   1.707503E-05  -6.038977E-07   5.774833E-07

               b2g   9        b2g  10        b2g  11        b2g  12        b2g  13        b2g  14        b2g  15        b2g  16
  b2g   1   3.084358E-03  -7.286684E-03   4.352472E-04   2.271957E-03   6.417102E-04  -2.521424E-03  -1.103352E-03  -3.004633E-05
  b2g   2  -4.640527E-05  -1.021050E-03   9.873097E-06  -4.078961E-04  -3.928613E-06  -5.047479E-05  -5.843256E-05   1.253304E-07
  b2g   3   1.127756E-03  -1.213309E-05   5.406829E-04   6.446689E-06  -5.747841E-05  -2.236045E-06  -6.554122E-06  -3.638244E-05
  b2g   4  -4.755750E-03   8.624864E-05  -2.037064E-03  -2.350729E-05  -2.232564E-04   1.272244E-05   2.542210E-05   9.492906E-05
  b2g   5   4.021946E-05   5.326241E-04   1.535511E-06   1.113598E-03  -2.123146E-07   8.762959E-05   9.487995E-05  -1.121497E-06
  b2g   6   2.702894E-03  -6.310967E-05   8.053205E-04   1.340669E-05   6.962367E-04  -9.179749E-06  -1.150970E-05   1.707503E-05
  b2g   7   3.606786E-05   8.636285E-04  -3.211015E-06   1.083260E-03  -1.871710E-06  -3.573423E-05   3.616337E-04  -6.038977E-07
  b2g   8   1.112289E-05   7.727105E-04  -9.246129E-06  -1.213049E-03  -1.380115E-06   2.759088E-04   2.789594E-04   5.774833E-07
  b2g   9   2.306809E-03  -1.338939E-05   5.069341E-04  -6.030365E-07   9.893316E-04  -1.364899E-06  -9.240239E-06   1.383012E-04
  b2g  10  -1.338939E-05   1.607509E-03  -2.133630E-05  -1.056209E-04  -3.123384E-05   2.463704E-04  -7.701314E-05  -2.749534E-06
  b2g  11   5.069341E-04  -2.133630E-05   5.885531E-04   1.918691E-06  -2.385077E-04  -5.117961E-06  -3.336074E-06  -6.358357E-05
  b2g  12  -6.030365E-07  -1.056209E-04   1.918691E-06   1.695961E-03  -3.343894E-06  -2.488781E-04  -4.664817E-05  -1.967112E-06
  b2g  13   9.893316E-04  -3.123384E-05  -2.385077E-04  -3.343894E-06   1.134416E-03  -2.368729E-06   1.042713E-06   2.892856E-04
  b2g  14  -1.364899E-06   2.463704E-04  -5.117961E-06  -2.488781E-04  -2.368729E-06   6.555690E-04   2.023898E-04   2.081020E-07
  b2g  15  -9.240239E-06  -7.701314E-05  -3.336074E-06  -4.664817E-05   1.042713E-06   2.023898E-04   9.309410E-04   5.413030E-07
  b2g  16   1.383012E-04  -2.749534E-06  -6.358357E-05  -1.967112E-06   2.892856E-04   2.081020E-07   5.413030E-07   3.574604E-04

Natural orbital populations,block 6
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.47691536     0.02318667     0.00504956     0.00316090     0.00229759     0.00113732     0.00060528     0.00042332
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00030560     0.00022606     0.00007184     0.00006262     0.00005820     0.00001201     0.00000556     0.00000147

          modens reordered block   1

               b3g   1        b3g   2        b3g   3        b3g   4        b3g   5        b3g   6        b3g   7        b3g   8
  b3g   1    1.47149       1.565969E-03   1.801675E-04  -6.807549E-04  -2.525132E-03  -4.678978E-03  -6.758455E-03   1.034557E-03
  b3g   2   1.565969E-03   1.053659E-03   8.987492E-04   2.672160E-05   1.123416E-03   1.566439E-04   9.655771E-04  -1.640535E-06
  b3g   3   1.801675E-04   8.987492E-04   1.311360E-03   1.314591E-05   1.356434E-03  -5.840134E-04   5.304897E-04   1.961308E-06
  b3g   4  -6.807549E-04   2.672160E-05   1.314591E-05   2.826937E-03   5.507689E-06   2.024096E-06   1.464474E-05  -2.781975E-03
  b3g   5  -2.525132E-03   1.123416E-03   1.356434E-03   5.507689E-06   1.726131E-03  -2.607673E-04   8.556233E-04   1.366197E-05
  b3g   6  -4.678978E-03   1.566439E-04  -5.840134E-04   2.024096E-06  -2.607673E-04   1.384045E-03   7.375763E-04   1.720726E-06
  b3g   7  -6.758455E-03   9.655771E-04   5.304897E-04   1.464474E-05   8.556233E-04   7.375763E-04   1.548102E-03  -9.586538E-07
  b3g   8   1.034557E-03  -1.640535E-06   1.961308E-06  -2.781975E-03   1.366197E-05   1.720726E-06  -9.586538E-07   3.510943E-03
  b3g   9   1.954843E-03   4.281136E-04   1.124761E-03   7.120358E-06   1.099865E-03  -1.205735E-03  -7.589028E-05  -1.947877E-06
  b3g  10  -2.362282E-03   2.721656E-05   7.896669E-05  -2.708403E-06  -5.544789E-05   2.576222E-04   2.207589E-04   2.312254E-06
  b3g  11  -1.513286E-03   4.737945E-05   9.291644E-05  -3.807090E-06   3.575537E-04   2.748583E-04  -9.024302E-05   7.807033E-06

               b3g   9        b3g  10        b3g  11
  b3g   1   1.954843E-03  -2.362282E-03  -1.513286E-03
  b3g   2   4.281136E-04   2.721656E-05   4.737945E-05
  b3g   3   1.124761E-03   7.896669E-05   9.291644E-05
  b3g   4   7.120358E-06  -2.708403E-06  -3.807090E-06
  b3g   5   1.099865E-03  -5.544789E-05   3.575537E-04
  b3g   6  -1.205735E-03   2.576222E-04   2.748583E-04
  b3g   7  -7.589028E-05   2.207589E-04  -9.024302E-05
  b3g   8  -1.947877E-06   2.312254E-06   7.807033E-06
  b3g   9   1.696052E-03  -2.443759E-04  -4.778710E-05
  b3g  10  -2.443759E-04   6.413007E-04   1.922627E-04
  b3g  11  -4.778710E-05   1.922627E-04   9.390455E-04

Natural orbital populations,block 7
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.47154606     0.00597110     0.00507930     0.00304691     0.00114211     0.00060657     0.00036603     0.00023097
              MO     9       MO    10       MO    11
  occ(*)=     0.00006531     0.00005609     0.00001217

          modens reordered block   1

               au    1        au    2        au    3        au    4        au    5        au    6        au    7        au    8
  au    1   0.517501       1.029100E-03   5.678634E-04  -4.185742E-04   5.769065E-04  -1.255513E-04  -3.074222E-03  -6.948470E-04
  au    2   1.029100E-03   8.785175E-04   7.756121E-04  -7.841316E-04  -6.905475E-04  -1.569588E-06   6.320114E-04   1.156080E-04
  au    3   5.678634E-04   7.756121E-04   1.060451E-03  -5.612827E-04  -1.204093E-03   4.190485E-08   2.856581E-04   4.016571E-04
  au    4  -4.185742E-04  -7.841316E-04  -5.612827E-04   1.583221E-03   1.226869E-04   2.086827E-06  -1.063640E-03   9.934360E-04
  au    5   5.769065E-04  -6.905475E-04  -1.204093E-03   1.226869E-04   1.721377E-03  -1.388175E-06   7.786920E-05  -1.133908E-03
  au    6  -1.255513E-04  -1.569588E-06   4.190485E-08   2.086827E-06  -1.388175E-06   8.851724E-05  -2.378952E-06   2.579786E-06
  au    7  -3.074222E-03   6.320114E-04   2.856581E-04  -1.063640E-03   7.786920E-05  -2.378952E-06   1.286919E-03  -7.392436E-04
  au    8  -6.948470E-04   1.156080E-04   4.016571E-04   9.934360E-04  -1.133908E-03   2.579786E-06  -7.392436E-04   2.139804E-03
  au    9  -9.481128E-04  -2.129628E-04  -7.226387E-04   3.959527E-04   9.318980E-04  -1.979770E-06   1.503722E-04   4.370821E-05
  au   10   7.420745E-04   3.502567E-05   6.197184E-05  -1.702006E-04  -1.819159E-04   2.829123E-08  -9.908196E-05   4.476898E-05
  au   11  -2.869323E-04  -2.083946E-06  -4.379142E-07   2.390395E-06  -6.306889E-07   1.675049E-04  -3.322508E-07   2.230870E-06

               au    9        au   10        au   11
  au    1  -9.481128E-04   7.420745E-04  -2.869323E-04
  au    2  -2.129628E-04   3.502567E-05  -2.083946E-06
  au    3  -7.226387E-04   6.197184E-05  -4.379142E-07
  au    4   3.959527E-04  -1.702006E-04   2.390395E-06
  au    5   9.318980E-04  -1.819159E-04  -6.306889E-07
  au    6  -1.979770E-06   2.829123E-08   1.675049E-04
  au    7   1.503722E-04  -9.908196E-05  -3.322508E-07
  au    8   4.370821E-05   4.476898E-05   2.230870E-06
  au    9   1.261733E-03  -2.979804E-04  -9.428804E-07
  au   10  -2.979804E-04   4.733652E-04  -5.544586E-07
  au   11  -9.428804E-07  -5.544586E-07   5.706162E-04

Natural orbital populations,block 8
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     0.51752711     0.00422710     0.00387521     0.00139279     0.00062290     0.00047854     0.00022434     0.00013103
              MO     9       MO    10       MO    11
  occ(*)=     0.00004247     0.00003602     0.00000824


 total number of electrons =   42.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        ag  partial gross atomic populations
   ao class       1ag        2ag        3ag        4ag        5ag        6ag 
    C1_ s       0.024115   1.965539   0.531284   0.472413  -0.012245   0.059382
    C1_ p       0.003468  -0.000704   0.038078   0.129293   0.360705   0.254612
    C1_ d       0.000470   0.000179   0.005236  -0.113717   0.039760  -0.025774
    C2_ s       1.965551   0.024164   1.072838   0.236668  -0.022350   0.033198
    C2_ p       0.000801   0.007425   0.078000   0.109048   0.702361   1.303907
    C2_ d      -0.000546   0.001050   0.010357  -0.053639   0.078994   0.011059
    H1_ s       0.001179   0.000814   0.076352   0.778020   0.284989   0.265182
    H1_ p       0.000351   0.000587   0.005670   0.031274   0.000203  -0.036489
    H2_ s       0.003058   0.000734   0.157680   0.389853   0.545910   0.139153
    H2_ p       0.001552   0.000211   0.011951   0.002199   0.000998  -0.030653
 
   ao class       7ag        8ag        9ag       10ag       11ag       12ag 
    C1_ s       0.004299   0.000525   0.000489   0.000001   0.000887   0.000568
    C1_ p       0.001625   0.000736   0.003881   0.001486   0.000057   0.000425
    C1_ d       0.000370   0.000945  -0.001158   0.000751   0.000792   0.000235
    C2_ s       0.002198   0.000894   0.000250  -0.000162   0.002395  -0.000102
    C2_ p       0.004297   0.001476   0.002438   0.002856   0.000278   0.002260
    C2_ d       0.001498   0.001985  -0.000106   0.001679   0.001876   0.000671
    H1_ s       0.000151   0.002218   0.002835  -0.000027  -0.000265  -0.000257
    H1_ p      -0.000050  -0.000124   0.000173   0.000052   0.000022   0.000068
    H2_ s       0.000101   0.004254   0.001462  -0.000015  -0.000689  -0.000071
    H2_ p       0.000047  -0.000239   0.000078   0.000098   0.000074   0.000860
 
   ao class      13ag       14ag       15ag       16ag       17ag       18ag 
    C1_ s      -0.000407   0.000644  -0.000142   0.000002   0.000074  -0.000083
    C1_ p       0.001330  -0.000043  -0.000125   0.000337   0.000231   0.000011
    C1_ d       0.000727   0.000289   0.000776   0.000080   0.000205   0.000158
    C2_ s      -0.000183   0.001220  -0.000049   0.000019   0.000035  -0.000043
    C2_ p       0.000616  -0.000082   0.000308   0.000788   0.000063   0.000699
    C2_ d       0.001350   0.000595   0.000638   0.000232   0.000156  -0.000106
    H1_ s      -0.001050  -0.000297  -0.000221  -0.000092  -0.000024   0.000000
    H1_ p       0.000286   0.000003   0.000268   0.000031   0.000048   0.000025
    H2_ s      -0.000506  -0.000523  -0.000136  -0.000178  -0.000010   0.000000
    H2_ p       0.000188  -0.000001   0.000183   0.000051   0.000145   0.000105
 
   ao class      19ag       20ag       21ag       22ag       23ag       24ag 
    C1_ s      -0.000097   0.000192   0.000047  -0.000086  -0.000009  -0.000003
    C1_ p      -0.000035   0.000195   0.000090   0.000018   0.000002  -0.000013
    C1_ d      -0.000073   0.000137   0.000064   0.000059   0.000074   0.000027
    C2_ s       0.000201  -0.000094   0.000024  -0.000328  -0.000011  -0.000046
    C2_ p       0.000162  -0.000089  -0.000036   0.000070  -0.000076  -0.000032
    C2_ d       0.000211   0.000005   0.000127   0.000221   0.000169   0.000047
    H1_ s      -0.000054   0.000030   0.000108   0.000099   0.000004   0.000048
    H1_ p       0.000040  -0.000023  -0.000119  -0.000008   0.000012  -0.000007
    H2_ s       0.000037  -0.000036   0.000031   0.000263  -0.000002   0.000123
    H2_ p       0.000078   0.000147  -0.000036  -0.000041  -0.000022  -0.000022
 
   ao class      25ag       26ag       27ag       28ag       29ag       30ag 
    C1_ s      -0.000213   0.000002   0.000046   0.000033  -0.000007   0.000000
    C1_ p       0.000072  -0.000009   0.000024  -0.000024  -0.000005  -0.000013
    C1_ d       0.000044  -0.000006   0.000035   0.000006  -0.000002   0.000004
    C2_ s      -0.000035  -0.000027   0.000014   0.000018   0.000012  -0.000001
    C2_ p       0.000141  -0.000006   0.000008  -0.000024  -0.000000   0.000003
    C2_ d       0.000096   0.000028  -0.000051   0.000069   0.000002  -0.000002
    H1_ s       0.000037   0.000038  -0.000027   0.000013  -0.000002   0.000028
    H1_ p      -0.000008  -0.000003   0.000027  -0.000000   0.000008  -0.000011
    H2_ s       0.000000   0.000081  -0.000008   0.000005  -0.000004   0.000015
    H2_ p      -0.000030  -0.000002   0.000020  -0.000044   0.000021  -0.000003
 
   ao class      31ag       32ag       33ag       34ag       35ag       36ag 
    C1_ s      -0.000005   0.000011  -0.000012  -0.000019   0.000025  -0.000001
    C1_ p       0.000002  -0.000007  -0.000003   0.000003  -0.000006   0.000004
    C1_ d      -0.000000  -0.000000   0.000004  -0.000002  -0.000000  -0.000000
    C2_ s       0.000004  -0.000009   0.000026   0.000031  -0.000009  -0.000000
    C2_ p       0.000020   0.000001  -0.000010  -0.000001  -0.000008   0.000001
    C2_ d      -0.000001   0.000000  -0.000011  -0.000000   0.000000  -0.000000
    H1_ s      -0.000000  -0.000002  -0.000005  -0.000004   0.000002  -0.000002
    H1_ p      -0.000001   0.000014   0.000002  -0.000001   0.000000   0.000000
    H2_ s      -0.000001  -0.000004   0.000004  -0.000010  -0.000001  -0.000001
    H2_ p      -0.000001   0.000004   0.000012   0.000010  -0.000000   0.000000
 
   ao class      37ag       38ag       39ag 
    C1_ s      -0.000001   0.000000   0.000000
    C1_ p       0.000000  -0.000000  -0.000000
    C1_ d      -0.000000  -0.000000   0.000000
    C2_ s       0.000000  -0.000000   0.000000
    C2_ p       0.000001   0.000000  -0.000000
    C2_ d       0.000000  -0.000000   0.000000
    H1_ s      -0.000000   0.000000   0.000000
    H1_ p      -0.000000   0.000000   0.000000
    H2_ s       0.000000   0.000000   0.000000
    H2_ p      -0.000000  -0.000000   0.000000

                        b3u partial gross atomic populations
   ao class       1b3u       2b3u       3b3u       4b3u       5b3u       6b3u
    C1_ s       0.022924   1.952003   0.723266   0.074066   0.082915  -0.000148
    C1_ p      -0.000455   0.001641  -0.038037   0.223761   0.717938   0.000930
    C1_ d      -0.000102   0.000126   0.042607  -0.187010  -0.016586   0.000983
    C2_ s       1.958708   0.025923   0.371430   0.176859   0.026998   0.000156
    C2_ p       0.013886   0.014358   0.207897   0.416594   0.567401   0.008735
    C2_ d       0.003923   0.004386  -0.031103  -0.368314  -0.040912   0.002358
    H1_ s       0.000617   0.000533   0.444565   0.559998   0.427795   0.000925
    H1_ p       0.000315   0.000354   0.040468  -0.003205  -0.015722  -0.000036
    H2_ s       0.000067   0.000410   0.218804   1.089802   0.231553   0.000364
    H2_ p       0.000116   0.000266   0.005094  -0.003844  -0.004769   0.000099
 
   ao class       7b3u       8b3u       9b3u      10b3u      11b3u      12b3u
    C1_ s       0.003513   0.000831   0.000186  -0.000066  -0.000573   0.000159
    C1_ p       0.000471   0.001947   0.002696   0.001858  -0.000159  -0.000093
    C1_ d      -0.000171   0.000110  -0.001094   0.000934   0.001300   0.000347
    C2_ s       0.006744   0.000383   0.000450  -0.000162  -0.000293   0.000075
    C2_ p       0.000639   0.002835   0.005520   0.004263   0.000816  -0.000577
    C2_ d      -0.000574   0.000230  -0.002207   0.000883   0.002167   0.001654
    H1_ s       0.000521   0.003103   0.000858  -0.002258  -0.000057   0.000048
    H1_ p       0.000195  -0.000066   0.000108   0.000015   0.000115  -0.000013
    H2_ s       0.001269   0.001483   0.001806  -0.001123  -0.000026   0.000023
    H2_ p       0.000403  -0.000006   0.000224   0.000394   0.000034   0.000490
 
   ao class      13b3u      14b3u      15b3u      16b3u      17b3u      18b3u
    C1_ s       0.000272   0.000519  -0.000167   0.000005   0.000063   0.000148
    C1_ p       0.001428   0.000592   0.000004  -0.000130   0.000213  -0.000031
    C1_ d      -0.000527   0.000292   0.000089   0.000361   0.000012   0.000095
    C2_ s       0.000616   0.000237  -0.000082  -0.000004   0.000101   0.000096
    C2_ p       0.002823  -0.000083   0.000645  -0.000252   0.000246  -0.000092
    C2_ d      -0.001074   0.000128   0.000242   0.000717  -0.000034   0.000159
    H1_ s      -0.000853  -0.000442   0.000016   0.000048  -0.000054   0.000038
    H1_ p       0.000336   0.000111   0.000012  -0.000026  -0.000008  -0.000007
    H2_ s      -0.001676  -0.000214   0.000008   0.000080  -0.000115   0.000021
    H2_ p       0.000663   0.000117   0.000036  -0.000047   0.000020   0.000015
 
   ao class      19b3u      20b3u      21b3u      22b3u      23b3u      24b3u
    C1_ s      -0.000199  -0.000080  -0.000211  -0.000039   0.000005  -0.000082
    C1_ p       0.000045   0.000100  -0.000007   0.000001   0.000005   0.000033
    C1_ d       0.000017   0.000023   0.000148   0.000000   0.000062   0.000062
    C2_ s      -0.000102  -0.000063  -0.000211  -0.000036   0.000003  -0.000131
    C2_ p       0.000419   0.000171   0.000509   0.000334  -0.000026   0.000057
    C2_ d       0.000051   0.000039  -0.000178  -0.000242   0.000044   0.000115
    H1_ s      -0.000031   0.000058   0.000083   0.000041  -0.000012   0.000023
    H1_ p       0.000014  -0.000048  -0.000011   0.000003   0.000017  -0.000019
    H2_ s      -0.000015   0.000110   0.000045   0.000022  -0.000007   0.000046
    H2_ p       0.000125  -0.000092   0.000013   0.000051   0.000007  -0.000037
 
   ao class      25b3u      26b3u      27b3u      28b3u      29b3u      30b3u
    C1_ s       0.000014   0.000007   0.000014   0.000006  -0.000007   0.000005
    C1_ p      -0.000006  -0.000011  -0.000033  -0.000001   0.000005  -0.000011
    C1_ d       0.000008  -0.000009   0.000008  -0.000000  -0.000001   0.000002
    C2_ s       0.000028   0.000005   0.000011   0.000012  -0.000004  -0.000003
    C2_ p      -0.000011  -0.000018  -0.000105  -0.000001   0.000014  -0.000006
    C2_ d       0.000014   0.000061   0.000106   0.000002  -0.000007   0.000006
    H1_ s      -0.000017   0.000012   0.000033  -0.000002  -0.000004   0.000008
    H1_ p       0.000016  -0.000001   0.000004   0.000003   0.000011  -0.000001
    H2_ s      -0.000035   0.000006   0.000016  -0.000004  -0.000002   0.000007
    H2_ p       0.000032  -0.000021  -0.000024   0.000005   0.000005  -0.000000
 
   ao class      31b3u      32b3u      33b3u      34b3u      35b3u      36b3u
    C1_ s      -0.000009   0.000037  -0.000016  -0.000000   0.000002  -0.000000
    C1_ p       0.000008  -0.000003   0.000006   0.000002  -0.000002   0.000001
    C1_ d       0.000001  -0.000000  -0.000001  -0.000000  -0.000000  -0.000000
    C2_ s      -0.000001  -0.000000   0.000019  -0.000000   0.000003  -0.000001
    C2_ p       0.000011  -0.000031  -0.000011   0.000002  -0.000003   0.000001
    C2_ d      -0.000006  -0.000002  -0.000000  -0.000001  -0.000000  -0.000000
    H1_ s      -0.000006   0.000002  -0.000001  -0.000001   0.000000  -0.000000
    H1_ p       0.000000   0.000000   0.000002   0.000000   0.000000  -0.000000
    H2_ s       0.000002   0.000001   0.000002  -0.000001   0.000000  -0.000000
    H2_ p       0.000007   0.000000   0.000005   0.000001   0.000000  -0.000000
 
   ao class      37b3u      38b3u      39b3u
    C1_ s      -0.000000  -0.000000   0.000000
    C1_ p       0.000000  -0.000000  -0.000000
    C1_ d      -0.000000  -0.000000   0.000000
    C2_ s      -0.000000  -0.000000   0.000000
    C2_ p      -0.000000   0.000000  -0.000000
    C2_ d       0.000000  -0.000000   0.000000
    H1_ s       0.000000   0.000000   0.000000
    H1_ p       0.000000   0.000000  -0.000000
    H2_ s       0.000000  -0.000000   0.000000
    H2_ p      -0.000000   0.000000  -0.000000

                        b2u partial gross atomic populations
   ao class       1b2u       2b2u       3b2u       4b2u       5b2u       6b2u
    C1_ p       0.016768   0.151870   0.144970   0.627720   0.005251   0.001315
    C1_ d       0.005055  -0.036042  -0.020761   0.017471   0.001190   0.000146
    C2_ s       1.965412   1.103487   0.114203   0.000626  -0.000103   0.001201
    C2_ p       0.007796   0.007812   1.143849   1.335535   0.004206   0.003445
    C2_ d       0.001844   0.046168  -0.038413   0.036741   0.002227   0.000236
    H1_ p       0.000103  -0.010382   0.002365  -0.014305   0.000077   0.000015
    H2_ s       0.001851   0.666249   0.653235   0.000506   0.001335   0.004536
    H2_ p       0.001171   0.055852  -0.022907  -0.029846  -0.000020  -0.000086
 
   ao class       7b2u       8b2u       9b2u      10b2u      11b2u      12b2u
    C1_ p       0.000742   0.002487   0.000531  -0.000347  -0.000247   0.000374
    C1_ d       0.000917   0.000297   0.000983   0.000964  -0.000010  -0.000041
    C2_ s      -0.000001  -0.000592  -0.000748   0.000238   0.000741   0.000000
    C2_ p       0.001716   0.003967   0.000083  -0.000306   0.000763   0.000668
    C2_ d       0.001762   0.001598   0.002455   0.001017   0.000436  -0.000047
    H1_ p       0.000459   0.000268  -0.000019   0.000334   0.000041   0.000057
    H2_ s      -0.000000  -0.003336  -0.000110   0.000071  -0.000644  -0.000000
    H2_ p       0.000944   0.000143   0.000170   0.000147   0.000181   0.000106
 
   ao class      13b2u      14b2u      15b2u      16b2u      17b2u      18b2u
    C1_ p       0.000408   0.000026  -0.000107   0.000277  -0.000039   0.000332
    C1_ d       0.000118   0.000142   0.000060   0.000033   0.000173  -0.000204
    C2_ s      -0.000235  -0.000001   0.000236  -0.000325  -0.000001  -0.000427
    C2_ p       0.000218   0.000031  -0.000028   0.000196  -0.000020   0.000192
    C2_ d       0.000214   0.000256   0.000213   0.000050   0.000276   0.000159
    H1_ p       0.000019   0.000079   0.000021   0.000078  -0.000060   0.000018
    H2_ s       0.000025  -0.000000   0.000057  -0.000045   0.000000   0.000128
    H2_ p       0.000030   0.000153  -0.000010   0.000063  -0.000113  -0.000015
 
   ao class      19b2u      20b2u      21b2u      22b2u      23b2u      24b2u
    C1_ p       0.000226  -0.000025  -0.000010  -0.000055  -0.000001   0.000008
    C1_ d      -0.000163   0.000007   0.000049   0.000063  -0.000006  -0.000003
    C2_ s      -0.000080   0.000017   0.000013   0.000026  -0.000000  -0.000010
    C2_ p       0.000104  -0.000005  -0.000023  -0.000079  -0.000003   0.000010
    C2_ d      -0.000075   0.000100   0.000010   0.000044  -0.000014  -0.000004
    H1_ p       0.000034  -0.000001  -0.000015  -0.000016   0.000015  -0.000001
    H2_ s       0.000076  -0.000020   0.000016   0.000049   0.000000  -0.000004
    H2_ p       0.000014   0.000027  -0.000008  -0.000002   0.000030   0.000016
 
   ao class      25b2u      26b2u      27b2u      28b2u      29b2u      30b2u
    C1_ p       0.000010  -0.000023  -0.000003   0.000001   0.000000   0.000000
    C1_ d      -0.000005  -0.000000   0.000000  -0.000001   0.000000  -0.000000
    C2_ s      -0.000013   0.000042   0.000003  -0.000001  -0.000000  -0.000000
    C2_ p       0.000017  -0.000016   0.000004   0.000003   0.000000  -0.000000
    C2_ d      -0.000002  -0.000002  -0.000000  -0.000000  -0.000000  -0.000000
    H1_ p       0.000006   0.000000  -0.000000   0.000001  -0.000000   0.000000
    H2_ s      -0.000010   0.000003   0.000000  -0.000002   0.000000   0.000000
    H2_ p       0.000002   0.000000  -0.000000   0.000000   0.000000   0.000000

                        b1g partial gross atomic populations
   ao class       1b1g       2b1g       3b1g       4b1g       5b1g       6b1g
    C1_ p       0.003407   0.034162   0.760366   0.002184   0.000350   0.002714
    C1_ d       0.000348   0.002090   0.015453   0.000964   0.000325   0.000634
    C2_ s       1.989749   0.722242   0.085691   0.006558   0.000750   0.000002
    C2_ p       0.004512   0.206612   0.789669   0.003881   0.005958   0.005244
    C2_ d       0.001538  -0.166458  -0.032009   0.000810  -0.001607   0.000975
    H1_ p       0.000044  -0.009472  -0.007544   0.000048  -0.000005   0.000055
    H2_ s       0.000202   1.149521   0.421191   0.000226   0.004306   0.000000
    H2_ p       0.000200   0.042742  -0.059158  -0.000020   0.000255   0.000117
 
   ao class       7b1g       8b1g       9b1g      10b1g      11b1g      12b1g
    C1_ p       0.001574  -0.000011  -0.000257   0.000251  -0.000033   0.000450
    C1_ d       0.000377   0.000682   0.000689   0.000193   0.000043  -0.000124
    C2_ s       0.000485  -0.000574   0.000000  -0.000197   0.000107  -0.000134
    C2_ p       0.001040   0.001936  -0.000518  -0.000068   0.000339   0.000268
    C2_ d       0.000547   0.001383   0.001367   0.001220   0.000311   0.000174
    H1_ p       0.000560   0.000034   0.000159   0.000025   0.000082   0.000064
    H2_ s      -0.000295  -0.001550   0.000000  -0.000352  -0.000033   0.000000
    H2_ p       0.000362   0.000445   0.000321   0.000422   0.000109   0.000073
 
   ao class      13b1g      14b1g      15b1g      16b1g      17b1g      18b1g
    C1_ p      -0.000008  -0.000041   0.000155  -0.000046   0.000065   0.000002
    C1_ d       0.000051   0.000062  -0.000148   0.000083   0.000068  -0.000043
    C2_ s       0.000118   0.000078  -0.000000  -0.000020  -0.000250   0.000064
    C2_ p       0.000219   0.000090   0.000314  -0.000024   0.000141   0.000041
    C2_ d       0.000040   0.000135  -0.000298   0.000155   0.000087   0.000013
    H1_ p       0.000096   0.000009   0.000065  -0.000018  -0.000023   0.000010
    H2_ s      -0.000048   0.000137   0.000000   0.000003   0.000040  -0.000043
    H2_ p       0.000001  -0.000168   0.000130   0.000009  -0.000021   0.000046
 
   ao class      19b1g      20b1g      21b1g      22b1g      23b1g      24b1g
    C1_ p       0.000179  -0.000013   0.000012  -0.000067  -0.000000  -0.000005
    C1_ d      -0.000183   0.000050  -0.000002   0.000093  -0.000000  -0.000008
    C2_ s       0.000000   0.000048  -0.000001  -0.000000   0.000000   0.000006
    C2_ p       0.000346  -0.000033  -0.000022  -0.000126  -0.000007  -0.000005
    C2_ d      -0.000350   0.000027   0.000005   0.000180   0.000000  -0.000001
    H1_ p       0.000025  -0.000031   0.000003  -0.000021   0.000002   0.000009
    H2_ s       0.000000   0.000019   0.000043   0.000000  -0.000004   0.000003
    H2_ p       0.000048  -0.000016  -0.000017  -0.000040   0.000018   0.000007
 
   ao class      25b1g      26b1g      27b1g      28b1g      29b1g      30b1g
    C1_ p      -0.000002  -0.000001   0.000001   0.000000   0.000000   0.000000
    C1_ d       0.000000  -0.000000  -0.000001   0.000000   0.000000   0.000000
    C2_ s       0.000017  -0.000002   0.000000  -0.000001   0.000000   0.000000
    C2_ p      -0.000014   0.000007   0.000001   0.000001  -0.000000   0.000000
    C2_ d      -0.000000  -0.000000  -0.000003  -0.000000   0.000000   0.000000
    H1_ p      -0.000000   0.000000   0.000001  -0.000000   0.000000  -0.000000
    H2_ s       0.000003  -0.000003  -0.000000  -0.000000   0.000000   0.000000
    H2_ p       0.000000   0.000000   0.000002  -0.000000   0.000000  -0.000000

                        b1u partial gross atomic populations
   ao class       1b1u       2b1u       3b1u       4b1u       5b1u       6b1u
    C1_ p       0.521048   0.352896  -0.000798   0.000010   0.001148   0.001911
    C1_ d       0.010426  -0.001386   0.000105   0.001647   0.000916   0.000003
    C2_ p       1.356118   0.145875   0.006270   0.000055  -0.000201   0.000018
    C2_ d       0.021414   0.031664   0.000185   0.003399   0.001601   0.001827
    H1_ p       0.003572   0.009933   0.000023   0.000120   0.000518   0.000045
    H2_ p       0.011762   0.005370   0.000346   0.000229   0.000134   0.000002
 
   ao class       7b1u       8b1u       9b1u      10b1u      11b1u      12b1u
    C1_ p       0.000172  -0.000004  -0.000023   0.000003   0.000005   0.000093
    C1_ d       0.000373  -0.000034  -0.000013   0.000086   0.000103   0.000105
    C2_ p       0.000154   0.000492   0.000314  -0.000000   0.000006  -0.000001
    C2_ d       0.000516   0.000051  -0.000002   0.000121   0.000345   0.000036
    H1_ p       0.000117   0.000087   0.000134   0.000180  -0.000031  -0.000018
    H2_ p       0.000088   0.000176   0.000223   0.000091  -0.000076  -0.000004
 
   ao class      13b1u      14b1u      15b1u      16b1u
    C1_ p       0.000024   0.000002   0.000016   0.000002
    C1_ d       0.000033  -0.000010   0.000019  -0.000002
    C2_ p      -0.000009   0.000001  -0.000005  -0.000000
    C2_ d       0.000079   0.000007  -0.000002  -0.000001
    H1_ p      -0.000002   0.000023  -0.000003   0.000007
    H2_ p      -0.000001   0.000026   0.000013   0.000003

                        b2g partial gross atomic populations
   ao class       1b2g       2b2g       3b2g       4b2g       5b2g       6b2g
    C1_ p       0.973528   0.007269   0.001111   0.001419   0.000189   0.000024
    C1_ d       0.000945  -0.000114   0.001647   0.000663   0.000282   0.000042
    C2_ p       0.450796   0.015739   0.000485   0.000587   0.000406   0.000012
    C2_ d       0.031690  -0.000209   0.001084   0.000380   0.000536   0.001052
    H1_ p       0.013781   0.000177   0.000486   0.000077   0.000298   0.000005
    H2_ p       0.006177   0.000325   0.000237   0.000035   0.000587   0.000003
 
   ao class       7b2g       8b2g       9b2g      10b2g      11b2g      12b2g
    C1_ p      -0.000003   0.000202   0.000031   0.000031   0.000027  -0.000032
    C1_ d       0.000049  -0.000102   0.000002   0.000150   0.000016   0.000055
    C2_ p      -0.000000   0.000432   0.000084   0.000017   0.000017  -0.000028
    C2_ d       0.000068  -0.000200   0.000026   0.000075   0.000008   0.000117
    H1_ p       0.000326   0.000030   0.000057  -0.000032   0.000001  -0.000016
    H2_ p       0.000164   0.000062   0.000106  -0.000015   0.000002  -0.000034
 
   ao class      13b2g      14b2g      15b2g      16b2g
    C1_ p       0.000013   0.000001  -0.000000   0.000001
    C1_ d       0.000005  -0.000002  -0.000001   0.000000
    C2_ p      -0.000016   0.000001  -0.000000   0.000001
    C2_ d       0.000057  -0.000001  -0.000002   0.000000
    H1_ p       0.000000   0.000008   0.000004  -0.000000
    H2_ p      -0.000002   0.000005   0.000006  -0.000000

                        b3g partial gross atomic populations
   ao class       1b3g       2b3g       3b3g       4b3g       5b3g       6b3g
    C1_ d       0.018587   0.002000   0.000123   0.000015   0.000691   0.000026
    C2_ p       1.421810   0.000000   0.001488   0.002031   0.000035  -0.000003
    C2_ d       0.011717   0.003971   0.002758   0.000898   0.000409   0.000092
    H2_ p       0.019433   0.000000   0.000711   0.000102   0.000007   0.000491
 
   ao class       7b3g       8b3g       9b3g      10b3g      11b3g
    C1_ d       0.000120   0.000000   0.000022   0.000018  -0.000000
    C2_ p       0.000000   0.000044   0.000021   0.000018   0.000003
    C2_ d       0.000246   0.000235   0.000021   0.000013  -0.000003
    H2_ p      -0.000000  -0.000049   0.000001   0.000008   0.000013

                        au  partial gross atomic populations
   ao class       1au        2au        3au        4au        5au        6au 
    C1_ d       0.023716   0.000511   0.001499   0.000248   0.000209   0.000050
    C2_ p       0.471095   0.001431   0.001479   0.000283   0.000000   0.000005
    C2_ d       0.007585   0.001592   0.000865   0.000652   0.000414   0.000155
    H2_ p       0.015132   0.000692   0.000033   0.000210   0.000000   0.000269
 
   ao class       7au        8au        9au       10au       11au 
    C1_ d       0.000012   0.000055   0.000001   0.000012  -0.000000
    C2_ p       0.000134   0.000013   0.000011   0.000000   0.000002
    C2_ d       0.000101   0.000063   0.000032   0.000024  -0.000004
    H2_ p      -0.000023   0.000001  -0.000003  -0.000000   0.000010


                        gross atomic populations
     ao           C1_        C2_        H1_        H2_
      s         5.906610  11.872518   2.845405   5.680632
      p         5.328563  10.872082   0.018527   0.042305
      d        -0.191171  -0.375472   0.000000   0.000000
    total      11.044002  22.369128   2.863932   5.722937
 

 Total number of electrons:   42.00000000

