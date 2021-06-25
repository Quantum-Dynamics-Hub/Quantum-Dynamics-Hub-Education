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

                CSFs      3909      465486    11047539    12628699    24145633
      internal walks     15280       34818       13545       15537       79180
valid internal walks      3909       32226        6039        6699       48873
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
 compressed index vector length=                 29763
 echo of the input for program ciudg:
 ------------------------------------------------------------------------
  &input
  NTYPE = 0,
  GSET = 0,
   DAVCOR =10,
  NCOREL = 30
  NROOT = 1
  IVMODE = 3
  NBKITR = 1
  NVBKMN = 1
  RTOLBK = 1e-3,
  NITER = 20
  NVCIMN = 3
  RTOLCI = 1e-3,
  NVCIMX = 6
  NVRFMX = 6
  NVBKMX = 6
  IDEN  = 1
  CSFPRN = 10,
 /&end
 ------------------------------------------------------------------------
lodens (list->root)=  1
invlodens (root->list)=  1
 bummer (warning):resetting fileloc for seriel operation0
 USING SEGMENTS OF EQUAL SIZE

****************  list of control variables  ****************
 lvlprt =    0      nroot  =    1      noldv  =   0      noldhv =   0
 nunitv =    3      nbkitr =    1      niter  =  20      davcor =  10
 csfprn =   10      ivmode =    3      istrt  =   0      vout   =   0
 iortls =    0      nvbkmx =    6      ibktv  =  -1      ibkthv =  -1
 nvcimx =    6      icitv  =   -1      icithv =  -1      frcsub =   0
 nvbkmn =    1      nvcimn =    3      maxseg = 300      nrfitr =  30
 ncorel =   30      nvrfmx =    6      nvrfmn =   3      iden   =   1
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
 Computing density:                    .drt1.state1
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
 Hermit Integral Program : SIFS version  compute-0-12      14:50:12.965 24-Jun-21
  cidrt_title                                                                    
 MO-coefficients from mcscf.x                                                    
  with dummy occupation 1.0 for active orbitals                                  
  total ao core energy =  201.595839628                                          
 MCSCF energy =    -230.596267333                                                
 SIFS file created by program tran.      compute-0-12      14:50:28.812 24-Jun-21

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
 !timer: first half-sort required        cpu_time=     0.965 walltime=     3.128

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
 !timer: second half-sort required       cpu_time=     1.052 walltime=     5.578
 !timer: cisrt complete                  cpu_time=     2.021 walltime=     8.716
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
 nrow  =   148 nsym  =     8 ssym  =     1 lenbuf=  1600
 nwalk,xbar:      79180    15280    34818    13545    15537
 nvalwt,nvalw:    48873     3909    32226     6039     6699
 ncsft:        24145633
 total number of valid internal walks:   48873
 nvalz,nvaly,nvalx,nvalw =     3909   32226    6039    6699

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
 calcthrxt: niot,maxw1=                    18                 16712
 block size     0
 pthz,pthy,pthx,pthw: 15280 34818 13545 15537 total internal walks:   79180
 maxlp3,n2lp,n1lp,n0lp 16712     0     0     0
 orbsym(*)= 1 1 1 1 2 2 2 3 3 3 4 4 5 5 5 6 7 8
setref: retained number of references =    55
 setref: total/valid number of walks=                 15280
                  3909
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
    threx            160285
    twoex             36190
    onex               8593
    allin              6144
    diagon            32097
               =======
   maximum           160285
 
  __ static summary __ 
   reflst              3909
   hrfspc              3909
               -------
   static->            3909
 
  __ core required  __ 
   totstc              3909
   max n-ex          160285
               -------
   totnec->          164194
 
  __ core available __ 
   totspc        2621439999
   totnec -          164194
               -------
   totvec->      2621275805

 number of external paths / symmetry
 vertex x    2038    2119    2056    2055    1469    1467    1413    1411
 vertex w    2206    2119    2056    2055    1469    1467    1413    1411
segment: free space=  2621275805
 reducing frespc by                214114 to             2621061691 
  for index/conft/indsym storage .
 resegmenting ...



                   segmentation summary for type all-internal
 -------------------------------------------------------------------------------
 seg.      no. of|    no. of|  starting|  internal|  starting|  starting|
  no.    internal|        ci|       csf|     walks|      walk|       DRT|
            paths|  elements|    number|     /seg.|    number|    record|
 -------------------------------------------------------------------------------
  Z 1       15280|      3909|         0|      3909|         0|         1|
 -------------------------------------------------------------------------------
  Y 2       34818|    465486|      3909|     32226|      3909|         2|
 -------------------------------------------------------------------------------
  X 3       12681|  11047539|    469395|      6039|     36135|         5|
 -------------------------------------------------------------------------------
  W 4       15231|  12628699|  11516934|      6699|     42174|         6|
 -------------------------------------------------------------------------------
max. additional memory requirements:index=       16828DP  conft+indsym=      128904DP  drtbuffer=       68382 DP

dimension of the ci-matrix ->>>  24145633

 executing brd_struct for civct
 gentasklist: ntask=                    20
                    TASKLIST
----------------------------------------------------------------------------------------------------
TASK# BRA# KET#  T-TYPE    DESCR.   SEGMENTTYPE    SEGEL              SEGCI          VWALKS   
----------------------------------------------------------------------------------------------------
     1  3   1    24      two-ext xz   2X  3 1   12681   15280   11047539       3909    6039    3909
     2  4   1    25      two-ext wz   2X  4 1   15231   15280   12628699       3909    6699    3909
     3  4   3    26      two-ext wx*  WX  4 3   15231   12681   12628699   11047539    6699    6039
     4  4   3    27      two-ext wx+  WX  4 3   15231   12681   12628699   11047539    6699    6039
     5  2   1    11      one-ext yz   1X  2 1   34818   15280     465486       3909   32226    3909
     6  3   2    15      1ex3ex yx    3X  3 2   12681   34818   11047539     465486    6039   32226
     7  4   2    16      1ex3ex yw    3X  4 2   15231   34818   12628699     465486    6699   32226
     8  1   1     1      allint zz    OX  1 1   15280   15280       3909       3909    3909    3909
     9  2   2     5      0ex2ex yy    OX  2 2   34818   34818     465486     465486   32226   32226
    10  3   3     6      0ex2ex xx*   OX  3 3   12681   12681   11047539   11047539    6039    6039
    11  3   3    18      0ex2ex xx+   OX  3 3   12681   12681   11047539   11047539    6039    6039
    12  4   4     7      0ex2ex ww*   OX  4 4   15231   15231   12628699   12628699    6699    6699
    13  4   4    19      0ex2ex ww+   OX  4 4   15231   15231   12628699   12628699    6699    6699
    14  2   2    42      four-ext y   4X  2 2   34818   34818     465486     465486   32226   32226
    15  3   3    43      four-ext x   4X  3 3   12681   12681   11047539   11047539    6039    6039
    16  4   4    44      four-ext w   4X  4 4   15231   15231   12628699   12628699    6699    6699
    17  1   1    75      dg-024ext z  OX  1 1   15280   15280       3909       3909    3909    3909
    18  2   2    76      dg-024ext y  OX  2 2   34818   34818     465486     465486   32226   32226
    19  3   3    77      dg-024ext x  OX  3 3   12681   12681   11047539   11047539    6039    6039
    20  4   4    78      dg-024ext w  OX  4 4   15231   15231   12628699   12628699    6699    6699
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
 initializing v-file: 1:              24145633

    ---------trial vector generation----------

    trial vectors will be created by: 

    (ivmode= 3) diagonalizing h in the reference space.                     

      3 vectors will be written to unit 11 beginning with logical record   1

            3 vectors will be created
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:      113824 2x:           0 4x:           0
All internal counts: zz :      581198 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:           0    task #     2:           0    task #     3:           0    task #     4:           0
task #     5:           0    task #     6:           0    task #     7:           0    task #     8:      492288
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:       80535    task #    18:           0    task #    19:           0    task #    20:           0
 reference space has dimension      55
 dsyevx: computed roots 1 to    6(converged:   6)

    root           eigenvalues
    ----           ------------
       1        -230.7649519907
       2        -230.4274558047
       3        -230.3774972608
       4        -230.3216094616
       5        -230.3002238035
       6        -230.2734751886

 strefv generated    3 initial ci vector(s).
    ---------end of vector generation---------

 ufvoutnew: ... writing  recamt=                  3909

         vector  1 from unit 11 written to unit 49 filename cirefv              

 ************************************************************************
 beginning the bk-type iterative procedure (nzcsf=  3909)...
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:          24145633
 number of initial trial vectors:                         3
 number of initial matrix-vector products:                0
 maximum dimension of the subspace vectors:               6
 number of roots to converge:                             1
 number of iterations:                                    1
 residual norm convergence criteria:               0.001000

          starting bk iteration   1

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:           0 xx:           0 ww:           0
One-external counts: yz :     2332820 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       15843 wz:       18722 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:           0    task #     4:           0
task #     5:     1946327    task #     6:           0    task #     7:           0    task #     8:      492288
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:           0 xx:           0 ww:           0
One-external counts: yz :     2332820 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       15843 wz:       18722 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:           0    task #     4:           0
task #     5:     1946327    task #     6:           0    task #     7:           0    task #     8:      492288
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:           0 xx:           0 ww:           0
One-external counts: yz :     2332820 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       15843 wz:       18722 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:           0    task #     4:           0
task #     5:     1946327    task #     6:           0    task #     7:           0    task #     8:      492288
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3
   ht   1    -6.31429113
   ht   2     0.00000000    -5.97679494
   ht   3     0.00000000    -0.00000000    -5.92683640

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3
 ref    1   -1.00000       8.612902E-15   8.673855E-15

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3
 ref    1    1.00000       7.418208E-29   7.523576E-29

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3
 ref:   1    -1.00000000     0.00000000     0.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.7649519907 -3.5527E-14  8.7920E-01  1.5702E+00  1.0000E-03   
 mr-sdci #  1  2   -230.4274558047  6.2172E-15  0.0000E+00  1.5760E+00  1.0000E-04   
 mr-sdci #  1  3   -230.3774972608  9.9476E-14  0.0000E+00  1.5424E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.180675
time for cinew                         2.284701
time for eigenvalue solver             0.000137
time for vector access                 0.000001

 mr-sdci  convergence not reached after  1 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.7649519907 -3.5527E-14  8.7920E-01  1.5702E+00  1.0000E-03   
 mr-sdci #  1  2   -230.4274558047  6.2172E-15  0.0000E+00  1.5760E+00  1.0000E-04   
 mr-sdci #  1  3   -230.3774972608  9.9476E-14  0.0000E+00  1.5424E+00  1.0000E-04   
 
    1 of the   4 expansion vectors are transformed.
    1 of the   3 matrix-vector products are transformed.

    1 expansion eigenvectors written to unit nvfile (= 11)
    1 matrix-vector products written to unit nhvfil (= 10)

 ************************************************************************
 beginning the ci iterative diagonalization procedure... 
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:          24145633
 number of initial trial vectors:                         1
 number of initial matrix-vector products:                1
 maximum dimension of the subspace vectors:               6
 number of roots to converge:                             1
 number of iterations:                                   20
 residual norm convergence criteria:               0.001000

          starting ci iteration   1

 Final subspace hamiltonian 

                ht   1
   ht   1    -6.31429113

          calcsovref: eigensolution overlap with references block   1

              v      1
 ref    1   -1.00000    

          calcsovref: reference weight per eigenvector block   1

              v      1
 ref    1    1.00000    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1

          calcsovref: aa weight per eigenvector block   1

              v      1

          reference overlap matrix  block   1

                ci   1
 ref:   1    -1.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.7649519907  1.7764E-15  8.7920E-01  1.5702E+00  1.0000E-03   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.022423
time for cinew                         2.002576
time for eigenvalue solver             0.000000
time for vector access                 0.000001

          starting ci iteration   2

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:     5571191 xx:      651346 ww:      725522
One-external counts: yz :     2332820 yx:     3683468 yw:     3794578
Two-external counts: yy :     2169709 ww:      357828 xx:      351534 xz:       15843 wz:       18722 wx:      577746
Three-ext.   counts: yx :      757585 yw:      797389

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:      240809    task #     4:      120404
task #     5:     1946327    task #     6:     2786664    task #     7:     2855666    task #     8:      492288
task #     9:     4799038    task #    10:      158133    task #    11:      158133    task #    12:      147051
task #    13:      147051    task #    14:       73525    task #    15:       73527    task #    16:       73527
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 Final subspace hamiltonian 

                ht   1         ht   2
   ht   1    -6.31429113
   ht   2    -0.87919702    -1.44829669

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2
 ref    1   0.911300      -0.411742    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2
 ref    1   0.830468       0.169532    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2

          reference overlap matrix  block   1

                ci   1         ci   2
 ref:   1     0.91130028    -0.41174240

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  2  1   -231.4087655143  6.4381E-01  3.1532E-02  3.1307E-01  1.0000E-03   
 mr-sdci #  2  2   -227.6111681957 -2.8163E+00  0.0000E+00  1.2770E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.065285
time for cinew                         1.612221
time for eigenvalue solver             0.000069
time for vector access                 0.000000

          starting ci iteration   3

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:     5571191 xx:      651346 ww:      725522
One-external counts: yz :     2332820 yx:     3683468 yw:     3794578
Two-external counts: yy :     2169709 ww:      357828 xx:      351534 xz:       15843 wz:       18722 wx:      577746
Three-ext.   counts: yx :      757585 yw:      797389

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:      240809    task #     4:      120404
task #     5:     1946327    task #     6:     2786664    task #     7:     2855666    task #     8:      492288
task #     9:     4799038    task #    10:      158133    task #    11:      158133    task #    12:      147051
task #    13:      147051    task #    14:       73525    task #    15:       73527    task #    16:       73527
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3
   ht   1    -6.31429113
   ht   2    -0.87919702    -1.44829669
   ht   3     0.05338729    -0.05653528    -0.04781701

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3
 ref    1  -0.910351       8.028195E-02  -0.405976    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3
 ref    1   0.828738       6.445192E-03   0.164817    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3
 ref:   1    -0.91035053     0.08028195    -0.40597627

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  3  1   -231.4346672832  2.5902E-02  2.3707E-03  8.3095E-02  1.0000E-03   
 mr-sdci #  3  2   -228.2886883879  6.7752E-01  0.0000E+00  1.2817E+00  1.0000E-04   
 mr-sdci #  3  3   -227.4942975008 -2.8832E+00  0.0000E+00  1.3313E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.106598
time for cinew                         1.773010
time for eigenvalue solver             0.000076
time for vector access                 0.000000

          starting ci iteration   4

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:     5571191 xx:      651346 ww:      725522
One-external counts: yz :     2332820 yx:     3683468 yw:     3794578
Two-external counts: yy :     2169709 ww:      357828 xx:      351534 xz:       15843 wz:       18722 wx:      577746
Three-ext.   counts: yx :      757585 yw:      797389

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:      240809    task #     4:      120404
task #     5:     1946327    task #     6:     2786664    task #     7:     2855666    task #     8:      492288
task #     9:     4799038    task #    10:      158133    task #    11:      158133    task #    12:      147051
task #    13:      147051    task #    14:       73525    task #    15:       73527    task #    16:       73527
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1    -6.31429113
   ht   2    -0.87919702    -1.44829669
   ht   3     0.05338729    -0.05653528    -0.04781701
   ht   4     0.02266828     0.00738004    -0.00160225    -0.00388195

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4
 ref    1  -0.910440       6.736120E-02   1.250301E-02  -0.407928    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4
 ref    1   0.828901       4.537531E-03   1.563252E-04   0.166405    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1    -0.91043997     0.06736120     0.01250301    -0.40792793

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  4  1   -231.4363578597  1.6906E-03  2.0691E-04  2.5256E-02  1.0000E-03   
 mr-sdci #  4  2   -228.4972530892  2.0856E-01  0.0000E+00  1.1462E+00  1.0000E-04   
 mr-sdci #  4  3   -227.9397430303  4.4545E-01  0.0000E+00  1.8330E+00  1.0000E-04   
 mr-sdci #  4  4   -227.4850210483  3.0344E+00  0.0000E+00  1.3132E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.138763
time for cinew                         1.928741
time for eigenvalue solver             0.000000
time for vector access                 0.000015

          starting ci iteration   5

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:     5571191 xx:      651346 ww:      725522
One-external counts: yz :     2332820 yx:     3683468 yw:     3794578
Two-external counts: yy :     2169709 ww:      357828 xx:      351534 xz:       15843 wz:       18722 wx:      577746
Three-ext.   counts: yx :      757585 yw:      797389

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:      240809    task #     4:      120404
task #     5:     1946327    task #     6:     2786664    task #     7:     2855666    task #     8:      492288
task #     9:     4799038    task #    10:      158133    task #    11:      158133    task #    12:      147051
task #    13:      147051    task #    14:       73525    task #    15:       73527    task #    16:       73527
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.31429113
   ht   2    -0.87919702    -1.44829669
   ht   3     0.05338729    -0.05653528    -0.04781701
   ht   4     0.02266828     0.00738004    -0.00160225    -0.00388195
   ht   5    -0.00443541    -0.00240718     0.00047654     0.00014615    -0.00031379

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.910245       1.102398E-02  -0.169987      -0.103794      -0.362854    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.828546       1.215281E-04   2.889575E-02   1.077325E-02   0.131663    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.91024517     0.01102398    -0.16998750    -0.10379426    -0.36285426

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  5  1   -231.4365270740  1.6921E-04  2.4254E-05  8.1247E-03  1.0000E-03   
 mr-sdci #  5  2   -228.7871272791  2.8987E-01  0.0000E+00  1.1791E+00  1.0000E-04   
 mr-sdci #  5  3   -228.1073661607  1.6762E-01  0.0000E+00  1.2162E+00  1.0000E-04   
 mr-sdci #  5  4   -227.8948000945  4.0978E-01  0.0000E+00  2.0140E+00  1.0000E-04   
 mr-sdci #  5  5   -227.3587091574  2.9080E+00  0.0000E+00  1.3411E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.175659
time for cinew                         2.032715
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   6

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:     5571191 xx:      651346 ww:      725522
One-external counts: yz :     2332820 yx:     3683468 yw:     3794578
Two-external counts: yy :     2169709 ww:      357828 xx:      351534 xz:       15843 wz:       18722 wx:      577746
Three-ext.   counts: yx :      757585 yw:      797389

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:      240809    task #     4:      120404
task #     5:     1946327    task #     6:     2786664    task #     7:     2855666    task #     8:      492288
task #     9:     4799038    task #    10:      158133    task #    11:      158133    task #    12:      147051
task #    13:      147051    task #    14:       73525    task #    15:       73527    task #    16:       73527
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.31429113
   ht   2    -0.87919702    -1.44829669
   ht   3     0.05338729    -0.05653528    -0.04781701
   ht   4     0.02266828     0.00738004    -0.00160225    -0.00388195
   ht   5    -0.00443541    -0.00240718     0.00047654     0.00014615    -0.00031379
   ht   6    -0.00030685     0.00050792    -0.00016212     0.00005900     0.00001033    -0.00004661

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.910144      -5.741254E-02   0.134930       3.513240E-02  -0.104301       0.371514    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.828362       3.296200E-03   1.820621E-02   1.234285E-03   1.087880E-02   0.138023    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.91014382    -0.05741254     0.13493039     0.03513240    -0.10430147     0.37151412

 trial vector basis is being transformed.  new dimension:   3

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  6  1   -231.4365486447  2.1571E-05  4.5954E-06  3.7983E-03  1.0000E-03   
 mr-sdci #  6  2   -229.5029411048  7.1581E-01  0.0000E+00  1.3196E+00  1.0000E-04   
 mr-sdci #  6  3   -228.3246765287  2.1731E-01  0.0000E+00  1.1297E+00  1.0000E-04   
 mr-sdci #  6  4   -227.9211518134  2.6352E-02  0.0000E+00  1.7919E+00  1.0000E-04   
 mr-sdci #  6  5   -227.7274017213  3.6869E-01  0.0000E+00  1.8021E+00  1.0000E-04   
 mr-sdci #  6  6   -227.3511519298  2.9005E+00  0.0000E+00  1.3344E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.224945
time for cinew                         3.383240
time for eigenvalue solver             0.000092
time for vector access                 0.000000

          starting ci iteration   7

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:     5571191 xx:      651346 ww:      725522
One-external counts: yz :     2332820 yx:     3683468 yw:     3794578
Two-external counts: yy :     2169709 ww:      357828 xx:      351534 xz:       15843 wz:       18722 wx:      577746
Three-ext.   counts: yx :      757585 yw:      797389

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:      240809    task #     4:      120404
task #     5:     1946327    task #     6:     2786664    task #     7:     2855666    task #     8:      492288
task #     9:     4799038    task #    10:      158133    task #    11:      158133    task #    12:      147051
task #    13:      147051    task #    14:       73525    task #    15:       73527    task #    16:       73527
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1    -6.98588778
   ht   2     0.00000000    -5.05228024
   ht   3     0.00000000     0.00000000    -3.87401566
   ht   4     0.00007422     0.00219712    -0.00047645    -0.00000847

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4
 ref    1   0.910058      -6.744388E-02  -9.197736E-02  -9.317484E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4
 ref    1   0.828206       4.548677E-03   8.459835E-03   8.681551E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1     0.91005834    -0.06744388    -0.09197736    -0.09317484

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  7  1   -231.4365539467  5.3020E-06  2.0854E-06  2.4087E-03  1.0000E-03   
 mr-sdci #  7  2   -230.2295755342  7.2663E-01  0.0000E+00  1.1782E+00  1.0000E-04   
 mr-sdci #  7  3   -228.5101859618  1.8551E-01  0.0000E+00  1.2062E+00  1.0000E-04   
 mr-sdci #  7  4   -227.7928059828 -1.2835E-01  0.0000E+00  1.6605E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.142212
time for cinew                         2.112671
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   8

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:     5571191 xx:      651346 ww:      725522
One-external counts: yz :     2332820 yx:     3683468 yw:     3794578
Two-external counts: yy :     2169709 ww:      357828 xx:      351534 xz:       15843 wz:       18722 wx:      577746
Three-ext.   counts: yx :      757585 yw:      797389

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:      240809    task #     4:      120404
task #     5:     1946327    task #     6:     2786664    task #     7:     2855666    task #     8:      492288
task #     9:     4799038    task #    10:      158133    task #    11:      158133    task #    12:      147051
task #    13:      147051    task #    14:       73525    task #    15:       73527    task #    16:       73527
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.98588778
   ht   2     0.00000000    -5.05228024
   ht   3     0.00000000     0.00000000    -3.87401566
   ht   4     0.00007422     0.00219712    -0.00047645    -0.00000847
   ht   5     0.00114377    -0.00117346     0.00032570     0.00000262    -0.00000464

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.910134      -4.470120E-03   6.771298E-02  -0.167669      -3.773538E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.828344       1.998197E-05   4.585047E-03   2.811286E-02   1.423959E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.91013380    -0.00447012     0.06771298    -0.16766890    -0.03773538

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  8  1   -231.4365568042  2.8575E-06  8.2291E-07  1.5954E-03  1.0000E-03   
 mr-sdci #  8  2   -230.7891079544  5.5953E-01  0.0000E+00  7.7231E-01  1.0000E-04   
 mr-sdci #  8  3   -229.0967132869  5.8653E-01  0.0000E+00  1.0961E+00  1.0000E-04   
 mr-sdci #  8  4   -228.2246301288  4.3182E-01  0.0000E+00  1.0414E+00  1.0000E-04   
 mr-sdci #  8  5   -227.1413258836 -5.8608E-01  0.0000E+00  1.8046E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.186523
time for cinew                         2.201233
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   9

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1161301 2x:      281616 4x:       44964
All internal counts: zz :      581198 yy:     5571191 xx:      651346 ww:      725522
One-external counts: yz :     2332820 yx:     3683468 yw:     3794578
Two-external counts: yy :     2169709 ww:      357828 xx:      351534 xz:       15843 wz:       18722 wx:      577746
Three-ext.   counts: yx :      757585 yw:      797389

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:      240809    task #     4:      120404
task #     5:     1946327    task #     6:     2786664    task #     7:     2855666    task #     8:      492288
task #     9:     4799038    task #    10:      158133    task #    11:      158133    task #    12:      147051
task #    13:      147051    task #    14:       73525    task #    15:       73527    task #    16:       73527
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       99576
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.98588778
   ht   2     0.00000000    -5.05228024
   ht   3     0.00000000     0.00000000    -3.87401566
   ht   4     0.00007422     0.00219712    -0.00047645    -0.00000847
   ht   5     0.00114377    -0.00117346     0.00032570     0.00000262    -0.00000464
   ht   6    -0.00101989    -0.00007188     0.00043348     0.00000093    -0.00000055    -0.00000157

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.910077       2.484839E-02  -2.409408E-02  -9.640530E-02  -0.174310      -7.365629E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.828240       6.174423E-04   5.805248E-04   9.293982E-03   3.038410E-02   5.425249E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.91007691     0.02484839    -0.02409408    -0.09640530    -0.17431036    -0.07365629

 trial vector basis is being transformed.  new dimension:   3

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  9  1   -231.4365577719  9.6764E-07  1.5333E-07  7.0506E-04  1.0000E-03   
 mr-sdci #  9  2   -230.9692764228  1.8017E-01  0.0000E+00  4.1904E-01  1.0000E-04   
 mr-sdci #  9  3   -229.4126719405  3.1596E-01  0.0000E+00  1.0346E+00  1.0000E-04   
 mr-sdci #  9  4   -228.2940197830  6.9390E-02  0.0000E+00  1.0750E+00  1.0000E-04   
 mr-sdci #  9  5   -227.9251820385  7.8386E-01  0.0000E+00  1.3745E+00  1.0000E-04   
 mr-sdci #  9  6   -227.0799841458 -2.7117E-01  0.0000E+00  1.8674E+00  1.0000E-04   
 
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.219055
time for cinew                         3.214417
time for eigenvalue solver             0.000000
time for vector access                 0.000000

 mr-sdci  convergence criteria satisfied after  9 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  9  1   -231.4365577719  9.6764E-07  1.5333E-07  7.0506E-04  1.0000E-03   
 mr-sdci #  9  2   -230.9692764228  1.8017E-01  0.0000E+00  4.1904E-01  1.0000E-04   
 mr-sdci #  9  3   -229.4126719405  3.1596E-01  0.0000E+00  1.0346E+00  1.0000E-04   

####################CIUDGINFO####################

   ci vector at position   1 energy= -231.436557771864

################END OF CIUDGINFO################

 
    1 of the   4 expansion vectors are transformed.
    1 of the   3 matrix-vector products are transformed.

    1 expansion eigenvectors written to unit nvfile (= 11)
    1 matrix-vector products written to unit nhvfil (= 10)


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 1) =      -231.4365577719

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
                                               18

                                          orbital     3    4    5    6   42   43   44   80   81   82  110  111  139  140  141  155  171
                                              182

                                         symmetry   ag   ag   ag   ag   b3u  b3u  b3u  b2u  b2u  b2u  b1g  b1g  b1u  b1u  b1u  b2g  b3g
                                              au 

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       2  0.031640                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-      
                                                
 z*  1  1       3  0.098142                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +- 
                                                
 z*  1  1       4 -0.018853                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                
                                             +- 
 z*  1  1       5  0.028621                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +-      
                                                
 z*  1  1       6  0.038763                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -        +- 
                                                
 z*  1  1       8 -0.063074                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -   +  
                                              - 
 z*  1  1       9  0.094944                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +     - 
                                              - 
 z*  1  1      10  0.015534                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-   +-      
                                                
 z*  1  1      11  0.010132                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-        +- 
                                                
 z*  1  1      13 -0.013514                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +     -   +  
                                              - 
 z*  1  1      14  0.054249                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +    +     - 
                                              - 
 z*  1  1      15 -0.885525                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +- 
                                                
 z*  1  1      16  0.104014                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-      
                                             +- 
 z*  1  1      17  0.042348                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +- 
                                             +- 
 z*  1  1      44  0.026630                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-        +-   +- 
                                                
 z*  1  1      49  0.012609                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +     -   +-   +- 
                                                
 z*  1  1      55  0.029974                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +-   +- 
                                             +- 
 y   1  1    4937 -0.010187              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -             +-   +  
                                              - 
 y   1  1    5212 -0.010356              3( b3g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -        +-    - 
                                                
 y   1  1    5231 -0.010253              2( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -         -   +- 
                                                
 y   1  1    5232  0.038298              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -         -   +- 
                                                
 y   1  1    5234 -0.017415              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -         -   +- 
                                                
 y   1  1    5237 -0.011328              8( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -         -   +- 
                                                
 y   1  1    5411  0.015053              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -    -   +- 
                                                
 y   1  1    5499  0.010621              2( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-    - 
                                              - 
 y   1  1    5500 -0.039545              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-    - 
                                              - 
 y   1  1    5502  0.017987              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-    - 
                                              - 
 y   1  1    5505  0.011964              8( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-    - 
                                              - 
 y   1  1    5515 -0.010770              3( b3g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +               -   +- 
                                              - 
 y   1  1    7165 -0.014413              9( b2u)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +          -   +- 
                                                
 y   1  1    7834 -0.010331             11( b3u)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -   +-             +    +- 
                                              - 
 y   1  1    8943 -0.010901             12( b1g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -   +     -        +-   +- 
                                                
 y   1  1   38805 -0.012826              8( b1g)   +-   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-             +-   +  
                                              - 
 y   1  1   60679 -0.011143             12( ag )   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-   +         +-    - 
                                                
 y   1  1   60742  0.013029              8( b1g)   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-   +          -   +- 
                                                
 y   1  1   60746 -0.010532             12( b1g)   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-   +          -   +- 
                                                
 y   1  1   61353 -0.013310              8( b1g)   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-             +-   +  
                                              - 
 y   1  1   61357 -0.011462             12( b1g)   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-             +-   +  
                                              - 
 y   1  1   61410 -0.012774             12( ag )   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-             +    +- 
                                              - 
 y   1  1   62512 -0.013770              9( b2u)   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +     -        +-   +- 
                                                
 y   1  1   63646  0.014289             11( b3u)   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +              +-   +- 
                                              - 
 y   1  1  102973 -0.011763              8( b1g)   +-   +-   +-   +-   +-   +-   +-   +    +-   +-   +-   +-   +-             +-    - 
                                              - 
 y   1  1  119740  0.012691             12( ag )   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-   +-   +-   +          -   +- 
                                                
 y   1  1  217746  0.011423             12( ag )   +-   +-   +-   +-   +    +-   +-   +-   +-   +-   +-   +-   +-    -         -   +- 
                                                
 y   1  1  242531 -0.014599             11( b3u)   +-   +-   +-    -   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +  
                                              - 
 y   1  1  244857  0.011288             12( b1g)   +-   +-   +-    -   +-   +-   +-   +-   +-   +-   +-   +-   +              +-   +- 
                                              - 
 w   1  111537994  0.011082    3( b2g)   3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-      
                                                
 w   1  111544311  0.011020    3( b2g)   3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +- 
                                                
 w   1  111728386  0.011157    3( b2g)   3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +-   +- 
                                                

 ci coefficient statistics:
           rq > 0.1                2
      0.1> rq > 0.01              47
     0.01> rq > 0.001          30264
    0.001> rq > 0.0001        424547
   0.0001> rq > 0.00001      3281819
  0.00001> rq > 0.000001     7831200
 0.000001> rq               12577754
           all              24145633
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:      113824 2x:           0 4x:           0
All internal counts: zz :      581198 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23127    task #     2:       26209    task #     3:      240809    task #     4:      120404
task #     5:     1946327    task #     6:     2786664    task #     7:     2855666    task #     8:      492288
task #     9:     4799038    task #    10:      158133    task #    11:      158133    task #    12:      147051
task #    13:      147051    task #    14:       73525    task #    15:       73527    task #    16:       73527
task #    17:       80535    task #    18:      584724    task #    19:       89882    task #    20:       14353
  iref  icsf         v(icsf)             hv(icsf)
     1     1     -0.000808880141      0.186823039272
     2     2      0.031640195671     -7.303197417722
     3     3      0.098141520145    -22.652490909072
     4     4     -0.018853079636      4.353792812947
     5     5      0.028620867423     -6.607487227660
     6     6      0.038762930560     -8.946032722822
     7     7     -0.009867143660      2.278711814537
     8     8     -0.063073780524     14.560500817392
     9     9      0.094943617666    -21.913020464766
    10    10      0.015533523439     -3.586321904424
    11    11      0.010131550381     -2.337990558968
    12    12     -0.003479838084      0.803637219535
    13    13     -0.013514372632      3.119281627957
    14    14      0.054249197709    -12.522707141528
    15    15     -0.885525493513    204.345160226824
    16    16      0.104014435506    -24.007660073303
    17    17      0.042347812777     -9.775202914800
    18    18     -0.001069565099      0.247032332640
    19    19      0.000002877322     -0.000080255250
    20    20      0.000264423092     -0.061312843115
    21    21     -0.006171109383      1.426496256298
    22    22     -0.007402744441      1.711183853570
    23    23     -0.000043216939      0.009851811801
    24    24     -0.000461258396      0.106874285516
    25    25      0.000122783752     -0.028456522982
    26    26      0.000028124053     -0.006169153153
    27    27     -0.002316629893      0.535836159646
    28    28     -0.001438101109      0.331362926785
    29    29      0.007520586317     -1.737927818644
    30    30     -0.007601281043      1.756772856615
    31    31     -0.000466771732      0.107933435673
    32    32     -0.005449078262      1.258830569575
    33    33      0.001914874225     -0.442504292877
    34    34      0.000387978670     -0.089592446919
    35    35     -0.002606536720      0.602298354250
    36    36     -0.001288080325      0.300910294037
    37    37      0.002903681357     -0.671718223732
    38    38     -0.001772855452      0.409947808729
    39    39     -0.000391884837      0.090514308917
    40    40     -0.000672577859      0.155355676329
    41    41      0.000195738039     -0.045259612413
    42    42     -0.000466082583      0.107700317350
    43    43     -0.001009084709      0.233135516766
    44    44      0.026630372267     -6.146252248863
    45    45     -0.004296177278      0.992047531936
    46    46     -0.005045033309      1.165011481441
    47    47     -0.000488068107      0.112726453361
    48    48      0.000155729619     -0.035938652033
    49    49      0.012609366360     -2.910193268627
    50    50     -0.002819682957      0.651202022991
    51    51     -0.001800878341      0.415802705741
    52    52      0.008080633050     -1.865264573425
    53    53     -0.001689582558      0.390259100441
    54    54     -0.000651880707      0.150516633632
    55    55      0.029973575819     -6.917741319235

 number of reference csfs (nref) is    55.  root number (iroot) is  1.
 c0**2 =   0.82881229  c**2 (all zwalks) =   0.82956317

 pople ci energy extrapolation is computed with 30 correlated electrons.

 eref      =   -230.764561570741   "relaxed" cnot**2         =   0.828812292110
 eci       =   -231.436557771864   deltae = eci - eref       =  -0.671996201123
 eci+dv1   =   -231.551595261245   dv1 = (1-cnot**2)*deltae  =  -0.115037489381
 eci+dv2   =   -231.575355772903   dv2 = dv1 / cnot**2       =  -0.138798001038
 eci+dv3   =   -231.611486587090   dv3 = dv1 / (2*cnot**2-1) =  -0.174928815226
 eci+pople =   -231.593731470990   ( 30e- scaled deltae )    =  -0.829169900248
maximum overlap with reference    1(overlap= 0.91008)

 information on vector: 1 from unit 11 written to unit 48 filename civout              
 passed aftci ... 
 readint2: molcas,dalton2=                     0                     0
 files%faoints=aoints              
lodens (list->root)=  1
sifcfg setup: record length 4096 DP
# d1 elements per record  3272
# d2 elements per record  2730
  The MR-CISD density will be calculated.
 item #                     1 suffix=:.drt1.state1:
 read_civout: repnuc=  -224.450660864067     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  1max overlap with ref# 91% root-following 0
 MR-CISD energy:  -231.43655777    -6.98589691
 residuum:     0.00070506
 deltae:     0.00000097
 apxde:     0.00000015

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.91007691     0.02484839    -0.02409408    -0.09640530    -0.17431036    -0.07365629     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1     0.91007691     0.02484839    -0.02409408     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

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
   DYZ=   55198  DYX=   80667  DYW=   85218
   D0Z=   13462  D0Y=  131165  D0X=   20182  D0W=   22620
  DDZI=   27756 DDYI=  211458 DDXI=   33432 DDWI=   36726
  DDZE=       0 DDYE=   32226 DDXE=    6039 DDWE=    6699
================================================================================
Trace of MO density:    30.000000
   30  correlated and    12  frozen core electrons

          modens reordered block   1

               ag    1        ag    2        ag    3        ag    4        ag    5        ag    6        ag    7        ag    8
  ag    1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    3    0.00000        0.00000        1.98749      -2.804007E-05   1.251394E-03  -2.108234E-05   1.251597E-03  -8.201719E-05
  ag    4    0.00000        0.00000      -2.804007E-05    1.98055       4.059821E-06  -2.328094E-03   1.340872E-04  -8.114057E-04
  ag    5    0.00000        0.00000       1.251394E-03   4.059821E-06    1.97936       5.771481E-05   1.878498E-03  -1.403938E-04
  ag    6    0.00000        0.00000      -2.108234E-05  -2.328094E-03   5.771481E-05    1.97414      -9.198288E-05   1.721665E-03
  ag    7    0.00000        0.00000       1.251597E-03   1.340872E-04   1.878498E-03  -9.198288E-05   6.220341E-04  -1.695918E-06
  ag    8    0.00000        0.00000      -8.201719E-05  -8.114057E-04  -1.403938E-04   1.721665E-03  -1.695918E-06   2.246569E-04
  ag    9    0.00000        0.00000      -1.320444E-03   3.572056E-04   2.319499E-03  -2.428024E-04   5.953663E-04   3.735487E-06
  ag   10    0.00000        0.00000       2.777018E-04   2.464191E-03   1.353785E-04  -7.608623E-04   2.431676E-05  -4.031752E-04
  ag   11    0.00000        0.00000      -3.407702E-03  -2.822968E-04  -1.060446E-03  -1.108465E-05  -9.999439E-04  -3.316895E-06
  ag   12    0.00000        0.00000       3.421698E-05  -5.087806E-04   3.673788E-04  -4.066081E-03  -2.521495E-06  -3.678555E-04
  ag   13    0.00000        0.00000       2.724715E-04   3.021908E-03  -1.742357E-04   2.456380E-03  -2.580524E-06  -1.828613E-04
  ag   14    0.00000        0.00000      -1.494369E-03  -3.149966E-04   1.140323E-03  -1.530915E-04  -1.221062E-03   3.537791E-06
  ag   15    0.00000        0.00000       4.907173E-04   2.920594E-03  -2.930154E-04   1.281112E-03  -2.564528E-06  -6.047907E-04
  ag   16    0.00000        0.00000      -1.075249E-04  -2.572584E-03   9.265879E-04  -7.318217E-03   3.666469E-07  -5.959566E-04
  ag   17    0.00000        0.00000      -4.057853E-03  -6.197769E-04   1.018307E-03  -3.182225E-04  -7.644323E-04   1.384079E-06
  ag   18    0.00000        0.00000      -5.718804E-04  -9.167258E-06   8.911899E-04  -2.147508E-03  -3.719660E-07   3.640486E-04
  ag   19    0.00000        0.00000      -7.166301E-04   1.422921E-04   1.492640E-03  -3.799377E-04   3.887587E-04   1.957760E-06
  ag   20    0.00000        0.00000       2.062467E-03   2.868587E-05  -5.219846E-03   7.485319E-04   4.019671E-04   2.487441E-07
  ag   21    0.00000        0.00000       3.944225E-04   2.847784E-03  -1.801462E-04   2.861487E-03  -6.891512E-07  -1.589337E-04
  ag   22    0.00000        0.00000       1.391051E-04   2.576435E-03  -5.526077E-04   2.442339E-03  -2.374472E-06   3.857142E-04
  ag   23    0.00000        0.00000      -1.528363E-04   3.520858E-04  -4.499576E-04   1.881533E-03  -7.139391E-06   5.056222E-04
  ag   24    0.00000        0.00000      -4.620364E-04   1.998864E-03   4.879887E-04   2.957251E-04   1.024796E-05   9.818416E-05
  ag   25    0.00000        0.00000      -1.238227E-03   5.645264E-05   1.013382E-03  -5.633192E-04   7.232018E-04   4.619572E-06
  ag   26    0.00000        0.00000       2.899910E-05   1.093667E-03  -1.875886E-04   2.778288E-04  -6.280688E-06   1.262453E-04
  ag   27    0.00000        0.00000      -1.802056E-03  -2.966008E-04   1.713681E-03  -2.728494E-04  -4.981848E-04  -5.915494E-07
  ag   28    0.00000        0.00000       4.391980E-04   2.089858E-03  -1.540809E-04  -3.972419E-04   1.129709E-04  -3.477045E-07
  ag   29    0.00000        0.00000      -1.142747E-03   2.035644E-04   7.780147E-05  -1.262795E-04  -4.129234E-04   1.283147E-06
  ag   30    0.00000        0.00000      -2.050805E-05   2.696374E-04   4.414941E-05   3.204251E-03  -9.131809E-07  -1.268047E-04
  ag   31    0.00000        0.00000       2.330783E-03   2.891870E-05  -2.966160E-04   2.128985E-04  -1.507917E-04  -1.896789E-06
  ag   32    0.00000        0.00000      -1.132191E-03  -3.726374E-05   3.958189E-03  -3.744615E-04   1.684395E-04  -6.115054E-07
  ag   33    0.00000        0.00000       5.124644E-06  -3.291127E-03   1.091766E-04  -4.607313E-03  -9.949505E-07   1.477433E-04
  ag   34    0.00000        0.00000      -2.865704E-04  -7.654810E-04   2.529379E-05   1.216782E-04   1.792830E-07   1.530195E-04
  ag   35    0.00000        0.00000      -1.223311E-04   2.570351E-03  -5.259598E-05  -5.609782E-04   4.679889E-07  -1.052091E-06
  ag   36    0.00000        0.00000       3.161503E-05   7.009686E-04   2.995428E-04  -3.888888E-03  -1.357653E-05   5.618141E-06
  ag   37    0.00000        0.00000      -7.998313E-04   5.543838E-05  -2.809728E-03  -5.336174E-04   1.469465E-04   8.693506E-07
  ag   38    0.00000        0.00000      -2.534865E-05  -3.228136E-03   1.492393E-04  -3.174761E-03  -3.118882E-07  -1.024910E-04
  ag   39    0.00000        0.00000      -4.163618E-05   6.196795E-04   1.532665E-05   9.152073E-04   4.182914E-08   1.111021E-05

               ag    9        ag   10        ag   11        ag   12        ag   13        ag   14        ag   15        ag   16
  ag    3  -1.320444E-03   2.777018E-04  -3.407702E-03   3.421698E-05   2.724715E-04  -1.494369E-03   4.907173E-04  -1.075249E-04
  ag    4   3.572056E-04   2.464191E-03  -2.822968E-04  -5.087806E-04   3.021908E-03  -3.149966E-04   2.920594E-03  -2.572584E-03
  ag    5   2.319499E-03   1.353785E-04  -1.060446E-03   3.673788E-04  -1.742357E-04   1.140323E-03  -2.930154E-04   9.265879E-04
  ag    6  -2.428024E-04  -7.608623E-04  -1.108465E-05  -4.066081E-03   2.456380E-03  -1.530915E-04   1.281112E-03  -7.318217E-03
  ag    7   5.953663E-04   2.431676E-05  -9.999439E-04  -2.521495E-06  -2.580524E-06  -1.221062E-03  -2.564528E-06   3.666469E-07
  ag    8   3.735487E-06  -4.031752E-04  -3.316895E-06  -3.678555E-04  -1.828613E-04   3.537791E-06  -6.047907E-04  -5.959566E-04
  ag    9   2.136340E-03   2.325355E-05  -4.575554E-04  -4.579020E-06  -1.031498E-05  -2.956385E-03  -2.773339E-05  -9.071790E-06
  ag   10   2.325355E-05   8.668490E-04  -2.305754E-05   5.351261E-04   5.714977E-04  -6.639963E-05   1.594375E-03   6.217009E-04
  ag   11  -4.575554E-04  -2.305754E-05   1.842867E-03   1.355385E-05   1.456401E-05   1.696963E-03   2.968180E-05   7.994829E-06
  ag   12  -4.579020E-06   5.351261E-04   1.355385E-05   7.660395E-04   1.793315E-05   1.437670E-06   6.419186E-04   1.574167E-03
  ag   13  -1.031498E-05   5.714977E-04   1.456401E-05   1.793315E-05   7.303732E-04   7.337137E-06   1.145474E-03  -5.431617E-04
  ag   14  -2.956385E-03  -6.639963E-05   1.696963E-03   1.437670E-06   7.337137E-06   5.814678E-03   1.410286E-05   1.036638E-05
  ag   15  -2.773339E-05   1.594375E-03   2.968180E-05   6.419186E-04   1.145474E-03   1.410286E-05   3.994781E-03   1.632276E-04
  ag   16  -9.071790E-06   6.217009E-04   7.994829E-06   1.574167E-03  -5.431617E-04   1.036638E-05   1.632276E-04   3.991868E-03
  ag   17   7.002312E-04  -2.140097E-05   1.864195E-03   1.073261E-05   4.231851E-06   5.743425E-04   3.092570E-06   7.361239E-06
  ag   18   9.985447E-06  -1.332067E-03  -2.212522E-05  -2.015778E-04  -1.062010E-03   1.660321E-05  -4.747440E-03   1.089463E-03
  ag   19   9.201601E-04   9.269970E-06  -3.055528E-04  -2.503111E-06  -4.873652E-06  -6.511686E-04  -1.987003E-05   3.197240E-06
  ag   20   1.463427E-03   2.651722E-05  -7.846148E-04  -2.324001E-06  -1.047431E-05  -4.193924E-03  -3.054337E-05  -1.108367E-05
  ag   21  -3.244394E-06   5.324520E-04   1.170062E-05  -1.065003E-04   9.364002E-04  -8.004174E-06   8.547160E-04  -1.001146E-03
  ag   22   1.297580E-06  -4.464864E-04  -3.410233E-06  -1.062455E-03   4.567259E-04   9.378921E-07  -4.391785E-04  -2.827381E-03
  ag   23  -6.616175E-06  -8.626650E-04  -6.696054E-06  -1.069945E-03  -1.486395E-04   1.678163E-05  -1.146789E-03  -2.607179E-03
  ag   24   2.103268E-05  -4.210360E-04  -1.904763E-05  -9.157239E-05  -3.359088E-04  -1.604603E-05  -2.169856E-03   6.433224E-04
  ag   25   1.703484E-03   2.611665E-05  -6.395163E-04  -1.444321E-05   3.232830E-06  -1.647586E-03   1.318296E-05  -4.281612E-05
  ag   26   3.777763E-07  -1.868123E-04   1.000450E-05  -2.833592E-04  -4.753127E-05   1.400711E-05   6.299994E-05  -9.024224E-04
  ag   27  -4.123736E-04  -2.022463E-05   1.027090E-03   7.453204E-06   5.276929E-06   1.793477E-03   7.057165E-07   1.499129E-05
  ag   28  -1.539353E-04   9.865302E-05  -2.593816E-04  -2.664078E-04   5.490638E-04   1.475306E-04  -1.201117E-04  -1.094957E-03
  ag   29   5.751747E-04   1.632611E-05   9.620156E-04  -6.605803E-05   1.474421E-04  -5.759476E-04  -3.845686E-05  -2.923079E-04
  ag   30  -1.844709E-07   1.325361E-04   2.721160E-06   2.380912E-04   4.977783E-05  -8.526791E-06  -8.994943E-05   3.190788E-04
  ag   31  -6.995504E-04  -4.371485E-06   6.356149E-05   2.770889E-06   2.963272E-06   8.880485E-04   1.479778E-05   7.827837E-06
  ag   32   1.190466E-04   3.387393E-06  -1.959709E-04   1.628661E-06  -2.923658E-06   1.113823E-04  -3.792641E-06   8.416540E-06
  ag   33   2.959226E-06  -4.355083E-04  -6.076270E-06  -7.867396E-05  -4.260941E-04   5.192712E-06  -7.831337E-04   2.155674E-05
  ag   34  -1.042850E-06  -3.147948E-04  -6.637207E-06  -3.030336E-04  -2.293553E-04   5.896435E-06  -6.309421E-04  -5.391288E-04
  ag   35   1.764851E-07  -3.740863E-05  -2.305552E-06   8.243515E-06  -4.379448E-05  -1.208395E-06  -1.847504E-04   4.868691E-05
  ag   36  -3.008782E-05  -8.206171E-05   1.818974E-05   1.742778E-05  -9.638949E-05   6.394608E-05  -3.111812E-04   2.322252E-05
  ag   37   3.353741E-04  -1.319485E-06  -2.108123E-04   6.831845E-09  -1.055568E-05  -6.691525E-04  -3.282296E-05  -2.050119E-06
  ag   38  -2.400144E-06   1.759412E-04   3.656094E-06   2.175574E-04   2.645100E-05   8.069555E-07   2.999820E-04   4.773724E-04
  ag   39   1.033572E-06  -4.755333E-05  -1.290154E-06   4.218691E-05  -9.029624E-05  -1.367577E-06  -1.203974E-04   2.445572E-04

               ag   17        ag   18        ag   19        ag   20        ag   21        ag   22        ag   23        ag   24
  ag    3  -4.057853E-03  -5.718804E-04  -7.166301E-04   2.062467E-03   3.944225E-04   1.391051E-04  -1.528363E-04  -4.620364E-04
  ag    4  -6.197769E-04  -9.167258E-06   1.422921E-04   2.868587E-05   2.847784E-03   2.576435E-03   3.520858E-04   1.998864E-03
  ag    5   1.018307E-03   8.911899E-04   1.492640E-03  -5.219846E-03  -1.801462E-04  -5.526077E-04  -4.499576E-04   4.879887E-04
  ag    6  -3.182225E-04  -2.147508E-03  -3.799377E-04   7.485319E-04   2.861487E-03   2.442339E-03   1.881533E-03   2.957251E-04
  ag    7  -7.644323E-04  -3.719660E-07   3.887587E-04   4.019671E-04  -6.891512E-07  -2.374472E-06  -7.139391E-06   1.024796E-05
  ag    8   1.384079E-06   3.640486E-04   1.957760E-06   2.487441E-07  -1.589337E-04   3.857142E-04   5.056222E-04   9.818416E-05
  ag    9   7.002312E-04   9.985447E-06   9.201601E-04   1.463427E-03  -3.244394E-06   1.297580E-06  -6.616175E-06   2.103268E-05
  ag   10  -2.140097E-05  -1.332067E-03   9.269970E-06   2.651722E-05   5.324520E-04  -4.464864E-04  -8.626650E-04  -4.210360E-04
  ag   11   1.864195E-03  -2.212522E-05  -3.055528E-04  -7.846148E-04   1.170062E-05  -3.410233E-06  -6.696054E-06  -1.904763E-05
  ag   12   1.073261E-05  -2.015778E-04  -2.503111E-06  -2.324001E-06  -1.065003E-04  -1.062455E-03  -1.069945E-03  -9.157239E-05
  ag   13   4.231851E-06  -1.062010E-03  -4.873652E-06  -1.047431E-05   9.364002E-04   4.567259E-04  -1.486395E-04  -3.359088E-04
  ag   14   5.743425E-04   1.660321E-05  -6.511686E-04  -4.193924E-03  -8.004174E-06   9.378921E-07   1.678163E-05  -1.604603E-05
  ag   15   3.092570E-06  -4.747440E-03  -1.987003E-05  -3.054337E-05   8.547160E-04  -4.391785E-04  -1.146789E-03  -2.169856E-03
  ag   16   7.361239E-06   1.089463E-03   3.197240E-06  -1.108367E-05  -1.001146E-03  -2.827381E-03  -2.607179E-03   6.433224E-04
  ag   17   2.750431E-03  -2.720097E-06   4.586107E-04  -6.894388E-04   3.521097E-06  -5.506473E-06  -1.038935E-05  -1.566101E-06
  ag   18  -2.720097E-06   7.510351E-03   2.678184E-05   1.059006E-05  -6.645209E-04  -2.458441E-04   9.394691E-05   4.612869E-03
  ag   19   4.586107E-04   2.678184E-05   1.098005E-03  -5.488271E-04  -3.790465E-06  -2.223802E-06  -9.055989E-06   2.827416E-05
  ag   20  -6.894388E-04   1.059006E-05  -5.488271E-04   4.802985E-03   5.172276E-06   6.702534E-06   8.333621E-06   5.352038E-06
  ag   21   3.521097E-06  -6.645209E-04  -3.790465E-06   5.172276E-06   1.484332E-03   1.022179E-03   2.302159E-04  -2.787620E-04
  ag   22  -5.506473E-06  -2.458441E-04  -2.223802E-06   6.702534E-06   1.022179E-03   2.272410E-03   2.030975E-03  -3.371132E-04
  ag   23  -1.038935E-05   9.394691E-05  -9.055989E-06   8.333621E-06   2.302159E-04   2.030975E-03   2.682012E-03  -4.694546E-04
  ag   24  -1.566101E-06   4.612869E-03   2.827416E-05   5.352038E-06  -2.787620E-04  -3.371132E-04  -4.694546E-04   4.000944E-03
  ag   25   5.687914E-04  -6.051303E-05   1.192840E-03  -7.882661E-04   7.917527E-06   2.180184E-05   7.725444E-06  -1.828726E-05
  ag   26   1.391729E-05  -8.751906E-04  -1.579956E-06  -1.393256E-05   2.260435E-05   6.668979E-04   9.877304E-04  -9.809982E-04
  ag   27   1.282050E-03   2.113119E-05   3.543887E-04  -1.740650E-03   2.646385E-06  -5.328189E-06  -7.006867E-06   9.427910E-06
  ag   28  -3.832182E-04   4.600863E-04   4.192128E-05  -2.301419E-04   1.151682E-03   1.237497E-03   7.398658E-04   3.563722E-04
  ag   29   1.408878E-03   1.269025E-04  -1.560152E-04   8.769956E-04   3.119965E-04   3.291009E-04   1.915328E-04   9.942178E-05
  ag   30   2.009074E-07   3.044272E-04  -4.762194E-06   1.706056E-05   2.508540E-04  -3.775004E-05   2.137019E-04   1.136627E-05
  ag   31  -5.067627E-04  -1.052273E-05  -7.249376E-04  -3.633951E-04  -1.783164E-06  -7.904175E-06  -6.734317E-06  -7.443679E-06
  ag   32   5.389339E-05   8.294660E-06   3.234101E-04  -1.201494E-03  -5.139209E-06  -8.015034E-06  -1.031688E-05   1.350429E-05
  ag   33  -4.627918E-07   8.323543E-05  -9.409224E-07  -1.043727E-06  -4.142721E-04   8.535321E-05   6.664410E-04  -9.742432E-04
  ag   34  -6.357523E-06   8.792257E-04   1.167830E-06  -1.006225E-06  -3.410954E-04   2.922153E-04   7.356156E-04   1.002829E-03
  ag   35  -3.434317E-06   4.292612E-04   6.463798E-07   4.545236E-06   1.170225E-06   9.434504E-05   1.526494E-04   5.626040E-04
  ag   36   4.247603E-06   3.947271E-04  -1.091354E-05  -2.470347E-05  -3.294833E-05   1.287243E-04   2.621723E-04   5.594322E-05
  ag   37  -4.583100E-05   3.559173E-05   1.381722E-04   2.626053E-04  -2.832369E-06   1.365006E-05   2.242506E-05   1.205033E-05
  ag   38   2.536254E-06  -2.406816E-04  -2.048621E-06  -8.388730E-07  -4.877204E-05  -3.535069E-04  -4.112417E-04  -1.407604E-04
  ag   39  -5.653516E-07   2.119780E-04   6.936937E-07   3.037445E-06  -1.413155E-04  -1.865010E-04  -2.619673E-04   1.810263E-04

               ag   25        ag   26        ag   27        ag   28        ag   29        ag   30        ag   31        ag   32
  ag    3  -1.238227E-03   2.899910E-05  -1.802056E-03   4.391980E-04  -1.142747E-03  -2.050805E-05   2.330783E-03  -1.132191E-03
  ag    4   5.645264E-05   1.093667E-03  -2.966008E-04   2.089858E-03   2.035644E-04   2.696374E-04   2.891870E-05  -3.726374E-05
  ag    5   1.013382E-03  -1.875886E-04   1.713681E-03  -1.540809E-04   7.780147E-05   4.414941E-05  -2.966160E-04   3.958189E-03
  ag    6  -5.633192E-04   2.778288E-04  -2.728494E-04  -3.972419E-04  -1.262795E-04   3.204251E-03   2.128985E-04  -3.744615E-04
  ag    7   7.232018E-04  -6.280688E-06  -4.981848E-04   1.129709E-04  -4.129234E-04  -9.131809E-07  -1.507917E-04   1.684395E-04
  ag    8   4.619572E-06   1.262453E-04  -5.915494E-07  -3.477045E-07   1.283147E-06  -1.268047E-04  -1.896789E-06  -6.115054E-07
  ag    9   1.703484E-03   3.777763E-07  -4.123736E-04  -1.539353E-04   5.751747E-04  -1.844709E-07  -6.995504E-04   1.190466E-04
  ag   10   2.611665E-05  -1.868123E-04  -2.022463E-05   9.865302E-05   1.632611E-05   1.325361E-04  -4.371485E-06   3.387393E-06
  ag   11  -6.395163E-04   1.000450E-05   1.027090E-03  -2.593816E-04   9.620156E-04   2.721160E-06   6.356149E-05  -1.959709E-04
  ag   12  -1.444321E-05  -2.833592E-04   7.453204E-06  -2.664078E-04  -6.605803E-05   2.380912E-04   2.770889E-06   1.628661E-06
  ag   13   3.232830E-06  -4.753127E-05   5.276929E-06   5.490638E-04   1.474421E-04   4.977783E-05   2.963272E-06  -2.923658E-06
  ag   14  -1.647586E-03   1.400711E-05   1.793477E-03   1.475306E-04  -5.759476E-04  -8.526791E-06   8.880485E-04   1.113823E-04
  ag   15   1.318296E-05   6.299994E-05   7.057165E-07  -1.201117E-04  -3.845686E-05  -8.994943E-05   1.479778E-05  -3.792641E-06
  ag   16  -4.281612E-05  -9.024224E-04   1.499129E-05  -1.094957E-03  -2.923079E-04   3.190788E-04   7.827837E-06   8.416540E-06
  ag   17   5.687914E-04   1.391729E-05   1.282050E-03  -3.832182E-04   1.408878E-03   2.009074E-07  -5.067627E-04   5.389339E-05
  ag   18  -6.051303E-05  -8.751906E-04   2.113119E-05   4.600863E-04   1.269025E-04   3.044272E-04  -1.052273E-05   8.294660E-06
  ag   19   1.192840E-03  -1.579956E-06   3.543887E-04   4.192128E-05  -1.560152E-04  -4.762194E-06  -7.249376E-04   3.234101E-04
  ag   20  -7.882661E-04  -1.393256E-05  -1.740650E-03  -2.301419E-04   8.769956E-04   1.706056E-05  -3.633951E-04  -1.201494E-03
  ag   21   7.917527E-06   2.260435E-05   2.646385E-06   1.151682E-03   3.119965E-04   2.508540E-04  -1.783164E-06  -5.139209E-06
  ag   22   2.180184E-05   6.668979E-04  -5.328189E-06   1.237497E-03   3.291009E-04  -3.775004E-05  -7.904175E-06  -8.015034E-06
  ag   23   7.725444E-06   9.877304E-04  -7.006867E-06   7.398658E-04   1.915328E-04   2.137019E-04  -6.734317E-06  -1.031688E-05
  ag   24  -1.828726E-05  -9.809982E-04   9.427910E-06   3.563722E-04   9.942178E-05   1.136627E-05  -7.443679E-06   1.350429E-05
  ag   25   3.624546E-03   2.307852E-05  -2.482644E-04  -1.806118E-05   1.118778E-04  -7.555191E-06  -2.743216E-04   1.259201E-03
  ag   26   2.307852E-05   8.526546E-04   3.728936E-06   8.345980E-05   2.508803E-05   4.866583E-05  -3.983747E-06   1.577117E-06
  ag   27  -2.482644E-04   3.728936E-06   1.559259E-03  -5.501463E-05   2.018990E-04  -2.219500E-06  -4.188558E-04   3.702793E-04
  ag   28  -1.806118E-05   8.345980E-05  -5.501463E-05   1.545845E-03   2.311953E-05   3.016817E-04   5.035927E-05   8.905067E-05
  ag   29   1.118778E-04   2.508803E-05   2.018990E-04   2.311953E-05   1.457526E-03   8.423891E-05  -1.973013E-04  -3.332279E-04
  ag   30  -7.555191E-06   4.866583E-05  -2.219500E-06   3.016817E-04   8.423891E-05   7.599618E-04  -1.200108E-07  -6.871187E-06
  ag   31  -2.743216E-04  -3.983747E-06  -4.188558E-04   5.035927E-05  -1.973013E-04  -1.200108E-07   1.151774E-03  -2.909064E-04
  ag   32   1.259201E-03   1.577117E-06   3.702793E-04   8.905067E-05  -3.332279E-04  -6.871187E-06  -2.909064E-04   1.744860E-03
  ag   33   1.667112E-05   3.249445E-04  -5.100592E-06  -3.328951E-04  -9.076806E-05   1.018960E-04   9.814456E-07   1.336721E-06
  ag   34  -4.504529E-06  -1.913222E-05  -1.582302E-06  -3.916565E-05  -1.414153E-05   4.175748E-05   1.447204E-07   1.739863E-06
  ag   35  -8.901070E-06  -2.421646E-04   2.771154E-06   3.974906E-04   1.049427E-04   2.742345E-04  -1.176527E-06  -2.832636E-06
  ag   36  -7.309025E-05   2.353411E-04   2.864484E-05   1.742747E-04   3.455102E-05   1.399160E-04   1.438873E-05  -4.273521E-05
  ag   37   7.905179E-04   2.163633E-05  -3.031078E-04  -7.283762E-06   1.080703E-04   1.268994E-05  -1.811263E-04   4.598041E-04
  ag   38  -5.015028E-06  -2.277859E-04   2.701437E-06  -9.060511E-05  -2.251582E-05   6.086334E-05   2.968808E-06  -1.495152E-06
  ag   39  -7.522775E-06  -1.374100E-04   1.376878E-06  -2.928201E-05  -7.410845E-06   3.736321E-06  -2.929516E-07  -7.457089E-07

               ag   33        ag   34        ag   35        ag   36        ag   37        ag   38        ag   39
  ag    3   5.124644E-06  -2.865704E-04  -1.223311E-04   3.161503E-05  -7.998313E-04  -2.534865E-05  -4.163618E-05
  ag    4  -3.291127E-03  -7.654810E-04   2.570351E-03   7.009686E-04   5.543838E-05  -3.228136E-03   6.196795E-04
  ag    5   1.091766E-04   2.529379E-05  -5.259598E-05   2.995428E-04  -2.809728E-03   1.492393E-04   1.532665E-05
  ag    6  -4.607313E-03   1.216782E-04  -5.609782E-04  -3.888888E-03  -5.336174E-04  -3.174761E-03   9.152073E-04
  ag    7  -9.949505E-07   1.792830E-07   4.679889E-07  -1.357653E-05   1.469465E-04  -3.118882E-07   4.182914E-08
  ag    8   1.477433E-04   1.530195E-04  -1.052091E-06   5.618141E-06   8.693506E-07  -1.024910E-04   1.111021E-05
  ag    9   2.959226E-06  -1.042850E-06   1.764851E-07  -3.008782E-05   3.353741E-04  -2.400144E-06   1.033572E-06
  ag   10  -4.355083E-04  -3.147948E-04  -3.740863E-05  -8.206171E-05  -1.319485E-06   1.759412E-04  -4.755333E-05
  ag   11  -6.076270E-06  -6.637207E-06  -2.305552E-06   1.818974E-05  -2.108123E-04   3.656094E-06  -1.290154E-06
  ag   12  -7.867396E-05  -3.030336E-04   8.243515E-06   1.742778E-05   6.831845E-09   2.175574E-04   4.218691E-05
  ag   13  -4.260941E-04  -2.293553E-04  -4.379448E-05  -9.638949E-05  -1.055568E-05   2.645100E-05  -9.029624E-05
  ag   14   5.192712E-06   5.896435E-06  -1.208395E-06   6.394608E-05  -6.691525E-04   8.069555E-07  -1.367577E-06
  ag   15  -7.831337E-04  -6.309421E-04  -1.847504E-04  -3.111812E-04  -3.282296E-05   2.999820E-04  -1.203974E-04
  ag   16   2.155674E-05  -5.391288E-04   4.868691E-05   2.322252E-05  -2.050119E-06   4.773724E-04   2.445572E-04
  ag   17  -4.627918E-07  -6.357523E-06  -3.434317E-06   4.247603E-06  -4.583100E-05   2.536254E-06  -5.653516E-07
  ag   18   8.323543E-05   8.792257E-04   4.292612E-04   3.947271E-04   3.559173E-05  -2.406816E-04   2.119780E-04
  ag   19  -9.409224E-07   1.167830E-06   6.463798E-07  -1.091354E-05   1.381722E-04  -2.048621E-06   6.936937E-07
  ag   20  -1.043727E-06  -1.006225E-06   4.545236E-06  -2.470347E-05   2.626053E-04  -8.388730E-07   3.037445E-06
  ag   21  -4.142721E-04  -3.410954E-04   1.170225E-06  -3.294833E-05  -2.832369E-06  -4.877204E-05  -1.413155E-04
  ag   22   8.535321E-05   2.922153E-04   9.434504E-05   1.287243E-04   1.365006E-05  -3.535069E-04  -1.865010E-04
  ag   23   6.664410E-04   7.356156E-04   1.526494E-04   2.621723E-04   2.242506E-05  -4.112417E-04  -2.619673E-04
  ag   24  -9.742432E-04   1.002829E-03   5.626040E-04   5.594322E-05   1.205033E-05  -1.407604E-04   1.810263E-04
  ag   25   1.667112E-05  -4.504529E-06  -8.901070E-06  -7.309025E-05   7.905179E-04  -5.015028E-06  -7.522775E-06
  ag   26   3.249445E-04  -1.913222E-05  -2.421646E-04   2.353411E-04   2.163633E-05  -2.277859E-04  -1.374100E-04
  ag   27  -5.100592E-06  -1.582302E-06   2.771154E-06   2.864484E-05  -3.031078E-04   2.701437E-06   1.376878E-06
  ag   28  -3.328951E-04  -3.916565E-05   3.974906E-04   1.742747E-04  -7.283762E-06  -9.060511E-05  -2.928201E-05
  ag   29  -9.076806E-05  -1.414153E-05   1.049427E-04   3.455102E-05   1.080703E-04  -2.251582E-05  -7.410845E-06
  ag   30   1.018960E-04   4.175748E-05   2.742345E-04   1.399160E-04   1.268994E-05   6.086334E-05   3.736321E-06
  ag   31   9.814456E-07   1.447204E-07  -1.176527E-06   1.438873E-05  -1.811263E-04   2.968808E-06  -2.929516E-07
  ag   32   1.336721E-06   1.739863E-06  -2.832636E-06  -4.273521E-05   4.598041E-04  -1.495152E-06  -7.457089E-07
  ag   33   1.481375E-03  -1.106207E-04  -3.111510E-04   2.983895E-04   2.943254E-05   1.987563E-04  -1.956757E-04
  ag   34  -1.106207E-04   1.035304E-03   3.046881E-04   2.178158E-05   2.677218E-06  -7.156270E-05  -1.166070E-04
  ag   35  -3.111510E-04   3.046881E-04   8.712417E-04   2.221625E-05   2.503790E-06   1.238508E-04   2.336600E-04
  ag   36   2.983895E-04   2.178158E-05   2.221625E-05   4.086371E-04  -1.309023E-05   9.329108E-05  -2.134222E-05
  ag   37   2.943254E-05   2.677218E-06   2.503790E-06  -1.309023E-05   5.587607E-04   8.249955E-06  -2.303087E-06
  ag   38   1.987563E-04  -7.156270E-05   1.238508E-04   9.329108E-05   8.249955E-06   5.000443E-04   1.182442E-04
  ag   39  -1.956757E-04  -1.166070E-04   2.336600E-04  -2.134222E-05  -2.303087E-06   1.182442E-04   3.271730E-04

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98769784     1.98133120     1.97921301     1.97347927     0.01459221     0.01268763
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.01034372     0.00676976     0.00558468     0.00465853     0.00226484     0.00179614     0.00146782     0.00126820
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00092796     0.00075894     0.00047455     0.00047080     0.00029283     0.00026148     0.00013990     0.00011541
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00009769     0.00009668     0.00008110     0.00005004     0.00002169     0.00002017     0.00001515     0.00000764
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000647     0.00000560     0.00000153     0.00000119     0.00000019     0.00000010     0.00000003

          modens reordered block   1

               b3u   1        b3u   2        b3u   3        b3u   4        b3u   5        b3u   6        b3u   7        b3u   8
  b3u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   3    0.00000        0.00000        1.98485      -3.370752E-05   1.749287E-03  -2.293351E-04   1.003960E-04   1.054017E-03
  b3u   4    0.00000        0.00000      -3.370752E-05    1.97863       6.043251E-05   3.404526E-04  -1.045762E-03  -4.257446E-04
  b3u   5    0.00000        0.00000       1.749287E-03   6.043251E-05    1.97680      -2.876143E-03   7.554307E-05   2.781079E-03
  b3u   6    0.00000        0.00000      -2.293351E-04   3.404526E-04  -2.876143E-03   3.947993E-04  -1.661185E-06  -5.101481E-04
  b3u   7    0.00000        0.00000       1.003960E-04  -1.045762E-03   7.554307E-05  -1.661185E-06   1.034820E-04   3.322856E-06
  b3u   8    0.00000        0.00000       1.054017E-03  -4.257446E-04   2.781079E-03  -5.101481E-04   3.322856E-06   7.320933E-04
  b3u   9    0.00000        0.00000      -6.682121E-04  -5.439771E-04   3.763991E-03  -6.037216E-04   4.721165E-06   6.746666E-04
  b3u  10    0.00000        0.00000      -3.176987E-04   2.470273E-03  -1.513832E-04  -3.139386E-06  -3.245997E-04   1.450078E-07
  b3u  11    0.00000        0.00000      -2.242795E-03   8.155345E-05   4.186256E-04   2.676598E-04  -1.426556E-06  -5.479311E-04
  b3u  12    0.00000        0.00000       3.778240E-04  -1.908569E-03   2.455293E-04  -3.454656E-05   4.376821E-04   3.920823E-05
  b3u  13    0.00000        0.00000       1.532414E-03   1.155797E-04  -1.634011E-03   9.283932E-04   8.398994E-06  -8.982999E-04
  b3u  14    0.00000        0.00000      -4.981521E-04   1.208240E-03  -3.823676E-05   7.608054E-07  -4.283693E-04  -9.510712E-06
  b3u  15    0.00000        0.00000       4.840437E-04  -8.909429E-04   4.349221E-03  -5.848168E-04  -4.508661E-06   8.120312E-04
  b3u  16    0.00000        0.00000      -2.041049E-03  -4.800327E-04   3.634650E-03   2.496497E-04  -2.484570E-07  -4.230976E-04
  b3u  17    0.00000        0.00000       1.513150E-03  -3.956486E-04   4.875690E-04  -2.097901E-04  -5.056880E-07   5.144006E-04
  b3u  18    0.00000        0.00000      -7.684187E-04  -2.220782E-03   6.366824E-04   7.830449E-06  -3.989339E-04  -1.480803E-05
  b3u  19    0.00000        0.00000      -1.587100E-03   5.817913E-05  -8.642732E-04  -3.186088E-04  -3.594487E-06   1.641970E-04
  b3u  20    0.00000        0.00000      -9.580209E-05   6.414826E-04  -4.311945E-04   2.281330E-06  -2.678941E-04  -4.942110E-06
  b3u  21    0.00000        0.00000      -5.569649E-04  -6.881082E-04   1.811354E-03  -7.117630E-04   2.805113E-08   9.327729E-04
  b3u  22    0.00000        0.00000      -4.672741E-04  -2.961822E-03   1.011315E-03   6.928220E-07   2.214045E-04   1.053650E-06
  b3u  23    0.00000        0.00000      -5.122305E-04  -1.462729E-03   5.384198E-04  -2.171278E-06  -1.170204E-04   1.941641E-06
  b3u  24    0.00000        0.00000       1.455806E-03   2.352026E-04  -1.718838E-03   1.442502E-05   4.885557E-07   4.377736E-05
  b3u  25    0.00000        0.00000      -2.132075E-04  -1.864342E-04  -2.160269E-04  -3.529529E-04   9.483318E-07   4.812644E-04
  b3u  26    0.00000        0.00000       2.523959E-04   6.795565E-05  -1.908198E-03  -1.743146E-04   1.409775E-06   1.292287E-04
  b3u  27    0.00000        0.00000      -6.415884E-05   1.983957E-03   3.377391E-04   1.745693E-06   2.017784E-04  -2.084362E-06
  b3u  28    0.00000        0.00000      -1.094725E-03   9.264903E-05  -1.518494E-03   1.568034E-05  -4.080385E-06  -1.931281E-04
  b3u  29    0.00000        0.00000      -2.164690E-03   5.710517E-05   2.190793E-03   2.573087E-04  -2.267111E-06  -4.747391E-04
  b3u  30    0.00000        0.00000      -8.635853E-05  -1.869512E-03  -4.760344E-05  -4.722294E-06  -7.601968E-05   7.025315E-06
  b3u  31    0.00000        0.00000       5.223058E-05   9.344220E-05  -3.554786E-03   9.569527E-05  -1.052997E-06  -6.309919E-05
  b3u  32    0.00000        0.00000       1.240933E-03  -8.011774E-05   8.171475E-06  -6.643473E-05   1.227461E-07   9.620199E-05
  b3u  33    0.00000        0.00000       1.759470E-04  -1.795828E-03  -1.598777E-04   3.386174E-07   1.040422E-04  -3.328425E-08
  b3u  34    0.00000        0.00000      -1.804222E-03  -1.168689E-04  -1.021232E-03  -1.480441E-04   2.260839E-07   2.033768E-04
  b3u  35    0.00000        0.00000       2.136798E-04  -2.922181E-03   5.274229E-04  -1.136990E-06  -1.172808E-05   1.054291E-05
  b3u  36    0.00000        0.00000      -6.917655E-04  -1.090085E-03  -1.581307E-03  -2.470507E-07  -4.685242E-06  -2.278728E-05
  b3u  37    0.00000        0.00000      -2.342104E-03  -5.986843E-05   2.835558E-03  -5.080953E-05   1.242880E-07   6.620646E-05
  b3u  38    0.00000        0.00000       3.908212E-05   2.430698E-04   2.351193E-05   3.460027E-07   5.980272E-05   2.661727E-07
  b3u  39    0.00000        0.00000       2.710447E-05   9.316432E-04   1.733859E-05   3.545496E-07   1.020662E-05  -3.357779E-07

               b3u   9        b3u  10        b3u  11        b3u  12        b3u  13        b3u  14        b3u  15        b3u  16
  b3u   3  -6.682121E-04  -3.176987E-04  -2.242795E-03   3.778240E-04   1.532414E-03  -4.981521E-04   4.840437E-04  -2.041049E-03
  b3u   4  -5.439771E-04   2.470273E-03   8.155345E-05  -1.908569E-03   1.155797E-04   1.208240E-03  -8.909429E-04  -4.800327E-04
  b3u   5   3.763991E-03  -1.513832E-04   4.186256E-04   2.455293E-04  -1.634011E-03  -3.823676E-05   4.349221E-03   3.634650E-03
  b3u   6  -6.037216E-04  -3.139386E-06   2.676598E-04  -3.454656E-05   9.283932E-04   7.608054E-07  -5.848168E-04   2.496497E-04
  b3u   7   4.721165E-06  -3.245997E-04  -1.426556E-06   4.376821E-04   8.398994E-06  -4.283693E-04  -4.508661E-06  -2.484570E-07
  b3u   8   6.746666E-04   1.450078E-07  -5.479311E-04   3.920823E-05  -8.982999E-04  -9.510712E-06   8.120312E-04  -4.230976E-04
  b3u   9   1.152876E-03  -1.548335E-06  -6.008578E-05   8.661387E-05  -2.101505E-03  -5.662764E-06   9.231364E-04  -2.320075E-04
  b3u  10  -1.548335E-06   1.037350E-03  -4.265398E-07  -1.423906E-03  -4.867497E-05   1.461674E-03   2.886485E-05  -4.561243E-06
  b3u  11  -6.008578E-05  -4.265398E-07   1.047923E-03  -8.054134E-06   3.048267E-04  -4.980353E-06   2.559262E-04   1.649617E-03
  b3u  12   8.661387E-05  -1.423906E-03  -8.054134E-06   2.004438E-03  -1.544853E-04  -2.124690E-03  -3.466070E-05  -9.564763E-05
  b3u  13  -2.101505E-03  -4.867497E-05   3.048267E-04  -1.544853E-04   6.205318E-03  -1.155560E-04   5.204601E-04   3.521949E-03
  b3u  14  -5.662764E-06   1.461674E-03  -4.980353E-06  -2.124690E-03  -1.155560E-04   2.561060E-03  -9.546750E-06  -8.084999E-05
  b3u  15   9.231364E-04   2.886485E-05   2.559262E-04  -3.466070E-05   5.204601E-04  -9.546750E-06   3.215718E-03   3.094927E-03
  b3u  16  -2.320075E-04  -4.561243E-06   1.649617E-03  -9.564763E-05   3.521949E-03  -8.084999E-05   3.094927E-03   6.462296E-03
  b3u  17  -1.386055E-04   6.557507E-06  -7.582374E-04  -5.456096E-05   1.381253E-03  -6.855998E-06   8.215006E-04   4.012565E-04
  b3u  18  -2.638024E-05   1.515588E-03   6.026311E-06  -2.414797E-03  -1.709327E-05   3.506141E-03   5.139225E-05   1.089999E-05
  b3u  19   9.947136E-04   2.265624E-05   3.088968E-04   7.623488E-05  -3.152122E-03   7.719380E-05  -4.186471E-04  -1.466580E-03
  b3u  20  -1.172234E-05   8.014259E-04  -4.929652E-06  -1.045847E-03  -2.634902E-05   7.386754E-04  -3.473509E-06  -1.977446E-05
  b3u  21   1.443571E-03   1.828337E-05   2.034699E-04   5.443525E-05  -1.607976E-03  -7.499383E-06   3.063056E-03   2.185163E-03
  b3u  22   6.561198E-07  -6.117812E-04  -2.448769E-06   6.839584E-04   2.809839E-05  -9.143823E-05  -1.182212E-05  -1.436927E-05
  b3u  23  -1.235157E-06   4.594743E-04  -5.033972E-06  -7.854734E-04  -2.579272E-05   1.180339E-03   1.686605E-05  -2.054403E-05
  b3u  24  -8.868445E-05  -2.831844E-06  -6.078346E-04   3.744355E-05  -1.225382E-03   2.212972E-05  -1.462601E-03  -3.191068E-03
  b3u  25   4.862932E-04   4.258253E-06  -5.436530E-04   3.758210E-05  -1.206749E-03   1.331956E-05  -1.514388E-04  -1.484923E-03
  b3u  26   3.134264E-04  -1.752783E-06   1.385744E-04   1.813827E-05  -3.355422E-04  -7.221009E-06   8.891851E-05   2.582744E-04
  b3u  27   8.415380E-06  -7.065650E-04   9.031510E-06   1.056813E-03   7.035396E-06  -1.203296E-03  -2.179081E-05   1.293643E-06
  b3u  28   3.991975E-04   1.501614E-05   6.318598E-04   3.037393E-05  -1.460370E-03   3.591530E-05   7.620112E-06   1.590197E-04
  b3u  29  -2.826714E-04   1.167827E-06   5.789498E-04  -2.039001E-05   4.267968E-04   6.752113E-06  -5.319555E-04   4.190262E-04
  b3u  30   4.388314E-06   1.824981E-04  -7.807772E-06  -2.022580E-04  -1.394238E-05  -6.110373E-05   8.317403E-06  -4.192891E-06
  b3u  31  -2.165487E-04   8.981396E-07  -2.333740E-05  -1.859875E-05   4.778629E-04  -3.596365E-06   1.834410E-05  -6.921465E-05
  b3u  32   1.166023E-04   1.311286E-06  -2.102330E-05   5.787562E-06  -1.816838E-04   4.178226E-07   1.931114E-04  -3.511836E-04
  b3u  33   3.265603E-06  -3.566015E-04   3.018521E-06   5.453253E-04   1.307105E-05  -5.662362E-04  -9.592188E-06   3.483399E-06
  b3u  34   2.929067E-04   2.746658E-06  -6.825750E-05   1.855349E-05  -5.284223E-04  -2.955793E-07   4.294694E-04   3.011814E-05
  b3u  35  -2.735604E-05   8.384348E-05  -2.658162E-05  -1.981667E-04   1.274914E-04   3.294722E-04   3.579549E-05   3.439969E-05
  b3u  36   7.251064E-05   3.344771E-05   6.470375E-05  -6.511624E-05  -3.410787E-04   1.340404E-04  -6.565211E-05  -8.801445E-05
  b3u  37   7.593327E-05   7.867228E-07   3.155065E-05  -2.644821E-06   1.085917E-04  -3.853294E-06   2.941181E-04   4.396745E-04
  b3u  38   7.460146E-07  -2.024835E-04  -6.037515E-08   2.925033E-04   8.497332E-06  -3.313075E-04  -6.451354E-06  -1.178812E-07
  b3u  39  -4.788498E-07  -2.848893E-05  -9.949357E-08   3.570095E-05   1.576456E-06  -4.534572E-05  -1.977590E-06  -4.098930E-07

               b3u  17        b3u  18        b3u  19        b3u  20        b3u  21        b3u  22        b3u  23        b3u  24
  b3u   3   1.513150E-03  -7.684187E-04  -1.587100E-03  -9.580209E-05  -5.569649E-04  -4.672741E-04  -5.122305E-04   1.455806E-03
  b3u   4  -3.956486E-04  -2.220782E-03   5.817913E-05   6.414826E-04  -6.881082E-04  -2.961822E-03  -1.462729E-03   2.352026E-04
  b3u   5   4.875690E-04   6.366824E-04  -8.642732E-04  -4.311945E-04   1.811354E-03   1.011315E-03   5.384198E-04  -1.718838E-03
  b3u   6  -2.097901E-04   7.830449E-06  -3.186088E-04   2.281330E-06  -7.117630E-04   6.928220E-07  -2.171278E-06   1.442502E-05
  b3u   7  -5.056880E-07  -3.989339E-04  -3.594487E-06  -2.678941E-04   2.805113E-08   2.214045E-04  -1.170204E-04   4.885557E-07
  b3u   8   5.144006E-04  -1.480803E-05   1.641970E-04  -4.942110E-06   9.327729E-04   1.053650E-06   1.941641E-06   4.377736E-05
  b3u   9  -1.386055E-04  -2.638024E-05   9.947136E-04  -1.172234E-05   1.443571E-03   6.561198E-07  -1.235157E-06  -8.868445E-05
  b3u  10   6.557507E-06   1.515588E-03   2.265624E-05   8.014259E-04   1.828337E-05  -6.117812E-04   4.594743E-04  -2.831844E-06
  b3u  11  -7.582374E-04   6.026311E-06   3.088968E-04  -4.929652E-06   2.034699E-04  -2.448769E-06  -5.033972E-06  -6.078346E-04
  b3u  12  -5.456096E-05  -2.414797E-03   7.623488E-05  -1.045847E-03   5.443525E-05   6.839584E-04  -7.854734E-04   3.744355E-05
  b3u  13   1.381253E-03  -1.709327E-05  -3.152122E-03  -2.634902E-05  -1.607976E-03   2.809839E-05  -2.579272E-05  -1.225382E-03
  b3u  14  -6.855998E-06   3.506141E-03   7.719380E-05   7.386754E-04  -7.499383E-06  -9.143823E-05   1.180339E-03   2.212972E-05
  b3u  15   8.215006E-04   5.139225E-05  -4.186471E-04  -3.473509E-06   3.063056E-03  -1.182212E-05   1.686605E-05  -1.462601E-03
  b3u  16   4.012565E-04   1.089999E-05  -1.466580E-03  -1.977446E-05   2.185163E-03  -1.436927E-05  -2.054403E-05  -3.191068E-03
  b3u  17   1.324568E-03   3.981249E-05  -1.165687E-03  -9.929589E-07   1.587212E-04   1.796902E-05   1.734195E-05  -3.851264E-04
  b3u  18   3.981249E-05   6.510303E-03   6.450557E-05  -1.053483E-04   2.002193E-05   2.595949E-03   2.950925E-03  -3.855766E-05
  b3u  19  -1.165687E-03   6.450557E-05   2.015240E-03  -5.418006E-06   4.567251E-04   3.702952E-05   4.228069E-05   7.453252E-04
  b3u  20  -9.929589E-07  -1.053483E-04  -5.418006E-06   1.134624E-03  -2.383600E-05  -1.751160E-03  -3.200902E-04   1.893698E-05
  b3u  21   1.587212E-04   2.002193E-05   4.567251E-04  -2.383600E-05   4.996748E-03  -5.345757E-06   1.692418E-05  -1.593885E-03
  b3u  22   1.796902E-05   2.595949E-03   3.702952E-05  -1.751160E-03  -5.345757E-06   5.400900E-03   2.690054E-03  -2.477682E-05
  b3u  23   1.734195E-05   2.950925E-03   4.228069E-05  -3.200902E-04   1.692418E-05   2.690054E-03   2.158517E-03  -1.124478E-05
  b3u  24  -3.851264E-04  -3.855766E-05   7.453252E-04   1.893698E-05  -1.593885E-03  -2.477682E-05  -1.124478E-05   2.534266E-03
  b3u  25   2.680940E-04  -1.320752E-05   5.819074E-04   9.443001E-06  -5.499362E-04  -4.185787E-06   8.577245E-07   7.518657E-04
  b3u  26   2.716906E-05  -1.517222E-05   3.204697E-04   2.362816E-06  -6.346761E-04   2.413256E-06  -6.194278E-06  -2.250967E-04
  b3u  27  -2.225268E-05  -1.014653E-03   5.930050E-06  -7.652204E-04   1.337626E-05   1.950571E-03   4.469235E-04  -6.390820E-06
  b3u  28  -1.050544E-03   1.140074E-05   1.167749E-03   6.440343E-06   1.077903E-03  -3.830679E-05  -7.355123E-06  -1.034617E-04
  b3u  29  -5.269162E-04   1.290257E-05   2.955619E-04   6.093973E-06  -1.348701E-03   1.962331E-05   9.108415E-06   1.175892E-04
  b3u  30   6.184503E-06  -3.103487E-04  -8.354555E-06   3.364947E-04   1.452972E-05   2.417477E-04   2.968020E-04  -5.239425E-06
  b3u  31   2.161887E-05  -9.710172E-07  -5.307948E-05   3.398528E-06  -4.136455E-05   8.378021E-07   4.026450E-06   4.775379E-04
  b3u  32  -6.310336E-05  -1.381952E-06   1.668015E-04   3.448449E-06  -2.215997E-04  -1.626099E-06   3.576481E-06   9.131161E-04
  b3u  33  -7.326033E-06  -8.840580E-04  -7.181693E-06  -2.760057E-04  -1.020682E-05  -3.329723E-04  -7.589312E-04   8.694286E-06
  b3u  34   6.194684E-05  -6.012037E-06   2.466375E-04  -1.506531E-06   7.224668E-04  -5.235635E-06  -6.051683E-07  -1.565766E-04
  b3u  35   8.035388E-05   7.664266E-04  -5.332292E-05  -6.035351E-06  -1.010750E-04   1.895126E-04   3.369347E-04   2.371897E-05
  b3u  36  -1.872020E-04   2.983234E-04   1.648266E-04  -3.609677E-06   2.859322E-04   7.032085E-05   1.308515E-04  -7.575843E-05
  b3u  37   1.168406E-04   5.031516E-06  -9.170613E-05  -3.464973E-06   3.172825E-04   5.245897E-06   2.302407E-06  -4.764370E-04
  b3u  38  -1.781699E-06  -3.293179E-04  -3.965880E-06  -1.612822E-04  -3.661866E-06   2.946231E-04   5.466807E-07   3.748043E-07
  b3u  39  -1.053965E-06  -1.650565E-04  -2.752742E-06   3.886160E-05  -2.786760E-06  -3.005768E-04  -1.496783E-04   2.396659E-06

               b3u  25        b3u  26        b3u  27        b3u  28        b3u  29        b3u  30        b3u  31        b3u  32
  b3u   3  -2.132075E-04   2.523959E-04  -6.415884E-05  -1.094725E-03  -2.164690E-03  -8.635853E-05   5.223058E-05   1.240933E-03
  b3u   4  -1.864342E-04   6.795565E-05   1.983957E-03   9.264903E-05   5.710517E-05  -1.869512E-03   9.344220E-05  -8.011774E-05
  b3u   5  -2.160269E-04  -1.908198E-03   3.377391E-04  -1.518494E-03   2.190793E-03  -4.760344E-05  -3.554786E-03   8.171475E-06
  b3u   6  -3.529529E-04  -1.743146E-04   1.745693E-06   1.568034E-05   2.573087E-04  -4.722294E-06   9.569527E-05  -6.643473E-05
  b3u   7   9.483318E-07   1.409775E-06   2.017784E-04  -4.080385E-06  -2.267111E-06  -7.601968E-05  -1.052997E-06   1.227461E-07
  b3u   8   4.812644E-04   1.292287E-04  -2.084362E-06  -1.931281E-04  -4.747391E-04   7.025315E-06  -6.309919E-05   9.620199E-05
  b3u   9   4.862932E-04   3.134264E-04   8.415380E-06   3.991975E-04  -2.826714E-04   4.388314E-06  -2.165487E-04   1.166023E-04
  b3u  10   4.258253E-06  -1.752783E-06  -7.065650E-04   1.501614E-05   1.167827E-06   1.824981E-04   8.981396E-07   1.311286E-06
  b3u  11  -5.436530E-04   1.385744E-04   9.031510E-06   6.318598E-04   5.789498E-04  -7.807772E-06  -2.333740E-05  -2.102330E-05
  b3u  12   3.758210E-05   1.813827E-05   1.056813E-03   3.037393E-05  -2.039001E-05  -2.022580E-04  -1.859875E-05   5.787562E-06
  b3u  13  -1.206749E-03  -3.355422E-04   7.035396E-06  -1.460370E-03   4.267968E-04  -1.394238E-05   4.778629E-04  -1.816838E-04
  b3u  14   1.331956E-05  -7.221009E-06  -1.203296E-03   3.591530E-05   6.752113E-06  -6.110373E-05  -3.596365E-06   4.178226E-07
  b3u  15  -1.514388E-04   8.891851E-05  -2.179081E-05   7.620112E-06  -5.319555E-04   8.317403E-06   1.834410E-05   1.931114E-04
  b3u  16  -1.484923E-03   2.582744E-04   1.293643E-06   1.590197E-04   4.190262E-04  -4.192891E-06  -6.921465E-05  -3.511836E-04
  b3u  17   2.680940E-04   2.716906E-05  -2.225268E-05  -1.050544E-03  -5.269162E-04   6.184503E-06   2.161887E-05  -6.310336E-05
  b3u  18  -1.320752E-05  -1.517222E-05  -1.014653E-03   1.140074E-05   1.290257E-05  -3.103487E-04  -9.710172E-07  -1.381952E-06
  b3u  19   5.819074E-04   3.204697E-04   5.930050E-06   1.167749E-03   2.955619E-04  -8.354555E-06  -5.307948E-05   1.668015E-04
  b3u  20   9.443001E-06   2.362816E-06  -7.652204E-04   6.440343E-06   6.093973E-06   3.364947E-04   3.398528E-06   3.448449E-06
  b3u  21  -5.499362E-04  -6.346761E-04   1.337626E-05   1.077903E-03  -1.348701E-03   1.452972E-05  -4.136455E-05  -2.215997E-04
  b3u  22  -4.185787E-06   2.413256E-06   1.950571E-03  -3.830679E-05   1.962331E-05   2.417477E-04   8.378021E-07  -1.626099E-06
  b3u  23   8.577245E-07  -6.194278E-06   4.469235E-04  -7.355123E-06   9.108415E-06   2.968020E-04   4.026450E-06   3.576481E-06
  b3u  24   7.518657E-04  -2.250967E-04  -6.390820E-06  -1.034617E-04   1.175892E-04  -5.239425E-06   4.775379E-04   9.131161E-04
  b3u  25   1.398707E-03   6.300148E-04  -1.564783E-05  -4.155137E-04   1.374729E-04   3.709352E-06  -4.690642E-04   4.782040E-04
  b3u  26   6.300148E-04   1.304202E-03  -7.737239E-06  -2.838121E-04   3.887235E-04   3.689258E-06  -6.153603E-04   2.964951E-04
  b3u  27  -1.564783E-05  -7.737239E-06   1.667462E-03  -2.879297E-06   8.666622E-06   4.083425E-04   1.297978E-05  -9.878941E-06
  b3u  28  -4.155137E-04  -2.838121E-04  -2.879297E-06   1.487804E-03   1.945860E-04  -1.415156E-05   3.335396E-04  -2.843345E-04
  b3u  29   1.374729E-04   3.887235E-04   8.666622E-06   1.945860E-04   1.274821E-03  -6.364369E-06  -4.702668E-05   1.351360E-04
  b3u  30   3.709352E-06   3.689258E-06   4.083425E-04  -1.415156E-05  -6.364369E-06   6.344124E-04  -3.679851E-06   8.427349E-07
  b3u  31  -4.690642E-04  -6.153603E-04   1.297978E-05   3.335396E-04  -4.702668E-05  -3.679851E-06   1.217655E-03  -1.689294E-04
  b3u  32   4.782040E-04   2.964951E-04  -9.878941E-06  -2.843345E-04   1.351360E-04   8.427349E-07  -1.689294E-04   1.132199E-03
  b3u  33   1.820353E-06   3.053113E-06  -4.366223E-05   2.186795E-06   1.839375E-07  -2.209530E-04  -1.001518E-06   1.876918E-06
  b3u  34   4.509805E-04   9.566055E-05  -1.643139E-06   8.301502E-05  -7.217951E-05   3.888466E-06  -3.876268E-04   2.099401E-04
  b3u  35   3.129545E-05  -2.941073E-05  -3.770708E-04  -7.161962E-05   5.939627E-05  -1.545316E-04   8.779431E-05   1.818668E-05
  b3u  36  -7.833036E-05   6.287726E-05  -1.466187E-04   2.044816E-04  -1.581770E-04  -5.821537E-05  -2.255642E-04  -4.790747E-05
  b3u  37  -6.850656E-05  -4.198536E-05   2.853879E-06   1.759134E-05   1.175884E-04  -4.248807E-07  -9.652388E-05  -3.722058E-04
  b3u  38  -1.289832E-06   4.488110E-07   2.400654E-04  -4.348589E-06   5.580429E-07  -3.597325E-06   1.015724E-06  -4.928352E-07
  b3u  39   4.743598E-07   5.716716E-07  -7.987706E-05   7.127623E-07  -1.094078E-06  -6.151375E-05  -5.082636E-07   8.016731E-07

               b3u  33        b3u  34        b3u  35        b3u  36        b3u  37        b3u  38        b3u  39
  b3u   3   1.759470E-04  -1.804222E-03   2.136798E-04  -6.917655E-04  -2.342104E-03   3.908212E-05   2.710447E-05
  b3u   4  -1.795828E-03  -1.168689E-04  -2.922181E-03  -1.090085E-03  -5.986843E-05   2.430698E-04   9.316432E-04
  b3u   5  -1.598777E-04  -1.021232E-03   5.274229E-04  -1.581307E-03   2.835558E-03   2.351193E-05   1.733859E-05
  b3u   6   3.386174E-07  -1.480441E-04  -1.136990E-06  -2.470507E-07  -5.080953E-05   3.460027E-07   3.545496E-07
  b3u   7   1.040422E-04   2.260839E-07  -1.172808E-05  -4.685242E-06   1.242880E-07   5.980272E-05   1.020662E-05
  b3u   8  -3.328425E-08   2.033768E-04   1.054291E-05  -2.278728E-05   6.620646E-05   2.661727E-07  -3.357779E-07
  b3u   9   3.265603E-06   2.929067E-04  -2.735604E-05   7.251064E-05   7.593327E-05   7.460146E-07  -4.788498E-07
  b3u  10  -3.566015E-04   2.746658E-06   8.384348E-05   3.344771E-05   7.867228E-07  -2.024835E-04  -2.848893E-05
  b3u  11   3.018521E-06  -6.825750E-05  -2.658162E-05   6.470375E-05   3.155065E-05  -6.037515E-08  -9.949357E-08
  b3u  12   5.453253E-04   1.855349E-05  -1.981667E-04  -6.511624E-05  -2.644821E-06   2.925033E-04   3.570095E-05
  b3u  13   1.307105E-05  -5.284223E-04   1.274914E-04  -3.410787E-04   1.085917E-04   8.497332E-06   1.576456E-06
  b3u  14  -5.662362E-04  -2.955793E-07   3.294722E-04   1.340404E-04  -3.853294E-06  -3.313075E-04  -4.534572E-05
  b3u  15  -9.592188E-06   4.294694E-04   3.579549E-05  -6.565211E-05   2.941181E-04  -6.451354E-06  -1.977590E-06
  b3u  16   3.483399E-06   3.011814E-05   3.439969E-05  -8.801445E-05   4.396745E-04  -1.178812E-07  -4.098930E-07
  b3u  17  -7.326033E-06   6.194684E-05   8.035388E-05  -1.872020E-04   1.168406E-04  -1.781699E-06  -1.053965E-06
  b3u  18  -8.840580E-04  -6.012037E-06   7.664266E-04   2.983234E-04   5.031516E-06  -3.293179E-04  -1.650565E-04
  b3u  19  -7.181693E-06   2.466375E-04  -5.332292E-05   1.648266E-04  -9.170613E-05  -3.965880E-06  -2.752742E-06
  b3u  20  -2.760057E-04  -1.506531E-06  -6.035351E-06  -3.609677E-06  -3.464973E-06  -1.612822E-04   3.886160E-05
  b3u  21  -1.020682E-05   7.224668E-04  -1.010750E-04   2.859322E-04   3.172825E-04  -3.661866E-06  -2.786760E-06
  b3u  22  -3.329723E-04  -5.235635E-06   1.895126E-04   7.032085E-05   5.245897E-06   2.946231E-04  -3.005768E-04
  b3u  23  -7.589312E-04  -6.051683E-07   3.369347E-04   1.308515E-04   2.302407E-06   5.466807E-07  -1.496783E-04
  b3u  24   8.694286E-06  -1.565766E-04   2.371897E-05  -7.575843E-05  -4.764370E-04   3.748043E-07   2.396659E-06
  b3u  25   1.820353E-06   4.509805E-04   3.129545E-05  -7.833036E-05  -6.850656E-05  -1.289832E-06   4.743598E-07
  b3u  26   3.053113E-06   9.566055E-05  -2.941073E-05   6.287726E-05  -4.198536E-05   4.488110E-07   5.716716E-07
  b3u  27  -4.366223E-05  -1.643139E-06  -3.770708E-04  -1.466187E-04   2.853879E-06   2.400654E-04  -7.987706E-05
  b3u  28   2.186795E-06   8.301502E-05  -7.161962E-05   2.044816E-04   1.759134E-05  -4.348589E-06   7.127623E-07
  b3u  29   1.839375E-07  -7.217951E-05   5.939627E-05  -1.581770E-04   1.175884E-04   5.580429E-07  -1.094078E-06
  b3u  30  -2.209530E-04   3.888466E-06  -1.545316E-04  -5.821537E-05  -4.248807E-07  -3.597325E-06  -6.151375E-05
  b3u  31  -1.001518E-06  -3.876268E-04   8.779431E-05  -2.255642E-04  -9.652388E-05   1.015724E-06  -5.082636E-07
  b3u  32   1.876918E-06   2.099401E-04   1.818668E-05  -4.790747E-05  -3.722058E-04  -4.928352E-07   8.016731E-07
  b3u  33   8.212113E-04   6.938979E-07  -1.388781E-04  -5.505962E-05  -1.521191E-06   4.201941E-05  -9.275589E-05
  b3u  34   6.938979E-07   5.855716E-04  -3.936957E-05   1.029712E-04   2.059025E-05  -1.548986E-06  -2.988497E-07
  b3u  35  -1.388781E-04  -3.936957E-05   4.607269E-04   4.384883E-05  -4.162023E-06   7.450036E-05  -2.740134E-05
  b3u  36  -5.505962E-05   1.029712E-04   4.384883E-05   3.689277E-04   1.263648E-05   2.904043E-05  -1.075980E-05
  b3u  37  -1.521191E-06   2.059025E-05  -4.162023E-06   1.263648E-05   4.495134E-04  -3.157431E-07  -6.092511E-07
  b3u  38   4.201941E-05  -1.548986E-06   7.450036E-05   2.904043E-05  -3.157431E-07   2.537550E-04  -7.994701E-05
  b3u  39  -9.275589E-05  -2.988497E-07  -2.740134E-05  -1.075980E-05  -6.092511E-07  -7.994701E-05   1.100673E-04

Natural orbital populations,block 2
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98522815     1.97865519     1.97649100     0.01442235     0.01282474     0.01088726
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00857989     0.00499418     0.00342785     0.00213462     0.00194110     0.00126871     0.00079570     0.00071041
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00044449     0.00044156     0.00032964     0.00021544     0.00019332     0.00013171     0.00009929     0.00006318
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00003828     0.00003164     0.00002855     0.00001842     0.00001220     0.00000648     0.00000619     0.00000446
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000384     0.00000143     0.00000041     0.00000012     0.00000012     0.00000006     0.00000000

          modens reordered block   1

               b2u   1        b2u   2        b2u   3        b2u   4        b2u   5        b2u   6        b2u   7        b2u   8
  b2u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b2u   2    0.00000        1.98486       8.715511E-05  -1.737082E-03   4.105867E-04   1.645478E-04   1.430994E-03  -5.726971E-04
  b2u   3    0.00000       8.715511E-05    1.97414      -4.907983E-05   3.394190E-05  -5.666025E-03   9.248246E-05   3.021091E-05
  b2u   4    0.00000      -1.737082E-03  -4.907983E-05    1.97681      -3.252571E-03   8.290046E-06  -3.145256E-03  -4.439562E-03
  b2u   5    0.00000       4.105867E-04   3.394190E-05  -3.252571E-03   3.994766E-04   8.498523E-06   5.168243E-04   6.078747E-04
  b2u   6    0.00000       1.645478E-04  -5.666025E-03   8.290046E-06   8.498523E-06   4.720071E-04   5.597157E-06   1.075218E-05
  b2u   7    0.00000       1.430994E-03   9.248246E-05  -3.145256E-03   5.168243E-04   5.597157E-06   7.409086E-04   6.808830E-04
  b2u   8    0.00000      -5.726971E-04   3.021091E-05  -4.439562E-03   6.078747E-04   1.075218E-05   6.808830E-04   1.157595E-03
  b2u   9    0.00000       2.839171E-03  -3.823164E-04   3.519177E-04   2.691014E-04   1.139431E-05   5.496709E-04   5.806527E-05
  b2u  10    0.00000       1.678617E-03   6.970103E-04   2.841890E-03  -9.417232E-04  -6.515849E-06  -9.211926E-04  -2.122178E-03
  b2u  11    0.00000      -2.291372E-04   9.577317E-03  -1.568300E-04  -1.162381E-05  -9.884306E-04   1.670874E-06  -2.418087E-05
  b2u  12    0.00000      -9.528976E-04  -9.576932E-04   4.394880E-03  -5.795324E-04  -5.846963E-06  -8.075680E-04  -9.022560E-04
  b2u  13    0.00000       2.362560E-03  -1.812811E-03   3.082996E-03   2.556277E-04   1.316951E-05   4.334382E-04   2.414300E-04
  b2u  14    0.00000      -2.195306E-03  -9.210466E-05   4.544081E-04  -2.070622E-04  -1.139421E-05  -5.102929E-04   1.510950E-04
  b2u  15    0.00000       2.053652E-03   4.122453E-04  -3.103410E-04  -3.240504E-04  -2.964528E-06  -1.758041E-04  -1.002871E-03
  b2u  16    0.00000      -2.070614E-04   7.661431E-03  -3.373853E-04   8.523382E-07  -9.991854E-04   1.620600E-05  -1.344714E-07
  b2u  17    0.00000      -1.886337E-04  -8.137971E-04   1.560026E-03  -7.115586E-04  -1.248911E-05  -9.359655E-04  -1.434068E-03
  b2u  18    0.00000       1.288666E-03  -1.036665E-03   1.584202E-03  -1.449147E-05   2.716850E-06   4.430007E-05  -8.833737E-05
  b2u  19    0.00000      -3.363356E-04  -3.410411E-04  -3.766269E-04   3.545087E-04   1.002175E-05   4.836352E-04   4.850480E-04
  b2u  20    0.00000       7.304551E-05   1.088762E-03   1.531615E-04  -5.863234E-06   4.012203E-04  -1.165650E-05  -9.436586E-06
  b2u  21    0.00000      -2.522769E-04   2.223308E-05   1.283492E-03   1.733404E-04   1.145076E-05   1.281081E-04   3.106236E-04
  b2u  22    0.00000      -1.287552E-03   5.508331E-05   1.565326E-03  -1.658634E-05  -3.712835E-06  -1.927750E-04   4.013872E-04
  b2u  23    0.00000       9.203535E-05  -8.772212E-03   2.283964E-04   6.400137E-07   5.579023E-04  -6.864912E-06   1.000011E-06
  b2u  24    0.00000       2.693063E-03  -1.130279E-04   2.322939E-03   2.581889E-04   8.560836E-06   4.760926E-04   2.811303E-04
  b2u  25    0.00000      -2.668666E-04   2.297770E-04  -3.925103E-03   9.728974E-05   1.134992E-06   6.495307E-05   2.194262E-04
  b2u  26    0.00000       9.215649E-04  -1.014372E-04  -1.439407E-04   6.663650E-05   1.058017E-07   9.705879E-05   1.153789E-04
  b2u  27    0.00000       7.939659E-05  -4.532359E-03   7.729818E-05   1.439078E-06   3.043750E-05   1.958418E-06   2.430566E-06
  b2u  28    0.00000       1.780013E-03  -4.025736E-05  -8.896290E-04  -1.483392E-04  -2.655226E-06  -2.039889E-04  -2.921885E-04
  b2u  29    0.00000      -7.386762E-04   1.947506E-05   1.519377E-03  -1.660116E-06  -1.247516E-07  -2.663413E-05   7.539868E-05
  b2u  30    0.00000       2.213035E-03  -2.424790E-04   2.844216E-03  -5.084458E-05  -7.880104E-08  -6.549378E-05  -7.613625E-05

               b2u   9        b2u  10        b2u  11        b2u  12        b2u  13        b2u  14        b2u  15        b2u  16
  b2u   2   2.839171E-03   1.678617E-03  -2.291372E-04  -9.528976E-04   2.362560E-03  -2.195306E-03   2.053652E-03  -2.070614E-04
  b2u   3  -3.823164E-04   6.970103E-04   9.577317E-03  -9.576932E-04  -1.812811E-03  -9.210466E-05   4.122453E-04   7.661431E-03
  b2u   4   3.519177E-04   2.841890E-03  -1.568300E-04   4.394880E-03   3.082996E-03   4.544081E-04  -3.103410E-04  -3.373853E-04
  b2u   5   2.691014E-04  -9.417232E-04  -1.162381E-05  -5.795324E-04   2.556277E-04  -2.070622E-04  -3.240504E-04   8.523382E-07
  b2u   6   1.139431E-05  -6.515849E-06  -9.884306E-04  -5.846963E-06   1.316951E-05  -1.139421E-05  -2.964528E-06  -9.991854E-04
  b2u   7   5.496709E-04  -9.211926E-04   1.670874E-06  -8.075680E-04   4.334382E-04  -5.102929E-04  -1.758041E-04   1.620600E-05
  b2u   8   5.806527E-05  -2.122178E-03  -2.418087E-05  -9.022560E-04   2.414300E-04   1.510950E-04  -1.002871E-03  -1.344714E-07
  b2u   9   1.045635E-03  -2.989828E-04  -5.344164E-06   2.500460E-04   1.645592E-03  -7.585942E-04   3.069771E-04  -8.782688E-07
  b2u  10  -2.989828E-04   6.212838E-03   4.005748E-05  -5.349181E-04  -3.508818E-03  -1.398067E-03   3.145850E-03  -1.325881E-05
  b2u  11  -5.344164E-06   4.005748E-05   2.160138E-03  -6.180665E-06  -3.098004E-05  -4.629287E-06   3.004462E-05   2.242559E-03
  b2u  12   2.500460E-04  -5.349181E-04  -6.180665E-06   3.195680E-03   3.093983E-03   8.299308E-04  -4.264027E-04  -1.873992E-05
  b2u  13   1.645592E-03  -3.508818E-03  -3.098004E-05   3.093983E-03   6.461463E-03   4.114898E-04  -1.458316E-03  -1.275715E-05
  b2u  14  -7.585942E-04  -1.398067E-03  -4.629287E-06   8.299308E-04   4.114898E-04   1.333854E-03  -1.169118E-03   7.431861E-07
  b2u  15   3.069771E-04   3.145850E-03   3.004462E-05  -4.264027E-04  -1.458316E-03  -1.169118E-03   1.999617E-03   6.039906E-06
  b2u  16  -8.782688E-07  -1.325881E-05   2.242559E-03  -1.873992E-05  -1.275715E-05   7.431861E-07   6.039906E-06   2.458689E-03
  b2u  17   2.055279E-04   1.591876E-03   2.432462E-05   3.054209E-03   2.210873E-03   1.639831E-04   4.477547E-04  -5.896904E-06
  b2u  18   6.065358E-04  -1.217475E-03  -8.636218E-06   1.468316E-03   3.197466E-03   3.922816E-04  -7.403516E-04  -5.189401E-06
  b2u  19   5.444387E-04  -1.212034E-03  -1.485084E-05   1.699560E-04   1.502144E-03  -2.601509E-04  -5.853211E-04  -4.316554E-06
  b2u  20  -1.617661E-06   2.718669E-05  -8.775334E-04   2.211904E-07  -1.161484E-05  -1.630696E-06   1.291398E-05  -1.077019E-03
  b2u  21  -1.406842E-04  -3.296845E-04  -2.405328E-05  -8.464427E-05  -2.641977E-04  -2.518315E-05  -3.157325E-04  -2.533614E-05
  b2u  22  -6.344058E-04  -1.455055E-03  -1.406795E-05   2.075822E-06  -1.657252E-04   1.055883E-03  -1.160392E-03  -3.861640E-06
  b2u  23   4.924572E-06   4.127026E-06  -1.456571E-03   1.399017E-05   1.798289E-05  -2.684009E-06  -5.149780E-06  -1.719103E-03
  b2u  24   5.758131E-04  -4.261744E-04  -6.802583E-06  -5.346355E-04   4.113399E-04  -5.260165E-04   2.912208E-04   2.312838E-06
  b2u  25  -2.380952E-05  -4.834782E-04  -5.341851E-06   2.474663E-05  -6.201552E-05   2.363645E-05  -5.692911E-05  -1.374922E-06
  b2u  26   2.242914E-05  -1.809816E-04   1.743515E-06  -1.911080E-04   3.538144E-04   6.377531E-05  -1.652452E-04   8.669073E-06
  b2u  27   2.599746E-06  -6.863369E-06  -1.052849E-04   1.580553E-07   1.765898E-05   1.541495E-06  -6.131500E-06  -3.022965E-04
  b2u  28  -6.804667E-05   5.323759E-04   5.839480E-06   4.209109E-04   2.623377E-05   5.849207E-05   2.486124E-04  -1.341793E-06
  b2u  29  -7.070556E-05  -3.595617E-04  -4.097931E-06   8.095988E-05   9.707281E-05   2.055692E-04  -1.705406E-04  -1.938449E-06
  b2u  30   3.217846E-05  -1.014683E-04  -1.445786E-06   2.889294E-04   4.351195E-04   1.139863E-04  -8.753239E-05  -2.012691E-06

               b2u  17        b2u  18        b2u  19        b2u  20        b2u  21        b2u  22        b2u  23        b2u  24
  b2u   2  -1.886337E-04   1.288666E-03  -3.363356E-04   7.304551E-05  -2.522769E-04  -1.287552E-03   9.203535E-05   2.693063E-03
  b2u   3  -8.137971E-04  -1.036665E-03  -3.410411E-04   1.088762E-03   2.223308E-05   5.508331E-05  -8.772212E-03  -1.130279E-04
  b2u   4   1.560026E-03   1.584202E-03  -3.766269E-04   1.531615E-04   1.283492E-03   1.565326E-03   2.283964E-04   2.322939E-03
  b2u   5  -7.115586E-04  -1.449147E-05   3.545087E-04  -5.863234E-06   1.733404E-04  -1.658634E-05   6.400137E-07   2.581889E-04
  b2u   6  -1.248911E-05   2.716850E-06   1.002175E-05   4.012203E-04   1.145076E-05  -3.712835E-06   5.579023E-04   8.560836E-06
  b2u   7  -9.359655E-04   4.430007E-05   4.836352E-04  -1.165650E-05   1.281081E-04  -1.927750E-04  -6.864912E-06   4.760926E-04
  b2u   8  -1.434068E-03  -8.833737E-05   4.850480E-04  -9.436586E-06   3.106236E-04   4.013872E-04   1.000011E-06   2.811303E-04
  b2u   9   2.055279E-04   6.065358E-04   5.444387E-04  -1.617661E-06  -1.406842E-04  -6.344058E-04   4.924572E-06   5.758131E-04
  b2u  10   1.591876E-03  -1.217475E-03  -1.212034E-03   2.718669E-05  -3.296845E-04  -1.455055E-03   4.127026E-06  -4.261744E-04
  b2u  11   2.432462E-05  -8.636218E-06  -1.485084E-05  -8.775334E-04  -2.405328E-05  -1.406795E-05  -1.456571E-03  -6.802583E-06
  b2u  12   3.054209E-03   1.468316E-03   1.699560E-04   2.211904E-07  -8.464427E-05   2.075822E-06   1.399017E-05  -5.346355E-04
  b2u  13   2.210873E-03   3.197466E-03   1.502144E-03  -1.161484E-05  -2.641977E-04  -1.657252E-04   1.798289E-05   4.113399E-04
  b2u  14   1.639831E-04   3.922816E-04  -2.601509E-04  -1.630696E-06  -2.518315E-05   1.055883E-03  -2.684009E-06  -5.260165E-04
  b2u  15   4.477547E-04  -7.403516E-04  -5.853211E-04   1.291398E-05  -3.157325E-04  -1.160392E-03  -5.149780E-06   2.912208E-04
  b2u  16  -5.896904E-06  -5.189401E-06  -4.316554E-06  -1.077019E-03  -2.533614E-05  -3.861640E-06  -1.719103E-03   2.312838E-06
  b2u  17   5.000255E-03   1.618892E-03   5.774704E-04  -2.121762E-05   6.373767E-04  -1.069706E-03   4.964846E-06  -1.355967E-03
  b2u  18   1.618892E-03   2.545524E-03   7.643835E-04  -5.913513E-06  -2.258038E-04  -1.077400E-04   8.397764E-06  -1.270279E-04
  b2u  19   5.774704E-04   7.643835E-04   1.409325E-03  -1.985415E-05   6.298665E-04  -4.208491E-04   3.885413E-06  -1.460294E-04
  b2u  20  -2.121762E-05  -5.913513E-06  -1.985415E-05   9.347348E-04  -9.710096E-06   9.514070E-06   5.025396E-04   1.016596E-05
  b2u  21   6.373767E-04  -2.258038E-04   6.298665E-04  -9.710096E-06   1.300211E-03  -2.842336E-04   1.098271E-05  -3.920522E-04
  b2u  22  -1.069706E-03  -1.077400E-04  -4.208491E-04   9.514070E-06  -2.842336E-04   1.487238E-03  -2.248087E-07  -1.931506E-04
  b2u  23   4.964846E-06   8.397764E-06   3.885413E-06   5.025396E-04   1.098271E-05  -2.248087E-07   2.075393E-03  -1.361262E-06
  b2u  24  -1.355967E-03  -1.270279E-04  -1.460294E-04   1.016596E-05  -3.920522E-04  -1.931506E-04  -1.361262E-06   1.274575E-03
  b2u  25  -3.267749E-05  -4.713184E-04   4.746895E-04  -1.450326E-05   6.182238E-04  -3.359425E-04   1.311026E-06  -5.047526E-05
  b2u  26   2.260829E-04   9.181634E-04   4.807095E-04  -1.649236E-05   2.930807E-04  -2.851055E-04  -1.288763E-05  -1.381919E-04
  b2u  27   9.692596E-06   2.737673E-05   1.427955E-05   2.261398E-04   1.149639E-05  -6.733068E-06   5.232293E-04  -4.302795E-06
  b2u  28   7.168875E-04   1.581426E-04  -4.484344E-04   4.960065E-06  -9.478149E-05  -8.198366E-05   1.546393E-07  -7.211708E-05
  b2u  29  -2.941270E-04  -7.557934E-05  -8.721274E-05  -2.736398E-07   6.839025E-05   2.148921E-04   1.067529E-06   1.696686E-04
  b2u  30   3.173034E-04   4.764906E-04   7.255146E-05  -1.864433E-06   4.012281E-05  -2.015055E-05   1.549013E-06   1.146203E-04

               b2u  25        b2u  26        b2u  27        b2u  28        b2u  29        b2u  30
  b2u   2  -2.668666E-04   9.215649E-04   7.939659E-05   1.780013E-03  -7.386762E-04   2.213035E-03
  b2u   3   2.297770E-04  -1.014372E-04  -4.532359E-03  -4.025736E-05   1.947506E-05  -2.424790E-04
  b2u   4  -3.925103E-03  -1.439407E-04   7.729818E-05  -8.896290E-04   1.519377E-03   2.844216E-03
  b2u   5   9.728974E-05   6.663650E-05   1.439078E-06  -1.483392E-04  -1.660116E-06  -5.084458E-05
  b2u   6   1.134992E-06   1.058017E-07   3.043750E-05  -2.655226E-06  -1.247516E-07  -7.880104E-08
  b2u   7   6.495307E-05   9.705879E-05   1.958418E-06  -2.039889E-04  -2.663413E-05  -6.549378E-05
  b2u   8   2.194262E-04   1.153789E-04   2.430566E-06  -2.921885E-04   7.539868E-05  -7.613625E-05
  b2u   9  -2.380952E-05   2.242914E-05   2.599746E-06  -6.804667E-05  -7.070556E-05   3.217846E-05
  b2u  10  -4.834782E-04  -1.809816E-04  -6.863369E-06   5.323759E-04  -3.595617E-04  -1.014683E-04
  b2u  11  -5.341851E-06   1.743515E-06  -1.052849E-04   5.839480E-06  -4.097931E-06  -1.445786E-06
  b2u  12   2.474663E-05  -1.911080E-04   1.580553E-07   4.209109E-04   8.095988E-05   2.889294E-04
  b2u  13  -6.201552E-05   3.538144E-04   1.765898E-05   2.623377E-05   9.707281E-05   4.351195E-04
  b2u  14   2.363645E-05   6.377531E-05   1.541495E-06   5.849207E-05   2.055692E-04   1.139863E-04
  b2u  15  -5.692911E-05  -1.652452E-04  -6.131500E-06   2.486124E-04  -1.705406E-04  -8.753239E-05
  b2u  16  -1.374922E-06   8.669073E-06  -3.022965E-04  -1.341793E-06  -1.938449E-06  -2.012691E-06
  b2u  17  -3.267749E-05   2.260829E-04   9.692596E-06   7.168875E-04  -2.941270E-04   3.173034E-04
  b2u  18  -4.713184E-04   9.181634E-04   2.737673E-05   1.581426E-04  -7.557934E-05   4.764906E-04
  b2u  19   4.746895E-04   4.807095E-04   1.427955E-05  -4.484344E-04  -8.721274E-05   7.255146E-05
  b2u  20  -1.450326E-05  -1.649236E-05   2.261398E-04   4.960065E-06  -2.736398E-07  -1.864433E-06
  b2u  21   6.182238E-04   2.930807E-04   1.149639E-05  -9.478149E-05   6.839025E-05   4.012281E-05
  b2u  22  -3.359425E-04  -2.851055E-04  -6.733068E-06  -8.198366E-05   2.148921E-04  -2.015055E-05
  b2u  23   1.311026E-06  -1.288763E-05   5.232293E-04   1.546393E-07   1.067529E-06   1.549013E-06
  b2u  24  -5.047526E-05  -1.381919E-04  -4.302795E-06  -7.211708E-05   1.696686E-04   1.146203E-04
  b2u  25   1.221178E-03   1.696781E-04   3.774812E-06  -3.901722E-04   2.389673E-04  -9.577015E-05
  b2u  26   1.696781E-04   1.132669E-03   1.198691E-05  -2.070353E-04  -5.119409E-05   3.736623E-04
  b2u  27   3.774812E-06   1.198691E-05   6.608320E-04  -4.577023E-06  -1.008853E-06   9.262961E-06
  b2u  28  -3.901722E-04  -2.070353E-04  -4.577023E-06   5.842463E-04  -1.081717E-04   2.008504E-05
  b2u  29   2.389673E-04  -5.119409E-05  -1.008853E-06  -1.081717E-04   3.495532E-04  -1.296368E-05
  b2u  30  -9.577015E-05   3.736623E-04   9.262961E-06   2.008504E-05  -1.296368E-05   4.482211E-04

Natural orbital populations,block 3
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98523923     1.97651748     1.97428565     0.01444704     0.01087243     0.00659837     0.00498697
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00342947     0.00213618     0.00126700     0.00110802     0.00079230     0.00067695     0.00044324     0.00032965
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00021101     0.00019312     0.00013117     0.00009917     0.00003161     0.00002853     0.00002068     0.00001215
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00000647     0.00000435     0.00000268     0.00000142     0.00000012     0.00000006

          modens reordered block   1

               b1g   1        b1g   2        b1g   3        b1g   4        b1g   5        b1g   6        b1g   7        b1g   8
  b1g   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b1g   2    0.00000        1.98060       2.272239E-03  -1.134518E-03  -3.043083E-03  -7.099825E-06  -6.565175E-06   3.255254E-03
  b1g   3    0.00000       2.272239E-03    1.97409      -1.904095E-03  -1.034828E-03  -4.308417E-03   2.303905E-05  -2.260484E-03
  b1g   4    0.00000      -1.134518E-03  -1.904095E-03   2.258584E-04   4.056601E-04   3.695663E-04   2.160657E-07  -1.833569E-04
  b1g   5    0.00000      -3.043083E-03  -1.034828E-03   4.056601E-04   8.693765E-04   5.397508E-04  -1.453670E-07  -5.697008E-04
  b1g   6    0.00000      -7.099825E-06  -4.308417E-03   3.695663E-04   5.397508E-04   7.687419E-04   9.220866E-07  -1.957478E-05
  b1g   7    0.00000      -6.565175E-06   2.303905E-05   2.160657E-07  -1.453670E-07   9.220866E-07   9.452892E-06   9.084338E-07
  b1g   8    0.00000       3.255254E-03  -2.260484E-03  -1.833569E-04  -5.697008E-04  -1.957478E-05   9.084338E-07   7.254535E-04
  b1g   9    0.00000       4.191619E-03  -1.511881E-03  -6.103170E-04  -1.602297E-03  -6.512706E-04   1.417948E-06   1.141966E-03
  b1g  10    0.00000       2.063800E-03  -7.851252E-03   5.945219E-04   6.197490E-04   1.574774E-03   2.870418E-06   5.459709E-04
  b1g  11    0.00000       1.645076E-03  -3.610209E-03  -3.647482E-04  -1.332038E-03  -2.010792E-04   4.250737E-06   1.057815E-03
  b1g  12    0.00000      -1.820299E-04   5.262745E-04   5.706225E-06   2.197013E-05   1.992162E-06   9.135322E-05  -1.805678E-05
  b1g  13    0.00000      -3.041705E-03   2.536852E-03   1.605473E-04   5.343673E-04  -1.043549E-04  -1.225372E-06  -9.344467E-04
  b1g  14    0.00000       2.265313E-03  -2.663180E-03   3.883796E-04   4.526789E-04   1.069371E-03   1.612039E-06   4.560788E-04
  b1g  15    0.00000      -5.124766E-04   1.318160E-03  -2.061887E-06  -1.156629E-06  -6.822715E-06   1.954106E-04  -2.076276E-06
  b1g  16    0.00000       8.682634E-06  -2.326008E-03   5.034512E-04   8.613269E-04   1.063013E-03   2.410777E-06  -1.525353E-04
  b1g  17    0.00000       9.084807E-04   1.169929E-03   9.740144E-05   4.194340E-04   9.092951E-05  -8.494482E-09  -3.317968E-04
  b1g  18    0.00000       4.165195E-04  -8.036296E-04  -3.643921E-07  -3.163742E-06   8.299943E-07  -6.610422E-05   2.987728E-06
  b1g  19    0.00000       1.183628E-03  -7.452361E-04   1.263922E-04   1.859287E-04   2.858978E-04   2.760777E-07  -4.350281E-05
  b1g  20    0.00000       2.011074E-03   7.138155E-04  -2.236782E-06  -1.005944E-04   2.697429E-04   9.607994E-07   5.671181E-04
  b1g  21    0.00000       3.193618E-04  -3.154119E-03  -1.260911E-04  -1.322058E-04  -2.378070E-04  -4.173484E-07   4.698564E-05
  b1g  22    0.00000       2.031331E-04  -2.161869E-04   1.328146E-07   4.210652E-07   2.800804E-07  -1.740159E-06  -9.724434E-07
  b1g  23    0.00000       3.149252E-03  -4.008697E-03  -1.489277E-04  -4.373814E-04  -7.965794E-05   3.157132E-07   4.263617E-04
  b1g  24    0.00000      -9.646270E-04   3.583690E-05   1.533316E-04   3.153739E-04   3.042977E-04   3.893850E-07  -2.275324E-04
  b1g  25    0.00000      -1.121058E-04   5.312770E-05   9.855827E-08  -1.808099E-08   8.197498E-07  -3.945675E-05   4.226274E-07
  b1g  26    0.00000       2.485545E-03   7.811250E-04  -2.232290E-06   3.540761E-05  -1.034805E-05   1.913345E-08  -4.262157E-05
  b1g  27    0.00000       7.094348E-04   3.717391E-03   5.363058E-06   8.094907E-05  -1.783748E-05  -1.331309E-07  -9.628472E-05
  b1g  28    0.00000       3.004980E-03  -3.064230E-03   1.016021E-04   1.745376E-04   2.156106E-04   7.512398E-08  -2.644225E-05
  b1g  29    0.00000      -1.068795E-04   1.912920E-05   1.100877E-07   1.569063E-07   2.880512E-07  -4.752285E-06   1.624989E-07
  b1g  30    0.00000      -5.742605E-04   8.779383E-04  -1.075749E-05  -4.647239E-05   4.291754E-05   2.054716E-07   8.974792E-05

               b1g   9        b1g  10        b1g  11        b1g  12        b1g  13        b1g  14        b1g  15        b1g  16
  b1g   2   4.191619E-03   2.063800E-03   1.645076E-03  -1.820299E-04  -3.041705E-03   2.265313E-03  -5.124766E-04   8.682634E-06
  b1g   3  -1.511881E-03  -7.851252E-03  -3.610209E-03   5.262745E-04   2.536852E-03  -2.663180E-03   1.318160E-03  -2.326008E-03
  b1g   4  -6.103170E-04   5.945219E-04  -3.647482E-04   5.706225E-06   1.605473E-04   3.883796E-04  -2.061887E-06   5.034512E-04
  b1g   5  -1.602297E-03   6.197490E-04  -1.332038E-03   2.197013E-05   5.343673E-04   4.526789E-04  -1.156629E-06   8.613269E-04
  b1g   6  -6.512706E-04   1.574774E-03  -2.010792E-04   1.992162E-06  -1.043549E-04   1.069371E-03  -6.822715E-06   1.063013E-03
  b1g   7   1.417948E-06   2.870418E-06   4.250737E-06   9.135322E-05  -1.225372E-06   1.612039E-06   1.954106E-04   2.410777E-06
  b1g   8   1.141966E-03   5.459709E-04   1.057815E-03  -1.805678E-05  -9.344467E-04   4.560788E-04  -2.076276E-06  -1.525353E-04
  b1g   9   4.011897E-03  -1.528472E-04   4.751641E-03  -7.679172E-05  -8.661646E-04  -4.475631E-04   2.553295E-06  -1.149349E-03
  b1g  10  -1.528472E-04   4.000866E-03   1.126128E-03  -2.079157E-05  -1.006360E-03   2.845193E-03  -1.490014E-05   2.581169E-03
  b1g  11   4.751641E-03   1.126128E-03   7.506104E-03  -1.005055E-04  -6.876695E-04   2.603119E-04   5.311543E-05  -8.695815E-05
  b1g  12  -7.679172E-05  -2.079157E-05  -1.005055E-04   9.926669E-04   1.131649E-05  -9.337395E-06   2.382423E-03   1.449673E-05
  b1g  13  -8.661646E-04  -1.006360E-03  -6.876695E-04   1.131649E-05   1.482134E-03  -1.024806E-03   1.119834E-08  -2.210715E-04
  b1g  14  -4.475631E-04   2.845193E-03   2.603119E-04  -9.337395E-06  -1.024806E-03   2.298412E-03  -1.926375E-05   2.026484E-03
  b1g  15   2.553295E-06  -1.490014E-05   5.311543E-05   2.382423E-03   1.119834E-08  -1.926375E-05   6.419992E-03   2.779396E-05
  b1g  16  -1.149349E-03   2.581169E-03  -8.695815E-05   1.449673E-05  -2.210715E-04   2.026484E-03   2.779396E-05   2.645399E-03
  b1g  17  -2.164714E-03  -6.570548E-04  -4.596882E-03   8.425148E-05   2.903741E-04  -3.451233E-04   2.098719E-05  -4.680679E-04
  b1g  18   1.238833E-05   5.843944E-06  -1.331987E-06  -1.124968E-03  -7.792617E-07   8.679421E-06  -3.770961E-03  -1.971870E-05
  b1g  19   6.242163E-05   9.130236E-04   8.758853E-04  -1.707961E-05  -3.401854E-05   6.828389E-04  -9.258310E-06   9.862334E-04
  b1g  20  -1.308640E-04   1.120580E-03  -4.860803E-04   8.898953E-06  -1.189579E-03   1.277829E-03   1.120670E-06   7.433003E-04
  b1g  21  -8.243424E-05  -3.191231E-04  -2.874255E-04   4.463933E-06  -2.430188E-04  -3.917167E-05   1.302526E-06   2.150020E-04
  b1g  22  -1.015894E-06   9.508089E-07  -6.443109E-06  -2.204061E-04   1.941093E-06   2.141761E-06  -1.114037E-03  -2.861376E-06
  b1g  23   7.930044E-04   2.923629E-05   1.017620E-04  -4.699644E-06  -4.149530E-04  -8.573534E-05  -9.469407E-06  -6.627722E-04
  b1g  24  -6.357449E-04   5.370155E-04  -8.850907E-04   1.726755E-05   3.408517E-04   2.957931E-04   4.212377E-06   7.301821E-04
  b1g  25  -9.508815E-07   1.577134E-06  -7.690688E-06  -2.357273E-04  -4.088368E-07   8.254708E-07  -6.158073E-05  -2.696445E-08
  b1g  26  -1.811267E-04  -5.207256E-05  -4.229209E-04   7.999468E-06  -1.134718E-06   9.389870E-05   9.634918E-07   1.506535E-04
  b1g  27  -3.130314E-04  -2.668765E-05  -3.982647E-04   7.859100E-06   3.426831E-05   1.281057E-04   4.586304E-06   2.589787E-04
  b1g  28  -2.949931E-04   4.727629E-04  -2.286460E-04   1.107190E-06  -4.680658E-05   3.521643E-04  -9.209175E-06   4.049604E-04
  b1g  29  -4.305395E-07   4.640720E-07  -3.221998E-07   2.599594E-05  -4.576619E-07  -1.208488E-07   2.505153E-04   1.430050E-06
  b1g  30   1.163253E-04   2.457263E-04   2.049006E-04  -3.256179E-06  -1.420820E-04   1.889109E-04  -3.884121E-07   2.598761E-04

               b1g  17        b1g  18        b1g  19        b1g  20        b1g  21        b1g  22        b1g  23        b1g  24
  b1g   2   9.084807E-04   4.165195E-04   1.183628E-03   2.011074E-03   3.193618E-04   2.031331E-04   3.149252E-03  -9.646270E-04
  b1g   3   1.169929E-03  -8.036296E-04  -7.452361E-04   7.138155E-04  -3.154119E-03  -2.161869E-04  -4.008697E-03   3.583690E-05
  b1g   4   9.740144E-05  -3.643921E-07   1.263922E-04  -2.236782E-06  -1.260911E-04   1.328146E-07  -1.489277E-04   1.533316E-04
  b1g   5   4.194340E-04  -3.163742E-06   1.859287E-04  -1.005944E-04  -1.322058E-04   4.210652E-07  -4.373814E-04   3.153739E-04
  b1g   6   9.092951E-05   8.299943E-07   2.858978E-04   2.697429E-04  -2.378070E-04   2.800804E-07  -7.965794E-05   3.042977E-04
  b1g   7  -8.494482E-09  -6.610422E-05   2.760777E-07   9.607994E-07  -4.173484E-07  -1.740159E-06   3.157132E-07   3.893850E-07
  b1g   8  -3.317968E-04   2.987728E-06  -4.350281E-05   5.671181E-04   4.698564E-05  -9.724434E-07   4.263617E-04  -2.275324E-04
  b1g   9  -2.164714E-03   1.238833E-05   6.242163E-05  -1.308640E-04  -8.243424E-05  -1.015894E-06   7.930044E-04  -6.357449E-04
  b1g  10  -6.570548E-04   5.843944E-06   9.130236E-04   1.120580E-03  -3.191231E-04   9.508089E-07   2.923629E-05   5.370155E-04
  b1g  11  -4.596882E-03  -1.331987E-06   8.758853E-04  -4.860803E-04  -2.874255E-04  -6.443109E-06   1.017620E-04  -8.850907E-04
  b1g  12   8.425148E-05  -1.124968E-03  -1.707961E-05   8.898953E-06   4.463933E-06  -2.204061E-04  -4.699644E-06   1.726755E-05
  b1g  13   2.903741E-04  -7.792617E-07  -3.401854E-05  -1.189579E-03  -2.430188E-04   1.941093E-06  -4.149530E-04   3.408517E-04
  b1g  14  -3.451233E-04   8.679421E-06   6.828389E-04   1.277829E-03  -3.917167E-05   2.141761E-06  -8.573534E-05   2.957931E-04
  b1g  15   2.098719E-05  -3.770961E-03  -9.258310E-06   1.120670E-06   1.302526E-06  -1.114037E-03  -9.469407E-06   4.212377E-06
  b1g  16  -4.680679E-04  -1.971870E-05   9.862334E-04   7.433003E-04   2.150020E-04  -2.861376E-06  -6.627722E-04   7.301821E-04
  b1g  17   3.979294E-03  -3.209376E-05  -9.773582E-04   3.803123E-04  -1.278005E-06  -2.094601E-06   9.575183E-04   1.008522E-03
  b1g  18  -3.209376E-05   2.978403E-03   8.903181E-06  -4.462924E-06  -3.434332E-06   1.207577E-03   1.390049E-06  -1.182885E-05
  b1g  19  -9.773582E-04   8.903181E-06   8.530312E-04   8.948768E-05   5.445211E-05   2.916100E-06  -3.175418E-04  -2.274218E-05
  b1g  20   3.803123E-04  -4.462924E-06   8.948768E-05   1.551044E-03   3.069415E-04  -1.571554E-06   3.468331E-04  -3.721695E-05
  b1g  21  -1.278005E-06  -3.434332E-06   5.445211E-05   3.069415E-04   7.552029E-04   7.262620E-07  -1.059123E-04   4.185251E-05
  b1g  22  -2.094601E-06   1.207577E-03   2.916100E-06  -1.571554E-06   7.262620E-07   6.595690E-04   1.317162E-06   3.319385E-07
  b1g  23   9.575183E-04   1.390049E-06  -3.175418E-04   3.468331E-04  -1.059123E-04   1.317162E-06   1.475825E-03   1.122672E-04
  b1g  24   1.008522E-03  -1.182885E-05  -2.274218E-05  -3.721695E-05   4.185251E-05   3.319385E-07   1.122672E-04   1.036276E-03
  b1g  25   3.658189E-06  -6.101733E-04  -1.119288E-06   1.142012E-06   1.694794E-06  -4.669291E-04   3.296561E-06   3.770647E-06
  b1g  26   5.549220E-04  -4.168717E-06  -2.376167E-04   4.151353E-04   2.715893E-04   7.831220E-07   3.062937E-04   3.046821E-04
  b1g  27   5.873841E-05  -4.368362E-06   2.369898E-04   1.743391E-04   1.387063E-04  -8.942852E-07  -2.985784E-04   2.090116E-05
  b1g  28   1.324798E-04   4.845491E-06   2.283228E-04   8.738254E-05  -6.013007E-05   2.563205E-06   1.989939E-04   7.020015E-05
  b1g  29   2.083220E-06  -4.089709E-04  -9.888881E-07   6.060128E-07   1.110081E-06  -3.102285E-04   9.073324E-07   1.549836E-06
  b1g  30  -1.741789E-04   3.632849E-07   1.358387E-04   3.014517E-05  -3.047631E-06   1.576913E-07  -1.933964E-04   1.161301E-04

               b1g  25        b1g  26        b1g  27        b1g  28        b1g  29        b1g  30
  b1g   2  -1.121058E-04   2.485545E-03   7.094348E-04   3.004980E-03  -1.068795E-04  -5.742605E-04
  b1g   3   5.312770E-05   7.811250E-04   3.717391E-03  -3.064230E-03   1.912920E-05   8.779383E-04
  b1g   4   9.855827E-08  -2.232290E-06   5.363058E-06   1.016021E-04   1.100877E-07  -1.075749E-05
  b1g   5  -1.808099E-08   3.540761E-05   8.094907E-05   1.745376E-04   1.569063E-07  -4.647239E-05
  b1g   6   8.197498E-07  -1.034805E-05  -1.783748E-05   2.156106E-04   2.880512E-07   4.291754E-05
  b1g   7  -3.945675E-05   1.913345E-08  -1.331309E-07   7.512398E-08  -4.752285E-06   2.054716E-07
  b1g   8   4.226274E-07  -4.262157E-05  -9.628472E-05  -2.644225E-05   1.624989E-07   8.974792E-05
  b1g   9  -9.508815E-07  -1.811267E-04  -3.130314E-04  -2.949931E-04  -4.305395E-07   1.163253E-04
  b1g  10   1.577134E-06  -5.207256E-05  -2.668765E-05   4.727629E-04   4.640720E-07   2.457263E-04
  b1g  11  -7.690688E-06  -4.229209E-04  -3.982647E-04  -2.286460E-04  -3.221998E-07   2.049006E-04
  b1g  12  -2.357273E-04   7.999468E-06   7.859100E-06   1.107190E-06   2.599594E-05  -3.256179E-06
  b1g  13  -4.088368E-07  -1.134718E-06   3.426831E-05  -4.680658E-05  -4.576619E-07  -1.420820E-04
  b1g  14   8.254708E-07   9.389870E-05   1.281057E-04   3.521643E-04  -1.208488E-07   1.889109E-04
  b1g  15  -6.158073E-05   9.634918E-07   4.586304E-06  -9.209175E-06   2.505153E-04  -3.884121E-07
  b1g  16  -2.696445E-08   1.506535E-04   2.589787E-04   4.049604E-04   1.430050E-06   2.598761E-04
  b1g  17   3.658189E-06   5.549220E-04   5.873841E-05   1.324798E-04   2.083220E-06  -1.741789E-04
  b1g  18  -6.101733E-04  -4.168717E-06  -4.368362E-06   4.845491E-06  -4.089709E-04   3.632849E-07
  b1g  19  -1.119288E-06  -2.376167E-04   2.369898E-04   2.283228E-04  -9.888881E-07   1.358387E-04
  b1g  20   1.142012E-06   4.151353E-04   1.743391E-04   8.738254E-05   6.060128E-07   3.014517E-05
  b1g  21   1.694794E-06   2.715893E-04   1.387063E-04  -6.013007E-05   1.110081E-06  -3.047631E-06
  b1g  22  -4.669291E-04   7.831220E-07  -8.942852E-07   2.563205E-06  -3.102285E-04   1.576913E-07
  b1g  23   3.296561E-06   3.062937E-04  -2.985784E-04   1.989939E-04   9.073324E-07  -1.933964E-04
  b1g  24   3.770647E-06   3.046821E-04   2.090116E-05   7.020015E-05   1.549836E-06   1.161301E-04
  b1g  25   6.244149E-04   1.837619E-06   1.646852E-07  -3.539872E-07   2.886423E-04  -1.913459E-07
  b1g  26   1.837619E-06   8.693599E-04   2.022379E-05  -1.257550E-04   7.709529E-07  -2.326979E-04
  b1g  27   1.646852E-07   2.022379E-05   4.053621E-04  -9.348386E-05   4.007150E-08   2.242833E-05
  b1g  28  -3.539872E-07  -1.257550E-04  -9.348386E-05   4.991294E-04  -3.817092E-07   1.181921E-04
  b1g  29   2.886423E-04   7.709529E-07   4.007150E-08  -3.817092E-07   3.361305E-04  -3.799575E-08
  b1g  30  -1.913459E-07  -2.326979E-04   2.242833E-05   1.181921E-04  -3.799575E-08   3.262462E-04

Natural orbital populations,block 4
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98134819     1.97347934     0.01458152     0.01033177     0.00994142     0.00465277     0.00226045
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00177381     0.00146680     0.00092678     0.00076073     0.00047076     0.00029216     0.00021873     0.00013954
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00009743     0.00008099     0.00006359     0.00005004     0.00002012     0.00001900     0.00000762     0.00000647
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00000154     0.00000118     0.00000107     0.00000019     0.00000003     0.00000002

          modens reordered block   1

               b1u   1        b1u   2        b1u   3        b1u   4        b1u   5        b1u   6        b1u   7        b1u   8
  b1u   1    1.96513       2.121384E-03   1.247412E-03   4.581088E-04   4.062886E-03   1.153608E-04  -6.028558E-03   1.782067E-04
  b1u   2   2.121384E-03   6.365757E-02   2.387448E-02  -4.729524E-03  -3.691716E-04  -3.678991E-03   1.305246E-05   3.556877E-03
  b1u   3   1.247412E-03   2.387448E-02   1.183934E-02  -1.660861E-03   3.923832E-04  -1.448960E-03   7.263627E-04   1.423557E-03
  b1u   4   4.581088E-04  -4.729524E-03  -1.660861E-03   6.266560E-04   9.542763E-05   5.893686E-04   1.480810E-05  -5.606660E-04
  b1u   5   4.062886E-03  -3.691716E-04   3.923832E-04   9.542763E-05   2.086923E-03   4.011502E-05  -1.644283E-03  -2.942260E-05
  b1u   6   1.153608E-04  -3.678991E-03  -1.448960E-03   5.893686E-04   4.011502E-05   9.105625E-04  -1.778712E-05  -4.265373E-04
  b1u   7  -6.028558E-03   1.305246E-05   7.263627E-04   1.480810E-05  -1.644283E-03  -1.778712E-05   2.219008E-03   1.525045E-05
  b1u   8   1.782067E-04   3.556877E-03   1.423557E-03  -5.606660E-04  -2.942260E-05  -4.265373E-04   1.525045E-05   1.349892E-03
  b1u   9   6.040046E-03   5.801065E-04  -2.145497E-03  -2.040129E-04  -6.300503E-04  -3.249233E-05  -7.437285E-04   2.571558E-05
  b1u  10   8.199296E-05  -3.097798E-03  -1.238838E-03   5.751635E-04   2.907985E-05   1.099352E-03  -1.217145E-05  -6.425255E-05
  b1u  11  -5.782070E-06   3.462379E-03   1.363143E-03  -5.280413E-04  -2.789839E-05  -2.374821E-04   1.415668E-05   9.475142E-04
  b1u  12   3.246811E-03  -1.964415E-04   3.615699E-04   5.916171E-05   1.772121E-03   2.057211E-05  -1.410444E-03  -1.470244E-05
  b1u  13  -1.459526E-04   7.351288E-04   2.822654E-04  -1.383310E-04  -9.863078E-06  -4.035687E-04   1.422504E-05  -9.289853E-04
  b1u  14   2.230380E-05   2.297675E-04   1.436649E-04  -1.690843E-04  -7.080726E-06  -6.827307E-04   2.173426E-06   3.606792E-04
  b1u  15   2.755680E-03  -2.822385E-05  -5.851782E-05   6.364843E-06   2.653173E-04   4.610468E-06  -6.413028E-04  -1.103860E-06
  b1u  16   1.573898E-05   7.013586E-05   1.512283E-05   2.165365E-05   3.084894E-06   4.736040E-05  -5.436917E-06  -1.681130E-04

               b1u   9        b1u  10        b1u  11        b1u  12        b1u  13        b1u  14        b1u  15        b1u  16
  b1u   1   6.040046E-03   8.199296E-05  -5.782070E-06   3.246811E-03  -1.459526E-04   2.230380E-05   2.755680E-03   1.573898E-05
  b1u   2   5.801065E-04  -3.097798E-03   3.462379E-03  -1.964415E-04   7.351288E-04   2.297675E-04  -2.822385E-05   7.013586E-05
  b1u   3  -2.145497E-03  -1.238838E-03   1.363143E-03   3.615699E-04   2.822654E-04   1.436649E-04  -5.851782E-05   1.512283E-05
  b1u   4  -2.040129E-04   5.751635E-04  -5.280413E-04   5.916171E-05  -1.383310E-04  -1.690843E-04   6.364843E-06   2.165365E-05
  b1u   5  -6.300503E-04   2.907985E-05  -2.789839E-05   1.772121E-03  -9.863078E-06  -7.080726E-06   2.653173E-04   3.084894E-06
  b1u   6  -3.249233E-05   1.099352E-03  -2.374821E-04   2.057211E-05  -4.035687E-04  -6.827307E-04   4.610468E-06   4.736040E-05
  b1u   7  -7.437285E-04  -1.217145E-05   1.415668E-05  -1.410444E-03   1.422504E-05   2.173426E-06  -6.413028E-04  -5.436917E-06
  b1u   8   2.571558E-05  -6.425255E-05   9.475142E-04  -1.470244E-05  -9.289853E-04   3.606792E-04  -1.103860E-06  -1.681130E-04
  b1u   9   2.880498E-03  -2.199050E-05   2.380840E-05  -4.886190E-04  -2.225333E-05   8.907970E-06   7.356037E-05  -2.914232E-08
  b1u  10  -2.199050E-05   1.634127E-03   8.226798E-05   1.109279E-05  -1.114942E-03  -9.006547E-04   3.465900E-06   1.621066E-04
  b1u  11   2.380840E-05   8.226798E-05   1.209723E-03  -1.220669E-05  -6.723006E-04  -1.561259E-04  -1.620779E-06   8.804642E-05
  b1u  12  -4.886190E-04   1.109279E-05  -1.220669E-05   2.029925E-03  -4.047203E-06  -6.687362E-06   2.658534E-04   3.018571E-06
  b1u  13  -2.225333E-05  -1.114942E-03  -6.723006E-04  -4.047203E-06   2.089234E-03  -4.677871E-05  -4.827006E-06  -2.892834E-05
  b1u  14   8.907970E-06  -9.006547E-04  -1.561259E-04  -6.687362E-06  -4.677871E-05   1.250353E-03   1.027177E-06  -2.887338E-04
  b1u  15   7.356037E-05   3.465900E-06  -1.620779E-06   2.658534E-04  -4.827006E-06   1.027177E-06   6.132201E-04   1.812117E-06
  b1u  16  -2.914232E-08   1.621066E-04   8.804642E-05   3.018571E-06  -2.892834E-05  -2.887338E-04   1.812117E-06   4.658644E-04

Natural orbital populations,block 5
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.96518556     0.07412123     0.00562952     0.00512613     0.00384875     0.00271508     0.00104791     0.00070034
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00046840     0.00046410     0.00034126     0.00016796     0.00010107     0.00004227     0.00002387     0.00000686

          modens reordered block   1

               b2g   1        b2g   2        b2g   3        b2g   4        b2g   5        b2g   6        b2g   7        b2g   8
  b2g   1    1.91564      -5.279886E-03  -5.456724E-04   1.788270E-03   3.363464E-03  -6.053172E-04   1.008331E-03  -4.529191E-03
  b2g   2  -5.279886E-03   1.143945E-03   9.725140E-06  -5.961924E-05  -9.762491E-04   3.003690E-05  -1.218182E-03  -1.082583E-04
  b2g   3  -5.456724E-04   9.725140E-06   9.208210E-04  -3.460831E-03   1.179096E-05   1.455980E-03   7.205464E-06  -8.974074E-06
  b2g   4   1.788270E-03  -5.961924E-05  -3.460831E-03   1.334054E-02  -2.749117E-05  -5.971394E-03  -4.330794E-06   4.284599E-05
  b2g   5   3.363464E-03  -9.762491E-04   1.179096E-05  -2.749117E-05   1.365581E-03   9.352576E-06   1.436858E-03  -6.034866E-04
  b2g   6  -6.053172E-04   3.003690E-05   1.455980E-03  -5.971394E-03   9.352576E-06   3.089520E-03  -1.783107E-06  -2.433849E-05
  b2g   7   1.008331E-03  -1.218182E-03   7.205464E-06  -4.330794E-06   1.436858E-03  -1.783107E-06   1.822744E-03  -3.155309E-04
  b2g   8  -4.529191E-03  -1.082583E-04  -8.974074E-06   4.284599E-05  -6.034866E-04  -2.433849E-05  -3.155309E-04   1.299019E-03
  b2g   9  -4.527702E-04  -5.601931E-06   9.093084E-04  -4.021313E-03   1.468280E-05   2.400139E-03   1.677826E-05   8.168585E-06
  b2g  10  -6.017996E-03  -9.914743E-04  -2.165416E-05   1.227712E-04   5.445605E-04  -7.866693E-05   8.738499E-04   7.032820E-04
  b2g  11  -2.450557E-04   1.321009E-05   4.905501E-04  -1.862999E-03  -2.487667E-07   7.313123E-04  -5.247565E-06  -6.426357E-06
  b2g  12   3.065641E-03  -4.866818E-04   3.660809E-06  -8.629063E-06   1.166891E-03   5.052277E-06   1.167516E-03  -1.178478E-03
  b2g  13   1.821577E-05   5.223901E-06  -7.979053E-05  -1.374234E-04  -6.215395E-06   6.532006E-04  -7.159343E-06  -1.629338E-06
  b2g  14  -2.755769E-03   3.351617E-06  -2.536721E-06   1.256337E-05   6.106658E-05  -9.082827E-06  -8.257014E-05   2.273081E-04
  b2g  15  -7.536290E-04  -5.468037E-05  -2.776816E-06   1.280521E-05   1.087890E-04  -7.421350E-06   3.678348E-04   2.238451E-04
  b2g  16   1.493311E-04  -6.152306E-08  -2.073230E-05   3.781098E-05  -1.301677E-06   4.057112E-05  -2.313927E-06   6.668566E-07

               b2g   9        b2g  10        b2g  11        b2g  12        b2g  13        b2g  14        b2g  15        b2g  16
  b2g   1  -4.527702E-04  -6.017996E-03  -2.450557E-04   3.065641E-03   1.821577E-05  -2.755769E-03  -7.536290E-04   1.493311E-04
  b2g   2  -5.601931E-06  -9.914743E-04   1.321009E-05  -4.866818E-04   5.223901E-06   3.351617E-06  -5.468037E-05  -6.152306E-08
  b2g   3   9.093084E-04  -2.165416E-05   4.905501E-04   3.660809E-06  -7.979053E-05  -2.536721E-06  -2.776816E-06  -2.073230E-05
  b2g   4  -4.021313E-03   1.227712E-04  -1.862999E-03  -8.629063E-06  -1.374234E-04   1.256337E-05   1.280521E-05   3.781098E-05
  b2g   5   1.468280E-05   5.445605E-04  -2.487667E-07   1.166891E-03  -6.215395E-06   6.106658E-05   1.087890E-04  -1.301677E-06
  b2g   6   2.400139E-03  -7.866693E-05   7.313123E-04   5.052277E-06   6.532006E-04  -9.082827E-06  -7.421350E-06   4.057112E-05
  b2g   7   1.677826E-05   8.738499E-04  -5.247565E-06   1.167516E-03  -7.159343E-06  -8.257014E-05   3.678348E-04  -2.313927E-06
  b2g   8   8.168585E-06   7.032820E-04  -6.426357E-06  -1.178478E-03  -1.629338E-06   2.273081E-04   2.238451E-04   6.668566E-07
  b2g   9   2.126840E-03  -2.369590E-05   4.589775E-04  -5.758696E-06   9.581449E-04  -1.491151E-06  -9.071749E-06   1.502193E-04
  b2g  10  -2.369590E-05   1.565146E-03  -2.001416E-05  -4.739105E-05  -3.314835E-05   1.957731E-04  -1.004846E-04  -4.230155E-06
  b2g  11   4.589775E-04  -2.001416E-05   5.676012E-04  -1.515723E-06  -2.453455E-04  -1.111793E-06  -1.542518E-06  -5.452320E-05
  b2g  12  -5.758696E-06  -4.739105E-05  -1.515723E-06   1.706640E-03  -4.844601E-06  -2.341209E-04  -9.909741E-06  -1.661321E-06
  b2g  13   9.581449E-04  -3.314835E-05  -2.453455E-04  -4.844601E-06   1.123929E-03  -2.600299E-06  -1.051827E-06   2.841176E-04
  b2g  14  -1.491151E-06   1.957731E-04  -1.111793E-06  -2.341209E-04  -2.600299E-06   6.310345E-04   1.688308E-04  -2.470010E-08
  b2g  15  -9.071749E-06  -1.004846E-04  -1.542518E-06  -9.909741E-06  -1.051827E-06   1.688308E-04   9.108979E-04  -1.456763E-06
  b2g  16   1.502193E-04  -4.230155E-06  -5.452320E-05  -1.661321E-06   2.841176E-04  -2.470010E-08  -1.456763E-06   3.495863E-04

Natural orbital populations,block 6
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.91570189     0.01861417     0.00535375     0.00294762     0.00221752     0.00111017     0.00060789     0.00035053
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00028169     0.00024289     0.00005696     0.00005369     0.00004952     0.00001086     0.00000418     0.00000036

          modens reordered block   1

               b3g   1        b3g   2        b3g   3        b3g   4        b3g   5        b3g   6        b3g   7        b3g   8
  b3g   1    1.91604       7.791965E-03   3.418190E-03  -7.592497E-04   1.074437E-03  -3.623930E-03  -5.219327E-03   9.975304E-04
  b3g   2   7.791965E-03   1.138554E-03   9.719472E-04  -3.968673E-07   1.208339E-03   1.089566E-04   9.632319E-04   2.675084E-06
  b3g   3   3.418190E-03   9.719472E-04   1.361716E-03   1.254143E-06   1.433219E-03  -5.972709E-04   5.436550E-04   6.808857E-08
  b3g   4  -7.592497E-04  -3.968673E-07   1.254143E-06   2.436529E-03  -5.443755E-06  -6.699527E-07   5.239000E-06  -2.588325E-03
  b3g   5   1.074437E-03   1.208339E-03   1.433219E-03  -5.443755E-06   1.819033E-03  -3.075650E-04   8.665215E-04   8.510508E-06
  b3g   6  -3.623930E-03   1.089566E-04  -5.972709E-04  -6.699527E-07  -3.075650E-04   1.305612E-03   6.937471E-04   7.282881E-07
  b3g   7  -5.219327E-03   9.632319E-04   5.436550E-04   5.239000E-06   8.665215E-04   6.937471E-04   1.523486E-03  -3.069171E-06
  b3g   8   9.975304E-04   2.675084E-06   6.808857E-08  -2.588325E-03   8.510508E-06   7.282881E-07  -3.069171E-06   3.392732E-03
  b3g   9   2.604802E-03   4.904490E-04   1.165995E-03   2.946814E-06   1.168322E-03  -1.172635E-03  -3.041267E-05  -2.614841E-06
  b3g  10  -2.298831E-03  -9.532652E-06   5.795706E-05  -2.543055E-06  -9.154800E-05   2.196287E-04   1.864323E-04   2.041120E-06
  b3g  11  -8.853147E-04   4.992066E-05   1.064180E-04  -6.104168E-06   3.617377E-04   2.257950E-04  -1.058234E-04   7.444716E-06

               b3g   9        b3g  10        b3g  11
  b3g   1   2.604802E-03  -2.298831E-03  -8.853147E-04
  b3g   2   4.904490E-04  -9.532652E-06   4.992066E-05
  b3g   3   1.165995E-03   5.795706E-05   1.064180E-04
  b3g   4   2.946814E-06  -2.543055E-06  -6.104168E-06
  b3g   5   1.168322E-03  -9.154800E-05   3.617377E-04
  b3g   6  -1.172635E-03   2.196287E-04   2.257950E-04
  b3g   7  -3.041267E-05   1.864323E-04  -1.058234E-04
  b3g   8  -2.614841E-06   2.041120E-06   7.444716E-06
  b3g   9   1.702002E-03  -2.308070E-04  -1.488116E-05
  b3g  10  -2.308070E-04   6.216931E-04   1.636217E-04
  b3g  11  -1.488116E-05   1.636217E-04   9.177665E-04

Natural orbital populations,block 7
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.91610220     0.00554599     0.00532864     0.00291317     0.00111213     0.00060670     0.00028252     0.00024203
              MO     9       MO    10       MO    11
  occ(*)=     0.00005672     0.00005339     0.00001079

          modens reordered block   1

               au    1        au    2        au    3        au    4        au    5        au    6        au    7        au    8
  au    1   7.149307E-02  -5.510183E-03  -4.122160E-03   4.149184E-03   3.364710E-03   9.784627E-06  -3.855365E-03  -6.824407E-04
  au    2  -5.510183E-03   7.152588E-04   6.463297E-04  -6.329649E-04  -6.142546E-04  -7.207676E-07   5.833569E-04   1.436418E-04
  au    3  -4.122160E-03   6.463297E-04   9.411524E-04  -4.672912E-04  -1.117688E-03   8.883499E-07   2.603393E-04   4.069292E-04
  au    4   4.149184E-03  -6.329649E-04  -4.672912E-04   1.401710E-03   9.559578E-05   2.726409E-06  -9.826288E-04   9.264823E-04
  au    5   3.364710E-03  -6.142546E-04  -1.117688E-03   9.559578E-05   1.638653E-03  -2.515691E-06   7.007896E-05  -1.112044E-03
  au    6   9.784627E-06  -7.207676E-07   8.883499E-07   2.726409E-06  -2.515691E-06   7.853242E-05  -2.569317E-06   3.394543E-06
  au    7  -3.855365E-03   5.833569E-04   2.603393E-04  -9.826288E-04   7.007896E-05  -2.569317E-06   1.229649E-03  -6.778234E-04
  au    8  -6.824407E-04   1.436418E-04   4.069292E-04   9.264823E-04  -1.112044E-03   3.394543E-06  -6.778234E-04   2.089126E-03
  au    9   2.967857E-04  -1.703530E-04  -6.833234E-04   3.689709E-04   8.983121E-04  -1.404108E-06   1.539617E-04   4.854720E-05
  au   10   7.303057E-05   2.005039E-05   4.565715E-05  -1.688653E-04  -1.592839E-04   2.463536E-07  -8.730505E-05   2.716765E-05
  au   11   1.584034E-05  -1.176426E-06  -3.374513E-07   4.187202E-06  -1.373149E-06   1.545269E-04  -2.529623E-06   5.067912E-06

               au    9        au   10        au   11
  au    1   2.967857E-04   7.303057E-05   1.584034E-05
  au    2  -1.703530E-04   2.005039E-05  -1.176426E-06
  au    3  -6.833234E-04   4.565715E-05  -3.374513E-07
  au    4   3.689709E-04  -1.688653E-04   4.187202E-06
  au    5   8.983121E-04  -1.592839E-04  -1.373149E-06
  au    6  -1.404108E-06   2.463536E-07   1.545269E-04
  au    7   1.539617E-04  -8.730505E-05  -2.529623E-06
  au    8   4.854720E-05   2.716765E-05   5.067912E-06
  au    9   1.248161E-03  -2.872156E-04   1.155534E-06
  au   10  -2.872156E-04   4.642419E-04  -3.868101E-07
  au   11   1.155534E-06  -3.868101E-07   5.515406E-04

Natural orbital populations,block 8
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     0.07280623     0.00384812     0.00272822     0.00106612     0.00059752     0.00046756     0.00016888     0.00010237
              MO     9       MO    10       MO    11
  occ(*)=     0.00003252     0.00002660     0.00000695


 total number of electrons =   42.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        ag  partial gross atomic populations
   ao class       1ag        2ag        3ag        4ag        5ag        6ag 
    C1_ s       0.024115   1.965539   0.526173   0.467569  -0.010767   0.061887
    C1_ p       0.003468  -0.000704   0.036777   0.131341   0.351952   0.261809
    C1_ d       0.000470   0.000179   0.004833  -0.113661   0.038915  -0.022982
    C2_ s       1.965551   0.024164   1.073193   0.233856  -0.019274   0.030864
    C2_ p       0.000801   0.007425   0.076595   0.116824   0.712136   1.291405
    C2_ d      -0.000546   0.001050   0.009970  -0.053412   0.078436   0.011475
    H1_ s       0.001179   0.000814   0.080091   0.776869   0.277225   0.272276
    H1_ p       0.000351   0.000587   0.005903   0.031267   0.000127  -0.037445
    H2_ s       0.003058   0.000734   0.162131   0.388378   0.550234   0.134600
    H2_ p       0.001552   0.000211   0.012032   0.002299   0.000227  -0.030410
 
   ao class       7ag        8ag        9ag       10ag       11ag       12ag 
    C1_ s       0.004386   0.000447   0.000562  -0.000041   0.001091   0.000367
    C1_ p       0.001713   0.000776   0.003871   0.001393   0.000166   0.000266
    C1_ d       0.000301   0.000968  -0.001207   0.000844   0.000903   0.000146
    C2_ s       0.002158   0.000927   0.000282  -0.000062   0.002139   0.000204
    C2_ p       0.004389   0.001485   0.002422   0.002828   0.000327   0.002105
    C2_ d       0.001422   0.001945  -0.000162   0.001703   0.001790   0.000898
    H1_ s       0.000129   0.002175   0.002851  -0.000010  -0.000322  -0.000147
    H1_ p      -0.000040  -0.000121   0.000189   0.000043   0.000040   0.000061
    H2_ s       0.000052   0.004333   0.001449  -0.000013  -0.000629  -0.000076
    H2_ p       0.000082  -0.000247   0.000085   0.000085   0.000078   0.000834
 
   ao class      13ag       14ag       15ag       16ag       17ag       18ag 
    C1_ s      -0.000325   0.000628  -0.000234   0.000003   0.000076  -0.000060
    C1_ p       0.001340  -0.000042  -0.000133   0.000360   0.000216   0.000013
    C1_ d       0.000684   0.000296   0.000714   0.000107   0.000205   0.000182
    C2_ s      -0.000163   0.001267  -0.000118   0.000009   0.000039  -0.000031
    C2_ p       0.000666  -0.000083   0.000413   0.000728   0.000064   0.000620
    C2_ d       0.001174   0.000586   0.000685   0.000204   0.000172  -0.000085
    H1_ s      -0.001127  -0.000287  -0.000169  -0.000072  -0.000028  -0.000000
    H1_ p       0.000337   0.000001   0.000232   0.000026   0.000051   0.000024
    H2_ s      -0.000566  -0.000572  -0.000083  -0.000151  -0.000016  -0.000000
    H2_ p       0.000245   0.000002   0.000160   0.000054   0.000150   0.000097
 
   ao class      19ag       20ag       21ag       22ag       23ag       24ag 
    C1_ s       0.000065   0.000034   0.000041  -0.000138  -0.000013  -0.000017
    C1_ p       0.000041   0.000113   0.000097   0.000026   0.000001  -0.000015
    C1_ d       0.000094  -0.000017   0.000056   0.000095   0.000060   0.000022
    C2_ s       0.000017   0.000073   0.000022  -0.000270  -0.000007  -0.000027
    C2_ p      -0.000019   0.000097  -0.000012   0.000058  -0.000079  -0.000029
    C2_ d       0.000108   0.000116   0.000110   0.000189   0.000196   0.000046
    H1_ s       0.000018  -0.000044   0.000094   0.000112  -0.000002   0.000056
    H1_ p       0.000039  -0.000017  -0.000115  -0.000014   0.000017  -0.000010
    H2_ s       0.000006  -0.000007   0.000042   0.000234  -0.000001   0.000108
    H2_ p       0.000108   0.000123  -0.000042  -0.000031  -0.000031  -0.000019
 
   ao class      25ag       26ag       27ag       28ag       29ag       30ag 
    C1_ s      -0.000175   0.000007   0.000026   0.000036   0.000001  -0.000000
    C1_ p       0.000103  -0.000008  -0.000001  -0.000023  -0.000001  -0.000016
    C1_ d       0.000009   0.000006   0.000054   0.000005   0.000000   0.000001
    C2_ s      -0.000065  -0.000039   0.000015   0.000018   0.000002   0.000000
    C2_ p       0.000170  -0.000009  -0.000013  -0.000025  -0.000003   0.000006
    C2_ d       0.000055   0.000016  -0.000008   0.000065  -0.000001   0.000000
    H1_ s       0.000016   0.000041  -0.000009   0.000009  -0.000002   0.000028
    H1_ p      -0.000003  -0.000002   0.000011   0.000003   0.000010  -0.000012
    H2_ s       0.000005   0.000091  -0.000004   0.000004  -0.000004   0.000014
    H2_ p      -0.000018  -0.000006   0.000010  -0.000042   0.000020  -0.000002
 
   ao class      31ag       32ag       33ag       34ag       35ag       36ag 
    C1_ s      -0.000001   0.000002   0.000003   0.000007   0.000004  -0.000000
    C1_ p       0.000008  -0.000007  -0.000002   0.000000  -0.000003   0.000004
    C1_ d      -0.000000   0.000000   0.000002  -0.000000  -0.000000  -0.000000
    C2_ s      -0.000001   0.000001   0.000002   0.000014   0.000002  -0.000000
    C2_ p       0.000016  -0.000005  -0.000007  -0.000000  -0.000002   0.000000
    C2_ d      -0.000001  -0.000001  -0.000010  -0.000000  -0.000000  -0.000000
    H1_ s      -0.000001  -0.000002   0.000002  -0.000006   0.000000  -0.000002
    H1_ p      -0.000000   0.000011   0.000002   0.000001  -0.000000   0.000000
    H2_ s      -0.000002  -0.000001   0.000001  -0.000012   0.000000  -0.000001
    H2_ p      -0.000001   0.000009   0.000014   0.000003  -0.000000   0.000000
 
   ao class      37ag       38ag       39ag 
    C1_ s      -0.000000  -0.000000   0.000000
    C1_ p       0.000000  -0.000000  -0.000000
    C1_ d      -0.000000  -0.000000   0.000000
    C2_ s      -0.000000  -0.000000   0.000000
    C2_ p       0.000001  -0.000000  -0.000000
    C2_ d      -0.000000  -0.000000   0.000000
    H1_ s      -0.000000   0.000000   0.000000
    H1_ p      -0.000000  -0.000000   0.000000
    H2_ s      -0.000000   0.000000   0.000000
    H2_ p      -0.000000  -0.000000   0.000000

                        b3u partial gross atomic populations
   ao class       1b3u       2b3u       3b3u       4b3u       5b3u       6b3u
    C1_ s       0.022924   1.952003   0.714956   0.079004   0.065666  -0.000085
    C1_ p      -0.000455   0.001641  -0.038320   0.214331   0.726241   0.001027
    C1_ d      -0.000102   0.000126   0.041886  -0.184282  -0.017122   0.001068
    C2_ s       1.958708   0.025923   0.363035   0.168380   0.034698  -0.000041
    C2_ p       0.013886   0.014358   0.235194   0.426492   0.554162   0.008681
    C2_ d       0.003923   0.004386  -0.025491  -0.368758  -0.040507   0.002473
    H1_ s       0.000617   0.000533   0.430874   0.551782   0.450204   0.000830
    H1_ p       0.000315   0.000354   0.039530  -0.002317  -0.016376  -0.000037
    H2_ s       0.000067   0.000410   0.217884   1.098902   0.223566   0.000414
    H2_ p       0.000116   0.000266   0.005681  -0.004878  -0.004040   0.000093
 
   ao class       7b3u       8b3u       9b3u      10b3u      11b3u      12b3u
    C1_ s       0.003474   0.000690   0.000211  -0.000952  -0.000288   0.000172
    C1_ p       0.000383   0.001897   0.002761   0.001800  -0.000073  -0.000093
    C1_ d      -0.000319   0.000110  -0.001121   0.000992   0.001319   0.000348
    C2_ s       0.006856   0.000356   0.000420  -0.000500  -0.000142   0.000090
    C2_ p       0.000725   0.002905   0.005537   0.005406   0.000493  -0.000609
    C2_ d      -0.000625   0.000353  -0.002263   0.001100   0.002144   0.001663
    H1_ s       0.000562   0.003102   0.000882  -0.002191  -0.000130   0.000048
    H1_ p       0.000226  -0.000061   0.000117  -0.000008   0.000144  -0.000013
    H2_ s       0.001095   0.001556   0.001804  -0.001087  -0.000065   0.000024
    H2_ p       0.000446  -0.000018   0.000232   0.000433   0.000026   0.000504
 
   ao class      13b3u      14b3u      15b3u      16b3u      17b3u      18b3u
    C1_ s       0.000317   0.000477  -0.000091  -0.000007   0.000149   0.000040
    C1_ p       0.001309   0.000618   0.000010  -0.000110   0.000088   0.000098
    C1_ d      -0.000495   0.000300   0.000116   0.000337   0.000136  -0.000021
    C2_ s       0.000620   0.000250  -0.000048  -0.000007   0.000054   0.000115
    C2_ p       0.002621  -0.000041   0.000533  -0.000221  -0.000167   0.000358
    C2_ d      -0.000979   0.000114   0.000204   0.000675   0.000125   0.000008
    H1_ s      -0.000825  -0.000457   0.000015   0.000039   0.000037  -0.000058
    H1_ p       0.000345   0.000113   0.000012  -0.000025  -0.000025   0.000014
    H2_ s      -0.001658  -0.000232   0.000008   0.000077   0.000013  -0.000110
    H2_ p       0.000686   0.000126   0.000036  -0.000048   0.000034  -0.000001
 
   ao class      19b3u      20b3u      21b3u      22b3u      23b3u      24b3u
    C1_ s      -0.000228  -0.000056  -0.000332  -0.000009   0.000001  -0.000064
    C1_ p       0.000043   0.000095   0.000009   0.000005   0.000003   0.000030
    C1_ d       0.000026   0.000017   0.000165   0.000002   0.000062   0.000057
    C2_ s      -0.000114  -0.000089  -0.000171  -0.000002   0.000000  -0.000134
    C2_ p       0.000447   0.000184   0.000547   0.000252  -0.000002   0.000061
    C2_ d       0.000061   0.000035  -0.000163  -0.000229   0.000028   0.000114
    H1_ s      -0.000031   0.000057   0.000096   0.000044  -0.000011   0.000017
    H1_ p       0.000015  -0.000048  -0.000019   0.000001   0.000016  -0.000017
    H2_ s      -0.000016   0.000114   0.000049   0.000022  -0.000006   0.000033
    H2_ p       0.000127  -0.000094   0.000012   0.000047   0.000007  -0.000033
 
   ao class      25b3u      26b3u      27b3u      28b3u      29b3u      30b3u
    C1_ s       0.000012   0.000009   0.000019   0.000005  -0.000007  -0.000008
    C1_ p      -0.000006  -0.000012  -0.000034  -0.000000   0.000004   0.000006
    C1_ d       0.000005  -0.000010   0.000008   0.000001  -0.000002   0.000000
    C2_ s       0.000026   0.000004   0.000009   0.000010  -0.000004  -0.000004
    C2_ p      -0.000012  -0.000026  -0.000094  -0.000001   0.000015   0.000017
    C2_ d       0.000009   0.000076   0.000087   0.000002  -0.000006  -0.000007
    H1_ s      -0.000012   0.000009   0.000034  -0.000002  -0.000003  -0.000005
    H1_ p       0.000014   0.000001   0.000003   0.000003   0.000012  -0.000001
    H2_ s      -0.000024   0.000004   0.000017  -0.000004  -0.000002  -0.000002
    H2_ p       0.000027  -0.000024  -0.000020   0.000005   0.000004   0.000009
 
   ao class      31b3u      32b3u      33b3u      34b3u      35b3u      36b3u
    C1_ s       0.000002   0.000031  -0.000001  -0.000000   0.000001  -0.000000
    C1_ p      -0.000008   0.000000   0.000000   0.000002  -0.000001   0.000000
    C1_ d       0.000003  -0.000001  -0.000001  -0.000000  -0.000000  -0.000000
    C2_ s       0.000004   0.000015  -0.000002  -0.000000   0.000002  -0.000000
    C2_ p      -0.000015  -0.000043   0.000000   0.000002  -0.000002   0.000001
    C2_ d       0.000007  -0.000001  -0.000001  -0.000001  -0.000000  -0.000000
    H1_ s       0.000006   0.000001   0.000000  -0.000001  -0.000000  -0.000000
    H1_ p      -0.000001   0.000000   0.000003   0.000000   0.000000  -0.000000
    H2_ s       0.000012   0.000001   0.000000  -0.000001  -0.000000  -0.000000
    H2_ p      -0.000003   0.000000   0.000005   0.000001   0.000000  -0.000000
 
   ao class      37b3u      38b3u      39b3u
    C1_ s      -0.000000  -0.000000   0.000000
    C1_ p      -0.000000  -0.000000  -0.000000
    C1_ d      -0.000000  -0.000000   0.000000
    C2_ s       0.000000  -0.000000   0.000000
    C2_ p       0.000000   0.000000  -0.000000
    C2_ d       0.000000  -0.000000   0.000000
    H1_ s       0.000000   0.000000   0.000000
    H1_ p       0.000000   0.000000  -0.000000
    H2_ s       0.000000   0.000000   0.000000
    H2_ p      -0.000000   0.000000  -0.000000

                        b2u partial gross atomic populations
   ao class       1b2u       2b2u       3b2u       4b2u       5b2u       6b2u
    C1_ p       0.016768   0.163056   0.123370   0.650245   0.005478   0.001304
    C1_ d       0.005055  -0.032971  -0.021455   0.018764   0.001293   0.000209
    C2_ s       1.965412   1.094704   0.100047  -0.000000  -0.000133   0.001038
    C2_ p       0.007796   0.022097   1.158723   1.313051   0.004245   0.003477
    C2_ d       0.001844   0.046768  -0.036284   0.037034   0.002265   0.000255
    H1_ p       0.000103  -0.009866   0.002600  -0.014905   0.000074   0.000007
    H2_ s       0.001851   0.646836   0.672385  -0.000000   0.001244   0.004673
    H2_ p       0.001171   0.054614  -0.022870  -0.029902  -0.000018  -0.000090
 
   ao class       7b2u       8b2u       9b2u      10b2u      11b2u      12b2u
    C1_ p       0.000814   0.003020   0.000354  -0.000376  -0.000240   0.000330
    C1_ d       0.000913   0.000409   0.000992   0.000994  -0.000023  -0.000003
    C2_ s       0.000000  -0.001486  -0.000435   0.000260   0.000735   0.000000
    C2_ p       0.001634   0.004213   0.000069  -0.000324   0.000813   0.000648
    C2_ d       0.001834   0.001682   0.002476   0.001019   0.000435  -0.000002
    H1_ p       0.000468   0.000292  -0.000031   0.000342   0.000045   0.000045
    H2_ s       0.000000  -0.003278  -0.000197   0.000071  -0.000694  -0.000000
    H2_ p       0.000936   0.000134   0.000201   0.000150   0.000195   0.000090
 
   ao class      13b2u      14b2u      15b2u      16b2u      17b2u      18b2u
    C1_ p       0.000345   0.000029  -0.000089   0.000283  -0.000020   0.000359
    C1_ d       0.000100   0.000106   0.000057   0.000032   0.000150  -0.000166
    C2_ s      -0.000133   0.000000   0.000201  -0.000342  -0.000000  -0.000497
    C2_ p       0.000191   0.000060   0.000002   0.000207  -0.000038   0.000193
    C2_ d       0.000219   0.000212   0.000210   0.000055   0.000297   0.000166
    H1_ p       0.000020   0.000089   0.000022   0.000080  -0.000059   0.000015
    H2_ s       0.000024   0.000000   0.000055  -0.000047   0.000000   0.000144
    H2_ p       0.000027   0.000181  -0.000013   0.000062  -0.000118  -0.000021
 
   ao class      19b2u      20b2u      21b2u      22b2u      23b2u      24b2u
    C1_ p       0.000170  -0.000002  -0.000014  -0.000051  -0.000002   0.000009
    C1_ d      -0.000153  -0.000002   0.000055   0.000055  -0.000007  -0.000003
    C2_ s      -0.000018   0.000001   0.000013   0.000028   0.000000  -0.000011
    C2_ p       0.000091   0.000003  -0.000024  -0.000077  -0.000003   0.000011
    C2_ d      -0.000073   0.000092   0.000012   0.000039  -0.000014  -0.000004
    H1_ p       0.000031  -0.000001  -0.000016  -0.000014   0.000015  -0.000001
    H2_ s       0.000066  -0.000017   0.000013   0.000051   0.000000  -0.000005
    H2_ p       0.000016   0.000024  -0.000007  -0.000003   0.000030   0.000017
 
   ao class      25b2u      26b2u      27b2u      28b2u      29b2u      30b2u
    C1_ p       0.000010  -0.000029   0.000001   0.000000  -0.000000   0.000000
    C1_ d      -0.000005  -0.000000  -0.000000  -0.000001   0.000000  -0.000000
    C2_ s      -0.000012   0.000046   0.000000  -0.000001  -0.000000  -0.000000
    C2_ p       0.000014  -0.000014   0.000002   0.000003   0.000000  -0.000000
    C2_ d      -0.000002  -0.000002  -0.000000  -0.000000   0.000000  -0.000000
    H1_ p       0.000006   0.000000  -0.000000   0.000001  -0.000000   0.000000
    H2_ s      -0.000007   0.000002  -0.000000  -0.000002   0.000000   0.000000
    H2_ p       0.000002   0.000001  -0.000000   0.000000   0.000000   0.000000

                        b1g partial gross atomic populations
   ao class       1b1g       2b1g       3b1g       4b1g       5b1g       6b1g
    C1_ p       0.003407   0.030606   0.765414   0.002368   0.000342   0.002736
    C1_ d       0.000348   0.001964   0.015246   0.000846   0.000299   0.000510
    C2_ s       1.989749   0.708440   0.094173   0.006520   0.000838   0.000000
    C2_ p       0.004512   0.212746   0.786610   0.003763   0.005945   0.005480
    C2_ d       0.001538  -0.169173  -0.026746   0.000862  -0.001676   0.001032
    H1_ p       0.000044  -0.008916  -0.007958   0.000067  -0.000005   0.000061
    H2_ s       0.000202   1.163185   0.406521   0.000178   0.004312   0.000000
    H2_ p       0.000200   0.042496  -0.059781  -0.000022   0.000277   0.000122
 
   ao class       7b1g       8b1g       9b1g      10b1g      11b1g      12b1g
    C1_ p       0.001315   0.000002  -0.000270   0.000317  -0.000033   0.000399
    C1_ d       0.000544   0.000554   0.000699   0.000214   0.000049  -0.000116
    C2_ s       0.000570  -0.000489   0.000000  -0.000355   0.000115  -0.000089
    C2_ p       0.001049   0.002007  -0.000541  -0.000034   0.000312   0.000231
    C2_ d       0.000493   0.001301   0.001399   0.001183   0.000328   0.000214
    H1_ p       0.000537   0.000052   0.000162   0.000031   0.000082   0.000055
    H2_ s      -0.000215  -0.001697  -0.000000  -0.000251  -0.000046  -0.000000
    H2_ p       0.000360   0.000530   0.000325   0.000361   0.000120   0.000065
 
   ao class      13b1g      14b1g      15b1g      16b1g      17b1g      18b1g
    C1_ p      -0.000007  -0.000039   0.000155  -0.000054   0.000081  -0.000008
    C1_ d       0.000056   0.000054  -0.000147   0.000113   0.000037  -0.000024
    C2_ s       0.000109   0.000063   0.000000  -0.000020  -0.000241   0.000041
    C2_ p       0.000219   0.000124   0.000310  -0.000025   0.000193  -0.000005
    C2_ d       0.000040   0.000112  -0.000295   0.000144   0.000027   0.000070
    H1_ p       0.000098   0.000009   0.000065  -0.000026  -0.000011   0.000003
    H2_ s      -0.000053   0.000135   0.000000  -0.000004   0.000021  -0.000013
    H2_ p       0.000007  -0.000167   0.000130   0.000012  -0.000010   0.000018
 
   ao class      19b1g      20b1g      21b1g      22b1g      23b1g      24b1g
    C1_ p       0.000176  -0.000007   0.000008  -0.000065  -0.000001  -0.000004
    C1_ d      -0.000178   0.000041  -0.000000   0.000092  -0.000001  -0.000007
    C2_ s       0.000000   0.000055  -0.000000  -0.000000   0.000003   0.000005
    C2_ p       0.000345  -0.000041  -0.000018  -0.000129  -0.000011  -0.000004
    C2_ d      -0.000352   0.000029   0.000002   0.000183   0.000000  -0.000001
    H1_ p       0.000024  -0.000029   0.000003  -0.000021   0.000002   0.000009
    H2_ s      -0.000000   0.000013   0.000042   0.000000  -0.000003   0.000003
    H2_ p       0.000048  -0.000011  -0.000016  -0.000041   0.000017   0.000007
 
   ao class      25b1g      26b1g      27b1g      28b1g      29b1g      30b1g
    C1_ p      -0.000000  -0.000001   0.000001   0.000000   0.000000   0.000000
    C1_ d       0.000000   0.000000  -0.000001   0.000000   0.000000   0.000000
    C2_ s       0.000007  -0.000000  -0.000000  -0.000001   0.000000   0.000000
    C2_ p      -0.000005   0.000005   0.000001   0.000001  -0.000000   0.000000
    C2_ d      -0.000000  -0.000000  -0.000002  -0.000000   0.000000   0.000000
    H1_ p      -0.000000  -0.000000   0.000001  -0.000000   0.000000  -0.000000
    H2_ s       0.000001  -0.000002  -0.000000  -0.000000   0.000000   0.000000
    H2_ p      -0.000000   0.000000   0.000002  -0.000000   0.000000  -0.000000

                        b1u partial gross atomic populations
   ao class       1b1u       2b1u       3b1u       4b1u       5b1u       6b1u
    C1_ p       0.641539   0.046704   0.000480   0.000877   0.000045   0.000880
    C1_ d       0.009120  -0.000835   0.001180   0.000600   0.000668   0.000338
    C2_ p       1.282689   0.023436   0.001106   0.002383   0.000025   0.000414
    C2_ d       0.018791   0.003595   0.002364   0.001113   0.002732   0.000706
    H1_ p       0.004391   0.000806   0.000169   0.000050   0.000250   0.000247
    H2_ p       0.008656   0.000415   0.000331   0.000103   0.000128   0.000130
 
   ao class       7b1u       8b1u       9b1u      10b1u      11b1u      12b1u
    C1_ p       0.000471   0.000008   0.000000   0.000037   0.000002   0.000001
    C1_ d       0.000069   0.000021   0.000084  -0.000008   0.000148   0.000060
    C2_ p       0.000219   0.000047   0.000009   0.000404   0.000006  -0.000001
    C2_ d       0.000181   0.000042   0.000111  -0.000016   0.000286   0.000125
    H1_ p       0.000073   0.000196   0.000166   0.000025  -0.000033  -0.000012
    H2_ p       0.000035   0.000387   0.000099   0.000022  -0.000067  -0.000006
 
   ao class      13b1u      14b1u      15b1u      16b1u
    C1_ p       0.000010  -0.000001   0.000012   0.000002
    C1_ d       0.000061  -0.000006   0.000003  -0.000002
    C2_ p       0.000000   0.000002  -0.000000   0.000000
    C2_ d       0.000058  -0.000005  -0.000000  -0.000001
    H1_ p      -0.000019   0.000019   0.000004   0.000005
    H2_ p      -0.000010   0.000033   0.000005   0.000002

                        b2g partial gross atomic populations
   ao class       1b2g       2b2g       3b2g       4b2g       5b2g       6b2g
    C1_ p       1.233502   0.006104   0.001135   0.001312   0.000146   0.000017
    C1_ d       0.001772  -0.000041   0.001852   0.000614   0.000298   0.000035
    C2_ p       0.610212   0.012261   0.000550   0.000632   0.000312   0.000008
    C2_ d       0.044207  -0.000092   0.001155   0.000270   0.000581   0.001040
    H1_ p       0.017348   0.000125   0.000447   0.000078   0.000293   0.000007
    H2_ p       0.008662   0.000258   0.000216   0.000042   0.000587   0.000003
 
   ao class       7b2g       8b2g       9b2g      10b2g      11b2g      12b2g
    C1_ p       0.000001   0.000101   0.000124   0.000032   0.000020   0.000003
    C1_ d       0.000047  -0.000014  -0.000073   0.000167   0.000010   0.000000
    C2_ p       0.000000   0.000205   0.000244   0.000019   0.000011   0.000001
    C2_ d       0.000062  -0.000028  -0.000159   0.000083   0.000012   0.000045
    H1_ p       0.000333   0.000028   0.000047  -0.000038   0.000003   0.000003
    H2_ p       0.000165   0.000058   0.000098  -0.000020   0.000002   0.000001
 
   ao class      13b2g      14b2g      15b2g      16b2g
    C1_ p      -0.000023   0.000001  -0.000000   0.000000
    C1_ d       0.000057  -0.000002  -0.000001   0.000000
    C2_ p      -0.000048   0.000001  -0.000000   0.000000
    C2_ d       0.000116  -0.000001  -0.000003   0.000000
    H1_ p      -0.000017   0.000008   0.000003  -0.000000
    H2_ p      -0.000035   0.000004   0.000006  -0.000000

                        b3g partial gross atomic populations
   ao class       1b3g       2b3g       3b3g       4b3g       5b3g       6b3g
    C1_ d       0.027597   0.001855   0.000154  -0.000033   0.000682   0.000026
    C2_ p       1.845802   0.000000   0.001623   0.001948   0.000025   0.000001
    C2_ d       0.017017   0.003691   0.002888   0.000880   0.000394   0.000083
    H2_ p       0.025686   0.000000   0.000664   0.000118   0.000010   0.000497
 
   ao class       7b3g       8b3g       9b3g      10b3g      11b3g
    C1_ d       0.000096  -0.000004   0.000009   0.000027  -0.000000
    C2_ p       0.000000   0.000050   0.000028   0.000008   0.000002
    C2_ d       0.000186   0.000255   0.000017   0.000014  -0.000003
    H2_ p      -0.000000  -0.000059   0.000004   0.000005   0.000012

                        au  partial gross atomic populations
   ao class       1au        2au        3au        4au        5au        6au 
    C1_ d       0.002765   0.001583   0.000353   0.000110   0.000198   0.000045
    C2_ p       0.068698   0.000070   0.001323   0.000676   0.000000  -0.000001
    C2_ d       0.000111   0.001814   0.000679   0.000167   0.000399   0.000155
    H2_ p       0.001232   0.000382   0.000373   0.000113   0.000000   0.000269
 
   ao class       7au        8au        9au       10au       11au 
    C1_ d       0.000062   0.000020   0.000011   0.000000  -0.000000
    C2_ p       0.000007   0.000012   0.000000   0.000013   0.000002
    C2_ d       0.000119   0.000097   0.000022   0.000006  -0.000003
    H2_ p      -0.000019  -0.000027   0.000000   0.000008   0.000008


                        gross atomic populations
     ao           C1_        C2_        H1_        H2_
      s         5.879334  11.831989   2.847819   5.681697
      p         5.413694  10.888473   0.012549   0.026842
      d        -0.196745  -0.385651   0.000000   0.000000
    total      11.096282  22.334811   2.860368   5.708539
 

 Total number of electrons:   42.00000000

