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

                CSFs      3760      437620     8722024     9868008    19031412
      internal walks     15207       34818       13467       15446       78938
valid internal walks      3760       29922        4857        5335       43874
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
 compressed index vector length=                 31479
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
 Hermit Integral Program : SIFS version  compute-0-55      14:52:49.849 24-Jun-21
  cidrt_title                                                                    
 MO-coefficients from mcscf.x                                                    
  with dummy occupation 1.0 for active orbitals                                  
  total ao core energy =  201.595839628                                          
 MCSCF energy =    -230.596267333                                                
 SIFS file created by program tran.      compute-0-55      14:53:00.701 24-Jun-21

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
 !timer: first half-sort required        cpu_time=     1.020 walltime=     1.280

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
 !timer: second half-sort required       cpu_time=     1.185 walltime=     4.400
 !timer: cisrt complete                  cpu_time=     2.212 walltime=     5.690
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
 nrow  =   148 nsym  =     8 ssym  =     3 lenbuf=  1600
 nwalk,xbar:      78938    15207    34818    13467    15446
 nvalwt,nvalw:    43874     3760    29922     4857     5335
 ncsft:        19031412
 total number of valid internal walks:   43874
 nvalz,nvaly,nvalx,nvalw =     3760   29922    4857    5335

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
 calcthrxt: niot,maxw1=                    18                 14192
 block size     0
 pthz,pthy,pthx,pthw: 15207 34818 13467 15446 total internal walks:   78938
 maxlp3,n2lp,n1lp,n0lp 14192     0     0     0
 orbsym(*)= 1 1 1 1 2 2 2 3 3 3 4 4 5 5 5 6 7 8
setref: retained number of references =    40
 setref: total/valid number of walks=                 15207
                  3760
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
    threx            145165
    twoex             36190
    onex               8593
    allin              6144
    diagon            32097
               =======
   maximum           145165
 
  __ static summary __ 
   reflst              3760
   hrfspc              3760
               -------
   static->            3760
 
  __ core required  __ 
   totstc              3760
   max n-ex          145165
               -------
   totnec->          148925
 
  __ core available __ 
   totspc        2621439999
   totnec -          148925
               -------
   totvec->      2621291074

 number of external paths / symmetry
 vertex x    2038    2119    2056    2055    1469    1467    1413    1411
 vertex w    2206    2119    2056    2055    1469    1467    1413    1411
segment: free space=  2621291074
 reducing frespc by                209386 to             2621081688 
  for index/conft/indsym storage .
 resegmenting ...



                   segmentation summary for type all-internal
 -------------------------------------------------------------------------------
 seg.      no. of|    no. of|  starting|  internal|  starting|  starting|
  no.    internal|        ci|       csf|     walks|      walk|       DRT|
            paths|  elements|    number|     /seg.|    number|    record|
 -------------------------------------------------------------------------------
  Z 1       15204|      3760|         0|      3760|         0|         1|
 -------------------------------------------------------------------------------
  Y 2       34818|    437620|      3760|     29922|      3760|         2|
 -------------------------------------------------------------------------------
  X 3       13467|   8722024|    441380|      4857|     33682|         5|
 -------------------------------------------------------------------------------
  W 4       15446|   9868008|   9163404|      5335|     38539|         6|
 -------------------------------------------------------------------------------
max. additional memory requirements:index=       19012DP  conft+indsym=      119688DP  drtbuffer=       70686 DP

dimension of the ci-matrix ->>>  19031412

 executing brd_struct for civct
 gentasklist: ntask=                    20
                    TASKLIST
----------------------------------------------------------------------------------------------------
TASK# BRA# KET#  T-TYPE    DESCR.   SEGMENTTYPE    SEGEL              SEGCI          VWALKS   
----------------------------------------------------------------------------------------------------
     1  3   1    24      two-ext xz   2X  3 1   13467   15204    8722024       3760    4857    3760
     2  4   1    25      two-ext wz   2X  4 1   15446   15204    9868008       3760    5335    3760
     3  4   3    26      two-ext wx*  WX  4 3   15446   13467    9868008    8722024    5335    4857
     4  4   3    27      two-ext wx+  WX  4 3   15446   13467    9868008    8722024    5335    4857
     5  2   1    11      one-ext yz   1X  2 1   34818   15204     437620       3760   29922    3760
     6  3   2    15      1ex3ex yx    3X  3 2   13467   34818    8722024     437620    4857   29922
     7  4   2    16      1ex3ex yw    3X  4 2   15446   34818    9868008     437620    5335   29922
     8  1   1     1      allint zz    OX  1 1   15204   15204       3760       3760    3760    3760
     9  2   2     5      0ex2ex yy    OX  2 2   34818   34818     437620     437620   29922   29922
    10  3   3     6      0ex2ex xx*   OX  3 3   13467   13467    8722024    8722024    4857    4857
    11  3   3    18      0ex2ex xx+   OX  3 3   13467   13467    8722024    8722024    4857    4857
    12  4   4     7      0ex2ex ww*   OX  4 4   15446   15446    9868008    9868008    5335    5335
    13  4   4    19      0ex2ex ww+   OX  4 4   15446   15446    9868008    9868008    5335    5335
    14  2   2    42      four-ext y   4X  2 2   34818   34818     437620     437620   29922   29922
    15  3   3    43      four-ext x   4X  3 3   13467   13467    8722024    8722024    4857    4857
    16  4   4    44      four-ext w   4X  4 4   15446   15446    9868008    9868008    5335    5335
    17  1   1    75      dg-024ext z  OX  1 1   15204   15204       3760       3760    3760    3760
    18  2   2    76      dg-024ext y  OX  2 2   34818   34818     437620     437620   29922   29922
    19  3   3    77      dg-024ext x  OX  3 3   13467   13467    8722024    8722024    4857    4857
    20  4   4    78      dg-024ext w  OX  4 4   15446   15446    9868008    9868008    5335    5335
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
 initializing v-file: 1:              19031412

    ---------trial vector generation----------

    trial vectors will be created by: 

    (ivmode= 3) diagonalizing h in the reference space.                     

      4 vectors will be written to unit 11 beginning with logical record   1

            4 vectors will be created
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:      111498 2x:           0 4x:           0
All internal counts: zz :      557690 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:           0    task #     2:           0    task #     3:           0    task #     4:           0
task #     5:           0    task #     6:           0    task #     7:           0    task #     8:      474076
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:       79366    task #    18:           0    task #    19:           0    task #    20:           0
 reference space has dimension      40
 dsyevx: computed roots 1 to    8(converged:   8)

    root           eigenvalues
    ----           ------------
       1        -230.5534777361
       2        -230.4831753627
       3        -230.4171087368
       4        -230.1503914526
       5        -230.0981520134
       6        -230.0795027436
       7        -230.0658011474
       8        -230.0626325355

 strefv generated    4 initial ci vector(s).
    ---------end of vector generation---------

 ufvoutnew: ... writing  recamt=                  3760

         vector  1 from unit 11 written to unit 49 filename cirefv              
 ufvoutnew: ... writing  recamt=                  3760

         vector  2 from unit 11 written to unit 49 filename cirefv              

 ************************************************************************
 beginning the bk-type iterative procedure (nzcsf=  3760)...
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:          19031412
 number of initial trial vectors:                         4
 number of initial matrix-vector products:                0
 maximum dimension of the subspace vectors:               7
 number of roots to converge:                             2
 number of iterations:                                    1
 residual norm convergence criteria:               0.001000  0.001000

          starting bk iteration   1

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:           0 xx:           0 ww:           0
One-external counts: yz :     2196660 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       12692 wz:       14788 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:           0    task #     4:           0
task #     5:     1940382    task #     6:           0    task #     7:           0    task #     8:      474076
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:           0 xx:           0 ww:           0
One-external counts: yz :     2196660 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       12692 wz:       14788 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:           0    task #     4:           0
task #     5:     1940382    task #     6:           0    task #     7:           0    task #     8:      474076
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:           0 xx:           0 ww:           0
One-external counts: yz :     2196660 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       12692 wz:       14788 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:           0    task #     4:           0
task #     5:     1940382    task #     6:           0    task #     7:           0    task #     8:      474076
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:           0 xx:           0 ww:           0
One-external counts: yz :     2196660 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:       12692 wz:       14788 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:           0    task #     4:           0
task #     5:     1940382    task #     6:           0    task #     7:           0    task #     8:      474076
task #     9:           0    task #    10:           0    task #    11:           0    task #    12:           0
task #    13:           0    task #    14:           0    task #    15:           0    task #    16:           0
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1    -6.10281687
   ht   2     0.00000000    -6.03251450
   ht   3     0.00000000    -0.00000000    -5.96644787
   ht   4    -0.00000000    -0.00000000    -0.00000000    -5.69973059

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4
 ref    1   -1.00000       4.275191E-13  -6.511458E-14  -3.239055E-14
 ref    2   4.271583E-13    1.00000       4.900802E-13  -8.132327E-14

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4
 ref    1    1.00000        1.00000       2.444185E-25   7.662622E-27

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1    -1.00000000     0.00000000    -0.00000000    -0.00000000
 ref:   2     0.00000000     1.00000000     0.00000000    -0.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.5534777361 -1.9540E-14  9.0379E-01  1.5765E+00  1.0000E-03   
 mr-sdci #  1  2   -230.4831753627 -5.6843E-14  0.0000E+00  1.5484E+00  1.0000E-03   
 mr-sdci #  1  3   -230.4171087368  4.2633E-14  0.0000E+00  1.5985E+00  1.0000E-04   
 mr-sdci #  1  4   -230.1503914526 -3.3751E-14  0.0000E+00  1.5516E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.207803
time for cinew                         1.900022
time for eigenvalue solver             0.000000
time for vector access                 0.000001

 mr-sdci  convergence not reached after  1 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.5534777361 -1.9540E-14  9.0379E-01  1.5765E+00  1.0000E-03   
 mr-sdci #  1  2   -230.4831753627 -5.6843E-14  0.0000E+00  1.5484E+00  1.0000E-03   
 mr-sdci #  1  3   -230.4171087368  4.2633E-14  0.0000E+00  1.5985E+00  1.0000E-04   
 mr-sdci #  1  4   -230.1503914526 -3.3751E-14  0.0000E+00  1.5516E+00  1.0000E-04   
 
    2 of the   5 expansion vectors are transformed.
    2 of the   4 matrix-vector products are transformed.

    2 expansion eigenvectors written to unit nvfile (= 11)
    2 matrix-vector products written to unit nhvfil (= 10)

 ************************************************************************
 beginning the ci iterative diagonalization procedure... 
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:          19031412
 number of initial trial vectors:                         2
 number of initial matrix-vector products:                2
 maximum dimension of the subspace vectors:               7
 number of roots to converge:                             2
 number of iterations:                                   40
 residual norm convergence criteria:               0.001000  0.001000

          starting ci iteration   1

 Final subspace hamiltonian 

                ht   1         ht   2
   ht   1    -6.10281687
   ht   2     0.00000000    -6.03251450

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2
 ref    1   -1.00000       4.271865E-13
 ref    2   4.271583E-13    1.00000    

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
 ref:   2     0.00000000     1.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.5534777361  0.0000E+00  9.0379E-01  1.5765E+00  1.0000E-03   
 mr-sdci #  1  2   -230.4831753627 -5.3291E-15  0.0000E+00  1.5484E+00  1.0000E-03   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.068747
time for cinew                         1.533045
time for eigenvalue solver             0.000092
time for vector access                 0.000002

          starting ci iteration   2

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3
   ht   1    -6.10281687
   ht   2     0.00000000    -6.03251450
   ht   3    -0.90378992    -0.00143938    -1.74953524

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3
 ref    1   0.903881       1.470832E-03   0.427781    
 ref    2  -1.296646E-03   0.999999      -6.985265E-04

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3
 ref    1   0.817003        1.00000       0.182997    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3
 ref:   1     0.90388100     0.00147083     0.42778145
 ref:   2    -0.00129665     0.99999892    -0.00069853

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  2  1   -231.1914858291  6.3801E-01  5.6712E-02  3.9019E-01  1.0000E-03   
 mr-sdci #  2  2   -230.4831755273  1.6468E-07  0.0000E+00  1.5484E+00  1.0000E-03   
 mr-sdci #  2  3   -227.7050534440 -2.7121E+00  0.0000E+00  1.4675E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.085686
time for cinew                         1.711418
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration   3

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1    -6.10281687
   ht   2     0.00000000    -6.03251450
   ht   3    -0.90378992    -0.00143938    -1.74953524
   ht   4     0.06463561     0.00134135    -0.16257698    -0.11510506

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4
 ref    1   0.905094       3.017289E-03  -0.119260      -0.408133    
 ref    2  -2.855599E-03   0.999985      -4.036802E-03   2.239682E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4
 ref    1   0.819204       0.999980       1.423925E-02   0.166577    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1     0.90509428     0.00301729    -0.11926003    -0.40813267
 ref:   2    -0.00285560     0.99998527    -0.00403680     0.00223968

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  3  1   -231.2355931411  4.4107E-02  3.1446E-03  9.7957E-02  1.0000E-03   
 mr-sdci #  3  2   -230.4832155113  3.9984E-05  0.0000E+00  1.5484E+00  1.0000E-03   
 mr-sdci #  3  3   -228.6026922517  8.9764E-01  0.0000E+00  1.3037E+00  1.0000E-04   
 mr-sdci #  3  4   -227.3654357906 -2.7850E+00  0.0000E+00  1.3672E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.092361
time for cinew                         1.890366
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration   4

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.10281687
   ht   2     0.00000000    -6.03251450
   ht   3    -0.90378992    -0.00143938    -1.74953524
   ht   4     0.06463561     0.00134135    -0.16257698    -0.11510506
   ht   5    -0.01249065     0.00224937    -0.01353850     0.00008730    -0.00459448

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.905656       4.706425E-03  -0.102337       2.177806E-04   0.411451    
 ref    2  -5.178678E-03   0.999580      -2.236071E-02   1.676813E-02  -5.605360E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.820239       0.999183       1.097291E-02   2.812175E-04   0.169324    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.90565588     0.00470643    -0.10233723     0.00021778     0.41145147
 ref:   2    -0.00517868     0.99958021    -0.02236071     0.01676813    -0.00560536

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  4  1   -231.2382542643  2.6611E-03  2.9522E-04  2.8718E-02  1.0000E-03   
 mr-sdci #  4  2   -230.4847922283  1.5767E-03  0.0000E+00  1.5463E+00  1.0000E-03   
 mr-sdci #  4  3   -228.8639674276  2.6128E-01  0.0000E+00  1.2103E+00  1.0000E-04   
 mr-sdci #  4  4   -227.8961168114  5.3068E-01  0.0000E+00  1.5059E+00  1.0000E-04   
 mr-sdci #  4  5   -227.3402983119  2.8896E+00  0.0000E+00  1.3636E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.123535
time for cinew                         1.970200
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration   5

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.10281687
   ht   2     0.00000000    -6.03251450
   ht   3    -0.90378992    -0.00143938    -1.74953524
   ht   4     0.06463561     0.00134135    -0.16257698    -0.11510506
   ht   5    -0.01249065     0.00224937    -0.01353850     0.00008730    -0.00459448
   ht   6     0.00580871     0.00176883     0.00219462    -0.00093772     0.00006420    -0.00053771

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.905469       8.392131E-03  -4.818890E-02   0.147736       9.042045E-02  -0.384359    
 ref    2  -7.144066E-03   0.995462      -7.908414E-02  -5.207758E-02  -2.321196E-03  -5.742817E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.819925       0.991015       8.576472E-03   2.453788E-02   8.181246E-03   0.147765    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.90546897     0.00839213    -0.04818890     0.14773561     0.09042045    -0.38435878
 ref:   2    -0.00714407     0.99546181    -0.07908414    -0.05207758    -0.00232120    -0.00574282

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  5  1   -231.2384927565  2.3849E-04  3.7548E-05  1.0420E-02  1.0000E-03   
 mr-sdci #  5  2   -230.4980384716  1.3246E-02  0.0000E+00  1.5269E+00  1.0000E-03   
 mr-sdci #  5  3   -229.0809534697  2.1699E-01  0.0000E+00  1.2567E+00  1.0000E-04   
 mr-sdci #  5  4   -228.3164168844  4.2030E-01  0.0000E+00  1.3198E+00  1.0000E-04   
 mr-sdci #  5  5   -227.8263343848  4.8604E-01  0.0000E+00  1.7491E+00  1.0000E-04   
 mr-sdci #  5  6   -227.2564171267  2.8058E+00  0.0000E+00  1.3753E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.143555
time for cinew                         2.114349
time for eigenvalue solver             0.000153
time for vector access                 0.000000

          starting ci iteration   6

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.10281687
   ht   2     0.00000000    -6.03251450
   ht   3    -0.90378992    -0.00143938    -1.74953524
   ht   4     0.06463561     0.00134135    -0.16257698    -0.11510506
   ht   5    -0.01249065     0.00224937    -0.01353850     0.00008730    -0.00459448
   ht   6     0.00580871     0.00176883     0.00219462    -0.00093772     0.00006420    -0.00053771
   ht   7    -0.00060833     0.00145023    -0.00042930    -0.00060855    -0.00006715    -0.00001419    -0.00006642

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.905375      -1.506757E-02   7.737163E-03  -0.166361       3.828385E-02   3.476236E-02  -0.386859    
 ref    2   9.303016E-03  -0.958217       0.261108       9.352914E-02  -4.286763E-02   5.105551E-02  -1.910371E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.819790       0.918406       6.823744E-02   3.642377E-02   3.303287E-03   3.815087E-03   0.150025    

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.90537462    -0.01506757     0.00773716    -0.16636128     0.03828385     0.03476236    -0.38685883
 ref:   2     0.00930302    -0.95821654     0.26110836     0.09352914    -0.04286763     0.05105551    -0.01910371

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  6  1   -231.2385289246  3.6168E-05  8.8555E-06  5.3349E-03  1.0000E-03   
 mr-sdci #  6  2   -230.5959575832  9.7919E-02  0.0000E+00  1.3610E+00  1.0000E-03   
 mr-sdci #  6  3   -229.4175230794  3.3657E-01  0.0000E+00  1.3076E+00  1.0000E-04   
 mr-sdci #  6  4   -228.4432972787  1.2688E-01  0.0000E+00  1.2166E+00  1.0000E-04   
 mr-sdci #  6  5   -227.9799186900  1.5358E-01  0.0000E+00  1.2146E+00  1.0000E-04   
 mr-sdci #  6  6   -227.6679629534  4.1155E-01  0.0000E+00  2.0132E+00  1.0000E-04   
 mr-sdci #  6  7   -227.2404967838  2.7898E+00  0.0000E+00  1.3602E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.162079
time for cinew                         3.049896
time for eigenvalue solver             0.000153
time for vector access                 0.000000

          starting ci iteration   7

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.78786806
   ht   2     0.00000000    -6.14529672
   ht   3     0.00000000     0.00000000    -4.96686222
   ht   4    -0.00000000     0.00000000     0.00000000    -3.99263641
   ht   5    -0.00037118     0.00157658     0.00080001     0.00029676    -0.00001315

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.905292       2.230904E-02   3.913904E-03   0.159841      -4.548765E-02
 ref    2  -1.192070E-02   0.905374       0.342721      -0.164158      -0.176336    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.819696       0.820199       0.117473       5.249718E-02   3.316348E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.90529205     0.02230904     0.00391390     0.15984115    -0.04548765
 ref:   2    -0.01192070     0.90537360     0.34272063    -0.16415843    -0.17633592

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  7  1   -231.2385395810  1.0656E-05  6.0395E-06  4.1774E-03  1.0000E-03   
 mr-sdci #  7  2   -230.8343195154  2.3836E-01  0.0000E+00  9.8308E-01  1.0000E-03   
 mr-sdci #  7  3   -229.4476011457  3.0078E-02  0.0000E+00  1.3696E+00  1.0000E-04   
 mr-sdci #  7  4   -228.4981684016  5.4871E-02  0.0000E+00  1.3286E+00  1.0000E-04   
 mr-sdci #  7  5   -227.8068481673 -1.7307E-01  0.0000E+00  1.6813E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.124451
time for cinew                         1.977509
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration   8

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.78786806
   ht   2     0.00000000    -6.14529672
   ht   3     0.00000000     0.00000000    -4.96686222
   ht   4    -0.00000000     0.00000000     0.00000000    -3.99263641
   ht   5    -0.00037118     0.00157658     0.00080001     0.00029676    -0.00001315
   ht   6     0.00026088    -0.00037188     0.00030437    -0.00044338    -0.00000096    -0.00001020

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.905221      -2.178611E-02  -1.276775E-03   0.128942      -0.108498      -3.175093E-03
 ref    2  -1.596840E-02  -0.873565      -0.203021      -0.330803      -0.286824      -9.190567E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.819680       0.763591       4.121899E-02   0.126056       9.403988E-02   9.454773E-05

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.90522097    -0.02178611    -0.00127678     0.12894184    -0.10849799    -0.00317509
 ref:   2    -0.01596840    -0.87356510    -0.20302059    -0.33080279    -0.28682410    -0.00919057

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  8  1   -231.2385487276  9.1466E-06  3.7149E-06  3.2718E-03  1.0000E-03   
 mr-sdci #  8  2   -231.0383257318  2.0401E-01  0.0000E+00  5.1121E-01  1.0000E-03   
 mr-sdci #  8  3   -229.5235137725  7.5913E-02  0.0000E+00  1.1894E+00  1.0000E-04   
 mr-sdci #  8  4   -228.5740650627  7.5897E-02  0.0000E+00  1.5885E+00  1.0000E-04   
 mr-sdci #  8  5   -228.1447578888  3.3791E-01  0.0000E+00  1.4800E+00  1.0000E-04   
 mr-sdci #  8  6   -227.6437803924 -2.4183E-02  0.0000E+00  1.8739E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.149353
time for cinew                         2.109314
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration   9

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.78786806
   ht   2     0.00000000    -6.14529672
   ht   3     0.00000000     0.00000000    -4.96686222
   ht   4    -0.00000000     0.00000000     0.00000000    -3.99263641
   ht   5    -0.00037118     0.00157658     0.00080001     0.00029676    -0.00001315
   ht   6     0.00026088    -0.00037188     0.00030437    -0.00044338    -0.00000096    -0.00001020
   ht   7     0.00000563    -0.00068686     0.00158657    -0.00015906     0.00000177     0.00000088    -0.00000749

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.905186       2.069792E-02  -1.101877E-02   0.130007       8.910277E-02  -6.208083E-02  -5.870981E-03
 ref    2  -1.871766E-02   0.874624       9.375155E-02  -0.324455       0.331405      -8.780440E-02   1.090202E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.819711       0.765395       8.910767E-03   0.122173       0.117769       1.156364E-02   1.533224E-04

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.90518556     0.02069792    -0.01101877     0.13000675     0.08910277    -0.06208083    -0.00587098
 ref:   2    -0.01871766     0.87462373     0.09375155    -0.32445543     0.33140519    -0.08780440     0.01090202

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  9  1   -231.2385529881  4.2605E-06  1.3693E-06  1.9051E-03  1.0000E-03   
 mr-sdci #  9  2   -231.1017478650  6.3422E-02  0.0000E+00  2.6201E-01  1.0000E-03   
 mr-sdci #  9  3   -229.7382433636  2.1473E-01  0.0000E+00  1.0058E+00  1.0000E-04   
 mr-sdci #  9  4   -228.5742511628  1.8610E-04  0.0000E+00  1.5760E+00  1.0000E-04   
 mr-sdci #  9  5   -228.2673357997  1.2258E-01  0.0000E+00  1.6593E+00  1.0000E-04   
 mr-sdci #  9  6   -228.0129331118  3.6915E-01  0.0000E+00  1.8041E+00  1.0000E-04   
 mr-sdci #  9  7   -227.5871099779  3.4661E-01  0.0000E+00  1.6370E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.202881
time for cinew                         3.040894
time for eigenvalue solver             0.000183
time for vector access                 0.000000

          starting ci iteration  10

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.78789212
   ht   2    -0.00000000    -6.65108700
   ht   3     0.00000000    -0.00000000    -5.28758250
   ht   4    -0.00000000     0.00000000     0.00000000    -4.12359030
   ht   5     0.00098653    -0.00018017     0.00104007    -0.00048825    -0.00000282

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.905125       2.585717E-02   1.477432E-03  -6.841748E-02  -0.138384    
 ref    2   1.981584E-02   0.885873      -0.101786       0.237802       0.179282    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.819644       0.785439       1.036267E-02   6.123067E-02   5.129243E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.90512511     0.02585717     0.00147743    -0.06841748    -0.13838435
 ref:   2     0.01981584     0.88587273    -0.10178650     0.23780184     0.17928247

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 10  1   -231.2385541758  1.1877E-06  4.4538E-07  1.1588E-03  1.0000E-03   
 mr-sdci # 10  2   -231.1198478836  1.8100E-02  0.0000E+00  1.9594E-01  1.0000E-03   
 mr-sdci # 10  3   -229.9336995975  1.9546E-01  0.0000E+00  9.4635E-01  1.0000E-04   
 mr-sdci # 10  4   -228.8835352498  3.0928E-01  0.0000E+00  1.6036E+00  1.0000E-04   
 mr-sdci # 10  5   -227.6878098408 -5.7953E-01  0.0000E+00  1.5588E+00  1.0000E-04   
 
 root number  1 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.114502
time for cinew                         1.992798
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  11

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.78789212
   ht   2    -0.00000000    -6.65108700
   ht   3     0.00000000    -0.00000000    -5.28758250
   ht   4    -0.00000000     0.00000000     0.00000000    -4.12359030
   ht   5     0.00098653    -0.00018017     0.00104007    -0.00048825    -0.00000282
   ht   6     0.00011116     0.00032772    -0.00035962     0.00038425     0.00000001    -0.00000095

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.905124      -2.388209E-02   2.327186E-02  -5.495678E-02   4.302037E-02  -0.136493    
 ref    2   2.041583E-02  -0.885923      -4.914872E-02   5.524587E-02  -0.268920       0.164058    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.819667       0.785429       2.957176E-03   6.072354E-03   7.416888E-02   4.554536E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.90512443    -0.02388209     0.02327186    -0.05495678     0.04302037    -0.13649341
 ref:   2     0.02041583    -0.88592257    -0.04914872     0.05524587    -0.26892030     0.16405764

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 11  1   -231.2385546918  5.1603E-07  0.0000E+00  7.1759E-04  1.0000E-03   
 mr-sdci # 11  2   -231.1312673333  1.1419E-02  7.4832E-03  1.5023E-01  1.0000E-03   
 mr-sdci # 11  3   -230.2058544109  2.7215E-01  0.0000E+00  8.9382E-01  1.0000E-04   
 mr-sdci # 11  4   -229.3644327942  4.8090E-01  0.0000E+00  1.2400E+00  1.0000E-04   
 mr-sdci # 11  5   -228.0037680341  3.1596E-01  0.0000E+00  1.8467E+00  1.0000E-04   
 mr-sdci # 11  6   -227.6864331155 -3.2650E-01  0.0000E+00  1.6029E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.133118
time for cinew                         4.166687
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  12

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.78789212
   ht   2    -0.00000000    -6.65108700
   ht   3     0.00000000    -0.00000000    -5.28758250
   ht   4    -0.00000000     0.00000000     0.00000000    -4.12359030
   ht   5     0.00098653    -0.00018017     0.00104007    -0.00048825    -0.00000282
   ht   6     0.00011116     0.00032772    -0.00035962     0.00038425     0.00000001    -0.00000095
   ht   7    -0.00840056     0.08404373    -0.01295352    -0.06960938    -0.00006253    -0.00000341    -0.01780469

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.905124      -2.148920E-02   3.400315E-02  -2.880082E-02  -5.112632E-02  -5.821865E-02  -0.127138    
 ref    2   2.094182E-02  -0.898185      -3.701111E-02   4.468221E-02  -5.736132E-02   0.266756       0.161968    

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.819688       0.807198       2.526037E-03   2.825987E-03   5.904221E-03   7.454794E-02   4.239774E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.90512376    -0.02148920     0.03400315    -0.02880082    -0.05112632    -0.05821865    -0.12713798
 ref:   2     0.02094182    -0.89818514    -0.03701111     0.04468221    -0.05736132     0.26675556     0.16196814

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 12  1   -231.2385549334  2.4154E-07  0.0000E+00  5.1775E-04  1.0000E-03   
 mr-sdci # 12  2   -231.1432074902  1.1940E-02  5.7528E-03  1.1307E-01  1.0000E-03   
 mr-sdci # 12  3   -230.5133983153  3.0754E-01  0.0000E+00  8.1908E-01  1.0000E-04   
 mr-sdci # 12  4   -229.7313208136  3.6689E-01  0.0000E+00  1.0632E+00  1.0000E-04   
 mr-sdci # 12  5   -228.3806905644  3.7692E-01  0.0000E+00  1.6669E+00  1.0000E-04   
 mr-sdci # 12  6   -227.9895847426  3.0315E-01  0.0000E+00  1.7261E+00  1.0000E-04   
 mr-sdci # 12  7   -227.6587773531  7.1667E-02  0.0000E+00  1.6089E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.158020
time for cinew                         3.030640
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  13

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.78789407
   ht   2     0.00000000    -6.69254663
   ht   3     0.00000000    -0.00000000    -6.06273745
   ht   4    -0.00000000     0.00000000    -0.00000000    -5.28065995
   ht   5    -0.04385866    -0.03361455    -0.12225575     0.06587300    -0.01641780

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.905114      -2.264841E-02  -1.914730E-02   5.452010E-04  -7.745768E-02
 ref    2  -2.115251E-02  -0.902628       7.512286E-03  -3.376918E-02  -2.627167E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.819679       0.815250       4.230537E-04   1.140655E-03   6.689893E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.90511431    -0.02264841    -0.01914730     0.00054520    -0.07745768
 ref:   2    -0.02115251    -0.90262788     0.00751229    -0.03376918    -0.02627167

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 13  1   -231.2385550134  8.0067E-08  0.0000E+00  4.1928E-04  1.0000E-03   
 mr-sdci # 13  2   -231.1482028845  4.9954E-03  1.7193E-03  6.9651E-02  1.0000E-03   
 mr-sdci # 13  3   -230.7076831821  1.9428E-01  0.0000E+00  4.6886E-01  1.0000E-04   
 mr-sdci # 13  4   -229.9102326625  1.7891E-01  0.0000E+00  8.4965E-01  1.0000E-04   
 mr-sdci # 13  5   -227.9291618575 -4.5153E-01  0.0000E+00  1.8203E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.119263
time for cinew                         1.978821
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  14

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.78789407
   ht   2     0.00000000    -6.69254663
   ht   3     0.00000000    -0.00000000    -6.06273745
   ht   4    -0.00000000     0.00000000    -0.00000000    -5.28065995
   ht   5    -0.04385866    -0.03361455    -0.12225575     0.06587300    -0.01641780
   ht   6     0.02930309     0.03356298     0.00714497     0.02355962    -0.00005644    -0.00366908

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.905116      -2.169594E-02  -2.215548E-02   5.366704E-04   1.503321E-02   8.142516E-02
 ref    2  -2.123715E-02  -0.901649      -2.266218E-02  -1.215568E-02  -3.943607E-02   4.538770E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.819686       0.813441       1.004440E-03   1.480485E-04   1.781201E-03   8.690099E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.90511620    -0.02169594    -0.02215548     0.00053667     0.01503321     0.08142516
 ref:   2    -0.02123715    -0.90164870    -0.02266218    -0.01215568    -0.03943607     0.04538770

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 14  1   -231.2385550515  3.8116E-08  0.0000E+00  2.7803E-04  1.0000E-03   
 mr-sdci # 14  2   -231.1502269292  2.0240E-03  5.7747E-04  3.4658E-02  1.0000E-03   
 mr-sdci # 14  3   -230.7754783628  6.7795E-02  0.0000E+00  2.8420E-01  1.0000E-04   
 mr-sdci # 14  4   -230.0186056697  1.0837E-01  0.0000E+00  8.2795E-01  1.0000E-04   
 mr-sdci # 14  5   -228.9378540500  1.0087E+00  0.0000E+00  1.3270E+00  1.0000E-04   
 mr-sdci # 14  6   -227.6013387500 -3.8825E-01  0.0000E+00  1.8926E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.159302
time for cinew                         2.086670
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  15

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.78789407
   ht   2     0.00000000    -6.69254663
   ht   3     0.00000000    -0.00000000    -6.06273745
   ht   4    -0.00000000     0.00000000    -0.00000000    -5.28065995
   ht   5    -0.04385866    -0.03361455    -0.12225575     0.06587300    -0.01641780
   ht   6     0.02930309     0.03356298     0.00714497     0.02355962    -0.00005644    -0.00366908
   ht   7     0.00989891    -0.02565954     0.01240025     0.00015484    -0.00136180     0.00004858    -0.00208287

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.905118       2.112927E-02   2.496352E-02   5.237337E-03  -2.021101E-03  -9.249152E-03   8.956306E-02
 ref    2  -2.125224E-02   0.901270       3.162576E-02   1.831251E-03  -3.345266E-02   3.584903E-02   3.941804E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.819690       0.812734       1.623366E-03   3.078318E-05   1.123165E-03   1.370699E-03   9.575323E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.90511760     0.02112927     0.02496352     0.00523734    -0.00202110    -0.00924915     0.08956306
 ref:   2    -0.02125224     0.90127013     0.03162576     0.00183125    -0.03345266     0.03584903     0.03941804

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 15  1   -231.2385550568  5.2431E-09  0.0000E+00  2.6287E-04  1.0000E-03   
 mr-sdci # 15  2   -231.1507031037  4.7617E-04  1.8202E-04  2.3865E-02  1.0000E-03   
 mr-sdci # 15  3   -230.7960949320  2.0617E-02  0.0000E+00  2.9779E-01  1.0000E-04   
 mr-sdci # 15  4   -230.1567077712  1.3810E-01  0.0000E+00  8.2796E-01  1.0000E-04   
 mr-sdci # 15  5   -229.5452550261  6.0740E-01  0.0000E+00  1.3275E+00  1.0000E-04   
 mr-sdci # 15  6   -228.3290870872  7.2775E-01  0.0000E+00  1.9197E+00  1.0000E-04   
 mr-sdci # 15  7   -227.4992816205 -1.5950E-01  0.0000E+00  1.7653E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.154907
time for cinew                         3.038330
time for eigenvalue solver             0.000244
time for vector access                 0.000000

          starting ci iteration  16

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.78789419
   ht   2    -0.00000000    -6.70004224
   ht   3     0.00000000    -0.00000000    -6.34543407
   ht   4     0.00000000     0.00000000     0.00000000    -5.70604691
   ht   5    -0.00825694    -0.00085932     0.00818456     0.00611796    -0.00035906

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.905116       2.208489E-02  -1.355851E-02   3.187823E-02   9.401174E-02
 ref    2   2.125444E-02   0.900374      -4.469048E-02  -2.377640E-02  -5.305133E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.819687       0.811160       2.181072E-03   1.581539E-03   1.165265E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.90511629     0.02208489    -0.01355851     0.03187823     0.09401174
 ref:   2     0.02125444     0.90037352    -0.04469048    -0.02377640    -0.05305133

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 16  1   -231.2385550572  4.1375E-10  0.0000E+00  2.6515E-04  1.0000E-03   
 mr-sdci # 16  2   -231.1509482751  2.4517E-04  1.8106E-04  1.7237E-02  1.0000E-03   
 mr-sdci # 16  3   -230.8173537944  2.1259E-02  0.0000E+00  2.9681E-01  1.0000E-04   
 mr-sdci # 16  4   -230.3606767621  2.0397E-01  0.0000E+00  7.1788E-01  1.0000E-04   
 mr-sdci # 16  5   -228.7277429943 -8.1751E-01  0.0000E+00  1.8751E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.119507
time for cinew                         1.984131
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  17

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.78789419
   ht   2    -0.00000000    -6.70004224
   ht   3     0.00000000    -0.00000000    -6.34543407
   ht   4     0.00000000     0.00000000     0.00000000    -5.70604691
   ht   5    -0.00825694    -0.00085932     0.00818456     0.00611796    -0.00035906
   ht   6     0.01432037    -0.00660852     0.00278029     0.00912766    -0.00019789    -0.00084639

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.905117       2.182972E-02  -7.977314E-03   2.272030E-02   1.018452E-02  -0.123443    
 ref    2   2.125807E-02   0.897976      -9.218925E-02  -3.210687E-02   2.741083E-02   2.147808E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.819689       0.806837       8.562496E-03   1.547063E-03   8.550784E-04   1.569956E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.90511699     0.02182972    -0.00797731     0.02272030     0.01018452    -0.12344333
 ref:   2     0.02125807     0.89797589    -0.09218925    -0.03210687     0.02741083     0.02147808

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 17  1   -231.2385550607  3.5340E-09  0.0000E+00  2.8663E-04  1.0000E-03   
 mr-sdci # 17  2   -231.1512448632  2.9659E-04  1.3498E-04  2.0443E-02  1.0000E-03   
 mr-sdci # 17  3   -230.9065201846  8.9166E-02  0.0000E+00  4.7413E-01  1.0000E-04   
 mr-sdci # 17  4   -230.6773780075  3.1670E-01  0.0000E+00  5.0507E-01  1.0000E-04   
 mr-sdci # 17  5   -230.0362917088  1.3085E+00  0.0000E+00  8.7075E-01  1.0000E-04   
 mr-sdci # 17  6   -227.5870034674 -7.4208E-01  0.0000E+00  1.7283E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.142944
time for cinew                         2.104736
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  18

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.78789419
   ht   2    -0.00000000    -6.70004224
   ht   3     0.00000000    -0.00000000    -6.34543407
   ht   4     0.00000000     0.00000000     0.00000000    -5.70604691
   ht   5    -0.00825694    -0.00085932     0.00818456     0.00611796    -0.00035906
   ht   6     0.01432037    -0.00660852     0.00278029     0.00912766    -0.00019789    -0.00084639
   ht   7    -0.00026087     0.01045045    -0.00313518    -0.00452998    -0.00005002    -0.00008901    -0.00023909

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.905116       2.195559E-02  -1.852263E-03   2.386606E-02   1.086705E-02  -3.824225E-02   0.118952    
 ref    2   2.125555E-02   0.894989      -0.130633      -2.286888E-02   3.112081E-02  -4.175097E-02  -4.724599E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.819686       0.801487       1.706841E-02   1.092575E-03   1.086597E-03   3.205613E-03   1.638172E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.90511578     0.02195559    -0.00185226     0.02386606     0.01086705    -0.03824225     0.11895183
 ref:   2     0.02125555     0.89498879    -0.13063299    -0.02286888     0.03112081    -0.04175097    -0.04724599

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 18  1   -231.2385550715  1.0720E-08  0.0000E+00  2.6797E-04  1.0000E-03   
 mr-sdci # 18  2   -231.1514342545  1.8939E-04  1.1800E-04  1.5831E-02  1.0000E-03   
 mr-sdci # 18  3   -230.9843052964  7.7785E-02  0.0000E+00  3.3107E-01  1.0000E-04   
 mr-sdci # 18  4   -230.7253279747  4.7950E-02  0.0000E+00  3.1873E-01  1.0000E-04   
 mr-sdci # 18  5   -230.1512456245  1.1495E-01  0.0000E+00  7.4227E-01  1.0000E-04   
 mr-sdci # 18  6   -228.0757357430  4.8873E-01  0.0000E+00  1.6611E+00  1.0000E-04   
 mr-sdci # 18  7   -227.4617296486 -3.7552E-02  0.0000E+00  1.7691E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.167114
time for cinew                         3.044800
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  19

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.78789421
   ht   2    -0.00000000    -6.70077339
   ht   3     0.00000000     0.00000000    -6.53364443
   ht   4    -0.00000000    -0.00000000    -0.00000000    -6.27466711
   ht   5     0.00443555    -0.00619440    -0.01870210     0.00036325    -0.00044604

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.905113       2.233111E-02  -5.157839E-03   2.685170E-02   2.591293E-02
 ref    2  -2.125286E-02   0.892345       0.150842      -1.256569E-02  -8.834296E-03

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.819681       0.796778       2.277996E-02   8.789105E-04   7.495246E-04

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.90511298     0.02233111    -0.00515784     0.02685170     0.02591293
 ref:   2    -0.02125286     0.89234464     0.15084216    -0.01256569    -0.00883430

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 19  1   -231.2385550794  7.9253E-09  0.0000E+00  2.7166E-04  1.0000E-03   
 mr-sdci # 19  2   -231.1515661861  1.3193E-04  1.1155E-04  1.7696E-02  1.0000E-03   
 mr-sdci # 19  3   -231.0299487225  4.5643E-02  0.0000E+00  3.1074E-01  1.0000E-04   
 mr-sdci # 19  4   -230.7459484763  2.0621E-02  0.0000E+00  3.0310E-01  1.0000E-04   
 mr-sdci # 19  5   -229.2315971921 -9.1965E-01  0.0000E+00  1.9865E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.118164
time for cinew                         1.977295
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  20

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.78789421
   ht   2    -0.00000000    -6.70077339
   ht   3     0.00000000     0.00000000    -6.53364443
   ht   4    -0.00000000    -0.00000000    -0.00000000    -6.27466711
   ht   5     0.00443555    -0.00619440    -0.01870210     0.00036325    -0.00044604
   ht   6    -0.00206088     0.00628851     0.00409005    -0.00333668    -0.00007297    -0.00023041

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.905115       2.203920E-02   6.840241E-04  -2.000851E-02   2.086121E-02  -6.027807E-02
 ref    2  -2.124076E-02   0.882463       0.209662       8.379187E-03   2.764814E-02  -5.353116E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.819684       0.779227       4.395883E-02   4.705514E-04   1.199610E-03   6.499031E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.90511461     0.02203920     0.00068402    -0.02000851     0.02086121    -0.06027807
 ref:   2    -0.02124076     0.88246319     0.20966249     0.00837919     0.02764814    -0.05353116

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 20  1   -231.2385550898  1.0414E-08  0.0000E+00  2.9189E-04  1.0000E-03   
 mr-sdci # 20  2   -231.1518860837  3.1990E-04  2.5961E-04  2.5990E-02  1.0000E-03   
 mr-sdci # 20  3   -231.0871052006  5.7156E-02  0.0000E+00  2.7926E-01  1.0000E-04   
 mr-sdci # 20  4   -230.7726148492  2.6666E-02  0.0000E+00  3.2024E-01  1.0000E-04   
 mr-sdci # 20  5   -230.3794884474  1.1479E+00  0.0000E+00  9.3255E-01  1.0000E-04   
 mr-sdci # 20  6   -227.4385929489 -6.3714E-01  0.0000E+00  1.8145E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.131958
time for cinew                         2.126709
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  21

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.78789421
   ht   2    -0.00000000    -6.70077339
   ht   3     0.00000000     0.00000000    -6.53364443
   ht   4    -0.00000000    -0.00000000    -0.00000000    -6.27466711
   ht   5     0.00443555    -0.00619440    -0.01870210     0.00036325    -0.00044604
   ht   6    -0.00206088     0.00628851     0.00409005    -0.00333668    -0.00007297    -0.00023041
   ht   7     0.00196021    -0.00869933    -0.02083015    -0.01532775    -0.00030023    -0.00005004    -0.00067104

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.905116       2.150084E-02   5.115692E-03  -1.279669E-02  -2.901157E-02   7.907822E-03  -6.185093E-02
 ref    2  -2.123215E-02   0.868028       0.263934       6.000043E-04  -7.174766E-03   3.629823E-02  -4.356982E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.819686       0.753935       6.968739E-02   1.641154E-04   8.931483E-04   1.380095E-03   5.723867E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.90511630     0.02150084     0.00511569    -0.01279669    -0.02901157     0.00790782    -0.06185093
 ref:   2    -0.02123215     0.86802804     0.26393412     0.00060000    -0.00717477     0.03629823    -0.04356982

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 21  1   -231.2385550993  9.4852E-09  0.0000E+00  2.7315E-04  1.0000E-03   
 mr-sdci # 21  2   -231.1522443655  3.5828E-04  1.4033E-04  2.1096E-02  1.0000E-03   
 mr-sdci # 21  3   -231.1129008034  2.5796E-02  0.0000E+00  1.4667E-01  1.0000E-04   
 mr-sdci # 21  4   -230.7940105671  2.1396E-02  0.0000E+00  2.7155E-01  1.0000E-04   
 mr-sdci # 21  5   -230.5836887697  2.0420E-01  0.0000E+00  4.6902E-01  1.0000E-04   
 mr-sdci # 21  6   -227.6232728251  1.8468E-01  0.0000E+00  1.7220E+00  1.0000E-04   
 mr-sdci # 21  7   -227.4095194930 -5.2210E-02  0.0000E+00  1.7489E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.177246
time for cinew                         3.048462
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  22

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.78789424
   ht   2    -0.00000000    -6.70158350
   ht   3    -0.00000000     0.00000000    -6.66223994
   ht   4     0.00000000    -0.00000000     0.00000000    -6.34334970
   ht   5    -0.00827360     0.00387803    -0.00314603    -0.00218987    -0.00024867

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.905113       2.195568E-02  -2.188712E-03   1.614545E-02   6.437606E-02
 ref    2   2.122586E-02   0.860076      -0.289638       5.677583E-04  -3.152802E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.819680       0.740212       8.389483E-02   2.609980E-04   5.138293E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.90511271     0.02195568    -0.00218871     0.01614545     0.06437606
 ref:   2     0.02122586     0.86007573    -0.28963777     0.00056776    -0.03152802

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 22  1   -231.2385551064  7.0841E-09  0.0000E+00  2.3708E-04  1.0000E-03   
 mr-sdci # 22  2   -231.1524029774  1.5861E-04  6.6473E-05  1.2798E-02  1.0000E-03   
 mr-sdci # 22  3   -231.1192311433  6.3303E-03  0.0000E+00  7.7460E-02  1.0000E-04   
 mr-sdci # 22  4   -230.8000007693  5.9902E-03  0.0000E+00  2.5053E-01  1.0000E-04   
 mr-sdci # 22  5   -228.5269991002 -2.0567E+00  0.0000E+00  1.7861E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.111694
time for cinew                         1.990845
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  23

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.78789424
   ht   2    -0.00000000    -6.70158350
   ht   3    -0.00000000     0.00000000    -6.66223994
   ht   4     0.00000000    -0.00000000     0.00000000    -6.34334970
   ht   5    -0.00827360     0.00387803    -0.00314603    -0.00218987    -0.00024867
   ht   6     0.00375367    -0.00708555     0.01043565     0.00422292    -0.00001577    -0.00016894

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1  -0.905114       2.157867E-02  -4.359952E-03   1.367676E-02  -2.176950E-02  -0.109746    
 ref    2   2.122510E-02   0.856865      -0.298682       2.508489E-03   5.226065E-03   4.614967E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.819683       0.734683       8.922978E-02   1.933463E-04   5.012227E-04   1.417392E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1    -0.90511448     0.02157867    -0.00435995     0.01367676    -0.02176950    -0.10974576
 ref:   2     0.02122510     0.85686510    -0.29868173     0.00250849     0.00522606     0.04614967

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 23  1   -231.2385551081  1.7740E-09  0.0000E+00  2.2551E-04  1.0000E-03   
 mr-sdci # 23  2   -231.1524673381  6.4361E-05  1.5133E-05  6.6519E-03  1.0000E-03   
 mr-sdci # 23  3   -231.1212889605  2.0578E-03  0.0000E+00  4.3934E-02  1.0000E-04   
 mr-sdci # 23  4   -230.8027425440  2.7418E-03  0.0000E+00  2.4512E-01  1.0000E-04   
 mr-sdci # 23  5   -229.6357358879  1.1087E+00  0.0000E+00  1.2935E+00  1.0000E-04   
 mr-sdci # 23  6   -227.5054307437 -1.1784E-01  0.0000E+00  1.7158E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.151611
time for cinew                         2.114258
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  24

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.78789424
   ht   2    -0.00000000    -6.70158350
   ht   3    -0.00000000     0.00000000    -6.66223994
   ht   4     0.00000000    -0.00000000     0.00000000    -6.34334970
   ht   5    -0.00827360     0.00387803    -0.00314603    -0.00218987    -0.00024867
   ht   6     0.00375367    -0.00708555     0.01043565     0.00422292    -0.00001577    -0.00016894
   ht   7    -0.00135626     0.00305624    -0.00419814     0.00055003    -0.00003751    -0.00000192    -0.00003590

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1  -0.905115       2.153312E-02  -4.609220E-03  -1.297781E-02  -1.864756E-02   7.288350E-03   0.110191    
 ref    2   2.122454E-02   0.855723      -0.302342      -1.244039E-03  -1.931552E-02  -4.350928E-02  -4.870669E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.819683       0.732726       9.143181E-02   1.699711E-04   7.208208E-04   1.946178E-03   1.451450E-02

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1    -0.90511455     0.02153312    -0.00460922    -0.01297781    -0.01864756     0.00728835     0.11019148
 ref:   2     0.02122454     0.85572303    -0.30234180    -0.00124404    -0.01931552    -0.04350928    -0.04870669

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 24  1   -231.2385551083  1.1491E-10  0.0000E+00  2.2524E-04  1.0000E-03   
 mr-sdci # 24  2   -231.1524855088  1.8171E-05  1.0485E-05  4.8572E-03  1.0000E-03   
 mr-sdci # 24  3   -231.1218708014  5.8184E-04  0.0000E+00  3.5200E-02  1.0000E-04   
 mr-sdci # 24  4   -230.8051126854  2.3701E-03  0.0000E+00  2.4247E-01  1.0000E-04   
 mr-sdci # 24  5   -230.2204325236  5.8470E-01  0.0000E+00  1.0450E+00  1.0000E-04   
 mr-sdci # 24  6   -228.0181828124  5.1275E-01  0.0000E+00  1.4964E+00  1.0000E-04   
 mr-sdci # 24  7   -227.5038360813  9.4317E-02  0.0000E+00  1.7289E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.196899
time for cinew                         3.044922
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  25

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.78789424
   ht   2    -0.00000000    -6.70182464
   ht   3     0.00000000     0.00000000    -6.67120994
   ht   4    -0.00000000     0.00000000     0.00000000    -6.35445182
   ht   5     0.00039643    -0.00269714    -0.00378829    -0.00129477    -0.00002790

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.905114       2.155349E-02   4.483707E-03   1.330377E-02   1.379106E-02
 ref    2  -2.122466E-02   0.855400       0.303038       2.573685E-03   4.582139E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.819683       0.732173       9.185222E-02   1.836140E-04   2.289793E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1     0.90511447     0.02155349     0.00448371     0.01330377     0.01379106
 ref:   2    -0.02122466     0.85539964     0.30303814     0.00257368     0.04582139

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 25  1   -231.2385551083  7.4283E-11  0.0000E+00  2.2521E-04  1.0000E-03   
 mr-sdci # 25  2   -231.1524924665  6.9577E-06  2.7779E-06  2.6845E-03  1.0000E-03   
 mr-sdci # 25  3   -231.1221007776  2.2998E-04  0.0000E+00  2.7401E-02  1.0000E-04   
 mr-sdci # 25  4   -230.8063549718  1.2423E-03  0.0000E+00  2.3785E-01  1.0000E-04   
 mr-sdci # 25  5   -228.4985964762 -1.7218E+00  0.0000E+00  2.0410E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.118652
time for cinew                         1.977417
time for eigenvalue solver             0.000122
time for vector access                 0.000000

          starting ci iteration  26

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6
   ht   1    -6.78789424
   ht   2    -0.00000000    -6.70182464
   ht   3     0.00000000     0.00000000    -6.67120994
   ht   4    -0.00000000     0.00000000     0.00000000    -6.35445182
   ht   5     0.00039643    -0.00269714    -0.00378829    -0.00129477    -0.00002790
   ht   6    -0.00028185     0.00147850     0.00184890     0.00015703    -0.00000175    -0.00000779

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.905115       2.147027E-02   4.969990E-03   1.187258E-02  -3.372704E-02   4.257109E-02
 ref    2  -2.122433E-02   0.855121       0.303941       2.174245E-03  -6.675632E-03   6.354781E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6
 ref    1   0.819683       0.731692       9.240486E-02   1.456856E-04   1.182078E-03   5.850622E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6
 ref:   1     0.90511489     0.02147027     0.00496999     0.01187258    -0.03372704     0.04257109
 ref:   2    -0.02122433     0.85512075     0.30394105     0.00217424    -0.00667563     0.06354781

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 26  1   -231.2385551084  1.2422E-10  0.0000E+00  2.2539E-04  1.0000E-03   
 mr-sdci # 26  2   -231.1524965807  4.1142E-06  1.6636E-06  2.2373E-03  1.0000E-03   
 mr-sdci # 26  3   -231.1222407978  1.4002E-04  0.0000E+00  2.6063E-02  1.0000E-04   
 mr-sdci # 26  4   -230.8072164210  8.6145E-04  0.0000E+00  2.4096E-01  1.0000E-04   
 mr-sdci # 26  5   -230.2211323381  1.7225E+00  0.0000E+00  1.1200E+00  1.0000E-04   
 mr-sdci # 26  6   -227.5104695383 -5.0771E-01  0.0000E+00  1.8229E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.148071
time for cinew                         2.120850
time for eigenvalue solver             0.000000
time for vector access                 0.000000

          starting ci iteration  27

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5         ht   6         ht   7
   ht   1    -6.78789424
   ht   2    -0.00000000    -6.70182464
   ht   3     0.00000000     0.00000000    -6.67120994
   ht   4    -0.00000000     0.00000000     0.00000000    -6.35445182
   ht   5     0.00039643    -0.00269714    -0.00378829    -0.00129477    -0.00002790
   ht   6    -0.00028185     0.00147850     0.00184890     0.00015703    -0.00000175    -0.00000779
   ht   7     0.00056943    -0.00107461    -0.00141247     0.00072529    -0.00000331     0.00000031    -0.00000307

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.905115       2.144204E-02   5.145399E-03   1.053185E-02   3.325168E-02  -3.057590E-02  -3.077740E-02
 ref    2  -2.122434E-02   0.855034       0.304148       2.570045E-03  -6.943409E-04  -5.281919E-02  -4.073863E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7
 ref    1   0.819684       0.731543       9.253251E-02   1.175251E-04   1.106156E-03   3.724752E-03   2.606885E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5        ev    6        ev    7

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5       v      6       v      7

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7
 ref:   1     0.90511502     0.02144204     0.00514540     0.01053185     0.03325168    -0.03057590    -0.03077740
 ref:   2    -0.02122434     0.85503402     0.30414804     0.00257004    -0.00069434    -0.05281919    -0.04073863

 trial vector basis is being transformed.  new dimension:   4

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 27  1   -231.2385551085  4.6596E-11  0.0000E+00  2.2523E-04  1.0000E-03   
 mr-sdci # 27  2   -231.1524981827  1.6020E-06  5.1863E-07  1.2193E-03  1.0000E-03   
 mr-sdci # 27  3   -231.1223000643  5.9266E-05  0.0000E+00  2.3812E-02  1.0000E-04   
 mr-sdci # 27  4   -230.8082395375  1.0231E-03  0.0000E+00  2.3999E-01  1.0000E-04   
 mr-sdci # 27  5   -230.5132035888  2.9207E-01  0.0000E+00  7.1921E-01  1.0000E-04   
 mr-sdci # 27  6   -227.5651757274  5.4706E-02  0.0000E+00  1.7656E+00  1.0000E-04   
 mr-sdci # 27  7   -227.4842076047 -1.9628E-02  0.0000E+00  1.6906E+00  1.0000E-04   
 
 root number  2 is used to define the new expansion vector.
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.167969
time for cinew                         3.047852
time for eigenvalue solver             0.000244
time for vector access                 0.000000

          starting ci iteration  28

 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:     1078425 2x:      256352 4x:       40114
All internal counts: zz :      557690 yy:     5173911 xx:      525198 ww:      577418
One-external counts: yz :     2196660 yx:     3245408 yw:     3330558
Two-external counts: yy :     2027211 ww:      296420 xx:      295854 xz:       12692 wz:       14788 wx:      488146
Three-ext.   counts: yx :      661113 yw:      693909

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       95779
 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4         ht   5
   ht   1    -6.78789424
   ht   2     0.00000000    -6.70183732
   ht   3     0.00000000    -0.00000000    -6.67163920
   ht   4     0.00000000    -0.00000000     0.00000000    -6.35757867
   ht   5    -0.00042518     0.00020028    -0.00043340     0.00022299    -0.00000111

          calcsovref: eigensolution overlap with references block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1  -0.905115       2.146402E-02  -5.002838E-03  -1.135164E-02   4.821206E-02
 ref    2   2.122424E-02   0.854996      -0.304282      -2.352919E-03  -1.622481E-02

          calcsovref: reference weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5
 ref    1   0.819683       0.731478       9.261280E-02   1.343960E-04   2.587647E-03

          calcsovref: sovlaa in eigenvector basis block   1

               ev    1        ev    2        ev    3        ev    4        ev    5

          calcsovref: aa weight per eigenvector block   1

              v      1       v      2       v      3       v      4       v      5

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5
 ref:   1    -0.90511482     0.02146402    -0.00500284    -0.01135164     0.04821206
 ref:   2     0.02122424     0.85499564    -0.30428239    -0.00235292    -0.01622481

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 28  1   -231.2385551085  4.2071E-11  0.0000E+00  2.2506E-04  1.0000E-03   
 mr-sdci # 28  2   -231.1524986650  4.8230E-07  1.7525E-07  7.1340E-04  1.0000E-03   
 mr-sdci # 28  3   -231.1223197219  1.9658E-05  0.0000E+00  2.2830E-02  1.0000E-04   
 mr-sdci # 28  4   -230.8088140316  5.7449E-04  0.0000E+00  2.3915E-01  1.0000E-04   
 mr-sdci # 28  5   -228.8093303807 -1.7039E+00  0.0000E+00  1.7768E+00  1.0000E-04   
 
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.124023
time for cinew                         1.992432
time for eigenvalue solver             0.000000
time for vector access                 0.000000

 mr-sdci  convergence criteria satisfied after 28 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci # 28  1   -231.2385551085  4.2071E-11  0.0000E+00  2.2506E-04  1.0000E-03   
 mr-sdci # 28  2   -231.1524986650  4.8230E-07  1.7525E-07  7.1340E-04  1.0000E-03   
 mr-sdci # 28  3   -231.1223197219  1.9658E-05  0.0000E+00  2.2830E-02  1.0000E-04   
 mr-sdci # 28  4   -230.8088140316  5.7449E-04  0.0000E+00  2.3915E-01  1.0000E-04   
 mr-sdci # 28  5   -228.8093303807 -1.7039E+00  0.0000E+00  1.7768E+00  1.0000E-04   

####################CIUDGINFO####################

   ci vector at position   1 energy= -231.238555108538
   ci vector at position   2 energy= -231.152498664989

################END OF CIUDGINFO################

 
    2 of the   6 expansion vectors are transformed.
    2 of the   5 matrix-vector products are transformed.

    2 expansion eigenvectors written to unit nvfile (= 11)
    2 matrix-vector products written to unit nhvfil (= 10)


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 1) =      -231.2385551085

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
                                               18

                                          orbital     3    4    5    6   42   43   44   80   81   82  110  111  139  140  141  155  171
                                              182

                                         symmetry   ag   ag   ag   ag   b3u  b3u  b3u  b2u  b2u  b2u  b1g  b1g  b1u  b1u  b1u  b2g  b3g
                                              au 

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       2  0.024799                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +       
                                              - 
 z*  1  1       6  0.583188                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +-    - 
                                                
 z*  1  1       7 -0.019509                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +               - 
                                             +- 
 z*  1  1       9  0.195433                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +    +-    - 
                                                
 z*  1  1      10 -0.017219                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +          - 
                                             +- 
 z*  1  1      11 -0.612385                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +    +- 
                                              - 
 z*  1  1      15 -0.103324                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-        +-    - 
                                                
 z*  1  1      18 -0.032813                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +    +-    - 
                                                
 z*  1  1      20 -0.083349                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -        +    +- 
                                              - 
 z*  1  1      22 -0.053882                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +     -   +-    - 
                                                
 z*  1  1      24  0.142503                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +          -   +- 
                                              - 
 z*  1  1      25 -0.017052                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +-   +-    - 
                                                
 z*  1  1      27 -0.024815                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -   +    +- 
                                              - 
 z*  1  1      28  0.045383                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +     -   +- 
                                              - 
 z*  1  1      29 -0.120402                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-    - 
                                             +- 
 z*  1  1      33  0.049696                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-        +    +- 
                                              - 
 z*  1  1      36  0.024698                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +     -   +    +- 
                                              - 
 z*  1  1      38 -0.052161                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +         +-    - 
                                             +- 
 z*  1  1      40 -0.019710                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +    +-    - 
                                             +- 
 y   1  1    3868  0.019883              2( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -        +     - 
                                                
 y   1  1    3869 -0.068109              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -        +     - 
                                                
 y   1  1    3871  0.028008              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -        +     - 
                                                
 y   1  1    3874  0.016720              8( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -        +     - 
                                                
 y   1  1    4051 -0.025404              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-         -   +     - 
                                                
 y   1  1    4053  0.010538              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-         -   +     - 
                                                
 y   1  1    4145 -0.012687              1( b1u)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-    - 
                                                
 y   1  1    4159 -0.013890              2( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-      
                                              - 
 y   1  1    4160  0.047393              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-      
                                              - 
 y   1  1    4162 -0.019493              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-      
                                              - 
 y   1  1    4165 -0.011665              8( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-      
                                              - 
 y   1  1    4173 -0.010945              1( au )   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-              -   +- 
                                                
 y   1  1    4214  0.014517              2( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +- 
                                              - 
 y   1  1    4215 -0.049764              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +- 
                                              - 
 y   1  1    4217  0.020554              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +- 
                                              - 
 y   1  1    4220  0.012308              8( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +- 
                                              - 
 y   1  1    4806 -0.010184              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -             +     - 
                                             +- 
 y   1  1    4869  0.011186              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-         -    - 
                                                
 y   1  1    5362  0.012817              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +               -    - 
                                             +- 
 y   1  1    5824 -0.010029              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-         -        +     - 
                                             +- 
 w   1  1 9171553  0.016215    3( b2g)   3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +               - 
                                                
 w   1  1 9188489 -0.017257    3( b2g)   3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +       
                                              - 
 w   1  1 9188496  0.010548    3( b2g)   5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +       
                                              - 

 ci coefficient statistics:
           rq > 0.1                6
      0.1> rq > 0.01              36
     0.01> rq > 0.001          26810
    0.001> rq > 0.0001        865094
   0.0001> rq > 0.00001      4008585
  0.00001> rq > 0.000001     6164504
 0.000001> rq                7966377
           all              19031412
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:      111498 2x:           0 4x:           0
All internal counts: zz :      557690 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       13954
  iref  icsf         v(icsf)             hv(icsf)
     1     1      0.004680476484     -1.080556105141
     2     2      0.024798982109     -5.722826589998
     3     3     -0.000517038444      0.118695609161
     4     4      0.007203463066     -1.662239284662
     5     5     -0.009348266860      2.157684186878
     6     6      0.583187538688   -134.454004513077
     7     7     -0.019508803299      4.502021325540
     8     8     -0.003156633090      0.728037018838
     9     9      0.195432733188    -45.055697594027
    10    10     -0.017218831034      3.973855655972
    11    11     -0.612384950203    141.183731488311
    12    12      0.000611285605     -0.141151151845
    13    13      0.000508272046     -0.117368930875
    14    14      0.002005818235     -0.462967820463
    15    15     -0.103323650836     23.825462479747
    16    16      0.000617455244     -0.142465152051
    17    17     -0.000645421634      0.148963685208
    18    18     -0.032812625592      7.567803349327
    19    19     -0.002204389257      0.508720503222
    20    20     -0.083348864315     19.219584109047
    21    21      0.001985110229     -0.458225572487
    22    22     -0.053882309644     12.423504043266
    23    23     -0.000126601499      0.029198322130
    24    24      0.142502800469    -32.860048622568
    25    25     -0.017051541359      3.931614528254
    26    26      0.000285483700     -0.065967986320
    27    27     -0.024814803980      5.723904882313
    28    28      0.045382827470    -10.464186497428
    29    29     -0.120401934441     27.763851748128
    30    30      0.000171224486     -0.039485884659
    31    31      0.001234510759     -0.285115900474
    32    32      0.000243429469     -0.056202335702
    33    33      0.049695638181    -11.463049443070
    34    34      0.004067174584     -0.938344383707
    35    35      0.000282122859     -0.065099784126
    36    36      0.024697931235     -5.696809224336
    37    37     -0.001287935902      0.297550104008
    38    38     -0.052160512916     12.031658324188
    39    39      0.009326644575     -2.151193379668
    40    40     -0.019709971347      4.546821349707

 number of reference csfs (nref) is    40.  root number (iroot) is  1.
 c0**2 =   0.82046314  c**2 (all zwalks) =   0.82131608

 pople ci energy extrapolation is computed with 30 correlated electrons.

 eref      =   -230.552885700372   "relaxed" cnot**2         =   0.820463137691
 eci       =   -231.238555108538   deltae = eci - eref       =  -0.685669408166
 eci+dv1   =   -231.361658042661   dv1 = (1-cnot**2)*deltae  =  -0.123102934123
 eci+dv2   =   -231.388595894611   dv2 = dv1 / cnot**2       =  -0.150040786073
 eci+dv3   =   -231.430625458672   dv3 = dv1 / (2*cnot**2-1) =  -0.192070350134
 eci+pople =   -231.410589320921   ( 30e- scaled deltae )    =  -0.857703620549


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 2) =      -231.1524986650

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
                                               18

                                          orbital     3    4    5    6   42   43   44   80   81   82  110  111  139  140  141  155  171
                                              182

                                         symmetry   ag   ag   ag   ag   b3u  b3u  b3u  b2u  b2u  b2u  b1g  b1g  b1u  b1u  b1u  b2g  b3g
                                              au 

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1  0.067969                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          - 
                                                
 z*  1  1       2 -0.085736                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +       
                                              - 
 z*  1  1       3  0.023985                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-         - 
                                                
 z*  1  1       6 -0.590338                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +-    - 
                                                
 z*  1  1       7 -0.081568                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +               - 
                                             +- 
 z*  1  1       8  0.015349                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-   +       
                                              - 
 z*  1  1       9  0.518691                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +    +-    - 
                                                
 z*  1  1      10  0.011579                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +          - 
                                             +- 
 z*  1  1      11 -0.394477                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +    +- 
                                              - 
 z*  1  1      15  0.095456                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-        +-    - 
                                                
 z*  1  1      16  0.023178                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-              - 
                                             +- 
 z*  1  1      18  0.032940                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +    +-    - 
                                                
 z*  1  1      20  0.042466                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -        +    +- 
                                              - 
 z*  1  1      22 -0.019857                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +     -   +-    - 
                                                
 z*  1  1      24  0.017977                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +          -   +- 
                                              - 
 z*  1  1      25 -0.057153                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +-   +-    - 
                                                
 z*  1  1      27 -0.054885                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -   +    +- 
                                              - 
 z*  1  1      28  0.056799                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +         +     -   +- 
                                              - 
 z*  1  1      29 -0.059143                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +              +-    - 
                                             +- 
 z*  1  1      31 -0.028163                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-   +    +-    - 
                                                
 z*  1  1      34 -0.018174                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +    +-   +-    - 
                                                
 z*  1  1      38  0.019052                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +         +-    - 
                                             +- 
 z*  1  1      39  0.012075                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +    +- 
                                              - 
 z*  1  1      40 -0.030687                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +    +-    - 
                                             +- 
 y   1  1    3857  0.013621              1( b3g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -        +-      
                                                
 y   1  1    3869  0.027756              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -        +     - 
                                                
 y   1  1    3871 -0.013648              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -        +     - 
                                                
 y   1  1    3963 -0.010534              1( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -    - 
                                                
 y   1  1    3965  0.015245              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -    - 
                                                
 y   1  1    4039 -0.010016              1( b3g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-         -   +-      
                                                
 y   1  1    4049 -0.012217              1( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-         -   +     - 
                                                
 y   1  1    4051 -0.024059              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-         -   +     - 
                                                
 y   1  1    4053  0.011321              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-         -   +     - 
                                                
 y   1  1    4158  0.010168              1( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-      
                                              - 
 y   1  1    4173 -0.017091              1( au )   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-              -   +- 
                                                
 y   1  1    4174 -0.010648              2( au )   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-              -   +- 
                                                
 y   1  1    4205  0.010593              3( b3g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +     - 
                                              - 
 y   1  1    4215 -0.027937              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +- 
                                              - 
 y   1  1    4217  0.011660              5( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-                  +- 
                                              - 
 y   1  1    4670  0.012457              6( b1u)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -        +    +-    - 
                                                
 y   1  1    4761 -0.017841              1( b3g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -             +-   +- 
                                                
 y   1  1    4869 -0.017680              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +-         -    - 
                                                
 y   1  1    4955 -0.015719              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -    -   +     - 
                                                
 y   1  1    5013  0.019414              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +     -    - 
                                                
 y   1  1    5064 -0.016573              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -        +-      
                                              - 
 y   1  1    5119  0.013402              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -             +- 
                                              - 
 y   1  1    5246  0.010753              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          -   +-      
                                              - 
 y   1  1    5362  0.012239              3( b2g)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +               -    - 
                                             +- 
 y   1  1    7615 -0.011080             11( b3u)   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-    -   +-             +-   +- 
                                                
 y   1  1   37774  0.010478              6( b2u)   +-   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-        +    +-    - 
                                                
 y   1  1   59134 -0.012258             12( ag )   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-   +         +     - 
                                              - 
 y   1  1   59646 -0.012852             12( ag )   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +-             +-   +- 
                                                
 y   1  1   60615  0.010570              9( b2u)   +-   +-   +-   +-   +-   +-   +-   +-    -   +-   +-   +-   +     -   +    +-    - 
                                                
 y   1  1  230181 -0.011020              9( b2u)   +-   +-   +-    -   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +- 
                                                

 ci coefficient statistics:
           rq > 0.1                3
      0.1> rq > 0.01              51
     0.01> rq > 0.001          26163
    0.001> rq > 0.0001        869689
   0.0001> rq > 0.00001      4687181
  0.00001> rq > 0.000001     8275483
 0.000001> rq                5172842
           all              19031412
 =========== Executing IN-CORE method ==========


====================================================================================================
Diagonal     counts:  0x:      111498 2x:           0 4x:           0
All internal counts: zz :      557690 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================




LOOPCOUNT per task:
task #     1:       23733    task #     2:       28437    task #     3:      238239    task #     4:      119119
task #     5:     1940382    task #     6:     2634426    task #     7:     2698545    task #     8:      474076
task #     9:     4691072    task #    10:      163317    task #    11:      163317    task #    12:      154171
task #    13:      154171    task #    14:       77085    task #    15:       77087    task #    16:       77087
task #    17:       79366    task #    18:      573900    task #    19:       86684    task #    20:       13954
  iref  icsf         v(icsf)             hv(icsf)
     1     1      0.067968746978    -15.669881020415
     2     2     -0.085736278829     19.765653552694
     3     3      0.023984640706     -5.529772665760
     4     4      0.000816165749     -0.188931094974
     5     5     -0.003052038393      0.703617170957
     6     6     -0.590337878465    136.051406316279
     7     7     -0.081568162221     18.804463222349
     8     8      0.015348732426     -3.538964593240
     9     9      0.518691071415   -119.557531532207
    10    10      0.011578980460     -2.668917817502
    11    11     -0.394476573595     90.904313983260
    12    12      0.004120699334     -0.950063920802
    13    13      0.009375326982     -2.162476155258
    14    14     -0.006978307371      1.610885943759
    15    15      0.095455917909    -22.003695916735
    16    16      0.023177993995     -5.347468104440
    17    17      0.006854963724     -1.580779484007
    18    18      0.032939795346     -7.593295255569
    19    19      0.006305692552     -1.453598179815
    20    20      0.042466107295     -9.787448498348
    21    21     -0.002430200377      0.560946492441
    22    22     -0.019857139278      4.580564119299
    23    23     -0.001160350611      0.267469644121
    24    24      0.017977239855     -4.141655180711
    25    25     -0.057153089031     13.176658256126
    26    26     -0.004769273033      1.100407529284
    27    27     -0.054884589246     12.655235191949
    28    28      0.056799132168    -13.094366119508
    29    29     -0.059143185649     13.631728253639
    30    30      0.000914089918     -0.210976554593
    31    31     -0.028162773228      6.493071940290
    32    32     -0.006748390742      1.557147933921
    33    33      0.004533469820     -1.044323101229
    34    34     -0.018173529991      4.190628571868
    35    35     -0.003542424257      0.817450850327
    36    36      0.009769755132     -2.252477319328
    37    37     -0.003318689285      0.764349657103
    38    38      0.019051804293     -4.392189130354
    39    39      0.012075197382     -2.784585833950
    40    40     -0.030687234462      7.076875948009

 number of reference csfs (nref) is    40.  root number (iroot) is  2.
 c0**2 =   0.82201769  c**2 (all zwalks) =   0.82327169

 pople ci energy extrapolation is computed with 30 correlated electrons.

 eref      =   -230.475446010070   "relaxed" cnot**2         =   0.822017692245
 eci       =   -231.152498664989   deltae = eci - eref       =  -0.677052654919
 eci+dv1   =   -231.273002058983   dv1 = (1-cnot**2)*deltae  =  -0.120503393994
 eci+dv2   =   -231.299093312684   dv2 = dv1 / cnot**2       =  -0.146594647695
 eci+dv3   =   -231.339605455851   dv3 = dv1 / (2*cnot**2-1) =  -0.187106790863
 eci+pople =   -231.320187605093   ( 30e- scaled deltae )    =  -0.844741595023
maximum overlap with reference    1(overlap= 0.90511)
weight of reference states=  0.8197

 information on vector: 1 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    2(overlap= 0.85500)
weight of reference states=  0.7315

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
 MR-CISD energy:  -231.23855511    -6.78789424
 residuum:     0.00022506
 deltae:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.90511482     0.02146402    -0.00500284    -0.01135164     0.04821206    -0.03057590    -0.03077740     0.00000000
 ref:   2     0.02122424     0.85499564    -0.30428239    -0.00235292    -0.01622481    -0.05281919    -0.04073863     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.90511482     0.02146402    -0.00500284    -0.01135164     0.04821206     0.00000000     0.00000000     0.00000000
 ref:   2     0.02122424     0.85499564    -0.30428239    -0.00235292    -0.01622481     0.00000000     0.00000000     0.00000000

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
   DYZ=   51682  DYX=   70131  DYW=   73894
   D0Z=   12972  D0Y=  119645  D0X=   15912  D0W=   17672
  DDZI=   27020 DDYI=  198402 DDXI=   27696 DDWI=   30254
  DDZE=       0 DDYE=   29922 DDXE=    4857 DDWE=    5335
================================================================================
Trace of MO density:    30.000000
   30  correlated and    12  frozen core electrons

          modens reordered block   1

               ag    1        ag    2        ag    3        ag    4        ag    5        ag    6        ag    7        ag    8
  ag    1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    3    0.00000        0.00000        1.98717      -1.991496E-05   1.160736E-03  -1.341489E-05   1.410750E-03  -8.973523E-05
  ag    4    0.00000        0.00000      -1.991496E-05    1.98078       4.143548E-06  -2.252080E-03   2.094868E-04  -8.207923E-04
  ag    5    0.00000        0.00000       1.160736E-03   4.143548E-06    1.97965       5.323003E-05   1.300648E-03  -1.539031E-04
  ag    6    0.00000        0.00000      -1.341489E-05  -2.252080E-03   5.323003E-05    1.97466      -9.406437E-05   1.650912E-03
  ag    7    0.00000        0.00000       1.410750E-03   2.094868E-04   1.300648E-03  -9.406437E-05   6.056367E-04  -2.103101E-06
  ag    8    0.00000        0.00000      -8.973523E-05  -8.207923E-04  -1.539031E-04   1.650912E-03  -2.103101E-06   2.211054E-04
  ag    9    0.00000        0.00000       1.797211E-04   2.702305E-04   8.953093E-04  -2.398525E-04   6.012951E-04   4.025056E-06
  ag   10    0.00000        0.00000       3.056486E-04   2.110424E-03   1.137270E-04  -9.779562E-04   2.435897E-05  -3.977535E-04
  ag   11    0.00000        0.00000      -3.232648E-03  -4.155008E-04  -3.444121E-04  -4.199077E-05  -9.686152E-04  -2.706241E-06
  ag   12    0.00000        0.00000       4.295285E-05  -2.715560E-04   3.934219E-04  -3.551503E-03  -2.225847E-06  -3.607229E-04
  ag   13    0.00000        0.00000       2.989841E-04   2.832038E-03  -1.702682E-04   1.573666E-03   3.228395E-07  -1.840105E-04
  ag   14    0.00000        0.00000      -3.246392E-03  -3.076383E-04   3.937717E-03  -2.338084E-04  -1.214718E-03   3.246749E-06
  ag   15    0.00000        0.00000       4.818612E-04   1.235007E-03  -3.107438E-04   5.972133E-04  -5.627553E-06  -5.960871E-04
  ag   16    0.00000        0.00000      -1.086044E-04  -1.747558E-03   9.732212E-04  -5.608782E-03   1.069180E-07  -5.844964E-04
  ag   17    0.00000        0.00000      -3.319279E-03  -6.974028E-04   1.160469E-03  -3.700166E-04  -7.310494E-04   1.608049E-06
  ag   18    0.00000        0.00000      -5.769609E-04   2.193467E-03   9.312593E-04  -1.053728E-03   4.266455E-06   3.611482E-04
  ag   19    0.00000        0.00000      -2.279261E-04   1.313401E-04   1.263269E-03  -4.388336E-04   3.870816E-04   1.842539E-06
  ag   20    0.00000        0.00000       3.530300E-03   2.838800E-05  -7.735165E-03   8.521678E-04   4.106630E-04   5.579682E-07
  ag   21    0.00000        0.00000       4.254702E-04   2.992728E-03  -1.858464E-04   1.498481E-03   4.035630E-06  -1.586429E-04
  ag   22    0.00000        0.00000       1.501298E-04   2.309140E-03  -5.692123E-04   1.191564E-03  -2.836693E-07   3.788573E-04
  ag   23    0.00000        0.00000      -1.599168E-04   4.793429E-05  -4.523368E-04   1.025269E-03  -8.571370E-06   5.002797E-04
  ag   24    0.00000        0.00000      -4.872823E-04   2.582757E-03   5.036002E-04   9.148504E-04   1.048593E-05   9.836854E-05
  ag   25    0.00000        0.00000      -1.759811E-03   3.615533E-05   2.091356E-03  -5.871193E-04   7.132395E-04   4.241973E-06
  ag   26    0.00000        0.00000       4.069359E-05   1.112390E-03  -1.827268E-04  -7.700432E-05  -6.229735E-06   1.270594E-04
  ag   27    0.00000        0.00000      -1.425563E-03  -3.049126E-04   1.826651E-03  -2.942830E-04  -4.873771E-04  -5.130740E-07
  ag   28    0.00000        0.00000       3.485878E-04   2.254525E-03  -5.008512E-05  -1.145765E-03   1.107957E-04   1.612648E-06
  ag   29    0.00000        0.00000      -8.197142E-04   2.482145E-04  -3.027806E-04  -3.204915E-04  -3.956815E-04   1.952749E-06
  ag   30    0.00000        0.00000      -2.263139E-05   2.416571E-04   4.045699E-05   2.317986E-03  -6.339160E-07  -1.214363E-04
  ag   31    0.00000        0.00000       1.773259E-03   3.650137E-05   3.779872E-04   2.466122E-04  -1.523936E-04  -2.072165E-06
  ag   32    0.00000        0.00000      -1.463834E-03  -3.405483E-05   4.493194E-03  -3.975186E-04   1.580453E-04  -5.670005E-07
  ag   33    0.00000        0.00000       3.636358E-06  -2.908767E-03   1.136570E-04  -3.836908E-03  -1.660070E-06   1.443747E-04
  ag   34    0.00000        0.00000      -2.939851E-04  -1.101085E-03   3.362518E-05   4.591530E-04  -2.160404E-06   1.477246E-04
  ag   35    0.00000        0.00000      -1.224843E-04   2.306280E-03  -4.404243E-05  -1.179466E-03  -6.142613E-08  -8.306750E-07
  ag   36    0.00000        0.00000       4.772757E-05   9.369309E-04   2.820257E-04  -3.884160E-03  -1.364935E-05   7.831193E-06
  ag   37    0.00000        0.00000      -9.933313E-04   7.752231E-05  -2.628097E-03  -5.297395E-04   1.474466E-04   1.063540E-06
  ag   38    0.00000        0.00000      -2.757057E-05  -2.882916E-03   1.459569E-04  -2.956563E-03  -5.929735E-07  -1.022755E-04
  ag   39    0.00000        0.00000      -4.233952E-05   5.435425E-04   1.164864E-05   5.725238E-04   2.942727E-07   1.162064E-05

               ag    9        ag   10        ag   11        ag   12        ag   13        ag   14        ag   15        ag   16
  ag    3   1.797211E-04   3.056486E-04  -3.232648E-03   4.295285E-05   2.989841E-04  -3.246392E-03   4.818612E-04  -1.086044E-04
  ag    4   2.702305E-04   2.110424E-03  -4.155008E-04  -2.715560E-04   2.832038E-03  -3.076383E-04   1.235007E-03  -1.747558E-03
  ag    5   8.953093E-04   1.137270E-04  -3.444121E-04   3.934219E-04  -1.702682E-04   3.937717E-03  -3.107438E-04   9.732212E-04
  ag    6  -2.398525E-04  -9.779562E-04  -4.199077E-05  -3.551503E-03   1.573666E-03  -2.338084E-04   5.972133E-04  -5.608782E-03
  ag    7   6.012951E-04   2.435897E-05  -9.686152E-04  -2.225847E-06   3.228395E-07  -1.214718E-03  -5.627553E-06   1.069180E-07
  ag    8   4.025056E-06  -3.977535E-04  -2.706241E-06  -3.607229E-04  -1.840105E-04   3.246749E-06  -5.960871E-04  -5.844964E-04
  ag    9   2.056512E-03   2.268162E-05  -4.972865E-04  -4.742032E-06  -1.388260E-05  -2.894353E-03  -2.320827E-05  -8.688273E-06
  ag   10   2.268162E-05   8.537129E-04  -2.333304E-05   5.298081E-04   5.600717E-04  -6.552858E-05   1.571629E-03   6.221557E-04
  ag   11  -4.972865E-04  -2.333304E-05   1.771503E-03   1.301678E-05   9.796083E-06   1.703568E-03   3.364780E-05   8.650288E-06
  ag   12  -4.742032E-06   5.298081E-04   1.301678E-05   7.483869E-04   2.309972E-05   2.312712E-06   6.478690E-04   1.536734E-03
  ag   13  -1.388260E-05   5.600717E-04   9.796083E-06   2.309972E-05   7.386670E-04   6.876489E-06   1.075004E-03  -5.209459E-04
  ag   14  -2.894353E-03  -6.552858E-05   1.703568E-03   2.312712E-06   6.876489E-06   5.736793E-03   1.391579E-05   1.200376E-05
  ag   15  -2.320827E-05   1.571629E-03   3.364780E-05   6.478690E-04   1.075004E-03   1.391579E-05   4.000906E-03   2.053312E-04
  ag   16  -8.688273E-06   6.221557E-04   8.650288E-06   1.536734E-03  -5.209459E-04   1.200376E-05   2.053312E-04   3.898068E-03
  ag   17   6.167625E-04  -2.151475E-05   1.771943E-03   9.944844E-06   2.383951E-06   6.112250E-04   3.439333E-06   7.252415E-06
  ag   18   3.862816E-06  -1.317790E-03  -2.827323E-05  -2.214983E-04  -9.766752E-04   1.655951E-05  -4.769736E-03   1.013665E-03
  ag   19   8.900610E-04   9.058912E-06  -3.158772E-04  -2.307753E-06  -5.084023E-06  -6.384891E-04  -2.017207E-05   3.551503E-06
  ag   20   1.438005E-03   2.652772E-05  -8.034372E-04  -3.452669E-06  -9.486157E-06  -4.143996E-03  -2.944223E-05  -1.331627E-05
  ag   21  -9.008567E-06   5.140946E-04   3.991726E-06  -1.041794E-04   9.545735E-04  -8.851389E-06   7.454971E-04  -9.827248E-04
  ag   22  -1.356315E-06  -4.485632E-04  -6.745607E-06  -1.044527E-03   4.613340E-04  -3.808767E-07  -4.919280E-04  -2.779535E-03
  ag   23  -3.827005E-06  -8.565219E-04  -4.118865E-06  -1.051345E-03  -1.651989E-04   1.579652E-05  -1.132001E-03  -2.561289E-03
  ag   24   2.008094E-05  -4.161754E-04  -1.894579E-05  -9.691155E-05  -3.169349E-04  -1.554359E-05  -2.150914E-03   6.171045E-04
  ag   25   1.662851E-03   2.558939E-05  -6.373309E-04  -1.344977E-05   2.471415E-06  -1.611847E-03   1.135623E-05  -4.010755E-05
  ag   26   4.362287E-07  -1.917246E-04   9.569389E-06  -2.809136E-04  -5.308539E-05   1.401670E-05   4.274746E-05  -8.914609E-04
  ag   27  -4.309089E-04  -2.014195E-05   9.974456E-04   7.026932E-06   5.431705E-06   1.790340E-03   2.470677E-07   1.428500E-05
  ag   28  -1.476909E-04   8.799319E-05  -2.519096E-04  -2.688419E-04   5.595881E-04   1.402618E-04  -1.724851E-04  -1.092980E-03
  ag   29   5.409807E-04   1.351238E-05   9.197791E-04  -6.746753E-05   1.503690E-04  -5.521493E-04  -5.306826E-05  -2.926891E-04
  ag   30   5.460946E-07   1.209110E-04   2.404988E-06   2.329187E-04   4.478458E-05  -8.499987E-06  -1.162691E-04   3.133417E-04
  ag   31  -6.921326E-04  -4.071121E-06   6.961041E-05   3.205002E-06   2.611407E-06   8.884105E-04   1.564802E-05   8.788694E-06
  ag   32   1.192035E-04   2.938871E-06  -1.762318E-04   1.697654E-06  -2.836481E-06   1.151663E-04  -4.348914E-06   8.204757E-06
  ag   33   3.374616E-06  -4.219853E-04  -4.863972E-06  -7.873553E-05  -4.182707E-04   5.475984E-06  -7.427946E-04   1.599018E-05
  ag   34   1.601561E-06  -2.988668E-04  -2.573156E-06  -2.909240E-04  -2.406436E-04   6.457266E-06  -5.632508E-04  -5.163755E-04
  ag   35   7.836916E-07  -3.663465E-05  -1.305299E-06   8.189114E-06  -4.741895E-05  -6.779947E-07  -1.736819E-04   4.610343E-05
  ag   36  -3.046383E-05  -8.449124E-05   1.826954E-05   1.154792E-05  -9.185024E-05   6.468915E-05  -3.145528E-04   8.829293E-06
  ag   37   3.391174E-04  -1.551314E-06  -2.101808E-04  -4.466080E-07  -1.016609E-05  -6.727653E-04  -3.348937E-05  -3.132183E-06
  ag   38  -2.415754E-06   1.787959E-04   4.090308E-06   2.148891E-04   2.706959E-05   1.361892E-06   3.171426E-04   4.688423E-04
  ag   39   9.177194E-07  -4.868905E-05  -1.598772E-06   4.003157E-05  -8.631433E-05  -1.370814E-06  -1.267668E-04   2.367805E-04

               ag   17        ag   18        ag   19        ag   20        ag   21        ag   22        ag   23        ag   24
  ag    3  -3.319279E-03  -5.769609E-04  -2.279261E-04   3.530300E-03   4.254702E-04   1.501298E-04  -1.599168E-04  -4.872823E-04
  ag    4  -6.974028E-04   2.193467E-03   1.313401E-04   2.838800E-05   2.992728E-03   2.309140E-03   4.793429E-05   2.582757E-03
  ag    5   1.160469E-03   9.312593E-04   1.263269E-03  -7.735165E-03  -1.858464E-04  -5.692123E-04  -4.523368E-04   5.036002E-04
  ag    6  -3.700166E-04  -1.053728E-03  -4.388336E-04   8.521678E-04   1.498481E-03   1.191564E-03   1.025269E-03   9.148504E-04
  ag    7  -7.310494E-04   4.266455E-06   3.870816E-04   4.106630E-04   4.035630E-06  -2.836693E-07  -8.571370E-06   1.048593E-05
  ag    8   1.608049E-06   3.611482E-04   1.842539E-06   5.579682E-07  -1.586429E-04   3.788573E-04   5.002797E-04   9.836854E-05
  ag    9   6.167625E-04   3.862816E-06   8.900610E-04   1.438005E-03  -9.008567E-06  -1.356315E-06  -3.827005E-06   2.008094E-05
  ag   10  -2.151475E-05  -1.317790E-03   9.058912E-06   2.652772E-05   5.140946E-04  -4.485632E-04  -8.565219E-04  -4.161754E-04
  ag   11   1.771943E-03  -2.827323E-05  -3.158772E-04  -8.034372E-04   3.991726E-06  -6.745607E-06  -4.118865E-06  -1.894579E-05
  ag   12   9.944844E-06  -2.214983E-04  -2.307753E-06  -3.452669E-06  -1.041794E-04  -1.044527E-03  -1.051345E-03  -9.691155E-05
  ag   13   2.383951E-06  -9.766752E-04  -5.084023E-06  -9.486157E-06   9.545735E-04   4.613340E-04  -1.651989E-04  -3.169349E-04
  ag   14   6.112250E-04   1.655951E-05  -6.384891E-04  -4.143996E-03  -8.851389E-06  -3.808767E-07   1.579652E-05  -1.554359E-05
  ag   15   3.439333E-06  -4.769736E-03  -2.017207E-05  -2.944223E-05   7.454971E-04  -4.919280E-04  -1.132001E-03  -2.150914E-03
  ag   16   7.252415E-06   1.013665E-03   3.551503E-06  -1.331627E-05  -9.827248E-04  -2.779535E-03  -2.561289E-03   6.171045E-04
  ag   17   2.607740E-03  -2.773375E-06   4.264440E-04  -7.173591E-04   7.614136E-07  -6.179950E-06  -8.817765E-06  -1.395445E-07
  ag   18  -2.773375E-06   7.520530E-03   2.699528E-05   9.112207E-06  -5.301604E-04  -1.665345E-04   9.468678E-05   4.559044E-03
  ag   19   4.264440E-04   2.699528E-05   1.073932E-03  -5.372294E-04  -3.952929E-06  -2.405565E-06  -8.623739E-06   2.793438E-05
  ag   20  -7.173591E-04   9.112207E-06  -5.372294E-04   4.742409E-03   6.364805E-06   8.185838E-06   8.787142E-06   4.512618E-06
  ag   21   7.614136E-07  -5.301604E-04  -3.952929E-06   6.364805E-06   1.525735E-03   1.037467E-03   2.115838E-04  -2.479564E-04
  ag   22  -6.179950E-06  -1.665345E-04  -2.405565E-06   8.185838E-06   1.037467E-03   2.259066E-03   1.996243E-03  -3.162522E-04
  ag   23  -8.817765E-06   9.468678E-05  -8.623739E-06   8.787142E-06   2.115838E-04   1.996243E-03   2.652003E-03  -4.618704E-04
  ag   24  -1.395445E-07   4.559044E-03   2.793438E-05   4.512618E-06  -2.479564E-04  -3.162522E-04  -4.618704E-04   3.938435E-03
  ag   25   5.450958E-04  -5.624278E-05   1.170353E-03  -7.790350E-04   7.043599E-06   2.027696E-05   7.284692E-06  -1.570789E-05
  ag   26   1.295517E-05  -8.371270E-04  -1.377536E-06  -1.385417E-05   2.506865E-05   6.626178E-04   9.830545E-04  -9.572392E-04
  ag   27   1.232198E-03   2.144454E-05   3.401397E-04  -1.735019E-03   3.007806E-06  -4.695972E-06  -6.695967E-06   9.461423E-06
  ag   28  -3.674768E-04   5.174895E-04   4.270144E-05  -2.199850E-04   1.174677E-03   1.246874E-03   7.315600E-04   3.616900E-04
  ag   29   1.347798E-03   1.437537E-04  -1.590698E-04   8.424668E-04   3.184795E-04   3.325444E-04   1.897924E-04   1.020461E-04
  ag   30   3.921531E-07   3.324656E-04  -4.057079E-06   1.644604E-05   2.559934E-04  -3.477900E-05   2.150413E-04   1.644231E-05
  ag   31  -4.932904E-04  -1.119044E-05  -7.194852E-04  -3.606740E-04  -2.794574E-06  -9.055309E-06  -7.081679E-06  -7.009307E-06
  ag   32   7.678777E-05   8.911645E-06   3.233871E-04  -1.194600E-03  -4.559128E-06  -7.526007E-06  -9.942795E-06   1.309433E-05
  ag   33  -3.653272E-07   4.401890E-05  -1.120531E-06  -1.179728E-06  -4.168743E-04   8.218454E-05   6.586197E-04  -9.718362E-04
  ag   34  -4.476676E-06   7.921101E-04   9.224178E-07  -1.617944E-06  -3.718214E-04   2.626189E-04   7.167604E-04   9.671434E-04
  ag   35  -2.663903E-06   4.098094E-04   7.810702E-07   3.848902E-06  -4.793646E-06   9.243634E-05   1.541505E-04   5.475839E-04
  ag   36   4.203963E-06   3.928967E-04  -1.086307E-05  -2.533667E-05  -2.329807E-05   1.392658E-04   2.700311E-04   5.102279E-05
  ag   37  -4.252729E-05   3.606402E-05   1.404733E-04   2.648351E-04  -1.762258E-06   1.459013E-05   2.303098E-05   1.177504E-05
  ag   38   2.479380E-06  -2.678971E-04  -2.193474E-06  -1.390562E-06  -5.405713E-05  -3.520430E-04  -4.043524E-04  -1.507643E-04
  ag   39  -4.272487E-07   2.170099E-04   9.569591E-07   2.747122E-06  -1.317567E-04  -1.780474E-04  -2.555457E-04   1.769713E-04

               ag   25        ag   26        ag   27        ag   28        ag   29        ag   30        ag   31        ag   32
  ag    3  -1.759811E-03   4.069359E-05  -1.425563E-03   3.485878E-04  -8.197142E-04  -2.263139E-05   1.773259E-03  -1.463834E-03
  ag    4   3.615533E-05   1.112390E-03  -3.049126E-04   2.254525E-03   2.482145E-04   2.416571E-04   3.650137E-05  -3.405483E-05
  ag    5   2.091356E-03  -1.827268E-04   1.826651E-03  -5.008512E-05  -3.027806E-04   4.045699E-05   3.779872E-04   4.493194E-03
  ag    6  -5.871193E-04  -7.700432E-05  -2.942830E-04  -1.145765E-03  -3.204915E-04   2.317986E-03   2.466122E-04  -3.975186E-04
  ag    7   7.132395E-04  -6.229735E-06  -4.873771E-04   1.107957E-04  -3.956815E-04  -6.339160E-07  -1.523936E-04   1.580453E-04
  ag    8   4.241973E-06   1.270594E-04  -5.130740E-07   1.612648E-06   1.952749E-06  -1.214363E-04  -2.072165E-06  -5.670005E-07
  ag    9   1.662851E-03   4.362287E-07  -4.309089E-04  -1.476909E-04   5.409807E-04   5.460946E-07  -6.921326E-04   1.192035E-04
  ag   10   2.558939E-05  -1.917246E-04  -2.014195E-05   8.799319E-05   1.351238E-05   1.209110E-04  -4.071121E-06   2.938871E-06
  ag   11  -6.373309E-04   9.569389E-06   9.974456E-04  -2.519096E-04   9.197791E-04   2.404988E-06   6.961041E-05  -1.762318E-04
  ag   12  -1.344977E-05  -2.809136E-04   7.026932E-06  -2.688419E-04  -6.746753E-05   2.329187E-04   3.205002E-06   1.697654E-06
  ag   13   2.471415E-06  -5.308539E-05   5.431705E-06   5.595881E-04   1.503690E-04   4.478458E-05   2.611407E-06  -2.836481E-06
  ag   14  -1.611847E-03   1.401670E-05   1.790340E-03   1.402618E-04  -5.521493E-04  -8.499987E-06   8.884105E-04   1.151663E-04
  ag   15   1.135623E-05   4.274746E-05   2.470677E-07  -1.724851E-04  -5.306826E-05  -1.162691E-04   1.564802E-05  -4.348914E-06
  ag   16  -4.010755E-05  -8.914609E-04   1.428500E-05  -1.092980E-03  -2.926891E-04   3.133417E-04   8.788694E-06   8.204757E-06
  ag   17   5.450958E-04   1.295517E-05   1.232198E-03  -3.674768E-04   1.347798E-03   3.921531E-07  -4.932904E-04   7.678777E-05
  ag   18  -5.624278E-05  -8.371270E-04   2.144454E-05   5.174895E-04   1.437537E-04   3.324656E-04  -1.119044E-05   8.911645E-06
  ag   19   1.170353E-03  -1.377536E-06   3.401397E-04   4.270144E-05  -1.590698E-04  -4.057079E-06  -7.194852E-04   3.233871E-04
  ag   20  -7.790350E-04  -1.385417E-05  -1.735019E-03  -2.199850E-04   8.424668E-04   1.644604E-05  -3.606740E-04  -1.194600E-03
  ag   21   7.043599E-06   2.506865E-05   3.007806E-06   1.174677E-03   3.184795E-04   2.559934E-04  -2.794574E-06  -4.559128E-06
  ag   22   2.027696E-05   6.626178E-04  -4.695972E-06   1.246874E-03   3.325444E-04  -3.477900E-05  -9.055309E-06  -7.526007E-06
  ag   23   7.284692E-06   9.830545E-04  -6.695967E-06   7.315600E-04   1.897924E-04   2.150413E-04  -7.081679E-06  -9.942795E-06
  ag   24  -1.570789E-05  -9.572392E-04   9.461423E-06   3.616900E-04   1.020461E-04   1.644231E-05  -7.009307E-06   1.309433E-05
  ag   25   3.569670E-03   2.215477E-05  -2.466548E-04  -1.718067E-05   1.071474E-04  -6.614029E-06  -2.696740E-04   1.235212E-03
  ag   26   2.215477E-05   8.463028E-04   3.572924E-06   9.360630E-05   2.713190E-05   5.691838E-05  -4.220660E-06   1.861056E-06
  ag   27  -2.466548E-04   3.572924E-06   1.525751E-03  -4.915702E-05   1.821665E-04  -2.224708E-06  -4.073831E-04   3.787771E-04
  ag   28  -1.718067E-05   9.360630E-05  -4.915702E-05   1.547298E-03   3.554261E-05   3.057499E-04   4.888058E-05   8.417266E-05
  ag   29   1.071474E-04   2.713190E-05   1.821665E-04   3.554261E-05   1.419579E-03   8.518866E-05  -1.956853E-04  -3.135671E-04
  ag   30  -6.614029E-06   5.691838E-05  -2.224708E-06   3.057499E-04   8.518866E-05   7.602438E-04  -3.802610E-07  -6.361248E-06
  ag   31  -2.696740E-04  -4.220660E-06  -4.073831E-04   4.888058E-05  -1.956853E-04  -3.802610E-07   1.157173E-03  -3.019335E-04
  ag   32   1.235212E-03   1.861056E-06   3.787771E-04   8.417266E-05  -3.135671E-04  -6.361248E-06  -3.019335E-04   1.730849E-03
  ag   33   1.607246E-05   3.176788E-04  -5.164235E-06  -3.375465E-04  -9.227757E-05   9.737519E-05   1.141102E-06   1.293923E-06
  ag   34  -3.789216E-06  -2.780949E-05  -1.726649E-06  -6.360572E-05  -2.044590E-05   2.556811E-05   1.253606E-06   1.203882E-06
  ag   35  -8.135773E-06  -2.340478E-04   2.813719E-06   3.901425E-04   1.031593E-04   2.640065E-04  -9.209988E-07  -2.947106E-06
  ag   36  -7.224873E-05   2.380272E-04   2.821850E-05   1.782875E-04   3.558541E-05   1.432662E-04   1.449291E-05  -4.213025E-05
  ag   37   7.833266E-04   2.189361E-05  -2.968123E-04  -7.079778E-06   1.093181E-04   1.319445E-05  -1.857687E-04   4.550886E-04
  ag   38  -4.732429E-06  -2.259161E-04   2.554362E-06  -9.712947E-05  -2.457741E-05   5.847971E-05   3.317702E-06  -1.460193E-06
  ag   39  -7.041445E-06  -1.310085E-04   1.536932E-06  -2.374665E-05  -5.835070E-06   6.272083E-06  -6.365271E-07  -5.408337E-07

               ag   33        ag   34        ag   35        ag   36        ag   37        ag   38        ag   39
  ag    3   3.636358E-06  -2.939851E-04  -1.224843E-04   4.772757E-05  -9.933313E-04  -2.757057E-05  -4.233952E-05
  ag    4  -2.908767E-03  -1.101085E-03   2.306280E-03   9.369309E-04   7.752231E-05  -2.882916E-03   5.435425E-04
  ag    5   1.136570E-04   3.362518E-05  -4.404243E-05   2.820257E-04  -2.628097E-03   1.459569E-04   1.164864E-05
  ag    6  -3.836908E-03   4.591530E-04  -1.179466E-03  -3.884160E-03  -5.297395E-04  -2.956563E-03   5.725238E-04
  ag    7  -1.660070E-06  -2.160404E-06  -6.142613E-08  -1.364935E-05   1.474466E-04  -5.929735E-07   2.942727E-07
  ag    8   1.443747E-04   1.477246E-04  -8.306750E-07   7.831193E-06   1.063540E-06  -1.022755E-04   1.162064E-05
  ag    9   3.374616E-06   1.601561E-06   7.836916E-07  -3.046383E-05   3.391174E-04  -2.415754E-06   9.177194E-07
  ag   10  -4.219853E-04  -2.988668E-04  -3.663465E-05  -8.449124E-05  -1.551314E-06   1.787959E-04  -4.868905E-05
  ag   11  -4.863972E-06  -2.573156E-06  -1.305299E-06   1.826954E-05  -2.101808E-04   4.090308E-06  -1.598772E-06
  ag   12  -7.873553E-05  -2.909240E-04   8.189114E-06   1.154792E-05  -4.466080E-07   2.148891E-04   4.003157E-05
  ag   13  -4.182707E-04  -2.406436E-04  -4.741895E-05  -9.185024E-05  -1.016609E-05   2.706959E-05  -8.631433E-05
  ag   14   5.475984E-06   6.457266E-06  -6.779947E-07   6.468915E-05  -6.727653E-04   1.361892E-06  -1.370814E-06
  ag   15  -7.427946E-04  -5.632508E-04  -1.736819E-04  -3.145528E-04  -3.348937E-05   3.171426E-04  -1.267668E-04
  ag   16   1.599018E-05  -5.163755E-04   4.610343E-05   8.829293E-06  -3.132183E-06   4.688423E-04   2.367805E-04
  ag   17  -3.653272E-07  -4.476676E-06  -2.663903E-06   4.203963E-06  -4.252729E-05   2.479380E-06  -4.272487E-07
  ag   18   4.401890E-05   7.921101E-04   4.098094E-04   3.928967E-04   3.606402E-05  -2.678971E-04   2.170099E-04
  ag   19  -1.120531E-06   9.224178E-07   7.810702E-07  -1.086307E-05   1.404733E-04  -2.193474E-06   9.569591E-07
  ag   20  -1.179728E-06  -1.617944E-06   3.848902E-06  -2.533667E-05   2.648351E-04  -1.390562E-06   2.747122E-06
  ag   21  -4.168743E-04  -3.718214E-04  -4.793646E-06  -2.329807E-05  -1.762258E-06  -5.405713E-05  -1.317567E-04
  ag   22   8.218454E-05   2.626189E-04   9.243634E-05   1.392658E-04   1.459013E-05  -3.520430E-04  -1.780474E-04
  ag   23   6.586197E-04   7.167604E-04   1.541505E-04   2.700311E-04   2.303098E-05  -4.043524E-04  -2.555457E-04
  ag   24  -9.718362E-04   9.671434E-04   5.475839E-04   5.102279E-05   1.177504E-05  -1.507643E-04   1.769713E-04
  ag   25   1.607246E-05  -3.789216E-06  -8.135773E-06  -7.224873E-05   7.833266E-04  -4.732429E-06  -7.041445E-06
  ag   26   3.176788E-04  -2.780949E-05  -2.340478E-04   2.380272E-04   2.189361E-05  -2.259161E-04  -1.310085E-04
  ag   27  -5.164235E-06  -1.726649E-06   2.813719E-06   2.821850E-05  -2.968123E-04   2.554362E-06   1.536932E-06
  ag   28  -3.375465E-04  -6.360572E-05   3.901425E-04   1.782875E-04  -7.079778E-06  -9.712947E-05  -2.374665E-05
  ag   29  -9.227757E-05  -2.044590E-05   1.031593E-04   3.558541E-05   1.093181E-04  -2.457741E-05  -5.835070E-06
  ag   30   9.737519E-05   2.556811E-05   2.640065E-04   1.432662E-04   1.319445E-05   5.847971E-05   6.272083E-06
  ag   31   1.141102E-06   1.253606E-06  -9.209988E-07   1.449291E-05  -1.857687E-04   3.317702E-06  -6.365271E-07
  ag   32   1.293923E-06   1.203882E-06  -2.947106E-06  -4.213025E-05   4.550886E-04  -1.460193E-06  -5.408337E-07
  ag   33   1.463067E-03  -9.965605E-05  -3.022184E-04   2.944167E-04   2.889183E-05   1.990493E-04  -1.927444E-04
  ag   34  -9.965605E-05   1.038376E-03   3.008245E-04   1.621431E-05   1.872009E-06  -6.397986E-05  -1.201391E-04
  ag   35  -3.022184E-04   3.008245E-04   8.550245E-04   2.397571E-05   2.664335E-06   1.246334E-04   2.282258E-04
  ag   36   2.944167E-04   1.621431E-05   2.397571E-05   4.052586E-04  -1.290797E-05   8.980241E-05  -2.010862E-05
  ag   37   2.889183E-05   1.872009E-06   2.664335E-06  -1.290797E-05   5.538447E-04   7.868396E-06  -2.062784E-06
  ag   38   1.990493E-04  -6.397986E-05   1.246334E-04   8.980241E-05   7.868396E-06   4.955173E-04   1.157896E-04
  ag   39  -1.927444E-04  -1.201391E-04   2.282258E-04  -2.010862E-05  -2.062784E-06   1.157896E-04   3.219951E-04

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98736576     1.98154522     1.97954395     1.97398104     0.01445780     0.01252392
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.01022412     0.00661539     0.00531394     0.00463613     0.00232976     0.00177791     0.00148544     0.00126012
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00091401     0.00077099     0.00046512     0.00046155     0.00029768     0.00024346     0.00013933     0.00012065
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00010568     0.00009445     0.00008784     0.00005108     0.00002170     0.00002046     0.00001472     0.00000835
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000639     0.00000537     0.00000304     0.00000129     0.00000021     0.00000010     0.00000003

          modens reordered block   1

               b3u   1        b3u   2        b3u   3        b3u   4        b3u   5        b3u   6        b3u   7        b3u   8
  b3u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   3    0.00000        0.00000        1.98459      -3.594985E-05   1.748097E-03  -3.953962E-04   1.099249E-04   9.764352E-04
  b3u   4    0.00000        0.00000      -3.594985E-05    1.97896       5.980638E-05   3.716078E-04  -9.871073E-04  -4.857994E-04
  b3u   5    0.00000        0.00000       1.748097E-03   5.980638E-05    1.97718      -2.421984E-03   8.242412E-05   2.548190E-03
  b3u   6    0.00000        0.00000      -3.953962E-04   3.716078E-04  -2.421984E-03   3.832409E-04  -1.728938E-06  -4.941725E-04
  b3u   7    0.00000        0.00000       1.099249E-04  -9.871073E-04   8.242412E-05  -1.728938E-06   1.024355E-04   3.422618E-06
  b3u   8    0.00000        0.00000       9.764352E-04  -4.857994E-04   2.548190E-03  -4.941725E-04   3.422618E-06   7.046149E-04
  b3u   9    0.00000        0.00000       8.918493E-06  -5.253712E-04   2.582386E-03  -5.915100E-04   4.738257E-06   6.666560E-04
  b3u  10    0.00000        0.00000      -3.380894E-04   2.166490E-03  -1.739335E-04  -2.884839E-06  -3.218697E-04  -1.824176E-07
  b3u  11    0.00000        0.00000      -1.568713E-03   1.380292E-04  -6.883030E-05   2.519505E-04  -1.540253E-06  -5.080111E-04
  b3u  12    0.00000        0.00000       4.568250E-04  -1.571787E-03   1.549642E-04  -3.414725E-05   4.341650E-04   3.927502E-05
  b3u  13    0.00000        0.00000       1.084264E-04  -1.345315E-05   1.212065E-03   9.111471E-04   8.405241E-06  -8.903784E-04
  b3u  14    0.00000        0.00000      -5.030584E-04  -2.906812E-05  -6.785502E-05   3.035080E-08  -4.287596E-04  -8.214040E-06
  b3u  15    0.00000        0.00000       3.732515E-04  -9.979130E-04   4.368616E-03  -5.756599E-04  -4.458381E-06   8.000835E-04
  b3u  16    0.00000        0.00000      -2.136509E-03  -5.373060E-04   4.564746E-03   2.335152E-04  -3.350858E-07  -3.869973E-04
  b3u  17    0.00000        0.00000       8.792962E-04  -4.753786E-04   1.642049E-03  -1.984531E-04  -3.505681E-07   4.848576E-04
  b3u  18    0.00000        0.00000      -8.126291E-04  -4.138363E-03   6.765650E-04   6.136241E-06  -4.066858E-04  -1.230254E-05
  b3u  19    0.00000        0.00000      -4.582449E-04   1.065773E-04  -2.539362E-03  -3.205857E-04  -3.790576E-06   1.848677E-04
  b3u  20    0.00000        0.00000      -1.007824E-04   1.389552E-03  -4.257008E-04   3.320035E-06  -2.604283E-04  -6.591760E-06
  b3u  21    0.00000        0.00000      -7.051867E-04  -7.056961E-04   1.129855E-03  -7.017118E-04  -1.438076E-07   9.154609E-04
  b3u  22    0.00000        0.00000      -4.972272E-04  -3.482419E-03   1.054762E-03  -1.443593E-07   2.112354E-04   2.311862E-06
  b3u  23    0.00000        0.00000      -5.371664E-04  -1.384010E-03   5.758737E-04  -1.969108E-06  -1.172351E-04   1.535015E-06
  b3u  24    0.00000        0.00000       1.655871E-03   2.466808E-04  -2.350545E-03   1.581439E-05   6.159952E-07   3.971290E-05
  b3u  25    0.00000        0.00000      -1.889562E-04  -2.032607E-04  -1.251942E-04  -3.455076E-04   1.051674E-06   4.703771E-04
  b3u  26    0.00000        0.00000       3.642071E-04   6.951511E-05  -1.011148E-03  -1.743937E-04   1.596757E-06   1.408454E-04
  b3u  27    0.00000        0.00000      -6.530543E-05   2.349222E-03   3.242248E-04   1.739083E-06   2.013690E-04  -2.038653E-06
  b3u  28    0.00000        0.00000      -7.598784E-04   1.145389E-04  -2.470164E-03   6.764153E-06  -4.262234E-06  -1.729517E-04
  b3u  29    0.00000        0.00000      -1.774554E-03   6.763596E-05   1.510137E-03   2.496424E-04  -2.413839E-06  -4.581619E-04
  b3u  30    0.00000        0.00000      -9.036105E-05  -1.388734E-03  -9.987551E-06  -4.088234E-06  -7.241148E-05   6.016851E-06
  b3u  31    0.00000        0.00000       2.447372E-04   1.162519E-04  -4.133927E-03   9.056717E-05  -1.018051E-06  -5.615573E-05
  b3u  32    0.00000        0.00000       1.132591E-03  -9.359733E-05   2.788487E-04  -6.591598E-05   2.049591E-07   9.639693E-05
  b3u  33    0.00000        0.00000       1.738024E-04  -2.368162E-03  -1.808759E-04  -4.040720E-07   9.617727E-05   1.162753E-06
  b3u  34    0.00000        0.00000      -1.839575E-03  -1.168470E-04  -9.759604E-04  -1.463323E-04   2.182730E-07   1.978177E-04
  b3u  35    0.00000        0.00000       2.681906E-04  -3.039177E-03   4.595821E-04  -7.455108E-07  -1.173774E-05   1.008834E-05
  b3u  36    0.00000        0.00000      -8.189941E-04  -1.139735E-03  -1.418293E-03  -1.059110E-06  -4.675178E-06  -2.208818E-05
  b3u  37    0.00000        0.00000      -2.035929E-03  -5.076429E-05   2.550551E-03  -4.756211E-05   5.522476E-08   6.066138E-05
  b3u  38    0.00000        0.00000       3.957258E-05   2.460657E-05   1.417751E-05   3.083784E-07   5.987433E-05   3.136650E-07
  b3u  39    0.00000        0.00000       2.684295E-05   1.048699E-03   2.323671E-05   4.497196E-07   1.172691E-05  -4.697951E-07

               b3u   9        b3u  10        b3u  11        b3u  12        b3u  13        b3u  14        b3u  15        b3u  16
  b3u   3   8.918493E-06  -3.380894E-04  -1.568713E-03   4.568250E-04   1.084264E-04  -5.030584E-04   3.732515E-04  -2.136509E-03
  b3u   4  -5.253712E-04   2.166490E-03   1.380292E-04  -1.571787E-03  -1.345315E-05  -2.906812E-05  -9.979130E-04  -5.373060E-04
  b3u   5   2.582386E-03  -1.739335E-04  -6.883030E-05   1.549642E-04   1.212065E-03  -6.785502E-05   4.368616E-03   4.564746E-03
  b3u   6  -5.915100E-04  -2.884839E-06   2.519505E-04  -3.414725E-05   9.111471E-04   3.035080E-08  -5.756599E-04   2.335152E-04
  b3u   7   4.738257E-06  -3.218697E-04  -1.540253E-06   4.341650E-04   8.405241E-06  -4.287596E-04  -4.458381E-06  -3.350858E-07
  b3u   8   6.666560E-04  -1.824176E-07  -5.080111E-04   3.927502E-05  -8.903784E-04  -8.214040E-06   8.000835E-04  -3.869973E-04
  b3u   9   1.124346E-03  -1.806359E-06  -7.429022E-05   8.499551E-05  -2.055310E-03  -5.853285E-06   9.093437E-04  -2.400418E-04
  b3u  10  -1.806359E-06   1.030840E-03  -5.028088E-08  -1.414640E-03  -4.825782E-05   1.468744E-03   2.888404E-05  -3.867777E-06
  b3u  11  -7.429022E-05  -5.028088E-08   9.570435E-04  -9.324230E-06   3.170506E-04  -5.710859E-06   2.430635E-04   1.546129E-03
  b3u  12   8.499551E-05  -1.414640E-03  -9.324230E-06   1.990736E-03  -1.505006E-04  -2.127625E-03  -3.421112E-05  -9.531325E-05
  b3u  13  -2.055310E-03  -4.825782E-05   3.170506E-04  -1.505006E-04   6.073791E-03  -1.135429E-04   4.875707E-04   3.455014E-03
  b3u  14  -5.853285E-06   1.468744E-03  -5.710859E-06  -2.127625E-03  -1.135429E-04   2.618345E-03  -6.019120E-06  -7.729537E-05
  b3u  15   9.093437E-04   2.888404E-05   2.430635E-04  -3.421112E-05   4.875707E-04  -6.019120E-06   3.150005E-03   3.011243E-03
  b3u  16  -2.400418E-04  -3.867777E-06   1.546129E-03  -9.531325E-05   3.455014E-03  -7.729537E-05   3.011243E-03   6.244957E-03
  b3u  17  -1.279075E-04   6.129572E-06  -7.005376E-04  -5.260585E-05   1.341792E-03  -5.525354E-06   8.000041E-04   4.303773E-04
  b3u  18  -2.601025E-05   1.547331E-03   4.564240E-06  -2.443832E-03  -1.686026E-05   3.616076E-03   5.683100E-05   1.536040E-05
  b3u  19   9.672458E-04   2.325033E-05   2.525802E-04   7.316164E-05  -3.081466E-03   7.770338E-05  -4.034301E-04  -1.483345E-03
  b3u  20  -1.171082E-05   7.746703E-04  -3.142610E-06  -1.014997E-03  -2.605620E-05   6.764697E-04  -5.500526E-06  -1.960516E-05
  b3u  21   1.433253E-03   1.893615E-05   2.128196E-04   5.342964E-05  -1.619592E-03  -3.721372E-06   3.014673E-03   2.138087E-03
  b3u  22   5.873948E-07  -5.778430E-04  -3.049383E-06   6.443881E-04   2.924179E-05  -2.448474E-05  -7.701093E-06  -9.923222E-06
  b3u  23  -1.222932E-06   4.582697E-04  -3.849383E-06  -7.812832E-04  -2.439705E-05   1.171403E-03   1.822350E-05  -1.656474E-05
  b3u  24  -8.482162E-05  -3.417370E-06  -5.840437E-04   3.772648E-05  -1.197109E-03   1.942665E-05  -1.423199E-03  -3.108277E-03
  b3u  25   4.799520E-04   3.770964E-06  -5.251105E-04   3.765820E-05  -1.185877E-03   1.214041E-05  -1.366898E-04  -1.436311E-03
  b3u  26   2.939508E-04  -2.342237E-06   9.526534E-05   1.767009E-05  -3.083172E-04  -7.904372E-06   7.124762E-05   2.116368E-04
  b3u  27   8.465782E-06  -7.063170E-04   8.898503E-06   1.054759E-03   7.399583E-06  -1.213012E-03  -2.161345E-05   1.463171E-06
  b3u  28   3.962380E-04   1.575061E-05   5.933967E-04   2.891284E-05  -1.447451E-03   3.706252E-05   1.433060E-05   1.215862E-04
  b3u  29  -2.796407E-04   1.669210E-06   5.537477E-04  -2.092477E-05   4.273226E-04   6.588734E-06  -5.084247E-04   4.125903E-04
  b3u  30   4.142730E-06   1.696633E-04  -6.778736E-06  -1.889111E-04  -1.315616E-05  -9.276084E-05   6.984082E-06  -3.819864E-06
  b3u  31  -2.086040E-04   8.577976E-07  -3.183670E-05  -1.768852E-05   4.536285E-04  -3.623030E-06   1.134690E-05  -9.454787E-05
  b3u  32   1.138082E-04   1.011885E-06  -2.299262E-05   5.881877E-06  -1.715912E-04  -3.083000E-07   1.963236E-04  -3.385804E-04
  b3u  33   3.210295E-06  -3.269539E-04   1.546688E-06   5.063231E-04   1.197519E-05  -4.916969E-04  -7.750130E-06   2.453597E-06
  b3u  34   2.964566E-04   2.772192E-06  -5.423577E-05   1.877810E-05  -5.338486E-04  -4.247422E-08   4.331738E-04   4.651561E-05
  b3u  35  -2.777406E-05   8.321387E-05  -2.632707E-05  -1.964721E-04   1.271825E-04   3.244097E-04   3.452806E-05   3.318865E-05
  b3u  36   7.361701E-05   3.317555E-05   6.519312E-05  -6.447414E-05  -3.394602E-04   1.319908E-04  -6.204921E-05  -8.264685E-05
  b3u  37   7.363837E-05   9.445517E-07   3.783474E-05  -2.845304E-06   1.077797E-04  -3.480287E-06   2.895038E-04   4.431007E-04
  b3u  38   8.185445E-07  -2.027497E-04  -3.527665E-08   2.918978E-04   8.666013E-06  -3.337263E-04  -6.209139E-06   2.342153E-07
  b3u  39  -4.964296E-07  -3.424562E-05  -8.038195E-08   4.326890E-05   1.608755E-06  -6.061212E-05  -2.617106E-06  -8.725018E-07

               b3u  17        b3u  18        b3u  19        b3u  20        b3u  21        b3u  22        b3u  23        b3u  24
  b3u   3   8.792962E-04  -8.126291E-04  -4.582449E-04  -1.007824E-04  -7.051867E-04  -4.972272E-04  -5.371664E-04   1.655871E-03
  b3u   4  -4.753786E-04  -4.138363E-03   1.065773E-04   1.389552E-03  -7.056961E-04  -3.482419E-03  -1.384010E-03   2.466808E-04
  b3u   5   1.642049E-03   6.765650E-04  -2.539362E-03  -4.257008E-04   1.129855E-03   1.054762E-03   5.758737E-04  -2.350545E-03
  b3u   6  -1.984531E-04   6.136241E-06  -3.205857E-04   3.320035E-06  -7.017118E-04  -1.443593E-07  -1.969108E-06   1.581439E-05
  b3u   7  -3.505681E-07  -4.066858E-04  -3.790576E-06  -2.604283E-04  -1.438076E-07   2.112354E-04  -1.172351E-04   6.159952E-07
  b3u   8   4.848576E-04  -1.230254E-05   1.848677E-04  -6.591760E-06   9.154609E-04   2.311862E-06   1.535015E-06   3.971290E-05
  b3u   9  -1.279075E-04  -2.601025E-05   9.672458E-04  -1.171082E-05   1.433253E-03   5.873948E-07  -1.222932E-06  -8.482162E-05
  b3u  10   6.129572E-06   1.547331E-03   2.325033E-05   7.746703E-04   1.893615E-05  -5.778430E-04   4.582697E-04  -3.417370E-06
  b3u  11  -7.005376E-04   4.564240E-06   2.525802E-04  -3.142610E-06   2.128196E-04  -3.049383E-06  -3.849383E-06  -5.840437E-04
  b3u  12  -5.260585E-05  -2.443832E-03   7.316164E-05  -1.014997E-03   5.342964E-05   6.443881E-04  -7.812832E-04   3.772648E-05
  b3u  13   1.341792E-03  -1.686026E-05  -3.081466E-03  -2.605620E-05  -1.619592E-03   2.924179E-05  -2.439705E-05  -1.197109E-03
  b3u  14  -5.525354E-06   3.616076E-03   7.770338E-05   6.764697E-04  -3.721372E-06  -2.448474E-05   1.171403E-03   1.942665E-05
  b3u  15   8.000041E-04   5.683100E-05  -4.034301E-04  -5.500526E-06   3.014673E-03  -7.701093E-06   1.822350E-05  -1.423199E-03
  b3u  16   4.303773E-04   1.536040E-05  -1.483345E-03  -1.960516E-05   2.138087E-03  -9.923222E-06  -1.656474E-05  -3.108277E-03
  b3u  17   1.272493E-03   4.132323E-05  -1.110311E-03  -2.529779E-06   1.339856E-04   1.945516E-05   1.699540E-05  -3.818032E-04
  b3u  18   4.132323E-05   6.679190E-03   6.608041E-05  -1.814885E-04   2.635150E-05   2.659979E-03   2.922506E-03  -4.260398E-05
  b3u  19  -1.110311E-03   6.608041E-05   1.946272E-03  -5.849828E-06   4.808274E-04   3.637099E-05   4.132958E-05   7.362508E-04
  b3u  20  -2.529779E-06  -1.814885E-04  -5.849828E-06   1.147979E-03  -2.526565E-05  -1.753651E-03  -3.163702E-04   1.940599E-05
  b3u  21   1.339856E-04   2.635150E-05   4.808274E-04  -2.526565E-05   4.918436E-03  -1.327796E-06   1.868237E-05  -1.550680E-03
  b3u  22   1.945516E-05   2.659979E-03   3.637099E-05  -1.753651E-03  -1.327796E-06   5.363259E-03   2.660113E-03  -2.794365E-05
  b3u  23   1.699540E-05   2.922506E-03   4.132958E-05  -3.163702E-04   1.868237E-05   2.660113E-03   2.129425E-03  -1.344237E-05
  b3u  24  -3.818032E-04  -4.260398E-05   7.362508E-04   1.940599E-05  -1.550680E-03  -2.794365E-05  -1.344237E-05   2.484453E-03
  b3u  25   2.579630E-04  -1.463702E-05   5.750848E-04   9.098746E-06  -5.175198E-04  -5.441829E-06  -4.272041E-07   7.251475E-04
  b3u  26   6.061746E-05  -1.576893E-05   2.829508E-04   2.314789E-06  -6.275227E-04   2.989200E-06  -5.564853E-06  -2.188250E-04
  b3u  27  -2.183806E-05  -1.044324E-03   4.971213E-06  -7.460762E-04   1.305154E-05   1.910021E-03   4.375226E-04  -6.310075E-06
  b3u  28  -1.019314E-03   1.411363E-05   1.138473E-03   6.230296E-06   1.079855E-03  -3.665281E-05  -6.219056E-06  -9.147400E-05
  b3u  29  -5.070330E-04   1.261007E-05   2.709962E-04   6.539096E-06  -1.299791E-03   1.827411E-05   8.649209E-06   1.044942E-04
  b3u  30   5.496005E-06  -3.546645E-04  -8.824594E-06   3.504219E-04   1.284417E-05   2.224786E-04   2.948282E-04  -4.944391E-06
  b3u  31   2.368236E-05  -1.730754E-06  -4.573091E-05   3.627032E-06  -4.191321E-05   6.889501E-07   4.044225E-06   4.858603E-04
  b3u  32  -5.759452E-05  -2.462137E-06   1.566983E-04   3.418839E-06  -2.072614E-04  -2.456653E-06   2.866579E-06   8.930410E-04
  b3u  33  -6.079705E-06  -7.660330E-04  -5.753895E-06  -2.953713E-04  -8.569616E-06  -2.950488E-04  -7.365721E-04   8.541063E-06
  b3u  34   5.017695E-05  -5.450968E-06   2.556813E-04  -1.907418E-06   7.253184E-04  -5.072374E-06  -8.469905E-07  -1.613949E-04
  b3u  35   7.962027E-05   7.509355E-04  -5.329561E-05  -2.364588E-06  -1.009672E-04   1.776758E-04   3.282321E-04   2.463308E-05
  b3u  36  -1.860754E-04   2.925804E-04   1.642048E-04  -2.314198E-06   2.856068E-04   6.616354E-05   1.277216E-04  -7.876898E-05
  b3u  37   1.105538E-04   5.368748E-06  -8.643975E-05  -3.432743E-06   3.092044E-04   5.102189E-06   2.116298E-06  -4.767616E-04
  b3u  38  -1.635186E-06  -3.363934E-04  -4.240191E-06  -1.559407E-04  -3.535638E-06   2.842537E-04  -1.760199E-06   2.104795E-07
  b3u  39  -1.231592E-06  -1.878812E-04  -3.028609E-06   4.334040E-05  -3.372604E-06  -3.028741E-04  -1.492177E-04   2.753886E-06

               b3u  25        b3u  26        b3u  27        b3u  28        b3u  29        b3u  30        b3u  31        b3u  32
  b3u   3  -1.889562E-04   3.642071E-04  -6.530543E-05  -7.598784E-04  -1.774554E-03  -9.036105E-05   2.447372E-04   1.132591E-03
  b3u   4  -2.032607E-04   6.951511E-05   2.349222E-03   1.145389E-04   6.763596E-05  -1.388734E-03   1.162519E-04  -9.359733E-05
  b3u   5  -1.251942E-04  -1.011148E-03   3.242248E-04  -2.470164E-03   1.510137E-03  -9.987551E-06  -4.133927E-03   2.788487E-04
  b3u   6  -3.455076E-04  -1.743937E-04   1.739083E-06   6.764153E-06   2.496424E-04  -4.088234E-06   9.056717E-05  -6.591598E-05
  b3u   7   1.051674E-06   1.596757E-06   2.013690E-04  -4.262234E-06  -2.413839E-06  -7.241148E-05  -1.018051E-06   2.049591E-07
  b3u   8   4.703771E-04   1.408454E-04  -2.038653E-06  -1.729517E-04  -4.581619E-04   6.016851E-06  -5.615573E-05   9.639693E-05
  b3u   9   4.799520E-04   2.939508E-04   8.465782E-06   3.962380E-04  -2.796407E-04   4.142730E-06  -2.086040E-04   1.138082E-04
  b3u  10   3.770964E-06  -2.342237E-06  -7.063170E-04   1.575061E-05   1.669210E-06   1.696633E-04   8.577976E-07   1.011885E-06
  b3u  11  -5.251105E-04   9.526534E-05   8.898503E-06   5.933967E-04   5.537477E-04  -6.778736E-06  -3.183670E-05  -2.299262E-05
  b3u  12   3.765820E-05   1.767009E-05   1.054759E-03   2.891284E-05  -2.092477E-05  -1.889111E-04  -1.768852E-05   5.881877E-06
  b3u  13  -1.185877E-03  -3.083172E-04   7.399583E-06  -1.447451E-03   4.273226E-04  -1.315616E-05   4.536285E-04  -1.715912E-04
  b3u  14   1.214041E-05  -7.904372E-06  -1.213012E-03   3.706252E-05   6.588734E-06  -9.276084E-05  -3.623030E-06  -3.083000E-07
  b3u  15  -1.366898E-04   7.124762E-05  -2.161345E-05   1.433060E-05  -5.084247E-04   6.984082E-06   1.134690E-05   1.963236E-04
  b3u  16  -1.436311E-03   2.116368E-04   1.463171E-06   1.215862E-04   4.125903E-04  -3.819864E-06  -9.454787E-05  -3.385804E-04
  b3u  17   2.579630E-04   6.061746E-05  -2.183806E-05  -1.019314E-03  -5.070330E-04   5.496005E-06   2.368236E-05  -5.759452E-05
  b3u  18  -1.463702E-05  -1.576893E-05  -1.044324E-03   1.411363E-05   1.261007E-05  -3.546645E-04  -1.730754E-06  -2.462137E-06
  b3u  19   5.750848E-04   2.829508E-04   4.971213E-06   1.138473E-03   2.709962E-04  -8.824594E-06  -4.573091E-05   1.566983E-04
  b3u  20   9.098746E-06   2.314789E-06  -7.460762E-04   6.230296E-06   6.539096E-06   3.504219E-04   3.627032E-06   3.418839E-06
  b3u  21  -5.175198E-04  -6.275227E-04   1.305154E-05   1.079855E-03  -1.299791E-03   1.284417E-05  -4.191321E-05  -2.072614E-04
  b3u  22  -5.441829E-06   2.989200E-06   1.910021E-03  -3.665281E-05   1.827411E-05   2.224786E-04   6.889501E-07  -2.456653E-06
  b3u  23  -4.272041E-07  -5.564853E-06   4.375226E-04  -6.219056E-06   8.649209E-06   2.948282E-04   4.044225E-06   2.866579E-06
  b3u  24   7.251475E-04  -2.188250E-04  -6.310075E-06  -9.147400E-05   1.044942E-04  -4.944391E-06   4.858603E-04   8.930410E-04
  b3u  25   1.366560E-03   6.117378E-04  -1.523124E-05  -3.946189E-04   1.413524E-04   3.328159E-06  -4.581353E-04   4.642487E-04
  b3u  26   6.117378E-04   1.272975E-03  -7.221002E-06  -3.060194E-04   3.562487E-04   4.355341E-06  -6.023874E-04   2.829602E-04
  b3u  27  -1.523124E-05  -7.221002E-06   1.656550E-03  -3.361464E-06   8.026843E-06   4.084135E-04   1.280680E-05  -9.693022E-06
  b3u  28  -3.946189E-04  -3.060194E-04  -3.361464E-06   1.460264E-03   1.892479E-04  -1.428737E-05   3.263785E-04  -2.807851E-04
  b3u  29   1.413524E-04   3.562487E-04   8.026843E-06   1.892479E-04   1.257539E-03  -6.346504E-06  -5.962044E-05   1.291004E-04
  b3u  30   3.328159E-06   4.355341E-06   4.084135E-04  -1.428737E-05  -6.346504E-06   6.391343E-04  -3.204485E-06   8.133840E-07
  b3u  31  -4.581353E-04  -6.023874E-04   1.280680E-05   3.263785E-04  -5.962044E-05  -3.204485E-06   1.205492E-03  -1.632602E-04
  b3u  32   4.642487E-04   2.829602E-04  -9.693022E-06  -2.807851E-04   1.291004E-04   8.133840E-07  -1.632602E-04   1.110776E-03
  b3u  33   2.393824E-06   1.949516E-06  -6.272546E-05   2.771535E-06   6.906631E-07  -2.373501E-04  -1.152626E-06   1.847349E-06
  b3u  34   4.466814E-04   9.435277E-05  -1.648221E-06   9.326973E-05  -6.006700E-05   3.330070E-06  -3.834470E-04   2.067134E-04
  b3u  35   3.039788E-05  -2.898500E-05  -3.756368E-04  -7.108946E-05   5.707043E-05  -1.505319E-04   8.709534E-05   1.820539E-05
  b3u  36  -7.658254E-05   6.233640E-05  -1.459188E-04   2.035509E-04  -1.519929E-04  -5.675925E-05  -2.237100E-04  -4.813208E-05
  b3u  37  -6.214769E-05  -3.719471E-05   2.646157E-06   2.196126E-05   1.257676E-04  -6.712628E-07  -1.017707E-04  -3.670720E-04
  b3u  38  -1.305597E-06   5.994000E-07   2.361058E-04  -4.424718E-06   3.326932E-07  -1.672625E-06   9.947428E-07  -4.995103E-07
  b3u  39   4.676824E-07   7.049893E-07  -7.233216E-05   3.304666E-07  -1.262049E-06  -5.591186E-05  -4.640262E-07   8.989474E-07

               b3u  33        b3u  34        b3u  35        b3u  36        b3u  37        b3u  38        b3u  39
  b3u   3   1.738024E-04  -1.839575E-03   2.681906E-04  -8.189941E-04  -2.035929E-03   3.957258E-05   2.684295E-05
  b3u   4  -2.368162E-03  -1.168470E-04  -3.039177E-03  -1.139735E-03  -5.076429E-05   2.460657E-05   1.048699E-03
  b3u   5  -1.808759E-04  -9.759604E-04   4.595821E-04  -1.418293E-03   2.550551E-03   1.417751E-05   2.323671E-05
  b3u   6  -4.040720E-07  -1.463323E-04  -7.455108E-07  -1.059110E-06  -4.756211E-05   3.083784E-07   4.497196E-07
  b3u   7   9.617727E-05   2.182730E-07  -1.173774E-05  -4.675178E-06   5.522476E-08   5.987433E-05   1.172691E-05
  b3u   8   1.162753E-06   1.978177E-04   1.008834E-05  -2.208818E-05   6.066138E-05   3.136650E-07  -4.697951E-07
  b3u   9   3.210295E-06   2.964566E-04  -2.777406E-05   7.361701E-05   7.363837E-05   8.185445E-07  -4.964296E-07
  b3u  10  -3.269539E-04   2.772192E-06   8.321387E-05   3.317555E-05   9.445517E-07  -2.027497E-04  -3.424562E-05
  b3u  11   1.546688E-06  -5.423577E-05  -2.632707E-05   6.519312E-05   3.783474E-05  -3.527665E-08  -8.038195E-08
  b3u  12   5.063231E-04   1.877810E-05  -1.964721E-04  -6.447414E-05  -2.845304E-06   2.918978E-04   4.326890E-05
  b3u  13   1.197519E-05  -5.338486E-04   1.271825E-04  -3.394602E-04   1.077797E-04   8.666013E-06   1.608755E-06
  b3u  14  -4.916969E-04  -4.247422E-08   3.244097E-04   1.319908E-04  -3.480287E-06  -3.337263E-04  -6.061212E-05
  b3u  15  -7.750130E-06   4.331738E-04   3.452806E-05  -6.204921E-05   2.895038E-04  -6.209139E-06  -2.617106E-06
  b3u  16   2.453597E-06   4.651561E-05   3.318865E-05  -8.264685E-05   4.431007E-04   2.342153E-07  -8.725018E-07
  b3u  17  -6.079705E-06   5.017695E-05   7.962027E-05  -1.860754E-04   1.105538E-04  -1.635186E-06  -1.231592E-06
  b3u  18  -7.660330E-04  -5.450968E-06   7.509355E-04   2.925804E-04   5.368748E-06  -3.363934E-04  -1.878812E-04
  b3u  19  -5.753895E-06   2.556813E-04  -5.329561E-05   1.642048E-04  -8.643975E-05  -4.240191E-06  -3.028609E-06
  b3u  20  -2.953713E-04  -1.907418E-06  -2.364588E-06  -2.314198E-06  -3.432743E-06  -1.559407E-04   4.334040E-05
  b3u  21  -8.569616E-06   7.253184E-04  -1.009672E-04   2.856068E-04   3.092044E-04  -3.535638E-06  -3.372604E-06
  b3u  22  -2.950488E-04  -5.072374E-06   1.776758E-04   6.616354E-05   5.102189E-06   2.842537E-04  -3.028741E-04
  b3u  23  -7.365721E-04  -8.469905E-07   3.282321E-04   1.277216E-04   2.116298E-06  -1.760199E-06  -1.492177E-04
  b3u  24   8.541063E-06  -1.613949E-04   2.463308E-05  -7.876898E-05  -4.767616E-04   2.104795E-07   2.753886E-06
  b3u  25   2.393824E-06   4.466814E-04   3.039788E-05  -7.658254E-05  -6.214769E-05  -1.305597E-06   4.676824E-07
  b3u  26   1.949516E-06   9.435277E-05  -2.898500E-05   6.233640E-05  -3.719471E-05   5.994000E-07   7.049893E-07
  b3u  27  -6.272546E-05  -1.648221E-06  -3.756368E-04  -1.459188E-04   2.646157E-06   2.361058E-04  -7.233216E-05
  b3u  28   2.771535E-06   9.326973E-05  -7.108946E-05   2.035509E-04   2.196126E-05  -4.424718E-06   3.304666E-07
  b3u  29   6.906631E-07  -6.006700E-05   5.707043E-05  -1.519929E-04   1.257676E-04   3.326932E-07  -1.262049E-06
  b3u  30  -2.373501E-04   3.330070E-06  -1.505319E-04  -5.675925E-05  -6.712628E-07  -1.672625E-06  -5.591186E-05
  b3u  31  -1.152626E-06  -3.834470E-04   8.709534E-05  -2.237100E-04  -1.017707E-04   9.947428E-07  -4.640262E-07
  b3u  32   1.847349E-06   2.067134E-04   1.820539E-05  -4.813208E-05  -3.670720E-04  -4.995103E-07   8.989474E-07
  b3u  33   8.380115E-04   1.374859E-06  -1.358733E-04  -5.395196E-05  -1.065903E-06   3.959517E-05  -9.993051E-05
  b3u  34   1.374859E-06   5.822677E-04  -3.961167E-05   1.034250E-04   2.429413E-05  -1.534147E-06  -4.617768E-07
  b3u  35  -1.358733E-04  -3.961167E-05   4.545409E-04   4.343919E-05  -4.578174E-06   7.359652E-05  -2.682762E-05
  b3u  36  -5.395196E-05   1.034250E-04   4.343919E-05   3.637304E-04   1.346437E-05   2.870248E-05  -1.054449E-05
  b3u  37  -1.065903E-06   2.429413E-05  -4.578174E-06   1.346437E-05   4.472537E-04  -3.825293E-07  -7.307516E-07
  b3u  38   3.959517E-05  -1.534147E-06   7.359652E-05   2.870248E-05  -3.825293E-07   2.510796E-04  -7.824458E-05
  b3u  39  -9.993051E-05  -4.617768E-07  -2.682762E-05  -1.054449E-05  -7.307516E-07  -7.824458E-05   1.102455E-04

Natural orbital populations,block 2
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98499830     1.97899301     1.97684393     0.01404657     0.01298406     0.01074558
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00845376     0.00472889     0.00332251     0.00209972     0.00199569     0.00124508     0.00080197     0.00076118
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00044200     0.00044047     0.00032138     0.00021548     0.00017269     0.00013282     0.00009846     0.00006392
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00004304     0.00003181     0.00002977     0.00001851     0.00001192     0.00000675     0.00000668     0.00000400
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000388     0.00000141     0.00000068     0.00000018     0.00000012     0.00000006     0.00000000

          modens reordered block   1

               b2u   1        b2u   2        b2u   3        b2u   4        b2u   5        b2u   6        b2u   7        b2u   8
  b2u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b2u   2    0.00000        1.98461       7.007674E-05  -1.730350E-03   6.143904E-04   1.639386E-04   1.424629E-03   9.914219E-05
  b2u   3    0.00000       7.007674E-05    1.97460      -5.064987E-05   2.134867E-05  -4.967244E-03   7.060504E-05   1.854489E-05
  b2u   4    0.00000      -1.730350E-03  -5.064987E-05    1.97721      -2.830540E-03   6.463429E-05  -2.937485E-03  -3.269010E-03
  b2u   5    0.00000       6.143904E-04   2.134867E-05  -2.830540E-03   3.884267E-04   8.270661E-06   5.018318E-04   5.955056E-04
  b2u   6    0.00000       1.639386E-04  -4.967244E-03   6.463429E-05   8.270661E-06   4.541500E-04   5.458618E-06   1.043023E-05
  b2u   7    0.00000       1.424629E-03   7.060504E-05  -2.937485E-03   5.018318E-04   5.458618E-06   7.154702E-04   6.724112E-04
  b2u   8    0.00000       9.914219E-05   1.854489E-05  -3.269010E-03   5.955056E-04   1.043023E-05   6.724112E-04   1.129216E-03
  b2u   9    0.00000       2.261002E-03  -4.186971E-04  -1.581158E-04   2.549976E-04   1.110066E-05   5.129291E-04   7.187183E-05
  b2u  10    0.00000       2.536732E-04   8.367176E-04   1.189212E-05  -9.238083E-04  -7.653850E-06  -9.120057E-04  -2.075576E-03
  b2u  11    0.00000      -2.393090E-04   8.039313E-03  -2.876183E-04  -1.130010E-05  -9.496193E-04   1.471478E-06  -2.342480E-05
  b2u  12    0.00000      -8.437591E-04  -1.040708E-03   4.402879E-03  -5.703821E-04  -4.827020E-06  -7.956871E-04  -8.888007E-04
  b2u  13    0.00000       2.507078E-03  -1.969646E-03   3.995537E-03   2.404321E-04   1.368994E-05   3.988528E-04   2.491944E-04
  b2u  14    0.00000      -1.589137E-03  -1.336341E-04   1.630222E-03  -1.970091E-04  -1.030377E-05  -4.828333E-04   1.399215E-04
  b2u  15    0.00000       9.592268E-04   4.691824E-04  -1.966457E-03  -3.251484E-04  -3.799828E-06  -1.948117E-04  -9.751566E-04
  b2u  16    0.00000      -2.047465E-04   6.089674E-03  -4.328395E-04   8.652039E-07  -9.669072E-04   1.572778E-05   1.833901E-07
  b2u  17    0.00000      -7.550575E-05  -8.637992E-04   8.050070E-04  -7.013454E-04  -1.175948E-05  -9.188376E-04  -1.423587E-03
  b2u  18    0.00000       1.481294E-03  -1.098416E-03   2.231552E-03  -1.587041E-05   2.998658E-06   4.017090E-05  -8.467662E-05
  b2u  19    0.00000      -3.356286E-04  -3.397253E-04  -5.097929E-04   3.471932E-04   9.758277E-06   4.727979E-04   4.787071E-04
  b2u  20    0.00000       5.637477E-05   9.702900E-04   1.723514E-04  -5.533821E-06   4.035928E-04  -1.126233E-05  -9.465979E-06
  b2u  21    0.00000      -1.774506E-04   2.539273E-05   3.547221E-04   1.741437E-04   1.151650E-05   1.405317E-04   2.916497E-04
  b2u  22    0.00000      -9.470676E-04   3.667110E-05   2.528448E-03  -8.280174E-06  -2.956537E-06  -1.733467E-04   3.979067E-04
  b2u  23    0.00000       1.026471E-04  -6.897248E-03   3.012823E-04   3.755161E-07   5.322621E-04  -6.950118E-06   1.079155E-06
  b2u  24    0.00000       2.304571E-03  -1.112452E-04   1.661727E-03   2.515983E-04   7.986269E-06   4.610501E-04   2.788632E-04
  b2u  25    0.00000      -4.708209E-04   2.533179E-04  -4.531770E-03   9.217920E-05   7.894864E-07   5.792878E-05   2.115230E-04
  b2u  26    0.00000       8.046806E-04  -1.329763E-04  -4.169568E-04   6.611405E-05   3.123362E-07   9.714072E-05   1.126392E-04
  b2u  27    0.00000       8.574619E-05  -3.553659E-03   1.068368E-04   1.313786E-06   2.345993E-05   1.881286E-06   2.422428E-06
  b2u  28    0.00000       1.818463E-03  -4.476019E-05  -8.447556E-04  -1.465184E-04  -2.518499E-06  -1.983053E-04  -2.957466E-04
  b2u  29    0.00000      -8.836098E-04   2.858522E-05   1.340365E-03  -6.490816E-07  -2.114494E-07  -2.569659E-05   7.682677E-05
  b2u  30    0.00000       1.908450E-03  -2.512730E-04   2.568564E-03  -4.740455E-05   9.707400E-08  -5.981561E-05  -7.351507E-05

               b2u   9        b2u  10        b2u  11        b2u  12        b2u  13        b2u  14        b2u  15        b2u  16
  b2u   2   2.261002E-03   2.536732E-04  -2.393090E-04  -8.437591E-04   2.507078E-03  -1.589137E-03   9.592268E-04  -2.047465E-04
  b2u   3  -4.186971E-04   8.367176E-04   8.039313E-03  -1.040708E-03  -1.969646E-03  -1.336341E-04   4.691824E-04   6.089674E-03
  b2u   4  -1.581158E-04   1.189212E-05  -2.876183E-04   4.402879E-03   3.995537E-03   1.630222E-03  -1.966457E-03  -4.328395E-04
  b2u   5   2.549976E-04  -9.238083E-04  -1.130010E-05  -5.703821E-04   2.404321E-04  -1.970091E-04  -3.251484E-04   8.652039E-07
  b2u   6   1.110066E-05  -7.653850E-06  -9.496193E-04  -4.827020E-06   1.368994E-05  -1.030377E-05  -3.799828E-06  -9.669072E-04
  b2u   7   5.129291E-04  -9.120057E-04   1.471478E-06  -7.956871E-04   3.988528E-04  -4.828333E-04  -1.948117E-04   1.572778E-05
  b2u   8   7.187183E-05  -2.075576E-03  -2.342480E-05  -8.888007E-04   2.491944E-04   1.399215E-04  -9.751566E-04   1.833901E-07
  b2u   9   9.591761E-04  -3.107797E-04  -5.695164E-06   2.371480E-04   1.543920E-03  -7.028448E-04   2.522537E-04  -1.308661E-06
  b2u  10  -3.107797E-04   6.080144E-03   4.115652E-05  -5.022547E-04  -3.442196E-03  -1.357170E-03   3.074766E-03  -1.153647E-05
  b2u  11  -5.695164E-06   4.115652E-05   2.082720E-03  -7.484402E-06  -3.224889E-05  -5.417987E-06   3.027879E-05   2.174097E-03
  b2u  12   2.371480E-04  -5.022547E-04  -7.484402E-06   3.130035E-03   3.009834E-03   8.086446E-04  -4.118173E-04  -1.964398E-05
  b2u  13   1.543920E-03  -3.442196E-03  -3.224889E-05   3.009834E-03   6.242943E-03   4.402091E-04  -1.475188E-03  -1.393064E-05
  b2u  14  -7.028448E-04  -1.357170E-03  -5.417987E-06   8.086446E-04   4.402091E-04   1.281401E-03  -1.113619E-03  -1.055944E-07
  b2u  15   2.522537E-04   3.074766E-03   3.027879E-05  -4.118173E-04  -1.475188E-03  -1.113619E-03   1.931248E-03   6.557065E-06
  b2u  16  -1.308661E-06  -1.153647E-05   2.174097E-03  -1.964398E-05  -1.393064E-05  -1.055944E-07   6.557065E-06   2.396914E-03
  b2u  17   2.139934E-04   1.602079E-03   2.353310E-05   3.006269E-03   2.163061E-03   1.405506E-04   4.701444E-04  -6.084170E-06
  b2u  18   5.824804E-04  -1.189316E-03  -8.922495E-06   1.428667E-03   3.113352E-03   3.894495E-04  -7.317790E-04  -5.372270E-06
  b2u  19   5.260625E-04  -1.191719E-03  -1.434000E-05   1.559116E-04   1.453608E-03  -2.498294E-04  -5.787705E-04  -3.871551E-06
  b2u  20  -9.022138E-07   2.543432E-05  -8.689678E-04   1.466293E-06  -9.627768E-06  -6.928535E-07   1.182477E-05  -1.073298E-03
  b2u  21  -9.643961E-05  -3.044672E-04  -2.316694E-05  -6.614090E-05  -2.148431E-04  -5.860739E-05  -2.795515E-04  -2.467084E-05
  b2u  22  -5.961599E-04  -1.440798E-03  -1.462801E-05  -4.435953E-06  -1.277460E-04   1.023762E-03  -1.130307E-03  -4.522087E-06
  b2u  23   3.847549E-06   4.397540E-06  -1.410515E-03   1.263759E-05   1.544526E-05  -2.411280E-06  -4.947772E-06  -1.670147E-03
  b2u  24   5.517696E-04  -4.285131E-04  -6.305603E-06  -5.120576E-04   4.050693E-04  -5.056615E-04   2.659306E-04   2.725872E-06
  b2u  25  -3.218930E-05  -4.592146E-04  -4.676346E-06   1.812444E-05  -8.687569E-05   2.536630E-05  -4.950646E-05  -8.973778E-07
  b2u  26   2.440109E-05  -1.715389E-04   1.649798E-06  -1.935588E-04   3.421531E-04   5.849035E-05  -1.556808E-04   8.434064E-06
  b2u  27   2.014214E-06  -6.093928E-06  -9.489817E-05  -9.526437E-07   1.559971E-05   1.489775E-06  -5.689173E-06  -2.872559E-04
  b2u  28  -5.403704E-05   5.376435E-04   5.777274E-06   4.243312E-04   4.256911E-05   4.697557E-05   2.573151E-04  -1.427323E-06
  b2u  29  -7.099456E-05  -3.583152E-04  -3.901033E-06   7.705475E-05   9.191101E-05   2.041013E-04  -1.700276E-04  -1.735253E-06
  b2u  30   3.850405E-05  -1.022512E-04  -1.595847E-06   2.849616E-04   4.396735E-04   1.083311E-04  -8.329155E-05  -2.046640E-06

               b2u  17        b2u  18        b2u  19        b2u  20        b2u  21        b2u  22        b2u  23        b2u  24
  b2u   2  -7.550575E-05   1.481294E-03  -3.356286E-04   5.637477E-05  -1.774506E-04  -9.470676E-04   1.026471E-04   2.304571E-03
  b2u   3  -8.637992E-04  -1.098416E-03  -3.397253E-04   9.702900E-04   2.539273E-05   3.667110E-05  -6.897248E-03  -1.112452E-04
  b2u   4   8.050070E-04   2.231552E-03  -5.097929E-04   1.723514E-04   3.547221E-04   2.528448E-03   3.012823E-04   1.661727E-03
  b2u   5  -7.013454E-04  -1.587041E-05   3.471932E-04  -5.533821E-06   1.741437E-04  -8.280174E-06   3.755161E-07   2.515983E-04
  b2u   6  -1.175948E-05   2.998658E-06   9.758277E-06   4.035928E-04   1.151650E-05  -2.956537E-06   5.322621E-04   7.986269E-06
  b2u   7  -9.188376E-04   4.017090E-05   4.727979E-04  -1.126233E-05   1.405317E-04  -1.733467E-04  -6.950118E-06   4.610501E-04
  b2u   8  -1.423587E-03  -8.467662E-05   4.787071E-04  -9.465979E-06   2.916497E-04   3.979067E-04   1.079155E-06   2.788632E-04
  b2u   9   2.139934E-04   5.824804E-04   5.260625E-04  -9.022138E-07  -9.643961E-05  -5.961599E-04   3.847549E-06   5.517696E-04
  b2u  10   1.602079E-03  -1.189316E-03  -1.191719E-03   2.543432E-05  -3.044672E-04  -1.440798E-03   4.397540E-06  -4.285131E-04
  b2u  11   2.353310E-05  -8.922495E-06  -1.434000E-05  -8.689678E-04  -2.316694E-05  -1.462801E-05  -1.410515E-03  -6.305603E-06
  b2u  12   3.006269E-03   1.428667E-03   1.559116E-04   1.466293E-06  -6.614090E-05  -4.435953E-06   1.263759E-05  -5.120576E-04
  b2u  13   2.163061E-03   3.113352E-03   1.453608E-03  -9.627768E-06  -2.148431E-04  -1.277460E-04   1.544526E-05   4.050693E-04
  b2u  14   1.405506E-04   3.894495E-04  -2.498294E-04  -6.928535E-07  -5.860739E-05   1.023762E-03  -2.411280E-06  -5.056615E-04
  b2u  15   4.701444E-04  -7.317790E-04  -5.787705E-04   1.182477E-05  -2.795515E-04  -1.130307E-03  -4.947772E-06   2.659306E-04
  b2u  16  -6.084170E-06  -5.372270E-06  -3.871551E-06  -1.073298E-03  -2.467084E-05  -4.522087E-06  -1.670147E-03   2.725872E-06
  b2u  17   4.922977E-03   1.575310E-03   5.455668E-04  -1.975504E-05   6.319504E-04  -1.070994E-03   3.302776E-06  -1.308620E-03
  b2u  18   1.575310E-03   2.495032E-03   7.371599E-04  -4.804189E-06  -2.183142E-04  -9.520251E-05   7.050006E-06  -1.142070E-04
  b2u  19   5.455668E-04   7.371599E-04   1.376998E-03  -1.887540E-05   6.131101E-04  -3.998347E-04   2.746558E-06  -1.496215E-04
  b2u  20  -1.975504E-05  -4.804189E-06  -1.887540E-05   9.497547E-04  -8.713178E-06   1.000940E-05   4.866104E-04   8.948350E-06
  b2u  21   6.319504E-04  -2.183142E-04   6.131101E-04  -8.713178E-06   1.270669E-03  -3.067312E-04   1.007097E-05  -3.598959E-04
  b2u  22  -1.070994E-03  -9.520251E-05  -3.998347E-04   1.000940E-05  -3.067312E-04   1.458709E-03   2.992203E-07  -1.868036E-04
  b2u  23   3.302776E-06   7.050006E-06   2.746558E-06   4.866104E-04   1.007097E-05   2.992203E-07   2.035630E-03  -1.075315E-06
  b2u  24  -1.308620E-03  -1.142070E-04  -1.496215E-04   8.948350E-06  -3.598959E-04  -1.868036E-04  -1.075315E-06   1.256334E-03
  b2u  25  -3.227184E-05  -4.795229E-04   4.642846E-04  -1.414409E-05   6.060146E-04  -3.290882E-04   7.519383E-07  -6.300499E-05
  b2u  26   2.128382E-04   8.986818E-04   4.670530E-04  -1.554259E-05   2.807486E-04  -2.816270E-04  -1.318021E-05  -1.321476E-04
  b2u  27   8.160551E-06   2.612451E-05   1.323357E-05   2.182777E-04   1.069685E-05  -6.348928E-06   5.097783E-04  -3.971253E-06
  b2u  28   7.195428E-04   1.630673E-04  -4.440435E-04   4.875294E-06  -9.385792E-05  -9.189749E-05   2.915820E-07  -6.062052E-05
  b2u  29  -2.937935E-04  -7.858647E-05  -8.499311E-05  -3.437864E-07   6.796550E-05   2.136538E-04   9.691771E-07   1.632479E-04
  b2u  30   3.097029E-04   4.773771E-04   6.653173E-05  -1.603255E-06   3.588231E-05  -2.410702E-05   1.386225E-06   1.223661E-04

               b2u  25        b2u  26        b2u  27        b2u  28        b2u  29        b2u  30
  b2u   2  -4.708209E-04   8.046806E-04   8.574619E-05   1.818463E-03  -8.836098E-04   1.908450E-03
  b2u   3   2.533179E-04  -1.329763E-04  -3.553659E-03  -4.476019E-05   2.858522E-05  -2.512730E-04
  b2u   4  -4.531770E-03  -4.169568E-04   1.068368E-04  -8.447556E-04   1.340365E-03   2.568564E-03
  b2u   5   9.217920E-05   6.611405E-05   1.313786E-06  -1.465184E-04  -6.490816E-07  -4.740455E-05
  b2u   6   7.894864E-07   3.123362E-07   2.345993E-05  -2.518499E-06  -2.114494E-07   9.707400E-08
  b2u   7   5.792878E-05   9.714072E-05   1.881286E-06  -1.983053E-04  -2.569659E-05  -5.981561E-05
  b2u   8   2.115230E-04   1.126392E-04   2.422428E-06  -2.957466E-04   7.682677E-05  -7.351507E-05
  b2u   9  -3.218930E-05   2.440109E-05   2.014214E-06  -5.403704E-05  -7.099456E-05   3.850405E-05
  b2u  10  -4.592146E-04  -1.715389E-04  -6.093928E-06   5.376435E-04  -3.583152E-04  -1.022512E-04
  b2u  11  -4.676346E-06   1.649798E-06  -9.489817E-05   5.777274E-06  -3.901033E-06  -1.595847E-06
  b2u  12   1.812444E-05  -1.935588E-04  -9.526437E-07   4.243312E-04   7.705475E-05   2.849616E-04
  b2u  13  -8.687569E-05   3.421531E-04   1.559971E-05   4.256911E-05   9.191101E-05   4.396735E-04
  b2u  14   2.536630E-05   5.849035E-05   1.489775E-06   4.697557E-05   2.041013E-04   1.083311E-04
  b2u  15  -4.950646E-05  -1.556808E-04  -5.689173E-06   2.573151E-04  -1.700276E-04  -8.329155E-05
  b2u  16  -8.973778E-07   8.434064E-06  -2.872559E-04  -1.427323E-06  -1.735253E-06  -2.046640E-06
  b2u  17  -3.227184E-05   2.128382E-04   8.160551E-06   7.195428E-04  -2.937935E-04   3.097029E-04
  b2u  18  -4.795229E-04   8.986818E-04   2.612451E-05   1.630673E-04  -7.858647E-05   4.773771E-04
  b2u  19   4.642846E-04   4.670530E-04   1.323357E-05  -4.440435E-04  -8.499311E-05   6.653173E-05
  b2u  20  -1.414409E-05  -1.554259E-05   2.182777E-04   4.875294E-06  -3.437864E-07  -1.603255E-06
  b2u  21   6.060146E-04   2.807486E-04   1.069685E-05  -9.385792E-05   6.796550E-05   3.588231E-05
  b2u  22  -3.290882E-04  -2.816270E-04  -6.348928E-06  -9.189749E-05   2.136538E-04  -2.410702E-05
  b2u  23   7.519383E-07  -1.318021E-05   5.097783E-04   2.915820E-07   9.691771E-07   1.386225E-06
  b2u  24  -6.300499E-05  -1.321476E-04  -3.971253E-06  -6.062052E-05   1.632479E-04   1.223661E-04
  b2u  25   1.209276E-03   1.642077E-04   3.305707E-06  -3.860496E-04   2.371026E-04  -1.009072E-04
  b2u  26   1.642077E-04   1.111843E-03   1.160067E-05  -2.037145E-04  -5.124748E-05   3.689597E-04
  b2u  27   3.305707E-06   1.160067E-05   6.431000E-04  -4.411989E-06  -1.012054E-06   9.009753E-06
  b2u  28  -3.860496E-04  -2.037145E-04  -4.411989E-06   5.807513E-04  -1.087457E-04   2.368643E-05
  b2u  29   2.371026E-04  -5.124748E-05  -1.012054E-06  -1.087457E-04   3.445252E-04  -1.377468E-05
  b2u  30  -1.009072E-04   3.689597E-04   9.009753E-06   2.368643E-05  -1.377468E-05   4.462556E-04

Natural orbital populations,block 3
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98501230     1.97687893     1.97470115     0.01407435     0.01073017     0.00645401     0.00472697
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00332626     0.00210169     0.00124334     0.00110672     0.00079821     0.00067576     0.00043956     0.00032176
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00020544     0.00017368     0.00013222     0.00009845     0.00003181     0.00002982     0.00002014     0.00001192
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00000672     0.00000401     0.00000367     0.00000140     0.00000013     0.00000006

          modens reordered block   1

               b1g   1        b1g   2        b1g   3        b1g   4        b1g   5        b1g   6        b1g   7        b1g   8
  b1g   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b1g   2    0.00000        1.98081       2.210620E-03  -1.182671E-03  -2.762806E-03  -2.803404E-04  -6.670195E-06   3.136108E-03
  b1g   3    0.00000       2.210620E-03    1.97460      -1.857693E-03  -1.294194E-03  -3.782020E-03   3.395701E-05  -1.330283E-03
  b1g   4    0.00000      -1.182671E-03  -1.857693E-03   2.225673E-04   4.005374E-04   3.626550E-04   1.979400E-07  -1.851186E-04
  b1g   5    0.00000      -2.762806E-03  -1.294194E-03   4.005374E-04   8.567274E-04   5.345545E-04  -1.201842E-07  -5.593131E-04
  b1g   6    0.00000      -2.803404E-04  -3.782020E-03   3.626550E-04   5.345545E-04   7.512571E-04   8.483787E-07  -2.522790E-05
  b1g   7    0.00000      -6.670195E-06   3.395701E-05   1.979400E-07  -1.201842E-07   8.483787E-07   9.282496E-06   8.287138E-07
  b1g   8    0.00000       3.136108E-03  -1.330283E-03  -1.851186E-04  -5.593131E-04  -2.522790E-05   8.287138E-07   7.352781E-04
  b1g   9    0.00000       2.616427E-03  -8.410599E-04  -6.008565E-04  -1.578601E-03  -6.558188E-04   1.299035E-06   1.071835E-03
  b1g  10    0.00000       1.231202E-03  -6.093350E-03   5.832693E-04   6.200908E-04   1.537535E-03   2.691439E-06   5.233787E-04
  b1g  11    0.00000      -4.439862E-04  -2.620421E-03  -3.600528E-04  -1.314815E-03  -2.184034E-04   4.066619E-06   9.716306E-04
  b1g  12    0.00000      -1.489892E-04   5.818763E-04   5.534490E-06   2.186215E-05   1.978449E-06   8.960599E-05  -1.698317E-05
  b1g  13    0.00000      -3.233520E-03   1.139400E-03   1.613550E-04   5.179520E-04  -1.009225E-04  -1.149717E-06  -9.543201E-04
  b1g  14    0.00000       2.011694E-03  -1.369375E-03   3.812205E-04   4.540143E-04   1.051105E-03   1.533568E-06   4.609476E-04
  b1g  15    0.00000      -5.264252E-04   1.446429E-03  -2.244887E-06  -7.139208E-07  -7.248401E-06   1.913395E-04  -2.513663E-06
  b1g  16    0.00000      -3.111034E-04  -1.458616E-03   4.984449E-04   8.555436E-04   1.044792E-03   2.373154E-06  -1.696125E-04
  b1g  17    0.00000       1.433308E-03   6.671852E-04   9.664268E-05   4.127382E-04   9.490783E-05   2.336416E-08  -3.121225E-04
  b1g  18    0.00000       4.343854E-04  -8.483153E-04  -2.829242E-07  -3.399946E-06   9.768264E-07  -6.424733E-05   3.185365E-06
  b1g  19    0.00000       1.190546E-03  -3.878314E-04   1.273649E-04   1.912920E-04   2.836718E-04   2.680345E-07  -4.914873E-05
  b1g  20    0.00000       2.198254E-03   1.501181E-03  -9.506171E-07  -9.120860E-05   2.717817E-04   9.471249E-07   5.794052E-04
  b1g  21    0.00000       2.650256E-04  -2.262348E-03  -1.211325E-04  -1.211550E-04  -2.331905E-04  -3.568869E-07   4.247202E-05
  b1g  22    0.00000       2.074636E-04  -2.042639E-04   1.106234E-07   2.980874E-07   2.253035E-07  -1.099606E-06  -9.155000E-07
  b1g  23    0.00000       2.746902E-03  -3.192730E-03  -1.459993E-04  -4.245801E-04  -8.023964E-05   3.241539E-07   4.192519E-04
  b1g  24    0.00000      -1.295671E-03  -2.790712E-04   1.485455E-04   3.003318E-04   2.926352E-04   3.889577E-07  -2.396970E-04
  b1g  25    0.00000      -1.183452E-04   5.818531E-05   9.949518E-08   4.282974E-08   7.950162E-07  -3.894121E-05   3.174209E-07
  b1g  26    0.00000       2.215858E-03   1.433355E-03  -2.066117E-06   3.461084E-05  -1.034684E-05   8.114812E-08  -4.622183E-05
  b1g  27    0.00000       9.538544E-04   3.699286E-03   7.482826E-06   8.314242E-05  -1.193788E-05  -1.192928E-07  -9.148821E-05
  b1g  28    0.00000       2.646688E-03  -2.820640E-03   1.014360E-04   1.775930E-04   2.129524E-04   8.660479E-08  -2.723738E-05
  b1g  29    0.00000      -1.095303E-04   1.869572E-05   1.387047E-07   2.819877E-07   2.717134E-07  -5.132134E-06   5.825316E-09
  b1g  30    0.00000      -4.950808E-04   5.325268E-04  -1.108675E-05  -4.726594E-05   4.096774E-05   1.739945E-07   8.553970E-05

               b1g   9        b1g  10        b1g  11        b1g  12        b1g  13        b1g  14        b1g  15        b1g  16
  b1g   2   2.616427E-03   1.231202E-03  -4.439862E-04  -1.489892E-04  -3.233520E-03   2.011694E-03  -5.264252E-04  -3.111034E-04
  b1g   3  -8.410599E-04  -6.093350E-03  -2.620421E-03   5.818763E-04   1.139400E-03  -1.369375E-03   1.446429E-03  -1.458616E-03
  b1g   4  -6.008565E-04   5.832693E-04  -3.600528E-04   5.534490E-06   1.613550E-04   3.812205E-04  -2.244887E-06   4.984449E-04
  b1g   5  -1.578601E-03   6.200908E-04  -1.314815E-03   2.186215E-05   5.179520E-04   4.540143E-04  -7.139208E-07   8.555436E-04
  b1g   6  -6.558188E-04   1.537535E-03  -2.184034E-04   1.978449E-06  -1.009225E-04   1.051105E-03  -7.248401E-06   1.044792E-03
  b1g   7   1.299035E-06   2.691439E-06   4.066619E-06   8.960599E-05  -1.149717E-06   1.533568E-06   1.913395E-04   2.373154E-06
  b1g   8   1.071835E-03   5.233787E-04   9.716306E-04  -1.698317E-05  -9.543201E-04   4.609476E-04  -2.513663E-06  -1.696125E-04
  b1g   9   4.014291E-03  -1.930320E-04   4.765413E-03  -7.792734E-05  -7.594683E-04  -4.982022E-04  -7.567593E-08  -1.134274E-03
  b1g  10  -1.930320E-04   3.906922E-03   1.052536E-03  -2.020750E-05  -9.864093E-04   2.796633E-03  -1.570736E-05   2.535488E-03
  b1g  11   4.765413E-03   1.052536E-03   7.500446E-03  -1.020052E-04  -5.553888E-04   1.834413E-04   4.831344E-05  -8.721714E-05
  b1g  12  -7.792734E-05  -2.020750E-05  -1.020052E-04   9.733683E-04   9.330394E-06  -8.141052E-06   2.331100E-03   1.457456E-05
  b1g  13  -7.594683E-04  -9.864093E-04  -5.553888E-04   9.330394E-06   1.523983E-03  -1.039150E-03  -4.450355E-08  -2.013605E-04
  b1g  14  -4.982022E-04   2.796633E-03   1.834413E-04  -8.141052E-06  -1.039150E-03   2.284156E-03  -1.884003E-05   1.991434E-03
  b1g  15  -7.567593E-08  -1.570736E-05   4.831344E-05   2.331100E-03  -4.450355E-08  -1.884003E-05   6.269494E-03   2.812744E-05
  b1g  16  -1.134274E-03   2.535488E-03  -8.721714E-05   1.457456E-05  -2.013605E-04   1.991434E-03   2.812744E-05   2.615731E-03
  b1g  17  -2.140001E-03  -6.310876E-04  -4.532960E-03   8.359292E-05   2.596382E-04  -3.247425E-04   2.256402E-05  -4.603143E-04
  b1g  18   1.379924E-05   6.050268E-06   1.406367E-06  -1.095987E-03  -6.945264E-07   8.305960E-06  -3.672577E-03  -1.995829E-05
  b1g  19   4.037005E-05   9.015985E-04   8.351239E-04  -1.625182E-05  -3.600157E-05   6.783361E-04  -8.730516E-06   9.813053E-04
  b1g  20  -1.817683E-04   1.118350E-03  -5.410804E-04   9.922874E-06  -1.213555E-03   1.287220E-03   1.688858E-06   7.341919E-04
  b1g  21  -1.088217E-04  -3.147363E-04  -3.164239E-04   5.783614E-06  -2.481698E-04  -3.690246E-05   3.604004E-06   2.153748E-04
  b1g  22  -5.122946E-07   8.931557E-07  -5.422068E-06  -2.086122E-04   1.867295E-06   2.002521E-06  -1.075045E-03  -2.928527E-06
  b1g  23   7.532643E-04   2.293532E-05   6.273148E-05  -3.468472E-06  -4.179943E-04  -8.273674E-05  -7.470196E-06  -6.553259E-04
  b1g  24  -5.689129E-04   5.153588E-04  -7.985186E-04   1.584036E-05   3.713439E-04   2.668473E-04   4.303043E-06   7.125400E-04
  b1g  25  -1.179556E-06   1.494930E-06  -8.057098E-06  -2.333059E-04  -2.998054E-07   8.005208E-07  -6.659856E-05   3.375897E-08
  b1g  26  -1.702334E-04  -4.928918E-05  -4.036942E-04   8.424513E-06   4.489720E-06   9.230223E-05   2.965708E-06   1.525934E-04
  b1g  27  -3.159310E-04  -1.212480E-05  -3.956276E-04   7.813260E-06   2.461225E-05   1.387054E-04   4.512796E-06   2.666978E-04
  b1g  28  -3.128060E-04   4.638064E-04  -2.568570E-04   1.937481E-06  -5.157564E-05   3.504187E-04  -8.072742E-06   3.980092E-04
  b1g  29  -7.760211E-07   3.071100E-07  -8.448246E-07   2.166969E-05  -2.813136E-07  -2.041229E-07   2.389174E-04   1.425318E-06
  b1g  30   1.221625E-04   2.382525E-04   2.094785E-04  -3.726445E-06  -1.325000E-04   1.806460E-04  -1.503524E-06   2.537652E-04

               b1g  17        b1g  18        b1g  19        b1g  20        b1g  21        b1g  22        b1g  23        b1g  24
  b1g   2   1.433308E-03   4.343854E-04   1.190546E-03   2.198254E-03   2.650256E-04   2.074636E-04   2.746902E-03  -1.295671E-03
  b1g   3   6.671852E-04  -8.483153E-04  -3.878314E-04   1.501181E-03  -2.262348E-03  -2.042639E-04  -3.192730E-03  -2.790712E-04
  b1g   4   9.664268E-05  -2.829242E-07   1.273649E-04  -9.506171E-07  -1.211325E-04   1.106234E-07  -1.459993E-04   1.485455E-04
  b1g   5   4.127382E-04  -3.399946E-06   1.912920E-04  -9.120860E-05  -1.211550E-04   2.980874E-07  -4.245801E-04   3.003318E-04
  b1g   6   9.490783E-05   9.768264E-07   2.836718E-04   2.717817E-04  -2.331905E-04   2.253035E-07  -8.023964E-05   2.926352E-04
  b1g   7   2.336416E-08  -6.424733E-05   2.680345E-07   9.471249E-07  -3.568869E-07  -1.099606E-06   3.241539E-07   3.889577E-07
  b1g   8  -3.121225E-04   3.185365E-06  -4.914873E-05   5.794052E-04   4.247202E-05  -9.155000E-07   4.192519E-04  -2.396970E-04
  b1g   9  -2.140001E-03   1.379924E-05   4.037005E-05  -1.817683E-04  -1.088217E-04  -5.122946E-07   7.532643E-04  -5.689129E-04
  b1g  10  -6.310876E-04   6.050268E-06   9.015985E-04   1.118350E-03  -3.147363E-04   8.931557E-07   2.293532E-05   5.153588E-04
  b1g  11  -4.532960E-03   1.406367E-06   8.351239E-04  -5.410804E-04  -3.164239E-04  -5.422068E-06   6.273148E-05  -7.985186E-04
  b1g  12   8.359292E-05  -1.095987E-03  -1.625182E-05   9.922874E-06   5.783614E-06  -2.086122E-04  -3.468472E-06   1.584036E-05
  b1g  13   2.596382E-04  -6.945264E-07  -3.600157E-05  -1.213555E-03  -2.481698E-04   1.867295E-06  -4.179943E-04   3.713439E-04
  b1g  14  -3.247425E-04   8.305960E-06   6.783361E-04   1.287220E-03  -3.690246E-05   2.002521E-06  -8.273674E-05   2.668473E-04
  b1g  15   2.256402E-05  -3.672577E-03  -8.730516E-06   1.688858E-06   3.604004E-06  -1.075045E-03  -7.470196E-06   4.303043E-06
  b1g  16  -4.603143E-04  -1.995829E-05   9.813053E-04   7.341919E-04   2.153748E-04  -2.928527E-06  -6.553259E-04   7.125400E-04
  b1g  17   3.909991E-03  -3.310555E-05  -9.518015E-04   3.838475E-04   4.937534E-06  -2.590303E-06   9.553965E-04   9.722808E-04
  b1g  18  -3.310555E-05   2.902093E-03   8.375321E-06  -4.769317E-06  -4.864692E-06   1.176006E-03  -5.568381E-08  -1.181429E-05
  b1g  19  -9.518015E-04   8.375321E-06   8.463447E-04   9.973922E-05   6.197731E-05   2.810071E-06  -3.104717E-04  -3.045631E-05
  b1g  20   3.838475E-04  -4.769317E-06   9.973922E-05   1.554991E-03   3.108320E-04  -1.536192E-06   3.516063E-04  -6.196277E-05
  b1g  21   4.937534E-06  -4.864692E-06   6.197731E-05   3.108320E-04   7.552823E-04   3.309270E-07  -1.010147E-04   2.637575E-05
  b1g  22  -2.590303E-06   1.176006E-03   2.810071E-06  -1.536192E-06   3.309270E-07   6.493684E-04   8.353078E-07   4.237457E-08
  b1g  23   9.553965E-04  -5.568381E-08  -3.104717E-04   3.516063E-04  -1.010147E-04   8.353078E-07   1.458012E-03   1.016762E-04
  b1g  24   9.722808E-04  -1.181429E-05  -3.045631E-05  -6.196277E-05   2.637575E-05   4.237457E-08   1.016762E-04   1.038130E-03
  b1g  25   3.988116E-06  -5.958108E-04  -1.008489E-06   1.135912E-06   1.860812E-06  -4.576247E-04   3.492131E-06   3.799939E-06
  b1g  26   5.399516E-04  -5.372020E-06  -2.293165E-04   4.075953E-04   2.620335E-04   3.536877E-07   2.978132E-04   3.006869E-04
  b1g  27   5.340369E-05  -4.314756E-06   2.396685E-04   1.781709E-04   1.418652E-04  -8.298606E-07  -2.946740E-04   1.560019E-05
  b1g  28   1.433738E-04   4.053088E-06   2.261670E-04   9.394818E-05  -5.808828E-05   2.299337E-06   1.992033E-04   6.350135E-05
  b1g  29   2.292855E-06  -4.010334E-04  -9.514708E-07   5.179037E-07   1.254122E-06  -3.068568E-04   1.012182E-06   1.623600E-06
  b1g  30  -1.701583E-04   1.097821E-06   1.296721E-04   2.468459E-05  -5.639573E-06   3.678608E-07  -1.905251E-04   1.192566E-04

               b1g  25        b1g  26        b1g  27        b1g  28        b1g  29        b1g  30
  b1g   2  -1.183452E-04   2.215858E-03   9.538544E-04   2.646688E-03  -1.095303E-04  -4.950808E-04
  b1g   3   5.818531E-05   1.433355E-03   3.699286E-03  -2.820640E-03   1.869572E-05   5.325268E-04
  b1g   4   9.949518E-08  -2.066117E-06   7.482826E-06   1.014360E-04   1.387047E-07  -1.108675E-05
  b1g   5   4.282974E-08   3.461084E-05   8.314242E-05   1.775930E-04   2.819877E-07  -4.726594E-05
  b1g   6   7.950162E-07  -1.034684E-05  -1.193788E-05   2.129524E-04   2.717134E-07   4.096774E-05
  b1g   7  -3.894121E-05   8.114812E-08  -1.192928E-07   8.660479E-08  -5.132134E-06   1.739945E-07
  b1g   8   3.174209E-07  -4.622183E-05  -9.148821E-05  -2.723738E-05   5.825316E-09   8.553970E-05
  b1g   9  -1.179556E-06  -1.702334E-04  -3.159310E-04  -3.128060E-04  -7.760211E-07   1.221625E-04
  b1g  10   1.494930E-06  -4.928918E-05  -1.212480E-05   4.638064E-04   3.071100E-07   2.382525E-04
  b1g  11  -8.057098E-06  -4.036942E-04  -3.956276E-04  -2.568570E-04  -8.448246E-07   2.094785E-04
  b1g  12  -2.333059E-04   8.424513E-06   7.813260E-06   1.937481E-06   2.166969E-05  -3.726445E-06
  b1g  13  -2.998054E-07   4.489720E-06   2.461225E-05  -5.157564E-05  -2.813136E-07  -1.325000E-04
  b1g  14   8.005208E-07   9.230223E-05   1.387054E-04   3.504187E-04  -2.041229E-07   1.806460E-04
  b1g  15  -6.659856E-05   2.965708E-06   4.512796E-06  -8.072742E-06   2.389174E-04  -1.503524E-06
  b1g  16   3.375897E-08   1.525934E-04   2.666978E-04   3.980092E-04   1.425318E-06   2.537652E-04
  b1g  17   3.988116E-06   5.399516E-04   5.340369E-05   1.433738E-04   2.292855E-06  -1.701583E-04
  b1g  18  -5.958108E-04  -5.372020E-06  -4.314756E-06   4.053088E-06  -4.010334E-04   1.097821E-06
  b1g  19  -1.008489E-06  -2.293165E-04   2.396685E-04   2.261670E-04  -9.514708E-07   1.296721E-04
  b1g  20   1.135912E-06   4.075953E-04   1.781709E-04   9.394818E-05   5.179037E-07   2.468459E-05
  b1g  21   1.860812E-06   2.620335E-04   1.418652E-04  -5.808828E-05   1.254122E-06  -5.639573E-06
  b1g  22  -4.576247E-04   3.536877E-07  -8.298606E-07   2.299337E-06  -3.068568E-04   3.678608E-07
  b1g  23   3.492131E-06   2.978132E-04  -2.946740E-04   1.992033E-04   1.012182E-06  -1.905251E-04
  b1g  24   3.799939E-06   3.006869E-04   1.560019E-05   6.350135E-05   1.623600E-06   1.192566E-04
  b1g  25   6.142056E-04   1.973594E-06   2.116637E-07  -2.479706E-07   2.862154E-04  -3.246205E-07
  b1g  26   1.973594E-06   8.536043E-04   2.194649E-05  -1.262358E-04   8.935726E-07  -2.274790E-04
  b1g  27   2.116637E-07   2.194649E-05   4.018027E-04  -9.004587E-05   2.950024E-08   2.131592E-05
  b1g  28  -2.479706E-07  -1.262358E-04  -9.004587E-05   4.945196E-04  -2.569934E-07   1.159120E-04
  b1g  29   2.862154E-04   8.935726E-07   2.950024E-08  -2.569934E-07   3.331860E-04  -1.375261E-07
  b1g  30  -3.246205E-07  -2.274790E-04   2.131592E-05   1.159120E-04  -1.375261E-07   3.209873E-04

Natural orbital populations,block 4
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98154911     1.97395654     0.01442896     0.01021310     0.00969454     0.00463151     0.00232492
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00175134     0.00148403     0.00091267     0.00077211     0.00046487     0.00029700     0.00021482     0.00013884
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00010548     0.00008762     0.00006676     0.00005109     0.00002041     0.00001908     0.00000838     0.00000639
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00000312     0.00000129     0.00000108     0.00000021     0.00000003     0.00000002

          modens reordered block   1

               b1u   1        b1u   2        b1u   3        b1u   4        b1u   5        b1u   6        b1u   7        b1u   8
  b1u   1    1.88869       2.832354E-03  -8.103390E-03   2.924181E-05   2.339145E-03  -3.763451E-05  -4.720978E-03   1.966038E-04
  b1u   2   2.832354E-03   0.492614       0.167691       4.010499E-03  -1.286460E-03   2.383854E-03   1.638635E-03  -6.641533E-04
  b1u   3  -8.103390E-03   0.167691       5.950504E-02   1.207667E-03   1.410166E-04   5.584994E-04   1.064732E-03   7.909338E-05
  b1u   4   2.924181E-05   4.010499E-03   1.207667E-03   9.138089E-04   7.139363E-05   7.973990E-04   4.757690E-05  -8.083016E-04
  b1u   5   2.339145E-03  -1.286460E-03   1.410166E-04   7.139363E-05   2.018694E-03   2.245769E-05  -1.558018E-03  -2.123896E-05
  b1u   6  -3.763451E-05   2.383854E-03   5.584994E-04   7.973990E-04   2.245769E-05   1.073574E-03   1.082379E-05  -5.962372E-04
  b1u   7  -4.720978E-03   1.638635E-03   1.064732E-03   4.757690E-05  -1.558018E-03   1.082379E-05   2.081764E-03  -1.808948E-06
  b1u   8   1.966038E-04  -6.641533E-04   7.909338E-05  -8.083016E-04  -2.123896E-05  -5.962372E-04  -1.808948E-06   1.575309E-03
  b1u   9   7.755738E-03  -8.031473E-04  -2.381897E-03  -2.261228E-04  -6.124016E-04  -5.707552E-05  -7.011096E-04   4.355126E-05
  b1u  10  -3.885956E-05   1.367344E-03   2.570902E-04   6.978369E-04   1.408497E-05   1.201572E-03   1.196973E-05  -1.655711E-04
  b1u  11   4.312145E-05   2.411098E-03   1.036755E-03  -6.148355E-04  -2.548877E-05  -2.980060E-04   8.070891E-06   1.037284E-03
  b1u  12   4.014145E-03  -7.297017E-04   1.482372E-04   4.092056E-05   1.740838E-03   1.021104E-05  -1.404513E-03  -7.631061E-06
  b1u  13  -1.038354E-04  -1.154021E-03  -3.851094E-04  -1.127067E-04  -2.267756E-06  -3.862075E-04   3.095424E-06  -9.576219E-04
  b1u  14   1.771894E-05  -4.462331E-04  -6.571113E-05  -2.252046E-04  -5.459583E-06  -7.292143E-04  -1.289357E-06   3.975904E-04
  b1u  15   2.974795E-03  -4.292177E-04  -1.808568E-04  -4.814939E-06   2.431895E-04  -3.220287E-06  -6.233820E-04   4.240182E-06
  b1u  16   2.980244E-05   1.673947E-04   4.267007E-05   3.279377E-05   2.369117E-06   5.983435E-05  -4.629986E-06  -1.711243E-04

               b1u   9        b1u  10        b1u  11        b1u  12        b1u  13        b1u  14        b1u  15        b1u  16
  b1u   1   7.755738E-03  -3.885956E-05   4.312145E-05   4.014145E-03  -1.038354E-04   1.771894E-05   2.974795E-03   2.980244E-05
  b1u   2  -8.031473E-04   1.367344E-03   2.411098E-03  -7.297017E-04  -1.154021E-03  -4.462331E-04  -4.292177E-04   1.673947E-04
  b1u   3  -2.381897E-03   2.570902E-04   1.036755E-03   1.482372E-04  -3.851094E-04  -6.571113E-05  -1.808568E-04   4.267007E-05
  b1u   4  -2.261228E-04   6.978369E-04  -6.148355E-04   4.092056E-05  -1.127067E-04  -2.252046E-04  -4.814939E-06   3.279377E-05
  b1u   5  -6.124016E-04   1.408497E-05  -2.548877E-05   1.740838E-03  -2.267756E-06  -5.459583E-06   2.431895E-04   2.369117E-06
  b1u   6  -5.707552E-05   1.201572E-03  -2.980060E-04   1.021104E-05  -3.862075E-04  -7.292143E-04  -3.220287E-06   5.983435E-05
  b1u   7  -7.011096E-04   1.196973E-05   8.070891E-06  -1.404513E-03   3.095424E-06  -1.289357E-06  -6.233820E-04  -4.629986E-06
  b1u   8   4.355126E-05  -1.655711E-04   1.037284E-03  -7.631061E-06  -9.576219E-04   3.975904E-04   4.240182E-06  -1.711243E-04
  b1u   9   2.782536E-03  -4.246822E-05   3.063708E-05  -4.098675E-04  -1.411903E-05   1.254750E-05   1.023274E-04  -1.482726E-07
  b1u  10  -4.246822E-05   1.698401E-03   5.251147E-05   3.393815E-06  -1.102411E-03  -9.388471E-04  -2.892940E-06   1.792504E-04
  b1u  11   3.063708E-05   5.251147E-05   1.266246E-03  -8.125496E-06  -7.113814E-04  -1.466782E-04   1.091785E-06   1.004911E-04
  b1u  12  -4.098675E-04   3.393815E-06  -8.125496E-06   2.029495E-03  -1.579701E-06  -6.466911E-06   2.751718E-04   3.699126E-06
  b1u  13  -1.411903E-05  -1.102411E-03  -7.113814E-04  -1.579701E-06   2.107195E-03  -3.969078E-05  -2.706138E-06  -4.397262E-05
  b1u  14   1.254750E-05  -9.388471E-04  -1.466782E-04  -6.466911E-06  -3.969078E-05   1.260167E-03   1.936863E-06  -2.988787E-04
  b1u  15   1.023274E-04  -2.892940E-06   1.091785E-06   2.751718E-04  -2.706138E-06   1.936863E-06   6.136980E-04   2.130830E-06
  b1u  16  -1.482726E-07   1.792504E-04   1.004911E-04   3.699126E-06  -4.397262E-05  -2.988787E-04   2.130830E-06   4.727561E-04

Natural orbital populations,block 5
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.88878446     0.55002708     0.00543444     0.00468663     0.00422845     0.00379068     0.00140316     0.00069089
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00048095     0.00038577     0.00034189     0.00022792     0.00012112     0.00004670     0.00003941     0.00000806

          modens reordered block   1

               b2g   1        b2g   2        b2g   3        b2g   4        b2g   5        b2g   6        b2g   7        b2g   8
  b2g   1    1.47451       5.587171E-04  -6.451209E-04   1.880681E-03   1.190935E-06  -5.992054E-04  -3.740442E-03  -6.837704E-03
  b2g   2   5.587171E-04   1.028939E-03   3.463154E-05  -1.392172E-04  -8.831000E-04   6.066610E-05  -1.115740E-03  -1.589268E-04
  b2g   3  -6.451209E-04   3.463154E-05   1.731246E-03  -6.135908E-03   1.405117E-05   2.501695E-03   4.347335E-06  -2.062495E-05
  b2g   4   1.880681E-03  -1.392172E-04  -6.135908E-03   2.223411E-02  -3.858545E-05  -9.492187E-03   3.655132E-06   8.632913E-05
  b2g   5   1.190935E-06  -8.831000E-04   1.405117E-05  -3.858545E-05   1.298066E-03   1.454433E-05   1.344103E-03  -5.805132E-04
  b2g   6  -5.992054E-04   6.066610E-05   2.501695E-03  -9.492187E-03   1.454433E-05   4.512823E-03  -5.201915E-06  -4.329119E-05
  b2g   7  -3.740442E-03  -1.115740E-03   4.347335E-06   3.655132E-06   1.344103E-03  -5.201915E-06   1.717662E-03  -2.455460E-04
  b2g   8  -6.837704E-03  -1.589268E-04  -2.062495E-05   8.632913E-05  -5.805132E-04  -4.329119E-05  -2.455460E-04   1.388071E-03
  b2g   9  -4.592332E-04   1.205879E-05   1.497457E-03  -6.023782E-03   1.727572E-05   3.224287E-03   1.425061E-05  -2.638703E-06
  b2g  10  -6.523983E-03  -9.773899E-04  -4.604609E-05   2.072604E-04   5.281932E-04  -1.144379E-04   8.653846E-04   7.483929E-04
  b2g  11  -1.638885E-04   1.626491E-05   5.789672E-04  -2.167026E-03   2.607377E-07   8.595608E-04  -6.413963E-06  -1.098992E-05
  b2g  12   2.671466E-03  -4.123429E-04   6.488438E-06  -2.101551E-05   1.111220E-03   1.071042E-05   1.077129E-03  -1.211900E-03
  b2g  13  -2.861182E-06   7.942199E-06   2.388879E-05  -5.013462E-04  -5.649246E-06   8.084244E-04  -7.479574E-06  -3.464806E-06
  b2g  14  -3.516877E-03  -3.708401E-05  -5.340706E-06   2.433362E-05   8.544768E-05  -1.485200E-05  -3.223308E-05   2.784727E-04
  b2g  15  -3.883415E-03  -5.790877E-05  -4.119888E-06   1.945059E-05   9.557484E-05  -1.119068E-05   3.714267E-04   2.996889E-04
  b2g  16   1.130701E-04  -4.820541E-08  -2.782103E-05   6.640605E-05  -1.461910E-06   2.813014E-05  -2.193243E-06   1.211384E-06

               b2g   9        b2g  10        b2g  11        b2g  12        b2g  13        b2g  14        b2g  15        b2g  16
  b2g   1  -4.592332E-04  -6.523983E-03  -1.638885E-04   2.671466E-03  -2.861182E-06  -3.516877E-03  -3.883415E-03   1.130701E-04
  b2g   2   1.205879E-05  -9.773899E-04   1.626491E-05  -4.123429E-04   7.942199E-06  -3.708401E-05  -5.790877E-05  -4.820541E-08
  b2g   3   1.497457E-03  -4.604609E-05   5.789672E-04   6.488438E-06   2.388879E-05  -5.340706E-06  -4.119888E-06  -2.782103E-05
  b2g   4  -6.023782E-03   2.072604E-04  -2.167026E-03  -2.101551E-05  -5.013462E-04   2.433362E-05   1.945059E-05   6.640605E-05
  b2g   5   1.727572E-05   5.281932E-04   2.607377E-07   1.111220E-03  -5.649246E-06   8.544768E-05   9.557484E-05  -1.461910E-06
  b2g   6   3.224287E-03  -1.144379E-04   8.595608E-04   1.071042E-05   8.084244E-04  -1.485200E-05  -1.119068E-05   2.813014E-05
  b2g   7   1.425061E-05   8.653846E-04  -6.413963E-06   1.077129E-03  -7.479574E-06  -3.223308E-05   3.714267E-04  -2.193243E-06
  b2g   8  -2.638703E-06   7.483929E-04  -1.098992E-05  -1.211900E-03  -3.464806E-06   2.784727E-04   2.996889E-04   1.211384E-06
  b2g   9   2.613421E-03  -4.575519E-05   5.387048E-04  -3.484290E-06   1.053240E-03  -4.670397E-06  -1.119914E-05   1.439280E-04
  b2g  10  -4.575519E-05   1.576044E-03  -2.579781E-05  -9.949494E-05  -3.750251E-05   2.358306E-04  -6.839008E-05  -3.532446E-06
  b2g  11   5.387048E-04  -2.579781E-05   5.848430E-04   6.472130E-07  -2.288891E-04  -3.229196E-06  -3.552291E-06  -6.160466E-05
  b2g  12  -3.484290E-06  -9.949494E-05   6.472130E-07   1.690864E-03  -4.517400E-06  -2.523186E-04  -5.633416E-05  -2.146376E-06
  b2g  13   1.053240E-03  -3.750251E-05  -2.288891E-04  -4.517400E-06   1.143126E-03  -3.574782E-06  -1.680060E-06   2.892633E-04
  b2g  14  -4.670397E-06   2.358306E-04  -3.229196E-06  -2.523186E-04  -3.574782E-06   6.576702E-04   2.159720E-04   2.464803E-07
  b2g  15  -1.119914E-05  -6.839008E-05  -3.552291E-06  -5.633416E-05  -1.680060E-06   2.159720E-04   9.458581E-04  -1.369595E-06
  b2g  16   1.439280E-04  -3.532446E-06  -6.160466E-05  -2.146376E-06   2.892633E-04   2.464803E-07  -1.369595E-06   3.553193E-04

Natural orbital populations,block 6
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.47461025     0.03005313     0.00502627     0.00307323     0.00228459     0.00114816     0.00060293     0.00045944
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00031074     0.00022470     0.00006507     0.00006041     0.00005359     0.00001193     0.00000550     0.00000132

          modens reordered block   1

               b3g   1        b3g   2        b3g   3        b3g   4        b3g   5        b3g   6        b3g   7        b3g   8
  b3g   1    1.47211       1.644169E-03   1.093327E-05  -6.297162E-04  -3.682072E-03  -6.013475E-03  -5.751158E-03   6.974265E-04
  b3g   2   1.644169E-03   1.009611E-03   8.768934E-04   6.500941E-06   1.098061E-03   1.450510E-04   9.393489E-04  -1.975436E-06
  b3g   3   1.093327E-05   8.768934E-04   1.293167E-03   2.140067E-06   1.337796E-03  -5.774040E-04   5.239755E-04  -1.395912E-06
  b3g   4  -6.297162E-04   6.500941E-06   2.140067E-06   2.471393E-03  -3.731090E-06   2.695186E-06   7.674309E-06  -2.528625E-03
  b3g   5  -3.682072E-03   1.098061E-03   1.337796E-03  -3.731090E-06   1.708925E-03  -2.443531E-04   8.497357E-04   5.952190E-06
  b3g   6  -6.013475E-03   1.450510E-04  -5.774040E-04   2.695186E-06  -2.443531E-04   1.384792E-03   7.274567E-04  -1.865431E-06
  b3g   7  -5.751158E-03   9.393489E-04   5.239755E-04   7.674309E-06   8.497357E-04   7.274567E-04   1.521104E-03  -5.193889E-06
  b3g   8   6.974265E-04  -1.975436E-06  -1.395912E-06  -2.528625E-03   5.952190E-06  -1.865431E-06  -5.193889E-06   3.347441E-03
  b3g   9   2.262850E-03   4.207973E-04   1.110612E-03   1.730702E-06   1.079488E-03  -1.202400E-03  -7.866139E-05  -2.275666E-06
  b3g  10  -3.071838E-03   2.509951E-05   8.101811E-05  -1.492940E-06  -4.436283E-05   2.658762E-04   2.209706E-04   1.261319E-06
  b3g  11  -3.973806E-03   4.726789E-05   9.283999E-05  -5.303841E-06   3.647567E-04   2.994914E-04  -7.716236E-05   6.515574E-06

               b3g   9        b3g  10        b3g  11
  b3g   1   2.262850E-03  -3.071838E-03  -3.973806E-03
  b3g   2   4.207973E-04   2.509951E-05   4.726789E-05
  b3g   3   1.110612E-03   8.101811E-05   9.283999E-05
  b3g   4   1.730702E-06  -1.492940E-06  -5.303841E-06
  b3g   5   1.079488E-03  -4.436283E-05   3.647567E-04
  b3g   6  -1.202400E-03   2.658762E-04   2.994914E-04
  b3g   7  -7.866139E-05   2.209706E-04  -7.716236E-05
  b3g   8  -2.275666E-06   1.261319E-06   6.515574E-06
  b3g   9   1.683596E-03  -2.467056E-04  -6.024900E-05
  b3g  10  -2.467056E-04   6.449454E-04   2.094929E-04
  b3g  11  -6.024900E-05   2.094929E-04   9.544000E-04

Natural orbital populations,block 7
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.47218817     0.00547517     0.00499979     0.00301538     0.00115018     0.00060155     0.00034314     0.00022531
              MO     9       MO    10       MO    11
  occ(*)=     0.00006464     0.00005298     0.00001181

          modens reordered block   1

               au    1        au    2        au    3        au    4        au    5        au    6        au    7        au    8
  au    1   0.544585       2.575557E-03   1.585249E-03   9.675515E-04  -1.144516E-03   4.953753E-05  -3.526991E-03   1.863515E-03
  au    2   2.575557E-03   9.560910E-04   8.217778E-04  -8.441737E-04  -7.170969E-04  -2.820986E-07   6.483793E-04   1.214294E-04
  au    3   1.585249E-03   8.217778E-04   1.081472E-03  -6.088534E-04  -1.206651E-03   1.227580E-06   3.038590E-04   3.949637E-04
  au    4   9.675515E-04  -8.441737E-04  -6.088534E-04   1.603350E-03   1.763299E-04   3.060774E-06  -1.058690E-03   9.594942E-04
  au    5  -1.144516E-03  -7.170969E-04  -1.206651E-03   1.763299E-04   1.696685E-03  -2.995490E-06   5.335951E-05  -1.107534E-03
  au    6   4.953753E-05  -2.820986E-07   1.227580E-06   3.060774E-06  -2.995490E-06   9.012180E-05  -3.042950E-06   4.213296E-06
  au    7  -3.526991E-03   6.483793E-04   3.038590E-04  -1.058690E-03   5.335951E-05  -3.042950E-06   1.277666E-03  -7.194041E-04
  au    8   1.863515E-03   1.214294E-04   3.949637E-04   9.594942E-04  -1.107534E-03   4.213296E-06  -7.194041E-04   2.116499E-03
  au    9  -1.923180E-04  -2.206903E-04  -7.243118E-04   4.022146E-04   9.315328E-04  -1.393157E-06   1.468857E-04   4.254328E-05
  au   10   2.023366E-04   3.204128E-05   5.850205E-05  -1.728327E-04  -1.763122E-04   3.266799E-07  -9.918489E-05   4.164905E-05
  au   11   1.164103E-04  -3.740121E-07   1.626043E-07   5.412525E-06  -2.209771E-06   1.701669E-04  -4.236697E-06   7.230693E-06

               au    9        au   10        au   11
  au    1  -1.923180E-04   2.023366E-04   1.164103E-04
  au    2  -2.206903E-04   3.204128E-05  -3.740121E-07
  au    3  -7.243118E-04   5.850205E-05   1.626043E-07
  au    4   4.022146E-04  -1.728327E-04   5.412525E-06
  au    5   9.315328E-04  -1.763122E-04  -2.209771E-06
  au    6  -1.393157E-06   3.266799E-07   1.701669E-04
  au    7   1.468857E-04  -9.918489E-05  -4.236697E-06
  au    8   4.254328E-05   4.164905E-05   7.230693E-06
  au    9   1.256082E-03  -2.970826E-04   1.508438E-06
  au   10  -2.970826E-04   4.706558E-04  -2.304385E-07
  au   11   1.508438E-06  -2.304385E-07   5.733707E-04

Natural orbital populations,block 8
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     0.54463532     0.00430345     0.00380309     0.00140926     0.00062721     0.00047894     0.00023852     0.00012285
              MO     9       MO    10       MO    11
  occ(*)=     0.00004386     0.00003621     0.00000811


 total number of electrons =   42.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        ag  partial gross atomic populations
   ao class       1ag        2ag        3ag        4ag        5ag        6ag 
    C1_ s       0.024115   1.965539   0.526288   0.470801  -0.011174   0.063456
    C1_ p       0.003468  -0.000704   0.036422   0.130027   0.352801   0.262920
    C1_ d       0.000470   0.000179   0.005171  -0.114374   0.039269  -0.024996
    C2_ s       1.965551   0.024164   1.074875   0.234740  -0.019634   0.031049
    C2_ p       0.000801   0.007425   0.076201   0.106956   0.713795   1.298503
    C2_ d      -0.000546   0.001050   0.010556  -0.053705   0.079023   0.010937
    H1_ s       0.001179   0.000814   0.079201   0.782337   0.276120   0.267159
    H1_ p       0.000351   0.000587   0.005915   0.031430   0.000280  -0.036891
    H2_ s       0.003058   0.000734   0.160625   0.391171   0.548605   0.131907
    H2_ p       0.001552   0.000211   0.012110   0.002162   0.000459  -0.030065
 
   ao class       7ag        8ag        9ag       10ag       11ag       12ag 
    C1_ s       0.004412   0.000450   0.000466  -0.000055   0.001112   0.000317
    C1_ p       0.001555   0.000762   0.003875   0.001399   0.000093   0.000344
    C1_ d       0.000397   0.000945  -0.001162   0.000811   0.000876   0.000185
    C2_ s       0.002171   0.000943   0.000235  -0.000106   0.002212   0.000160
    C2_ p       0.004311   0.001448   0.002422   0.002821   0.000188   0.002250
    C2_ d       0.001326   0.001897  -0.000081   0.001648   0.001763   0.000775
    H1_ s       0.000187   0.002149   0.002797  -0.000018  -0.000337  -0.000203
    H1_ p      -0.000052  -0.000116   0.000174   0.000047   0.000029   0.000061
    H2_ s       0.000077   0.004283   0.001422  -0.000026  -0.000680  -0.000102
    H2_ p       0.000073  -0.000237   0.000077   0.000092   0.000059   0.000851
 
   ao class      13ag       14ag       15ag       16ag       17ag       18ag 
    C1_ s      -0.000355   0.000620  -0.000125   0.000004   0.000070  -0.000101
    C1_ p       0.001299  -0.000044  -0.000126   0.000375   0.000238   0.000011
    C1_ d       0.000695   0.000293   0.000746   0.000094   0.000202   0.000158
    C2_ s      -0.000172   0.001222  -0.000056   0.000014   0.000036  -0.000051
    C2_ p       0.000615  -0.000088   0.000297   0.000771   0.000098   0.000727
    C2_ d       0.001320   0.000574   0.000660   0.000191   0.000117  -0.000117
    H1_ s      -0.001033  -0.000272  -0.000227  -0.000089  -0.000028   0.000000
    H1_ p       0.000287   0.000001   0.000266   0.000027   0.000045   0.000027
    H2_ s      -0.000521  -0.000527  -0.000119  -0.000181  -0.000016   0.000000
    H2_ p       0.000196  -0.000002   0.000170   0.000055   0.000152   0.000117
 
   ao class      19ag       20ag       21ag       22ag       23ag       24ag 
    C1_ s       0.000093   0.000003   0.000047  -0.000126  -0.000015  -0.000017
    C1_ p       0.000175  -0.000018   0.000087   0.000025   0.000004  -0.000017
    C1_ d       0.000031   0.000039   0.000059   0.000081   0.000074   0.000026
    C2_ s       0.000012   0.000085   0.000025  -0.000281  -0.000007  -0.000027
    C2_ p       0.000041   0.000038  -0.000027   0.000059  -0.000074  -0.000033
    C2_ d       0.000066   0.000142   0.000122   0.000184   0.000167   0.000052
    H1_ s      -0.000024   0.000002   0.000093   0.000110   0.000002   0.000056
    H1_ p      -0.000036   0.000054  -0.000111  -0.000013   0.000011  -0.000009
    H2_ s      -0.000025   0.000029   0.000041   0.000235  -0.000000   0.000109
    H2_ p       0.000133   0.000087  -0.000038  -0.000031  -0.000023  -0.000018
 
   ao class      25ag       26ag       27ag       28ag       29ag       30ag 
    C1_ s      -0.000152  -0.000007   0.000023   0.000039   0.000000  -0.000001
    C1_ p       0.000065  -0.000006   0.000028  -0.000025  -0.000002  -0.000016
    C1_ d       0.000049   0.000008   0.000014   0.000005   0.000000   0.000001
    C2_ s      -0.000076  -0.000020   0.000011   0.000019   0.000003  -0.000000
    C2_ p       0.000123  -0.000013   0.000030  -0.000027  -0.000003   0.000005
    C2_ d       0.000102   0.000014  -0.000021   0.000064  -0.000001   0.000002
    H1_ s       0.000018   0.000041  -0.000022   0.000012  -0.000002   0.000029
    H1_ p      -0.000004  -0.000002   0.000019   0.000002   0.000010  -0.000012
    H2_ s       0.000007   0.000084  -0.000010   0.000006  -0.000004   0.000015
    H2_ p      -0.000028  -0.000005   0.000015  -0.000043   0.000020  -0.000003
 
   ao class      31ag       32ag       33ag       34ag       35ag       36ag 
    C1_ s      -0.000002   0.000003   0.000004  -0.000012   0.000020  -0.000001
    C1_ p       0.000006  -0.000006  -0.000001  -0.000000  -0.000008   0.000005
    C1_ d      -0.000000   0.000000   0.000002   0.000000  -0.000000  -0.000000
    C2_ s      -0.000000  -0.000003   0.000001   0.000031  -0.000002  -0.000000
    C2_ p       0.000017  -0.000002  -0.000007  -0.000001  -0.000008   0.000001
    C2_ d      -0.000001  -0.000000  -0.000010  -0.000001   0.000000  -0.000000
    H1_ s      -0.000001  -0.000002   0.000002  -0.000008   0.000002  -0.000002
    H1_ p      -0.000000   0.000012   0.000001   0.000001   0.000000   0.000000
    H2_ s      -0.000003  -0.000002   0.000000  -0.000009   0.000000  -0.000001
    H2_ p      -0.000001   0.000007   0.000015   0.000004  -0.000000   0.000000
 
   ao class      37ag       38ag       39ag 
    C1_ s      -0.000000   0.000000   0.000000
    C1_ p       0.000000  -0.000000  -0.000000
    C1_ d      -0.000000   0.000000   0.000000
    C2_ s      -0.000000  -0.000000   0.000000
    C2_ p       0.000001  -0.000000  -0.000000
    C2_ d       0.000000  -0.000000   0.000000
    H1_ s      -0.000000   0.000000   0.000000
    H1_ p      -0.000000  -0.000000   0.000000
    H2_ s      -0.000000   0.000000   0.000000
    H2_ p      -0.000000  -0.000000   0.000000

                        b3u partial gross atomic populations
   ao class       1b3u       2b3u       3b3u       4b3u       5b3u       6b3u
    C1_ s       0.022924   1.952003   0.718798   0.080749   0.078978  -0.000013
    C1_ p      -0.000455   0.001641  -0.039050   0.214387   0.727938   0.001060
    C1_ d      -0.000102   0.000126   0.042511  -0.185403  -0.018480   0.001061
    C2_ s       1.958708   0.025923   0.364486   0.169443   0.041759  -0.000034
    C2_ p       0.013886   0.014358   0.208064   0.428309   0.559094   0.008392
    C2_ d       0.003923   0.004386  -0.030070  -0.371175  -0.040116   0.002248
    H1_ s       0.000617   0.000533   0.447845   0.550984   0.433753   0.000839
    H1_ p       0.000315   0.000354   0.040406  -0.002243  -0.016559  -0.000036
    H2_ s       0.000067   0.000410   0.226390   1.098734   0.214772   0.000432
    H2_ p       0.000116   0.000266   0.005620  -0.004792  -0.004295   0.000098
 
   ao class       7b3u       8b3u       9b3u      10b3u      11b3u      12b3u
    C1_ s       0.003424   0.000812   0.000208  -0.000097  -0.000485   0.000138
    C1_ p       0.000359   0.001914   0.002711   0.001786  -0.000118  -0.000095
    C1_ d      -0.000229   0.000137  -0.001084   0.000930   0.001319   0.000344
    C2_ s       0.006772   0.000439   0.000415  -0.000065  -0.000241   0.000075
    C2_ p       0.000705   0.002840   0.005426   0.004126   0.000738  -0.000541
    C2_ d      -0.000411   0.000164  -0.002185   0.000893   0.002065   0.001642
    H1_ s       0.000621   0.002986   0.000868  -0.002166  -0.000077   0.000045
    H1_ p       0.000197  -0.000061   0.000111   0.000015   0.000118  -0.000011
    H2_ s       0.001160   0.001513   0.001763  -0.001077  -0.000039   0.000023
    H2_ p       0.000386   0.000000   0.000220   0.000386   0.000042   0.000479
 
   ao class      13b3u      14b3u      15b3u      16b3u      17b3u      18b3u
    C1_ s       0.000290   0.000478  -0.000213  -0.000006   0.000092   0.000122
    C1_ p       0.001421   0.000595  -0.000001  -0.000139   0.000274  -0.000092
    C1_ d      -0.000537   0.000295   0.000094   0.000355   0.000043   0.000062
    C2_ s       0.000582   0.000241  -0.000127   0.000030   0.000070   0.000114
    C2_ p       0.002835  -0.000047   0.000718  -0.000255   0.000161   0.000013
    C2_ d      -0.001065   0.000115   0.000262   0.000717  -0.000057   0.000175
    H1_ s      -0.000836  -0.000440   0.000015   0.000043  -0.000039   0.000023
    H1_ p       0.000330   0.000109   0.000008  -0.000021  -0.000023   0.000011
    H2_ s      -0.001682  -0.000221   0.000012   0.000075  -0.000107   0.000014
    H2_ p       0.000657   0.000119   0.000035  -0.000037   0.000030  -0.000001
 
   ao class      19b3u      20b3u      21b3u      22b3u      23b3u      24b3u
    C1_ s      -0.000202  -0.000059  -0.000302  -0.000031   0.000014  -0.000067
    C1_ p       0.000042   0.000091   0.000013   0.000009   0.000005   0.000029
    C1_ d       0.000019   0.000021   0.000153   0.000001   0.000063   0.000058
    C2_ s      -0.000101  -0.000084  -0.000164  -0.000017   0.000007  -0.000140
    C2_ p       0.000418   0.000174   0.000505   0.000318  -0.000034   0.000059
    C2_ d       0.000051   0.000042  -0.000159  -0.000259   0.000038   0.000115
    H1_ s      -0.000030   0.000056   0.000086   0.000036  -0.000012   0.000022
    H1_ p       0.000014  -0.000046  -0.000014   0.000003   0.000017  -0.000018
    H2_ s      -0.000015   0.000111   0.000044   0.000018  -0.000006   0.000042
    H2_ p       0.000126  -0.000092   0.000012   0.000054   0.000007  -0.000036
 
   ao class      25b3u      26b3u      27b3u      28b3u      29b3u      30b3u
    C1_ s       0.000014   0.000008   0.000017   0.000007  -0.000007   0.000004
    C1_ p      -0.000005  -0.000012  -0.000035  -0.000002   0.000004  -0.000010
    C1_ d       0.000008  -0.000009   0.000007   0.000000  -0.000001   0.000002
    C2_ s       0.000023   0.000004   0.000010   0.000013  -0.000004  -0.000001
    C2_ p      -0.000009  -0.000015  -0.000110  -0.000003   0.000014  -0.000008
    C2_ d       0.000016   0.000057   0.000114   0.000001  -0.000006   0.000007
    H1_ s      -0.000016   0.000013   0.000032  -0.000002  -0.000003   0.000008
    H1_ p       0.000015  -0.000002   0.000005   0.000003   0.000012  -0.000001
    H2_ s      -0.000032   0.000007   0.000016  -0.000005  -0.000001   0.000008
    H2_ p       0.000030  -0.000020  -0.000026   0.000006   0.000004  -0.000001
 
   ao class      31b3u      32b3u      33b3u      34b3u      35b3u      36b3u
    C1_ s      -0.000011   0.000029  -0.000004  -0.000000   0.000002  -0.000000
    C1_ p       0.000009   0.000000   0.000001   0.000002  -0.000002   0.000001
    C1_ d       0.000001  -0.000001  -0.000001  -0.000000  -0.000000  -0.000000
    C2_ s      -0.000001   0.000014   0.000001  -0.000000   0.000004  -0.000001
    C2_ p       0.000013  -0.000040  -0.000000   0.000002  -0.000004   0.000001
    C2_ d      -0.000007  -0.000001  -0.000001  -0.000001  -0.000000  -0.000000
    H1_ s      -0.000007   0.000001   0.000000  -0.000001   0.000000  -0.000000
    H1_ p      -0.000000   0.000000   0.000002   0.000000   0.000000  -0.000000
    H2_ s       0.000001   0.000001   0.000001  -0.000001   0.000000  -0.000000
    H2_ p       0.000008   0.000000   0.000004   0.000001   0.000000  -0.000000
 
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
    C1_ p       0.016768   0.145249   0.124929   0.654325   0.005290   0.001237
    C1_ d       0.005055  -0.035967  -0.020938   0.018001   0.001144   0.000082
    C2_ s       1.965412   1.098983   0.119960   0.000017  -0.000047   0.001242
    C2_ p       0.007796   0.013141   1.163653   1.310508   0.004187   0.003502
    C2_ d       0.001844   0.046079  -0.037905   0.035723   0.002166   0.000219
    H1_ p       0.000103  -0.010146   0.002515  -0.014634   0.000078   0.000018
    H2_ s       0.001851   0.671913   0.647883   0.000004   0.001272   0.004511
    H2_ p       0.001171   0.055760  -0.023218  -0.029242  -0.000017  -0.000080
 
   ao class       7b2u       8b2u       9b2u      10b2u      11b2u      12b2u
    C1_ p       0.000817   0.002168   0.000530  -0.000332  -0.000237   0.000353
    C1_ d       0.000861   0.000291   0.000941   0.000982  -0.000020  -0.000031
    C2_ s       0.000000  -0.000185  -0.000727   0.000212   0.000728  -0.000000
    C2_ p       0.001630   0.003759   0.000089  -0.000303   0.000778   0.000683
    C2_ d       0.001731   0.001535   0.002452   0.001007   0.000430  -0.000056
    H1_ p       0.000471   0.000252  -0.000012   0.000325   0.000043   0.000053
    H2_ s       0.000000  -0.003240  -0.000118   0.000068  -0.000664   0.000000
    H2_ p       0.000943   0.000147   0.000171   0.000144   0.000185   0.000105
 
   ao class      13b2u      14b2u      15b2u      16b2u      17b2u      18b2u
    C1_ p       0.000480   0.000019  -0.000108   0.000265  -0.000022   0.000329
    C1_ d       0.000145   0.000125   0.000060   0.000029   0.000148  -0.000163
    C2_ s      -0.000333   0.000000   0.000241  -0.000306  -0.000000  -0.000460
    C2_ p       0.000230   0.000037  -0.000032   0.000196  -0.000039   0.000188
    C2_ d       0.000208   0.000250   0.000209   0.000043   0.000290   0.000154
    H1_ p       0.000021   0.000081   0.000021   0.000080  -0.000058   0.000014
    H2_ s       0.000027   0.000000   0.000059  -0.000045   0.000000   0.000129
    H2_ p       0.000020   0.000163  -0.000009   0.000060  -0.000115  -0.000015
 
   ao class      19b2u      20b2u      21b2u      22b2u      23b2u      24b2u
    C1_ p       0.000213  -0.000025  -0.000006  -0.000062  -0.000001   0.000008
    C1_ d      -0.000172   0.000004   0.000042   0.000073  -0.000007  -0.000003
    C2_ s      -0.000056   0.000023   0.000012   0.000027  -0.000000  -0.000011
    C2_ p       0.000119  -0.000005  -0.000022  -0.000084  -0.000003   0.000010
    C2_ d      -0.000085   0.000095   0.000008   0.000048  -0.000013  -0.000004
    H1_ p       0.000034  -0.000001  -0.000013  -0.000018   0.000015  -0.000001
    H2_ s       0.000056  -0.000017   0.000020   0.000049   0.000000  -0.000004
    H2_ p       0.000021   0.000025  -0.000009  -0.000002   0.000030   0.000017
 
   ao class      25b2u      26b2u      27b2u      28b2u      29b2u      30b2u
    C1_ p       0.000010  -0.000022  -0.000003   0.000001   0.000000   0.000000
    C1_ d      -0.000005  -0.000000   0.000000  -0.000001   0.000000  -0.000000
    C2_ s      -0.000013   0.000039   0.000003  -0.000001  -0.000000  -0.000000
    C2_ p       0.000016  -0.000014   0.000004   0.000003   0.000000  -0.000000
    C2_ d      -0.000002  -0.000001  -0.000000  -0.000000  -0.000000  -0.000000
    H1_ p       0.000006   0.000000  -0.000000   0.000001  -0.000000   0.000000
    H2_ s      -0.000009   0.000002   0.000000  -0.000002   0.000000   0.000000
    H2_ p       0.000002   0.000000  -0.000000   0.000000   0.000000   0.000000

                        b1g partial gross atomic populations
   ao class       1b1g       2b1g       3b1g       4b1g       5b1g       6b1g
    C1_ p       0.003407   0.022569   0.773911   0.002372   0.000321   0.002687
    C1_ d       0.000348   0.001935   0.015655   0.000743   0.000337   0.000471
    C2_ s       1.989749   0.710104   0.097070   0.006552   0.000702   0.000000
    C2_ p       0.004512   0.209261   0.787551   0.003521   0.005967   0.005370
    C2_ d       0.001538  -0.170588  -0.029149   0.000958  -0.001595   0.000978
    H1_ p       0.000044  -0.008947  -0.007969   0.000067  -0.000006   0.000063
    H2_ s       0.000202   1.174586   0.395773   0.000257   0.004233   0.000000
    H2_ p       0.000200   0.042631  -0.058885  -0.000041   0.000254   0.000126
 
   ao class       7b1g       8b1g       9b1g      10b1g      11b1g      12b1g
    C1_ p       0.001383  -0.000014  -0.000253   0.000241  -0.000018   0.000483
    C1_ d       0.000454   0.000649   0.000678   0.000189   0.000015  -0.000133
    C2_ s       0.000478  -0.000525   0.000000  -0.000187   0.000106  -0.000150
    C2_ p       0.001207   0.001929  -0.000508  -0.000065   0.000349   0.000254
    C2_ d       0.000500   0.001360   0.001358   0.001214   0.000307   0.000175
    H1_ p       0.000547   0.000037   0.000158   0.000024   0.000086   0.000068
    H2_ s      -0.000301  -0.001559  -0.000000  -0.000342  -0.000045   0.000000
    H2_ p       0.000364   0.000447   0.000317   0.000410   0.000113   0.000075
 
   ao class      13b1g      14b1g      15b1g      16b1g      17b1g      18b1g
    C1_ p      -0.000006  -0.000047   0.000152  -0.000052   0.000060   0.000011
    C1_ d       0.000056   0.000062  -0.000145   0.000090   0.000055  -0.000018
    C2_ s       0.000105   0.000071   0.000000  -0.000022  -0.000228   0.000034
    C2_ p       0.000223   0.000107   0.000306  -0.000018   0.000129   0.000048
    C2_ d       0.000041   0.000119  -0.000291   0.000152   0.000096   0.000010
    H1_ p       0.000094   0.000010   0.000064  -0.000019  -0.000018   0.000004
    H2_ s      -0.000050   0.000134   0.000000   0.000002   0.000027  -0.000032
    H2_ p       0.000003  -0.000159   0.000129   0.000007  -0.000014   0.000031
 
   ao class      19b1g      20b1g      21b1g      22b1g      23b1g      24b1g
    C1_ p       0.000185  -0.000009   0.000008  -0.000065  -0.000000  -0.000004
    C1_ d      -0.000189   0.000041   0.000001   0.000092  -0.000000  -0.000008
    C2_ s       0.000000   0.000058  -0.000001   0.000000   0.000000   0.000004
    C2_ p       0.000365  -0.000044  -0.000018  -0.000131  -0.000007  -0.000004
    C2_ d      -0.000370   0.000028   0.000002   0.000185   0.000000  -0.000001
    H1_ p       0.000026  -0.000029   0.000002  -0.000021   0.000002   0.000009
    H2_ s      -0.000000   0.000018   0.000044   0.000000  -0.000004   0.000003
    H2_ p       0.000050  -0.000012  -0.000017  -0.000041   0.000018   0.000007
 
   ao class      25b1g      26b1g      27b1g      28b1g      29b1g      30b1g
    C1_ p      -0.000002  -0.000001   0.000001   0.000000   0.000000   0.000000
    C1_ d       0.000000  -0.000000  -0.000001   0.000000   0.000000   0.000000
    C2_ s       0.000019  -0.000002  -0.000000  -0.000001   0.000000   0.000000
    C2_ p      -0.000016   0.000007   0.000001   0.000001  -0.000000   0.000000
    C2_ d      -0.000000  -0.000000  -0.000003  -0.000000   0.000000   0.000000
    H1_ p      -0.000000  -0.000000   0.000001  -0.000000   0.000000  -0.000000
    H2_ s       0.000003  -0.000003  -0.000000  -0.000000   0.000000   0.000000
    H2_ p       0.000000   0.000000   0.000002  -0.000000   0.000000  -0.000000

                        b1u partial gross atomic populations
   ao class       1b1u       2b1u       3b1u       4b1u       5b1u       6b1u
    C1_ p       0.613143   0.335261   0.000031  -0.000495   0.002606   0.000790
    C1_ d       0.009569  -0.001816   0.001330   0.000565   0.000539   0.000210
    C2_ p       1.231643   0.171175   0.000583   0.004186  -0.000442   0.000116
    C2_ d       0.019925   0.030037   0.003002   0.000248   0.000927   0.002619
    H1_ p       0.004811   0.010177   0.000134   0.000043   0.000470   0.000032
    H2_ p       0.009695   0.005192   0.000355   0.000139   0.000128   0.000023
 
   ao class       7b1u       8b1u       9b1u      10b1u      11b1u      12b1u
    C1_ p       0.000166   0.000011   0.000001  -0.000100   0.000002   0.000144
    C1_ d       0.000364   0.000015   0.000082  -0.000054   0.000153   0.000091
    C2_ p       0.000107   0.000045   0.000014   0.000454   0.000017  -0.000021
    C2_ d       0.000558   0.000042   0.000115   0.000043   0.000276   0.000030
    H1_ p       0.000130   0.000191   0.000187   0.000014  -0.000036  -0.000013
    H2_ p       0.000078   0.000387   0.000081   0.000028  -0.000070  -0.000004
 
   ao class      13b1u      14b1u      15b1u      16b1u
    C1_ p       0.000022   0.000001   0.000016   0.000002
    C1_ d       0.000024  -0.000006   0.000020  -0.000002
    C2_ p      -0.000008   0.000003  -0.000007  -0.000000
    C2_ d       0.000082   0.000014  -0.000007  -0.000001
    H1_ p       0.000001   0.000019  -0.000003   0.000007
    H2_ p      -0.000000   0.000015   0.000020   0.000003

                        b2g partial gross atomic populations
   ao class       1b2g       2b2g       3b2g       4b2g       5b2g       6b2g
    C1_ p       0.956339   0.010253   0.001018   0.001383   0.000197   0.000025
    C1_ d       0.001337  -0.000499   0.001747   0.000597   0.000267   0.000042
    C2_ p       0.474551   0.020548   0.000495   0.000650   0.000419   0.000014
    C2_ d       0.022231  -0.001029   0.001063   0.000349   0.000519   0.001061
    H1_ p       0.013403   0.000252   0.000474   0.000063   0.000294   0.000004
    H2_ p       0.006750   0.000527   0.000229   0.000033   0.000588   0.000002
 
   ao class       7b2g       8b2g       9b2g      10b2g      11b2g      12b2g
    C1_ p      -0.000002   0.000253   0.000019   0.000024   0.000007  -0.000024
    C1_ d       0.000052  -0.000133   0.000036   0.000153   0.000008   0.000058
    C2_ p      -0.000001   0.000501   0.000035   0.000016   0.000004  -0.000038
    C2_ d       0.000067  -0.000272   0.000063   0.000080   0.000045   0.000118
    H1_ p       0.000326   0.000036   0.000051  -0.000031  -0.000000  -0.000017
    H2_ p       0.000162   0.000075   0.000106  -0.000017   0.000001  -0.000036
 
   ao class      13b2g      14b2g      15b2g      16b2g
    C1_ p       0.000020   0.000002  -0.000000   0.000000
    C1_ d       0.000001  -0.000002  -0.000001   0.000000
    C2_ p       0.000006   0.000001  -0.000000   0.000001
    C2_ d       0.000018  -0.000001  -0.000002   0.000000
    H1_ p       0.000006   0.000008   0.000003  -0.000000
    H2_ p       0.000003   0.000004   0.000007  -0.000000

                        b3g partial gross atomic populations
   ao class       1b3g       2b3g       3b3g       4b3g       5b3g       6b3g
    C1_ d       0.013511   0.001832   0.000127   0.000025   0.000692   0.000028
    C2_ p       1.429654   0.000000   0.001435   0.002036   0.000038  -0.000003
    C2_ d       0.009113   0.003643   0.002733   0.000862   0.000413   0.000091
    H2_ p       0.019911   0.000000   0.000705   0.000092   0.000006   0.000486
 
   ao class       7b3g       8b3g       9b3g      10b3g      11b3g
    C1_ d       0.000115  -0.000002   0.000032   0.000008  -0.000000
    C2_ p       0.000000   0.000040   0.000009   0.000027   0.000003
    C2_ d       0.000229   0.000236   0.000024   0.000008  -0.000003
    H2_ p      -0.000000  -0.000049   0.000000   0.000010   0.000012

                        au  partial gross atomic populations
   ao class       1au        2au        3au        4au        5au        6au 
    C1_ d       0.022189   0.000286   0.001717   0.000259   0.000208   0.000054
    C2_ p       0.498872   0.002128   0.000809   0.000264  -0.000000   0.000008
    C2_ d       0.007838   0.001237   0.001207   0.000673   0.000419   0.000148
    H2_ p       0.015737   0.000652   0.000070   0.000213  -0.000000   0.000269
 
   ao class       7au        8au        9au       10au       11au 
    C1_ d       0.000009   0.000053   0.000001   0.000012  -0.000000
    C2_ p       0.000153   0.000013   0.000011   0.000000   0.000001
    C2_ d       0.000098   0.000054   0.000037   0.000024  -0.000003
    H2_ p      -0.000022   0.000003  -0.000004  -0.000000   0.000010


                        gross atomic populations
     ao           C1_        C2_        H1_        H2_
      s         5.903355  11.873974   2.845842   5.679302
      p         5.389896  10.840647   0.019980   0.042407
      d        -0.199934  -0.395468   0.000000   0.000000
    total      11.093317  22.319153   2.865822   5.721708
 

 Total number of electrons:   42.00000000

 item #                     2 suffix=:.drt1.state2:
 read_civout: repnuc=  -224.450660864067     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref# 91% root-following 0
 MR-CISD energy:  -231.23855511    -6.78789424
 residuum:     0.00022506
 deltae:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  1max overlap with ref# 85% root-following 0
 MR-CISD energy:  -231.15249866    -6.70183780
 residuum:     0.00071340
 deltae:     0.00000048
 apxde:     0.00000018

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.90511482     0.02146402    -0.00500284    -0.01135164     0.04821206    -0.03057590    -0.03077740     0.00000000
 ref:   2     0.02122424     0.85499564    -0.30428239    -0.00235292    -0.01622481    -0.05281919    -0.04073863     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -0.90511482     0.02146402    -0.00500284    -0.01135164     0.04821206     0.00000000     0.00000000     0.00000000
 ref:   2     0.02122424     0.85499564    -0.30428239    -0.00235292    -0.01622481     0.00000000     0.00000000     0.00000000

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
   DYZ=   51682  DYX=   70131  DYW=   73894
   D0Z=   12972  D0Y=  119645  D0X=   15912  D0W=   17672
  DDZI=   27020 DDYI=  198402 DDXI=   27696 DDWI=   30254
  DDZE=       0 DDYE=   29922 DDXE=    4857 DDWE=    5335
================================================================================
Trace of MO density:    30.000000
   30  correlated and    12  frozen core electrons

          modens reordered block   1

               ag    1        ag    2        ag    3        ag    4        ag    5        ag    6        ag    7        ag    8
  ag    1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  ag    3    0.00000        0.00000        1.98725      -5.156907E-05   1.102022E-03   1.590330E-04  -3.529553E-05   1.706724E-04
  ag    4    0.00000        0.00000      -5.156907E-05    1.98064      -1.474311E-04  -2.183645E-03  -2.658749E-04   6.661191E-04
  ag    5    0.00000        0.00000       1.102022E-03  -1.474311E-04    1.97925      -1.499208E-04  -4.307557E-03   3.772401E-04
  ag    6    0.00000        0.00000       1.590330E-04  -2.183645E-03  -1.499208E-04    1.97384       1.323983E-04  -1.597919E-03
  ag    7    0.00000        0.00000      -3.529553E-05  -2.658749E-04  -4.307557E-03   1.323983E-04   8.899320E-04  -3.242703E-05
  ag    8    0.00000        0.00000       1.706724E-04   6.661191E-04   3.772401E-04  -1.597919E-03  -3.242703E-05   2.210037E-04
  ag    9    0.00000        0.00000       1.471842E-03  -1.333830E-03  -7.344389E-03   7.164333E-04   2.211289E-04   3.995883E-05
  ag   10    0.00000        0.00000      -5.632491E-04  -1.046174E-03  -6.335871E-04   2.210763E-03   5.854286E-05  -3.991937E-04
  ag   11    0.00000        0.00000       5.103802E-04   9.304835E-04   6.170016E-03   1.957408E-04  -1.354599E-03   3.675608E-05
  ag   12    0.00000        0.00000      -1.203774E-04  -1.494583E-03  -8.468582E-04   3.234977E-03   1.242398E-05  -3.553358E-04
  ag   13    0.00000        0.00000      -6.661317E-04  -7.625051E-04   3.561298E-04  -7.487175E-05   1.095505E-04  -1.999186E-04
  ag   14    0.00000        0.00000      -2.816547E-03   1.553299E-03   1.218107E-02   5.380629E-04  -1.227564E-03   8.219935E-07
  ag   15    0.00000        0.00000      -1.296662E-03  -1.262389E-04   2.577427E-04   3.555527E-03  -1.010613E-04  -5.967974E-04
  ag   16    0.00000        0.00000       1.918770E-04  -3.343899E-03  -2.208861E-03   4.751786E-03   1.294331E-05  -5.713167E-04
  ag   17    0.00000        0.00000       1.466291E-03   1.954139E-03   3.345499E-03   1.033759E-03  -7.401764E-04   6.098009E-06
  ag   18    0.00000        0.00000       1.384759E-03  -1.460202E-03  -1.942486E-03  -4.037158E-03   1.193107E-04   3.849602E-04
  ag   19    0.00000        0.00000      -4.913328E-04  -4.079798E-04  -2.791800E-03   1.245027E-03   3.976535E-04   2.272676E-06
  ag   20    0.00000        0.00000       3.342993E-03  -7.495377E-04  -7.963751E-03  -2.397172E-03   4.651190E-04  -3.903416E-06
  ag   21    0.00000        0.00000      -1.028966E-03  -1.176591E-03   4.428919E-04  -3.168021E-04   1.236382E-04  -1.694536E-04
  ag   22    0.00000        0.00000      -2.653477E-04   2.349855E-03   1.521484E-03  -4.285752E-03   3.727601E-05   3.756888E-04
  ag   23    0.00000        0.00000       5.593211E-04   3.246476E-03   1.262281E-03  -1.151174E-03  -6.115407E-05   5.022563E-04
  ag   24    0.00000        0.00000       1.077090E-03  -5.711216E-04  -1.443546E-03  -3.452156E-03   1.124308E-05   1.232965E-04
  ag   25    0.00000        0.00000      -1.438356E-03  -7.464601E-06  -3.473550E-03   1.740149E-03   7.421420E-04  -1.846062E-06
  ag   26    0.00000        0.00000      -2.678178E-05   2.775361E-03   5.917072E-04   1.314346E-03  -1.391594E-05   1.279821E-04
  ag   27    0.00000        0.00000      -8.712053E-04   8.918772E-04   2.726138E-03   8.898963E-04  -4.683156E-04   1.301369E-06
  ag   28    0.00000        0.00000      -4.438300E-04   8.758795E-04   2.097681E-04  -2.685376E-03   1.328962E-04   1.127060E-05
  ag   29    0.00000        0.00000       6.640570E-04   1.180004E-03   8.244178E-04  -8.308987E-04  -3.383381E-04   3.185251E-06
  ag   30    0.00000        0.00000       1.131513E-04   1.101663E-04  -3.451326E-04   4.133862E-03   2.862906E-06  -1.252601E-04
  ag   31    0.00000        0.00000       2.801261E-03   1.067890E-04   5.502507E-04  -6.888650E-04  -2.053071E-04  -4.952300E-06
  ag   32    0.00000        0.00000      -1.799614E-03   2.094550E-04   4.575093E-03   1.252885E-03   1.299136E-04   5.331950E-06
  ag   33    0.00000        0.00000       1.317997E-05  -2.057261E-03  -4.464757E-06  -3.066560E-03  -7.054075E-06   1.478887E-04
  ag   34    0.00000        0.00000       7.536583E-04  -7.513279E-04  -1.000679E-04  -6.193736E-04  -4.451830E-05   1.562409E-04
  ag   35    0.00000        0.00000       3.919896E-04   2.332355E-03   3.736882E-05  -1.612397E-03  -1.538306E-05   8.038940E-06
  ag   36    0.00000        0.00000       2.843545E-04   1.325883E-03   3.452140E-04  -3.512356E-03  -1.852292E-05   2.233994E-05
  ag   37    0.00000        0.00000      -1.113850E-03  -1.063464E-05  -3.139040E-03   5.603431E-05   1.690587E-04   6.934053E-07
  ag   38    0.00000        0.00000       2.000366E-05  -3.175031E-03  -2.200978E-04  -2.808212E-03   1.918827E-06  -1.011318E-04
  ag   39    0.00000        0.00000       1.669002E-04   8.082754E-04  -9.677901E-05   8.785858E-04   4.935589E-06   5.344907E-06

               ag    9        ag   10        ag   11        ag   12        ag   13        ag   14        ag   15        ag   16
  ag    3   1.471842E-03  -5.632491E-04   5.103802E-04  -1.203774E-04  -6.661317E-04  -2.816547E-03  -1.296662E-03   1.918770E-04
  ag    4  -1.333830E-03  -1.046174E-03   9.304835E-04  -1.494583E-03  -7.625051E-04   1.553299E-03  -1.262389E-04  -3.343899E-03
  ag    5  -7.344389E-03  -6.335871E-04   6.170016E-03  -8.468582E-04   3.561298E-04   1.218107E-02   2.577427E-04  -2.208861E-03
  ag    6   7.164333E-04   2.210763E-03   1.957408E-04   3.234977E-03  -7.487175E-05   5.380629E-04   3.555527E-03   4.751786E-03
  ag    7   2.211289E-04   5.854286E-05  -1.354599E-03   1.242398E-05   1.095505E-04  -1.227564E-03  -1.010613E-04   1.294331E-05
  ag    8   3.995883E-05  -3.991937E-04   3.675608E-05  -3.553358E-04  -1.999186E-04   8.219935E-07  -5.967974E-04  -5.713167E-04
  ag    9   2.530921E-03  -2.501229E-05   1.839854E-05  -1.476706E-05  -1.539538E-04  -2.856843E-03   5.831927E-05  -5.666709E-07
  ag   10  -2.501229E-05   8.584829E-04  -7.389350E-05   5.263960E-04   5.826326E-04  -7.185045E-05   1.589612E-03   6.039126E-04
  ag   11   1.839854E-05  -7.389350E-05   2.321446E-03  -2.479224E-06  -1.383429E-04   1.747250E-03   1.375023E-04  -2.073041E-07
  ag   12  -1.476706E-05   5.263960E-04  -2.479224E-06   7.337934E-04   4.136184E-05   1.458884E-05   6.427709E-04   1.508096E-03
  ag   13  -1.539538E-04   5.826326E-04  -1.383429E-04   4.136184E-05   7.541431E-04  -1.241740E-05   1.119968E-03  -4.946716E-04
  ag   14  -2.856843E-03  -7.185045E-05   1.747250E-03   1.458884E-05  -1.241740E-05   5.757773E-03   1.452363E-05   3.497951E-05
  ag   15   5.831927E-05   1.589612E-03   1.375023E-04   6.427709E-04   1.119968E-03   1.452363E-05   4.094931E-03   1.290923E-04
  ag   16  -5.666709E-07   6.039126E-04  -2.073041E-07   1.508096E-03  -4.946716E-04   3.497951E-05   1.290923E-04   3.867540E-03
  ag   17   6.080080E-04  -3.830866E-05   1.819223E-03   6.347975E-06  -1.308993E-06   6.750230E-04  -5.307713E-05   4.133505E-06
  ag   18  -8.115123E-05  -1.395016E-03  -1.489300E-04  -2.241977E-04  -1.095793E-03  -1.567463E-05  -4.979021E-03   1.146614E-03
  ag   19   8.895262E-04  -1.225465E-07  -3.201034E-04   7.278985E-07  -6.996194E-06  -6.477884E-04  -6.108124E-05   1.765163E-05
  ag   20   1.452548E-03   5.605748E-05  -8.930446E-04  -1.842407E-05   3.370346E-05  -4.257578E-03   1.483992E-05  -5.610773E-05
  ag   21  -1.580483E-04   5.303668E-04  -1.497834E-04  -8.160609E-05   9.233763E-04  -4.395014E-05   8.373742E-04  -9.396824E-04
  ag   22  -6.533359E-05  -4.431888E-04  -5.571717E-05  -1.037693E-03   4.351266E-04  -1.859015E-05  -4.313541E-04  -2.778340E-03
  ag   23   5.153749E-05  -8.525140E-04   6.411299E-05  -1.058048E-03  -1.663822E-04   2.574606E-05  -1.097233E-03  -2.593538E-03
  ag   24   4.904950E-05  -4.940754E-04  -6.466111E-06  -1.013650E-04  -4.037101E-04  -3.600879E-05  -2.366690E-03   7.094422E-04
  ag   25   1.588905E-03   1.601674E-05  -6.702087E-04   6.423810E-06  -9.135629E-07  -1.538918E-03  -5.759556E-05   1.378441E-05
  ag   26  -5.662685E-06  -1.831656E-04   1.713117E-05  -2.893330E-04  -4.055374E-05   3.040031E-05   9.764795E-05  -9.316692E-04
  ag   27  -4.924318E-04  -2.594478E-05   9.704345E-04  -6.896347E-07   1.811973E-05   1.818074E-03  -3.145361E-05  -2.839883E-06
  ag   28  -1.705477E-04   7.439022E-05  -2.707599E-04  -2.766040E-04   5.137034E-04   1.064772E-04  -1.271694E-04  -1.101335E-03
  ag   29   4.459729E-04   7.488454E-06   8.399368E-04  -7.370577E-05   1.578206E-04  -5.155862E-04  -8.735220E-05  -3.043461E-04
  ag   30   1.007671E-05   1.333342E-04  -1.228978E-06   2.399130E-04   4.091889E-05  -1.754258E-05  -6.949345E-05   3.296123E-04
  ag   31  -6.308471E-04   2.114709E-06   1.542239E-04   1.365653E-05  -1.337898E-05   8.994386E-04   5.544672E-05   2.668649E-05
  ag   32   7.922993E-05  -1.130457E-05  -1.502186E-04  -1.009183E-06  -1.790205E-05   2.090903E-04  -1.889951E-05   7.298567E-06
  ag   33   7.018083E-06  -4.153688E-04  -1.652659E-06  -1.006961E-04  -3.878263E-04   8.183349E-06  -7.348953E-04  -3.778059E-05
  ag   34   5.260795E-05  -3.240397E-04   5.056350E-05  -3.031787E-04  -2.326630E-04   2.384348E-05  -6.716305E-04  -5.206786E-04
  ag   35   1.003808E-05  -5.082026E-05   1.984364E-05  -9.081473E-06  -4.496163E-05   1.230383E-05  -2.062052E-04   1.595991E-05
  ag   36  -3.557993E-05  -1.008526E-04   2.405254E-05  -1.809257E-05  -9.007619E-05   6.754587E-05  -3.226636E-04  -4.444359E-05
  ag   37   3.727765E-04  -1.860507E-06  -2.414737E-04   4.055031E-06  -1.568442E-05  -7.028001E-04  -4.028967E-05   1.293143E-05
  ag   38  -3.772148E-06   1.814783E-04  -2.488921E-06   2.084300E-04   4.318852E-05  -5.591433E-06   3.079372E-04   4.607399E-04
  ag   39   5.951132E-06  -3.805608E-05  -6.304178E-06   5.144059E-05  -8.715572E-05  -9.377566E-06  -1.069925E-04   2.570716E-04

               ag   17        ag   18        ag   19        ag   20        ag   21        ag   22        ag   23        ag   24
  ag    3   1.466291E-03   1.384759E-03  -4.913328E-04   3.342993E-03  -1.028966E-03  -2.653477E-04   5.593211E-04   1.077090E-03
  ag    4   1.954139E-03  -1.460202E-03  -4.079798E-04  -7.495377E-04  -1.176591E-03   2.349855E-03   3.246476E-03  -5.711216E-04
  ag    5   3.345499E-03  -1.942486E-03  -2.791800E-03  -7.963751E-03   4.428919E-04   1.521484E-03   1.262281E-03  -1.443546E-03
  ag    6   1.033759E-03  -4.037158E-03   1.245027E-03  -2.397172E-03  -3.168021E-04  -4.285752E-03  -1.151174E-03  -3.452156E-03
  ag    7  -7.401764E-04   1.193107E-04   3.976535E-04   4.651190E-04   1.236382E-04   3.727601E-05  -6.115407E-05   1.124308E-05
  ag    8   6.098009E-06   3.849602E-04   2.272676E-06  -3.903416E-06  -1.694536E-04   3.756888E-04   5.022563E-04   1.232965E-04
  ag    9   6.080080E-04  -8.115123E-05   8.895262E-04   1.452548E-03  -1.580483E-04  -6.533359E-05   5.153749E-05   4.904950E-05
  ag   10  -3.830866E-05  -1.395016E-03  -1.225465E-07   5.605748E-05   5.303668E-04  -4.431888E-04  -8.525140E-04  -4.940754E-04
  ag   11   1.819223E-03  -1.489300E-04  -3.201034E-04  -8.930446E-04  -1.497834E-04  -5.571717E-05   6.411299E-05  -6.466111E-06
  ag   12   6.347975E-06  -2.241977E-04   7.278985E-07  -1.842407E-05  -8.160609E-05  -1.037693E-03  -1.058048E-03  -1.013650E-04
  ag   13  -1.308993E-06  -1.095793E-03  -6.996194E-06   3.370346E-05   9.233763E-04   4.351266E-04  -1.663822E-04  -4.037101E-04
  ag   14   6.750230E-04  -1.567463E-05  -6.477884E-04  -4.257578E-03  -4.395014E-05  -1.859015E-05   2.574606E-05  -3.600879E-05
  ag   15  -5.307713E-05  -4.979021E-03  -6.108124E-05   1.483992E-05   8.373742E-04  -4.313541E-04  -1.097233E-03  -2.366690E-03
  ag   16   4.133505E-06   1.146614E-03   1.765163E-05  -5.610773E-05  -9.396824E-04  -2.778340E-03  -2.593538E-03   7.094422E-04
  ag   17   2.654696E-03   8.949626E-05   4.352552E-04  -7.735072E-04   3.037415E-05   1.808468E-05  -1.476844E-05   4.300091E-05
  ag   18   8.949626E-05   7.880368E-03   7.737366E-05  -1.262440E-05  -7.069033E-04  -2.875820E-04   4.838154E-05   4.890608E-03
  ag   19   4.352552E-04   7.737366E-05   1.079878E-03  -5.244516E-04  -5.746741E-07  -7.416251E-06  -1.821777E-05   4.970061E-05
  ag   20  -7.735072E-04  -1.262440E-05  -5.244516E-04   4.889457E-03   5.446551E-05   3.835278E-05   4.765751E-06   1.311131E-06
  ag   21   3.037415E-05  -7.069033E-04  -5.746741E-07   5.446551E-05   1.432176E-03   9.778164E-04   2.170926E-04  -3.340512E-04
  ag   22   1.808468E-05  -2.875820E-04  -7.416251E-06   3.835278E-05   9.778164E-04   2.253751E-03   2.035933E-03  -3.788112E-04
  ag   23  -1.476844E-05   4.838154E-05  -1.821777E-05   4.765751E-06   2.170926E-04   2.035933E-03   2.677426E-03  -4.885352E-04
  ag   24   4.300091E-05   4.890608E-03   4.970061E-05   1.311131E-06  -3.340512E-04  -3.788112E-04  -4.885352E-04   4.177316E-03
  ag   25   5.713205E-04   4.037336E-05   1.220741E-03  -8.735657E-04   1.052251E-05  -5.367816E-06  -2.656525E-05   2.750519E-05
  ag   26   6.210743E-06  -9.457400E-04  -6.764483E-06  -2.448037E-05   3.812155E-05   6.940177E-04   1.008892E-03  -1.040852E-03
  ag   27   1.228179E-03   6.092376E-05   3.123654E-04  -1.738485E-03   3.001372E-05   2.010459E-05  -3.643902E-06   1.982809E-05
  ag   28  -3.214095E-04   4.248035E-04   4.715071E-05  -1.899907E-04   1.096262E-03   1.223824E-03   7.665939E-04   3.422085E-04
  ag   29   1.342202E-03   1.958546E-04  -1.574248E-04   8.295549E-04   3.415915E-04   3.569522E-04   1.916831E-04   1.243483E-04
  ag   30  -4.643197E-06   2.816074E-04  -8.912957E-07   2.028516E-05   2.384913E-04  -5.492199E-05   1.874509E-04   1.984513E-06
  ag   31  -4.537695E-04  -6.506098E-05  -7.103501E-04  -3.903123E-04  -3.538748E-05  -3.461813E-05  -4.544847E-06  -2.456045E-05
  ag   32   5.236495E-05   1.723340E-05   3.174658E-04  -1.260780E-03  -1.432261E-05  -9.441496E-06  -2.294665E-06   1.472565E-05
  ag   33  -1.807034E-05   4.258310E-05  -1.096039E-05   5.321690E-06  -3.835262E-04   1.213613E-04   6.699640E-04  -9.728755E-04
  ag   34  -2.369756E-05   9.543333E-04  -7.436637E-06  -1.518190E-05  -3.365204E-04   2.933966E-04   7.148176E-04   1.050313E-03
  ag   35  -6.223927E-08   4.544835E-04  -2.358265E-06  -2.829820E-06   5.240842E-06   1.126460E-04   1.673625E-04   5.788214E-04
  ag   36   4.003984E-06   3.823081E-04  -1.811508E-05  -1.367588E-05  -2.196507E-05   1.689442E-04   2.979936E-04   4.716934E-05
  ag   37  -6.504320E-05   5.456845E-05   1.666971E-04   2.477172E-04  -6.326762E-06   1.932927E-06   1.083929E-05   2.634350E-05
  ag   38  -6.853550E-06  -2.389578E-04  -5.708747E-06   8.775598E-06  -3.074281E-05  -3.470423E-04  -4.192888E-04  -1.286243E-04
  ag   39   3.030976E-06   2.022715E-04   8.287333E-06   2.174427E-06  -1.385385E-04  -1.961223E-04  -2.687862E-04   1.760715E-04

               ag   25        ag   26        ag   27        ag   28        ag   29        ag   30        ag   31        ag   32
  ag    3  -1.438356E-03  -2.678178E-05  -8.712053E-04  -4.438300E-04   6.640570E-04   1.131513E-04   2.801261E-03  -1.799614E-03
  ag    4  -7.464601E-06   2.775361E-03   8.918772E-04   8.758795E-04   1.180004E-03   1.101663E-04   1.067890E-04   2.094550E-04
  ag    5  -3.473550E-03   5.917072E-04   2.726138E-03   2.097681E-04   8.244178E-04  -3.451326E-04   5.502507E-04   4.575093E-03
  ag    6   1.740149E-03   1.314346E-03   8.898963E-04  -2.685376E-03  -8.308987E-04   4.133862E-03  -6.888650E-04   1.252885E-03
  ag    7   7.421420E-04  -1.391594E-05  -4.683156E-04   1.328962E-04  -3.383381E-04   2.862906E-06  -2.053071E-04   1.299136E-04
  ag    8  -1.846062E-06   1.279821E-04   1.301369E-06   1.127060E-05   3.185251E-06  -1.252601E-04  -4.952300E-06   5.331950E-06
  ag    9   1.588905E-03  -5.662685E-06  -4.924318E-04  -1.705477E-04   4.459729E-04   1.007671E-05  -6.308471E-04   7.922993E-05
  ag   10   1.601674E-05  -1.831656E-04  -2.594478E-05   7.439022E-05   7.488454E-06   1.333342E-04   2.114709E-06  -1.130457E-05
  ag   11  -6.702087E-04   1.713117E-05   9.704345E-04  -2.707599E-04   8.399368E-04  -1.228978E-06   1.542239E-04  -1.502186E-04
  ag   12   6.423810E-06  -2.893330E-04  -6.896347E-07  -2.766040E-04  -7.370577E-05   2.399130E-04   1.365653E-05  -1.009183E-06
  ag   13  -9.135629E-07  -4.055374E-05   1.811973E-05   5.137034E-04   1.578206E-04   4.091889E-05  -1.337898E-05  -1.790205E-05
  ag   14  -1.538918E-03   3.040031E-05   1.818074E-03   1.064772E-04  -5.155862E-04  -1.754258E-05   8.994386E-04   2.090903E-04
  ag   15  -5.759556E-05   9.764795E-05  -3.145361E-05  -1.271694E-04  -8.735220E-05  -6.949345E-05   5.544672E-05  -1.889951E-05
  ag   16   1.378441E-05  -9.316692E-04  -2.839883E-06  -1.101335E-03  -3.043461E-04   3.296123E-04   2.668649E-05   7.298567E-06
  ag   17   5.713205E-04   6.210743E-06   1.228179E-03  -3.214095E-04   1.342202E-03  -4.643197E-06  -4.537695E-04   5.236495E-05
  ag   18   4.037336E-05  -9.457400E-04   6.092376E-05   4.248035E-04   1.958546E-04   2.816074E-04  -6.506098E-05   1.723340E-05
  ag   19   1.220741E-03  -6.764483E-06   3.123654E-04   4.715071E-05  -1.574248E-04  -8.912957E-07  -7.103501E-04   3.174658E-04
  ag   20  -8.735657E-04  -2.448037E-05  -1.738485E-03  -1.899907E-04   8.295549E-04   2.028516E-05  -3.903123E-04  -1.260780E-03
  ag   21   1.052251E-05   3.812155E-05   3.001372E-05   1.096262E-03   3.415915E-04   2.384913E-04  -3.538748E-05  -1.432261E-05
  ag   22  -5.367816E-06   6.940177E-04   2.010459E-05   1.223824E-03   3.569522E-04  -5.492199E-05  -3.461813E-05  -9.441496E-06
  ag   23  -2.656525E-05   1.008892E-03  -3.643902E-06   7.665939E-04   1.916831E-04   1.874509E-04  -4.544847E-06  -2.294665E-06
  ag   24   2.750519E-05  -1.040852E-03   1.982809E-05   3.422085E-04   1.243483E-04   1.984513E-06  -2.456045E-05   1.472565E-05
  ag   25   3.687153E-03   5.223947E-06  -2.577490E-04  -1.922783E-05   1.065665E-04   3.433169E-07  -2.939620E-04   1.258831E-03
  ag   26   5.223947E-06   8.655507E-04   8.442247E-06   9.998552E-05   2.330098E-05   4.569285E-05  -5.586425E-06   6.407925E-06
  ag   27  -2.577490E-04   8.442247E-06   1.526641E-03  -2.564679E-05   2.010139E-04  -5.837034E-06  -3.846928E-04   3.719025E-04
  ag   28  -1.922783E-05   9.998552E-05  -2.564679E-05   1.503083E-03   6.498330E-05   2.899842E-04   1.638282E-05   9.031325E-05
  ag   29   1.065665E-04   2.330098E-05   2.010139E-04   6.498330E-05   1.430670E-03   7.901187E-05  -1.913410E-04  -3.461935E-04
  ag   30   3.433169E-07   4.569285E-05  -5.837034E-06   2.899842E-04   7.901187E-05   7.366703E-04  -3.396706E-06   9.051555E-07
  ag   31  -2.939620E-04  -5.586425E-06  -3.846928E-04   1.638282E-05  -1.913410E-04  -3.396706E-06   1.130555E-03  -2.683234E-04
  ag   32   1.258831E-03   6.407925E-06   3.719025E-04   9.031325E-05  -3.461935E-04   9.051555E-07  -2.683234E-04   1.727500E-03
  ag   33   2.180890E-07   3.390472E-04  -1.414764E-05  -3.217812E-04  -9.610537E-05   9.491968E-05   1.372274E-05  -3.977443E-06
  ag   34  -7.903709E-06  -3.945182E-05  -1.852288E-05  -1.825137E-05  -2.987613E-05   2.550248E-05   2.959364E-05  -1.247413E-06
  ag   35  -1.051002E-05  -2.371275E-04   5.069342E-06   4.057066E-04   1.054890E-04   2.624797E-04   4.930017E-06  -5.041225E-06
  ag   36  -9.038992E-05   2.404550E-04   3.015212E-05   1.809747E-04   4.181654E-05   1.335397E-04   1.355295E-05  -4.893391E-05
  ag   37   8.476488E-04   1.456839E-05  -3.153758E-04  -4.138035E-06   8.801373E-05   1.727938E-05  -1.857987E-04   4.733476E-04
  ag   38  -3.781212E-06  -2.375666E-04  -4.861703E-06  -8.856605E-05  -2.836177E-05   6.619088E-05   1.191892E-05  -6.404642E-06
  ag   39   3.100772E-06  -1.394266E-04   5.147095E-06  -3.156690E-05  -4.584717E-06   8.271211E-06  -9.262611E-06   3.200241E-06

               ag   33        ag   34        ag   35        ag   36        ag   37        ag   38        ag   39
  ag    3   1.317997E-05   7.536583E-04   3.919896E-04   2.843545E-04  -1.113850E-03   2.000366E-05   1.669002E-04
  ag    4  -2.057261E-03  -7.513279E-04   2.332355E-03   1.325883E-03  -1.063464E-05  -3.175031E-03   8.082754E-04
  ag    5  -4.464757E-06  -1.000679E-04   3.736882E-05   3.452140E-04  -3.139040E-03  -2.200978E-04  -9.677901E-05
  ag    6  -3.066560E-03  -6.193736E-04  -1.612397E-03  -3.512356E-03   5.603431E-05  -2.808212E-03   8.785858E-04
  ag    7  -7.054075E-06  -4.451830E-05  -1.538306E-05  -1.852292E-05   1.690587E-04   1.918827E-06   4.935589E-06
  ag    8   1.478887E-04   1.562409E-04   8.038940E-06   2.233994E-05   6.934053E-07  -1.011318E-04   5.344907E-06
  ag    9   7.018083E-06   5.260795E-05   1.003808E-05  -3.557993E-05   3.727765E-04  -3.772148E-06   5.951132E-06
  ag   10  -4.153688E-04  -3.240397E-04  -5.082026E-05  -1.008526E-04  -1.860507E-06   1.814783E-04  -3.805608E-05
  ag   11  -1.652659E-06   5.056350E-05   1.984364E-05   2.405254E-05  -2.414737E-04  -2.488921E-06  -6.304178E-06
  ag   12  -1.006961E-04  -3.031787E-04  -9.081473E-06  -1.809257E-05   4.055031E-06   2.084300E-04   5.144059E-05
  ag   13  -3.878263E-04  -2.326630E-04  -4.496163E-05  -9.007619E-05  -1.568442E-05   4.318852E-05  -8.715572E-05
  ag   14   8.183349E-06   2.384348E-05   1.230383E-05   6.754587E-05  -7.028001E-04  -5.591433E-06  -9.377566E-06
  ag   15  -7.348953E-04  -6.716305E-04  -2.062052E-04  -3.226636E-04  -4.028967E-05   3.079372E-04  -1.069925E-04
  ag   16  -3.778059E-05  -5.206786E-04   1.595991E-05  -4.444359E-05   1.293143E-05   4.607399E-04   2.570716E-04
  ag   17  -1.807034E-05  -2.369756E-05  -6.223927E-08   4.003984E-06  -6.504320E-05  -6.853550E-06   3.030976E-06
  ag   18   4.258310E-05   9.543333E-04   4.544835E-04   3.823081E-04   5.456845E-05  -2.389578E-04   2.022715E-04
  ag   19  -1.096039E-05  -7.436637E-06  -2.358265E-06  -1.811508E-05   1.666971E-04  -5.708747E-06   8.287333E-06
  ag   20   5.321690E-06  -1.518190E-05  -2.829820E-06  -1.367588E-05   2.477172E-04   8.775598E-06   2.174427E-06
  ag   21  -3.835262E-04  -3.365204E-04   5.240842E-06  -2.196507E-05  -6.326762E-06  -3.074281E-05  -1.385385E-04
  ag   22   1.213613E-04   2.933966E-04   1.126460E-04   1.689442E-04   1.932927E-06  -3.470423E-04  -1.961223E-04
  ag   23   6.699640E-04   7.148176E-04   1.673625E-04   2.979936E-04   1.083929E-05  -4.192888E-04  -2.687862E-04
  ag   24  -9.728755E-04   1.050313E-03   5.788214E-04   4.716934E-05   2.634350E-05  -1.286243E-04   1.760715E-04
  ag   25   2.180890E-07  -7.903709E-06  -1.051002E-05  -9.038992E-05   8.476488E-04  -3.781212E-06   3.100772E-06
  ag   26   3.390472E-04  -3.945182E-05  -2.371275E-04   2.404550E-04   1.456839E-05  -2.375666E-04  -1.394266E-04
  ag   27  -1.414764E-05  -1.852288E-05   5.069342E-06   3.015212E-05  -3.153758E-04  -4.861703E-06   5.147095E-06
  ag   28  -3.217812E-04  -1.825137E-05   4.057066E-04   1.809747E-04  -4.138035E-06  -8.856605E-05  -3.156690E-05
  ag   29  -9.610537E-05  -2.987613E-05   1.054890E-04   4.181654E-05   8.801373E-05  -2.836177E-05  -4.584717E-06
  ag   30   9.491968E-05   2.550248E-05   2.624797E-04   1.335397E-04   1.727938E-05   6.619088E-05   8.271211E-06
  ag   31   1.372274E-05   2.959364E-05   4.930017E-06   1.355295E-05  -1.857987E-04   1.191892E-05  -9.262611E-06
  ag   32  -3.977443E-06  -1.247413E-06  -5.041225E-06  -4.893391E-05   4.733476E-04  -6.404642E-06   3.200241E-06
  ag   33   1.458854E-03  -1.151223E-04  -3.112098E-04   3.060240E-04   2.392353E-05   1.865680E-04  -1.935036E-04
  ag   34  -1.151223E-04   1.033924E-03   3.139033E-04   2.441348E-05  -7.049473E-07  -7.362836E-05  -1.171724E-04
  ag   35  -3.112098E-04   3.139033E-04   8.545224E-04   2.288929E-05   1.718437E-06   1.243104E-04   2.284029E-04
  ag   36   3.060240E-04   2.441348E-05   2.288929E-05   4.037492E-04  -1.763635E-05   8.393534E-05  -2.497789E-05
  ag   37   2.392353E-05  -7.049473E-07   1.718437E-06  -1.763635E-05   5.756671E-04   8.970564E-06   2.200365E-06
  ag   38   1.865680E-04  -7.362836E-05   1.243104E-04   8.393534E-05   8.970564E-06   4.958328E-04   1.204013E-04
  ag   39  -1.935036E-04  -1.171724E-04   2.284029E-04  -2.497789E-05   2.200365E-06   1.204013E-04   3.237680E-04

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98741828     1.98132971     1.97932141     1.97327200     0.01532853     0.01255706
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.01021523     0.00678534     0.00562580     0.00452819     0.00228608     0.00179037     0.00147998     0.00122898
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00102274     0.00089054     0.00074856     0.00044914     0.00042201     0.00028158     0.00016610     0.00013306
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00011405     0.00010406     0.00009499     0.00008078     0.00004821     0.00002615     0.00001969     0.00001457
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000973     0.00000669     0.00000608     0.00000157     0.00000127     0.00000029     0.00000003

          modens reordered block   1

               b3u   1        b3u   2        b3u   3        b3u   4        b3u   5        b3u   6        b3u   7        b3u   8
  b3u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   2    0.00000        4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b3u   3    0.00000        0.00000        1.98469       1.046004E-04   1.798425E-03   7.324549E-04  -2.391746E-04  -1.318660E-03
  b3u   4    0.00000        0.00000       1.046004E-04    1.97835      -7.631229E-05  -6.878020E-04   1.810581E-03   8.736239E-04
  b3u   5    0.00000        0.00000       1.798425E-03  -7.631229E-05    1.97656       2.321693E-03  -2.117676E-04  -2.357238E-03
  b3u   6    0.00000        0.00000       7.324549E-04  -6.878020E-04   2.321693E-03   3.882470E-04  -5.771636E-06  -5.160076E-04
  b3u   7    0.00000        0.00000      -2.391746E-04   1.810581E-03  -2.117676E-04  -5.771636E-06   1.008946E-04   9.977804E-06
  b3u   8    0.00000        0.00000      -1.318660E-03   8.736239E-04  -2.357238E-03  -5.160076E-04   9.977804E-06   7.571053E-04
  b3u   9    0.00000        0.00000      -5.287556E-04   1.435127E-03  -4.181815E-03  -5.622885E-04   5.314131E-06   6.309182E-04
  b3u  10    0.00000        0.00000       8.000121E-04  -5.486322E-03   5.111236E-04   5.008380E-06  -3.164099E-04  -1.249730E-05
  b3u  11    0.00000        0.00000       1.333762E-03  -1.981404E-04   1.202198E-03   2.922126E-04  -9.312420E-06  -5.877958E-04
  b3u  12    0.00000        0.00000      -1.178073E-03   7.283488E-03  -7.919718E-04  -4.479865E-05   4.303773E-04   5.790681E-05
  b3u  13    0.00000        0.00000       4.319689E-04  -8.219103E-04   9.729917E-03   8.855615E-04   4.302667E-06  -8.679703E-04
  b3u  14    0.00000        0.00000       1.335533E-03  -6.428556E-03  -3.818379E-05  -1.063127E-05  -4.160848E-04   1.299935E-05
  b3u  15    0.00000        0.00000      -1.115171E-03   2.032770E-03   2.708282E-04  -5.743197E-04  -5.839359E-06   8.072138E-04
  b3u  16    0.00000        0.00000       7.482151E-05   8.230492E-04   4.649929E-03   2.441310E-04  -1.085360E-05  -4.049444E-04
  b3u  17    0.00000        0.00000      -1.245760E-03   9.382592E-04   7.944327E-04  -1.998243E-04   1.449619E-06   4.874891E-04
  b3u  18    0.00000        0.00000       2.188297E-03  -7.141907E-03  -1.854963E-03  -1.094513E-05  -3.978157E-04   1.974045E-05
  b3u  19    0.00000        0.00000      -5.006642E-05   6.827873E-05  -5.133004E-03  -3.184958E-04  -1.579348E-06   1.786030E-04
  b3u  20    0.00000        0.00000       3.442001E-04  -4.524260E-03   1.097764E-03   2.469023E-05  -2.670932E-04  -4.544208E-05
  b3u  21    0.00000        0.00000      -1.740151E-03   1.946760E-03  -2.991440E-03  -7.167002E-04  -2.021871E-06   9.503071E-04
  b3u  22    0.00000        0.00000       1.203527E-03   1.652477E-04  -3.161051E-03  -1.182167E-05   2.127515E-04   2.139486E-05
  b3u  23    0.00000        0.00000       1.441543E-03  -2.886770E-03  -1.778573E-03   4.620682E-06  -1.212920E-04  -1.373972E-05
  b3u  24    0.00000        0.00000       7.230010E-04  -4.851770E-04  -1.008994E-03   1.365571E-05   6.956189E-06   3.256937E-05
  b3u  25    0.00000        0.00000      -1.806524E-03   7.378766E-04   9.913927E-06  -3.476349E-04   2.682239E-06   4.747477E-04
  b3u  26    0.00000        0.00000      -3.644832E-05   9.912351E-05  -2.291857E-03  -1.672376E-04   2.421268E-06   1.269908E-04
  b3u  27    0.00000        0.00000      -1.708586E-05   4.030796E-03  -8.974536E-04  -3.110679E-06   2.092246E-04   3.261260E-06
  b3u  28    0.00000        0.00000      -4.857381E-04  -1.819520E-04  -2.746962E-03  -1.354879E-05  -5.650590E-06  -1.449422E-04
  b3u  29    0.00000        0.00000      -1.907282E-03  -3.979602E-04   2.912968E-03   2.474589E-04  -5.544076E-06  -4.452249E-04
  b3u  30    0.00000        0.00000       3.373987E-04  -2.678674E-03  -6.195744E-05   1.435834E-06  -7.700779E-05  -3.974194E-06
  b3u  31    0.00000        0.00000       3.009747E-04   6.992412E-05  -4.527561E-03   7.888932E-05   1.355131E-06  -5.396119E-05
  b3u  32    0.00000        0.00000       4.936741E-04   2.070297E-04   6.016710E-04  -6.743225E-05   2.643165E-06   9.950835E-05
  b3u  33    0.00000        0.00000      -4.506236E-04  -1.671315E-03   3.428363E-04  -6.019332E-06   9.968627E-05   1.487418E-05
  b3u  34    0.00000        0.00000      -2.599065E-03   6.001945E-04  -7.870143E-04  -1.649706E-04   1.534279E-06   2.267787E-04
  b3u  35    0.00000        0.00000       3.811616E-04  -3.253945E-03   5.200890E-04   8.470157E-06  -2.444240E-05  -3.563381E-06
  b3u  36    0.00000        0.00000      -4.592082E-04  -1.303841E-03  -1.550870E-03  -8.852922E-06  -8.847021E-06  -1.184691E-05
  b3u  37    0.00000        0.00000      -2.170647E-03  -3.382347E-05   2.794493E-03  -4.036438E-05  -2.714901E-06   5.678673E-05
  b3u  38    0.00000        0.00000      -9.647700E-05   4.747741E-04  -9.193600E-05  -3.432744E-07   6.093040E-05   1.816261E-06
  b3u  39    0.00000        0.00000      -1.257632E-04   1.283253E-03   7.362910E-05   6.101022E-07   1.052694E-05  -7.677644E-07

               b3u   9        b3u  10        b3u  11        b3u  12        b3u  13        b3u  14        b3u  15        b3u  16
  b3u   3  -5.287556E-04   8.000121E-04   1.333762E-03  -1.178073E-03   4.319689E-04   1.335533E-03  -1.115171E-03   7.482151E-05
  b3u   4   1.435127E-03  -5.486322E-03  -1.981404E-04   7.283488E-03  -8.219103E-04  -6.428556E-03   2.032770E-03   8.230492E-04
  b3u   5  -4.181815E-03   5.111236E-04   1.202198E-03  -7.919718E-04   9.729917E-03  -3.818379E-05   2.708282E-04   4.649929E-03
  b3u   6  -5.622885E-04   5.008380E-06   2.922126E-04  -4.479865E-05   8.855615E-04  -1.063127E-05  -5.743197E-04   2.441310E-04
  b3u   7   5.314131E-06  -3.164099E-04  -9.312420E-06   4.303773E-04   4.302667E-06  -4.160848E-04  -5.839359E-06  -1.085360E-05
  b3u   8   6.309182E-04  -1.249730E-05  -5.877958E-04   5.790681E-05  -8.679703E-04   1.299935E-05   8.072138E-04  -4.049444E-04
  b3u   9   1.101581E-03  -4.666611E-06  -5.133569E-05   8.525770E-05  -2.020579E-03  -1.548267E-05   8.874169E-04  -2.533056E-04
  b3u  10  -4.666611E-06   1.012376E-03   1.883423E-05  -1.400813E-03  -2.887348E-05   1.427365E-03   4.253804E-05   3.858709E-05
  b3u  11  -5.133569E-05   1.883423E-05   1.085197E-03  -4.446546E-05   3.464414E-04  -6.705069E-06   2.543120E-04   1.636186E-03
  b3u  12   8.525770E-05  -1.400813E-03  -4.446546E-05   1.987459E-03  -1.848895E-04  -2.086581E-03  -6.620150E-05  -1.730612E-04
  b3u  13  -2.020579E-03  -2.887348E-05   3.464414E-04  -1.848895E-04   6.126393E-03  -5.703096E-05   6.047756E-04   3.635089E-03
  b3u  14  -1.548267E-05   1.427365E-03  -6.705069E-06  -2.086581E-03  -5.703096E-05   2.525103E-03   5.607278E-05   3.017271E-05
  b3u  15   8.874169E-04   4.253804E-05   2.543120E-04  -6.620150E-05   6.047756E-04   5.607278E-05   3.234008E-03   3.160677E-03
  b3u  16  -2.533056E-04   3.858709E-05   1.636186E-03  -1.730612E-04   3.635089E-03   3.017271E-05   3.160677E-03   6.550427E-03
  b3u  17  -1.201742E-04   1.962290E-06  -7.026337E-04  -4.717869E-05   1.361919E-03  -4.917677E-06   8.441822E-04   4.646127E-04
  b3u  18  -3.445170E-05   1.513687E-03   2.687944E-05  -2.417588E-03   1.009438E-04   3.503271E-03   2.052261E-04   2.571096E-04
  b3u  19   9.748564E-04   1.463597E-05   2.634573E-04   8.903871E-05  -3.151771E-03   5.664736E-05  -4.483824E-04  -1.548979E-03
  b3u  20  -7.763975E-06   7.970149E-04   2.982561E-05  -1.050524E-03  -4.213512E-05   7.213465E-04  -4.136717E-05  -5.068314E-05
  b3u  21   1.414898E-03   3.153185E-05   1.725678E-04   2.954342E-05  -1.541405E-03   5.174832E-05   3.075129E-03   2.214354E-03
  b3u  22   8.662908E-06  -5.824803E-04   2.389157E-05   6.477296E-04   9.875074E-05  -3.375467E-05   1.111712E-04   1.626299E-04
  b3u  23   1.051460E-05   4.730529E-04   4.826608E-05  -8.107186E-04   1.985723E-05   1.201971E-03   8.886714E-05   1.048928E-04
  b3u  24  -6.600326E-05  -2.715871E-05  -6.002798E-04   7.694550E-05  -1.316582E-03  -3.805086E-05  -1.526353E-03  -3.265332E-03
  b3u  25   4.889711E-04  -2.695214E-06  -5.333444E-04   4.940527E-05  -1.225709E-03  -5.664777E-06  -1.507598E-04  -1.503502E-03
  b3u  26   3.018672E-04  -4.547464E-06   1.366438E-04   2.240996E-05  -3.106020E-04  -7.228361E-06   1.089515E-04   2.741255E-04
  b3u  27   2.182526E-05  -7.262267E-04   1.180710E-05   1.082271E-03  -4.605783E-06  -1.206069E-03  -3.515391E-06   8.740463E-06
  b3u  28   4.072969E-04   2.281067E-05   5.732039E-04   1.939628E-05  -1.486691E-03   6.292927E-05   1.553207E-05   1.132395E-04
  b3u  29  -3.074978E-04   1.362238E-05   5.188578E-04  -3.895261E-05   4.329192E-04   3.470990E-05  -5.560958E-04   3.547815E-04
  b3u  30   4.634061E-06   1.895859E-04   4.021943E-06  -2.174668E-04  -2.678279E-06  -3.168145E-05   7.728123E-06   1.697818E-06
  b3u  31  -1.774525E-04  -6.284459E-06  -2.087094E-05  -5.418461E-06   3.741315E-04  -1.206824E-05  -2.017889E-06  -1.096154E-04
  b3u  32   1.120385E-04  -6.283440E-06  -3.330976E-05   1.669309E-05  -1.873536E-04  -1.091640E-05   1.728500E-04  -3.720865E-04
  b3u  33  -5.957584E-06  -3.449977E-04  -1.664928E-05   5.359799E-04   2.136375E-05  -5.642850E-04  -2.304447E-06   9.938253E-06
  b3u  34   3.146083E-04  -3.466779E-07  -8.120718E-05   2.445116E-05  -5.596468E-04   1.188768E-07   4.555485E-04   2.076073E-05
  b3u  35  -3.903089E-05   1.210366E-04  -1.690185E-05  -2.454753E-04   1.379718E-04   3.707567E-04   1.719930E-05   2.993046E-05
  b3u  36   8.389272E-05   4.571224E-05   6.089355E-05  -8.058547E-05  -3.446908E-04   1.479702E-04  -4.536990E-05  -7.413051E-05
  b3u  37   5.610051E-05   9.680878E-06   2.993830E-05  -1.558088E-05   1.371211E-04   1.055233E-05   2.838260E-04   4.448084E-04
  b3u  38   3.339264E-07  -2.059129E-04  -3.409576E-06   2.987979E-04   7.781132E-06  -3.336861E-04  -8.205685E-06  -4.995212E-06
  b3u  39  -1.583362E-06  -2.912990E-05  -3.248527E-06   3.623210E-05  -1.888884E-06  -4.358314E-05  -1.057090E-05  -1.341349E-05

               b3u  17        b3u  18        b3u  19        b3u  20        b3u  21        b3u  22        b3u  23        b3u  24
  b3u   3  -1.245760E-03   2.188297E-03  -5.006642E-05   3.442001E-04  -1.740151E-03   1.203527E-03   1.441543E-03   7.230010E-04
  b3u   4   9.382592E-04  -7.141907E-03   6.827873E-05  -4.524260E-03   1.946760E-03   1.652477E-04  -2.886770E-03  -4.851770E-04
  b3u   5   7.944327E-04  -1.854963E-03  -5.133004E-03   1.097764E-03  -2.991440E-03  -3.161051E-03  -1.778573E-03  -1.008994E-03
  b3u   6  -1.998243E-04  -1.094513E-05  -3.184958E-04   2.469023E-05  -7.167002E-04  -1.182167E-05   4.620682E-06   1.365571E-05
  b3u   7   1.449619E-06  -3.978157E-04  -1.579348E-06  -2.670932E-04  -2.021871E-06   2.127515E-04  -1.212920E-04   6.956189E-06
  b3u   8   4.874891E-04   1.974045E-05   1.786030E-04  -4.544208E-05   9.503071E-04   2.139486E-05  -1.373972E-05   3.256937E-05
  b3u   9  -1.201742E-04  -3.445170E-05   9.748564E-04  -7.763975E-06   1.414898E-03   8.662908E-06   1.051460E-05  -6.600326E-05
  b3u  10   1.962290E-06   1.513687E-03   1.463597E-05   7.970149E-04   3.153185E-05  -5.824803E-04   4.730529E-04  -2.715871E-05
  b3u  11  -7.026337E-04   2.687944E-05   2.634573E-04   2.982561E-05   1.725678E-04   2.389157E-05   4.826608E-05  -6.002798E-04
  b3u  12  -4.717869E-05  -2.417588E-03   8.903871E-05  -1.050524E-03   2.954342E-05   6.477296E-04  -8.107186E-04   7.694550E-05
  b3u  13   1.361919E-03   1.009438E-04  -3.151771E-03  -4.213512E-05  -1.541405E-03   9.875074E-05   1.985723E-05  -1.316582E-03
  b3u  14  -4.917677E-06   3.503271E-03   5.664736E-05   7.213465E-04   5.174832E-05  -3.375467E-05   1.201971E-03  -3.805086E-05
  b3u  15   8.441822E-04   2.052261E-04  -4.483824E-04  -4.136717E-05   3.075129E-03   1.111712E-04   8.886714E-05  -1.526353E-03
  b3u  16   4.646127E-04   2.571096E-04  -1.548979E-03  -5.068314E-05   2.214354E-03   1.626299E-04   1.048928E-04  -3.265332E-03
  b3u  17   1.288637E-03   4.457142E-05  -1.143735E-03  -1.041196E-05   1.912002E-04   2.883164E-05   1.756053E-05  -4.124458E-04
  b3u  18   4.457142E-05   6.552683E-03   2.688658E-05  -1.129997E-04   1.587457E-04   2.707235E-03   3.025463E-03  -1.635306E-04
  b3u  19  -1.143735E-03   2.688658E-05   2.014561E-03  -2.264619E-06   4.482738E-04   2.042418E-05   3.232805E-05   7.920603E-04
  b3u  20  -1.041196E-05  -1.129997E-04  -2.264619E-06   1.159257E-03  -6.732081E-05  -1.786219E-03  -3.272551E-04   4.001621E-05
  b3u  21   1.912002E-04   1.587457E-04   4.482738E-04  -6.732081E-05   4.979131E-03   1.126549E-04   7.742012E-05  -1.634662E-03
  b3u  22   2.883164E-05   2.707235E-03   2.042418E-05  -1.786219E-03   1.126549E-04   5.526233E-03   2.776048E-03  -1.135049E-04
  b3u  23   1.756053E-05   3.025463E-03   3.232805E-05  -3.272551E-04   7.742012E-05   2.776048E-03   2.227870E-03  -6.815510E-05
  b3u  24  -4.124458E-04  -1.635306E-04   7.920603E-04   4.001621E-05  -1.634662E-03  -1.135049E-04  -6.815510E-05   2.561081E-03
  b3u  25   2.510318E-04  -5.125817E-05   6.016299E-04   1.347900E-05  -5.352312E-04  -2.474865E-05  -1.477954E-05   7.828624E-04
  b3u  26   3.307038E-05  -1.028613E-05   3.060116E-04  -5.555368E-06  -5.927757E-04   2.068857E-05   2.853177E-06  -2.316929E-04
  b3u  27  -1.477905E-05  -9.728972E-04   1.277805E-05  -8.161546E-04   3.573421E-05   2.037961E-03   4.787044E-04  -1.182554E-05
  b3u  28  -1.028294E-03   7.200758E-05   1.171992E-03  -9.939229E-06   1.085150E-03   6.658290E-06   1.580871E-05  -9.393339E-05
  b3u  29  -5.162392E-04   5.189297E-05   2.766405E-04   7.406905E-06  -1.362637E-03   1.768058E-05   1.085556E-05   1.493963E-04
  b3u  30   1.063189E-05  -2.531246E-04  -1.402643E-05   3.299784E-04   3.161720E-06   2.450645E-04   3.084112E-04  -4.584114E-06
  b3u  31  -2.616277E-06  -9.708745E-06  -1.454507E-05  -7.417106E-07  -2.698841E-05   7.591936E-06   5.035511E-06   4.774827E-04
  b3u  32  -5.823049E-05  -2.339701E-05   1.710908E-04   2.713407E-06  -2.350673E-04  -1.486102E-05  -7.621510E-06   9.115769E-04
  b3u  33  -7.817140E-06  -9.289291E-04  -1.042072E-05  -2.542412E-04   6.768146E-06  -4.072262E-04  -7.954548E-04  -1.791425E-07
  b3u  34   6.895520E-05   2.647341E-06   2.611915E-04  -1.207061E-05   7.522733E-04   1.208456E-05   3.049636E-06  -1.506844E-04
  b3u  35   7.085917E-05   7.922266E-04  -5.398854E-05   2.489937E-05  -1.243266E-04   1.554187E-04   3.398013E-04   3.241447E-05
  b3u  36  -1.753652E-04   3.137343E-04   1.628934E-04   3.092229E-06   2.998311E-04   7.081429E-05   1.381524E-04  -8.983519E-05
  b3u  37   1.192028E-04   2.106216E-05  -1.028822E-04   1.676147E-07   3.060139E-04   7.848993E-08   1.334023E-06  -4.823164E-04
  b3u  38  -8.751100E-08  -3.334465E-04  -3.853543E-06  -1.698681E-04  -6.013642E-06   2.973341E-04  -5.063487E-06   2.189760E-06
  b3u  39  -1.026310E-06  -1.570339E-04  -2.751546E-06   3.571312E-05  -1.271209E-05  -2.888340E-04  -1.457337E-04   9.827673E-06

               b3u  25        b3u  26        b3u  27        b3u  28        b3u  29        b3u  30        b3u  31        b3u  32
  b3u   3  -1.806524E-03  -3.644832E-05  -1.708586E-05  -4.857381E-04  -1.907282E-03   3.373987E-04   3.009747E-04   4.936741E-04
  b3u   4   7.378766E-04   9.912351E-05   4.030796E-03  -1.819520E-04  -3.979602E-04  -2.678674E-03   6.992412E-05   2.070297E-04
  b3u   5   9.913927E-06  -2.291857E-03  -8.974536E-04  -2.746962E-03   2.912968E-03  -6.195744E-05  -4.527561E-03   6.016710E-04
  b3u   6  -3.476349E-04  -1.672376E-04  -3.110679E-06  -1.354879E-05   2.474589E-04   1.435834E-06   7.888932E-05  -6.743225E-05
  b3u   7   2.682239E-06   2.421268E-06   2.092246E-04  -5.650590E-06  -5.544076E-06  -7.700779E-05   1.355131E-06   2.643165E-06
  b3u   8   4.747477E-04   1.269908E-04   3.261260E-06  -1.449422E-04  -4.452249E-04  -3.974194E-06  -5.396119E-05   9.950835E-05
  b3u   9   4.889711E-04   3.018672E-04   2.182526E-05   4.072969E-04  -3.074978E-04   4.634061E-06  -1.774525E-04   1.120385E-04
  b3u  10  -2.695214E-06  -4.547464E-06  -7.262267E-04   2.281067E-05   1.362238E-05   1.895859E-04  -6.284459E-06  -6.283440E-06
  b3u  11  -5.333444E-04   1.366438E-04   1.180710E-05   5.732039E-04   5.188578E-04   4.021943E-06  -2.087094E-05  -3.330976E-05
  b3u  12   4.940527E-05   2.240996E-05   1.082271E-03   1.939628E-05  -3.895261E-05  -2.174668E-04  -5.418461E-06   1.669309E-05
  b3u  13  -1.225709E-03  -3.106020E-04  -4.605783E-06  -1.486691E-03   4.329192E-04  -2.678279E-06   3.741315E-04  -1.873536E-04
  b3u  14  -5.664777E-06  -7.228361E-06  -1.206069E-03   6.292927E-05   3.470990E-05  -3.168145E-05  -1.206824E-05  -1.091640E-05
  b3u  15  -1.507598E-04   1.089515E-04  -3.515391E-06   1.553207E-05  -5.560958E-04   7.728123E-06  -2.017889E-06   1.728500E-04
  b3u  16  -1.503502E-03   2.741255E-04   8.740463E-06   1.132395E-04   3.547815E-04   1.697818E-06  -1.096154E-04  -3.720865E-04
  b3u  17   2.510318E-04   3.307038E-05  -1.477905E-05  -1.028294E-03  -5.162392E-04   1.063189E-05  -2.616277E-06  -5.823049E-05
  b3u  18  -5.125817E-05  -1.028613E-05  -9.728972E-04   7.200758E-05   5.189297E-05  -2.531246E-04  -9.708745E-06  -2.339701E-05
  b3u  19   6.016299E-04   3.060116E-04   1.277805E-05   1.171992E-03   2.766405E-04  -1.402643E-05  -1.454507E-05   1.710908E-04
  b3u  20   1.347900E-05  -5.555368E-06  -8.161546E-04  -9.939229E-06   7.406905E-06   3.299784E-04  -7.417106E-07   2.713407E-06
  b3u  21  -5.352312E-04  -5.927757E-04   3.573421E-05   1.085150E-03  -1.362637E-03   3.161720E-06  -2.698841E-05  -2.350673E-04
  b3u  22  -2.474865E-05   2.068857E-05   2.037961E-03   6.658290E-06   1.768058E-05   2.450645E-04   7.591936E-06  -1.486102E-05
  b3u  23  -1.477954E-05   2.853177E-06   4.787044E-04   1.580871E-05   1.085556E-05   3.084112E-04   5.035511E-06  -7.621510E-06
  b3u  24   7.828624E-04  -2.316929E-04  -1.182554E-05  -9.393339E-05   1.493963E-04  -4.584114E-06   4.774827E-04   9.115769E-04
  b3u  25   1.406080E-03   6.229923E-04  -9.756786E-06  -3.996251E-04   1.428393E-04   1.556117E-06  -4.645626E-04   4.941480E-04
  b3u  26   6.229923E-04   1.277385E-03   4.235771E-06  -2.720495E-04   3.790096E-04   8.247352E-06  -6.131953E-04   2.909850E-04
  b3u  27  -9.756786E-06   4.235771E-06   1.717664E-03  -2.533048E-06  -1.310938E-05   3.875183E-04   1.646827E-05  -7.347604E-06
  b3u  28  -3.996251E-04  -2.720495E-04  -2.533048E-06   1.484403E-03   1.874524E-04  -2.130174E-05   3.452732E-04  -2.821182E-04
  b3u  29   1.428393E-04   3.790096E-04  -1.310938E-05   1.874524E-04   1.272763E-03  -1.320094E-05  -5.226699E-05   1.465359E-04
  b3u  30   1.556117E-06   8.247352E-06   3.875183E-04  -2.130174E-05  -1.320094E-05   6.138521E-04  -5.610101E-08  -6.257307E-07
  b3u  31  -4.645626E-04  -6.131953E-04   1.646827E-05   3.452732E-04  -5.226699E-05  -5.610101E-08   1.201028E-03  -1.659662E-04
  b3u  32   4.941480E-04   2.909850E-04  -7.347604E-06  -2.821182E-04   1.465359E-04  -6.257307E-07  -1.659662E-04   1.125747E-03
  b3u  33   7.737746E-07  -8.548818E-06  -5.678416E-05   1.584475E-05   1.785221E-05  -2.272586E-04   1.753896E-06  -3.628875E-06
  b3u  34   4.643569E-04   1.065801E-04   5.539721E-06   8.752407E-05  -8.191656E-05  -1.727024E-06  -3.863474E-04   2.185421E-04
  b3u  35   2.439148E-05  -2.999650E-05  -3.985914E-04  -7.156918E-05   6.599824E-05  -1.468927E-04   8.899384E-05   2.091011E-05
  b3u  36  -7.020548E-05   7.124751E-05  -1.489482E-04   1.991136E-04  -1.624702E-04  -5.585117E-05  -2.308730E-04  -5.148823E-05
  b3u  37  -7.303032E-05  -4.050082E-05  -7.435850E-06   1.832542E-05   1.147920E-04  -2.968049E-06  -9.874126E-05  -3.787763E-04
  b3u  38   2.070462E-06   5.267681E-06   2.473745E-04  -7.675644E-06  -1.911838E-06  -1.018773E-05  -1.221399E-06   2.255629E-06
  b3u  39   8.686787E-07   1.396952E-06  -7.845661E-05  -4.607497E-06  -6.669690E-06  -5.898540E-05  -7.087373E-07   2.470652E-06

               b3u  33        b3u  34        b3u  35        b3u  36        b3u  37        b3u  38        b3u  39
  b3u   3  -4.506236E-04  -2.599065E-03   3.811616E-04  -4.592082E-04  -2.170647E-03  -9.647700E-05  -1.257632E-04
  b3u   4  -1.671315E-03   6.001945E-04  -3.253945E-03  -1.303841E-03  -3.382347E-05   4.747741E-04   1.283253E-03
  b3u   5   3.428363E-04  -7.870143E-04   5.200890E-04  -1.550870E-03   2.794493E-03  -9.193600E-05   7.362910E-05
  b3u   6  -6.019332E-06  -1.649706E-04   8.470157E-06  -8.852922E-06  -4.036438E-05  -3.432744E-07   6.101022E-07
  b3u   7   9.968627E-05   1.534279E-06  -2.444240E-05  -8.847021E-06  -2.714901E-06   6.093040E-05   1.052694E-05
  b3u   8   1.487418E-05   2.267787E-04  -3.563381E-06  -1.184691E-05   5.678673E-05   1.816261E-06  -7.677644E-07
  b3u   9  -5.957584E-06   3.146083E-04  -3.903089E-05   8.389272E-05   5.610051E-05   3.339264E-07  -1.583362E-06
  b3u  10  -3.449977E-04  -3.466779E-07   1.210366E-04   4.571224E-05   9.680878E-06  -2.059129E-04  -2.912990E-05
  b3u  11  -1.664928E-05  -8.120718E-05  -1.690185E-05   6.089355E-05   2.993830E-05  -3.409576E-06  -3.248527E-06
  b3u  12   5.359799E-04   2.445116E-05  -2.454753E-04  -8.058547E-05  -1.558088E-05   2.987979E-04   3.623210E-05
  b3u  13   2.136375E-05  -5.596468E-04   1.379718E-04  -3.446908E-04   1.371211E-04   7.781132E-06  -1.888884E-06
  b3u  14  -5.642850E-04   1.188768E-07   3.707567E-04   1.479702E-04   1.055233E-05  -3.336861E-04  -4.358314E-05
  b3u  15  -2.304447E-06   4.555485E-04   1.719930E-05  -4.536990E-05   2.838260E-04  -8.205685E-06  -1.057090E-05
  b3u  16   9.938253E-06   2.076073E-05   2.993046E-05  -7.413051E-05   4.448084E-04  -4.995212E-06  -1.341349E-05
  b3u  17  -7.817140E-06   6.895520E-05   7.085917E-05  -1.753652E-04   1.192028E-04  -8.751100E-08  -1.026310E-06
  b3u  18  -9.289291E-04   2.647341E-06   7.922266E-04   3.137343E-04   2.106216E-05  -3.334465E-04  -1.570339E-04
  b3u  19  -1.042072E-05   2.611915E-04  -5.398854E-05   1.628934E-04  -1.028822E-04  -3.853543E-06  -2.751546E-06
  b3u  20  -2.542412E-04  -1.207061E-05   2.489937E-05   3.092229E-06   1.676147E-07  -1.698681E-04   3.571312E-05
  b3u  21   6.768146E-06   7.522733E-04  -1.243266E-04   2.998311E-04   3.060139E-04  -6.013642E-06  -1.271209E-05
  b3u  22  -4.072262E-04   1.208456E-05   1.554187E-04   7.081429E-05   7.848993E-08   2.973341E-04  -2.888340E-04
  b3u  23  -7.954548E-04   3.049636E-06   3.398013E-04   1.381524E-04   1.334023E-06  -5.063487E-06  -1.457337E-04
  b3u  24  -1.791425E-07  -1.506844E-04   3.241447E-05  -8.983519E-05  -4.823164E-04   2.189760E-06   9.827673E-06
  b3u  25   7.737746E-07   4.643569E-04   2.439148E-05  -7.020548E-05  -7.303032E-05   2.070462E-06   8.686787E-07
  b3u  26  -8.548818E-06   1.065801E-04  -2.999650E-05   7.124751E-05  -4.050082E-05   5.267681E-06   1.396952E-06
  b3u  27  -5.678416E-05   5.539721E-06  -3.985914E-04  -1.489482E-04  -7.435850E-06   2.473745E-04  -7.845661E-05
  b3u  28   1.584475E-05   8.752407E-05  -7.156918E-05   1.991136E-04   1.832542E-05  -7.675644E-06  -4.607497E-06
  b3u  29   1.785221E-05  -8.191656E-05   6.599824E-05  -1.624702E-04   1.147920E-04  -1.911838E-06  -6.669690E-06
  b3u  30  -2.272586E-04  -1.727024E-06  -1.468927E-04  -5.585117E-05  -2.968049E-06  -1.018773E-05  -5.898540E-05
  b3u  31   1.753896E-06  -3.863474E-04   8.899384E-05  -2.308730E-04  -9.874126E-05  -1.221399E-06  -7.087373E-07
  b3u  32  -3.628875E-06   2.185421E-04   2.091011E-05  -5.148823E-05  -3.787763E-04   2.255629E-06   2.470652E-06
  b3u  33   8.117716E-04   8.101044E-06  -1.428883E-04  -5.828418E-05   9.372216E-06   4.456591E-05  -8.624562E-05
  b3u  34   8.101044E-06   5.962158E-04  -4.513333E-05   1.073687E-04   1.830803E-05   1.105529E-06  -3.981156E-06
  b3u  35  -1.428883E-04  -4.513333E-05   4.603152E-04   4.381254E-05  -6.222134E-06   6.780745E-05  -2.773172E-05
  b3u  36  -5.828418E-05   1.073687E-04   4.381254E-05   3.670219E-04   1.294336E-05   2.739291E-05  -1.090625E-05
  b3u  37   9.372216E-06   1.830803E-05  -6.222134E-06   1.294336E-05   4.471556E-04  -2.787568E-06  -3.328054E-06
  b3u  38   4.456591E-05   1.105529E-06   6.780745E-05   2.739291E-05  -2.787568E-06   2.542896E-04  -7.958824E-05
  b3u  39  -8.624562E-05  -3.981156E-06  -2.773172E-05  -1.090625E-05  -3.328054E-06  -7.958824E-05   1.066260E-04

Natural orbital populations,block 2
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     1.98509405     1.97848905     1.97629978     0.01483023     0.01281238     0.01068945
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00863793     0.00484928     0.00336972     0.00208965     0.00188776     0.00123506     0.00078421     0.00067525
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00043150     0.00042411     0.00034400     0.00026656     0.00020953     0.00013663     0.00009844     0.00007430
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30       MO    31       MO    32
  occ(*)=     0.00004997     0.00003203     0.00002926     0.00002682     0.00001801     0.00001067     0.00000733     0.00000633
              MO    33       MO    34       MO    35       MO    36       MO    37       MO    38       MO    39
  occ(*)=     0.00000433     0.00000175     0.00000155     0.00000078     0.00000022     0.00000010     0.00000000

          modens reordered block   1

               b2u   1        b2u   2        b2u   3        b2u   4        b2u   5        b2u   6        b2u   7        b2u   8
  b2u   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b2u   2    0.00000        1.98427      -2.319120E-04  -1.777465E-03  -1.053748E-03  -4.212432E-04  -1.982461E-03  -1.061383E-03
  b2u   3    0.00000      -2.319120E-04    1.97390      -4.486697E-04  -7.950766E-05  -6.858825E-04  -2.708251E-06  -9.002689E-05
  b2u   4    0.00000      -1.777465E-03  -4.486697E-04    1.97653       3.153882E-03  -3.730079E-04   3.180318E-03   5.965780E-03
  b2u   5    0.00000      -1.053748E-03  -7.950766E-05   3.153882E-03   4.215572E-04   1.075898E-05   5.678412E-04   5.671831E-04
  b2u   6    0.00000      -4.212432E-04  -6.858825E-04  -3.730079E-04   1.075898E-05   4.427654E-04   1.481231E-05  -3.456102E-06
  b2u   7    0.00000      -1.982461E-03  -2.708251E-06   3.180318E-03   5.678412E-04   1.481231E-05   8.411891E-04   6.291512E-04
  b2u   8    0.00000      -1.061383E-03  -9.002689E-05   5.965780E-03   5.671831E-04  -3.456102E-06   6.291512E-04   1.125084E-03
  b2u   9    0.00000      -2.301293E-03   9.063231E-04   1.279837E-03   3.451949E-04   4.365989E-05   6.674271E-04   4.745129E-05
  b2u  10    0.00000       6.514684E-04  -1.599408E-03  -1.265328E-02  -9.222068E-04  -4.603588E-05  -9.078563E-04  -2.091819E-03
  b2u  11    0.00000       7.839578E-04   1.791541E-03   1.006366E-03  -1.697052E-07  -9.211392E-04   1.474697E-05  -7.415484E-06
  b2u  12    0.00000       2.396846E-03   2.105390E-03  -1.701952E-04  -5.533322E-04   4.860887E-05  -7.960016E-04  -8.317637E-04
  b2u  13    0.00000      -7.923202E-04   4.579432E-03   6.036758E-03   2.907489E-04   7.926222E-05   4.412618E-04   3.237966E-04
  b2u  14    0.00000       2.649938E-03   1.517703E-04   7.994463E-04  -2.089737E-04   4.200089E-06  -4.954476E-04   1.349751E-04
  b2u  15    0.00000      -6.349539E-04  -8.293597E-04  -6.537281E-03  -3.342695E-04  -2.329367E-05  -2.083261E-04  -9.963925E-04
  b2u  16    0.00000       7.283365E-04   9.710658E-04   1.528904E-03   5.505621E-06  -9.418442E-04   1.852531E-05   1.653752E-05
  b2u  17    0.00000       3.895980E-03   2.462615E-03  -2.259217E-03  -7.123986E-04   1.052593E-05  -9.614116E-04  -1.384951E-03
  b2u  18    0.00000       1.124781E-03   2.970554E-03   1.516989E-03  -9.323413E-06   2.066988E-05   2.647901E-05  -4.256014E-05
  b2u  19    0.00000      -1.503450E-03   1.234607E-03   1.955206E-03   3.547054E-04   1.372471E-05   4.795710E-04   4.984433E-04
  b2u  20    0.00000      -3.280683E-04   3.738971E-03  -6.340147E-04   3.603476E-06   3.917195E-04   7.019376E-06  -2.140995E-05
  b2u  21    0.00000       1.321559E-03   4.700038E-05   4.196225E-03   1.858246E-04   9.559422E-06   1.606712E-04   3.013912E-04
  b2u  22    0.00000      -3.481595E-04  -4.104645E-04   2.461506E-03   1.616047E-05   6.240598E-06  -1.253529E-04   3.964538E-04
  b2u  23    0.00000      -2.770251E-04  -7.640862E-03  -1.323822E-03  -1.225308E-05   5.188992E-04  -2.829701E-05  -2.008532E-06
  b2u  24    0.00000       6.780426E-04   8.169311E-05   2.617271E-03   2.575752E-04  -5.388759E-06   4.503852E-04   3.258143E-04
  b2u  25    0.00000       4.437780E-04  -2.763591E-04  -3.423347E-03   7.820625E-05  -1.508484E-06   5.030288E-05   1.806032E-04
  b2u  26    0.00000       1.378988E-03   9.339205E-04   8.005273E-05   6.696050E-05  -9.195553E-07   9.796270E-05   1.129164E-04
  b2u  27    0.00000      -1.273707E-04  -4.174120E-03  -6.397948E-04  -2.604965E-06   2.493050E-05  -4.715748E-06   1.034593E-06
  b2u  28    0.00000       2.754811E-03   1.220649E-04  -1.299772E-03  -1.668557E-04  -2.829297E-06  -2.281472E-04  -3.203523E-04
  b2u  29    0.00000      -3.685340E-04  -2.258962E-04   2.129331E-03   1.622851E-05  -2.434078E-06  -4.829402E-06   1.042394E-04
  b2u  30    0.00000       2.582496E-03   5.563515E-04   2.842830E-03  -3.576052E-05  -2.307834E-06  -5.273765E-05  -4.619949E-05

               b2u   9        b2u  10        b2u  11        b2u  12        b2u  13        b2u  14        b2u  15        b2u  16
  b2u   2  -2.301293E-03   6.514684E-04   7.839578E-04   2.396846E-03  -7.923202E-04   2.649938E-03  -6.349539E-04   7.283365E-04
  b2u   3   9.063231E-04  -1.599408E-03   1.791541E-03   2.105390E-03   4.579432E-03   1.517703E-04  -8.293597E-04   9.710658E-04
  b2u   4   1.279837E-03  -1.265328E-02   1.006366E-03  -1.701952E-04   6.036758E-03   7.994463E-04  -6.537281E-03   1.528904E-03
  b2u   5   3.451949E-04  -9.222068E-04  -1.697052E-07  -5.533322E-04   2.907489E-04  -2.089737E-04  -3.342695E-04   5.505621E-06
  b2u   6   4.365989E-05  -4.603588E-05  -9.211392E-04   4.860887E-05   7.926222E-05   4.200089E-06  -2.329367E-05  -9.418442E-04
  b2u   7   6.674271E-04  -9.078563E-04   1.474697E-05  -7.960016E-04   4.412618E-04  -4.954476E-04  -2.083261E-04   1.852531E-05
  b2u   8   4.745129E-05  -2.091819E-03  -7.415484E-06  -8.317637E-04   3.237966E-04   1.349751E-04  -9.963925E-04   1.653752E-05
  b2u   9   1.129852E-03  -3.794383E-04  -2.424161E-05   2.331245E-04   1.586893E-03  -6.702237E-04   1.979775E-04  -2.223039E-05
  b2u  10  -3.794383E-04   6.259143E-03   1.030869E-04  -6.916569E-04  -3.795496E-03  -1.354348E-03   3.179364E-03   2.779454E-05
  b2u  11  -2.424161E-05   1.030869E-04   2.028085E-03  -9.615183E-05  -1.235178E-04  -3.377622E-05   6.196082E-05   2.130463E-03
  b2u  12   2.331245E-04  -6.916569E-04  -9.615183E-05   3.215496E-03   3.168557E-03   8.762179E-04  -5.058015E-04  -8.364222E-05
  b2u  13   1.586893E-03  -3.795496E-03  -1.235178E-04   3.168557E-03   6.514820E-03   5.765697E-04  -1.694238E-03  -7.742926E-05
  b2u  14  -6.702237E-04  -1.354348E-03  -3.377622E-05   8.762179E-04   5.765697E-04   1.227981E-03  -1.085828E-03  -2.156525E-05
  b2u  15   1.979775E-04   3.179364E-03   6.196082E-05  -5.058015E-04  -1.694238E-03  -1.085828E-03   1.969261E-03   2.845586E-05
  b2u  16  -2.223039E-05   2.779454E-05   2.130463E-03  -8.364222E-05  -7.742926E-05  -2.156525E-05   2.845586E-05   2.363189E-03
  b2u  17   1.582887E-04   1.458018E-03  -2.061543E-05   3.091748E-03   2.265871E-03   2.289675E-04   4.001679E-04  -4.097261E-05
  b2u  18   5.651618E-04  -1.361690E-03  -3.385583E-05   1.525712E-03   3.245276E-03   4.566061E-04  -8.355796E-04  -2.298083E-05
  b2u  19   5.373478E-04  -1.289199E-03  -2.064081E-05   1.962965E-04   1.555181E-03  -2.224186E-04  -6.364958E-04  -9.646641E-06
  b2u  20   2.721404E-05   1.363887E-05  -8.556134E-04   2.067306E-05   1.844283E-05   1.084829E-05  -1.883262E-08  -1.057405E-03
  b2u  21  -7.344150E-05  -3.500544E-04  -2.115184E-05  -6.931629E-05  -1.678756E-04  -5.296605E-05  -3.005557E-04  -2.329328E-05
  b2u  22  -5.147957E-04  -1.444310E-03  -2.346399E-05   9.180982E-06  -3.691657E-05   9.683474E-04  -1.105557E-03  -1.055486E-05
  b2u  23  -1.859956E-05   1.582446E-05  -1.393379E-03   1.794761E-05  -1.213178E-05  -1.565172E-07   3.164951E-06  -1.668887E-03
  b2u  24   4.779818E-04  -4.792647E-04   1.311289E-05  -5.895728E-04   2.841358E-04  -4.612475E-04   2.079893E-04   2.443566E-05
  b2u  25  -2.097287E-05  -3.944720E-04  -3.301807E-06   2.230519E-05  -8.344207E-05   2.977165E-06  -2.562977E-05  -3.067328E-06
  b2u  26   3.522249E-05  -1.995853E-04   3.944242E-06  -1.609184E-04   3.926616E-04   6.189403E-05  -1.759944E-04   7.275980E-06
  b2u  27  -8.839007E-06   1.004705E-05  -1.030226E-04  -1.020435E-05  -1.133051E-05  -8.425365E-07   3.523553E-06  -2.933158E-04
  b2u  28  -8.004643E-05   5.807520E-04   6.271262E-06   4.437679E-04   1.005657E-05   6.473710E-05   2.709450E-04  -1.379563E-06
  b2u  29  -6.140293E-05  -3.947734E-04   3.825877E-07   5.116949E-05   9.692092E-05   1.914582E-04  -1.812470E-04   3.362178E-06
  b2u  30   2.756190E-05  -1.523868E-04   2.700048E-06   2.786265E-04   4.483978E-04   1.235829E-04  -1.139823E-04   2.417381E-06

               b2u  17        b2u  18        b2u  19        b2u  20        b2u  21        b2u  22        b2u  23        b2u  24
  b2u   2   3.895980E-03   1.124781E-03  -1.503450E-03  -3.280683E-04   1.321559E-03  -3.481595E-04  -2.770251E-04   6.780426E-04
  b2u   3   2.462615E-03   2.970554E-03   1.234607E-03   3.738971E-03   4.700038E-05  -4.104645E-04  -7.640862E-03   8.169311E-05
  b2u   4  -2.259217E-03   1.516989E-03   1.955206E-03  -6.340147E-04   4.196225E-03   2.461506E-03  -1.323822E-03   2.617271E-03
  b2u   5  -7.123986E-04  -9.323413E-06   3.547054E-04   3.603476E-06   1.858246E-04   1.616047E-05  -1.225308E-05   2.575752E-04
  b2u   6   1.052593E-05   2.066988E-05   1.372471E-05   3.917195E-04   9.559422E-06   6.240598E-06   5.188992E-04  -5.388759E-06
  b2u   7  -9.614116E-04   2.647901E-05   4.795710E-04   7.019376E-06   1.606712E-04  -1.253529E-04  -2.829701E-05   4.503852E-04
  b2u   8  -1.384951E-03  -4.256014E-05   4.984433E-04  -2.140995E-05   3.013912E-04   3.964538E-04  -2.008532E-06   3.258143E-04
  b2u   9   1.582887E-04   5.651618E-04   5.373478E-04   2.721404E-05  -7.344150E-05  -5.147957E-04  -1.859956E-05   4.779818E-04
  b2u  10   1.458018E-03  -1.361690E-03  -1.289199E-03   1.363887E-05  -3.500544E-04  -1.444310E-03   1.582446E-05  -4.792647E-04
  b2u  11  -2.061543E-05  -3.385583E-05  -2.064081E-05  -8.556134E-04  -2.115184E-05  -2.346399E-05  -1.393379E-03   1.311289E-05
  b2u  12   3.091748E-03   1.525712E-03   1.962965E-04   2.067306E-05  -6.931629E-05   9.180982E-06   1.794761E-05  -5.895728E-04
  b2u  13   2.265871E-03   3.245276E-03   1.555181E-03   1.844283E-05  -1.678756E-04  -3.691657E-05  -1.213178E-05   2.841358E-04
  b2u  14   2.289675E-04   4.566061E-04  -2.224186E-04   1.084829E-05  -5.296605E-05   9.683474E-04  -1.565172E-07  -4.612475E-04
  b2u  15   4.001679E-04  -8.355796E-04  -6.364958E-04  -1.883262E-08  -3.005557E-04  -1.105557E-03   3.164951E-06   2.079893E-04
  b2u  16  -4.097261E-05  -2.298083E-05  -9.646641E-06  -1.057405E-03  -2.329328E-05  -1.055486E-05  -1.668887E-03   2.443566E-05
  b2u  17   5.004961E-03   1.660476E-03   5.567427E-04  -3.325144E-06   6.035667E-04  -1.060057E-03   2.243983E-06  -1.389851E-03
  b2u  18   1.660476E-03   2.560022E-03   7.974615E-04   1.119567E-05  -1.956921E-04  -6.885048E-05  -1.675179E-05  -1.720968E-04
  b2u  19   5.567427E-04   7.974615E-04   1.405145E-03  -1.280026E-05   6.349554E-04  -3.879892E-04  -9.067977E-06  -1.321669E-04
  b2u  20  -3.325144E-06   1.119567E-05  -1.280026E-05   9.161441E-04  -3.585157E-06   1.764765E-05   4.976712E-04  -1.069203E-05
  b2u  21   6.035667E-04  -1.956921E-04   6.349554E-04  -3.585157E-06   1.264211E-03  -2.910260E-04   3.857401E-06  -3.548517E-04
  b2u  22  -1.060057E-03  -6.885048E-05  -3.879892E-04   1.764765E-05  -2.910260E-04   1.429231E-03  -4.043016E-06  -1.335163E-04
  b2u  23   2.243983E-06  -1.675179E-05  -9.067977E-06   4.976712E-04   3.857401E-06  -4.043016E-06   2.058185E-03  -2.987670E-06
  b2u  24  -1.389851E-03  -1.720968E-04  -1.321669E-04  -1.069203E-05  -3.548517E-04  -1.335163E-04  -2.987670E-06   1.210970E-03
  b2u  25  -1.480144E-05  -4.682851E-04   4.601341E-04  -1.308857E-05   6.157336E-04  -3.452634E-04   3.349526E-06  -4.640157E-05
  b2u  26   2.466518E-04   9.245466E-04   4.934304E-04  -4.114472E-06   2.977667E-04  -2.829122E-04  -2.347929E-05  -1.378377E-04
  b2u  27  -1.214632E-07   1.040219E-05   3.084572E-06   2.170684E-04   6.294956E-06  -9.530854E-06   5.278632E-04  -3.992720E-06
  b2u  28   7.525060E-04   1.539149E-04  -4.621617E-04   4.770189E-06  -1.166501E-04  -8.476278E-05   1.130398E-06  -9.926925E-05
  b2u  29  -3.250742E-04  -8.678764E-05  -6.862407E-05  -3.872870E-06   8.137970E-05   2.125314E-04  -1.856933E-06   1.806198E-04
  b2u  30   3.128145E-04   4.877685E-04   8.941905E-05   2.527864E-06   5.115779E-05  -1.446915E-05  -6.963419E-06   1.013508E-04

               b2u  25        b2u  26        b2u  27        b2u  28        b2u  29        b2u  30
  b2u   2   4.437780E-04   1.378988E-03  -1.273707E-04   2.754811E-03  -3.685340E-04   2.582496E-03
  b2u   3  -2.763591E-04   9.339205E-04  -4.174120E-03   1.220649E-04  -2.258962E-04   5.563515E-04
  b2u   4  -3.423347E-03   8.005273E-05  -6.397948E-04  -1.299772E-03   2.129331E-03   2.842830E-03
  b2u   5   7.820625E-05   6.696050E-05  -2.604965E-06  -1.668557E-04   1.622851E-05  -3.576052E-05
  b2u   6  -1.508484E-06  -9.195553E-07   2.493050E-05  -2.829297E-06  -2.434078E-06  -2.307834E-06
  b2u   7   5.030288E-05   9.796270E-05  -4.715748E-06  -2.281472E-04  -4.829402E-06  -5.273765E-05
  b2u   8   1.806032E-04   1.129164E-04   1.034593E-06  -3.203523E-04   1.042394E-04  -4.619949E-05
  b2u   9  -2.097287E-05   3.522249E-05  -8.839007E-06  -8.004643E-05  -6.140293E-05   2.756190E-05
  b2u  10  -3.944720E-04  -1.995853E-04   1.004705E-05   5.807520E-04  -3.947734E-04  -1.523868E-04
  b2u  11  -3.301807E-06   3.944242E-06  -1.030226E-04   6.271262E-06   3.825877E-07   2.700048E-06
  b2u  12   2.230519E-05  -1.609184E-04  -1.020435E-05   4.437679E-04   5.116949E-05   2.786265E-04
  b2u  13  -8.344207E-05   3.926616E-04  -1.133051E-05   1.005657E-05   9.692092E-05   4.483978E-04
  b2u  14   2.977165E-06   6.189403E-05  -8.425365E-07   6.473710E-05   1.914582E-04   1.235829E-04
  b2u  15  -2.562977E-05  -1.759944E-04   3.523553E-06   2.709450E-04  -1.812470E-04  -1.139823E-04
  b2u  16  -3.067328E-06   7.275980E-06  -2.933158E-04  -1.379563E-06   3.362178E-06   2.417381E-06
  b2u  17  -1.480144E-05   2.466518E-04  -1.214632E-07   7.525060E-04  -3.250742E-04   3.128145E-04
  b2u  18  -4.682851E-04   9.245466E-04   1.040219E-05   1.539149E-04  -8.678764E-05   4.877685E-04
  b2u  19   4.601341E-04   4.934304E-04   3.084572E-06  -4.621617E-04  -6.862407E-05   8.941905E-05
  b2u  20  -1.308857E-05  -4.114472E-06   2.170684E-04   4.770189E-06  -3.872870E-06   2.527864E-06
  b2u  21   6.157336E-04   2.977667E-04   6.294956E-06  -1.166501E-04   8.137970E-05   5.115779E-05
  b2u  22  -3.452634E-04  -2.829122E-04  -9.530854E-06  -8.476278E-05   2.125314E-04  -1.446915E-05
  b2u  23   3.349526E-06  -2.347929E-05   5.278632E-04   1.130398E-06  -1.856933E-06  -6.963419E-06
  b2u  24  -4.640157E-05  -1.378377E-04  -3.992720E-06  -9.926925E-05   1.806198E-04   1.013508E-04
  b2u  25   1.191106E-03   1.619024E-04   4.229522E-06  -3.883314E-04   2.446734E-04  -9.169453E-05
  b2u  26   1.619024E-04   1.127319E-03   6.302153E-06  -2.126268E-04  -5.157555E-05   3.861774E-04
  b2u  27   4.229522E-06   6.302153E-06   6.396852E-04  -3.142707E-06  -1.494179E-06   2.934928E-06
  b2u  28  -3.883314E-04  -2.126268E-04  -3.142707E-06   5.950953E-04  -1.199255E-04   1.447896E-05
  b2u  29   2.446734E-04  -5.157555E-05  -1.494179E-06  -1.199255E-04   3.507160E-04  -1.142037E-05
  b2u  30  -9.169453E-05   3.861774E-04   2.934928E-06   1.447896E-05  -1.142037E-05   4.473648E-04

Natural orbital populations,block 3
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98470871     1.97640519     1.97388597     0.01497696     0.01060731     0.00642796     0.00463081
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00327523     0.00209351     0.00123723     0.00108675     0.00074705     0.00064658     0.00046507     0.00037884
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00028160     0.00020824     0.00013511     0.00010154     0.00004254     0.00003009     0.00002657     0.00002035
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00001066     0.00000837     0.00000398     0.00000187     0.00000065     0.00000014

          modens reordered block   1

               b1g   1        b1g   2        b1g   3        b1g   4        b1g   5        b1g   6        b1g   7        b1g   8
  b1g   1    4.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000    
  b1g   2    0.00000        1.98027       2.458518E-03   1.343955E-03   2.289019E-03   2.822087E-03   3.214810E-05  -1.192778E-03
  b1g   3    0.00000       2.458518E-03    1.97353       2.041028E-03   2.993109E-03   3.728775E-03  -7.580660E-05  -5.661037E-04
  b1g   4    0.00000       1.343955E-03   2.041028E-03   2.279646E-04   4.115131E-04   3.623153E-04   7.249339E-08  -2.066701E-04
  b1g   5    0.00000       2.289019E-03   2.993109E-03   4.115131E-04   8.770675E-04   5.379240E-04   1.225693E-06  -5.930506E-04
  b1g   6    0.00000       2.822087E-03   3.728775E-03   3.623153E-04   5.379240E-04   7.400072E-04  -5.394469E-07  -5.455041E-05
  b1g   7    0.00000       3.214810E-05  -7.580660E-05   7.249339E-08   1.225693E-06  -5.394469E-07   9.065616E-06  -1.212820E-06
  b1g   8    0.00000      -1.192778E-03  -5.661037E-04  -2.066701E-04  -5.930506E-04  -5.455041E-05  -1.212820E-06   7.520520E-04
  b1g   9    0.00000      -3.173191E-03  -3.275942E-03  -6.030665E-04  -1.587046E-03  -6.349946E-04  -3.721132E-06   1.127592E-03
  b1g  10    0.00000       5.108926E-03   5.738730E-03   5.729922E-04   6.057962E-04   1.514565E-03  -1.477415E-06   4.928233E-04
  b1g  11    0.00000      -2.793260E-03  -2.491940E-04  -3.555627E-04  -1.316134E-03  -1.777903E-04  -2.349836E-06   1.054573E-03
  b1g  12    0.00000       4.627166E-04  -1.333386E-03   6.109411E-06   3.710120E-05  -8.439311E-06   8.880031E-05  -3.811515E-05
  b1g  13    0.00000       1.518414E-03   8.124847E-04   1.816083E-04   5.558819E-04  -6.735420E-05   3.757478E-07  -9.466232E-04
  b1g  14    0.00000       3.588431E-03   4.585929E-03   3.809455E-04   4.485928E-04   1.049303E-03  -3.065728E-07   4.320025E-04
  b1g  15    0.00000       1.382715E-03  -3.934631E-03   8.492665E-07   3.492944E-05  -2.658140E-05   1.931696E-04  -4.919150E-05
  b1g  16    0.00000       4.379883E-03   2.249812E-03   5.048924E-04   8.622909E-04   1.060258E-03   1.720384E-06  -1.753521E-04
  b1g  17    0.00000       2.544113E-03  -9.293355E-04   9.173308E-05   4.194558E-04   6.462819E-05   2.367528E-06  -3.650891E-04
  b1g  18    0.00000      -1.211530E-03   2.557501E-03  -3.162660E-06  -2.293256E-05   8.536982E-06  -6.897617E-05   2.932113E-05
  b1g  19    0.00000       2.601879E-03   3.209644E-04   1.332076E-04   1.940603E-04   3.023462E-04  -4.354163E-08  -3.310128E-05
  b1g  20    0.00000       1.570853E-03   1.642231E-03   3.470905E-06  -9.514594E-05   2.790655E-04   6.200491E-07   5.522357E-04
  b1g  21    0.00000      -1.208792E-04  -4.549432E-03  -1.284244E-04  -1.379031E-04  -2.443093E-04   2.137231E-07   4.890012E-05
  b1g  22    0.00000      -6.005193E-04   6.683967E-04  -1.921646E-07  -2.848404E-06   2.858069E-06  -3.647486E-06   5.956021E-06
  b1g  23    0.00000       2.462747E-03  -4.903062E-03  -1.643866E-04  -4.498140E-04  -1.168055E-04   3.584750E-07   4.153975E-04
  b1g  24    0.00000      -5.981434E-05  -3.901243E-05   1.553336E-04   3.263167E-04   3.010251E-04   1.196319E-07  -2.428547E-04
  b1g  25    0.00000       3.489271E-04  -2.998334E-04   1.199615E-06   2.744226E-06   2.423346E-06  -3.672184E-05  -2.772164E-06
  b1g  26    0.00000       2.782152E-03   6.284950E-04   4.965222E-06   4.517737E-05   6.568466E-06   4.724273E-07  -4.130204E-05
  b1g  27    0.00000       1.468095E-03   4.008657E-03   2.672058E-05   1.073609E-04   2.506046E-05  -3.617463E-07  -9.477619E-05
  b1g  28    0.00000       3.743893E-03  -3.059447E-03   9.913858E-05   1.763034E-04   2.066129E-04   5.756980E-07  -4.058554E-05
  b1g  29    0.00000       3.405207E-04  -1.656242E-04   1.258020E-06   3.231384E-06   1.360975E-06  -4.868060E-06  -3.433215E-06
  b1g  30    0.00000      -9.856928E-04   1.119757E-03  -2.259402E-06  -3.006754E-05   5.640438E-05  -2.593157E-07   8.302839E-05

               b1g   9        b1g  10        b1g  11        b1g  12        b1g  13        b1g  14        b1g  15        b1g  16
  b1g   2  -3.173191E-03   5.108926E-03  -2.793260E-03   4.627166E-04   1.518414E-03   3.588431E-03   1.382715E-03   4.379883E-03
  b1g   3  -3.275942E-03   5.738730E-03  -2.491940E-04  -1.333386E-03   8.124847E-04   4.585929E-03  -3.934631E-03   2.249812E-03
  b1g   4  -6.030665E-04   5.729922E-04  -3.555627E-04   6.109411E-06   1.816083E-04   3.809455E-04   8.492665E-07   5.048924E-04
  b1g   5  -1.587046E-03   6.057962E-04  -1.316134E-03   3.710120E-05   5.558819E-04   4.485928E-04   3.492944E-05   8.622909E-04
  b1g   6  -6.349946E-04   1.514565E-03  -1.777903E-04  -8.439311E-06  -6.735420E-05   1.049303E-03  -2.658140E-05   1.060258E-03
  b1g   7  -3.721132E-06  -1.477415E-06  -2.349836E-06   8.880031E-05   3.757478E-07  -3.065728E-07   1.931696E-04   1.720384E-06
  b1g   8   1.127592E-03   4.928233E-04   1.054573E-03  -3.811515E-05  -9.466232E-04   4.320025E-04  -4.919150E-05  -1.753521E-04
  b1g   9   3.973599E-03  -1.193975E-04   4.693142E-03  -1.292823E-04  -8.792993E-04  -4.174825E-04  -1.170310E-04  -1.121709E-03
  b1g  10  -1.193975E-04   3.887262E-03   1.179904E-03  -5.746038E-05  -9.434028E-04   2.813542E-03  -9.196306E-05   2.576513E-03
  b1g  11   4.693142E-03   1.179904E-03   7.398652E-03  -1.697253E-04  -7.355051E-04   3.344562E-04  -1.089272E-04  -4.732998E-05
  b1g  12  -1.292823E-04  -5.746038E-05  -1.697253E-04   9.764015E-04   2.634362E-05  -2.715003E-05   2.372607E-03   9.597812E-06
  b1g  13  -8.792993E-04  -9.434028E-04  -7.355051E-04   2.634362E-05   1.464624E-03  -9.827092E-04   3.389561E-05  -1.919799E-04
  b1g  14  -4.174825E-04   2.813542E-03   3.344562E-04  -2.715003E-05  -9.827092E-04   2.296169E-03  -5.861216E-05   2.042855E-03
  b1g  15  -1.170310E-04  -9.196306E-05  -1.089272E-04   2.372607E-03   3.389561E-05  -5.861216E-05   6.439355E-03   1.942224E-05
  b1g  16  -1.121709E-03   2.576513E-03  -4.732998E-05   9.597812E-06  -1.919799E-04   2.042855E-03   1.942224E-05   2.655982E-03
  b1g  17  -2.144245E-03  -7.224531E-04  -4.557713E-03   1.140730E-04   3.458855E-04  -4.110048E-04   9.580074E-05  -4.900052E-04
  b1g  18   7.378324E-05   4.255136E-05   8.577572E-05  -1.149504E-03  -2.078103E-05   3.023311E-05  -3.825957E-03  -1.699679E-05
  b1g  19   5.728529E-05   9.543366E-04   9.025476E-04  -1.849689E-05  -4.947219E-05   7.235521E-04  -5.664655E-06   1.017442E-03
  b1g  20  -9.012898E-05   1.143145E-03  -3.822645E-04   4.795777E-06  -1.149900E-03   1.275479E-03  -6.373756E-06   7.658138E-04
  b1g  21  -6.971203E-05  -3.454418E-04  -2.792784E-04   1.375547E-05  -2.294012E-04  -6.948897E-05   2.881666E-05   1.724812E-04
  b1g  22   6.242658E-06   9.642018E-06   2.940059E-06  -2.396518E-04  -4.360532E-06   9.494348E-06  -1.158973E-03  -2.527134E-06
  b1g  23   7.955350E-04  -4.430480E-05   9.942826E-05   3.929677E-06  -3.983684E-04  -1.378874E-04   2.478126E-05  -6.911255E-04
  b1g  24  -6.780394E-04   5.155934E-04  -9.622659E-04   1.407449E-05   3.464969E-04   2.954547E-04  -1.018651E-05   7.240443E-04
  b1g  25  -8.388210E-06   1.602823E-06  -2.059264E-05  -2.175131E-04   3.125697E-06   1.785595E-07  -3.706620E-05   1.819334E-06
  b1g  26  -1.954525E-04  -1.588467E-05  -4.344609E-04   1.307110E-05  -3.583344E-06   1.115970E-04   1.204991E-05   1.711356E-04
  b1g  27  -3.267628E-04   5.678050E-05  -3.668022E-04   2.578384E-06   2.978443E-05   1.778582E-04  -1.062781E-05   3.100883E-04
  b1g  28  -2.962519E-04   4.559071E-04  -2.281717E-04   1.126730E-05  -2.723479E-05   3.432348E-04   2.429033E-05   4.105070E-04
  b1g  29  -6.919825E-06   1.131336E-07  -7.359056E-06   2.688880E-05   3.253881E-06  -6.617815E-07   2.615987E-04   4.729485E-06
  b1g  30   8.564841E-05   2.649897E-04   1.771732E-04  -9.702504E-06  -1.388826E-04   2.071223E-04  -1.952080E-05   2.754911E-04

               b1g  17        b1g  18        b1g  19        b1g  20        b1g  21        b1g  22        b1g  23        b1g  24
  b1g   2   2.544113E-03  -1.211530E-03   2.601879E-03   1.570853E-03  -1.208792E-04  -6.005193E-04   2.462747E-03  -5.981434E-05
  b1g   3  -9.293355E-04   2.557501E-03   3.209644E-04   1.642231E-03  -4.549432E-03   6.683967E-04  -4.903062E-03  -3.901243E-05
  b1g   4   9.173308E-05  -3.162660E-06   1.332076E-04   3.470905E-06  -1.284244E-04  -1.921646E-07  -1.643866E-04   1.553336E-04
  b1g   5   4.194558E-04  -2.293256E-05   1.940603E-04  -9.514594E-05  -1.379031E-04  -2.848404E-06  -4.498140E-04   3.263167E-04
  b1g   6   6.462819E-05   8.536982E-06   3.023462E-04   2.790655E-04  -2.443093E-04   2.858069E-06  -1.168055E-04   3.010251E-04
  b1g   7   2.367528E-06  -6.897617E-05  -4.354163E-08   6.200491E-07   2.137231E-07  -3.647486E-06   3.584750E-07   1.196319E-07
  b1g   8  -3.650891E-04   2.932113E-05  -3.310128E-05   5.522357E-04   4.890012E-05   5.956021E-06   4.153975E-04  -2.428547E-04
  b1g   9  -2.144245E-03   7.378324E-05   5.728529E-05  -9.012898E-05  -6.971203E-05   6.242658E-06   7.955350E-04  -6.780394E-04
  b1g  10  -7.224531E-04   4.255136E-05   9.543366E-04   1.143145E-03  -3.454418E-04   9.642018E-06  -4.430480E-05   5.155934E-04
  b1g  11  -4.557713E-03   8.577572E-05   9.025476E-04  -3.822645E-04  -2.792784E-04   2.940059E-06   9.942826E-05  -9.622659E-04
  b1g  12   1.140730E-04  -1.149504E-03  -1.849689E-05   4.795777E-06   1.375547E-05  -2.396518E-04   3.929677E-06   1.407449E-05
  b1g  13   3.458855E-04  -2.078103E-05  -4.947219E-05  -1.149900E-03  -2.294012E-04  -4.360532E-06  -3.983684E-04   3.464969E-04
  b1g  14  -4.110048E-04   3.023311E-05   7.235521E-04   1.275479E-03  -6.948897E-05   9.494348E-06  -1.378874E-04   2.954547E-04
  b1g  15   9.580074E-05  -3.825957E-03  -5.664655E-06  -6.373756E-06   2.881666E-05  -1.158973E-03   2.478126E-05  -1.018651E-05
  b1g  16  -4.900052E-04  -1.699679E-05   1.017442E-03   7.658138E-04   1.724812E-04  -2.527134E-06  -6.911255E-04   7.240443E-04
  b1g  17   3.966503E-03  -7.795862E-05  -1.012909E-03   3.163714E-04   6.737377E-06  -5.810711E-06   9.523527E-04   1.057816E-03
  b1g  18  -7.795862E-05   3.005950E-03   1.495444E-06   9.087210E-07  -2.277108E-05   1.231474E-03  -2.779214E-05  -8.770125E-07
  b1g  19  -1.012909E-03   1.495444E-06   8.694288E-04   1.084729E-04   3.953340E-05  -2.533705E-06  -3.441326E-04  -2.725345E-05
  b1g  20   3.163714E-04   9.087210E-07   1.084729E-04   1.512883E-03   2.818324E-04   3.438144E-06   3.194467E-04  -1.653326E-05
  b1g  21   6.737377E-06  -2.277108E-05   3.953340E-05   2.818324E-04   7.310912E-04  -5.263597E-06  -8.905401E-05   3.964960E-05
  b1g  22  -5.810711E-06   1.231474E-03  -2.533705E-06   3.438144E-06  -5.263597E-06   6.707798E-04  -7.729185E-06   4.836631E-06
  b1g  23   9.523527E-04  -2.779214E-05  -3.441326E-04   3.194467E-04  -8.905401E-05  -7.729185E-06   1.471598E-03   1.262380E-04
  b1g  24   1.057816E-03  -8.770125E-07  -2.725345E-05  -1.653326E-05   3.964960E-05   4.836631E-06   1.262380E-04   1.026696E-03
  b1g  25   1.365995E-05  -6.065195E-04   2.331910E-06   8.431449E-08   5.443081E-06  -4.686286E-04   8.428417E-06   1.787804E-06
  b1g  26   5.629544E-04  -1.033379E-05  -2.293343E-04   4.133931E-04   2.713899E-04   6.045128E-07   3.093135E-04   3.206941E-04
  b1g  27   2.915281E-05   5.388140E-06   2.525850E-04   1.797528E-04   1.259062E-04   2.942705E-06  -3.193767E-04   2.938243E-05
  b1g  28   1.339620E-04  -2.084758E-05   2.307949E-04   7.962305E-05  -6.731832E-05  -7.485269E-06   1.905687E-04   8.680374E-05
  b1g  29   4.440337E-06  -4.241180E-04   2.657812E-06  -2.005573E-06   3.714292E-06  -3.200315E-04   9.381317E-07  -1.407482E-07
  b1g  30  -1.638490E-04   1.331308E-05   1.455310E-04   4.000057E-05  -1.065386E-05   3.028353E-06  -1.937125E-04   1.080763E-04

               b1g  25        b1g  26        b1g  27        b1g  28        b1g  29        b1g  30
  b1g   2   3.489271E-04   2.782152E-03   1.468095E-03   3.743893E-03   3.405207E-04  -9.856928E-04
  b1g   3  -2.998334E-04   6.284950E-04   4.008657E-03  -3.059447E-03  -1.656242E-04   1.119757E-03
  b1g   4   1.199615E-06   4.965222E-06   2.672058E-05   9.913858E-05   1.258020E-06  -2.259402E-06
  b1g   5   2.744226E-06   4.517737E-05   1.073609E-04   1.763034E-04   3.231384E-06  -3.006754E-05
  b1g   6   2.423346E-06   6.568466E-06   2.506046E-05   2.066129E-04   1.360975E-06   5.640438E-05
  b1g   7  -3.672184E-05   4.724273E-07  -3.617463E-07   5.756980E-07  -4.868060E-06  -2.593157E-07
  b1g   8  -2.772164E-06  -4.130204E-05  -9.477619E-05  -4.058554E-05  -3.433215E-06   8.302839E-05
  b1g   9  -8.388210E-06  -1.954525E-04  -3.267628E-04  -2.962519E-04  -6.919825E-06   8.564841E-05
  b1g  10   1.602823E-06  -1.588467E-05   5.678050E-05   4.559071E-04   1.131336E-07   2.649897E-04
  b1g  11  -2.059264E-05  -4.344609E-04  -3.668022E-04  -2.281717E-04  -7.359056E-06   1.771732E-04
  b1g  12  -2.175131E-04   1.307110E-05   2.578384E-06   1.126730E-05   2.688880E-05  -9.702504E-06
  b1g  13   3.125697E-06  -3.583344E-06   2.978443E-05  -2.723479E-05   3.253881E-06  -1.388826E-04
  b1g  14   1.785595E-07   1.115970E-04   1.778582E-04   3.432348E-04  -6.617815E-07   2.071223E-04
  b1g  15  -3.706620E-05   1.204991E-05  -1.062781E-05   2.429033E-05   2.615987E-04  -1.952080E-05
  b1g  16   1.819334E-06   1.711356E-04   3.100883E-04   4.105070E-04   4.729485E-06   2.754911E-04
  b1g  17   1.365995E-05   5.629544E-04   2.915281E-05   1.339620E-04   4.440337E-06  -1.638490E-04
  b1g  18  -6.065195E-04  -1.033379E-05   5.388140E-06  -2.084758E-05  -4.241180E-04   1.331308E-05
  b1g  19   2.331910E-06  -2.293343E-04   2.525850E-04   2.307949E-04   2.657812E-06   1.455310E-04
  b1g  20   8.431449E-08   4.133931E-04   1.797528E-04   7.962305E-05  -2.005573E-06   4.000057E-05
  b1g  21   5.443081E-06   2.713899E-04   1.259062E-04  -6.731832E-05   3.714292E-06  -1.065386E-05
  b1g  22  -4.686286E-04   6.045128E-07   2.942705E-06  -7.485269E-06  -3.200315E-04   3.028353E-06
  b1g  23   8.428417E-06   3.093135E-04  -3.193767E-04   1.905687E-04   9.381317E-07  -1.937125E-04
  b1g  24   1.787804E-06   3.206941E-04   2.938243E-05   8.680374E-05  -1.407482E-07   1.080763E-04
  b1g  25   6.143781E-04   1.998584E-06   3.911719E-07   4.511894E-06   2.976547E-04  -2.747475E-06
  b1g  26   1.998584E-06   8.580018E-04   1.786465E-05  -1.215625E-04   2.744375E-07  -2.303013E-04
  b1g  27   3.911719E-07   1.786465E-05   4.036929E-04  -8.753703E-05   2.308441E-07   3.147069E-05
  b1g  28   4.511894E-06  -1.215625E-04  -8.753703E-05   4.979727E-04   4.503826E-06   1.250425E-04
  b1g  29   2.976547E-04   2.744375E-07   2.308441E-07   4.503826E-06   3.411392E-04  -9.019321E-07
  b1g  30  -2.747475E-06  -2.303013E-04   3.147069E-05   1.250425E-04  -9.019321E-07   3.212068E-04

Natural orbital populations,block 4
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     1.98116236     1.97280513     0.01453128     0.01023865     0.00998463     0.00458390     0.00221206
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00173560     0.00142171     0.00089896     0.00072347     0.00045231     0.00028293     0.00020978     0.00014851
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO    23       MO    24
  occ(*)=     0.00012666     0.00009463     0.00006007     0.00004622     0.00001943     0.00001766     0.00001215     0.00000693
              MO    25       MO    26       MO    27       MO    28       MO    29       MO    30
  occ(*)=     0.00000493     0.00000142     0.00000100     0.00000031     0.00000004     0.00000003

          modens reordered block   1

               b1u   1        b1u   2        b1u   3        b1u   4        b1u   5        b1u   6        b1u   7        b1u   8
  b1u   1    1.93109      -6.894680E-02   4.171020E-02  -1.688664E-03  -5.033487E-03   2.581489E-05  -2.357125E-03  -3.640860E-04
  b1u   2  -6.894680E-02   0.498930      -0.362185      -8.282072E-04  -1.976185E-03  -4.200317E-04  -1.410929E-03   1.467413E-03
  b1u   3   4.171020E-02  -0.362185       0.357490      -5.524818E-03   1.847090E-03  -4.700353E-03   1.947917E-03   3.451492E-03
  b1u   4  -1.688664E-03  -8.282072E-04  -5.524818E-03   7.677158E-04   6.121597E-06   6.732199E-04   5.790721E-05  -6.768810E-04
  b1u   5  -5.033487E-03  -1.976185E-03   1.847090E-03   6.121597E-06   2.041216E-03  -1.998918E-05  -1.582856E-03   1.284756E-05
  b1u   6   2.581489E-05  -4.200317E-04  -4.700353E-03   6.732199E-04  -1.998918E-05   9.540258E-04   1.691040E-05  -4.850964E-04
  b1u   7  -2.357125E-03  -1.410929E-03   1.947917E-03   5.790721E-05  -1.582856E-03   1.691040E-05   2.126979E-03  -3.047575E-05
  b1u   8  -3.640860E-04   1.467413E-03   3.451492E-03  -6.768810E-04   1.284756E-05  -4.850964E-04  -3.047575E-05   1.377971E-03
  b1u   9   1.584586E-02   1.976688E-03  -3.048589E-03  -2.940693E-04  -5.639156E-04  -1.409781E-04  -6.952140E-04   1.289407E-04
  b1u  10   3.257675E-04  -1.475346E-03  -3.247197E-03   6.223789E-04  -1.305573E-06   1.116821E-03   2.287689E-05  -1.290510E-04
  b1u  11  -1.505774E-04   3.550586E-03   1.219037E-03  -5.636185E-04  -1.017501E-05  -2.578738E-04  -6.632769E-05   9.735696E-04
  b1u  12   9.737946E-04   1.689278E-04   4.288051E-05   5.459051E-05   1.685298E-03   3.058505E-05  -1.378052E-03  -2.079390E-05
  b1u  13  -6.149196E-05   5.531069E-04   7.479808E-04  -1.303509E-04  -6.128499E-06  -3.909815E-04   8.968719E-06  -8.619192E-04
  b1u  14  -2.400254E-04  -2.511633E-04   6.936983E-04  -1.862153E-04  -1.106482E-05  -6.830372E-04   2.419788E-06   3.400163E-04
  b1u  15   1.596739E-03   2.880130E-04  -3.168497E-04  -2.224106E-07   2.449426E-04   3.097793E-06  -6.074069E-04   4.032585E-06
  b1u  16   5.738619E-05   1.181785E-04  -2.062577E-04   3.599715E-05   1.456261E-05   6.832973E-05  -2.736467E-05  -1.532777E-04

               b1u   9        b1u  10        b1u  11        b1u  12        b1u  13        b1u  14        b1u  15        b1u  16
  b1u   1   1.584586E-02   3.257675E-04  -1.505774E-04   9.737946E-04  -6.149196E-05  -2.400254E-04   1.596739E-03   5.738619E-05
  b1u   2   1.976688E-03  -1.475346E-03   3.550586E-03   1.689278E-04   5.531069E-04  -2.511633E-04   2.880130E-04   1.181785E-04
  b1u   3  -3.048589E-03  -3.247197E-03   1.219037E-03   4.288051E-05   7.479808E-04   6.936983E-04  -3.168497E-04  -2.062577E-04
  b1u   4  -2.940693E-04   6.223789E-04  -5.636185E-04   5.459051E-05  -1.303509E-04  -1.862153E-04  -2.224106E-07   3.599715E-05
  b1u   5  -5.639156E-04  -1.305573E-06  -1.017501E-05   1.685298E-03  -6.128499E-06  -1.106482E-05   2.449426E-04   1.456261E-05
  b1u   6  -1.409781E-04   1.116821E-03  -2.578738E-04   3.058505E-05  -3.909815E-04  -6.830372E-04   3.097793E-06   6.832973E-05
  b1u   7  -6.952140E-04   2.287689E-05  -6.632769E-05  -1.378052E-03   8.968719E-06   2.419788E-06  -6.074069E-04  -2.736467E-05
  b1u   8   1.289407E-04  -1.290510E-04   9.735696E-04  -2.079390E-05  -8.619192E-04   3.400163E-04   4.032585E-06  -1.532777E-04
  b1u   9   2.786288E-03  -1.358394E-04   1.324808E-04  -4.711183E-04  -7.735454E-06   6.827650E-05   6.047473E-05  -1.042670E-05
  b1u  10  -1.358394E-04   1.620296E-03   5.522482E-05   2.811557E-05  -1.063491E-03  -9.127831E-04   2.304782E-06   1.850001E-04
  b1u  11   1.324808E-04   5.522482E-05   1.223457E-03   1.365683E-05  -6.729441E-04  -1.538076E-04   2.574227E-05   8.854867E-05
  b1u  12  -4.711183E-04   2.811557E-05   1.365683E-05   2.008270E-03  -8.316130E-06  -2.871421E-05   2.730152E-04   2.776821E-05
  b1u  13  -7.735454E-06  -1.063491E-03  -6.729441E-04  -8.316130E-06   2.005609E-03  -1.629811E-05  -4.965272E-06  -4.936267E-05
  b1u  14   6.827650E-05  -9.127831E-04  -1.538076E-04  -2.871421E-05  -1.629811E-05   1.233729E-03  -8.591595E-06  -2.941593E-04
  b1u  15   6.047473E-05   2.304782E-06   2.574227E-05   2.730152E-04  -4.965272E-06  -8.591595E-06   5.835670E-04   2.141319E-05
  b1u  16  -1.042670E-05   1.850001E-04   8.854867E-05   2.776821E-05  -4.936267E-05  -2.941593E-04   2.141319E-05   4.480343E-04

Natural orbital populations,block 5
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.93688450     0.79172872     0.06043451     0.00526308     0.00372516     0.00335597     0.00240538     0.00101984
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00064381     0.00044737     0.00037391     0.00019076     0.00010601     0.00006396     0.00003458     0.00000843

          modens reordered block   1

               b2g   1        b2g   2        b2g   3        b2g   4        b2g   5        b2g   6        b2g   7        b2g   8
  b2g   1    1.72740       1.149870E-02   1.028344E-03   8.386215E-04  -5.941985E-03  -2.302695E-03  -9.278537E-03  -3.555333E-03
  b2g   2   1.149870E-02   1.510091E-03   1.715342E-04  -4.482123E-04  -9.489450E-04   1.425448E-04  -1.262379E-03  -2.860128E-04
  b2g   3   1.028344E-03   1.715342E-04   1.010364E-03  -3.647906E-03  -3.793931E-05   1.486340E-03  -5.809341E-05  -8.298351E-05
  b2g   4   8.386215E-04  -4.482123E-04  -3.647906E-03   1.383837E-02   1.056114E-04  -6.104731E-03   2.285881E-04   3.526380E-04
  b2g   5  -5.941985E-03  -9.489450E-04  -3.793931E-05   1.056114E-04   1.284863E-03  -3.685915E-05   1.377779E-03  -5.497547E-04
  b2g   6  -2.302695E-03   1.425448E-04   1.486340E-03  -6.104731E-03  -3.685915E-05   3.133461E-03  -1.177454E-04  -1.870417E-04
  b2g   7  -9.278537E-03  -1.262379E-03  -5.809341E-05   2.285881E-04   1.377779E-03  -1.177454E-04   1.799096E-03  -2.360659E-04
  b2g   8  -3.555333E-03  -2.860128E-04  -8.298351E-05   3.526380E-04  -5.497547E-04  -1.870417E-04  -2.360659E-04   1.356350E-03
  b2g   9  -3.023167E-03   4.741531E-05   9.112094E-04  -4.072451E-03  -2.110694E-05   2.413192E-03  -7.315140E-05  -1.162795E-04
  b2g  10  -1.156271E-02  -1.187519E-03  -1.171001E-04   4.477147E-04   5.765380E-04  -2.319259E-04   9.547205E-04   7.789616E-04
  b2g  11  -7.109364E-04   3.940276E-05   5.007365E-04  -1.915314E-03  -1.735155E-05   7.622598E-04  -3.363987E-05  -4.297973E-05
  b2g  12  -6.580489E-05  -4.082578E-04   2.686337E-05  -9.628509E-05   1.118792E-03   4.412851E-05   1.125494E-03  -1.184034E-03
  b2g  13  -8.416311E-04   1.721666E-05  -8.363575E-05  -1.274909E-04  -2.144631E-05   6.338715E-04  -4.381853E-05  -3.749295E-05
  b2g  14  -1.262132E-03  -4.256112E-05  -4.190236E-05   1.742214E-04   6.213068E-05  -9.294871E-05  -6.108203E-05   2.495068E-04
  b2g  15  -1.079725E-03  -6.472630E-05  -3.363630E-05   1.670325E-04   1.064222E-04  -9.805418E-05   3.617018E-04   2.323736E-04
  b2g  16  -2.682019E-04  -4.519152E-06  -2.868330E-05   5.172823E-05  -2.268401E-06   4.553552E-05  -8.052547E-06  -1.024500E-06

               b2g   9        b2g  10        b2g  11        b2g  12        b2g  13        b2g  14        b2g  15        b2g  16
  b2g   1  -3.023167E-03  -1.156271E-02  -7.109364E-04  -6.580489E-05  -8.416311E-04  -1.262132E-03  -1.079725E-03  -2.682019E-04
  b2g   2   4.741531E-05  -1.187519E-03   3.940276E-05  -4.082578E-04   1.721666E-05  -4.256112E-05  -6.472630E-05  -4.519152E-06
  b2g   3   9.112094E-04  -1.171001E-04   5.007365E-04   2.686337E-05  -8.363575E-05  -4.190236E-05  -3.363630E-05  -2.868330E-05
  b2g   4  -4.072451E-03   4.477147E-04  -1.915314E-03  -9.628509E-05  -1.274909E-04   1.742214E-04   1.670325E-04   5.172823E-05
  b2g   5  -2.110694E-05   5.765380E-04  -1.735155E-05   1.118792E-03  -2.144631E-05   6.213068E-05   1.064222E-04  -2.268401E-06
  b2g   6   2.413192E-03  -2.319259E-04   7.622598E-04   4.412851E-05   6.338715E-04  -9.294871E-05  -9.805418E-05   4.553552E-05
  b2g   7  -7.315140E-05   9.547205E-04  -3.363987E-05   1.125494E-03  -4.381853E-05  -6.108203E-05   3.617018E-04  -8.052547E-06
  b2g   8  -1.162795E-04   7.789616E-04  -4.297973E-05  -1.184034E-03  -3.749295E-05   2.495068E-04   2.323736E-04  -1.024500E-06
  b2g   9   2.123830E-03  -1.414698E-04   4.808602E-04   1.931284E-05   9.353971E-04  -6.911871E-05  -8.668329E-05   1.581696E-04
  b2g  10  -1.414698E-04   1.669188E-03  -6.579149E-05  -6.441005E-05  -7.381767E-05   2.142452E-04  -7.816835E-05  -4.262208E-06
  b2g  11   4.808602E-04  -6.579149E-05   5.606173E-04   8.431390E-06  -2.382584E-04  -2.134580E-05  -1.928883E-05  -5.706298E-05
  b2g  12   1.931284E-05  -6.441005E-05   8.431390E-06   1.706412E-03  -3.798496E-06  -2.452888E-04  -1.332100E-05  -9.220835E-06
  b2g  13   9.353971E-04  -7.381767E-05  -2.382584E-04  -3.798496E-06   1.106268E-03  -2.309234E-05  -3.349001E-05   2.914160E-04
  b2g  14  -6.911871E-05   2.142452E-04  -2.134580E-05  -2.452888E-04  -2.309234E-05   6.105084E-04   1.670568E-04   3.243127E-06
  b2g  15  -8.668329E-05  -7.816835E-05  -1.928883E-05  -1.332100E-05  -3.349001E-05   1.670568E-04   8.755394E-04  -1.316605E-05
  b2g  16   1.581696E-04  -4.262208E-06  -5.706298E-05  -9.220835E-06   2.914160E-04   3.243127E-06  -1.316605E-05   3.381653E-04

Natural orbital populations,block 6
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.72764697     0.01923788     0.00518300     0.00316780     0.00221418     0.00106976     0.00058296     0.00038886
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.00029595     0.00025840     0.00013043     0.00007029     0.00005249     0.00001677     0.00000808     0.00000240

          modens reordered block   1

               b3g   1        b3g   2        b3g   3        b3g   4        b3g   5        b3g   6        b3g   7        b3g   8
  b3g   1    1.18015      -2.019087E-02  -7.815709E-03  -3.517389E-03  -1.088612E-02  -4.474975E-03  -1.297234E-02  -1.617461E-05
  b3g   2  -2.019087E-02   1.658630E-03   1.054467E-03   2.143877E-04   1.206889E-03   1.388420E-04   1.014962E-03  -7.262990E-05
  b3g   3  -7.815709E-03   1.054467E-03   1.293511E-03   5.481145E-05   1.297950E-03  -5.539314E-04   5.232600E-04  -2.078905E-05
  b3g   4  -3.517389E-03   2.143877E-04   5.481145E-05   2.564394E-03   3.617158E-05   6.801429E-06   7.959694E-05  -2.614625E-03
  b3g   5  -1.088612E-02   1.206889E-03   1.297950E-03   3.617158E-05   1.639049E-03  -2.564996E-04   8.267097E-04  -9.721292E-06
  b3g   6  -4.474975E-03   1.388420E-04  -5.539314E-04   6.801429E-06  -2.564996E-04   1.262268E-03   6.367275E-04   9.759090E-06
  b3g   7  -1.297234E-02   1.014962E-03   5.232600E-04   7.959694E-05   8.267097E-04   6.367275E-04   1.404684E-03  -6.908453E-05
  b3g   8  -1.617461E-05  -7.262990E-05  -2.078905E-05  -2.614625E-03  -9.721292E-06   9.759090E-06  -6.908453E-05   3.408288E-03
  b3g   9  -3.896111E-04   4.079040E-04   1.065855E-03   6.326423E-06   1.051696E-03  -1.136064E-03  -4.123943E-05  -1.526999E-05
  b3g  10  -9.702725E-04   1.355561E-05   7.025085E-05  -9.617250E-06  -6.248419E-05   2.083912E-04   1.614067E-04   1.263856E-05
  b3g  11  -2.282614E-04   1.129445E-05   6.754615E-05  -4.906958E-05   3.098116E-04   2.535130E-04  -9.696908E-05   5.224014E-05

               b3g   9        b3g  10        b3g  11
  b3g   1  -3.896111E-04  -9.702725E-04  -2.282614E-04
  b3g   2   4.079040E-04   1.355561E-05   1.129445E-05
  b3g   3   1.065855E-03   7.025085E-05   6.754615E-05
  b3g   4   6.326423E-06  -9.617250E-06  -4.906958E-05
  b3g   5   1.051696E-03  -6.248419E-05   3.098116E-04
  b3g   6  -1.136064E-03   2.083912E-04   2.535130E-04
  b3g   7  -4.123943E-05   1.614067E-04  -9.696908E-05
  b3g   8  -1.526999E-05   1.263856E-05   5.224014E-05
  b3g   9   1.638905E-03  -2.149004E-04  -4.876428E-05
  b3g  10  -2.149004E-04   5.843639E-04   1.852995E-04
  b3g  11  -4.876428E-05   1.852995E-04   8.893873E-04

Natural orbital populations,block 7
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     1.18082538     0.00564998     0.00476831     0.00264839     0.00108216     0.00056251     0.00046849     0.00030604
              MO     9       MO    10       MO    11
  occ(*)=     0.00010620     0.00005975     0.00002102

          modens reordered block   1

               au    1        au    2        au    3        au    4        au    5        au    6        au    7        au    8
  au    1   0.249960       2.769586E-03   9.800658E-04  -7.842998E-04   9.840037E-04   1.641316E-04  -2.163480E-03  -1.585437E-03
  au    2   2.769586E-03   1.079446E-03   8.853866E-04  -6.915806E-04  -7.987710E-04   1.037579E-05   5.424755E-04   2.586199E-04
  au    3   9.800658E-04   8.853866E-04   1.093314E-03  -4.810462E-04  -1.230528E-03   3.735436E-06   2.215341E-04   4.962006E-04
  au    4  -7.842998E-04  -6.915806E-04  -4.810462E-04   1.280037E-03   1.327025E-04  -1.784578E-06  -8.685455E-04   7.916172E-04
  au    5   9.840037E-04  -7.987710E-04  -1.230528E-03   1.327025E-04   1.723900E-03  -4.767585E-06   7.623066E-05  -1.163705E-03
  au    6   1.641316E-04   1.037579E-05   3.735436E-06  -1.784578E-06  -4.767585E-06   7.839602E-05  -2.929767E-06   1.064819E-05
  au    7  -2.163480E-03   5.424755E-04   2.215341E-04  -8.685455E-04   7.623066E-05  -2.929767E-06   1.112836E-03  -5.952540E-04
  au    8  -1.585437E-03   2.586199E-04   4.962006E-04   7.916172E-04  -1.163705E-03   1.064819E-05  -5.952540E-04   2.028186E-03
  au    9  -6.681568E-04  -1.702512E-04  -6.632950E-04   3.077800E-04   8.907310E-04   1.611921E-06   1.975358E-04  -7.985564E-06
  au   10   2.701227E-05   3.568501E-05   6.572444E-05  -1.650987E-04  -1.771259E-04  -6.303596E-07  -8.279401E-05   3.654897E-05
  au   11   4.167997E-04   1.486744E-05   3.994278E-06   8.068503E-06  -4.407384E-06   1.520296E-04  -1.183775E-05   2.501056E-05

               au    9        au   10        au   11
  au    1  -6.681568E-04   2.701227E-05   4.167997E-04
  au    2  -1.702512E-04   3.568501E-05   1.486744E-05
  au    3  -6.632950E-04   6.572444E-05   3.994278E-06
  au    4   3.077800E-04  -1.650987E-04   8.068503E-06
  au    5   8.907310E-04  -1.771259E-04  -4.407384E-06
  au    6   1.611921E-06  -6.303596E-07   1.520296E-04
  au    7   1.975358E-04  -8.279401E-05  -1.183775E-05
  au    8  -7.985564E-06   3.654897E-05   2.501056E-05
  au    9   1.207125E-03  -2.886689E-04   1.244726E-05
  au   10  -2.886689E-04   4.378274E-04  -4.611418E-06
  au   11   1.244726E-05  -4.611418E-06   5.334484E-04

Natural orbital populations,block 8
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     0.25003291     0.00436742     0.00323529     0.00141410     0.00057808     0.00044013     0.00024348     0.00013424
              MO     9       MO    10       MO    11
  occ(*)=     0.00004575     0.00003210     0.00001150


 total number of electrons =   42.0000000000



          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        ag  partial gross atomic populations
   ao class       1ag        2ag        3ag        4ag        5ag        6ag 
    C1_ s       0.024115   1.965539   0.505563   0.485423  -0.009445   0.053012
    C1_ p       0.003468  -0.000704   0.029021   0.125352   0.314771   0.303780
    C1_ d       0.000470   0.000179   0.004608  -0.109464   0.035364  -0.027948
    C2_ s       1.965551   0.024164   1.114118   0.262360  -0.030171   0.019351
    C2_ p       0.000801   0.007425   0.089874   0.118303   0.754427   1.246798
    C2_ d      -0.000546   0.001050   0.012788  -0.054690   0.086225   0.008381
    H1_ s       0.001179   0.000814   0.083538   0.759762   0.235915   0.310373
    H1_ p       0.000351   0.000587   0.006624   0.031102   0.001354  -0.037382
    H2_ s       0.003058   0.000734   0.132790   0.360937   0.589262   0.124120
    H2_ p       0.001552   0.000211   0.008494   0.002246   0.001618  -0.027215
 
   ao class       7ag        8ag        9ag       10ag       11ag       12ag 
    C1_ s       0.004576   0.000477   0.000629  -0.000133  -0.000045   0.000469
    C1_ p       0.001972   0.000779   0.003737   0.001609   0.000013   0.000145
    C1_ d       0.000370   0.001128  -0.001269   0.000751   0.000698   0.000202
    C2_ s       0.002247   0.001026   0.000343  -0.000005   0.003810   0.000269
    C2_ p       0.004701   0.001301   0.002545   0.002831   0.000494   0.001978
    C2_ d       0.001245   0.002036  -0.000304   0.001558   0.001815   0.000876
    H1_ s       0.000084   0.002215   0.002740  -0.000000  -0.000298  -0.000114
    H1_ p      -0.000018  -0.000135   0.000166   0.000050   0.000022   0.000061
    H2_ s       0.000051   0.003934   0.001562   0.000078  -0.000904  -0.000127
    H2_ p       0.000101  -0.000204   0.000066   0.000046   0.000021   0.000769
 
   ao class      13ag       14ag       15ag       16ag       17ag       18ag 
    C1_ s      -0.000727   0.000919  -0.000234  -0.000016  -0.000491   0.000091
    C1_ p       0.001246   0.000010  -0.000041   0.000332  -0.000044   0.000193
    C1_ d       0.000867   0.000405   0.000534   0.000082   0.000015   0.000202
    C2_ s       0.000364   0.000880   0.000224   0.000048   0.000484   0.000083
    C2_ p       0.000497  -0.000046   0.000327   0.000715   0.000381   0.000042
    C2_ d       0.001079   0.000535   0.000583   0.000238   0.000446   0.000168
    H1_ s      -0.000766  -0.000681  -0.000095  -0.000099  -0.000012  -0.000034
    H1_ p       0.000241   0.000088   0.000187   0.000035   0.000012   0.000055
    H2_ s      -0.000709  -0.000293  -0.000231  -0.000147   0.000188  -0.000027
    H2_ p       0.000192  -0.000026   0.000226   0.000042   0.000043   0.000119
 
   ao class      19ag       20ag       21ag       22ag       23ag       24ag 
    C1_ s      -0.000065   0.000102   0.000076   0.000047  -0.000052  -0.000017
    C1_ p       0.000007   0.000170  -0.000036   0.000080   0.000009   0.000000
    C1_ d       0.000159   0.000061  -0.000010   0.000049   0.000137   0.000042
    C2_ s      -0.000052   0.000004   0.000002   0.000015  -0.000020  -0.000011
    C2_ p       0.000655  -0.000006   0.000064  -0.000003  -0.000000  -0.000061
    C2_ d      -0.000084   0.000047   0.000135   0.000111   0.000011   0.000187
    H1_ s      -0.000001  -0.000008   0.000013   0.000082   0.000061  -0.000001
    H1_ p       0.000032  -0.000037   0.000046  -0.000112  -0.000013   0.000019
    H2_ s       0.000003  -0.000032   0.000082   0.000062   0.000037   0.000016
    H2_ p       0.000095   0.000150   0.000050  -0.000051  -0.000003  -0.000039
 
   ao class      25ag       26ag       27ag       28ag       29ag       30ag 
    C1_ s      -0.000068  -0.000034  -0.000054   0.000027   0.000028  -0.000020
    C1_ p      -0.000008   0.000033   0.000004   0.000017  -0.000026   0.000016
    C1_ d       0.000007   0.000013   0.000018   0.000022   0.000001  -0.000001
    C2_ s       0.000013  -0.000106   0.000022   0.000020   0.000018   0.000033
    C2_ p       0.000017   0.000115  -0.000008  -0.000004  -0.000028  -0.000005
    C2_ d       0.000077   0.000013   0.000003  -0.000002   0.000080   0.000001
    H1_ s       0.000056  -0.000023   0.000043   0.000001   0.000014  -0.000010
    H1_ p      -0.000011   0.000006   0.000002   0.000001   0.000001   0.000000
    H2_ s       0.000028   0.000115   0.000072  -0.000014   0.000011  -0.000001
    H2_ p       0.000004  -0.000028  -0.000007   0.000013  -0.000051   0.000012
 
   ao class      31ag       32ag       33ag       34ag       35ag       36ag 
    C1_ s      -0.000002  -0.000017   0.000058   0.000013  -0.000001  -0.000000
    C1_ p      -0.000008   0.000010  -0.000026  -0.000016   0.000002   0.000000
    C1_ d       0.000008  -0.000003  -0.000001   0.000001   0.000001   0.000000
    C2_ s      -0.000003   0.000024  -0.000030   0.000010  -0.000001   0.000000
    C2_ p       0.000004  -0.000007   0.000009  -0.000017  -0.000002   0.000001
    C2_ d      -0.000005   0.000004   0.000000  -0.000001  -0.000009  -0.000000
    H1_ s       0.000023  -0.000009   0.000011   0.000001   0.000000   0.000001
    H1_ p      -0.000009   0.000011  -0.000000   0.000006   0.000003   0.000000
    H2_ s       0.000015  -0.000009  -0.000013   0.000001   0.000000  -0.000000
    H2_ p      -0.000003   0.000011   0.000003   0.000008   0.000012   0.000000
 
   ao class      37ag       38ag       39ag 
    C1_ s      -0.000001  -0.000001   0.000000
    C1_ p       0.000005   0.000001  -0.000000
    C1_ d      -0.000000  -0.000000   0.000000
    C2_ s      -0.000000  -0.000000   0.000000
    C2_ p      -0.000000   0.000001  -0.000000
    C2_ d      -0.000000  -0.000000   0.000000
    H1_ s      -0.000002  -0.000000   0.000000
    H1_ p       0.000000  -0.000000   0.000000
    H2_ s       0.000001  -0.000000   0.000000
    H2_ p       0.000000  -0.000000   0.000000

                        b3u partial gross atomic populations
   ao class       1b3u       2b3u       3b3u       4b3u       5b3u       6b3u
    C1_ s       0.022924   1.952003   0.714934   0.146636   0.013492   0.000997
    C1_ p      -0.000455   0.001641  -0.031294   0.152584   0.783053   0.001667
    C1_ d      -0.000102   0.000126   0.043749  -0.179947  -0.027249   0.001495
    C2_ s       1.958708   0.025923   0.400414   0.133372   0.096628  -0.000204
    C2_ p       0.013886   0.014358   0.189842   0.485507   0.509977   0.007789
    C2_ d       0.003923   0.004386  -0.030194  -0.387115  -0.036673   0.001935
    H1_ s       0.000617   0.000533   0.444798   0.467646   0.503914   0.000410
    H1_ p       0.000315   0.000354   0.040273  -0.000269  -0.018557  -0.000015
    H2_ s       0.000067   0.000410   0.207471   1.165144   0.151869   0.000614
    H2_ p       0.000116   0.000266   0.005100  -0.005070  -0.000153   0.000143
 
   ao class       7b3u       8b3u       9b3u      10b3u      11b3u      12b3u
    C1_ s       0.002950   0.000293   0.000244  -0.000852  -0.000398   0.000154
    C1_ p       0.000114   0.001490   0.002958   0.001788  -0.000074  -0.000085
    C1_ d      -0.000649   0.000274  -0.001299   0.000916   0.001288   0.000366
    C2_ s       0.006629   0.000859   0.000480  -0.000449  -0.000223   0.000081
    C2_ p       0.001636   0.003130   0.005635   0.005229   0.000628  -0.000583
    C2_ d      -0.000023   0.000057  -0.002539   0.001096   0.002179   0.001614
    H1_ s       0.001023   0.002856   0.001010  -0.002169  -0.000117   0.000046
    H1_ p       0.000254  -0.000121   0.000119  -0.000019   0.000128  -0.000011
    H2_ s       0.000512   0.001824   0.001806  -0.001108  -0.000049   0.000027
    H2_ p       0.000366   0.000028   0.000223   0.000418   0.000008   0.000481
 
   ao class      13b3u      14b3u      15b3u      16b3u      17b3u      18b3u
    C1_ s       0.000036   0.000703   0.000017  -0.000082  -0.000132   0.000167
    C1_ p       0.001314   0.000587   0.000075  -0.000148   0.000035   0.000184
    C1_ d      -0.000449   0.000272   0.000156   0.000253   0.000020   0.000156
    C2_ s       0.000877   0.000138  -0.000131   0.000050   0.000284   0.000004
    C2_ p       0.002484  -0.000126   0.000376  -0.000084   0.000227  -0.000164
    C2_ d      -0.000951   0.000121   0.000215   0.000641   0.000106   0.000057
    H1_ s      -0.000764  -0.000494  -0.000001   0.000051  -0.000059   0.000023
    H1_ p       0.000322   0.000123   0.000020  -0.000028   0.000027  -0.000032
    H2_ s      -0.001615  -0.000192   0.000023   0.000057  -0.000072  -0.000016
    H2_ p       0.000633   0.000104   0.000035  -0.000035  -0.000003   0.000046
 
   ao class      19b3u      20b3u      21b3u      22b3u      23b3u      24b3u
    C1_ s      -0.000271  -0.000060   0.000094   0.000165  -0.000009  -0.000087
    C1_ p      -0.000044  -0.000023   0.000041  -0.000015   0.000022   0.000026
    C1_ d       0.000104   0.000029   0.000007   0.000016   0.000080   0.000061
    C2_ s      -0.000353  -0.000007  -0.000254  -0.000000   0.000026  -0.000097
    C2_ p       0.000705   0.000266   0.000247   0.000192  -0.000054   0.000039
    C2_ d       0.000113  -0.000120   0.000039  -0.000318   0.000021   0.000101
    H1_ s       0.000002   0.000104   0.000053   0.000014   0.000002   0.000055
    H1_ p       0.000006  -0.000017  -0.000038   0.000013   0.000005  -0.000025
    H2_ s      -0.000000   0.000045   0.000133   0.000011  -0.000005   0.000020
    H2_ p       0.000083   0.000050  -0.000112   0.000059   0.000010  -0.000018
 
   ao class      25b3u      26b3u      27b3u      28b3u      29b3u      30b3u
    C1_ s       0.000052  -0.000008   0.000007  -0.000014   0.000011  -0.000009
    C1_ p       0.000001   0.000012  -0.000011  -0.000022  -0.000003   0.000006
    C1_ d       0.000016  -0.000007  -0.000008   0.000003   0.000001   0.000001
    C2_ s       0.000018   0.000133   0.000004   0.000014   0.000008   0.000001
    C2_ p      -0.000054  -0.000129  -0.000030  -0.000049  -0.000001   0.000007
    C2_ d       0.000023   0.000017   0.000084   0.000057   0.000002  -0.000005
    H1_ s      -0.000021   0.000004   0.000002   0.000038  -0.000001  -0.000006
    H1_ p       0.000014  -0.000000   0.000004  -0.000003   0.000001   0.000012
    H2_ s      -0.000009   0.000006   0.000004   0.000019  -0.000003  -0.000002
    H2_ p       0.000009   0.000004  -0.000027  -0.000016   0.000001   0.000007
 
   ao class      31b3u      32b3u      33b3u      34b3u      35b3u      36b3u
    C1_ s      -0.000007  -0.000006  -0.000003  -0.000010   0.000019  -0.000000
    C1_ p      -0.000000   0.000002   0.000001   0.000007  -0.000010   0.000000
    C1_ d      -0.000000   0.000004  -0.000001   0.000000  -0.000000  -0.000000
    C2_ s      -0.000003   0.000003   0.000001   0.000015  -0.000012  -0.000000
    C2_ p       0.000013  -0.000008  -0.000000  -0.000010   0.000005   0.000001
    C2_ d      -0.000001   0.000002  -0.000001  -0.000001  -0.000001  -0.000000
    H1_ s      -0.000001  -0.000000   0.000002  -0.000003   0.000002  -0.000000
    H1_ p       0.000001  -0.000001   0.000001   0.000000   0.000000   0.000000
    H2_ s      -0.000003   0.000013   0.000001   0.000002  -0.000003  -0.000000
    H2_ p       0.000010  -0.000003   0.000005   0.000001   0.000001   0.000000
 
   ao class      37b3u      38b3u      39b3u
    C1_ s      -0.000000   0.000000   0.000000
    C1_ p       0.000001  -0.000000  -0.000000
    C1_ d      -0.000000  -0.000000   0.000000
    C2_ s      -0.000001   0.000000   0.000000
    C2_ p       0.000001  -0.000000  -0.000000
    C2_ d      -0.000000  -0.000000   0.000000
    H1_ s      -0.000000   0.000000   0.000000
    H1_ p       0.000000  -0.000000  -0.000000
    H2_ s      -0.000001   0.000000   0.000000
    H2_ p       0.000000   0.000000  -0.000000

                        b2u partial gross atomic populations
   ao class       1b2u       2b2u       3b2u       4b2u       5b2u       6b2u
    C1_ p       0.016768   0.099746   0.055938   0.716085   0.005734   0.001148
    C1_ d       0.005055  -0.042057  -0.023711   0.019929   0.001223   0.000123
    C2_ s       1.965412   1.172075   0.113774   0.009971   0.000100   0.001235
    C2_ p       0.007796   0.004643   1.262508   1.228835   0.004444   0.003643
    C2_ d       0.001844   0.045224  -0.032859   0.030515   0.002146   0.000030
    H1_ p       0.000103  -0.010156   0.002791  -0.015484   0.000094   0.000017
    H2_ s       0.001851   0.661122   0.620798   0.012262   0.001259   0.004500
    H2_ p       0.001171   0.054112  -0.022833  -0.028228  -0.000023  -0.000089
 
   ao class       7b2u       8b2u       9b2u      10b2u      11b2u      12b2u
    C1_ p       0.001034   0.001690   0.000829  -0.000342  -0.000293   0.000404
    C1_ d       0.000828   0.000314   0.001115   0.001009  -0.000002  -0.000045
    C2_ s      -0.000008   0.000752  -0.001289   0.000164   0.000871  -0.000001
    C2_ p       0.001315   0.003650   0.000050  -0.000291   0.000619   0.000584
    C2_ d       0.001891   0.001318   0.002487   0.001035   0.000427   0.000020
    H1_ p       0.000478   0.000228  -0.000009   0.000317   0.000029   0.000054
    H2_ s      -0.000001  -0.003434  -0.000029   0.000050  -0.000582  -0.000001
    H2_ p       0.000891   0.000113   0.000121   0.000150   0.000169   0.000072
 
   ao class      13b2u      14b2u      15b2u      16b2u      17b2u      18b2u
    C1_ p       0.000248   0.000053   0.000150   0.000271   0.000125   0.000026
    C1_ d       0.000059   0.000130   0.000063   0.000025  -0.000114   0.000098
    C2_ s      -0.000002  -0.000008  -0.000280  -0.000308  -0.000003   0.000005
    C2_ p       0.000160   0.000025   0.000085   0.000120   0.000160  -0.000082
    C2_ d       0.000216   0.000203   0.000405   0.000138  -0.000052   0.000326
    H1_ p       0.000010   0.000083  -0.000001   0.000043   0.000064  -0.000046
    H2_ s       0.000023   0.000001  -0.000022   0.000126   0.000093   0.000001
    H2_ p       0.000032   0.000159   0.000065  -0.000036   0.000009  -0.000120
 
   ao class      19b2u      20b2u      21b2u      22b2u      23b2u      24b2u
    C1_ p       0.000080  -0.000014  -0.000082  -0.000027  -0.000009  -0.000006
    C1_ d      -0.000164  -0.000005   0.000005   0.000071   0.000035  -0.000015
    C2_ s       0.000116   0.000020   0.000144   0.000015   0.000000   0.000001
    C2_ p       0.000129   0.000004  -0.000072  -0.000028  -0.000042  -0.000003
    C2_ d      -0.000111   0.000081   0.000015   0.000021   0.000005   0.000001
    H1_ p       0.000032   0.000001   0.000000  -0.000021  -0.000011   0.000016
    H2_ s       0.000023  -0.000001   0.000032   0.000003   0.000057   0.000003
    H2_ p       0.000029   0.000017  -0.000001  -0.000003  -0.000009   0.000024
 
   ao class      25b2u      26b2u      27b2u      28b2u      29b2u      30b2u
    C1_ p       0.000003   0.000014  -0.000000   0.000001   0.000000  -0.000000
    C1_ d      -0.000003  -0.000005  -0.000000  -0.000001  -0.000000  -0.000000
    C2_ s      -0.000005  -0.000017   0.000001  -0.000001   0.000000   0.000000
    C2_ p       0.000004   0.000021   0.000004   0.000003   0.000000  -0.000000
    C2_ d       0.000000  -0.000002  -0.000000  -0.000001  -0.000000  -0.000000
    H1_ p      -0.000000   0.000007   0.000000   0.000001  -0.000000   0.000000
    H2_ s      -0.000008  -0.000011   0.000000  -0.000002  -0.000000   0.000000
    H2_ p       0.000019   0.000002  -0.000000   0.000001   0.000000  -0.000000

                        b1g partial gross atomic populations
   ao class       1b1g       2b1g       3b1g       4b1g       5b1g       6b1g
    C1_ p       0.003407   0.008964   0.760400   0.003011  -0.000270   0.002785
    C1_ d       0.000348   0.001790   0.015891   0.000537   0.000208   0.000351
    C2_ s       1.989749   0.746745   0.105762   0.006183   0.001113   0.000003
    C2_ p       0.004512   0.190054   0.811237   0.003535   0.006531   0.005293
    C2_ d       0.001538  -0.177806  -0.027339   0.001124  -0.001763   0.001233
    H1_ p       0.000044  -0.008886  -0.008614   0.000103   0.000001   0.000053
    H2_ s       0.000202   1.176919   0.372535   0.000096   0.004202   0.000163
    H2_ p       0.000200   0.043383  -0.057066  -0.000058   0.000217   0.000105
 
   ao class       7b1g       8b1g       9b1g      10b1g      11b1g      12b1g
    C1_ p       0.001181   0.000033  -0.000274   0.000285  -0.000023   0.000466
    C1_ d       0.000659   0.000553   0.000676   0.000160   0.000056  -0.000130
    C2_ s       0.000663  -0.000502  -0.000000  -0.000350   0.000123  -0.000073
    C2_ p       0.000961   0.001933  -0.000555   0.000031   0.000207   0.000126
    C2_ d       0.000450   0.001304   0.001435   0.001159   0.000359   0.000232
    H1_ p       0.000498   0.000060   0.000149   0.000041   0.000063   0.000045
    H2_ s      -0.000167  -0.001622  -0.000001  -0.000274   0.000006   0.000001
    H2_ p       0.000339   0.000453   0.000305   0.000370   0.000108   0.000057
 
   ao class      13b1g      14b1g      15b1g      16b1g      17b1g      18b1g
    C1_ p      -0.000008  -0.000047   0.000153  -0.000018   0.000024   0.000010
    C1_ d       0.000049   0.000065  -0.000148   0.000008   0.000152  -0.000021
    C2_ s       0.000138   0.000065  -0.000000  -0.000066  -0.000166   0.000034
    C2_ p       0.000177   0.000108   0.000337   0.000033   0.000056   0.000069
    C2_ d       0.000040   0.000113  -0.000329   0.000141   0.000104  -0.000013
    H1_ p       0.000099   0.000011   0.000064   0.000004  -0.000048   0.000009
    H2_ s      -0.000048   0.000129  -0.000000   0.000059   0.000026  -0.000043
    H2_ p       0.000005  -0.000160   0.000134  -0.000012  -0.000021   0.000050
 
   ao class      19b1g      20b1g      21b1g      22b1g      23b1g      24b1g
    C1_ p       0.000147   0.000001  -0.000024  -0.000021  -0.000003  -0.000009
    C1_ d      -0.000159   0.000035   0.000023   0.000063  -0.000001  -0.000008
    C2_ s       0.000001   0.000043   0.000001   0.000000   0.000005   0.000020
    C2_ p       0.000332  -0.000036   0.000006  -0.000152  -0.000002  -0.000018
    C2_ d      -0.000328   0.000029  -0.000006   0.000181   0.000000   0.000002
    H1_ p       0.000023  -0.000026  -0.000003  -0.000014   0.000000   0.000011
    H2_ s       0.000000   0.000012   0.000035   0.000002  -0.000003   0.000006
    H2_ p       0.000044  -0.000011  -0.000012  -0.000042   0.000016   0.000005
 
   ao class      25b1g      26b1g      27b1g      28b1g      29b1g      30b1g
    C1_ p      -0.000002   0.000000  -0.000000   0.000000  -0.000000   0.000000
    C1_ d      -0.000000  -0.000000  -0.000001  -0.000000   0.000000   0.000000
    C2_ s       0.000014  -0.000002  -0.000000  -0.000001   0.000000   0.000000
    C2_ p      -0.000012   0.000006   0.000002   0.000001   0.000000  -0.000000
    C2_ d      -0.000001  -0.000000  -0.000003  -0.000000   0.000000   0.000000
    H1_ p       0.000000   0.000000   0.000001  -0.000000  -0.000000  -0.000000
    H2_ s       0.000001  -0.000003  -0.000000  -0.000000   0.000000   0.000000
    H2_ p       0.000005   0.000000   0.000002  -0.000000   0.000000  -0.000000

                        b1u partial gross atomic populations
   ao class       1b1u       2b1u       3b1u       4b1u       5b1u       6b1u
    C1_ p       0.494826   0.004136   0.026433   0.000035   0.000089   0.002546
    C1_ d       0.010462   0.000007  -0.000961   0.001729   0.000778   0.000067
    C2_ p       1.398831   0.764093   0.032363   0.000066   0.000092   0.000086
    C2_ d       0.020611   0.023571   0.001955   0.003094   0.002330   0.000498
    H1_ p       0.002487   0.001963   0.000572   0.000113   0.000321   0.000147
    H2_ p       0.009666  -0.002042   0.000071   0.000227   0.000115   0.000012
 
   ao class       7b1u       8b1u       9b1u      10b1u      11b1u      12b1u
    C1_ p       0.000027   0.000322   0.000007   0.000006   0.000021   0.000028
    C1_ d       0.000244   0.000034   0.000013   0.000100   0.000047   0.000094
    C2_ p       0.001301   0.000283   0.000033  -0.000002   0.000010   0.000005
    C2_ d       0.000521   0.000275   0.000023   0.000101   0.000385   0.000079
    H1_ p       0.000090   0.000072   0.000216   0.000145  -0.000017  -0.000016
    H2_ p       0.000223   0.000035   0.000353   0.000097  -0.000071   0.000001
 
   ao class      13b1u      14b1u      15b1u      16b1u
    C1_ p       0.000015  -0.000001   0.000022   0.000002
    C1_ d       0.000048   0.000015   0.000002  -0.000002
    C2_ p      -0.000003  -0.000000  -0.000005  -0.000000
    C2_ d       0.000051  -0.000019   0.000010  -0.000001
    H1_ p      -0.000017   0.000019   0.000007   0.000007
    H2_ p       0.000012   0.000050  -0.000001   0.000003

                        b2g partial gross atomic populations
   ao class       1b2g       2b2g       3b2g       4b2g       5b2g       6b2g
    C1_ p       1.111492   0.006770   0.001225   0.001252   0.000142   0.000011
    C1_ d       0.001773   0.000028   0.001607   0.000684   0.000344   0.000027
    C2_ p       0.549269   0.012077   0.000661   0.000619   0.000370   0.000017
    C2_ d       0.041056  -0.000052   0.001046   0.000473   0.000514   0.001005
    H1_ p       0.015857   0.000109   0.000458   0.000069   0.000287   0.000008
    H2_ p       0.008201   0.000306   0.000186   0.000071   0.000559   0.000001
 
   ao class       7b2g       8b2g       9b2g      10b2g      11b2g      12b2g
    C1_ p       0.000008   0.000239   0.000051   0.000001  -0.000002  -0.000009
    C1_ d       0.000050  -0.000003   0.000142  -0.000152   0.000096   0.000025
    C2_ p       0.000007   0.000193   0.000024   0.000307   0.000051  -0.000026
    C2_ d       0.000050  -0.000090   0.000012   0.000027   0.000006   0.000109
    H1_ p       0.000317   0.000007  -0.000002   0.000048  -0.000017  -0.000007
    H2_ p       0.000151   0.000043   0.000069   0.000027  -0.000003  -0.000022
 
   ao class      13b2g      14b2g      15b2g      16b2g
    C1_ p       0.000002   0.000000   0.000000   0.000001
    C1_ d       0.000032  -0.000000  -0.000002  -0.000000
    C2_ p      -0.000028   0.000000   0.000000   0.000001
    C2_ d       0.000065  -0.000004  -0.000001  -0.000000
    H1_ p      -0.000011   0.000009   0.000006   0.000000
    H2_ p      -0.000008   0.000011   0.000005   0.000000

                        b3g partial gross atomic populations
   ao class       1b3g       2b3g       3b3g       4b3g       5b3g       6b3g
    C1_ d       0.018712   0.002037   0.000125   0.000048   0.000619   0.000017
    C2_ p       1.137272   0.000051   0.001059   0.001854   0.000041   0.000003
    C2_ d       0.009224   0.003548   0.002837   0.000654   0.000415   0.000090
    H2_ p       0.015617   0.000013   0.000748   0.000092   0.000007   0.000452
 
   ao class       7b3g       8b3g       9b3g      10b3g      11b3g
    C1_ d       0.000008   0.000087   0.000001   0.000040  -0.000000
    C2_ p       0.000258   0.000021   0.000033   0.000002   0.000002
    C2_ d       0.000209   0.000202   0.000088   0.000018  -0.000005
    H2_ p      -0.000007  -0.000004  -0.000016  -0.000000   0.000024

                        au  partial gross atomic populations
   ao class       1au        2au        3au        4au        5au        6au 
    C1_ d       0.012050   0.000809   0.001025   0.000246   0.000181   0.000043
    C2_ p       0.224014   0.001137   0.001599   0.000230   0.000000   0.000000
    C2_ d       0.005976   0.001769   0.000589   0.000707   0.000397   0.000155
    H2_ p       0.007992   0.000653   0.000022   0.000232   0.000000   0.000242
 
   ao class       7au        8au        9au       10au       11au 
    C1_ d      -0.000011   0.000070   0.000002   0.000011   0.000000
    C2_ p       0.000199   0.000003   0.000008   0.000000   0.000002
    C2_ d       0.000082   0.000046   0.000044   0.000021  -0.000005
    H2_ p      -0.000027   0.000016  -0.000008   0.000000   0.000015


                        gross atomic populations
     ao           C1_        C2_        H1_        H2_
      s         5.883691  12.100252   2.814344   5.591988
      p         5.030735  11.146547   0.011777   0.029889
      d        -0.214063  -0.395159   0.000000   0.000000
    total      10.700363  22.851639   2.826121   5.621877
 

 Total number of electrons:   42.00000000

