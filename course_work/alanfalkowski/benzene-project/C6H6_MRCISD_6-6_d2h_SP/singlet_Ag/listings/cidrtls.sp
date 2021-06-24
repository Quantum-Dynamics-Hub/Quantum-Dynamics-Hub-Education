 
 program cidrt 7.0  

 distinct row table construction, reference csf selection, and internal
 walk selection for multireference single- and double-excitation
configuration interaction.

 references:  r. shepard, i. shavitt, r. m. pitzer, d. c. comeau, m. pepper
                  h. lischka, p. g. szalay, r. ahlrichs, f. b. brown, and
                  j.-g. zhao, int. j. quantum chem. symp. 22, 149 (1988).
              h. lischka, r. shepard, f. b. brown, and i. shavitt,
                  int. j. quantum chem. symp. 15, 91 (1981).

 based on the initial version by  Ron Shepard

 extended for spin-orbit CI calculations ( Russ Pitzer, OSU)

 and large active spaces (Thomas Müller, FZ(21 Juelich)

 This Version of Program CIDRT is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de

*********************** File revision status: ***********************
* cidrt1.F9 Revision: 1.1.6.2           Date: 2013/04/11 14:37:29   * 
* cidrt2.F9 Revision: 1.1.6.6           Date: 2015/02/26 17:04:32   * 
* cidrt3.F9 Revision: 1.1.6.2           Date: 2013/04/11 14:37:29   * 
* cidrt4.F9 Revision: 1.1.6.2           Date: 2013/04/11 14:37:29   * 
********************************************************************

 workspace allocation parameters: lencor=2621440000 mem1=         0 ifirst=         1
 expanded "keystrokes" are being written to file:
 /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/MRCI_6-6/singlets_Ag
 Spin-Orbit CI Calculation?(y,[n])
 Spin-Free Calculation
 
 input the spin multiplicity [  0]:
 spin multiplicity, smult            :   1    singlet 
 input the total number of electrons [  0]:
 total number of electrons, nelt     :    42
 input the number of irreps (1:8) [  0]:
 point group dimension, nsym         :     8
 enter symmetry labels:(y,[n])
 enter 8 labels (a4):
 enter symmetry label, default=   1
 enter symmetry label, default=   2
 enter symmetry label, default=   3
 enter symmetry label, default=   4
 enter symmetry label, default=   5
 enter symmetry label, default=   6
 enter symmetry label, default=   7
 enter symmetry label, default=   8
 symmetry labels: (symmetry, slabel)
 ( 1,  ag ) ( 2,  b3u) ( 3,  b2u) ( 4,  b1g) ( 5,  b1u) ( 6,  b2g) ( 7,  b3g) ( 8,  au ) 
 input nmpsy(*):
 nmpsy(*)=        39  39  30  30  16  16  11  11
 
   symmetry block summary
 block(*)=         1   2   3   4   5   6   7   8
 slabel(*)=      ag  b3u b2u b1g b1u b2g b3g au 
 nmpsy(*)=        39  39  30  30  16  16  11  11
 
 total molecular orbitals            :   192
 input the molecular spatial symmetry (irrep 1:nsym) [  0]:
 state spatial symmetry label        :  ag 
 
 input the frozen core orbitals (sym(i),rmo(i),i=1,nfct):
 total frozen core orbitals, nfct    :     6
 
 fcorb(*)=         1   2  40  41  79 109
 slabel(*)=      ag  ag  b3u b3u b2u b1g
 
 number of frozen core orbitals      :     6
 number of frozen core electrons     :    12
 number of internal electrons        :    30
 
 input the frozen virtual orbitals (sym(i),rmo(i),i=1,nfvt):
 total frozen virtual orbitals, nfvt :     0

 no frozen virtual orbitals entered
 
 input the internal orbitals (sym(i),rmo(i),i=1,niot):
 niot                                :    18
 
 modrt(*)=         3   4   5   6  42  43  44  80  81  82 110 111 139 140 141
                 155 171 182
 slabel(*)=      ag  ag  ag  ag  b3u b3u b3u b2u b2u b2u b1g b1g b1u b1u b1u
                 b2g b3g au 
 
 total number of orbitals            :   192
 number of frozen core orbitals      :     6
 number of frozen virtual orbitals   :     0
 number of internal orbitals         :    18
 number of external orbitals         :   168
 
 orbital-to-level mapping vector
 map(*)=          -1  -1 169 170 171 172   1   2   3   4   5   6   7   8   9
                  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24
                  25  26  27  28  29  30  31  32  33  -1  -1 173 174 175  34
                  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49
                  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64
                  65  66  67  -1 176 177 178  68  69  70  71  72  73  74  75
                  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90
                  91  92  93  -1 179 180  94  95  96  97  98  99 100 101 102
                 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117
                 118 119 120 181 182 183 121 122 123 124 125 126 127 128 129
                 130 131 132 133 184 134 135 136 137 138 139 140 141 142 143
                 144 145 146 147 148 185 149 150 151 152 153 154 155 156 157
                 158 186 159 160 161 162 163 164 165 166 167 168
 
 input the number of ref-csf doubly-occupied orbitals [  0]:
 (ref) doubly-occupied orbitals      :    12
 
 no. of internal orbitals            :    18
 no. of doubly-occ. (ref) orbitals   :    12
 no. active (ref) orbitals           :     6
 no. of active electrons             :     6
 
 input the active-orbital, active-electron occmnr(*):
 139140141155171182
 input the active-orbital, active-electron occmxr(*):
 139140141155171182
 
 actmo(*) =      139 140 141 155 171 182
 occmnr(*)=        0   0   0   0   0   6
 occmxr(*)=        6   6   6   6   6   6
 reference csf cumulative electron occupations:
 modrt(*)=         3   4   5   6  42  43  44  80  81  82 110 111 139 140 141
                 155 171 182
 occmnr(*)=        2   4   6   8  10  12  14  16  18  20  22  24  24  24  24
                  24  24  30
 occmxr(*)=        2   4   6   8  10  12  14  16  18  20  22  24  30  30  30
                  30  30  30
 
 input the active-orbital bminr(*):
 139140141155171182
 input the active-orbital bmaxr(*):
 139140141155171182
 reference csf b-value constraints:
 modrt(*)=         3   4   5   6  42  43  44  80  81  82 110 111 139 140 141
                 155 171 182
 bminr(*)=         0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
                   0   0   0
 bmaxr(*)=         0   0   0   0   0   0   0   0   0   0   0   0   6   6   6
                   6   6   6
 input the active orbital smaskr(*):
 139140141155171182
 modrt:smaskr=
   3:1000   4:1000   5:1000   6:1000  42:1000  43:1000  44:1000  80:1000
  81:1000  82:1000 110:1000 111:1000 139:1111 140:1111 141:1111 155:1111
 171:1111 182:1111
 
 input the maximum excitation level from the reference csfs [  2]:
 maximum excitation from ref. csfs:  :     2
 number of internal electrons:       :    30
 
 input the internal-orbital mrsdci occmin(*):
   3  4  5  6 42 43 44 80 81 82110111139140141155171182
 input the internal-orbital mrsdci occmax(*):
   3  4  5  6 42 43 44 80 81 82110111139140141155171182
 mrsdci csf cumulative electron occupations:
 modrt(*)=         3   4   5   6  42  43  44  80  81  82 110 111 139 140 141
                 155 171 182
 occmin(*)=        0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
                   0   0  28
 occmax(*)=       30  30  30  30  30  30  30  30  30  30  30  30  30  30  30
                  30  30  30
 
 input the internal-orbital mrsdci bmin(*):
   3  4  5  6 42 43 44 80 81 82110111139140141155171182
 input the internal-orbital mrsdci bmax(*):
   3  4  5  6 42 43 44 80 81 82110111139140141155171182
 mrsdci b-value constraints:
 modrt(*)=         3   4   5   6  42  43  44  80  81  82 110 111 139 140 141
                 155 171 182
 bmin(*)=          0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
                   0   0   0
 bmax(*)=         30  30  30  30  30  30  30  30  30  30  30  30  30  30  30
                  30  30  30
 
 input the internal-orbital smask(*):
   3  4  5  6 42 43 44 80 81 82110111139140141155171182
 modrt:smask=
   3:1111   4:1111   5:1111   6:1111  42:1111  43:1111  44:1111  80:1111
  81:1111  82:1111 110:1111 111:1111 139:1111 140:1111 141:1111 155:1111
 171:1111 182:1111
 
 internal orbital summary:
 block(*)=         1   1   1   1   2   2   2   3   3   3   4   4   5   5   5
                   6   7   8
 slabel(*)=      ag  ag  ag  ag  b3u b3u b3u b2u b2u b2u b1g b1g b1u b1u b1u
                 b2g b3g au 
 rmo(*)=           3   4   5   6   3   4   5   2   3   4   2   3   1   2   3
                   1   1   1
 modrt(*)=         3   4   5   6  42  43  44  80  81  82 110 111 139 140 141
                 155 171 182
 
 reference csf info:
 occmnr(*)=        2   4   6   8  10  12  14  16  18  20  22  24  24  24  24
                  24  24  30
 occmxr(*)=        2   4   6   8  10  12  14  16  18  20  22  24  30  30  30
                  30  30  30
 
 bminr(*)=         0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
                   0   0   0
 bmaxr(*)=         0   0   0   0   0   0   0   0   0   0   0   0   6   6   6
                   6   6   6
 
 
 mrsdci csf info:
 occmin(*)=        0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
                   0   0  28
 occmax(*)=       30  30  30  30  30  30  30  30  30  30  30  30  30  30  30
                  30  30  30
 
 bmin(*)=          0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
                   0   0   0
 bmax(*)=         30  30  30  30  30  30  30  30  30  30  30  30  30  30  30
                  30  30  30
 

 a priori removal of distinct rows:

 input the level, a, and b values for the vertices 
 to be removed (-1/ to end).

 input level, a, and b (-1/ to end):
 no vertices marked for removal
 
 impose generalized interacting space restrictions?(y,[n])
 generalized interacting space restrictions will be imposed.
 multp(*)=
  hmult                     0
 lxyzir   0   0   0
 symmetry of spin functions (spnir)
       --------------------------Ms ----------------------------
   S     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
   1     1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   2     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   3     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   4     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   5     1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   6     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   7     0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   8     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   9     1  0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0
  10     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  11     0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0
  12     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  13     1  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  0  0  0
  14     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  15     0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  0  0  0  0
  16     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  17     1  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  1  0  0
  18     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  19     0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0

 number of rows in the drt : 148
    23 arcs removed due to generalized interacting space restrictions.

 manual arc removal step:


 input the level, a, b, and step values 
 for the arcs to be removed (-1/ to end).

 input the level, a, b, and step (-1/ to end):
 remarc:   0 arcs removed out of   0 specified.

 xbarz=       17815
 xbary=      148008
 xbarx=      410109
 xbarw=      479199
        --------
 nwalk=     1055131
 input the range of drt levels to print (l1,l2):
 levprt(*)        -1   0

 reference-csf selection step 1:
 total number of z-walks in the drt, nzwalk=   17815

 input the list of allowed reference symmetries:
 allowed reference symmetries:             1
 allowed reference symmetry labels:      ag 
 keep all of the z-walks as references?(y,[n])
 all z-walks are initially deleted.
 
 generate walks while applying reference drt restrictions?([y],n)
 reference drt restrictions will be imposed on the z-walks.
 
 impose additional orbital-group occupation restrictions?(y,[n])
 
 apply primary reference occupation restrictions?(y,[n])
 
 manually select individual walks?(y,[n])

 step 1 reference csf selection complete.
       55 csfs initially selected from   17815 total walks.

 beginning step-vector based selection.
 enter [internal_orbital_step_vector/disposition] pairs:

 enter internal orbital step vector, (-1/ to end):
   3  4  5  6 42 43 44 80 81 82110111139140141155171182

 step 2 reference csf selection complete.
       55 csfs currently selected from   17815 total walks.

 beginning numerical walk based selection.
 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end:

 input reference walk number (0 to end) [  0]:

 numerical walk-number based selection complete.
       55 reference csfs selected from   17815 total z-walks.
 
 input the reference occupations, mu(*):
 reference occupations:
 mu(*)=            2   2   2   2   2   2   2   2   2   2   2   2   0   0   0
                   0   0   0
 
 interacting space determination:
 checking diagonal loops...
 checking 2-internal loops...
 checking 3-internal loops...
 checking 4-internal loops...
 limint: nvalw=                  3909                 32226
                  6039                  6699
 post-limint icd(*)=                     1                     1
                 17818               1072951                     0
                     0
 
 this is an obsolete prompt.(y,[n])

 final mrsdci walk selection step:

 nvalw(*)=    3909   32226    6039    6699 nvalwt=   48873

 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end.

 input mrsdci walk number (0 to end) [  0]:

 end of manual mrsdci walk selection.
 number added=   0 number removed=   0

 nvalw(*)=    3909   32226    6039    6699 nvalwt=   48873

 lprune input numv1,nwalk=                 48873               1055131
 lprune input xbar(1,1),nref=                 17815                    55

 lprune: l(*,*,*) pruned with nwalk= 1055131 nvalwt=   48873=3909****60396699
 lprune:  z-drt, nprune=   183
 lprune:  y-drt, nprune=   147
 lprune: wx-drt, nprune=   171

 xbarz=       15280
 xbary=       34818
 xbarx=       13545
 xbarw=       15537
        --------
 nwalk=       79180
 levprt(*)        -1   0

 beginning the reference csf index recomputation...

     iref   iwalk  step-vector
   ------  ------  ------------
        1       1  333333333333333000
        2       2  333333333333330300
        3       5  333333333333330030
        4       7  333333333333330003
        5       8  333333333333312300
        6      11  333333333333312030
        7      13  333333333333312003
        8      17  333333333333310212
        9      19  333333333333310122
       10      22  333333333333303300
       11      25  333333333333303030
       12      27  333333333333303003
       13      31  333333333333301212
       14      33  333333333333301122
       15      36  333333333333300330
       16      38  333333333333300303
       17      41  333333333333300033
       18      42  333333333333132300
       19      45  333333333333132030
       20      47  333333333333132003
       21      51  333333333333130212
       22      53  333333333333130122
       23      56  333333333333123300
       24      59  333333333333123030
       25      61  333333333333123003
       26      65  333333333333121212
       27      67  333333333333121122
       28      70  333333333333120330
       29      72  333333333333120303
       30      75  333333333333120033
       31      82  333333333333112212
       32      84  333333333333112122
       33      87  333333333333111222
       34      94  333333333333103212
       35      96  333333333333103122
       36      99  333333333333102330
       37     101  333333333333102303
       38     104  333333333333102033
       39     111  333333333333033300
       40     114  333333333333033030
       41     116  333333333333033003
       42     120  333333333333031212
       43     122  333333333333031122
       44     125  333333333333030330
       45     127  333333333333030303
       46     130  333333333333030033
       47     134  333333333333013212
       48     136  333333333333013122
       49     139  333333333333012330
       50     141  333333333333012303
       51     144  333333333333012033
       52     151  333333333333003330
       53     153  333333333333003303
       54     156  333333333333003033
       55     160  333333333333000333
 indx01:    55 elements set in vec01(*)

 beginning the valid upper walk index recomputation...
 indx01: 48873 elements set in vec01(*)

 beginning the final csym(*) computation...

  number of valid internal walks of each symmetry:

       ag      b3u     b2u     b1g     b1u     b2g     b3g     au 
      ----    ----    ----    ----    ----    ----    ----    ----
 z            3909               0               0               0               0               0               0               0
 y            1104            1092            1092            1080            6993            6963            6963            6939
 x             736            1018            1018             963             588             576             576             564
 w            1408            1014            1014             959             588             576             576             564

 csfs grouped by internal walk symmetry:

       ag      b3u     b2u     b1g     b1u     b2g     b3g     au 
      ----    ----    ----    ----    ----    ----    ----    ----
 z            3909               0               0               0               0               0               0               0
 y           36432           37128           28392           29160           90909          104445           69630           69390
 x         1499968         2157142         2093008         1978965          863772          844992          813888          795804
 w         3106048         2148666         2084784         1970745          863772          844992          813888          795804

 total csf counts:
 z-vertex:            3909
 y-vertex:          465486
 x-vertex:        11047539
 w-vertex:        12628699
           --------
 total:        24145633
 
 input a title card, default=cidrt_title
 title card:
  cidrt_title                                                                   
  
 
 input a drt file name, default=cidrtfl
 drt and indexing arrays will be written to file:
 /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/MRCI_6-6/singlets_Ag
 
 write the drt file?([y],n)
 drt file is being written...
 wrtstr:  ag  b3u b2u b1g b1u b2g b3g au 
nwalk=   79180 cpos=   29412 maxval=    9 cmprfactor=   62.85 %.
nwalk=   79180 cpos=   25605 maxval=   99 cmprfactor=   35.32 %.
 compressed with: nwalk=   79180 cpos=   29763 maxval=    9 cmprfactor=   62.41 %.
initial index vector length:     79180
compressed index vector length:     29763reduction:  62.41%
nwalk=   15280 cpos=    1993 maxval=    9 cmprfactor=   86.96 %.
nwalk=   15280 cpos=     258 maxval=   99 cmprfactor=   96.62 %.
nwalk=   15280 cpos=     119 maxval=  999 cmprfactor=   97.66 %.
nwalk=   15280 cpos=     105 maxval= 9999 cmprfactor=   97.25 %.
 compressed with: nwalk=   15280 cpos=     119 maxval=  999 cmprfactor=   97.66 %.
initial ref vector length:     15280
compressed ref vector length:       119reduction:  99.22%
