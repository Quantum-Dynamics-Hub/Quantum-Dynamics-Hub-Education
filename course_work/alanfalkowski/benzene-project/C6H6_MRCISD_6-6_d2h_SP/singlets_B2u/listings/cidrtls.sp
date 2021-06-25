 
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
 /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/MRCI_6-6/singlets/WO
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
 state spatial symmetry label        :  b2u
 
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
 allowed reference symmetries:             3
 allowed reference symmetry labels:      b2u
 keep all of the z-walks as references?(y,[n])
 all z-walks are initially deleted.
 
 generate walks while applying reference drt restrictions?([y],n)
 reference drt restrictions will be imposed on the z-walks.
 
 impose additional orbital-group occupation restrictions?(y,[n])
 
 apply primary reference occupation restrictions?(y,[n])
 
 manually select individual walks?(y,[n])

 step 1 reference csf selection complete.
       40 csfs initially selected from   17815 total walks.

 beginning step-vector based selection.
 enter [internal_orbital_step_vector/disposition] pairs:

 enter internal orbital step vector, (-1/ to end):
   3  4  5  6 42 43 44 80 81 82110111139140141155171182

 step 2 reference csf selection complete.
       40 csfs currently selected from   17815 total walks.

 beginning numerical walk based selection.
 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end:

 input reference walk number (0 to end) [  0]:

 numerical walk-number based selection complete.
       40 reference csfs selected from   17815 total z-walks.
 
 input the reference occupations, mu(*):
 reference occupations:
 mu(*)=            2   2   2   2   2   2   2   2   2   2   2   2   0   0   0
                   0   0   0
 
 interacting space determination:
 checking diagonal loops...
 checking 2-internal loops...
 checking 3-internal loops...
 checking 4-internal loops...
 limint: nvalw=                  3760                 29922
                  4857                  5335
 post-limint icd(*)=                     1                     1
                 17818               1072951                     0
                     0
 
 this is an obsolete prompt.(y,[n])

 final mrsdci walk selection step:

 nvalw(*)=    3760   29922    4857    5335 nvalwt=   43874

 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end.

 input mrsdci walk number (0 to end) [  0]:

 end of manual mrsdci walk selection.
 number added=   0 number removed=   0

 nvalw(*)=    3760   29922    4857    5335 nvalwt=   43874

 lprune input numv1,nwalk=                 43874               1055131
 lprune input xbar(1,1),nref=                 17815                    40

 lprune: l(*,*,*) pruned with nwalk= 1055131 nvalwt=   43874=3760****48575335
 lprune:  z-drt, nprune=   183
 lprune:  y-drt, nprune=   147
 lprune: wx-drt, nprune=   174

 xbarz=       15207
 xbary=       34818
 xbarx=       13467
 xbarw=       15446
        --------
 nwalk=       78938
 levprt(*)        -1   0

 beginning the reference csf index recomputation...

     iref   iwalk  step-vector
   ------  ------  ------------
        1       1  333333333333331020
        2       4  333333333333330102
        3       8  333333333333313020
        4      11  333333333333312102
        5      16  333333333333311202
        6      18  333333333333310320
        7      25  333333333333310023
        8      27  333333333333303102
        9      31  333333333333301320
       10      38  333333333333301023
       11      42  333333333333300132
       12      45  333333333333133020
       13      48  333333333333132102
       14      53  333333333333131202
       15      55  333333333333130320
       16      62  333333333333130023
       17      64  333333333333123102
       18      68  333333333333121320
       19      75  333333333333121023
       20      79  333333333333120132
       21      83  333333333333113202
       22      85  333333333333112320
       23      92  333333333333112023
       24      95  333333333333110232
       25      97  333333333333103320
       26     104  333333333333103023
       27     108  333333333333102132
       28     112  333333333333101232
       29     115  333333333333100323
       30     118  333333333333033102
       31     122  333333333333031320
       32     129  333333333333031023
       33     133  333333333333030132
       34     136  333333333333013320
       35     143  333333333333013023
       36     147  333333333333012132
       37     151  333333333333011232
       38     154  333333333333010323
       39     159  333333333333003132
       40     163  333333333333001323
 indx01:    40 elements set in vec01(*)

 beginning the valid upper walk index recomputation...
 indx01: 43874 elements set in vec01(*)

 beginning the final csym(*) computation...

  number of valid internal walks of each symmetry:

       ag      b3u     b2u     b1g     b1u     b2g     b3g     au 
      ----    ----    ----    ----    ----    ----    ----    ----
 z               0               0            3760               0               0               0               0               0
 y            1104            1092            1092            1080            6397            6367            6407            6383
 x             741             708             548             748             530             518             538             526
 w             751             704            1024             744             530             518             538             526

 csfs grouped by internal walk symmetry:

       ag      b3u     b2u     b1g     b1u     b2g     b3g     au 
      ----    ----    ----    ----    ----    ----    ----    ----
 z               0               0            3760               0               0               0               0               0
 y           28704           29484           36036           36720           63970           63670           83291           95745
 x         1523496         1454940         1116824         1585012          748890          730898          790322          771642
 w         1544056         1446720         2258944         1576536          748890          730898          790322          771642

 total csf counts:
 z-vertex:            3760
 y-vertex:          437620
 x-vertex:         8722024
 w-vertex:         9868008
           --------
 total:        19031412
 
 input a title card, default=cidrt_title
 title card:
  cidrt_title                                                                   
  
 
 input a drt file name, default=cidrtfl
 drt and indexing arrays will be written to file:
 /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/MRCI_6-6/singlets/WO
 
 write the drt file?([y],n)
 drt file is being written...
 wrtstr:  ag  b3u b2u b1g b1u b2g b3g au 
nwalk=   78938 cpos=   31178 maxval=    9 cmprfactor=   60.50 %.
nwalk=   78938 cpos=   29553 maxval=   99 cmprfactor=   25.12 %.
 compressed with: nwalk=   78938 cpos=   31479 maxval=    9 cmprfactor=   60.12 %.
initial index vector length:     78938
compressed index vector length:     31479reduction:  60.12%
nwalk=   15207 cpos=    1960 maxval=    9 cmprfactor=   87.11 %.
nwalk=   15207 cpos=     233 maxval=   99 cmprfactor=   96.94 %.
nwalk=   15207 cpos=      95 maxval=  999 cmprfactor=   98.13 %.
nwalk=   15207 cpos=      81 maxval= 9999 cmprfactor=   97.87 %.
 compressed with: nwalk=   15207 cpos=      95 maxval=  999 cmprfactor=   98.13 %.
initial ref vector length:     15207
compressed ref vector length:        95reduction:  99.38%
