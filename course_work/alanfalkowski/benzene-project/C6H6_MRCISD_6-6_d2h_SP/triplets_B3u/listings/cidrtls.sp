 
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
 /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/MRCI_6-6/triplets/WO
 Spin-Orbit CI Calculation?(y,[n])
 Spin-Free Calculation
 
 input the spin multiplicity [  0]:
 spin multiplicity, smult            :   3    triplet 
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
 state spatial symmetry label        :  b3u
 
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

 number of rows in the drt : 161
    23 arcs removed due to generalized interacting space restrictions.

 manual arc removal step:


 input the level, a, b, and step values 
 for the arcs to be removed (-1/ to end).

 input the level, a, b, and step (-1/ to end):
 remarc:   0 arcs removed out of   0 specified.

 xbarz=       26757
 xbary=      223596
 xbarx=      613725
 xbarw=      714807
        --------
 nwalk=     1578885
 input the range of drt levels to print (l1,l2):
 levprt(*)        -1   0

 reference-csf selection step 1:
 total number of z-walks in the drt, nzwalk=   26757

 input the list of allowed reference symmetries:
 allowed reference symmetries:             2
 allowed reference symmetry labels:      b3u
 keep all of the z-walks as references?(y,[n])
 all z-walks are initially deleted.
 
 generate walks while applying reference drt restrictions?([y],n)
 reference drt restrictions will be imposed on the z-walks.
 
 impose additional orbital-group occupation restrictions?(y,[n])
 
 apply primary reference occupation restrictions?(y,[n])
 
 manually select individual walks?(y,[n])

 step 1 reference csf selection complete.
       48 csfs initially selected from   26757 total walks.

 beginning step-vector based selection.
 enter [internal_orbital_step_vector/disposition] pairs:

 enter internal orbital step vector, (-1/ to end):
   3  4  5  6 42 43 44 80 81 82110111139140141155171182

 step 2 reference csf selection complete.
       48 csfs currently selected from   26757 total walks.

 beginning numerical walk based selection.
 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end:

 input reference walk number (0 to end) [  0]:

 numerical walk-number based selection complete.
       48 reference csfs selected from   26757 total z-walks.
 
 input the reference occupations, mu(*):
 reference occupations:
 mu(*)=            2   2   2   2   2   2   2   2   2   2   2   2   0   0   0
                   0   0   0
 
 interacting space determination:
 checking diagonal loops...
 checking 2-internal loops...
 checking 3-internal loops...
 checking 4-internal loops...
 limint: nvalw=                  5628                 41898
                  6259                  6729
 post-limint icd(*)=                     1                     1
                 26760               1605647                     0
                     0
 
 this is an obsolete prompt.(y,[n])

 final mrsdci walk selection step:

 nvalw(*)=    5628   41898    6259    6729 nvalwt=   60514

 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end.

 input mrsdci walk number (0 to end) [  0]:

 end of manual mrsdci walk selection.
 number added=   0 number removed=   0

 nvalw(*)=    5628   41898    6259    6729 nvalwt=   60514

 lprune input numv1,nwalk=                 60514               1578885
 lprune input xbar(1,1),nref=                 26757                    48

 lprune: l(*,*,*) pruned with nwalk= 1578885 nvalwt=   60514=5628****62596729
 lprune:  z-drt, nprune=   200
 lprune:  y-drt, nprune=   161
 lprune: wx-drt, nprune=   191

 xbarz=       23025
 xbary=       49686
 xbarx=       15792
 xbarw=       17871
        --------
 nwalk=      106374
 levprt(*)        -1   0

 beginning the reference csf index recomputation...

     iref   iwalk  step-vector
   ------  ------  ------------
        1       1  333333333333331100
        2       4  333333333333330011
        3       5  333333333333313100
        4       8  333333333333312011
        5      15  333333333333311021
        6      16  333333333333311012
        7      21  333333333333310130
        8      24  333333333333310103
        9      29  333333333333303011
       10      33  333333333333301130
       11      36  333333333333301103
       12      39  333333333333300311
       13      42  333333333333133100
       14      45  333333333333132011
       15      52  333333333333131021
       16      53  333333333333131012
       17      58  333333333333130130
       18      61  333333333333130103
       19      66  333333333333123011
       20      70  333333333333121130
       21      73  333333333333121103
       22      76  333333333333120311
       23      85  333333333333113021
       24      86  333333333333113012
       25      91  333333333333112130
       26      94  333333333333112103
       27      99  333333333333111230
       28     102  333333333333111203
       29     107  333333333333110321
       30     108  333333333333110312
       31     118  333333333333103130
       32     121  333333333333103103
       33     124  333333333333102311
       34     128  333333333333101321
       35     129  333333333333101312
       36     138  333333333333100133
       37     141  333333333333033011
       38     145  333333333333031130
       39     148  333333333333031103
       40     151  333333333333030311
       41     157  333333333333013130
       42     160  333333333333013103
       43     163  333333333333012311
       44     167  333333333333011321
       45     168  333333333333011312
       46     177  333333333333010133
       47     178  333333333333003311
       48     183  333333333333001133
 indx01:    48 elements set in vec01(*)

 beginning the valid upper walk index recomputation...
 indx01: 60514 elements set in vec01(*)

 beginning the final csym(*) computation...

  number of valid internal walks of each symmetry:

       ag      b3u     b2u     b1g     b1u     b2g     b3g     au 
      ----    ----    ----    ----    ----    ----    ----    ----
 z               0            5628               0               0               0               0               0               0
 y            1772            1761            1761            1750            8697            8731            8695            8731
 x             919             668             872             920             711             729             711             729
 w             885            1228             844             892             711             729             711             729

 csfs grouped by internal walk symmetry:

       ag      b3u     b2u     b1g     b1u     b2g     b3g     au 
      ----    ----    ----    ----    ----    ----    ----    ----
 z               0            5628               0               0               0               0               0               0
 y           60248           58113           47547           45500          130455          113503           86950           87310
 x         1947361         1361384         1791960         1891520         1043037         1070901         1003221         1030077
 w         1875315         2708968         1734420         1833952         1043037         1070901         1003221         1030077

 total csf counts:
 z-vertex:            5628
 y-vertex:          629626
 x-vertex:        11139461
 w-vertex:        12299891
           --------
 total:        24074606
 
 input a title card, default=cidrt_title
 title card:
  cidrt_title                                                                   
  
 
 input a drt file name, default=cidrtfl
 drt and indexing arrays will be written to file:
 /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/MRCI_6-6/triplets/WO
 
 write the drt file?([y],n)
 drt file is being written...
 wrtstr:  ag  b3u b2u b1g b1u b2g b3g au 
nwalk=  106374 cpos=   35092 maxval=    9 cmprfactor=   67.01 %.
nwalk=  106374 cpos=   31937 maxval=   99 cmprfactor=   39.95 %.
 compressed with: nwalk=  106374 cpos=   35500 maxval=    9 cmprfactor=   66.63 %.
initial index vector length:    106374
compressed index vector length:     35500reduction:  66.63%
nwalk=   23025 cpos=    2935 maxval=    9 cmprfactor=   87.25 %.
nwalk=   23025 cpos=     313 maxval=   99 cmprfactor=   97.28 %.
nwalk=   23025 cpos=     102 maxval=  999 cmprfactor=   98.67 %.
nwalk=   23025 cpos=      82 maxval= 9999 cmprfactor=   98.58 %.
 compressed with: nwalk=   23025 cpos=     102 maxval=  999 cmprfactor=   98.67 %.
initial ref vector length:     23025
compressed ref vector length:       102reduction:  99.56%
