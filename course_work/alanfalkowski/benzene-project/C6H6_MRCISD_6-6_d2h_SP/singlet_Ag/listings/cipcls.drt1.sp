
 program cipc      

 print the csf info for mrsdci wave functions

 written by: ron shepard

 version date: 06-jun-96

 This Version of Program cipc is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              CIPC        **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************


 workspace allocation parameters: lencor=2621440000 mem1=         0 ifirst=         1

 drt header information:
  cidrt_title                                                                    
 nmot  =   192 niot  =    18 nfct  =     6 nfvt  =     0
 nrow  =   148 nsym  =     8 ssym  =     1 lenbuf=  1600
 spnorb=     F spnodd=     F lxyzir(1:3)= 0 0 0
 nwalk,xbar:         79180       15280       34818       13545       15537
 nvalwt,nvalw:       48873        3909       32226        6039        6699
 ncsft:           24145633
 map(*)=    -1 -1169170171172  1  2  3  4  5  6  7  8  9 10 11 12 13 14
            15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 -1
            -1173174175 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
            50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 -1176
           177178 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85
            86 87 88 89 90 91 92 93 -1179180 94 95 96 97 98 99100101102
           103104105106107108109110111112113114115116117118119120181182
           183121122123124125126127128129130131132133184134135136137138
           139140141142143144145146147148185149150151152153154155156157
           158186159160161162163164165166167168
 mu(*)=      2  2  2  2  2  2  2  2  2  2  2  2  0  0  0  0  0  0
 syml(*) =   1  1  1  1  2  2  2  3  3  3  4  4  5  5  5  6  7  8
 rmo(*)=     3  4  5  6  3  4  5  2  3  4  2  3  1  2  3  1  1  1

 indx01: 48873 indices saved in indxv(*)
 test nroots froot                      1                     1
===================================ROOT # 1===================================

 rdhciv: CI vector file information:
  cidrt_title                                                                    
 energy computed by program ciudg.       compute-0-12      15:05:03.069 24-Jun-21

 lenrec =   32768 lenci =  24145633 ninfo =  6 nenrgy =  7 ntitle =  2

 Max. overlap with ref vector #        1
 Valid ci vector #        1
 Method:        0       91% overlap
 energy( 1)=  2.015958396276E+02, ietype=   -1,    core energy of type: Nuc.Rep.
 energy( 2)= -2.935437811656E+02, ietype=    6,   fcore energy of type: H1(*)   
 energy( 3)= -1.325027193261E+02, ietype=    5,   fcore energy of type: Vref(*) 
 energy( 4)= -2.314365577719E+02, ietype=-1026,   total energy of type: MRSDCI  
 energy( 5)=  7.050612000940E-04, ietype=-2055, cnvginf energy of type: CI-Resid
 energy( 6)=  9.676383614377E-07, ietype=-2056, cnvginf energy of type: CI-D.E. 
 energy( 7)=  1.533324097188E-07, ietype=-2057, cnvginf energy of type: CI-ApxDE
==================================================================================
space sufficient for valid walk range           1      48873
               respectively csf range           1   24145633

 space is available for 873742976 coefficients.

 updated histogram parameters:
 csfmn = 0.0000E+00 csfmx = 1.0000E+00 fhist = 5.0000E-01 nhist =  20

 this program will print the csfs generated from
 the drt according to the following print options :

 1) run in batch mode: all valid roots are automatically
    analysed and csf info is printed by default contribution
    threshold 0.01 
 2) run in interactive mode
 3) generate files for cioverlap without symmetry
 4) generate files for cioverlap with symmetry

 input menu number [  1]:
================================================================================
===================================VECTOR # 1===================================
================================================================================


 rdcivnew:          49 coefficients were selected.
 workspace: ncsfmx=    24145633
 ncsfmx=              24145633

 histogram parameters:
 csfmn = 1.0000E-02 csfmx = 1.0000E+00 fhist = 5.0000E-01
 nhist =  20 icsfmn =           1 icsfmx =    24145633 ncsft =    24145633 ncsf =          49
 nhist =  20 fhist = 0.50000

    cmin                cmax        num  '*'=     1 csfs.
 ----------          ----------   ----- ---------|---------|---------|---------|
 5.0000E-01 <= |c| < 1.0000E+00       1 *
 2.5000E-01 <= |c| < 5.0000E-01       0
 1.2500E-01 <= |c| < 2.5000E-01       0
 6.2500E-02 <= |c| < 1.2500E-01       4 ****
 3.1250E-02 <= |c| < 6.2500E-02       6 ******
 1.5625E-02 <= |c| < 3.1250E-02       6 ******
 7.8125E-03 <= |c| < 1.5625E-02      32 ********************************
 3.9062E-03 <= |c| < 7.8125E-03       0
 1.9531E-03 <= |c| < 3.9062E-03       0
 9.7656E-04 <= |c| < 1.9531E-03       0
 4.8828E-04 <= |c| < 9.7656E-04       0
 2.4414E-04 <= |c| < 4.8828E-04       0
 1.2207E-04 <= |c| < 2.4414E-04       0
 6.1035E-05 <= |c| < 1.2207E-04       0
 3.0518E-05 <= |c| < 6.1035E-05       0
 1.5259E-05 <= |c| < 3.0518E-05       0
 7.6294E-06 <= |c| < 1.5259E-05       0
 3.8147E-06 <= |c| < 7.6294E-06       0
 0.0000E+00 <= |c| < 3.8147E-06       0
                                  ----- ---------|---------|---------|---------|
                  total read =           49 total stored =          49

 from the selected csfs,
 min(|csfvec(:)|) = 1.0132E-02    max(|csfvec(:)|) = 8.8553E-01
 norm=  0.999999999994306     
 csfs will be printed based on coefficient magnitudes.

 current csfvec(*) selection parameters:
 csfmn = 1.0000E-02 csfmx = 1.0000E+00 fhist = 5.0000E-01
 nhist =  20 icsfmn =           1 icsfmx =    24145633 ncsft =    24145633 ncsf =          49

 i:slabel(i) =  1: ag   2: b3u  3: b2u  4: b1g  5: b1u  6: b2g  7: b3g  8: au 
 
 frozen orbital =    1    2    3    4    5    6
 symfc(*)       =    1    1    2    2    3    4
 label          =  ag   ag   b3u  b3u  b2u  b1g
 rmo(*)         =    1    2    1    2    1    1
 
 internal level =    1    2    3    4    5    6    7    8    9   10
   11   12   13   14   15   16   17   18
 syml(*)        =    1    1    1    1    2    2    2    3    3    3
    4    4    5    5    5    6    7    8
 label          =  ag   ag   ag   ag   b3u  b3u  b3u  b2u  b2u  b2u
  b1g  b1g  b1u  b1u  b1u  b2g  b3g  au 
 rmo(*)         =    3    4    5    6    3    4    5    2    3    4
    2    3    1    2    3    1    1    1

 printing selected csfs in sorted order from cmin = 0.00000 to cmax = 1.00000

   indcsf     c     c**2   v  lab:rmo  lab:rmo   step(*)
  ------- -------- ------- - ---- --- ---- --- ------------
         15 -0.88553 0.78416 z*                    333333333333300330
         16  0.10401 0.01082 z*                    333333333333300303
          3  0.09814 0.00963 z*                    333333333333330030
          9  0.09494 0.00901 z*                    333333333333310122
          8 -0.06307 0.00398 z*                    333333333333310212
         14  0.05425 0.00294 z*                    333333333333301122
         17  0.04235 0.00179 z*                    333333333333300033
       5500 -0.03954 0.00156 y           b2g:  4  1333333333333100322
          6  0.03876 0.00150 z*                    333333333333312030
       5232  0.03830 0.00147 y           b2g:  4  1333333333333120230
          2  0.03164 0.00100 z*                    333333333333330300
         55  0.02997 0.00090 z*                    333333333333000333
          5  0.02862 0.00082 z*                    333333333333312300
         44  0.02663 0.00071 z*                    333333333333030330
          4 -0.01885 0.00036 z*                    333333333333330003
       5502  0.01799 0.00032 y           b2g:  6  1333333333333100322
       5234 -0.01742 0.00030 y           b2g:  6  1333333333333120230
         10  0.01553 0.00024 z*                    333333333333303300
       5411  0.01505 0.00023 y           b2g:  4  1333333333333102230
     242531 -0.01460 0.00021 y           b3u: 16  1333233333333300312
       7165 -0.01441 0.00021 y           b2u: 13  1333333333332310230
      63646  0.01429 0.00020 y           b3u: 16  1333333332333100332
      62512 -0.01377 0.00019 y           b2u: 13  1333333332333120330
         13 -0.01351 0.00018 z*                    333333333333301212
      61353 -0.01331 0.00018 y           b1g: 11  1333333332333300312
      60742  0.01303 0.00017 y           b1g: 11  1333333332333310230
      38805 -0.01283 0.00016 y           b1g: 11  1333333333233300312
      61410 -0.01277 0.00016 y           ag : 18  1333333332333300132
     119740  0.01269 0.00016 y           ag : 18  1333333233333310230
         49  0.01261 0.00016 z*                    333333333333012330
       5505  0.01196 0.00014 y           b2g:  9  1333333333333100322
     102973 -0.01176 0.00014 y           b1g: 11  1333333313333300322
      61357 -0.01146 0.00013 y           b1g: 15  1333333332333300312
     217746  0.01142 0.00013 y           ag : 18  1333313333333320230
       5237 -0.01133 0.00013 y           b2g:  9  1333333333333120230
     244857  0.01129 0.00013 y           b1g: 15  1333233333333100332
   11728386  0.01116 0.00012 w           b2g:  4  3333333333333000330
      60679 -0.01114 0.00012 y           ag : 18  1333333332333310320
   11537994  0.01108 0.00012 w           b2g:  4  3333333333333300300
   11544311  0.01102 0.00012 w           b2g:  4  3333333333333300030
       8943 -0.01090 0.00012 y           b1g: 15  1333333333332120330
       5515 -0.01077 0.00012 y           b3g:  4  1333333333333100232
       5499  0.01062 0.00011 y           b2g:  3  1333333333333100322
      60746 -0.01053 0.00011 y           b1g: 15  1333333332333310230
       5212 -0.01036 0.00011 y           b3g:  4  1333333333333120320
       7834 -0.01033 0.00011 y           b3u: 16  1333333333332300132
       5231 -0.01025 0.00011 y           b2g:  3  1333333333333120230
       4937 -0.01019 0.00010 y           b2g:  4  1333333333333200312
         11  0.01013 0.00010 z*                    333333333333303030
           49 csfs were printed in this range.
