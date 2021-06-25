
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
 nrow  =   161 nsym  =     8 ssym  =     2 lenbuf=  1600
 spnorb=     F spnodd=     F lxyzir(1:3)= 0 0 0
 nwalk,xbar:        106374       23025       49686       15792       17871
 nvalwt,nvalw:       60514        5628       41898        6259        6729
 ncsft:           24074606
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

 indx01: 60514 indices saved in indxv(*)
 test nroots froot                      1                    -1
===================================ROOT # 1===================================

 rdhciv: CI vector file information:
  cidrt_title                                                                    
 energy computed by program ciudg.       compute-0-12      15:45:02.531 24-Jun-21

 lenrec =   32768 lenci =  24074606 ninfo =  6 nenrgy =  6 ntitle =  2

 Max. overlap with ref vector #        1
 Valid ci vector #        1
 Method:        0       91% overlap
 energy( 1)=  2.015958396276E+02, ietype=   -1,    core energy of type: Nuc.Rep.
 energy( 2)= -2.935437811656E+02, ietype=    6,   fcore energy of type: H1(*)   
 energy( 3)= -1.325027193261E+02, ietype=    5,   fcore energy of type: Vref(*) 
 energy( 4)= -2.312777567153E+02, ietype=-1026,   total energy of type: MRSDCI  
 energy( 5)=  1.498852581201E-04, ietype=-2055, cnvginf energy of type: CI-Resid
 energy( 6)=  9.114042853753E-11, ietype=-2056, cnvginf energy of type: CI-D.E. 
==================================================================================
 test nroots froot                      2                    -1
===================================ROOT # 2===================================

 rdhciv: CI vector file information:
  cidrt_title                                                                    
 energy computed by program ciudg.       compute-0-12      15:45:02.558 24-Jun-21

 lenrec =   32768 lenci =  24074606 ninfo =  6 nenrgy =  7 ntitle =  2

 Max. overlap with ref vector #        2
 Valid ci vector #        2
 Method:        0       91% overlap
 energy( 1)=  2.015958396276E+02, ietype=   -1,    core energy of type: Nuc.Rep.
 energy( 2)= -2.935437811656E+02, ietype=    6,   fcore energy of type: H1(*)   
 energy( 3)= -1.325027193261E+02, ietype=    5,   fcore energy of type: Vref(*) 
 energy( 4)= -2.312525798997E+02, ietype=-1026,   total energy of type: MRSDCI  
 energy( 5)=  7.943632483301E-04, ietype=-2055, cnvginf energy of type: CI-Resid
 energy( 6)=  6.402463030852E-07, ietype=-2056, cnvginf energy of type: CI-D.E. 
 energy( 7)=  2.604125635624E-07, ietype=-2057, cnvginf energy of type: CI-ApxDE
==================================================================================
space sufficient for valid walk range           1      60514
               respectively csf range           1   24074606
space sufficient for valid walk range           1      60514
               respectively csf range           1   24074606

 space is available for 873728753 coefficients.

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


 rdcivnew:          39 coefficients were selected.
 workspace: ncsfmx=    24074606
 ncsfmx=              24074606

 histogram parameters:
 csfmn = 1.0000E-02 csfmx = 1.0000E+00 fhist = 5.0000E-01
 nhist =  20 icsfmn =           1 icsfmx =    24074606 ncsft =    24074606 ncsf =          39
 nhist =  20 fhist = 0.50000

    cmin                cmax        num  '*'=     1 csfs.
 ----------          ----------   ----- ---------|---------|---------|---------|
 5.0000E-01 <= |c| < 1.0000E+00       2 **
 2.5000E-01 <= |c| < 5.0000E-01       0
 1.2500E-01 <= |c| < 2.5000E-01       1 *
 6.2500E-02 <= |c| < 1.2500E-01       2 **
 3.1250E-02 <= |c| < 6.2500E-02       8 ********
 1.5625E-02 <= |c| < 3.1250E-02      15 ***************
 7.8125E-03 <= |c| < 1.5625E-02      11 ***********
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
                  total read =           39 total stored =          39

 from the selected csfs,
 min(|csfvec(:)|) = 1.0382E-02    max(|csfvec(:)|) = 6.3471E-01
 norm=  0.999999999996277     
 csfs will be printed based on coefficient magnitudes.

 current csfvec(*) selection parameters:
 csfmn = 1.0000E-02 csfmx = 1.0000E+00 fhist = 5.0000E-01
 nhist =  20 icsfmn =           1 icsfmx =    24074606 ncsft =    24074606 ncsf =          39

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
         12 -0.63471 0.40286 z*                    333333333333300311
          7  0.58802 0.34576 z*                    333333333333310130
         10  0.21059 0.04435 z*                    333333333333301130
          8 -0.08948 0.00801 z*                    333333333333310103
          2  0.08709 0.00758 z*                    333333333333330011
       7857 -0.05943 0.00353 y           b2g:  4  1333333333333100330
         30  0.05261 0.00277 z*                    333333333333110312
         36 -0.05068 0.00257 z*                    333333333333100133
         17 -0.04974 0.00247 z*                    333333333333130130
          4  0.04299 0.00185 z*                    333333333333312011
         11 -0.03720 0.00138 z*                    333333333333301103
       7145  0.03685 0.00136 y           b2g:  4  1333333333333120211
         22 -0.03497 0.00122 z*                    333333333333120311
         35  0.02814 0.00079 z*                    333333333333101312
       6162  0.02704 0.00073 y           b2g:  4  1333333333333300121
         46 -0.02616 0.00068 z*                    333333333333010133
         34 -0.02435 0.00059 z*                    333333333333101321
       7489  0.02412 0.00058 y           b2g:  4  1333333333333110122
         40  0.02406 0.00058 z*                    333333333333030311
       7859  0.02362 0.00056 y           b2g:  6  1333333333333100330
       5881 -0.02103 0.00044 y           b2g:  4  1333333333333310030
       5826  0.01999 0.00040 y           b2g:  4  1333333333333310300
       7898  0.01922 0.00037 y           b2g:  4  1333333333333100303
       7856  0.01799 0.00032 y           b2g:  3  1333333333333100330
       6892  0.01770 0.00031 y           b2g:  4  1333333333333130030
       6134  0.01722 0.00030 y           b2g:  4  1333333333333300211
       7147 -0.01665 0.00028 y           b2g:  6  1333333333333120211
         25 -0.01616 0.00026 z*                    333333333333112130
         29 -0.01509 0.00023 z*                    333333333333110321
       7655  0.01437 0.00021 y           b2g:  4  1333333333333102211
       7862  0.01365 0.00019 y           b2g:  9  1333333333333100330
          9  0.01300 0.00017 z*                    333333333333303011
         33 -0.01299 0.00017 z*                    333333333333102311
         43  0.01286 0.00017 z*                    333333333333012311
       6164 -0.01106 0.00012 y           b2g:  6  1333333333333300121
       6488 -0.01095 0.00012 y           b2g:  4  1333333333333210211
       7150 -0.01087 0.00012 y           b2g:  9  1333333333333120211
       7491 -0.01060 0.00011 y           b2g:  6  1333333333333110122
     591437 -0.01038 0.00011 y           b3u: 18  1133333333333300330
           39 csfs were printed in this range.

================================================================================
===================================VECTOR # 2===================================
================================================================================


 rdcivnew:          37 coefficients were selected.
 workspace: ncsfmx=    24074606
 ncsfmx=              24074606

 histogram parameters:
 csfmn = 1.0000E-02 csfmx = 1.0000E+00 fhist = 5.0000E-01
 nhist =  20 icsfmn =           1 icsfmx =    24074606 ncsft =    24074606 ncsf =          37
 nhist =  20 fhist = 0.50000

    cmin                cmax        num  '*'=     1 csfs.
 ----------          ----------   ----- ---------|---------|---------|---------|
 5.0000E-01 <= |c| < 1.0000E+00       2 **
 2.5000E-01 <= |c| < 5.0000E-01       0
 1.2500E-01 <= |c| < 2.5000E-01       2 **
 6.2500E-02 <= |c| < 1.2500E-01       2 **
 3.1250E-02 <= |c| < 6.2500E-02       4 ****
 1.5625E-02 <= |c| < 3.1250E-02       9 *********
 7.8125E-03 <= |c| < 1.5625E-02      18 ******************
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
                  total read =           37 total stored =          37

 from the selected csfs,
 min(|csfvec(:)|) = 1.0096E-02    max(|csfvec(:)|) = 6.1959E-01
 norm=  0.999999999996266     
 csfs will be printed based on coefficient magnitudes.

 current csfvec(*) selection parameters:
 csfmn = 1.0000E-02 csfmx = 1.0000E+00 fhist = 5.0000E-01
 nhist =  20 icsfmn =           1 icsfmx =    24074606 ncsft =    24074606 ncsf =          37

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
         12  0.61959 0.38389 z*                    333333333333300311
          7  0.60849 0.37026 z*                    333333333333310130
         10  0.16918 0.02862 z*                    333333333333301130
         29  0.15508 0.02405 z*                    333333333333110321
         17 -0.07489 0.00561 z*                    333333333333130130
         34  0.07049 0.00497 z*                    333333333333101321
       6177  0.05260 0.00277 y           b2g:  4  1333333333333300112
          2 -0.03933 0.00155 z*                    333333333333330011
          8 -0.03835 0.00147 z*                    333333333333310103
         46 -0.03586 0.00129 z*                    333333333333010133
         40 -0.03024 0.00091 z*                    333333333333030311
       5881 -0.02982 0.00089 y           b2g:  4  1333333333333310030
       6162 -0.02710 0.00073 y           b2g:  4  1333333333333300121
         25 -0.02584 0.00067 z*                    333333333333112130
       6179 -0.02045 0.00042 y           b2g:  6  1333333333333300112
         43 -0.01955 0.00038 z*                    333333333333012311
       7145 -0.01772 0.00031 y           b2g:  4  1333333333333120211
       6176 -0.01614 0.00026 y           b2g:  3  1333333333333300112
       7173 -0.01595 0.00025 y           b2g:  4  1333333333333120121
         35 -0.01533 0.00023 z*                    333333333333101312
       7489  0.01424 0.00020 y           b2g:  4  1333333333333110122
       5841  0.01364 0.00019 y           b3g:  4  1333333333333310210
         18  0.01248 0.00016 z*                    333333333333130103
          5  0.01229 0.00015 z*                    333333333333311021
       5883  0.01210 0.00015 y           b2g:  6  1333333333333310030
       6164  0.01202 0.00014 y           b2g:  6  1333333333333300121
   11795561 -0.01178 0.00014 w           b2g:  4  3333333333333300011
       6182 -0.01149 0.00013 y           b2g:  9  1333333333333300112
       6892  0.01143 0.00013 y           b2g:  4  1333333333333130030
   11778816 -0.01124 0.00013 w           b2g:  4  3333333333333310100
         48 -0.01121 0.00013 z*                    333333333333001133
         20 -0.01100 0.00012 z*                    333333333333121130
       9880  0.01091 0.00012 y           b2u: 13  1333333333332310211
       6205 -0.01063 0.00011 y           b3g:  4  1333333333333300031
      88487 -0.01032 0.00011 y           b1g: 11  1333333332333310211
         47 -0.01027 0.00011 z*                    333333333333003311
       6488  0.01010 0.00010 y           b2g:  4  1333333333333210211
           37 csfs were printed in this range.
