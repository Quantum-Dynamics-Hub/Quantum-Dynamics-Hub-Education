1
 program ciudg      
 multireference single and double excitation configuration
 interaction based on the graphical unitary group approach.


 ************************************************************************
 beginning the bk-type iterative procedure (nzcsf=  3909)...
 ************************************************************************

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.7649519907 -3.5527E-14  8.7920E-01  1.5702E+00  1.0000E-03   

 mr-sdci  convergence not reached after  1 iterations.

 final mr-sdci  convergence information:
 mr-sdci #  1  1   -230.7649519907 -3.5527E-14  8.7920E-01  1.5702E+00  1.0000E-03   

 from bk iterations: iconv=   1

 ************************************************************************
 beginning the ci iterative diagonalization procedure... 
 ************************************************************************

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -230.7649519907  1.7764E-15  8.7920E-01  1.5702E+00  1.0000E-03   
 mr-sdci #  2  1   -231.4087655143  6.4381E-01  3.1532E-02  3.1307E-01  1.0000E-03   
 mr-sdci #  3  1   -231.4346672832  2.5902E-02  2.3707E-03  8.3095E-02  1.0000E-03   
 mr-sdci #  4  1   -231.4363578597  1.6906E-03  2.0691E-04  2.5256E-02  1.0000E-03   
 mr-sdci #  5  1   -231.4365270740  1.6921E-04  2.4254E-05  8.1247E-03  1.0000E-03   
 mr-sdci #  6  1   -231.4365486447  2.1571E-05  4.5954E-06  3.7983E-03  1.0000E-03   
 mr-sdci #  7  1   -231.4365539467  5.3020E-06  2.0854E-06  2.4087E-03  1.0000E-03   
 mr-sdci #  8  1   -231.4365568042  2.8575E-06  8.2291E-07  1.5954E-03  1.0000E-03   
 mr-sdci #  9  1   -231.4365577719  9.6764E-07  1.5333E-07  7.0506E-04  1.0000E-03   

 mr-sdci  convergence criteria satisfied after  9 iterations.

 final mr-sdci  convergence information:
 mr-sdci #  9  1   -231.4365577719  9.6764E-07  1.5333E-07  7.0506E-04  1.0000E-03   

 number of reference csfs (nref) is    55.  root number (iroot) is  1.
 c0**2 =   0.82881229  c**2 (all zwalks) =   0.82956317

 eref      =   -230.764561570741   "relaxed" cnot**2         =   0.828812292110
 eci       =   -231.436557771864   deltae = eci - eref       =  -0.671996201123
 eci+dv1   =   -231.551595261245   dv1 = (1-cnot**2)*deltae  =  -0.115037489381
 eci+dv2   =   -231.575355772903   dv2 = dv1 / cnot**2       =  -0.138798001038
 eci+dv3   =   -231.611486587090   dv3 = dv1 / (2*cnot**2-1) =  -0.174928815226
 eci+pople =   -231.593731470990   ( 30e- scaled deltae )    =  -0.829169900248
