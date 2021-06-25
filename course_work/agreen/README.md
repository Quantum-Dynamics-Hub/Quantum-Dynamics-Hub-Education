Author: Austin Green
e: austig2@uci.edu
6/24/2021

The DVR wavepacket calculations attached use Libra's gaussian wavepacket module to 
calculate density matrix elements (the populations and coherences as a function of
time) with intended applications for spectroscopy. The python scripts and outputs
follow the directory hierarchy below
………………………………………………………………………………………………………………………………………………………………….
DVR_WAVEPACKET_CALCULATIONS
	Morse
		de_0.3
			morse.ipynb
			Holstein-morse
				Dyn.txt
de_0.5
			morse.ipynb
			Holstein-morse
				Dyn.txt
de_1.0
			morse.ipynb
			Holstein-morse
				Dyn.txt
	Harmonic
		de_0.3
			harmonic.ipynb
			Holstein-2
				Dyn.txt
de_0.5
			harmonic.ipynb
			Holstein-2
				Dyn.txt
de_1.0
			harmonic.ipynb
			Holstein-2
				Dyn.txt
………………………………………………………………………………………………………………………………………………………………….
where de_X.X refers to the dissociation energy of the excited state of a morse potential
or the potential derived from the harmonic approximation to a morse potential of that
disocciation energy. The morse potential function was created from the Holstein-2 model
of coupled oscillators in Libra’s model library by changing the Holstein-2 function
(input parameters, function, and derivatives). Dyn.txt outputs are tab delimited outputs
of the elements of the form: time, rho_{11}, rho_{22}, and Re[rho_{12}]. The morse
potential should be added and merge with the rest of the Hamiltons in the model library.
