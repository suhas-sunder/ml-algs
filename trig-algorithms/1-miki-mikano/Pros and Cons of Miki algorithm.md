### Assumptions
- The waveform of the signal is a sinusoid of the nominal  frequency  
- The frequency of the signal is invariant.
### WHAT WE NEED TO MODEL EQNS IN MATLAB

![[20250514205238.png]]
![[20250514205332.png]]
![[20250514205357.png]]
- This last slide gives us all the equations we need to plot on MATLAB. 
- The imaginary and real parts can help us plot the filter itself.
- The Estimated phasor is what we will use to model the filter. These equations will be applied to input frequencies in order to get a filtered signal/output. 
### WHAT WE WILL BE GIVEN TO PLUG IN
![[20250514205830.png]]
- We will probably only be given freq = 60 Hz and sf = 720Hz
- Delta T = 1/sf
- Omega Delta T is calculated using nominal frequency and sampling frequency
	- omega delta T = 720/60 = 30 degrees
