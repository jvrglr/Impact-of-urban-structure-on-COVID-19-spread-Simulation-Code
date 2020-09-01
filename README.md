# Impact-of-urban-structure-on-COVID-19-spread-Simulation-Code
Code associated to simulations of *Aguilar, Javier, et al. "Impact of urban structure on COVID-19 spread."
  arXiv preprint arXiv:2007.15367 (2020).*  
  
  Due to privacy agreements with google, it is not allowed to share directly the data used in the manuscript. A random fully-connected network is provided (filled with random populations), so it is possible to run tests and grasp the data format.  
  
  # Run code
  
  In order to run the test code, just execute the bash script *Super_script_P.sh*, this will:
  1. Compile FORTRAN 95 code.
  2. Create folders to store output data.
  3. run code. By default, the code is run one time. The user can change the number of times to execute the code with the variable *realizations* of the *Super_script_P.sh* script.
  
  # Aknowledgements
  
  Richard Chandler and Paul Northrop, authors of RANDGEN.F (source of Random2.f)  
  https://www.ucl.ac.uk/~ucakarc/work/software/randgen.txt
