# Collapse-of-Classical-Electron

Project on C++ simulations of electron motion under Classical Collapse due to theoretical radiation.

Please read the REPORT I wrote in the PDF file.

----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------

# Objectives and Process

I made this project in my Computational Physics class at Fordham University. This was a 
research/simulation project that I made to "investigate" Bohr's postulates before Quantum Mechanics
was formulated. In those latters, the electron would lose energy due to its radiation and thus 
collapse.The project was meant to simulate the motion of such a "Classical" Electron and calculate 
the time for which it would exist.

  My goal was to implement the equations of motions related to this case, as well as compare two 
algorithms, the Velocity Verlet and the Runge-Kutta 4. I thus wrote four pieces of code: two for 
normal centripetal motion; two for the electron motion.

  I then wrote a five pages report on what I found and what I was able to accomplished in the time 
frame I was given. (This is the separate PDF one can find with my name on it)

-----------------------------------------------------------------------------------------------------

# A couple of things about this project

- This project was written in C++ on a Linux platform and uses 
  Gnuplot for all graph related matters.  

- The "orbit_rk.cpp" / "orbit_vl.cpp" files contains the code for 
  normal, classical, centripetal motion for the Runge-Kutta 4 (rk)
  and Velocity Verlet (vl) algorithms.
  
- The "orbit_electron_rk.cpp" / "orbit_electron_vl.cpp" files 
  contains the electron collapse code for both algorithms.

- In the "Results" folder, there are two sub-folders containing the 
  raw data as DAT files and orbital graphs as PS and PDF files for
  respectively the comparison of the two algorithms and the motion
  of the collapsing electron.

-----------------------------------------------------------------------------------------------------

# A couple of things about the code

- If you want to compile it, please make sure you refer to the usage
  given. You should not need any particular library downloaded in 
  order to make the program run.

- The program is not the most efficient if you want to use realistic
  scaling. What I mean by that is that if you want to use small 
  enough time step while initially placing the electron at Bohr's 
  radius, you most likely will not be able to obtain any real results.
  This is what I explain in my paper. 

- In "Results," there are my results as previously mentioned.
  However, do not be alarmed if you do not see such a file 
  when you run the code. "DatFile" files should appear; they
  are the results.
  
-----------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------

Fall 2021.

R. Van Laer.
