# UMAREL


 
<img src="umarel_logo.png" alt="alt text">  
 
 Parallel code in Jula language to simulate the propagation of UHECRs in cosmological simulations, developed by F.Vazza, A.Firinu (University of Bologna) and C. Evoli (GSSI). 


UMAREL (Ultra-high-energy cosmic rays in Magnetic fields Affected by Rigidity diffusion and Energy Losses) injects large sets
of cosmic rays into a simulated volume and self-consistently evolve their spatial trajectories and energies in time. See Firinu, Vazza & Evoli 2026 (submitted) for details.

Key features

* particle propagator: Borish busher;
* loss terms: continuous loss terms from tabulated tables including interaction with the EBL and the CMB;
* sources selected from a galaxy catalog;
* cosmological effects;
* particles are evolved while the background simulation also is evolved, by combining differnt timesteps;
* multiple generation epochs of particles are allowed;
* the codes is parallelised using Julia and has been tested up to 128 processors. 

This is the simulated propagation of 1e5 UHECR protons from z=1 to z=0.


<img src="map_t.png" alt="" width="600" height="600">


This is a movie showing multiple injections of UHECRs, where each color represents a different generation of particles

 
<img src="movie.gif" alt="" width="600" height="600">


Umarel gives a fresh view on old problems! 
 
<img src="cover.png" alt="alt text" width="600" height="580">
