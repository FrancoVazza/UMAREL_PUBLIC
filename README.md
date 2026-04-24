# UMAREL


 
<img src="umarel_logo.png" alt="alt text">  
 
 Parallel code in Jula language to simulate the propagation of UHECRs in cosmological simulations, developed by F.Vazza, A.Firinu (University of Bologna) and C. Evoli (GSSI). 


UMAREL (Ultra-high-energy cosmic rays in Magnetic fields Affected by Rigidity diffusion and Energy Losses) injects large sets
of cosmic rays into a simulated volume and self-consistently evolve their spatial trajectories and energies in time. See Firinu, Vazza & Evoli 2026 (submitted) for details.

## Key features

* particle propagator: Borish busher;
* loss terms: continuous loss terms from tabulated tables including interaction with the EBL and the CMB;
* sources selected from a galaxy catalog;
* cosmological effects;
* particles are evolved while the background simulation also is evolved, by combining differnt timesteps;
* multiple generation epochs of particles are allowed;
* the codes is parallelised using Julia and has been tested up to 128 processors. 

## Examples of results
This is a movie showing multiple injections of UHECRs, where each color represents a different generation of particles
 
<img src="movie.gif" alt="" width="600" height="600">

This is the simulated propagation of 100,000 UHECR protons injected at z=1 and evolved until z=0.

<img src="map_t.png" alt="" width="600" height="600">

This is the evolution of the energy of each simulated UHECR proton, as a function of time. Multiple spikes mark the epochs of generation of new families of UHECR protons.

<img src="map_t.png" alt="" width="600" height="600">



This is one of the statistics which can be produced with UMAREL: the average distance covered by UHECR protons since their injection, with (solid lines) or without (dashed) the effects of extragalactic magnetic fields.

 
<img src="fidvsnolosses.png" alt="" width="600" height="600">


## What is the public version of UMAREL
The public version of UMAREL shared here is meant to work on a laptop, using a sequence of 3D snapshots of cosmological simulations 

Umarel gives a fresh view on old problems! 
 
<img src="cover.png" alt="alt text" width="600" height="580">
