#using Distributed, SlurmClusterManager
#addprocs(SlurmManager())
#@everywhere println("hello from $(myid()):$(gethostname())")

#.......START OF UMAREL 

@everywhere using LaTeXStrings
@everywhere using Plots
@everywhere Plots.PyPlotBackend()
@everywhere using Statistics
@everywhere using SpecialFunctions
@everywhere using DelimitedFiles
@everywhere using HDF5
@everywhere using Unitful
@everywhere using Random
@everywhere using DelimitedFiles
@everywhere using Interpolations

#...just two parameters to give here 
@everywhere tag = "test"    #....will be attached to all output file names to differentiate them if necessary 
@everywhere np = 60000   #...total number of UHECR     
#..PATH TO FOLDERS OF ROUTINES AND FILES
#@everywhere      main="/leonardo_scratch/fast/IscrC_UMAREL/Julia/UMAREL_P/"  #..main folder containing UMAREL functions
@everywhere main = "/Users/francovazza/Library/CloudStorage/Dropbox/Julia_prog/UHECRpropa/UMAREL/UMAREL_P/PUBLIC/"
@everywhere include(string(main, "/constants.jl"))
@everywhere include(string(main, "parameters_UMAREL.jl"))
@everywhere include(string(main, "functions_UMAREL_new0.jl"))   #...external module with all relevant functions used for the transport of CRs    
#@everywhere      energy,dEdt,interp_losses=lossesC(Z,main)   #....loading the appropriate tabulated loss function 
@everywhere interp_losses, Emin, Emax = load_losses()
@everywhere times, zed, dt, max_it = define_times(courant, scale, zfin, zin, ch, cOmegaM)   #...computes the maximum iterations and the array of redshift to cover the entire evolution


@everywhere ngen = 1   #....ref: ngen=25
@everywhere inj_times = Array{Float64}(undef, ngen)
for i in 1:ngen
  inj_times[i] = convert(Int64, (i - 1) * trunc(max_it / ngen) + 1)
end
println(inj_times)


#....this array will store the trajectories and energy evolution of all particles 
@everywhere p = Array{Float64}(undef, 11, np)
@everywhere p .= 0.0

@everywhere npath = convert(Int64, 1 + trunc(max_it / skip_path))
@everywhere path = Array{Float64}(undef, np, 10, npath)
@everywhere path .= -100


#....MAIN SIMULATION CODE 
@everywhere nw = nworkers()
proc_0 = 2 #...
pid = collect(proc_0:proc_0+nw-1)

@everywhere np_proc = convert(Int64, trunc(np / nw))
@everywhere idp = Array{Int64}(undef, nw, np_proc)
idp .= 0
@inbounds for j in 1:np_proc
  @inbounds for i in 1:nw
    idp[i, j] = i + nw * (j - 1)
  end
end

@everywhere iz0 = 0   #...useful counters 
@everywhere inj = 0
@everywhere inj0 = 0
@everywhere i0 = 1

@inbounds for t in i0:max_it-1   #...main time loop 
  global iz0, inj, np, inj0, path, p, cosmo, use_syntheticB

  println("iteration ", t, "of ", max_it, "z=", zed[t])

  iz = findz(zeds, zed[t])  #...from the list of available snapshots, we find the one closest to the epoch of this t-iteration

  #...INJECTION OF NEW UHECR 
  iw = findall(inj_times .== t)
  nij = size(iw)
  if nij[1] == 1 #...injection  
    global inj += 1
    inj0 = inj
    println(root0, root_halos, snap[iz], i1, i2, j1, j2, l1, l2, mass_source, E_initial, np, " ", inj, " ", ngen, " ", n)
    global p = inject_new_UHECR(root0, root_halos, snap[iz], p, i1, i2, j1, j2, l1, l2, mass_source, E_initial, np, inj, ngen, n, sim, Lbox, zed[t])
    println("NEW INJECTION at z=", zed[t])
  end
  iz0 = iz[1]



  #...MAIN PROPAGATION & LOSSES ROUTINE
  t0 = time()
  println(procs())
  #....SENDING THE ADVECTION ROUTINE TO nw PROCESSORS, EACH OF THEM PROCESS THE p[:,cc1[i]:cc2[i]] PORTION OF THE ARRAY 
  @sync for i in 1:nw
    @async p[:, idp[i, :]] = remotecall_fetch(move_CR, pid[i], p[:, idp[i, :]], t, skip_path, scale, courant, dt, dx, i1, i2, j1, j2, l1, l2, cdd, cv, cb, Emin, Emax, interp_losses, Z, zed[t], zed[t+1], ngen, inj, root0, snap[iz0], sim)#,bx,by,bz) #...evolve CR in time 

  end

  println("dt=", time() - t0)

  #....we write in the path[] file (to be written on disk) only one step every skip_path, to save memory 
  it_path = convert(Int64, round(1 + t / skip_path))
  if trunc(1 + t / skip_path) == it_path
    t0 = time()
    path[:, 1, it_path] = p[1, :]   #x position with periodic BC 
    path[:, 2, it_path] = p[2, :]   #y position with periodic BC 
    path[:, 3, it_path] = p[3, :]   #z position with periodic BC 
    path[:, 4, it_path] = p[10, :]  #energy [eV] 
    path[:, 5, it_path] = p[7, :]   # x position without periodic BC - coordinates go beyond the box size 
    path[:, 6, it_path] = p[8, :]   # y position without periodic BC - coordinates go beyond the box size 
    path[:, 7, it_path] = p[9, :]   # z position without periodic BC - coordinates go beyond the box size 
    path[:, 8, it_path] = p[11, :]  #...physical magnetic field amplitude
    path[:, 9, it_path] .= zed[t]   #...redshift
    path[:, 10, it_path] .= t * dt  #...time since the start of the simulation, in seconds

    a = write_path(root_out, E_initial, Z, cosmo, tag, path, t, it_path, inj)
    a = plot_distance(path, it_path, tag, root_out)

  end
end
println("the propagation of UHECRs is done, now plotting")

#....plotting and writing of data on disk 
it_path = convert(Int64, trunc(max_it / skip_path))
t = max_it
a = write_path(root_out, E_initial, Z, cosmo, tag, path, t, it_path, max_it)
a = plot_map_WBC(dx, np, path, E_initial, Z, root_out, cosmo, tag)   #...plot trajectories without periodic BC 
a = make_gif(skip_path, npath, E_initial, path, ngen, Z, root_out, n, cosmo, tag)    #...gif with UHECR propagation
a = plot_spectrum(path, np, root_out, E_initial, Z, cosmo, tag)
a = plot_energy(path, np, root_out, tag)


println("end of run")
