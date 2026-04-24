
@everywhere function cross_prod(a, b)
    c = similar(a)
    c[1] = a[2] * b[3] - a[3] * b[2]
    c[2] = a[3] * b[1] - a[1] * b[3]
    c[3] = a[1] * b[2] - a[2] * b[1]
    return c
end

@everywhere function findz(a, x)
    na = size(a)
    i0 = 1#na[1]
    @inbounds for i in eachindex(a)
        if i == 1
            continue
        end
        if x >= a[i] && x < a[i-1]
            i0 = i
        end
    end

    return i0
end

@everywhere function assign_CR_halo_z(p, i1, i2, j1, j2, l1, l2, filecat, mass_source, E_initial, np, inj, ngen, n, sim, Lbox, zed)
    #...this function assigns the initial position of UHECRs using the position of matter halos given by a catalog

    a = readdlm(filecat)

    imass = findall(a[:, 5] .>= mass_source)
    mass = a[imass, 5]
    radius = a[imass, 4]
    x = a[imass, 1]
    y = a[imass, 2]
    z = a[imass, 3]

    if sim == "MAKITRA"

        x .*= n   #..notice this is instead needed for MAKITRA runs
        y .*= n
        z .*= n
    end

    x .+= (-i1 + 1)
    y .+= (-j1 + 1)
    z .+= (-l1 + 1)

    im = sortperm(mass, rev=true)

    npt = convert(Int64, trunc(np / ngen))   #....number of new UHECR to be injected at this step 
    nim = size(im)
    println("sources with >Mass= ", nim, " ", filecat)

    if nim[1] <= 1
        imass = findall(a[:, 5] .>= mass_source * 0.1)
        mass = a[imass, 5]
        radius = a[imass, 4]
        x = a[imass, 1]
        y = a[imass, 2]
        z = a[imass, 3]
        if sim == "MAKITRA"
            x .*= n
            y .*= n
            z .*= n
        end
       
        x .+= (-i1 + 1)
        y .+= (-j1 + 1)
        z .+= (-l1 + 1)
        im = sortperm(mass, rev=true)
        nim = size(im)
        println("sources with >Mass= ", nim, " ", filecat)

    end

    tpi = 2 * pi

    @inbounds for i in 1:npt  #...loop over the number of new UHECRs 
        nE = length(E_initial)
        ip = npt * (inj - 1) + i        #...counter to begin the generation starting from the last UHECR evolved this far       
        if i <= nim[1]
            x1 = x[im[i]]
            y1 = y[im[i]]
            z1 = z[im[i]]
        end
        if i > nim[1]
            ni = convert(Int64, trunc(i / nim[1]))
            ii = i - ni * nim[1] + 1
            x1 = x[im[ii]]
            y1 = y[im[ii]]
            z1 = z[im[ii]]
        end

        ix = x1 + rand()
        iy = y1 + rand()
        iz = z1 + rand()
        θ = pi * rand()
        φ = tpi * rand()
        vx = vc * sin(θ) * cos(φ)     #...random initial velocity vector (v=c)
        vy = vc * sin(θ) * sin(φ)
        vz = vc * cos(θ)

        p[1, ip] = ix
        p[2, ip] = iy
        p[3, ip] = iz
        p[4, ip] = vx
        p[5, ip] = vy
        p[6, ip] = vz

        p[7, ip] = ix
        p[8, ip] = iy
        p[9, ip] = iz

        if nE == 1 && E_initial[1] > 0
            p[10, ip] = 10^E_initial[1]  #...initial energy in eV 
        end
        if nE > 1 && E_initial[1] > 0
            ii = convert(Int64, round(nE * rand()))
            if ii < 1
                ii = 1
            end
            if ii > nE
                ii = nE
            end
            p[10, ip] = 10^(E_initial[ii])  #....this generates are random distribution of initial energies picked from the given energy bins.
        end

        if nE > 1 && E_initial[1] < 0
            Er = (E_initial[3]) * rand()
            Ein_random = E_initial[2] + Er[1]
            p[10, ip] = 10^(Ein_random)  #....this generates a random distribution of initial energies picked from the given energy bins.

        end
    end

    return p
end


#...Boris pusher   
@everywhere function integ_kdk(p1, dvx, pb, dt, qm, scale, ng, courant, iqm, interp_losses, z1, z2)

    vp2 = p1[4:6]
    modv = (p1[4]^2 + p1[5]^2 + p1[6]^2)  #...velocity module before kick (must be=c)
    vp2 .+= dvx .* dt * 0.5  #..kick   - velocity is updated for half a step 

    #...enforce energy conservation
    diss = (vp2[1]^2 + vp2[2]^2 + vp2[3]^2)
    eps2 = (modv - diss) / (modv) #..energy dissipation by the numerical scheme
    vp2 .*= sqrt(modv / diss)   #..ensures the velocity module after kick is =c

    #...positions with enforced periodic BC (see below)
    p1[1] += vp2[1] * dt / (scale)  #...drift
    p1[2] += vp2[2] * dt / (scale)  #...drift
    p1[3] += vp2[3] * dt / (scale)  #...drift

    #...positions without enforced periodic BC
    p1[7] += vp2[1] * dt / (scale) #...drift
    p1[8] += vp2[2] * dt / (scale) #...drift
    p1[9] += vp2[3] * dt / (scale)  #...drift


    #...periodic boundary conditions
    if p1[1] < 1
        p1[1] = ng + p1[1] - 1
    end
    if p1[2] < 1
        p1[2] = ng + p1[2] - 1
    end
    if p1[3] < 1
        p1[3] = ng + p1[3] - 1
    end
    if p1[1] > ng
        p1[1] = p1[1] - ng + 1
    end
    if p1[2] > ng
        p1[2] = p1[2] - ng + 1
    end
    if p1[3] > ng
        p1[3] = p1[3] - ng + 1
    end

    i1 = convert(Int64, trunc(p1[1]))
    i2 = convert(Int64, trunc(p1[2]))
    i3 = convert(Int64, trunc(p1[3]))

    dvx2 = qm * cross_prod(vp2, pb) / vc
    if isnan(dvx2[1]) || isnan(dvx2[2]) || isnan(dvx2[3])   #...in case of some odd Nan numbers
        dvx2 .= 0.0
    end
    p1[4:6] .= vp2 + dvx2 * dt * 0.5   #...kick

    #...enforce energy conservation
    diss = (p1[4]^2 + p1[5]^2 + p1[6]^2)
    eps2 = abs((modv - diss) / (modv)) #...energy dissipation by the numerical scheme
    p1[4:6] .*= sqrt(modv / (diss))

    p1[10] = integrate_energy(p1[10], z1, z2, interp_losses)

    return p1[1:11]
end


# -------------------------
# RK4 integrator
# -------------------------

function integrate_energyB(E0, z1, z2, interp; Nsteps=30)   #...z2<z1
    z = z2
    E = E0
    dz = (z1 - z2) / Nsteps
    @inbounds @simd for i in 1:Nsteps
        k1 = dEdz(E, z, interp)
        k2 = dEdz(E + 0.5 * dz * k1, z + 0.5 * dz, interp)
        k3 = dEdz(E + 0.5 * dz * k2, z + 0.5 * dz, interp)
        k4 = dEdz(E + dz * k3, z + dz, interp)
        E -= dz * (k1 + 2k2 + 2k3 + k4) / 6
        z += dz
    end
    return E
end


function integrate_energy(E0, z1, z2, interp; Nsteps=65)   #...symmetric in z, centred at (z2+z1)/0.5
    z = (z2 + z1) / 2.0
    Δz = (z1 - z2) / 2.0
    E = E0
    dz = (z1 - z2) / Nsteps
    @inbounds for i in 1:Nsteps
        k1 = dEdz(E, z + Δz, interp)
        k2 = dEdz(E + 0.5 * dz * k1, z + 0.5 * dz + Δz, interp)
        k3 = dEdz(E + 0.5 * dz * k2, z + 0.5 * dz + Δz, interp)
        k4 = dEdz(E + dz * k3, z + dz + Δz, interp)
        E -= dz * (k1 + 2k2 + 2k3 + k4) / 6
        z += dz
    end
    return E
end

function dEdz(E, z, interp)
    Hz = H_z(z)
    b = beta_losses(E, z, interp)
    return (b / ((1 + z) * Hz)) * E
end


function beta_losses(E, z, interp)

    E_now = min(1e21, E * (1 + z))
    x = log10(E_now)

    return H_z(z) + (1 + z)^3 * 10.0^(interp(x))
end

function H_z(z; H0_km_s_Mpc=ch * 100, Ωm=cOmegaM, ΩΛ=1 - cOmegaM, Ωr=0.0, Ωk=nothing)
    if Ωk === nothing
        Ωk = 1.0 - (Ωm + ΩΛ + Ωr)
    end
    Mpc_in_km = 3.0856775814913673e19
    H0 = H0_km_s_Mpc / Mpc_in_km
    Ez2 = Ωr * (1 + z)^4 + Ωm * (1 + z)^3 + Ωk * (1 + z)^2 + ΩΛ
    return H0 * sqrt(Ez2)
end


@everywhere function move_CR(p::Array{Float64,2}, t::Int64, skip_path::Int64, scale::Float64, courant::Float64, dt::Float64, dx::Float64, i1::Int64, i2::Int64, j1::Int64, j2::Int64, l1::Int64, l2::Int64, cd::Float64, cv::Float64, cb::Float64, Emin, Emax, interp_losses, Z::Int64, zed1::Float64, zed2::Float64, ngen::Int64, inj::Int64, root0::String, snapiz::Int64, sim::String)

    #...selectio of a specific MAKITRA 512^3 RUN. 
    if sim == "MAKITRA"
        model = "minus1.0_08"
        file1 = string(root0, model, "RD00", snapiz, ".cpu0000")   #...ENZO snapshots with 3D data
        filecat = string(root_halos, snapiz, "_", model, "_halof_200_new.dat")   #...halo catalogs
    end

    #....we apply or not cosmological correction factors related to z 
    if cosmo == 0   #...the redshift is fixed to 0 everywhere 
        zed = 0.0
    end
    if cosmo == 1    #...the comoving scale is made proper at the given z  
        zed = zed2
        scale = scale / (1 + zed)
    end

    #...selection of atomic mass number based on the nuclear charge number
    #...this may be set in a more general way for all particles at once, but it does not cost time and it allows us to assign different particles a different Z if we want 
    if Z == 1  #proton 
        A = 1
    end
    if Z > 1
        println("currently not supported Z")
        error()
    end
    np1 = size(p)
    ng = i2 - i1 + 1

    npt = convert(Int64, trunc(np1[2] / (ngen)))
    nev = convert(Int64, npt * inj)


    @inbounds for i in 1:nev  #...only evolve the UHECRs injected so far - the other remain idle 
        zz = zed  #...we assume the redshift as constant in this interval 

        i1 = convert(Int64, trunc(p[1, i]))
        i2 = convert(Int64, trunc(p[2, i]))
        i3 = convert(Int64, trunc(p[3, i]))

        ng0 = ng

        #...periodic boundary conditions are enforced here 

        if i1 < 1
            i1 = ng0 + i1
        end
        if i2 < 1
            i2 = ng0 + i2
        end
        if i3 < 1
            i3 = ng0 + i3
        end
        if i1 > ng0
            i1 = i1 - ng0
        end
        if i2 > ng0
            i2 = i2 - ng0
        end
        if i3 > ng0
            i3 = i3 - ng0
        end
        vp = p[4:6, i]

        #....reading of B-fields only where particles are - necessary to spare memory for big data, but slower

        pbx = h5read(file1, "Grid00000001/Bx", (i1, i2, i3))
        pby = h5read(file1, "Grid00000001/By", (i1, i2, i3))
        pbz = h5read(file1, "Grid00000001/Bz", (i1, i2, i3))

        norm = 5.095941297750065e-6
        pb = [norm * pbx, norm * pby, norm * pbz] * ((1 + zed)^1.5) #/0.37   #.....magnetic field in the particle reference frame .   (1+z)^1.5 is to make it properly scale with the physical (1+z)^2

        γ = p[10, i] * evtoerg / (A * prest)     #....Lorentz factor of particles.

        qm = Z * qe / (A * mp * γ)  #...mass, charge and gamma of the particle to be set here 
        iqm = 1 / (qm)

        dvx = qm * cross_prod(vp, pb)#/vc   #....acceleration from Lorentz force
        @fastmath p[11, i] = sqrt(pb[1]^2 + pb[2]^2 + pb[3]^2)
        rl = iqm / p[11, i]

        p1 = p[1:11, i]
        pnew = integ_kdk(p1, dvx, pb, dt, qm, scale, ng, courant, iqm, interp_losses, zed1, zed2)  #...integrator of particle motion with kick-drift-kick method in the Borish pusher
        p[1:11, i] .= pnew


    end

    return p
end


function load_losses(filename="/Users/francovazza/Dropbox/Julia_prog/UHECRpropa/UMAREL/UMAREL_P/SimProp_losses_proton_new.txt")

    data = readdlm(filename)

    E_losses = data[:, 1]
    beta_pair = data[:, 3]
    beta_pion = data[:, 5]

    tiny = floatmin(Float64)

    beta_pair .= max.(beta_pair, tiny)
    beta_pion .= max.(beta_pion, tiny)

    beta0 = beta_pair .+ beta_pion

    logE = log10.(E_losses)
    logb = log10.(beta0) .- log10(3.1e16)

    interp = interpolate((logE,), logb, Gridded(Linear()))

    return interp, logE[1], logE[end]
end

@everywhere function lossesC(Z, main)  #...new function by C. Evoli, only for protons

    #....tabulated loss function 

    file_losses = string(main, "/SimProp_losses_proton_new.txt")
    a = readdlm(file_losses)
    energy = a[:, 1]
    nE = length(energy)
    loss_adiabatic = a[:, 2]  #adiabatic
    loss_pair_CMB = a[:, 3]  #pair production  CMB
    loss_pair_EBL = a[:, 4]  #pair production - EBL
    loss_pp_CMB = a[:, 5]  # photo pion production - CMB 
    loss_pp_EBL = a[:, 6]  #photo pion production - EBL 
    dEdt = Array{Float64}(undef, nE, 2)
    logE = similar(energy)
    timesE = similar(energy)
    for e in eachindex(energy)
        dEdt[e, 1] = (loss_adiabatic[e]) / (1e9 * yr)
        dEdt[e, 2] = (floatmin(Float64) + loss_pp_CMB[e] + loss_pp_EBL[e] + loss_pair_CMB[e] + loss_pair_EBL[e]) / (1e9 * yr)
        timesE[e] = 1 / (dEdt[e, 1] + dEdt[e, 2]) / (1e9 * yr)
        logE[e] = log10(energy[e] + floatmin(Float64))
    end

    #...cooling losses interpolated on a finer Energy grid in log(E) space 
    #         logE = log10.(energy)
    println(size(logE), " ", typeof(logE))


    logb = log10.(dEdt[:, 2])
    println(size(logb), " ", typeof(logb))
    interp = interpolate((logE,), logb, Gridded(Linear()))
    println(typeof(interp))
    @inbounds for i in 1:1000
        interp[i] = 10^interp[i]
    end

    return energy, dEdt, interp
end


function define_times(courant, scale, zfin, zin, ch, cOmegaM)

    #....DEFINITION OF THE FINAL ARRAY WITH ALL PARTICLE INFORMATION
    dt = courant * scale / vc
    max_it = convert(Int64, 1 + trunc(time_tot / dt))      #...maximum number of iterations
    println("going to evolve UHECR for ", max_it, " iterations, to cover ", time_tot / Gyr, "Gyr of evolution")

    times = Array{Float64}(undef, max_it)
    zed = Array{Float64}(undef, max_it)

    #...DEFINITION OF TIMESTEPS
    #...we use the approximated inversion between time and scale factor as in https://academic.oup.com/mnras/article/505/2/2764/6289943, eq.15-17 
    dz = (zin - zfin) / max_it
    time_tilde = (6.519 / ch) * (1 / sqrt(1 - cOmegaM))
    @inbounds for i in 1:max_it #...we define an appropriate array of times and redshifts 
        times[i] = dt * i / Gyr
        a = (sqrt(cOmegaM / (1 - cOmegaM)) * sinh((Δt / Gyr + times[i]) / time_tilde))^0.6666
        zed[i] = 1 / a - 1
    end

    return times, zed, dt, max_it

end


@everywhere function inject_new_UHECR(root0, root_halos, snapiz, p, i1, i2, j1, j2, l1, l2, mass_source, E_initial, np, inj, ngen, n, sim, LBox, z)

    #...MAKITRA 512^3 RUNs

    if sim == "MAKITRA"
        model = "minus1.0_08"
        file1 = string(root0, model, "RD00", snapiz, ".cpu0000")   #...ENZO snapshots with 3D data
        filecat = string(root_halos, snapiz, "_minus1.0_08_halof_200_new.dat")   #...halo catalogs 
    end
    p = assign_CR_halo_z(p, i1, i2, j1, j2, l1, l2, filecat, mass_source, E_initial, np, inj, ngen, n, sim, Lbox, z)   #...assign UHECR based on density  
    inj0 = inj
    ntr = size(p)    #....number of UHECR - it depends on the dsource density choosen above 
    np = ntr[2]
    npt = convert(Int64, trunc(np / ngen))   #....number of new UHECR to be injected at this step 
    nj1 = (inj - 1) * npt + 1
    nj2 = nj1 + npt - 1

    #....we plot the initial location of UHECR
    plo = plot(p[1, nj1:nj2], p[2, nj1:nj2], seriestype=:scatter, ms=1, label="", grid=false, aspect_ratio=1.0)  #...plot initial location of CRs
    xlims!(i1, i2)
    ylims!(j1, j2)
    filep1 = string(root_out, "_initial_map_newU", snapiz, "_cosmo", cosmo, ".png")
    savefig(filep1)

    return p
end




@everywhere function plot_map_WBC(dx, np, path, E_initial, Z, root_out, cosmo, tag)
    #1) static plot without boundary conditions
    xs = dx * 1e-3  #...to convert positions from cells to Mpc
    jump = 10
    #the following plot can use a lot of memory, do for i in 1:jump:np with jump = any integer number to plot only 1 every jump particles
    @inbounds for i in 1:jump:np-1   #....all particle trajectories are plotted without periodic BC 
        if i == 1
            plo = plot(path[i, 5, :] * xs, path[i, 6, :] * xs, label="", aspect_ratio=:equal, dpi=1000, lw=0.5, grid=false)
        end
        io = findall((path[i, 5, :] .> 0) .| (path[i, 5, :] .< 0))
        #          io=findall(abs.(path[i,5,:]). > 0) 
        plot!(path[i, 5, io] * xs, path[i, 6, io] * xs, label="", lw=0.5, alpha=0.5)#,seriestype=:scatter,ms=0.4)
    end
    tit = string(" E0=", E_initial)
    if E_initial[1] == -0
        tit = L" N(E) \propto E^{-1}"
    end
    title!(string("Z=", Z, " ,", tit, " eV"), fonts=20)
    yaxis!("[Mpc]", fonts=20)
    xaxis!("[Mpc]", fonts=20)
    filep1 = string(root_out, "_UHECR_path_color_map_", E_initial[1], "_Z_", Z, "_noBC_cosmo", cosmo, "_", tag, ".png")
    savefig(filep1)

end



@everywhere function make_gif(skip_path, npath, E_initial, path, nz, Z, root_out, n, cosmo, tag)
    #2) animated gif with periodic boundary conditions
    #...plotting evolution of UHECR in the periodic BC as an animated gif 
    it = 0
    xs = dx * 1e-3  #...to convert positions from cells to Mpc
    #...this plotting function in general produces way too many frames and it wastes some time
    #   nframe=npath  if one wants to produce the full movie with all steps in the path file
    #   nframe<npath  if one wants to limit the gif to a smaller number of frames, like npath=100 for the first 100 steps only  
    nframe = npath
    it = 0
    anim = @animate for t in 1:skip_path:nframe*skip_path   #....plotting positions 
        it += 1
        println(t)
        tit = string(" E0=", E_initial)
        if E_initial[1] == -1
            tit = L" N(E) \propto E^{-1}"
        end
        npt = convert(Int64, trunc(np / ngen))
        plot(path[1:2, 1, it] * xs, path[1:2, 2, it] * xs, label="", aspect_ratio=:equal, line=0, dpi=300, lw=0.0, marker=:circle, ms=0.1, grid=false)

        @inbounds for gg in 1:nz[1]
            g1 = (gg - 1) * npt + 1
            g2 = g1 + npt - 1
            plot!(path[g1:g2, 1, it] * xs, path[g1:g2, 2, it] * xs, label="", line=0, lw=0.5, marker=:circle, ms=1.2, seriestype=:scatter, alpha=0.5)
        end
        title!(string("Z=", Z, " ,", tit, " eV, time= ", trunc(times[1+(it-1)*skip_path], digits=4), "Gyr"), fonts=20)
        yaxis!("[Mpc]", (0, dx * n * 1e-3), fonts=20)
        xaxis!("[Mpc]", (0, dx * n * 1e-3), fonts=20)
    end
    file_gif = string(root_out, "_UHECR_path_color_map_", E_initial[1], "_Z_", Z, "_BC_redshift_newU_cosmo", cosmo, "_", tag, ".gif")
    gif(anim, file_gif, fps=10)

end


@everywhere function plot_spectrum(path, np, root_out, E_initial, Z, cosmo, tag)

    #....plotting the final spectrum 
    nspec = 10
    mie = 16.0   #...logE
    mae = 21.0
    bie = (mae - mie) / nspec
    xe = Array{Float64}(undef, nspec)
    spec = Array{Float64}(undef, nspec, npath)
    spec .= 0.0
    iobs = npath  #...snapshot where we produce the final spectrum 
    xe .= 0.0
    @inbounds for i in 1:nspec
        xe[i] = 10.0^(mie + bie * (i - 0.5))
    end

    @inbounds for j in 1:npath
        @inbounds @simd for i in 1:np
            en = path[i, 4, j]
            if en > 10.0^mie && en < 10.0^mae
                ie = convert(Int64, trunc((log10(en) - mie) / bie)) + 1
                spec[ie, j] += 1.0
            end
        end
    end

    e1 = 17.0
    e2 = 21.0
    plot(xe, spec[:, 1], line=:solid, dpi=300, lw=1, grid=false, alpha=0.5, label=string("z=", trunc(path[1, 9, 1], digits=6)), color="black")
    plot!(xe, spec[:, iobs], line=:solid, lw=2, alpha=0.5, label=string("z=0"), color="red")
    yaxis!("N(E)", :log10, (1, np), fonts=20)
    xaxis!(L"log_{10}E[eV]", :log10, (10^e1, 3 * 10^e2), fonts=20, xticks=[1e17, 1e18, 1e19, 1e20, 1e21])

    file = string(root_out, "_UHECR_spectum_", E_initial[1], "_Z_", Z, "_newU_cosmo", cosmo, "_", tag, ".png")
    savefig(file)


end


#....routines which writes the trajectory and energy evolution of all particles (once every 'skip_path snapshots) on an HDF5 file
@everywhere function write_path(root_out, E_initial, Z, cosmo, tag, path, t, it_path, inj)


    filep1 = string(root_out, "path_", E_initial[1], "Z", Z, "_spec_redshift_cosmo", cosmo, "_", tag, ".hdf5")
    #...we delete the file hdf5 if existing already
    command = `rm -r -f `
    run(`$command $filep1 `)
    h5write(filep1, "injections", inj)
    h5write(filep1, "snapshot", t)  #...timestep of the simulation
    h5write(filep1, "path step", it_path)   #..step in the path file (which uses a reduced number of steps to spare memory)
    h5write(filep1, "px_pbc", path[:, 1, :])  #...X position with periodic BC 
    h5write(filep1, "py_pbc", path[:, 2, :])  #...Y position with periodic BC 
    h5write(filep1, "pz_pbc", path[:, 3, :])  #...Z position with periodic BC 
    h5write(filep1, "px", path[:, 5, :])  #...X position without periodic BC 
    h5write(filep1, "py", path[:, 6, :])  #...Y position without periodic BC 
    h5write(filep1, "pz", path[:, 7, :])  #...Z position without periodic BC 
    h5write(filep1, "E[eV]", path[:, 4, :])  #...energy 
    h5write(filep1, "B[G]", path[:, 8, :])  #...physical mag.field in Gauss
    h5write(filep1, "redshift", path[:, 9, :])  #...redshift
    h5write(filep1, "time", path[:, 10, :])  #...elapsed time  
end


#...energy evolution of alla particles 
@everywhere function plot_energy(path, np, root_out, tag)

    plot(path[1, 10, :], path[1, 4, :], line=:solid, dpi=300, lw=1, grid=false, alpha=0.5, label=string("Energy"), color="black")
    @inbounds for i in 1:np
        plot!(path[i, 10, :], path[i, 4, :], lw=0.5, label="", alpha=0.6)
    end
    yaxis!("E[eV]", :log10, (1e17, 1e22), fonts=12)
    xaxis!(L"time", (0, 4e17), fonts=12)

    file = string(root_out, "_Energy_time_", tag, ".png")
    savefig(file)

end


#...plot of average distance vs time for all particles 
@everywhere function plot_distance(path, it_path, tag, root_out)
    ns = size(path)
    Dist = Array{Float64}(undef, ns[3], 2)
    Dist .= 0
    @inbounds for t in 1:ns[3]
        @inbounds for i in 1:ns[1]
            xd = (path[i, 5, 1] - path[i, 5, t])^2
            yd = (path[i, 6, 1] - path[i, 6, t])^2
            zd = (path[i, 7, 1] - path[i, 7, t])^2
            Dist[t, 1] += sqrt(xd + yd + zd)
            Dist[t, 2] += 1.0
        end
    end

    Dist[:, 1] ./= Dist[:, 2]
    plot(Dist[:, 1], line=:solid, dpi=300, grid=false, alpha=0.5, label=string("average distance"), color="black")
    file = string(root_out, "_Distance_time_", tag, ".png")
    savefig(file)

end