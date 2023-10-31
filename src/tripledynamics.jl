
using Syzygy, WGLMakie, Unitful, UnitfulAstro, LinearAlgebra
using DataStructures: CircularBuffer
WGLMakie.set_theme!(theme_black())

function animate_2d()
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 2000
    es = [[0.1, 0.1], [0.1, 0.5], [0.1, 0.9]]
    smas = [[0.1, 2.0]u"AU", [0.3, 2.0]u"AU", [0.6, 2.0]u"AU", [0.8, 2.0]u"AU"]
    params = [[es[i], smas[j]] for i in eachindex(es), j in eachindex(smas)]
    positions = zeros(size(params)..., 3, 3, N)

    for i in axes(params, 1)
        for j in axes(params, 2)
            triple = multibodysystem(masses, a=params[i, j][2], e=params[i, j][1])
            res = simulate(triple, t_sim=2, npoints=N, callbacks=[])
            sol = analyse_simulation(res)
            positions[i, j, :, :, :] = ustrip.(u"AU", sol.r)
        end
    end

    fig = Figure(resolution=(1920, 1080))
    axs = [Axis(fig[i, j], aspect=1) for i in axes(params, 1), j in axes(params, 2)]
    # axs = Axis(fig[1, 1])

    for ax in axs
        hidedecorations!(ax)
        hidespines!(ax)
        xlims!(ax, -2, 2)
        ylims!(ax, -2, 2)
    end

    frames = 2:N

    # r1s = [Observable(Point2f[]) for i in axes(params, 1), j in axes(params, 2)]
    # r2s = [Observable(Point2f[]) for i in axes(params, 1), j in axes(params, 2)]
    # r3s = [Observable(Point2f[]) for i in axes(params, 1), j in axes(params, 2)]

    sc1s = [Observable{Point2f}() for i in axes(params, 1), j in axes(params, 2)]
    sc2s = [Observable{Point2f}() for i in axes(params, 1), j in axes(params, 2)]
    sc3s = [Observable{Point2f}() for i in axes(params, 1), j in axes(params, 2)]

    
    nt = N÷100*20
    r1s = [CircularBuffer{Point2f}(nt) for i in axes(params, 1), j in axes(params, 2)]
    r2s = [CircularBuffer{Point2f}(nt) for i in axes(params, 1), j in axes(params, 2)]
    r3s = [CircularBuffer{Point2f}(nt) for i in axes(params, 1), j in axes(params, 2)]

    for i in CartesianIndices(axs)
        fill!(r1s[i], positions[i, 1:2, 1, 1])
        fill!(r2s[i], positions[i, 1:2, 2, 1])
        fill!(r3s[i], positions[i, 1:2, 3, 1])
    end

    r1s = Observable.(r1s)
    r2s = Observable.(r2s)
    r3s = Observable.(r3s)

    colors = Makie.wong_colors()[[1, 2, 3]]
    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]

    for i in CartesianIndices(axs)
        lines!(axs[i], r1s[i], color=trailcolors[1], linewidth=2.5)
        lines!(axs[i], r2s[i], color=trailcolors[2], linewidth=2.5)
        lines!(axs[i], r3s[i], color=trailcolors[3], linewidth=2.5)
        scatter!(axs[i], sc1s[i], color=colors[1], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))
        scatter!(axs[i], sc2s[i], color=colors[2], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))
        scatter!(axs[i], sc3s[i], color=colors[3], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))

    end

    record(fig, "figures/triple_animation.mp4", frames; framerate = N÷25) do frame

        for i in CartesianIndices(axs)
            
            push!(r1s[i][], positions[i, 1:2, 1, frame])
            push!(r2s[i][], positions[i, 1:2, 2, frame])
            push!(r3s[i][], positions[i, 1:2, 3, frame])

            sc1s[i][] = positions[i, 1:2, 1, frame] |> Point2f
            sc2s[i][] = positions[i, 1:2, 2, frame] |> Point2f
            sc3s[i][] = positions[i, 1:2, 3, frame] |> Point2f

            
            notify(r1s[i])
            notify(r2s[i])
            notify(r3s[i])


            notify(sc1s[i])
            notify(sc2s[i])
            notify(sc3s[i])
        end

    end
end


function animate_3d()
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 2000
    e = 0.1
    smas = [0.3, 2.0]u"AU"
    incs = [[0.0, 0.0]u"rad", [π/4, 0.0]u"rad", [π/2, 0.0]u"rad", [π, 0.0]u"rad"]
    incs = reshape(incs, 2, 2)
    positions = zeros(size(incs)..., 3, 3, N)

    for i ∈ CartesianIndices(incs)
        triple = multibodysystem(masses, a=smas, e=e, i = incs[i])
        res = simulate(triple, t_sim=2, npoints=N, callbacks=[])
        sol = analyse_simulation(res)
        positions[i, :, :, :] = ustrip.(u"AU", sol.r)
    end

    fig = Figure(resolution=(1920, 1080))
    colors = [:red, :cyan, :yellow]
    axs = axs = [Axis3(fig[i, j], xticklabelsvisible=false, 
                                  yticklabelsvisible=false, 
                                  zticklabelsvisible=false) for i in axes(incs, 1), j in axes(incs, 2)]

    for ax in axs
        # hidedecorations!(ax)
        hidespines!(ax)
        xlims!(ax, -2, 2)
        ylims!(ax, -2, 2)
        zlims!(ax, -2, 2)
    end

    frames = 2:N

    nt = N÷100*20
    r1s = [CircularBuffer{Point3f}(nt) for i in axes(incs, 1), j in axes(incs, 2)]
    r2s = [CircularBuffer{Point3f}(nt) for i in axes(incs, 1), j in axes(incs, 2)]
    r3s = [CircularBuffer{Point3f}(nt) for i in axes(incs, 1), j in axes(incs, 2)]

    for i in CartesianIndices(axs)
        fill!(r1s[i], positions[i, :, 1, 1])
        fill!(r2s[i], positions[i, :, 2, 1])
        fill!(r3s[i], positions[i, :, 3, 1])
    end

    r1s = Observable.(r1s)
    r2s = Observable.(r2s)
    r3s = Observable.(r3s)

    colors = Makie.wong_colors()[[1, 2, 3]]
    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]

    sc1s = [Observable{Point3f}() for i in axes(incs, 1), j in axes(incs, 2)]
    sc2s = [Observable{Point3f}() for i in axes(incs, 1), j in axes(incs, 2)]
    sc3s = [Observable{Point3f}() for i in axes(incs, 1), j in axes(incs, 2)]


    for i in CartesianIndices(axs)
        lines!(axs[i], r1s[i], color=trailcolors[1], linewidth=2.5)
        lines!(axs[i], r2s[i], color=trailcolors[2], linewidth=2.5)
        lines!(axs[i], r3s[i], color=trailcolors[3], linewidth=2.5)
        scatter!(axs[i], sc1s[i], color=colors[1], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))
        scatter!(axs[i], sc2s[i], color=colors[2], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))
        scatter!(axs[i], sc3s[i], color=colors[3], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))

    end

    record(fig, "figures/triple_animation_inclination.mp4", frames; framerate = N÷25) do frame

        for i in CartesianIndices(axs)
            push!(r1s[i][], positions[i, :, 1, frame])
            push!(r2s[i][], positions[i, :, 2, frame])
            push!(r3s[i][], positions[i, :, 3, frame])

            sc1s[i][] = positions[i, :, 1, frame] |> Point3f
            sc2s[i][] = positions[i, :, 2, frame] |> Point3f
            sc3s[i][] = positions[i, :, 3, frame] |> Point3f

            notify(r1s[i])
            notify(r2s[i])
            notify(r3s[i])

            notify(sc1s[i])
            notify(sc2s[i])
            notify(sc3s[i])

            axs[i].azimuth[] = 2π*sin(2frame/N)
            axs[i].elevation[] = 0.1π*sin(2frame/N)
        end
    end
end

function animate_kozailidov()
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 10000
    e = 0.1
    i = [π/2, 0.0]u"rad"
    smas = [0.3, 2.0]u"AU"
    positions = zeros(3, 3, N)


    triple = multibodysystem(masses, a=smas, e=e, i=i)
    res = simulate(triple, t_sim=100, npoints=N, callbacks=[])
    sol = analyse_simulation(res)
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    

    fig = Figure(resolution=(1080, 1080))
    colors = [:red, :cyan, :yellow]
    ax = Axis3(fig[1, 1], xticklabelsvisible=false, 
                                 yticklabelsvisible=false, 
                                 zticklabelsvisible=false)

    hidespines!(ax)
    xlims!(ax, -2, 2)
    ylims!(ax, -2, 2)
    zlims!(ax, -2, 2)

    frames = 1:N

    r1s = Observable(Point3f[])
    r2s = Observable(Point3f[])
    r3s = Observable(Point3f[])


    lines!(ax, r1s, color=colors[1], linewidth=1)
    lines!(ax, r2s, color=colors[2], linewidth=1)
    lines!(ax, r3s, color=colors[3], linewidth=1)
    
    # lines!(axs[1, 1], seco[1], seco[2], color=colors[2], linewidth=1)
    # lines!(axs[1, 1], tert[1], tert[2], color=colors[3], linewidth=1)

    record(fig, "figures/triple_animation_kozailidov.mp4", frames; framerate = N÷15) do frame


        push!(r1s[], positions[:, 1, frame])
        push!(r2s[], positions[:, 2, frame])
        push!(r3s[], positions[:, 3, frame])
        notify(r1s)
        notify(r2s)
        notify(r3s)
            # axs[i].azimuth[] = 2π*sin(frame/N)

    end
end

function animate_kozailidov_inner_bin()
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 10000
    e = 0.1
    i = [π/2, 0.0]u"rad"
    smas = [0.3, 2.0]u"AU"
    positions = zeros(3, 3, N)


    triple = multibodysystem(masses, a=smas, e=e, i=i)
    res = simulate(triple, t_sim=30, npoints=N, callbacks=[])
    sol = analyse_simulation(res)
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    com_in = map(i -> Syzygy.centre_of_mass(sol.r[:, 1:2, i], masses[1:2]), eachindex(sol.t))
    com_in = map(x -> ustrip.(u"AU", x), com_in)
    fig = Figure(resolution=(1080, 1080))
    colors = [(:red, 0.2), (:cyan, 0.2), :yellow]
    ax = Axis(fig[1, 1], xticklabelsvisible=false, 
                                 yticklabelsvisible=false)

    hidespines!(ax)
    xlims!(ax, -0.5, 0.5)
    ylims!(ax, -0.5, 0.5)
    # zlims!(ax, -0.5, 0.5)

    frames = 1:N

    r1s = Observable(Point2f[])
    r2s = Observable(Point2f[])
    # r3s = Observable(Point3f[])


    lines!(ax, r1s, color=colors[1], linewidth=1)
    lines!(ax, r2s, color=colors[2], linewidth=1)
    # lines!(ax, r3s, color=colors[3], linewidth=1)
    
    # lines!(axs[1, 1], seco[1], seco[2], color=colors[2], linewidth=1)
    # lines!(axs[1, 1], tert[1], tert[2], color=colors[3], linewidth=1)

    record(fig, "figures/triple_animation_kozailidov_inner_2.mp4", frames; framerate = N÷20) do frame


        push!(r1s[], positions[[1,3], 1, frame] .- com_in[frame][[1,3]])
        push!(r2s[], positions[[1,3], 2, frame] .- com_in[frame][[1,3]])
        # push!(r3s[], positions[:, 3, frame])
        notify(r1s)
        notify(r2s)
        # notify(r3s)
        # ax.azimuth[] = 2π*sin(frame/N)

    end
end

function plot_orbital_elements()
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 10000
    e = 0.1
    i = [π/2, 0.0]u"rad"
    smas = [0.3, 2.0]u"AU"
    step = 10
    # positions = zeros(3, 3, N)


    triple = multibodysystem(masses, a=smas, e=e, i=i)
    res = simulate(triple, t_sim=200, npoints=N, callbacks=[])
    sol = analyse_simulation(res)
    # return sol
    r = sol.r[:,:,1:step:end]
    v = sol.v[:,:,1:step:end]
    t = sol.t[1:step:end] ./ triple.binaries[2].elements.P .|> upreferred
    # positions[:, :, :] = ustrip.(u"AU", sol.r)
    

    fig = Figure()#resolution=(1080, 1080)
    colors = [:red, :cyan, :yellow]
    
    # ax_aout = Axis(fig[1, 1], title=L"a_{{\text{out}}} \quad [\text{R}_\odot]", 
    #                titlesize=30, xticklabelsvisible=false)
    # ax_ain  = Axis(fig[1, 2], title=L"a_{{\text{in}}} \quad [\text{R}_\odot]", 
    #                titlesize=30, xticklabelsvisible=false)
    ax_eout = Axis(fig[1, 1], title=L"e_{{\text{out}}}", 
                   titlesize=30, xticklabelsvisible=false)
    ax_ein  = Axis(fig[1, 2], title=L"e_{{\text{in}}}", 
                   titlesize=30, xticklabelsvisible=false)
    ax_imut = Axis(fig[2,1:2], title=L"i_{{\text{mut}}} \quad [°]", 
                   titlesize=30, xlabel="# outer orbits")

    # hidespines!(ax)
    # xlims!(ax, -2, 2)
    ylims!(ax_eout, 0, 1)
    # zlims!(ax, -2, 2)

    com_in = map(i -> Syzygy.centre_of_mass(r[:, 1:2, i], masses[1:2]), eachindex(t))
    v_com_in = map(i -> Syzygy.centre_of_mass_velocity(v[:, 1:2, i], masses[1:2]), eachindex(t))
    
    com_in = reduce(hcat, com_in)
    v_com_in = reduce(hcat, v_com_in)

    r_rel = r[particle=3] .- com_in
    r_in = r[particle=2] .- r[particle=1]

    v_rel = v[particle=3] .- v_com_in
    v_in = v[particle=2] .- v[particle=1]
    
    d_rel = norm.(eachcol(r_rel))
    v²_rel = norm.(eachcol(v_rel)) .^ 2

    d_in = norm.(eachcol(r_in))
    v²_in = norm.(eachcol(v_in)) .^ 2

    a_out = Syzygy.semi_major_axis.(d_rel, v²_rel, sum(masses)) 
    a_in = Syzygy.semi_major_axis.(d_in, v²_in, sum(masses)) 

    e_out = Syzygy.eccentricity.(eachcol(r_rel), eachcol(v_rel), a_out, sum(masses))
    e_in = Syzygy.eccentricity.(eachcol(r_in), eachcol(v_in), a_in, sum(masses[1:2]))

    h_in = Syzygy.angular_momentum.(eachcol(r_in), eachcol(v_in))

    i_in = Syzygy.inclination.(h_in)

    h_out = Syzygy.angular_momentum.(eachcol(r_rel), eachcol(v_rel))
    i_out = Syzygy.inclination.(h_out)


    i_mut = i_out + i_in

    # scatter!(ax_aout, t, ustrip.(u"Rsun", a_out), markersize=4)
    # scatter!(ax_ain, t, ustrip.(u"Rsun", a_in), markersize=4)
    scatter!(ax_eout, t, ustrip.(e_out), markersize=4)
    scatter!(ax_ein, t, ustrip.(e_in), markersize=4)
    scatter!(ax_imut, t, ustrip.(u"°", i_mut), markersize=4)
    save("figures/orbital_elements.png", fig)
    # fig

end


function animate_destablization(;showplot=false, kwargs...)
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 5000
    e = 0.1
    i = [0.0, 0.0]u"rad"
    a = [0.8, 2.0]u"AU"
    positions = zeros(3, 3, N)


    triple = multibodysystem(masses, a=a, e=e, i=i)
    res = simulate(triple, t_sim=10, npoints=N, callbacks=[])
    sol = analyse_simulation(res)
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    
    if showplot
        figg = Figure()
        axx = Axis(figg[1, 1], aspect=1)
        lines!(axx, positions[1,1,:],  positions[2,1,:])
        lines!(axx, positions[1,2,:],  positions[2,2,:])
        lines!(axx, positions[1,3,:],  positions[2,3,:])
        return figg
    end

    fig = Figure(resolution=(1080, 1080))
    colors = [:red, :cyan, :yellow]
    ax = Axis(fig[1, 1], xticklabelsvisible=false, 
                         yticklabelsvisible=false)

    hidespines!(ax)
    xlims!(ax, -3, 3)
    ylims!(ax, -3, 3)

    frames = 2:N

    nt = N÷100*20
    r1s = CircularBuffer{Point2f}(nt)
    r2s = CircularBuffer{Point2f}(nt)
    r3s = CircularBuffer{Point2f}(nt)

    fill!(r1s, positions[1:2, 1, 1])
    fill!(r2s, positions[1:2, 2, 1])
    fill!(r3s, positions[1:2, 3, 1])

    r1s = Observable(r1s)
    r2s = Observable(r2s)
    r3s = Observable(r3s)

    colors = Makie.wong_colors()[[1, 2, 3]]
    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]

    sc1s = Observable{Point2f}()
    sc2s = Observable{Point2f}()
    sc3s = Observable{Point2f}()
    
    lines!(ax, r1s, color=trailcolors[1], linewidth=2.5)
    lines!(ax, r2s, color=trailcolors[2], linewidth=2.5)
    lines!(ax, r3s, color=trailcolors[3], linewidth=2.5)
    scatter!(ax, sc1s, color=colors[1], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))
    scatter!(ax, sc2s, color=colors[2], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))
    scatter!(ax, sc3s, color=colors[3], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))


    record(fig, "figures/triple_animation_destabilize.mp4", frames; framerate = N÷25) do frame


        push!(r1s[], positions[1:2, 1, frame])
        push!(r2s[], positions[1:2, 2, frame])
        push!(r3s[], positions[1:2, 3, frame])

        sc1s[] = positions[1:2, 1, frame] |> Point2f
        sc2s[] = positions[1:2, 2, frame] |> Point2f
        sc3s[] = positions[1:2, 3, frame] |> Point2f

        notify(r1s)
        notify(r2s)
        notify(r3s)

        notify(sc1s)
        notify(sc2s)
        notify(sc3s)

            # axs[i].azimuth[] = 2π*sin(frame/N)

    end
end
