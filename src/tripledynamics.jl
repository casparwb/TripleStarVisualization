
using Syzygy, LinearAlgebra, ProgressMeter
using DataStructures: CircularBuffer
using JLD2
using GLMakie
# using CairoMakie
# CairoMakie.set_theme!(theme_black())


const FIGPATH = joinpath(@__DIR__, "..", "figures")
function many_triples_2d(;outname="many_triples_2d")
    masses = [1.0, 2.0, 1.0]u"Msun"
    N = 2000
    es = [[0.1, 0.1], [0.1, 0.5], [0.1, 0.9]]
    smas = [[0.1, 3.0]u"AU", [0.3, 3.0]u"AU", [0.6, 3.0]u"AU", [0.8, 2.0]u"AU"]
    params = [[es[i], smas[j]] for i in eachindex(es), j in eachindex(smas)]
    positions = zeros(size(params)..., 3, 3, N)

    for i in axes(params, 1)
        for j in axes(params, 2)
            triple = multibodysystem(masses, a=params[i, j][2], e=params[i, j][1])
            res = simulate(triple, t_sim=2, npoints=N, callbacks=[])
            sol = to_solution(res)
            positions[i, j, :, :, :] = ustrip.(u"AU", sol.r)
        end
    end

    fig = Figure(size=(1920, 1080))
    axs = [Axis(fig[i, j], aspect=1) for i in axes(params, 1), j in axes(params, 2)]
    # axs = Axis(fig[1, 1])
    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)
    for ax in axs
        hidedecorations!(ax)
        hidespines!(ax)
        xlims!(ax, -3, 3)
        ylims!(ax, -3, 3)
    end

    frames = 2:N

    # r1s = [Observable(Point2f[]) for i in axes(params, 1), j in axes(params, 2)]
    # r2s = [Observable(Point2f[]) for i in axes(params, 1), j in axes(params, 2)]
    # r3s = [Observable(Point2f[]) for i in axes(params, 1), j in axes(params, 2)]

    sc1s = [Observable{Point2f}() for i in axes(params, 1), j in axes(params, 2)]
    sc2s = [Observable{Point2f}() for i in axes(params, 1), j in axes(params, 2)]
    sc3s = [Observable{Point2f}() for i in axes(params, 1), j in axes(params, 2)]

    
    nt = N÷100*40
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
        lines!(axs[i], r1s[i], color=trailcolors[1], linewidth=4.0)
        lines!(axs[i], r2s[i], color=trailcolors[2], linewidth=4.0)
        lines!(axs[i], r3s[i], color=trailcolors[3], linewidth=4.0)
        scatter!(axs[i], sc1s[i], color=colors[1], marker=:star5, markersize=20,
                 colormap=Makie.wong_colors(), colorrange=(1, 3))
        scatter!(axs[i], sc2s[i], color=colors[2], marker=:star5, markersize=20,
                 colormap=Makie.wong_colors(), colorrange=(1, 3))
        scatter!(axs[i], sc3s[i], color=colors[3], marker=:star5, markersize=20,
                 colormap=Makie.wong_colors(), colorrange=(1, 3))
    end

    p = Progress(N)
    savepath = joinpath(FIGPATH, outname)*".mp4"
    record(fig, savepath, frames; framerate = N÷25) do frame

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

        next!(p)
    end
end


function mutual_inclination(;outname = "mutual_inclination")
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 2000
    e = 0.1
    smas = [0.5, 2.0]u"AU"
    incs = [[0.0, 0.0]u"rad", [π/4, 0.0]u"rad", [π/2, 0.0]u"rad", [π, 0.0]u"rad"]
    incs = reshape(incs, 2, 2)
    positions = zeros(size(incs)..., 3, 3, N)

    for i ∈ CartesianIndices(incs)
        triple = multibodysystem(masses, a=smas, e=e, i = incs[i])
        res = simulate(triple, t_sim=2, npoints=N, callbacks=[])
        sol = to_solution(res)
        positions[i, :, :, :] = ustrip.(u"AU", sol.r)
    end

    fig = Figure(size=(1920, 1080))
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

    nt = N÷100*30
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
        lines!(axs[i], r1s[i], color=trailcolors[1], linewidth=4.0)
        lines!(axs[i], r2s[i], color=trailcolors[2], linewidth=4.0)
        lines!(axs[i], r3s[i], color=trailcolors[3], linewidth=4.0)
        scatter!(axs[i], sc1s[i], color=colors[1], marker=:star5, markersize=20,
                 colormap=Makie.wong_colors(), colorrange=(1, 3))
        scatter!(axs[i], sc2s[i], color=colors[2], marker=:star5, markersize=20,
                 colormap=Makie.wong_colors(), colorrange=(1, 3))
        scatter!(axs[i], sc3s[i], color=colors[3], marker=:star5, markersize=20,
                 colormap=Makie.wong_colors(), colorrange=(1, 3))

    end

    p = Progress(N)
    savepath = joinpath(FIGPATH, outname)*".mp4"
    record(fig, savepath, frames; framerate = N÷25) do frame

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

            axs[i].azimuth[] = 2π*cos(2frame/N)
            axs[i].elevation[] = 0.1π*sin(2frame/N)
        end
        next!(p)
    end
end

function kozai_lidov(;outname = "kozai_lidov")
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 30000
    e = 0.1
    i = [π/2, 0.0]u"rad"
    smas = [0.3, 2.0]u"AU"
    positions = zeros(3, 3, N)


    triple = multibodysystem(masses, a=smas, e=e, i=i)
    res = simulate(triple, t_sim=60, npoints=N, callbacks=[])
    sol = to_solution(res)
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    

    fig = Figure(size=(1080, 1080))
    colors = [:red, :cyan, :yellow]
    # colors = Makie.wong_colors()
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

    p = Progress(N)
    savepath = joinpath(FIGPATH, outname)*".mp4"
    record(fig, savepath, frames; framerate = N÷20) do frame


        push!(r1s[], positions[:, 1, frame])
        push!(r2s[], positions[:, 2, frame])
        push!(r3s[], positions[:, 3, frame])
        notify(r1s)
        notify(r2s)
        notify(r3s)

        next!(p)
            # axs[i].azimuth[] = 2π*sin(frame/N)

    end
end

function kozai_lidov_inner_binary(;outname="kozai_lidov_inner_binary")
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 10000
    e = 0.1
    i = [π/2, 0.0]u"rad"
    smas = [0.3, 2.0]u"AU"
    positions = zeros(3, 3, N)

    triple = multibodysystem(masses, a=smas, e=e, i=i)
    res = simulate(triple, t_sim=30, npoints=N, callbacks=[])
    sol = to_solution(res)
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    com_in = map(i -> Syzygy.centre_of_mass(sol.r[:, 1:2, i], masses[1:2]), eachindex(sol.t))
    com_in = map(x -> ustrip.(u"AU", x), com_in)

    P_outs = sol.t ./ triple.binaries[2].elements.P .|> upreferred

    fig = Figure(size=(1920, 1080))
    colors = [(:red, 0.2), (:cyan, 0.2), :yellow]
    ax = Axis(fig[1:2, 2], xticklabelsvisible=false, 
                           yticklabelsvisible=false)
    ax_ecc = Axis(fig[1, 1])
    ax_inc = Axis(fig[2, 1])

    xlims!(ax_ecc, 0, P_outs[end])
    ylims!(ax_ecc, 0, 1)

    xlims!(ax_inc, 0, P_outs[end])
    ylims!(ax_inc, 0, 100)

    hidespines!(ax)
    xlims!(ax, -0.5, 0.5)
    ylims!(ax, -0.5, 0.5)

    frames = 2:N

    nt = N÷100*10
    r1s = CircularBuffer{Point2f}(nt)
    r2s = CircularBuffer{Point2f}(nt)

    fill!(r1s, positions[[1,3], 1, 1] .- com_in[1][[1, 3]])
    fill!(r2s, positions[[1,3], 2, 1] .- com_in[1][[1, 3]])
    
    r1s = Observable(r1s)
    r2s = Observable(r2s)  


    e_ins = Observable(Point2f[])
    i_muts = Observable(Point2f[])
    
    e_in, i_mut = get_orbital_elements(sol.r[:,:,1], sol.v[:,:,1], masses)

    i_mut = ustrip.(u"°", i_mut)

    colors = Makie.wong_colors()[[1, 2]]
    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]

    lines!(ax, r1s, color=trailcolors[1], linewidth=1)
    lines!(ax, r2s, color=trailcolors[2], linewidth=1)

    scatter!(ax_ecc, Point2f(P_outs[1], e_in))
    scatter!(ax_inc, Point2f(P_outs[1], i_mut))

    savepath = joinpath(FIGPATH, outname)*".mp4"

    p = Progress(N)
    record(fig, savepath, frames; framerate = N÷20) do frame

        push!(r1s[], positions[[1,3], 1, frame] .- com_in[frame][[1,3]])
        push!(r2s[], positions[[1,3], 2, frame] .- com_in[frame][[1,3]])
        
        notify(r1s)
        notify(r2s)

        e_in, i_mut = get_orbital_elements(sol.r[:,:,frame], sol.v[:,:,frame], masses)
        push!(e_ins[], Point2f(P_outs[frame], e_in))
        push!(i_muts[], Point2f(P_outs[frame], i_mut))

        notify(e_ins)
        notify(i_muts)
        next!(p)
    end
end

function get_orbital_elements(r, v, masses)
    com_in = Syzygy.centre_of_mass(r[:, 1:2], masses[1:2])
    v_com_in = Syzygy.centre_of_mass_velocity(v[:, 1:2], masses[1:2])
    
    # com_in = reduce(hcat, com_in)
    # v_com_in = reduce(hcat, v_com_in)

    M = sum(masses)
    M12 = masses[1] + masses[2]

    r_rel = r[particle=3] - com_in
    r_in = r[particle=2] - r[particle=1]

    v_rel = v[particle=3] - v_com_in
    v_in = v[particle=2] - v[particle=1]
    
    # d_rel = norm(r_rel)
    # v²_rel = norm(v_rel) ^ 2

    d_in = norm(r_in)
    v²_in = norm(v_in) ^ 2

    # a_out = Syzygy.semi_major_axis(d_rel, v²_rel, M) 
    a_in = Syzygy.semi_major_axis(d_in, v²_in, M) 

    # e_out = Syzygy.eccentricity(r_rel, v_rel, a_out, M)
    e_in = Syzygy.eccentricity(r_in, v_in, a_in, M12)

    h_in = Syzygy.angular_momentum(r_in, v_in)

    h_out = Syzygy.angular_momentum(r_rel, v_rel)

    i_mut =  Syzygy.mutual_inclination(h_in, h_out)

    return e_in, i_mut
end

function get_inner_eccentricity_from_eccentricity_vector(sol)

    r, v = sol.r, sol.v

    r_in = r[particle=2] .- r[particle=1]
    v_in = v[particle=2] .- v[particle=1]

    d_in = norm.(eachcol(r_in))

    n = eachcol(r_in) ./ d_in
    μ = GRAVCONST*(sol.ic.particles.mass[1] + sol.ic.particles.mass[2])


    e_in = zeros(Float64, length(sol.t))
    for i in eachindex(sol.t)
        e_vec = (v_in[:,i] × (r_in[:,i] × v_in[:,i]))/μ - n[i] |> norm
        e_in[i] = e_vec
    end

    return e_in
    # return Syzygy.eccentricity_vector.(eachcol(r_in), eachcol(d_in), eachcol(v_in), μ)
end

function get_inner_eccentricity_from_eccentricity_vector(r, v, M12)

    r_in = r[particle=2] - r[particle=1]
    v_in = v[particle=2] - v[particle=1]

    d_in = norm(r_in)

    n = r_in / d_in
    μ = GRAVCONST*M12


    (v_in × (r_in × v_in))/μ - n |> norm
    # return Syzygy.eccentricity_vector.(eachcol(r_in), eachcol(d_in), eachcol(v_in), μ)
end

function get_imut_and_ein(r, v, masses)
    com_in = Syzygy.centre_of_mass(r[:, 1:2], masses[1:2])
    v_com_in = Syzygy.centre_of_mass_velocity(v[:, 1:2], masses[1:2])
    
    # com_in = reduce(hcat, com_in)
    # v_com_in = reduce(hcat, v_com_in)

    # M = sum(masses)
    M12 = masses[1] + masses[2]

    r_rel = r[particle=3] - com_in
    r_in = r[particle=2] - r[particle=1]

    v_rel = v[particle=3] - v_com_in
    v_in = v[particle=2] - v[particle=1]
    
    # d_rel = norm(r_rel)
    # v²_rel = norm(v_rel) ^ 2

    # d_in = norm(r_in)
    # v²_in = norm(v_in) ^ 2

    # a_out = Syzygy.semi_major_axis(d_rel, v²_rel, M) 
    # a_in = Syzygy.semi_major_axis(d_in, v²_in, M) 

    # e_out = Syzygy.eccentricity(r_rel, v_rel, a_out, M)
    # e_in = Syzygy.eccentricity(r_in, v_in, a_in, M12)
    e_in = get_inner_eccentricity_from_eccentricity_vector(r, v, M12)

    h_in = Syzygy.angular_momentum(r_in, v_in)

    h_out = Syzygy.angular_momentum(r_rel, v_rel)

    i_mut =  Syzygy.mutual_inclination(h_in, h_out)

    return e_in, i_mut
end

function kozai_lidov_full(;outname="kozai_lidov_full")
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 2000#8000#20000
    e = 0.1
    i = [π/2, 0.0]u"rad"
    smas = [0.3, 2.0]u"AU"
    positions = zeros(3, 3, N)

    T = 20#60


    triple = multibodysystem(masses, a=smas, e=e, i=i)
    res = simulate(triple, t_sim=T, npoints=N, callbacks=[])
    sol = to_solution(res)
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    com_in = map(i -> Syzygy.centre_of_mass(sol.r[:, 1:2, i], masses[1:2]), eachindex(sol.t))
    com_in = map(x -> ustrip.(u"AU", x), com_in)
    P_outs = sol.t ./ triple.binaries[2].elements.P .|> upreferred


    fig = Figure(size=(1920, 1080))


    # colors = [:red, :cyan, :yellow]
    colors = Makie.wong_colors()

    # ax = Axis3(fig[1, 1],  xticklabelsvisible=false, 
    #                          yticklabelsvisible=false, 
    #                          zticklabelsvisible=false)
    ax_ecc = Axis(fig[1, 2], xticklabelsvisible=false, xgridvisible=false, ygridvisible=false, yticklabelsize=20, ylabel="Inner eccentricity", ylabelsize=30)
    ax_inc = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false, yticklabelsize=20, xticksvisible=false, xticklabelsvisible=false, ylabel="Inclination [°]", ylabelsize=30)
    ax_inn = Axis(fig[1:2, 1], xticklabelsvisible=false, yticksvisible=false, xticksvisible=false,
                             yticklabelsvisible=false, xgridvisible=false, ygridvisible=false, title="Inner orbit", titlesize=40)
    
    xlims!(ax_ecc, 0, P_outs[end])
    xlims!(ax_inc, 0, P_outs[end])
    ylims!(ax_ecc, -0.1, 1.1)
    ylims!(ax_inc, 20, 90)
    
    xlims!(ax_inn, -0.5, 0.5)
    ylims!(ax_inn, -0.5, 0.5)

    hidespines!(ax_inn)

    # return fig
    # hidespines!(ax)
    # xlims!(ax, -2, 2)
    # ylims!(ax, -2, 2)
    # zlims!(ax, -2, 2)

    frames = 2:N
    nt = N÷100*15

    # r1s = Observable(Point3f[])
    # r2s = Observable(Point3f[])
    # r3s = Observable(Point3f[])

    # lines!(ax, r1s, color=colors[1], linewidth=1)
    # lines!(ax, r2s, color=colors[2], linewidth=1)
    # lines!(ax, r3s, color=colors[3], linewidth=1)

    r1s_inn = CircularBuffer{Point2f}(nt)
    r2s_inn = CircularBuffer{Point2f}(nt)

    fill!(r1s_inn, positions[[1,3], 1, 1] .- com_in[1][[1, 3]])
    fill!(r2s_inn, positions[[1,3], 2, 1] .- com_in[1][[1, 3]])
    
    r1s_inn = Observable(r1s_inn)
    r2s_inn = Observable(r2s_inn)  
    
    e_ins = Observable(Point2f[(0.0, 0.0)])
    i_muts = Observable(Point2f[(0.0, 0.0)])
    
    e_in, i_mut = get_imut_and_ein(sol.r[:,:,1], sol.v[:,:,1], masses)#get_orbital_elements(sol.r[:,:,1], sol.v[:,:,1], masses)
    # e_prev = e_in

    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors[1:2]]

    lines!(ax_inn, r1s_inn, color=trailcolors[1], linewidth=1)
    lines!(ax_inn, r2s_inn, color=trailcolors[2], linewidth=1)

    # scatter!(ax_ecc, Point2f(P_outs[1], e_in))
    # scatter!(ax_inc, Point2f(P_outs[1], ustrip(u"°", i_mut)))

    # scatter!(ax_ecc, e_ins, color=:pink)
    # scatter!(ax_inc, i_muts, color=:pink)

    lines!(ax_ecc, e_ins, color=:pink)
    lines!(ax_inc, i_muts, color=:pink)

    savepath = joinpath(FIGPATH, outname)*".mp4"
    p = Progress(N)
    record(fig, savepath, frames; framerate = N÷25) do frame

        # push!(r1s[], positions[:, 1, frame])
        # push!(r2s[], positions[:, 2, frame])
        # push!(r3s[], positions[:, 3, frame])
        # notify(r1s)
        # notify(r2s)
        # notify(r3s)

        push!(r1s_inn[], positions[[1,3], 1, frame] .- com_in[frame][[1,3]])
        push!(r2s_inn[], positions[[1,3], 2, frame] .- com_in[frame][[1,3]])
        
        notify(r1s_inn)
        notify(r2s_inn)

        # if iszero(mod(frame, 10))
            e_in, i_mut = get_imut_and_ein(sol.r[:,:,frame], sol.v[:,:,frame], masses)
            # if abs(e_in - e_prev)
            # println(P_outs[frame], " ", e_in, " ", i_mut)
            e_ins[] = push!(e_ins[], Point2f(P_outs[frame], e_in))
            i_muts[] = push!(i_muts[], Point2f(P_outs[frame], ustrip(u"°", i_mut)))
        # end
        # notify(e_ins)
        # notify(i_muts)
        next!(p) 
    end

end

function plot_orbital_elements()
    masses = [1.0, 1.0, 1.0]u"Msun"
    # N = 10000
    e = [0.1, 0.4]
    i = [π/2, 0.0]u"rad"
    smas = [0.3, 2.0]u"AU"
    # positions = zeros(3, 3, N)

    triple = multibodysystem(masses, a=smas, e=e, i=i)
    P = triple.binaries[2].elements.P
    res = simulate(triple, t_sim=1000, saveat=upreferred(P).val, callbacks=[])
    sol = to_solution(res)
    # return sol
    r = sol.r[:,:,:]
    v = sol.v[:,:,:]
    t = sol.t ./ triple.binaries[2].elements.P .|> upreferred
    # positions[:, :, :] = ustrip.(u"AU", sol.r)
    

    fig = Figure(size=(1980, 1080))
    colors = [:red, :cyan, :yellow]
    
    # ax_aout = Axis(fig[1, 1], title=L"a_{{\text{out}}} \quad [\text{R}_\odot]", 
    #                titlesize=30, xticklabelsvisible=false)
    # ax_ain  = Axis(fig[1, 2], title=L"a_{{\text{in}}} \quad [\text{R}_\odot]", 
    #                titlesize=30, xticklabelsvisible=false)
    # ax_eout = Axis(fig[1, 1], title=L"e_{{\text{out}}}", 
    #                titlesize=30, xticklabelsvisible=false)
    # ax_ein  = Axis(fig[1, 2], title=L"e_{{\text{in}}}", 
    #                titlesize=30, xticklabelsvisible=false)
    # ax_imut = Axis(fig[2,1:2], title=L"i_{{\text{mut}}} \quad [°]", 
    #                titlesize=30, xlabel="# outer orbits")

    axe = Axis(fig[1, 1], yscale = log10, ylabel = "1 - eᵢₙ", ylabelsize=30)
    axi = Axis(fig[2, 1], ylabel = L"i_$\text{mut}$ [°]", ylabelsize=30,
                          xlabel="# outer orbits", xlabelsize=20, xticks=Makie.WilkinsonTicks(10))

    hidexdecorations!(axe)

    # hidespines!(ax)
    # xlims!(ax, -2, 2)

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

    h_out = Syzygy.angular_momentum.(eachcol(r_rel), eachcol(v_rel))

    i_mut =  Syzygy.mutual_inclination.(h_in, h_out)

    # scatter!(ax_eout, t, ustrip.(e_out), markersize=4)
    # scatter!(ax_ein, t, ustrip.(e_in), markersize=4)
    # scatter!(ax_imut, t, ustrip.(u"°", i_mut), markersize=4)
    # lines!(ax_eout, t, ustrip.(e_out))
    # lines!(ax_ein, t, ustrip.(e_in))
    # lines!(ax_imut, t, ustrip.(u"°", i_mut))

    # ylims!(axe, -0.1, 1.1)
    ylims!(axi, 20, 100)

    lines!(axe, t, 1 .- ustrip.(e_in))
    lines!(axi, t, ustrip.(u"°", i_mut), color=:red)


    save("figures/orbital_elements.png", fig)
    fig

end


function unstable_triple(;outname="unstable_triple", showplot=false, kwargs...)
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 10000
    # e = 0.2
    # i = [0.0, 0.0]u"rad"
    # a = [0.6, 2.0]u"AU"

    e = 0.4
    i = [0.0, 0.0]u"rad"
    a = [0.6, 2.0]u"AU"

    positions = zeros(3, 3, N)


    triple = multibodysystem(masses, a=a, e=e, i=i)
    res = simulate(triple, t_sim=25, npoints=N, callbacks=[])
    sol = to_solution(res)
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    
    if showplot
        figg = Figure()
        axx = Axis(figg[1, 1], aspect=1)
        lines!(axx, positions[1,1,:],  positions[2,1,:])
        lines!(axx, positions[1,2,:],  positions[2,2,:])
        lines!(axx, positions[1,3,:],  positions[2,3,:])
        return figg
    end

    fig = Figure(size=(1080, 1080))
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

    savepath = joinpath(FIGPATH, outname)*".mp4"
    p = Progress(N)
    record(fig, savepath, frames; framerate = N÷25) do frame


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
        next!(p)
    end
end


function collision_mesh_animation(T=6.242690855674927e8u"s"; showplot=false)
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 5_000

    R = [1.0, 1.0, 1.0]u"Rsun"
    """ Very quickly destabilizing system:"""
    # e = 0.1
    # i = [0.0, 0.0]u"rad"
    # a = [0.8, 2.0]u"AU"


    """ Collision system """
    e = [0.1, 0.6]
    i = [120.0, 0.0]u"°"
    a = [0.4, 2.0]u"AU"


    
    
    triple = multibodysystem(masses, a=a, e=e, i=i, R=R)
    res = simulate(triple, t_sim=T, npoints=N, callbacks=[])
    @show res.retcode
    sol = to_solution(res)
    N = length(sol.t)
    positions = zeros(3, 3, length(sol.t))
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    
    Rs = ustrip.(u"AU", R)
    com_12 = Syzygy.centre_of_mass(sol, [1, 2]) 
    com_12 = ustrip.(u"AU", com_12)
    if showplot
        figg = Figure()
        axx = Axis3(figg[1, 1], aspect=:equal)
    	# xlims!(axx, -a[1].val*5, a[1].val*5)
    	# ylims!(axx, -a[1].val*5, a[1].val*5)
    	# zlims!(axx, -a[1].val*5, a[1].val*5)

        # xlims!(axx, -0.05, 0.05)
        # ylims!(axx, -0.05, 0.05)
        # zlims!(axx, -0.05, 0.05)
        # lines!(axx, positions[1,1,end-100:end] .- com_12[1,end-100:end],  positions[2,1,end-100:end] .- com_12[2,end-100:end], positions[3,1,end-100:end] .- com_12[3,end-100:end])
        # lines!(axx, positions[1,2,end-100:end] .- com_12[1,end-100:end],  positions[2,2,end-100:end] .- com_12[2,end-100:end], positions[3,2,end-100:end] .- com_12[3,end-100:end])
        # lines!(axx, positions[1,3,:],  positions[2,3,:], positions[3,3,:])
        # save(figg, "figures/test.png")

        mesh!(axx, Sphere(Point3f(positions[:, 1, end] .- com_12[:,end]), Rs[1]))#color=colors[1], colormap=Makie.wong_colors(), colorrange=(1, 3))
        mesh!(axx, Sphere(Point3f(positions[:, 2, end] .- com_12[:,end]), Rs[2]))#color=colors[2], colormap=Makie.wong_colors(), colorrange=(1, 3))
        # mesh!(axx, Sphere(Point3f(positions[:, 3, end]), Rs[3]))#
        return figg
    end


    fig = Figure(size=(1920, 1080))
    colors = [:red, :cyan, :yellow]
    ax = Axis3(fig[1, 1], xticklabelsvisible=false, 
                          yticklabelsvisible=false,
                          zticklabelsvisible=false)
    ax2 = Axis3(fig[1, 2], xticklabelsvisible=false, 
                          yticklabelsvisible=false,
                          zticklabelsvisible=false)
    hidespines!(ax)
    hidespines!(ax2)
    xlims!(ax, -a[2].val, a[2].val)
    ylims!(ax, -a[2].val, a[2].val)
    zlims!(ax, -a[2].val, a[2].val)

    xlims!(ax2, -a[1].val/3, a[1].val/3)
    ylims!(ax2, -a[1].val/3, a[1].val/3)
    zlims!(ax2, -a[1].val/3, a[1].val/3)

    extra_frames = (N÷100)*10
    frames = 2:(N + extra_frames)

    nt = N÷100*20
    r1s = CircularBuffer{Point3f}(nt)
    r2s = CircularBuffer{Point3f}(nt)
    r3s = CircularBuffer{Point3f}(nt)

    fill!(r1s, positions[:, 1, 1])
    fill!(r2s, positions[:, 2, 1])
    fill!(r3s, positions[:, 3, 1])

    r1s = Observable(r1s)
    r2s = Observable(r2s)
    r3s = Observable(r3s)

    r1s_in = CircularBuffer{Point3f}(nt)
    r2s_in = CircularBuffer{Point3f}(nt)

    fill!(r1s_in, positions[:, 1, 1] .- com_12[:,1])
    fill!(r2s_in, positions[:, 2, 1] .- com_12[:,1])

    r1s_in = Observable(r1s_in)
    r2s_in = Observable(r2s_in)

    colors = Makie.wong_colors()[[1, 2, 3]]
    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]

    sc1s = Observable{Point3f}()
    sc2s = Observable{Point3f}()
    # sc3s = Observable{Point3f}()
    
    lines!(ax, r1s, color=trailcolors[1], linewidth=2.5)
    lines!(ax, r2s, color=trailcolors[2], linewidth=2.5)
    lines!(ax, r3s, color=trailcolors[3], linewidth=2.5)

    lines!(ax2, r1s_in, color=trailcolors[1], linewidth=2.5)
    lines!(ax2, r2s_in, color=trailcolors[2], linewidth=2.5) 
    
    mesh!(ax, Sphere(sc1s[], Rs[1]), color=Makie.wong_colors()[1], colorrange=(1, 3))
    mesh!(ax, Sphere(sc2s[], Rs[2]), color=Makie.wong_colors()[2], colorrange=(1, 3))
    # mesh!(ax, Sphere(sc3s[], Rs[3]), color=Makie.wong_colors()[3], colorrange=(1, 3))


    p = Progress(N)
    record(fig, "figures/triple_animation_collide_mesh.mp4", frames; framerate = N÷25) do frame

        if frame <= N
            push!(r1s[], positions[:, 1, frame])
            push!(r2s[], positions[:, 2, frame])
            push!(r3s[], positions[:, 3, frame])
            

            r1_12 = positions[:, 1, frame] .- com_12[:,frame]
            r2_12 = positions[:, 2, frame] .- com_12[:,frame]

            # d12 = norm(r1_12 .- r2_12)
            push!(r1s_in[], r1_12)
            push!(r2s_in[], r2_12)

            sc1s[] = positions[:, 1, frame] |> Point3f
            sc2s[] = positions[:, 2, frame] |> Point3f
            # sc3s[] = positions[:, 3, frame] |> Point3f

            notify(r1s)
            notify(r2s)
            notify(r3s)

            notify(r1s_in)
            notify(r2s_in)

            notify(sc1s)
            notify(sc2s)
        # notify(sc3s)
        else
            # println((frame - N)/((N÷100)*10))
            ax2.azimuth[] = 2π*sin((N - frame)/((N÷100)*10))
        end
        next!(p)
    end

# 	row = df[1222,:]
# 	# row = df[1622,:]
# 	colliders = row.retcodes[:Collision][1]

# 	r1 = getproperty(row, Symbol("r$(colliders[1])_final")) .|> u"Rsun" .|> ustrip
# 	r2 = getproperty(row, Symbol("r$(colliders[2])_final")) .|> u"Rsun" .|> ustrip

# 	v1 = getproperty(row, Symbol("v$(colliders[1])_final")) .|> u"km/cs" .|> ustrip
# 	v2 = getproperty(row, Symbol("v$(colliders[2])_final")) .|> u"km/cs" .|> ustrip
	
# 	r1 = r1 .- r2
# 	r2 = zeros(3)

# 	R1 = getproperty(row, Symbol("R$(colliders[1])_final")) |> u"Rsun" |> ustrip
# 	R2 = getproperty(row, Symbol("R$(colliders[2])_final")) |> u"Rsun" |> ustrip
# 	Rmax = max(R1, R2)

# 	# sl_x = Slider(fig[2, 1], range = 0:0.01:6, startvalue = 0.1)
# 	# azimuth = lift(sl_x.value) do az
#  #    	az
# 	# end
# 	azimuth=0.4π
# 	ax = Axis3(fig[1, 1], aspect=:equal, azimuth=azimuth)

# 	xlims!(ax, -Rmax*2, Rmax*2)
# 	ylims!(ax, -Rmax*2, Rmax*2)
# 	zlims!(ax, -Rmax*2, Rmax*2)

# 	mesh!(ax, Sphere(Point3f(r1), R1))
# 	mesh!(ax, Sphere(Point3f(r2), R2))

# 	arrows!(ax, [[r] for r in r1]..., [[v] for v in v1]...)
# 	arrows!(ax, [[r] for r in r2]..., [[v] for v in v2]...)

# 	fig

end

function simple_chaotic_threebody_vis(T=1; outname="simple_threebody_scatter")
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 5000
    e = 0.1
    i = [0.0, 0.0]u"rad"
    a = [0.1, 0.1]u"AU"
    positions = zeros(3, 3, N)


    triple = multibodysystem(masses, a=a, e=e, i=i)
    res = simulate(triple, t_sim=T, npoints=N, callbacks=[])
    sol = to_solution(res)
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    
    if showplot
        figg = Figure()
        axx = Axis(figg[1, 1], aspect=1)
        lines!(axx, positions[1,1,:],  positions[2,1,:])
        lines!(axx, positions[1,2,:],  positions[2,2,:])
        lines!(axx, positions[1,3,:],  positions[2,3,:])
        return figg
    end

    fig = Figure(size=(1080, 1080))
    colors = [:red, :cyan, :yellow]
    ax = Axis(fig[1, 1], xticklabelsvisible=false, 
                         yticklabelsvisible=false)
    hidedecorations!(ax)
    hidespines!(ax)
    xlims!(ax, -a[2].val, a[2].val)
    ylims!(ax, -a[2].val, a[2].val)

    frames = 2:N

    colors = Makie.wong_colors()[[1, 2, 3]]

    sc1s = Observable{Point2f}()
    sc2s = Observable{Point2f}()
    sc3s = Observable{Point2f}()
    
    scatter!(ax, sc1s, color=colors[1], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))
    scatter!(ax, sc2s, color=colors[2], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))
    scatter!(ax, sc3s, color=colors[3], marker=:star5, colormap=Makie.wong_colors(), colorrange=(1, 3))


    savepath = joinpath(FIGPATH, outname)*".mp4"
    record(fig, savepath, frames; framerate = N÷25) do frame


        sc1s[] = positions[1:2, 1, frame] |> Point2f
        sc2s[] = positions[1:2, 2, frame] |> Point2f
        sc3s[] = positions[1:2, 3, frame] |> Point2f

        notify(sc1s)
        notify(sc2s)
        notify(sc3s)

    end

end

function chaotic_3body(;outname="chaotic_3body", showplot=false)
    masses = [1.0, 1.0, 1.0]u"Msun"
    N = 3000
    e = [0.1, 0.5]
    i = [0.0, 0.0]u"rad"
    a = [0.8, 1.1]u"AU"
    positions = zeros(3, 3, N)

    triple = multibodysystem(masses, a=a, e=e, i=i)
    res = simulate(triple, t_sim=15.0, npoints=N, callbacks=[])
    sol = to_solution(res)
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    
    if showplot
        figg = Figure()
        axx = Axis(figg[1, 1], aspect=1)
        lines!(axx, positions[1,1,:],  positions[2,1,:])
        lines!(axx, positions[1,2,:],  positions[2,2,:])
        lines!(axx, positions[1,3,:],  positions[2,3,:])
        return figg
    end

    fig = Figure(size=(1080, 1080))
    colors = [:red, :cyan, :yellow]
    ax = Axis(fig[1, 1], xticklabelsvisible=false, 
                         yticklabelsvisible=false,
                         xticksvisible=false,
                         yticksvisible=false,
                         xgridvisible=false, ygridvisible=false, backgroundcolor=:transparent)
    
    hidespines!(ax)
    xlims!(ax, -2.5, 2.5)
    ylims!(ax, -2.5, 2.5)

    frames = 2:N
    colors = Makie.wong_colors()[[1, 2, 3]]


    sc1s = Observable{Point2f}()
    sc2s = Observable{Point2f}()
    sc3s = Observable{Point2f}()
    

    scatter!(ax, sc1s, color=colors[1], colormap=Makie.wong_colors(), colorrange=(1, 3), markersize=100, markerspace=:pixel)
    scatter!(ax, sc2s, color=colors[2], colormap=Makie.wong_colors(), colorrange=(1, 3), markersize=100, markerspace=:pixel)
    scatter!(ax, sc3s, color=colors[3], colormap=Makie.wong_colors(), colorrange=(1, 3), markersize=100, markerspace=:pixel)

    savepath = joinpath(FIGPATH, outname)*".gif"

    p = Progress(N)
    record(fig, savepath, frames; framerate = N÷25) do frame



        sc1s[] = positions[1:2, 1, frame] |> Point2f
        sc2s[] = positions[1:2, 2, frame] |> Point2f
        sc3s[] = positions[1:2, 3, frame] |> Point2f

        notify(sc1s)
        notify(sc2s)
        notify(sc3s)

            # axs[i].azimuth[] = 2π*sin(frame/N)
        next!(p)
    end

end

function many_binaries_2d(;outname="many_binaries_2d")
    N = 2000
    es = [0.1, 0.3, 0.6, 0.9]
    m1 = 1.0u"Msun"
    a = 1.0u"AU"
    qs = [1.0, 0.5, 0.1]
    params = [[es[i], qs[j]] for j in eachindex(qs), i in eachindex(es)]
    positions = zeros(size(params)..., 3, 2, N)

    for i in axes(params, 1)
        for j in axes(params, 2)
            e, q = params[i, j]
            masses = [m1, m1*q]
            triple = multibodysystem(masses, a=a, e=e)
            res = simulate(triple, t_sim=2, npoints=N, callbacks=[])
            sol = to_solution(res)
            positions[i, j, :, :, :] = ustrip.(u"AU", sol.r)
        end
    end

    fig = Figure(size=(1920, 1080))
    axs = [Axis(fig[i, j], aspect=1) for i in axes(params, 1), j in axes(params, 2)]
    # axs = Axis(fig[1, 1])
    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)
    for ax in axs
        hidedecorations!(ax)
        hidespines!(ax)
        xlims!(ax, -a.val*1.1, a.val*1.1)
        ylims!(ax, -a.val*1.1, a.val*1.1)
    end

    frames = 2:N

    sc1s = [Observable{Point2f}() for i in axes(params, 1), j in axes(params, 2)]
    sc2s = [Observable{Point2f}() for i in axes(params, 1), j in axes(params, 2)]

    
    nt = N÷100*40
    r1s = [CircularBuffer{Point2f}(nt) for i in axes(params, 1), j in axes(params, 2)]
    r2s = [CircularBuffer{Point2f}(nt) for i in axes(params, 1), j in axes(params, 2)]

    for i in CartesianIndices(axs)
        fill!(r1s[i], positions[i, 1:2, 1, 1])
        fill!(r2s[i], positions[i, 1:2, 2, 1])
    end

    r1s = Observable.(r1s)
    r2s = Observable.(r2s)

    colors = Makie.wong_colors()[[2, 3]]
    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]


    for i in CartesianIndices(axs)
        lines!(axs[i], r1s[i], color=trailcolors[1], linewidth=4)
        lines!(axs[i], r2s[i], color=trailcolors[2], linewidth=4)
        scatter!(axs[i], sc1s[i], color=colors[1], marker=:star5, markersize=20,
                 colormap=Makie.wong_colors(), colorrange=(1, 3))
        scatter!(axs[i], sc2s[i], color=colors[2], marker=:star5, markersize=20,
                 colormap=Makie.wong_colors(), colorrange=(1, 3))
    end

    p = Progress(N)
    savepath = joinpath(FIGPATH, outname)*".mp4"
    record(fig, savepath, frames; framerate = N÷18) do frame

        for i in CartesianIndices(axs)
            
            push!(r1s[i][], positions[i, 1:2, 1, frame])
            push!(r2s[i][], positions[i, 1:2, 2, frame])

            sc1s[i][] = positions[i, 1:2, 1, frame] |> Point2f
            sc2s[i][] = positions[i, 1:2, 2, frame] |> Point2f

            
            notify(r1s[i])
            notify(r2s[i])


            notify(sc1s[i])
            notify(sc2s[i])
        end

        next!(p)
    end
end

function real_unstable_triples(T=20)

    N = 10_000

    files = readdir(joinpath(@__DIR__, "..", "data"), join=true)

    ejections = filter(x -> occursin("ejection", x), files)
    collisions = filter(x -> occursin("collision", x), files)

    triple_ejection = JLD2.load(ejections[1], "triple")
    triple_collision = JLD2.load(collisions[1], "triple")

    animate_1_triple(triple_ejection, T, N, "real_triple_ejection_white")
    animate_1_triple(triple_collision, 14, N, "real_triple_collision_white")

end

function animate_1_triple(triple, T, N, outname)
    a_out = triple.binaries[2].elements.a |> u"AU" |> ustrip

    radii = [p.structure.R for (k, p) in sort(triple.particles)]

    triple = multibodysystem(triple.particles.mass, a=t.binaries.a, e=t.binaries.e, 
                             ν=t.binaries.ν, ω=t.binaries.ω, Ω=t.binaries.Ω, 
                             R=radii)

    res = simulate(triple, t_sim=T, npoints=N, callbacks=[EscapeCB(), CollisionCB()])
    sol = to_solution(res)
    positions = zeros(3, 3, length(sol.t))
    positions[:, :, :] = ustrip.(u"AU", sol.r)
    # return res.retcode, length(sol.t)
    N = length(sol.t)

    fig = Figure(size=(1080, 1080))
    colors = [:red, :cyan, :yellow]
    # colors = Makie.wong_colors()
    ax = Axis3(fig[1, 1], xticklabelsvisible=false, 
                                 yticklabelsvisible=false, 
                                 zticklabelsvisible=false)
    # ax.elevation = 0.25π
    hidespines!(ax)
    xlims!(ax, -a_out*1.1, a_out*1.1)
    ylims!(ax, -a_out*1.1, a_out*1.1)
    zlims!(ax, -a_out*1.1, a_out*1.1)

    frames = 1:N

    # r1s = Observable(Point3f[])
    # r2s = Observable(Point3f[])
    # r3s = Observable(Point3f[])

    nt = N÷100*10
    r1s = CircularBuffer{Point3f}(nt)
    r2s = CircularBuffer{Point3f}(nt)
    r3s = CircularBuffer{Point3f}(nt)

    fill!(r1s, positions[:, 1, 1])
    fill!(r2s, positions[:, 2, 1])
    fill!(r3s, positions[:, 3, 1])

    r1s = Observable(r1s)
    r2s = Observable(r2s)
    r3s = Observable(r3s)
    colors = Makie.wong_colors()[[1, 2, 3]]
    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]


    lines!(ax, r1s, color=trailcolors[1], linewidth=3)
    lines!(ax, r2s, color=trailcolors[2], linewidth=3)
    lines!(ax, r3s, color=trailcolors[3], linewidth=3)

    p = Progress(N)
    savepath = joinpath(FIGPATH, outname)*".mp4"
    record(fig, savepath, frames; framerate = N÷30) do frame


        push!(r1s[], positions[:, 1, frame])
        push!(r2s[], positions[:, 2, frame])
        push!(r3s[], positions[:, 3, frame])
        notify(r1s)
        notify(r2s)
        notify(r3s)

        next!(p)
        ax.azimuth[] = 2π*sin(frame/N) - π

    end
end

function chaotic_triple_varying_nu(T=20; outname="varying_nu")

    N = 3500#10_000

    files = readdir(joinpath(@__DIR__, "..", "data"), join=true)

    collisions = filter(x -> occursin("14838", x), files)
    collisions = filter(x -> !occursin("collision", x), collisions)
    n = length(collisions)

    positions = zeros(n, 3, 3, N)

    for i ∈ 1:n
        triple = JLD2.load(collisions[i], "triple")
        radii = [p.structure.R for (k, p) in sort(triple.particles)]

        triple = multibodysystem(triple.particles.mass, a=t.binaries.a, e=t.binaries.e, 
                                 ν=t.binaries.ν, ω=t.binaries.ω, Ω=t.binaries.Ω, 
                                 R=radii)
        res = simulate(triple, t_sim=T, npoints=N, callbacks=[])
        sol = to_solution(res)
        positions[i, :, :, :] = ustrip.(u"AU", sol.r)
    end

    fig = Figure(size=(1920, 1080))
    colors = [:red, :cyan, :yellow]
    ax = Axis3(fig[1, 1], xticklabelsvisible=false, 
                          yticklabelsvisible=false, 
                          zticklabelsvisible=false)

    hidespines!(ax)
    xlims!(ax, -10, 10)
    ylims!(ax, -10, 10)
    zlims!(ax, -10, 10)

    frames = 2:N

    nt = N÷100*10
    r1s = [CircularBuffer{Point3f}(nt) for i = 1:n]
    r2s = [CircularBuffer{Point3f}(nt) for i = 1:n]
    r3s = [CircularBuffer{Point3f}(nt) for i = 1:n]

    for i = 1:n
        fill!(r1s[i], positions[i, :, 1, 1])
        fill!(r2s[i], positions[i, :, 2, 1])
        fill!(r3s[i], positions[i, :, 3, 1])
    end

    r1s = Observable.(r1s)
    r2s = Observable.(r2s)
    r3s = Observable.(r3s)

    colors = Makie.ColorSchemes.viridis[1:end÷(n-1):end]#Makie.wong_colors()
    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]

    sc1s = [Observable{Point3f}() for i = 1:n]
    sc2s = [Observable{Point3f}() for i = 1:n]
    sc3s = [Observable{Point3f}() for i = 1:n]



    for i = 1:n
        lines!(ax, r1s[i], color=trailcolors[i], linewidth=3, linestyle=:solid)
        lines!(ax, r2s[i], color=trailcolors[i], linewidth=3, linestyle=:dash)
        lines!(ax, r3s[i], color=trailcolors[i], linewidth=3, linestyle=:dashdot)
        scatter!(ax, sc1s[i], color=colors[i], marker=:star5, colormap=:viridis, colorrange=(1, 256))
        scatter!(ax, sc2s[i], color=colors[i], marker=:star6, colormap=:viridis, colorrange=(1, 256))
        scatter!(ax, sc3s[i], color=colors[i], marker=:star8, colormap=:viridis, colorrange=(1, 256))

    end

    p = Progress(N)
    savepath = joinpath(FIGPATH, outname)*".mp4"
    record(fig, savepath, frames; framerate = N÷40) do frame

        for i = 1:n
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

        end
        ax.azimuth[] = 2π*sin(2frame/N) - π
        # ax.elevation[] = 0.1π*sin(2frame/N)
        next!(p)
    end
end


function making_logo()

    jl_red   = Makie.RGBf(([203, 60 ,61] ./ 255)...)
    jl_green = Makie.RGBf(([56, 152 ,38] ./ 255)...)
    jl_purp  = Makie.RGBf(([149, 88 ,178] ./ 255)...)


    fig = Figure(size=(1080, 1080))#, backgroundcolor=:transparent)
    # colors = [jl_red, jl_green, jl_purp]
    colors = Makie.wong_colors()[2:5]
    ax = Axis(fig[1, 1], xticklabelsvisible=false,
                         yticklabelsvisible=false,
                         xticksvisible=false,
                         yticksvisible=false,
                         xgridvisible=false, ygridvisible=false)#, backgroundcolor=:transparent)
    
    hidespines!(ax)

    t0 = 1925.0u"d"
    triple = multibodysystem([2.0, 1.0, 1.5]u"Msun", time=t0,
                             a=[1.0, 8.0]u"AU", e=[0.0, 0.0])
    # p = triple.binaries[2].elements.P

    # p1 = [-0.7, -0.7, 0.0]u"AU"
    # p2 = [0.7, -0.7, 0.0]u"AU"
    # p3 = [0.0, 1.0, 0.0]u"AU"

    # ps = [p1, p2, p3]
    # vs = [zeros(3)u"km/s" for i = 1:3]

    # particles = [Syzygy.Particle(p.key, p.parent, p.sibling, p.mass, pos, vel, p.structure, nothing) for (p, pos, vel) in zip(values(sort(triple.particles)), ps, vs)]
    # particles = Dict(i => p for (i, p) in enumerate(particles))
    # triple = Syzygy.MultiBodySystem(triple.n, triple.time, particles, triple.pairs,
    #                                 triple.binaries, triple.levels, triple.root, triple.hierarchy, nothing)
    sol = simulate(triple, t_sim=1.3, npoints=5_000) |> to_solution

    N = length(sol.t)
    positions = zeros(3, 3, N)

    positions[:, :, :] = ustrip.(u"AU", sol.r)

    nt = N#÷100*20

    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]

    lines!(ax, positions[1, 1, :], positions[2, 1, :], color=trailcolors[1], linewidth=5)
    lines!(ax, positions[1, 2, :], positions[2, 2, :], color=trailcolors[2], linewidth=5)
    lines!(ax, positions[1, 3, :], positions[2, 3, :], color=trailcolors[3], linewidth=5)

    scatter!(ax, positions[1, 1, end], positions[2, 1, end], color=colors[1], markersize=100)
    scatter!(ax, positions[1, 2, end], positions[2, 2, end], color=colors[2], markersize=100)
    scatter!(ax, positions[1, 3, end], positions[2, 3, end], color=colors[3], markersize=100)

    save(joinpath(@__DIR__, "..", "grouplogo.png"), fig)
    # fig
end


function simple_hierarchical_triple_figure()


    fig = Figure(size=(1080, 700), backgroundcolor=Makie.RGBf(([41,42,53] ./ 255)...) )
    # colors = [jl_red, jl_green, jl_purp]
    colors = Makie.wong_colors()[[3, 2, 4]]
    ax = Axis(fig[1, 1], xticklabelsvisible=false,
                         yticklabelsvisible=false,
                         xticksvisible=false,
                         yticksvisible=false,
                         xgridvisible=false, ygridvisible=false, backgroundcolor=Makie.RGBf(([41,42,53] ./ 255)...))
    
    hidespines!(ax)

    triple = multibodysystem([2.0, 1.0, 1.5]u"Msun",
                             a=[0.5, 8.0]u"AU", e=[0.1, 0.5])

    sol = simulate(triple, t_sim=1.2) |> to_solution

    N = length(sol.t)
    positions = zeros(3, 3, N)

    positions[:, :, :] = ustrip.(u"AU", sol.r)

    nt = N#÷100*20

    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]

    lines!(ax, positions[1, 1, :], positions[2, 1, :], color=trailcolors[1], linewidth=5)
    lines!(ax, positions[1, 2, :], positions[2, 2, :], color=trailcolors[2], linewidth=5)
    lines!(ax, positions[1, 3, :], positions[2, 3, :], color=trailcolors[3], linewidth=5)

    scatter!(ax, positions[1, 1, end], positions[2, 1, end], color=colors[1], markersize=30)
    scatter!(ax, positions[1, 2, end], positions[2, 2, end], color=colors[2], markersize=30)
    scatter!(ax, positions[1, 3, end], positions[2, 3, end], color=colors[3], markersize=30)

    save(joinpath(@__DIR__, "..", "simple-ht.png"), fig)
    fig
end

function simple_unstable_triple_figure()


    fig = Figure(size=(1080, 1080), backgroundcolor=Makie.RGBf(([41,42,53] ./ 255)...) )
    # colors = [jl_red, jl_green, jl_purp]
    colors = Makie.wong_colors()[[3, 2, 4]]
    ax = Axis(fig[1, 1], xticklabelsvisible=false,
                         yticklabelsvisible=false,
                         xticksvisible=false,
                         yticksvisible=false,
                         xgridvisible=false, ygridvisible=false, backgroundcolor=Makie.RGBf(([41,42,53] ./ 255)...))
    
    hidespines!(ax)

    triple = multibodysystem([2.0, 1.0, 1.5]u"Msun",
                             a=[0.35, 1.0]u"AU", e=[0.1, 0.2])

    sol = simulate(triple, t_sim=5.8) |> to_solution

    N = length(sol.t)
    positions = zeros(3, 3, N)

    positions[:, :, :] = ustrip.(u"AU", sol.r)

    Nn = N÷100*80
    nt = Nn

    trailcolors = [[RGBAf(c.r, c.g, c.b, (i/nt)^2.5) for i in 1:nt] for c in colors]

    lines!(ax, positions[1, 1, end-Nn+1:end], positions[2, 1, end-Nn+1:end], color=trailcolors[1], linewidth=5)
    lines!(ax, positions[1, 2, end-Nn+1:end], positions[2, 2, end-Nn+1:end], color=trailcolors[2], linewidth=5)
    lines!(ax, positions[1, 3, end-Nn+1:end], positions[2, 3, end-Nn+1:end], color=trailcolors[3], linewidth=5)

    scatter!(ax, positions[1, 1, end], positions[2, 1, end], color=colors[1], markersize=30)
    scatter!(ax, positions[1, 2, end], positions[2, 2, end], color=colors[2], markersize=30)
    scatter!(ax, positions[1, 3, end], positions[2, 3, end], color=colors[3], markersize=30)

    save(joinpath(@__DIR__, "..", "simple-ut.png"), fig)
    fig
end