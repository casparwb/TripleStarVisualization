### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ c17bfbb0-520b-11ee-1ad1-89753b6ee0a8
begin
	using Pkg
	Pkg.add(url="https://github.com/casparwb/Syzygy.jl")
	Pkg.add(["WGLMakie", "Unitful", "UnitfulAstro"])
end

# ╔═╡ a5ca170c-63f6-4659-aaf7-c37abc138457
begin
	using Syzygy, WGLMakie, Unitful, UnitfulAstro
	WGLMakie.set_theme!(theme_black())
end

# ╔═╡ a787c632-f59d-4c54-8cfd-85120997188d
html"""<style>
main {
    max-width: 800px;
}
"""

# ╔═╡ d78a0863-1113-4e85-814c-b99b1e3fa991
html"<button onclick=present()>Present</button<"

# ╔═╡ de693fc5-c4ba-44d1-82bb-d2ee3d9ef78d
md"""
# Test
"""

# ╔═╡ afb48fef-4f34-4e3b-96f9-739e6640fa67
md"""
# Figs
"""

# ╔═╡ 2ff5852c-6c0c-455e-ab52-a284837c5687
let
	masses = [1.0, 1.0, 1.0]u"Msun"

	es = [[0.1, 0.1], [0.1, 0.5], [0.1, 0.9]]
	smas = [[0.1, 2.0]u"AU", [0.3, 2.0]u"AU", [0.6, 2.0]u"AU", [0.8, 2.0]u"AU"]
	params = [[es[i], smas[j]] for i in eachindex(es), j in eachindex(smas)]
	# params = hcat(es, smas)
	fig = Figure(resolution=(1280, 720))
	colors = [:red, :cyan, :yellow]
	for i in axes(params, 1)
		for j in axes(params, 2)
			triple = multibodysystem(masses, a=params[i, j][2], e=params[i, j][1])
			res = simulate(triple, t_sim=1, npoints=1000, callbacks=[])
			sol = analyse_simulation(res)

			a = ustrip.(u"AU", params[i, j][2])
			e = params[i, j][1]
			ax = Axis(fig[i, j], aspect=1)#, title="a = $a, e = $e")
			hidedecorations!(ax)
			hidespines!(ax)
			for k = 1:3
				lines!(ax, ustrip.(u"AU", sol.r[1, k, :]), 
						  ustrip.(u"AU", sol.r[2,k,:]), color=colors[k],
						  linewidth=1)
			end
		end
	end
	# save("triples_with_varying_e_and_a.png", fig)
	fig
end

# ╔═╡ 036a1c27-6712-4661-a599-04347e01350b
md"""
## Animation
"""

# ╔═╡ f477f9fd-6f72-4890-9e2c-c220612f6f80
let
	masses = [1.0, 1.0, 1.0]u"Msun"
	N = 1000
	es = [[0.1, 0.1], [0.1, 0.5], [0.1, 0.9]]
	smas = [[0.1, 2.0]u"AU", [0.3, 2.0]u"AU", [0.6, 2.0]u"AU", [0.8, 2.0]u"AU"]
	params = [[es[i], smas[j]] for i in eachindex(es), j in eachindex(smas)]
	positions = zeros(size(params)..., 3, 3, N)

	for i in axes(params, 1)
		for j in axes(params, 2)
			triple = multibodysystem(masses, a=params[i, j][2], e=params[i, j][1])
			res = simulate(triple, t_sim=1, npoints=N, callbacks=[])
			sol = analyse_simulation(res)
			positions[i, j, :, :, :] = ustrip.(u"AU", sol.r)
		end
	end

	fig = Figure(resolution=(1280, 720))
	colors = [:red, :cyan, :yellow]
	axs = [Axis(fig[i, j], aspect=1) for i in axes(params, 1), j in axes(params, 2)]
	
	for ax in axs
		hidedecorations!(ax)
		hidespines!(ax)
	end

	frames = 1:N
	# prim = @lift((positions[1, 1, 1,1,1:$frame], positions[1, 1, 2,1,1:$frame]))
	# seco = @lift((positions[1, 1, 1,2,1:$frame], positions[1, 1, 2,2,1:$frame]))
	# tert = @lift((positions[1, 1, 1,3,1:$frame], positions[1, 1, 2,3,1:$frame]))

	r1 = Observable(Point2f[])
	r2 = Observable(Point2f[])
	r3 = Observable(Point2f[])

	lines!(axs[1, 1], r1, color=colors[1], linewidth=1)
	lines!(axs[1, 1], r2, color=colors[1], linewidth=1)
	lines!(axs[1, 1], r3, color=colors[1], linewidth=1)
	# lines!(axs[1, 1], seco[1], seco[2], color=colors[2], linewidth=1)
	# lines!(axs[1, 1], tert[1], tert[2], color=colors[3], linewidth=1)

	record(fig, "triple_animation.mp4", frames; framerate = N÷10) do frame
		# for i in axes(params, 1)
		# 	for j in axes(params, 2)
		# 			r = positions[i, j, :, :, :]
		# 			for  k = 1:3
						
		# 			end
		# 		end
		# 	end
		push!(r1[], positions[1, 1, 1:2, 1, frame])
		push!(r2[], positions[1, 1, 1:2, 2, frame])
		push!(r3[], positions[1, 1, 1:2, 3, frame])
	end;
	# save("triples_with_varying_e_and_a.png", fig)
	# fig
end;

# ╔═╡ 429f83a0-c5d5-45dd-8fbe-d4cafda6a759
[[0.0, 0.0]u"rad" [0.0, π/6]u"rad"; [0.0, π/2]u"rad" [0.0, π]u"rad"]

# ╔═╡ Cell order:
# ╠═a787c632-f59d-4c54-8cfd-85120997188d
# ╠═d78a0863-1113-4e85-814c-b99b1e3fa991
# ╠═de693fc5-c4ba-44d1-82bb-d2ee3d9ef78d
# ╠═c17bfbb0-520b-11ee-1ad1-89753b6ee0a8
# ╠═a5ca170c-63f6-4659-aaf7-c37abc138457
# ╟─afb48fef-4f34-4e3b-96f9-739e6640fa67
# ╠═2ff5852c-6c0c-455e-ab52-a284837c5687
# ╠═036a1c27-6712-4661-a599-04347e01350b
# ╠═f477f9fd-6f72-4890-9e2c-c220612f6f80
# ╠═91b96b00-ab6e-41e7-9dff-551f95826041
# ╠═429f83a0-c5d5-45dd-8fbe-d4cafda6a759
