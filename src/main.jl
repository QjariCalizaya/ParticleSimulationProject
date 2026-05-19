using GLMakie

x = 0:0.1:10
y = sin.(x)

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, x, y)

display(fig)

println("Presiona Enter para cerrar...")
readline()