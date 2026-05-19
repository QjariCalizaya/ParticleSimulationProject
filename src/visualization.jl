using GLMakie

function plot_trajectories_3d(trajectories, output_path)
    fig = Figure(size = (900, 700))

    ax = Axis3(
        fig[1, 1],
        xlabel = "x",
        ylabel = "y",
        zlabel = "z",
        title = "Simulación 3D de partículas cargadas con gravedad"
    )

    for i in eachindex(trajectories)
        xs, ys, zs = trajectories[i]

        lines!(ax, xs, ys, zs, linewidth = 2)
        scatter!(ax, [xs[end]], [ys[end]], [zs[end]], markersize = 12)
    end

    save(output_path, fig)
    display(fig)
end