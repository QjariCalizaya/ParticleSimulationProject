function plot_energy_error(method_name, times, energies)
    E0 = energies[1]

    relative_errors = [
        abs(E - E0) / abs(E0) for E in energies
    ]

    fig = Figure(size = (900, 500))

    ax = Axis(
        fig[1, 1],
        xlabel = "Tiempo (s)",
        ylabel = "Error relativo de energía",
        title = "Conservación de energía usando método $(method_name)"
    )

    lines!(ax, times, relative_errors, linewidth = 2)

    save("plots/energy_error_$(lowercase(method_name)).png", fig)
    display(fig)
end