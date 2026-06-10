function translate_method_name(method_name)
    if method_name == "Euler" || method_name == "Эйлер"
        return "Эйлера"
    elseif method_name == "Verlet" || method_name == "Верле"
        return "Верле"
    else
        return method_name
    end
end

function normalize_method_filename(method_name)
    if method_name == "Euler" || method_name == "Эйлер"
        return "euler"
    elseif method_name == "Verlet" || method_name == "Верле"
        return "verlet"
    else
        return lowercase(replace(method_name, r"\s+" => "_"))
    end
end

function plot_energy_error(method_name, times, energies)
    E0 = energies[1]
    method_ru = translate_method_name(method_name)

    relative_errors = [
        abs(E - E0) / max(abs(E0), 1e-12) for E in energies
    ]

    fig = Figure(size = (900, 500))

    ax = Axis(
        fig[1, 1],
        xlabel = "Время (с)",
        ylabel = "Относительная ошибка энергии",
        title = "Сохранение энергии методом $(method_ru)"
    )

    lines!(ax, times, relative_errors, linewidth = 2)

    save("plots/energy_error_$(normalize_method_filename(method_name)).png", fig)
    display(fig)
end

function plot_energy_comparison(times_euler, energies_euler, times_verlet, energies_verlet)
    E0_euler = energies_euler[1]
    E0_verlet = energies_verlet[1]

    error_euler = [
        abs(E - E0_euler) / max(abs(E0_euler), 1e-12) for E in energies_euler
    ]

    error_verlet = [
        abs(E - E0_verlet) / max(abs(E0_verlet), 1e-12) for E in energies_verlet
    ]

    fig = Figure(size = (1000, 550))

    ax = Axis(
        fig[1, 1],
        xlabel = "Время (с)",
        ylabel = "Относительная ошибка энергии",
        title = "Сравнение сохранения энергии: метод Эйлера и метод Верле"
    )

    lines!(ax, times_euler, error_euler, linewidth = 2, label = "Метод Эйлера")
    lines!(ax, times_verlet, error_verlet, linewidth = 2, label = "Метод Верле")

    axislegend(ax)

    save("plots/energy_comparison_euler_vs_verlet.png", fig)
    display(fig)
end
