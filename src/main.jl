include("particles.jl")
include("physics.jl")
include("integrators.jl")
include("visualization.jl")

using GLMakie

function create_particles()
    return [
        Particle3D(1000000000000000,  0, [-1.0, 0.0, 2.0], [2.0, 0.5, 4.0], [0.0, 0.0, 0.0]),
        Particle3D(1, 0, [ 1.0, 0.0, 2.5], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        Particle3D(1, 0, [ 0.0, 1.0, 3.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    ]
end

function plot_energy(method_name, times, energies)
    fig = Figure(size = (900, 500))

    ax = Axis(
        fig[1, 1],
        xlabel = "Tiempo (s)",
        ylabel = "Energía total",
        title = "Energía total usando método $(method_name)"
    )

    lines!(ax, times, energies, linewidth = 2)

    save("plots/energy_$(method_name).png", fig)
    display(fig)
end

function simulate(method_name)
    dt = 0.001
    total_time = 10.0
    steps = Int(total_time / dt)

    particles = create_particles()

    trajectories = [
        (Float64[], Float64[], Float64[]) for _ in particles
    ]

    times = Float64[]
    energies = Float64[]

    reset_forces!(particles)
    apply_pairwise_forces!(particles)

    for step in 1:steps
        current_time = step * dt

        push!(times, current_time)
        push!(energies, total_energy(particles))

        for i in eachindex(particles)
            push!(trajectories[i][1], particles[i].position[1])
            push!(trajectories[i][2], particles[i].position[2])
            push!(trajectories[i][3], particles[i].position[3])
        end

        if method_name == "euler"
            reset_forces!(particles)
            apply_pairwise_forces!(particles)
            euler_step!(particles, dt)

        elseif method_name == "verlet"
            verlet_step!(particles, dt)

        else
            error("Método desconocido: $method_name")
        end
    end

    plot_trajectories_3d(
        trajectories,
        "plots/$(method_name)_3d.png"
    )

    plot_energy(method_name, times, energies)
end

simulate("euler")
simulate("verlet")