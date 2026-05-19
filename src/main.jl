include("particles.jl")
include("physics.jl")
include("integrators.jl")
include("visualization.jl")

function create_particles()
    return [
        Particle3D(1.0,  1e-6, [-1.0, 0.0, 2.0], [2.0, 0.5, 4.0], [0.0, 0.0, 0.0]),
        Particle3D(1.0, -1e-6, [ 1.0, 0.0, 2.5], [-2.0, -0.3, 3.0], [0.0, 0.0, 0.0]),
        Particle3D(1.0,  1e-6, [ 0.0, 1.0, 3.0], [0.5, -1.5, 2.0], [0.0, 0.0, 0.0])
    ]
end

function simulate(method_name)
    dt = 0.001
    total_time = 2.0
    steps = Int(total_time / dt)

    particles = create_particles()

    trajectories = [
        (Float64[], Float64[], Float64[]) for _ in particles
    ]

    reset_forces!(particles)
    apply_coulomb_forces!(particles)

    for step in 1:steps
        for i in eachindex(particles)
            push!(trajectories[i][1], particles[i].position[1])
            push!(trajectories[i][2], particles[i].position[2])
            push!(trajectories[i][3], particles[i].position[3])
        end

        if method_name == "euler"
            reset_forces!(particles)
            apply_coulomb_forces!(particles)
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
end

simulate("euler")
simulate("verlet")