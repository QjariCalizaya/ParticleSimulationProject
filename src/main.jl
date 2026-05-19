include("particles.jl")
include("physics.jl")
include("integrators.jl")
include("visualization.jl")

function simulate()
    dt = 0.001
    total_time = 2.0
    steps = Int(total_time / dt)

    particles = [
        Particle3D(1.0,  1e-6, [-1.0, 0.0, 2.0], [2.0, 0.5, 4.0], [0.0, 0.0, 0.0]),
        Particle3D(1.0, -1e-6, [ 1.0, 0.0, 2.5], [-2.0, -0.3, 3.0], [0.0, 0.0, 0.0]),
        Particle3D(1.0,  1e-6, [ 0.0, 1.0, 3.0], [0.5, -1.5, 2.0], [0.0, 0.0, 0.0])
    ]

    trajectories = [
        (Float64[], Float64[], Float64[]) for _ in particles
    ]

    for step in 1:steps
        for i in eachindex(particles)
            push!(trajectories[i][1], particles[i].position[1])
            push!(trajectories[i][2], particles[i].position[2])
            push!(trajectories[i][3], particles[i].position[3])
        end

        reset_forces!(particles)
        apply_coulomb_forces!(particles)
        euler_step!(particles, dt)
    end

    plot_trajectories_3d(
        trajectories,
        "plots/particles_3d.png"
    )
end

simulate()