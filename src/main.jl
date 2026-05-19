using GLMakie
using LinearAlgebra

const GROUND_Z = 0.0
const g = 9.81
const k = 8.9875517923e9

mutable struct Particle3D
    mass::Float64
    charge::Float64
    position::Vector{Float64}
    velocity::Vector{Float64}
    force::Vector{Float64}
end

function reset_forces!(particles)
    for p in particles
        p.force = [0.0, 0.0, -p.mass * g]
    end
end

function apply_coulomb_forces!(particles)
    n = length(particles)

    for i in 1:n-1
        for j in i+1:n
            pi = particles[i]
            pj = particles[j]

            r_vec = pj.position - pi.position
            distance = norm(r_vec)

            if distance < 1e-6
                continue
            end

            direction = r_vec / distance
            force_magnitude = k * pi.charge * pj.charge / distance^2
            force = force_magnitude * direction

            pi.force += force
            pj.force -= force
        end
    end
end

function euler_step!(particles, dt)
    for p in particles
        acceleration = p.force / p.mass
        p.velocity += acceleration * dt
        p.position += p.velocity * dt

        if p.position[3] < GROUND_Z
            p.position[3] = GROUND_Z
            p.velocity[3] = -0.6 * p.velocity[3]
        end
    end
end

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

    save("plots/particles_3d.png", fig)
    display(fig)
end

simulate()