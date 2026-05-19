using GLMakie

include("Particles.jl")

function apply_gravity!(p::Particle, g::Float64)
    p.force = [0.0, -p.mass * g]
end

function euler_step!(p::Particle, dt::Float64)
    acceleration = p.force / p.mass
    p.velocity = p.velocity + acceleration * dt
    p.position = p.position + p.velocity * dt
end

function simulate()
    dt = 0.01
    total_time = 5.0
    steps = Int(total_time / dt)

    particle = Particle(
        1.0,
        [0.0, 0.0],
        [5.0, 10.0],
        [0.0, 0.0]
    )

    g = 9.81

    xs = Float64[]
    ys = Float64[]

    for i in 1:steps
        push!(xs, particle.position[1])
        push!(ys, particle.position[2])

        apply_gravity!(particle, g)
        euler_step!(particle, dt)
    end

    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "x",
        ylabel = "y",
        title = "Movimiento de una partícula con gravedad"
    )

    lines!(ax, xs, ys, linewidth = 2)

    save("plots/trajectory.png", fig)
    display(fig)
end

simulate()