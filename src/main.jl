using GLMakie

include("particles.jl")
include("physics.jl")
include("integrators.jl")
include("visualization.jl")
include("energyError.jl")

mkpath("plots")


const N_PARTICLES = 3


function simulate_particles(particles, method_name, dt, total_time)
    steps = Int(total_time / dt)

    trajectories = [
        (Float64[], Float64[], Float64[]) for _ in particles
    ]

    times = Float64[]
    energies = Float64[]

    reset_forces!(particles)
    apply_pairwise_forces!(particles)

    for step in 1:steps
        t = step * dt

        push!(times, t)
        push!(energies, total_energy(particles))

        for i in eachindex(particles)
            push!(trajectories[i][1], particles[i].position[1])
            push!(trajectories[i][2], particles[i].position[2])
            push!(trajectories[i][3], particles[i].position[3])
        end

        if method_name == "Euler"
            reset_forces!(particles)
            apply_pairwise_forces!(particles)
            euler_step!(particles, dt)
        elseif method_name == "Verlet"
            verlet_step!(particles, dt)
        else
            error("Método desconocido: $method_name")
        end
    end

    return trajectories, times, energies
end



fig = Figure(size = (1400, 850))

ax = Axis3(
    fig[1:8, 1:3],
    xlabel = "x",
    ylabel = "y",
    zlabel = "z",
    title = "Simulación 3D interactiva de partículas"
)

xs_obs = [Observable(Float64[]) for _ in 1:N_PARTICLES]
ys_obs = [Observable(Float64[]) for _ in 1:N_PARTICLES]
zs_obs = [Observable(Float64[]) for _ in 1:N_PARTICLES]

x_end_obs = [Observable(Float64[]) for _ in 1:N_PARTICLES]
y_end_obs = [Observable(Float64[]) for _ in 1:N_PARTICLES]
z_end_obs = [Observable(Float64[]) for _ in 1:N_PARTICLES]

for i in 1:N_PARTICLES
    lines!(ax, xs_obs[i], ys_obs[i], zs_obs[i], linewidth = 2)
    scatter!(ax, x_end_obs[i], y_end_obs[i], z_end_obs[i], markersize = 12)
end

Label(fig[1, 4], "Parámetros generales", tellwidth = false)

method_menu = Menu(
    fig[2, 4],
    options = ["Verlet", "Euler"],
    default = "Verlet"
)

dt_slider = Slider(
    fig[3, 4],
    range = 0.001:0.001:0.05,
    startvalue = 0.01
)

time_slider = Slider(
    fig[4, 4],
    range = 1.0:1.0:100.0,
    startvalue = 20.0
)

Label(fig[3, 5], lift(x -> "dt = $(round(x, digits=4))", dt_slider.value))
Label(fig[4, 5], lift(x -> "tiempo = $(round(x, digits=2))", time_slider.value))

run_button = Button(fig[5, 4], label = "Simular")
energy_button = Button(fig[6, 4], label = "Guardar error de energía")

Label(fig[1, 6], "Partícula 1", tellwidth = false)
Label(fig[1, 7], "Partícula 2", tellwidth = false)
Label(fig[1, 8], "Partícula 3", tellwidth = false)


m1 = Slider(fig[2, 6], range = 1.0:10.0:2000.0, startvalue = 1000.0)
m2 = Slider(fig[2, 7], range = 0.1:0.1:20.0, startvalue = 1.0)
m3 = Slider(fig[2, 8], range = 0.1:0.1:20.0, startvalue = 0.5)


x1 = Slider(fig[3, 6], range = -30.0:0.5:30.0, startvalue = 0.0)
x2 = Slider(fig[3, 7], range = -30.0:0.5:30.0, startvalue = 10.0)
x3 = Slider(fig[3, 8], range = -30.0:0.5:30.0, startvalue = 16.0)

y1 = Slider(fig[4, 6], range = -30.0:0.5:30.0, startvalue = 0.0)
y2 = Slider(fig[4, 7], range = -30.0:0.5:30.0, startvalue = 0.0)
y3 = Slider(fig[4, 8], range = -30.0:0.5:30.0, startvalue = 0.0)

z1 = Slider(fig[5, 6], range = -30.0:0.5:30.0, startvalue = 0.0)
z2 = Slider(fig[5, 7], range = -30.0:0.5:30.0, startvalue = 0.0)
z3 = Slider(fig[5, 8], range = -30.0:0.5:30.0, startvalue = 3.0)

vx1 = Slider(fig[6, 6], range = -15.0:0.5:15.0, startvalue = 0.0)
vx2 = Slider(fig[6, 7], range = -15.0:0.5:15.0, startvalue = 0.0)
vx3 = Slider(fig[6, 8], range = -15.0:0.5:15.0, startvalue = 0.0)

vy1 = Slider(fig[7, 6], range = -15.0:0.5:15.0, startvalue = 0.0)
vy2 = Slider(fig[7, 7], range = -15.0:0.5:15.0, startvalue = 9.0)
vy3 = Slider(fig[7, 8], range = -15.0:0.5:15.0, startvalue = 7.0)

vz1 = Slider(fig[8, 6], range = -15.0:0.5:15.0, startvalue = 0.0)
vz2 = Slider(fig[8, 7], range = -15.0:0.5:15.0, startvalue = 0.0)
vz3 = Slider(fig[8, 8], range = -15.0:0.5:15.0, startvalue = 1.0)

Label(fig[9, 6], "m, x, y, z, vx, vy, vz")
Label(fig[9, 7], "m, x, y, z, vx, vy, vz")
Label(fig[9, 8], "m, x, y, z, vx, vy, vz")


last_times = Ref(Float64[])
last_energies = Ref(Float64[])
last_method = Ref("Verlet")

function create_particles_from_sliders()
    return [
        Particle3D(
            m1.value[],
            0.0,
            [x1.value[], y1.value[], z1.value[]],
            [vx1.value[], vy1.value[], vz1.value[]],
            [0.0, 0.0, 0.0]
        ),

        Particle3D(
            m2.value[],
            0.0,
            [x2.value[], y2.value[], z2.value[]],
            [vx2.value[], vy2.value[], vz2.value[]],
            [0.0, 0.0, 0.0]
        ),

        Particle3D(
            m3.value[],
            0.0,
            [x3.value[], y3.value[], z3.value[]],
            [vx3.value[], vy3.value[], vz3.value[]],
            [0.0, 0.0, 0.0]
        )
    ]
end

function update_simulation!()
    particles = create_particles_from_sliders()

    method_name = method_menu.selection[]
    dt = dt_slider.value[]
    total_time = time_slider.value[]

    trajectories, times, energies = simulate_particles(
        particles,
        method_name,
        dt,
        total_time
    )

    for i in 1:N_PARTICLES
        xs, ys, zs = trajectories[i]

        xs_obs[i][] = xs
        ys_obs[i][] = ys
        zs_obs[i][] = zs

        x_end_obs[i][] = [xs[end]]
        y_end_obs[i][] = [ys[end]]
        z_end_obs[i][] = [zs[end]]
    end

    last_times[] = times
    last_energies[] = energies
    last_method[] = method_name
end

on(run_button.clicks) do _
    update_simulation!()
end

on(energy_button.clicks) do _
    if !isempty(last_times[])
        plot_energy_error(last_method[], last_times[], last_energies[])
    end
end

update_simulation!()

screen = display(fig)
wait(screen)