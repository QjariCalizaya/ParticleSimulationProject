using GLMakie

include("particles.jl")
include("physics.jl")
include("integrators.jl")
include("energyError.jl")

mkpath("plots")

const N_PARTICLES = 3

# Simulación

function simulate_particles(particles, method_name, dt, total_time)
    steps = max(1, floor(Int, total_time / dt))

    trajectories = [
        (Float64[], Float64[], Float64[]) for _ in particles
    ]

    times = Float64[]
    energies = Float64[]

    reset_forces!(particles)
    apply_pairwise_forces!(particles)

    for step in 1:steps
        t = (step - 1) * dt

        push!(times, t)
        push!(energies, total_energy(particles))

        for i in eachindex(particles)
            push!(trajectories[i][1], particles[i].position[1])
            push!(trajectories[i][2], particles[i].position[2])
            push!(trajectories[i][3], particles[i].position[3])
        end

        if method_name == "Эйлер" || method_name == "Euler"
            reset_forces!(particles)
            apply_pairwise_forces!(particles)
            euler_step!(particles, dt)

        elseif method_name == "Верле" || method_name == "Verlet"
            verlet_step!(particles, dt)

        else
            error("Неизвестный метод: $method_name")
        end
    end

    return trajectories, times, energies
end

# Interfaz principal

fig = Figure(size = (1500, 900))

main_grid = fig[1, 1] = GridLayout()

plot_panel = main_grid[1, 1] = GridLayout()
control_panel = main_grid[1, 2] = GridLayout()
time_panel = main_grid[2, 1:2] = GridLayout()

colsize!(main_grid, 1, Relative(0.72))
colsize!(main_grid, 2, Relative(0.28))
rowsize!(main_grid, 1, Relative(0.88))
rowsize!(main_grid, 2, Relative(0.12))

ax = Axis3(
    plot_panel[1, 1],
    xlabel = "x",
    ylabel = "y",
    zlabel = "z",
    title = "Трёхмерное моделирование частиц",
)

# Trayectorias completas
xs_path = [Observable(Float64[]) for _ in 1:N_PARTICLES]
ys_path = [Observable(Float64[]) for _ in 1:N_PARTICLES]
zs_path = [Observable(Float64[]) for _ in 1:N_PARTICLES]

# Posición actual según slider de tiempo
xs_now = [Observable(Float64[]) for _ in 1:N_PARTICLES]
ys_now = [Observable(Float64[]) for _ in 1:N_PARTICLES]
zs_now = [Observable(Float64[]) for _ in 1:N_PARTICLES]

for i in 1:N_PARTICLES
    lines!(ax, xs_path[i], ys_path[i], zs_path[i], linewidth = 2)
    scatter!(ax, xs_now[i], ys_now[i], zs_now[i], markersize = 16)
end

# Panel de controles
Label(control_panel[1, 1:3], "Управление симуляцией", fontsize = 20)

preset_options = [
    "Гравитационная орбита",
    "Неустойчивый метод Эйлера",
    "Электрическое притяжение",
    "Электрическое отталкивание",
    "Смешанная система"
]

method_menu = Menu(
    control_panel[2, 1:3],
    options = ["Верле", "Эйлер"],
    default = "Верле"
)

dt_slider = Slider(
    control_panel[3, 1:2],
    range = 0.001:0.001:0.05,
    startvalue = 0.01
)

Label(
    control_panel[3, 3],
    lift(x -> "dt = $(round(x, digits = 4))", dt_slider.value)
)

duration_slider = Slider(
    control_panel[4, 1:2],
    range = 5.0:5.0:10000.0,
    startvalue = 20.0
)

Label(
    control_panel[4, 3],
    lift(x -> "Длительность = $(round(x, digits = 1)) с", duration_slider.value)
)

run_button = Button(
    control_panel[5, 1:3],
    label = "Запустить симуляцию"
)

energy_button = Button(
    control_panel[6, 1:3],
    label = "Сохранить график ошибки"
)

compare_button = Button(
    control_panel[7, 1:3],
    label = "Сравнить Эйлер и Верле"
)

Label(control_panel[8, 1:3], "Частицы", fontsize = 18)

# Sliders de partículas

# Masa
m1 = Slider(control_panel[9, 1], range = 100.0:100.0:3000.0, startvalue = 1000.0)
m2 = Slider(control_panel[9, 2], range = 0.1:0.1:20.0, startvalue = 1.0)
m3 = Slider(control_panel[9, 3], range = 0.1:0.1:20.0, startvalue = 0.5)

Label(control_panel[10, 1], lift(x -> "m1=$(round(x, digits=2))", m1.value))
Label(control_panel[10, 2], lift(x -> "m2=$(round(x, digits=2))", m2.value))
Label(control_panel[10, 3], lift(x -> "m3=$(round(x, digits=2))", m3.value))

# Carga eléctrica
q1 = Slider(control_panel[11, 1], range = -0.002:0.0001:0.002, startvalue = 0.0)
q2 = Slider(control_panel[11, 2], range = -0.002:0.0001:0.002, startvalue = 0.0)
q3 = Slider(control_panel[11, 3], range = -0.002:0.0001:0.002, startvalue = 0.0)

Label(control_panel[12, 1], lift(q -> "q1=$(round(q, sigdigits=3))", q1.value))
Label(control_panel[12, 2], lift(q -> "q2=$(round(q, sigdigits=3))", q2.value))
Label(control_panel[12, 3], lift(q -> "q3=$(round(q, sigdigits=3))", q3.value))

# Posición X
x1 = Slider(control_panel[13, 1], range = -30.0:0.5:30.0, startvalue = 0.0)
x2 = Slider(control_panel[13, 2], range = -30.0:0.5:30.0, startvalue = 10.0)
x3 = Slider(control_panel[13, 3], range = -30.0:0.5:30.0, startvalue = 16.0)

Label(control_panel[14, 1], lift(x -> "x1=$(round(x,digits=1))", x1.value))
Label(control_panel[14, 2], lift(x -> "x2=$(round(x,digits=1))", x2.value))
Label(control_panel[14, 3], lift(x -> "x3=$(round(x,digits=1))", x3.value))

# Posición Y
y1 = Slider(control_panel[15, 1], range = -30.0:0.5:30.0, startvalue = 0.0)
y2 = Slider(control_panel[15, 2], range = -30.0:0.5:30.0, startvalue = 0.0)
y3 = Slider(control_panel[15, 3], range = -30.0:0.5:30.0, startvalue = 0.0)

Label(control_panel[16, 1], lift(y -> "y1=$(round(y,digits=1))", y1.value))
Label(control_panel[16, 2], lift(y -> "y2=$(round(y,digits=1))", y2.value))
Label(control_panel[16, 3], lift(y -> "y3=$(round(y,digits=1))", y3.value))

# Posición Z
z1 = Slider(control_panel[17, 1], range = -30.0:0.5:30.0, startvalue = 0.0)
z2 = Slider(control_panel[17, 2], range = -30.0:0.5:30.0, startvalue = 0.0)
z3 = Slider(control_panel[17, 3], range = -30.0:0.5:30.0, startvalue = 3.0)

Label(control_panel[18, 1], lift(z -> "z1=$(round(z,digits=1))", z1.value))
Label(control_panel[18, 2], lift(z -> "z2=$(round(z,digits=1))", z2.value))
Label(control_panel[18, 3], lift(z -> "z3=$(round(z,digits=1))", z3.value))

# Velocidad X
vx1 = Slider(control_panel[19, 1], range = -15.0:0.5:15.0, startvalue = 0.0)
vx2 = Slider(control_panel[19, 2], range = -15.0:0.5:15.0, startvalue = 0.0)
vx3 = Slider(control_panel[19, 3], range = -15.0:0.5:15.0, startvalue = 0.0)

Label(control_panel[20, 1], lift(v -> "vx1=$(round(v,digits=1))", vx1.value))
Label(control_panel[20, 2], lift(v -> "vx2=$(round(v,digits=1))", vx2.value))
Label(control_panel[20, 3], lift(v -> "vx3=$(round(v,digits=1))", vx3.value))

# Velocidad Y
vy1 = Slider(control_panel[21, 1], range = -15.0:0.5:15.0, startvalue = 0.0)
vy2 = Slider(control_panel[21, 2], range = -15.0:0.5:15.0, startvalue = 9.0)
vy3 = Slider(control_panel[21, 3], range = -15.0:0.5:15.0, startvalue = 7.0)

Label(control_panel[22, 1], lift(v -> "vy1=$(round(v,digits=1))", vy1.value))
Label(control_panel[22, 2], lift(v -> "vy2=$(round(v,digits=1))", vy2.value))
Label(control_panel[22, 3], lift(v -> "vy3=$(round(v,digits=1))", vy3.value))

# Velocidad Z
vz1 = Slider(control_panel[23, 1], range = -15.0:0.5:15.0, startvalue = 0.0)
vz2 = Slider(control_panel[23, 2], range = -15.0:0.5:15.0, startvalue = 0.0)
vz3 = Slider(control_panel[23, 3], range = -15.0:0.5:15.0, startvalue = 1.0)

Label(control_panel[24, 1], lift(v -> "vz1=$(round(v,digits=1))", vz1.value))
Label(control_panel[24, 2], lift(v -> "vz2=$(round(v,digits=1))", vz2.value))
Label(control_panel[24, 3], lift(v -> "vz3=$(round(v,digits=1))", vz3.value))


# Slider de tiempo real
Label(time_panel[1, 1], "Время визуализации", fontsize = 18)

time_slider = Slider(
    time_panel[1, 2],
    range = 0.0:0.001:1.0,
    startvalue = 0.0
)

time_label = Label(time_panel[1, 3], "t = 0.0 с")

preset_menu = Menu(
    time_panel[1, 4],
    options = preset_options,
    default = "Гравитационная орбита"
)

preset_button = Button(
    time_panel[1, 5],
    label = "Применить сценарий"
)

colsize!(time_panel, 4, Relative(0.15))
colsize!(time_panel, 5, Relative(0.15))

colsize!(time_panel, 1, Relative(0.18))
colsize!(time_panel, 2, Relative(0.65))
colsize!(time_panel, 3, Relative(0.17))


# Estado interno
last_trajectories = Ref([
    (Float64[], Float64[], Float64[]) for _ in 1:N_PARTICLES
])

last_times = Ref(Float64[])
last_energies = Ref(Float64[])
last_method = Ref("Верле")

function create_particles_from_sliders()
    return [
        Particle3D(
            m1.value[],
            q1.value[],
            [x1.value[], y1.value[], z1.value[]],
            [vx1.value[], vy1.value[], vz1.value[]],
            [0.0, 0.0, 0.0]
        ),

        Particle3D(
            m2.value[],
            q2.value[],
            [x2.value[], y2.value[], z2.value[]],
            [vx2.value[], vy2.value[], vz2.value[]],
            [0.0, 0.0, 0.0]
        ),

        Particle3D(
            m3.value[],
            q3.value[],
            [x3.value[], y3.value[], z3.value[]],
            [vx3.value[], vy3.value[], vz3.value[]],
            [0.0, 0.0, 0.0]
        )
    ]
end

function set_slider!(slider, value)
    set_close_to!(slider, value)
    notify(slider.value)
end

function apply_preset!(preset_name)
    if preset_name == "Гравитационная орбита"
        method_menu.selection[] = "Верле"

        set_slider!(dt_slider, 0.01)
        set_slider!(duration_slider, 40.0)

        set_slider!(m1, 1000.0)
        set_slider!(m2, 1.0)
        set_slider!(m3, 0.5)

        set_slider!(q1, 0.0)
        set_slider!(q2, 0.0)
        set_slider!(q3, 0.0)

        set_slider!(x1, 0.0)
        set_slider!(y1, 0.0)
        set_slider!(z1, 0.0)

        set_slider!(x2, 10.0)
        set_slider!(y2, 0.0)
        set_slider!(z2, 0.0)

        set_slider!(x3, 16.0)
        set_slider!(y3, 0.0)
        set_slider!(z3, 3.0)

        set_slider!(vx1, 0.0)
        set_slider!(vy1, 0.0)
        set_slider!(vz1, 0.0)

        set_slider!(vx2, 0.0)
        set_slider!(vy2, 9.0)
        set_slider!(vz2, 0.0)

        set_slider!(vx3, 0.0)
        set_slider!(vy3, 7.0)
        set_slider!(vz3, 1.0)

    elseif preset_name == "Неустойчивый метод Эйлера"
        method_menu.selection[] = "Эйлер"

        set_slider!(dt_slider, 0.03)
        set_slider!(duration_slider, 60.0)

        set_slider!(m1, 1000.0)
        set_slider!(m2, 1.0)
        set_slider!(m3, 0.5)

        set_slider!(q1, 0.0)
        set_slider!(q2, 0.0)
        set_slider!(q3, 0.0)

        set_slider!(x1, 0.0)
        set_slider!(y1, 0.0)
        set_slider!(z1, 0.0)

        set_slider!(x2, 10.0)
        set_slider!(y2, 0.0)
        set_slider!(z2, 0.0)

        set_slider!(x3, 16.0)
        set_slider!(y3, 0.0)
        set_slider!(z3, 3.0)

        set_slider!(vx1, 0.0)
        set_slider!(vy1, 0.0)
        set_slider!(vz1, 0.0)

        set_slider!(vx2, 0.0)
        set_slider!(vy2, 9.0)
        set_slider!(vz2, 0.0)

        set_slider!(vx3, 0.0)
        set_slider!(vy3, 7.0)
        set_slider!(vz3, 1.0)

    elseif preset_name == "Электрическое притяжение"
        method_menu.selection[] = "Верле"

        set_slider!(dt_slider, 0.005)
        set_slider!(duration_slider, 100.0)

        set_slider!(m1, 10.0)
        set_slider!(m2, 10.0)
        set_slider!(m3, 10.0)

        set_slider!(q1, 0.001)
        set_slider!(q2, -0.001)
        set_slider!(q3, -0.0002)

        set_slider!(x1, -10.0)
        set_slider!(y1, 0.0)
        set_slider!(z1, 0.0)

        set_slider!(x2, 10.0)
        set_slider!(y2, 0.0)
        set_slider!(z2, 0.0)

        set_slider!(x3, 0.0)
        set_slider!(y3, 10.0)
        set_slider!(z3, 0.0)

        set_slider!(vx1, 0.0)
        set_slider!(vy1, 2.5)
        set_slider!(vz1, 0.0)

        set_slider!(vx2, 0.0)
        set_slider!(vy2, -5.0)
        set_slider!(vz2, 0.0)

        set_slider!(vx3, 0.0)
        set_slider!(vy3, 0.0)
        set_slider!(vz3, 4.0)

    elseif preset_name == "Электрическое отталкивание"
        method_menu.selection[] = "Верле"

        set_slider!(dt_slider, 0.005)
        set_slider!(duration_slider, 15.0)

        set_slider!(m1, 10.0)
        set_slider!(m2, 10.0)
        set_slider!(m3, 10.0)

        set_slider!(q1, 0.001)
        set_slider!(q2, 0.001)
        set_slider!(q3, 0.001)

        set_slider!(x1, -5.0)
        set_slider!(y1, 0.0)
        set_slider!(z1, 0.0)

        set_slider!(x2, 5.0)
        set_slider!(y2, 0.0)
        set_slider!(z2, 0.0)

        set_slider!(x3, 0.0)
        set_slider!(y3, 8.0)
        set_slider!(z3, 0.0)

        set_slider!(vx1, 0.0)
        set_slider!(vy1, 0.0)
        set_slider!(vz1, 0.0)

        set_slider!(vx2, 0.0)
        set_slider!(vy2, 0.0)
        set_slider!(vz2, 0.0)

        set_slider!(vx3, 0.0)
        set_slider!(vy3, 0.0)
        set_slider!(vz3, 0.0)

    elseif preset_name == "Смешанная система"
        method_menu.selection[] = "Верле"

        set_slider!(dt_slider, 0.005)
        set_slider!(duration_slider, 30.0)

        set_slider!(m1, 1000.0)
        set_slider!(m2, 5.0)
        set_slider!(m3, 5.0)

        set_slider!(q1, 0.0005)
        set_slider!(q2, -0.0003)
        set_slider!(q3, 0.0003)

        set_slider!(x1, 0.0)
        set_slider!(y1, 0.0)
        set_slider!(z1, 0.0)

        set_slider!(x2, 12.0)
        set_slider!(y2, 0.0)
        set_slider!(z2, 0.0)

        set_slider!(x3, -14.0)
        set_slider!(y3, 0.0)
        set_slider!(z3, 4.0)

        set_slider!(vx1, 0.0)
        set_slider!(vy1, 0.0)
        set_slider!(vz1, 0.0)

        set_slider!(vx2, 0.0)
        set_slider!(vy2, 8.0)
        set_slider!(vz2, 0.0)

        set_slider!(vx3, 0.0)
        set_slider!(vy3, -7.0)
        set_slider!(vz3, 1.0)
    end

    time_slider.value[] = 0.0
    update_simulation!()
end

function update_time_view!()
    times = last_times[]
    trajectories = last_trajectories[]

    if isempty(times)
        return
    end

    alpha = time_slider.value[]
    index = clamp(round(Int, 1 + alpha * (length(times) - 1)), 1, length(times))

    current_time = times[index]
    time_label.text[] = "t = $(round(current_time, digits = 2)) с"

    for i in 1:N_PARTICLES
        xs, ys, zs = trajectories[i]

        xs_now[i][] = [xs[index]]
        ys_now[i][] = [ys[index]]
        zs_now[i][] = [zs[index]]
    end
end

function update_simulation!()
    particles = create_particles_from_sliders()

    method_name = method_menu.selection[]
    dt = dt_slider.value[]
    total_time = duration_slider.value[]

    trajectories, times, energies = simulate_particles(
        particles,
        method_name,
        dt,
        total_time
    )

    for i in 1:N_PARTICLES
        xs, ys, zs = trajectories[i]

        xs_path[i][] = xs
        ys_path[i][] = ys
        zs_path[i][] = zs
    end

    last_trajectories[] = trajectories
    last_times[] = times
    last_energies[] = energies
    last_method[] = method_name

    update_time_view!()
end

on(run_button.clicks) do _
    update_simulation!()
end

on(energy_button.clicks) do _
    if !isempty(last_times[])
        plot_energy_error(last_method[], last_times[], last_energies[])
    end
end

on(time_slider.value) do _
    update_time_view!()
end

on(compare_button.clicks) do _
    dt = dt_slider.value[]
    total_time = duration_slider.value[]

    particles_euler = create_particles_from_sliders()
    particles_verlet = create_particles_from_sliders()

    _, times_euler, energies_euler = simulate_particles(
        particles_euler,
        "Эйлер",
        dt,
        total_time
    )

    _, times_verlet, energies_verlet = simulate_particles(
        particles_verlet,
        "Верле",
        dt,
        total_time
    )

    plot_energy_comparison(
        times_euler,
        energies_euler,
        times_verlet,
        energies_verlet
    )
end

function get_selected_preset()
    selected = preset_menu.selection[]

    if selected isa String
        return selected
    elseif selected isa Int
        return preset_options[selected]
    else
        try
            return preset_options[preset_menu.i_selected[]]
        catch
            return string(selected)
        end
    end
end

function selected_preset_name()
    selected = preset_menu.selection[]

    if selected isa String
        return selected
    elseif selected isa Int
        return preset_options[selected]
    else
        return string(selected)
    end
end

on(preset_menu.selection) do _
    preset_name = selected_preset_name()
    println("Применение сценария из меню: ", preset_name)
    apply_preset!(preset_name)
end


# Primera simulación automática
update_simulation!()

screen = display(fig)
wait(screen)

