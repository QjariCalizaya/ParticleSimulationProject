function euler_step!(particles, dt)
    for p in particles
        acceleration = p.force / p.mass

        p.velocity += acceleration * dt
        p.position += p.velocity * dt

        handle_ground_collision!(p)
    end
end

function verlet_step!(particles, dt)
    old_accelerations = [p.force / p.mass for p in particles]

    for i in eachindex(particles)
        p = particles[i]
        a = old_accelerations[i]

        p.position += p.velocity * dt + 0.5 * a * dt^2
        handle_ground_collision!(p)
    end

    reset_forces!(particles)
    apply_coulomb_forces!(particles)

    new_accelerations = [p.force / p.mass for p in particles]

    for i in eachindex(particles)
        p = particles[i]
        a_old = old_accelerations[i]
        a_new = new_accelerations[i]

        p.velocity += 0.5 * (a_old + a_new) * dt
        handle_ground_collision!(p)
    end
end

function handle_ground_collision!(p)
    if p.position[3] < GROUND_Z
        p.position[3] = GROUND_Z
        p.velocity[3] = -0.6 * p.velocity[3]
    end
end