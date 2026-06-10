function euler_step!(particles, dt)
    for p in particles
        acceleration = p.force / p.mass

        p.velocity += acceleration * dt
        p.position += p.velocity * dt
    end
end

function verlet_step!(particles, dt)
    old_accelerations = [p.force / p.mass for p in particles]

    for i in eachindex(particles)
        p = particles[i]
        a = old_accelerations[i]

        p.position += p.velocity * dt + 0.5 * a * dt^2
    end

    reset_forces!(particles)
    apply_pairwise_forces!(particles)

    new_accelerations = [p.force / p.mass for p in particles]

    for i in eachindex(particles)
        p = particles[i]

        p.velocity += 0.5 * (old_accelerations[i] + new_accelerations[i]) * dt
    end
end
