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