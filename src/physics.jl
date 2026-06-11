using LinearAlgebra

const G = 6.67430e-11
const k = 8.9875517923e9
const MIN_DISTANCE = 1e-6

function reset_forces!(particles)
    for p in particles
        p.force = [0.0, 0.0, 0.0]
    end
end

function apply_pairwise_forces!(particles)
    n = length(particles)

    for i in 1:n-1
        for j in i+1:n
            pi = particles[i]
            pj = particles[j]

            r_vec = pj.position - pi.position
            distance = norm(r_vec)

            if distance < MIN_DISTANCE
                continue
            end

            direction = r_vec / distance

            coulomb_magnitude = k * pi.charge * pj.charge / distance^2
            gravity_magnitude = G * pi.mass * pj.mass / distance^2

            coulomb_force = -coulomb_magnitude * direction
            gravity_force = gravity_magnitude * direction

            total_force = coulomb_force + gravity_force

            pi.force += total_force
            pj.force -= total_force
        end
    end
end

function kinetic_energy(particles)
    return sum(0.5 * p.mass * dot(p.velocity, p.velocity) for p in particles)
end



function potential_energy_coulomb(particles)
    n = length(particles)
    energy = 0.0

    for i in 1:n-1
        for j in i+1:n
            pi = particles[i]
            pj = particles[j]

            r = norm(pj.position - pi.position)

            if r < MIN_DISTANCE
                continue
            end

            energy += k * pi.charge * pj.charge / r
        end
    end

    return energy
end

function potential_energy_gravity(particles)
    n = length(particles)
    energy = 0.0

    for i in 1:n-1
        for j in i+1:n
            pi = particles[i]
            pj = particles[j]

            r = norm(pj.position - pi.position)

            if r < MIN_DISTANCE
                continue
            end

            energy += -G * pi.mass * pj.mass / r
        end
    end

    return energy
end

function total_energy(particles)
    return kinetic_energy(particles) +
           potential_energy_coulomb(particles) +
           potential_energy_gravity(particles)
end
