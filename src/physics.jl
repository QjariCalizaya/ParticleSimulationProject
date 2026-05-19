using LinearAlgebra

const GROUND_Z = 0.0
const g = 9.81
const k = 8.9875517923e9

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