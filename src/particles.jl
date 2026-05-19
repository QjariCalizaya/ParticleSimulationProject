mutable struct Particle3D
    mass::Float64
    charge::Float64
    position::Vector{Float64}
    velocity::Vector{Float64}
    force::Vector{Float64}
end