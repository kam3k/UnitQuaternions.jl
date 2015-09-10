module UnitQuaternions

import Base: +, -, ^, ==, ×, angle, inv, log, mean, show, getindex

export UnitQuaternion, ⊕, ⊞, ⊟, axis, covariance, rotatevector, rotateframe, rotationmatrix, slerp, vector

const EPS = 1e-9

####################################################################################################
# Type Definition
####################################################################################################

immutable UnitQuaternion
    ϵ::Vector{Float64}
    η::Float64

    function UnitQuaternion(ϵ, η) 
        length(ϵ) == 3 || error("Must be a 3-vector.")
        η = abs(η) < EPS ? zero(η) : η # small scalars are exactly zero
        # normalization algorithm taken from http://stackoverflow.com/a/12934750/1053656
        squared_mag = η * η + sum(abs2(ϵ))
        if abs(one(squared_mag) - squared_mag) < 2.107342e-08
            scale = 2.0 / (1.0 + squared_mag)
        else
            scale = 1.0 / sqrt(squared_mag)
        end
        η >= zero(η) ? new(ϵ * scale, η * scale) : new(-ϵ * scale, -η * scale)
    end
end

####################################################################################################
# Constructors
####################################################################################################

UnitQuaternion() = UnitQuaternion([0, 0, 0], 1)

UnitQuaternion(ϵ1::Real, ϵ2::Real, ϵ3::Real, η::Real) = UnitQuaternion([ϵ1, ϵ2, ϵ3], η)

function UnitQuaternion(v::AbstractVector) 
    if length(v) == 3
        θ = norm(v)
        return θ > EPS ? UnitQuaternion(sin(θ/2.0) * v / θ, cos(θ/2.0)) : UnitQuaternion()
    elseif length(v) == 4
        return UnitQuaternion([v[1], v[2], v[3]], v[4])
    else
        error("Must be a 3-vector (rotation vector) or a 4-vector (quaternion).")
    end
end

####################################################################################################
# Operators
####################################################################################################

+(q::UnitQuaternion) = [q.η * eye(3) - ×(q.ϵ) q.ϵ; -q.ϵ' q.η]

+(p::UnitQuaternion, q::UnitQuaternion) = UnitQuaternion(+(p) * vector(q))

-(q::UnitQuaternion) = UnitQuaternion(-q.ϵ, -q.η)

^(q::UnitQuaternion, n::Integer) = UnitQuaternion(2n * log(q))
^(q::UnitQuaternion, n::Real) = UnitQuaternion(2n * log(q))

⊕(q::UnitQuaternion) = [q.η * eye(3) + ×(q.ϵ) q.ϵ; -q.ϵ' q.η]

⊕(p::UnitQuaternion, q::UnitQuaternion) = UnitQuaternion(⊕(p) * vector(q))

⊟(p::UnitQuaternion, q::UnitQuaternion) = 2.0 * log(inv(q) + p)

function ⊞(q::UnitQuaternion, rotationvector::AbstractVector) 
    length(rotationvector) == 3 || error("Must be 3-vector.")
    q + UnitQuaternion(rotationvector)
end

==(p::UnitQuaternion, q::UnitQuaternion) = abs(p.η - q.η) < EPS && all(i->(abs(i) < EPS), p.ϵ - q.ϵ)

getindex(q::UnitQuaternion, i::Integer) = i == 4 ? q.η : q.ϵ[i]

####################################################################################################
# Methods
####################################################################################################

angle(q::UnitQuaternion) = 2.0 * atan2(norm(q.ϵ), q.η)

axis(q::UnitQuaternion) = norm(q.ϵ) > EPS ? q.ϵ / norm(q.ϵ) : [0.0, 0.0, 1.0]

function covariance{T<:Real}(v::AbstractVector{UnitQuaternion}, qmean::UnitQuaternion,
                             weights::Vector{T} = ones(length(v)))
    normweights = map(w -> w/sum(weights), weights)
    cov = zeros(3, 3)
    for i = 1:length(v)
        diff = v[i] ⊟ qmean
        cov += normweights[i] * diff * diff'
    end
    return cov
end

function covariance{T<:Real}(v::AbstractVector{UnitQuaternion}, weights::Vector{T} = ones(length(v)))
    qmean = mean(v, weights)
    return covariance(v, qmean, weights)
end

inv(q::UnitQuaternion) = UnitQuaternion(-q.ϵ, q.η)

function mean{T<:Real}(v::AbstractVector{UnitQuaternion}, weights::Vector{T} = ones(length(v)))
    normweights = map(w -> w/sum(weights), weights)
    M = sum([w * vector(q) * vector(q)' for (w, q) in zip(normweights, v)])
    evals, evecs = eig(M)
    largestevalindex = sortperm(evals)[end]
    return UnitQuaternion(evecs[:,largestevalindex])
end

function rotatevector(q::UnitQuaternion, v::AbstractVector)
    # Note, the "natural order" is used here (see Trawny or Shuster), which is different
    # from the textbook by Kuipers, who uses the "historical" (Hamiltonian) order.
    length(v) == 3 || error("Must be a 3-vector.")
    (2.0 * q.η*q.η - 1.0)*v + 2.0 * dot(v, q.ϵ)*q.ϵ + 2.0 * q.η*cross(v, q.ϵ)
end

function rotateframe(q::UnitQuaternion, v::AbstractVector)
    # Note, the "natural order" is used here (see Trawny or Shuster), which is different
    # from the textbook by Kuipers, who uses the "historical" (Hamiltonian) order.
    length(v) == 3 || error("Must be a 3-vector.")
    (q.η*q.η - sum(abs2(q.ϵ)))*v + 2.0 * dot(q.ϵ, v)*q.ϵ + 2.0 * q.η*cross(q.ϵ, v)
end

"""
`rotationmatrix(q)`

Creates a rotation matrix (direction cosine matrix) of the rotation parameterized by the unit quaternion `q`.
"""
function rotationmatrix(q::UnitQuaternion)
    # Returned rotation matrix implements rotatevector (i.e., rotatevector(q, r) = rotationmatrix(q) * r)
    return (+(q) * ⊕(inv(q)))[1:3, 1:3]
end

"""
`log(q)`

Takes the log of the unit quaternion `q`.
"""
function log(q::UnitQuaternion)
    if q == UnitQuaternion()
        return [0.0, 0.0, 0.0]
    elseif abs(q.η) < EPS
        return π/2 * axis(q)
    else
        return angle(q)/2.0 * axis(q)
    end
end

"""
`show(io, q)`

Prints the unit quaternion `q` to the stream `io`.
"""
function show(io::IO, q::UnitQuaternion)
    out = @sprintf("ϵ = [%.3f, %.3f, %.3f], η = %.3f", q.ϵ[1], q.ϵ[2], q.ϵ[3], q.η)
    print(io, out)
end

"""
`slerp(p, q)`

Performs spherical linear interpolation from unit quaternion `p` to unit quaternion `q`. The
interpolation parameter `t` specifies the fraction of arc to traverse. The returned unit
quaternion parameterizes a rotation from `p` to the end point specified by `t`.
"""
slerp(p::UnitQuaternion, q::UnitQuaternion, t::Real) = (q + inv(p))^t + p

"""
`vector(q)`

Returns the fields of the unit quaternion `q` as a 4-vector.
"""
vector(q::UnitQuaternion) = [q.ϵ[1], q.ϵ[2], q.ϵ[3], q.η]

####################################################################################################
# Utility Functions
####################################################################################################

function ×(a::AbstractVector)  
    length(a) == 3 || error("Must be a 3-vector.")
    return [0.0 -a[3] a[2]; a[3] 0.0 -a[1]; -a[2] a[1] 0.0]
end

end # module
