module UnitQuaternions

import Base: +, -, angle, inv, log, show, getindex

export UnitQuaternion, ⊕, ⊞, ⊟, ι, axis, rotatevector, vector

const EPS = 1e-9

####################################################################################################
# Type Definition
####################################################################################################

immutable UnitQuaternion{T<:FloatingPoint}
    ϵ::Vector{T}
    η::T

    function UnitQuaternion(ϵ, η) 
        η = abs(η) < EPS ? 0.0 : η # small scalars are exactly zero
        # normalization algorithm taken from http://stackoverflow.com/a/12934750/1053656
        squared_mag = η * η + sum(abs2(ϵ))
        if abs(1 - squared_mag) < 2.107342e-08
            scale = 2 / (1 + squared_mag)
        else
            scale = 1 / sqrt(squared_mag)
        end
        new(ϵ * scale, η * scale)
    end
end

####################################################################################################
# Constructors
####################################################################################################

function UnitQuaternion{T<:Real}(ϵ::AbstractVector, η::T) 
    length(ϵ) == 3 || error("Must be a 3-vector.")
    UnitQuaternion{Float64}(float(ϵ), float(η))
end

const ι = UnitQuaternion([0, 0, 0], 1)

UnitQuaternion(ϵ1::Real, ϵ2::Real, ϵ3::Real, η::Real) = UnitQuaternion([ϵ1, ϵ2, ϵ3], η)

function UnitQuaternion(v::AbstractVector) 
    if length(v) == 3
        θ = norm(v)
        return θ > EPS ? UnitQuaternion(sin(θ/2) * v / θ, cos(θ/2)) : ι
    elseif length(v) == 4
        return UnitQuaternion([v[1], v[2], v[3]], v[4])
    else
        error("Must be a 3-vector (rotation vector) or a 4-vector (quaternion).")
    end
end

####################################################################################################
# Operators
####################################################################################################

+(q::UnitQuaternion) = [q.η * eye(3) - _cross(q.ϵ) q.ϵ; -q.ϵ' q.η]

+(p::UnitQuaternion, q::UnitQuaternion) = UnitQuaternion(+(p) * vector(q))

-(q::UnitQuaternion) = UnitQuaternion(-q.ϵ, -q.η)

⊕(q::UnitQuaternion) = [q.η * eye(3) + _cross(q.ϵ) q.ϵ; -q.ϵ' q.η]

⊕(p::UnitQuaternion, q::UnitQuaternion) = UnitQuaternion(⊕(p) * vector(q))

⊟(p::UnitQuaternion, q::UnitQuaternion) = 2log(inv(q) + p)

function ⊞(q::UnitQuaternion, rotationvector::AbstractArray) 
    length(rotationvector) == 3 || error("Must be 3-vector.")
    q + UnitQuaternion(rotationvector)
end

function ==(p::UnitQuaternion, q::UnitQuaternion) 
    # Because -q == q, compare p with q based on signs of their scalar parts
    sign(p.η) == sign(q.η) ? abs(p.η - q.η) < EPS && all(i->(abs(i) < EPS), p.ϵ - q.ϵ) :
                             abs(p.η + q.η) < EPS && all(i->(abs(i) < EPS), p.ϵ + q.ϵ)
end

getindex(q::UnitQuaternion, i::Integer) = i == 4 ? q.η : q.ϵ[i]

####################################################################################################
# Methods
####################################################################################################

angle(q::UnitQuaternion) = 2atan2(norm(q.ϵ), q.η)

axis(q::UnitQuaternion) = norm(q.ϵ) > EPS ? q.ϵ / norm(q.ϵ) : [0.0, 0.0, 1.0]

inv(q::UnitQuaternion) = UnitQuaternion(-q.ϵ, q.η)

function rotatevector(q::UnitQuaternion, v::AbstractVector)
    length(v) == 3 || error("Must be a 3-vector.")
    (2q.η*q.η - 1)*v + 2dot(q.ϵ, v)*q.ϵ + 2q.η*cross(q.ϵ, v)
end

function log(q::UnitQuaternion)
    if q == ι
        return [0.0, 0.0, 0.0]
    elseif abs(q.η) < EPS
        return π/2 * axis(q)
    else
        return angle(q)/2 * axis(q)
    end
end

function show(io::IO, q::UnitQuaternion)
    print(io, "ϵ = [$(q.ϵ[1]) $(q.ϵ[2]) $(q.ϵ[3])], η = $(q.η)")
end

vector(q::UnitQuaternion) = [q.ϵ[1], q.ϵ[2], q.ϵ[3], q.η]

####################################################################################################
# Utility Functions
####################################################################################################

function _cross{T<:Real}(a::Vector{T})  
    length(a) == 3 || error("Must be a 3-vector.")
    [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
end

end # module
