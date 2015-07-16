using Base.Test
using UnitQuaternions

####################################################################################################
# Constructors
####################################################################################################

print("Testing constructors...")

q = UnitQuaternion([1, 1, 1], 1)
for i = 1:4
    @test_approx_eq q[i] 0.5
end

@test_throws ErrorException UnitQuaternion([1, 1, 1, 1], 1)

@test UnitQuaternion() == UnitQuaternion([0, 0, 0], 1)

for i = 1:100
    a, b, c, d = 5randn(4)
    p = UnitQuaternion(a, b, c, d)
    q = UnitQuaternion([a, b, c, d])
    @test p == UnitQuaternion([a, b, c], d)
    @test q == UnitQuaternion([a, b, c], d)
    @test p == q
    @test_approx_eq sqrt(q.η^2 + norm(q.ϵ)^2) 1
end

@test UnitQuaternion([0, 0, 0]) == UnitQuaternion()
@test UnitQuaternion([π, 0, 0]) == UnitQuaternion(1, 0, 0, 0)
@test UnitQuaternion([2π, 0, 0]) == UnitQuaternion(0, 0, 0, -1)
@test UnitQuaternion([0, π/2, 0]) == UnitQuaternion(0, sqrt(2)/2, 0, sqrt(2)/2)
@test UnitQuaternion([0, 0, π/3]) == UnitQuaternion(0, 0, 0.5, sqrt(3)/2)

println(" done.")

####################################################################################################
# Operators
####################################################################################################

print("Testing operators...")

@test_approx_eq +(UnitQuaternion()) eye(4)
@test_approx_eq +(UnitQuaternion(1, 0, 0, 0)) [0 0 0 1; 0 0 1 0; 0 -1 0 0; -1 0 0 0]
@test_approx_eq +(UnitQuaternion(1, 1, 1, 1)) [ 0.5  0.5 -0.5  0.5; 
                                               -0.5  0.5  0.5  0.5;
                                                0.5 -0.5  0.5  0.5;
                                               -0.5 -0.5 -0.5  0.5;]

@test UnitQuaternion(1, 2, 3, 4) + UnitQuaternion() == UnitQuaternion(1, 2, 3, 4)
@test UnitQuaternion() + UnitQuaternion(1, 2, 3, 4) == UnitQuaternion(1, 2, 3, 4)
@test UnitQuaternion(1, 2, 3, 4) ⊕ UnitQuaternion() == UnitQuaternion(1, 2, 3, 4)
@test UnitQuaternion() ⊕ UnitQuaternion(1, 2, 3, 4) == UnitQuaternion(1, 2, 3, 4)

p = UnitQuaternion(sqrt(2)/2, 0, 0, sqrt(2)/2)
q = UnitQuaternion(0, 0, sqrt(2)/2, sqrt(2)/2)
@test p + q == UnitQuaternion(0.5, 0.5, 0.5, 0.5)
@test q ⊕ p == UnitQuaternion(0.5, 0.5, 0.5, 0.5)
@test q + p == UnitQuaternion(0.5, -0.5, 0.5, 0.5)
@test p ⊕ q == UnitQuaternion(0.5, -0.5, 0.5, 0.5)

r = UnitQuaternion(4, 3, -1, -3)
@test -r == UnitQuaternion(-4, -3, 1, 3)
@test r == -r

@test q ⊞ [0, 0, 0] == q
@test p ⊞ (q ⊟ p) == q
@test_approx_eq (q ⊞ [0.1, 0.2, 0.3]) ⊟ q [0.1, 0.2, 0.3]

r = inv(p) + q
@test r == UnitQuaternion(-0.5, -0.5, 0.5, 0.5)
@test_approx_eq q ⊟ p 2atan2(norm(r.ϵ), r.η) * r.ϵ / norm(r.ϵ)

@test p ⊞ [0, 0, π/2] == UnitQuaternion(0.5, 0.5, 0.5, 0.5)

@test_throws BoundsError q[5]
@test_approx_eq q[1] 0
@test_approx_eq q[4] sqrt(2)/2

println(" done.")

####################################################################################################
# Methods
####################################################################################################

print("Testing methods...")

q = UnitQuaternion(0, 0, sqrt(2)/2, sqrt(2)/2)

@test_approx_eq angle(q) π/2

@test_approx_eq axis(q) [0, 0, 1]

qmean = UnitQuaternion(sind(90/2), 0, 0, cosd(90/2))
v = [UnitQuaternion(sind(60/2), 0, 0, cosd(60/2)), UnitQuaternion(sind(120/2), 0, 0, cosd(120/2))]
cov = covariance(v, qmean, [1, 1])
@test_approx_eq rad2deg(sqrt(cov[1])) 30

@test inv(q) == UnitQuaternion(0, 0, -sqrt(2)/2, sqrt(2)/2)

qs = [UnitQuaternion(), UnitQuaternion()]
weights = [rand(), rand()]
@test mean(qs, weights) == UnitQuaternion()
@test mean(qs) == UnitQuaternion()
qs = [UnitQuaternion(1, 0, 0, 0), UnitQuaternion(1, 0, 0, 0)]
weights = [rand(), rand()]
@test mean(qs, weights) == UnitQuaternion(1, 0, 0, 0)
@test mean(qs) == UnitQuaternion(1, 0, 0, 0)
qs = [UnitQuaternion(1.1, 0, 0, 0), UnitQuaternion(0.9, 0, 0, 0)]
weights = [1, 1]
@test mean(qs, weights) == UnitQuaternion(1, 0, 0, 0)
@test mean(qs) == UnitQuaternion(1, 0, 0, 0)
qs = [UnitQuaternion(0, 1.2, 0, 0), UnitQuaternion(0, 0.8, 0, 0)]
weights = [1, 1]
@test mean(qs, weights) == UnitQuaternion(0, 1, 0, 0)
@test mean(qs) == UnitQuaternion(0, 1, 0, 0)
qs = [UnitQuaternion(0, 0, 1.3, 0), UnitQuaternion(0, 0, 0.7, 0)]
weights = [1, 1]
@test mean(qs, weights) == UnitQuaternion(0, 0, 1, 0)
@test mean(qs) == UnitQuaternion(0, 0, 1, 0)
qs = [UnitQuaternion(0, 0, 0, 1.4), UnitQuaternion(0, 0, 0, 0.6)]
weights = [1, 1]
@test mean(qs, weights) == UnitQuaternion(0, 0, 0, 1)
@test mean(qs) == UnitQuaternion(0, 0, 0, 1)

@test_approx_eq log(UnitQuaternion()) [0, 0, 0]
@test_approx_eq log(UnitQuaternion(1, 1, 1, 0)) π/2 * [sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]
@test_approx_eq log(q) [0, 0, π/4]

@test_approx_eq rotatevector(q, [1, 0, 0]) [0, 1, 0]
p = UnitQuaternion(0.5, 0.5, 0.5, 0.5)
@test_approx_eq rotatevector(p, [1, 0, 0]) [0, 1, 0]

@test_approx_eq vector(UnitQuaternion()) [0, 0, 0, 1]
@test_approx_eq vector(q) [0, 0, sqrt(2)/2, sqrt(2)/2]

println(" done.")
