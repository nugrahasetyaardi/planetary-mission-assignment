function [coe] =car2kep(R,V,mu)
%{
% This function computes the classical orbital elements (coe)
% from the state vector (R,V).
%
mu - gravitational parameter (km^3/s^2)
R - position vector in the geocentric equatorial frame (km)
V - velocity vector in the geocentric equatorial frame (km)
r, v - the magnitudes of R and V
vr - radial velocity component (km/s)
H - the angular momentum vector (km^2/s)
h - the magnitude of H (km^2/s)
incl - inclination of the orbit (rad)
N - the node line vector (km^2/s)
n - the magnitude of N
cp - cross product of N and R
RA - right ascension of the ascending node (rad)
E - eccentricity vector
e - eccentricity (magnitude of E)
eps - a small number below which the eccentricity is considered
to be zero
w - argument of perigee (rad)
TA - true anomaly (rad)
a - semimajor axis (km)
coe - vector of orbital elements [a e incl RA w f h]
%}

eps = 1.e-10;
r = norm(R);
v = norm(V);
vr = dot(R,V)/r;
H = cross(R,V);
h = norm(H);

incl = acos(H(3)/h);

N = cross([0 0 1],H);
n = norm(N);

if n ~=  0
RA = acos(N(1)/n);
if N(2) < 0
RA = 2*pi - RA;
end
else
RA = 0;
end

E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);

if n ~= 0
if e > eps
w = acos(dot(N,E)/n/e);
if E(3) < 0
w = 2*pi - w;
end
else
w = 0;
end
else
w = 0;
end

if e > eps
f = acos(dot(E,R)/e/r);
if vr < 0
f = 2*pi - f;
end
else
cp = cross(N,R);
if cp(3) >= 0
f = acos(dot(N,R)/n/r);
else
    f= 2*pi - acos(dot(N,R)/n/r);
end
end

a = h^2/mu/(1 - e^2);
coe = [a e incl RA w f];
end 