function [r, v] = kep2car(coe,mu)
%kep2car.m
% 
% DESCRIPTION: 
%             This function computes the state vector (r,v) from the
%             classical orbital elements (coe).
%             
%             
% INPUT:
 
% mu       gravitational parameter (km^3/s^2)
% coe      orbital elements [a e i RA w f]
% where 
%
%        - a = semi-major axis (km^2/s)
%        - e = eccentricity
%        - i = inclination of the orbit (rad)
%        - RA = right ascension of the ascending node (rad)
%        - w = argument of perigee (rad)
%        - f = true anomaly (rad)


% OUTPUT:

% r     position vector in the geocentric equatorial frame (km)
% v     velocity vector in the geocentric equatorial frame (km/s)
%
% AUTHORS: Biagetti Giorgia, Galzignato Sara, Ardi Nugraha Setya, Tagnin Emanuele



a = coe(1);
e = coe(2);
i = coe(3);
RA = coe(4);
w = coe(5);
f = coe(6);

h=sqrt(mu*a*(1-e^2));

rp = (h^2/mu) * (1/(1 + e*cos(f))) * (cos(f)*[1;0;0] + sin(f)*[0;1;0]); %position vector in the perifocal frame (km)
vp = (mu/h) * (-sin(f)*[1;0;0] + (e + cos(f))*[0;1;0]); %velocity vector in the perifocal frame (km/s)

R3_W = [ cos(RA) sin(RA) 0; -sin(RA) cos(RA) 0; 0 0 1]; %Rotation matrix about the z-axis through the angle RA

R1_i= [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];       %Rotation matrix about the x-axis through the angle i

R3_w = [ cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1];     %Rotation matrix about the z-axis through the angle w
 
Q_pX = (R3_w*R1_i*R3_W);                                %Matrix of the transformation from perifocal to geocentric equatorial frame

r = Q_pX'*rp;
v = Q_pX'*vp;
% Convert r and v into row vectors:
r = r';
v = v';
end