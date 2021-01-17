function dy=perturbed_propagator(y,mu)
%perturbed_propagator.m
% 
% DESCRIPTION: 
%             This function gives as output dy, the system differential
%             equations, governing the dynamics of the satellite in the
%             perturbed case due to the presence of J2.
%             
%             
% INPUT:
%     y             Vector of position and velocity
%     mu            Planetary Constant
% OUTPUT: 
%     dy            System differential equations for the perturbed case
%     
%
% AUTHORS: Biagetti Giorgia, Galzignato Sara, Ardi Nugraha Setya, Tagnin Emanuele

r=sqrt(y(1)^2+y(2)^2+y(3)^2);
R_E=astroConstants(23);
J2=astroConstants(9);
dy=[y(4); y(5); y(6); -mu/r^3*y(1)+3/2*((J2*mu*R_E^2)/(r)^4)*((y(1)/r)*(5*(y(3)^2/r^2)-1)); -mu/r^3*y(2)+3/2*((J2*mu*R_E^2)/(r)^4)*((y(2)/r)*(5*(y(3)^2/r^2)-1)); -mu/r^3*y(3)+3/2*((J2*mu*R_E^2)/(r)^4)*((y(3)/r)*(5*(y(3)^2/r^2)-3))];
end