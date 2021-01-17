function dy=propagator(y,mu)
%propagator.m
% 
% DESCRIPTION: 
%             This function gives as output dy, the system differential
%             equations, governing the dynamics of the satellite for the unperturbed case.
%             
%             
% INPUT:
%     y             Vector of position and velocity
%     mu            Planetary Constant [km^3/s^2]
% OUTPUT: 
%     dy            System differential equations
%     
%
% AUTHORS: Biagetti Giorgia, Galzignato Sara, Ardi Nugraha Setya, Tagnin Emanuele

r=sqrt(y(1)^2+y(2)^2+y(3)^2);

dy=[y(4); y(5); y(6); -(mu/r^3)*y(1); -(mu/r^3)*y(2); -(mu/r^3)*y(3)];
end 

    