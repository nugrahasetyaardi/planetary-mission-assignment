function dy = perturbed_J2_D(y,mu)
%perturbed_J2_D.m
% 
%  
% DESCRIPTION: 
%             This function gives as output dy, the system differential
%             equations, governing the dynamics of the satellite in the
%             perturbed case due to the presence of both J2 and drag perturbations.
%             It also takes into account the exponential model of density.
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
om_E = [0,0,2*pi/(24*3600)];
AM_ratio = 0.080;
Cd = 2.2;

density=Nominal_Density;

altitude=r-R_E;
rho_zero=0;
altitude_0=0;
H=0;
for q=1:27
    if  altitude>=density(q,1) &&  altitude<density(q+1,1)
    rho_zero=density(q,2);
    altitude_0=density(q,1);
    H=density(q,3);
    elseif altitude>=1000
        rho_zero=density(28,2);
        altitude_0=density(28,1);
        H=density(28,3);
    end 
end
rho=rho_zero*exp(-(altitude-altitude_0)/H);
 
K = 0.5*AM_ratio*Cd*rho;
V = [y(4),y(5),y(6)];
R = [y(1),y(2),y(3)];

v_rel = V - cross (om_E,R);

a_drag = K.*(norm(v_rel)^2).*v_rel./norm(v_rel);

dy=[y(4); y(5); y(6); -mu/r^3*y(1)+3/2*((J2*mu*R_E^2)/(r)^4)*((y(1)/r)*(5*(y(3)^2/r^2)-1)) - a_drag(1); -mu/r^3*y(2)+3/2*((J2*mu*R_E^2)/(r)^4)*((y(2)/r)*(5*(y(3)^2/r^2)-1))- a_drag(2); -mu/r^3*y(3)+3/2*((J2*mu*R_E^2)/(r)^4)*((y(3)/r)*(5*(y(3)^2/r^2)-3))-a_drag(3)];

% dy = [y(4); y(5); y(6); -mu/r^3*y(1) - a_drag(1); -mu/r^3*y(2) - a_drag(2); -mu/r^3*y(3) - a_drag(3)]; 
end