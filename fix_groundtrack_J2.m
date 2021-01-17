function a_J2 = fix_groundtrackJ(w_e,i,J2,e,mu,R_e,m,k,a_guess)
% fix_groundtrackJ.m
% 
% DESCRIPTION: 
%             This function computes the semi major-axis in order to obtain a match of the initial and final position
%             of the groundtrack within k orbits. The function also takes takes into account the J2 perturbation term.
%             
%             
% INPUT:
%     w_e[1]        Earth angular velocity [rad/s]
%     i[1]          Orbit inclination [rad]
%     J2[1]         Perturbation term
%     e[1]          Orbit eccentricity 
%     mu[1]         Planetary constant[km^3/s^2]
%     R_e[1]        Radius of the Earth[km]
%     m[1]          Number of Earth's rotation
%     k[1]          Number of satellite's orbit
%     a_guess[1]    Initial guess used in fzero (it can be used the unperturbed semi major axis)
%
%
% AUTHORS: Biagetti Giorgia, Galzignato Sara, Ardi Nugraha Setya, Tagnin Emanuele
    

K = (3/2 *sqrt(mu)*J2*R_e^2)/((1-e^2)^2);

funz = @(a) (w_e + (K.*cos(i))./(a.^(7/2)))./(sqrt(mu./(a.^3))- K.*(5/2 .*sin(i).^2 -2)./(a.^(7/2))+ 3/2*(sqrt(mu)*J2*R_e^2)/( (1-e^2)^1.5 * a^(7/2) )*(1-1.5*sin(i)^2)) - m/k;

a_J2 = fzero(funz, a_guess);



end