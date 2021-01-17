function [alpha_grad, delta_grad, lat_grad, lon_grad] = groundtrack (r,t,w_e,theta_g0,t0)

%groundtrack.m
% 
% DESCRIPTION: 
%             This function computes right ascension and declination in ECEI reference frame, 
%             and longitute and latitude in ECEF reference frame.
%             
%             
% INPUT:
%     r             Position vector [km]
%     t             Time form integration [s]
%     w_e[1]        Earth angular velocity [rad/s]
%     theta_g0      Greenwich theta [rad]
%     t0            Initial time [s]
%     
%
% OUTPUT:
%     alpha_grad    Right ascension [deg]
%     delta_grad    Declination     [deg]
%     lat_grad      Latitude        [deg]
%     lon_grad      Longitute       [deg]
     

% AUTHORS: Biagetti Giorgia, Galzignato Sara, Ardi Nugraha Setya, Tagnin Emanuele
    
for i = 1:length(t)
    norm_r(i) = norm(r(i,:));
    delta(i) = asin(r(i,3)/norm_r(i));
    
    delta_grad(i) = rad2deg(delta(i));
    
    if (r(i,2))/norm_r(i) >0
        alpha(i) = acos((r(i,1)/norm_r(i))/cos(delta(i)));
        alpha_grad(i) = rad2deg(alpha(i));
    else
        alpha(i) = 2*pi - acos((r(i,1)/norm_r(i))/cos(delta(i)));
        alpha_grad(i) = rad2deg(alpha(i));
    end

lat_grad(i) = delta_grad(i);
theta_g(i) =  theta_g0 + w_e.*(t(i)-t0);
theta_g_grad(i) = rad2deg(theta_g(i));
lon(i) = alpha_grad(i) - theta_g_grad(i);
lon_grad(i) = wrapTo180(lon(i));


end
end 