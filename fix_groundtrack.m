function a_fixed = fix_groundtrack (m,k,w_e,mu)
%fix_groundtrack.m
% 
% DESCRIPTION: 
%             This function gives as output the semi major axis for the 
%             unperturbed  case.
%             
%             
% INPUT:
%    m
%    k
%    w_e           Earth's angular velocity [rad/s]
%    mu            Planetary Constant       [km^3/s^2]

% OUTPUT: 
%     a_fixed      Semi major axis [km]
%
% AUTHORS: Biagetti Giorgia, Galzignato Sara, Ardi Nugraha Setya, Tagnin Emanuele

n = w_e *(k/m);             % mean motion

a_fixed = (mu/(n)^2)^(1/3); %semi major axis

end

