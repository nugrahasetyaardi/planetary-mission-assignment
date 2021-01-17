function [dkep_dt]=propagate_Gauss_J2_D(t,kep) 
%propagate_Gauss_J2_D.m
% 

% DESCRIPTION: 
%             This function gives as output dkep_dt, the system
%             differential Gauss equations, representing the variation of 
%             the keplerian elements of the satellite in time in the perturbed
%             case due to the presence of both J2 and drag perturbations.
%             It also takes into account the exponential model of density.
%             
% INPUT:
%     t             Time vector
%     kep           Keplerian elements [a,e,i,RAAN, omega,f]
% OUTPUT: 
%     dkep_dt       System differential Gauss equations for the perturbed case
%     
%
% AUTHORS: Biagetti Giorgia, Galzignato Sara, Ardi Nugraha Setya, Tagnin Emanuele



muE=astroConstants(13);
R_E=astroConstants(23);
J2=astroConstants(9);


a=kep(1);
e=kep(2);
i=kep(3);
RAAN=kep(4);
omega=kep(5);
f=kep(6);


b=a*sqrt(1-e^2);
n=sqrt(muE/a^3);
h=n*a*b;
p=b^2/a;
r=p/(1+e*cos(f));


[r_inertial,v_inertial]=kep2car(kep,muE);

om_E = [0,0,2*pi/(24*3600)];
AM_ratio = 0.080;
Cd = 2.2;

density=Nominal_Density;

altitude=norm(r_inertial)-R_E;
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
    else
    end
end
rho=rho_zero*exp(-(altitude-altitude_0)/H);

K = 0.5*AM_ratio*Cd*rho;
v_rel = v_inertial - cross (om_E,r_inertial);


a_drag_inertial = K.*(norm(v_rel)^2).*v_rel./norm(v_rel);



%Creation of rotation matrix

R_3_RAAN=[cos(RAAN), sin(RAAN), 0; -sin(RAAN), cos(RAAN), 0; 0,0,1];


R_1_i=[1,0,0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];


R_3_omegaf=[cos(omega+f), sin(omega+f), 0; -sin(omega+f), cos(omega+f), 0; 0,0,1];


a_drag_rsw=R_3_omegaf*R_1_i*R_3_RAAN*a_drag_inertial';
a_drag_r=a_drag_rsw(1);
a_drag_s=a_drag_rsw(2);
a_drag_w=a_drag_rsw(3);


ar=-3/2*((J2*muE*R_E^2)/(r^4))*(1-3*sin(i)^2*sin(f+omega)^2)-a_drag_r;
as=-3/2*((J2*muE*R_E^2)/(r^4))*(sin(i)^2*sin(2*(f+omega)))-a_drag_s;
aw=-3/2*((J2*muE*R_E^2)/(r^4))*(sin(2*i)*sin(f+omega))-a_drag_w;



da=(2*a^2/h)*(e*sin(f)*ar+(p/r)*as);
de=(1/h)*(p*sin(f)*ar+((p+r)*cos(f)+r*e)*as);
di=(r*cos(f+omega))*aw/h;
dRAAN=(r*sin(f+omega))*aw/(h*sin(i));
domega=1/(h*e)*(-p*cos(f)*ar+(p+r)*sin(f)*as)-(r*sin(f+omega)*cos(i))*aw/(h*sin(i));
df=h/r^2+(1/(e*h))*(p*cos(f)*ar-(p+r)*sin(f)*as);


dkep_dt=[da; de; di; dRAAN; domega; df];
end 
