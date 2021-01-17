% PlanetaryMission_group_39.m
% 
% DESCRIPTION:
%             This is the main script that computes the Ground Track(in
%             unperturbed and perturbed conditions)for a chosen orbit and 
%             propagate the keplerian elements comparing two different 
%             methods, applying J2 and drag prsturbations:
%                   - Integration of the equation of motion in Cartesian 
%                     coordinates and 
%                   - Gauss planetary equations.
%              
%             Then a low-pass filter is applied. 
%             
%             In the last sections it compares the two method with real
%             data coming from DSX satellite. 
%               
%     
% AUTHORS: Giorgia Biagetti, Sara Galzignato, Nugraha Setya Ardi, Emanuele Tagnin 

   clear 
   close all
   clc 
%% DATA
muE = astroConstants(13) ;      %Earth's gravitational constant [km^3/s^2]
R_E=astroConstants(23);         %Radius of the Earth [km]
J2=astroConstants(9);           %Disturbance term due to Earth's oblateness
theta_g0 =0 ;                   %Greenwich theta [rad]
w_e = deg2rad((15.04)/(3600));  %Earth's angular velocity [rad/s]
t0  = 0;                        %initial time [s]
k = 11;
m = 2;
day = 24*3600;                  %just plotting for on day [s]

%% PROJECT ORBITAL ELEMENTS
a=13336;                        %semi major axis [km]
e=0.4785;                       %eccentricity [-]
i=deg2rad(32.4346);             %orbit inclination [rad]
RAAN=deg2rad(23.4355);          %orbit RAAN [rad]
omega=deg2rad(10.32);           %pericenter anomaly [rad]
f=deg2rad(90);                  %true anomaly [rad]

%% Kep2car and orbit propagation with 2BP: UNPERTURBED CASE 
n = 10000;                        %number of iteration
T_orb=2*pi*sqrt(a^3/muE);         %orbital period
tspan=linspace(0,T_orb,n);        %time vector
kep0=[a, e, i, RAAN, omega, f ];  %initial vector of keplerian element
[r0, v0] = kep2car(kep0,muE);     %initial state vector: position and velocity

y0=[r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];   %initial conditions

options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[t, Y_unp] = ode113( @(t,y)propagator(y, muE), tspan, y0, options);

R_unp=[Y_unp(:,1), Y_unp(:,2), Y_unp(:,3)]; % position vector in unperturbed case

%% Orbit plot: UNPERTURBED CASE 
% plots one orbit in the unperturbed case
figure (1)
title('Orbit Plot')
xlabel ('x [{\itkm}]')
ylabel ('y [{\itkm}]')
zlabel ('z [{\itkm}]')
hold on 
grid on 
plot3(R_unp(:,1),R_unp(:,2),R_unp(:,3),'y','Linewidth',2);
mArrow3([0 0 0],[12000 0 0], 'facealpha', 0.5, 'color', 'red', 'stemWidth', 100);
mArrow3([0 0 0],[0 12000 0], 'facealpha', 0.5, 'color', 'red', 'stemWidth', 100);
mArrow3([0 0 0],[0 0 12000], 'facealpha', 0.5, 'color', 'red', 'stemWidth', 100);
[x, y, z] = sphere (180);
 
X1 = x * R_E;
Y_unp = y * R_E;
Z1 = z * R_E;
Earth = surf(X1,Y_unp,Z1);
cdata = imread('EarthTexture.jpg'); 
set(Earth, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
set(gca,'GridColor',[1 1 1])
set(gca,'Color','k')
set(gca,'GridColor',[1 1 1])
view (45,45)
axis equal

figure (2)
title ('Orbit Plot GIF')
xlabel ('x [{\itkm}]')
ylabel ('y [{\itkm}]')
zlabel ('z [{\itkm}]')
hold on 
grid on 
h = [];
earth = [];
mArrow3([0 0 0],[12000 0 0], 'facealpha', 0.5, 'color', 'red', 'stemWidth', 100);
mArrow3([0 0 0],[0 12000 0], 'facealpha', 0.5, 'color', 'red', 'stemWidth', 100);
mArrow3([0 0 0],[0 0 12000], 'facealpha', 0.5, 'color', 'red', 'stemWidth', 100);
for ii = 1:200:n
    delete(h)
    h = plot3(R_unp(ii,1),R_unp(ii,2),R_unp(ii,3),'mo','MarkerFaceColor','m');
    plot3(2e4,2e4,2e4,'w',-2e4,-2e4,-2e4,'w',2e4,-2e4,2e4,'w');
    hold on 
    h1 = plot3(R_unp(ii,1),R_unp(ii,2),R_unp(ii,3),'y');
    drawnow
    delete(earth) 
    [x, y, z] = sphere (180);
   
    X1 = x * R_E;
    Y_unp = y * R_E;
    Z1 = z * R_E;
    earth = surf(X1,Y_unp,Z1);
    rot = ii*w_e*200;
    rotate(earth,[0 0 1],rot)
    cdata = imread('EarthTexture.jpg'); 
    set(earth, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    hold on 
    axis equal
    set(gca,'GridColor',[1 1 1])
    set(gca,'Color','k')
    set(gca,'GridColor',[1 1 1])
    view (45,45)
    drawnow
   
   
end 
hold off

%% Ground Track - UNPERTURBED CASE
% plot the Ground Track for one orbit in the unperturbed case
t1_orbit = linspace(0,T_orb,n);
[alpha_grad, delta_grad, lat_grad, lon_grad] = groundtrack(R_unp,t,w_e,theta_g0,t0);

img = imread('planisphererid5.jpg');
min_x = -180;
max_x = 180;
min_y = -90;
max_y = 90;

figure (3)

imagesc([min_x max_x], [min_y max_y], flipud(img));
hold on 

for jj = 1:length(lon_grad)-1 %this cicle is to avoid the line at the discontinuity 
    tol = 100;
    if abs(lon_grad(jj)-lon_grad(jj+1))> tol
        lon_grad(jj) = NaN;
    end
end    

%plotting groundtrack on the planysphere
plot (lon_grad,lat_grad,'y','Linewidth',2)
hold on 
plot(lon_grad(1),lat_grad(1),'r o',lon_grad(end),lat_grad(end),'m o','MarkerSize',6,'LineWidth',2);
legend ('Ground track','Initial point', 'Final point')
xlabel('longitude [°]');
ylabel('lattitude [°]');
title ('Ground track - unperturbed case')
% set the y-axis back to normal.
set(gca,'ydir','normal');


%% Fix Ground track - UNPERTURBED CASE
% plot the Ground Track in the unperturbed case for 11 full orbits every 2 Eart's rotations 
a_fix = fix_groundtrack (m,k,w_e,muE); %correction of semimajor-axis

kep_fix = [a_fix,e,i,RAAN,omega,f];
[r0_fix,v0_fix]= kep2car(kep_fix,muE);
T_orb_fix = 2*pi*sqrt(a_fix^3/muE);
t2_orbit = linspace(0,T_orb_fix*k,n); %plotting for k  orbital period

y0_fix=[r0_fix(1); r0_fix(2); r0_fix(3); v0_fix(1); v0_fix(2); v0_fix(3)];

options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[~, Y_fix] = ode113( @(t,y)propagator(y, muE), t2_orbit, y0_fix, options);

R_fix=[Y_fix(:,1), Y_fix(:,2), Y_fix(:,3)]; % state vectors
V_fix=[Y_fix(:,4), Y_fix(:,5), Y_fix(:,6)];

[alpha_grad_fix, delta_grad_fix, lat_grad_fix, lon_grad_fix] = groundtrack(R_fix,t2_orbit,w_e,theta_g0,t0);

img = imread('planisphererid5.jpg');
figure (4)

imagesc([min_x max_x], [min_y max_y], flipud(img));
hold on 

for ll = 1:length(lon_grad_fix)-1 %this cicle is to avoid the line at the discontinuity 
    tol = 100;
    if abs(lon_grad_fix(ll)-lon_grad_fix(ll+1))> tol
        lon_grad_fix(ll) = NaN;
    end
end    

plot (lon_grad_fix,lat_grad_fix,'y','Linewidth',2)
hold on 
plot(lon_grad_fix(1),lat_grad_fix(1),'r o',lon_grad_fix(end),lat_grad_fix(end),'m o','MarkerSize',6,'LineWidth',2);
legend ('Fixed Ground track','Initial point', 'Final point')
xlabel('longitude [°]');
ylabel('lattitude [°]');
title ('Fixed Ground track - unperturbed case')
% set the y-axis back to normal
set(gca,'ydir','normal');


%% Fix Ground track - J2 case
% add J2 perturbation
a_guess = a; % strarting point for fzero
a_J2 = fix_groundtrack_J2(w_e,i,J2,e,muE,R_E,m,k,a_guess); % fixed semi-major axis

kep_J2 = [a_J2,e,i,RAAN,omega,f];
[r0_J2,v0_J2]= kep2car(kep_J2,muE);
T_orb_J2 = 2*pi*sqrt(a_J2^3/muE);
t3_orbit = linspace(0,T_orb_J2*k,n);
y0_J2 =[r0_J2(1); r0_J2(2); r0_J2(3); v0_J2(1); v0_J2(2); v0_J2(3)];

options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[t, Y_J2] = ode113( @(t,y)perturbed_propagator(y,muE), t3_orbit, y0_J2, options);

R_J2=[Y_J2(:,1), Y_J2(:,2), Y_J2(:,3)]; % state vectors 
V_J2=[Y_J2(:,4), Y_J2(:,5), Y_J2(:,6)];

[alpha_grad_J2, delta_grad_J2, lat_grad_J2, lon_grad_J2] = groundtrack(R_J2,t3_orbit,w_e,theta_g0,t0);

img = imread('planisphererid5.jpg');
figure (5)

imagesc([min_x max_x], [min_y max_y], flipud(img));
hold on

for ff = 1:length(lon_grad_J2)-1 %this cicle is to avoid the line at the discontinuity 
    tol = 100;
    if abs(lon_grad_J2(ff)-lon_grad_J2(ff+1))> tol
        lon_grad_J2(ff) = NaN;
    end
end    
 % plots the ground track in the unperturbed case for 11 full orbits 
plot (lon_grad_J2,lat_grad_J2,'y','Linewidth',2)
hold on 
plot(lon_grad_J2(1),lat_grad_J2(1),'r o',lon_grad_J2(end),lat_grad_J2(end),'m o','MarkerSize',6,'LineWidth',2);
legend ('Fixed Ground track','Initial point', 'Final point')
xlabel('longitude [°]');
ylabel('lattitude [°]');
title ('Ground track - J2')

% set the y-axis back to normal
set(gca,'ydir','normal');


%% J2 and Drag perturbation
% plotting for a larger number of orbital periods
n=10000;
k=100;
T_orb=2*pi*sqrt(a^3/muE);
tspan=linspace(0,k*T_orb,n);
 
kep0=[a, e, i, RAAN, omega, f ];
[r0, v0] = kep2car(kep0,muE);
y0=[r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];
 
options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

%Propagator in cartesian coordinates UNPERTURBED case
[t, Y_unp_1] = ode113( @(t,y)propagator(y, muE), tspan, y0, options);
%Propagation with Gauss planetary equations PERTURBED case
[T,kep_Gauss] = ode113( 'propagate_Gauss_J2_D', tspan, kep0, options);
%Propagator in cartesian coordinates PERTURBED case
[t, Y] = ode113( @(t,y)perturbed_J2_D(y, muE), tspan, y0, options);
%unperturbed position vector
R_unp_1=[Y_unp_1(:,1), Y_unp_1(:,2), Y_unp_1(:,3)];

%perturbed position and velocity vectors
R=[Y(:,1), Y(:,2), Y(:,3)];
V=[Y(:,4), Y(:,5), Y(:,6)];

kep=zeros(n,6);
for i=1:n
[kep(i,:)] = car2kep(R(i,:),V(i,:),muE);
end

%% Definig ERRORS on keplerian elements 
% error in the keplerian elements for the two methods 
a_error=abs((kep(:,1)-kep_Gauss(:,1))/kep0(1));
e_error=abs((kep(:,2)-kep_Gauss(:,2)));
i_error=abs((unwrap(kep(:,3))-kep_Gauss(:,3))/2*pi);
RAAN_error=abs((unwrap(kep(:,4))-kep_Gauss(:,4))/2*pi);
omega_error=abs((unwrap(kep(:,5))-kep_Gauss(:,5))/2*pi);
f_error=zeros(n,1);
for i=1:n
    f_error(i)=abs((kep(i,6)-wrapTo2Pi(kep_Gauss(i,6)))); 
end 

%Creation of state matrix from kep_Gauss 
r_Gauss=zeros(n,3);
v_Gauss=zeros(n,3);
for i=1:n
[r_Gauss(i,:), v_Gauss(i,:)] = kep2car(kep_Gauss(i,:),muE);
end 

% plots perturbed vs unperturbed orbits
figure (6)
title ('Perturbed Orbit Plot vs Unperturbed GIF')
xlabel ('x [{\itkm}]')
ylabel ('y [{\itkm}]')
zlabel ('z [{\itkm}]')
hold on 
grid on 
h = [];
h1 = [];
earth = [];
for ii = 1:250

    h = plot3(R(ii,1),R(ii,2),R(ii,3),'mo','MarkerFaceColor','m','MarkerSize',5);
      


    hold on 
    h1 = plot3(R_unp_1(ii,1),R_unp_1(ii,2),R_unp_1(ii,3),'go','MarkerFaceColor','g','MarkerSize',5);

hold on 
    h1 = plot3(R_unp(ii,1),R_unp(ii,2),R_unp(ii,3),'y');
    drawnow
    delete(earth)
    
      [x, y, z] = sphere (180);
   
    X1 = x * R_E;
    Y_unp = y * R_E;
    Z1 = z * R_E;
    earth = surf(X1,Y_unp,Z1);
    rot = ii*w_e*200;
    rotate(earth,[0 0 1],rot)
    cdata = imread('EarthTexture.jpg'); 
    set(earth, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');

    axis equal
    set(gca,'GridColor',[1 1 1])
    set(gca,'Color','k')
    set(gca,'GridColor',[1 1 1])
    view (45,45)
    drawnow
   
end 
legend('Perturbed orbit','Unperturbed orbit')
legend('boxoff')
hold off

figure (7)
%Unperturbed orbit
plot3(R_unp_1(:,1),R_unp_1(:,2),R_unp_1(:,3),'r','LineWidth',2)
hold on
%Perturbed orbit (Cartesian)
plot3 (R(:,1),R(:,2),R(:,3),'g');
hold on
%Perturbed orbit (Gauss)
plot3(r_Gauss(:,1),r_Gauss(:,2),r_Gauss(:,3),'b')
xlabel ('x [{\itkm}]')
ylabel ('y [{\itkm}]')
zlabel ('z [{\itkm}]')

hold on 

grid on
X1 = x * R_E;
Y1 = y * R_E;
Z1 = z * R_E;
Earth = surf(X1,Y1,Z1);
cdata = imread('EarthTexture.jpg'); 
set(Earth, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
axis equal
title ('Orbit perturbation - Long time window')
set(gca,'GridColor',[1 1 1])
set(gca,'Color','k')
set(gca,'GridColor',[1 1 1])
view (45,45)
legend('Unperturbed orbit','Perturbed orbit (Cartesian)','Perturbed orbit (Gauss)' )
legend('boxoff')

%% ERROR GAUSS-CARTESIAN
% plot the errors between the two methods in semilogarithmic scale
n=10000;
figure (8)
semilogy(linspace(0,k,n),a_error)
hold on 
grid
xlabel( 'time [{\itT}]');
ylabel('|a_C_a_r - a_G_a_u_s_s|/a_0 [{\it-}] ');
hold off
figure (9)
semilogy(linspace(0,k,n),e_error)
hold on 
grid
xlabel( 'time [{\itT}]');
ylabel('|e_C_a_r - e_G_a_u_s_s| [{\it-}] ');
hold off
figure (10)
semilogy(linspace(0,k,n),i_error)
hold on 
grid
xlabel( 'time [{\itT}]');
ylabel('|i_C_a_r - i_G_a_u_s_s|/2\pi [{\it-}] ');
hold off
figure (11)
semilogy(linspace(0,k,n),RAAN_error)
hold on 
grid
xlabel( 'time [{\itT}]');
ylabel('|\Omega_C_a_r - \Omega_G_a_u_s_s|/2\pi [{\it-}] ');
hold off
figure (12)
semilogy(linspace(0,k,n),omega_error)
hold on 
grid
xlabel( 'time [{\itT}]');
ylabel('|\omega_C_a_r - \omega_G_a_u_s_s|/2\pi [{\it-}] ');
hold off
figure (13)
semilogy(linspace(0,k,n),rad2deg(f_error))
hold on 
grid
xlabel( 'time [{\itT}]');
ylabel('|{\itf}_C_a_r - {\itf}_G_a_u_s_s| [{\itdeg}] ');
hold off

%% Filter of the keplerian elements
% Cut-off period
Tfilter=3*T_orb;
% Number of points for the filtering window
nwindow=nearest(Tfilter/(sum(diff(tspan))/(length(tspan)-1))); % Tfilter/(average length of the time step)
% Filter elements
kep_Gauss_filtered=movmean(kep_Gauss, nwindow,1);

% plot keplerian element evolutions in both methods + filter
figure (14) 
plot(linspace(0,k,n),kep_Gauss(:,1))
hold on
plot(linspace(0,k,n),kep(:,1))
plot(linspace(0,k,n),kep_Gauss_filtered(:,1),'k')
grid
xlabel( 'time [{\itT}]');
ylabel('a [{\itkm}]');
legend('Gauss equations','Cartesian','Secular evolution (filtered)');
hold off

figure (15) 
plot(linspace(0,k,n),kep_Gauss(:,2))
hold on
plot(linspace(0,k,n),kep(:,2))
plot(linspace(0,k,n),kep_Gauss_filtered(:,2),'k')
grid
xlabel( 'time [{\itT}]');
ylabel('e [{\it-}]');
legend('Gauss equations','Cartesian','Secular evolution (filtered)');
hold off

figure (16)
plot(linspace(0,k,n),rad2deg(kep_Gauss(:,3)))
hold on
plot(linspace(0,k,n),rad2deg(unwrap(kep(:,3))))
plot(linspace(0,k,n),rad2deg(kep_Gauss_filtered(:,3)),'k')
grid
xlabel( 'time [{\itT}]');
ylabel('i [{\itdeg}]');
legend('Gauss equations','Cartesian','Secular evolution (filtered)');
hold off

figure (17) 
plot(linspace(0,k,n),rad2deg(kep_Gauss(:,4)),'LineWidth',2)
hold on
plot(linspace(0,k,n),rad2deg(unwrap(kep(:,4))),'LineWidth',2)
plot(linspace(0,k,n),rad2deg(kep_Gauss_filtered(:,4)),'k')
grid
xlabel( 'time [{\itT}]');
ylabel('\Omega [{\itdeg}]');
legend('Gauss equations','Cartesian','Secular evolution (filtered)');
hold off

figure (18) 
plot(linspace(0,k,n),rad2deg(kep_Gauss(:,5)),'LineWidth',2)
hold on
plot(linspace(0,k,n),rad2deg(unwrap(kep(:,5))),'LineWidth',2)
plot(linspace(0,k,n),rad2deg(kep_Gauss_filtered(:,5)),'k')
grid
xlabel( 'time [{\itT}]');
ylabel('\omega [{\itdeg}]');
legend('Gauss equations','Cartesian','Secular evolution (filtered)');
hold off

figure (19) 
plot(linspace(0,k,n),rad2deg(kep_Gauss(:,6)),'LineWidth',2)
hold on
plot(linspace(0,k,n),rad2deg(unwrap(kep(:,6))),'LineWidth',2)
plot(linspace(0,k,n),rad2deg(kep_Gauss_filtered(:,6)),'k')
grid
xlabel( 'time [{\itT}]');
ylabel('{\itf} [{\itdeg}]');
legend('Gauss equations','Cartesian','Secular evolution (filtered)');
hold off

%% COMPARISON WITH REAL SATELLITE

% ******************************************************************************* 
%                         DSX satellite
%
% Start time      : A.D. 2019-Jul-01 00:00:00.0000 TDB
% Stop  time      : A.D. 2021-Jan-01 00:00:00.0000 TDB
% Step-size       : 10000 steps
% *******************************************************************************

%reading text file
horizon_timetable=readtimetable('horizons_results.txt');
horizon_table=timetable2table(horizon_timetable);
horizon_table = removevars(horizon_table, {'Time','Var2','Var3','Var5','Var9','Var10','Var11','Var14','Var15'});
horizon_cell=table2array(horizon_table);

% extracting keplerian elements at the starting point
e_sat=str2num(cell2mat(horizon_cell(1,2)));
i_sat=str2num(cell2mat(horizon_cell(1,3)));
i_sat=deg2rad(i_sat);
RAAN_sat=str2num(cell2mat(horizon_cell(1,4)));
RAAN_sat=deg2rad(RAAN_sat);
omega_sat=str2num(cell2mat(horizon_cell(1,5)));
omega_sat=deg2rad(omega_sat);
f_sat=str2num(cell2mat(horizon_cell(1,6)));
f_sat=deg2rad(f_sat);
a_sat=str2num(cell2mat(horizon_cell(1,7)));

% define time window
t_start=24*60*60*str2num(cell2mat(horizon_cell(1,1)));
t_end=24*60*60*str2num(cell2mat(horizon_cell(end,1)));
dt=t_end-t_start;
n_it=max(size(horizon_cell));
tspan=linspace(0,dt,n_it);

% retrieve keplerian elements with numerical methods
kep0_sat=[a_sat, e_sat, i_sat, RAAN_sat, omega_sat, f_sat ];
[r0_sat, v0_sat] = kep2car(kep0_sat,muE);
y0_sat=[r0_sat(1); r0_sat(2); r0_sat(3); v0_sat(1); v0_sat(2); v0_sat(3)];
 
options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

%Propagation with Gauss planetary equations PERTURBED case
[T,kep_Gauss_sat] = ode113( 'propagate_Gauss_J2_D', tspan, kep0_sat, options);
%Propagator in cartesian coordinates PERTURBED case
[t, Y_sat] = ode113( @(t,y)perturbed_J2_D(y, muE), tspan, y0_sat, options);

%PERTURBED POSITION VECTOR (CARTESIAN)
R_sat=[Y_sat(:,1), Y_sat(:,2), Y_sat(:,3)];
V_sat=[Y_sat(:,4), Y_sat(:,5), Y_sat(:,6)];

kep_sat=zeros(n_it,6);

for i=1:n_it
[kep_sat(i,:)] = car2kep(R_sat(i,:),V_sat(i,:),muE);
end

% retrieving keplerian elemnts from ephemerides and TLEs
e_horizon_char=cell2mat(horizon_cell(:,2));
i_horizon_char=cell2mat(horizon_cell(:,3));
RAAN_horizon_char=cell2mat(horizon_cell(:,4));
omega_horizon_char=cell2mat(horizon_cell(:,5));
f_horizon_char=cell2mat(horizon_cell(:,6));
a_horizon_char=cell2mat(horizon_cell(:,7));


e_horizon=str2num(e_horizon_char);
i_horizon=str2num(i_horizon_char);
RAAN_horizon=str2num(RAAN_horizon_char);
omega_horizon=str2num(omega_horizon_char);
f_horizon=str2num(f_horizon_char);
a_horizon=str2num(a_horizon_char);

i_horizon=deg2rad(i_horizon);
RAAN_horizon=deg2rad(RAAN_horizon);
omega_horizon=deg2rad(omega_horizon);
f_horizon=deg2rad(f_horizon);

kep_horizon=[a_horizon, e_horizon, i_horizon, RAAN_horizon, omega_horizon, f_horizon];
%% COMPARING KEPLERIAN ELEMENTS FOUND BY HORIZON AND kep_sat
tol=deg2rad(280); % to avoid picks in plots due to step changes 

% error between real and numerical keplerian elements (Cartesian)
a_error_sat_Car=abs(a_horizon-kep_sat(:,1))/kep0_sat(1);
e_error_sat_Car=abs(e_horizon-kep_sat(:,2));

i_error_sat_Car=zeros(n_it,1);
for i=1:n_it
   if abs(i_horizon(i)-kep_sat(i,3)) > tol 
      i_error_sat_Car(i)=(2*pi-abs(i_horizon(i)-kep_sat(i,3)))/2*pi;
   else
      i_error_sat_Car(i)=abs(i_horizon(i)-kep_sat(i,3))/2*pi;
   end 
end

RAAN_error_sat_Car=zeros(n_it,1);
for i=1:n_it
   if abs(RAAN_horizon(i)-kep_sat(i,4)) > tol
      RAAN_error_sat_Car(i)=(2*pi-abs(RAAN_horizon(i)-kep_sat(i,4)))/2*pi;
   else
      RAAN_error_sat_Car(i)=abs(RAAN_horizon(i)-kep_sat(i,4))/2*pi;
   end 
end

omega_error_sat_Car=zeros(n_it,1);
for i=1:n_it
   if abs(omega_horizon(i)-kep_sat(i,5)) > tol
      omega_error_sat_Car(i)=(2*pi-abs(omega_horizon(i)-kep_sat(i,5)))/2*pi;
   else
      omega_error_sat_Car(i)=abs(omega_horizon(i)-kep_sat(i,5))/2*pi;
   end 
end

f_error_sat_Car=zeros(n_it,1);
for i=1:n_it
   if abs(f_horizon(i)-wrapTo2Pi(kep_sat(i,6))) > tol
      f_error_sat_Car(i)=2*pi-abs(f_horizon(i)-wrapTo2Pi(kep_sat(i,6)));
   else
      f_error_sat_Car(i)=abs(f_horizon(i)-wrapTo2Pi(kep_sat(i,6)));
   end 
end 


%% COMPARING KEPLERIAN ELEMENTS FOUND BY HORIZON AND kep_Gauss_sat
% error between real and numerical keplerian elements (Gauss)
a_error_sat_Gauss=abs((a_horizon-kep_Gauss_sat(:,1))/kep0_sat(1));
e_error_sat_Gauss=abs(e_horizon-kep_Gauss_sat(:,2));
i_error_sat_Gauss=abs(unwrap(i_horizon)-kep_Gauss_sat(:,3))/2*pi;
RAAN_error_sat_Gauss=abs(unwrap(RAAN_horizon)-kep_Gauss_sat(:,4))/2*pi;
omega_error_sat_Gauss=abs(unwrap(omega_horizon)-kep_Gauss_sat(:,5))/2*pi;

f_error_sat_Gauss=zeros(n_it,1);
for i=1:n_it
   if abs(f_horizon(i)-wrapTo2Pi(kep_Gauss_sat(i,6))) > tol
     f_error_sat_Gauss(i)=2*pi-abs(f_horizon(i)-wrapTo2Pi(kep_Gauss_sat(i,6)));
   else
     f_error_sat_Gauss(i)=abs(f_horizon(i)-wrapTo2Pi(kep_Gauss_sat(i,6)));
   end 
end 

tspan=tspan/(3600*24);

%% PLOTTING KEPLERIAN ELEMENTS ERRORS FOUND BY HORIZON AND kep_sat
figure (20)
semilogy(tspan,a_error_sat_Car)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|a_H_o_r_i_z_o_n - a_C_a_r|/a_0 [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (21)
semilogy(tspan,e_error_sat_Car)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|e_H_o_r_i_z_o_n - e_C_a_r| [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (22)
semilogy(tspan,i_error_sat_Car)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|i_H_o_r_i_z_o_n - i_C_a_r|/2\pi [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (23)
semilogy(tspan,RAAN_error_sat_Car)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|\Omega_H_o_r_i_z_o_n - \Omega_C_a_r|/2\pi [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (24)
semilogy(tspan,omega_error_sat_Car)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|\omega_H_o_r_i_z_o_n - \omega_C_a_r|/2\pi [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (25)
semilogy(tspan,rad2deg(f_error_sat_Car))
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|{\itf}_H_o_r_i_z_o_n - {\itf}_C_a_r| [{\itdeg}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

%% PLOTTING KEPLERIAN ELEMENTS ERRORS FOUND BY HORIZON AND kep_Gauss_sat
figure (26)
semilogy(tspan,a_error_sat_Gauss)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|a_H_o_r_i_z_o_n - a_G_a_u_s_s|/a_0 [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (27)
semilogy(tspan,e_error_sat_Gauss)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|e_H_o_r_i_z_o_n - e_G_a_u_s_s| [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (28)
semilogy(tspan,i_error_sat_Gauss)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|i_H_o_r_i_z_o_n - i_G_a_u_s_s|/2\pi [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (29)
semilogy(tspan,RAAN_error_sat_Gauss)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|\Omega_H_o_r_i_z_o_n - \Omega_G_a_u_s_s|/2\pi [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (30)
semilogy(tspan,omega_error_sat_Gauss)
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|\omega_H_o_r_i_z_o_n - \omega_G_a_u_s_s|/2\pi [{\it-}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (31)
semilogy(tspan,rad2deg(f_error_sat_Gauss))
hold on 
grid
xlabel( 'time [{\itdays}]');
ylabel('|{\itf}_H_o_r_i_z_o_n - {\itf}_G_a_u_s_s| [{\itdeg}] ');
hold off
set(gca, 'XLim', [0, tspan(end)])

%% PLOTTING THE KEPLERIAN ELEMENTS
figure (32) 
plot(tspan,kep_Gauss_sat(:,1),'g')
hold on
plot(tspan,kep_sat(:,1),'r')
plot(tspan,a_horizon,'b')
grid
xlabel( 'time [{\itdays}]');
ylabel('a [{\itkm}]');
legend('Gauss equations','Cartesian','Real data');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (33) 
plot(tspan,kep_Gauss_sat(:,2),'g')
hold on
plot(tspan,kep_sat(:,2),'r')
plot(tspan,e_horizon,'b')
grid
xlabel( 'time [{\itdays}]');
ylabel('e [{\it-}]');
legend('Gauss equations','Cartesian','Real data');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (34)
plot(tspan,rad2deg(kep_Gauss_sat(:,3)),'g')
hold on
plot(tspan,rad2deg(unwrap(kep_sat(:,3))),'r')
plot(tspan,rad2deg(unwrap(i_horizon)),'b')
grid
xlabel( 'time [{\itdays}]');
ylabel('i [{\itdeg}]');
legend('Gauss equations','Cartesian','Real data');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (35) 
plot(tspan,rad2deg(kep_Gauss_sat(:,4)),'g')
hold on
plot(tspan,rad2deg(unwrap(kep_sat(:,4))),'r')
plot(tspan,rad2deg(unwrap(RAAN_horizon)),'b')
grid
xlabel( 'time [{\itdays}]');
ylabel('\Omega [{\itdeg}]');
legend('Gauss equations','Cartesian','Real data');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (36) 
plot(tspan,rad2deg(kep_Gauss_sat(:,5)),'g')
hold on
plot(tspan,rad2deg(unwrap(kep_sat(:,5))),'r')
plot(tspan,rad2deg(unwrap(omega_horizon)),'b')
grid
xlabel( 'time [{\itdays}]');
ylabel('\omega [{\itdeg}]');
legend('Gauss equations','Cartesian','Real data');
hold off
set(gca, 'XLim', [0, tspan(end)])

figure (37) 
plot(tspan,rad2deg(kep_Gauss_sat(:,6)),'g')
hold on
plot(tspan,rad2deg(unwrap(kep_sat(:,6))),'r')
plot(tspan,rad2deg(unwrap(f_horizon)),'b')
grid
xlabel( 'time [{\itdays}]');
ylabel('{\itf} [{\itdeg}]');
legend('Gauss equations','Cartesian','Real data');
set(gca, 'XLim', [0, tspan(end)])
hold off

tspan=linspace(0,dt,n_it);

%% PLOT ORBIT
 figure (38)
%Perturbed orbit from horizons data
 for i=1:n_it
[R_horizon(i,:),V_horizon(i,:)] = kep2car(kep_horizon(i,:),muE);
end
plot3(R_horizon(:,1),R_horizon(:,2),R_horizon(:,3),'r','LineWidth',2)
hold on
%Perturbed orbit (Cartesian)
plot3 (R_sat(:,1),R_sat(:,2),R_sat(:,3),'g');
hold on
%Perturbed orbit (Gauss)
 for i=1:length(tspan)
     [R_Gauss(i,:),V_Gauss(i,:)] = kep2car(kep_Gauss_sat(i,:),muE);
end
plot3(R_Gauss(:,1),R_Gauss(:,2),R_Gauss(:,3),'b')
xlabel ('x [{​​​​\itkm}​​​​]')
ylabel ('y [{​​​​\itkm}​​​​]')
zlabel ('z [{​​​​\itkm}​​​​]')
legend('Real data','Cartesian','Gauss equations');
hold on 
grid on
[x, y, z] = sphere (180);
X1 = x * R_E;
Y1 = y * R_E;
Z1 = z * R_E;
Earth = surf(X1,Y1,Z1);
cdata = imread('EarthTexture.jpg'); 
set(Earth, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
axis equal
title ('Real orbit vs modelled orbit')
set(gca,'GridColor',[1 1 1])
set(gca,'Color','k')
set(gca,'GridColor',[1 1 1])
view (45,45)
legend('Real orbit','Modelled orbit (Cartesian)','Modelled orbit (Gauss)' )
legend('boxoff')
