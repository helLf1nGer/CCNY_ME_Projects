clear
close all
%experimental data import
E_data = readtable('T_v_Th_Exp_1.txt');
E_data = E_data {:,:};
E_data = sortrows(E_data);
Th_exp = E_data (:,2);
%LEGO parts dimensions
L=10^-3*(90); %Lego length !!!!NEED TO VERIFY!!!
r=10^-3*(40/2); %radius of the gear !!!NEED TO VERIFY!!!
h=7.38*10^-3; %height of the beam

mass_pin=0.2*10^-3;
mass_gear = 3.7*10^-3;
mass_beam = 3.7*10^-3;
mass1 = (mass_pin+mass_gear+mass_beam);
mass2 = (2*mass_pin+mass_beam);
mass3 = (mass_gear);
g=9.82; %acc of gravity
I_ab =((mass_beam)*L^2/12+(mass_gear+mass_pin)*r^2/2); %Combinet inertia 1 of the beam+gear
I_bd = ((2*mass_pin+mass_beam)*L^2/12);
W1 = mass1 * g; %weights
W2 = mass2 * g;
W3 = mass_gear * g;

%%%Initial conditions%%%
phi = 90*pi/180;
ksi = 0;
om1 = 0;
om2 = 0;
dt=10^-3;
c1 = 0.00005; %damping coefficient
c2 = 0.00005; %damping coefficient
t=2; %sample time
for i=1:1:t*1000
    A = [
        1 0 1 0 0 0 -mass1*L/2*cos(phi) 0;
        0 1 0 1 0 0 -mass1*L/2*sin(phi) 0;
        -L/2*cos(phi) -L/2*sin(phi) L/2*cos(phi) L/2*sin(phi) 0 0 -I_ab 0;
        0 0 -1 0 1 0 -mass2*L*cos(phi) -mass2*L/2*cos(ksi);
        0 0 0 -1 0 1 -mass2*L*sin(phi) -mass2*L/2*sin(ksi);
        0 0 L/2*cos(ksi) L/2*sin(ksi) L/2*cos(ksi) L/2*sin(ksi) 0 -I_bd;
        0 0 0 0 -1 0 -mass3*L*cos(phi) -mass3*L*cos(ksi);
        0 0 0 0 0 -1 -mass3*L*sin(phi) -mass3*L*sin(ksi)
        ];
    B = [
        -mass1*L/2*om1^2*sin(phi)
        mass1*L/2*om1^2*cos(phi)+W1
        c1*om1
        -mass2*L*om1^2*sin(phi)-mass2*L/2*om2^2*sin(ksi)
        mass2*L*om1^2*cos(phi)+mass2*L/2*om2^2*cos(ksi)+W2
        c2*om2
        -mass3*L*om1^2*sin(phi)-mass3*L*om2^2*sin(ksi)
        W3+mass3*L*om1^2*cos(phi)+mass3*L*om2^2*cos(ksi)
        ];
    X = inv(A)*B;
    x1 (i) = L/2*sin(phi); y1 (i) = -L/2*cos(phi); x2 (i) = L*sin(phi)+ L/2*sin(ksi) ; y2 (i) = -L*cos(phi)-L/2*cos(ksi); %x and y positions of the G1 G2
    x3 (i) = L*sin(phi)+ L*sin(ksi) ; y3 (i) = -L*cos(phi)-L*cos(ksi); xb (i) = L*sin(phi); yb (i) = -L*cos(phi); %x and y positions of the G3
    om1 = om1 + X(7)*dt; om2 = om2 + X(8)*dt;  phi = phi + om1*dt; ksi = ksi + om2*dt; %omega and phi, ksi step
    Theta_AB(i) = phi; Theta_BD(i) = ksi;  time(i) = i*dt;    omega_1(i)  = om1; omega_2(i) = om2; %iteration steps
    Ax(i) = X(1); Ay (i) = X(2); Bx(i) = X(3); By(i) = X(4); Dx(i) = X(5); Dy(i) = X(6); %output data
    alpha_1(i) = X(7); alpha_2(i) = X(8); %output data
    if x3(i)>0
        Theta(i) = atand(y3(i)/x3(i)); %angle between x axis and G3 !!!!NOT WORKING PROPERLY!!!
    elseif x3(i) == 0
        Theta(i) = -90;
    else
        Theta (i) = -90 -(90-abs(atand(y3(i)/x3(i))));
    end
end

%Plots
figure
xlabel('x-Position, mm')
ylabel('y-Position, mm')
title('Positions of G_1, G_2, G_3')
grid minor
for i=1:5:t*1000
    clf
    hold on
    plot (x1*1000, y1*1000, 'DisplayName', 'G_1','Color','b')
    plot (x2*1000, y2*1000, 'DisplayName', 'G_2','Color','r')
    plot (x3*1000, y3*1000, 'DisplayName', 'G_3','Color','g')
    plot (xb*1000, yb*1000, 'DisplayName', 'B','Color','y')
    line ([0, xb(i)*1000], [0, yb(i)*1000],'LineWidth', 2,'Color','black','DisplayName', 'Link AB')
    plot (x1(i)*1000,y1(i)*1000,'o','MarkerSize',50,'LineWidth', 2,'DisplayName', 'Gear C')
    line ([xb(i)*1000, x3(i)*1000], [yb(i)*1000, y3(i)*1000],'LineWidth', 2,'Color','black','DisplayName', 'Link BD')
    plot (x3(i)*1000,y3(i)*1000,'o','MarkerSize',50,'LineWidth', 2,'DisplayName', 'Gear D')    
    pause(0.0001)
    hold off
end
xlabel('x-Position, mm')
ylabel('y-Position, mm')
title('Positions of G_1, G_2, G_3')
grid minor
legend('Location','northwest')

%%%
figure
title('Data vs Time')

subplot (4,1,1)
grid minor
hold on
plot (time, Theta_AB*180/pi)

ylabel('phi, °')


subplot (4,1,2)
plot (time, omega_1, 'DisplayName', 'ω_1, rad/s')
ylabel('ω_1, rad/s')
grid minor
hold on

subplot (4,1,3)
plot (time, alpha_1, 'DisplayName', 'α_1, rad/s^{2}')
ylabel('α_1, rad/s^{2}')
grid minor
hold on

subplot (4,1,4)
plot (time, Ax, 'b',time, Ay, 'r')
ylabel('A_x A_y, N')
legend('A_x', 'A_y','Location','northwest')
grid minor
hold on


%%%_____%%%
figure
subplot (5,1,1)
grid minor
hold on
plot (time, Theta_BD*180/pi, 'DisplayName', 'ksi, °')
ylabel('ksi, °')

subplot (5,1,2)
plot (time, omega_2, 'DisplayName', 'ω_2, rad/s')
ylabel('ω_2, rad/s')
grid minor
hold on

subplot (5,1,3)
plot (time, alpha_2, 'DisplayName', 'α_2, rad/s^{2}')
ylabel('α_2, rad/s^{2}')
grid minor
hold on

subplot (5,1,4)
plot (time, Bx, 'b', time, By, 'r')
ylabel('B_x B_y, N')
grid minor
hold on
legend('B_x', 'B_y','Location','northwest')

subplot (5,1,5)
plot (time, Dx, 'b', time, Dy, 'r')
ylabel('D_x D_y, N')
grid minor
hold on
legend('D_x', 'D_y','Location','northwest')

xlabel('Time, s')

figure 
hold on
plot(E_data(:,1), E_data(:,2),'o','DisplayName','θ_{experiment}')
plot(time, Theta,'DisplayName','θ_{theory}')
legend('Location','northeast')
xlabel('Time, s')
ylabel('θ, °')