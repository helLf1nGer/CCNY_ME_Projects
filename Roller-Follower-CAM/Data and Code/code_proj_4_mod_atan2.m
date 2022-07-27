clear
close all

%Initial data
rb = 2*10^-2;
rr = 1.5*10^-2;
L = 11*10^-2;
X = 7*10^-2;
Y = 8*10^-2;
%design constraints
Theta = linspace(0, 2*pi, 450);
beta = 110*pi/180;
startAlpha = -31*pi/180; %starting point in rad
endAlpha = -11*pi/180; %endpoint in rad
delAlpha = abs(startAlpha - endAlpha);

%Functions
flift = delAlpha/2*(1-cos(pi*Theta/beta));
fliftDot = sin (pi*Theta/beta)*pi/beta*delAlpha/2;
freturn = delAlpha/2*(1+cos(pi*(Theta-2*pi+beta)/beta));
freturnDot = -delAlpha/2*pi/beta*sin(pi*(Theta-2*pi+beta)/beta);

%initial follower angle 
alpha0 = -atan2(Y,X)+acos((X^2+Y^2+L^2-(rr+rb)^2)/(2*L*sqrt(X^2+Y^2)));
%%%%follower angle%%%%
for i = 1:1:length(Theta)
    if Theta(i)<beta
        alpha(i) = alpha0 + flift(i);%gate angle
        Chi (i) = atan2((Y+L*(1+fliftDot(i))*sin(alpha(i))),(X-L*(1+fliftDot(i))*cos(alpha(i))));%contact angle on wheel
    elseif Theta(i)>(2*pi-beta)
        alpha(i) = alpha0 + freturn(i);
        Chi (i)= atan2((Y+L*(1+freturnDot(i))*sin(alpha(i))),(X-L*(1+freturnDot(i))*cos(alpha(i))));
    else
        alpha(i) = alpha0 + delAlpha;
        Chi (i)= atan2((Y+L*sin(alpha(i))),(X-L*cos(alpha(i))));
    end
end

%contact point
Xc = X - L*cos(alpha)-rr*cos(Chi);
Yc = Y + L*sin(alpha)-rr*sin(Chi);
rc = sqrt(Xc.^2+Yc.^2);
Gammac = atan2(Yc,Xc);
Xp = rc.*cos(Theta-Gammac); %X-coord of CAM
Yp = rc.*sin(Theta-Gammac); %Y-Coord of CAM

%pressure angle
phi = pi/2-alpha-Chi;

%Motion
figure
xlabel('x-Position, m')
ylabel('y-Position, m')
axis equal; grid minor;
for i = 1:1:length(Theta)
    clf
    grid minor
    hold on
    axis equal
    xlim ([-0.1 0.1])
    ylim ([-0.1 0.1])
    %cam center
    plot (0,0, 'o','DisplayName', 'Pin CAM')
    %bar pin
    plot (X,Y, 'o','DisplayName', 'Pin Bar')
    %point of contact
    plot (Xc(i),Yc(i),'o','DisplayName', 'Point of Contact')
    
    %Cam shape
    plot(Xp,Yp,'-.','Color','blue','DisplayName', 'CAM Shape')
    %CAM redrawing-rotation
    plot (Xp*cos(Theta(i)-31*pi/180+3/2*pi)-Yp*sin(Theta(i)-31*pi/180+3/2*pi),...
    Xp*sin(Theta(i)-31*pi/180+3/2*pi)+Yp*cos(Theta(i)-31*pi/180+3/2*pi),'Color',[0, 0.4470, 0.7410],'LineWidth', 3,...
    'DisplayName', 'CAM') 
     
    %cam base circle
    plot (rb*cos(Theta), rb*sin(Theta),'DisplayName', 'Base circle','Color','red')
 
    %bar
    line ([X, Xc(i)], [Y, Yc(i)+rr],'LineWidth', 5,'Color','green','DisplayName', 'gate bar')
    %follower roller
    plot(Xc(i)+rr*cos(Theta),Yc(i)+rr+rr*sin(Theta),'Color', [0.4940, 0.1840, 0.5560],'DisplayName', 'Roller follower')
  
    pause(0.001)
    hold off
end
legend('Location','northeast')

%%%IMPORT DATA%%%
EXP_data = readtable('Tracker_data_Alpha_Theta.xlsx');
EXP_data = EXP_data {:,:};
EXP_alpha = EXP_data (:,3)-32;
EXP_Theta = EXP_data (:,6)+180;
SW_data = readtable('SW_Aplha_vs_Theta.csv');
SW_data = SW_data {:,:};
SW_alpha = SW_data (:,3);
SW_Theta = SW_data (:,1)+180;

%plots
figure
grid minor; hold on; 
plot (Theta*180/pi,alpha*180/pi, 'DisplayName', 'Analytical')
plot (EXP_Theta, EXP_alpha, '-.', 'DisplayName', 'Experimental')
plot (SW_Theta, SW_alpha, '--', 'DisplayName', 'Computational')
ylabel('α, °'); xlabel('θ, °')
hold off; legend('Location','northeast')

figure
grid minor; hold on; axis equal
plot (0,0,'o','DisplayName', 'Center') 
plot (Xp,Yp, 'DisplayName', 'CAM','Color','black')
plot (rb*cos(Theta), rb*sin(Theta),'-.','DisplayName', 'Base circle','Color','red')
ylabel('y, m'); xlabel('x, m'); legend('Location','northeast'); hold off



