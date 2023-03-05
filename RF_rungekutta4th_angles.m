clear all, close all;

k=8.987551e+9; %Coulomb allando [Nm^2/C^2]
q1=2*1.6021766e-19; %alpha (helium) [C]
q2=79*1.6021766e-19; %arany [C]
m=6.624424e-27; %[kg]
km=k/m;
D=4.065e-10;
theta_=[];
y0=[];
th030=0;
th3060=0;
th6090=0;
th90120=0;
th120150=0;
th150180=0;


v0=5.2e+5; %kezdeti sebesseg [m/s]

dt=1e-18;
figure(4);
for k = -3:3
    rc=5.29177e-11;
    th=0:pi/60:2*pi;
    xunit=rc*cos(th);
    yunit=rc*sin(th)+k*D;
    hold on;
    plot(xunit, yunit)
end
hold on

for j=-180:180
    vx=v0;
    vy=0;
    rx=4e-10;
    ry=j*8.1e-12;
    x=-rx;
    y=ry;

    x_=[x];
    y_=[y];
    vx_=[vx];
    vy_=[vy];
    
    for i=0:10000
       
        r=sqrt(x^2+y^2);
        r1=sqrt(x^2+(y-D)^2);
        r2=sqrt(x^2+(y-2*D)^2);
        r3=sqrt(x^2+(y+D)^2);
        r4=sqrt(x^2+(y+2*D)^2);
        r5=sqrt(x^2+(y+3*D)^2);
        r6=sqrt(x^2+(y-3*D)^2);
        
        v_x(1)=dt*km*q1*q2*x*(1/r^3+1/r1^3+1/r2^3+1/r3^3+1/r4^3+1/r5^3+1/r6^3);
        v_y(1)=dt*km*q1*q2*(y/r^3+(y-D)/r1^3+(y-2*D)/r2^3+(y+D)/r3^3+(y+2*D)/r4^3+(y+3*D)/r5^3+(y-3*D)/r6^3);
        x_k(1)=dt*vx;
        y_k(1)=dt*vy;
        
        v_x(2)=dt*km*q1*q2*(x+x_k(1)/2)*(1/r^3+1/r1^3+1/r2^3+1/r3^3+1/r4^3+1/r5^3+1/r6^3);
        v_y(2)=dt*km*q1*q2*((y+y_k(1)/2)/r^3+((y+y_k(1)/2)-D)/r1^3+((y+y_k(1)/2)-2*D)/r2^3+((y+y_k(1)/2)+D)/r3^3+((y+y_k(1)/2)+2*D)/r4^3+((y+y_k(1)/2)+3*D)/r5^3+((y+y_k(1)/2)-3*D)/r6^3);
        x_k(2)=dt*(vx+v_x(1)/2);
        y_k(2)=dt*(vy+v_y(1)/2);
        
        v_x(3)=dt*km*q1*q2*(x+x_k(2)/2)*(1/r^3+1/r1^3+1/r2^3+1/r3^3+1/r4^3+1/r5^3+1/r6^3);
        v_y(3)=dt*km*q1*q2*((y+y_k(2)/2)/r^3+((y+y_k(2)/2)-D)/r1^3+((y+y_k(2)/2)-2*D)/r2^3+((y+y_k(2)/2)+D)/r3^3+((y+y_k(2)/2)+2*D)/r4^3+((y+y_k(2)/2)+3*D)/r5^3+((y+y_k(2)/2)-3*D)/r6^3);
        x_k(3)=dt*(vx+v_x(2)/2);
        y_k(3)=dt*(vy+v_y(2)/2);
        
        v_x(4)=dt*km*q1*q2*(x+x_k(3))*(1/r^3+1/r1^3+1/r2^3+1/r3^3+1/r4^3+1/r5^3+1/r6^3);
        v_y(4)=dt*km*q1*q2*((y+y_k(3))/r^3+((y+y_k(3))-D)/r1^3+((y+y_k(3))-2*D)/r2^3+((y+y_k(3))+D)/r3^3+((y+y_k(3))+2*D)/r4^3+((y+y_k(3))+3*D)/r5^3+((y+y_k(3))-3*D)/r6^3);
        x_k(4)=dt*(vx+v_x(3));
        y_k(4)=dt*(vy+v_y(3));
        
        vx=vx+(v_x(1)+2*v_x(2)+2*v_x(3)+v_x(4))/6;
        vy=vy+(v_y(1)+2*v_y(2)+2*v_y(3)+v_y(4))/6;
        
        vx_=[vx_ vx];
        vy_=[vy_ vy];
        
        x=x+(x_k(1)+2*x_k(2)+2*x_k(3)+x_k(4))/6;
        y=y+(y_k(1)+2*y_k(2)+2*y_k(3)+y_k(4))/6;
        
        x_=[x_ x];
        y_=[y_ y];

    end
    
    theta=atan2((y_(end)-y_(1)), (x_(end)-x_(1)))*180/pi;
    theta=abs(theta);
    
    if theta>=0 && theta<30
        th030=th030+1;
    end
    if theta>=30 && theta<60
        th3060=th3060+1;
    end
    if theta>=60 && theta<=90
        th6090=th6090+1;
    end
    if theta>90 && theta<=120
        th90120=th90120+1;
    end
    if theta>120 && theta<=150
        th120150=th120150+1;
    end
    if theta>150 && theta<=180
        th150180=th150180+1;
    end
    
    theta_=[theta_ theta];
    y0=[y0 ry];
    plot(x_, y_)
    hold on;
end
pbaspect([1 1 1])
xlim([-7e-10, 7e-10]);
ylim([-7e-10, 7e-10]);

xlabel('x coordinates [m]');
ylabel('y coordinates [m]');
title('Orbitals of alpha particles');

figure(5);
plot(y0, theta_)
xlim([-1.5e-9, 1.5e-9]);
ylim([0, 180]);
xlabel('y coordinates [m]');
ylabel('Angle [degree]');

figure(6);
X = categorical({'0°<theta<30°','30°<theta<60°','60°<theta<90°','90°<theta<120°', '120°<theta<150°','150°<theta<180°'});
X = reordercats(X,{'0°<theta<30°','30°<theta<60°','60°<theta<90°','90°<theta<120°', '120°<theta<150°','150°<theta<180°'});
Y = [th030 th3060 th6090 th90120 th120150 th150180];
bar(X,Y);
ylabel('Number of alpha particles');