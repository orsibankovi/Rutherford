clear all, close all;

k=8.987551e+9; %Coulomb allando [Nm^2/C^2]
q1=2*1.6021766e-19; %alpha (helium) [C]
q2=79*1.6021766e-19; %arany [C]
m=6.624424e-27; %[kg]
km=k/m;
D=4.065e-10;

v0=5.2e+5; %kezdeti sebesseg [m/s]

dt=1e-18;
figure(3);
for k = -3:3
    rc=5.29177e-11;
    th=0:pi/60:2*pi;
    xunit=rc*cos(th);
    yunit=rc*sin(th)+k*D;
    hold on;
    plot(xunit, yunit)
end
hold on
    
for j=-150:150
    vx=v0;
    vy=0;
    rx=4e-10;
    ry=j*1e-11;
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
        
        ax=km*q1*q2*x*(1/r^3+1/r1^3+1/r2^3+1/r3^3+1/r4^3+1/r5^3+1/r6^3);
        ay=km*q1*q2*(y/r^3+(y-D)/r1^3+(y-2*D)/r2^3+(y+D)/r3^3+(y+2*D)/r4^3+(y+3*D)/r5^3+(y-3*D)/r6^3);

        vx=vx+ax*dt;
        vy=vy+ay*dt;

        vx_=[vx_ vx];
        vy_=[vy_ vy];

        x=x+vx*dt;
        y=y+vy*dt;

        x_=[x_ x];
        y_=[y_ y];
        
    end
    plot(x_, y_)
    hold on;
end


pbaspect([1 1 1])
xlim([-7e-10, 7e-10]);
ylim([-7e-10, 7e-10]);
xlabel('x coordinates [m]');
ylabel('y coordinates [m]');
title('Orbitals of alpha particles');