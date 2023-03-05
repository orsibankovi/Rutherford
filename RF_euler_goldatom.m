clear all, close all;

k=8.987551e+9; %Coulomb allando [Nm^2/C^2]
q1=2*1.6021766e-19; %alpha (helium) [C]
q2=79*1.6021766e-19; %arany [C]
m=6.624424e-27; %[kg]
km=k/m;
D=4.065e-10;

v0=1e+6; %kezdeti sebesseg [m/s]

dt=1e-19;
figure(1);
rc=5.29177e-11;
th=0:pi/60:2*pi;
xunit=rc*cos(th);
yunit=rc*sin(th);
hold on;
plot(xunit, yunit)
hold on
for j=-50:50
    vx=v0;
    vy=0;
    rx=1e-10;
    ry=j*1e-12;
    x=-rx;
    y=ry;

    x_=[x];
    y_=[y];
    vx_=[vx];
    vy_=[vy];
    for i=0:10000
       
        r=sqrt(x^2+y^2);

        ax=km*q1*q2*x /r^3;
        ay=km*q1*q2*y /r^3;

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
 xlim([-10e-11, 6e-11]);
 ylim([-8e-11, 8e-11]);
 xlabel('x coordinates [m]');
ylabel('y coordinates [m]');
title('Orbitals of alpha particles');