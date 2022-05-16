%Ayush Bisen 21105025
%y''=-w^2*y=u'
%w=1
%u=y'
t0=0;
y0=1;
u0=0;
h=2.1;
tEnd=10;
N=floor((tEnd-t0)/h);

%% Initializing solutions
T=[t0:h:tEnd]';
Y=zeros(N+1,1);
U=zeros(N+1,1);
U(1)=u0;
Y(1)=y0;

%% Analytical Solution

Yact=cos(T);
plot(T,Yact,'-y');hold on;
%% Explicit Euler
for i=1:N
    Y(i+1)=Y(i)+h*U(i);
    U(i+1)=U(i)+h*(-Y(i));
end
plot(T,Y,'-r');hold on;

%% Implicit Euler
for i=1:N
    Y(i+1)=Y(i)+h*U(i);
    U(i+1)=U(i)+h*(-Y(i+1));
end
plot(T,Y,'-b');hold on;

%% Crank Nicolson
for i=1:N
    Y(i+1)=Y(i)+h*U(i);
    U(i+1)=U(i)+(h/2)*(-Y(i)-Y(i+1));
end
plot(T,Y,'-g');