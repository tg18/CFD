% Ayush Bisen 21105025
% Initializing constants
%y'=-y --> equation
%y(0)=1
t0=0;
y0=1;
tEnd=10;
h=2.1;
N=floor((tEnd-t0)/h);



%% Initializing solutions
T=[t0:h:tEnd]';
Y=zeros(N+1,1);
Y(1)=y0;

%% Solving using Euler's Explicit Method
for i=1:N
    fi= -Y(i);
    Y(i+1)=Y(i) + h*fi;
end

%% Plotting the results 
plot(T,Y,'-r');hold on;
Yact=exp(-T);
plot(T,Yact,'-y'); hold on;


%% Solving using Euler's Implicit Method
for i=1:N
    t=T(i)+h;
    y=fsolve(@(y) y-Y(i)+h*(y), Y(i));
    
    T(i+1)=t;
    Y(i+1)=y;
end

%% Plotiing the results

plot(T,Y,'-g');hold on;

%% Solving using Crank Nicolson Method
for i=1:N
    t=T(i)+h;
    fi= -Y(i);
    y=fsolve(@(y) y-Y(i)-h/2*(fi+y), Y(i));
    
    T(i+1)=t;
    Y(i+1)=y;
end

%% Plotting the results

plot(T,Y,'-b');

%% Solving using wRK3 method
for i=1:N
    Y1(i)=Y(i);
    k1=derivative(T(i),Y1(i));
    Y1(i)= Y1(i)+(h/3)*k1;
    k1=(-5/9)*k1+derivative(T(i)+h/3,Y1(i));
    Y1(i)=Y1(i)+(15/16)*h*k1;
    k1= (-153/128)*k1 + derivative(T(i)+3/4*h,Y1(i));
    Y(i+1) = Y1(i) + 8/15*h*k1;
end

%% Plotting the results
plot(T,Y,'-k');



