%% Initialization

xmin=-1;
xmax=1;
N=100;
dt=0.001;
t=0;
tmax=0.47;
Re=50;
nsteps=tmax/dt;
dx=(xmax-xmin)/N;
x=xmin-dx:dx:xmax+dx;
u0=zeros(1,length(x));
a=zeros(1,N);b=zeros(1,N);
c=zeros(1,N);d=zeros(1,N);


%% initial conditions

u0(x<=0)=1.0;

u=u0;
unpl=u0;
s=dt/(Re*(dx^2));



%% Explicit method

for n=1:nsteps
    u(2)=1;
    u(N+2)=0;
    
    for i=2:N+2
        
        unpl(i)=u(i)-dt/(2*dx)*(u(i+1)^2/2-u(i-1)^2/2)+dt/(Re*dx^2)*(u(i-1)-2*u(i)+u(i+1));
        
    end
    
    u=unpl;
    plot(x,u);
    pause(dt);
    
end


%% Implicit method

for n=1:nsteps

    for i=2:N+2
        a(i)=-0.25*(dt/dx)*u(i-1)-0.5*s;
        b(i)=1+s;
        c(i)=0.25*(dt/dx)*u(i+1)-0.5*s;
        d(i)=0.5*s*u(i-1)+(1-s)*u(i)+0.5*s*u(i+1);       
    end
    
    unpl=TDMA(c(2:end),b(2:end),a(2:end),d(2:end));
    u(2:end-1)=unpl;
    u(2)=1;
    u(end-1)=0;  
    plot(x(2:end-1),u(2:end-1));
    pause(dt);
    
end

% %% TDMA Solver
% 
% function X = TDMA(a,b,c,d)
% n=length(d);
% c(1)=c(1)/b(1);
% d(1)=d(1)/b(1);
% % for i=2:n-1
% %     c(i)=c(i)/(b(i)-a(i)*c(i-1));
% % end
% 
% for i=2:n
%     c(i)=c(i)/(b(i)-a(i)*c(i-1));
%     d(i)=(d(i)-a(i)*d(i-1))/(b(i)-a(i)*c(i-1));
% end
% X(n)=d(n);
% for i=n-1:-1:1
%     X(i)=d(i)-X(i+1)*c(i);
% end


        
        
    
    