%% Initializing constants

Nx=64;Nxp2=Nx+2;
Ny=64;Nyp2=Ny+2;
dx=1/Nx;dy=1/Ny;
Lx=1;Ly=1;
x=0:dx:1;
y=0:dy:1;
Re=10000;
tstart=0;
tend=0.5;
dt=0.001;
nsteps=(tend-tstart)/dt;


%% Initializing phi

rng(0,'twister');
a=-1;
b=1;
phi = (b-a).*rand(Nxp2,Nyp2) + a;
w=zeros(Nxp2,Nyp2);
u=zeros(Nxp2,Nyp2);v=zeros(Nxp2,Nyp2);
phi_old=phi;

%% Initializing u,v,w
for i=2:Nx+1
    for j=2:Ny+1
        u(i,j)=(phi(i,j+1)-phi(i,j-1))/2*dy;
        v(i,j)=-(phi(i+1,j)-phi(i-1,j))/2*dx;
    end
end
for i=2:Nx+1
    for j=2:Ny+1
        w(i,j)=(v(i+1,j)-v(i-1,j))/2*dx-(u(i,j+1)-u(i,j-1))/2*dy;
    end
end

%% RK3 Method
for n=1:nsteps    
    for i=2:Nx+1
        for j=2:Ny+1
            k1(i,j)=dt*(-u(i,j)*(w(i+1,j)-w(i-1,j))/(2*dx)-v(i,j)*(w(i,j+1)-w(i,j-1))/(2*dy)+1/Re*((w(i+1,j)-2*w(i,j)+w(i-1,j))/(dx^2)+(w(i,j+1)-2*w(i,j)+w(i,j-1))/(dy^2)));
            w(i,j)=w(i,j)+(8/15)*k1(i,j);
            k2(i,j)=dt*(-u(i,j)*(w(i+1,j)-w(i-1,j))/(2*dx)-v(i,j)*(w(i,j+1)-w(i,j-1))/(2*dy)+1/Re*((w(i+1,j)-2*w(i,j)+w(i-1,j))/(dx^2)+(w(i,j+1)-2*w(i,j)+w(i,j-1))/(dy^2)));
            w(i,j)=w(i,j)+0.25*k1(i,j)+(5/12)*k2(i,j);
            k3(i,j)=dt*(-u(i,j)*(w(i+1,j)-w(i-1,j))/(2*dx)-v(i,j)*(w(i,j+1)-w(i,j-1))/(2*dy)+1/Re*((w(i+1,j)-2*w(i,j)+w(i-1,j))/(dx^2)+(w(i,j+1)-2*w(i,j)+w(i,j-1))/(dy^2)));
            w(i,j)=w(i,j)+0.25*k1(i,j)+0.75*k3(i,j);            
        end
    end
    
    %% Gauss Seidel
    
    while true
        for i=2:(Nx+1)
            for j=2:(Ny+1)
                phi(i,j)=0.25*(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+w(i,j)*(dx^2));
            end
        end
    
        l2=norm(phi-phi_old,2);
    
    
        if (l2>10^-5)
           phi_old=phi;
           continue;
        end
    
        if(l2<10^-5)
           break;
        end
    
    end
    for i=2:Nx+1
        for j=2:Ny+1
            u(i,j)=(phi(i,j+1)-phi(i,j-1))/2*dy;
            v(i,j)=-(phi(i+1,j)-phi(i-1,j))/2*dx;
        end
    end
    
    %% Updating Boundary conditions
     
    for j=2:(Ny+1)
        u(1,j)=u(Nx+1,j);
        u(Nxp2,j)=u(2,j);
        v(1,j)=v(Nx+1,j);
        v(Nxp2,j)=v(2,j);
        w(1,j)=w(Nx+1,j);
        w(Nxp2,j)=w(2,j);
        phi(1,j)=phi(Nx+1,j);
        phi(Nxp2,j)=phi(2,j);
    end
    for i=2:(Nx+1);
        u(i,1)=u(i,Nx+1);
        u(i,Nyp2)=u(i,2);
        v(i,1)=v(i,Nx+1);
        v(i,Nyp2)=v(i,2);
        w(i,1)=w(i,Nx+1);
        w(i,Nyp2)=w(i,2);
        phi(i,1)=phi(i,Nx+1);
        phi(i,Nyp2)=phi(i,2);
    end
    
    %% Plotting results
    
    [X,Y]=meshgrid(x,y);
    contourf(x(2:end),y(2:end),w(2:(end-1),2:(end-1)));
    hold on;
    quiver(x(2:end),y(2:end),u(2:(end-1),2:(end-1)),v((2:end-1),(2:end-1)),'b');
    hold off;
    colorbar;
    pause(dt);
end
    
          
    
            
            
    
    