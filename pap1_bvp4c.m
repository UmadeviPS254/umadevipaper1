function pap1_bvp4c
clc
close all
%%%%%%%% 
Pr=6.7; 
Sc=1; 
Fr=0.2;
E=0.5;
n=0.5; 
beta=0.5; 
sigma=0.5; 
BiT=0.3; 
BiC=0.3; 
BiN=0.3; 
Hg=0.1; 
lambda=0.5; 
delta=0.5; 
LeB=1; 
Nb=0.5; 
Nt=0.2;   
Sr=0.2; 
Df=0.2; 
Rd=0.5; 
S=0.3; 
%%%%%%%%%%%%%%%%-----
s3 = 0;  %initial guess for f''(0)
s5 = 0;  %initial guess for g'(0)
s6 = 0;  %initial guess for theta(0)
s8 = 0;  %initial guess for chi(0)
s10 = 0; %initial guess for psi(0)   
%%%%%%----

% beta1=[0, 0.2, 0.4, 0.6, 0.8];
lambda1=[0.5, 0.75, 1, 1.25, 1.5];
% Rd1 = [0.4, 0.8, 1.2, 1.6, 2.0];
% Nb1 =[0.5, 0.75, 1, 1.25, 1.5]; 
% Nt1 =[0, 0.05, 0.1, 0.15, 0.2]; 
% E1 = [0, 1, 2, 3, 4];
% sigma1= [0.5, 1.25, 2, 2.75, 3.5];
% S1 =[0.3, 0.4, 0.5, 0.6, 0.7];

% Fr1=[0.1 0.7 1.3 1.9 2.5];
% BiT1=[0.1, 0.25, 0.5, 0.75, 1.0];
% BiC1=[0.1, 0.25, 0.5, 0.75, 1.0];
% BiN1=[0.1, 0.25, 0.5, 0.75, 1.0];
% Hg1 = [-0.15, -0.1, 0, 0.1, 0.15];
% delta1 = [0, 1.5, 3, 4.5, 6];
% n1=[-1, -0.5, 0, 0.5, 1];
% Sr1 = [0, 0.15, 0.3, 0.45, 0.6];
% Df1 = [0, 0.05, 0.1, 0.15, 0.2];



for j=1:5

% beta=beta1(j);
lambda=lambda1(j);
% Rd = Rd1(j);  
% Nb = Nb1(j);
% Nt = Nt1(j);
% E = E1(j);
% sigma = sigma1(j);
% S = S1(j);

% Fr=Fr1(j);
% BiT=BiT1(j);
% BiC=BiC1(j);
% BiN=BiN1(j);
% Hg = Hg1(j);
% delta = delta1(j);
% n = n1(j);
% Sr = Sr1(j);
% Df = Df1(j);


guess = [0 1+S*s3 s3 0 s5 s6 -BiT*(1-s6) s8 -BiC*(1-s8) s10 -BiN*(1-s10)];

options = bvpset('AbsTol', 1e-6, 'RelTol',1e-6);
infinity = 1;
maxinfinity = 10;

solinit = bvpinit(linspace(0,infinity,200), guess);
sol= bvp4c(@ddc_cone_ode,@ddc_cone_bc,solinit,options);

     for Bnew = infinity+1:maxinfinity

     solinit = bvpxtend(sol, Bnew);

     sol= bvp4c(@ddc_cone_ode, @ddc_cone_bc,solinit,options);

     end

    
      x   = sol.x;
      y1  = sol.y(1,:);
      y2  = sol.y(2,:);
      y3  = sol.y(3,:);
      y4  = sol.y(4,:);
      y5  = sol.y(5,:);
      y6  = sol.y(6,:);
      y7  = sol.y(7,:);
      y8  = sol.y(8,:);
      y9  = sol.y(9,:);
      y10 = sol.y(10,:);
      y11 = sol.y(11,:);

      xChi{j} = x;
      chiSol{j} = y8;


%.....velocity profile---
figure(1)
plot(x,y2,'linewidth',1.5);
set(gca,'linewidth',1.5);
xlabel('$\eta$','Interpreter','latex','Fontsize',20,'FontWeight','bold')
ylabel('$f''(\eta)$','Interpreter','latex','Fontsize',20,'FontWeight','bold')
hold on

% ---temperature profile----
figure(2)
plot(x, y6, 'linewidth', 1.5);
set(gca,'linewidth',1.5);
xlabel('$\eta$','Interpreter','latex', 'Fontsize', 20, 'FontWeight', 'bold');
ylabel('$\theta(\eta)$','Interpreter','latex','Fontsize', 20, 'FontWeight', 'bold');
hold on

%.....dissolved-solute concentration profile (chi)---
figure(3)
plot(x, y8,'linewidth',1.5);
set(gca,'linewidth',1.5);
xlabel('$\eta$','Interpreter','latex','Fontsize',20,'FontWeight','bold');
ylabel('$\chi(\eta)$','Interpreter','latex','Fontsize',20,'FontWeight','bold');
hold on


%.....nanoparticle volume-fraction profile (psi)  
figure(4)
plot(x, y10,'linewidth',1.5);
set(gca,'linewidth',1.5);
xlabel('$\eta$','Interpreter','latex','Fontsize',20,'FontWeight','bold');
ylabel('$\psi(\eta)$','Interpreter','latex','Fontsize',20,'FontWeight','bold');
hold on




% fprintf('beta=%.2f\n', beta);
% fprintf('f''=%f \n',-y3(1));
% fprintf('g''=%f \n', -y5(1));
% fprintf('theta''=%f \n', -(1 + (4/3)*Rd)*y7(1));
% fprintf('chi''=%f \n', -y9(1));
% fprintf('psi''=%f \n', -y11(1));
% fprintf('%f, %f, %f\n', -y3(1), -(1 + (4/3)*Rd)*y7(1), -y9(1));%eta=0 at 1 (checkin workspace)
% 
hold on
end
% 
% 



function dydx = ddc_cone_ode (~,y)
A1 = (1/Pr)*(1+((4/3)*Rd));
B1 = Df;           % Dufour: driven by solute chi''
A2 = Sr;           % Soret : driven by theta''
B2 = (1/Sc);
% Energy: (1/Pr)(1+4Rd/3)theta'' + f*theta' + Nb*psi'*theta' + Nt*(theta')^2 + Df*chi'' + Hg*theta = 0
R1 = -(y(1)*y(7) + Nb*y(11)*y(7) + Nt*y(7)^2 + Hg*y(6));
react= sigma*(1+delta*y(6))^n*exp(-E/(1+delta*y(6)))*y(8);  % reaction acts on chi
R2 = -(y(1)*y(9)-react);
sol = [A1 B1; A2 B2] \ [R1; R2];
thetapp = sol(1);   % theta''
chipp   = sol(2);   % chi''
psipp = -LeB*y(1)*y(11) - (Nt/Nb)*thetapp;


dydx = [y(2)
    y(3)
    -y(1)*y(3)+lambda*y(2)-2*beta*y(4)+(1+Fr)*(y(2)^2)
    y(5)
    -y(1)*y(5)+y(2)*y(4)+2*beta*y(2)+lambda*y(4)+Fr*(y(4)^2)
    y(7)
    thetapp
    y(9)
    chipp
    y(11)
    psipp];
end

%Boundary conditions at eta=0 & infinity (11 conditions for 11 unknowns)
function res = ddc_cone_bc (ya, yb)
    res = [ya(1)
          ya(2)-(1+S*ya(3))
          ya(4)
          ya(7)+(BiT*(1-ya(6)))     % thermal convective
          ya(9)+(BiC*(1-ya(8)))     % chi  convective
          ya(11)+(BiN*(1-ya(10)))   % psi  convective  <-- NEW
          % For a prescribed wall value instead, replace the line above with:
          % ya(10)-1
          yb(2)
          yb(4)
          yb(6)
          yb(8)
          yb(10)];                  % psi(inf)=0        <-- NEW
end
%-------
end
