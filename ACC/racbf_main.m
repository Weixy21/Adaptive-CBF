clc;
x0 = [0 20 100];%(x, v, z), initial state
g = 9.81; m = 1650; f0 = 0.1; f1 = 5; f2 = 0.25;%dynamics param.
vd = 24;%desired speed
eps = 10;%param. for desired speed with CLF
ca = 0.4; cd = 0.23; %max/min acceleration coe. (tentative)
psc = 1;
result1 = zeros(1,5); data1 = zeros(1,8); y_history = []; u_min = []; %for recording
c = 10; %minimum dis.
set(0,'DefaultTextInterpreter','latex')
warning('off');
global u v0;
v0 = 13.89; % speed of the preceding vehicle

fea = 1;  % feasibility of the QP
te = 300; % terminal time

dd = 6; % desired value for r/r1
p = [0.1,2,1,1]; %p1, q1, p2, q2; all are constants in the racbf
y = [dd,0]; % (r,r2) state for the auxiliary dynamics

for i = 1:te
    i

% CLF for desired speed
Fr = f0 + f1*x0(2) + f2*x0(2)^2; %Fr = 0;
phi0 = -2*(x0(2) - vd)*Fr/m + eps*(x0(2) - vd)^2;
phi1 = 2*(x0(2) - vd)/m;

b = x0(3) - c - y(1); %racbf
b_dot = v0 - x0(2) - y(2);
psi_1 = b_dot + p(1)*b^(p(2));
b_acbf = Fr/m + p(1)*b^(p(2))*p(2)*b_dot/b + p(3)*psi_1^(p(4)); %racbf constraint
A_acbf = [1/m, 0, 1, 0]; 
C2 = b_dot + p(1)*b^(p(2));b = x0(3) - c; C1 = b; 
data1(i,1:6) = [b, b_acbf, phi1, 1/m, -phi0, b_acbf];
u_min = [u_min;[b,-cd*g]];
    lb = 1; %ra in (51) in the paper for considering the noise
    b = y(1) - lb;  b_dot = y(2);
    b_cbf = 2*b_dot + b; %cbf constraint (51) in the paper
    A_cbf = [0,0,-1,0];
    k1 = 1;
    y2d = -k1*(y(1) - dd);
    V = (y(2) - y2d)^2;
    b_clf = -eps*V - 2*k1*y(2)*(y(2) - y2d); % CLF (52) in the paper for stabilizing r
    A_clf = [0,0,2*(y(2) - y2d),-1];
    QQ = 20;
    %(u, \delta for desired speed, \nu, \delta for stabilizing r)
    A = [phi1 -1 0 0;A_acbf;A_cbf; A_clf; 1 0 0 0;-1 0 0 0];%;LgBf 0;1 0
    b = [-phi0;b_acbf;b_cbf;b_clf; ca*m*g; cd*m*g];%;LfBf;ca*m*g
    H = [2/(m^2) 0 0 0; 0 2*psc 0 0;0 0 QQ 0; 0 0 0 QQ*psc];
    F = [-2*Fr/(m^2); 0; 0; 0];
if(fea) % if the QP is feasible
  %  tic
    options = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','off');
    [u1,fval,exitflag,output,lambda] = ...
       quadprog(H,F,A,b,[],[],[],[],[],options);
 %  toc
    t=[0 0.1];
    if(numel(u1) == 0)
        break;
        fea = 0; te = 100; infea = i
    else
        u = u1;
    end
end
if(i > te)% infeasibile, end the loop
    break;
end
data1(i,7:8) = [u(1)/m, u(2)];
[tt,xx]=ode45('acc_racbf',t,[x0, y]);% integrate the dynamics
x0 = [xx(end, 1) xx(end, 2) xx(end, 3)]; %update state
y = xx(end, 4:5);%update the auxiliary dynamics state
y_history = [y_history;[0.1*i,y]];
result1(i,:) = [0.1*i xx(end,2) u(1)/m C1 C2];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot for adaptive method
 plot( data1(:,1), data1(:,7),'r','linewidth',1.6)
axis([0 90 -5 5]);
xlabel('$b(\textbf{\textit x}(t))$','fontsize',15)
ylabel('$u(t)/M$','fontsize',15)
grid on
hold on
plot([0,90],[ca*g, ca*g],'k:',u_min(:,1),u_min(:,2),'k:','linewidth',1.6)
text(20,0,'$c_ag$', 'fontsize',15)
text(20,0,'$-c_dg$', 'fontsize',15)
annotation('arrow',[0.2,0.2],[0.4,0.2])


figure(2)
subplot(3, 1, 1)
plot(result1(:, 1),result1(:, 3), 'b-',[0,30],[0,0], 'k--', 'linewidth',1.5)
xlabel('$t/s$','fontsize',15)
ylabel('$b(u_i(t)/M)$','fontsize',15)
lg = legend('Linear','Quadratic');
set(lg,'box','off')
set(lg,'fontsize',12)
% axis([0 30 -10 100]);
subplot(3, 1, 2)
plot(result1(:, 1),y_history(:, 2), 'b-',[0,30],[0,0], 'k--', 'linewidth',1.5)
xlabel('$t/s$','fontsize',15)
ylabel('$p_1(t)$','fontsize',15)
lg = legend('Linear','Quadratic');
set(lg,'box','off')
set(lg,'fontsize',12)
subplot(3, 1, 3)
plot(result1(:, 1),y_history(:, 3), 'b-',[0,30],[0,0], 'k--', 'linewidth',1.5)
xlabel('$t/s$','fontsize',15)
ylabel('$p_2(t)$','fontsize',15)
lg = legend('Linear','Quadratic');
set(lg,'box','off')
set(lg,'fontsize',12)

 figure(3)
subplot(2, 1, 1)
plot(result1(:, 1),result1(:, 4),'g',[0,30],[0,0], 'k--', 'linewidth',1.5)
xlabel('$t/s$','fontsize',15)
ylabel('$b(\textbf{\textit x}(t))$','fontsize',15)
lg = legend('Noise $\times 1$','Noise $\times 2$');
set(lg,'box','off')
set(lg,'Interpreter','latex')
set(lg,'fontsize',12)
axis([0 30 -1 10]);
subplot(2, 1, 2)
plot(result1(:, 1),result1(:, 5),'g',[0,30],[0,0], 'k--', 'linewidth',1.5)
xlabel('$t/s$','fontsize',15)
ylabel('$\psi_1(\textbf{\textit x}(t))$','fontsize',15)
lg = legend('Noise $\times 1$','Noise $\times 2$');
set(lg,'box','off')
set(lg,'Interpreter','latex')
set(lg,'fontsize',12)
axis([0 30 -1 10]);
% axis([0 30 -10 160]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end plot





 