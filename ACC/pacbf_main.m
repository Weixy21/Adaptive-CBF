clc;
x0 = [0 20 100]; %(x, v, z), initial state
g = 9.81; m = 1650; f0 = 0.1; f1 = 5; f2 = 0.25; %dynamics param.
vd = 24; %desired speed
eps = 10; %param. for desired speed with CLF
ca = 0.4; cd = 0.37; %max/min acceleration coe. (tentative)
psc = 1;
result1 = zeros(1,5); data1 = zeros(1,8); p_history = [];%for recording
c = 10; %minimum dis.
set(0,'DefaultTextInterpreter','latex')
warning('off');
global u v0;
v0 = 13.89; % speed of the preceding vehicle

without_p2 = 0; % take p2 as a penalty function or not
fea = 1; % feasibility of the QP
te = 600; % terminal time

p = [0.103,2,1,1]; %p1, q1, p2, q2; q1 and q2 are set as constants 2 and 1, resp., i.e., class k function 1 is quadratic, while class k function 2 is linear.
for i = 1:te
    i
    if(i >=69 && i <= 135)  % linearly increase the lower control bound from -0.37 to -0.2 starting at time 6.9s.
        cd = -0.37*i/66 + 69*0.37/66 + 0.37;
    end
    if (cd < 0.2)
        cd = 0.2; 
    end

% CLF for desired speed
Fr = f0 + f1*x0(2) + f2*x0(2)^2; 
phi0 = -2*(x0(2) - vd)*Fr/m + eps*(x0(2) - vd)^2;
phi1 = 2*(x0(2) - vd)/m;

b = x0(3) - c; %pacbf
b_dot = v0 - x0(2);
psi_1 = b_dot + p(1)*b^(p(2));

if(without_p2) 
    b_acbf = Fr/m + p(1)*b^(p(2))*p(2)*b_dot/b + p(3)*psi_1^(p(4)); %pacbf constraint
    A_acbf = [1/m, 0, -b^(p(2)), 0];  %(u, \delta for desired speed, \nu, \delta for stabilizing p1)
    C2 = b_dot + p(1)*b^(p(2));

    b_cbf1 = p(1) - 0; % cbf for p1, 0 <= p1 <= 3, upper bound of p1 is for the purpose of stablizing the system (not used in the case study)
    A_cbf1 = [0,0,-1,0];
    b_cbf2 = 3 - p(1);
    A_cbf2 = [0,0,1,0];
    % CLF for stabilizing p1
    V = (p(1) - 0.1)^2;
    b_clf = -eps*V;
    A_clf = [0,0,2*(p(1) - 0.1),-0.000000001]; %leveraging the cost and the stablizing mentioned in the paper, not necessary.
    
    data1(i,1:6) = [b, b_acbf, phi1, 1/m, -phi0, b_acbf];
    QQ = 200000; %pealty for relaxation, stabilizing p1
    %(u, \delta for desired speed, \nu, \delta for stabilizing p1)
    A = [phi1 -1 0 0;A_acbf;A_cbf1; A_cbf2; A_clf; 1 0 0 0;-1 0 0 0];
    b = [-phi0;b_acbf;b_cbf1;b_cbf2; b_clf; ca*m*g; cd*m*g];
    H = [2/(m^2) 0 0 0; 0 2*psc 0 0; 0 0 0 0; 0 0 0 QQ*psc];
    F = [-2*Fr/(m^2); 0; 2; 0];
else
    b_acbf = Fr/m + p(1)*b^(p(2))*p(2)*b_dot/b;
    A_acbf = [1/m, 0, -b^(p(2)), 0, -psi_1^(p(4))]; %(u, \delta for desired speed, \nu, \delta for stabilizing p1, p2)
    C2 = b_dot + p(1)*b^(p(2));

    b_cbf1 = p(1) - 0;
    A_cbf1 = [0,0,-1,0,0];
    b_cbf2 = 3 - p(1);
    A_cbf2 = [0,0,1,0,0];
    
    V = (p(1) - 0.103)^2;
    b_clf = -eps*V;
    A_clf = [0,0,2*(p(1) - 0.103),-1, 0]; %wihout leveraging, replace -1 by -0.00001 for leveraging
    data1(i,1:6) = [b, b_acbf, phi1, 1/m, -phi0, b_acbf];
    p3_ref = 1;
    QQ = 2000000000000; %without leveraging, set QQ to 200000 for leveraging;
    QU = 1;
    %(u, \delta for desired speed, \nu, \delta for stabilizing p1, p2)
    A = [phi1 -1 0 0 0;A_acbf;A_cbf1; A_cbf2; A_clf; 1 0 0 0 0;-1 0 0 0 0; 0 0 0 0 -1];%;LgBf 0;1 0
    b = [-phi0;b_acbf;b_cbf1;b_cbf2; b_clf; ca*m*g; cd*m*g; 0];%;LfBf;ca*m*g
    H = [2*QU/(m^2) 0 0 0 0; 0 2*psc 0 0 0; 0 0 0 0 0; 0 0 0 QQ*psc 0;0 0 0 0 QQ*psc];
    F = [-2*QU*Fr/(m^2); 0; 2; 0; -QQ*psc*p3_ref]; %100000
end
if(fea) %if the QP is feasible
  %  tic
    options = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','off');
    [u1,fval,exitflag,output,lambda] = ...
       quadprog(H,F,A,b,[],[],[],[],[],options);
 %  toc
    t=[0 0.1];
    if(numel(u1) == 0)
        fea = 0; te = 100; infea = i
    else
        u = u1;
    end
end
if(without_p2 == 0)
    p(3) = u(5); %update p2
end
if(i > te) % infeasibile, end the loop
    break;
end
data1(i,7:8) = [u(1)/m, u(2)];
[tt,xx]=ode45('acc',t,[x0, p(1)]); % integrate the dynamics
x0 = [xx(end, 1) xx(end, 2) xx(end, 3)]; %update state
p(1) = xx(end, 4); %update p1
p_history = [p_history;p];
result1(i,:) = [0.1*i xx(end,2) u(1)/m x0(3)-c C2];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot for adaptive method
 plot( data1(:,1), data1(:,7),'r','linewidth',1.6)
axis([0 90 -5 5]);
xlabel('$b(\textbf{\textit x}(t))$','fontsize',15)
ylabel('$u(t)/M$','fontsize',15)
grid on
hold on
plot([0,90],[ca*g, ca*g],'k:',[0,90],[-ca*g, -ca*g],'k:','linewidth',1.6)
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
plot(result1(:, 1),p_history(:, 1), 'b-',[0,30],[0,0], 'k--', 'linewidth',1.5)
xlabel('$t/s$','fontsize',15)
ylabel('$p_1(t)$','fontsize',15)
lg = legend('Linear','Quadratic');
set(lg,'box','off')
set(lg,'fontsize',12)
subplot(3, 1, 3)
plot(result1(:, 1),p_history(:, 3), 'b-',[0,30],[0,0], 'k--', 'linewidth',1.5)
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



 



 