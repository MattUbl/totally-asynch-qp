clear all
close all 
[E5,E20,E50,Q,A5,A20,A50,q] = DDQGenA(.85);
iter = 999;
n = 20;
%A = 1*diag(alpha);
QA5 = Q+A5;
QA20 = Q+A20;
QA50 = Q+A50;
r = Q*1000*ones(n,1);
Sol = -Q\r;
nSol = norm(Sol);
SolA5 = -QA5\r;
SolA20 = -QA20\r;
SolA50 = -QA50\r;
%nSolA = norm(SolA);
f = 0.5*r'*Sol;
fA5 = -0.5*r'*(inv(QA5)*Q-2*eye(n))*SolA5;
fA20 = -0.5*r'*(inv(QA20)*Q-2*eye(n))*SolA20;
fA50 = -0.5*r'*(inv(QA50)*Q-2*eye(n))*SolA50;

% D = norm(Sol,Inf); %The worst performing block at time zero
% DA = norm(SolA,Inf); %The worst performing block at time zero
U = zeros(n,1); %Keeps track which agents have computed an update this cycle
S = eye(n); %If agents have computed an update, keeps track of which agents have recieved this update (eye(25) used becaue agents already "sent" it to themselves)
C = 0; %Number of cycles that have been completed

X = zeros(n);
XA5 = zeros(n);
XA20 = zeros(n);
XA50 = zeros(n);

G = 1./diag(diag(Q));
GA5 = 1./diag(diag(QA5));
GA20 = 1./diag(diag(QA20));
GA50 = 1./diag(diag(QA50));

b = NaN(iter,1); %Theoretical error bound for dHom
e = NaN(iter,1); %Error of Hom at each timestep
R = NaN(n,1); %Actual state of the Hom network
d = NaN(iter,1); %Worst-performing state variable on any agent, used for convergence rate verification

bA5 = NaN(iter,1);
eA5 = NaN(iter,1);
RA5 = NaN(n,1);
dA5 = NaN(iter,1);

bA20 = NaN(iter,1);
eA20 = NaN(iter,1);
RA20 = NaN(n,1);
dA20 = NaN(iter,1);

bA50 = NaN(iter,1);
eA50 = NaN(iter,1);
RA50 = NaN(n,1);
dA50 = NaN(iter,1);

last = 0;

for k = 1:iter
    for i = 1:n
        up = rand(1);
        if up < .1 %Probability agent i computes an update
            X(i,i) = X(i,i)-G(i,i)*(Q(i,:)*X(:,i)+r(i));
            XA5(i,i) = XA5(i,i)-GA5(i,i)*(QA5(i,:)*XA5(:,i)+r(i));
            XA20(i,i) = XA20(i,i)-GA20(i,i)*(QA20(i,:)*XA20(:,i)+r(i));
            XA50(i,i) = XA50(i,i)-GA50(i,i)*(QA50(i,:)*XA50(:,i)+r(i));
            U(i) = 1; %If agent i updates, record it in U
        end
    end
    for i = 1:n
        for j = 1:n
            com = rand(1);
            probcom = exp(log(.01)*(1-last/200));
            if com < .01 %Probability agent j sends its state to agent i
                X(j,i) = X(j,j); %Agent j shares with agent i
                XA5(j,i) = XA5(j,j);
                XA20(j,i) = XA20(j,j);
                XA50(j,i) = XA50(j,j);
                if U(j) == 1 %ONLY if agent j has updated prior to this transmission during this cycle
                    S(i,j) = 1; %IF agent j shares with agent i, record it in S
                end
           % elseif last == 200
               % X(j,i) = X(j,j); %Agent j shares with agent i
               % XA5(j,i) = XA5(j,j);
              %  XA20(j,i) = XA20(j,j);
               % XA50(j,i) = XA50(j,j);
               % if U(j) == 1 %ONLY if agent j has updated prior to this transmission during this cycle
               %     S(i,j) = 1; %IF agent j shares with agent i, record it in S
               % end
            end
        end
    end
        R = diag(X); %The actual state of the system, i.e. each agent's copy of its own state
        RA5 = diag(XA5);
        RA20 = diag(XA20);
        RA50 = diag(XA50);
    e(k) = norm(R-Sol)/nSol; %Normalized system error, distance to Sol divided by norm of Sol
    eA5(k) = norm(RA5-Sol)/nSol;
    eA20(k) = norm(RA20-Sol)/nSol;
    eA50(k) = norm(RA50-Sol)/nSol;
%     d(k) = max(vecnorm(X-Sol,Inf)); %The norm fo the worst performing block in the system
%     dA(k) = max(vecnorm(XA-Sol,Inf));
%     b(k) = D*max(q)^C; %Convergence bound
%     bA(k) = DA*.79^C;
    if norm(U,1) == length(U)
        if norm(S,'fro') == n %If every agent has updated and shared this update with every other agent (i.e. U and S are all ones)
            C = C+1; %Move to next cycle
            U = zeros(n,1); %Reset U and S
            S = eye(n);
            last = 0;
            C
            k
        else
            last = last+1;
        end
    end
end

qA5 = NaN(n,1);
for i = 1:n
    qA5(i) = (norm(QA5(i,:),1)-QA5(i,i))/QA5(i,i);
end

qA20 = NaN(n,1);
for i = 1:n
    qA20(i) = (norm(QA20(i,:),1)-QA20(i,i))/QA20(i,i);
end

qA50 = NaN(n,1);
for i = 1:n
    qA50(i) = (norm(QA50(i,:),1)-QA50(i,i))/QA50(i,i);
end

error = [(f-fA5)/f (f-fA20)/f (f-fA50)/f]';

% figure(1)
% plot(bA,'LineWidth',2)
% hold on
% plot(dA,'LineWidth',2)
% title('Independently Regularized Convergence')
% xlabel('Iteration Number')
% ylabel('Worst-Performing Block (Variable)')
% legend('Error Bound','Error')
% hold off
% 
% figure(2)
% plot(d,'LineWidth',2)
% hold on
% plot(b,'LineWidth',2)
% title('Homogenous Convergence')
% xlabel('Iteration Number')
% ylabel('Worst-Performing Block (Variable)')
% legend('Error','Error Bound')
% hold off

% figure(3)
% plot(d,'LineWidth',2)
% hold on
% plot(dA,'LineWidth',2)
% title('Convergence Comparison (Worst Block)')
% xlabel('Iteration Number')
% ylabel('Worst-Performing Block (Variable)')
% legend('Homogenous Stepsize','Independent Stepsizes')
% hold off

figure(4)
plot([1; e],'LineWidth',2)
hold on
plot([1; eA5],'--','LineWidth',2)
plot([1; eA20],':','LineWidth',2)
plot([1; eA50],'-.','LineWidth',2)
title('Convergence Comparison (System Norm)')
xlabel('Iteration Number','FontWeight','Bold')
ylabel('Norm of System Error (From Unregularized Solution)','FontWeight','Bold')
legend({'Unregularized','5% Reduction','15% Reduction','45% Reduction'},'FontWeight','Bold')
axes('position', [.25 .5 .35 .35])
box on;
plot(100:1:200,[e(100:1:200)],'LineWidth',2) 
hold on;
plot(100:1:200,[eA5(100:1:200)], '--','LineWidth',2) 
plot(100:1:200,[eA20(100:1:200)], ':','LineWidth',2) 
plot(100:1:200,[eA50(100:1:200)], '-.','LineWidth',2) 
ylim([0.05 0.25])
hold off