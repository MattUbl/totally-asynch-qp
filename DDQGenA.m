function [E5,E20,E50,Q,A5,A20,A50,q] = DDQGenA(qi)
good = 0;
z = 0;
%qr = .1
%qi = .99;
%while good == 0
% size of matrix
n = 20;
%Generate initial random Q
Q = rand(n);
q5 = qi*(1-.05);
q20 = qi*(1-.15);
q50 = qi*(1-.45);
e = .05;
%scale each row such that the matrix is (not strictly) diagonally dominant
Q = Q+Q';
for i = 1:n
    T = Q(i,:);
    T(i) = 0;
    Q(i,i) = norm(T,1)/qi;
end
%take off-diagonal pairs (e.g. Q(i,j) and Q(j,i)) and set each to the small
%value of the two, preserves diagonal dominance and ensures symmetry
% for i = 1:n
%     for j = 1:n
%         if Q(i,j) > Q(j,i)
%             Q(i,j) = Q(j,i);
%         else
%             Q(i,j) = Q(i,j);
%         end
%     end
% end
%above process usually results in a matrix with a really small q, in order
%to initially increase q, we calculate the beta diagonal dominance of each
%row, choose a fraction of the smallest beta_i, and subtract that from each
%diagonal element. Preserves diagonal dominance, and scales q up to
%something closer to 1
% Beta = NaN(n,1);
% for i = 1:n
%     Beta(i) = 2*Q(i,i)-norm(Q(i,:),1);
% end
% Beta;
% 
% peturb = .1*min(Beta);
% Q = Q - peturb*eye(n);
Q = Q*10;
%recalculate Beta for the new matrix, to be used later
Beta = NaN(n,1);
for i = 1:n
    Beta(i) = 2*Q(i,i)-norm(Q(i,:),1);
end
%calculate the q term for each row, the largest is q for the system
q = NaN(n,1);
for i = 1:n
    q(i) = (norm(Q(i,:),1)-Q(i,i))/Q(i,i);
end
%calculate the required regularization to hit our desired q* for each row
alpha5 = NaN(n,1);
for i = 1:n
    alpha5(i) = max([Q(i,i)*(q(i)/q5-1) 0]);
end
A5 = diag(alpha5);

alpha20 = NaN(n,1);
for i = 1:n
    alpha20(i) = max([Q(i,i)*(q(i)/q20-1) 0]);
end
A20 = diag(alpha20);

alpha50 = NaN(n,1);
for i = 1:n
    alpha50(i) = max([Q(i,i)*(q(i)/q50-1) 0]);
end
A50 = diag(alpha50);
%Based on these regularizations, calculate our error bound

E5 = max((alpha5./Beta).^2./(alpha5./Beta+1).^2);
E20 = max((alpha20./Beta).^2./(alpha20./Beta+1).^2);
E50 = max((alpha50./Beta).^2./(alpha50./Beta+1).^2);

qi

%q

%if cond(Q) > 10
         %   good = 1
%end
%z = z+1;
%cond(Q)
%end