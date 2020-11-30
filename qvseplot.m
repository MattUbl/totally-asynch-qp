qhist = [0:1/1000:1];
qistack = [.99 .95 .85 .7 .5 .3 .01];
ehist = NaN(length(qhist),length(qistack));
for j = 1:length(qistack)
for i = 1:length(qhist)
    qr = qhist(i);
    E = DDQGen(qr,qistack(j));
    ehist(i,j) = E;
end
end

t = [0:.1:100];

figure(1)
plot(t,ehist(:,1),t,ehist(:,2),'--',t,ehist(:,3),':',t,ehist(:,4),'-.',t,ehist(:,5),':',t,ehist(:,6),'--',t,ehist(:,7),'LineWidth',2)
title('Regularization vs Relative Cost Error')
xlabel('% Reduction in q','FontWeight','Bold')
ylabel('Error Bound \epsilon','FontWeight','Bold')
legend({'q_{initial} = 0.99','q_{initial} = 0.95','q_{initial} = 0.85','q_{initial} = 0.70','q_{initial} = 0.50','q_{initial} = 0.30','q_{initial} = 0.01'},'Location','SouthEast','FontWeight','Bold')
% hold off