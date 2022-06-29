function PRTable = helperPrecisionRecall(confmat)
% This function is only in support of XpwWaveletMLExample. It may change or
% be removed in a future release.
precisionA = confmat(1,1)/sum(confmat(:,1))*100;
precisionB = confmat(2,2)/sum(confmat(:,2))*100;
precisionC = confmat(3,3)/sum(confmat(:,3))*100;
precisionD = confmat(4,4)/sum(confmat(:,4))*100;
precisionE = confmat(5,5)/sum(confmat(:,5))*100;
recallA = confmat(1,1)/sum(confmat(1,:))*100;
recallB = confmat(2,2)/sum(confmat(2,:))*100;
recallC = confmat(3,3)/sum(confmat(3,:))*100;
recallD = confmat(4,4)/sum(confmat(4,:))*100;
recallE = confmat(5,5)/sum(confmat(5,:))*100;
F1A = 2*precisionA*recallA/(precisionA+recallA);
F1B = 2*precisionB*recallB/(precisionB+recallB);
F1C = 2*precisionC*recallC/(precisionC+recallC);
F1D = 2*precisionD*recallD/(precisionD+recallD);
F1E = 2*precisionE*recallE/(precisionE+recallE);
% Construct a MATLAB Table to display the results.
PRTable = array2table([precisionA recallA F1A;...
    precisionB recallB F1B; precisionC recallC...
    F1C; precisionD recallD F1D; precisionE recallE F1E;...
    ],'VariableNames',{'Precision','Recall','F1_Score'},'RowNames',...
    {'0%','0.01%','0.1%','1%','2%'});

end