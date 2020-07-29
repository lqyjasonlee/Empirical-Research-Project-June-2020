%(2) naive estimator
%(9) post-double-selection

% design
% 1 exact sparsity 
% 22 approximate sparsity

clear;

NUM_SIM = 10000;

p = 200; 
n = 100; 

alpha0 = 1/2;

design = 1;

rho = 0.5;


NumValC1 = 1;
NumValC2 = 1;
ALL_COVERAGE = zeros(1,1,16);
ALL_BIAS = zeros(1,1,16);
ALL_SD = zeros(1,1,16);

for i1 = 1 : 1 : 1
    for i2 = 1 : 1 : 1


R21 = 0.5;
R22 = 0.5;

fprintf('R-square First Stage %f\n', R21);                
fprintf('R-square Second Stage %f\n', R22);                
        
[ ALL_StdErr, ALL_ALPHA ] = MC_TE_FixedDesign_Heteroskedastic_Lasso_RedForm ( NUM_SIM, rho, alpha0, R21, R22, design, p, n );

%%% (2) naive estimator 
Bias_postLASSO = mean(ALL_ALPHA(:,2) - alpha0);
SD_postLASSO = sqrt(var(ALL_ALPHA(:,2)));
Zvalue = (ALL_ALPHA(:,2) - alpha0)./ALL_StdErr(:,2);

%%% (9) post-double-selection
Bias_9 = mean(ALL_ALPHA(:,9) - alpha0);
SD_9 = sqrt(var(ALL_ALPHA(:,9)));
Zvalue_pdse = (ALL_ALPHA(:,9) - mean(ALL_ALPHA(:,9)))./ALL_StdErr(:,9);
    end
end


