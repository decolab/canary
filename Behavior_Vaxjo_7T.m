clear all
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Tenet/TENET/'];
addpath(genpath(path3));
path4=[ '../../TaskReservoir/DataHCP100ordered'];
addpath(genpath(path4));
path5=[ '../../TaskReservoir/DMF_Receptors'];
addpath(genpath(path5));
path6=[ '../../QuantumDiffMap/CHARM'];
addpath(genpath(path6));
path7=[ '../../EntropyProdandFlow/HCP'];
addpath(genpath(path7));
path8=[ '../../TurbuThermo/HCPrest'];
addpath(genpath(path8));

NSUB=182;

load results_model_vortex_rest_7T.mat;
load results_Empirical_Vaxjo_7T.mat;

%% Behav G-factor

load('hcpbehaviouraldata.mat')
load (['hcp7t_rfMRI_REST1_PA_schaefer1000.mat']);

behav=hcpbehaviouraldata;
for jj=1:1113
    behavsub(jj)=behav{jj,1};
end

for nsub=1:NSUB
    idsub=subject{nsub}.id;
    idx=find(idsub==behavsub);
    behav_PMAT24_A_CR(nsub)=behav{idx,122}; %% fluid int
    behav_ProcSpeed_Unadj(nsub)=behav{idx,129}; %% speed
    behav_CardSort_Unadj(nsub)=behav{idx,118};  %% executive
    behav_ListSort_Unadj(nsub)=behav{idx,158};
    behav_Flanker_Unadj(nsub)=behav{idx,120};
    behav_VSPLOT_TC(nsub)=behav{idx,145};
    behav_IWRD_TOT(nsub)=behav{idx,156};
    behav_PicSeq_Unadj(nsub)=behav{idx,116};  
    behav_RedEng_Unadj(nsub)=behav{idx,125};
    behav_PicVocab_Unadj(nsub)=behav{idx,127};  %% crystallized int

    behav_all(1,nsub)=behav{idx,122};
    behav_all(2,nsub)=behav{idx,129};
    behav_all(3,nsub)=behav{idx,118};
    behav_all(4,nsub)=behav{idx,158};
    behav_all(5,nsub)=behav{idx,120};
    behav_all(6,nsub)=behav{idx,145};
    behav_all(7,nsub)=behav{idx,156};
    behav_all(8,nsub)=behav{idx,116};
    behav_all(9,nsub)=behav{idx,125};
    behav_all(10,nsub)=behav{idx,127};
end

% behav_y=behav_PMAT24_A_CR';  %% 
% behav_y=behav_ProcSpeed_Unadj'; %%
 % behav_y=behav_CardSort_Unadj'; %% 
%  behav_y=behav_PicVocab_Unadj'; %% 
% 
% [~ , indx]=find(isnan(behav_y'));


%% G Factor
[~, indx]=find(isnan(behav_all));
behav_all(:,indx)=[];
[lambda,psi,T,stats,F] = factoran(behav_all',1);
behav_y=F;
%%

for nsub=1:NSUB
    % C=squeeze(Ceff(nsub,:,:));
     % C=squeeze(FCempbold(nsub,:,:));
     C=squeeze(Interferencevaxjo_sub(nsub,:,:));
    Cvec(nsub,:)=C(:);
end
Cvec(indx,:)=[];

Threshold2=0.14; %%0.14
alpha_EN=0.01;

%%
%% ===================================================
%  1. Feature selection on entire dataset
% ====================================================

% y = behav_y(:);
% Xall = Cvec;
% 
% for i = 1:size(Xall,2)
%     co2(i) = corr(Xall(:,i), y, 'rows','pairwise');
% end
% 
% selected = find(abs(co2) > Threshold2);
% 
% fprintf("Selected %d features.\n", length(selected));
% 
% X = Xall(:, selected);
% 
% % 2) Standardize predictors on full data (CPM convention)
% X = zscore(X);    

% 3) Outer CV to get honest per-subject predictions
for repe=1:100
    %%
    y = behav_y(:);
    Xall = Cvec;

    for i = 1:size(Xall,2)
        co2(i) = corr(Xall(:,i), y, 'rows','pairwise');
    end

    selected = find(abs(co2) > Threshold2);

    fprintf("Selected %d features.\n", length(selected));

    X = Xall(:, selected);

    % 2) Standardize predictors on full data (CPM convention)
    X = zscore(X);    
    %%

    K = 10;
    cv = cvpartition(length(y),'KFold',K);
    yhat = zeros(size(y));

    for fold = 1:K
        % fprintf('Fold %d/%d\n', fold, K);
        train_idx = training(cv, fold);
        test_idx  = test(cv, fold);

        Xtrain = X(train_idx, :);
        ytrain = y(train_idx);

        Xtest  = X(test_idx, :);

        % Inside-fold: pick lambda with internal CV on the training set
        % NOTE: 'Standardize',false because we already z-scored globally by convention
        [B, FitInfo] = lasso(Xtrain, ytrain, 'Alpha', alpha_EN, 'CV', 10, 'Standardize', false);

        % Get lambda index that minimized MSE *on the training fold's internal CV*
        idxMin = FitInfo.IndexMinMSE;
        lambda_min = FitInfo.Lambda(idxMin);

        % Refit on the entire training set at the selected lambda
        % (lasso allows specifying Lambda to fit)
        [B_retrain, FitInfo2] = lasso( ...
            Xtrain, ytrain, ...
            'Alpha', alpha_EN, ...
            'Lambda', lambda_min, ...
            'Standardize', false);
        intercept = FitInfo2.Intercept;

        % Predict test set
        yhat(test_idx) = Xtest * B_retrain + intercept;
    end

    %% ===================================================
    %  4. Metrics
    % ===================================================

    SST = sum((y - mean(y)).^2);
    SSE = sum((y - yhat).^2);

    R2 = 1 - SSE/SST;
    Corr = corr(y, yhat);

    fprintf('\n========== RESULTS ==========\n');
    fprintf('Cross-validated R2   = %.4f\n', R2);
    fprintf('Cross-validated Corr = %.4f\n', Corr);
    R2all(repe)=R2;
    Corrall(repe)=Corr;
    yhatall(repe,:)=yhat;
end

yhatall=squeeze(mean(yhatall));

figure;
violinplot([R2all' Corrall']);

figure;
scatter(y,yhatall,'filled'); 
hold on;
p = polyfit(yhatall, y, 1);
xvals = linspace(min(yhatall), max(yhatall));
plot(xvals, polyval(p, xvals), 'r', 'LineWidth', 2)
hold off;
[cc pp]=corrcoef(y,yhatall)

save results_Behavior_Vaxjo_7T_Gfactor.mat y yhatall R2all Corrall;

%% EXTRA RENDERING
sel=zeros(1,10000);
sel(selected)=1;
msel=reshape(sel,100,100);
[aux idx]=sort(F,'descend');
topbad=idx(end-15:end);
xb=mean(Xall(topbad,:));
topgood=idx(1:15);
xg=mean(Xall(topgood,:));
mg=reshape(xg,100,100).*msel;
mb=reshape(xb,100,100).*msel;
gbcgood=mean(mg);
gbcbad=mean(mb);

save results_rendering_Vaxjo_7T_Gfactor.mat gbcgood gbcbad;