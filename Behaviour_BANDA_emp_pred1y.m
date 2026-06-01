clear all;
addpath(genpath('../../Turbulence/Basics'));
addpath(genpath('../../Turbulence/SC_longrange'));
addpath(genpath('../../Nonequilibrium'));
addpath(genpath('../../BANDAdata'));
addpath(genpath('../../BANDAdata/sFC_Schaefer2018_1000parcels_7Networks_order_REST_AP'));

CASOemp    = 1;      % 1 = Vaxjo interference, 2 = FC
alpha_EN   = 0.01;   % Elastic Net Alpha (Low value = Ridge-dominant for stability)
nRepeats   = 100;    % Number of resampling repetitions for robustness analysis

if CASOemp==1
    load results_BANDA_Hopf_Ceff_Vaxjo.mat;
    clear yhatfinal R2final Corrfinal;
end
if CASOemp==2
    load results_BANDA_Vaxjo_emp.mat;
end
T = readtable('rcads01.xlsx');

%% Loading IDs
load ("Original_groups_ID.mat")
load ("subjects_BANDA_char.mat")

%% Extract participants ids from fmri file order
subjNames_full = strtrim(string(cellstr(subjects)));
subjNames_base = extractBefore(subjNames_full, "_");

%% Group ID
groupIDs = struct();
groupIDs.CNT = strtrim(string(CNT_ID(:)));
groupIDs.ANX = strtrim(string(ANX_ID(:)));
groupIDs.DEP = strtrim(string(DEP_ID(:)));

%% Get indices in subject{} for each group
groupIdx = struct();
fn = fieldnames(groupIDs);
for i = 1:numel(fn)
    g = fn{i};
    ids = groupIDs.(g);
    [tf, loc] = ismember(ids, subjNames_base);
    groupIdx.(g) = loc(tf);
end

%% Loop and read timeseries per group
mindur=340; maxdur=530;
badguys={'BANDA005','BANDA028','BANDA195','BANDA205'}; % High motion exclusion

nn=1; nn1=0; nn2=0; nn3=0;
for i = 1:numel(fn)
    g = fn{i};
    ids = groupIDs.(g);

    if CASOemp == 1
        InterferenceMatrix = Interferencevaxjo_sub.(g);
    else
        InterferenceMatrix = FCemp_sub.(g);
    end

    for n = 1:size(ids,1)
        loc = find(ids(n)==string(T{:,5}));

        % Data extraction logic based on visit availability
        if size(loc,1)==4
            dur=days(datetime(T{loc(3),6})-datetime(T{loc(1),6}));
            if dur>mindur && dur<maxdur && ~ismember(ids(n),badguys)
                rcads_base(nn)=T{loc(1),59};
                rcads_1year(nn)=T{loc(3),59};
                C=squeeze(InterferenceMatrix(n,:,:));
                Cvec(nn,:)=C(:);
                finalID(nn)=ids(n);
                nn=nn+1; nn1=nn1+1;
            end
        elseif size(loc,1)==3
            dur=days(datetime(T{loc(3),6})-datetime(T{loc(1),6}));
            if dur>mindur && dur<maxdur && ~ismember(ids(n),badguys)
                rcads_base(nn)=T{loc(1),59};
                rcads_1year(nn)=T{loc(3),59};
                C=squeeze(InterferenceMatrix(n,:,:));
                Cvec(nn,:)=C(:);
                finalID(nn)=ids(n);
                nn=nn+1; nn2=nn2+1;
            end
        elseif size(loc,1)==2
            dur=days(datetime(T{loc(2),6})-datetime(T{loc(1),6}));
            if dur>mindur && dur<maxdur && ~ismember(ids(n),badguys)
                rcads_base(nn)=T{loc(1),59};
                rcads_1year(nn)=T{loc(2),59};
                C=squeeze(InterferenceMatrix(n,:,:));
                Cvec(nn,:)=C(:);
                finalID(nn)=ids(n);
                nn=nn+1; nn3=nn3+1;
            end
        end
    end
end

%% =========================================================
%  INFERENCE FRAMEWORK (NESTED THRESHOLD SELECTION)
% =========================================================
behav_y = str2double(rcads_base);
y1 = behav_y(:);
behav_y = str2double(rcads_1year);
y2 = behav_y(:);

Xall = Cvec;

% 1. Global selection based on baseline y1 (No y2 leakage here!)
co2 = corr(Xall, y1, 'rows', 'pairwise');
thrange = 0.08:0.01:max(abs(co2))-0.03;

%% 2. ROBUSTNESS ANALYSIS (NESTED CV FOR THRESHOLD)
R2all = zeros(nRepeats, 1);
Corrall = zeros(nRepeats, 1);
y_recon_accum = zeros(nRepeats, length(y2));

fprintf('\nRunning %d Nested Resampling Repetitions...\n', nRepeats);

for repe = 1:nRepeats
    cv_outer = cvpartition(length(y2), 'KFold', 10);
    y_recon = zeros(size(y2));
    
    for fold = 1:10
        tr_out = training(cv_outer, fold);
        te_out  = test(cv_outer, fold);
        
        % --- INNER LOOP: Find best threshold using ONLY outer training data ---
        cv_inner = cvpartition(sum(tr_out), 'KFold', 5);
        X_tr_outer = Xall(tr_out, :);
        y2_tr_outer = y2(tr_out);
        
        best_th = thrange(1);
        best_inner_corr = -Inf;
        
        for th = thrange
            sel_inner = find(abs(co2) > th);
            if length(sel_inner) < 2
                continue; % skip if too few features
            end
            
            y_recon_in = zeros(size(y2_tr_outer));
            for infold = 1:5
                tr_in = training(cv_inner, infold);
                te_in = test(cv_inner, infold);
                
                X_in_tr = X_tr_outer(tr_in, sel_inner);
                X_in_te = X_tr_outer(te_in, sel_inner);
                
                % Inner Z-scoring
                mu_in = mean(X_in_tr, 1); sig_in = std(X_in_tr, 0, 1); sig_in(sig_in==0)=1;
                X_in_tr_z = (X_in_tr - mu_in) ./ sig_in;
                X_in_te_z = (X_in_te - mu_in) ./ sig_in;
                
                % Fast inner lasso
                [B_in, FitInfo_in] = lasso(X_in_tr_z, y2_tr_outer(tr_in), 'Alpha', alpha_EN, 'CV', 3);
                y_recon_in(te_in) = X_in_te_z * B_in(:, FitInfo_in.IndexMinMSE) + FitInfo_in.Intercept(FitInfo_in.IndexMinMSE);
            end
            
            % Check if this threshold is the best for this specific training fold
            in_corr = corr(y2_tr_outer, y_recon_in);
            if in_corr > best_inner_corr
                best_inner_corr = in_corr;
                best_th = th;
            end
        end
        
        % --- OUTER LOOP: Apply the best nested threshold to unseen test data ---
        sel_final = find(abs(co2) > best_th);
        
        X_tr_final = Xall(tr_out, sel_final);
        X_te_final = Xall(te_out, sel_final);
        
        % Outer Z-scoring
        mu_out = mean(X_tr_final, 1); sig_out = std(X_tr_final, 0, 1); sig_out(sig_out==0)=1;
        X_tr_z = (X_tr_final - mu_out) ./ sig_out;
        X_te_z = (X_te_final - mu_out) ./ sig_out;

        % Outer Lasso Optimization
        [B, FitInfo] = lasso(X_tr_z, y2(tr_out), 'Alpha', alpha_EN, 'CV', 5, 'Standardize', false);
        idx_opt = FitInfo.IndexMinMSE;
        y_recon(te_out) = X_te_z * B(:, idx_opt) + FitInfo.Intercept(idx_opt);
    end
    
    % Calculate Metrics for this repetition
    SST = sum((y2 - mean(y2)).^2);
    SSE = sum((y2 - y_recon).^2);
    R2all(repe) = 1 - SSE/SST;
    Corrall(repe) = corr(y2, y_recon);
    y_recon_accum(repe,:) = y_recon;
    
    if mod(repe, 10) == 0, fprintf('Completed %d/%d repeats...\n', repe, nRepeats); end
end

fprintf('\n========== INFERENCE RESULTS ==========\n');
fprintf('Mean Reconstruction R2   = %.4f\n', mean(R2all));
fprintf('Mean Reconstruction Corr = %.4f\n', mean(Corrall));
R2final   = R2all;
Corrfinal = Corrall;
yhatfinal = mean(y_recon_accum, 1)'; % Transposed to column vector

%% SAVE
if CASOemp == 1
    save results_Behaviour_BANDA_emp_pred1y_Vaxjo.mat R2final Corrfinal yhatfinal y1 y2
else
    save results_Behaviour_BANDA_emp_pred1y_FC.mat R2final Corrfinal yhatfinal y1 y2
end

%% PLOT

figure
violinplot([R2all' Corrall']);

figure
scatter(yhatfinal, y2, 'filled');
hold on;
p = polyfit(yhatfinal, y2, 1);
xvals = linspace(min(yhatfinal), max(yhatfinal));
plot(xvals, polyval(p, xvals), 'r', 'LineWidth', 2)
hold off;
[cc pp] = corrcoef(y2, yhatfinal)

%% EXTRA RENDERING

sel=zeros(1,10000);
sel(sel_final)=1;
msel=reshape(sel,100,100);
[aux idx]=sort(str2double(rcads_1year),'descend');
topbad=idx(1:15);
topgood=idx(end-15:end);
xb=mean(Xall(topbad,:));
xg=mean(Xall(topgood,:));
mg=reshape(xg',100,100).*msel;
mb=reshape(xb',100,100).*msel;
gbcgood_1y=mean(mg);
gbcbad_1y=mean(mb);

[aux idx]=sort(str2double(rcads_base),'descend');
topbad=idx(1:15);
topgood=idx(end-15:end);
xb=mean(Xall(topbad,:));
xg=mean(Xall(topgood,:));
mg=reshape(xg,100,100).*msel;
mb=reshape(xb,100,100).*msel;
gbcgood_base=mean(mg);
gbcbad_base=mean(mb);

if CASOemp == 1
    save results_rendering_BANDA_Vaxjo.mat gbcgood_1y gbcbad_1y gbcgood_base gbcbad_base;
else
    save results_rendering_BANDA_FC.mat gbcgood_1y gbcbad_1y gbcgood_base gbcbad_base;
end






