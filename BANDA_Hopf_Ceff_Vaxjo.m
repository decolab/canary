clear all;
path2=[ '../../Turbulence/Basics'];
addpath(genpath(path2));
path3=[ '../../Turbulence/SC_longrange'];
addpath(genpath(path3));
path4=[ '../../Nonequilibrium'];
addpath(genpath(path4));
path5=[ '../../BANDAdata'];
addpath(genpath(path5));
path6=[ '../../BANDAdata/sFC_Schaefer2018_1000parcels_7Networks_order_REST_AP'];
addpath(genpath(path6));

load ('schaefer1000to100.mat');
partitiondef=schaefer1000to100;
NP=100;
N=1000;
Tmax=2000;
Tmax1=2500;

Isubdiag = find(tril(ones(N),-1));

%% Loading the SC, anatomy of HCP in schaefer 1000

load sc_schaefer_1000.mat;
C=sc_schaefer;
C=C-diag(diag(C));
C = C/max(max(C));

epsFC=0.0001;
maxC=0.2;

C=C*maxC;

for i=1:NP
    partition{i}=find(partitiondef==i);
end

%%
% Parameters of the data
TR=0.8;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

load hpcdata1003_f_diff_fce.mat;
f_diff=f_diff';

sig = 0.01;
dt=0.1*TR/2;
dsig = sqrt(dt)*sig;

AA=-0.02;

a=AA*ones(N,1);
a=repmat(a,1,2);
wo = f_diff*(2*pi);
omega = repmat(wo,1,2);
omega(:,1) = -omega(:,1);


%% Loading IDs and timeseries
load ("Original_groups_ID.mat")
load ("subjects_BANDA_char.mat")
load ("sFC_Schaefer2018_1000Parcels_7Networks_order_REST1_AP.mat")
%% Extracting parcitipants ids from fmri file order
% subjects is 202x11 char -> convert to clean string array (202x1)
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
    [tf, loc] = ismember(ids, subjNames_base);   % loc are indices into subject{}
    groupIdx.(g) = loc(tf);                % valid indices
end
%% Example: loop and read timeseries per group
for i = 1:numel(fn)
    g = fn{i};
    clear vaxjo2 Interferencevaxjo_submat;
    vaxjo2NULLall=[];
    idxs = groupIdx.(g);
    for nsub = 1:numel(idxs)
        nsub
        s = subject{idxs(nsub)};          % 1x1 struct
        ts = s.schaeferts;            % 100 x T (from your screenshot)
        clear signal_filt Phases;
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            if sum(isnan(ts(seed,:)))==0
                signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
            end
        end

        ts2=signal_filt(:,20:end-20);
        FCemp=corrcoef(ts2');

        %% GEC

        Cnew=C;
        olderror=0.001;
        errornow=0;
        nume=0;
        for iter=1:5000 %% Optimization iteration
            wC=Cnew;
            sumC = repmat(sum(wC,2),1,2);
            xs=zeros(Tmax,N);
            z = 0.1*ones(N,2);
            nn=0;
            % discard first 2000 time steps
            for t=0:dt:1000
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
            end
            % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
            for t=0:dt:((Tmax-1)*TR)
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
                if abs(mod(t,TR))<0.01
                    nn=nn+1;
                    xs(nn,:)=z(:,1)';
                end
            end
            ts=xs';
            clear signal_filt Phases;
            for seed=1:N
                ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
                signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
            end
            ts2=signal_filt(:,100:end-100);
            FCsim=corrcoef(ts2');

            aux=corrcoef(FCsim(Isubdiag),FCemp(Isubdiag),'Rows','complete');
            corrFC(iter)=aux(2);

            if mod(iter,10)<0.1
                errornow=errornow/nume;
                if  abs(errornow-olderror)/olderror<0.01
                    break;
                end
                % if  olderror>errornow
                %     break;
                % end
                olderror=errornow;
                errornow=0;
                nume=0;
            else
                aux=corrcoef(FCsim(Isubdiag),FCemp(Isubdiag),'Rows','complete');
                errornow=errornow+aux(2);
                nume=nume+1;
            end

            for i=1:N  %% learning
                for j=1:N
                    if (C(i,j)>0 || abs(j-i)==N/2)
                        if ~isnan(FCemp(i,j))
                            Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j));
                        end
                    end
                end
            end
            Cnew = Cnew/max(max(Cnew))*maxC;
        end
        FCfitting(nsub)=mean(corrFC(iter-10:iter));

        %% Spectral gap and entangelment
        [V D]=eig(Cnew);
        D=diag(real(D));
        Dall=D;
        [dmax]=sort(D,'descend');
        spacings = abs(dmax);
        gapThreshold = mean(spacings) + std(spacings);
        gapIdx = find(spacings > gapThreshold);
        sumgap=0;
        if ~isempty(gapIdx)
            for gap = 1:numel(gapIdx)
                sumgap=sumgap+spacings(gapIdx(gap));
            end
        end
        spectral_gap.(g)(nsub)=sumgap;

        %%
        wC=Cnew;
        sumC = repmat(sum(wC,2),1,2);
        xs=zeros(Tmax1,N);
        z = 0.1*ones(N,2);
        nn=0;
        % discard first 2000 time steps
        for t=0:dt:1000
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        end
        % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
        for t=0:dt:((Tmax1-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end
        ts=xs';
        clear signal_filt Phases;
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            if sum(isnan(ts(seed,:)))==0
                signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
                Xanalytic = hilbert(demean(signal_filt(seed,:)));
                Phases(seed,:) = angle(Xanalytic);
            end
        end

        Phases=Phases(:,100:end-100);

        % Define Partition

        for i=1:NP
            Vorticity(i,:)=abs(sum(complex(cos(Phases(partition{i},:)),sin(Phases(partition{i},:))))/length(partition{i}));
        end
        Turbulence.(g)(nsub)=std(Vorticity(:));

        Tmax2=size(Vorticity,2);
        NullVorticity=Vorticity(randperm(NP),randperm(Tmax2));

        %%%%
        %% Vaxjol
        for i=1:NP
            for j=1:Tmax2
                if Vorticity(i,j)<0.2
                    binvortex(i,j)=0;
                end
                if Vorticity(i,j)>=0.2 && Vorticity(i,j)<0.4
                    binvortex(i,j)=1;
                end
                if Vorticity(i,j)>=0.4 && Vorticity(i,j)<0.6
                    binvortex(i,j)=2;
                end
                if Vorticity(i,j)>=0.6 && Vorticity(i,j)<0.8
                    binvortex(i,j)=3;
                end
                if Vorticity(i,j)>=0.8
                    binvortex(i,j)=4;
                end

                %% surro
                if NullVorticity(i,j)<0.2
                    binvortexNULL(i,j)=0;
                end
                if NullVorticity(i,j)>=0.2 && NullVorticity(i,j)<0.4
                    binvortexNULL(i,j)=1;
                end
                if NullVorticity(i,j)>=0.4 && NullVorticity(i,j)<0.6
                    binvortexNULL(i,j)=2;
                end
                if NullVorticity(i,j)>=0.6 && NullVorticity(i,j)<0.8
                    binvortexNULL(i,j)=3;
                end
                if NullVorticity(i,j)>=0.8
                    binvortexNULL(i,j)=4;
                end
            end
        end

        npa=1;
        for i=1:NP
            for j=1:NP
                if i~=j
                    bin1=binvortex(i,:);
                    bin2=binvortex(j,:);

                    P1_0=length(find(bin1==0))/length(bin1);
                    P1_1=length(find(bin1==1))/length(bin1);
                    P1_2=length(find(bin1==2))/length(bin1);
                    P1_3=length(find(bin1==3))/length(bin1);
                    P1_4=length(find(bin1==4))/length(bin1);
                    P2_0=length(find(bin2==0))/length(bin1);
                    P2_1=length(find(bin2==1))/length(bin1);
                    P2_2=length(find(bin2==2))/length(bin1);
                    P2_3=length(find(bin2==3))/length(bin1);
                    P2_4=length(find(bin2==4))/length(bin1);
                    P12_00=length(intersect(find(bin1==0),find(bin2==0)))/length(find(bin2==0));
                    P12_01=length(intersect(find(bin1==0),find(bin2==1)))/length(find(bin2==1));
                    P12_02=length(intersect(find(bin1==0),find(bin2==2)))/length(find(bin2==2));
                    P12_03=length(intersect(find(bin1==0),find(bin2==3)))/length(find(bin2==3));
                    P12_04=length(intersect(find(bin1==0),find(bin2==4)))/length(find(bin2==4));
                    P12_10=length(intersect(find(bin1==1),find(bin2==0)))/length(find(bin2==0));
                    P12_11=length(intersect(find(bin1==1),find(bin2==1)))/length(find(bin2==1));
                    P12_12=length(intersect(find(bin1==1),find(bin2==2)))/length(find(bin2==2));
                    P12_13=length(intersect(find(bin1==1),find(bin2==3)))/length(find(bin2==3));
                    P12_14=length(intersect(find(bin1==1),find(bin2==4)))/length(find(bin2==4));
                    P12_20=length(intersect(find(bin1==2),find(bin2==0)))/length(find(bin2==0));
                    P12_21=length(intersect(find(bin1==2),find(bin2==1)))/length(find(bin2==1));
                    P12_22=length(intersect(find(bin1==2),find(bin2==2)))/length(find(bin2==2));
                    P12_23=length(intersect(find(bin1==2),find(bin2==3)))/length(find(bin2==3));
                    P12_24=length(intersect(find(bin1==2),find(bin2==4)))/length(find(bin2==4));
                    P12_30=length(intersect(find(bin1==3),find(bin2==0)))/length(find(bin2==0));
                    P12_31=length(intersect(find(bin1==3),find(bin2==1)))/length(find(bin2==1));
                    P12_32=length(intersect(find(bin1==3),find(bin2==2)))/length(find(bin2==2));
                    P12_33=length(intersect(find(bin1==3),find(bin2==3)))/length(find(bin2==3));
                    P12_34=length(intersect(find(bin1==3),find(bin2==4)))/length(find(bin2==4));
                    P12_40=length(intersect(find(bin1==4),find(bin2==0)))/length(find(bin2==0));
                    P12_41=length(intersect(find(bin1==4),find(bin2==1)))/length(find(bin2==1));
                    P12_42=length(intersect(find(bin1==4),find(bin2==2)))/length(find(bin2==2));
                    P12_43=length(intersect(find(bin1==4),find(bin2==3)))/length(find(bin2==3));
                    P12_44=length(intersect(find(bin1==4),find(bin2==4)))/length(find(bin2==4));

                    %% Vaxjo

                    deltapair(npa)=abs(P1_0-P2_0*P12_00-P2_1*P12_01-P2_2*P12_02-P2_3*P12_03-P2_4*P12_04) ...
                        /2/sqrt(P2_0*P12_00*P2_1*P12_01*P2_2*P12_02*P2_3*P12_03*P2_4*P12_04) ...
                        +abs(P1_1-P2_0*P12_10-P2_1*P12_11-P2_2*P12_12-P2_3*P12_13-P2_4*P12_14) ...
                        /2/sqrt(P2_0*P12_10*P2_1*P12_11*P2_2*P12_12*P2_3*P12_13*P2_4*P12_14) ...
                        +abs(P1_2-P2_0*P12_20-P2_1*P12_21-P2_2*P12_22-P2_3*P12_23-P2_4*P12_24) ...
                        /2/sqrt(P2_0*P12_20*P2_1*P12_21*P2_2*P12_22*P2_3*P12_23*P2_4*P12_24) ...
                        +abs(P1_3-P2_0*P12_30-P2_1*P12_31-P2_2*P12_32-P2_3*P12_33-P2_4*P12_34) ...
                        /2/sqrt(P2_0*P12_30*P2_1*P12_31*P2_2*P12_32*P2_3*P12_33*P2_4*P12_34) ...
                        +abs(P1_4-P2_0*P12_40-P2_1*P12_41-P2_2*P12_42-P2_3*P12_43-P2_4*P12_44) ...
                        /2/sqrt(P2_0*P12_40*P2_1*P12_41*P2_2*P12_42*P2_3*P12_43*P2_4*P12_44);

                    vaxjo2(nsub,i,j)=deltapair(npa);

                    %% Uncertainty

                    Prob=[P12_00 P12_10 P12_20 P12_30 P12_40];
                    Prob(Prob==0)=[];
                    Entropy0=-sum(Prob.*log2(Prob));
                    Prob=[P12_01 P12_11 P12_21 P12_31 P12_41];
                    Prob(Prob==0)=[];
                    Entropy1=-sum(Prob.*log2(Prob));
                    Prob=[P12_02 P12_12 P12_22 P12_32 P12_42];
                    Prob(Prob==0)=[];
                    Entropy2=-sum(Prob.*log2(Prob));
                    Prob=[P12_03 P12_13 P12_23 P12_33 P12_43];
                    Prob(Prob==0)=[];
                    Entropy3=-sum(Prob.*log2(Prob));
                    Prob=[P12_04 P12_14 P12_24 P12_34 P12_44];
                    Prob(Prob==0)=[];
                    Entropy4=-sum(Prob.*log2(Prob));

                    Entropy(npa)=(Entropy0*P2_0+Entropy1*P2_1+Entropy2*P2_2+Entropy3*P2_3+Entropy4*P2_4);

                    %% Ignition

                    MI(npa)=log2(5)-Entropy(npa);

                    %% Surrogates
                    bin1=binvortexNULL(i,:);
                    bin2=binvortexNULL(j,:);

                    P1_0=length(find(bin1==0))/length(bin1);
                    P1_1=length(find(bin1==1))/length(bin1);
                    P1_2=length(find(bin1==2))/length(bin1);
                    P1_3=length(find(bin1==3))/length(bin1);
                    P1_4=length(find(bin1==4))/length(bin1);
                    P2_0=length(find(bin2==0))/length(bin1);
                    P2_1=length(find(bin2==1))/length(bin1);
                    P2_2=length(find(bin2==2))/length(bin1);
                    P2_3=length(find(bin2==3))/length(bin1);
                    P2_4=length(find(bin2==4))/length(bin1);
                    P12_00=length(intersect(find(bin1==0),find(bin2==0)))/length(find(bin2==0));
                    P12_01=length(intersect(find(bin1==0),find(bin2==1)))/length(find(bin2==1));
                    P12_02=length(intersect(find(bin1==0),find(bin2==2)))/length(find(bin2==2));
                    P12_03=length(intersect(find(bin1==0),find(bin2==3)))/length(find(bin2==3));
                    P12_04=length(intersect(find(bin1==0),find(bin2==4)))/length(find(bin2==4));
                    P12_10=length(intersect(find(bin1==1),find(bin2==0)))/length(find(bin2==0));
                    P12_11=length(intersect(find(bin1==1),find(bin2==1)))/length(find(bin2==1));
                    P12_12=length(intersect(find(bin1==1),find(bin2==2)))/length(find(bin2==2));
                    P12_13=length(intersect(find(bin1==1),find(bin2==3)))/length(find(bin2==3));
                    P12_14=length(intersect(find(bin1==1),find(bin2==4)))/length(find(bin2==4));
                    P12_20=length(intersect(find(bin1==2),find(bin2==0)))/length(find(bin2==0));
                    P12_21=length(intersect(find(bin1==2),find(bin2==1)))/length(find(bin2==1));
                    P12_22=length(intersect(find(bin1==2),find(bin2==2)))/length(find(bin2==2));
                    P12_23=length(intersect(find(bin1==2),find(bin2==3)))/length(find(bin2==3));
                    P12_24=length(intersect(find(bin1==2),find(bin2==4)))/length(find(bin2==4));
                    P12_30=length(intersect(find(bin1==3),find(bin2==0)))/length(find(bin2==0));
                    P12_31=length(intersect(find(bin1==3),find(bin2==1)))/length(find(bin2==1));
                    P12_32=length(intersect(find(bin1==3),find(bin2==2)))/length(find(bin2==2));
                    P12_33=length(intersect(find(bin1==3),find(bin2==3)))/length(find(bin2==3));
                    P12_34=length(intersect(find(bin1==3),find(bin2==4)))/length(find(bin2==4));
                    P12_40=length(intersect(find(bin1==4),find(bin2==0)))/length(find(bin2==0));
                    P12_41=length(intersect(find(bin1==4),find(bin2==1)))/length(find(bin2==1));
                    P12_42=length(intersect(find(bin1==4),find(bin2==2)))/length(find(bin2==2));
                    P12_43=length(intersect(find(bin1==4),find(bin2==3)))/length(find(bin2==3));
                    P12_44=length(intersect(find(bin1==4),find(bin2==4)))/length(find(bin2==4));

                    %% Vaxjo

                    vaxjo2NULL(npa)=abs(P1_0-P2_0*P12_00-P2_1*P12_01-P2_2*P12_02-P2_3*P12_03-P2_4*P12_04) ...
                        /2/sqrt(P2_0*P12_00*P2_1*P12_01*P2_2*P12_02*P2_3*P12_03*P2_4*P12_04) ...
                        +abs(P1_1-P2_0*P12_10-P2_1*P12_11-P2_2*P12_12-P2_3*P12_13-P2_4*P12_14) ...
                        /2/sqrt(P2_0*P12_10*P2_1*P12_11*P2_2*P12_12*P2_3*P12_13*P2_4*P12_14) ...
                        +abs(P1_2-P2_0*P12_20-P2_1*P12_21-P2_2*P12_22-P2_3*P12_23-P2_4*P12_24) ...
                        /2/sqrt(P2_0*P12_20*P2_1*P12_21*P2_2*P12_22*P2_3*P12_23*P2_4*P12_24) ...
                        +abs(P1_3-P2_0*P12_30-P2_1*P12_31-P2_2*P12_32-P2_3*P12_33-P2_4*P12_34) ...
                        /2/sqrt(P2_0*P12_30*P2_1*P12_31*P2_2*P12_32*P2_3*P12_33*P2_4*P12_34) ...
                        +abs(P1_4-P2_0*P12_40-P2_1*P12_41-P2_2*P12_42-P2_3*P12_43-P2_4*P12_44) ...
                        /2/sqrt(P2_0*P12_40*P2_1*P12_41*P2_2*P12_42*P2_3*P12_43*P2_4*P12_44);

                    %%

                    npa=npa+1;

                end
            end
        end

        vaxjo2NULL(isinf(vaxjo2NULL))=[];
        vaxjo2NULL(isnan(vaxjo2NULL))=[];
        vaxjo2NULLall=[vaxjo2NULLall vaxjo2NULL];

        deltapair(deltapair==Inf)=[];
        interferencevaxjo.(g)(nsub)=nanmean(deltapair);
        uncertainty.(g)(nsub)=nanmean(Entropy);
        Ignition.(g)(nsub)=nanmean(MI);
    end

    mNULL=nanmean(vaxjo2NULLall);
    sNULL=nanstd(vaxjo2NULLall);

    for nsub=1:numel(idxs)
        for i=1:NP
            for j=1:NP
                if i~=j
                    Interferencevaxjo_submat(nsub,i,j)=(vaxjo2(nsub,i,j)-mNULL)/sNULL;
                    if isinf(Interferencevaxjo_submat(nsub,i,j))
                        Interferencevaxjo_submat(nsub,i,j)=randn*0.001;
                    end
                    if isnan(Interferencevaxjo_submat(nsub,i,j))
                        Interferencevaxjo_submat(nsub,i,j)=randn*0.001;
                    end
                end
            end
        end
    end

    Interferencevaxjo_sub.(g)=Interferencevaxjo_submat;
end

save results_BANDA_GECfitting.mat FCfitting;

figure(1)
a=(Interferencevaxjo_sub.CNT);
b=(Interferencevaxjo_sub.ANX);
c=(Interferencevaxjo_sub.DEP);
a=nansum(nansum(a,3),2);
b=nansum(nansum(b,3),2);
c=nansum(nansum(c,3),2);
boxplot([a; b; c],[zeros(size(a)); ones(size(b)); 2*ones(size(c))]);
stats=permutation_htest2_np([a',b'],[ones(1,numel(a')), 2*ones(1,numel(b'))],10000,0.01,'ranksum');
min(stats.pvals)
stats=permutation_htest2_np([a',c'],[ones(1,numel(a')), 2*ones(1,numel(c'))],10000,0.01,'ranksum');
min(stats.pvals)
stats=permutation_htest2_np([b',c'],[ones(1,numel(b')), 2*ones(1,numel(c'))],10000,0.01,'ranksum');
min(stats.pvals)

figure(2)
a=(spectral_gap.CNT)';
b=(spectral_gap.ANX)';
c=(spectral_gap.DEP)';
boxplot([a; b; c],[zeros(size(a)); ones(size(b)); 2*ones(size(c))]);
stats=permutation_htest2_np([a',b'],[ones(1,numel(a')), 2*ones(1,numel(b'))],10000,0.01,'ranksum');
min(stats.pvals)
stats=permutation_htest2_np([a',c'],[ones(1,numel(a')), 2*ones(1,numel(c'))],10000,0.01,'ranksum');
min(stats.pvals)
stats=permutation_htest2_np([b',c'],[ones(1,numel(b')), 2*ones(1,numel(c'))],10000,0.01,'ranksum');
min(stats.pvals)

figure(3)
a=(Turbulence.CNT)';
b=(Turbulence.ANX)';
c=(Turbulence.DEP)';
boxplot([a; b; c],[zeros(size(a)); ones(size(b)); 2*ones(size(c))]);
stats=permutation_htest2_np([a',b'],[ones(1,numel(a')), 2*ones(1,numel(b'))],10000,0.01,'ranksum');
min(stats.pvals)
stats=permutation_htest2_np([a',c'],[ones(1,numel(a')), 2*ones(1,numel(c'))],10000,0.01,'ranksum');
min(stats.pvals)
stats=permutation_htest2_np([b',c'],[ones(1,numel(b')), 2*ones(1,numel(c'))],10000,0.01,'ranksum');
min(stats.pvals)


%% Scores

Threshold2=0.18;

alpha_EN=0.01;

scoresdata=readtable('scores_02.csv');

for i = 2:3
    g = fn{i};
    ids = groupIDs.(g);
    InterferenceMatrix=Interferencevaxjo_sub.(g);
    nn=1;
    for n=1:size(ids,1)
        [tf, loc] = ismember(ids(n), string(scoresdata{:,1}));   % loc are indices into subject{}
        if tf==1
            score1(nn)=scoresdata{loc(tf),2};
            score2(nn)=scoresdata{loc(tf),3};
            C=squeeze(InterferenceMatrix(n,:,:));
            Cvec(nn,:)=C(:);
            nn=nn+1;
        end
    end
    Cvec(isnan(Cvec))=randn*0.001;

    for sc=1:2
        if sc==1
            behav_y=score1;
        end
        if sc==2
            behav_y=score2;
        end

        %%
        %% ===================================================
        %  1. Feature selection on entire dataset
        % ====================================================

        % 3) Outer CV to get honest per-subject predictions
        for repe=1:100
            %%
            y = behav_y(:);
            Xall = Cvec;

            for i = 1:size(Xall,2)
                co2(i) = corr(Xall(:,i), y, 'rows','pairwise');
            end

            selected = find(abs(co2) > Threshold2);

            % fprintf("Selected %d features.\n", length(selected));

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

            R2all(repe)=R2;
            Corrall(repe)=Corr;
            yhatall(repe,:)=yhat;
        end
        fprintf('\n========== RESULTS ==========\n');
        fprintf('Case   = %s\n', g);
        fprintf('Factor   = %i\n', sc);
        fprintf('Cross-validated R2   = %.4f\n', mean(R2all));
        fprintf('Cross-validated Corr = %.4f\n', mean(Corrall));
        R2final.g(sc,:)=R2all;
        Corrfinal.g(sc,:)=Corrall;
        yhatfinal.g(sc,:)=mean(yhatall);
    end
end

save results_BANDA_Hopf_Ceff_Vaxjo.mat spectral_gap ...
    Turbulence interferencevaxjo Interferencevaxjo_sub uncertainty Ignition ...
    R2final Corrfinal yhatfinal;
