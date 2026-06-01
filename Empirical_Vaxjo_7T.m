clear all;
path2=[ '../../Turbulence/Basics'];
addpath(genpath(path2));
path3=[ '../../Turbulence/SC_longrange'];
addpath(genpath(path3));
path4=[ '../../TNonequilibrium'];
addpath(genpath(path3));

load (['hcp7t_rfMRI_REST1_PA_schaefer1000.mat']);
load results_model_vortex_rest_7T.mat;

load ('schaefer1000to100.mat');
partitiondef=schaefer1000to100;
NP=100;

NSUB=182;
N=1000;
Tmax=900;

Isubdiag = find(tril(ones(NP),-1));

for i=1:NP
    partition{i}=find(partitiondef==i);
end

%%
% Parameters of the data
TR=1;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

%
for nsub=1:NSUB
    wC=squeeze(Ceff(nsub,:,:));
    [V D]=eig(wC);
    D=diag(real(D));
    Dall=D;
    [dmax]=sort(D,'descend');
    spacings = abs(dmax);
    gapThreshold = mean(spacings) + std(spacings);
    gapIdx = find(spacings > gapThreshold);
    nugap(nsub)=numel(gapIdx);
end
NumGap=ceil(mean(nugap));

for nsub=1:NSUB
    nsub
    clear signal_filt Phases;
    ts=subject{nsub}.schaeferts;
    ts=ts(:,1:Tmax);
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        if sum(isnan(ts(seed,:)))==0
            signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end
    end
    ts=signal_filt(:,50:end-50);
    Phases=Phases(:,50:end-50);

    % Define Partition

    for i=1:NP
        Vorticity(i,:)=abs(sum(complex(cos(Phases(partition{i},:)),sin(Phases(partition{i},:))))/length(partition{i}));
    end
    Turbulence(nsub)=std(Vorticity(:));

    Tmax2=size(Vorticity,2);
    NullVorticity=Vorticity(randperm(NP),randperm(Tmax2));

    %% Spectralgap
    wC=squeeze(Ceff(nsub,:,:));
    [V D]=eig(wC);
    D=diag(real(D));
    Dall=D;
    [dmax]=sort(D,'descend');
    spacings = abs(dmax);
    spacingsord=sort(spacings,'descend');
    sumgap=0;
    for g = 1:NumGap
        sumgap=sumgap+spacingsord(g);
    end
    spectral_gap(nsub)=sumgap;

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

                Entropy_sub(nsub,i,j)=Entropy(npa);

                %% Ignition

                MI(npa)=log2(5)-Entropy(npa);
                Ignition_sub(nsub,i,j)=MI(npa);

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
    
    vaxjo2NULL(vaxjo2NULL==Inf)=[];
    mvaxjoNULL(nsub)=nanmean(vaxjo2NULL);
    svaxjoNULL(nsub)=nanstd(vaxjo2NULL);

    deltapair(deltapair==Inf)=[];
    interferencevaxjo(nsub)=nanmean(deltapair);
    uncertainty(nsub)=nanmean(Entropy);
    Ignition(nsub)=nanmean(MI);
end


mNULL=nanmean(mvaxjoNULL);
sNULL=nanmean(svaxjoNULL);


for nsub=1:NSUB
    for i=1:NP
        for j=1:NP
            if i~=j
                Interferencevaxjo_sub(nsub,i,j)=(vaxjo2(nsub,i,j)-mNULL)/sNULL;
                if isinf(Interferencevaxjo_sub(nsub,i,j))
                    Interferencevaxjo_sub(nsub,i,j)=randn*0.001;
                end
                if isnan(Interferencevaxjo_sub(nsub,i,j))
                    Interferencevaxjo_sub(nsub,i,j)=randn*0.001;
                end
            end
        end
    end
end

figure(1)
scatter(interferencevaxjo,Ignition);
[cc pp]=corrcoef(interferencevaxjo,Ignition);
cc(2)
pp(2)

figure(2)
scatter(spectral_gap,Ignition);
[cc pp]=corrcoef(spectral_gap,Ignition);
cc(2)
pp(2)

figure(3)
scatter(interferencevaxjo,Turbulence);
[cc pp]=corrcoef(interferencevaxjo,Turbulence);
cc(2)
pp(2)

figure(4)
scatter(interferencevaxjo,spectral_gap);
[cc pp]=corrcoef(interferencevaxjo,spectral_gap);
cc(2)
pp(2)

save results_Empirical_Vaxjo_7T.mat interferencevaxjo uncertainty Turbulence spectral_gap Ignition ...
    Ignition_sub Entropy_sub Interferencevaxjo_sub vaxjo2;


%% EXTRA
for nsub=1:NSUB
    for i=1:NP
        for j=1:NP
            if isinf(vaxjo2(nsub,i,j))
                    vaxjo2(nsub,i,j)=0;
                end
                if isnan(vaxjo2(nsub,i,j))
                    vaxjo2(nsub,i,j)=0;
                end
        end
    end
end

load results_model_vortex_rest_7T.mat;
VAXTOT=squeeze(mean(vaxjo2));
for nsub=1:NSUB
    C=squeeze(Ceff(nsub,:,:));
    vax=squeeze(vaxjo2(nsub,:,:));    
    corrsub(nsub)=corr2(C(Isubdiag),vax(Isubdiag));
    corrvaxsub(nsub)=corr2(VAXTOT(Isubdiag),vax(Isubdiag));
end

figure
boxplot([corrsub' corrvaxsub']);


for nsub=1:NSUB
    nsub
    clear signal_filt Phases;
    ts=subject{nsub}.schaeferts;
    ts=ts(:,1:Tmax);
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        if sum(isnan(ts(seed,:)))==0
            signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end
    end
    ts=signal_filt(:,50:end-50);
    Phases=Phases(:,50:end-50);

    % Define Partition

    for i=1:NP
        Vorticity(i,:)=abs(sum(complex(cos(Phases(partition{i},:)),sin(Phases(partition{i},:))))/length(partition{i}));
    end
    FCV(nsub,:,:)=corrcoef(Vorticity');
end

FCVTOT=squeeze(mean(FCV));
for nsub=1:NSUB
    C=squeeze(Ceff(nsub,:,:));
    fc=squeeze(FCV(nsub,:,:));    
    corrsubfcv(nsub)=corr2(C(Isubdiag),fc(Isubdiag));
    corrvaxsubfcv(nsub)=corr2(FCVTOT(Isubdiag),fc(Isubdiag));
end

figure
boxplot([corrsubfcv' corrvaxsubfcv']);

for t=1:size(Vorticity,2)
    co(t)=corr2(Vorticity(:,t),squeeze(mean(squeeze(Ceff(nsub,:,:))))');
end