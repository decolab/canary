clear all;
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Tenet/TENET/'];
addpath(genpath(path3));
path4=[ '../../TaskResrvoir/DataHCP100ordered'];
addpath(genpath(path4));
path5=[ '../../Turbulence/Basics'];
addpath(genpath(path5));
path6=[ '../../Turbulence/SC_longrange'];
addpath(genpath(path6));
path7=[ '../../VortexModel/PartitionModel'];
addpath(genpath(path7));

N=1000;

load results_SC_surro.mat;

for trials=1:20
    Csurro=squeeze(C_surro2(trials,:,:));
    %% Spectral gap and entangelment
    [V D]=eig(Csurro);
    D=diag(real(D));
    Dall=D;
    [dmax]=sort(D,'descend');
    spacings = abs(dmax);
    gapThreshold = mean(spacings) + std(spacings);
    gapIdx = find(spacings > gapThreshold); 
    sumgap=0;
    if ~isempty(gapIdx)
        for g = 1:numel(gapIdx)
            sumgap=sumgap+spacings(gapIdx(g));
        end
    end
    spectral_gap(trials)=sumgap;
end

boxplot(spectral_gap);