clear all;
load results_Behaviour_BANDA_emp_pred1y_Vaxjo.mat;
R2final_Vaxjo=R2final;
Corrfinal_Vaxjo=Corrfinal;

load results_Behaviour_BANDA_emp_pred1y_FC.mat;
R2final_FC=R2final;
Corrfinal_FC=Corrfinal;

a=R2final_Vaxjo;
b=R2final_FC;
figure;
violinplot([b a]);
ranksum(a,b)

a=Corrfinal_Vaxjo;
b=Corrfinal_FC;
figure;
violinplot([b a]);
ranksum(a,b)

%%%%%

clear all;
load results_Behaviour_BANDA_emp_Vaxjo.mat;
R2final_Vaxjo=R2final;
Corrfinal_Vaxjo=Corrfinal;

load results_Behaviour_BANDA_emp_FC.mat;
R2final_FC=R2final;
Corrfinal_FC=Corrfinal;

a=R2final_Vaxjo(1,:)';
b=R2final_FC(1,:)';
figure;
violinplot([b a]);
ranksum(a,b)

a=R2final_Vaxjo(2,:)';
b=R2final_FC(2,:)';
figure;
violinplot([b a]);
ranksum(a,b)

a=Corrfinal_Vaxjo(1,:)';
b=Corrfinal_FC(1,:)';
figure;
violinplot([b a]);
ranksum(a,b)

a=Corrfinal_Vaxjo(2,:)';
b=Corrfinal_FC(2,:)';
figure;
violinplot([b a]);
ranksum(a,b)
