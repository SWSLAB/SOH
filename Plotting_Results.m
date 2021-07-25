%% ####################################################################################################################
% Code for the paper:
% Optimization of Multi-quality Water Networks: Can simple optimization heuristics compete with nonlinear solvers?
% By Mashor Housh, PhD
% University of Haifa, mhoush@univ.haifa.ac.il
%% ####################################################################################################################
% This code plots Figures 3-11 in the paper
clc
clear
close all

load('IPOPT_Results.mat')
Tvec=1:8;
NonlinOpt=Objopt;
NonlinTime=TotSolTime;
NonlinAllObj=Obj;
NonlinAllTime=SolTime;
NonlinAllFlag=SolFlag;

load('SOH_Results.mat')
LinOpt=Objopt;
LinTime=TotSolTime;
LinAllObj=Objval;
LinAllTime=SolTime;

load('ProblemsData.mat')

%% Figure 3 in the paper
figure
Tplot=1;
K=reqkOpt(Tplot);
objectiveiter=[];
for i=1:K
    objectiveiter(i)=sum(sum(f{Tplot}.*Q0{id(Tplot),Tplot,i+1}));
end
yyaxis left
plot(objectiveiter,'LineWidth',1.5)
ylabel('Objective Value (M$)')
yyaxis right
plot(100*squeeze(RelInfeasiblity(id(Tplot),Tplot,1:K)),'LineWidth',1.5)
hold all
plot(100*squeeze(eQ(id(Tplot),Tplot,1:K)),'LineStyle',':','LineWidth',1.5)
ylabel('Relative Error (%)')
xlabel('Iteration')
xlim([0 120])
set(gca,'YScale','log')
legend({'Objective Function','Relative Infeasibility','Stopping Criterion'},'FontSize',12)
set(gca,'FontSize',12)

%% Figure 4 in the paper

% Figure 4
figure
plot(Tvec(1:5),NonlinTime(1:5)/25,'LineWidth',1.5)
hold all
plot(Tvec(1:5),LinTime(1:5)/25,'LineWidth',1.5)
set(gca,'FontSize',12)
set(gca,'XTick',[1:5])
xlabel('T (years)')
ylabel('Average Solver Time (sec)')
legend({'Nonlinear Solver (IPOPT)','Optimization Heuristic'},'FontSize',12)

%% Figure 5 in the paper
% Figure 5a
figure
plot(Tvec,NonlinTime/25,'LineWidth',1.5)
hold all
plot(Tvec,LinTime/25,'LineWidth',1.5)
set(gca,'FontSize',12)
xlabel('T (years)')
ylabel('Average Solver Time (sec)')
plot([0 9], [30*60 30*60],'-g','LineWidth',1.5)
xlim([1 8])
ylim([0 1900])
legend({'Nonlinear Solver (IPOPT)','Optimization Heuristic','Solver Maximum Time Setting'},'FontSize',12)

% Figure 5b
figure
bar(Tvec,sum(NonlinAllFlag==0))
set(gca,'FontSize',12)
xlabel('T (years)')
ylabel('Declared Successful Runs')

%% Figure 6 in the paper
% Figure 6a
figure
plot([0 9], [30*60 30*60],'-g')
hold all
boxplot(min(NonlinAllTime,30*60),'Widths',0.6,'Whisker',inf)
plot(Tvec,mean(NonlinAllTime),'or')
text(0.5,1870,'Solver Maximum Time Setting')
set(gca,'FontSize',12)
ylabel('Solver Time (sec)')
xlabel('T (years)')

% Figure 6b
figure
plot([0 9], [30*60 30*60],'-g')
hold all
boxplot(min(LinAllTime,30*60),'Widths',0.6,'Whisker',inf)
plot(Tvec,mean(LinAllTime),'or')
set(gca,'FontSize',12)
ylabel('Solver Time (sec)')
xlabel('T (years)')

%% Figure 7 in the paper
figure
subplot(1,2,1)
NonlinAllTimeplot=min(NonlinAllTime,30*60);
histogram(NonlinAllTimeplot(:,1),10);
xlabel('Solver Time (sec)')
ylabel('Count (out of 25)')
ylim([0 25])
set(gca,'FontSize',12)
subplot(1,2,2)
histogram(NonlinAllTimeplot(:,8),10);
xlabel('Solver Time (sec)')
ylim([0 25])
set(gca,'FontSize',12)

%% Figure 8 in the paper
figure
subplot(1,2,1)
histogram(LinAllTime(:,1),10);
xlabel('Solver Time (sec)')
ylabel('Count (out of 25)')
ylim([0 25])
set(gca,'FontSize',12)
subplot(1,2,2)
histogram(LinAllTime(:,8),10);
xlabel('Solver Time (sec)')
ylim([0 25])
set(gca,'FontSize',12)

%% Figure 9 in the paper
figure
plot(Tvec,(LinOpt-NonlinOpt)./NonlinOpt*100,'LineWidth',1.5)
set(gca,'FontSize',12)
xlabel('T (years)')
ylabel('Objective Error (%)')

%% Figure 10 in the paper
figure
plot(Tvec,InfeasablityOpt*100,'LineWidth',1.5)
hold all
plot(Tvec,ones(1,8)*eps*100,'--r')
ylabel('Relative Solute Balance Infeasibility (%)')
xlabel('T (years)')
ylim([0.0999, 0.1001])
legend({'Resulted Infeasibility','Infeasibility Tolerance'},'FontSize',12)
set(gca,'FontSize',12)

%% Figure 11 in the paper
figure
for i=1:8
    MaxInfeas(i)=max(max(abs(A*(Qopt{i}.*Copt{i}))));
end

plot(Tvec,MaxInfeas,'LineWidth',1.5)
set(gca,'FontSize',12)
ylabel('Maximum Solute Balance Infeasiblity (kg)')
xlabel('T (years)')



