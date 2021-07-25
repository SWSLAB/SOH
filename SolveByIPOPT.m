%% ####################################################################################################################
% Code for the paper:
% Optimization of Multi-quality Water Networks: Can simple optimization heuristics compete with nonlinear solvers?
% By Mashor Housh, PhD
% University of Haifa, mhoush@univ.haifa.ac.il
%% ####################################################################################################################
% This code requires:
% YALMIP toolbox: https://yalmip.github.io/
% Optitoolbox: https://inverseproblem.co.nz/OPTI/
% Developed under Matlab 2018b
%% ####################################################################################################################
clc
clear
close all
%%
fprintf('\n################################################\n')
fprintf('######Solving nonlinear problem with IPOPT######\n')
fprintf('################################################\n\n\n')

% rand('seed',10)  % Change the seed for randomization
load('ProblemsData.mat');

Tmax=8;
for T=2:2
    Q=sdpvar(Ntot,T,'full');
    C=sdpvar(Ntot,T,'full');
    M=sdpvar(Ntot,T,'full');
    
    % Defining the constraints
    C1=[A*Q==0];
    C2=[A*M==0];
    C3=[B*C==0];
    C4=[M==Q.*C];
    C5=[Bs*C==Cs{T}];
    C6=[Bd*Q==Qd{T}];
    C7=[Qmin{T}<=Q<=Qmax{T}];
    C8=[Cmin{T}<=C<=Cmax{T}];
    C9=[sum(Bs*Q,2)<=Qsmax_total{T}];
    Const=[C1;C2;C3;C4;C5;C6;C7;C8;C9];
    % Objective function
    cost = sum(sum(f{T}.*Q));

    % Multistart solver
    options=sdpsettings('solver','ipopt');
    options.savesolverinput=1;
    options.ipopt.max_iter=1000000;
    options.ipopt.max_cpu_time = 30*60;
    options.usex0=1;
    loopS = 25;
    for i=1:loopS
        clear mex
        Q0{i,T} = unifrnd(Qmin{T},Qmax{T});
        C0{i,T} = unifrnd(Cmin{T},Cmax{T});
        M0{i,T} = unifrnd(Qmin{T}.*Cmin{T},Qmax{T}.*Cmax{T});
        assign(Q,Q0{i,T});
        assign(C,C0{i,T});
        assign(M,M0{i,T});
        solution{i,T}=optimize(Const,cost,options);
        SolTime(i,T) = solution{i,T}.solvertime;
        SolFlag(i,T) = solution{i,T}.problem;
        Obj(i,T)=double(cost);
        Qval{i,T}=double(Q);
        Cval{i,T}=double(C);
        Mval{i,T}=double(M);
    end
    [Objopt(T), id]=min(Obj(:,T));
    TotSolTime(T) = sum(SolTime(:,T));
    Qopt{T} = Qval{id,T};
    Copt{T} = Cval{id,T};
    Mopt{T} = Mval{id,T};
end
