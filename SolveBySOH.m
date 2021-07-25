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
fprintf('\n#############################################\n')
fprintf('#########Simple Optimization Heuristic####### \n')
fprintf('#############################################\n\n\n')

% rand('seed',10)  % Change the seed for randomization  
load('ProblemsData.mat');

loopS = 25;
beta=1e6;
kmax=500;
eps=0.001;
Tmax=8;

Q0=cell(loopS,Tmax,kmax+1);
C0=cell(loopS,Tmax,kmax+1);
Obj0=nan*ones(loopS,Tmax,kmax+1);
RelInfeasiblity=nan*ones(loopS,Tmax,kmax+1);
iter_time=zeros(loopS,Tmax,kmax+1);

QL=cell(loopS,Tmax,kmax+1);
CL=cell(loopS,Tmax,kmax+1);
dQ=cell(loopS,Tmax,kmax+1);
dC=cell(loopS,Tmax,kmax+1);
alphaopt=nan*ones(loopS,Tmax,kmax+1);
eObj=nan*ones(loopS,Tmax,kmax+1);
eQ=nan*ones(loopS,Tmax,kmax+1);
BalMaxLinInfeas=nan*ones(loopS,Tmax,kmax+1);
BalMaxNonlinInfeas=nan*ones(loopS,Tmax,kmax+1);


for T=1:Tmax
              
    % Define Fixed Parts of Optimization Problem
    Q=sdpvar(Ntot,T,'full');
    C=sdpvar(Ntot,T,'full');
    v=sdpvar(1,1);
    
    C1=[A*Q==0];
    C2=[];%C2=[A*(Q.*C)==0] This constrain is ommited, since it will be linearized
    C3=[B*C==0];
    C4=[Bs*C==Cs{T}];
    C5=[Bd*Q==Qd{T}];
    C6=[Qmin{T}<=Q<=Qmax{T}];
    C7=[Cmin{T}<=C<=Cmax{T}];
    C8=[sum(Bs*Q,2)<=Qsmax_total{T}];
    Const0=[C1;C2;C3;C4;C5;C6;C7;C8;];
    Obj=sum(sum(f{T}.*Q))+beta*v;
    
    options=sdpsettings('solver','cplex');
    options.verbose=0;
    
    % Multistart Optimizaion
    for i=1:loopS
        % Start Successive Linearization
        k=1;
        while 1
            clear mex
            % Solve Linear Approximation
            if k==1
                % Initial Guess
                Q0{i,T,k}=unifrnd(Qmin{T},Qmax{T});
                C0{i,T,k}=unifrnd(Cmin{T},Cmax{T});
                Obj0(i,T,k)=inf;
                RelInfeasiblity(i,T,k)=inf;                
            end
            Const=[Const0; -v<=A*(C0{i,T,k}.*Q+Q0{i,T,k}.*C-Q0{i,T,k}.*C0{i,T,k})<=v;];
            sol=solvesdp(Const,Obj,options);
            
            % Minimize the nonlinear infeasiblity while keeping linear feasiblity of
            % the problem. Using Golden Section single-variable optimization.
            QL{i,T,k}=double(Q);
            CL{i,T,k}=double(C);
            dQ{i,T,k}=QL{i,T,k}-Q0{i,T,k};
            dC{i,T,k}=CL{i,T,k}-C0{i,T,k};
            tic
            [alphaopt(i,T,k), fopt]=fminbnd(@(alpha)funeval(alpha,Q0{i,T,k},C0{i,T,k},dQ{i,T,k},dC{i,T,k},f{T},eps,beta,A),0,1);
            GS_time=toc;
            
            Q0{i,T,k+1}=Q0{i,T,k}+alphaopt(i,T,k)*dQ{i,T,k};
            C0{i,T,k+1}=C0{i,T,k}+alphaopt(i,T,k)*dC{i,T,k};
            Obj0(i,T,k+1)=fopt;
            
            % Store time
            iter_time(i,T,k)=sol.solvertime+GS_time;
            
            % Errors
            % Error in Objective
            eObj(i,T,k)=norm(Obj0(i,T,k+1)-Obj0(i,T,k))/Obj0(i,T,k);
            % Error in flow compared to previous iterations
            ex_tmp=zeros(k,1);
            for j=1:k
                ex_tmp(j)=norm(Q0{i,T,k+1}(:)-Q0{i,T,j}(:))/norm(Q0{i,T,j}(:));
            end
            eQ(i,T,k)=min(ex_tmp);
            % Maximum Linear Infeasiblity in Solute Balance
            BalMaxLinInfeas(i,T,k)=double(v);
            % Maximum Nonlinear Infeasiblity in Solute Balance
            BalMaxNonlinInfeas(i,T,k)=max(max(abs(A*(Q0{i,T,k+1}.*C0{i,T,k+1}))));
            % Relative Nonlinear Infeasiblity in Solute balance
            [~,RelInfeasiblity(i,T,k)]=funeval(alphaopt(i,T,k),Q0{i,T,k},C0{i,T,k},dQ{i,T,k},dC{i,T,k},f{T},eps,beta,A);
            % Report
            fprintf('#################### Iteration %d ########################\n',k)
            fprintf('Current Penalized Objective Value is: %f \n', Obj0(i,T,k+1))
            fprintf('Current nonPenalized Objective Value is: %f \n', sum(sum(f{T}.*Q0{i,T,k+1})))
            fprintf('Current Max Linear Infeasiblity is: %f \n',BalMaxLinInfeas(i,T,k))
            fprintf('Current Max NonLin Infeasiblity is: %f \n', BalMaxNonlinInfeas(i,T,k))
            fprintf('Current Relative Error in Q is: %f \n', eQ(i,T,k))
            fprintf('Current Relative Error in Objective is: %f \n', eObj(i,T,k))
            fprintf('Current Relative Infeasiblity in salt balance: %f \n', RelInfeasiblity(i,T,k))
            
            % Check Stopping Criteria
            if  eQ(i,T,k)<eps || k>kmax
                lag=k-find(ex_tmp(1:k)<eps,1);
                if ~isempty(lag)
                    Lags(i,T)=lag;
                else
                    Lags(i,T)=nan;
                end
                fprintf('The stopping criteria is obtained in Lag %d \n',Lags(i,T))
                break
            else
                k=k+1;
            end
        end
        % Store the solution after convergence
        Objval(i,T)=Obj0(i,T,k+1);
        ObjvalOrig(i,T)=sum(sum(f{T}.*Q0{i,T,k+1}));
        Infeasablity(i,T)=RelInfeasiblity(i,T,k);
        artvars(i,T)=BalMaxLinInfeas(i,T,k);
        MaxInfeasablity(i,T)=BalMaxNonlinInfeas(i,T,k);
        reqk(i,T)=k;
        SolTime(i,T)=sum(iter_time(i,T,:));
        Qval{i,T}=Q0{i,T,k+1};
        Cval{i,T}=C0{i,T,k+1};
    end
    % Choose best solution from multiple runs
    [Objopt(T), id(T)]=min(Objval(:,T));
    TotSolTime(T) = sum(SolTime(:,T));
    Qopt{T} = Qval{id(T),T};
    Copt{T} = Cval{id(T),T};
    InfeasablityOpt(T)=Infeasablity(id(T),T);
    MaxInfeasablityOpt(T)=MaxInfeasablity(id(T),T);
    artvarsOpt(T)=artvars((T),T);
    reqkOpt(T)=reqk(id(T),T);
    ObjvalOrigOpt(T)=ObjvalOrig(id(T),T);
end

function [f,Infeasiblity]=funeval(alpha,Q0,C0,dQ,dC,f,eps,beta,A)
Q=Q0+alpha*dQ;
C=C0+alpha*dC;
tmp0=A*(Q.*C);
tmp1=Q.*C;
Infeasiblity=norm(tmp0(:))/norm(tmp1(:));
f=sum(sum(f.*Q))+beta*max(Infeasiblity-eps,0);
end
