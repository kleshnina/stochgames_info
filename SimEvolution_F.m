function [pi,coop,Str]=SimEvolution_F(q,b1,b2,c,N,beta,epsi)

% Simulating the evolutionary dynamics for the population with FULL information
% [W,w1,w2,pi,coop]=SimEvolution_F(q,b1,b2,c,N1,N2,beta,epsi);
% INPUT:
% q=[q12, q11, q10, q02, q01, q00], where qij is the transition
% probability to move to State 1 if currently in State i and if j players
% have cooperated. 
% b1, b2, c .. game parameters
% N .. population size
% beta .. strength of selection
% epsi .. error rate for implementation errors%
% OUTPUT:
% pi=[pi1,pi2] .. average payoffs in the two populations
% coop=[coop1,coop2] .. respective average cooperation rates
% f1 .. frequency with which players are in state 1

%% Preparations
nStr=2^8; % Initiating all memory-1 strategies
Str=zeros(nStr,8); 
for i=1:nStr
    Str(i,:)=sscanf(dec2bin(i-1,8), '%1d' )';
end

%% Calculating the payoff matrices
Pi=zeros(1,nStr); % Pairwise payoff matrix 
Co=zeros(1,nStr); % Cooperation matrix
for i=1:nStr
    [T,v,Pai,Coop,F1]=CalcPayoff(Str(i,:),Str(i,:),q,epsi,b1,b2,c);
    % Calculating players' payoffs
    Pi(1,i)=Pai(1); 
    Co(1,i)=Coop(1);
end

%% Constructing the transition matrix
%% Off-diagonal entries
M=zeros(nStr,nStr); 
for iOld=1:nStr
     for iNew=1:nStr 
         if iNew~=iOld
             M(iOld,iNew) = 1/2*1/nStr*CalcRho(Pi(iOld),Pi(iNew),N,beta); 
         else
             M(iOld,iNew) =0; 
         end
     end
end

%% Diagonal entries
for i=1:nStr
    M(i,i)=1-sum(M(i,:)); 
end
            
%% Invariant distribution
w=null(M'-eye(nStr));
w=w'/sum(w);

%% Other outputs
pi=sum(sum(w.*Pi)); 
coop=sum(sum(w.*Co)); 
end


function rho=CalcRho(PiOld,PiNew,N,beta)
% Calculates the fixation probability of one S1 mutant in an S2 population
alpha=exp(-beta*(PiNew-PiOld));
rho=1/(1+sum(cumprod(alpha*ones(1,N-1)))); 
% Calculating the fixation probability according to formula given in SI
end
