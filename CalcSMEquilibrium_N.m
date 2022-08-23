function [coop,pay,s1,freq,Str]=CalcSMEquilibrium_N(qvec,piRound,beta,eps)

% Simulating the evolutionary dynamics for the population with NO information
% [coop,pay,s1,freq,Str]=CalcSMEquilibrium_N(qvec,piRound,beta,eps);
% Calculates the cooperation rate coop, payoff pay, frequency to be in state 1 s1,
% frequency of all memory-1 strategies freq, in the selection-mutation
% equilibrium of the dynamical process, using the method of Fudenberg and Imhof. 
% Input: qvec=[q12, q11, q10, q02, q01, q00], where qij is the transition
% probability to move to State 1 if currently in State i and if j players
% have cooperated. 
% piRound=[u1CC,u1CD,u1DC,u1DD,u2CC,u2CD,u2DC,u2DD] ... One-shot payoffs
% depending on current state and on players' actions
% beta .. strength of selection
% eps .. error rate for implementation errors

%% Setting up all objects
N=100; % Population size
pv1=piRound; % payoff vector from the perspective of player 1
pv2=piRound; pv2(2:3)=piRound(3:-1:2); pv2(6:7)=piRound(7:-1:6); % Creating the payoff vector from the perspective of player 2
Str=zeros(2^4,4); 
for i=1:2^4, Str(i,:)=sscanf(dec2bin(i-1,4), '%1d' )'; end % List of all memory-1 strategies
Str8=[Str(:,1)'; Str(:,2)'; Str(:,3)'; Str(:,4)'; Str(:,1)'; Str(:,2)'; Str(:,3)'; Str(:,4)']'; 
% Converting the 4-dim strategy vectors into 8-dim strategy vectors 
PayM=zeros(2^4,2^4); C=zeros(2^4,2^4); S1=zeros(2^4,2^4); % Initializing the pairwise payoff matrix, the cooperation matrix, and a matrix for the frequency of being in State 1
for i=1:2^4
    for j=i:2^4
        [pi1,pi2,cop1,cop2,s1]=payoff(Str8(i,:),Str8(j,:),qvec,pv1,pv2,eps);
        PayM(i,j)=pi1; PayM(j,i)=pi2; C(i,j)=cop1; C(j,i)=cop2;  S1(i,j)=s1;
        % Calculating and storing all pairwise quantities
    end
end

%% Setting up the transition matrix according to Fudenberg and Imhof, 
%% and calculating the invariant distribution
T=zeros(2^4,2^4); % Initializing the transition matrix
for i=1:2^4
    for j=1:2^4
        T(i,j)=1/(2^4-1)*CalcRho(j,i,PayM,N,beta);
        % Off-diagonal entries are the product of the respective mutant
        % strategy being chosen, and the respective strategy reaching
        % fixation. 
    end
    T(i,i)=0; T(i,i)=1-sum(T(i,:)); 
    %On-diagonal entries are chosen such that each row sums up to 1
end
v=null(T'-eye(2^4)); freq=v'/sum(v); % Calculating the unique normalized Eigenvector with respect to EV 1
coop=freq*diag(C); % Calculating the average cooperation rate
pay=freq*diag(PayM); % Calculating the average payoff
s1=freq*diag(S1); % Calculating the average frequency to be in State 1. 
end



function Rho=CalcRho(S1,S2,PayM,N,beta)
% Calculates the fixation probability of one S1 mutant in an S2 population
alpha=zeros(1,N-1);
for j=1:N-1 % j.. Number of mutants in the population
    pi1=(j-1)/(N-1)*PayM(S1,S1)+(N-j)/(N-1)*PayM(S1,S2); % Payoff mutant
    pi2=j/(N-1)*PayM(S2,S1)+(N-j-1)/(N-1)*PayM(S2,S2); % Payoff resident 
    alpha(j)=exp(-beta*(pi1-pi2));
end
Rho=1/(1+sum(cumprod(alpha))); % Calculating the fixation probability according to formula given in SI
end



function [pi1,pi2,cop1,cop2,s1]=payoff(p,q,qvec,piv1,piv2,eps)
p=p*(1-eps)+(1-p)*eps; q=q*(1-eps)+(1-q)*eps; % Adding errors to the players' strategies
M=[qvec(1)*p(1)*q(1), qvec(1)*p(1)*(1-q(1)), qvec(1)*(1-p(1))*q(1), qvec(1)*(1-p(1))*(1-q(1)), (1-qvec(1))*p(5)*q(5), (1-qvec(1))*p(5)*(1-q(5)), (1-qvec(1))*(1-p(5))*q(5), (1-qvec(1))*(1-p(5))*(1-q(5));
   qvec(2)*p(2)*q(3), qvec(2)*p(2)*(1-q(3)), qvec(2)*(1-p(2))*q(3), qvec(2)*(1-p(2))*(1-q(3)), (1-qvec(2))*p(6)*q(7), (1-qvec(2))*p(6)*(1-q(7)), (1-qvec(2))*(1-p(6))*q(7), (1-qvec(2))*(1-p(6))*(1-q(7));
   qvec(2)*p(3)*q(2), qvec(2)*p(3)*(1-q(2)), qvec(2)*(1-p(3))*q(2), qvec(2)*(1-p(3))*(1-q(2)), (1-qvec(2))*p(7)*q(6), (1-qvec(2))*p(7)*(1-q(6)), (1-qvec(2))*(1-p(7))*q(6), (1-qvec(2))*(1-p(7))*(1-q(6));
   qvec(3)*p(4)*q(4), qvec(3)*p(4)*(1-q(4)), qvec(3)*(1-p(4))*q(4), qvec(3)*(1-p(4))*(1-q(4)), (1-qvec(3))*p(8)*q(8), (1-qvec(3))*p(8)*(1-q(8)), (1-qvec(3))*(1-p(8))*q(8), (1-qvec(3))*(1-p(8))*(1-q(8));
   qvec(4)*p(1)*q(1), qvec(4)*p(1)*(1-q(1)), qvec(4)*(1-p(1))*q(1), qvec(4)*(1-p(1))*(1-q(1)), (1-qvec(4))*p(5)*q(5), (1-qvec(4))*p(5)*(1-q(5)), (1-qvec(4))*(1-p(5))*q(5), (1-qvec(4))*(1-p(5))*(1-q(5));
   qvec(5)*p(2)*q(3), qvec(5)*p(2)*(1-q(3)), qvec(5)*(1-p(2))*q(3), qvec(5)*(1-p(2))*(1-q(3)), (1-qvec(5))*p(6)*q(7), (1-qvec(5))*p(6)*(1-q(7)), (1-qvec(5))*(1-p(6))*q(7), (1-qvec(5))*(1-p(6))*(1-q(7));
   qvec(5)*p(3)*q(2), qvec(5)*p(3)*(1-q(2)), qvec(5)*(1-p(3))*q(2), qvec(5)*(1-p(3))*(1-q(2)), (1-qvec(5))*p(7)*q(6), (1-qvec(5))*p(7)*(1-q(6)), (1-qvec(5))*(1-p(7))*q(6), (1-qvec(5))*(1-p(7))*(1-q(6));
   qvec(6)*p(4)*q(4), qvec(6)*p(4)*(1-q(4)), qvec(6)*(1-p(4))*q(4), qvec(6)*(1-p(4))*(1-q(4)), (1-qvec(6))*p(8)*q(8), (1-qvec(6))*p(8)*(1-q(8)), (1-qvec(6))*(1-p(8))*q(8), (1-qvec(6))*(1-p(8))*(1-q(8))];
% Constructing the transition matrix M for the Markov chain of the game dynamics
v=null(M'-eye(8)); v=v/sum(v); % Calculating the normalized left eigenvector with respect to EV 1. 
pi1=piv1*v; pi2=piv2*v; % Calculating expected payoffs of the two players
cop1=v(1)+v(2)+v(5)+v(6); cop2=v(1)+v(3)+v(5)+v(7); % Calculating the cooperation frequency of the two players
s1=sum(v(1:4)); % Calculating how often players are in the first state
end
