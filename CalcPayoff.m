function [T,v,pi,coop,f1]=CalcPayoff(p,pt,q,epsi,b1,b2,c)

% [T,v,pi,coop,f1]=CalcPayoff(p,pt,q,epsi,b1,b2,c)
% INPUT:
% p=[p1CC,p1CD,...,p2DD] .. strategy of player 1 (fully informed player)
% pt=[pt1CC,...,pt2DD] .. strategy of player 2 (non-informed player:
% pt1CC=pt2CC, ... , pt1DD=pt2DD). 
% q=[q1CC,q1CD,...,q2DD] .. state transitions (probability to move to state
% 1 given previous state and outcome)
% epsi .. error rate
% b1,b2,c .. payoff parameters
% OUTPUT:
% T .. transition matrix for the respective game
% v .. invariant distribution 
% pi=[pi1,pi2] .. Resulting payoffs for the two players
% coop=[coop1,coop2] .. Respective cooperation rates
% f1 .. frequency that the two players will be in state 1



%% (1) Constructing the transition matrix %%
p=(1-epsi)*p+epsi*(1-p); 
pt=(1-epsi)*pt+epsi*(1-pt);
q=[q(1) q(2) q(2) q(3) q(4) q(5) q(5) q(6)]; 
T=[q(1)*p(1)*pt(1), q(1)*p(1)*(1-pt(1)), q(1)*(1-p(1))*pt(1), q(1)*(1-p(1))*(1-pt(1)), (1-q(1))*p(5)*pt(5), (1-q(1))*p(5)*(1-pt(5)), (1-q(1))*(1-p(5))*pt(5), (1-q(1))*(1-p(5))*(1-pt(5)); 
   q(2)*p(2)*pt(3), q(2)*p(2)*(1-pt(3)), q(2)*(1-p(2))*pt(3), q(2)*(1-p(2))*(1-pt(3)), (1-q(2))*p(6)*pt(7), (1-q(2))*p(6)*(1-pt(7)), (1-q(2))*(1-p(6))*pt(7), (1-q(2))*(1-p(6))*(1-pt(7));
    q(3)*p(3)*pt(2), q(3)*p(3)*(1-pt(2)), q(3)*(1-p(3))*pt(2), q(3)*(1-p(3))*(1-pt(2)), (1-q(3))*p(7)*pt(6), (1-q(3))*p(7)*(1-pt(6)), (1-q(3))*(1-p(7))*pt(6), (1-q(3))*(1-p(7))*(1-pt(6));
    q(4)*p(4)*pt(4), q(4)*p(4)*(1-pt(4)), q(4)*(1-p(4))*pt(4), q(4)*(1-p(4))*(1-pt(4)), (1-q(4))*p(8)*pt(8), (1-q(4))*p(8)*(1-pt(8)), (1-q(4))*(1-p(8))*pt(8), (1-q(4))*(1-p(8))*(1-pt(8));
    q(5)*p(1)*pt(1), q(5)*p(1)*(1-pt(1)), q(5)*(1-p(1))*pt(1), q(5)*(1-p(1))*(1-pt(1)), (1-q(5))*p(5)*pt(5), (1-q(5))*p(5)*(1-pt(5)), (1-q(5))*(1-p(5))*pt(5), (1-q(5))*(1-p(5))*(1-pt(5));
    q(6)*p(2)*pt(3), q(6)*p(2)*(1-pt(3)), q(6)*(1-p(2))*pt(3), q(6)*(1-p(2))*(1-pt(3)), (1-q(6))*p(6)*pt(7), (1-q(6))*p(6)*(1-pt(7)), (1-q(6))*(1-p(6))*pt(7), (1-q(6))*(1-p(6))*(1-pt(7));
    q(7)*p(3)*pt(2), q(7)*p(3)*(1-pt(2)), q(7)*(1-p(3))*pt(2), q(7)*(1-p(3))*(1-pt(2)), (1-q(7))*p(7)*pt(6), (1-q(7))*p(7)*(1-pt(6)), (1-q(7))*(1-p(7))*pt(6), (1-q(7))*(1-p(7))*(1-pt(6));
    q(8)*p(4)*pt(4), q(8)*p(4)*(1-pt(4)), q(8)*(1-p(4))*pt(4), q(8)*(1-p(4))*(1-pt(4)), (1-q(8))*p(8)*pt(8), (1-q(8))*p(8)*(1-pt(8)), (1-q(8))*(1-p(8))*pt(8), (1-q(8))*(1-p(8))*(1-pt(8))];

%% (2) Find the invariant distribution
v=null(T'-eye(8));
v=v(:,1)'/sum(v(:,1));

%% (3) Compute payoffs, cooperation rates, frequency of state 1
f1=sum(v(1:4)); 
coop1=v(1)+v(5)+v(2)+v(6); 
coop2=v(1)+v(5)+v(3)+v(7); 
coop=[coop1, coop2]; 
u1=[b1-c -c b1 0 b2-c -c b2 0]; 
u2=[b1-c b1 -c 0 b2-c b2 -c 0]; 
pi1=v*u1'; 
pi2=v*u2'; 
pi=[pi1,pi2]; 
end