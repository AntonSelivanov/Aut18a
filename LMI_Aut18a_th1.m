function flag=LMI_Aut18a_th1(A,B,C,Kbar,h,alpha)
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper 
% A. Selivanov and E. Fridman, "An improved time-delay implementation of derivative-dependent feedback,"  
% Automatica, vol. 98, pp. 269-276, 2018. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A, B, C   - parameters of (1) 
% Kbar      - cell array of gains from (6) 
% h         - delay from (6)
% alpha     - decay rate

% Output: 
% flag =1 if feasible, =0 otherwise
%% Decision variables and notations
r=length(Kbar); % relative degree 
[n,m]=size(B); 
P=sdpvar(n); 

D=A+B*Kbar{1}*C; 
for i=1:r-1 
    R{i}=sdpvar(m); %#ok<AGROW> 
    D=D+B*Kbar{i+1}*C*A^i; 
end 
%% LMIs 
M=blkvar; 
M(1,1)=D'*P+P*D+2*alpha*P; 
for i=1:r-2
    M(1,1)=M(1,1)+(i*h/2)^2*(Kbar{i+1}*C*A^(i+1))'*R{i}*(Kbar{i+1}*C*A^(i+1)); 
end
for i=1:r-1
    M(1,1+i)=P*B; 
    M(1+i,1+i)=-exp(-2*alpha*i*h)*R{i}; 
    M(1+i,r+1)=(r-1)*h/2*(Kbar{r}*C*A^(r-1)*B)'*R{r-1}; 
end
M(1,r+1)=(r-1)*h/2*(Kbar{r}*C*A^(r-1)*D)'*R{r-1}; 
M(r+1,r+1)=-R{r-1}; 
M=sdpvar(M); 
%% Solution of LMIs
LMIs=[P>=0,M<=0]; 

options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

flag=0; 
if sol.problem==0
    [primal,~]=check(LMIs); 
    flag=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end
