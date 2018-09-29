function flag=LMI_Aut18a_th2(A,B,C,Kbar,h,alpha)
% This MATLAB program checks the feasibility of LMIs from Theorem 2 of the paper 
% A. Selivanov and E. Fridman, "An improved time-delay implementation of derivative-dependent feedback," 
% Automatica, vol. 98, pp. 269-276, 2018. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A, B, C   - parameters of (1) 
% Kbar      - cell array of gains from (25) 
% h         - sampling period 
% alpha     - decay rate

% Output: 
% flag =1 if feasible, =0 otherwise
%% Decision variables and notations
[n,m]=size(B); 
r=length(Kbar); 
P=sdpvar(n); 
W0=sdpvar(m); 

D=A+B*Kbar{1}*C; 
for i=1:r-1
    W{i}=sdpvar(m);	%#ok<AGROW>
    R{i}=sdpvar(m);	%#ok<AGROW>    
    D=D+B*Kbar{i+1}*C*A^i; 
end

H=exp(2*alpha*(r-1)*h)*W{r-1}+((r-1)/2)^2*R{r-1}; 
%% LMIs 
N=blkvar; 
N(1,1)=D'*P+P*D+2*alpha*P+h^2*(Kbar{1}*C*A)'*W0*(Kbar{1}*C*A); 
for i=1:r-2
    N(1,1)=N(1,1)+h^2*exp(2*alpha*i*h)*(Kbar{i+1}*C*A^(i+1))'*W{i}*(Kbar{i+1}*C*A^(i+1))...
        +(i*h/2)^2*(Kbar{i+1}*C*A^(i+1))'*R{i}*(Kbar{i+1}*C*A^(i+1)); 
end
N(1,2)=P*B; 
N(2,2)=-pi^2/4*exp(-2*alpha*h)*W0; 
N(1,2*r+1)=h*(Kbar{r}*C*A^(r-1)*D)'*H; 
N(2,2*r+1)=h*(Kbar{r}*C*A^(r-1)*B)'*H; 
for i=1:r-1
    N(1,2+i)=P*B; 
    N(1,r+1+i)=P*B; 
    N(2+i,2+i)=-pi^2/4*exp(-2*alpha*h)*W{i}; 
    N(r+1+i,r+1+i)=-exp(-2*alpha*i*h)*R{i}; 
    N(2+i,2*r+1)=h*(Kbar{r}*C*A^(r-1)*B)'*H; 
    N(r+1+i,2*r+1)=h*(Kbar{r}*C*A^(r-1)*B)'*H; 
end
N(2*r+1,2*r+1)=-H; 
N=sdpvar(N); 
%% Solution of LMIs
LMIs=[P>=0,N<=0]; 

options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

flag=0; 
if sol.problem==0
    [primal,~]=check(LMIs); 
    flag=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end
