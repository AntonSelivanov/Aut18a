% This MATLAB program checks the feasibility of LMIs from Theorems 1 and 2 of the paper 
% A. Selivanov and E. Fridman, "An improved time-delay implementation of derivative-dependent feedback," Automatica, 2018. 

%% Example 1 
disp('Chain of three integrators') 
A=[0 1 0; 0 0 1; 0 0 0];    % parameters from (21)
B=[0; 0; 1];                % 
C=[1 0 0];                  % 
Kbar={-2e-4, -.06, -.342};  % parameters from (22)

% Continuous-time control 
h=2.529; % delay 
alpha=0; % decay rate
if LMI_Aut18a_th1(A,B,C,Kbar,h,alpha)
    disp(['Theorem 1: Feasible for h=' num2str(h)]); 
else
    disp(['Theorem 1: Not feasible for h=' num2str(h)]); 
end

% Sampled-data control 
h=1.436;	% sampling period 
alpha=1e-3;	% decay rate
if LMI_Aut18a_th2(A,B,C,Kbar,h,alpha)
    disp(['Theorem 2: Feasible for h=' num2str(h)]); 
else
    disp(['Theorem 2: Not feasible for h=' num2str(h)]); 
end
%% Example 2 
disp('Chain of four integrators') 
A=[0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 0]; % parameters from (23)
B=[0; 0; 0; 1];                         % 
C=[1 0 0 0];                            % 
Kbar={-.0208,-.32,-1.18,-.7};           % parameters from (24)

% Continuous-time control 
h=.169;  % delay 
alpha=0; % decay rate
if LMI_Aut18a_th1(A,B,C,Kbar,h,alpha)
    disp(['Theorem 1: Feasible for h=' num2str(h)]); 
else
    disp(['Theorem 1: Not feasible for h=' num2str(h)]); 
end

% Sampled-data control 
h=.1;       % sampling period 
alpha=.01;	% decay rate
if LMI_Aut18a_th2(A,B,C,Kbar,h,alpha)
    disp(['Theorem 2: Feasible for h=' num2str(h)]); 
else
    disp(['Theorem 2: Not feasible for h=' num2str(h)]); 
end
%% Example 3 
disp('Furuta pendulum') 
A=[0 1 0 0; 37.377 -.515 0 .142; 0 0 0 1; -8.228 .113 0 -.173]; % parameters from (30)
B=[0; -35.42; 0; 43.28];                                        % 
C=[1 0 0 0; 0 0 1 0];                                           % 
Kbar={[1.2826 1.3e-3],[.1209 8.6e-3]};                          % controller gains

% Sampled-data control 
h=.104;  % sampling period 
alpha=0; % decay rate
if LMI_Aut18a_th2(A,B,C,Kbar,h,alpha)
    disp(['Theorem 2: Feasible for h=' num2str(h)]); 
else
    disp(['Theorem 2: Not feasible for h=' num2str(h)]); 
end