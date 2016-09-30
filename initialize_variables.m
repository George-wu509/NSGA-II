function f = initialize_variables(N, M, V, SUL, SDL, I, D,s0,qu,dayflow_r)

%% function f = initialize_variables(N, M, V, SUL, SDL) 

% N - Population size
% M - Number of objective functions
% V - Number of decision variables
% SUL - M5 Upper Line
% SDL - M5 Down Line
% K is the total number of array elements. For ease of computation decision
% variables and objective functions are concatenated to form a single
% array. For crossover and mutation only the decision variables are used
% while for selection, only the objective variable are utilized.

K = V+M;
f=zeros(N,V+M);
%% Initialize each chromosome
% For each chromosome perform the following (N is the population size)

%% initial generation ( 1 )
for i = 1 : 0.3*N                  % ��l����1: �J�y�qI ����
    for j=1:V
        f(i,j)=I(j)-10+20*rand(1);
    end
    f(i,V + 1: K) = evaluate_objective(f(i,1:V),SUL, SDL, I, D,s0,qu,dayflow_r);
end

for i = 0.3*N+1 : 0.6*N           % ��l����2: ����Ъ�D ����
    for j=1:V
        f(i,j)=D(j)-10+20*rand(1);
    end
    f(i,V + 1: K) = evaluate_objective(f(i,1:V),SUL, SDL, I, D,s0,qu,dayflow_r);
end

for i=0.6*N+1:N                   % ��l����3: ����q0~60 ��
    for j=1:V
        f(i,j)=rand(1)*60;
    end
    f(i,V + 1: K) = evaluate_objective(f(i,1:V),SUL, SDL, I, D,s0,qu,dayflow_r);
end
