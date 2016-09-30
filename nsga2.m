function [solution,his,h10,GSI]=nsga2(I,s0,qu,dayflow_r)
%   [solution,his,h10,GSI]=nsga2(I,s0,qu,dayflow_r);
%% function nsga_2(pop,gen)

%% Set pop and gen
pop = 1000; % pop - Population size
gen = 400;  % gen - Total number of generations
M=2;  % M is the dimension of the objective space 

%% Set variables dimension and min-max range
V=36; % V is the dimension of decision variable space
% min_range and max_range are the range for the variables in the decision variable space. 
min_range=zeros(1,36);
max_range=ones(1,36)*254;

% SUL - M5 Upper Line, SDL - M5 Down Line
SUL = [254.00 254.00 254.00 248.00 242.09 237.11 228.94 220.96 213.16 207.05 201.06 194.45 188.70 183.05 176.84 176.84 176.84 176.84 178.89 180.97 183.75 193.00 202.55 213.16 213.16 213.16 213.16 226.53 240.42 254.00 254.00 254.00 254.00 254.00 254.00 254.00];
SDL = [213.16 213.16 213.16 201.05 188.94 176.84 166.15 155.46 144.77 127.32 109.87 92.41 85.48 78.54 71.61 71.61 71.61 71.61 78.54 85.48 92.41 100.52 108.62 116.72 116.72 116.72 116.72 126.07 135.42 144.77 155.46 166.15 176.84 188.94 201.05 213.16];
D=[15.94 15.94 17.53 27.41 27.96 22.70 29.76 30.34 34.60 29.19 28.32 28.50 26.67 25.81 28.35 27.49 27.92 28.79 31.09 31.27 35.03 31.05 30.95 33.80 28.10 27.25 28.11 29.13 29.13 31.09 28.83 25.81 23.63 15.95 15.95 17.54];

%% Initialize the population
% Population is initialized with random values which are within the
% specified range. Each chromosome consists of the decision variables. Also
% the value of the objective functions, rank and crowding distance
% information is also added to the chromosome vector but only the elements
% of the vector which has the decision variables are operated upon to
% perform the genetic operations like corssover and mutation.
chromosome = initialize_variables(pop, M, V, SUL, SDL, I, D, s0,qu,dayflow_r);

%% Sort the initialized population
% Sort the population using non-domination-sort. This returns two columns
% for each individual which are the rank and the crowding distance
% corresponding to their position in the front they belong. At this stage
% the rank and the crowding distance for each chromosome is added to the
% chromosome vector for easy of computation.
chromosome = non_domination_sort_mod(chromosome, M, V);

%% Start the evolution process
% The following are performed in each generation
% * Select the parents which are fit for reproduction
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.
his=zeros(gen,2*M);
h10=zeros(pop,gen/10*M);
for i = 1 : gen
    tic;
    % Select the parents
    % Parents are selected for reproduction to generate offspring. The
    % original NSGA-II uses a binary tournament selection based on the
    % crowded-comparision operator. The arguments are 
    % pool - size of the mating pool. It is common to have this to be half the
    %        population size.
    % tour - Tournament size. Original NSGA-II uses a binary tournament
    %        selection, but to see the effect of tournament size this is kept
    %        arbitary, to be choosen by the user.
    pool = round(pop/2);
    tour = 2;
    % Selection process
    % A binary tournament selection is employed in NSGA-II. In a binary
    % tournament selection process two individuals are selected at random
    % and their fitness is compared. The individual with better fitness is
    % selcted as a parent. Tournament selection is carried out until the
    % pool size is filled. Basically a pool size is the number of parents
    % to be selected. The input arguments to the function
    % tournament_selection are chromosome, pool, tour. The function uses
    % only the information from last two elements in the chromosome vector.
    % The last element has the crowding distance information while the
    % penultimate element has the rank information. Selection is based on
    % rank and if individuals with same rank are encountered, crowding
    % distance is compared. A lower rank and higher crowding distance is
    % the selection criteria.
    parent_chromosome = tournament_selection(chromosome, pool, tour);

    % Perfrom crossover and Mutation operator
    % The original NSGA-II algorithm uses Simulated Binary Crossover (SBX) and
    % Polynomial  mutation. Crossover probability pc = 0.9 and mutation
    % probability is pm = 1/n, where n is the number of decision variables.
    % Both real-coded GA and binary-coded GA are implemented in the original
    % algorithm, while in this program only the real-coded GA is considered.
    % The distribution indeices for crossover and mutation operators as mu = 20
    % and mum = 20 respectively.
    mu = 20;
    mum = 20;
    offspring_chromosome = ...
        genetic_operator(parent_chromosome, ...
        M, V, mu, mum, min_range, max_range, SUL, SDL, I, D,s0,qu,dayflow_r);

    % Intermediate population
    % Intermediate population is the combined population of parents and
    % offsprings of the current generation. The population size is two
    % times the initial population.
    
    [main_pop,temp] = size(chromosome);
    [offspring_pop,temp] = size(offspring_chromosome);
    % temp is a dummy variable.
    clear temp
    % intermediate_chromosome is a concatenation of current population and
    % the offspring population.
    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = ...
        offspring_chromosome;

    % Non-domination-sort of intermediate population
    % The intermediate population is sorted again based on non-domination sort
    % before the replacement operator is performed on the intermediate
    % population.
    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, M, V);
    % Perform Selection
    % Once the intermediate population is sorted only the best solution is
    % selected based on it rank and crowding distance. Each front is filled in
    % ascending order until the addition of population size is reached. The
    % last front is included in the population based on the individuals with
    % least crowding distance
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
    fprintf('%d/%d generation completed\n',i,gen);
    toc
    if ~mod(i,100)
        clc
        fprintf('%d generations completed\n',i);
    end
    
    if rem(i,10)==0
        [his(i,:),h10(:,(i/10-1)*M+1:i/10*M)]=average_fitness(chromosome,V,M);
    else
        his(i,:)=average_fitness(chromosome,V,M);
    end
end % end for

%% Result
% Save the result in ASCII text format.                          % S(t)為每旬水庫蓄水量(百萬力方公尺)
S=zeros(pop,36);
for i=1:pop
    S(i,1)=s0+I(1)-chromosome(i,1);
        for t=2:36
            S(i,t)=S(i,t-1)+I(t)-chromosome(i,t);
        end
end
chromosome(:,1:36)=S;
solution=chromosome;
save solution.txt chromosome -ASCII
save history.txt his -ASCII

%% Visualize
% The following is used to visualize the result if objective space
% dimension is visualizable.
if M == 2
    plot(-chromosome(:,V + 1),-chromosome(:,V + 2),'*');
elseif M ==3
    plot3(-chromosome(:,V + 1),chromosome(:,V + 2),-chromosome(:,V + 3),'*');
    grid on;
end
GSI=final_objective(chromosome,pop, I, D,s0);
%save GSI.txt GSI -ASCII
    
