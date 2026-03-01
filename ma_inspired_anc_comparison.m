clc;
clear all;
close all;

%% Simulation preparation
% Path settings
Pw = [0 0 0 0 0 0 0 0 0 0.8 0.6 -0.2 -0.5 -0.1 0.4 -0.05]; % primary path
Sw = [0 0 0 0 0 1 2.5 1.76 0.15 -0.4825 -0.18625 -0.005 -0.001875]; % secondary path
Tap_Sw = length(Sw);  

% Primary noise and desired signal
global DataLong; 
DataLong = 160000;      
snr = 30;
[X_noise, ~] = add_awgn(sine_wave_generator(200)+sine_wave_generator(400), snr); % primary noise
Yp = filter(Pw, 1, X_noise); % desired signal

% Key parameters
block = 50; % block length
dim = 20; % control filter length 
pop = 20; % population size
x_max = 0.1;
ub = x_max*ones(1, dim); % position upper bound
lb = -x_max*ones(1, dim); % position lower bound

maxIter1 = DataLong/(pop*block);
maxIter2 = DataLong/(pop*block*3); 
maxIter3 = DataLong/(pop*block*3/2);

q = 20; % number of simulation repetitions
ANR1 = zeros(q, DataLong);
ANR2 = zeros(q, DataLong);
ANR3 = zeros(q, DataLong);
ANR4 = zeros(q, DataLong);

%% Simulation start
for n = 1:q
    Me1 = zeros(1, DataLong);
    Me2 = zeros(1, DataLong);
    Me3 = zeros(1, DataLong);
    Me4 = zeros(1, DataLong);
    Md = zeros(1, DataLong);

    %% QPSO algorithm
    % Initialization
    alpha = 0.5; 
    ANCCx1 = zeros(1, dim);
    ActualSPx1 = zeros(1, Tap_Sw);  
    e_cont1 = zeros(1, DataLong);
    X = initialize_population(pop, ub, lb, dim);

    % Iteration start
    for t = 1:maxIter1
        fitness1 = zeros(1, pop);
        for p = 1:pop

            % calculate the index range of the current particle
            start_index = pop*(t-1)*block + (p-1)*block + 1; 
            end_index = pop*(t-1)*block + p*block;

            % traverse the data block of the current particle
            for j = start_index:end_index 
                ANCCx1 = [X_noise(j) ANCCx1(1:(dim-1))];
                ANCCy = X(p, :)*ANCCx1';   
                ActualSPx1 = [ANCCy ActualSPx1(:, (1:(Tap_Sw-1)))];
                e_cont1(j) = Yp(j) - ActualSPx1*Sw';
                fitness1(p) = fitness1(p) + e_cont1(j)^2/block;
            end
        end

        if t == 1

            % use the initial population as the historical optimum
            pBest = X; 
            pBestFitness = fitness1;

            % record the initial global optimum
            [~, index] = min(fitness1);
            gBestFitness = fitness1(index);
            gBest = X(index, :);
        else
            for i = 1:pop

                % update the historical optimum
                if fitness1(i) < pBestFitness(i)    
                    pBest(i, :) = X(i, :);
                    pBestFitness(i) = fitness1(i);
                end

                % update the global optimum
                if fitness1(i) < gBestFitness
                    gBestFitness = fitness1(i);
                    gBest = X(i, :);
                end
            end
        end

        C = sum(pBest, 1)/pop; 
        for i = 1:pop
            for j = 1:dim
                phi = rand;
                p = phi*pBest(i, j) + (1-phi)*gBest(j); 
                if rand < 0.5
                    X(i, j) = p + alpha*abs(X(i, j)-C(j))*log(1/rand);
                else
                    X(i, j) = p - alpha*abs(X(i, j)-C(j))*log(1/rand);
                end
            end
            X(i,:) = boundary_constraints(X(i, :), ub, lb, dim); % boundary constraint
        end
    end

    %% ABC algorithm
    % Initialization
    VarSize = [1 dim];
    empty_bee.Position = [];
    empty_bee.Cost = 0;
    npop = repmat(empty_bee, pop, 1);
    X = initialize_population(pop, ub, lb, dim);
    for i = 1:pop
        npop(i).Position = X(i, :); 
    end   
    a = 1; % acceleration coefficient upper bound                   
    L = round(0.6*dim*pop); % abandonment limit parameter
    ANCCx2 = zeros(1, dim);
    ActualSPx2 = zeros(1, Tap_Sw);  
    e_cont2 = zeros(1, DataLong);

    % Iteration start
    for t = 1:maxIter2
        for p = 1:pop
            start_index = 3*pop*(t-1)*block + (p-1)*block+1;
            end_index = 3*pop*(t-1)*block + p*block;
            for j = start_index:end_index
                ANCCx2 = [X_noise(j) ANCCx2(1:(dim-1))];
                ANCCy = npop(p).Position*ANCCx2';   
                ActualSPx2 = [ANCCy ActualSPx2(:, (1:(Tap_Sw-1)))]; 
                e_cont2(j) = Yp(j) - ActualSPx2*Sw';
                npop(p).Cost = npop(p).Cost + e_cont2(j)^2/block;
            end
        end

        BestSol.Cost = inf;
        for i = 1 : pop
            if npop(i).Cost <= BestSol.Cost
                BestSol = npop(i);
            end
        end

        C = zeros(pop, 1);
        % employed bee
        npop1 = repmat(empty_bee, pop, 1);
        for i = 1:pop                  
            K = [1:i-1 i+1:pop];
            k = K(randi([1 numel(K)])); % choose k randomly, not equal to i               
            phi = a*unifrnd(-1, +1, VarSize); % define acceleration coeff
            npop1(i).Position = npop(i).Position + phi.*(npop(i).Position-npop(k).Position); % new bee position             
            npop1(i).Position = boundary_constraints(npop1(i).Position, ub, lb, dim); % boundary constraint
        end

        for p = 1:pop
            start_index = pop*block + 3*pop*(t-1)*block + (p-1)*block + 1;
            end_index = pop*block + 3*pop*(t-1)*block + p*block;
            for j = start_index:end_index
                ANCCx2 = [X_noise(j) ANCCx2(1:(dim-1))];
                ANCCy = npop1(p).Position*ANCCx2';  
                ActualSPx2 = [ANCCy ActualSPx2(:, (1:(Tap_Sw-1)))];
                e_cont2(j) = Yp(j) - ActualSPx2*Sw';
                npop1(p).Cost = npop1(p).Cost + e_cont2(j)^2/block;
            end
        end

        for i = 1:pop
            if npop1(i).Cost <= npop(i).Cost
                npop(i) = npop1(i);
            else
                C(i) = C(i)+1;
            end
        end

        F = zeros(pop, 1);
        MeanCost = mean([npop.Cost]);
        for i = 1:pop
            F(i) = exp(-npop(i).Cost/MeanCost); % Convert Cost to Fitness
        end
        P = F/sum(F);

        % onlooker bee
        npop2 = repmat(empty_bee, pop, 1);
        for m = 1:pop
                
            % Select Source Site
            i = roulette_wheel_selection(P);
                
            % Choose k randomly, not equal to i
            K = [1:i-1 i+1:pop];
            k = K(randi([1 numel(K)]));
                
            % Define Acceleration Coeff.
            phi = a*unifrnd(-1, +1, VarSize);
                
            % New Bee Position
            npop2(m).Position = npop(i).Position + phi.*(npop(i).Position-npop(k).Position);
                
            % Apply Bounds
            npop2(m).Position = boundary_constraints(npop2(m).Position, ub, lb, dim); 
        end

        for p = 1:pop
            start_index = 2*pop*block + 3*pop*(t-1)*block + (p-1)*block + 1;
            end_index = 2*pop*block + 3*pop*(t-1)*block + p*block;       
            for j = start_index:end_index
                ANCCx2 = [X_noise(j) ANCCx2(1:(dim-1))];
                ANCCy = npop2(p).Position*ANCCx2';   
                ActualSPx2 = [ANCCy ActualSPx2(:, (1:(Tap_Sw-1)))]; % 0.66 * tanh(1.5 * ANCCy)
                e_cont2(j) = Yp(j) - ActualSPx2*Sw';
                npop2(p).Cost = npop2(p).Cost + e_cont2(j)^2/block;
            end
        end

        for i = 1:pop
            if npop2(i).Cost <= npop(i).Cost
                npop(i) = npop2(i);
            else
                C(i) = C(i) + 1;
            end
        end
    end

    %% RGA algorithm
    % Initialization       
    mut = 0.1; % mutation probability
    acr = 0.9; % crossover probability
    eta_c = 20;
    eta_m = 20;
    X = initialize_population(3*pop/2, ub, lb, dim); 
    ANCCx3 = zeros(1, dim);
    ActualSPx3 = zeros(1, Tap_Sw);  
    e_cont3 = zeros(1, DataLong);    

    % Iteration start
    for t = 1: maxIter3
        fitness3 = zeros(1, 3*pop/2); 
        for p = 1:pop
            start_index = 3/2*pop*(t-1)*block + (p-1)*block + 1;
            end_index = 3/2*pop*(t-1)*block + p* block;
            for j = start_index:end_index
                ANCCx3 = [X_noise(j) ANCCx3(1:(dim-1))];
                ANCCy = X(p, :)*ANCCx3';   
                ActualSPx3 = [ANCCy ActualSPx3(:, (1:(Tap_Sw-1)))]; % 0.66 * tanh(1.5 * ANCCy)
                e_cont3(j) = Yp(j) - ActualSPx3*Sw';
                fitness3(p) = fitness3(p) + e_cont3(j)^2/block;
            end
        end
    
        % Calculate the best fitness value
        chrome_pos = zeros(1, dim+1);
        chrome_pos(end) = inf;
        for i = 1:pop
            if fitness3(i) < chrome_pos(end) 
                chrome_pos(end) = fitness3(i); 
                chrome_pos(1:dim) = X(i, :); 
            end
        end
    
        % Selection
        for i = pop+1 : 3*pop/2
            idx1 = randi(pop);
            idx2 = randi(pop);
            if fitness3(idx1) < fitness3(idx2)
                X(i, :) = X(idx1, :);
            else
                X(i, :) = X(idx2, :);
            end
        end
    
        % Crossover
        parent1 = zeros(1, dim);
        parent2 = zeros(1, dim);
        child1 = zeros(1, dim);
        child2 = zeros(1, dim);
        for i = (pop+1):2:(3*pop/2)
            parent1 = X(i, :);
            parent2 = X(i+1, :);
            for j = 1:dim
                acr_rand = rand; 
                if acr_rand < acr 
                    r_h = rand;      
                    if r_h <= 0.5
                        alpha_h = (2*r_h)^(1/(eta_c+1));
                    else
                        alpha_h = (1/(2*(1-r_h)))^(1/(eta_c+1));
                    end
                    child1(j) = 0.5*((1-alpha_h)*parent1(j) + (1+alpha_h)*parent2(j));
                    child2(j) = 0.5*((1 + alpha_h)*parent1(j) + (1-alpha_h)*parent2(j));
                else
                    child1(j) = parent1(j);
                    child2(j) = parent2(j);
                end
            end  
            X(i, :) = child1;
            if i+1 <= 3*pop/2
                X(i+1, :) = child2;
            end
        end
        
        % Mutation
        for i = (pop+1):(3*pop/2)
            for j = 1:dim
                mut_rand = rand; 
                if mut_rand < mut  
                    r_h = rand;        
                    if r_h <= 0.5
                        delta_h = (2*r_h)^(1/(eta_m+1)) - 1;
                    else
                        delta_h = 1 - (2*(1-r_h))^(1/(eta_m+1));
                    end
                    X(i, j) = X(i, j) + delta_h;
                end
            end
            X(i, :) = boundary_constraints(X(i, :), ub, lb, dim);
        end
        
        for p = pop:(3*pop/2-1)
            start_index = 3/2*pop*(t-1)*block + p*block + 1;
            end_index = 3/2*pop*(t-1)*block + (p+1)*block;
            for j = start_index:end_index
                ANCCx3 = [X_noise(j) ANCCx3(1:(dim-1))];
                ANCCy = X(p+1, :)*ANCCx3';   
                ActualSPx3 = [ANCCy ActualSPx3(:, (1:(Tap_Sw-1)))]; 
                e_cont3(j) = Yp(j) - ActualSPx3*Sw';
                fitness3(p+1) = fitness3(p+1) + e_cont3(j)^2/block;
            end
        end
        
        [sorted_fitness, idx] = sort(fitness3);    
        top_20_vectors = X(idx(1:20), :);
        for i = 1:pop
            X(i, :) = top_20_vectors(i, :);
        end
    end

    %% GWO algorithm
    % Initialization 
    eta = 0.2;
    b = 0.1;
    X = initialize_population(pop, ub, lb, dim);    
    ANCCx4 = zeros(1, dim);
    ActualSPx4 = zeros(1, Tap_Sw);  
    e_cont4 = zeros(1, DataLong);
    D = zeros(1, maxIter1);
    a1 = zeros(1, maxIter1);

    % Iteration start
    for t = 1:maxIter1
        fitness4 = zeros(1, pop); 
        for p = 1:pop
            start_index = pop*(t-1)*block + (p-1)*block + 1;
            end_index = pop*(t-1)*block + p*block;
            for j = start_index:end_index
                ANCCx4 = [X_noise(j) ANCCx4(1:(dim-1))];
                ANCCy = X(p, :)*ANCCx4';
                ActualSPx4 = [ANCCy ActualSPx4(:, (1:(Tap_Sw-1)))]; 
                e_cont4(j) = Yp(j) - ActualSPx4*Sw';
                fitness4(p) = fitness4(p) + e_cont4(j)^2/block;
            end
        end

        % Initialize alpha, beta, and delta
        Alpha_pos = zeros(1, dim);
        Beta_pos = zeros(1, dim);
        Delta_pos = zeros(1, dim);
        Alpha_score = inf; 
        Beta_score = inf; 
        Delta_score = inf; 

        % Update alpha, beta, and delta
        for i = 1:pop

            % Update alpha
            if fitness4(i) < Alpha_score 
                Alpha_score = fitness4(i); 
                Alpha_pos = X(i, :);
            end

            % Update beta
            if (fitness4(i) > Alpha_score) && (fitness4(i) < Beta_score) 
                Beta_score = fitness4(i); 
                Beta_pos = X(i, :);
            end

            % Update delta
            if (fitness4(i) > Alpha_score) && (fitness4(i) > Beta_score) && (fitness4(i) < Delta_score)
                Delta_score = fitness4(i); 
                Delta_pos = X(i, :);
            end
        end

        X_avg = mean(X, 1); 
        distances = zeros(pop, 1); 
        for i = 1:pop
            distances(i) = norm(X(i, :) - X_avg); 
        end
        D(t) = mean(distances);
        
        if D(t) > eta 
            a1(t) = 2*(1-t/maxIter1);
            T = t;
            a_mid = a1(t);
        else
            a1(t) = a_mid*exp(-b*(t-T));
        end
        
        for i = 1:pop
            if isequal(X(i,:), Alpha_pos) || isequal(X(i,:), Beta_pos) || isequal(X(i,:), Delta_pos)
                continue
            end

            for j = 1:dim                
                r1 = rand(); % r1 is a random number in [0,1]
                r2 = rand(); % r2 is a random number in [0,1]                        
                A1 = 2*a1(t)*r1 - a1(t); 
                C1 = 2*r2;                       
                D_alpha = abs(C1*Alpha_pos(j)-X(i,j)); 
                X1 = Alpha_pos(j) - A1*D_alpha; 
                           
                r1 = rand();
                r2 = rand();
                A2 = 2*a1(t)*r1 - a1(t); 
                C2 = 2*r2;
                D_beta = abs(C2*Beta_pos(j)-X(i,j)); 
                X2 = Beta_pos(j) - A2*D_beta;       
                            
                r1 = rand();
                r2 = rand();           
                A3 = 2*a1(t)*r1 - a1(t); 
                C3 = 2*r2;  
                D_delta = abs(C3*Delta_pos(j)-X(i,j)); 
                X3 = Delta_pos(j) - A3*D_delta;          

                wa = 1/3;
                wb = 1/3;
                wc = 1/3;
                X(i, j) = wa*X1 + wb*X2 + wc*X3;
            end

            X(i,:) = boundary_constraints(X(i, :), ub, lb, dim);
        end
    end

    %% Performance evaluation
    lambda = 0.999;
    for i = 1:DataLong 
        if i == 1
            Me1(i) = (1-lambda)*abs(e_cont1(i));
            Me2(i) = (1-lambda)*abs(e_cont2(i));
            Me3(i) = (1-lambda)*abs(e_cont3(i));
            Me4(i) = (1-lambda)*abs(e_cont4(i));
            Md(i) = (1-lambda)*abs(Yp(i));
        else
            Me1(i) = lambda*Me1(i-1) + (1-lambda)*abs(e_cont1(i));
            Me2(i) = lambda*Me2(i-1) + (1-lambda)*abs(e_cont2(i));
            Me3(i) = lambda*Me3(i-1) + (1-lambda)*abs(e_cont3(i));
            Me4(i) = lambda*Me4(i-1) + (1-lambda)*abs(e_cont4(i));
            Md(i) = lambda*Md(i-1) + (1-lambda)*abs(Yp(i));
        end

        ANR1(n, i) = 20*log10(Md(i)/Me1(i));
        ANR2(n, i) = 20*log10(Md(i)/Me2(i));
        ANR3(n, i) = 20*log10(Md(i)/Me3(i));
        ANR4(n, i) = 20*log10(Md(i)/Me4(i));
    end
end

%% Data processing
% Average
AANR1 = mean(ANR1, 1);
AANR2 = mean(ANR2, 1);
AANR3 = mean(ANR3, 1);
AANR4 = mean(ANR4, 1);

% Gaussian filtering
windowSize = 5000; 
weights = gausswin(windowSize);
AANR11 = conv(AANR1, weights, 'same')/sum(weights);
AANR22 = conv(AANR2, weights, 'same')/sum(weights);
AANR33 = conv(AANR3, weights, 'same')/sum(weights);
AANR44 = conv(AANR4, weights, 'same')/sum(weights);

%% Plot
plot(AANR11, 'LineWidth', 1.5);
hold on;
plot(AANR22, 'LineWidth', 1.5);
hold on;
plot(AANR33, 'LineWidth', 1.5);
hold on;
plot(AANR44, 'LineWidth', 1.5);

ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xtickformat('%.0f');    
ytickformat('%.0f');
xlim([windowSize/2, DataLong-windowSize]);
ylim([-15, 35]);        
xlabel('{\fontname{Times New Roman}Sample Number}');
ylabel('{\fontname{Times New Roman}ANR/dB}');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
legend('{\fontname{Times New Roman}QPSO}','{\fontname{Times New Roman}ABC}','{\fontname{Times New Roman}RGA}','{\fontname{Times New Roman}IGWO}',...%    
     'NumColumns',1,'FontSize', 12, 'Location', 'southeast');
grid on
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 1);