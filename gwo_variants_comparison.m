clc;
clear all;
close all;

%% Simulation preparation
% Path settings
Pw = [0 0 0 0 0 0 0 0 0 0.8 0.6 -0.2 -0.5 -0.1 0.4 -0.05];    
Sw = [0 0 0 0 0 1 2.5 1.76 0.15 -0.4825 -0.18625 -0.005 -0.001875]; 
Tap_Sw = length(Sw);

% Key parameters
block = 50;
dim = 20; 
pop = 20;
x_max = 0.1;
ub = x_max * ones(1, dim);
lb = -x_max * ones(1, dim);

% Primary noise and desired signal
global DataLong; 
DataLong = 60000;       
snr = 30;    
[X_noise, ~] = add_awgn(sine_wave_generator(200)+sine_wave_generator(400), snr);
Yp = filter(Pw, 1, X_noise);  

q = 50;
maxIter = DataLong/(pop*block);

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

    %% Standard GWO 
    % Initialization
    ANCCx1 = zeros(1, dim);
    ActualSPx1 = zeros(1, Tap_Sw);  
    e_cont1 = zeros(1, DataLong);
    X = initialize_population(pop, ub, lb, dim); 

    %% Iteration start
    for t = 1: maxIter
        fitness1 = zeros(1,pop); 
        for p = 1:pop
            start_index = pop*(t-1)*block + (p-1)*block + 1;
            end_index = pop*(t-1)*block + p*block;
        
            for j = start_index:end_index
                ANCCx1 = [X_noise(j) ANCCx1(1:(dim-1))];
                ANCCy = X(p, :)*ANCCx1';  
                ActualSPx1 = [ANCCy ActualSPx1(:, (1:(Tap_Sw-1)))];
                e_cont1(j) = Yp(j) - ActualSPx1*Sw';
                fitness1(p) = fitness1(p) + e_cont1(j)^2/block;
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
            if fitness1(i) < Alpha_score 
                Alpha_score = fitness1(i); 
                Alpha_pos = X(i, :);
            end
    
            % Update beta
            if (fitness1(i) > Alpha_score) && (fitness1(i) < Beta_score) 
                Beta_score = fitness1(i); 
                Beta_pos = X(i, :);
            end
    
            % Update delta
            if (fitness1(i) > Alpha_score) && (fitness1(i) > Beta_score) && (fitness1(i) < Delta_score) 
                Delta_score = fitness1(i); 
                Delta_pos = X(i, :);
            end
        end
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop
            if isequal(X(i, :), Alpha_pos) || isequal(X(i, :), Beta_pos) || isequal(X(i, :), Delta_pos)
                continue
            end
    
            for j = 1:dim                           
                r1 = rand(); % r1 is a random number in [0,1]
                r2 = rand(); % r2 is a random number in [0,1]                        
                A1 = 2*a*r1 - a; 
                C1 = 2*r2;                         
                D_alpha = abs(C1*Alpha_pos(j)-X(i, j)); 
                X1 = Alpha_pos(j) - A1*D_alpha;
                          
                r1 = rand();
                r2 = rand();         
                A2 = 2*a*r1 - a; 
                C2 = 2*r2;             
                D_beta = abs(C2*Beta_pos(j)-X(i, j)); 
                X2 = Beta_pos(j) - A2*D_beta;       
               
                r1 = rand();
                r2 = rand();           
                A3 = 2*a*r1 - a; 
                C3 = 2*r2;           
                D_delta = abs(C3*Delta_pos(j)-X(i, j)); 
                X3 = Delta_pos(j) - A3*D_delta;              
    
                wa = 1/3;
                wb = 1/3;
                wc = 1/3;
                X(i, j) = wa*X1 + wb*X2 + wc*X3;
            end
    
            X(i, :) = boundary_constraints(X(i, :), ub, lb, dim);
        end
    end

    %% Improved GWO (sigma = 0.1)
    % Initialization    
    ANCCx2 = zeros(1, dim);
    ActualSPx2 = zeros(1, Tap_Sw);  
    e_cont2 = zeros(1, DataLong);
    D1 = zeros(1, maxIter);
    a1 = zeros(1, maxIter);
    X = initialize_population(pop, ub, lb, dim);

    % Iteration start
    for t = 1:maxIter
        fitness2 = zeros(1, pop);
        for p = 1:pop
        
            % calculate the index range of the current particle
            start_index = pop*(t-1)*block + (p-1)*block + 1; 
            end_index = pop*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
            for j = start_index:end_index 
                ANCCx2 = [X_noise(j) ANCCx2(1:(dim-1))];
                ANCCy = X(p, :)*ANCCx2';   
                ActualSPx2 = [ANCCy ActualSPx2(:, (1:(Tap_Sw-1)))];
                e_cont2(j) = Yp(j) - ActualSPx2*Sw';
                fitness2(p) = fitness2(p) + e_cont2(j)^2/block;
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
            if fitness2(i) < Alpha_score 
                Alpha_score = fitness2(i); 
                Alpha_pos = X(i, :);
            end
    
            % Update beta
            if (fitness2(i) > Alpha_score) && (fitness2(i) < Beta_score) 
                Beta_score = fitness2(i); 
                Beta_pos = X(i, :);
            end
    
            % Update delta
            if (fitness2(i) > Alpha_score) && (fitness2(i) > Beta_score) && (fitness2(i) < Delta_score) 
                Delta_score = fitness2(i); 
                Delta_pos = X(i, :);
            end
        end

        X_avg = mean(X, 1);
        distances = zeros(pop, 1);
        for i = 1:pop
            distances(i) = norm(X(i, :) - X_avg);
        end
        D1(t) = mean(distances);
        
        eta = 0.1;
        b = 0.1;
        if D1(t) > eta
            a1(t) = 2*(1 - t/maxIter);
            T = t;
            a_mid = a1(t);
        else
            a1(t) = a_mid*exp(-b*(t-T));
        end
        
        for i = 1:pop
            if isequal(X(i, :), Alpha_pos) || isequal(X(i, :), Beta_pos) || isequal(X(i, :), Delta_pos)
                continue
            end
    
            for j = 1:dim                            
                r1 = rand(); % r1 is a random number in [0,1]
                r2 = rand(); % r2 is a random number in [0,1]                        
                A1 = 2*a1(t)*r1 - a1(t); 
                C1 = 2*r2;                         
                D_alpha = abs(C1*Alpha_pos(j)-X(i, j)); 
                X1 = Alpha_pos(j) - A1*D_alpha;
                          
                r1 = rand();
                r2 = rand();         
                A2 = 2*a1(t)*r1 - a1(t); 
                C2 = 2*r2;             
                D_beta = abs(C2*Beta_pos(j)-X(i, j)); 
                X2 = Beta_pos(j) - A2*D_beta;       
               
                r1 = rand();
                r2 = rand();           
                A3 = 2*a1(t)*r1 - a1(t); 
                C3 = 2*r2;           
                D_delta = abs(C3*Delta_pos(j)-X(i, j)); 
                X3 = Delta_pos(j) - A3*D_delta;              
    
                wa = 1/3;
                wb = 1/3;
                wc = 1/3;
                X(i, j) = wa*X1 + wb*X2 + wc*X3;
            end
    
            X(i, :) = boundary_constraints(X(i, :), ub, lb, dim);
        end      
    end

    %% Improved GWO (sigma = 0.5)
    % Initialization    
    ANCCx3 = zeros(1, dim);
    ActualSPx3 = zeros(1, Tap_Sw);  
    e_cont3 = zeros(1, DataLong);
    D2 = zeros(1, maxIter);
    a2 = zeros(1, maxIter);
    X = initialize_population(pop, ub, lb, dim);
    
    % Iteration start
    for t = 1:maxIter
        fitness3 = zeros(1, pop);
        for p = 1:pop
        
            % calculate the index range of the current particle
            start_index = pop*(t-1)*block + (p-1)*block + 1; 
            end_index = pop*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
            for j = start_index:end_index 
                ANCCx3 = [X_noise(j) ANCCx3(1:(dim-1))];
                ANCCy = X(p, :)*ANCCx3';   
                ActualSPx3 = [ANCCy ActualSPx3(:, (1:(Tap_Sw-1)))];
                e_cont3(j) = Yp(j) - ActualSPx3*Sw';
                fitness3(p) = fitness3(p) + e_cont3(j)^2/block;
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
            if fitness3(i) < Alpha_score 
                Alpha_score = fitness3(i); 
                Alpha_pos = X(i, :);
            end
    
            % Update beta
            if (fitness3(i) > Alpha_score) && (fitness3(i) < Beta_score) 
                Beta_score = fitness3(i); 
                Beta_pos = X(i, :);
            end
    
            % Update delta
            if (fitness3(i) > Alpha_score) && (fitness3(i) > Beta_score) && (fitness3(i) < Delta_score) 
                Delta_score = fitness3(i); 
                Delta_pos = X(i, :);
            end
        end
    
        X_avg = mean(X, 1);
        distances = zeros(pop, 1);
        for i = 1:pop
            distances(i) = norm(X(i, :) - X_avg);
        end
        D2(t) = mean(distances);
        
        b = 0.5;
        if D2(t) > eta
            a2(t) = 2*(1 - t/maxIter);
            T = t;
            a_mid = a2(t);
        else
            a2(t) = a_mid*exp(-b*(t-T));
        end
        
        for i = 1:pop
            if isequal(X(i, :), Alpha_pos) || isequal(X(i, :), Beta_pos) || isequal(X(i, :), Delta_pos)
                continue
            end
    
            for j = 1:dim                            
                r1 = rand(); % r1 is a random number in [0,1]
                r2 = rand(); % r2 is a random number in [0,1]                        
                A1 = 2*a2(t)*r1 - a2(t); 
                C1 = 2*r2;                         
                D_alpha = abs(C1*Alpha_pos(j)-X(i, j)); 
                X1 = Alpha_pos(j) - A1*D_alpha;
                          
                r1 = rand();
                r2 = rand();         
                A2 = 2*a2(t)*r1 - a2(t); 
                C2 = 2*r2;             
                D_beta = abs(C2*Beta_pos(j)-X(i, j)); 
                X2 = Beta_pos(j) - A2*D_beta;       
               
                r1 = rand();
                r2 = rand();           
                A3 = 2*a2(t)*r1 - a2(t); 
                C3 = 2*r2;           
                D_delta = abs(C3*Delta_pos(j)-X(i, j)); 
                X3 = Delta_pos(j) - A3*D_delta;              
    
                wa = 1/3;
                wb = 1/3;
                wc = 1/3;
                X(i, j) = wa*X1 + wb*X2 + wc*X3;
            end
    
            X(i, :) = boundary_constraints(X(i, :), ub, lb, dim);
        end
    end

    %% Improved GWO (sigma = 1)
    ANCCx4 = zeros(1, dim);
    ActualSPx4 = zeros(1, Tap_Sw);  
    e_cont4 = zeros(1, DataLong);
    D3 = zeros(1, maxIter);
    a3 = zeros(1, maxIter);
    X = initialize_population(pop, ub, lb, dim);
    
    % Iteration start
    for t = 1: maxIter
        fitness4 = zeros(1, pop);
        for p = 1:pop
        
            % calculate the index range of the current particle
            start_index = pop*(t-1)*block + (p-1)*block + 1; 
            end_index = pop*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
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
        D3(t) = mean(distances);
        
        b = 1;
        if D3(t) > eta
            a3(t) = 2*(1 - t/maxIter);
            T = t;
            a_mid = a3(t);
        else
            a3(t) = a_mid*exp(-b*(t-T));
        end
        
        for i = 1:pop
            if isequal(X(i, :), Alpha_pos) || isequal(X(i, :), Beta_pos) || isequal(X(i, :), Delta_pos)
                continue
            end
    
            for j = 1:dim                            
                r1 = rand(); % r1 is a random number in [0,1]
                r2 = rand(); % r2 is a random number in [0,1]                        
                A1 = 2*a3(t)*r1 - a3(t); 
                C1 = 2*r2;                         
                D_alpha = abs(C1*Alpha_pos(j)-X(i, j)); 
                X1 = Alpha_pos(j) - A1*D_alpha;
                          
                r1 = rand();
                r2 = rand();         
                A2 = 2*a3(t)*r1 - a3(t); 
                C2 = 2*r2;             
                D_beta = abs(C2*Beta_pos(j)-X(i, j)); 
                X2 = Beta_pos(j) - A2*D_beta;       
               
                r1 = rand();
                r2 = rand();           
                A3 = 2*a3(t)*r1 - a3(t); 
                C3 = 2*r2;           
                D_delta = abs(C3*Delta_pos(j)-X(i, j)); 
                X3 = Delta_pos(j) - A3*D_delta;              
    
                wa = 1/3;
                wb = 1/3;
                wc = 1/3;
                X(i, j) = wa*X1 + wb*X2 + wc*X3;
            end
    
            X(i, :) = boundary_constraints(X(i, :), ub, lb, dim);
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
legend('{\fontname{Times New Roman}Standard GWO}', '{\fontname{Times New Roman}Improved GWO (\sigma=0.1)}', '{\fontname{Times New Roman}Improved GWO (\sigma=0.5)}', '{\fontname{Times New Roman}Improved GWO (\sigma=1)}',...%    
     'NumColumns', 1, 'FontSize', 12, 'Location', 'northeast');
grid on
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 1);