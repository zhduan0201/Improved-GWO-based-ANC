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

pop1 = 10;  
pop2 = 20;
pop3 = 30;
pop4 = 40;
pop5 = 50;

x_max = 0.1;
ub = x_max * ones(1, dim);
lb = -x_max * ones(1, dim);

% Primary noise and desired signal
maxIter = 30;
global DataLong;
DataLong = pop5*block*maxIter;   
snr = 30;    
[X_noise, ~] = add_awgn(sine_wave_generator(200)+sine_wave_generator(400), snr);
Yp=filter(Pw, 1, X_noise); 

DataLong1 = pop1*block*maxIter;
DataLong2 = pop2*block*maxIter;
DataLong3 = pop3*block*maxIter; 
DataLong4 = pop4*block*maxIter; 
DataLong5 = pop5*block*maxIter; 

IterCurve1 = zeros(1, maxIter); 
IterCurve2 = zeros(1, maxIter);
IterCurve3 = zeros(1, maxIter);
IterCurve4 = zeros(1, maxIter);
IterCurve5 = zeros(1, maxIter);

%% Simulation start
q = 50; 
for n = 1:q
    %% GWO algorithm (population size = 10)
    % Initialization    
    ANCCx1 = zeros(1, dim);
    ActualSPx1 = zeros(1, Tap_Sw);  
    e_cont1 = zeros(1, DataLong1);
    X = initialize_population(pop1, ub, lb, dim);

    % Iteration start
    for t = 1:maxIter
        fitness1 = zeros(1, pop1);
        for p = 1:pop1
        
            % calculate the index range of the current particle
            start_index = pop1*(t-1)*block + (p-1)*block + 1; 
            end_index = pop1*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
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
        for i = 1:pop1
    
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
    
        IterCurve1(t) = IterCurve1(t) + Alpha_score;
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop1
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

    %% GWO algorithm (population size = 20)
    % Initialization    
    ANCCx2 = zeros(1, dim);
    ActualSPx2 = zeros(1, Tap_Sw);  
    e_cont2 = zeros(1, DataLong2);
    X = initialize_population(pop2, ub, lb, dim);

    % Iteration start
    for t = 1:maxIter
        fitness2  =zeros(1, pop2); 
        for p = 1:pop2
        
            % calculate the index range of the current particle
            start_index = pop2*(t-1)*block + (p-1)*block + 1; 
            end_index = pop2*(t-1)*block + p*block;
        
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
        for i = 1:pop2
    
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
    
        IterCurve2(t) = IterCurve2(t) + Alpha_score;
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop2
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

    %% GWO algorithm (population size = 30)
    % Initialization    
    ANCCx3 = zeros(1, dim);
    ActualSPx3 = zeros(1, Tap_Sw);  
    e_cont3 = zeros(1, DataLong3);
    X = initialize_population(pop3, ub, lb, dim);

    % Iteration start
    for t = 1:maxIter
        fitness3 = zeros(1, pop3);
        for p = 1:pop3
        
            % calculate the index range of the current particle
            start_index = pop3*(t-1)*block + (p-1)*block + 1; 
            end_index = pop3*(t-1)*block + p*block;
        
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
        for i = 1:pop3
    
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
    
        IterCurve3(t) = IterCurve3(t) + Alpha_score;
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop3
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

    %% GWO algorithm (population size = 40)
    % Initialization   
    ANCCx4 = zeros(1, dim);
    ActualSPx4 = zeros(1, Tap_Sw);  
    e_cont4 = zeros(1, DataLong4);
    X = initialize_population(pop4, ub, lb, dim);

    % Iteration start
    for t = 1: maxIter
        fitness4 = zeros(1, pop4);
        for p = 1:pop4
        
            % calculate the index range of the current particle
            start_index = pop4*(t-1)*block + (p-1)*block + 1; 
            end_index = pop4*(t-1)*block + p*block;
        
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
        for i = 1:pop4
    
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
    
        IterCurve4(t) = IterCurve4(t) + Alpha_score;
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop4
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

    %% GWO algorithm (population size = 50)
    % Initialization    
    ANCCx5 = zeros(1, dim);
    ActualSPx5 = zeros(1, Tap_Sw);  
    e_cont5 = zeros(1, DataLong5);
    X = initialize_population(pop5, ub, lb, dim);

    % Iteration start
    for t = 1:maxIter
        fitness5 = zeros(1, pop5);
        for p = 1:pop5
        
            % calculate the index range of the current particle
            start_index = pop5*(t-1)*block + (p-1)*block + 1; 
            end_index = pop5*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
            for j = start_index:end_index 
                ANCCx5 = [X_noise(j) ANCCx5(1:(dim-1))];
                ANCCy = X(p, :)*ANCCx5';   
                ActualSPx5 = [ANCCy ActualSPx5(:, (1:(Tap_Sw-1)))];
                e_cont5(j) = Yp(j) - ActualSPx5*Sw';
                fitness5(p) = fitness5(p) + e_cont5(j)^2/block;
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
        for i = 1:pop5
    
            % Update alpha
            if fitness5(i) < Alpha_score 
                Alpha_score = fitness5(i); 
                Alpha_pos = X(i, :);
            end
    
            % Update beta
            if (fitness5(i) > Alpha_score) && (fitness5(i) < Beta_score) 
                Beta_score = fitness5(i); 
                Beta_pos = X(i, :);
            end
    
            % Update delta
            if (fitness5(i) > Alpha_score) && (fitness5(i) > Beta_score) && (fitness5(i) < Delta_score) 
                Delta_score = fitness5(i); 
                Delta_pos = X(i, :);
            end
        end
    
        IterCurve5(t) = IterCurve5(t) + Alpha_score;
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop5
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
end

%% Plot
plot(IterCurve1/q, 'LineWidth', 1.5);
hold on;
plot(IterCurve2/q, 'LineWidth', 1.5);
hold on;
plot(IterCurve3/q, 'LineWidth', 1.5);
hold on;
plot(IterCurve4/q, 'LineWidth', 1.5);
hold on;
plot(IterCurve5/q, 'LineWidth', 1.5);
xlabel('{\fontname{Times New Roman}Iteration Number}');
ylabel('{\fontname{Times New Roman}Best Fitness Value}');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
legend('{\fontname{Times New Roman}P=10}', '{\fontname{Times New Roman}P=20}', '{\fontname{Times New Roman}P=30}', '{\fontname{Times New Roman}P=40}', '{\fontname{Times New Roman}P=50}',...%    
     'NumColumns', 1, 'FontSize', 12, 'Location', 'northeast');
grid on
set(gca, 'GridLineStyle', ':');  % 设置为虚线 
set(gca, 'GridAlpha', 1);  % 设置透明度