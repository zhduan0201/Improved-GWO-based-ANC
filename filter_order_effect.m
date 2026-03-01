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

dim1 = 10;
dim2 = 20;
dim3 = 30;
dim4 = 40;
dim5 = 50;

pop = 20;
x_max = 0.1;
ub1 = x_max*ones(1, dim1);
lb1 = -x_max*ones(1, dim1);
ub2 = x_max*ones(1, dim2);
lb2 = -x_max*ones(1, dim2);
ub3 = x_max*ones(1, dim3);
lb3 = -x_max*ones(1, dim3);
ub4 = x_max*ones(1, dim4);
lb4 = -x_max*ones(1, dim4);
ub5 = x_max*ones(1, dim5);
lb5 = -x_max*ones(1, dim5);

% Primary noise and desired signal
maxIter = 30;
global DataLong; 
DataLong = pop*block*maxIter;      
snr = 30;    
[X_noise, ~] = add_awgn(sine_wave_generator(200)+sine_wave_generator(400), snr);
Yp=filter(Pw,1,X_noise); 

IterCurve1 = zeros(1, maxIter); 
IterCurve2 = zeros(1, maxIter);
IterCurve3 = zeros(1, maxIter);
IterCurve4 = zeros(1, maxIter);
IterCurve5 = zeros(1, maxIter);

%% Simulation start
q = 50; 
for n = 1:q
    %% GWO algorithm (filter order = 10)
    % Initialization    
    ANCCx1 = zeros(1, dim1);
    ActualSPx1 = zeros(1, Tap_Sw);  
    e_cont1 = zeros(1, DataLong);
    X = initialize_population(pop, ub1, lb1, dim1);
    
    %% Iteration start
    for t = 1:maxIter
        fitness1 = zeros(1, pop);
        for p = 1:pop
        
            % calculate the index range of the current particle
            start_index = pop*(t-1)*block + (p-1)*block + 1; 
            end_index = pop*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
            for j = start_index:end_index 
                ANCCx1 = [X_noise(j) ANCCx1(1:(dim1-1))];
                ANCCy = X(p, :)*ANCCx1';   
                ActualSPx1 = [ANCCy ActualSPx1(:, (1:(Tap_Sw-1)))];
                e_cont1(j) = Yp(j) - ActualSPx1*Sw';
                fitness1(p) = fitness1(p) + e_cont1(j)^2/block;
            end
        end
        
        % Initialize alpha, beta, and delta
        Alpha_pos = zeros(1, dim1);
        Beta_pos = zeros(1, dim1);
        Delta_pos = zeros(1, dim1);
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
    
        IterCurve1(t) = IterCurve1(t) + Alpha_score;
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop
            if isequal(X(i, :), Alpha_pos) || isequal(X(i, :), Beta_pos) || isequal(X(i, :), Delta_pos)
                continue
            end
    
            for j = 1:dim1                            
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
    
            X(i, :) = boundary_constraints(X(i, :), ub1, lb1, dim1);
        end
    end

    %% GWO algorithm (filter order = 20)
    % Initialization    
    ANCCx2 = zeros(1, dim2);
    ActualSPx2 = zeros(1, Tap_Sw);  
    e_cont2 = zeros(1, DataLong);
    X = initialize_population(pop, ub2, lb2, dim2);

    % Iteration start
    for t = 1:maxIter
        fitness2 = zeros(1, pop);
        for p = 1:pop
        
            % calculate the index range of the current particle
            start_index = pop*(t-1)*block + (p-1)*block + 1; 
            end_index = pop*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
            for j = start_index:end_index 
                ANCCx2 = [X_noise(j) ANCCx2(1:(dim2-1))];
                ANCCy = X(p, :)*ANCCx2';   
                ActualSPx2 = [ANCCy ActualSPx2(:, (1:(Tap_Sw-1)))];
                e_cont2(j) = Yp(j) - ActualSPx2*Sw';
                fitness2(p) = fitness2(p) + e_cont2(j)^2/block;
            end
        end
        
        % Initialize alpha, beta, and delta
        Alpha_pos = zeros(1, dim2);
        Beta_pos = zeros(1, dim2);
        Delta_pos = zeros(1, dim2);
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
    
        IterCurve2(t) = IterCurve2(t) + Alpha_score;
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop
            if isequal(X(i, :), Alpha_pos) || isequal(X(i, :), Beta_pos) || isequal(X(i, :), Delta_pos)
                continue
            end
    
            for j = 1:dim2                            
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
    
            X(i, :) = boundary_constraints(X(i, :), ub2, lb2, dim2);
        end      
    end

    %% GWO algorithm (filter order = 30)
    % Initialization    
    ANCCx3 = zeros(1, dim3);
    ActualSPx3 = zeros(1, Tap_Sw);  
    e_cont3 = zeros(1, DataLong);
    X = initialize_population(pop, ub3, lb3, dim3);

    % Iteration start
    for t = 1:maxIter
        fitness3 = zeros(1, pop);
        for p = 1:pop
        
            % calculate the index range of the current particle
            start_index = pop*(t-1)*block + (p-1)*block + 1; 
            end_index = pop*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
            for j = start_index:end_index 
                ANCCx3 = [X_noise(j) ANCCx3(1:(dim3-1))];
                ANCCy = X(p, :)*ANCCx3';   
                ActualSPx3 = [ANCCy ActualSPx3(:, (1:(Tap_Sw-1)))];
                e_cont3(j) = Yp(j) - ActualSPx3*Sw';
                fitness3(p) = fitness3(p) + e_cont3(j)^2/block;
            end
        end

        % Initialize alpha, beta, and delta
        Alpha_pos = zeros(1, dim3);
        Beta_pos = zeros(1, dim3);
        Delta_pos = zeros(1, dim3);
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
    
        IterCurve3(t) = IterCurve3(t) + Alpha_score;
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop
            if isequal(X(i, :), Alpha_pos) || isequal(X(i, :), Beta_pos) || isequal(X(i, :), Delta_pos)
                continue
            end
    
            for j = 1:dim3                            
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
    
            X(i, :) = boundary_constraints(X(i, :), ub3, lb3, dim3);
        end
    end

    %% GWO algorithm (filter order = 40)
    % Initialization   
    ANCCx4 = zeros(1, dim4);
    ActualSPx4 = zeros(1, Tap_Sw);  
    e_cont4 = zeros(1, DataLong);
    X = initialize_population(pop, ub4, lb4, dim4);

    % Iteration start
    for t = 1: maxIter
        fitness4 = zeros(1, pop);
        for p = 1:pop
        
            % calculate the index range of the current particle
            start_index = pop*(t-1)*block + (p-1)*block + 1; 
            end_index = pop*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
            for j = start_index:end_index 
                ANCCx4 = [X_noise(j) ANCCx4(1:(dim4-1))];
                ANCCy = X(p, :)*ANCCx4';   
                ActualSPx4 = [ANCCy ActualSPx4(:, (1:(Tap_Sw-1)))];
                e_cont4(j) = Yp(j) - ActualSPx4*Sw';
                fitness4(p) = fitness4(p) + e_cont4(j)^2/block;
            end
        end
       
        % Initialize alpha, beta, and delta
        Alpha_pos = zeros(1, dim4);
        Beta_pos = zeros(1, dim4);
        Delta_pos = zeros(1, dim4);
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
    
        IterCurve4(t) = IterCurve4(t) + Alpha_score;
    
        a = 2*(1-t/maxIter);
        
        for i = 1:pop
            if isequal(X(i, :), Alpha_pos) || isequal(X(i, :), Beta_pos) || isequal(X(i, :), Delta_pos)
                continue
            end
    
            for j = 1:dim4                            
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
    
            X(i, :) = boundary_constraints(X(i, :), ub4, lb4, dim4);
        end
    end

    %% GWO algorithm (filter order = 50)
    % Initialization    
    ANCCx5 = zeros(1, dim5);
    ActualSPx5 = zeros(1, Tap_Sw);  
    e_cont5 = zeros(1, DataLong);
    X = initialize_population(pop, ub5, lb5, dim5);

    % Iteration start
    for t = 1:maxIter
        fitness5 = zeros(1, pop);
        for p = 1:pop
        
            % calculate the index range of the current particle
            start_index = pop*(t-1)*block + (p-1)*block + 1; 
            end_index = pop*(t-1)*block + p*block;
        
            % traverse the data block of the current particle
            for j = start_index:end_index 
                ANCCx5 = [X_noise(j) ANCCx5(1:(dim5-1))];
                ANCCy = X(p, :)*ANCCx5';   
                ActualSPx5 = [ANCCy ActualSPx5(:, (1:(Tap_Sw-1)))];
                e_cont5(j) = Yp(j) - ActualSPx5*Sw';
                fitness5(p) = fitness5(p) + e_cont5(j)^2/block;
            end
        end


        % Initialize alpha, beta, and delta
        Alpha_pos = zeros(1, dim5);
        Beta_pos = zeros(1, dim5);
        Delta_pos = zeros(1, dim5);
        Alpha_score = inf; 
        Beta_score = inf; 
        Delta_score = inf; 
    
        % Update alpha, beta, and delta
        for i = 1:pop
    
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
        
        for i = 1:pop
            if isequal(X(i, :), Alpha_pos) || isequal(X(i, :), Beta_pos) || isequal(X(i, :), Delta_pos)
                continue
            end
    
            for j = 1:dim5                            
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
    
            X(i, :) = boundary_constraints(X(i, :), ub5, lb5, dim5);
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
legend('{\fontname{Times New Roman}N=10}', '{\fontname{Times New Roman}N=20}', '{\fontname{Times New Roman}N=30}', '{\fontname{Times New Roman}N=40}', '{\fontname{Times New Roman}N=50}',...%    
     'NumColumns', 1, 'FontSize', 12, 'Location', 'northeast');
grid on
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha', 1);