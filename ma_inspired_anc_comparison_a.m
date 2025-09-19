clc;clear all;close all;
global DataLong; 
%% 理想数据
%% 理想路径
Pw = [0 0 0 0 0 0 0 0 0 0.8 0.6 -0.2 -0.5 -0.1 0.4 -0.05];    
Sw = [0 0 0 0 0 1 2.5 1.76 0.15 -0.4825 -0.18625 -0.005 -0.001875]; 
Tap_Sw = length(Sw);  
q = 50;
%% 初级噪声
DataLong = 120000;      
snr = 30; 
[X_noise, ~] = add_awgn(sine_wave_generator(200) + sine_wave_generator(400), snr);
%% 期望信号
Yp=filter(Pw,1,X_noise);  
%% 关键参数
block=50;
dim=20; % 控制滤波器阶数 
pop=20; % 种群数量
%%
maxIter1 = DataLong / (pop * block); % 1000
maxIter2 = DataLong / (pop * block * 3); % 3000
maxIter3 = DataLong / (pop * block * 3 / 2); % 1500
a1 = zeros(1, maxIter1);
IterCurve1=zeros(1, maxIter1); 
IterCurve2=zeros(1, maxIter1);
IterCurve3=zeros(1, maxIter1);
IterCurve4=zeros(1, maxIter1);
ANR1 = zeros(q, DataLong);
ANR2 = zeros(q, DataLong);
ANR3 = zeros(q, DataLong);
ANR4 = zeros(q, DataLong);
Pe = zeros(q, DataLong);
%% 重复实验
for epoch = 1:q
%% 群智能优化算法参数设定
x_max = 0.1;
ub = x_max * ones(1, dim); % 位置范围
lb = -x_max * ones(1, dim);
%% ANR 
Me1 = zeros(1, DataLong);
Me2 = zeros(1, DataLong);
Me3 = zeros(1, DataLong);
Me4 = zeros(1, DataLong);
Md = zeros(1, DataLong);
%%
%% QPSO算法
%% 初始化  
alpha = 0.5; 
ANCCx1=zeros(1,dim);
ActualSPx1=zeros(1,Tap_Sw);  
e_cont1=zeros(1,DataLong);
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 迭代开始
for t = 1: maxIter1
fitness1=zeros(1,pop); 
%% 粒子1 
for j = (20 * (t - 1) * block + 1) : (20 * (t - 1) * block + block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(1)=fitness1(1)+e_cont1(j)^2/block;
end
%% 粒子2
for j = (20 * (t - 1) * block + block + 1) : (20 * (t - 1) * block + 2 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(2)=fitness1(2)+e_cont1(j)^2/block;
end
%% 粒子3
for j = (20 * (t - 1) * block + 2 * block + 1) : (20 * (t - 1) * block + 3 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(3)=fitness1(3)+e_cont1(j)^2/block;
end
%% 粒子4
for j = (20 * (t - 1) * block + 3 * block + 1) : (20 * (t - 1) * block + 4 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(4)=fitness1(4)+e_cont1(j)^2/block;
end
%% 粒子5
for j = (20 * (t - 1) * block + 4 * block + 1) : (20 * (t - 1) * block + 5 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(5)=fitness1(5)+e_cont1(j)^2/block;
end
%% 粒子6
for j = (20 * (t - 1) * block + 5 * block + 1) : (20 * (t - 1) * block + 6 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(6)=fitness1(6)+e_cont1(j)^2/block;
end
%% 粒子7
for j = (20 * (t - 1) * block + 6 * block + 1) : (20 * (t - 1) * block + 7 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(7)=fitness1(7)+e_cont1(j)^2/block;
end
%% 粒子8
for j = (20 * (t - 1) * block + 7 * block + 1) : (20 * (t - 1) * block + 8 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(8)=fitness1(8)+e_cont1(j)^2/block;
end
%% 粒子9
for j = (20 * (t - 1) * block + 8 * block + 1) : (20 * (t - 1) * block + 9 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(9)=fitness1(9)+e_cont1(j)^2/block;
end
%% 粒子10
for j = (20 * (t - 1) * block + 9 * block + 1) : (20 * (t - 1) * block + 10 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(10)=fitness1(10)+e_cont1(j)^2/block;
end
%% 粒子11
for j = (20 * (t - 1) * block + 10 * block + 1) : (20 * (t - 1) * block + 11 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(11)=fitness1(11)+e_cont1(j)^2/block;
end
%% 粒子12
for j = (20 * (t - 1) * block + 11 * block + 1) : (20 * (t - 1) * block + 12 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(12)=fitness1(12)+e_cont1(j)^2/block;
end
%% 粒子13
for j = (20 * (t - 1) * block + 12 * block + 1) : (20 * (t - 1) * block + 13 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(13)=fitness1(13)+e_cont1(j)^2/block;
end
%% 粒子14
for j = (20 * (t - 1) * block + 13 * block + 1) : (20 * (t - 1) * block + 14 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(14)=fitness1(14)+e_cont1(j)^2/block;
end
%% 粒子15
for j = (20 * (t - 1) * block + 14 * block + 1)  :(20 * (t - 1) * block + 15 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(15)=fitness1(15)+e_cont1(j)^2/block;
end
%% 粒子16
for j = (20 * (t - 1) * block + 15 * block + 1) : (20 * (t - 1) * block + 16 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(16)=fitness1(16)+e_cont1(j)^2/block;
end
%% 粒子17
for j = (20 * (t - 1) * block + 16 * block + 1) : (20 * (t - 1) * block + 17 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(17)=fitness1(17)+e_cont1(j)^2/block;
end
%% 粒子18
for j = (20 * (t - 1) * block + 17 * block + 1) : (20 * (t - 1) * block + 18 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(18)=fitness1(18)+e_cont1(j)^2/block;
end
%% 粒子19
for j = (20 * (t - 1) * block + 18 * block + 1) : (20 * (t - 1) * block + 19 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(19)=fitness1(19)+e_cont1(j)^2/block;
end
%% 粒子20
for j = (20 * (t - 1) * block + 19 * block + 1) : (20 * (t - 1) * block + 20 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(20)=fitness1(20)+e_cont1(j)^2/block;
end

if t == 1
%% 将初始种群作为历史最优值
    pBest=X;
    pBestFitness=fitness1;
%% 记录初始全局最优值
    [~,index]=min(fitness1);
    gBestFitness=fitness1(index);
    gBest=X(index,:);
else
    for i=1:pop
        % 更新历史最优值
        if fitness1(i)<pBestFitness(i)    
            pBest(i,:)=X(i,:);
            pBestFitness(i)=fitness1(i);
        end
        % 更新全局最优值
        if fitness1(i)<gBestFitness
            gBestFitness=fitness1(i);
            gBest=X(i,:);
        end
    end
end
C = sum(pBest, 1) / pop; 
%% 适应度函数
IterCurve1(t) = IterCurve1(t) + gBestFitness;
%% 对每个粒子进行更新
for i=1:pop
    for j = 1:dim
        phi = rand;
        p = phi * pBest(i, j) +(1 - phi) * gBest(j); 
        if rand < 0.5
            X(i, j) = p + alpha * abs(X(i, j) - C(j)) * log(1 / rand);
        else
            X(i, j) = p - alpha * abs(X(i, j) - C(j)) * log(1 / rand);
        end
    end
    % 位置边界检查及约束
    X(i,:)=boundary_constraints(X(i,:),ub,lb,dim);
end
end
%%
%% ABC算法
%% 初始化种群位置
VarSize = [1 dim];
empty_bee.Position = [];
empty_bee.Cost = 0;
npop = repmat(empty_bee, pop, 1);
X = initialize_population(pop,ub,lb,dim);
for i = 1:pop
    npop(i).Position = X(i, :); % 初始化蜜源
end
%% 初始化    
a = 1;                    % Acceleration Coefficient Upper Bound
L = round(0.6 * dim * pop); % Abandonment Limit Parameter (Trial Limit)
ANCCx2 = zeros(1, dim);
ActualSPx2 = zeros(1, Tap_Sw);  
e_cont2 = zeros(1, DataLong);
%% 迭代开始
for t = 1: maxIter2
%% 粒子1 
for j = (3 * pop * block * (t - 1) + 1) : (3 * pop * block * (t - 1) + block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(1).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(1).Cost=npop(1).Cost+e_cont2(j)^2/block;
end
%% 粒子2
for j = (3 * pop * block * (t - 1) + block + 1) : (3 * pop * block * (t - 1) + 2 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(2).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(2).Cost=npop(2).Cost+e_cont2(j)^2/block;
end
%% 粒子3
for j = (3 * pop * block * (t - 1) + 2 * block + 1) : (3 * pop * block * (t - 1) + 3 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(3).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(3).Cost=npop(3).Cost+e_cont2(j)^2/block;
end
%% 粒子4
for j = (3 * pop * block * (t - 1) + 3 * block + 1) : (3 * pop * block * (t - 1) + 4 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(4).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(4).Cost=npop(4).Cost+e_cont2(j)^2/block;
end
%% 粒子5
for j = (3 * pop * block * (t - 1) + 4 * block + 1) : (3 * pop * block * (t - 1) + 5 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(5).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(5).Cost=npop(5).Cost+e_cont2(j)^2/block;
end
%% 粒子6
for j = (3 * pop * block * (t - 1) + 5 * block + 1) : (3 * pop * block * (t - 1) + 6 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(6).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(6).Cost=npop(6).Cost+e_cont2(j)^2/block;
end
%% 粒子7
for j = (3 * pop * block * (t - 1) + 6 * block + 1) : (3 * pop * block * (t - 1) + 7 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(7).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(7).Cost=npop(7).Cost+e_cont2(j)^2/block;
end
%% 粒子8
for j = (3 * pop * block * (t - 1) + 7 * block + 1) : (3 * pop * block * (t - 1) + 8 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(8).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(8).Cost=npop(8).Cost+e_cont2(j)^2/block;
end
%% 粒子9
for j = (3 * pop * block * (t - 1) + 8 * block + 1) : (3 * pop * block * (t - 1) + 9 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(9).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(9).Cost=npop(9).Cost+e_cont2(j)^2/block;
end
%% 粒子10
for j = (3 * pop * block * (t - 1) + 9 * block + 1) : (3 * pop * block * (t - 1) + 10 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(10).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(10).Cost=npop(10).Cost+e_cont2(j)^2/block;
end
%% 粒子11
for j = (3 * pop * block * (t - 1) + 10 * block + 1) : (3 * pop * block * (t - 1) + 11 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(11).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(11).Cost=npop(11).Cost+e_cont2(j)^2/block;
end
%% 粒子12
for j = (3 * pop * block * (t - 1) + 11 * block + 1) : (3 * pop * block * (t - 1) + 12 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(12).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(12).Cost=npop(12).Cost+e_cont2(j)^2/block;
end
%% 粒子13
for j = (3 * pop * block * (t - 1) + 12 * block + 1) : (3 * pop * block * (t - 1) + 13 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(13).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(13).Cost=npop(13).Cost+e_cont2(j)^2/block;
end
%% 粒子14
for j = (3 * pop * block * (t - 1) + 13 * block + 1) : (3 * pop * block * (t - 1) + 14 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(14).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(14).Cost=npop(14).Cost+e_cont2(j)^2/block;
end
%% 粒子15
for j = (3 * pop * block * (t - 1) + 14 * block + 1) : (3 * pop * block * (t - 1) + 15 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(15).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(15).Cost=npop(15).Cost+e_cont2(j)^2/block;
end
%% 粒子16
for j = (3 * pop * block * (t - 1) + 15 * block + 1) : (3 * pop * block * (t - 1) + 16 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(16).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(16).Cost=npop(16).Cost+e_cont2(j)^2/block;
end
%% 粒子17
for j = (3 * pop * block * (t - 1) + 16 * block + 1) : (3 * pop * block * (t - 1) + 17 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(17).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(17).Cost=npop(17).Cost+e_cont2(j)^2/block;
end
%% 粒子18
for j = (3 * pop * block * (t - 1) + 17 * block + 1) : (3 * pop * block * (t - 1) + 18 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(18).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(18).Cost=npop(18).Cost+e_cont2(j)^2/block;
end
%% 粒子19
for j = (3 * pop * block * (t - 1) + 18 * block + 1) : (3 * pop * block * (t - 1) + 19 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(19).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(19).Cost=npop(19).Cost+e_cont2(j)^2/block;
end
%% 粒子20
for j = (3 * pop * block * (t - 1) + 19 * block + 1) : (3 * pop * block * (t - 1) + 20 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop(20).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop(20).Cost=npop(20).Cost+e_cont2(j)^2/block;
end
%% 寻优
BestSol.Cost = inf;
for i = 1:pop
    if npop(i).Cost <= BestSol.Cost
        BestSol = npop(i);
    end
end
IterCurve2(t) = IterCurve2(t) + BestSol.Cost;
%% 更新
% Abandonment Counter
C = zeros(pop, 1);
%% 雇佣蜂
npop1 = repmat(empty_bee, pop, 1);
for i = 1:pop        
    % Choose k randomly, not equal to i
    K = [1:i-1 i+1:pop];
    k = K(randi([1 numel(K)]));
       
    % Define Acceleration Coeff.
    phi = a*unifrnd(-1, +1, VarSize);
        
    % New Bee Position 
    npop1(i).Position = npop(i).Position + phi .* (npop(i).Position - npop(k).Position);
        
    % Apply Bounds
    npop1(i).Position = boundary_constraints(npop1(i).Position, ub, lb, dim); 
end

%% 粒子1 
for j = (pop * block + 3 * pop * block * (t - 1) + 1) : (pop * block + 3 * pop * block * (t - 1) + block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(1).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(1).Cost=npop1(1).Cost+e_cont2(j)^2/block;
end
%% 粒子2
for j = (pop * block + 3 * pop * block * (t - 1) + block + 1) : (pop * block + 3 * pop * block * (t - 1) + 2 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(2).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(2).Cost=npop1(2).Cost+e_cont2(j)^2/block;
end
%% 粒子3
for j = (pop * block + 3 * pop * block * (t - 1) + 2 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 3 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(3).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(3).Cost=npop1(3).Cost+e_cont2(j)^2/block;
end
%% 粒子4
for j = (pop * block + 3 * pop * block * (t - 1) + 3 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 4 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(4).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(4).Cost=npop1(4).Cost+e_cont2(j)^2/block;
end
%% 粒子5
for j = (pop * block + 3 * pop * block * (t - 1) + 4 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 5 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(5).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(5).Cost=npop1(5).Cost+e_cont2(j)^2/block;
end
%% 粒子6
for j = (pop * block + 3 * pop * block * (t - 1) + 5 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 6 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(6).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(6).Cost=npop1(6).Cost+e_cont2(j)^2/block;
end
%% 粒子7
for j = (pop * block + 3 * pop * block * (t - 1) + 6 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 7 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(7).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(7).Cost=npop1(7).Cost+e_cont2(j)^2/block;
end
%% 粒子8
for j = (pop * block + 3 * pop * block * (t - 1) + 7 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 8 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(8).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(8).Cost=npop1(8).Cost+e_cont2(j)^2/block;
end
%% 粒子9
for j = (pop * block + 3 * pop * block * (t - 1) + 8 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 9 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(9).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(9).Cost=npop1(9).Cost+e_cont2(j)^2/block;
end
%% 粒子10
for j = (pop * block + 3 * pop * block * (t - 1) + 9 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 10 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(10).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(10).Cost=npop1(10).Cost+e_cont2(j)^2/block;
end
%% 粒子11
for j = (pop * block + 3 * pop * block * (t - 1) + 10 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 11 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(11).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(11).Cost=npop1(11).Cost+e_cont2(j)^2/block;
end
%% 粒子12
for j = (pop * block + 3 * pop * block * (t - 1) + 11 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 12 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(12).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(12).Cost=npop1(12).Cost+e_cont2(j)^2/block;
end
%% 粒子13
for j = (pop * block + 3 * pop * block * (t - 1) + 12 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 13 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(13).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(13).Cost=npop1(13).Cost+e_cont2(j)^2/block;
end
%% 粒子14
for j = (pop * block + 3 * pop * block * (t - 1) + 13 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 14 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(14).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(14).Cost=npop1(14).Cost+e_cont2(j)^2/block;
end
%% 粒子15
for j = (pop * block + 3 * pop * block * (t - 1) + 14 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 15 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(15).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(15).Cost=npop1(15).Cost+e_cont2(j)^2/block;
end
%% 粒子16
for j = (pop * block + 3 * pop * block * (t - 1) + 15 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 16 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(16).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(16).Cost=npop1(16).Cost+e_cont2(j)^2/block;
end
%% 粒子17
for j = (pop * block + 3 * pop * block * (t - 1) + 16 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 17 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(17).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(17).Cost=npop1(17).Cost+e_cont2(j)^2/block;
end
%% 粒子18
for j = (pop * block + 3 * pop * block * (t - 1) + 17 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 18 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(18).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(18).Cost=npop1(18).Cost+e_cont2(j)^2/block;
end
%% 粒子19
for j = (pop * block + 3 * pop * block * (t - 1) + 18 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 19 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(19).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(19).Cost=npop1(19).Cost+e_cont2(j)^2/block;
end
%% 粒子20
for j = (pop * block + 3 * pop * block * (t - 1) + 19 * block + 1) : (pop * block + 3 * pop * block * (t - 1) + 20 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop1(20).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop1(20).Cost=npop1(20).Cost+e_cont2(j)^2/block;
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
%% 跟随蜂
npop2 = repmat(empty_bee, pop, 1);
for m = 1:pop
        
    % Select Source Site
    i = RouletteWheelSelection(P);
        
    % Choose k randomly, not equal to i
    K = [1:i-1 i+1:pop];
    k = K(randi([1 numel(K)]));
        
    % Define Acceleration Coeff.
    phi = a*unifrnd(-1, +1, VarSize);
        
    % New Bee Position
    npop2(m).Position = npop(i).Position + phi .* (npop(i).Position - npop(k).Position);
        
    % Apply Bounds
    npop2(m).Position = boundary_constraints(npop2(m).Position, ub, lb, dim); 
end

%% 粒子1 
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(1).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(1).Cost=npop2(1).Cost+e_cont2(j)^2/block;
end
%% 粒子2
for j = (2 * pop * block + 3 * pop * block * (t - 1) + block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 2 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(2).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(2).Cost=npop2(2).Cost+e_cont2(j)^2/block;
end
%% 粒子3
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 2 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 3 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(3).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(3).Cost=npop2(3).Cost+e_cont2(j)^2/block;
end
%% 粒子4
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 3 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 4 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(4).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(4).Cost=npop2(4).Cost+e_cont2(j)^2/block;
end
%% 粒子5
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 4 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 5 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(5).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(5).Cost=npop2(5).Cost+e_cont2(j)^2/block;
end
%% 粒子6
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 5 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 6 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(6).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(6).Cost=npop2(6).Cost+e_cont2(j)^2/block;
end
%% 粒子7
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 6 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 7 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(7).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(7).Cost=npop2(7).Cost+e_cont2(j)^2/block;
end
%% 粒子8
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 7 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 8 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(8).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(8).Cost=npop2(8).Cost+e_cont2(j)^2/block;
end
%% 粒子9
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 8 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 9 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(9).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(9).Cost=npop2(9).Cost+e_cont2(j)^2/block;
end
%% 粒子10
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 9 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 10 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(10).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(10).Cost=npop2(10).Cost+e_cont2(j)^2/block;
end
%% 粒子11
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 10 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 11 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(11).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(11).Cost=npop2(11).Cost+e_cont2(j)^2/block;
end
%% 粒子12
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 11 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 12 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(12).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(12).Cost=npop2(12).Cost+e_cont2(j)^2/block;
end
%% 粒子13
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 12 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 13 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(13).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(13).Cost=npop2(13).Cost+e_cont2(j)^2/block;
end
%% 粒子14
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 13 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 14 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(14).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(14).Cost=npop2(14).Cost+e_cont2(j)^2/block;
end
%% 粒子15
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 14 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 15 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(15).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(15).Cost=npop2(15).Cost+e_cont2(j)^2/block;
end
%% 粒子16
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 15 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 16 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(16).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(16).Cost=npop2(16).Cost+e_cont2(j)^2/block;
end
%% 粒子17
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 16 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 17 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(17).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(17).Cost=npop2(17).Cost+e_cont2(j)^2/block;
end
%% 粒子18
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 17 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 18 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(18).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(18).Cost=npop2(18).Cost+e_cont2(j)^2/block;
end
%% 粒子19
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 18 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 19 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(19).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(19).Cost=npop2(19).Cost+e_cont2(j)^2/block;
end
%% 粒子20
for j = (2 * pop * block + 3 * pop * block * (t - 1) + 19 * block + 1) : (2 * pop * block + 3 * pop * block * (t - 1) + 20 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=npop2(20).Position*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    npop2(20).Cost=npop2(20).Cost+e_cont2(j)^2/block;
end

for i = 1:pop
    if npop2(i).Cost <= npop(i).Cost
        npop(i) = npop2(i);
    else
        C(i) = C(i)+1;
    end
end

end
%%
%% RGA
%% 初始化种群位置
X = initialize_population(3 * pop / 2, ub, lb, dim);
%% 初始化    
mut = 0.1; % 突变概率
acr = 0.9; % 交叉概率
eta_c = 20;
eta_m = 20;
ANCCx3 = zeros(1, dim);
ActualSPx3 = zeros(1, Tap_Sw);  
e_cont3 = zeros(1, DataLong);
%% 迭代开始
for t = 1: maxIter3
fitness3=zeros(1, 3 * pop / 2); 
%% 粒子1 
for j = (30 * (t - 1) * block + 1) : (30 * (t - 1) * block + block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(1)=fitness3(1)+e_cont3(j)^2/block;
end
%% 粒子2
for j = (30 * (t - 1) * block + block + 1) : (30 * (t - 1) * block + 2 * block) 
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(2)=fitness3(2)+e_cont3(j)^2/block;
end
%% 粒子3
for j = (30 * (t - 1) * block + 2 * block + 1) : (30 * (t - 1) * block + 3 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(3)=fitness3(3)+e_cont3(j)^2/block;
end
%% 粒子4
for j = (30 * (t - 1) * block + 3 * block + 1) : (30 * (t - 1) * block + 4 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(4)=fitness3(4)+e_cont3(j)^2/block;
end
%% 粒子5
for j = (30 * (t - 1) * block + 4 * block + 1) : (30 * (t - 1) * block + 5 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(5)=fitness3(5)+e_cont3(j)^2/block;
end
%% 粒子6
for j = (30 * (t - 1) * block + 5 * block + 1) : (30 * (t - 1) * block + 6 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(6)=fitness3(6)+e_cont3(j)^2/block;
end
%% 粒子7
for j = (30 * (t - 1) * block + 6 * block + 1) : (30 * (t - 1) * block + 7 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(7)=fitness3(7)+e_cont3(j)^2/block;
end
%% 粒子8
for j = (30 * (t - 1) * block + 7 * block + 1) : (30 * (t - 1) * block + 8 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(8)=fitness3(8)+e_cont3(j)^2/block;
end
%% 粒子9
for j = (30 * (t - 1) * block + 8 * block + 1) : (30 * (t - 1) * block + 9 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(9)=fitness3(9)+e_cont3(j)^2/block;
end
%% 粒子10
for j = (30 * (t - 1) * block + 9 * block + 1) : (30 * (t - 1) * block + 10 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(10)=fitness3(10)+e_cont3(j)^2/block;
end
%% 粒子11
for j = (30 * (t - 1) * block + 10 * block + 1) : (30 * (t - 1) * block + 11 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(11)=fitness3(11)+e_cont3(j)^2/block;
end
%% 粒子12
for j = (30 * (t - 1) * block + 11 * block + 1) : (30 * (t - 1) * block + 12 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(12)=fitness3(12)+e_cont3(j)^2/block;
end
%% 粒子13
for j = (30 * (t - 1) * block + 12 * block + 1) : (30 * (t - 1) * block + 13 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(13)=fitness3(13)+e_cont3(j)^2/block;
end
%% 粒子14
for j = (30 * (t - 1) * block + 13 * block + 1) : (30 * (t - 1) * block + 14 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(14)=fitness3(14)+e_cont3(j)^2/block;
end
%% 粒子15
for j = (30 * (t - 1) * block + 14 * block + 1) : (30 * (t - 1) * block + 15 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(15)=fitness3(15)+e_cont3(j)^2/block;
end
%% 粒子16
for j = (30 * (t - 1) * block + 15 * block + 1) : (30 * (t - 1) * block + 16 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(16)=fitness3(16)+e_cont3(j)^2/block;
end
%% 粒子17
for j = (30 * (t - 1) * block + 16 * block + 1) : (30 * (t - 1) * block + 17 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(17)=fitness3(17)+e_cont3(j)^2/block;
end
%% 粒子18
for j = (30 * (t - 1) * block + 17 * block + 1) : (30 * (t - 1) * block + 18 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(18)=fitness3(18)+e_cont3(j)^2/block;
end
%% 粒子19
for j = (30 * (t - 1) * block + 18 * block + 1) : (30 * (t - 1) * block + 19 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(19)=fitness3(19)+e_cont3(j)^2/block;
end
%% 粒子20
for j = (30 * (t - 1) * block + 19 * block + 1) : (30 * (t - 1) * block + 20 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(20)=fitness3(20)+e_cont3(j)^2/block;
end
%% 计算最佳适应度值
chrome_pos = zeros(1, dim + 1);
chrome_pos(end) = inf;
for i = 1 : pop
if fitness3(i)<chrome_pos(end) 
    chrome_pos(end)=fitness3(i); 
    chrome_pos(1:dim)=X(i,:); 
end
end
IterCurve3(t) = IterCurve3(t) + chrome_pos(end); 
%% 选择
for i = pop + 1 : 3 * pop / 2
    idx1 = randi(pop);
    idx2 = randi(pop);

    if fitness3(idx1) < fitness3(idx2)
        X(i, :) = X(idx1, :);
    else
        X(i, :) = X(idx2, :);
    end
end
%% 交叉
parent1 = zeros(1, dim);
parent2 = zeros(1, dim);
child1 = zeros(1, dim);
child2 = zeros(1, dim);

for i = pop + 1 : 2 : 3 * pop / 2
    parent1 = X(i, :);
    parent2 = X(i + 1, :);

    for j = 1:dim

        acr_rand = rand; % 生成一个代表该个体是否产生交叉的概率大小，用于判别是否进行交叉处理
        if acr_rand < acr % 如果该个体的交叉概率值大于产生交叉处理的阈值，则对该个体的染色体（两条，因为此案例中有两个自变量）进行交叉处理

            r_h = rand;      
            if r_h <= 0.5
                alpha_h = (2 * r_h)^(1 / (eta_c + 1));
            else
                alpha_h = (1 / (2 * (1 - r_h)))^(1 / (eta_c + 1));
            end

            child1(j) = 0.5 * ((1 - alpha_h) * parent1(j) + (1 + alpha_h) * parent2(j));
            child2(j) = 0.5 * ((1 + alpha_h) * parent1(j) + (1 - alpha_h) * parent2(j));

        else
            child1(j) = parent1(j);
            child2(j) = parent2(j);
        end
    end

    X(i, :) = child1;
    if i + 1 <= 3 * pop / 2
        X(i + 1, :) = child2;
    end

end

%% 变异
for i = pop + 1 : 3 * pop / 2 %%N是个体总数，也就是每一代有多少头袋鼠
    for j = 1 : dim
        mut_rand = rand; %随机生成一个数，代表自然里的基因突变，然后用改值来决定是否产生突变。
        if mut_rand < mut  %mut代表突变概率，即产生突变的阈值，如果小于0.2的基因突变概率阈值才进行基因突变处理，否者不进行突变处理
            r_h = rand;        
            if r_h <= 0.5
                delta_h = (2 * r_h)^(1 / (eta_m + 1)) - 1;
            else
                delta_h = 1 - (2 * (1 - r_h))^(1 / (eta_m + 1));
            end

            X(i, j) = X(i, j) + delta_h;
        end
    end
    X(i, :) = boundary_constraints(X(i, :), ub, lb, dim);
end

%% 子代1
for j = (30 * (t - 1) * block + 20 * block + 1) : (30 * (t - 1) * block + 21 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(21,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(21)=fitness3(21)+e_cont3(j)^2/block;
end
%% 子代2
for j = (30 * (t - 1) * block + 21 * block + 1) : (30 * (t - 1) * block + 22 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(22,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(22)=fitness3(22)+e_cont3(j)^2/block;
end
%% 子代3
for j = (30 * (t - 1) * block + 22 * block + 1) : (30 * (t - 1) * block + 23 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(23,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(23)=fitness3(23)+e_cont3(j)^2/block;
end
%% 子代4
for j = (30 * (t - 1) * block + 23 * block + 1) : (30 * (t - 1) * block + 24 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(24,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(24)=fitness3(24)+e_cont3(j)^2/block;
end
%% 子代5
for j = (30 * (t - 1) * block + 24 * block + 1) : (30 * (t - 1) * block + 25 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(25,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(25)=fitness3(25)+e_cont3(j)^2/block;
end
%% 子代6
for j = (30 * (t - 1) * block + 25 * block + 1) : (30 * (t - 1) * block + 26 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(26,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(26)=fitness3(26)+e_cont3(j)^2/block;
end
%% 子代7
for j = (30 * (t - 1) * block + 26 * block + 1) : (30 * (t - 1) * block + 27 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(27,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(27)=fitness3(27)+e_cont3(j)^2/block;
end
%% 子代8
for j = (30 * (t - 1) * block + 27 * block + 1) : (30 * (t - 1) * block + 28 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(28,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(8)=fitness3(8)+e_cont3(j)^2/block;
end
%% 子代9
for j = (30 * (t - 1) * block + 28 * block + 1) : (30 * (t - 1) * block + 29 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(29,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(29)=fitness3(29)+e_cont3(j)^2/block;
end
%% 子代10 
for j = (30 * (t - 1) * block + 29 * block + 1) : (30 * (t - 1) * block + 30 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(30,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(30)=fitness3(30)+e_cont3(j)^2/block;
end

% 步骤 1: 根据适应度值排序
[sorted_fitness, idx] = sort(fitness3);

% 步骤 2: 选取前20个适应度最小的向量
top_20_vectors = X(idx(1:20), :);

for i = 1 : pop
    X(i, :) = top_20_vectors(i, :);
end
end
%%
%% GWO算法
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 初始化    
ANCCx4=zeros(1,dim);
ActualSPx4=zeros(1,Tap_Sw);  
e_cont4=zeros(1,DataLong);
D = zeros(1, maxIter1);
%% 迭代开始
for t = 1: maxIter1
fitness4=zeros(1,pop); 
%% 粒子1 
for j = (20 * (t - 1) * block + 1) : (20 * (t - 1) * block + block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(1)=fitness4(1)+e_cont4(j)^2/block;
end
%% 粒子2
for j = (20 * (t - 1) * block + block + 1) : (20 * (t - 1) * block + 2 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(2)=fitness4(2)+e_cont4(j)^2/block;
end
%% 粒子3
for j = (20 * (t - 1) * block + 2 * block + 1) : (20 * (t - 1) * block + 3 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(3)=fitness4(3)+e_cont4(j)^2/block;
end
%% 粒子4
for j = (20 * (t - 1) * block + 3 * block + 1) : (20 * (t - 1) * block + 4 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(4)=fitness4(4)+e_cont4(j)^2/block;
end
%% 粒子5
for j = (20 * (t - 1) * block + 4 * block + 1) : (20 * (t - 1) * block + 5 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(5)=fitness4(5)+e_cont4(j)^2/block;
end
%% 粒子6
for j = (20 * (t - 1) * block + 5 * block + 1) : (20 * (t - 1) * block + 6 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(6)=fitness4(6)+e_cont4(j)^2/block;
end
%% 粒子7
for j = (20 * (t - 1) * block + 6 * block + 1) : (20 * (t - 1) * block + 7 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(7)=fitness4(7)+e_cont4(j)^2/block;
end
%% 粒子8
for j = (20 * (t - 1) * block + 7 * block + 1) : (20 * (t - 1) * block + 8 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(8)=fitness4(8)+e_cont4(j)^2/block;
end
%% 粒子9
for j = (20 * (t - 1) * block + 8 * block + 1) : (20 * (t - 1) * block + 9 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(9)=fitness4(9)+e_cont4(j)^2/block;
end
%% 粒子10
for j = (20 * (t - 1) * block + 9 * block + 1) : (20 * (t - 1) * block + 10 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(10)=fitness4(10)+e_cont4(j)^2/block;
end
%% 粒子11
for j = (20 * (t - 1) * block + 10 * block + 1) : (20 * (t - 1) * block + 11 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(11)=fitness4(11)+e_cont4(j)^2/block;
end
%% 粒子12
for j = (20 * (t - 1) * block + 11 * block + 1) : (20 * (t - 1) * block + 12 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(12)=fitness4(12)+e_cont4(j)^2/block;
end
%% 粒子13
for j = (20 * (t - 1) * block + 12 * block + 1) : (20 * (t - 1) * block + 13 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(13)=fitness4(13)+e_cont4(j)^2/block;
end
%% 粒子14
for j = (20 * (t - 1) * block + 13 * block + 1) : (20 * (t - 1) * block + 14 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(14)=fitness4(14)+e_cont4(j)^2/block;
end
%% 粒子15
for j = (20 * (t - 1) * block + 14 * block + 1) : (20 * (t - 1) * block + 15 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(15)=fitness4(15)+e_cont4(j)^2/block;
end
%% 粒子16
for j = (20 * (t - 1) * block + 15 * block + 1) : (20 * (t - 1) * block + 16 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(16)=fitness4(16)+e_cont4(j)^2/block;
end
%% 粒子17
for j = (20 * (t - 1) * block + 16 * block + 1) : (20 * (t - 1) * block + 17 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(17)=fitness4(17)+e_cont4(j)^2/block;
end
%% 粒子18
for j = (20 * (t - 1) * block + 17 * block + 1) : (20 * (t - 1) * block + 18 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(18)=fitness4(18)+e_cont4(j)^2/block;
end
%% 粒子19
for j = (20 * (t - 1) * block + 18 * block + 1) : (20 * (t - 1) * block + 19 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(19)=fitness4(19)+e_cont4(j)^2/block;
end
%% 粒子20
for j = (20 * (t - 1) * block + 19 * block + 1) : (20 * (t - 1) * block + 20 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(20)=fitness4(20)+e_cont4(j)^2/block;
end
%% 初始化alpha，beta和delta
Alpha_pos=zeros(1,dim);
Beta_pos=zeros(1,dim);
Delta_pos=zeros(1,dim);
Alpha_score=inf; 
Beta_score=inf; 
Delta_score=inf; 
%% 更新Alpha，Beta和Delta
% 更新alpha
for i=1:pop
if fitness4(i)<Alpha_score 
    Alpha_score=fitness4(i); 
    Alpha_pos=X(i,:);
end
% 更新beta
if fitness4(i)>Alpha_score && fitness4(i)<Beta_score 
    Beta_score=fitness4(i); 
    Beta_pos=X(i,:);
end
% 更新delta
if fitness4(i)>Alpha_score && fitness4(i)>Beta_score && fitness4(i)<Delta_score 
    Delta_score=fitness4(i); 
    Delta_pos=X(i,:);
end
end
max_score = max(fitness4);
IterCurve4(t) = IterCurve4(t) + Alpha_score;
%% 更新
X_avg = mean(X, 1); % 逐列求平均，得到 1×d 的向量
distances = zeros(pop, 1); % 初始化距离存储
for i = 1:pop
    distances(i) = norm(X(i, :) - X_avg); % 欧氏距离
end
D(t) = mean(distances);

eta = 0.2;
b = 0.1;
if D(t) > eta 
    a1(t) = 2 * (1 - t / maxIter1);
    T = t;
    a_mid = a1(t);
else
    a1(t) = a_mid * exp(-b * (t - T));
end

for i = 1:pop
if isequal(X(i,:), Alpha_pos) || isequal(X(i,:), Beta_pos) || isequal(X(i,:), Delta_pos)
    continue
end
for j=1:dim     
%%                       
r1=rand(); % r1 is a random number in [0,1]
r2=rand(); % r2 is a random number in [0,1]
            
A1=2*a1(t)*r1-a1(t); % Equation (3.3)
C1=2*r2; % Equation (3.4)
            
D_alpha=abs(C1*Alpha_pos(j)-X(i,j)); % Equation (3.5)-part 1
X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
%%                       
r1=rand();
r2=rand();
            
A2=2*a1(t)*r1-a1(t); % Equation (3.3)
C2=2*r2; % Equation (3.4)
            
D_beta=abs(C2*Beta_pos(j)-X(i,j)); % Equation (3.5)-part 2
X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
%%            
r1=rand();
r2=rand(); 
            
A3=2*a1(t)*r1-a1(t); % Equation (3.3)
C3=2*r2; % Equation (3.4)
            
D_delta=abs(C3*Delta_pos(j)-X(i,j)); % Equation (3.5)-part 3
X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
%%
wa = 1 / 3;
wb = 1 / 3;
wc = 1 / 3;
X(i, j) = wa * X1 + wb * X2 + wc * X3;% Equation (3.7)
end
X(i,:)=boundary_constraints(X(i,:),ub,lb,dim);
end
end
%% ANR
lambda = 0.999;
for i = 1:DataLong 
    if i == 1
        Me1(i) = (1 - lambda) * abs(e_cont1(i));
        Me2(i) = (1 - lambda) * abs(e_cont2(i));
        Me3(i) = (1 - lambda) * abs(e_cont3(i));
        Me4(i) = (1 - lambda) * abs(e_cont4(i));
        Md(i) = (1 - lambda) * abs(Yp(i));
    else
        Me1(i) = lambda * Me1(i - 1) + (1 - lambda) * abs(e_cont1(i));
        Me2(i) = lambda * Me2(i - 1) + (1 - lambda) * abs(e_cont2(i));
        Me3(i) = lambda * Me3(i - 1) + (1 - lambda) * abs(e_cont3(i));
        Me4(i) = lambda * Me4(i - 1) + (1 - lambda) * abs(e_cont4(i));
        Md(i) = lambda * Md(i - 1) +(1 - lambda) * abs(Yp(i));
    end
    ANR1(epoch, i) = 20 * log10(Md(i) / Me1(i));
    ANR2(epoch, i) = 20 * log10(Md(i) / Me2(i));
    ANR3(epoch, i) = 20 * log10(Md(i) / Me3(i));
    ANR4(epoch, i) = 20 * log10(Md(i) / Me4(i));
end
end
%% 平均
AANR1 = mean(ANR1, 1);
AANR2 = mean(ANR2, 1);
AANR3 = mean(ANR3, 1);
AANR4 = mean(ANR4, 1);
% 高斯滤波法
windowSize = 5000; % 窗口大小
weights = gausswin(windowSize);
AANR11 = conv(AANR1, weights, 'same') / sum(weights);
AANR22 = conv(AANR2, weights, 'same') / sum(weights);
AANR33 = conv(AANR3, weights, 'same') / sum(weights);
AANR44 = conv(AANR4, weights, 'same') / sum(weights);
%% 画图
%% 降噪量曲线图
figure(1);
plot(AANR11,'LineWidth',1.5);
hold on;
plot(AANR22,'LineWidth',1.5);
hold on;
plot(AANR33,'LineWidth',1.5);
hold on;
plot(AANR44,'LineWidth',1.5);
% 获取当前轴
ax = gca;
% 设置 X 和 Y 轴的刻度标签格式
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
% 或者使用 xtickformat 和 ytickformat
xtickformat('%.0f');    
ytickformat('%.0f');
xlim([windowSize/2, DataLong-windowSize]);
ylim([-15, 35]);        
xlabel('{\fontname{Times New Roman}Sample Number}');
ylabel('{\fontname{Times New Roman}ANR}{\fontname{Times New Roman}/dB}');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
legend('{\fontname{Times New Roman}QPSO}','{\fontname{Times New Roman}ABC}','{\fontname{Times New Roman}RGA}','{\fontname{Times New Roman}IGWO}',...%    
     'NumColumns',1,'FontSize', 12, 'Location', 'southeast');
grid on
set(gca, 'GridLineStyle', ':');  % 设置为虚线 
set(gca, 'GridAlpha', 1);  % 设置透明度