clc;clear all;close all;
global DataLong; 
%% 理想数据
%% 理想路径
Pw = [0 0 0 0 0 0 0 0 0 0.8 0.6 -0.2 -0.5 -0.1 0.4 -0.05];    
Sw1 = [0 0 0 0 0 1 2.5 1.76 0.15 -0.4825 -0.18625 -0.005 -0.001875]; 
Sw2 = -Sw1; 
Tap_Sw=length(Sw1);
q=50;
%% 初级噪声    
DataLong = 60100;       
snr = 30;     
[X_noise,~]=add_awgn(sine_wave_generator(200)+sine_wave_generator(400),snr);
%% 期望信号
Yp=filter(Pw,1,X_noise);  
%% 关键参数
block=50;
dim=20; % 控制滤波器阶数 
pop=20; % 种群数量
%%
maxIter = (DataLong - 100) / (pop * block);
IterCurve1=zeros(1, maxIter); 
IterCurve2=zeros(1, maxIter);
IterCurve22=zeros(1, maxIter);
IterCurve3=zeros(1, maxIter);
a1 = zeros(1, maxIter);
a2 = zeros(1, maxIter);
a3 = zeros(1, maxIter);
lambda = 0.999;
ANR1 = zeros(q, DataLong);
ANR2 = zeros(q, DataLong);
ANR3 = zeros(q, DataLong);
%% 重复实验
for epoch=1:q
%% 群智能优化算法参数设定
x_max = 0.1;
ub = x_max * ones(1, dim); % 位置范围
lb = -x_max * ones(1, dim);
%% 
Me1 = zeros(1, DataLong);
Me2 = zeros(1, DataLong);
Me3 = zeros(1, DataLong);
Md = zeros(1, DataLong);
%%
%% GWO算法
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 初始化  
ff = 0;
alpha = 0.98;
beta = 3;
ANCCx1=zeros(1,dim);
ActualSPx1=zeros(1,Tap_Sw);  
e_cont1=zeros(1,DataLong);
D1 = zeros(1, maxIter);
Pe = zeros(1, DataLong);
flag = zeros(1, DataLong);
%% 迭代开始
for t = 1: maxIter
fitness4=zeros(1,pop); 
%% 粒子1 
for j = (20 * (t - 1) * block + 1) : (20 * (t - 1) * block + block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(1)=fitness4(1)+e_cont1(j)^2/block;
    if j == 1
        Pe(j) = alpha + (1 - alpha) * (e_cont1(j)^2);
    else
        Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    end
    if j > 1
        if Pe(j) > beta * Pe(j - 1)
            flag(j) = 1;
            J = j;
            ff = 1;
        end
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子2
for j = (20 * (t - 1) * block + block + 1) : (20 * (t - 1) * block + 2 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(2)=fitness4(2)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子3
for j = (20 * (t - 1) * block + 2 * block + 1) : (20 * (t - 1) * block + 3 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(3)=fitness4(3)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子4
for j = (20 * (t - 1) * block + 3 * block + 1) : (20 * (t - 1) * block + 4 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(4)=fitness4(4)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子5
for j = (20 * (t - 1) * block + 4 * block + 1) : (20 * (t - 1) * block + 5 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(5)=fitness4(5)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子6
for j = (20 * (t - 1) * block + 5 * block + 1) : (20 * (t - 1) * block + 6 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(6)=fitness4(6)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子7
for j = (20 * (t - 1) * block + 6 * block + 1) : (20 * (t - 1) * block + 7 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(7)=fitness4(7)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子8
for j = (20 * (t - 1) * block + 7 * block + 1) : (20 * (t - 1) * block + 8 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(8)=fitness4(8)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子9
for j = (20 * (t - 1) * block + 8 * block + 1) : (20 * (t - 1) * block + 9 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(9)=fitness4(9)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子10
for j = (20 * (t - 1) * block + 9 * block + 1) : (20 * (t - 1) * block + 10 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(10)=fitness4(10)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子11
for j = (20 * (t - 1) * block + 10 * block + 1) : (20 * (t - 1) * block + 11 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(11)=fitness4(11)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子12
for j = (20 * (t - 1) * block + 11 * block + 1) : (20 * (t - 1) * block + 12 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(12)=fitness4(12)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子13
for j = (20 * (t - 1) * block + 12 * block + 1) : (20 * (t - 1) * block + 13 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(13)=fitness4(13)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子14
for j = (20 * (t - 1) * block + 13 * block + 1) : (20 * (t - 1) * block + 14 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(14)=fitness4(14)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子15
for j = (20 * (t - 1) * block + 14 * block + 1) : (20 * (t - 1) * block + 15 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(15)=fitness4(15)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子16
for j = (20 * (t - 1) * block + 15 * block + 1) : (20 * (t - 1) * block + 16 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(16)=fitness4(16)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子17
for j = (20 * (t - 1) * block + 16 * block + 1) : (20 * (t - 1) * block + 17 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(17)=fitness4(17)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子18
for j = (20 * (t - 1) * block + 17 * block + 1) : (20 * (t - 1) * block + 18 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(18)=fitness4(18)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子19
for j = (20 * (t - 1) * block + 18 * block + 1) : (20 * (t - 1) * block + 19 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(19)=fitness4(19)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
end
%% 粒子20
for j = (20 * (t - 1) * block + 19 * block + 1) : (20 * (t - 1) * block + 20 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness4(20)=fitness4(20)+e_cont1(j)^2/block;
    Pe(j) = alpha * Pe(j - 1) + (1 - alpha) * (e_cont1(j)^2);
    if Pe(j) > beta * Pe(j - 1)
        flag(j) = 1;
        J = j;
        ff = 1;
    end
    if ff == 1
        break
    end
end
if ff == 1
    TT = t;
    break
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
IterCurve1(t)=IterCurve1(t)+Alpha_score;
%% 更新
% 计算种群的平均位置向量
X_avg = mean(X, 1); % 逐列求平均，得到 1×d 的向量
% 计算每个个体到种群平均位置的欧氏距离
distances = zeros(pop, 1); % 初始化距离存储
for i = 1:pop
    distances(i) = norm(X(i, :) - X_avg); % 欧氏距离
end
% 计算种群多样性（平均距离）
D1(t) = mean(distances);

eta = 0.1;
b = 0.1;
if D1(t) > eta 
    a1(t) = 2 * (1 - t / maxIter);
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
%%
%%
X(i,:)=boundary_constraints(X(i,:),ub,lb,dim);
end
end

%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 迭代开始
for t = TT: maxIter
fitness1=zeros(1,pop); 
%% 粒子1 
for j = (J + 20 * (t - TT) * block + 1) : (J + 20 * (t - TT) * block + block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(1)=fitness1(1)+e_cont1(j)^2/block;
end
%% 粒子2
for j = (J + 20 * (t - TT) * block + block + 1) : (J + 20 * (t - TT) * block + 2 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(2)=fitness1(2)+e_cont1(j)^2/block;
end
%% 粒子3
for j = (J + 20 * (t - TT) * block + 2 * block + 1) : (J + 20 * (t - TT) * block + 3 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(3)=fitness1(3)+e_cont1(j)^2/block;
end
%% 粒子4
for j = (J + 20 * (t - TT) * block + 3 * block + 1) : (J + 20 * (t - TT) * block + 4 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(4)=fitness1(4)+e_cont1(j)^2/block;
end
%% 粒子5
for j = (J + 20 * (t - TT) * block + 4 * block + 1) : (J + 20 * (t - TT) * block + 5 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(5)=fitness1(5)+e_cont1(j)^2/block;
end
%% 粒子6
for j = (J + 20 * (t - TT) * block + 5 * block + 1) : (J + 20 * (t - TT) * block + 6 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(6)=fitness1(6)+e_cont1(j)^2/block;
end
%% 粒子7
for j = (J + 20 * (t - TT) * block + 6 * block + 1) : (J + 20 * (t - TT) * block + 7 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(7)=fitness1(7)+e_cont1(j)^2/block;
end
%% 粒子8
for j = (J + 20 * (t - TT) * block + 7 * block + 1) : (J + 20 * (t - TT) * block + 8 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(8)=fitness1(8)+e_cont1(j)^2/block;
end
%% 粒子9
for j = (J + 20 * (t - TT) * block + 8 * block + 1) : (J + 20 * (t - TT) * block + 9 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(9)=fitness1(9)+e_cont1(j)^2/block;
end
%% 粒子10
for j = (J + 20 * (t - TT) * block + 9 * block + 1) : (J + 20 * (t - TT) * block + 10 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(10)=fitness1(10)+e_cont1(j)^2/block;
end
%% 粒子11
for j = (J + 20 * (t - TT) * block + 10 * block + 1) : (J + 20 * (t - TT) * block + 11 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(11)=fitness1(11)+e_cont1(j)^2/block;
end
%% 粒子12
for j = (J + 20 * (t - TT) * block + 11 * block + 1) : (J + 20 * (t - TT) * block + 12 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(12)=fitness1(12)+e_cont1(j)^2/block;
end
%% 粒子13
for j = (J + 20 * (t - TT) * block + 12 * block + 1) : (J + 20 * (t - TT) * block + 13 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(13)=fitness1(13)+e_cont1(j)^2/block;
end
%% 粒子14
for j = (J + 20 * (t - TT) * block + 13 * block + 1) : (J + 20 * (t - TT) * block + 14 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(14)=fitness1(14)+e_cont1(j)^2/block;
end
%% 粒子15
for j = (J + 20 * (t - TT) * block + 14 * block + 1) : (J + 20 * (t - TT) * block + 15 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(15)=fitness1(15)+e_cont1(j)^2/block;
end
%% 粒子16
for j = (J + 20 * (t - TT) * block + 15 * block + 1) : (J + 20 * (t - TT) * block + 16 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(16)=fitness1(16)+e_cont1(j)^2/block;
end
%% 粒子17
for j = (J + 20 * (t - TT) * block + 16 * block + 1) : (J + 20 * (t - TT) * block + 17 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(17)=fitness1(17)+e_cont1(j)^2/block;
end
%% 粒子18
for j = (J + 20 * (t - TT) * block + 17 * block + 1) : (J + 20 * (t - TT) * block + 18 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(18)=fitness1(18)+e_cont1(j)^2/block;
end
%% 粒子19
for j = (J + 20 * (t - TT) * block + 18 * block + 1) : (J + 20 * (t - TT) * block + 19 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(19)=fitness1(19)+e_cont1(j)^2/block;
end
%% 粒子20
for j = (J + 20 * (t - TT) * block + 19 * block + 1) : (J + 20 * (t - TT) * block + 20 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(20)=fitness1(20)+e_cont1(j)^2/block;
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
if fitness1(i)<Alpha_score 
    Alpha_score=fitness1(i); 
    Alpha_pos=X(i,:);
end
% 更新beta
if fitness1(i)>Alpha_score && fitness1(i)<Beta_score 
    Beta_score=fitness1(i); 
    Beta_pos=X(i,:);
end
% 更新delta
if fitness1(i)>Alpha_score && fitness1(i)>Beta_score && fitness1(i)<Delta_score 
    Delta_score=fitness1(i); 
    Delta_pos=X(i,:);
end
end
IterCurve1(t)=IterCurve1(t)+Alpha_score;
% IterCurve1(t)=IterCurve1(t) + mean(fitness1);
%% 更新
% 计算种群的平均位置向量
X_avg = mean(X, 1); % 逐列求平均，得到 1×d 的向量
% 计算每个个体到种群平均位置的欧氏距离
distances = zeros(pop, 1); % 初始化距离存储
for i = 1:pop
    distances(i) = norm(X(i, :) - X_avg); % 欧氏距离
end
% 计算种群多样性（平均距离）
D1(t) = mean(distances);

if D1(t) > eta 
    a1(t) = 2 * (1 - (t - TT + 1) / maxIter);
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
%%
%%
X(i,:)=boundary_constraints(X(i,:),ub,lb,dim);
end
end

%%
%% GWO算法
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 初始化    
ff = 0;
beta = 0.1; % 指数衰减速率
ANCCx2=zeros(1,dim);
ActualSPx2=zeros(1,Tap_Sw);  
e_cont2=zeros(1,DataLong);
D2 = zeros(1, maxIter);
%% 迭代开始
for t = 1: maxIter
fitness2=zeros(1,pop); 
% Sw = Sw1;
%% 粒子1 
for j = (20 * (t - 1) * block + 1) : (20 * (t - 1) * block + block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(1)=fitness2(1)+e_cont2(j)^2/block;
end
%% 粒子2
for j = (20 * (t - 1) * block + block + 1) : (20 * (t - 1) * block + 2 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(2)=fitness2(2)+e_cont2(j)^2/block;
end
%% 粒子3
for j = (20 * (t - 1) * block + 2 * block + 1) : (20 * (t - 1) * block + 3 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(3)=fitness2(3)+e_cont2(j)^2/block;
end
%% 粒子4
for j = (20 * (t - 1) * block + 3 * block + 1) : (20 * (t - 1) * block + 4 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(4)=fitness2(4)+e_cont2(j)^2/block;
end
%% 粒子5
for j = (20 * (t - 1) * block + 4 * block + 1) : (20 * (t - 1) * block + 5 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(5)=fitness2(5)+e_cont2(j)^2/block;
end
%% 粒子6
for j = (20 * (t - 1) * block + 5 * block + 1) : (20 * (t - 1) * block + 6 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(6)=fitness2(6)+e_cont2(j)^2/block;
end
%% 粒子7
for j = (20 * (t - 1) * block + 6 * block + 1) : (20 * (t - 1) * block + 7 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(7)=fitness2(7)+e_cont2(j)^2/block;
end
%% 粒子8
for j = (20 * (t - 1) * block + 7 * block + 1) : (20 * (t - 1) * block + 8 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(8)=fitness2(8)+e_cont2(j)^2/block;
end
%% 粒子9
for j = (20 * (t - 1) * block + 8 * block + 1) : (20 * (t - 1) * block + 9 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(9)=fitness2(9)+e_cont2(j)^2/block;
end
%% 粒子10
for j = (20 * (t - 1) * block + 9 * block + 1) : (20 * (t - 1) * block + 10 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(10)=fitness2(10)+e_cont2(j)^2/block;
end
%% 粒子11
for j = (20 * (t - 1) * block + 10 * block + 1) : (20 * (t - 1) * block + 11 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(11)=fitness2(11)+e_cont2(j)^2/block;
end
%% 粒子12
for j = (20 * (t - 1) * block + 11 * block + 1) : (20 * (t - 1) * block + 12 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(12)=fitness2(12)+e_cont2(j)^2/block;
end
%% 粒子13
for j = (20 * (t - 1) * block + 12 * block + 1) : (20 * (t - 1) * block + 13 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(13)=fitness2(13)+e_cont2(j)^2/block;
end
%% 粒子14
for j = (20 * (t - 1) * block + 13 * block + 1) : (20 * (t - 1) * block + 14 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(14)=fitness2(14)+e_cont2(j)^2/block;
end
%% 粒子15
for j = (20 * (t - 1) * block + 14 * block + 1) : (20 * (t - 1) * block + 15 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(15)=fitness2(15)+e_cont2(j)^2/block;
end
%% 粒子16
for j = (20 * (t - 1) * block + 15 * block + 1) : (20 * (t - 1) * block + 16 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(16)=fitness2(16)+e_cont2(j)^2/block;
end
%% 粒子17
for j = (20 * (t - 1) * block + 16 * block + 1) : (20 * (t - 1) * block + 17 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(17)=fitness2(17)+e_cont2(j)^2/block;
end
%% 粒子18
for j = (20 * (t - 1) * block + 17 * block + 1) : (20 * (t - 1) * block + 18 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(18)=fitness2(18)+e_cont2(j)^2/block;
end
%% 粒子19
for j = (20 * (t - 1) * block + 18 * block + 1) : (20 * (t - 1) * block + 19 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(19)=fitness2(19)+e_cont2(j)^2/block;
end
%% 粒子20
for j = (20 * (t - 1) * block + 19 * block + 1) : (20 * (t - 1) * block + 20 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(20)=fitness2(20)+e_cont2(j)^2/block;
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
if fitness2(i)<Alpha_score 
    Alpha_score=fitness2(i); 
    Alpha_pos=X(i,:);
end
% 更新beta
if fitness2(i)>Alpha_score && fitness2(i)<Beta_score 
    Beta_score=fitness2(i); 
    Beta_pos=X(i,:);
end
% 更新delta
if fitness2(i)>Alpha_score && fitness2(i)>Beta_score && fitness2(i)<Delta_score 
    Delta_score=fitness2(i); 
    Delta_pos=X(i,:);
end
end
max_score = max(fitness2);
IterCurve2(t) = IterCurve2(t) + Alpha_score;

MM = 0;
for i = 1:pop
    MM = MM + (fitness2(i) - mean(fitness2))^2;
end
IterCurve22(t) = sqrt(MM / (pop - 1));
%% 更新
if t > 2 && abs(IterCurve22(t) - IterCurve22(t - 1)) > 1 && abs(IterCurve22(t - 1) - IterCurve22(t - 2)) < 0.1
    ff = 1;
    TT = t;
    X=initialize_population(pop,ub,lb,dim);
else
    X_avg = mean(X, 1); % 逐列求平均，得到 1×d 的向量
    distances = zeros(pop, 1); % 初始化距离存储
    for i = 1:pop
        distances(i) = norm(X(i, :) - X_avg); % 欧氏距离
    end
    D2(t) = mean(distances);

    if ff == 0
        if D2(t) > eta 
            a2(t) = 2 * (1 - t / maxIter);
            T = t;
            a_mid = a2(t);
        else
            a2(t) = a_mid * exp(-b * (t - T));
        end
    else
        if D2(t) > eta 
            a2(t) = 2 * (1 - (t - TT) / maxIter);
            T = t;
            a_mid = a2(t);
        else
            a2(t) = a_mid * exp(-b * (t - T));
        end
    end
    
    for i = 1:pop
    if isequal(X(i,:), Alpha_pos) || isequal(X(i,:), Beta_pos) || isequal(X(i,:), Delta_pos)
        continue
    end
    for j=1:dim     
    %%                       
    r1=rand(); % r1 is a random number in [0,1]
    r2=rand(); % r2 is a random number in [0,1]
            
    A1=2*a2(t)*r1-a2(t); % Equation (3.3)
    C1=2*r2; % Equation (3.4)
            
    D_alpha=abs(C1*Alpha_pos(j)-X(i,j)); % Equation (3.5)-part 1
    X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
    %%                       
    r1=rand();
    r2=rand();
            
    A2=2*a2(t)*r1-a2(t); % Equation (3.3)
    C2=2*r2; % Equation (3.4)
            
    D_beta=abs(C2*Beta_pos(j)-X(i,j)); % Equation (3.5)-part 2
    X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
    %%            
    r1=rand();
    r2=rand(); 
            
    A3=2*a2(t)*r1-a2(t); % Equation (3.3)
    C3=2*r2; % Equation (3.4)
            
    D_delta=abs(C3*Delta_pos(j)-X(i,j)); % Equation (3.5)-part 3
    X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
    %%
    wa = 1 / 3;
    wb = 1 / 3;
    wc = 1 / 3;
    X(i, j) = wa * X1 + wb * X2 + wc * X3;% Equation (3.7)
            
    end
    %%
    %%
    X(i,:)=boundary_constraints(X(i,:),ub,lb,dim);
    end
end

end
%%
%% GWO算法
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 初始化    
ANCCx3=zeros(1,dim);
ActualSPx3=zeros(1,Tap_Sw);  
e_cont3=zeros(1,DataLong);
D3 = zeros(1, maxIter);
%% 迭代开始
for t = 1: maxIter
fitness3=zeros(1,pop); 
%% 粒子1 
for j = (20 * (t - 1) * block + 1) : (20 * (t - 1) * block + block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(1)=fitness3(1)+e_cont3(j)^2/block;
end
%% 粒子2
for j = (20 * (t - 1) * block + block + 1) : (20 * (t - 1) * block + 2 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(2)=fitness3(2)+e_cont3(j)^2/block;
end
%% 粒子3
for j = (20 * (t - 1) * block + 2 * block + 1) : (20 * (t - 1) * block + 3 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(3)=fitness3(3)+e_cont3(j)^2/block;
end
%% 粒子4
for j = (20 * (t - 1) * block + 3 * block + 1) : (20 * (t - 1) * block + 4 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(4)=fitness3(4)+e_cont3(j)^2/block;
end
%% 粒子5
for j = (20 * (t - 1) * block + 4 * block + 1) : (20 * (t - 1) * block + 5 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(5)=fitness3(5)+e_cont3(j)^2/block;
end
%% 粒子6
for j = (20 * (t - 1) * block + 5 * block + 1) : (20 * (t - 1) * block + 6 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(6)=fitness3(6)+e_cont3(j)^2/block;
end
%% 粒子7
for j = (20 * (t - 1) * block + 6 * block + 1) : (20 * (t - 1) * block + 7 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(7)=fitness3(7)+e_cont3(j)^2/block;
end
%% 粒子8
for j = (20 * (t - 1) * block + 7 * block + 1) : (20 * (t - 1) * block + 8 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(8)=fitness3(8)+e_cont3(j)^2/block;
end
%% 粒子9
for j = (20 * (t - 1) * block + 8 * block + 1) : (20 * (t - 1) * block + 9 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(9)=fitness3(9)+e_cont3(j)^2/block;
end
%% 粒子10
for j = (20 * (t - 1) * block + 9 * block + 1) : (20 * (t - 1) * block + 10 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(10)=fitness3(10)+e_cont3(j)^2/block;
end
%% 粒子11
for j = (20 * (t - 1) * block + 10 * block + 1) : (20 * (t - 1) * block + 11 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(11)=fitness3(11)+e_cont3(j)^2/block;
end
%% 粒子12
for j = (20 * (t - 1) * block + 11 * block + 1) : (20 * (t - 1) * block + 12 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(12)=fitness3(12)+e_cont3(j)^2/block;
end
%% 粒子13
for j = (20 * (t - 1) * block + 12 * block + 1) : (20 * (t - 1) * block + 13 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(13)=fitness3(13)+e_cont3(j)^2/block;
end
%% 粒子14
for j = (20 * (t - 1) * block + 13 * block + 1) : (20 * (t - 1) * block + 14 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(14)=fitness3(14)+e_cont3(j)^2/block;
end
%% 粒子15
for j = (20 * (t - 1) * block + 14 * block + 1) : (20 * (t - 1) * block + 15 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(15)=fitness3(15)+e_cont3(j)^2/block;
end
%% 粒子16
for j = (20 * (t - 1) * block + 15 * block + 1) : (20 * (t - 1) * block + 16 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(16)=fitness3(16)+e_cont3(j)^2/block;
end
%% 粒子17
for j = (20 * (t - 1) * block + 16 * block + 1) : (20 * (t - 1) * block + 17 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(17)=fitness3(17)+e_cont3(j)^2/block;
end
%% 粒子18
for j = (20 * (t - 1) * block + 17 * block + 1) : (20 * (t - 1) * block + 18 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(18)=fitness3(18)+e_cont3(j)^2/block;
end
%% 粒子19
for j = (20 * (t - 1) * block + 18 * block + 1) : (20 * (t - 1) * block + 19 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(19)=fitness3(19)+e_cont3(j)^2/block;
end
%% 粒子20
for j = (20 * (t - 1) * block + 19 * block + 1) : (20 * (t - 1) * block + 20 * block)
    if j <= (DataLong - 100) / 2 + 99
        Sw = Sw1;
    else
        Sw = Sw2;
    end
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(20)=fitness3(20)+e_cont3(j)^2/block;
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
if fitness3(i)<Alpha_score 
    Alpha_score=fitness3(i); 
    Alpha_pos=X(i,:);
end
% 更新beta
if fitness3(i)>Alpha_score && fitness3(i)<Beta_score 
    Beta_score=fitness3(i); 
    Beta_pos=X(i,:);
end
% 更新delta
if fitness3(i)>Alpha_score && fitness3(i)>Beta_score && fitness3(i)<Delta_score 
    Delta_score=fitness3(i); 
    Delta_pos=X(i,:);
end
end
IterCurve3(t)=IterCurve3(t)+Alpha_score;
%% 更新
% 计算种群的平均位置向量
X_avg = mean(X, 1); % 逐列求平均，得到 1×d 的向量
% 计算每个个体到种群平均位置的欧氏距离
distances = zeros(pop, 1); % 初始化距离存储
for i = 1:pop
    distances(i) = norm(X(i, :) - X_avg); % 欧氏距离
end
% 计算种群多样性（平均距离）
D3(t) = mean(distances);

eta = 0.1;
b = 0.1;
if D3(t) > eta
    a3(t) = 2 * (1 - t / maxIter);
    T = t;
    a_mid = a3(t);
else
    a3(t) = a_mid * exp(-b * (t - T));
end

for i = 1:pop
if isequal(X(i,:), Alpha_pos) || isequal(X(i,:), Beta_pos) || isequal(X(i,:), Delta_pos)
    continue
end
for j=1:dim     
%%                       
r1=rand(); % r1 is a random number in [0,1]
r2=rand(); % r2 is a random number in [0,1]
            
A1=2*a3(t)*r1-a3(t); % Equation (3.3)
C1=2*r2; % Equation (3.4)
            
D_alpha=abs(C1*Alpha_pos(j)-X(i,j)); % Equation (3.5)-part 1
X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
%%                       
r1=rand();
r2=rand();
            
A2=2*a3(t)*r1-a3(t); % Equation (3.3)
C2=2*r2; % Equation (3.4)
            
D_beta=abs(C2*Beta_pos(j)-X(i,j)); % Equation (3.5)-part 2
X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
%%            
r1=rand();
r2=rand(); 
            
A3=2*a3(t)*r1-a3(t); % Equation (3.3)
C3=2*r2; % Equation (3.4)
            
D_delta=abs(C3*Delta_pos(j)-X(i,j)); % Equation (3.5)-part 3
X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
%%
wa = 1 / 3;
wb = 1 / 3;
wc = 1 / 3;
X(i, j) = wa * X1 + wb * X2 + wc * X3;% Equation (3.7)
            
end
%%
%%
X(i,:)=boundary_constraints(X(i,:),ub,lb,dim);
end
end
%% ANR
for i = 1:DataLong
    if i == 1
        Me1(i) = (1 - lambda) * abs(e_cont1(i));
        Me2(i) = (1 - lambda) * abs(e_cont2(i));
        Me3(i) = (1 - lambda) * abs(e_cont3(i));
        Md(i) = (1 - lambda) * abs(Yp(i));
    else
        Me1(i) = lambda * Me1(i - 1) + (1 - lambda) * abs(e_cont1(i));
        Me2(i) = lambda * Me2(i - 1) + (1 - lambda) * abs(e_cont2(i));
        Me3(i) = lambda * Me3(i - 1) + (1 - lambda) * abs(e_cont3(i));
        Md(i) = lambda * Md(i - 1) +(1 - lambda) * abs(Yp(i));
    end
    ANR1(epoch, i) = 20 * log10(Md(i) / Me1(i));
    ANR2(epoch, i) = 20 * log10(Md(i) / Me2(i));
    ANR3(epoch, i) = 20 * log10(Md(i) / Me3(i));
end
end
%% 平均
AANR1 = mean(ANR1, 1);
AANR2 = mean(ANR2, 1);
AANR3 = mean(ANR3, 1);
% 高斯滤波法
windowSize = 5000; % 窗口大小
weights = gausswin(windowSize);
AANR11 = conv(AANR1, weights, 'same') / sum(weights);
AANR22 = conv(AANR2, weights, 'same') / sum(weights);
AANR33 = conv(AANR3, weights, 'same') / sum(weights);
%% 画图
%% 降噪量曲线图
figure(1);
plot(AANR33,'LineWidth',1.5);
hold on;
plot(AANR22,'LineWidth',1.5);
hold on;
plot(AANR11,'LineWidth',1.5);
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
ylabel('{\fontname{Times New Roman}ANR/dB}');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
legend('{\fontname{Times New Roman}Without reinitialization strategy}','{\fontname{Times New Roman}With traditional reinitialization strategy}','{\fontname{Times New Roman}With proposed reinitialization strategy}',...%    
     'NumColumns',1,'FontSize', 12, 'Location', 'southeast');
grid on
set(gca, 'GridLineStyle', ':');  % 设置为虚线 
set(gca, 'GridAlpha', 1);  % 设置透明度

figure(2);
plot(IterCurve3/q,'LineWidth',1.5);
hold on;
plot(IterCurve2/q,'LineWidth',1.5);
hold on;
plot(IterCurve1/q,'LineWidth',1.5);
xlabel('{\fontname{Times New Roman}Iteration Number}');
ylabel('{\fontname{Times New Roman}Best Fitness Value}');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
legend('{\fontname{Times New Roman}Without reinitialization strategy}','{\fontname{Times New Roman}With traditional reinitialization strategy}','{\fontname{Times New Roman}With proposed reinitialization strategy}',...%    
     'NumColumns',1,'FontSize', 12, 'Location', 'southeast');
grid on
set(gca, 'GridLineStyle', ':');  % 设置为虚线 
set(gca, 'GridAlpha', 1);  % 设置透明度
