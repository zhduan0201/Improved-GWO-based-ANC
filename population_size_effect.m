clc;clear all;close all;
global DataLong; 
%% 理想数据
%% 理想路径
Pw = [0 0 0 0 0 0 0 0 0 0.8 0.6 -0.2 -0.5 -0.1 0.4 -0.05];    
Sw = [0 0 0 0 0 1 2.5 1.76 0.15 -0.4825 -0.18625 -0.005 -0.001875]; 
Tap_Sw = length(Sw);     
q = 50; 
%% 关键参数
block = 50;
dim = 20; % 控制滤波器阶数
pop1 = 10;  
pop2 = 20;
pop3 = 30;
pop4 = 40;
pop5 = 50;
%% 初级噪声
maxIter = 30; 
DataLong1 = pop1 * block * maxIter;
DataLong2 = pop2 * block * maxIter;
DataLong3 = pop3 * block * maxIter; 
DataLong4 = pop4 * block * maxIter; 
DataLong5 = pop5 * block * maxIter; 
DataLong = pop5 * block * maxIter;      
snr = 30;    
[X_noise, ~] = add_awgn(sine_wave_generator(200) + sine_wave_generator(400), snr);
%% 期望信号
Yp=filter(Pw,1,X_noise); 
%% 
IterCurve1=zeros(1, maxIter); 
IterCurve2=zeros(1, maxIter);
IterCurve3=zeros(1, maxIter);
IterCurve4=zeros(1, maxIter);
IterCurve5=zeros(1, maxIter);
%% 重复实验
for epoch = 1:q
%% 群智能优化算法参数设定
x_max = 0.1;
ub = x_max * ones(1, dim); % 位置范围
lb = -x_max * ones(1, dim);
%%
%% GWO算法 P = 10
%% 初始化种群位置
X=initialize_population(pop1,ub,lb,dim);
%% 初始化    
ANCCx1=zeros(1,dim);
ActualSPx1=zeros(1,Tap_Sw);  
e_cont1=zeros(1,DataLong1);
%% 迭代开始
for t = 1: maxIter
fitness1=zeros(1,pop1); 
%% 粒子1 
for j = (10 * (t - 1) * block + 1) : (10 * (t - 1) * block + block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(1)=fitness1(1)+e_cont1(j)^2/block;
end
%% 粒子2
for j = (10 * (t - 1) * block + block + 1) : (10 * (t - 1) * block + 2 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(2)=fitness1(2)+e_cont1(j)^2/block;
end
%% 粒子3
for j = (10 * (t - 1) * block + 2 * block + 1) : (10 * (t - 1) * block + 3 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(3)=fitness1(3)+e_cont1(j)^2/block;
end
%% 粒子4
for j = (10 * (t - 1) * block + 3 * block + 1) : (10 * (t - 1) * block + 4 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(4)=fitness1(4)+e_cont1(j)^2/block;
end
%% 粒子5
for j = (10 * (t - 1) * block + 4 * block + 1) : (10 * (t - 1) * block + 5 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(5)=fitness1(5)+e_cont1(j)^2/block;
end
%% 粒子6
for j = (10 * (t - 1) * block + 5 * block + 1) : (10 * (t - 1) * block + 6 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(6)=fitness1(6)+e_cont1(j)^2/block;
end
%% 粒子7
for j = (10 * (t - 1) * block + 6 * block + 1) : (10 * (t - 1) * block + 7 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(7)=fitness1(7)+e_cont1(j)^2/block;
end
%% 粒子8
for j = (10 * (t - 1) * block + 7 * block + 1) : (10 * (t - 1) * block + 8 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(8)=fitness1(8)+e_cont1(j)^2/block;
end
%% 粒子9
for j = (10 * (t - 1) * block + 8 * block + 1) : (10 * (t - 1) * block + 9 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(9)=fitness1(9)+e_cont1(j)^2/block;
end
%% 粒子10
for j = (10 * (t - 1) * block + 9 * block + 1) : (10 * (t - 1) * block + 10 * block)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(10)=fitness1(10)+e_cont1(j)^2/block;
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
for i=1:pop1
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
IterCurve1(t) = IterCurve1(t) + Alpha_score;
%% 更新
a = 2 * (1 - t / maxIter );

for i = 1:pop1
if isequal(X(i,:), Alpha_pos) || isequal(X(i,:), Beta_pos) || isequal(X(i,:), Delta_pos)
    continue
end
for j=1:dim     
%%                       
r1=rand(); % r1 is a random number in [0,1]
r2=rand(); % r2 is a random number in [0,1]
            
A1=2*a*r1-a; % Equation (3.3)
C1=2*r2; % Equation (3.4)
            
D_alpha=abs(C1*Alpha_pos(j)-X(i,j)); % Equation (3.5)-part 1
X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
%%                       
r1=rand();
r2=rand();
            
A2=2*a*r1-a; % Equation (3.3)
C2=2*r2; % Equation (3.4)
            
D_beta=abs(C2*Beta_pos(j)-X(i,j)); % Equation (3.5)-part 2
X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
%%            
r1=rand();
r2=rand(); 
            
A3=2*a*r1-a; % Equation (3.3)
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
%%
%% GWO算法 P = 20
%% 初始化种群位置
X=initialize_population(pop2,ub,lb,dim);
%% 初始化    
ANCCx2=zeros(1,dim);
ActualSPx2=zeros(1,Tap_Sw);  
e_cont2=zeros(1,DataLong2);
%% 迭代开始
for t = 1: maxIter
fitness2=zeros(1,pop2); 
%% 粒子1 
for j = (20 * (t - 1) * block + 1) : (20 * (t - 1) * block + block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(1)=fitness2(1)+e_cont2(j)^2/block;
end
%% 粒子2
for j = (20 * (t - 1) * block + block + 1) : (20 * (t - 1) * block + 2 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(2)=fitness2(2)+e_cont2(j)^2/block;
end
%% 粒子3
for j = (20 * (t - 1) * block + 2 * block + 1) : (20 * (t - 1) * block + 3 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(3)=fitness2(3)+e_cont2(j)^2/block;
end
%% 粒子4
for j = (20 * (t - 1) * block + 3 * block + 1) : (20 * (t - 1) * block + 4 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(4)=fitness2(4)+e_cont2(j)^2/block;
end
%% 粒子5
for j = (20 * (t - 1) * block + 4 * block + 1) : (20 * (t - 1) * block + 5 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(5)=fitness2(5)+e_cont2(j)^2/block;
end
%% 粒子6
for j = (20 * (t - 1) * block + 5 * block + 1) : (20 * (t - 1) * block + 6 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(6)=fitness2(6)+e_cont2(j)^2/block;
end
%% 粒子7
for j = (20 * (t - 1) * block + 6 * block + 1) : (20 * (t - 1) * block + 7 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(7)=fitness2(7)+e_cont2(j)^2/block;
end
%% 粒子8
for j = (20 * (t - 1) * block + 7 * block + 1) : (20 * (t - 1) * block + 8 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(8)=fitness2(8)+e_cont2(j)^2/block;
end
%% 粒子9
for j = (20 * (t - 1) * block + 8 * block + 1) : (20 * (t - 1) * block + 9 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(9)=fitness2(9)+e_cont2(j)^2/block;
end
%% 粒子10
for j = (20 * (t - 1) * block + 9 * block + 1) : (20 * (t - 1) * block + 10 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(10)=fitness2(10)+e_cont2(j)^2/block;
end
%% 粒子11
for j = (20 * (t - 1) * block + 10 * block + 1) : (20 * (t - 1) * block + 11 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(11)=fitness2(11)+e_cont2(j)^2/block;
end
%% 粒子12
for j = (20 * (t - 1) * block + 11 * block + 1) : (20 * (t - 1) * block + 12 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(12)=fitness2(12)+e_cont2(j)^2/block;
end
%% 粒子13
for j = (20 * (t - 1) * block + 12 * block + 1) : (20 * (t - 1) * block + 13 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(13)=fitness2(13)+e_cont2(j)^2/block;
end
%% 粒子14
for j = (20 * (t - 1) * block + 13 * block + 1) : (20 * (t - 1) * block + 14 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(14)=fitness2(14)+e_cont2(j)^2/block;
end
%% 粒子15
for j = (20 * (t - 1) * block + 14 * block + 1) : (20 * (t - 1) * block + 15 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(15)=fitness2(15)+e_cont2(j)^2/block;
end
%% 粒子16
for j = (20 * (t - 1) * block + 15 * block + 1) : (20 * (t - 1) * block + 16 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(16)=fitness2(16)+e_cont2(j)^2/block;
end
%% 粒子17
for j = (20 * (t - 1) * block + 16 * block + 1) : (20 * (t - 1) * block + 17 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(17)=fitness2(17)+e_cont2(j)^2/block;
end
%% 粒子18
for j = (20 * (t - 1) * block + 17 * block + 1) : (20 * (t - 1) * block + 18 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(18)=fitness2(18)+e_cont2(j)^2/block;
end
%% 粒子19
for j = (20 * (t - 1) * block + 18 * block + 1) : (20 * (t - 1) * block + 19 * block)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(19)=fitness2(19)+e_cont2(j)^2/block;
end
%% 粒子20
for j = (20 * (t - 1) * block + 19 * block + 1) : (20 * (t - 1) * block + 20 * block)
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
for i=1:pop2
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
IterCurve2(t) = IterCurve2(t) + Alpha_score;
%% 更新
a = 2 * (1 - t / maxIter );

for i = 1:pop2
if isequal(X(i,:), Alpha_pos) || isequal(X(i,:), Beta_pos) || isequal(X(i,:), Delta_pos)
    continue
end
for j=1:dim     
%%                       
r1=rand(); % r1 is a random number in [0,1]
r2=rand(); % r2 is a random number in [0,1]
            
A1=2*a*r1-a; % Equation (3.3)
C1=2*r2; % Equation (3.4)
            
D_alpha=abs(C1*Alpha_pos(j)-X(i,j)); % Equation (3.5)-part 1
X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
%%                       
r1=rand();
r2=rand();
            
A2=2*a*r1-a; % Equation (3.3)
C2=2*r2; % Equation (3.4)
            
D_beta=abs(C2*Beta_pos(j)-X(i,j)); % Equation (3.5)-part 2
X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
%%            
r1=rand();
r2=rand(); 
            
A3=2*a*r1-a; % Equation (3.3)
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
%%
%% GWO算法 P = 30
%% 初始化种群位置
X=initialize_population(pop3,ub,lb,dim);
%% 初始化    
ANCCx3=zeros(1,dim);
ActualSPx3=zeros(1,Tap_Sw);  
e_cont3=zeros(1,DataLong3);
%% 迭代开始
for t = 1: maxIter
fitness3=zeros(1,pop3); 
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
%% 粒子21
for j = (30 * (t - 1) * block + 20 * block + 1) : (30 * (t - 1) * block + 21 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(21,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(21)=fitness3(21)+e_cont3(j)^2/block;
end
%% 粒子22
for j = (30 * (t - 1) * block + 21 * block + 1) : (30 * (t - 1) * block + 22 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(22,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(22)=fitness3(22)+e_cont3(j)^2/block;
end
%% 粒子23
for j = (30 * (t - 1) * block + 22 * block + 1) : (30 * (t - 1) * block + 23 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(23,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(23)=fitness3(23)+e_cont3(j)^2/block;
end
%% 粒子24
for j = (30 * (t - 1) * block + 23 * block + 1) : (30 * (t - 1) * block + 24 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(24,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(24)=fitness3(24)+e_cont3(j)^2/block;
end
%% 粒子25
for j = (30 * (t - 1) * block + 24 * block + 1) : (30 * (t - 1) * block + 25 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(25,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(25)=fitness3(25)+e_cont3(j)^2/block;
end
%% 粒子26
for j = (30 * (t - 1) * block + 25 * block + 1) : (30 * (t - 1) * block + 26 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(26,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(26)=fitness3(26)+e_cont3(j)^2/block;
end
%% 粒子27
for j = (30 * (t - 1) * block + 26 * block + 1) : (30 * (t - 1) * block + 27 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(27,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(27)=fitness3(27)+e_cont3(j)^2/block;
end
%% 粒子28
for j = (30 * (t - 1) * block + 27 * block + 1) : (30 * (t - 1) * block + 28 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(28,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(28)=fitness3(28)+e_cont3(j)^2/block;
end
%% 粒子29
for j = (30 * (t - 1) * block + 28 * block + 1) : (30 * (t - 1) * block + 29 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(29,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(29)=fitness3(29)+e_cont3(j)^2/block;
end
%% 粒子30
for j = (30 * (t - 1) * block + 29 * block + 1) : (30 * (t - 1) * block + 30 * block)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(30,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(30)=fitness3(30)+e_cont3(j)^2/block;
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
for i=1:pop3
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
IterCurve3(t) = IterCurve3(t) + Alpha_score;
%% 更新
a = 2 * (1 - t / maxIter );

for i = 1:pop3
if isequal(X(i,:), Alpha_pos) || isequal(X(i,:), Beta_pos) || isequal(X(i,:), Delta_pos)
    continue
end
for j=1:dim     
%%                       
r1=rand(); % r1 is a random number in [0,1]
r2=rand(); % r2 is a random number in [0,1]
            
A1=2*a*r1-a; % Equation (3.3)
C1=2*r2; % Equation (3.4)
            
D_alpha=abs(C1*Alpha_pos(j)-X(i,j)); % Equation (3.5)-part 1
X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
%%                       
r1=rand();
r2=rand();
            
A2=2*a*r1-a; % Equation (3.3)
C2=2*r2; % Equation (3.4)
            
D_beta=abs(C2*Beta_pos(j)-X(i,j)); % Equation (3.5)-part 2
X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
%%            
r1=rand();
r2=rand(); 
            
A3=2*a*r1-a; % Equation (3.3)
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
%%
%% GWO算法 P = 40
%% 初始化种群位置
X=initialize_population(pop4,ub,lb,dim);
%% 初始化    
ANCCx4=zeros(1,dim);
ActualSPx4=zeros(1,Tap_Sw);  
e_cont4=zeros(1,DataLong4);
%% 迭代开始
for t = 1: maxIter
fitness4=zeros(1,pop4); 
%% 粒子1 
for j = (40 * (t - 1) * block + 1) : (40 * (t - 1) * block + block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(1)=fitness4(1)+e_cont4(j)^2/block;
end
%% 粒子2
for j = (40 * (t - 1) * block + block + 1) : (40 * (t - 1) * block + 2 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(2)=fitness4(2)+e_cont4(j)^2/block;
end
%% 粒子3
for j = (40 * (t - 1) * block + 2 * block + 1) : (40 * (t - 1) * block + 3 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(3)=fitness4(3)+e_cont4(j)^2/block;
end
%% 粒子4
for j = (40 * (t - 1) * block + 3 * block + 1) : (40 * (t - 1) * block + 4 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(4)=fitness4(4)+e_cont4(j)^2/block;
end
%% 粒子5
for j = (40 * (t - 1) * block + 4 * block + 1) : (40 * (t - 1) * block + 5 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(5)=fitness4(5)+e_cont4(j)^2/block;
end
%% 粒子6
for j = (40 * (t - 1) * block + 5 * block + 1) : (40 * (t - 1) * block + 6 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(6)=fitness4(6)+e_cont4(j)^2/block;
end
%% 粒子7
for j = (40 * (t - 1) * block + 6 * block + 1) : (40 * (t - 1) * block + 7 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(7)=fitness4(7)+e_cont4(j)^2/block;
end
%% 粒子8
for j = (40 * (t - 1) * block + 7 * block + 1) : (40 * (t - 1) * block + 8 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(8)=fitness4(8)+e_cont4(j)^2/block;
end
%% 粒子9
for j = (40 * (t - 1) * block + 8 * block + 1) : (40 * (t - 1) * block + 9 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(9)=fitness4(9)+e_cont4(j)^2/block;
end
%% 粒子10
for j = (40 * (t - 1) * block + 9 * block + 1) : (40 * (t - 1) * block + 10 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(10)=fitness4(10)+e_cont4(j)^2/block;
end
%% 粒子11
for j = (40 * (t - 1) * block + 10 * block + 1) : (40 * (t - 1) * block + 11 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(11)=fitness4(11)+e_cont4(j)^2/block;
end
%% 粒子12
for j = (40 * (t - 1) * block + 11 * block + 1) : (40 * (t - 1) * block + 12 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(12)=fitness4(12)+e_cont4(j)^2/block;
end
%% 粒子13
for j = (40 * (t - 1) * block + 12 * block + 1) : (40 * (t - 1) * block + 13 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(13)=fitness4(13)+e_cont4(j)^2/block;
end
%% 粒子14
for j = (40 * (t - 1) * block + 13 * block + 1) : (40 * (t - 1) * block + 14 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(14)=fitness4(14)+e_cont4(j)^2/block;
end
%% 粒子15
for j = (40 * (t - 1) * block + 14 * block + 1) : (40 * (t - 1) * block + 15 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(15)=fitness4(15)+e_cont4(j)^2/block;
end
%% 粒子16
for j = (40 * (t - 1) * block + 15 * block + 1) : (40 * (t - 1) * block + 16 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(16)=fitness4(16)+e_cont4(j)^2/block;
end
%% 粒子17
for j = (40 * (t - 1) * block + 16 * block + 1) : (40 * (t - 1) * block + 17 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(17)=fitness4(17)+e_cont4(j)^2/block;
end
%% 粒子18
for j = (40 * (t - 1) * block + 17 * block + 1) : (40 * (t - 1) * block + 18 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(18)=fitness4(18)+e_cont4(j)^2/block;
end
%% 粒子19
for j = (40 * (t - 1) * block + 18 * block + 1) : (40 * (t - 1) * block + 19 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(19)=fitness4(19)+e_cont4(j)^2/block;
end
%% 粒子20
for j = (40 * (t - 1) * block + 19 * block + 1) : (40 * (t - 1) * block + 20 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(20)=fitness4(20)+e_cont4(j)^2/block;
end
%% 粒子21
for j = (40 * (t - 1) * block + 20 * block + 1) : (40 * (t - 1) * block + 21 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(21,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(21)=fitness4(21)+e_cont4(j)^2/block;
end
%% 粒子22
for j = (40 * (t - 1) * block + 21 * block + 1) : (40 * (t - 1) * block + 22 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(22,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(22)=fitness4(22)+e_cont4(j)^2/block;
end
%% 粒子23
for j = (40 * (t - 1) * block + 22 * block + 1) : (40 * (t - 1) * block + 23 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(23,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(23)=fitness4(23)+e_cont4(j)^2/block;
end
%% 粒子24
for j = (40 * (t - 1) * block + 23 * block + 1) : (40 * (t - 1) * block + 24 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(24,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(24)=fitness4(24)+e_cont4(j)^2/block;
end
%% 粒子25
for j = (40 * (t - 1) * block + 24 * block + 1) : (40 * (t - 1) * block + 25 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(25,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(25)=fitness4(25)+e_cont4(j)^2/block;
end
%% 粒子26
for j = (40 * (t - 1) * block + 25 * block + 1) : (40 * (t - 1) * block + 26 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(26,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(26)=fitness4(26)+e_cont4(j)^2/block;
end
%% 粒子27
for j = (40 * (t - 1) * block + 26 * block + 1) : (40 * (t - 1) * block + 27 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(27,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(27)=fitness4(27)+e_cont4(j)^2/block;
end
%% 粒子28
for j = (40 * (t - 1) * block + 27 * block + 1) : (40 * (t - 1) * block + 28 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(28,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(28)=fitness4(28)+e_cont4(j)^2/block;
end
%% 粒子29
for j = (40 * (t - 1) * block + 28 * block + 1) : (40 * (t - 1) * block + 29 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(29,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(29)=fitness4(29)+e_cont4(j)^2/block;
end
%% 粒子30
for j = (40 * (t - 1) * block + 29 * block + 1) : (40 * (t - 1) * block + 30 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(30,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(30)=fitness4(30)+e_cont4(j)^2/block;
end
%% 粒子31
for j = (40 * (t - 1) * block + 30 * block + 1) : (40 * (t - 1) * block + 31 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(31,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(31)=fitness4(31)+e_cont4(j)^2/block;
end
%% 粒子32
for j = (40 * (t - 1) * block + 31 * block + 1) : (40 * (t - 1) * block + 32 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(32,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(32)=fitness4(32)+e_cont4(j)^2/block;
end
%% 粒子33
for j = (40 * (t - 1) * block + 32 * block + 1) : (40 * (t - 1) * block + 33 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(33,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(33)=fitness4(33)+e_cont4(j)^2/block;
end
%% 粒子34
for j = (40 * (t - 1) * block + 33 * block + 1) : (40 * (t - 1) * block + 34 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(34,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(34)=fitness4(34)+e_cont4(j)^2/block;
end
%% 粒子35
for j = (40 * (t - 1) * block + 34 * block + 1) : (40 * (t - 1) * block + 35 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(35,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(35)=fitness4(35)+e_cont4(j)^2/block;
end
%% 粒子36
for j = (40 * (t - 1) * block + 35 * block + 1) : (40 * (t - 1) * block + 36 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(36,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(36)=fitness4(36)+e_cont4(j)^2/block;
end
%% 粒子37
for j = (40 * (t - 1) * block + 36 * block + 1) : (40 * (t - 1) * block + 37 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(37,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(37)=fitness4(37)+e_cont4(j)^2/block;
end
%% 粒子38
for j = (40 * (t - 1) * block + 37 * block + 1) : (40 * (t - 1) * block + 38 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(38,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(38)=fitness4(38)+e_cont4(j)^2/block;
end
%% 粒子39
for j = (40 * (t - 1) * block + 38 * block + 1) : (40 * (t - 1) * block + 39 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(39,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(39)=fitness4(39)+e_cont4(j)^2/block;
end
%% 粒子40
for j = (40 * (t - 1) * block + 39 * block + 1) : (40 * (t - 1) * block + 40 * block)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(40,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(40)=fitness4(40)+e_cont4(j)^2/block;
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
for i=1:pop4
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
IterCurve4(t) = IterCurve4(t) + Alpha_score;
%% 更新
a = 2 * (1 - t / maxIter );

for i = 1:pop4
if isequal(X(i,:), Alpha_pos) || isequal(X(i,:), Beta_pos) || isequal(X(i,:), Delta_pos)
    continue
end
for j=1:dim     
%%                       
r1=rand(); % r1 is a random number in [0,1]
r2=rand(); % r2 is a random number in [0,1]
            
A1=2*a*r1-a; % Equation (3.3)
C1=2*r2; % Equation (3.4)
            
D_alpha=abs(C1*Alpha_pos(j)-X(i,j)); % Equation (3.5)-part 1
X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
%%                       
r1=rand();
r2=rand();
            
A2=2*a*r1-a; % Equation (3.3)
C2=2*r2; % Equation (3.4)
            
D_beta=abs(C2*Beta_pos(j)-X(i,j)); % Equation (3.5)-part 2
X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
%%            
r1=rand();
r2=rand(); 
            
A3=2*a*r1-a; % Equation (3.3)
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
%%
%% GWO算法 P = 50
%% 初始化种群位置
X=initialize_population(pop5,ub,lb,dim);
%% 初始化    
ANCCx5=zeros(1,dim);
ActualSPx5=zeros(1,Tap_Sw);  
e_cont5=zeros(1,DataLong5);
%% 迭代开始
for t = 1: maxIter
fitness5=zeros(1,pop5); 
%% 粒子1 
for j = (50 * (t - 1) * block + 1) : (50 * (t - 1) * block + block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(1)=fitness5(1)+e_cont5(j)^2/block;
end
%% 粒子2
for j = (50 * (t - 1) * block + block + 1) : (50 * (t - 1) * block + 2 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(2)=fitness5(2)+e_cont5(j)^2/block;
end
%% 粒子3
for j = (50 * (t - 1) * block + 2 * block + 1) : (50 * (t - 1) * block + 3 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(3)=fitness5(3)+e_cont5(j)^2/block;
end
%% 粒子4
for j = (50 * (t - 1) * block + 3 * block + 1) : (50 * (t - 1) * block + 4 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(4)=fitness5(4)+e_cont5(j)^2/block;
end
%% 粒子5
for j = (50 * (t - 1) * block + 4 * block + 1) : (50 * (t - 1) * block + 5 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(5)=fitness5(5)+e_cont5(j)^2/block;
end
%% 粒子6
for j = (50 * (t - 1) * block + 5 * block + 1) : (50 * (t - 1) * block + 6 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(6)=fitness5(6)+e_cont5(j)^2/block;
end
%% 粒子7
for j = (50 * (t - 1) * block + 6 * block + 1) : (50 * (t - 1) * block + 7 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(7)=fitness5(7)+e_cont5(j)^2/block;
end
%% 粒子8
for j = (50 * (t - 1) * block + 7 * block + 1) : (50 * (t - 1) * block + 8 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(8)=fitness5(8)+e_cont5(j)^2/block;
end
%% 粒子9
for j = (50 * (t - 1) * block + 8 * block + 1) : (50 * (t - 1) * block + 9 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(9)=fitness5(9)+e_cont5(j)^2/block;
end
%% 粒子10
for j = (50 * (t - 1) * block + 9 * block + 1) : (50 * (t - 1) * block + 10 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(10)=fitness5(10)+e_cont5(j)^2/block;
end
%% 粒子11
for j = (50 * (t - 1) * block + 10 * block + 1) : (50 * (t - 1) * block + 11 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(11)=fitness5(11)+e_cont5(j)^2/block;
end
%% 粒子12
for j = (50 * (t - 1) * block + 11 * block + 1) : (50 * (t - 1) * block + 12 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(12)=fitness5(12)+e_cont5(j)^2/block;
end
%% 粒子13
for j = (50 * (t - 1) * block + 12 * block + 1) : (50 * (t - 1) * block + 13 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(13)=fitness5(13)+e_cont5(j)^2/block;
end
%% 粒子14
for j = (50 * (t - 1) * block + 13 * block + 1) : (50 * (t - 1) * block + 14 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(14)=fitness5(14)+e_cont5(j)^2/block;
end
%% 粒子15
for j = (50 * (t - 1) * block + 14 * block + 1) : (50 * (t - 1) * block + 15 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(15)=fitness5(15)+e_cont5(j)^2/block;
end
%% 粒子16
for j = (50 * (t - 1) * block + 15 * block + 1) : (50 * (t - 1) * block + 16 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(16)=fitness5(16)+e_cont5(j)^2/block;
end
%% 粒子17
for j = (50 * (t - 1) * block + 16 * block + 1) : (50 * (t - 1) * block + 17 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(17)=fitness5(17)+e_cont5(j)^2/block;
end
%% 粒子18
for j = (50 * (t - 1) * block + 17 * block + 1) : (50 * (t - 1) * block + 18 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(18)=fitness5(18)+e_cont5(j)^2/block;
end
%% 粒子19
for j = (50 * (t - 1) * block + 18 * block + 1) : (50 * (t - 1) * block + 19 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(19)=fitness5(19)+e_cont5(j)^2/block;
end
%% 粒子20
for j = (50 * (t - 1) * block + 19 * block + 1) : (50 * (t - 1) * block + 20 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(20)=fitness5(20)+e_cont5(j)^2/block;
end
%% 粒子21
for j = (50 * (t - 1) * block + 20 * block + 1) : (50 * (t - 1) * block + 21 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(21,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(21)=fitness5(21)+e_cont5(j)^2/block;
end
%% 粒子22
for j = (50 * (t - 1) * block + 21 * block + 1) : (50 * (t - 1) * block + 22 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(22,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(22)=fitness5(22)+e_cont5(j)^2/block;
end
%% 粒子23
for j = (50 * (t - 1) * block + 22 * block + 1) : (50 * (t - 1) * block + 23 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(23,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(23)=fitness5(23)+e_cont5(j)^2/block;
end
%% 粒子24
for j = (50 * (t - 1) * block + 23 * block + 1) : (50 * (t - 1) * block + 24 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(24,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(24)=fitness5(24)+e_cont5(j)^2/block;
end
%% 粒子25
for j = (50 * (t - 1) * block + 24 * block + 1) : (50 * (t - 1) * block + 25 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(25,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(25)=fitness5(25)+e_cont5(j)^2/block;
end
%% 粒子26
for j = (50 * (t - 1) * block + 25 * block + 1) : (50 * (t - 1) * block + 26 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(26,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(26)=fitness5(26)+e_cont5(j)^2/block;
end
%% 粒子27
for j = (50 * (t - 1) * block + 26 * block + 1) : (50 * (t - 1) * block + 27 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(27,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(27)=fitness5(27)+e_cont5(j)^2/block;
end
%% 粒子28
for j = (50 * (t - 1) * block + 27 * block + 1) : (50 * (t - 1) * block + 28 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(28,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(28)=fitness5(28)+e_cont5(j)^2/block;
end
%% 粒子29
for j = (50 * (t - 1) * block + 28 * block + 1) : (50 * (t - 1) * block + 29 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(29,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(29)=fitness5(29)+e_cont5(j)^2/block;
end
%% 粒子30
for j = (50 * (t - 1) * block + 29 * block + 1) : (50 * (t - 1) * block + 30 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(30,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(30)=fitness5(30)+e_cont5(j)^2/block;
end
%% 粒子31
for j = (50 * (t - 1) * block + 30 * block + 1) : (50 * (t - 1) * block + 31 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(31,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(31)=fitness5(31)+e_cont5(j)^2/block;
end
%% 粒子32
for j = (50 * (t - 1) * block + 31 * block + 1) : (50 * (t - 1) * block + 32 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(32,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(32)=fitness5(32)+e_cont5(j)^2/block;
end
%% 粒子33
for j = (50 * (t - 1) * block + 32 * block + 1) : (50 * (t - 1) * block + 33 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(33,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(33)=fitness5(33)+e_cont5(j)^2/block;
end
%% 粒子34
for j = (50 * (t - 1) * block + 33 * block + 1) : (50 * (t - 1) * block + 34 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(34,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(34)=fitness5(34)+e_cont5(j)^2/block;
end
%% 粒子35
for j = (50 * (t - 1) * block + 34 * block + 1) : (50 * (t - 1) * block + 35 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(35,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(35)=fitness5(35)+e_cont5(j)^2/block;
end
%% 粒子36
for j = (50 * (t - 1) * block + 35 * block + 1) : (50 * (t - 1) * block + 36 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(36,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(36)=fitness5(36)+e_cont5(j)^2/block;
end
%% 粒子37
for j = (50 * (t - 1) * block + 36 * block + 1) : (50 * (t - 1) * block + 37 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(37,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(37)=fitness5(37)+e_cont5(j)^2/block;
end
%% 粒子38
for j = (50 * (t - 1) * block + 37 * block + 1) : (50 * (t - 1) * block + 38 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(38,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(38)=fitness5(38)+e_cont5(j)^2/block;
end
%% 粒子39
for j = (50 * (t - 1) * block + 38 * block + 1) : (50 * (t - 1) * block + 39 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(39,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(39)=fitness5(39)+e_cont5(j)^2/block;
end
%% 粒子40
for j = (50 * (t - 1) * block + 39 * block + 1) : (50 * (t - 1) * block + 40 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(40,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(40)=fitness5(40)+e_cont5(j)^2/block;
end
%% 粒子41
for j = (50 * (t - 1) * block + 40 * block + 1) : (50 * (t - 1) * block + 41 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(41,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(41)=fitness5(41)+e_cont5(j)^2/block;
end
%% 粒子42
for j = (50 * (t - 1) * block + 41 * block + 1) : (50 * (t - 1) * block + 42 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(42,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(42)=fitness5(42)+e_cont5(j)^2/block;
end
%% 粒子43
for j = (50 * (t - 1) * block + 42 * block + 1) : (50 * (t - 1) * block + 43 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(43,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(43)=fitness5(43)+e_cont5(j)^2/block;
end
%% 粒子44
for j = (50 * (t - 1) * block + 43 * block + 1) : (50 * (t - 1) * block + 44 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(44,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(44)=fitness5(44)+e_cont5(j)^2/block;
end
%% 粒子45
for j = (50 * (t - 1) * block + 44 * block + 1) : (50 * (t - 1) * block + 45 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(45,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(45)=fitness5(45)+e_cont5(j)^2/block;
end
%% 粒子46
for j = (50 * (t - 1) * block + 45 * block + 1) : (50 * (t - 1) * block + 46 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(46,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(46)=fitness5(46)+e_cont5(j)^2/block;
end
%% 粒子47
for j = (50 * (t - 1) * block + 46 * block + 1) : (50 * (t - 1) * block + 47 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(47,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(47)=fitness5(47)+e_cont5(j)^2/block;
end
%% 粒子48
for j = (50 * (t - 1) * block + 47 * block + 1) : (50 * (t - 1) * block + 48 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(48,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(48)=fitness5(48)+e_cont5(j)^2/block;
end
%% 粒子49
for j = (50 * (t - 1) * block + 48 * block + 1) : (50 * (t - 1) * block + 49 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(49,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(49)=fitness5(49)+e_cont5(j)^2/block;
end
%% 粒子50
for j = (50 * (t - 1) * block + 49 * block + 1) : (50 * (t - 1) * block + 50 * block)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(50,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(50)=fitness5(50)+e_cont5(j)^2/block;
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
for i=1:pop5
if fitness5(i)<Alpha_score 
    Alpha_score=fitness5(i); 
    Alpha_pos=X(i,:);
end
% 更新beta
if fitness5(i)>Alpha_score && fitness5(i)<Beta_score 
    Beta_score=fitness5(i); 
    Beta_pos=X(i,:);
end
% 更新delta
if fitness5(i)>Alpha_score && fitness5(i)>Beta_score && fitness5(i)<Delta_score 
    Delta_score=fitness5(i); 
    Delta_pos=X(i,:);
end
end
IterCurve5(t) = IterCurve5(t) + Alpha_score;
%% 更新
a = 2 * (1 - t / maxIter );

for i = 1:pop5
if isequal(X(i,:), Alpha_pos) || isequal(X(i,:), Beta_pos) || isequal(X(i,:), Delta_pos)
    continue
end
for j=1:dim     
%%                       
r1=rand(); % r1 is a random number in [0,1]
r2=rand(); % r2 is a random number in [0,1]
            
A1=2*a*r1-a; % Equation (3.3)
C1=2*r2; % Equation (3.4)
            
D_alpha=abs(C1*Alpha_pos(j)-X(i,j)); % Equation (3.5)-part 1
X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
%%                       
r1=rand();
r2=rand();
            
A2=2*a*r1-a; % Equation (3.3)
C2=2*r2; % Equation (3.4)
            
D_beta=abs(C2*Beta_pos(j)-X(i,j)); % Equation (3.5)-part 2
X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
%%            
r1=rand();
r2=rand(); 
            
A3=2*a*r1-a; % Equation (3.3)
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
end
%% 画图
%% 适应度函数
figure(1);
plot(IterCurve1/q,'LineWidth',1.5);
hold on;
plot(IterCurve2/q,'LineWidth',1.5);
hold on;
plot(IterCurve3/q,'LineWidth',1.5);
hold on;
plot(IterCurve4/q,'LineWidth',1.5);
hold on;
plot(IterCurve5/q,'LineWidth',1.5);
xlabel('{\fontname{Times New Roman}Iteration Number}');
ylabel('{\fontname{Times New Roman}Best Fitness Value}');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
legend('{\fontname{Times New Roman}P=10}','{\fontname{Times New Roman}P=20}','{\fontname{Times New Roman}P=30}','{\fontname{Times New Roman}P=40}','{\fontname{Times New Roman}P=50}',...%    
     'NumColumns',1,'FontSize', 12, 'Location', 'northeast');
grid on
set(gca, 'GridLineStyle', ':');  % 设置为虚线 
set(gca, 'GridAlpha', 1);  % 设置透明度