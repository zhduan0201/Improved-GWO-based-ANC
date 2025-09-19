clc;clear all;close all;
global DataLong; 
%% 理想数据
%% 理想路径
Pw = [0 0 0 0 0 0 0 0 0 0.8 0.6 -0.2 -0.5 -0.1 0.4 -0.05];    
Sw = [0 0 0 0 0 1 2.5 1.76 0.15 -0.4825 -0.18625 -0.005 -0.001875]; 
Tap_Sw = length(Sw);     
q = 50; 
%% 关键参数 
block1 = 10;
block2 = 20;   
block3 = 30;
block4 = 50;    
block5 = 100;      
dim = 20; % 控制滤波器阶数
pop = 20;
%% 初级噪声
maxIter = 30; 
DataLong1 = pop * block1 * maxIter;
DataLong2 = pop * block2 * maxIter;
DataLong3 = pop * block3 * maxIter; 
DataLong4 = pop * block4 * maxIter; 
DataLong5 = pop * block5 * maxIter; 
DataLong = pop * block5 * maxIter;      
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
%% GWO算法 lambda = 10
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 初始化    
ANCCx1=zeros(1,dim);
ActualSPx1=zeros(1,Tap_Sw);  
e_cont1=zeros(1,DataLong2);
%% 迭代开始
for t = 1: maxIter
fitness1=zeros(1,pop); 
%% 粒子1 
for j = (20 * (t - 1) * block1 + 1) : (20 * (t - 1) * block1 + block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(1)=fitness1(1)+e_cont1(j)^2/block1;
end
%% 粒子2
for j = (20 * (t - 1) * block1 + block1 + 1) : (20 * (t - 1) * block1 + 2 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(2)=fitness1(2)+e_cont1(j)^2/block1;
end
%% 粒子3
for j = (20 * (t - 1) * block1 + 2 * block1 + 1) : (20 * (t - 1) * block1 + 3 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(3)=fitness1(3)+e_cont1(j)^2/block1;
end
%% 粒子4
for j = (20 * (t - 1) * block1 + 3 * block1 + 1) : (20 * (t - 1) * block1 + 4 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(4)=fitness1(4)+e_cont1(j)^2/block1;
end
%% 粒子5
for j = (20 * (t - 1) * block1 + 4 * block1 + 1) : (20 * (t - 1) * block1 + 5 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(5)=fitness1(5)+e_cont1(j)^2/block1;
end
%% 粒子6
for j = (20 * (t - 1) * block1 + 5 * block1 + 1) : (20 * (t - 1) * block1 + 6 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(6)=fitness1(6)+e_cont1(j)^2/block1;
end
%% 粒子7
for j = (20 * (t - 1) * block1 + 6 * block1 + 1) : (20 * (t - 1) * block1 + 7 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(7)=fitness1(7)+e_cont1(j)^2/block1;
end
%% 粒子8
for j = (20 * (t - 1) * block1 + 7 * block1 + 1) : (20 * (t - 1) * block1 + 8 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(8)=fitness1(8)+e_cont1(j)^2/block1;
end
%% 粒子9
for j = (20 * (t - 1) * block1 + 8 * block1 + 1) : (20 * (t - 1) * block1 + 9 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(9)=fitness1(9)+e_cont1(j)^2/block1;
end
%% 粒子10
for j = (20 * (t - 1) * block1 + 9 * block1 + 1) : (20 * (t - 1) * block1 + 10 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(10)=fitness1(10)+e_cont1(j)^2/block1;
end
%% 粒子11
for j = (20 * (t - 1) * block1 + 10 * block1 + 1) : (20 * (t - 1) * block1 + 11 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(11)=fitness1(11)+e_cont1(j)^2/block1;
end
%% 粒子12
for j = (20 * (t - 1) * block1 + 11 * block1 + 1) : (20 * (t - 1) * block1 + 12 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(12)=fitness1(12)+e_cont1(j)^2/block1;
end
%% 粒子13
for j = (20 * (t - 1) * block1 + 12 * block1 + 1) : (20 * (t - 1) * block1 + 13 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(13)=fitness1(13)+e_cont1(j)^2/block1;
end
%% 粒子14
for j = (20 * (t - 1) * block1 + 13 * block1 + 1) : (20 * (t - 1) * block1 + 14 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(14)=fitness1(14)+e_cont1(j)^2/block1;
end
%% 粒子15
for j = (20 * (t - 1) * block1 + 14 * block1 + 1) : (20 * (t - 1) * block1 + 15 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(15)=fitness1(15)+e_cont1(j)^2/block1;
end
%% 粒子16
for j = (20 * (t - 1) * block1 + 15 * block1 + 1) : (20 * (t - 1) * block1 + 16 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(16)=fitness1(16)+e_cont1(j)^2/block1;
end
%% 粒子17
for j = (20 * (t - 1) * block1 + 16 * block1 + 1) : (20 * (t - 1) * block1 + 17 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(17)=fitness1(17)+e_cont1(j)^2/block1;
end
%% 粒子18
for j = (20 * (t - 1) * block1 + 17 * block1 + 1) : (20 * (t - 1) * block1 + 18 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(18)=fitness1(18)+e_cont1(j)^2/block1;
end
%% 粒子19
for j = (20 * (t - 1) * block1 + 18 * block1 + 1) : (20 * (t - 1) * block1 + 19 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(19)=fitness1(19)+e_cont1(j)^2/block1;
end
%% 粒子20
for j = (20 * (t - 1) * block1 + 19 * block1 + 1) : (20 * (t - 1) * block1 + 20 * block1)
    ANCCx1=[X_noise(j) ANCCx1(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx1';   
    ActualSPx1=[ANCCy ActualSPx1(:,(1:(Tap_Sw-1)))]; 
    e_cont1(j)=Yp(j)-ActualSPx1*Sw';
    fitness1(20)=fitness1(20)+e_cont1(j)^2/block1;
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
IterCurve1(t) = IterCurve1(t) + Alpha_score;
%% 更新
a = 2 * (1 - t / maxIter );

for i = 1:pop
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
%% GWO算法 lambda = 20
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 初始化    
ANCCx2=zeros(1,dim);
ActualSPx2=zeros(1,Tap_Sw);  
e_cont2=zeros(1,DataLong2);
%% 迭代开始
for t = 1: maxIter
fitness2=zeros(1,pop); 
%% 粒子1 
for j = (20 * (t - 1) * block2 + 1) : (20 * (t - 1) * block2 + block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(1)=fitness2(1)+e_cont2(j)^2/block2;
end
%% 粒子2
for j = (20 * (t - 1) * block2 + block2 + 1) : (20 * (t - 1) * block2 + 2 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(2)=fitness2(2)+e_cont2(j)^2/block2;
end
%% 粒子3
for j = (20 * (t - 1) * block2 + 2 * block2 + 1) : (20 * (t - 1) * block2 + 3 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(3)=fitness2(3)+e_cont2(j)^2/block2;
end
%% 粒子4
for j = (20 * (t - 1) * block2 + 3 * block2 + 1) : (20 * (t - 1) * block2 + 4 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(4)=fitness2(4)+e_cont2(j)^2/block2;
end
%% 粒子5
for j = (20 * (t - 1) * block2 + 4 * block2 + 1) : (20 * (t - 1) * block2 + 5 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(5)=fitness2(5)+e_cont2(j)^2/block2;
end
%% 粒子6
for j = (20 * (t - 1) * block2 + 5 * block2 + 1) : (20 * (t - 1) * block2 + 6 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(6)=fitness2(6)+e_cont2(j)^2/block2;
end
%% 粒子7
for j = (20 * (t - 1) * block2 + 6 * block2 + 1) : (20 * (t - 1) * block2 + 7 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(7)=fitness2(7)+e_cont2(j)^2/block2;
end
%% 粒子8
for j = (20 * (t - 1) * block2 + 7 * block2 + 1) : (20 * (t - 1) * block2 + 8 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(8)=fitness2(8)+e_cont2(j)^2/block2;
end
%% 粒子9
for j = (20 * (t - 1) * block2 + 8 * block2 + 1) : (20 * (t - 1) * block2 + 9 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(9)=fitness2(9)+e_cont2(j)^2/block2;
end
%% 粒子10
for j = (20 * (t - 1) * block2 + 9 * block2 + 1) : (20 * (t - 1) * block2 + 10 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(10)=fitness2(10)+e_cont2(j)^2/block2;
end
%% 粒子11
for j = (20 * (t - 1) * block2 + 10 * block2 + 1) : (20 * (t - 1) * block2 + 11 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(11)=fitness2(11)+e_cont2(j)^2/block2;
end
%% 粒子12
for j = (20 * (t - 1) * block2 + 11 * block2 + 1) : (20 * (t - 1) * block2 + 12 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(12)=fitness2(12)+e_cont2(j)^2/block2;
end
%% 粒子13
for j = (20 * (t - 1) * block2 + 12 * block2 + 1) : (20 * (t - 1) * block2 + 13 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(13)=fitness2(13)+e_cont2(j)^2/block2;
end
%% 粒子14
for j = (20 * (t - 1) * block2 + 13 * block2 + 1) : (20 * (t - 1) * block2 + 14 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(14)=fitness2(14)+e_cont2(j)^2/block2;
end
%% 粒子15
for j = (20 * (t - 1) * block2 + 14 * block2 + 1) : (20 * (t - 1) * block2 + 15 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(15)=fitness2(15)+e_cont2(j)^2/block2;
end
%% 粒子16
for j = (20 * (t - 1) * block2 + 15 * block2 + 1) : (20 * (t - 1) * block2 + 16 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(16)=fitness2(16)+e_cont2(j)^2/block2;
end
%% 粒子17
for j = (20 * (t - 1) * block2 + 16 * block2 + 1) : (20 * (t - 1) * block2 + 17 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(17)=fitness2(17)+e_cont2(j)^2/block2;
end
%% 粒子18
for j = (20 * (t - 1) * block2 + 17 * block2 + 1) : (20 * (t - 1) * block2 + 18 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(18)=fitness2(18)+e_cont2(j)^2/block2;
end
%% 粒子19
for j = (20 * (t - 1) * block2 + 18 * block2 + 1) : (20 * (t - 1) * block2 + 19 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(19)=fitness2(19)+e_cont2(j)^2/block2;
end
%% 粒子20
for j = (20 * (t - 1) * block2 + 19 * block2 + 1) : (20 * (t - 1) * block2 + 20 * block2)
    ANCCx2=[X_noise(j) ANCCx2(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx2';   
    ActualSPx2=[ANCCy ActualSPx2(:,(1:(Tap_Sw-1)))]; 
    e_cont2(j)=Yp(j)-ActualSPx2*Sw';
    fitness2(20)=fitness2(20)+e_cont2(j)^2/block2;
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
IterCurve2(t) = IterCurve2(t) + Alpha_score;
%% 更新
a = 2 * (1 - t / maxIter );

for i = 1:pop
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
%% GWO算法 lambda = 30
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 初始化    
ANCCx3=zeros(1,dim);
ActualSPx3=zeros(1,Tap_Sw);  
e_cont3=zeros(1,DataLong3);
%% 迭代开始
for t = 1: maxIter
fitness3=zeros(1,pop); 
%% 粒子1 
for j = (20 * (t - 1) * block3 + 1) : (20 * (t - 1) * block3 + block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(1)=fitness3(1)+e_cont3(j)^2/block3;
end
%% 粒子2
for j = (20 * (t - 1) * block3 + block3 + 1) : (20 * (t - 1) * block3 + 2 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(2)=fitness3(2)+e_cont3(j)^2/block3;
end
%% 粒子3
for j = (20 * (t - 1) * block3 + 2 * block3 + 1) : (20 * (t - 1) * block3 + 3 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(3)=fitness3(3)+e_cont3(j)^2/block3;
end
%% 粒子4
for j = (20 * (t - 1) * block3 + 3 * block3 + 1) : (20 * (t - 1) * block3 + 4 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(4)=fitness3(4)+e_cont3(j)^2/block3;
end
%% 粒子5
for j = (20 * (t - 1) * block3 + 4 * block3 + 1) : (20 * (t - 1) * block3 + 5 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(5)=fitness3(5)+e_cont3(j)^2/block3;
end
%% 粒子6
for j = (20 * (t - 1) * block3 + 5 * block3 + 1) : (20 * (t - 1) * block3 + 6 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(6)=fitness3(6)+e_cont3(j)^2/block3;
end
%% 粒子7
for j = (20 * (t - 1) * block3 + 6 * block3 + 1) : (20 * (t - 1) * block3 + 7 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(7)=fitness3(7)+e_cont3(j)^2/block3;
end
%% 粒子8
for j = (20 * (t - 1) * block3 + 7 * block3 + 1) : (20 * (t - 1) * block3 + 8 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(8)=fitness3(8)+e_cont3(j)^2/block3;
end
%% 粒子9
for j = (20 * (t - 1) * block3 + 8 * block3 + 1) : (20 * (t - 1) * block3 + 9 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(9)=fitness3(9)+e_cont3(j)^2/block3;
end
%% 粒子10
for j = (20 * (t - 1) * block3 + 9 * block3 + 1) : (20 * (t - 1) * block3 + 10 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(10)=fitness3(10)+e_cont3(j)^2/block3;
end
%% 粒子11
for j = (20 * (t - 1) * block3 + 10 * block3 + 1) : (20 * (t - 1) * block3 + 11 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(11)=fitness3(11)+e_cont3(j)^2/block3;
end
%% 粒子12
for j = (20 * (t - 1) * block3 + 11 * block3 + 1) : (20 * (t - 1) * block3 + 12 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(12)=fitness3(12)+e_cont3(j)^2/block3;
end
%% 粒子13
for j = (20 * (t - 1) * block3 + 12 * block3 + 1) : (20 * (t - 1) * block3 + 13 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(13)=fitness3(13)+e_cont3(j)^2/block3;
end
%% 粒子14
for j = (20 * (t - 1) * block3 + 13 * block3 + 1) : (20 * (t - 1) * block3 + 14 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(14)=fitness3(14)+e_cont3(j)^2/block3;
end
%% 粒子15
for j = (20 * (t - 1) * block3 + 14 * block3 + 1) : (20 * (t - 1) * block3 + 15 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(15)=fitness3(15)+e_cont3(j)^2/block3;
end
%% 粒子16
for j = (20 * (t - 1) * block3 + 15 * block3 + 1) : (20 * (t - 1) * block3 + 16 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(16)=fitness3(16)+e_cont3(j)^2/block3;
end
%% 粒子17
for j = (20 * (t - 1) * block3 + 16 * block3 + 1) : (20 * (t - 1) * block3 + 17 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(17)=fitness3(17)+e_cont3(j)^2/block3;
end
%% 粒子18
for j = (20 * (t - 1) * block3 + 17 * block3 + 1) : (20 * (t - 1) * block3 + 18 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(18)=fitness3(18)+e_cont3(j)^2/block3;
end
%% 粒子19
for j = (20 * (t - 1) * block3 + 18 * block3 + 1) : (20 * (t - 1) * block3 + 19 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(19)=fitness3(19)+e_cont3(j)^2/block3;
end
%% 粒子20
for j = (20 * (t - 1) * block3 + 19 * block3 + 1) : (20 * (t - 1) * block3 + 20 * block3)
    ANCCx3=[X_noise(j) ANCCx3(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx3';   
    ActualSPx3=[ANCCy ActualSPx3(:,(1:(Tap_Sw-1)))]; 
    e_cont3(j)=Yp(j)-ActualSPx3*Sw';
    fitness3(20)=fitness3(20)+e_cont3(j)^2/block3;
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
IterCurve3(t) = IterCurve3(t) + Alpha_score;
%% 更新
a = 2 * (1 - t / maxIter );

for i = 1:pop
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
%% GWO算法 lambda = 50
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 初始化    
ANCCx4=zeros(1,dim);
ActualSPx4=zeros(1,Tap_Sw);  
e_cont4=zeros(1,DataLong4);
%% 迭代开始
for t = 1: maxIter
fitness4=zeros(1,pop); 
%% 粒子1 
for j = (20 * (t - 1) * block4 + 1) : (20 * (t - 1) * block4 + block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(1)=fitness4(1)+e_cont4(j)^2/block4;
end
%% 粒子2
for j = (20 * (t - 1) * block4 + block4 + 1) : (20 * (t - 1) * block4 + 2 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(2)=fitness4(2)+e_cont4(j)^2/block4;
end
%% 粒子3
for j = (20 * (t - 1) * block4 + 2 * block4 + 1) : (20 * (t - 1) * block4 + 3 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(3)=fitness4(3)+e_cont4(j)^2/block4;
end
%% 粒子4
for j = (20 * (t - 1) * block4 + 3 * block4 + 1) : (20 * (t - 1) * block4 + 4 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(4)=fitness4(4)+e_cont4(j)^2/block4;
end
%% 粒子5
for j = (20 * (t - 1) * block4 + 4 * block4 + 1) : (20 * (t - 1) * block4 + 5 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(5)=fitness4(5)+e_cont4(j)^2/block4;
end
%% 粒子6
for j = (20 * (t - 1) * block4 + 5 * block4 + 1) : (20 * (t - 1) * block4 + 6 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(6)=fitness4(6)+e_cont4(j)^2/block4;
end
%% 粒子7
for j = (20 * (t - 1) * block4 + 6 * block4 + 1) : (20 * (t - 1) * block4 + 7 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(7)=fitness4(7)+e_cont4(j)^2/block4;
end
%% 粒子8
for j = (20 * (t - 1) * block4 + 7 * block4 + 1) : (20 * (t - 1) * block4 + 8 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(8)=fitness4(8)+e_cont4(j)^2/block4;
end
%% 粒子9
for j = (20 * (t - 1) * block4 + 8 * block4 + 1) : (20 * (t - 1) * block4 + 9 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(9)=fitness4(9)+e_cont4(j)^2/block4;
end
%% 粒子10
for j = (20 * (t - 1) * block4 + 9 * block4 + 1) : (20 * (t - 1) * block4 + 10 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(10)=fitness4(10)+e_cont4(j)^2/block4;
end
%% 粒子11
for j = (20 * (t - 1) * block4 + 10 * block4 + 1) : (20 * (t - 1) * block4 + 11 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(11)=fitness4(11)+e_cont4(j)^2/block4;
end
%% 粒子12
for j = (20 * (t - 1) * block4 + 11 * block4 + 1) : (20 * (t - 1) * block4 + 12 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(12)=fitness4(12)+e_cont4(j)^2/block4;
end
%% 粒子13
for j = (20 * (t - 1) * block4 + 12 * block4 + 1) : (20 * (t - 1) * block4 + 13 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(13)=fitness4(13)+e_cont4(j)^2/block4;
end
%% 粒子14
for j = (20 * (t - 1) * block4 + 13 * block4 + 1) : (20 * (t - 1) * block4 + 14 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(14)=fitness4(14)+e_cont4(j)^2/block4;
end
%% 粒子15
for j = (20 * (t - 1) * block4 + 14 * block4 + 1) : (20 * (t - 1) * block4 + 15 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(15)=fitness4(15)+e_cont4(j)^2/block4;
end
%% 粒子16
for j = (20 * (t - 1) * block4 + 15 * block4 + 1) : (20 * (t - 1) * block4 + 16 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(16)=fitness4(16)+e_cont4(j)^2/block4;
end
%% 粒子17
for j = (20 * (t - 1) * block4 + 16 * block4 + 1) : (20 * (t - 1) * block4 + 17 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(17)=fitness4(17)+e_cont4(j)^2/block4;
end
%% 粒子18
for j = (20 * (t - 1) * block4 + 17 * block4 + 1) : (20 * (t - 1) * block4 + 18 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(18)=fitness4(18)+e_cont4(j)^2/block4;
end
%% 粒子19
for j = (20 * (t - 1) * block4 + 18 * block4 + 1) : (20 * (t - 1) * block4 + 19 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(19)=fitness4(19)+e_cont4(j)^2/block4;
end
%% 粒子20
for j = (20 * (t - 1) * block4 + 19 * block4 + 1) : (20 * (t - 1) * block4 + 20 * block4)
    ANCCx4=[X_noise(j) ANCCx4(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx4';   
    ActualSPx4=[ANCCy ActualSPx4(:,(1:(Tap_Sw-1)))]; 
    e_cont4(j)=Yp(j)-ActualSPx4*Sw';
    fitness4(20)=fitness4(20)+e_cont4(j)^2/block4;
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
IterCurve4(t) = IterCurve4(t) + Alpha_score;
%% 更新
a = 2 * (1 - t / maxIter );

for i = 1:pop
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
%% GWO算法 lambda = 100
%% 初始化种群位置
X=initialize_population(pop,ub,lb,dim);
%% 初始化    
ANCCx5=zeros(1,dim);
ActualSPx5=zeros(1,Tap_Sw);  
e_cont5=zeros(1,DataLong5);
%% 迭代开始
for t = 1: maxIter
fitness5=zeros(1,pop); 
%% 粒子1 
for j = (20 * (t - 1) * block5 + 1) : (20 * (t - 1) * block5 + block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(1,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(1)=fitness5(1)+e_cont5(j)^2/block5;
end
%% 粒子2
for j = (20 * (t - 1) * block5 + block5 + 1) : (20 * (t - 1) * block5 + 2 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(2,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(2)=fitness5(2)+e_cont5(j)^2/block5;
end
%% 粒子3
for j = (20 * (t - 1) * block5 + 2 * block5 + 1) : (20 * (t - 1) * block5 + 3 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(3,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(3)=fitness5(3)+e_cont5(j)^2/block5;
end
%% 粒子4
for j = (20 * (t - 1) * block5 + 3 * block5 + 1) : (20 * (t - 1) * block5 + 4 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(4,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(4)=fitness5(4)+e_cont5(j)^2/block5;
end
%% 粒子5
for j = (20 * (t - 1) * block5 + 4 * block5 + 1) : (20 * (t - 1) * block5 + 5 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(5,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(5)=fitness5(5)+e_cont5(j)^2/block5;
end
%% 粒子6
for j = (20 * (t - 1) * block5 + 5 * block5 + 1) : (20 * (t - 1) * block5 + 6 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(6,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(6)=fitness5(6)+e_cont5(j)^2/block5;
end
%% 粒子7
for j = (20 * (t - 1) * block5 + 6 * block5 + 1) : (20 * (t - 1) * block5 + 7 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(7,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(7)=fitness5(7)+e_cont5(j)^2/block5;
end
%% 粒子8
for j = (20 * (t - 1) * block5 + 7 * block5 + 1) : (20 * (t - 1) * block5 + 8 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(8,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(8)=fitness5(8)+e_cont5(j)^2/block5;
end
%% 粒子9
for j = (20 * (t - 1) * block5 + 8 * block5 + 1) : (20 * (t - 1) * block5 + 9 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(9,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(9)=fitness5(9)+e_cont5(j)^2/block5;
end
%% 粒子10
for j = (20 * (t - 1) * block5 + 9 * block5 + 1) : (20 * (t - 1) * block5 + 10 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(10,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(10)=fitness5(10)+e_cont5(j)^2/block5;
end
%% 粒子11
for j = (20 * (t - 1) * block5 + 10 * block5 + 1) : (20 * (t - 1) * block5 + 11 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(11,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(11)=fitness5(11)+e_cont5(j)^2/block5;
end
%% 粒子12
for j = (20 * (t - 1) * block5 + 11 * block5 + 1) : (20 * (t - 1) * block5 + 12 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(12,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(12)=fitness5(12)+e_cont5(j)^2/block5;
end
%% 粒子13
for j = (20 * (t - 1) * block5 + 12 * block5 + 1) : (20 * (t - 1) * block5 + 13 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(13,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(13)=fitness5(13)+e_cont5(j)^2/block5;
end
%% 粒子14
for j = (20 * (t - 1) * block5 + 13 * block5 + 1) : (20 * (t - 1) * block5 + 14 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(14,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(14)=fitness5(14)+e_cont5(j)^2/block5;
end
%% 粒子15
for j = (20 * (t - 1) * block5 + 14 * block5 + 1) : (20 * (t - 1) * block5 + 15 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(15,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(15)=fitness5(15)+e_cont5(j)^2/block5;
end
%% 粒子16
for j = (20 * (t - 1) * block5 + 15 * block5 + 1) : (20 * (t - 1) * block5 + 16 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(16,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(16)=fitness5(16)+e_cont5(j)^2/block5;
end
%% 粒子17
for j = (20 * (t - 1) * block5 + 16 * block5 + 1) : (20 * (t - 1) * block5 + 17 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(17,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(17)=fitness5(17)+e_cont5(j)^2/block5;
end
%% 粒子18
for j = (20 * (t - 1) * block5 + 17 * block5 + 1) : (20 * (t - 1) * block5 + 18 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(18,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(18)=fitness5(18)+e_cont5(j)^2/block5;
end
%% 粒子19
for j = (20 * (t - 1) * block5 + 18 * block5 + 1) : (20 * (t - 1) * block5 + 19 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(19,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(19)=fitness5(19)+e_cont5(j)^2/block5;
end
%% 粒子20
for j = (20 * (t - 1) * block5 + 19 * block5 + 1) : (20 * (t - 1) * block5 + 20 * block5)
    ANCCx5=[X_noise(j) ANCCx5(1:(dim-1))];
    ANCCy=X(20,:)*ANCCx5';   
    ActualSPx5=[ANCCy ActualSPx5(:,(1:(Tap_Sw-1)))]; 
    e_cont5(j)=Yp(j)-ActualSPx5*Sw';
    fitness5(20)=fitness5(20)+e_cont5(j)^2/block5;
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

for i = 1:pop
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
legend('{\fontname{Times New Roman}\lambda=10}','{\fontname{Times New Roman}\lambda=20}','{\fontname{Times New Roman}\lambda=30}','{\fontname{Times New Roman}\lambda=50}','{\fontname{Times New Roman}\lambda=100}',...%    
     'NumColumns',1,'FontSize', 12, 'Location', 'northeast');
grid on
set(gca, 'GridLineStyle', ':');  % 设置为虚线 
set(gca, 'GridAlpha', 1);  % 设置透明度