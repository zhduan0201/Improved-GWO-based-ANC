%% 粒子群初始化函数
function X=initialize_population(pop,ub,lb,dim)
% pop种群数量
% dim每个粒子的维度
% ub每个维度的变量上边界
% lb每个维度的变量下边界，维度为[1,dim]
% X输出的种群，维度为[pop,dim]
for i=1:pop
    for j=1:dim
        X(i,j)=lb(j)+(ub(j)-lb(j))*randn();
    end
end