%% 边界检查函数
function[X]=boundary_constraints(x,ub,lb,dim)
for i=1:dim
    if x(i)>ub(i)
        x(i)=ub(i);
    end
    if x(i)<lb(i)
        x(i)=lb(i);
    end
end
X=x;
end
