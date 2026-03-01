function X = initialize_population(pop, ub, lb, dim)
% Initialize population
% pop: population size
% dim: dimension of each particle
% ub: upper bound of variables for each dimension
% lb: lower bound of variables for each dimension
    for i = 1:pop
        for j = 1:dim
            X(i, j) = lb(j) + (ub(j)-lb(j))*randn();
        end
    end
end