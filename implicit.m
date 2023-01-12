function [T] = implicit(CFL,N,M,T)
A = zeros(N-2,N-2);
v = (1+2*CFL) * ones(N-2,1);
A = diag(v);

for i = 1:N-2
    for j = 1:N-2
        if j == i+1 || i == j+1
            A(i,j) = -CFL;
        end
    end
end

for k = 1:M-1
    for i = 2:N-1
        if i == 2
            b(i-1,k) = T(i,k) + CFL * T(i-1,k+1);
        elseif i == N-1
            %i
            b(i-1,k) = T(i,k) + CFL * T(i+1,k+1);
        else
            b(i-1,k) = T(i,k);
        end
    end
    T(2:N-1,k+1) = mldivide(A,b(:,k));
end
end

