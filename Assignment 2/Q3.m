%enter the value of n where n=2,4,6,8

display = 'enter the value of n';
n = input(display);
A = hilb(n)
e = ones(n,1);
b = A*e;
L = eye(n,n);

for j = 1:n
    if (A(j,j) - sum(L(j,1:j-1) .^ 2)) < 0
        disp('taking sqaure root of a negative number. Quitting')
    end
    
    L(j,j) = sqrt(A(j,j) - sum(L(j,1:j-1) .^ 2));
    for i = j+1:n
        if L(j,j) == 0
            disp('Dividing by zero. Quitting.')
        end
        
        L(i,j) = (A(i,j) - sum(L(i,1:j-1) .* L(j,1:j-1)))/ L(j,j);
         
    end
end

% L transpose
L_t = transpose(L);
A_estimated = L * L_t;
%solving linear equation

y = linsolve(L,b);
x_c = linsolve(L_t,y);

%calculating the results
disp('the metrics: relative error, relative residual and reltive matrix residual are calculated below')
relative_error = norm(x_c - e)/norm(e)
relative_residual = (norm(b - A*x_c)) /  ( norm(A,'fro')* norm(x_c))
relative_matrix_residual = norm(A - (L*L_t),'fro')/ norm(A,'fro')

%condition numbers
disp('condition number')
cond(A)
cond(A_estimated)







