% Question 4 Part B is solved first which is then followed by Question 4
% part C

% n equals 30 and matrix A0 is set below
% x is the calculated x which is x_c in the question
%actual exact solution is denoted by x_actual
n = 30;
x_actual = ones(n,1);
AO = eye(n,n);
for i = 2:n
    AO(i,1:i-1) = -1;
end
AO(:,n) = 1;

% Marix A and vector b are set below

A = AO + 10.^-8 .* (eye(n,n));
e = ones(n,1);
b = A * e;
L = eye(n,n);

% store matrix A as AA for later calulations in the end
AA = A;

% Guassian Elimination with Partial Pivoting

for j = 1:n-1
    %calculation of indexes and swapping of Matrix A
    [maxval, index] = max(abs(A(j:n,j)));
    temp = A(j,:);
    if j > 1
        index = index + j -1;
    end
    A(j,:) = A(index,:);
    A(index,:) = temp;
    
    if j>1
        tk = L(j,1:j-1);
        L(j,1:j-1) = L(index,1:j-1);
        L(index,1:j-1) = tk;
    end
    
        % GEPP:
    
    for i = j+1:n
              
        if A(j,j) == 0
            disp ('Division by zero to occur. Quitting')
        end        
        L(i,j) = A(i,j) / A(j,j);
        A(i,j) = 0;
        for k = j+1:n
            A(i,k) = A(i,k) - (L(i,j) * A(j,k));
        end
        
    end
    
end

% U matrix
U = A;

%solving the linear system of equations

y = linsolve(L,b);
x = linsolve(U,y);

%calucation of metrics

disp('fraction:relative error')
fraction = norm(x_actual - x) / norm(x_actual)

disp('growth factor')
growth_factor = max(abs(U(:)))/ max(abs(AA(:)))

% Question 4 part C

%iterative refinement

iter = 0;
stop_value = 10 ^ (-15);
condition = true;
while condition
     r = b - A*x;
     s = linsolve(L,r);
     d = linsolve(U,s);
     y = x+d;
     x = y;
     test = norm(d)/norm(x);
     iter = iter + 1;
     condition = (test <= stop_value);
end

disp('final_relative_error')
final_relative_error = norm(x-x_actual)/norm(x)

disp('number of iterations')
iter
