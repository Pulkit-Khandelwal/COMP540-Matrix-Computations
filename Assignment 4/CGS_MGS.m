% Define m and n for Vandermonde Matrices
m = 20;
n = 10;
% m = 30;
% n = 20;

A = zeros(m,n);
for i = 1:m
    for j =1:n
        A(i,j) = (j/n) ^ (i-1);
    end
end

% MATLAB function for QR factorization

[Q_m,R_m] = qr(A)
norm(A - (Q_m * R_m)) / norm(A)
norm((Q_m' * Q_m) - eye(m,m))

V = zeros(m,n);
Q = zeros(m,n);
R = zeros(n,n);

%CGS

for j = 1:n
    V(:,j) = A(:,j);
    for i = 1:j-1
        R(i,j) = (Q(:,i)') * A(:,j);
        V(:,j) = V(:,j) - R(i,j) * Q(:,i);
    end
    R(j,j) = norm(V(:,j));
    Q(:,j) = V(:,j)/R(j,j);
end


%MGS

for j = 1:n
    V(:,j) = A(:,j);
end

for i = 1:n
    R(i,i) = norm(V(:,i));
    Q(:,i) = V(:,i) / R(i,i);
    for j= i+1:n
        R(i,j) = Q(:,i)' * V(:,j);
        V(:,j) = V(:,j) - R(i,j)*Q(:,i);
    end
end

% Compute the results:

Q_1 = Q(:,1:n);
norm(A - Q_1 * R) / norm(A);
norm(((Q_1)' * Q_1) - eye(n,n));
   
        




  
        