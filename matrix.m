% below is the general formula to create an nxn matrix with random values
% between a and b
% r = a + (b-a).*rand(N,N) 
% in this case, a = -0.4 and b = 0.8
A = -4 + (0.8+0.4)*rand(500,500);
z=ones(500,1);
b=A*z;
% Solve Equation by partial pivoting
tic 
% Find the solution 
x = partial_pivot(A,b);
% Get the SMSE using the formula given in the question
smse = sqrt(sum((x-z).^2));
fprintf('\n         ==================================================================================');
fprintf('\n                Solve Ax = b using Gaussian Elimination Method with partial pivoting ');
fprintf('\n         ==================================================================================');
fprintf('\n         ----------------------------------------------------------------------------------');
fprintf('\nSMSE obtained through partial pivoting = %.16f\n',smse)
toc
% Solve Equation without partial pivoting
tic 
% Find the solution 
x = no_partial_pivot(A,b);
% Get the SMSE
smse = sqrt(sum((x-z).^2));
fprintf('\n         ==================================================================================');
fprintf('\n             Solve Ax = b using Gaussian Elimination Method without partial pivoting ');
fprintf('\n         ==================================================================================');
fprintf(' \n        ----------------------------------------------------------------------------------');
fprintf('\nSMSE obtained without pivoting = %.16f\n',smse)
toc
%
%---------------------main_2.m----------------------------------------------------------------------------

% Solve Equation without partial pivoting
%clc;clear;
A = -4 + (0.8+0.4)*rand(500,500);
z=ones(500,1);
b=A*z;
tic 
% Find the solution 
x = linsolve(A,b);
% Get the SMSE
smse = sqrt(sum((x-z).^2));
fprintf('\n         ==================================================================================');
fprintf('\n                               Solve Ax = b MATLAB Linsolve ');
fprintf('\n         ==================================================================================');
fprintf('\nSMSE obtained using Linsolve = %.16f\n',smse)
toc
%


function x = partial_pivot(A, b)
    % Program to Solve Ax = b using Gaussian Elimination Method with partial-
    % pivoting (Row-exchange)
    %==========================================================================
    % INPUTS:
    %==========================================================================
    % Right handside vector, b (n-by-1)
    % Solution vector, x
    %==========================================================================
    [m,n] = size(A);
    %==========================================================================
    % Initialization
    x = zeros(m,1);
    l = zeros(m,m-1);
    %==========================================================================
    % Main Program
    %==========================================================================
    % Reducing Matrix A to upper triangular form
    for k = 1:m-1
        % =========Performing Partial-pivoting=================================
            for p = k+1:m
                if (abs(A(k,k)) < abs(A(p,k)))
                    A([k p],:) = A([p k],:);
                      b([k p]) = b([p k]);
                end
            end
        % =====================================================================
        for i = k+1:m
            l(i,k) = A(i,k)/A(k,k);
            for j = k+1:n
                A(i,j) = A(i,j)-l(i,k)*A(k,j);
            end  
            b(i) = b(i)-l(i,k)*b(k);
        end
    end
    for k = 1:m-1
        for i = k+1:m
            A(i,k) = 0;
        end
    end
    %==========================================================================
    % Backward substitution
    %==========================================================================
    x(m) = b(m)/A(m,m);
    for i = m-1:-1:1
        sum = 0;
        for j = i+1:m
            sum = sum + A(i,j)*x(j);
        end
        x(i) = (b(i)- sum)/A(i,i);
    end
end

function x = no_partial_pivot(A, b)
    % Solve linear system Ax = b
    % using Gaussian elimination without pivoting
    % A is an n by n matrix
    % b is an n by k matrix (k copies of n-vectors)
    % x is an n by k matrix (k copies of solution vectors)
    [n, n] = size(A);     % Find size of matrix A
    [n, k] = size(b);     % Find size of matrix b
    x = zeros(n,k);      % Initialize x
    for i = 1:n-1
        m = -A(i+1:n,i)/A(i,i); % multipliers
        A(i+1:n,:) = A(i+1:n,:) + m*A(i,:);
        b(i+1:n,:) = b(i+1:n,:) + m*b(i,:);
    end
    % Use back substitution to find unknowns
    x(n,:) = b(n,:)/A(n,n);
    for i = n-1:-1:1
        x(i,:) = (b(i,:) - A(i,i+1:n)*x(i+1:n,:))/A(i,i);
    end
end
