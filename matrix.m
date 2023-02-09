%Luis Waldo 
%below is the general formula to create an nxn matrix with random values
% between a and b
% r = a + (b-a).*rand(N,N) 
% in this case, a = -0.4 and b = 0.8
A1 = -4 + (0.8+0.4)*rand(500,500);
z1=ones(500,1);
b1=A1*z1;
% Solve Equation by partial pivoting
fprintf('SOLVING WITH PARTIAL PIVOT, double precision \n');
tic
x1 = partial_pivot(A1,b1);
toc
smse_1 = sqrt(sum((x1-z1).^2));
fprintf('SMSE = %.32f\n',smse_1)
fprintf('-------------------------------------------------------------------------------\n');

A_S1 = single(-4 + (0.8+0.4)*rand(500,500));
z_S1=single(ones(500,1));
b_S1=single(A_S1*z_S1);
fprintf('SOLVING WITH PARTIAL PIVOT, single precision \n');
tic
x_S1 = single(partial_pivot(A_S1,b_S1));
toc
smse_S1 = single(sqrt(sum((x_S1-z_S1).^2)));
fprintf('SMSE = %.32f\n',smse_S1)

fprintf('\n==================================================================================\n\n');


% Solving without partial pivoting
A2 = -4 + (0.8+0.4)*rand(500,500);
z2=ones(500,1);
b2=A2*z2;
fprintf('SOLVING WITHOUT PARTIAL PIVOT \n')
tic
x2 = standard_gauss(A2,b2);
toc
smse_2 = sqrt(sum((x2-z2).^2));
fprintf('SMSE = %.32f\n',smse_2)
fprintf('-------------------------------------------------------------------------------\n');

A_S2 = single(-4 + (0.8+0.4)*rand(500,500));
z_S2=single(ones(500,1));
b_S2=single(A_S2*z_S2);
fprintf('SOLVING WITHOUT PARTIAL PIVOT, single precision \n');
tic
x_S2 = single(standard_gauss(A_S2,b_S2));
toc
smse_S2 = single(sqrt(sum((x_S2-z_S2).^2)));
fprintf('SMSE = %.32f\n',smse_S2)
fprintf('\n==================================================================================\n\n');



%Solving equation with linsolve
A3 = -4 + (0.8+0.4)*rand(500,500);
z3=ones(500,1);
b3=A3*z3;
fprintf('SOLVING WITH LINSOLVE, double precision \n')
tic;
x3 = linsolve(A3,b3);
toc;
smse_3 = sqrt(sum((x3-z3).^2));
fprintf('SMSE = %.32f\n',smse_3)

fprintf('-------------------------------------------------------------------------------\n');

fprintf('SOLVING WITH LINSOLVE, single precision \n')
A_S3 = single(-4 + (0.8+0.4)*rand(500,500));
z_S3=single(ones(500,1));
b_S3=single(A_S3*z_S3);
tic
x_S3 = single(linsolve(A_S3,b_S3));
toc
smse_S3 = single(sqrt(sum((x_S3-z_S3).^2)));
fprintf('SMSE = %.32f\n',smse_S3)




function x = partial_pivot(A, b)
    [m,n] = size(A);
    x = zeros(m,1);
    l = zeros(m,m-1);
    % Reducing Matrix A to upper triangular form
    for k = 1:m-1
        %finding element with largest absolute value, make that the pivot (partial pivoting)
            for p = k+1:m
                if (abs(A(k,k)) < abs(A(p,k)))
                    A([k p],:) = A([p k],:);
                      b([k p]) = b([p k]);
                end
            end
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
    %backward substitution
    x(m) = b(m)/A(m,m);
    for i = m-1:-1:1
        sum = 0;
        for j = i+1:m
            sum = sum + A(i,j)*x(j);
        end
        x(i) = (b(i)- sum)/A(i,i);
    end
end

function x = standard_gauss(A,b)
    n = length(b); 
    x = zeros(n,1);
    for k=1:n-1 
          for i=k+1:n
               xmult = A(i,k)/A(k,k);
               for j=k+1:n
                 A(i,j) = A(i,j)-xmult*A(k,j);
               end
               b(i) = b(i)-xmult*b(k);
          end
    end
     % back substitution
    x(n) = b(n)/A(n,n);
    for i=n-1:-1:1
       sum = b(i);
       for j=i+1:n
         sum = sum-A(i,j)*x(j);
       end
       x(i) = sum/A(i,i);
    end
end
