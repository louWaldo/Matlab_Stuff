fprintf('SOLVING WITH PARTIAL PIVOT, double precision \n');
A1 = zeros(500, 500);
z1=ones(500,1);
b1=A1*z1;
%creating ill-conditioned matrix
for i = 1 : 500

    for j = 1: 500

        A1(i, j) =1/( i + j-1);

    end

end

tic
x1 = partial_pivot(A1, b1);
toc
smse_1 = sqrt(sum((x1-z1).^2));
fprintf('SMSE = %.32f\n',smse_1)
fprintf('-------------------------------------------------------------------------------\n');


fprintf('SOLVING WITH PARTIAL PIVOT, single precision \n');
A_S1 = single(zeros(500, 500));
z_S1=single(ones(500,1));
b_S1=single(A_S1*z_S1);

for i = 1 : 500

    for j = 1: 500

        A_S1(i, j) =1/( i + j-1);

    end

end
tic
x_S1 = single(partial_pivot(A_S1, b_S1));
toc
smse_S1 = single(sqrt(sum((x_S1-z_S1).^2)));
fprintf('SMSE = %.32f\n',smse_S1)
fprintf('========================================================================================\n');



fprintf('SOLVING WITH LINSOLVE, double precision \n')
A2 = zeros(500, 500);
z2=ones(500,1);
b2=A2*z2;
for i = 1 : 500

    for j = 1: 500

        A2(i, j) =1/( i + j-1);

    end

end


tic;
x2 = linsolve(A2,b2);
toc;
smse_2 = sqrt(sum((x2-z2).^2));
fprintf('SMSE = %.32f\n',smse_2)

fprintf('-------------------------------------------------------------------------------\n');

fprintf('SOLVING WITH LINSOLVE, single precision \n')
A_S2 = single(zeros(500, 500));
z_S2=single(ones(500,1));
b_S2=single(A_S2*z_S2);
for i = 1 : 500

    for j = 1: 500

        A_S2(i, j) =1/( i + j-1);

    end

end

tic
x_S2 = single(linsolve(A_S2,b_S2));
toc
smse_S2 = single(sqrt(sum((x_S2-z_S2).^2)));
fprintf('SMSE = %.32f\n',smse_S2)


function x = partial_pivot(A, b)
    [m,n] = size(A);
    x = zeros(m,1);
    l = zeros(m,m-1);
    % Reducing Matrix A to upper triangular form
    for k = 1:m-1
        %finding element with largest absolut value, make that the pivot
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
    % Backward substitution
    x(m) = b(m)/A(m,m);
    for i = m-1:-1:1
        sum = 0;
        for j = i+1:m
            sum = sum + A(i,j)*x(j);
        end
        x(i) = (b(i)- sum)/A(i,i);
    end
end
