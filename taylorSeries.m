%creating the range
x = linspace(-1*pi,1*pi);

%defining sin(x) as well as taylor series approximations
trueSinX = sin(x);
s1 = x;
s3 = x - (x.^3)/factorial(3);
s5 = x - (x.^3)/factorial(3) + (x.^5)/factorial(5);
s7 = x - (x.^3)/factorial(3) + (x.^5)/factorial(5)- (x.^7)/factorial(7);
s9 = x - (x.^3)/factorial(3) + (x.^5)/factorial(5)- (x.^7)/factorial(7) + (x.^9)/factorial(9);

%plotting approximations
figure
hold on
p1 = plot(x, s1, "red");
p3 = plot(x, s3, "green");
p5 = plot(x, s5, "blue");
p7 = plot(x, s7, "cyan");
p9 = plot(x, s9, "magenta");
p0 = plot(x, trueSinX, "black");

legend("s1", "s3", "s5", "s7", "s9", "trueSinX")

