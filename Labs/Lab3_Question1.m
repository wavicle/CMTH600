clear all;
close all;
clc;
a = -2;
b = 2;
N = 10000000;
z = zeros(1 , N);
x = zeros(1 , N);
f = zeros(1 , N);
for i = 1:N
    z (i) = rand;
    x (i) = a + (b - a) * z(i);
    f (i) = exp (x(i));
end 
I = (1/N)*sum (f*(b-a))
fun = @(x) exp(x);
q = integral(fun, -2, 2)
