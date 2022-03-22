clear all;
close all;
clc;
a = -1;
b = 1;
c = 1;
d = 2;
N = 1000000;
x = zeros(1 , N);
y = zeros(1 , N);
f = zeros(1 , N);
for i = 1:N
    x (i) = a + (b - a) * rand;
    y (i) = c + (d - c) * rand;
    f (i) = x(i) * y(i);
end 
Ef = (1/N)*sum (f);
I = (b - a)*(d - c) * Ef
fun = @(x,y) x.*y;
q = integral2(fun, -1,1,1,2)
