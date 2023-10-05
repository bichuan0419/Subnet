clc;clear;close all
%%
A = 10*rand(100,100);
A(1:20,1:20) = 2*A(1:20,1:20);
A(21:30, 21:30) = 3*A(21:30, 21:30);
% A(21:30,21:30) = 2.5*A(21:30,21:30) ;
A = A - diag(diag(A));
A = A+A';
A_vec = squareform(A);
order_shuffled = randperm(length(A));
A = A(order_shuffled, order_shuffled);
figure;imagesc(A);colorbar;
writematrix(A,'A.csv');
save('A.mat','A');