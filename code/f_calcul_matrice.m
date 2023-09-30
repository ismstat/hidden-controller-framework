% on calcule le couple (A,B) et K
function [A,B,K] = f_calcul_matrice(y)
load tmp.mat m n % scaleDead
na = n*n;
nb = n*m;
A=reshape(y(1:na),n,n);
B=reshape(y((na+1):(na+nb)),n,m);
Ru = diag(y((end - 5):(end -3)));
Qx = diag(y((end - 2):end));
[K,PP]=LQ_regulator(A,B,Qx,Ru);
end
