% les contraintes de l'optimisation
function [c,ceq] = f_de_contrainte_k(y)
load M.mat M;
load tmp.mat np T n nsteps1 %scaleD
[A,B,K]= f_calcul_matrice(y);
nstep = floor(nsteps1 / T);
c = [];
for i = 1:nstep
    %c=[M(1:n,np,i)-ones(n,1)*1000;-M(1:n,np,i)-zeros(n,1)];
    %c=[M(1:n,np,i)-ones(n,1);-M(1:n,np,i)-zeros(n,1)];
    j = (i-1)*T+1;
    for t=2:T   % number of days in total
        % c=[c;(A+B*K)^(t-1)*M(1:n,np,i)-ones(3,1)*1000;-(A+B*K)^(t-1)*M(1:n,np,i)-zeros(n,1)];
        c=[c;(A+B*K)^(t-1)*M(1:n,np,j)-ones(n,1);-(A+B*K)^(t-1)*M(1:n,np,j)-zeros(n,1)];
    end
end
c=[c;-y((end-5):end)];
ceq=[];
end 

