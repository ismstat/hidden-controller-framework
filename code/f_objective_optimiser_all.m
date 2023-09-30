% fonction optimiser
function J = f_objective_optimiser_all(y)
load M.mat M;
load tmp.mat np T n nsteps1


%%% This implementation may have a problem.  {true, estimate, estimate}*
nstep = floor(nsteps1 / T);
for t = 1:nstep
    xx=M(1:n,np,(t-1)*T+1);   % initialization
    
    [A,B,K] = f_calcul_matrice(y);
    
    x(:,(t-1)*T+1)=xx;
    for i=((t-1)*T+2):((t-1)*T+T)
        x(:,i)=(A+B*K)*x(:,i-1);   % prediction
    end
end

for i = 1:(nstep*T)
    z(:,i)=norm(x(1:n,i)-M(1:n,np,i),2);  % for all sigmals, calculate norm between true observations and predicted ones.
end
J=sum(z);  % sum for all norms ... minimize J in fmincon.
end
