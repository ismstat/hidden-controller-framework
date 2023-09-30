%--------------------------------------------------------------------------
% Function for Linear-Quadratic Regulator (LQR)
% calculates the optimal gain matrix K for a driscrete time system.
% the state-feedback law u = Kx minimizes the quadratic cost function:
%        k=inf
% J =     SUM    x(k)'*Q*x(k) + u(k)'*R*u(k)
%         k=1
%--------------------------------------------------------------------------
%  Author: Edwin Alonso González Querubín
%          https://www.researchgate.net/profile/Edwin_Gonzalez_Querubin
%          https://es.mathworks.com/matlabcentral/profile/authors/15149689
%  Research Group: Predictive Control and Heuristic Optimization (CPOH)
%                  http://cpoh.upv.es
%  Unversity: Universidad Politécnica de Valencia
%             http://www.upv.es
%  Version: Beta
%  For new releases and bug fixing of this tool please visit:
%  https://es.mathworks.com/matlabcentral/fileexchange/75803
%--------------------------------------------------------------------------

function [K, P] = LQ_regulator(A ,B, Q, R)

%--------------------------------------------------------------------------
Pprev = ones(size(Q));
P = zeros(size(Q));
% Pprev1 = ones(size(Q));
flag = 1;
iflag = 0;
c = 0;
while (P ~= Pprev) % & (P ~= Pprev1)
    % Pprev1 = Pprev;
    Pprev = P;
    for i=1:100
        if (det(R + B'*P*B) ~= 0)
            P = Q + A'*P*A - A'*P*B*inv(R + B'*P*B)*B'*P*A;
            iflag = 1;
        end
    end
    % if (eigs((A - B*inv(R)*B'*P), 1) > 0)
    if (eigs((A - B*inv(R + B'*P*B)*B'*P*A), 1) > 1)
        flag = 0;
        break;
    end
    c = c + 1;
    if (c > 10000)
        disp('aho');
        break;
    end
end
if flag && iflag
    K = -inv(R + B'*P*B)*B'*P*A;
else
    K = NaN;
end