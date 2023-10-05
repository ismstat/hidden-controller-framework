%%%%%%%%%%%%%%%%%%%%
% Hidden Controller Framework for (ii) data from nine countries sourced from
%   a "worldwide epidemiological database for COVID-19" Guidotti (2022) and
%   the "COVID-19 data hub" Guidotti and Ardia (2020).
%%%%%%%%%%%%%%%%%%%%
addpath('code');
rng(123);
countrynames = ['AUS'; 'COL'; 'DEU'; 'JPN';'BRA'; 'CHL'; 'CZE'; 'LTU'; 'ZAF'];
for country = 1:9
    cname = countrynames(country,:);
    %%%%%%%%%%%%%%%%%%%%
    % Parameters for Matlab fmincon
    % Find minimum of constrained nonlinear multivariable function
    %%%%%%%%%%%%%%%%%%%%
    scaleY = 1;
    scaleD = 1;
    perturb = randn(1);
    ymin = 1e-4;
    options=optimoptions('fmincon');
    options.Algorithm = 'interior-point'; % 'sqp'
    options.Display = 'iter';
    options.OptimalityTolerance = 1e-4;
    options.FunctionTolerance = 1e-4;
    options.StepTolerance = 1e-4;
    options.UseParallel=true;
    options.MaxFunctionEvaluations = 1.000000e+04;
    options.MaxIterations = 1e+04;
    %%%%%%%%%%%%%%%%%%%%
    % Read data
    %%%%%%%%%%%%%%%%%%%%
    [MO, nsteps, npref, n, m] = readdata_pop_country([cname, '.csv']);
    xtype = 'original';
    %%%%%%%%%%%%%%%%%%%
    % For a goal of reducing the counts in both categories, e.g., by half.
    %%%%%%%%%%%%%%%%%%%
    xinfected = 1;
    xrecovered = 1;
    xdead = 1;
    % xtype = 'rate'
    % xinfected = 0.5;
    % xrecovered = 1;
    % xdead = 0.5;
    %%%%%%%%%%%%%%%%%%%%
    % Output
    %%%%%%%%%%%%%%%%%%%%
    logID = 1;
    odir = append('output', '_', options.Algorithm, '_', xtype, '_i', num2str(xinfected), '_d', num2str(xdead),'_', cname); 
    [status, msg, msgID] = mkdir(odir);
    %%%%%%%%%%%%%%%%%%%
    % Target data
    %%%%%%%%%%%%%%%%%%%
    M = target_x(MO, nsteps, npref, n, xtype, xinfected, xrecovered, xdead, odir);
    %%%%%%%%%%%%%%%%%%%%
    % fast_mpc: code for fast model predictive control
    % https://stanford.edu/~boyd/fast_mpc/
    %%%%%%%%%%%%%%%%%%%%
    % x(t+1) = A x(t) + B u(t) + G w(t)
    % y(t+1) = x(t+1)  or  y(t+1) = C x(t) + v(t)
    % x: m x nsteps
    % A: m x m
    % B: n x m
    % u: n x nsteps
    prefset = [1];
    m = m - 1;
    n = n - 1;
    save M.mat M
    nsteps1 = nsteps;
    T = 2;

    i = 1;
    for np = prefset
        fprintf(logID, 'prefecture: %i\n', np);
        yoptimal = rand(n*n+n*m+3,1) * scaleY;
        save yopt.mat yoptimal;
        for it = 1:1
            load yopt.mat yoptimal;
            save tmp.mat m n np T nsteps1 i scaleD
            if (abs(perturb) > 1e-1)
                yo = yoptimal + rand(n*n+n*m+3,1)/10;
            else
                yo = yoptimal;
            end
            nonlinco=@f_de_contrainte_k;
            [yoptimal,fval] = fmincon(@f_objective_optimiser_all,yo,[],[],[],[],[],[],nonlinco,options);
            save yopt.mat yoptimal
            [A,B,K] = f_calcul_matrice(yoptimal);
            nstep = floor(nsteps1 / T);
            for t = 1:nstep
                xx=M(1:n,np,(t-1)*T+1); % initialization
                x(:,(t-1)*T+1)=xx;
                for i=((t-1)*T+2):((t-1)*T+T)
                    x(:,i)=(A+B*K)*x(:,i-1);
                    disp(x(:,i-1))
                    disp(xx) % prediction
                end
            end
            % disp(A)
            % disp(B)
            % disp(K)
            % disp(yoptimal(end-2))
            % disp(yoptimal(end-1))
            % disp(yoptimal(end))
            Aname = append(odir, '/A', int2str(it), '_', int2str(np), '.csv');
            csvwrite(Aname, A);
            Bname = append(odir, '/B', int2str(it), '_', int2str(np), '.csv');
            csvwrite(Bname, B);
            Kname = append(odir, '/K', int2str(it), '_', int2str(np), '.csv');
            csvwrite(Kname, K);
            SRname = append(odir, '/SR', int2str(it), '_',int2str(np), '.csv');
            csvwrite(SRname, yoptimal(end-2));
            Qname = append(odir, '/Q', int2str(it), '_', int2str(np), '.csv');
            csvwrite(Qname, yoptimal(end-1));
            Rname = append(odir, '/R', int2str(it), '_', int2str(np), '.csv');
            csvwrite(Rname, yoptimal(end));
        end
    end
end




