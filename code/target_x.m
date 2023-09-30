function MTarget = target_x(M, nday, npref, nnum, xtype, xinfected, xrecovered, xdead, odir)
%%%%%%%
% X(:,1) =
% [M(1,1,1),M(2,1,1),M(3,1,1),M(4,1,1),M(1,2,1),M(2,2,1),...,M(nnum,npref,1),
%  M(1,1,2),M(2,1,2),...,M(nnum,npref,2),
%  ...
%  M(1,1,nday),...,M(nnum,npref,nday)]
%%%%%%%
MTarget = zeros(nnum, npref, nday);
for np = 1:npref
    M0 = zeros(nnum, nday);
    MT = zeros(nnum, nday);
    switch(xtype)
        case 'rate'
            nn = 1;
            for nd = 1:nday
                M0(nn, nd) = M(nn,np,nd);
                MT(nn, nd) = M(nn,np,nd) * xinfected;
                MTarget(nn, np, nd) = MT(nn, nd);
                
            end
            nn = 2;
            for nd = 1:nday
                M0(nn, nd) = M(nn,np,nd);
                MT(nn, nd) = M(nn,np,nd) * xrecovered;
                MTarget(nn, np, nd) = MT(nn, nd);
            end
            nn = 3;
            for nd = 1:nday
                M0(nn, nd) = M(nn,np,nd);
                MT(nn, nd) = M(nn,np,nd) * xdead;
                MTarget(nn, np, nd) = MT(nn, nd);
            end
            nn = 4;
            for nd = 1:nday
                M0(nn, nd) = M(nn,np,nd);
                MT(nn, nd) = 1 - MT(1, nd) - MT(2, nd) - MT(3, nd);
                MTarget(nn, np, nd) = MT(nn, nd);
            end
        case 'epsilon'
            nn = 1;
            for nd = 1:nday
                M0(nn, nd) = M(nn,np,nd);
                if M(nn,np,nd) > xinfected
                    MT(nn, nd) = xinfected;
                else
                    MT(nn, nd) = M(nn, np, nd);
                end
                MTarget(nn, np, nd) = MT(nn, nd);
            end
            nn = 2;
            for nd = 1:nday
                M0(nn, nd) = M(nn,np,nd);
                MT(nn, nd) = M(nn, np, nd);
                MTarget(nn, np, nd) = MT(nn, nd);
            end
            nn = 3;
            for nd = 1:nday
                M0(nn, nd) = M(nn,np,nd);
                if M(nn,np,nd) > xdead
                    MT(nn, nd) = xdead;
                else
                    MT(nn, nd) = M(nn, np, nd);
                end
                MTarget(nn, np, nd) = MT(nn, nd);
            end
            nn = 4;
            for nd = 1:nday
                M0(nn, nd) = M(nn,np,nd);
                MT(nn, nd) = 1 - MT(1, nd) - MT(2, nd) - MT(3, nd);
                MTarget(nn, np, nd) = MT(nn, nd);
            end
        otherwise
            for nn = 1:nnum
                for nd = 1:nday
                    M0(nn, nd) = M(nn,np,nd);
                    MT(nn, nd) = M(nn,np,nd);
                    MTarget(nn, np, nd) = M(nn, np, nd);
                end
            end
    end
    fname = append(odir, '/MO_',int2str(np), '.csv');
    csvwrite(fname,M0);
    fname = append(odir, '/MT_',int2str(np), '.csv');
    csvwrite(fname,MT);
end
end
