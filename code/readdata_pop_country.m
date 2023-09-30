function [MW, nsteps, npref, n, m] = readdata_pop_country(fname)

x = readtable(fname);
dead = x(:,4);
infected = x(:,3);
recover = x(:,5);
pop = x(:,37);

N = pop.population(1); 
npref = 1;

mindate = datenum('2021-4-1');
maxdate = datenum('2021-11-22');

datetime(mindate, 'ConvertFrom','datenum')
datetime(maxdate, 'ConvertFrom','datenum')

nday = maxdate - mindate + 1;
days = mindate:maxdate;


nnum = 4;  
nu = nnum; 

M = zeros(nnum,npref,nday);
nwin = 7;
nweek = floor(nday / nwin);
MW = zeros(nnum, npref, nweek);

[ndead, dummy] = size(dead);
[ninfected, dummy] = size(infected);
[nrecover, dummy] = size(recover);

if ndead ~= nday + 1
    msg = ['Warning: ndead != nday + 1'];
    disp(msg);
end
if ninfected ~= nday + 1
    msg = ['Warning: ninfected != nday + 1'];
    disp(msg);
end
if nrecover ~= nday + 1
    msg = ['Warning: nrecover != nday + 1'];
    disp(msg);
end


M(1,1,:) = diff(infected.confirmed);
M(2,1,:) = diff(dead.deaths);
M(3,1,:) = diff(recover.recovered);

for i = 1:1
    for j = 1:nweek
        for n = 1:(nnum - 1)
            MW(n, i, j) = sum(M(n, i, ((j - 1) * nwin + 1):(j * nwin)));
            if MW(n, i, j) < 0
                msg = ['Warning: MW ', num2str(n), ', ', num2str(i), ', ', num2str(j), ' is ', num2str(MW(n, i, j)), ' .'];
                disp(msg);
                MW(n, i, j) = 0;
            end
            MW(n, i, j) = MW(n, i, j) / N;
        end
    end
end


disp(std(MW(1,1,1:nweek)));
disp(std(MW(2, 1,1:nweek)));
disp(std(MW(3, 1,1:nweek)));

MW(1, 1, :) = (MW(1, 1, :) -min(MW(1, 1, :))) / max(MW(1, 1, :));
MW(2, 1, :) = (MW(2, 1, :) -min(MW(2, 1, :))) / max(MW(2, 1, :));
MW(3, 1, :) = (MW(3, 1, :) -min(MW(3, 1, :))) / max(MW(3, 1, :));

maxsum = max(MW(1,1,:)+MW(2,1,:)+ MW(3,1,:)) + 1;

MW(1, 1, :) = MW(1, 1, :) / maxsum;
MW(2, 1, :) = MW(2, 1, :) / maxsum;
MW(3, 1, :) = MW(3, 1, :) / maxsum;

n = nnum;
m = nu;
nsteps = nweek;




