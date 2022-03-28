f1 = datetime(1948,1,1,0,0,0);
f2 = datetime(2015,12,1,0,0,0);
ft = f1:calmonths(1):f2;

d1 = datetime(1955,11,1,0,0,0);
d2 = datetime(1955+10,12,1,0,0,0);
dt = d1:calmonths(1):d2;

fosi_date = ft;
dple_date = dt;

idst = find(fosi_date==dple_date(1))
iden = find(fosi_date==dple_date(end))

