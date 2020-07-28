%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tests of IntegrateWSimpsonsIrregular
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Test 1 - 
%test simple cubic function, regularly spaced.
%Both even and odd number of interval
h = [0.25; 0.25; 0.25; 0.25];
x = [0;cumsum(h)];
y = x.^3;
%integral from 0 to 1 of x^3  = 1/4 (analytically)
rslt = IntegrateWSimpsonsIrregular(y,h);
er = rslt - (1/4)

h = (1/5)*ones([5,1]);
x = [0;cumsum(h)];
y = x.^3;
%integral from 0 to 1 of x^3  = 1/3 (analytically)
rslt = IntegrateWSimpsonsIrregular(y,h);
er = rslt - (1/4)


%Test 2 - Unevenly spaced data
h = [0.2;0.17;0.33;0.3];
x = [0;cumsum(h)];
y = x.^3;
%integral from 0 to 1 of x^3  = 1/3 (analytically)
rslt = IntegrateWSimpsonsIrregular(y,h);
er = rslt - (1/4)

h = [0.17;0.21;0.33;0.2; 0.09];
x = [0;cumsum(h)];
y = x.^3;
%integral from 0 to 1 of x^3  = 1/3 (analytically)
rslt = IntegrateWSimpsonsIrregular(y,h);
er = rslt - (1/4)