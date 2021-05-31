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
er1 = rslt - (1/4)

n = 5;
h = (1/n)*ones([n,1]);
x = [0;cumsum(h)];
y = x.^3;
%integral from 0 to 1 of x^3  = 1/4 (analytically)
rslt = IntegrateWSimpsonsIrregular(y,h);
er2 = rslt - (1/4)


%Test 2 - Unevenly spaced data
%even number
h = [0.2;0.17;0.33;0.3];
x = [0;cumsum(h)];
y = x.^3;
%integral from 0 to 1 of x^3  = 1/4 (analytically)
rslt = IntegrateWSimpsonsIrregular(y,h);
er3 = rslt - (1/4)

%odd number
h = [0.17;0.21;0.33;0.2; 0.09];
x = [0;cumsum(h)];
y = x.^3;
%integral from 0 to 1 of x^3  = 1/4 (analytically)
rslt = IntegrateWSimpsonsIrregular(y,h);
er4 = rslt - (1/4)