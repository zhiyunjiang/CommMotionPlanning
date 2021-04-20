%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
n = 2;
lambdas = [1,3];
beta = 0.2;
S = [[0,1];[1,0]];

ps2 = PollSys(n, lambdas, beta, S);

ps2.plotWvsFreq()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
n = 3;
lambdas = [10,10,10];
beta = 0.005;

S = [[0, 2, 5];[2,0, sqrt(29)]; [5, sqrt(29), 0]];

ps3 = PollSys(n, lambdas, beta, S);

ps3.plotWvsFreq()

0.ps3.findOptiamlVisitFreqs()