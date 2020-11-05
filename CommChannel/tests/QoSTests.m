%QoS Tests
%%%%%%%%%%%%%%
%Initilization
%%%%%%%%%%%%%%
num_tests = 0;
num_fails = 0;
%values not necessarily realistic, but easy to verify that computation is
%correct

%%%%%%%%%%%%%%%%%
%Test Constructor
%%%%%%%%%%%%%%%%%
num_tests = num_tests + 1;
BER = exp(-1)/5; r = 0; rx_noise = -20;
did_catch = 0;
try 
    QoSParams(BER, r, rx_noise)
catch 
    did_catch = 1;
end
if ~did_catch
    fprintf('Failed QoSParams constructor test 1\n');
    num_fails = num_fails + 1;
end

num_tests = num_tests + 1;
BER = exp(-1)/5; r = -1; rx_noise = 0;
did_catch = 0;
try 
    QoSParams(BER, r, rx_noise)
catch 
    did_catch = 1;
end
if ~did_catch
    fprintf('Failed QoSParams constructor test 2\n');
    num_fails = num_fails + 1;
end

num_tests = num_tests + 1;
BER = 4; r = 1; rx_noise = 0;
did_catch = 0;
try 
    QoSParams(BER, r, rx_noise)
catch 
    did_catch = 1;
end
if ~did_catch
    fprintf('Failed QoSParams constructor test 3\n');
    num_fails = num_fails + 1;
end

%%%%%%%%%%%%%%%%
%Test reqTXPower
%%%%%%%%%%%%%%%%
num_tests = num_tests + 1;
BER = exp(-1)/5; r = 2; rx_noise = 1e-10;
qp = QoSParams(BER, r, rx_noise);

if 1 ~= qp.reqTXPower(10*log10(rx_noise/500))
    fprintf('Failed QoSParams reqTXPower test 1\n');
    num_fails = num_fails + 1;
end




if num_fails == 0
    fprintf('All tests passed\n');
end