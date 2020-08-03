%Unit Tests for RDTree
tcount = 0; fcount = 0;
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%Testing RDTree.getApproxGridSLPath
%%%%%%%%%%%%%%%%%%%%%%%%
%Straight up
start = [0, 0.5]; dest = [0, 2.5];
expected_path = [start; [0, 1.5]; dest];
[tcount, fcount, success] = testGetApproxGridSLPath(start, dest, expected_path, tcount, fcount);
if ~success
    fprintf('failed straight up test');
end

%straight left
start = [-1, -2]; dest = [-2, -2];
expected_path = [start; dest];
[tcount, fcount, success] = testGetApproxGridSLPath(start, dest, expected_path, tcount, fcount);
if ~success
    fprintf('failed left test');
end

%right two, down three
start = [5, 1]; dest = [7, -2];
expected_path = [start; [6, 0]; [6, -1]; dest];
[tcount, fcount, success] = testGetApproxGridSLPath(start, dest, expected_path, tcount, fcount);
if ~success
    fprintf('failed right two, down three test');
end
%start and destination the same
start = [5, 1];
[tcount, fcount, success] = testGetApproxGridSLPath(start, start, start, tcount, fcount);
if ~success
    fprintf('failed same start, destination test');
end




fprintf('Failed %d of %d tests\n', fcount, tcount);

function [tcount, fcount, success] = testGetApproxGridSLPath(start, dest, expected_path, tcount, fcount);
    path = RDTree.getApproxGridSLPath(start, dest);
    success = sum(abs(path - expected_path), 'all') == 0;
    
    if ~success
        fcount = fcount + 1;
    end
    
    tcount = tcount + 1;
end