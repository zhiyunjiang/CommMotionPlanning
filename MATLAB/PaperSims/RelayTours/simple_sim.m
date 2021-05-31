%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = [2,3,6];
D =[[0,1,2];[1,0,1];[2, 1, 0]];%for now, we assume this is symmetric
v = 1;
m = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup simple probabilsitic policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simple policy, make transition proportional to data rates
%p = [[0,1/4, 2/5];[1/3, 0, 3/5];[2/3, 3/4,0]];
% Example transition porbabilities for gated servicing strategy.
% Here, probabilites are simply proportional to Poisson rates
P = [[2/11, 2/11, 2/11];[3/11, 3/11, 3/11];[6/11, 6/11,6/11]];

%%%%%%%%%%%%%%%%%%%%%%
% Compute Mean Values
%%%%%%%%%%%%%%%%%%%%%%
%Compute Tij's
[A, b] = build_A_matrix(P, D/v, lambda, m);
T_calc = reshape(A\b,3,3)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate Policy to Validate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iterations = 50000;
data = zeros(iterations+1, 3);
states = -1*ones(iterations+1, 1);
service_times = ones(iterations, 1);
travel_times = ones(iterations, 1);
states(1) = 1;
for i=1:iterations
    states(i+1) = transition(states(i),P);
    travel_times(i) = D(states(i+1), states(i))/v;
    data(i+1,:) = data(i,:)+poissrnd(lambda*travel_times(i));
    
    %%%%%%%%%%%%%%%%%%%%
    % Service to Existing (Gated)
    %%%%%%%%%%%%%%%%%%%%
    service_times(i) = data(i+1,states(i+1))/m;
    data(i+1,states(i+1)) = 0;     
    data(i+1,:) = data(i+1,:) + poissrnd(lambda*service_times(i));
    
    %%%%%%%%%%%%%%%%%%%%
    % Service to Etmpty
    %%%%%%%%%%%%%%%%%%%%
%     service_times(i) = data(i+1,states(i+1))/(m-lambda(states(i)));
%     data(i+1,:) = data(i+1,:)+lambda*service_times(i);
%     data(i+1,states(i+1)) = 0;
end
%%
data_times = [0; cumsum(service_times) + cumsum(travel_times)];
plot(data_times, data)
legend
%%
T_emp = calc_T_emp(states, service_times, travel_times, 3);
[R,S,H] = calc_R_S_emp(states, service_times, travel_times, 3);

%%
x = unique(states);
counts = histc(states,x);
freqs = counts/iterations;

function x_next = transition(x, p)
    pdf = p(:,x);
    r_sample = rand();
    
    CDF = 0;
    for x_next = 1:length(pdf)
        CDF = CDF + pdf(x_next);
        if r_sample <= CDF
            break;
        end
    end
    
end

function [A,b] = build_A_matrix(P, D, lambda, m)
    % Setup as system of linear euqations
    % 3 states - 9 Tij variables
    % Tij -> T((i-1)*3 + j)
    % first calculate the average hop times from each node/region
    t_hop_avg_scaled = diag(D*P);
    n = size(P,1);
    scale2 = lambda/m;

    A = zeros(9,9);
    b = zeros(9,1);
    for i = 1:n
        for j = 1:n
           eq_num = (i-1)*n + j;

           b(eq_num) = -t_hop_avg_scaled(i);

           for q=1:n
               if q~=j
                   coef = P(q, i);%probability of moving to q from i
                   qj = (q-1)*n + j;
                   A(eq_num,qj) =  coef;
               end
           end
           ii = (i-1)*n + i;
           A(eq_num,ii) =  A(eq_num,ii)+scale2(i);
           %lastly, subtract what was on the left hand side
           A(eq_num,eq_num) = A(eq_num,eq_num) - 1;
        end
        
    end
end

function [R,S,H] = calc_R_S_emp(states, service_times, travel_times, n)
    R = zeros(n,1); S = zeros(n,1); H = zeros(n,1);
    for i=1:n
       visits = find(states == i);
       intervals = zeros(length(visits)-1, 1);
       for j = 1:length(visits)-1
           q = visits(j);
           q_n = visits(j+1);
           if q > 1
                intervals(j) = sum(service_times(q-1:q_n-2)) + sum(travel_times(q:q_n-1));
           end
       end
       R(i) = mean(intervals);
     
       
       if visits(1) == 1
           visits(1)=[]; %it's ok to throw one away for hop times, too
       end

       S(i) = mean(service_times(visits -1) );
      
       if visits(end) == length(states)
           visits(end)=[];
       end
       H(i) = mean(travel_times(visits) );
    end
end

function T = calc_T_emp(states, service_times, travel_times, n)
    T = zeros(n,n);
    T_state = inf*ones(n,n);
    T_counts = zeros(n,n);
    T_partial_sums = zeros(n,n);
    for i=1:length(states)
        state = states(i);
        if i~=1 && i~=length(states)
            %everything that was on is added
            T_partial_sums = T_partial_sums + (service_times(i-1) + travel_times(i))*(T_state == 1);
        end
        
        %if we've just found the ending token
        %if was on and is now off...
        % increment count
        %TODO some checking here
        T_counts(:,state) = T_counts(:,state)+ 1*(T_state(:,state) == 1);
        % add to total count
        T(:,state) = T(:,state) + T_partial_sums(:,state);
        T_partial_sums(:,state) = 0;
        T_state(:,state) = 0;%if we've detected the end token, stop
        
        T_state(state,:) = 1;%if we've detected the start token, start counting everything going forward
    end
    
    T = T./T_counts;
    
end


