classdef PollSys < handle
    %PollSys
    
    properties
        n; % # of 
        lambdas; % Time averaged arrival rate at queue
        beta; % average service time
        S; %matrix of switching times, with S(i,j) time form i to j
        rhos; % traffic at each queue
        
    end
    
    methods 
        function this = PollSys(n, lambdas, beta, S)
            this.n = n;
            this.lambdas = lambdas;
            this. beta = beta;
            this.S = S;
            this.rhos = lambdas*beta;
            
            if this.rhos > 1
                warning("Unstable (non-ergodic) polling system created. System traffic is greater than 1")
            end
        end
        
        function freqs = findOptiamlVisitFreqs(this)
            x0 = ones([this.n, 1])/this.n;%To begin, assume equal frequency
            f = @(x) calcAvgWait(this, x);
            Aeq = ones([1, this.n]); Beq = 1;%constrain frequencies to be in simplex
            A = -1*eye(this.n); B = zeros([this.n, 1]);
            [freqs, ~] = fmincon(f, x0, A, B, Aeq, Beq);
        end
        
        % Based on: Boxma, "Efficient Polling Orders for Polling Systems"
        % 1993
        function W = calcAvgWait(this, freqs)
            Sis = this.S*freqs;
            Si2s = (this.S.^2)*freqs;
            rho = this.rhos';
            rhoS = this.rhoS();
            W = (1/(this.lambdaS()*this.beta))*(...
                (rhoS*this.lambdaS*this.beta^2)/(2*(1-rhoS)) +...
                (freqs'*Sis/(1-rhoS))*(rho-rho.^2)'*(1./freqs) - ...
                rho'*Sis + rhoS*(freqs'*Si2s)/(2*freqs'*Sis));
        end
        
        function ls = lambdaS(this)
            ls = sum(this.lambdas);
        end
        
        function rs = rhoS(this)
            rs = sum(this.rhos); 
        end
        
        function plotWvsFreq(this)
           if this.n >3
               warning("This function is only valid for 2- ro 3-queue systems.\n");
           elseif this.n == 3
              f1 = 0.001:0.001:1;
              f2 = 0.001:0.001:1;
              d = length(f1);
              Z = -1*ones(d);
              for i=1:d
                  for j=1:d
                    if f1(i)+f2(j) < 1
                        Z(j,i) = this.calcAvgWait([f1(i);f2(j); 1 - f1(i)-f2(j)]);
                    end
                  end
              end
              %h = surf(f1, f2, Z);
              %set(h,'LineStyle','none')
              h = imagesc(f1, f2, Z);
              
           elseif this.n == 2 
              % create simplex
              f1 = 0.001:0.001:1;
              d = length(f1);
              Z = zeros([d,1]);
              
              for i=1:d
                    Z(i) = this.calcAvgWait([f1(i);1 - f1(i)]);
              end
              
               plot(f1, Z);
               
           end
        end
        
            
    end
end

