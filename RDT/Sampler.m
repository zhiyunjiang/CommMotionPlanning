classdef Sampler < handle
    
    properties (Constant)
        UNIRAND_DISCRETE = 1;
        DETERMINISTIC = 2;
        INFORMED = 3;
        UNIRAND_CONT = 4;
    end
    
    properties (Access = private)
        samplerType;
        destFreq;
        sequence;
        ranges;
        offsets;
        
        C_inf;
        L_inf;
        best_inf;
        
        center;
        destCount;
        
        halfCount;
    end
    
    methods (Access = public)
        function this = Sampler(sampler_type, dest_frequency, sequence)
            this.samplerType = sampler_type;
            if sampler_type == Sampler.INFORMED
               this.best_inf = Inf; 
            end
            this.destFreq = dest_frequency;
            this.sequence = sequence;
            
            this.halfCount = 0;
        end
        
        function setParams(this, pppi)
            gridRegion = pppi.getGridRegion();
            x_range = gridRegion(1) - gridRegion(2) + 1;
            y_range = gridRegion(3) - gridRegion(4) + 1;
            
            this.ranges = [x_range, y_range];
            this.offsets = [gridRegion(2), gridRegion(4)];
            
            goalRegion = pppi.getGoalRegion();
            dests = goalRegion.goalGridPoints;
            [this.destCount, ~] = size(dests);
            
            if this.destCount == 1
                this.center = (pppi.getSourceGrid() + dests)/2;
            end
        end
        
        function setHalfCount(this, hc)
           this.halfCount = hc; 
        end
        
        function setDestFreq(this, freq)
           this.destFreq = freq; 
        end
        
        function freq = getDestFreq(this)
           freq = this.destFreq;
        end
        
        function new_sample = sample(this, count, pppi, bsf)
            if this.halfCount
               count = floor(count/2); 
            end
            if mod(count, this.destFreq)==0 && ...
                    (this.samplerType == Sampler.UNIRAND_DISCRETE || this.samplerType == Sampler.UNIRAND_CONT)
                
                new_sample = this.sampleDestination(pppi);
                
            elseif this.samplerType == Sampler.UNIRAND_DISCRETE
                
                new_sample = this.sampleUniRandDiscrete();
                
            elseif this.samplerType == Sampler.DETERMINISTIC
                
                new_sample = this.deterministicSequence(count);
                
            elseif this.samplerType == Sampler.INFORMED
                
                new_sample = this.sampleInformedSet(this, count, pppi, bsf);
                
            elseif this.samplerType == Sampler.UNIRAND_CONT
                new_sample = this.sampleUniRandCont();
            end
            this.sequence(count + 1,1:2) = new_sample;
            
            %sampling always done in (1,1) shifted space, possibly
            %discretized. Map back to a point in true coordinates
            new_sample = pppi.toRawCoordinate(new_sample);

        end
    end
    
    methods (Access = private)
        
        function sample = sampleUniRandDiscrete(this)
            x_rand = randi(this.ranges(1),1) + this.offsets(1) - 1; 
            y_rand = randi(this.ranges(2),1) + this.offsets(2) - 1; 
            sample = [x_rand, y_rand];
        end
        
        function sample = sampleUniRandCont(this)
            %account for different behavior between rand and randi
            x_rand = (this.ranges(1)-1)*rand(1) + this.offsets(1); 
            y_rand = (this.ranges(2)-1)*rand(1) + this.offsets(2); 
            sample = [x_rand, y_rand];
        end
        
        function sample = sampleDestination(this, pppi)
             
            dest = pppi.getGoalRegion().goalPoints;   
            dest_pppi_grid = pppi.toGridCoordinate(dest);
            %randomly choose from the destination set
            [dest_count, ~] = size(dest_pppi_grid);
            rand_i = randi(dest_count,1);
            sample = dest_pppi_grid(rand_i,:);
        end
        
        %based on "Informed RRT*Optimal sampling-based path planning 
        % focused via direct sampling of an admissible ellipsoidal 
        % heuristic", Gammell, Srinivasa, et al. 2014 
        function sample = sampleInformedSet(this, count, pppi, bsf)
            if bsf == Inf
                sample = this.randSample(count, pppi);
            else
               %sample from the infromed set
               %if there's just one destination, sample directly from the
               %hyperellipsoid.
                if this.destCount == 1
                    
                    sample = this.sampleSingleDestinationInformed(pppi, bsf);
                    
                else
                    %otherwise, we will do rejection sampling, dependingo n
                    %the number of destinations, this may not in fact help
                    %all that much
                    sample = this.sampleMultiDestinationInformed(pppi, bsf);
                    
                end
            end
        end
        
        function sample = sampleMultiDestinationInformed(this, pppi, bsf)
            
            dests = pppi.getDestGrid();
            src = pppi.getSourceGrid();
            
            keep = 0;
            while ~keep
                sample = this.sampleUniRand();
                for i=1:this.destCount
                    dest = dests(i);
                   if norm(src - sample) + norm(dest - sample) < bsf 
                       %if it's in the ellipsoid of at least one, then we
                       %can keep it
                       keep = 1;
                       break
                   end
                end
            end
            
        end
        
        function sample =  sampleSingleDestinationInformed(this, pppi, bsf)
            
            if bsf < o.best_inf
                  %reset the ellipsoid parameters
                    
                    c_min = norm(pppi.getSourceGrid() + pppi.getDestGrid());
                    this.C_inf = Sampler.RotateToWorldFrame(start,goal, c_min);
               
                    c_max = bsf;
                    this.L_inf = [[c_max/2, 0]; [0, sqrt(c_max^2 - c_min^2)/2]];
            end

           circle_sample = Sampler.SampleFromUnitCircle();
           sample = this.C_inf*this.L_inf*circle_sample + this.center;
           while ~pppi.ptInRegion(sample)
               circle_sample = Sampler.SampleFromUnitCircle();
                sample = this.C_inf*this.L_inf*circle_sample + this.center;
                sample = pppi.toGridCoordinate(sample);
           end
        end
        
        function next = deterministicSequence(o, count)
            index = count + 1;
            if index <= length(o.sequence)
                next = o.sequence(index); 
            else
                error('Trying to sample beyond given sequence. \nFor same problem with the same termination criterion, no additional samples should be needed');
            end
        end

    end
    
    methods (Access = private, Static = true)
        function C = RotateToWorldFrame(start, goal, c_min)
            a1 = (goal-start)/c_min;
         
            M = [a1;0 0]';
            [U,~,V] = svd(M);
            C = U*diag([1, det(U)*dev(V)])*V';
        end
        
        function p = SampleFromUnitCircle()
           theta = 2*pi*rand(1);
           r = sqrt(rand(1));
           p = [r*cos(theta), r*sin(theta)];
        end
    end
end

