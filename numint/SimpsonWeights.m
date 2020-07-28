function w = SimpsonWeights(n, h)
    if n == 1
        %only a single interval, just average the two values, then multiply
        %by interval length
        w = ones([1,2]) * h(1)/2;
    else
   
        w = zeros([1, n+1]);
        is_odd = 0;
        if mod(n,2) ~= 0
            %we're odd
            is_odd = 1;
            %make it even, use odd flag to tack on last part
            n = n-1;
        end

        for i = 1:n/2
            ai = ( 2*h(2*i)^3 - h(2*i-1)^3 + 3*h(2*i)^2*h(2*i-1) )/...
                (6*h(2*i)*(h(2*i) + h(2*i - 1)));
            bi = ( h(2*i)^3 + h(2*i -1)^3 + 3*h(2*i)*h(2*i-1)*(h(2*i) + h(2*i-1)) )/...
                (6*h(2*i)*h(2*i-1));
            etai = (2*h(2*i - 1)^3 - h(2*i)^3 + 3*h(2*i) * h(2*i - 1)^2)/(6*h(2*i - 1)*(h(2*i) + h(2*i - 1) ));

            w(2*i - 1) = w(2*i - 1) + etai;
            w(2*i) = w(2*i) + bi;
            w(2*i+1) = w(2*i + 1) + ai;

        end

        %if odd, tack on the last bit
        if is_odd
           a = (2*h(end-1)^2 + 3*h(end-1)*h(end-2))/(6*(h(end-1) + h(end-2)));
           b = ( h(end-1)^2 + 3*h(end-1)*h(end-2) )/( 6*h(end-2) );
           eta = h(end - 1)^3/( 6*h(end-2)*(h(end-1) + h(end-2)) );

           w(end) = w(end) + a;
           w(end-1) = w(end-1) + b;
           w(end-2) = w(end-2) - eta;
        end
    end
    
end


