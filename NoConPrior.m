function [p_no_connection] = NoConPrior(ric_MP, sigma_SH, gamma_PL, gamma_TH)
    %use linear (not dB) power when working with MP
    J = @(w) normpdf(w, 0, sigma_SH).*cdf(ric_MP, toLin(gamma_TH - gamma_PL - w));
    
    %denominator to be used in a number of calculations
    p_no_connection = integral(J, -Inf, Inf);
end

function lin = toLin(val_db)
    lin = 10.^(val_db/10);
end