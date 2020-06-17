function [no_con_to_here, J] = NoConToHere(ric_MP, sigma_SH_c, gamma_PL, gamma_TH, J_prev, p_Y0_lessthan_TH)

    J = @(w) (cdf(ric_MP, toLin(gamma_TH - gamma_PL - w)) / ro) .* ...
        arrayfun( @(wk) integral(@(s) normpdf(s, wk,sigma_SH_c) .* J_prev(s/ro), -Inf , Inf), w );

    no_con_to_here = integral(J, -Inf, Inf)/p_Y0_lessthan_TH;

end

