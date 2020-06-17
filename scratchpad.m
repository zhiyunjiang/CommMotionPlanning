
%%
for K_ric = 1:10
    %values
    logx = [-20:20];
    linx = toLin(logx);
    nu_MP = sqrt(K_ric/(K_ric + 1));%non-centrality parameter
    sigma_MP = 1/sqrt(2*(1+K_ric));%spread parameter
    ric_MP = makedist('Rician','s',nu_MP,'sigma',sigma_MP);
    p = cdf(ric_MP,linx);
    plot(logx, p);
    hold on;
end


%%
function lin = toLin(val_db)
    lin = 10.^(val_db/10);
end