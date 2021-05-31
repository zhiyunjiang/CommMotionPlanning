function [G_v, G_v_pl, sigma_v] = estimate_channel(x_s, y_s, z_s, x_v, y_v, alpha_h, beta_h, eta_h, K_pl_h, n_pl_h, q_b)

x_s = x_s(:);
y_s = y_s(:); 
z_s = z_s(:); 

N = length(z_s);

eps = 0.1;
H = [ones(N,1) -10*log10(sqrt((x_s(:) - q_b(1)).^2 + (y_s - q_b(2)).^2 + eps^2))];

U_h = zeros(N);

for i=1:N
    U_h(i,i) = alpha_h + eta_h;
    for j = i+1:N
        U_h(i,j) = alpha_h*exp(- sqrt((x_s(i) - x_s(j))^2 +  (y_s(i) - y_s(j))^2) /beta_h);
        U_h(j,i) = U_h(i,j);
    end;
end;

U_h_inv = U_h^-1;

theta_h = [K_pl_h; n_pl_h];

x_v = x_v(:);
y_v = y_v(:);
             
E_x = repmat(x_v(:), 1, N) - repmat(x_s(:)', length(x_v), 1);
E_y = repmat(y_v(:), 1, N) - repmat(y_s(:)', length(x_v), 1);

phi = alpha_h*exp(-sqrt(E_x.^2 + E_y.^2)/beta_h);
h = [ones(length(x_v),1) -10*log10(sqrt((x_v(:) - q_b(1)).^2 + (y_v(:) - q_b(2)).^2 + eps^2))];

G_v = h*theta_h + phi*U_h_inv*( z_s(:) - H*theta_h);
G_v_pl = h*theta_h;

sigma_v =  sqrt(alpha_h + eta_h - diag(phi*U_h_inv*phi'));

if ~isreal(G_v) || ~isreal(sigma_v)
    error('')
end