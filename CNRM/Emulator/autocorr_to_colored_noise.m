function [Y_t] = autocorr_to_colored_noise(autocorr,n_seq)
%Given auto correlation, output colored noise

var_t = zeros(n_seq,1);
Y_t = zeros(n_seq,1);
phi = zeros(n_seq,n_seq);

var_t(1) = autocorr(1);
Y_t(1) = sqrt(var_t(1))*randn;
phi(2,2) = autocorr(2) / var_t(1);
var_t(2) = var_t(1) * (1 - phi(2,2)^2);
Y_t(2) = phi(2,2)*Y_t(1) + sqrt(var_t(2)) * randn;
for i=3:n_seq
    phi(i,i) = autocorr(i);
    for j=2:i-1
        phi(i,i) = phi(i,i) - phi(j,i-1)*autocorr(i-j+1);
    end
    phi(i,i) = phi(i,i) / var_t(i-1);
    for j=2:i-1
        phi(j,i) = phi(j,i-1) - phi(i,i) * phi(i-j+1,i-1);
    end
    var_t(i) = var_t(i-1) * (1 - phi(i,i)^2);
    for j=2:i
        Y_t(i) = Y_t(i) + phi(j,i) * Y_t(i-j+1);
    end
    Y_t(i) = Y_t(i) + sqrt(var_t(i))*randn;
end
end

