function bp = get_PLEC_BP(alpha,lambda,xmin)

    if xmin > 0
        x = plcutoff_dist(alpha, lambda, xmin);
        r = x.cv;
    else
        r = gamma(3-alpha).*gamma(1-alpha)./(gamma(2-alpha).^2) - 1;

    end
    bp = (r-1)./(r+1);
end