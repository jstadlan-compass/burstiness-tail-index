function bti = get_BTI_PLEC(alpha,lambda,xmin,phi,plo)

    % can't handle poles–
    if(mod(alpha, 1) == 0)
        alpha = alpha+1E-6;
    end

    x = plcutoff_dist(alpha, lambda, xmin);

    if ~exist("phi","var")
        plo = .75;
        phi = .99;
    end
    p_ref = 0.5;

    Qlo_EXP = -log(1-plo);
    Qhi_EXP = -log(1-phi);
    Qmed_EXP = log(2);
    
    ref_quantiles_lo_hi = [Qlo_EXP Qhi_EXP];
    RATIO_EXP = (Qhi_EXP-Qmed_EXP)/(Qlo_EXP-Qmed_EXP);
    
    ratio_data = (x.ppf(phi)-x.ppf(p_ref))./max(x.ppf(plo)-x.ppf(p_ref),1E-100);
    
    % the thicker the tail, the higher the ratio of 99:lo
    r = ratio_data / RATIO_EXP;
    
    % map to [-1,1]:
    bti = (r-1)./(r+1);

end