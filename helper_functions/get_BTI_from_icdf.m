function bti = get_BTI_from_icdf(icdf,phi,plo)

    if ~exist("phi","var")
        plo = .75;
        phi = .99;
    end
    p_ref = 0.5;

    Qlo_EXP = -log(1-plo);
    Qhi_EXP = -log(1-phi);
    Qmed_EXP = log(2);
    
    RATIO_EXP = (Qhi_EXP-Qmed_EXP)/(Qlo_EXP-Qmed_EXP);
    
    ratio_data = (icdf(phi)-icdf(p_ref))./max(icdf(plo)-icdf(p_ref),1E-100);
    
    % the thicker the tail, the higher the ratio of 99:lo
    r = ratio_data / RATIO_EXP;
    
    % map to [-1,1]:
    bti = (r-1)./(r+1);

end