function [remapped_ratio_data_to_ref, quantiles_lo_hi, ref_quantiles_lo_hi] = get_sample_BTI(data)
    pmedian = .5;
    plo = .75;
    phi = .99;

    quantiles_lo_hi = quantile(data, [plo phi pmedian]);
    Qlo_EXP = -log(1-plo);
    Qhi_EXP = -log(1-phi);
    Qmed_EXP = log(2);

    ref_quantiles_lo_hi = [Qlo_EXP Qhi_EXP];
    RATIO_EXP = (Qhi_EXP-Qmed_EXP)/(Qlo_EXP-Qmed_EXP);

    if isvector(quantiles_lo_hi)
        ratio_data = (quantiles_lo_hi(2)-quantiles_lo_hi(3))./max(quantiles_lo_hi(1)-quantiles_lo_hi(3),1E-100);
    else
        ratio_data = (quantiles_lo_hi(2,:)-quantiles_lo_hi(3,:))./max(quantiles_lo_hi(1,:)-quantiles_lo_hi(3,:),1E-100);    
    end
    % the thicker the tail, the higher the ratio of 99:lo
    r = ratio_data / RATIO_EXP;

    % map to [-1,1]:
    remapped_ratio_data_to_ref = (r-1)./(r+1);
end

