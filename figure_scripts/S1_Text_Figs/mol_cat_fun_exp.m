
function [k_p, k_c, kl,kh_p,kh_c] = mol_cat_fun_exp(alpha,gamma,epsilon,hp_score,l)
 
        kh_p = exp(hp_score*alpha)/1e3;%(exp((hp_score*k_hs)))/1e3;%
        kh_c = exp(-(hp_score*alpha))/1e3;%(exp((-hp_score*k_hs)))/1e3;%
        kl = (1/(1+exp(-l+gamma)))*epsilon; %sigmoid function, gamma shifts function horizontally and epsilon changes the plateau value.
        k_p = kh_p * kl;
        k_c = kh_c * kl;

end



