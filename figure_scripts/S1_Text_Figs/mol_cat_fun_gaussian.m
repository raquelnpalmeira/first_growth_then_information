function [k_p, k_c, kl,kh_p,kh_c] = mol_cat_fun_gaussian(beta,gamma,epsilon,hp_score,l,omega)
 

        kh_p = exp(-((hp_score-beta)^2)/omega)*0.1;
        kh_c = exp(-((hp_score-beta)^2)/omega)*0.1;%0.0546;%exp(-((hp_score-k_hs)^2)/0.02)*0.1;%0.0546;
        kl = (1/(1+exp(-l+gamma)))*epsilon; %sigmoid function, gamma shifts function horizontally and epsilon changes the plateau value.

        k_p = kh_p * kl;
        k_c = kh_c * kl;

end
