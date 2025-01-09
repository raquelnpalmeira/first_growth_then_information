
function [k_p, k_c, kl,kh_p,kh_c] = mol_cat_fun(k_hs,gamma,epsilon,hp_score,l)
 
        kh_p = exp(-((hp_score+k_hs)^2)/0.02)*0.0546;%(exp((hp_score*k_hs)))/1e3;%
        kh_c = exp(-((hp_score+k_hs)^2)/0.02)*0.0546;%(exp((-hp_score*k_hs)))/1e3;%
        kl = (1/(1+exp(-l+gamma)))*epsilon; %sigmoid function, gamma shifts function horizontally and epsilon changes the plateau value.
        k_p = kh_p * kl;
        k_c = kh_c * kl;

end










% function [k_p, k_c, kl,kh_p,kh_c] = mol_cat_fun(k_hs,gamma,epsilon,hp_score,l,constant)
%  
%         kh_p = constant * (exp(hp_score*k_hs) - exp(-k_hs));
%         kh_c = constant *exp(hp_score*(-k_hs)) - constant *exp(-k_hs));
%         kl = (1/(1+exp(-l+gamma)))*epsilon; %sigmoid function, gamma shifts function horizontally and epsilon changes the plateau value.
% 
%         k_p = kh_p * kl;
%         k_c = kh_c * kl;
% 
% end
