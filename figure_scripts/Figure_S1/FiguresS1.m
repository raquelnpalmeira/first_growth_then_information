%% Figure S1 A

hp_score = 1;
lengths = 0:0.1:15;
alpha = 1;
gamma = 4:2:10;% shifts curve to the right
kl = zeros(length(lengths),length(gamma));
epsilon = 1.5; % changes plateau 

for i = 1:length(gamma)
for j = 1:length(lengths)
    
    [~, ~, kl(j,i),~,~] = mol_cat_fun(alpha,gamma(i),epsilon,hp_score,lengths(j));
    
end
end

colourpalvar=[213 97 95;87 120 164;229 146 67;88 151 146;106 159 89;151 120 110;168 203 229;246 192 135;147 188 186;242 163 159;156 206 133;208 183 168;]/256;
figure('Name','catalysis based on length')
ax1 = subplot(1,2,1);
plot(lengths,[kl],'LineWidth',2);
legend({'\gamma = 4', '\gamma = 6', '\gamma = 8', '\gamma = 10'},'Location', 'southeast')
xlim(ax1,[1,15]);
xlabel(ax1, 'length');
ylabel(ax1, '\it{k_l}');
title(ax1, 'Effect of \gamma' );
set(ax1, 'FontSize', 18,'FontName', 'Times');
ax1.ColorOrder = colourpalvar;

hp_score = 1;
lengths = 0:0.1:15;
alpha = 1;
gamma = 7;% shifts curve to the right
epsilon = 1.5:1.5:6; % changes plateau 
kl = zeros(length(lengths),length(epsilon));


for i = 1:length(epsilon)
for j = 1:length(lengths)
    
    [~, ~, kl(j,i),~,~] = mol_cat_fun(alpha,gamma,epsilon(i),hp_score,lengths(j));
    
end
end

ax2 = subplot(1,2,2);
plot(lengths,[kl],'LineWidth',2);
legend({'\it{k_{l}max} = 1.5', '\it{k_{l}max} = 3', '\it{k_{l}max} = 4.5', '\it{k_{l}max} = 6'},'Location', 'west')
xlim(ax2,[1,15]);
xlabel(ax2, 'length');
ylabel(ax2, '\it{k_l}');
title(ax2, 'Effect of \it{k_{l}max}' );
set(ax2, 'FontSize', 18,'FontName', 'Times');
ax2.ColorOrder = colourpalvar;



%% Figure S1 B
colourpalvar=[0.3010 0.7450 0.9330;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;];
beta = [-0.8:0.4:0.8];%-1:0.2:1;%1;%0:0.2:1;% 0:0.2:2; % catalytic constant for catalysis based on hydrophobicity
hp_score = -1:0.01:1;
lengths = 8;%10;
k_p = zeros(length(hp_score),length(beta));
k_c = zeros(length(hp_score),length(beta));
k_s = 1;
kl_max =  1.5; % "vmax" for catalysis based on length
omega = 0.02;
tiledlayout(1,2)

for j = 1:length(hp_score)
for i = 1:length(beta)
    
    [k_p(j,i),k_c(j,i),] = mol_cat_fun_gaussian(beta(i),kl_max,k_s,hp_score(j),lengths,omega);
    
end
end

ax1 = nexttile;
plot(hp_score,[k_c],'LineWidth',2);
xlabel(ax1, 'hydrophobicity score');
ylabel(ax1, '\it{k_h}');
title(ax1, 'Effect of \beta' );
set(ax1, 'FontSize', 16,'FontName', 'Times');
ax1.ColorOrder = colourpalvar;
legend(['\beta= ',num2str(beta(1))],['\beta = ',num2str(beta(2))],['\beta = ',num2str(beta(3))],['\beta = ',num2str(beta(4))],['\beta = ',num2str(beta(5))], 'Location', 'southeast');


%%%

beta = 0;%-1:0.2:1;%1;%0:0.2:1;% 0:0.2:2; % catalytic constant for catalysis based on hydrophobicity
hp_score = -1:0.01:1;
lengths = 8;%10;
k_p = zeros(length(hp_score),length(beta));
k_c = zeros(length(hp_score),length(beta));
k_s = 1;
kl_max =  1.5; % "vmax" for catalysis based on length
constant = 0;%[1,2,3,4,5];% 0.1:0.1:1; 
omega = 0.01:0.02:0.1;
for j = 1:length(hp_score)
for i = 1:length(omega)
    
    [k_p(j,i),k_c(j,i),] = mol_cat_fun_gaussian(beta,kl_max,k_s,hp_score(j),lengths,omega(i));
    
end
end

ax2 = nexttile;
plot(hp_score,[k_c],'LineWidth',2);
xlabel(ax2, 'hydrophobicity score');
ylabel(ax2, '\it{k_h}');
title(ax2, 'Effect of \omega' );
set(ax2, 'FontSize', 16,'FontName', 'Times');
ax2.ColorOrder = colourpalvar;
legend(['\omega = ',num2str(omega(1))],['\omega = ',num2str(omega(2))],['\omega = ',num2str(omega(3))],['\omega = ',num2str(omega(4))],['\omega = ',num2str(omega(5))], 'Location', 'southeast');



%% Figure S1 C

colourpalvar=[0.3010 0.7450 0.9330;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;];
figure()
tiledlayout(1,2)

alpha = 1:1:5;%0:0.2:1;%1;%0:0.2:1;% 0:0.2:2; % catalytic constant for catalysis based on hydrophobicity
hp_score = -1:0.01:1;
l = 8;%10;
k_p = zeros(length(hp_score),length(alpha));
k_c = zeros(length(hp_score),length(alpha));
gamma = 7;
epsilon = 1;
for j = 1:length(hp_score)
for i = 1:length(alpha)
    
    [k_p(j,i),k_c(j,i),] = mol_cat_fun_exp(alpha(i),gamma,epsilon,hp_score(j),l);
    
end
end

ax1 = nexttile;
plot(hp_score,k_c,'LineWidth',2);
title(ax1, 'Effect of \alpha_{fix}' );
set(ax1, 'FontSize', 18,'FontName', 'Times');
xlabel(ax1, 'hydrophobicity score');
ylabel(ax1, '\it{k_{fix}^h}');
legend(['\alpha_{fix} =',num2str(k_hs(1))],['\alpha_{fix} =',num2str(alpha(2))],['\alpha_{fix} =',num2str(alpha(3))],['\alpha_{fix} =',num2str(alpha(4))],['\alpha_{fix} =',num2str(alpha(5))]);
ax1.ColorOrder = colourpalvar;

%%%%
k_hs = 1:1:5;%0:0.2:1;%1;%0:0.2:1;% 0:0.2:2; % catalytic constant for catalysis based on hydrophobicity
hp_score = -1:0.01:1;
l = 8;%10;
k_p = zeros(length(hp_score),length(k_hs));
k_c = zeros(length(hp_score),length(k_hs));
gamma = 7;
epsilon = 1;
for j = 1:length(hp_score)
for i = 1:length(k_hs)
    
    [k_p(j,i),k_c(j,i),] = mol_cat_fun_exp(k_hs(i),gamma,epsilon,hp_score(j),l);
    
end
end

ax2 = nexttile;
plot(hp_score,k_p,'LineWidth',2);
title(ax2, 'Effect of \alpha_t =' );
set(ax2, 'FontSize', 18,'FontName', 'Times');
xlabel(ax2, 'hydrophobicity score');
ylabel(ax2, '\it{k_t^h}');
legend(['\alpha_t =',num2str(k_hs(1))],['\alpha_t =',num2str(k_hs(2))],['\alpha_t =',num2str(k_hs(3))],['\alpha_t =',num2str(k_hs(4))],['\alpha_t =',num2str(k_hs(5))], 'Location', 'northwest');
ax2.ColorOrder = colourpalvar;


