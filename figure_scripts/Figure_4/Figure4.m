
%% Figure 4

figure()
tiledlayout(1,2);


%%  Timecourse Easy Growth

data = {data_pm0_01_6828,data_pm0_001_6828,data_pm0_0001_6828,data_pm0_00001_6828};
%%%% FIRST PLOT TIMESERIES
legends = {'p_{n}=p_{aa}=0.01','p_{n}=p_{aa}=0.001','p_{n}=p_{aa}=0.0001','p_{n}=p_{aa}=0.00001'};%

ax1 = nexttile([1 1]); %create tile that spams 1 row by 3 columns

for j = 1:length(data)
div_per_gen_reps = data{j}{1,1}; %use {1,1} if the data has been cleaned because there are many reps, {1,2} for smaller datasets


n_gen = length(div_per_gen_reps{1,1});
reps = length(div_per_gen_reps);
interval_for_moving_sum = 50;

%colourpalvar=[87 120 164;229 146 67;88 151 146;213 97 95;106 159 89;151 120 110;168 203 229;246 192 135;147 188 186;242 163 159;156 206 133;208 183 168;]/256;
%colourpalvar=[213 97 95;87 120 164;229 146 67;88 151 146;106 159 89;151 120 110;168 203 229;246 192 135;147 188 186;242 163 159;156 206 133;208 183 168;]/256;
%colourpalvar=[0.3010 0.7450 0.9330;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840;];
colourpalvar = [
    0.3010 0.7450 0.9330;  % light blue
    0.8500 0.3250 0.0980;  % orange
    0.9290 0.6940 0.1250;  % yellow
    0.4940 0.1840 0.5560;  % purple
    0.4660 0.6740 0.1880;  % green
    0.6350 0.0780 0.1840;  % dark red
    0.0000 0.4470 0.7410;  % classic MATLAB blue
    0.3010 0.3000 0.3000   % neutral dark grey
];

y = colourpalvar;
colourpalvar2 = y.^.3; %any colour palette but lighter
colour = j;

%subplot(2,3,[1:3]);

for i = 1:reps 
x = movsum(div_per_gen_reps{1,i},interval_for_moving_sum, "Endpoints","discard");
plot(x, 'Color', colourpalvar2(colour,:), 'LineStyle','--', 'HandleVisibility','off'); % last argument so that these lines don't appear in the legend 
drawnow
hold on
end

div_per_gen_mat = cell2mat(div_per_gen_reps');
div_per_gen_avg = mean(div_per_gen_mat,1);

y = movsum(div_per_gen_avg,interval_for_moving_sum, "Endpoints","discard");
plot(y, 'Color', colourpalvar(colour,:), 'LineStyle','-','LineWidth',2);
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times';
ylabel({'Division rate'}, 'Fontsize',20,'FontName', 'Times', 'DisplayName',legends{j});
xlabel(['Time steps'],'Fontsize',20,'FontName', 'Times');
hold on
end
xline(ax1,8000,'--','Color','red', 'LineWidth', 2); % add vertical line
legend({'p_{n}=p_{aa}=0.01','p_{n}=p_{aa}=0.001','p_{n}=p_{aa}=0.0001','p_{n}=p_{aa}=0.00001'}, 'Fontsize', 20,'FontName', 'Times', 'Location', 'west');
%legend({'\gamma_{fix} < \gamma_t and |\alpha_{fix}| < |\alpha_t|','\gamma_{fix} < \gamma_t and |\alpha_{fix}| > |\alpha_t|','\gamma_{fix} > \gamma_t and |\alpha_{fix}| < |\alpha_t|','\gamma_{fix} > \gamma_t and |\alpha_{fix}| > |\alpha_t|'}, 'Fontsize', 20,'FontName', 'Times', 'Location', 'west');
xlim([0,10000-interval_for_moving_sum]);
%text(115, 245,'\it{(c)}', 'FontSize',20,'FontName', 'Times'); % This adds the "A" to the top left hand corner of the plot

%%  Timecourse Easy templated polymerisation

data = {data_pm0_01_8682,data_pm0_001_8682,data_pm0_0001_8682,data_pm0_00001_8682};
%%%% FIRST PLOT TIMESERIES
legends = {'p_{n}=p_{aa}=0.01','p_{n}=p_{aa}=0.001','p_{n}=p_{aa}=0.0001','p_{n}=p_{aa}=0.00001'};%

ax1 = nexttile([1 1]); %create tile that spams 1 row by 3 columns

for j = 1:length(data)
div_per_gen_reps = data{j}{1,1}; %use {1,1} if the data has been cleaned because there are many reps, {1,2} for smaller datasets


n_gen = length(div_per_gen_reps{1,1});
reps = length(div_per_gen_reps);
interval_for_moving_sum = 50;

%colourpalvar=[87 120 164;229 146 67;88 151 146;213 97 95;106 159 89;151 120 110;168 203 229;246 192 135;147 188 186;242 163 159;156 206 133;208 183 168;]/256;
%colourpalvar=[213 97 95;87 120 164;229 146 67;88 151 146;106 159 89;151 120 110;168 203 229;246 192 135;147 188 186;242 163 159;156 206 133;208 183 168;]/256;
%colourpalvar=[0.3010 0.7450 0.9330;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840;];
colourpalvar = [
    0.3010 0.7450 0.9330;  % light blue
    0.8500 0.3250 0.0980;  % orange
    0.9290 0.6940 0.1250;  % yellow
    0.4940 0.1840 0.5560;  % purple
    0.4660 0.6740 0.1880;  % green
    0.6350 0.0780 0.1840;  % dark red
    0.0000 0.4470 0.7410;  % classic MATLAB blue
    0.3010 0.3000 0.3000   % neutral dark grey
];

y = colourpalvar;
colourpalvar2 = y.^.3; %any colour palette but lighter
colour = j;

%subplot(2,3,[1:3]);

for i = 1:reps 
x = movsum(div_per_gen_reps{1,i},interval_for_moving_sum, "Endpoints","discard");
plot(x, 'Color', colourpalvar2(colour,:), 'LineStyle','--', 'HandleVisibility','off'); % last argument so that these lines don't appear in the legend 
drawnow
hold on
end

div_per_gen_mat = cell2mat(div_per_gen_reps');
div_per_gen_avg = mean(div_per_gen_mat,1);

y = movsum(div_per_gen_avg,interval_for_moving_sum, "Endpoints","discard");
plot(y, 'Color', colourpalvar(colour,:), 'LineStyle','-','LineWidth',2);
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times';
ylabel({'Division rate'}, 'Fontsize',20,'FontName', 'Times', 'DisplayName',legends{j});
xlabel(['Time steps'],'Fontsize',20,'FontName', 'Times');
hold on
end
xline(ax1,8000,'--','Color','red', 'LineWidth', 2); % add vertical line
legend({'p_{n}=p_{aa}=0.01','p_{n}=p_{aa}=0.001','p_{n}=p_{aa}=0.0001','p_{n}=p_{aa}=0.00001'}, 'Fontsize', 20,'FontName', 'Times', 'Location', 'west');
%legend({'\gamma_{fix} < \gamma_t and |\alpha_{fix}| < |\alpha_t|','\gamma_{fix} < \gamma_t and |\alpha_{fix}| > |\alpha_t|','\gamma_{fix} > \gamma_t and |\alpha_{fix}| < |\alpha_t|','\gamma_{fix} > \gamma_t and |\alpha_{fix}| > |\alpha_t|'}, 'Fontsize', 20,'FontName', 'Times', 'Location', 'west');
xlim([0,10000-interval_for_moving_sum]);
%text(115, 245,'\it{(c)}', 'FontSize',20,'FontName', 'Times'); % This adds the "A" to the top left hand corner of the plot


