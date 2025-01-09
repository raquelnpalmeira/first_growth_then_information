%% Figure 3



figure()
tiledlayout(2,3);


%% FIRST PLOT TIMESERIES
legends = {'p_p = 0.0001','p_p = 0.001', 'p_p = 0.02'};
%{'\gamma_c = 6, \gamma_p = 10, \alpha = 0.2, \beta = -0.8','\gamma_c = 6, \gamma_p = 10, \alpha = 0.8, \beta = -0.2','\gamma_c = 10, \gamma_p = 6, \alpha = 0.2, \beta = -0.8', '\gamma_c = 10, \gamma_p = 6, \alpha = 0.8, \beta = -0.2'};%{'p_t = 0.0001', 'p_t = 0.001'};%, 'p_t = 0.01' 'p_t = 0.1'};%{'p_t = 0','p_t = 0.001', 'p_t = 0.01', 'p_t = 0.1'};
outputs1 = {output, output2,output6};

ax1 = nexttile([1 2]); %create tile that spams 1 row by 3 columns

for j = 1:length(outputs1)
div_per_gen_reps = outputs1{j}{1,2};


n_gen = length(div_per_gen_reps{1,1});
reps = length(div_per_gen_reps);
interval_for_moving_sum = 50;

%colourpalvar=[87 120 164;229 146 67;88 151 146;213 97 95;106 159 89;151 120 110;168 203 229;246 192 135;147 188 186;242 163 159;156 206 133;208 183 168;]/256;
%colourpalvar=[213 97 95;87 120 164;229 146 67;88 151 146;106 159 89;151 120 110;168 203 229;246 192 135;147 188 186;242 163 159;156 206 133;208 183 168;]/256;
%colourpalvar=[230 159 0;86 180 233;0 158 115;240 228 66;0 114 178;204 121 167;168 203 229;246 192 135;147 188 186;242 163 159;156 206 133;208 183 168;]/256;
colourpalvar=[0.3010 0.7450 0.9330;0.8500 0.3250 0.0980;0.6350 0.0780 0.1840;];

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
ylabel({['divisions per ',num2str(interval_for_moving_sum),' time steps']}, 'Fontsize',20,'FontName', 'Times', 'DisplayName',legends{j});
xlabel(['Time steps'],'Fontsize',20,'FontName', 'Times');
hold on
end
xline(ax1,8000,'--','Color','red', 'LineWidth', 2); % add vertical line
legend({'p_p = 0.0001','p_p = 0.001', 'p_p = 0.02'}, 'Fontsize', 18,'FontName', 'Times', 'Location', 'west', 'NumColumns',1);
xlim([0,10000-interval_for_moving_sum]);
%text(115, 430,'\it{(a)}', 'FontSize',24,'FontName', 'Times'); % This adds the "A" to the top left hand corner of the plot


%% boxplots 

nexttile([1 1]); %create tile that spams 1 row by 21 column

colourpalvar=[0.3010 0.7450 0.9330;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840;];

total_time = 8000;
outputs = {output, output2, output3, output4, output5, output6};
avg_div_per_gen = zeros(reps,2000);
reps = 10;
avg_div = zeros(reps,length(outputs));

for i = 1:length(outputs)
    for j = 1:reps
    x = outputs{i}{1,2}{1,j};
    avg_div_per_gen(j,:) = x(6001:8000);
    avg_div(j,i) = mean(avg_div_per_gen(j,:));
    end  
end


avg_div_per_50 = avg_div .* 50;

boxplot(avg_div_per_50, 'Colors', colourpalvar, 'Widths',1, 'Symbol', '+');
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times';
bx = findobj('Tag','boxplot');
set(bx.Children,'LineWidth',1.5)
xticklabels({'0.0001','0.001', '0.0025', '0.005','0.01','0.02'});
xlabel("p_p",'Fontsize',20,'FontName', 'Times');
ylabel([['divisions per ',num2str(interval_for_moving_sum),' time steps'],  "at equilibrium"],'Fontsize',20,'FontName', 'Times');
%text(0.1, 360,'\it(b)', 'FontSize',24,'FontName', 'Times'); % This adds the "A" to the top left hand corner of the plot


%% THEN PLOT HEATMAPS
titles = {'p_p = 0.0001','p_p = 0.001', 'p_p = 0.02'};%{'p_d = 0.0001', 'p_d = 0.001', 'p_d = 0.01'};
%letters = {'\it{(c)}','\it{(d)}','\it{(e)}','\it{(f)}','\it{(g)}'};
snps = 3; %which snaptshot do I wanna look at?
reps = 10; %how many reps are there in the data
cell_matrix_size = 60;
avg_cells_at_nth_snapshot_reps = cell(1,reps);
avg_cells_at_nth_snapshot = cell(1,length(outputs1));


for i = 1:length(outputs1)
    nucs_snapshot_cells_reps = outputs1{1,i}{1,3}; % First curly brackets says which bit of data (i.e. output of which parameter), and the second index (set as 3) tells us to get the snapshot cell compositions
    for j = 1:reps
        nucs_snapshot_cells_nth_rep = nucs_snapshot_cells_reps{j};
        nucs_snapshot_cells_nth_rep_nth_snapshot = nucs_snapshot_cells_nth_rep(:,snps);
        tri_matrix = cat(3, nucs_snapshot_cells_nth_rep_nth_snapshot{:});
        x = mean(tri_matrix,3);
        x_corrected = pop4heatmap_fun({x});
        if length(x_corrected{:}) < cell_matrix_size
            x_corrected{:}(cell_matrix_size,cell_matrix_size)=0;
        elseif length(x_corrected{:}) > cell_matrix_size
            disp("Warning! Matrix larger than maximum size, pick cell_matrix_size again!");
        end
        avg_cells_at_nth_snapshot_reps{j} = x_corrected{:};
    end
    tri_matrix2 = cat(3, avg_cells_at_nth_snapshot_reps{:});
    x2 = mean(tri_matrix2,3);
    avg_cells_at_nth_snapshot{i} = squeeze(x2);
end


for k = 1:length(avg_cells_at_nth_snapshot)
    nexttile;
    surf(avg_cells_at_nth_snapshot{k});  view(2);
    xlim([1,21]);
    ylim([1,21]);
    xticks(6:5:21);
    yticks(6:5:21);
    xticklabels(5:5:20);
    yticklabels(5:5:20);
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times';
    xlabel('number of purines','FontName', 'Times', 'Fontsize',20);
    ylabel('number of pyrimidines','FontName', 'Times','Fontsize',20);  
    
    set(gca,'ColorScale','log');
    title(titles{k},'Fontsize',20,'FontName', 'Times');
    %text(1.5, 20,letters{k}, 'Color','white','FontSize',24,'FontName', 'Times'); % This adds the letter to the top left hand corner of the plot
    caxis([0.00051, 51.67]);%caxis([0, 51.6730]); %renamed clim from 2022a 
    %caxis([0.0004,330]); %renamed clim from 2022a 
    %colorbar;
    lims = caxis;
    disp(lims)
    
end
 colorbar; 
 