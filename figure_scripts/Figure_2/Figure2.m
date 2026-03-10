
%% Figure 2

figure()
tiledlayout(2,2);


%%  Timecourse

data = {data_6828,data_8682};

%%%% FIRST PLOT TIMESERIES
legends = {'Easy CO_{2} Fixation','Easy Templated Polymerisation'};%

ax1 = nexttile([1 2]); %create tile that spams 1 row by 3 columns

for j = 1:length(data)
div_per_gen_reps = data{j}{1,1}; %use {1,1} if the data has been cleaned because there are many reps, {1,2} for smaller datasets


n_gen = length(div_per_gen_reps{1,1});
reps = length(div_per_gen_reps);
interval_for_moving_sum = 50;


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
legend({'Easy CO_{2} Fixation','Easy Templated Polymerisation'}, 'Fontsize', 20,'FontName', 'Times', 'Location', 'west');
xlim([0,10000-interval_for_moving_sum]);



%% Plot heatmaps

titles = {'Easy CO_{2} Fixation','Easy Templated Polymerisation'};%{'p_p = 0.0025','p_p = 0.005', 'p_p = 0.01', 'p_p = 0.02'};%{'\gamma_{fix} < \gamma_t and |\alpha_{fix}| > |\alpha_t|','\gamma_{fix} > \gamma_t and |\alpha_{fix}| < |\alpha_t|'};%{'p_d = 0.0001', 'p_d = 0.001', 'p_d = 0.01'};
snps = 4; %which snaptshot do I wanna look at?
reps = 10; %how many reps are there in the data
cell_matrix_size = 60;
avg_cells_at_nth_snapshot_reps = cell(1,reps);
avg_cells_at_nth_snapshot = cell(1,length(data));


for i = 1:length(data)
    nucs_snapshot_cells_reps = data{1,i}{1,2}; % %use {1,2} if the data has been cleaned because there are many reps, {1,3} for smaller datasets. First curly brackets says which bit of data (i.e. output of which parameter), and the second index (set as 3) tells us to get the snapshot cell compositions
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
    nexttile([1,1]);
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
    xlabel('Number of purines','FontName', 'Times', 'Fontsize',20);
    ylabel('Number of pyrimidines','FontName', 'Times','Fontsize',20);  
    
    set(gca,'ColorScale','log');
    title(titles{k},'Fontsize',20,'FontName', 'Times');
    caxis([0.00051   51.67]);%caxis([0, 51.6730]); %renamed clim from 2022a 
    lims = caxis;
    disp(lims)
    
end
 colorbar; 

