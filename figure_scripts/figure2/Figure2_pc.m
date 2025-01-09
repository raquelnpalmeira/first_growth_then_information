

 %% Timecourse plus heatmaps for successive snapshots in for all parameters



outputs = {output,output1};%, output1, output2, output3};%, output3, output4};%{job61_output, job60_output, job59_output};%{output, output2, output3};%

%%%% FIRST PLOT TIMESERIES
legends = {'p_c = 0.001', 'p_c = 0.1'};
figure()
tiledlayout(3,3);
ax1 = nexttile([1 3]); %create tile that spams 1 row by 3 columns

for j = 1:length(outputs)
div_per_gen_reps = outputs{j}{1,2};


n_gen = length(div_per_gen_reps{1,1});
reps = length(div_per_gen_reps);
interval_for_moving_sum = 50;

colourpalvar=[0.3010 0.7450 0.9330;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;];
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
ylabel({['divisions per ',num2str(interval_for_moving_sum),' time steps']},'Fontsize', 24,'FontName', 'Times', 'DisplayName',legends{j});
xlabel(['Time steps'],'Fontsize', 24,'FontName', 'Times');
hold on
end
xline(ax1,[1000 4000 8000],'--', 'Color', 'r', 'LineWidth', 2); % add vertical line use "-interval_for_moving_sum" if plotting the last time step
legend({'p_c = 0.001', 'p_c = 0.1'},'Fontsize', 20,'FontName', 'Times', 'Location', 'northeast');
xlim([0,10000-interval_for_moving_sum]);
%text(115, 130,'\it(a)', 'FontSize',20,'FontName', 'Times'); % This adds the "A" to the top left hand corner of the plot
% title('a) ','Fontsize',16, 'FontWeight','Normal');
% ax = gca;
% ax.TitleHorizontalAlignment = 'left';


%%%%% THEN PLOT HEATMAPS
 titles = {'t = 1000', 't = 4000','t = 8000'};%{'t = 1000', 't = 2000', 't = 3000','t = 4000','t = 8000'};
 %letters = {'\it(b)','\it(c)','\it(d)','\it(e)','\it(f)','\it(g)','\it(h)','\it(i)','\it(j)','\it(k)','\it(l)','\it(m)','\it(n)','\it(o)','\it(p)','\it(q)','\it(r)','\it(s)','\it(t)','\it(u)'};
snps = 5; %number of snapshots
reps = 10; %how many reps are there in the data
cell_matrix_size = 60;
avg_cells_at_nth_snapshot_reps = cell(1,reps);
avg_cells_at_nth_snapshot = cell(snps,length(outputs));


for i = 1:length(outputs)
    nucs_snapshot_cells_reps = outputs{1,i}{1,3}; % First curly brackets says which bit of data (i.e. output of which parameter), and the second index (set as 3) tells us to get the snapshot cell compositions
    for k = 1:snps
    for j = 1:reps
        nucs_snapshot_cells_nth_rep = nucs_snapshot_cells_reps{j};
        nucs_snapshot_cells_nth_rep_nth_snapshot = nucs_snapshot_cells_nth_rep(:,k);
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
    avg_cells_at_nth_snapshot{k,i} = squeeze(x2);
    end
end
i = 0;

j_index = [1,4,5];
for k = 1:length(outputs)
    for j = 1:3%1:5%snps
        nexttile;
        i = i+1;
        surf(avg_cells_at_nth_snapshot{j,k});  view(2);
        ax = gca;
        ax.FontSize = 14;
        ax.FontName = 'Times';
        xlim([1,20]);
        ylim([1,20]);
        %if j == 1 && k == 2
            ylabel('number of pyrimidines','Fontsize', 20,'FontName', 'Times');
        %end
        %if j == 3 && k == 3
            xlabel('number of purines','Fontsize', 20,'FontName', 'Times');
        %end
        %dim = [0.2 0.5 0.3 0.3];
        set(gca,'ColorScale','log');
        %if k == 1
        title(titles{j},'Fontsize', 20,'FontName', 'Times');
        %end
        %ax = gca;
        %ax.TitleHorizontalAlignment = 'left';
        %caxis([0.0004,1.3e3]); %renamed clim from 2022a
        %annotation('textbox',dim,'String',letters{i},'FitBoxToText','on');
        %text(1.5, 19.5,letters{i}, 'Color','white','FontSize',20,'FontName', 'Times', 'Horiz','left', 'Vert','top'); % This adds the letter to the top left hand corner of the plot
        %colorbar;
        lims = caxis;
        disp(lims)
    end
    colorbar;
end
