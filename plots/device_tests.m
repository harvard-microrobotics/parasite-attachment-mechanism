clc;
close all;
clear;

%% PARAMS

% Colors
paper_main_blue = [5, 166, 166]/255;
paper_main_magenta = [166, 51, 93]/255;

% Plot parameters
line_width = 2;
legend_font_size = 20;
ticks_font_size = 15;
axes_font_size = 25;
epsil_display = 1e-3;


% pull test plot params
disp_step = 1e-4;
diff_smooth_factor = 10;

% trigger plot params
peak_window_width = 80;
diff_treshold = 0.08;


%% PULL TEST PLOT
F = [];

x_min= 0;
x_max = 0;

nb_tests = 3;
nb_devices = 4;

figure
for i = 1:nb_devices
    for j = 1:nb_tests

        %read csv
        file_path = strcat('exp_data/TAM_device_',int2str(i),'.is_tens_Exports/TAM_device_',int2str(i),'_',int2str(j),'.csv');
        data = readPullTestData(file_path);
        displacement = data(:,2);
        force = data(:,3);

        % compute offsets
        [~,detach_idx] = min(diff(smooth(force,diff_smooth_factor,'loess')));
        detach_idx = detach_idx + 1;
        disp_offset = displacement(detach_idx);
        force_offset = mean(force(detach_idx:end));
        
        % align curves
        displacement = displacement - disp_offset;
        force = force - force_offset;
        
        % limits of x-axis      
        if i == 1 && j == 1
            x_min = displacement(1);
            x_max = displacement(end);
        else
            if displacement(1) > x_min
                x_min_prev = x_min;
                x_min = displacement(1);
                idx_min = round((x_min-x_min_prev)/disp_step) + 1;
                F = F(idx_min:end,:);
            end
            if displacement(end) < x_max
                x_max_prev = x_max;
                x_max = displacement(end);
                idx_max = size(F,1) - round((x_max_prev-x_max)/disp_step);
                F = F(1:idx_max,:);
            end   
        end
        
        [disp_unique,ia,~] = unique(displacement,'legacy');
        force_unique = force(ia);
        
        force_inter = interp1(disp_unique,force_unique, x_min:disp_step:x_max);
       
        F(:,(i-1)*nb_tests + j) = force_inter'; %#ok<SAGROW>
        
        % plot inidvidual curves
        plot(x_min:disp_step:x_max, force_inter);
        hold on
    end
end
set(gca,'FontSize',ticks_font_size)
xlabel('Distance to detachment point [mm]','Interpreter','latex', 'FontSize', axes_font_size);
ylabel('Pulling force [mN]', 'Interpreter','latex', 'FontSize', axes_font_size);

% plot
figure
F_mean = mean(F,2);
F_std = std(F,0,2);
X = x_min:disp_step:x_max;

% Plot bounded line. Kelly Kearney, boundedline.m, https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m
boundedline(X,F_mean,F_std,'Color',paper_main_blue,'LineWidth',line_width)


% plotting settings
xlim([X(1)-epsil_display Inf]); % mm
ylim([-Inf 70]); % mN


lgd = legend('$\pm 1\sigma$','Mean force','Interpreter','latex','Location','northwest');
lgd.FontSize = legend_font_size;

labels = get(legend(), 'String');
plots = flipud(get(gca, 'children'));
neworder = [2 1];
legend(plots(neworder), labels(neworder))

set(gca,'FontSize',ticks_font_size)
xlabel('Distance to detachment point [mm]','Interpreter','latex', 'FontSize', axes_font_size);
ylabel('Pulling force [mN]', 'Interpreter','latex', 'FontSize', axes_font_size);


%% TRIGGER TEST PLOT

F = [];

x_min= 0;
x_max = 0;
x_min_prev = 0;
x_max_prev = 0;

nb_tests = 5;

figure
for i = 1:nb_tests

    %read csv
    file_path = strcat('exp_data/trigger.is_comp_Exports/trigger_1_',int2str(i),'.csv');
    data = readPullTestData(file_path);
    time = data(:,1);
    force = data(:,2);

    %find last peak
    [~,pks_idx] = findpeaks(force,'MinPeakHeight',1.5,'MinPeakWidth',5);
    last_pk_idx = pks_idx(end);
    time_last_peak = time(last_pk_idx-peak_window_width/2:last_pk_idx+peak_window_width/2);
    force_last_peak = force(last_pk_idx-peak_window_width/2:last_pk_idx+peak_window_width/2);


    % compute offsets
    diff_force = diff(force_last_peak);
    [~,min_idx] = min(diff_force);
    post_min_pos_diff = find(diff_force(min_idx:end) > diff_treshold);
    end_idx = post_min_pos_diff(1) + min_idx + 1;
    time_offset = time_last_peak(end_idx);

    % align curves
    time_last_peak = time_last_peak - time_offset;
    
    % limits of x-axis      
    if i == 1
        x_min = time_last_peak(1);
        x_max = time_last_peak(end);
    else
        if time_last_peak(1) > x_min
            x_min_prev = x_min;
            x_min = time_last_peak(1);
            idx_min = round((x_min-x_min_prev)/disp_step) + 1;
            F = F(idx_min:end,:);
        end
        if time_last_peak(end) < x_max
            x_max_prev = x_max;
            x_max = time_last_peak(end);
            idx_max = size(F,1) - round((x_max_prev-x_max)/disp_step)-1;
            F = F(1:idx_max,:);
        end   
    end
    
    [time_unique,ia,~] = unique(time_last_peak,'legacy');
    force_unique = force_last_peak(ia);

    force_inter = interp1(time_unique,force_unique, x_min:disp_step:x_max);

    F(:,i) = force_inter'; %#ok<SAGROW>

    % plot inidvidual curves
    plot(x_min:disp_step:x_max, force_inter);
    hold on
end
set(gca,'FontSize',ticks_font_size)
xlabel('Time [s]','Interpreter','latex', 'FontSize', axes_font_size);
ylabel('Triggering force [N]', 'Interpreter','latex', 'FontSize', axes_font_size);

% plot
figure
F_mean = mean(F,2);
F_std = std(F,0,2);
T = x_min:disp_step:x_max;
[~,mean_max_idx] = max(F_mean);
T = T - T(mean_max_idx);

% Plot bounded line. Kelly Kearney, boundedline.m, https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m
boundedline(T,F_mean,F_std,'Color',paper_main_magenta,'LineWidth',2)


% plotting settings
xlim([T(1)-epsil_display Inf]); % mm
ylim([-Inf 7]); % mN


lgd = legend('$\pm 1\sigma$','Mean force','Interpreter','latex','Location','northwest');
lgd.FontSize = legend_font_size;

labels = get(legend(), 'String');
plots = flipud(get(gca, 'children'));
neworder = [2 1];
legend(plots(neworder), labels(neworder))

set(gca,'FontSize',ticks_font_size)
xlabel('Time [s]','Interpreter','latex', 'FontSize', axes_font_size);
ylabel('Triggering force [N]', 'Interpreter','latex', 'FontSize', axes_font_size);
