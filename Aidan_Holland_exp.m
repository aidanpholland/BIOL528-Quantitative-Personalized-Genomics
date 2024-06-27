function Aidan_Holland_exp(starting_number, growth_rate, percent_noise)
    % This function simulates cell growth with no noise
    [cell_counter, steps] = sim_growth(starting_number, growth_rate);
    % This function simulates 9 doublings and plots the effect of starting cell count on final cell count
    ninedoublings(growth_rate);
    % This function simulates cell growth with noise
    [sampling_times, cell_count_noise] = sim_growth_noise(starting_number, growth_rate, percent_noise);
    % This function simulates exponential growth for modeling
    [sampling_times_model, cell_count_model] = exp_growth(starting_number, growth_rate, false);
    
    % The following code recognizes that the noisy cell count has error and
    % moves to optimize it using an inital guess and fminsearch so that it
    % better represents an exponential growth model.
    objective_function = @(parameters) compute_error(parameters, cell_count_noise);
    initial_guess = [starting_number, growth_rate];
    optimal_parameters = fminsearch(objective_function, initial_guess);
    optimal_starting_number = optimal_parameters(1);
    optimal_growth_rate = optimal_parameters(2);
    [~, optimal_cell_count_model] = exp_growth(optimal_starting_number, optimal_growth_rate, false);
    disp(['Optimized parameters: Starting Number = ', num2str(optimal_starting_number), ', Growth Rate = ', num2str(optimal_growth_rate)]);

    
    % The following code produces a plot with specified colors and line
    % width to the respective lines. Plotted is cell growth with noise, an
    % exponential growth model, and cell growth optimized using the model.
    figure;
    plot(sampling_times, cell_count_noise, '-g', 'LineWidth', 0.5);
    hold on;
    plot(sampling_times_model, cell_count_model, '-b', 'LineWidth', 0.5);
    plot(sampling_times_model, optimal_cell_count_model, '-r', 'LineWidth', 0.5);
    xlabel('Time');
    ylabel('Cell Count');
    title('Simulated Data + Noise and Modeled Data');
    legend('Simulated Data + Noise', 'Modeled Data', 'Optimized Fit');
    grid on;
    hold off;
end

function [intervals, cell_count] = sim_growth(starting_number, growth_rate)
    % Using the rule of 70, we compute doubling time using growth rate.
    doubling_time = 69.314718056 / growth_rate;
    six_doublings = 6 * doubling_time;
    disp(['The doubling time for six doublings is ', num2str(six_doublings), ' units of time.']);
    % 100 intervals over the course of six_doublings time are made, and
    % zeroes for each of those positions to be assigned a value depending
    % on the for loop below that represents growth
    intervals = linspace(0, six_doublings, 100);
    cell_count = zeros(1, length(intervals));
    cell_count(1) = starting_number;
    for t = 1:length(intervals) - 1
        dt = intervals(t + 1) - intervals(t);
        total_cells = cell_count(t) * (2^(dt / doubling_time));
        cell_count(t + 1) = total_cells;
    end
    disp(['Final cell count: ', num2str(total_cells)]);
    
    % This plots the data just produced
    plot(intervals, cell_count, 'b.-', 'LineWidth', 0.5);
    xlabel('Time');
    ylabel('Cell Count');
    title('Cell Count Over Time');
    grid on;
    legend(['Growth Rate: ', num2str(growth_rate)], 'Location', 'northwest');
end

function ninedoublings(growth_rate)
    % This code is similar to sim_growth but instead produces a linear
    % plot of final number vs. starting number after 9 doublings.
    growth_rate_ac = growth_rate / 100;
    doubling_time = 69.314718056 / growth_rate_ac;
    nine_doublings = doubling_time * 9;
    starting_numbers = 1:10000;
    final_cell_count = zeros(size(starting_numbers));
    for i = 1:length(starting_numbers)
        final_cell_count(i) = starting_numbers(i) * (2^(nine_doublings / doubling_time));
    end
    % This code plots our graph
    plot(starting_numbers, final_cell_count, '-');
    xlabel('Starting Number of Cells');
    ylabel('Final Number of Cells after 9 Doublings');
    title('Effect of Starting Cell Count on Final Cell Count');
    slope = (final_cell_count(end) - final_cell_count(1)) / (starting_numbers(end) - starting_numbers(1));
    legendText = sprintf('Slope: %.2f', slope);
    legend(legendText);
    saveas(gcf, 'Figure2_Cell_Count_Plot.pdf');
end

function [sampling_times, cell_count] = sim_growth_noise(starting_number, growth_rate, percent_noise)
    [intervals, cell_count_base] = sim_growth(starting_number, growth_rate);
    noise_amplitude = percent_noise / 100 * max(cell_count_base); 
    % The randn() function simulates the amplitude that will represent our
    % noise in the graph.
    noise = noise_amplitude * randn(size(cell_count_base));
    % The new cell count will be modified by whatever the noise was for
    % that given point.
    cell_count = cell_count_base + noise;
    sampling_times = intervals;
end

function [sampling_times, cell_count] = exp_growth(starting_number, growth_rate, gen_plot)
    % Rather than the more arduous task that we used with sim_growth, this
    % code creates a model of exponential growth using exp() and plots it.
    sampling_times = linspace(0, 6 * (70 / growth_rate), 100);
    cell_count = starting_number * exp((growth_rate / 100) * sampling_times);
    if gen_plot
        plot(sampling_times_exp, cell_count_exp, '-');
        xlabel('Time');
        ylabel('Cell Count');
        title('Exponential Growth Curve');
        legend('Cell Count');
    end
end

function [error] = compute_error(parameters, cell_count_noise)
    % In order to optimize our noisy cell growth, the error between that
    % and the model (exp_growth function) must be calculated.
    starting_number = parameters(1);
    growth_rate = parameters(2);
    [sampling_times, cell_count] = exp_growth(starting_number, growth_rate, false);
    error = sum((cell_count_noise - cell_count).^2);
    error = sqrt(error);
end
