function aidan_holland_clust(FileName, group_size)
    data_table = readtable(FileName);
    num_data = table2array(data_table(:, 2:end));
    txt_data = data_table(:, 1);
    [number_of_people, number_of_questions] = size(num_data);
    n_groups = ceil(number_of_people / group_size);
    disp(['Number of people: ', num2str(number_of_people)]);
    disp(['Number of questions: ', num2str(number_of_questions)]);
    disp(['Number of groups: ', num2str(n_groups)]);
    Identifiers = cell(number_of_people, 1);
    for i = 1:number_of_people
        if isnan(num_data(i, 1))
            Identifiers{i} = '';
        else
            Identifiers{i} = txt_data{i, 1};
        end
    end
    disp('All Identifiers:');
    disp(Identifiers);
    % Identifiers{3}
    grouping = generate_random_grouping(number_of_people, n_groups);
    disp('Randomly Generated Grouping:');
    disp(grouping);
    % grouping{2}
    score = score_grouping(grouping, num_data, n_groups)
    % disp(score)

    min_score = Inf;
    best_grouping = {};
    
    for trial = 1:1000000
        grouping = generate_random_grouping(number_of_people, n_groups);
        score = score_grouping(grouping, num_data, n_groups);
        if score < min_score
            min_score = score;
            best_grouping = grouping;
        end
    end
    
    disp('Best Grouping:');
    disp(best_grouping);
    disp('Minimum Score:');
    disp(min_score);

end

function grouping = generate_random_grouping(number_of_people, n_groups)
    random_indices = randperm(number_of_people);
    group_size = ceil(number_of_people / n_groups);
    grouping = cell(n_groups, 1);
    for i = 1:n_groups
        start_index = (i - 1) * group_size + 1;
        end_index = min(i * group_size, number_of_people);
        grouping{i} = random_indices(start_index:end_index);
    end
end

function score = score_grouping(grouping, num_data, n_groups)
    n_qgroups = 5;
    mean_scores = zeros(n_groups, n_qgroups);
    aggregate_scores = zeros(n_groups, 1);
    for i = 1:n_groups
        indexed_groups = grouping{i};
        group_data = num_data(indexed_groups, :);
        for j = 1:n_qgroups
            group_columns = (j - 1) * 3 + 1 : j * 3;
            group_question_data = group_data(:, group_columns);
            mean_scores(i, j) = mean(group_question_data, 'all');
        end
        aggregate_scores(i) = sum(mean_scores(i, :));
    end
    total_aggregate_score = sum(aggregate_scores);
    target_aggregate_score = total_aggregate_score / n_groups;
    deviation = abs(aggregate_scores - target_aggregate_score);
    score = mean(deviation);
end