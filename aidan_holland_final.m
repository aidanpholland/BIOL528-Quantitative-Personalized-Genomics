function aidan_holland_final
    % The purpose of this exercise is to generate a population of individuals with two independent SNPs and determine the most common genotype.    
    % set seed for reproducibility
    rng(123);

    % define SNP frequencies matrix as provided in instructions
    % oddly, the code suggests that 
    SNP_frequencies = [0.2 0.3 0.5; 0.1 0.4 0.5];

    % initialize a cell array to store genotypes
    genotypes = cell(1, 10000);

    % generate 10000 genotypes
    for i = 1:10000
        genotypes{i} = generate_individual(SNP_frequencies);
    end

    % display all 10000 genotypes
    for i = 1:10000
        disp(['Individual ', num2str(i), ': ', genotypes{i}{1}, ', ', genotypes{i}{2}]);
    end

    % determine unique genotypes
    unique_genotypes = unique_genotype_combinations(genotypes);
    disp('Unique Genotypes:');
    for i = 1:numel(unique_genotypes)
        disp(['Genotype ', num2str(i), ': ', unique_genotypes{i}{1}, ', ', unique_genotypes{i}{2}]);
    end

    % count frequency of each genotype of the 10000
    genotype_frequency = count_genotype_frequency(genotypes, unique_genotypes);
    disp('Genotype Frequency:');
    disp_genotype_frequency(genotype_frequency, 10000);

    % single out the most common genotype 
    identify_most_common_genotype(genotype_frequency);
end

function individual_genotype = generate_individual(SNP_frequencies)
    % generate two SNP's random frequencies
    tmp_freq_SNP1 = rand(1,1);
    tmp_freq_SNP2 = rand(1,1);

    % determine the genotype for SNP1 given the random number produced
    % using the first 3 numbers of SNP_frequencies
    if tmp_freq_SNP1 < SNP_frequencies(1)
        genotype_SNP1 = 'AA';
    elseif tmp_freq_SNP1 >= SNP_frequencies(3) && tmp_freq_SNP1 < SNP_frequencies(5)
        genotype_SNP1 = 'AB';
    else
        genotype_SNP1 = 'BB';
    end

    % determine the genotype for SNP2 given the random number produced
    % using the second 3 numbers of SNP_frequencies
    if tmp_freq_SNP2 < SNP_frequencies(2)
        genotype_SNP2 = 'AA';
    elseif tmp_freq_SNP2 >= SNP_frequencies(4) && tmp_freq_SNP2 < SNP_frequencies(6)
        genotype_SNP2 = 'AB';
    else
        genotype_SNP2 = 'BB';
    end

    % return an array of specific genotypes
    individual_genotype = {genotype_SNP1, genotype_SNP2};
end

function unique_genotypes = unique_genotype_combinations(genotypes)
    % unique genotype combinations identifier
    unique_genotypes = {};
    for i = 1:numel(genotypes)
        if ~any(cellfun(@(x) isequal(x, genotypes{i}), unique_genotypes))
            unique_genotypes{end+1} = genotypes{i};
        end
    end
end

function genotype_frequency = count_genotype_frequency(genotypes, unique_genotypes)
    % count the frequency of each genotype
    genotype_frequency = containers.Map();
    for i = 1:numel(unique_genotypes)
        genotype_frequency(strjoin(unique_genotypes{i}, ',')) = 0;
    end
    for i = 1:numel(genotypes)
        genotype_str = strjoin(genotypes{i}, ',');
        genotype_frequency(genotype_str) = genotype_frequency(genotype_str) + 1;
    end
end

function disp_genotype_frequency(genotype_frequency, total_individuals)
    % display genotype frequencies in an easy to understand fashion
    keys = genotype_frequency.keys;
    values = genotype_frequency.values;
    for i = 1:numel(keys)
        disp([keys{i}, ': ', num2str(values{i}), ' (', num2str(values{i}/total_individuals*100), '%)']);
    end
end

function most_common_genotype = identify_most_common_genotype(genotype_frequency)
    % identify the genotype with highest prevalence
    values = cell2mat(genotype_frequency.values);
    keys = genotype_frequency.keys;
    [~, idx] = max(values);
    most_common_genotype = keys{idx};
    disp(['Most common genotype: ', most_common_genotype, ' with frequency ', num2str(values(idx)), ' (', num2str(values(idx)/sum(values)*100), '%)']);
end