function aidan_holland_ribosnitches
% RB1 5' UTR Sequence
sequence = "GCUCAGUUGCCGGGCGGGGGAGGGCGCGUCCGGUUUUUCUCAGGGGACGUUGAAAUUAUUUUUGUAACGGGAGUCGGGAGAGGACGGGGCGUGCCCCGACGUGCGCGCGCGUCGUCCUCCCCGGCGCUCCUCCACAGCUCGCUGGCUCCCGCCGCGGAAAGGCGUC";
% minimum loop size
minLoopSize = 7;
% RNAfold
rna_structure=rnafold(sequence,'MinLoopSize', minLoopSize);
% plot RNAfold
rnaplot(rna_structure);;
predict_structure
% ribonucs

data = readtable('RB1_aidanph.txt', 'Delimiter', '\t');
[min_score, min_index] = min(data{:, 4});
mutated_position = data{min_index, 1};
original_nucleotide = data{min_index, 2};
mutated_nucleotide = data{min_index, 3};
disp(['Lowest score comparison: ', num2str(min_score), ', Position: ', num2str(mutated_position), ', Original Nucleotide: ', original_nucleotide,', Mutated Nucleotide: ', mutated_nucleotide]);
mutated_sequence = "GCUCGGUUGCCGGGCGGGGGAGGGCGCGUCCGGUUUUUCUCAGGGGACGUUGAAAUUAUUUUUGUAACGGGAGUCGGGAGAGGACGGGGCGUGCCCCGACGUGCGCGCGCGUCGUCCUCCCCGGCGCUCCUCCACAGCUCGCUGGCUCCCGCCGCGGAAAGGCGUC";
mutated_structure = rnafold(mutated_sequence,'MinLoopSize', minLoopSize);
rnaplot(mutated_structure);

analyze_similarity_scores(data);
end

function[score]=score_comp(str1,str2)
% score_comparison compares 2 RNA structures in dot bracket form
str1=strrep(str1,'.','A'); str1 = strrep(str1,'(','C'); str1 = strrep(str1,')','G');
str2=strrep(str2,'.','A'); str2 = strrep(str2,'(','C'); str2 = strrep(str2,')','G');
score=nwalign(str1,str2,'Alphabet','NT','ScoringMatrix',eye(4),'GAPOPEN',.5,'EXTENDGAP',0.25);
score=score / max(length(str1),length(str2));
end

function predict_structure
% input the predicted structure
predicted_structure = '(((.....((((((((.(......).)))))))).((((...(....)...))))...(((((((..((((((((.((.((((((((((....))))((((((....)))))))))))).))((((..((........))..))))))))))))..))))))))))...';
% input wildtype structures
wt1='(((...((((..((((((((..(((((((((.............))))).............(((...(((.((((((.((((((((.((((((..........)))))).))))))))))))))..)))...)))...))))..))))))))))))...)))...';
wt2='(((....(.((((((((((((((((((((((.............))))))......................((((((.((((((((.((((((..........)))))).)))))))))))))))).)))))..((((....)))))))))).))))..)))...';
wt3='........(((.(((((((((((.(((((((.............)))))).)....................((((((.((((((((.((((((..........)))))).))))))))))))))...)))))....)))))))))...((((.......))))..';
% compare wildtype structures to each other and to input predicted, with
% score output
score1=score_comp(predicted_structure, wt1);
score2=score_comp(predicted_structure, wt2);
score3=score_comp(predicted_structure, wt3); 
score4=score_comp(wt1, wt2); 
score5=score_comp(wt1, wt3); 
score6=score_comp(wt2, wt3); 
disp(['WT1/Predict Comparison Score: ', num2str(score1)]);
disp(['WT2/Predict Comparison Score: ', num2str(score2)]);
disp(['WT3/Predict Comparison Score: ', num2str(score3)]);
disp(['WT1/WT2 Comparison Score: ', num2str(score4)]);
disp(['WT1/WT3 Comparison Score: ', num2str(score5)]);
disp(['WT2/WT3 Comparison Score: ', num2str(score6)]);
end

% function ribonucs
%     sequence = 'GCUCAGUUGCCGGGCGGGGGAGGGCGCGUCCGGUUUUUCUCAGGGGACGUUGAAAUUAUUUUUGUAACGGGAGUCGGGAGAGGACGGGGCGUGCCCCGACGUGCGCGCGCGUCGUCCUCCCCGGCGCUCUUCCUCCACAGCUCGCUGGCUCCCGCCGCGGAAAGGCGUC';
%     [structure] = rnafold(sequence);
%     disp('Dot-Bracket Notation:')
%     disp(structure)
%     nucleotides = 'ACGU';
%     fileID = fopen('RB1_aidanph.txt', 'w');
%     fprintf(fileID, 'Mutated Position\tOriginal Nucleotide\tMutated Nucleotide\tComparison Score (higher is more similar)\n');
%     for i = 1:length(sequence)
%         original_nucleotide = sequence(i);
%         for j = 1:length(nucleotides)
%             if nucleotides(j) == original_nucleotide
%                 continue;
%             end
%             mutated_sequence = sequence;
%             mutated_sequence(i) = nucleotides(j);
%             mutated_structure = rnafold(mutated_sequence);
%             score = score_comp(rnafold(sequence), mutated_structure);
%             fprintf(fileID, '%d\t%s\t%s\t%f\n', i, original_nucleotide, nucleotides(j), score);
%         end
%     end
%     fclose(fileID);
%     disp('Mutation analysis completed. Results saved in RB1_onyen.txt');
% end

function analyze_similarity_scores(data)
    similarity_scores = data{:, 4};
    figure;
    histogram(similarity_scores, 'BinWidth', 0.1);
    xlabel('Similarity Score');
    ylabel('Frequency');
    title('Histogram of Similarity Scores for Induced Mutations');
    hold on;
    wild_type_score = 0.82;
    y_limits = ylim;
    line([wild_type_score, wild_type_score], y_limits, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
    hold off;
    hold on;
    wt_predict_score = 0.60;
    y_limits = ylim;
    line([wt_predict_score, wt_predict_score], y_limits, 'Color', 'g', 'LineWidth', 1.5, 'LineStyle', '--');
    hold off;
    legend('Similarity Scores', 'WT-WT Comparison Score', 'WT-Predict Comparison Score');
    fprintf('Mean similarity score: %.4f\n', mean(similarity_scores));
    fprintf('Median similarity score: %.4f\n', median(similarity_scores));
    fprintf('Standard deviation of similarity scores: %.4f\n', std(similarity_scores));
end