from global_sequence_aligner import GlobalSequenceAligner

aligner = GlobalSequenceAligner("GATTACA", "GTCGACGCA", 1, -1, -1)

print()

print(aligner.get_paths_and_aligned_sequences_with_aop())
#print(aligner.get_optimal_path_and_aligned_sequences_with_aop())
def print_alignment_matrix(matrix):
    for row in matrix:
        formatted_row = []
        for score, parents in row:
            p_str = ",".join(f"{i},{j}" for i, j in parents)
            formatted_row.append(f"{score:>3} [{p_str}]")
        print(" | ".join(formatted_row))
    print()

print_alignment_matrix(aligner.get_matrix_with_aop())
