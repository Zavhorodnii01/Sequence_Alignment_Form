import numpy as np

class GlobalSequenceAligner2:
    def __init__(self, seq1, seq2, match, mismatch, gap):
        self.first_seq = seq1
        self.second_seq = seq2
        self.match_score = match
        self.mismatch_score = mismatch
        self.gap_penalty = gap

    def get_matrix_with_aop(self):
        return self._fill_matrix_with_aop()

    def _fill_matrix_with_aop(self):
        matrix = np.empty((len(self.first_seq) + 1, len(self.second_seq) + 1), dtype=object)

        score = 0
        for i in range(len(self.second_seq) + 1):
            matrix[0][i] = (score, [(0, i - 1)] if i > 0 else [(-1, 0)])
            score += self.gap_penalty

        score = 0
        for i in range(len(self.first_seq) + 1):
            matrix[i][0] = (score, [(i - 1, 0)] if i > 0 else [(0, -1)])
            score += self.gap_penalty

        for i in range(1, len(self.first_seq) + 1):
            for j in range(1, len(self.second_seq) + 1):
                vertical = matrix[i - 1][j][0] + self.gap_penalty
                horizontal = matrix[i][j - 1][0] + self.gap_penalty
                diagonal = matrix[i - 1][j - 1][0] + (
                    self.match_score if self.first_seq[i - 1] == self.second_seq[j - 1]
                    else self.mismatch_score
                )

                max_score = max(vertical, horizontal, diagonal)
                backtrack = []
                if vertical == max_score:
                    backtrack.append((i - 1, j))
                if horizontal == max_score:
                    backtrack.append((i, j - 1))
                if diagonal == max_score:
                    backtrack.append((i - 1, j - 1))

                matrix[i][j] = (max_score, backtrack)
        return matrix

# Run it
aligner = GlobalSequenceAligner2("GATT", "GCTT", 1, -1, -1)
matrix = aligner.get_matrix_with_aop()
print(matrix)
"""# Print nicely
for row in matrix:
    print([f"{cell[0]} {cell[1]}" for cell in row])"""



