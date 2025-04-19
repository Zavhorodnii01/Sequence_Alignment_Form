import numpy as np

class GlobalSequenceAligner:

    def __init__(self, first_seq, second_seq, match_score, mismatch_score, gap_penalty):
        self.first_seq = first_seq
        self.second_seq = second_seq
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        self.__matrix = self._fill_matrix()
        self.__matrix_with_aop = self._fill_matrix_with_aop()

    def get_matrix(self):
        return self.__matrix

    def get_matrix_with_aop(self):
        return self.__matrix_with_aop

    # fill the matrix with all optimal paths
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

    def _fill_matrix(self):
        matrix = np.empty((len(self.first_seq) + 1, len(self.second_seq) + 1),dtype=object)  # Creates an m x n matrix initialized to 0

        firstNum = 0
        for i in range(len(self.second_seq) + 1):
            matrix[0][i] = (firstNum,(0,i-1))
            firstNum += self.gap_penalty

        firstNum = 0
        for i in range(len(self.first_seq) + 1):
            matrix[i][0] = (firstNum,(i-1,0))
            firstNum += self.gap_penalty

        for i in range(1, len(self.first_seq) + 1):
            for j in range(1, len(self.second_seq) + 1):
                vertical = matrix[i - 1][j][0] + self.gap_penalty
                horizontal = matrix[i][j - 1][0] + self.gap_penalty

                if self.first_seq[i - 1] == self.second_seq[j - 1]:
                    diagonal = matrix[i - 1][j - 1][0] + self.match_score
                else:
                    diagonal = matrix[i - 1][j - 1][0] + self.mismatch_score

                if vertical == max(vertical,horizontal,diagonal):
                    matrix[i][j] = (vertical,(i - 1,j))

                elif horizontal == max(horizontal,vertical,diagonal):
                    matrix[i][j] = (horizontal,(i,j - 1))

                else:
                    matrix[i][j] = (diagonal,(i - 1,j - 1))

        return matrix


    def get_paths_and_aligned_sequences_with_aop(self):
        align1 = []
        align2 = []

        def get_all_paths(path, i , j):
            path = path + [(i,j)]
            if (i,j) == (0,0):
                optimal_paths.append(path)
            else:
                for i1,j1 in self.__matrix_with_aop[i][j][1]:
                    get_all_paths(path, i1 ,j1)

        i = len(self.first_seq)
        j = len(self.second_seq)
        optimal_paths = []
        get_all_paths([], i, j)
        """trace = [(i,j)]
        for i,j in trace:
            trace_for_one = []
            while (0, 0) not in self.__matrix_with_aop[i][j][1]:
                optimal_path.append((i, j))
                trace_for_one = self.__matrix_with_aop[i][j][1]
            trace = trace_for_one"""

        """optimal_path.append((i, j))
        optimal_path.append((0, 0))
        optimal_path = optimal_path[::-1]"""

        optimal_path = optimal_paths[0]
        for i in range(len(optimal_path) - 1):
            if optimal_path[i][0] == optimal_path[i + 1][0]:
                align2.append(self.second_seq[optimal_path[i + 1][1] - 1])
                align1.append("-")
            elif optimal_path[i][1] == optimal_path[i + 1][1]:
                align1.append(self.first_seq[optimal_path[i + 1][0] - 1])
                align2.append("-")
            else:
                align2.append(self.second_seq[optimal_path[i + 1][1] - 1])
                align1.append(self.first_seq[optimal_path[i + 1][0] - 1])

        return optimal_paths, align1, align2

    def get_optimal_path_and_aligned_sequences(self):
        optimal_path = []
        align1 = []
        align2 = []

        i = len(self.__matrix) - 1
        j = len(self.__matrix[0]) - 1
        while self.__matrix[i][j][1] != (0, 0):
            optimal_path.append((i, j))
            i, j = self.__matrix[i][j][1]
        optimal_path.append((i, j))
        optimal_path.append((0, 0))
        optimal_path = optimal_path[::-1]

        for i in range(len(optimal_path) - 1):
            if optimal_path[i][0] == optimal_path[i + 1][0]:
                align2.append(self.second_seq[optimal_path[i + 1][1] - 1])
                align1.append("-")
            elif optimal_path[i][1] == optimal_path[i + 1][1]:
                align1.append(self.first_seq[optimal_path[i + 1][0] - 1])
                align2.append("-")
            else:
                align2.append(self.second_seq[optimal_path[i + 1][1] - 1])
                align1.append(self.first_seq[optimal_path[i + 1][0] - 1])

        return optimal_path, align1, align2

    def get_score(self):
        score = str(self.__matrix[len(self.first_seq)][len(self.second_seq)][0])
        return score