import numpy as np


class GlobalSequenceAligner:
    """
    A class to perform global sequence alignment using the Needleman-Wunsch algorithm.

    Attributes:
        first_seq (str): The first sequence to align.
        second_seq (str): The second sequence to align.
        match_score (int): The score given for a character match.
        mismatch_score (int): The penalty given for a character mismatch.
        gap_penalty (int): The penalty given for introducing a gap.
        __matrix (numpy.ndarray): The scoring matrix for a single optimal path.
        __matrix_with_aop (numpy.ndarray): The scoring matrix storing all possible optimal paths.
    """

    def __init__(self, first_seq, second_seq, match_score, mismatch_score, gap_penalty):
        """
        Initializes the aligner with sequences and scoring parameters.

        Args:
            first_seq (str): The first sequence to be aligned.
            second_seq (str): The second sequence to be aligned.
            match_score (int): The score for a match between characters.
            mismatch_score (int): The penalty for a mismatch between characters.
            gap_penalty (int): The penalty for introducing a gap in alignment.
        """
        self.first_seq = first_seq
        self.second_seq = second_seq
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        self.__matrix = self._fill_matrix()
        self.__matrix_with_aop = self._fill_matrix_with_aop()

    def get_matrix(self):
        """
        Get the scoring matrix for a single optimal alignment path with backtracking pointers for each cell in the matrix,.

        Returns:
            numpy.ndarray: The filled scoring matrix with backtracking pointers.
        """
        return self.__matrix

    def _fill_matrix(self):
        """
        Constructs the scoring matrix for a single optimal path using dynamic programming.

        Returns:
            numpy.ndarray: The filled scoring matrix with backtracking pointers (int, (int, int)).
        """
        matrix = np.empty((len(self.first_seq) + 1, len(self.second_seq) + 1), dtype=object)
        first_num = 0

        # Initialize first row
        for i in range(len(self.second_seq) + 1):
            matrix[0][i] = (first_num, (0, i - 1))
            first_num += self.gap_penalty

        # Initialize first column
        first_num = 0
        for i in range(len(self.first_seq) + 1):
            matrix[i][0] = (first_num, (i - 1, 0))
            first_num += self.gap_penalty

        # Fill rest of the matrix
        for i in range(1, len(self.first_seq) + 1):
            for j in range(1, len(self.second_seq) + 1):
                vertical = matrix[i - 1][j][0] + self.gap_penalty
                horizontal = matrix[i][j - 1][0] + self.gap_penalty
                if self.first_seq[i - 1] == self.second_seq[j - 1]:
                    diagonal = matrix[i - 1][j - 1][0] + self.match_score
                else:
                    diagonal = matrix[i - 1][j - 1][0] + self.mismatch_score

                if vertical == max(vertical, horizontal, diagonal):
                    matrix[i][j] = (vertical, (i - 1, j))
                elif horizontal == max(horizontal, vertical, diagonal):
                    matrix[i][j] = (horizontal, (i, j - 1))
                else:
                    matrix[i][j] = (diagonal, (i - 1, j - 1))

        return matrix

    def get_optimal_path_and_aligned_sequences(self):
        """
        Retrieve the optimal alignment path and aligned sequences.

        Returns:
            tuple: A tuple containing:
                - list: A list of coordinates for the optimal alignment path.
                - list: A list of tuples (aligned sequences).
        """
        optimal_path = []
        align1 = []
        align2 = []
        i = len(self.__matrix) - 1
        j = len(self.__matrix[0]) - 1

        # Traceback path
        while self.__matrix[i][j][1] != (0, 0):
            optimal_path.append((i, j))
            i, j = self.__matrix[i][j][1]
        optimal_path.append((i, j))
        optimal_path.append((0, 0))
        optimal_path = optimal_path[::-1]

        # Create aligned sequences
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

        return [optimal_path], [(align1, align2)]

    def get_matrix_with_aop(self):
        """
        Get the scoring matrix containing all optimal paths.

        Returns:
            numpy.ndarray: The scoring matrix with all optimal path pointers.
        """
        return self.__matrix_with_aop

    def _fill_matrix_with_aop(self):
        """
        Constructs the scoring matrix storing all optimal paths using dynamic programming.

        Returns:
            numpy.ndarray: The filled scoring matrix with multiple backtracking pointers.
        """
        matrix = np.empty((len(self.first_seq) + 1, len(self.second_seq) + 1), dtype=object)
        score = 0

        # Initialize first row
        for i in range(len(self.second_seq) + 1):
            matrix[0][i] = (score, [(0, i - 1)] if i > 0 else [(-1, 0)])
            score += self.gap_penalty

        # Initialize first column
        score = 0
        for i in range(len(self.first_seq) + 1):
            matrix[i][0] = (score, [(i - 1, 0)] if i > 0 else [(0, -1)])
            score += self.gap_penalty

        # Fill rest of the matrix
        for i in range(1, len(self.first_seq) + 1):
            for j in range(1, len(self.second_seq) + 1):
                vertical = matrix[i - 1][j][0] + self.gap_penalty
                horizontal = matrix[i][j - 1][0] + self.gap_penalty
                diagonal = matrix[i - 1][j - 1][0] + (
                    self.match_score if self.first_seq[i - 1] == self.second_seq[j - 1] else self.mismatch_score
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

    def get_paths_and_aligned_sequences_with_aop(self):
        """
        Retrieve all optimal paths and their corresponding aligned sequences.

        Returns:
            tuple: A tuple containing:
                - list: A list of all optimal paths (each path as a list of coordinates).
                - list: A list of aligned sequences for each optimal path.
        """

        def get_all_paths(path, i, j):
            """
            Recursively find all paths.

            Args:
                path (list): Current path being traced.
                i (int): Current row index.
                j (int): Current column index.
            """
            path = path + [(i, j)]
            if (i, j) == (0, 0):
                optimal_paths.append(path)
            else:
                for i1, j1 in self.__matrix_with_aop[i][j][1]:
                    get_all_paths(path, i1, j1)

        i = len(self.first_seq)
        j = len(self.second_seq)
        optimal_paths = []
        get_all_paths([], i, j)

        optimal_paths = [path[::-1] for path in optimal_paths]
        # Create aligned sequences for each path
        aligned_sequences = []
        for optimal_path in optimal_paths:
            align1 = []
            align2 = []
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
            aligned_sequences.append([align1, align2])

        return optimal_paths, aligned_sequences

    def get_score(self):
        """
        Retrieve the highest alignment score from the matrix.

        Returns:
            str: The alignment score as a string.
        """
        return str(self.__matrix[len(self.first_seq)][len(self.second_seq)][0])
