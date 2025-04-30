import os
import re
import tkinter as tk
from tkinter import filedialog, font
import pandas as pd
import numpy as np
from global_sequence_aligner import GlobalSequenceAligner

# Number of results frames
NUM_FRAMES = 0
root = tk.Tk()
root.geometry("800x700")
root.title("Sequence Alignment Form")

# Introductory Text
intro_label = tk.Label(root, text="This is an interactive example of"
    " Needleman-Wunsch algorithm used for Global Alignment",
    font=("Arial", 11), justify="center", wraplength=600)
intro_label.pack(pady=2)
message_label = tk.Label(root, text="", fg="red",
    font=("Arial", 11), justify="center", wraplength=600)
message_label.pack(pady=2)

# Set default font size globally
default_font = font.nametofont("TkDefaultFont")
default_font.configure(size=10)

# Scrollable canvas container
canvas_container = tk.Frame(root)
canvas_container.pack(fill="both", expand=True)

# Horizontal scrollbar (at top)
scrollbar_x = tk.Scrollbar(canvas_container, orient="horizontal")
scrollbar_x.pack(side="top", fill="x")

# Vertical scrollbar (at right)
scrollbar_y = tk.Scrollbar(canvas_container, orient="vertical")
scrollbar_y.pack(side="right", fill="y")

# Canvas
main_canvas = tk.Canvas(canvas_container, yscrollcommand=scrollbar_y.set,
                        xscrollcommand=scrollbar_x.set)
main_canvas.pack(side="left", fill="both", expand=True)

scrollbar_y.config(command=main_canvas.yview)
scrollbar_x.config(command=main_canvas.xview)

# Scrollable frame inside canvas
scrollable_frame = tk.Frame(main_canvas)
main_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

# Update scroll region dynamically
def on_frame_configure(event):
    """
    :param event: The event object containing information about the
     frame configuration,
     typically a resizing or layout-related event.
    :return: None
    """
    main_canvas.configure(scrollregion=main_canvas.bbox("all"))
scrollable_frame.bind("<Configure>", on_frame_configure)

# Mousewheel Support
def _on_mousewheel(event):
    """
    Handles mouse wheel events to provide scroll functionality.
    A horizontal scroll is initiated if the Shift key is held during
     the event;
    otherwise, vertical scroll is performed.

    :param event: The mouse wheel event containing information such
     as delta and state.
    :return: None
    """
    if event.state & 0x0001:  # Shift held -> horizontal scroll
        main_canvas.xview_scroll(-1 * int(event.delta / 120), "units")
    else:
        main_canvas.yview_scroll(-1 * int(event.delta / 120), "units")

def _bind_mousewheel(widget):
    """
    Binds mouse wheel events to the specified widget for scrolling
     functionality, including support for horizontal scrolling with
      Shift key and Linux-specific mouse wheel buttons.

    :param widget: The widget to bind mouse wheel events to.
    :return: None
    """
    widget.bind_all("<MouseWheel>", _on_mousewheel)
    widget.bind_all("<Shift-MouseWheel>", _on_mousewheel)
    widget.bind_all("<Button-4>", lambda e: main_canvas.yview_scroll
        (-1, "units")) # Linux scroll up
    widget.bind_all("<Button-5>", lambda e: main_canvas.yview_scroll
        (1, "units")) # Linux scroll down

_bind_mousewheel(root)

# GUI Content
main_frame = tk.Frame(scrollable_frame)
main_frame.pack(fill="both", expand=True)

# Top Frame
top_frame = tk.Frame(main_frame)
top_frame.grid(row=0, column=0, columnspan=3, pady=10)

# First Sequence entry
tk.Label(top_frame, text="Sequence 1").grid(row=0, column=0, padx=10,
                                            pady=10, sticky="e")
sequence1_entry = tk.Entry(top_frame, width=30)
sequence1_entry.grid(row=0, column=1, columnspan=2, padx=10, pady=10)

# Second Sequence entry
tk.Label(top_frame, text="Sequence 2").grid(row=1, column=0, padx=10,
                                            pady=10, sticky="e")
sequence2_entry = tk.Entry(top_frame, width=30)
sequence2_entry.grid(row=1, column=1, columnspan=2, padx=10, pady=10)

loaded_first_fasta = ""
loaded_second_fasta = ""
def load_first_fasta_file():
    """
    Loads the first FASTA file selected by the user via a file dialog
    window. If a valid file is selected, its content is read and stored
    in the global variable `loaded_first_fasta`. Additionally, updates
    a label with the file name if successfully loaded or an error message
    if the file couldn't be loaded.

    :return: None
    """
    global loaded_first_fasta
    file_path = filedialog.askopenfilename(
        filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])

    if file_path:
        with open(file_path, 'r') as file:
            fasta_content = file.read()
        loaded_first_fasta = fasta_content
        file_name = os.path.basename(file_path)
        first_fasta_label.config(text=f"Loaded: {file_name}", fg="green")
    else:
        first_fasta_label.config(text="Can't load this file", fg="red")

def load_second_fasta_file():
    """
    Prompts the user to select a second FASTA file using a file dialog.
    If a valid file is selected, its content is read and stored in the
    global variable `loaded_second_fasta`. The label `second_fasta_label` is
    updated to display the file name or an error message if the file cannot
    be loaded.

    :return: None
    """
    global loaded_second_fasta
    file_path = filedialog.askopenfilename(
        filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")]
    )
    if file_path:
        with open(file_path, 'r') as file:
            fasta_content = file.read()
        loaded_second_fasta = fasta_content
        file_name = os.path.basename(file_path)
        second_fasta_label.config(text=f"Loaded: {file_name}", fg="green")
    else:
        second_fasta_label.config(text="Can't load this file", fg="red")

def extract_sequence_from_fasta_file(loaded_fasta):
    """
    :param loaded_fasta: A string representation of a loaded FASTA file.
                         It includes header lines and sequences, where
                         valid sequences contain amino acid codes.
    :return: A concatenated string of amino acid sequences extracted
             from the input FASTA content.
    """
    extracted_sequence = ""
    for line in loaded_fasta.split("\n"):
        line = line.strip()
        if re.fullmatch(
            r'[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+', line):
            extracted_sequence += line
    return extracted_sequence

# First load button
download_button1 = tk.Button(top_frame, text="load first fasta file",
                            command=load_first_fasta_file)
download_button1.grid(row=2, column=0, pady=10)
first_fasta_label = tk.Label(top_frame, text="No file loaded", fg="gray")
first_fasta_label.grid(row=2, column=1, padx=1, sticky="w")

# Second load button
download_button2 = tk.Button(top_frame, text="load second fasta file",
                            command=load_second_fasta_file)
download_button2.grid(row=2, column=2, pady=10)
second_fasta_label = tk.Label(top_frame, text="No file loaded", fg="gray")
second_fasta_label.grid(row=2, column=3, padx=1, sticky="w")

# Scoring Labels
tk.Label(top_frame, text="Match Score").grid(row=3, column=0, padx=10,
                                             pady=5, sticky="e")
tk.Label(top_frame, text="Mismatch Score").grid(row=3, column=1, padx=10,
                                                pady=5, sticky="e")
tk.Label(top_frame, text="Gap Penalty").grid(row=3, column=2, padx=10,
                                             pady=5, sticky="e")

# Scoring Inputs
match_score_entry = tk.Entry(top_frame, width=10)
match_score_entry.insert(0, "1")
match_score_entry.grid(row=4, column=0, padx=10, pady=5)

mismatch_score_entry = tk.Entry(top_frame, width=10)
mismatch_score_entry.insert(0, "-1")
mismatch_score_entry.grid(row=4, column=1, padx=10, pady=5)

gap_penalty_entry = tk.Entry(top_frame, width=10)
gap_penalty_entry.insert(0, "-2")
gap_penalty_entry.grid(row=4, column=2, padx=10, pady=5)

# Canvas frame (for the matrix)
canvas_frame = tk.Frame(main_frame)
canvas_frame.grid(row=3, column=0, columnspan=3, pady=10)

canvas = tk.Canvas(canvas_frame, width=200, height=200)
canvas.pack()

# Submit Button to compute Optimal Alignment
submit_button = tk.Button(main_frame)
submit_button.grid(row=1, column=0, columnspan=1, pady=10)
submit_button.config(text="Compute Optimal Alignment",
                     command=lambda: compute_optimal_alignment(False))

# Submit Button to compute Optimal Alignments
submit_button = tk.Button(main_frame)
submit_button.grid(row=1, column=1, columnspan=1, pady=10)
submit_button.config(text="Compute All Alignments",
                     command=lambda: compute_optimal_alignment(True))

# Clear Button to clear resul
clear_button = tk.Button(main_frame, text="Clear Result",
                         command=lambda: clear_result(True))
clear_button.grid(row=1, column=2, columnspan=1, pady=10)

# Matrix, Optimal Paths and Widget References variables
matrix = []
optimal_paths = []
widget_references = {}
def compute_optimal_alignment(aop: bool):
    """
    Computes the optimal sequence alignment based on the provided
    input parameters and updates the relevant UI components.

    :param aop: A boolean value indicating if the matrix should
    have All Optimal Paths (AOP).
    :return: None
    """

    # Clears all previous results as well as warning messages
    message_label.config(text="")
    clear_result(False)

    # Checks the input correctness
    global loaded_first_fasta
    global loaded_second_fasta
    if not valid_user_input(sequence1_entry.get().upper().strip(), sequence2_entry.get().upper().strip(),
            extract_sequence_from_fasta_file(loaded_first_fasta), extract_sequence_from_fasta_file(loaded_second_fasta)):
        message_label.config(text="Sequences contain invalid nucleotides, try again")
    else:
        # Extra method to validate user's inputs
        def is_valid_number(s):
            return (len(s) > 0 and s.lstrip('-').isdigit() and
                    (s.count('-') <= 1 and (s.startswith('-') or '-' not in s)))
        if (is_valid_number(match_score_entry.get()) == False
                or is_valid_number(mismatch_score_entry.get()) == False
                or is_valid_number(gap_penalty_entry.get()) == False):
            message_label.config(text="Fill match, mismatch and gap penalty scores")

        # if everything is correct - start computing optimal alignment
        else:
            # Disables the user from changing the inputs until the 'Clear Result' button is pressed,
            # ensuring that the initial data cannot be changed if the user saves the results as a TXT file.
            sequence1_entry.config(state="disabled")
            sequence2_entry.config(state="disabled")
            match_score_entry.config(state="disabled")
            mismatch_score_entry.config(state="disabled")
            gap_penalty_entry.config(state="disabled")
            download_button1.config(state="disabled")
            download_button2.config(state="disabled")


            global optimal_paths
            global widget_references
            global NUM_FRAMES
            global matrix


            # Ensures user enters only one sequence for one entry
            sequence1_loaded = extract_sequence_from_fasta_file(loaded_first_fasta)
            sequence2_loaded = extract_sequence_from_fasta_file(loaded_second_fasta)
            if ((bool(sequence1_entry.get().strip()) and bool(sequence1_loaded))
                    or (not sequence1_entry.get().strip() and not sequence1_loaded)
                    or (bool(sequence2_entry.get().strip()) and bool(sequence2_loaded))
                    or (not sequence2_entry.get().strip() and not sequence2_loaded)):
                message_label.config(
                    text="You loaded two sequences or none for one entry, try again")

            else:
                sequence1 = sequence1_entry.get().upper().strip()\
                    if sequence1_entry.get() else sequence1_loaded
                sequence2 = sequence2_entry.get().upper().strip()\
                    if sequence2_entry.get() else sequence2_loaded
                match_score = int(match_score_entry.get())
                mismatch_score = int(mismatch_score_entry.get())
                gap_penalty = int(gap_penalty_entry.get())

                # Uses GlobalSequenceAligner Class to get all needed results
                aligner = GlobalSequenceAligner(sequence1, sequence2, match_score,
                                                mismatch_score, gap_penalty)
                # Checks whether the results should contain All Optimal Paths(AOP)
                if aop :
                    matrix = aligner.get_matrix_with_aop()
                    optimal_paths, aligned_sequences = (
                        aligner.get_paths_and_aligned_sequences_with_aop())
                    NUM_FRAMES = len(aligned_sequences) # equals to number of aligned sequences
                else:
                    # For one optimal path
                    NUM_FRAMES = 1
                    matrix = aligner.get_matrix()
                    optimal_paths, aligned_sequences = (
                        aligner.get_optimal_path_and_aligned_sequences())
                score = aligner.get_score()

                # Calls method to draw matrix with all alignments
                draw_matrix(optimal_paths, matrix, sequence1, sequence2, aop)
                # Creates result frames
                widget_references = create_result_frames(NUM_FRAMES)
                for idx, widgets in widget_references.items():
                    widgets['first_alignment_label'].config(
                        text=" ".join(aligned_sequences[idx][0]))
                    widgets['second_alignment_label'].config(
                        text=" ".join(aligned_sequences[idx][1]))
                    widgets['score_label'].config(text="".join(score))
                    widgets['alignment_length_label'].config(
                        text="".join(str(len(aligned_sequences[idx][0]))))
                    widgets['identity_percentage_label'].config(
                        text=f"{calculate_identety(aligned_sequences[idx][0], aligned_sequences[idx][1])}"
                             f" to {len(aligned_sequences[idx][0])} / "
                             f"{round(calculate_identety(aligned_sequences[idx][0], aligned_sequences[idx][1])
                                      / len(aligned_sequences[idx][0]) * 100, 2)}%")
                    widgets['gaps_label'].config(
                        text=f"{calculate_gaps(aligned_sequences[idx][0], aligned_sequences[idx][1])}"
                             f" to {len(aligned_sequences[idx][0])} / "
                             f"{round(calculate_gaps(aligned_sequences[idx][0], aligned_sequences[idx][1])
                                      / len(aligned_sequences[idx][0]) * 100, 2)}%")


def create_result_frames(num_frames=NUM_FRAMES):
    """
    Creates a specified number of result frames, each containing labels and widgets to
    display alignment details,
     and stores their references in a dictionary.

    :param num_frames: Number of result frames to create. Default is `NUM_FRAMES`.
    :return: A dictionary containing widget references for each created result frame,
     keyed by the frame index.
    """
    widget_references = {}  # Dictionary to store widget references by frame index

    # Create the specified number of result frames
    for idx in range(num_frames):
        result_frame = tk.Frame(main_frame)
        result_frame.grid(row=0, column=3 + idx, rowspan=3, padx=6, pady=10, sticky="nw")

        # Section title
        section_title = tk.Label(result_frame, text=f"Result {idx + 1}",
                                 font=("Arial", 12, "bold"))
        section_title.grid(row=0, column=0, columnspan=2, sticky="n", pady=(0, 10))

        # Alignment 1
        alignment1_label = tk.Label(result_frame, text="Alignment 1:")
        alignment1_label.grid(row=1, column=0, sticky="w")
        first_alignment_label = tk.Label(result_frame, text="", wraplength=300,
                                         font=("Consolas", 10, "bold"))
        first_alignment_label.grid(row=1, column=1, columnspan=2, sticky="w")

        # Alignment 2
        alignment2_label = tk.Label(result_frame, text="Alignment 2:")
        alignment2_label.grid(row=2, column=0, sticky="w")
        second_alignment_label = tk.Label(result_frame, text="", wraplength=300,
                                          font=("Consolas", 10, "bold"))
        second_alignment_label.grid(row=2, column=1, columnspan=2, sticky="w")

        # Score
        score_label_text = tk.Label(result_frame, text="Alignment Score:")
        score_label_text.grid(row=3, column=0, sticky="w")
        score_label = tk.Label(result_frame, text="", font=("Arial", 10, "bold"))
        score_label.grid(row=3, column=1, sticky="w")

        # Alignment Length
        alignment_length_label_text = tk.Label(result_frame,
                                               text="Alignment Length:")
        alignment_length_label_text.grid(row=4, column=0, sticky="w")
        alignment_length_label = tk.Label(result_frame, text="",
                                          font=("Arial", 10, "bold"))
        alignment_length_label.grid(row=4, column=1, sticky="w")

        # Identity Percentage
        identity_percentage_label_text = tk.Label(result_frame,
                                                  text="Identity num/%:")
        identity_percentage_label_text.grid(row=5, column=0, sticky="w")
        identity_percentage_label = tk.Label(result_frame, text="",
                                             font=("Arial", 10, "bold"))
        identity_percentage_label.grid(row=5, column=1, sticky="w")

        # Gaps
        gaps_label_text = tk.Label(result_frame, text="Gaps num/%:")
        gaps_label_text.grid(row=6, column=0, sticky="w")
        gaps_label = tk.Label(result_frame, text="", font=("Arial", 10, "bold"))
        gaps_label.grid(row=6, column=1, sticky="w")

        # Button to save results as .txt
        save_button = tk.Button(result_frame, text="Save result as .txt",
                                command=lambda i=idx: save_result_to_text_file(i))
        save_button.grid(row=7, column=0, columnspan=2, pady=5)

        # Store widget references by index
        widget_references[idx] = {
            'result_frame': result_frame,
            'section_title': section_title,
            'first_alignment_label': first_alignment_label,
            'second_alignment_label': second_alignment_label,
            'score_label': score_label,
            'alignment_length_label': alignment_length_label,
            'identity_percentage_label': identity_percentage_label,
            'gaps_label': gaps_label,
            'save_button': save_button
        }

    return widget_references  # Return the dictionary of widget references by index

def clear_result(with_scores : bool):
    """
    Clears all result frames, resets canvas dimensions and content,
     and optionally sets default score values.

    :param with_scores: A boolean flag to determine whether to reset
     scoring-related entry fields to default values.
    :return: None
    """
    # Enables user to change inputs after clearing a result
    sequence1_entry.config(state="normal")
    sequence2_entry.config(state="normal")
    match_score_entry.config(state="normal")
    mismatch_score_entry.config(state="normal")
    gap_penalty_entry.config(state="normal")
    download_button1.config(state="normal")
    download_button2.config(state="normal")

    for idx, widgets in widget_references.items():
        result_frame = widgets['result_frame']
        result_frame.destroy()  # Destroy the result frame

    widget_references.clear()
    canvas.config(width=0,height=0, bg="white")
    canvas.delete("all")
    if with_scores:
        match_score_entry.delete(0, tk.END)  # Deletes old value
        match_score_entry.insert(0, "1")  # Inserts new value
        mismatch_score_entry.delete(0, tk.END)  # Deletes old value
        mismatch_score_entry.insert(0, "-1")  # Inserts new value
        gap_penalty_entry.delete(0, tk.END)  # Deletes old value
        gap_penalty_entry.insert(0, "-2")  # Inserts new value

    # Removes loaded files
    loaded_first_fasta = ""
    loaded_second_fasta = ""
    first_fasta_label.config(text="No file loaded", fg="gray")
    second_fasta_label.config(text="No file loaded", fg="gray")


# Arrow images in .png format for visualization (diagonal, top, left)
arrow_left_img = tk.PhotoImage(file="left_arrow.png")
arrow_top_img = tk.PhotoImage(file="top_arrow.png")
arrow_diag_img = tk.PhotoImage(file="diagonal_arrow.png")
# Maps direction vectors to corresponding arrow images for visualization
# (diagonal, top, left)
arrow_map = {
    (-1, -1): arrow_diag_img,  # diagonal
    (-1, 0): arrow_top_img,  # from top
    (0, -1): arrow_left_img,  # from left
}

def draw_matrix(optimal_paths: list, matrix, first_seq, second_seq,
                draw_with_aop):
    """
    Draws the specified matrix, with the specified alignment,
    and optionally, with all alignments

    :param optimal_paths: List of lists or tuples representing the
    optimal paths to be visualized on the matrix.
    :param matrix: Two-dimensional list or array where each element contains
     a tuple representing a value and its corresponding directions
    or links to other cells in the matrix.
    :param first_seq: The first sequence (typically a string or a list)
    used to label the rows of the matrix.
    :param second_seq: The second sequence (typically a string or a list)
    used to label the columns of the matrix.
    :param draw_with_aop: Boolean flag indicating whether to draw all
    available optimal paths (if True) or just one (if False).
    :return: None
    """

    # Configures canvas for matrix
    canvas.config(width=50 * len(matrix[0])+50, height=50 * len(matrix) + 50, bg="white")

    # Fills the first column of the matrix
    for i in range(len(matrix[0]) + 1):
        if i < 2:
            canvas.create_rectangle(i * 50, 0, i * 50 + 50, 50,fill="blue")
        else:
            canvas.create_rectangle(i * 50, 0, i * 50 + 50, 50, fill="lightblue")
            if i - 2 < len(second_seq):
                canvas.create_text( i * 50 + 25, 25, text=second_seq[i - 2],
                                    font=("Arial", 19), fill="#335b00")

    # Fills the second column of the matrix
    for j in range(len(matrix) + 1):
        if j < 2:
            canvas.create_rectangle(0, j * 50, 50, j * 50 + 50, fill="blue")
        else:
            canvas.create_rectangle(0, j * 50, 50, j * 50 + 50, fill="lightblue")
            if j - 2 < len(first_seq):
                canvas.create_text( 25, j * 50 + 25 , text=first_seq[j-2],
                                    font=("Arial", 19), fill="#335b00")

    # Creates empty grey rectangles
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            (canvas.create_rectangle
             (j * 50 + 50, i * 50 + 50, 50 + j * 50 + 50, 50 + i * 50 + 50, fill="lightgrey",
                                    outline="black"))
    # Checks whether matrix should have All Optimal Alignments (AOP)
    if draw_with_aop:
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                for optimal_path in optimal_paths:
                    if optimal_path.__contains__((i, j)):
                        # Creates red rectangles indicating alignment
                        (canvas.create_rectangle
                         (j * 50 + 50, i * 50 + 50, 50 + j * 50 + 50, 50 + i * 50 + 50,
                                                fill="red", outline="black"))
                # Compares current coordinates with backtrack coordinates to determine arrow direction
                for ii, jj in matrix[i][j][1]:
                    dif_i = ii - i
                    dif_j = jj - j
                    if arrow_map.__contains__((dif_i, dif_j)):
                        # Draws appropriate arrow on matrix
                        (canvas.create_image
                         (j * 50 + 50, i * 50 + 50, image=arrow_map[(dif_i, dif_j)], anchor="nw"))
                    canvas.create_text(j * 50 + 25 + 50, i * 50 + 25 + 50, text=matrix[i][j][0],
                                       font=("Arial", 16, "bold"), fill="black")
    else:
        optimal_path = optimal_paths[0]
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if optimal_path.__contains__((i, j)):
                    # Creates red rectangles indicating alignment
                    (canvas.create_rectangle
                     (j * 50 + 50, i * 50 + 50, 50 + j * 50 + 50, 50 + i * 50 + 50, fill="red",
                                            outline="black"))
                dif_i = matrix[i][j][1][0] - i
                dif_j = matrix[i][j][1][1] - j
                if arrow_map.__contains__((dif_i, dif_j)):
                    # Draws appropriate arrow on matrix
                    (canvas.create_image
                     (j * 50 + 50, i * 50 + 50, image=arrow_map[(dif_i, dif_j)], anchor="nw"))
                (canvas.create_text
                 (j * 50 + 25 + 50, i * 50 + 25 + 50, text=matrix[i][j][0], font=("Arial", 16, "bold"),
                                   fill="black"))

    canvas.create_rectangle(50, 50, 100, 100, fill="lightgreen", outline="black")
    canvas.create_text(50 + 25, 50 + 25, text="0",font=("Arial", 16, "bold"), fill="black")


def valid_user_input(seq1, seq2, seq1_loaded, seq2_loaded):
    """
    Validates user input by checking whether it is a valid sequence.

    :param seq1: A string representing the first amino acid sequence.
    :param seq2: A string representing the second amino acid sequence.

    :return: Returns True if both sequences contain only valid amino acid characters
     (ACDEFGHIKLMNPQRSTVWY, case-insensitive), otherwise returns False.
    """
    return True if re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+', seq1 + seq2 + seq1_loaded + seq2_loaded) else False

def save_result_to_text_file(idx):
    """
    Saves the alignment results and matrix visualization to a specified text file.
     If no alignment has been computed, an error message is displayed.

    :param idx: Index of the alignment for which the results are to be saved to
     frame with index (widget_references.get(idx)).
     Used to retrieve the relevant optimal path and widget references.
    :return: None
    """


    global optimal_paths
    optimal_path = optimal_paths[idx]
    global widget_references
    widget = widget_references.get(idx)

    sequence1_loaded = extract_sequence_from_fasta_file(loaded_first_fasta)
    sequence2_loaded = extract_sequence_from_fasta_file(loaded_second_fasta)

    sequence1 = sequence1_entry.get().upper().strip() \
        if sequence1_entry.get() else sequence1_loaded
    sequence2 = sequence2_entry.get().upper().strip() \
        if sequence2_entry.get() else sequence2_loaded
    print(sequence1)

    # Ensures all required data has been entered by the user
    if len(widget['first_alignment_label'].cget("text")) == 0:
        message_label.config(text="Press Compute Alignment Button!", fg="red")
    else:
        # Creates visual matrix
        matrix_txt = np.empty(
            (len(sequence1) + 2, len(sequence2) + 2),
            dtype=object)

        for i in range(len(sequence1) + 2):
            if i > 1:
                matrix_txt[i][0] = sequence1[i - 2]

        for j in range(len(sequence2) + 2):
            if j > 1:
                matrix_txt[0][j] = sequence2[j - 2]

        # Fills headers
        matrix_txt[0][0] = "."
        matrix_txt[1][0] = "."
        matrix_txt[0][1] = "."

        arrow_map = {
            (-1, -1): "↖",  # diagonal
            (-1, 0): "↑",  # from top
            (0, -1): "←",  # from left
        }
        # Fills matrix_txt with score + arrow
        for i in range(1, len(sequence1) + 2):
            for j in range(1, len(sequence2) + 2):
                arrow = ""
                if optimal_path.__contains__((i - 1, j - 1)):
                    index = optimal_path.index((i - 1, j - 1))
                    (next_i, next_j) = optimal_path[index - 1]
                    # Calculates direction vector
                    di = next_i - (i - 1)
                    dj = next_j - (j - 1)
                    arrow = arrow_map.get((di, dj), " ")  # fallback if not matched
                matrix_txt[i][j] = f"{matrix[i - 1][j - 1][0]}{arrow}"

        # Formats the matrix using pandas
        matrix_txt = np.where(matrix_txt == None, ".", matrix_txt)
        df = pd.DataFrame(matrix_txt)

        # Saves data to the selected file path
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if file_path:
            with open(file_path, "w", encoding="utf-8") as file:
                file.write("Sequence 1:\n")
                file.write(f"{sequence1}\n")
                file.write("Sequence 2:\n")
                file.write(f"{sequence2}\n")
                file.write(f"Match Score:   {match_score_entry.get()}\n")
                file.write(f"Mismatch Score:   {mismatch_score_entry.get()}\n")
                file.write(f"Gap Penalty:   {gap_penalty_entry.get()}\n")
                file.write(f"Alignment 1:   {widget['first_alignment_label'].cget("text")}\n")
                file.write(f"Alignment 2:   {widget['second_alignment_label'].cget("text")}\n")
                file.write(f"Score:   {widget['score_label'].cget("text")}\n")
                file.write(f"Alignment Length:  {widget['alignment_length_label'].cget("text")}\n")
                file.write(f"Identity  num/%:    {widget['identity_percentage_label'].cget("text")}\n")
                file.write(f"Gaps  num/%:    {widget['gaps_label'].cget("text")}\n")
                file.write(df.to_string(index=False, header=False))

def calculate_gaps(aligned_seq1, aligned_seq2):
    """
    Calculates the number of gaps in the alignment.

    :param aligned_seq1: The first aligned sequence containing gaps represented as '-'.
    :param aligned_seq2: The second aligned sequence containing gaps represented as '-'.

    :return: The total number of gaps ('-') present in both aligned sequences combined.
    """
    gaps = aligned_seq1.count("-") + aligned_seq2.count("-")
    return gaps

def calculate_identety(aligned_seq1, aligned_seq2):
    """
    Calculates the number of identity in the alignment.

    :param aligned_seq1: The first aligned sequence to compare.
    :param aligned_seq2: The second aligned sequence to compare.

    :return: The count of identical positions between the two aligned sequences,
     excluding gaps ("-").
    """
    identical_positions = 0
    for i in range(len(aligned_seq1)):
        if (aligned_seq1[i] == aligned_seq2[i]) and aligned_seq1 != "-":
            identical_positions += 1
    return identical_positions

root.mainloop()