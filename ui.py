import os
import re
import tkinter as tk
from tkinter import filedialog, font, messagebox

import pandas as pd
import numpy as np
from global_sequence_aligner import GlobalSequenceAligner

root = tk.Tk()
root.geometry("700x600")
root.title("Sequence Alignment Form")

# === Introductory Text ===
intro_label = tk.Label(root, text="This is an interactive example of Needleman-Wunsch"
    " algorithm used for Global Alignment",
    font=("Arial", 11), justify="center", wraplength=600)
intro_label.pack(pady=2)
message_label = tk.Label(root, text="", fg="red",
    font=("Arial", 11), justify="center", wraplength=600)
message_label.pack(pady=2)



# Set default font size globally
default_font = font.nametofont("TkDefaultFont")
default_font.configure(size=10)

# === Scrollable canvas container ===
canvas_container = tk.Frame(root)
canvas_container.pack(fill="both", expand=True)

# === Horizontal scrollbar (at top) ===
scrollbar_x = tk.Scrollbar(canvas_container, orient="horizontal")
scrollbar_x.pack(side="top", fill="x")

# === Vertical scrollbar (at right) ===
scrollbar_y = tk.Scrollbar(canvas_container, orient="vertical")
scrollbar_y.pack(side="right", fill="y")

# === Canvas ===
main_canvas = tk.Canvas(canvas_container, yscrollcommand=scrollbar_y.set, xscrollcommand=scrollbar_x.set)
main_canvas.pack(side="left", fill="both", expand=True)

scrollbar_y.config(command=main_canvas.yview)
scrollbar_x.config(command=main_canvas.xview)

# === Scrollable frame inside canvas ===
scrollable_frame = tk.Frame(main_canvas)
main_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

# === Update scroll region dynamically ===
def on_frame_configure(event):
    main_canvas.configure(scrollregion=main_canvas.bbox("all"))

scrollable_frame.bind("<Configure>", on_frame_configure)
# === Mousewheel Support ===
def _on_mousewheel(event):
    if event.state & 0x0001:  # Shift held -> horizontal scroll
        main_canvas.xview_scroll(-1 * int(event.delta / 120), "units")
    else:
        main_canvas.yview_scroll(-1 * int(event.delta / 120), "units")

def _bind_mousewheel(widget):
    widget.bind_all("<MouseWheel>", _on_mousewheel)
    widget.bind_all("<Shift-MouseWheel>", _on_mousewheel)
    widget.bind_all("<Button-4>", lambda e: main_canvas.yview_scroll(-1, "units"))  # Linux scroll up
    widget.bind_all("<Button-5>", lambda e: main_canvas.yview_scroll(1, "units"))   # Linux scroll down

_bind_mousewheel(root)

# === GUI Content ===
main_frame = tk.Frame(scrollable_frame)
main_frame.pack(fill="both", expand=True)

top_frame = tk.Frame(main_frame)
top_frame.grid(row=0, column=0, columnspan=3, pady=10)

tk.Label(top_frame, text="Sequence 1").grid(row=0, column=0, padx=10, pady=10, sticky="e")
sequence1_entry = tk.Entry(top_frame, width=30)
sequence1_entry.grid(row=0, column=1, columnspan=2, padx=10, pady=10)

tk.Label(top_frame, text="Sequence 2").grid(row=1, column=0, padx=10, pady=10, sticky="e")
sequence2_entry = tk.Entry(top_frame, width=30)
sequence2_entry.grid(row=1, column=1, columnspan=2, padx=10, pady=10)

loaded_first_fasta = ""
loaded_second_fasta = ""


def valid_user_input(seq1, seq2):
    return True if re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+', seq1 + seq2) else False


def extract_sequence_from_fasta_file(loaded_fasta):
    extracted_sequence = ""
    for line in loaded_fasta.split("\n"):
        line = line.strip()
        if re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+', line):
            extracted_sequence += line
    return extracted_sequence


def load_first_fasta_file():
    global loaded_first_fasta
    file_path = filedialog.askopenfilename(
        filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")]
    )
    if file_path:
        with open(file_path, 'r') as file:
            fasta_content = file.read()
        loaded_first_fasta = fasta_content

        file_name = os.path.basename(file_path)
        first_fasta_label.config(text=f"Loaded: {file_name}", fg="green")
    else:
        first_fasta_label.config(text="Can't load this file", fg="red")

def load_second_fasta_file():
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

def save_result_to_text_file():
    if len(first_alignment_label.cget("text")) == 0:
        message_label.config(text="Press Compute Alignment Button!", fg="red")
    else:
        matrix_txt = np.empty(
            (len(sequence1_entry.get().strip()) + 2, len(sequence2_entry.get().strip()) + 2), dtype=object)

        for i in range(len(sequence1_entry.get().strip()) + 2):
            if i > 1:
                matrix_txt[i][0] = sequence1_entry.get().upper().strip()[i - 2]

        for j in range (len(sequence2_entry.get().strip()) + 2):
            if j > 1:
                matrix_txt[0][j] = sequence2_entry.get().upper().strip()[j - 2]

        arrow_map = {
            (-1, -1): "↖",  # diagonal
            (-1, 0): "↑",  # from top
            (0, -1): "←",  # from left
        }

        # Create visual matrix
        matrix_txt = np.empty(
            (len(sequence1_entry.get().strip()) + 2, len(sequence2_entry.get().strip()) + 2),
            dtype=object
        )

        # Fill headers
        matrix_txt[0][0] = "."
        matrix_txt[1][0] = "."
        matrix_txt[0][1] = "."

        for i, ch in enumerate(sequence1_entry.get().upper().strip()):
            matrix_txt[i + 2][0] = ch
        for j, ch in enumerate(sequence2_entry.get().upper().strip()):
            matrix_txt[0][j + 2] = ch

        # Now fill matrix_txt with score + arrow
        for i in range(1, len(sequence1_entry.get().strip()) + 2):
            for j in range(1, len(sequence2_entry.get().strip()) + 2):
                score, (prev_i, prev_j) = matrix[i - 1][j - 1]

                # Calculate direction vector
                di = prev_i - (i - 1)
                dj = prev_j - (j - 1)

                arrow = arrow_map.get((di, dj), " ")  # fallback if not matched
                matrix_txt[i][j] = f"{score}{arrow}"

        matrix_txt = np.where(matrix_txt == None, ".", matrix_txt)
        df = pd.DataFrame(matrix_txt)

        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if file_path:
            with open(file_path, "w", encoding="utf-8") as file:
                file.write("Sequence 1:\n")
                file.write(f"{sequence1_entry.get().upper().strip()}\n")
                file.write("Sequence 2:\n")
                file.write(f"{sequence2_entry.get().upper().strip()}\n")
                file.write(f"Match Score:   {match_score_entry.get()}\n")
                file.write(f"Mismatch Score:   {mismatch_score_entry.get()}\n")
                file.write(f"Gap Penalty:   {gap_penalty_entry.get()}\n")
                file.write(f"Alignment 1:   {first_alignment_label.cget('text')}\n")
                file.write(f"Alignment 2:   {second_alignment_label.cget('text')}\n")
                file.write(f"Score:   {score_label.cget('text')}\n")
                file.write(f"Alignment Length:  {alignment_length_label.cget('text')}\n")
                file.write(f"Identity  num/%:    {identity_percentage_label.cget('text')}\n")
                file.write(f"Gaps  num/%:    {gaps_label.cget('text')}\n")
                file.write(df.to_string(index=False, header=False))


#row 2
download_button = tk.Button(top_frame, text="load first fasta file", command=load_first_fasta_file)
download_button.grid(row=2, column=0, pady=10)
first_fasta_label = tk.Label(top_frame, text="No file loaded", fg="gray")
first_fasta_label.grid(row=2, column=1, padx=1, sticky="w")

download_button = tk.Button(top_frame, text="load second fasta file", command=load_second_fasta_file)
download_button.grid(row=2, column=2, pady=10)
second_fasta_label = tk.Label(top_frame, text="No file loaded", fg="gray")
second_fasta_label.grid(row=2, column=3, padx=1, sticky="w")





# ==== Row 2: Scoring Labels ====
tk.Label(top_frame, text="Match Score").grid(row=3, column=0, padx=10, pady=5, sticky="e")
tk.Label(top_frame, text="Mismatch Score").grid(row=3, column=1, padx=10, pady=5, sticky="e")
tk.Label(top_frame, text="Gap Penalty").grid(row=3, column=2, padx=10, pady=5, sticky="e")

# ==== Row 3: Scoring Inputs ====
match_score_entry = tk.Entry(top_frame, width=10)
match_score_entry.insert(0, "1")
match_score_entry.grid(row=4, column=0, padx=10, pady=5)

mismatch_score_entry = tk.Entry(top_frame, width=10)
mismatch_score_entry.insert(0, "-1")
mismatch_score_entry.grid(row=4, column=1, padx=10, pady=5)

gap_penalty_entry = tk.Entry(top_frame, width=10)
gap_penalty_entry.insert(0, "-2")
gap_penalty_entry.grid(row=4, column=2, padx=10, pady=5)


# === Canvas frame (for the matrix) ===
canvas_frame = tk.Frame(main_frame)
canvas_frame.grid(row=3, column=0, columnspan=3, pady=10)

canvas = tk.Canvas(canvas_frame, width=200, height=200)
canvas.pack()

matrix = []
# ==== Compute alignment button ====
def compute_optimal_alignment(aop: bool):
    message_label.config(text="")

    if valid_user_input(sequence1_entry.get().upper().strip(), sequence2_entry.get().upper().strip()) == False:
        message_label.config(text="Sequences contain invalid nucleotides, try again")
    else:
        global matrix
        global loaded_first_fasta
        global loaded_second_fasta
        sequence1_loaded = extract_sequence_from_fasta_file(loaded_first_fasta)
        sequence2_loaded = extract_sequence_from_fasta_file(loaded_second_fasta)
        if ((bool(sequence1_entry.get().strip()) and bool(sequence1_loaded)) or (not sequence1_entry.get().strip() and not sequence1_loaded)
                or (bool(sequence2_entry.get().strip()) and bool(sequence2_loaded)) or (not sequence2_entry.get().strip() and not sequence2_loaded)):
            message_label.config(text="You loaded two sequences or none for one entry, try again")


        else:
            sequence1 = sequence1_entry.get().upper().strip() if sequence1_entry.get() else sequence1_loaded
            sequence2 = sequence2_entry.get().upper().strip() if sequence2_entry.get() else sequence2_loaded
            match_score = int(match_score_entry.get())
            mismatch_score = int(mismatch_score_entry.get())
            gap_penalty = int(gap_penalty_entry.get())

            aligner = GlobalSequenceAligner(sequence1, sequence2, match_score, mismatch_score, gap_penalty)
            if aop :
                matrix = aligner.get_matrix_with_aop()
                optimal_path, first_alignment, second_alignment = aligner.get_paths_and_aligned_sequences_with_aop()
                draw_matrix(optimal_path, matrix, sequence1, sequence2, True)
            else:
                matrix = aligner.get_matrix()
                optimal_path, first_alignment, second_alignment = aligner.get_optimal_path_and_aligned_sequences()
                draw_matrix(optimal_path, matrix, sequence1, sequence2, False)
            #chnages
            score = aligner.get_score()
            score_label.config(text="".join(score))
            second_alignment_label.config(text=" ".join(list(second_alignment)))
            first_alignment_label.config(text=" ".join(list(first_alignment)))
            alignment_length_label.config(text="".join(str(len(first_alignment))))
            identity_percentage_label.config(
                text=f"{calculate_identety(first_alignment, second_alignment)} to {len(first_alignment)} / "
                     f"{round(calculate_identety(first_alignment, second_alignment) / len(first_alignment) * 100, 2)}%"
            )
            gaps_label.config(text=f"{calculate_gaps(first_alignment, second_alignment)} to {len(first_alignment)} / {round(calculate_gaps(first_alignment, second_alignment) / len(first_alignment) * 100, 2)}%")



        loaded_first_fasta = ""
        loaded_second_fasta = ""
        first_fasta_label.config(text="No file loaded", fg="gray")
        second_fasta_label.config(text="No file loaded", fg="gray")



def clear_result():
    first_alignment_label.config(text="")
    second_alignment_label.config(text="")
    score_label.config(text="")
    alignment_length_label.config(text="")
    identity_percentage_label.config(text="")
    gaps_label.config(text="")
    canvas.delete("all")
    match_score_entry.delete(0, tk.END)  # Deletes old value
    match_score_entry.insert(0, "1")  # Inserts new value
    mismatch_score_entry.delete(0, tk.END)  # Deletes old value
    mismatch_score_entry.insert(0, "-1")  # Inserts new value
    gap_penalty_entry.delete(0, tk.END)  # Deletes old value
    gap_penalty_entry.insert(0, "-2")  # Inserts new value

def calculate_gaps(aligned_seq1, aligned_seq2):
    gaps = aligned_seq1.count("-") + aligned_seq2.count("-")
    return gaps

def calculate_identety(aligned_seq1, aligned_seq2):
    identical_positions = 0
    for i in range(len(aligned_seq1)):
        if (aligned_seq1[i] == aligned_seq2[i]) and aligned_seq1 != "-":
            identical_positions += 1
    return identical_positions

arrow_left_img = tk.PhotoImage(file="left_arrow.png")
arrow_top_img = tk.PhotoImage(file="top_arrow.png")
arrow_diag_img = tk.PhotoImage(file="diagonal_arrow.png")
arrow_map = {
    (-1, -1): arrow_diag_img,  # diagonal
    (-1, 0): arrow_top_img,  # from top
    (0, -1): arrow_left_img,  # from left
}

def draw_matrix(optimal_paths: list, matrix, first_seq, second_seq, draw_with_aop):

    canvas.config(width=50 * len(matrix[0])+50, height=50 * len(matrix) + 50, bg="white")

    for i in range(len(matrix[0]) + 1):
        if i < 2:
            canvas.create_rectangle(i * 50, 0, i * 50 + 50, 50,fill="blue")
        else:
            canvas.create_rectangle(i * 50, 0, i * 50 + 50, 50, fill="lightblue")
            if i - 2 < len(second_seq):
                canvas.create_text( i * 50 + 25, 25, text=second_seq[i - 2], font=("Arial", 19), fill="#335b00")

    for j in range(len(matrix) + 1):
        if j < 2:
            canvas.create_rectangle(0, j * 50, 50, j * 50 + 50, fill="blue")
        else:
            canvas.create_rectangle(0, j * 50, 50, j * 50 + 50, fill="lightblue")
            if j - 2 < len(first_seq):
                canvas.create_text( 25, j * 50 + 25 , text=first_seq[j-2], font=("Arial", 19), fill="#335b00")

    if draw_with_aop:
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):

                    canvas.create_rectangle(j * 50+50 , i * 50+50, 50 + j * 50+50, 50 + i * 50+50, fill="lightgrey", outline="black")
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                for optimal_path in optimal_paths:
                    if optimal_path.__contains__((i, j)):
                        canvas.create_rectangle(j * 50 + 50, i * 50 + 50, 50 + j * 50 + 50, 50 + i * 50 + 50,
                                                fill="red", outline="black")
                for ii, jj in matrix[i][j][1]:
                    dif_i = ii - i
                    dif_j = jj - j
                    if arrow_map.__contains__((dif_i, dif_j)):
                        canvas.create_image(j * 50 + 50, i * 50 + 50, image=arrow_map[(dif_i, dif_j)], anchor="nw")
                    canvas.create_text(j * 50 + 25 + 50, i * 50 + 25 + 50, text=matrix[i][j][0],
                                       font=("Arial", 16, "bold"), fill="black")

    else:
        optimal_path = optimal_paths
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if optimal_path.__contains__((i,j)):
                    canvas.create_rectangle(j * 50 + 50 , i * 50+50, 50 + j * 50+50, 50 + i * 50+50, fill="red", outline="black")
                else:
                    canvas.create_rectangle(j * 50+50 , i * 50+50, 50 + j * 50+50, 50 + i * 50+50, fill="lightgrey", outline="black")

                dif_i = matrix[i][j][1][0] - i
                dif_j = matrix[i][j][1][1] - j
                if arrow_map.__contains__((dif_i, dif_j)):
                    canvas.create_image(j * 50 + 50, i * 50 + 50, image=arrow_map[(dif_i, dif_j)], anchor="nw")
                canvas.create_text(j * 50 + 25 + 50, i * 50 + 25 + 50 , text=matrix[i][j][0], font=("Arial", 16, "bold"), fill="black")

    canvas.create_rectangle(50, 50, 100, 100, fill="lightgreen", outline="black")
    canvas.create_text(50 + 25, 50 + 25, text="0",font=("Arial", 16, "bold"), fill="black")



submit_button = tk.Button(main_frame)
submit_button.grid(row=1, column=0, columnspan=1, pady=10)
submit_button.config(text="Compute Optimal Alignment", command=lambda: compute_optimal_alignment(False))


submit_button = tk.Button(main_frame)
submit_button.grid(row=1, column=1, columnspan=1, pady=10)
submit_button.config(text="Compute All Alignments", command=lambda: compute_optimal_alignment(True))


clear_button = tk.Button(main_frame)
clear_button.grid(row=1, column=2, columnspan=1, pady=10)
clear_button.config(text="Clear Result", command=clear_result)



# ==== Result Section ====
result_frame = tk.Frame(main_frame)
result_frame.grid(row=0, column=3, rowspan=3, padx=6, pady=10, sticky="nw")

# Section title
tk.Label(result_frame, text="Result", font=("Arial", 12, "bold")).grid(row=0, column=0, columnspan=2, sticky="n", pady=(0, 10))
#
# Alignment 1
tk.Label(result_frame, text="Alignment 1:").grid(row=1, column=0, sticky="w")
first_alignment_label = tk.Label(result_frame, text="", wraplength=300, font=("Consolas", 10, "bold"))
first_alignment_label.grid(row=1, column=1, columnspan=2, sticky="w")

# Alignment 2
tk.Label(result_frame, text="Alignment 2:").grid(row=2, column=0, sticky="w")
second_alignment_label = tk.Label(result_frame, text="", wraplength=300, font=("Consolas", 10, "bold"))
second_alignment_label.grid(row=2, column=1, columnspan=2, sticky="w")

# Score
tk.Label(result_frame, text="Alignment Score:").grid(row=3, column=0, sticky="w")
score_label = tk.Label(result_frame, text="", font=("Arial", 10, "bold"))
score_label.grid(row=3, column=1, sticky="w")

# Alignment Length
tk.Label(result_frame, text="Alignment Length:").grid(row=4, column=0, sticky="w")
alignment_length_label = tk.Label(result_frame, text="", font=("Arial", 10, "bold"))
alignment_length_label.grid(row=4, column=1, sticky="w")

# Identity Percentage
tk.Label(result_frame, text="Identity num/%:").grid(row=5, column=0, sticky="w")
identity_percentage_label = tk.Label(result_frame, text="", font=("Arial", 10, "bold"))
identity_percentage_label.grid(row=5, column=1, sticky="w")

# Gaps
tk.Label(result_frame, text="Gaps num/%:").grid(row=6, column=0, sticky="w")
gaps_label = tk.Label(result_frame, text="", font=("Arial", 10, "bold"))
gaps_label.grid(row=6, column=1, sticky="w")




save_button = tk.Button(result_frame, text="Save result as .txt", command=save_result_to_text_file)
save_button.grid(row=8, column=0, columnspan=2, pady=5)





optimal_path_label = tk.Label(result_frame, text="")

root.mainloop()
