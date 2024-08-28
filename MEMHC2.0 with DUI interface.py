
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText
import pandas as pd
from mhcflurry import Class1AffinityPredictor
import numpy as np
from PIL import Image, ImageTk
import threading
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os

# Initialize global variables for real-time saving
session_save_path = "realtime_results.csv"
session_df = pd.DataFrame()

def log_message(message):
    log_text.insert(tk.END, message + "\n")
    log_text.see(tk.END)

def save_realtime_results(df):
    """Save the results in real-time to avoid data loss."""
    global session_df
    session_df = pd.concat([session_df, df], ignore_index=True)
    session_df.to_csv(session_save_path, index=False)

def generate_peptides(protein_sequence, min_length, max_length):
    peptides = []
    start_positions = []
    peptide_length = []
    for i in range(len(protein_sequence)):
        for j in range(i + min_length, min(i + max_length + 1, len(protein_sequence) + 1)):
            peptides.append(protein_sequence[i:j])
            peptide_length += [len(peptides)]
            start_positions.append(i)
    return peptides, start_positions, peptide_length


def update_chart(x_data, y_data):
    ax.clear()
    ax.plot(x_data, y_data, color='green')
    ax.fill_between(x_data, y_data, color='green', alpha=0.3)

    ax.set_xlabel('Percentage of Protein Area Covered', labelpad=15)
    ax.set_ylabel('Percentage of HLA Coverage', labelpad=15)

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    # Adjust the figure to increase the white space around the plot
    fig.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.2)

    # Increase the font size of the labels to make them more readable
    ax.xaxis.label.set_size(12)
    ax.yaxis.label.set_size(12)

    canvas.draw()

    # Force the canvas to update, ensuring labels are drawn on top
    ax.figure.canvas.draw_idle()



def run_script():
    log_text.delete(1.0, tk.END)
    log_message("Script started.")

    protein_sequence = protein_entry.get().upper()
    min_length = min_length_entry.get()
    max_length = max_length_entry.get()
    hla_file = hla_file_entry.get()
    kd_threshold = int(kd_threshold_entry.get())
    hla_coverage_threshold = float(hla_coverage_threshold_entry.get())

    # Validate inputs
    if not protein_sequence or not min_length.isdigit() or not max_length.isdigit() or not hla_file or not kd_threshold:
        messagebox.showerror("Input Error", "Please ensure all fields are filled correctly.")
        return

    min_length = int(min_length)
    max_length = int(max_length)

    accepted_aas = set("ACDEFGHIKLMNPQRSTVWY")
    if set(protein_sequence) - accepted_aas:
        messagebox.showerror("Input Error", "Protein sequence contains invalid amino acids.")
        return

    log_message("Generating peptides...")
    peptides, start_positions, peptide_length = generate_peptides(protein_sequence, min_length, max_length)

    log_message("Predicting binding affinities...")
    predictor = Class1AffinityPredictor.load()
    HLA_type_I = pd.read_csv(hla_file)
    HLA_type_I_list = list(HLA_type_I.columns)
    peptide_list = list(peptides)

    kd_predictions = []
    hla_hits_for_peptide = []
    total_alleles = len(HLA_type_I_list)
    covered_alleles = set()
    progress_bar['value'] = 0
    coverage_bar['value'] = 0
    progress_step = 100 / (len(HLA_type_I_list) * len(peptide_list))

    x_data = []
    y_data = []

    for pep in peptide_list:
        hits = []
        for mhc_allele in HLA_type_I_list:
            try:
                predicted_affinity = predictor.predict([pep], [mhc_allele])
                kd_predictions += list(predicted_affinity)
                if predicted_affinity[0] <= kd_threshold:
                    hits.append(mhc_allele)
                    covered_alleles.add(mhc_allele)
            except ValueError:
                kd_predictions += ["NA"]

            progress_bar['value'] += progress_step
            progress_percent.set(f"{int(progress_bar['value'])}%")
            current_coverage = (len(covered_alleles) / total_alleles) * 100
            coverage_bar['value'] = current_coverage
            coverage_percent.set(f"{current_coverage:.2f}%")

            protein_area_covered = (start_positions[peptide_list.index(pep)] + len(pep)) / len(protein_sequence) * 100
            x_data.append(protein_area_covered)
            y_data.append(current_coverage)
            update_chart(x_data, y_data)

            root.update_idletasks()

            # Stop processing if coverage threshold is achieved
            if current_coverage >= hla_coverage_threshold:
                log_message(f"{hla_coverage_threshold}% HLA coverage threshold achieved. Stopping early.")
                messagebox.showinfo("Complete", f"{hla_coverage_threshold}% HLA coverage achieved. Process completed.")
                return

        hla_hits_for_peptide.append(", ".join(hits) if hits else "None")

    log_message("Creating DataFrame and processing results...")
    df = pd.DataFrame(kd_predictions, index=pd.MultiIndex.from_product([peptide_list, HLA_type_I_list]),
                      columns=["Binding Affinity"])
    df.index.names = ["Peptide", "MHC Allele"]
    df = df.reset_index()
    df = df.pivot(index="Peptide", columns="MHC Allele", values="Binding Affinity")
    df["peptide"] = peptides
    df["start"] = start_positions
    df["HLA Hits"] = hla_hits_for_peptide
    df1 = df.reindex(columns=["peptide", "start", "HLA Hits"] + list(df.columns[:-3]))

    # Save the results in real-time
    save_realtime_results(df1)

    df1.to_csv("binding_affinities.csv", index=False)

    df1 = pd.read_csv("binding_affinities.csv")
    grouped = df1.groupby(['start'])
    filtered_df = pd.DataFrame(columns=df1.columns)

    for group_name, group_data in grouped:
        min_values = group_data.drop(['start', 'peptide', 'HLA Hits'], axis=1).min()
        filtered_group = pd.DataFrame(columns=df1.columns)
        filtered_group['start'] = group_data['start']
        filtered_group['peptide'] = group_data['peptide']
        filtered_group['HLA Hits'] = group_data['HLA Hits']
        for col in min_values.index:
            filtered_group[col] = group_data[col].apply(lambda x: x if x == min_values[col] else 1000)
        filtered_df = pd.concat([filtered_df, filtered_group], ignore_index=True)
    df1 = filtered_df

    Kd_treshold_nM = 500
    exclude_cols = ["start", "peptide", "HLA Hits"]
    numeric_cols = df1.columns.drop(exclude_cols)
    df1[numeric_cols] = df1[numeric_cols].apply(pd.to_numeric, errors="coerce")
    df1[numeric_cols] = np.where(~np.isnan(df1[numeric_cols]), np.where(df1[numeric_cols] <= Kd_treshold_nM, 1, 0),
                                 np.nan)

    df1.to_csv("filterd_for_kd.csv", index=False)

    relative_HLA_coverage = 0
    exclude_cols = ["start", "peptide", "HLA Hits"]
    covered_hlas = []
    for col in df1.columns[3:]:
        if col not in exclude_cols:
            sum_value = df1[col].sum()
            if sum_value > 0:
                sum_value = 1
                covered_hlas.append(col)
            else:
                sum_value = 0
            relative_HLA_coverage += sum_value

    df1.to_csv('filtered_e7_netmhcpan_2915_alleles_8_to_11_mers.csv', index=False)
    df1['sum_weight'] = df1.sum(axis=1, skipna=True, numeric_only=True)
    df1.to_csv('E7_clement_result_table.csv', index=False)
    total1 = df1.shape[1] - 1
    df1 = df1.sort_values(by='sum_weight', axis=0, ascending=False, inplace=False, kind='quicksort', na_position='last',
                          ignore_index=False, key=None)
    maximum1 = max(df1['sum_weight'])
    lead_peptide = df1.iloc[0, 0]
    Lead_p = len(df1.iloc[0, 0])

    leadcoverage = (maximum1 / relative_HLA_coverage) * 100
    row_number = df1[df1['sum_weight'] == maximum1].index
    df1['sum_weight'] = df1.sum(axis=1, skipna=True, numeric_only=True)

    totalcoverage = leadcoverage
    total = total1
    df1.to_csv('jojojo.csv', index=False)

    col_name = 'start'
    col_to_move = df1.pop(col_name)
    df1.insert(1, col_name, col_to_move)

    max_num = maximum1

    totalscore = [total1]
    lead_pep_score = [max_num]
    lead_pep = [lead_peptide]
    Lead_peptide_lenght = [int(Lead_p)]
    rowindex = [row_number]
    accuulativcoverage = [leadcoverage]
    start = [df1.iloc[0, 1]]
    start_data = df1['start']
    df1['start'] = df1['start'].astype(str)
    while total > 0:
        df = df1.iloc[0]
        s = df.index
        g = []
        for i in range(len(s)):
            h = df.iloc[i]
            if type(h) != str:
                if h == 0:
                    g = [s[i]] + g
        ali = ['peptide', 'start', 'HLA Hits'] + g
        df1 = df1[ali]
        df1 = df1[1:(df1.shape[0])]
        df1['sum_weight'] = df1.sum(axis=1, skipna=True, numeric_only=True)
        total = df1['sum_weight'].sum()
        df1 = df1.sort_values(by='sum_weight', axis=0, ascending=False, inplace=False, kind='quicksort',
                              na_position='last', ignore_index=False, key=None)
        maximum = max(df1['sum_weight'])
        max_num += maximum
        lead_peptide = df1.iloc[0, 0]

        Lead_peptide_lenght += [len(df1.iloc[0, 0])]
        leadcoverage = [(maximum / total1) * 100]
        row_number = df1[df1['sum_weight'] == maximum].index
        absolut_acumulativelcoverage = (max_num / total1) * 100
        acumulativelcoverage = (max_num / relative_HLA_coverage) * 100

        accuulativcoverage += [acumulativelcoverage]
        lead_pep_score += [df1.iloc[0]['sum_weight']]
        lead_pep += [lead_peptide]
        rowindex += [(row_number[0]) + 1]
        start += df1.iloc[0, 1]

    data = {
        "start position in protein": start,
        "Peptide length": Lead_peptide_lenght,
        'peptide sequence': lead_pep,
        '# HLA allele hits': lead_pep_score,
        'accumulative coverage (%)': accuulativcoverage
    }
    global df3
    df3 = pd.DataFrame(data)

    update_treeview(df3)

    # Display HLA hits and coverage percentage
    hla_hits_text.set(", ".join(covered_hlas))
    hla_coverage_percentage = (relative_HLA_coverage / len(HLA_type_I_list)) * 100
    hla_coverage_text.set(f"{hla_coverage_percentage:.2f}%")

    if hla_coverage_percentage < 100:
        messagebox.showwarning("Incomplete Coverage",
                               "The generated peptides could not cover the entire list of HLA alleles. "
                               "Consider increasing the length of your input protein sequence, widening range of peptide length and increasing the Kd threshold for MHC binding.")

    log_message("Script finished successfully.")
    messagebox.showinfo("Success", "The script has been run successfully.")

def update_treeview(df):
    tree.delete(*tree.get_children())
    tree["column"] = list(df.columns)
    tree["show"] = "headings"

    for column in df.columns:
        tree.heading(column, text=column)

    df_rows = df.to_numpy().tolist()
    for row in df_rows:
        tree.insert("", "end", values=row)

def browse_file():
    filename = filedialog.askopenfilename()
    hla_file_entry.delete(0, tk.END)
    hla_file_entry.insert(0, filename)

def save_file():
    filetypes = [('CSV Files', '*.csv'), ('Excel Files', '*.xlsx')]
    filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=filetypes)
    if not filepath:
        return
    if filepath.endswith('.csv'):
        df3.to_csv(filepath, index=False)
    elif filepath.endswith('.xlsx'):
        df3.to_excel(filepath, index=False)
    messagebox.showinfo("Save Successful", f"File saved as {filepath}")

def toggle_log():
    if log_frame.winfo_ismapped():
        log_frame.grid_remove()
        toggle_log_button.config(text="Show Log")
    else:
        log_frame.grid(row=14, column=0, columnspan=2, sticky="nsew")
        toggle_log_button.config(text="Hide Log")

def start_long_running_task():
    task_thread = threading.Thread(target=run_script)
    task_thread.start()

# Setting up the main window
root = tk.Tk()
root.title("Minimal Epitiope for Maximum MHC Coverage 2.0")
root.geometry("1024x800")

# Creating the input fields
tk.Label(root, text="Protein Sequence:").grid(row=0, column=0, sticky="e")
protein_entry = tk.Entry(root, width=50)
protein_entry.grid(row=0, column=1, padx=10, pady=10)

tk.Label(root, text="Min Peptide Length:").grid(row=1, column=0, sticky="e")
min_length_entry = tk.Entry(root)
min_length_entry.grid(row=1, column=1, padx=10, pady=10)
min_length_entry.insert(0, "8")

tk.Label(root, text="Max Peptide Length:").grid(row=2, column=0, sticky="e")
max_length_entry = tk.Entry(root)
max_length_entry.grid(row=2, column=1, padx=10, pady=10)
max_length_entry.insert(0, "9")

tk.Label(root, text="HLA File:").grid(row=3, column=0, sticky="e")
hla_file_entry = tk.Entry(root, width=50)
hla_file_entry.grid(row=3, column=1, padx=10, pady=10)
tk.Button(root, text="Browse", command=browse_file).grid(row=3, column=2, padx=10, pady=10)

tk.Label(root, text="Kd Threshold (nM):").grid(row=4, column=0, sticky="e")
kd_threshold_entry = tk.Entry(root)
kd_threshold_entry.grid(row=4, column=1, padx=10, pady=10)
kd_threshold_entry.insert(0, "500")

tk.Label(root, text="HLA Coverage Threshold (%):").grid(row=5, column=0, sticky="e")
hla_coverage_threshold_entry = tk.Entry(root)
hla_coverage_threshold_entry.grid(row=5, column=1, padx=10, pady=10)
hla_coverage_threshold_entry.insert(0, "100")

# Progress bar with percentage label
tk.Label(root, text="Processing Progress:").grid(row=6, column=0, sticky="e")

progress_percent = tk.StringVar()
progress_percent.set("0%")
progress_bar = ttk.Progressbar(root, orient="horizontal", length=400, mode="determinate",
                               style="green.Horizontal.TProgressbar")
progress_bar.grid(row=6, column=1, padx=10, pady=10)
progress_label = tk.Label(root, textvariable=progress_percent)
progress_label.grid(row=6, column=2, padx=10, pady=10)

# Coverage bar with percentage label
tk.Label(root, text="HLA Allele Hit Coverage:").grid(row=7, column=0, sticky="e")

coverage_percent = tk.StringVar()
coverage_percent.set("0%")
coverage_bar = ttk.Progressbar(root, orient="horizontal", length=400, mode="determinate",
                               style="green.Horizontal.TProgressbar")
coverage_bar.grid(row=7, column=1, padx=10, pady=10)
coverage_label = tk.Label(root, textvariable=coverage_percent)
coverage_label.grid(row=7, column=2, padx=10, pady=10)

# Display HLA Hits and Coverage
tk.Label(root, text="HLA Alleles Hit:").grid(row=8, column=0, sticky="e")
hla_hits_text = tk.StringVar()
hla_hits_label = tk.Label(root, textvariable=hla_hits_text, wraplength=800, justify="left")
hla_hits_label.grid(row=8, column=1, padx=10, pady=10)

tk.Label(root, text="Percentage of HLA Coverage:").grid(row=9, column=0, sticky="e")
hla_coverage_text = tk.StringVar()
hla_coverage_label = tk.Label(root, textvariable=hla_coverage_text)
hla_coverage_label.grid(row=9, column=1, padx=10, pady=10)

# Treeview to display the results
tree = ttk.Treeview(root, height=10)
tree.grid(row=10, column=0, columnspan=3, padx=10, pady=10, sticky="nsew")

# Bottom row with buttons
run_button = tk.Button(root, text="Run", command=start_long_running_task)
run_button.grid(row=11, column=0, sticky="e", padx=10, pady=10)

save_button = tk.Button(root, text="Download", command=save_file)
save_button.grid(row=11, column=1, sticky="w", padx=10, pady=10)

abort_button = tk.Button(root, text="Abort", command=root.quit)
abort_button.grid(row=11, column=2, sticky="e", padx=10, pady=10)

toggle_log_button = tk.Button(root, text="Show Log", command=toggle_log)
toggle_log_button.grid(row=11, column=3, sticky="w", padx=10, pady=10)

log_frame = tk.Frame(root)
log_text = ScrolledText(log_frame, height=10, state='normal')
log_text.pack(fill=tk.BOTH, expand=True)

root.grid_rowconfigure(10, weight=1)
root.grid_columnconfigure(1, weight=1)

# Adding Logo and MHC-Peptide Image
logo_label = tk.Label(root, text="MEMHC2.0", font=("Helvetica", 24, "bold"))
logo_label.grid(row=0, column=3, rowspan=2, padx=10, pady=10, sticky="nsew")

# Live updating chart for HLA coverage vs Protein Area
fig, ax = plt.subplots(figsize=(5, 5))
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=2, column=3, rowspan=8, padx=10, pady=10, sticky="nsew")

# Start the Tkinter event loop
root.mainloop()
