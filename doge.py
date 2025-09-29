#!/usr/bin/env python3

# doge.py - A script to discover organelle genome variants.
# Copyright (C) 2025  Zhi-Jie Xu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__version__ = "1.0.1"

import sys
import os
import re
import argparse
import subprocess
import glob
from collections import defaultdict
from Bio import SeqIO
# from Bio.Align import MultipleSeqAlignment
# from Bio.Phylo.Applications import _Fasttree
import concurrent.futures
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                            QLineEdit, QPushButton, QCheckBox, QSpinBox,
                            QGroupBox, QTextEdit, QFileDialog, QMessageBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

class DOGE_UI(QWidget):
    def __init__(self, args=False):
        super().__init__()
        self.setWindowTitle("DOGE: Discover Organelle Genomes")
        self.setGeometry(300, 300, 600, 500)
        self.init_ui()

        if args:
            self.prefill_from_args(args)
        
        # 存储工具路径
        self.mauve_path = ""
        self.mafft_path = ""
        self.trimal_path = ""
    
    def prefill_from_args(self, args):
        if args.input:
            self.input_edit.setText(os.path.abspath(args.input))
        if args.output:
            self.output_edit.setText(os.path.abspath(args.output))
        
        
        self.mauve_edit.setText(os.path.abspath(args.mauve))
        if args.mafft:
            self.mafft_edit.setText(os.path.abspath(args.mafft))
        if args.trimal:
            self.trimal_edit.setText(os.path.abspath(args.trimal))
        
        self.filter_check.setChecked(args.filter == 'True')
        self.trimal_check.setChecked(args.trim == 'True')
        self.mafft_check.setChecked(args.refine_alignment == 'True')
        
        self.coverage_spin.setValue(args.lowest_coverage)

        if args.filter == 'True':
            self.coverage_spin.setEnabled(True)
        if args.refine_alignment == 'True':
            self.mafft_edit.setEnabled(True)
        if args.trim == 'True':
            self.trimal_edit.setEnabled(True)
        
        self.highlight_prefilled_fields(args)

    def highlight_prefilled_fields(self, args):
        highlighted_style = "QLineEdit { background-color: #eeffee; border: 1px solid #00ff00; }"
        
        fields_to_highlight = [ ]
        if args.input:
            fields_to_highlight.append(self.input_edit)
        if args.output:
            fields_to_highlight.append(self.output_edit)
        if args.mauve:
            fields_to_highlight.append(self.mauve_edit)
        
        if self.mafft_check.isChecked():
            fields_to_highlight.append(self.mafft_edit)
        if self.trimal_check.isChecked():
            fields_to_highlight.append(self.trimal_edit)
        if self.filter_check.isChecked():
            self.coverage_spin.setStyleSheet("QSpinBox { background-color: #eeffee  ; }")
        
        for field in fields_to_highlight:
            field.setStyleSheet(highlighted_style)
            field.setToolTip("This argument comes from preset.")

    def init_ui(self):
        main_layout = QVBoxLayout()
        
        folder_group = QGroupBox("Input/Output Directories")
        folder_layout = QVBoxLayout()
        
        input_layout = QHBoxLayout()
        input_layout.addWidget(QLabel("Input Directory:"))
        self.input_edit = QLineEdit()
        self.input_edit.setPlaceholderText("Select folder containing genome FASTA files")
        input_layout.addWidget(self.input_edit)
        input_btn = QPushButton("Browse...")
        input_btn.clicked.connect(self.select_input_dir)
        input_layout.addWidget(input_btn)
        folder_layout.addLayout(input_layout)
        
        output_layout = QHBoxLayout()
        output_layout.addWidget(QLabel("Output Directory:"))
        self.output_edit = QLineEdit()
        self.output_edit.setPlaceholderText("Select folder for analysis results")
        output_layout.addWidget(self.output_edit)
        output_btn = QPushButton("Browse...")
        output_btn.clicked.connect(self.select_output_dir)
        output_layout.addWidget(output_btn)
        folder_layout.addLayout(output_layout)
        
        folder_group.setLayout(folder_layout)
        main_layout.addWidget(folder_group)
        
        param_group = QGroupBox("Analysis Parameters")
        param_layout = QVBoxLayout()
        
        mauve_layout = QHBoxLayout()
        mauve_layout.addWidget(QLabel("ProgressiveMauve Path:"))
        self.mauve_edit = QLineEdit()
        self.mauve_edit.setPlaceholderText("Path to progressiveMauve executable")
        mauve_layout.addWidget(self.mauve_edit)
        mauve_btn = QPushButton("Locate...")
        mauve_btn.clicked.connect(self.locate_mauve)
        mauve_layout.addWidget(mauve_btn)
        param_layout.addLayout(mauve_layout)
        
        self.mafft_check = QCheckBox("Refine alignment with MAFFT (May be Harmful!)")
        self.mafft_check.setChecked(False)
        self.mafft_check.stateChanged.connect(self.toggle_mafft)
        
        mafft_layout = QHBoxLayout()
        mafft_layout.addWidget(QLabel("MAFFT Path:"))
        self.mafft_edit = QLineEdit()
        self.mafft_edit.setPlaceholderText("Path to MAFFT executable")
        self.mafft_edit.setEnabled(False)
        mafft_layout.addWidget(self.mafft_edit)
        self.mafft_btn = QPushButton("Locate...")
        self.mafft_btn.clicked.connect(self.locate_mafft)
        self.mafft_btn.setEnabled(False)
        mafft_layout.addWidget(self.mafft_btn)
        
        param_layout.addWidget(self.mafft_check)
        param_layout.addLayout(mafft_layout)
        
        self.filter_check = QCheckBox("Filter by sequence coverage")
        self.filter_check.setChecked(True)
        self.filter_check.stateChanged.connect(self.toggle_filter)
        
        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("Minimum coverage:"))
        self.coverage_spin = QSpinBox()
        self.coverage_spin.setRange(2, 100)
        self.coverage_spin.setValue(2)
        filter_layout.addWidget(self.coverage_spin)
        filter_layout.addWidget(QLabel("sequences per LCB"))
        filter_layout.addStretch()
        
        param_layout.addWidget(self.filter_check)
        param_layout.addLayout(filter_layout)
        
        self.trimal_check = QCheckBox("Trim alignments with TrimAl")
        self.trimal_check.setChecked(True)
        self.trimal_check.stateChanged.connect(self.toggle_trimal)
        
        trimal_layout = QHBoxLayout()
        trimal_layout.addWidget(QLabel("TrimAl Path:"))
        self.trimal_edit = QLineEdit()
        self.trimal_edit.setPlaceholderText("Path to trimAl executable")
        trimal_layout.addWidget(self.trimal_edit)
        self.trimal_btn = QPushButton("Locate...")
        self.trimal_btn.clicked.connect(self.locate_trimal)
        trimal_layout.addWidget(self.trimal_btn)
        
        param_layout.addWidget(self.trimal_check)
        param_layout.addLayout(trimal_layout)
        
        param_group.setLayout(param_layout)
        main_layout.addWidget(param_group)
        
        run_layout = QHBoxLayout()
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold;")
        self.run_btn.clicked.connect(self.run_analysis)
        self.run_btn.setFixedHeight(40)
        run_layout.addWidget(self.run_btn)
        
        self.log_btn = QPushButton("Show Logs")
        self.log_btn.setEnabled(False)
        self.log_btn.setFixedHeight(40)
        run_layout.addWidget(self.log_btn)
        
        main_layout.addLayout(run_layout)
        
        # Logger
        log_group = QGroupBox("Analysis Log")
        log_layout = QVBoxLayout()
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        log_layout.addWidget(self.log_text)
        log_group.setLayout(log_layout)
        main_layout.addWidget(log_group)
        
        self.setLayout(main_layout)
    
    def select_input_dir(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Input Directory")
        if directory:
            self.input_edit.setText(directory)
    
    def select_output_dir(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if directory:
            self.output_edit.setText(directory)
    
    def locate_mauve(self):
        path, _ = QFileDialog.getOpenFileName(self, "Locate ProgressiveMauve")
        if path:
            self.mauve_edit.setText(path)
            self.mauve_path = path
    
    def locate_mafft(self):
        path, _ = QFileDialog.getOpenFileName(self, "Locate MAFFT")
        if path:
            self.mafft_edit.setText(path)
            self.mafft_path = path
    
    def locate_trimal(self):
        path, _ = QFileDialog.getOpenFileName(self, "Locate TrimAl")
        if path:
            self.trimal_edit.setText(path)
            self.trimal_path = path
    
    def toggle_mafft(self, state):
        self.mafft_edit.setEnabled(state)
        self.mafft_btn.setEnabled(state)
        # for i in range(self.mafft_edit.parent().layout().count()):
        #     widget = self.mafft_edit.parent().layout().itemAt(i).widget()
        #     if widget and widget != self.mafft_edit:
        #         widget.setEnabled(state)
    
    def toggle_filter(self, state):
        self.coverage_spin.setEnabled(state)
    
    def toggle_trimal(self, state):
        self.trimal_edit.setEnabled(state)
        self.trimal_btn.setEnabled(state)
        # for i in range(self.trimal_edit.parent().layout().count()):
        #     widget = self.trimal_edit.parent().layout().itemAt(i).widget()
        #     if widget and widget != self.trimal_edit:
        #         widget.setEnabled(state)
    
    def validate_inputs(self):
        errors = []
        
        if not self.input_edit.text():
            errors.append("Input directory is required")
        elif not os.path.isdir(self.input_edit.text()):
            errors.append("Input directory does not exist")
            
        if not self.output_edit.text():
            errors.append("Output directory is required")
            
        if not self.mauve_edit.text():
            errors.append("ProgressiveMauve path is required")
        elif not os.path.isfile(self.mauve_edit.text()):
            errors.append("ProgressiveMauve executable not found")
            
        if self.mafft_check.isChecked() and not self.mafft_edit.text():
            errors.append("MAFFT path is required when refinement is enabled")
        elif self.mafft_check.isChecked() and not os.path.isfile(self.mafft_edit.text()):
            errors.append("MAFFT executable not found")
        if self.trimal_check.isChecked() and not self.trimal_edit.text():
            errors.append("TrimAl path is required when trimming is enabled")
        elif self.trimal_check.isChecked() and not os.path.isfile(self.trimal_edit.text()):
            errors.append("TrimAl executable not found")
            
        if self.input_edit.text():
            fasta_files = [f for f in os.listdir(self.input_edit.text()) 
                          if f.endswith(('.fasta', '.fa', '.fna'))]
            if not fasta_files:
                errors.append("No FASTA files found in input directory")
        
        return errors
    
    def run_analysis(self):
        errors = self.validate_inputs()
        if errors:
            QMessageBox.critical(self, "Input Errors", "\n".join(errors))
            return
        
        self.run_btn.setEnabled(False)
        self.run_btn.setText("Running...")
        
        params = {
            'input_dir': self.input_edit.text(),
            'output_dir': self.output_edit.text(),
            'mauve_path': self.mauve_edit.text(),
            'mafft_path': self.mafft_edit.text(),
            'trimal_path': self.trimal_edit.text(),
            'lowest_coverage': self.coverage_spin.value(),
            'refine_alignment': self.mafft_check.isChecked(),
            'filter': self.filter_check.isChecked(),
            'trim': self.trimal_check.isChecked()
        }
        
        self.log_text.clear()
        self.log_text.append("=== Starting DOGE Analysis ===")
        self.log_text.append(f"Input directory: {params['input_dir']}")
        self.log_text.append(f"Output directory: {params['output_dir']}")
        self.log_text.append("\nParameters:")
        self.log_text.append(f"  ProgressiveMauve path: {params['mauve_path']}")
        self.log_text.append(f"  Refine with MAFFT: {params['refine_alignment']}")
        self.log_text.append(f"  Filter coverage: {params['filter']} (min: {params['lowest_coverage']})")
        self.log_text.append(f"  Trim with TrimAl: {params['trim']}")
        self.log_text.append("\nRunning ProgressiveMauve alignment...")
        
        self.analysis_thread = AnalysisWorker(params)
        
        self.analysis_thread.log_message.connect(self.append_log)
        self.analysis_thread.progress_update.connect(self.update_progress)
        self.analysis_thread.finished.connect(self.analysis_completed)
        self.analysis_thread.error_occurred.connect(self.analysis_failed)

        self.analysis_thread.start()

    def append_log(self, message):
        self.log_text.append(message)
    
    def update_progress(self, current, total):
        percent = int((current / total) * 100)
        self.log_text.append(f"Progress: {current}/{total} ({percent}%)")
    
    def analysis_completed(self, results):
        self.log_text.append("\n✔️ Analysis Completed!")
        self.log_text.append(f"Final Output: {results['nexus']}")
        self.log_text.append(f"Proceeded {results['lcb_count']} LCBs")
        
        self.run_btn.setEnabled(True)
        self.log_btn.setEnabled(True)
        QMessageBox.information(self, "Completed", "Analysis Completed Successfully!")
        
    def analysis_failed(self, error_msg):
        self.log_text.append(f"\n❌ Error: {error_msg}")
        self.set_ui_enabled(True)
        QMessageBox.critical(self, "Analysis Failed", error_msg)


class AnalysisWorker(QThread):
    log_message = pyqtSignal(str)
    progress_update = pyqtSignal(int, int)
    finished = pyqtSignal(dict)
    error_occurred = pyqtSignal(str)

    def __init__(self, params, parent=None):
        super().__init__(parent)
        self.params = params
        self.is_running = True

    def run(self):
        try:
            total_steps = 8
            current_step = 1
            
            self.log_message.emit("===== Restart =====")
            self.progress_update.emit(current_step, total_steps)
            current_step += 1
            
            self.log_message.emit("\n=== Running ProgressiveMauve ===")
            mauve_output = run_mauve(
                self.params['mauve_path'],
                self.params['input_dir'],
                os.path.join(self.params['output_dir'], "mauve_alignment.xmfa")
            )
            self.progress_update.emit(current_step, total_steps)
            current_step += 1
            
            self.log_message.emit("\n=== Splitting xmfa into LCBs ===")
            lcb_dir = os.path.join(self.params['output_dir'], "lcbs")
            lcb_files = split_xmfa(mauve_output, lcb_dir)
            self.progress_update.emit(current_step, total_steps)
            current_step += 1
            
            if self.params['refine_alignment']:
                self.log_message.emit("\n=== Refining Alignments with MAFFT ===")
                aligned_dir = os.path.join(self.params['output_dir'], "aligned_lcbs")
                lcb_files = align_all_lcbs(
                    self.params['mafft_path'], 
                    lcb_files, 
                    aligned_dir
                )
            self.progress_update.emit(current_step, total_steps)
            current_step += 1
            
            if self.params['filter']:
                self.log_message.emit("\n=== Filtering out Low-coverage LCBs ===")
                lcb_files = filter_lcbs(
                    self.params['lowest_coverage'],
                    lcb_files
                )
            self.progress_update.emit(current_step, total_steps)
            current_step += 1
            
            if self.params['trim']:
                self.log_message.emit("\n=== Trimming alignments with TrimAl ===")
                trimmed_dir = os.path.join(self.params['output_dir'], "trimmed_lcbs")
                lcb_files = trim_lcbs(
                    self.params['trimal_path'],
                    lcb_files,
                    trimmed_dir 
                )
            self.progress_update.emit(current_step, total_steps)
            current_step += 1
            
            self.log_message.emit("\n=== Constructing Supermatrix ===")
            supermatrix_file = os.path.join(self.params['output_dir'], "supermatrix.fasta")
            supermatrix, partitions = create_supermatrix(lcb_files, supermatrix_file)
            self.progress_update.emit(current_step, total_steps)
            current_step += 1
            
            self.log_message.emit("\n=== Calling SNPs; Generating Partition File ===")
            nexus_file = call_snps(
                supermatrix, 
                partitions, 
                os.path.join(self.params['output_dir'], "snps")
            )
            
            self.log_message.emit("\n=== Analysis Completed! ===")
            self.finished.emit({
                "supermatrix": supermatrix_file,
                "nexus": nexus_file,
                "lcb_count": len(lcb_files)
            })
            
        except Exception as e:
            self.error_occurred.emit(f"Analysis Failed: {str(e)}")
        finally:
            self.is_running = False

    def stop(self):
        self.is_running = False
        self.terminate()

def run_mauve(mauve_path, input_dir, output_file):
    """Run progressiveMauve on all fasta files in input directory"""
    input_files = glob.glob(os.path.join(input_dir, "*.fasta"))
    input_files += glob.glob(os.path.join(input_dir, "*.fa"))
    input_files += glob.glob(os.path.join(input_dir, "*.fna"))
    
    if not input_files:
        print("Error: No input fasta files found in directory")
        sys.exit(1)
    
    cmd = [mauve_path, f"--output={output_file}"] + input_files
    print(f"Running Mauve: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Mauve failed with error:\n{result.stderr}")
        sys.exit(1)
    
    print("Mauve alignment completed successfully")
    return output_file

def split_xmfa(xmfa_file, output_dir):
    """Split .xmfa file into separate LCB FASTA files"""
    if not os.path.isfile(xmfa_file):
        print(f"Error: File '{xmfa_file}' not found!")
        sys.exit(1)
    
    os.makedirs(output_dir, exist_ok=True)
    
    with open(xmfa_file, 'r') as f:
        lines = f.readlines()
    
    lcb_blocks = []
    current_block = []
    in_block = False
    block_count = 0
    
    separator_pattern = re.compile(r'^[=#]+$')
    
    for line in lines:
        stripped_line = line.strip()
        
        if separator_pattern.match(stripped_line):
            if in_block and current_block:
                lcb_blocks.append(current_block)
                current_block = []
                in_block = False
            continue
        
        if line.startswith('>'):
            in_block = True
            clean_header = re.sub(r'>\s*(\d+)\s*:\s*(\d+)\s*-\s*(\d+)\s*', r'>\1:\2-\3', line)
            current_block.append(clean_header)
        elif in_block:
            current_block.append(line)
    
    if current_block:
        lcb_blocks.append(current_block)
    
    output_files = []
    for i, block in enumerate(lcb_blocks):
        output_file = os.path.join(output_dir, f"LCB_{i+1}.fasta")
        with open(output_file, 'w') as f_out:
            f_out.writelines(block)
        output_files.append(output_file)
        print(f"Created {output_file} with {len(block)//2} sequences")
    
    print(f"\nSuccessfully split {len(lcb_blocks)} LCBs from {xmfa_file}")
    return output_files

def run_mafft(mafft_path, input_file, output_file):
    """Run MAFFT on a single LCB file"""
    cmd = [mafft_path, '--auto', '--inputorder', input_file]
    print(f"Running MAFFT on {os.path.basename(input_file)}")
    
    try:
        with open(output_file, 'w') as f_out:
            result = subprocess.run(cmd, stdout=f_out, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            print(f"MAFFT failed for {input_file}:\n{result.stderr}")
            return None
        
        print(f"MAFFT completed for {os.path.basename(input_file)}")
        return output_file
    except Exception as e:
        print(f"Error running MAFFT: {str(e)}")
        return None

def align_all_lcbs(mafft_path, lcb_files, output_dir):
    """Align all LCB files using MAFFT in parallel"""
    os.makedirs(output_dir, exist_ok=True)
    aligned_files = []
    
    output_paths = [os.path.join(output_dir, f"aligned_{os.path.basename(f)}") for f in lcb_files]
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for input_file, output_file in zip(lcb_files, output_paths):
            futures.append(executor.submit(run_mafft, mafft_path, input_file, output_file))
        
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                aligned_files.append(result)
    
    return aligned_files

def filter_lcbs(lowest_coverage, lcb_files):
    """Filter out low-coverage LCBs"""
    valid_lcbs = []
    
    for lcb_file in lcb_files:

        seq_count = 0
        total_length = 0

        try:
            for record in SeqIO.parse(lcb_file, "fasta"):
                seq_count += 1
        except:
            print(f"Error reading {lcb_file}, skipping")
            continue
        
        if seq_count >= lowest_coverage:
            valid_lcbs.append(lcb_file)
        else:
            print(f"Discarded LCB: {os.path.basename(lcb_file)} "
                    f"(Taxa: {seq_count})")
    
    print(f"Retained {len(valid_lcbs)}/{len(lcb_files)} LCBs after filtering")
    return valid_lcbs

def trim_lcbs(trimal_path, aligned_files, output_dir):
    """Apply TrimAl to each aligned LCB in parallel"""
    os.makedirs(output_dir, exist_ok=True)
    trimmed_files = []
    
    def run_trimal(aligned_file):
        output_file = os.path.join(output_dir, f"trimmed_{os.path.basename(aligned_file)}")
        cmd = [trimal_path, "-in", aligned_file, "-out", output_file, "-automated1"]
        
        print(f"Running TrimAl on {os.path.basename(aligned_file)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"TrimAl failed for {aligned_file}:\n{result.stderr}")
            return None
        else:
            print(f"TrimAl completed for {os.path.basename(aligned_file)}")
            return output_file
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(run_trimal, file) for file in aligned_files]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                trimmed_files.append(result)
    
    return trimmed_files

def create_supermatrix(aligned_files, output_file):
    """Combine aligned LCBs into a supermatrix with proper gap handling"""
    all_sample_ids = set()
    lcb_lengths = []  
    
    for aligned_file in aligned_files:
        for record in SeqIO.parse(aligned_file, "fasta"):
            sample_id = record.id.split(':')[0]
            all_sample_ids.add(sample_id)
    
    sequences = {sample_id: [] for sample_id in all_sample_ids}
    partitions = []
    current_position = 1
    
    for i, aligned_file in enumerate(aligned_files):
        lcb_id = i + 1
        
        lcb_sequences = {}
        lcb_length = 0
        for record in SeqIO.parse(aligned_file, "fasta"):
            sample_id = record.id.split(':')[0]
            seq_str = str(record.seq)
            lcb_sequences[sample_id] = seq_str
            
            if lcb_length == 0:
                lcb_length = len(seq_str)
        
        partitions.append((lcb_id, current_position, current_position + lcb_length - 1))
        current_position += lcb_length
        
        for sample_id in all_sample_ids:
            if sample_id in lcb_sequences:

                sequences[sample_id].append(lcb_sequences[sample_id])
            else:

                sequences[sample_id].append('?' * lcb_length)
    
    concatenated_sequences = {}
    for sample_id, seq_parts in sequences.items():
        concatenated_sequences[sample_id] = ''.join(seq_parts)
    
    with open(output_file, 'w') as f_out:
        for sample_id, seq in concatenated_sequences.items():
            f_out.write(f">{sample_id}\n{seq}\n")
    
    total_length = len(next(iter(concatenated_sequences.values())))
    print(f"Created supermatrix with {len(concatenated_sequences)} samples and {total_length} positions")
    return output_file, partitions

def extract_core_snps(supermatrix_file, partitions, output_prefix):
    """Extract core SNPs (present in all samples without gaps/missing) and write NEXUS.

    Core SNP definition here: positions where every sample has a canonical base
    (A/C/G/T, case-insensitive) and there is variation among samples.
    """

    valid_bases = set(['A', 'C', 'G', 'T'])

    sequences_dict = {}
    sample_ids = []
    for record in SeqIO.parse(supermatrix_file, "fasta"):
        seq = str(record.seq).upper()
        sequences_dict[record.id] = seq
        sample_ids.append(record.id)

    if not sequences_dict:
        return None

    seq_length = len(next(iter(sequences_dict.values())))
    if any(len(seq) != seq_length for seq in sequences_dict.values()):
        print("Error: Sequences in supermatrix have different lengths!")
        return None

    # Global position -> partition id
    position_to_partition = {}
    for part_id, start, end in partitions:
        for pos in range(start - 1, end):  # convert to 0-based
            position_to_partition[pos] = part_id

    # Identify core variable positions
    core_variable_positions = []  # 0-based indices
    for pos in range(seq_length):
        bases = []
        is_core = True
        for seq in sequences_dict.values():
            base = seq[pos]
            if base not in valid_bases:
                is_core = False
                break
            bases.append(base)
        if not is_core:
            continue
        if len(set(bases)) > 1:
            core_variable_positions.append(pos)

    # Build SNP matrix
    snp_matrix = []
    for sample_id in sample_ids:
        seq = sequences_dict[sample_id]
        snp_seq = [seq[pos] for pos in core_variable_positions]
        snp_matrix.append(''.join(snp_seq))

    nexus_file = f"{output_prefix}\\core_snps.nex"
    with open(nexus_file, 'w') as f:
        f.write("#NEXUS\n\n")
        f.write("BEGIN DATA;\n")
        f.write(f"DIMENSIONS NTAX={len(sample_ids)} NCHAR={len(core_variable_positions)};\n")
        f.write("FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        f.write("MATRIX\n")

        for sample_id, seq in zip(sample_ids, snp_matrix):
            f.write(f"{sample_id:<20} {seq}\n")

        f.write(";\nEND;\n\n")

        # Partition scheme for core SNPs (map original positions to SNP columns)
        f.write("BEGIN SETS;\n")
        partition_positions = defaultdict(list)
        for i, pos in enumerate(core_variable_positions):
            part_id = position_to_partition.get(pos)
            if part_id is not None:
                partition_positions[part_id].append(i + 1)  # 1-based in NEXUS

        for part_id, positions in partition_positions.items():
            if positions:
                start_i = min(positions)
                end_i = max(positions)
                f.write(f"CHARSET LCB_{part_id} = {start_i}-{end_i};\n")
        f.write("END;\n")

    # Also write an IQ-TREE partition file limited to core SNP columns
    iqpartition_file = f"{output_prefix}\\core_snps_iqpartition.nex"
    with open(iqpartition_file, 'w') as fp:
        fp.write('#NEXUS\n')
        fp.write('begin sets;\n')
        for part_id, positions in partition_positions.items():
            if positions:
                start_i = min(positions)
                end_i = max(positions)
                fp.write(f"CHARSET LCB_{part_id} = {start_i}-{end_i};\n")
        fp.write(';\nend;\n')

    print(f"Created core SNP NEXUS file: {nexus_file}")
    return nexus_file

def call_snps(supermatrix_file, partitions, output_prefix):
    """Call SNPs from supermatrix and create partitioned NEXUS file"""

    sequences_dict = {}
    sample_ids = []
    for record in SeqIO.parse(supermatrix_file, "fasta"):
        sequence = str(record.seq)
        if set(sequence) == {'-'}: # process missing sequence in a single LCB
            sequence_dict[record.id] = '?'*len(sequence)
        else:
            sequences_dict[record.id] = str(record.seq)
        sample_ids.append(record.id)
    
    seq_length = len(next(iter(sequences_dict.values())))
    if any(len(seq) != seq_length for seq in sequences_dict.values()):
        print("Error: Sequences in supermatrix have different lengths!")
        return None
    
    # partitioning
    position_to_partition = {}
    for part_id, start, end in partitions:
        for pos in range(start, end + 1):
            position_to_partition[pos] = part_id
    
    # SNP calling
    variable_positions = []
    for pos in range(seq_length):
        bases = set()
        for seq in sequences_dict.values():
            base = seq[pos]
            if base not in ['-', '?']:
                bases.add(base)
        
        if len(bases) > 1:
            variable_positions.append(pos)
    
    snp_matrix = []
    for sample_id in sample_ids:
        snp_seq = []
        for pos in variable_positions:
            base = sequences_dict[sample_id][pos]
            snp_seq.append(base)
        snp_matrix.append(''.join(snp_seq))
    
    nexus_file = f"{output_prefix}.nex"
    with open(nexus_file, 'w') as f:
        f.write("#NEXUS\n\n")
        f.write("BEGIN DATA;\n")
        f.write(f"DIMENSIONS NTAX={len(sample_ids)} NCHAR={len(variable_positions)};\n")
        f.write("FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        f.write("MATRIX\n")
        
        for sample_id, seq in zip(sample_ids, snp_matrix):
            f.write(f"{sample_id:<20} {seq}\n")
        
        f.write(";\nEND;\n\n")
        
        # partition scheme
        f.write("BEGIN SETS;\n")
        
        partition_positions = defaultdict(list)
        for i, pos in enumerate(variable_positions):
            part_id = position_to_partition.get(pos)
            if part_id:
                partition_positions[part_id].append(str(i+1))
        
        for part_id, positions in partition_positions.items():
            if positions:
                f.write(f"CHARSET LCB_{part_id} = {min(positions)}-{max(positions)};\n")
        
        f.write("END;\n")

        # create another partition NEXUS file for modelfinder / other programs
        with open(f'{output_prefix}_iqpartition.nex', 'w') as fp:
            fp.write('#NEXUS\n')
            fp.write('begin sets;\n')
            for part_id, positions in partition_positions.items():
                if positions:
                    f.write(f"CHARSET LCB_{part_id} = {min(positions)}-{max(positions)};\n")
            fp.write(';\nend;\n')
        

    
    print(f"Created partitioned NEXUS file: {nexus_file}")
    return nexus_file

# def run_fasttree(fasttree_path, nexus_file, output_tree):
#     """Run FastTree on the SNP matrix"""
#     cmd = [fasttree_path, "-nt", "-gtr", nexus_file]
    
#     try:
#         with open(output_tree, 'w') as f_out:
#             result = subprocess.run(cmd, stdout=f_out, stderr=subprocess.PIPE, text=True)
        
#         if result.returncode != 0:
#             print(f"FastTree failed:\n{result.stderr}")
#             return None
        
#         print(f"Phylogenetic tree saved to {output_tree}")
#         return output_tree
#     except Exception as e:
#         print(f"Error running FastTree: {str(e)}")
#         return None

def launch_gui_with_args(args):
    app = QApplication(sys.argv)
    ui = DOGE_UI(args)
    ui.show()
    sys.exit(app.exec_())

def main():
    parser = argparse.ArgumentParser(description='DOGE: Discover Organelle Genome Variants')
    parser.add_argument('-i', '--input', help='Input directory with genome FASTA files [required]')
    parser.add_argument('-o', '--output', help='Output directory for results [required]')
    parser.add_argument('--mauve', default='./mauve/progressiveMauve.exe', help='Path to Mauve executable')
    parser.add_argument('--mafft', default='./mafft/mafft.bat', help='Path to MAFFT executable')
    parser.add_argument('--trimal', default='./trimAl/bin/trimal.exe', help='Path to TrimAl executable')
    parser.add_argument('-l', '--lowest_coverage', type=int, default='3', help='Minimum number of sequences in an LCB file')
    parser.add_argument('-r', '--refine_alignment', default='False', help='Refine MAUVE alignment of LCBs with MAFFT (Warning: perhaps harmful for distantly related taxa!)')
    parser.add_argument('-f', '--filter', default='True', help='Filter out low-coverage LCBs')
    parser.add_argument('-c', '--core', action='store_true', help='Only extract core SNPs (resemble parSNP)')
    parser.add_argument('-t', '--trim', default='True', help='Trim gap regions with TrimAl before SNP calling')
    parser.add_argument('-w', '--window', action='store_true', help='Launch GUI with pre-filled parameters')
    parser.add_argument('-v', '--version', action='store_true', help='Show version information')
    # parser.add_argument('--fasttree', default='FastTree', help='Path to FastTree executable')
    
    args = parser.parse_args()

    if args.version:
        print(f"DOGE version {__version__}")
        return

    if args.window:
        launch_gui_with_args(args)
    
    if not (args.input and args.output):
        raise Excption("input and output directorirs required!")
        return
    
    os.makedirs(args.output, exist_ok=True)
    mauve_output = os.path.join(args.output, "mauve_alignment.xmfa")
    lcb_dir = os.path.join(args.output, "lcbs")
    aligned_dir = os.path.join(args.output, "aligned_lcbs")
    trimmed_dir = os.path.join(args.output, "trimmed_lcbs")
    
    print("\n=== Running Mauve alignment ===")
    mauve_output_file = run_mauve(args.mauve, args.input, mauve_output)
    
    print("\n=== Splitting XMFA into LCBs ===")
    lcb_files = split_xmfa(mauve_output_file, lcb_dir)
    
    if args.refine_alignment == 'True':
        print("\n=== Refining LCB Alignments with MAFFT ===")
        aligned_files = align_all_lcbs(args.mafft, lcb_files, aligned_dir)
        lcb_files = aligned_files

    if args.filter == 'True' or args.core == 'True':
        print("\n=== Filtering LCB Coverages ===")
        filtered_lcb_files = filter_lcbs(args.lowest_coverage, lcb_files)
        lcb_files = filtered_lcb_files

    if args.trim == 'True' and args.core == 'False':
        print("\n=== Trimming Valid LCBs ===")
        trimmed_filtered_lcb_files = trim_lcbs(args.trimal, lcb_files, trimmed_dir)
        lcb_files = trimmed_filtered_lcb_files
    
    print("\n=== Creating supermatrix ===")
    supermatrix_file = os.path.join(args.output, "supermatrix.fasta")
    supermatrix, partitions = create_supermatrix(lcb_files, supermatrix_file)

    if args.core:
        print("\n=== Extracting Core SNPs ===")
        core_snps = extract_core_snps(supermatrix, partitions, args.output)
        nexus_file = core_snps
    else:
        print("\n=== Calling SNPs and creating partitioned NEXUS ===")
        nexus_file = call_snps(supermatrix, partitions, os.path.join(args.output, "snps"))

    
    print("\n=== DiscoverOrganelle pipeline completed successfully! ===")

if __name__ == "__main__":
    main()
    # Next time I may add a "--continue" function for breakpoints

class Plugin:
    # This is the plugin entry point for YRTools.
    def run(self):
        return DOGE_UI()