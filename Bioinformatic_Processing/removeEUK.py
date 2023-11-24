#!/usr/bin/env python3

import argparse
import os
import subprocess

def run_analysis(input_dir, output_dir, sample_name, kraken_db, threads):
    # Load necessary modules
    subprocess.run(["module", "load", "trimmomatic", "kraken2", "pigz"], check=True)

    # Create temporary directory
    if os.path.exists("./" + os.getlogin() + "/tmp") and os.path.isdir("./" + os.getlogin() + "/tmp"):
        tmp_dir = "./" + os.getlogin() + "/tmp"
    else:
        subprocess.run(["rm", "-rf", "./" + os.getlogin() + "/tmp"], check=True, stderr=subprocess.DEVNULL)
        os.makedirs("./" + os.getlogin() + "/tmp")
        tmp_dir = subprocess.check_output(["mktemp", "-d", "./" + os.getlogin() + "/tmp/XXXX"], universal_newlines=True).strip()

    os.environ["TMPDIR"] = tmp_dir
    os.environ["TMP"] = tmp_dir
    os.environ["TEMP"] = tmp_dir

    # Set input and output file names
    input_file_1 = os.path.join(input_dir, f"{sample_name}_1.fq.gz")
    input_file_2 = os.path.join(input_dir, f"{sample_name}_2.fq.gz")
    output_file_1 = os.path.join(output_dir, f"{sample_name}_filtered_1.fastq.gz")
    output_file_2 = os.path.join(output_dir, f"{sample_name}_filtered_2.fastq.gz")

    # Run Trimmomatic to trim low quality reads
    trim_cmd = ["trimmomatic", "PE", "-threads", str(threads), "-phred33", input_file_1, input_file_2, output_file_1, "/dev/null", output_file_2, "/dev/null", "LEADING:20", "TRAILING:20", "SLIDINGWINDOW:4:20", "MINLEN:50"]
    subprocess.run(trim_cmd, check=True)

    # Remove eukaryotic data from the metagenome using Kraken2
    kraken_cmd = ["/projects/mjolnir1/apps/conda/kraken2-2.1.2/bin/kraken2", "--db", kraken_db, "--threads", str(threads), "--gzip-compressed", "--paired", output_file_1, output_file_2]
    kraken_out = os.path.join(output_dir, f"{sample_name}_kraken_output.txt")
    with open(kraken_out, "w") as kraken_file:
        subprocess.run(kraken_cmd, check=True, stdout=kraken_file)

    # Filter out eukaryotic reads from the metagenome using Kraken2 output
    extract_cmd = ["/projects/mjolnir1/apps/bin/extract_kraken_readsBB.py", "-k", os.path.join(output_dir, f"{sample_name}_12_krk2_500.kraken2.gz"), "-s1", output_file_1, "-s2", output_file_2, "-o", os.path.join(output_dir, f"noEUK_{sample_name}_1.fq"), "-o2", os.path.join(output_dir, f"noEUK_{sample_name}_2.fq"), "--fastq-output", "--exclude", "--taxid", "2759", "--include-children", "-r", os.path.join(output_dir, f"{sample_name}_12_krk2_500.kraken2.report.gz")]
    subprocess.run(extract_cmd, check=True)

    # Compress noEUK files
    pigz_cmd = ["pigz", "-p", str(threads), os.path.join(output_dir, f"noEUK_{sample_name}_1.fq"), os.path.join(output_dir, f"noEUK_{sample_name}_2.fq")]
    subprocess.run(pigz_cmd, check=True)

    # Check if noEUK files exist before removing filtered files
    if os.path.exists(os.path.join(output_dir, f"noEUK_{sample_name}_1.fq.gz")) and os.path.exists(os.path.join(output_dir, f"noEUK_{sample_name}_2.fq.gz")):
        # Remove filtered files
        os.remove(output_file_1)
        os.remove(output_file_2)

        # Create done file
        with open(os.path.join(output_dir, f"done.{sample_name}"), "w") as done_file:
            done_file.write(f"{sample_name} is done\n")
    else:
        print("Error: noEUK files do not exist.")
        exit(1)

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run metagenomic analysis pipeline")
    parser.add_argument("--input", required=True, help="Input directory containing paired-end FASTQ files")
    parser.add_argument("--output", required=True, help="Output directory for filtered FASTQ files and Kraken2 output")
    parser.add_argument("--sample", required=True, help="Sample name (prefix of FASTQ file names)")
    parser.add_argument("--kraken-db", required=True, help="Path to Kraken2 database")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use (default: 1)")
    args = parser.parse_args()

    # Run analysis
    run_analysis(args.input, args.output, args.sample, args.kraken_db, args.threads)
