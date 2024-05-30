import subprocess

def main(input_fasta , output_file ):
    cd_hit_executable = "/home/dina/cd-hit-v4.8.1-2019-0228/cd-hit"

    # Define CD-HIT command
    cd_hit_command = [cd_hit_executable, "-i", input_fasta, "-o", output_file, "-c", "0.7", "-n", "5"]

    # Execute CD-HIT command
    subprocess.run(cd_hit_command)

    print("CD-HIT finished successfully.")

if __name__ == "__main__":
    main(input_fasta , output_file )