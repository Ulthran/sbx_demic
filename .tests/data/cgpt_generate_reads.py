import numpy as np

COMPLEMENT_BASES = {
    "T": "A",
    "C": "G",
    "A": "T",
    "G": "C",
}

def generate_cyclical_read(genome, read_length, start_position):
    genome_length = len(genome)
    
    # Ensure start_position is within the genome boundaries
    start_position %= genome_length
    
    # Calculate the end position of the read
    end_position = (start_position + read_length) % genome_length
    
    # Handle the case where the read wraps around the genome
    if end_position <= start_position:
        read = genome[start_position:] + genome[:end_position]
    else:
        read = genome[start_position:end_position]
    
    return read

def generate_reverse_cyclical_read(genome, read_length, start_position):
    fwd_read = generate_cyclical_read(genome, read_length, start_position)
    transposed_read = [COMPLEMENT_BASES[x] for x in fwd_read]
    transposed_read.reverse()
    return "".join(transposed_read)

def write_reads(genome, fp, num_reads, amplitude):
    # Genome length and number of reads
    genome_length = len(genome)

    # Parameters for the sine wave
    amplitude = amplitude
    period = genome_length
    phase = -1 * np.pi / 2
    vertical_shift = 2 * amplitude

    # Generate the sine wave for coverage
    x = np.arange(0, genome_length)
    coverage = amplitude * np.sin(2 * np.pi * x / period + phase) + vertical_shift
    coverage = np.round(coverage).astype(int)

    # Generate random positions for paired reads
    read_positions = np.random.choice(genome_length, num_reads, p=coverage / np.sum(coverage))

    # Write reads to file
    read_length = 100
    #assert num_reads * read_length > 6 * genome_length, f"{num_reads * read_length} must be greater than {6 * genome_length}"

    with open(fp.format("1"), "a") as f_r1, open(fp.format("2"), "a") as f_r2:
        for i, start in enumerate(read_positions):
            f_r1.write(f"@{i}\n{generate_cyclical_read(genome, read_length, start)}\n+\n{'~' * read_length}\n")
            f_r2.write(f"@{i}\n{generate_reverse_cyclical_read(genome, read_length, start)}\n+\n{'~' * read_length}\n")

def do_the_thing(ref, out, num_reads, amplitude):
    genome = ""
    with open(ref) as f_Ecoli:
        f_Ecoli.readline() # Header
        genome = "".join([line.strip() for line in f_Ecoli.readlines()])

    write_reads(genome, out, num_reads, amplitude)

do_the_thing("/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/reference/Ecoli.fasta", "/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/readsFULLSIM/3.1_R{}.fastq", 5000, 5)
do_the_thing("/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/reference/Bfragilis.fasta", "/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/readsFULLSIM/3.1_R{}.fastq", 5000, 5)
#do_the_thing("/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/reference/akk-genome.fasta", "/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/readsFULLSIM/3.1_R{}.fastq", 20000, 5)
do_the_thing("/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/reference/Ecoli.fasta", "/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/readsFULLSIM/3.2_R{}.fastq", 5000, 2)
do_the_thing("/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/reference/Bfragilis.fasta", "/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/readsFULLSIM/3.2_R{}.fastq", 5000, 2)
#do_the_thing("/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/reference/akk-genome.fasta", "/mnt/d/Penn/sunbeam/extensions/sbx_demic/.tests/data/readsFULLSIM/3.2_R{}.fastq", 20000, 2)