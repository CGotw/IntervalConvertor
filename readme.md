# IntervalConvertor

This code repository contains the following three tools：

- **Interval Conversion**

  The **Interval Conversion** is a tool which could convert the location of fragments that from other genomes into the location of‘Fuji’unphased genome using the message in the gaf file created by minigraph with‘Fuji’unphased genome and the other genome.
- **Species Passing Information**

  Run the Species Passing Information to get the species information through the gene sequences in the input file. For the reference genome 'Fuji', program queries directly. For the non-reference genome, the program first performs Interval Conversion, and then queries.
- **SV Calling**

  SV Calling is using for classifying the structure variations in the graph file by four types,  (1) Bi-allelic (insertions/deletions): with two paths in the bubble, and the length of the shorter path should be below 50bp, insertion or deletion was defined based on the result of comparison with reference path; (2) Divergent: with two paths in the bubble, and the length of the shorter path should be above 50bp; (3) Multi-allelic: with more than two paths in the bubble after excluding those non-reference paths below 50bp.

## Installation

Some input files of the tools are generated using minigraph，therefore, minigraph needs to be installed first.

```bash
git clone https://github.com/lh3/minigraph
cd minigraph && make
```

Running environment required to use the **conda** installer.

```bash
# The following commands install ”xxx“ in a new conda environment called `xxx`
conda create -n intervalConvertor python=3.9
conda activate intervalConvertor
git clone https://github.com/CGotw/IntervalConvertor.git
cd IntervalConvertor
pip install -r requirements.txt
```

## Usage

### 1. Interval Conversion

Five required parameters need to be set when running Interval Conversion：

- gaf\_dir: This is a folder with gaf files of all the genomes in the pan-genome graph. Each gaf file is created with a query genome and the pan-genome graph file, and the paramater is: `minigraph -x lr pangenome.gfa query.fasta > query.gaf`.
- gfa\_path: The gfa file is created by minigraph with the parameter 'xggs', with one genome as reference and all other query genomes, including two types of messages: the description of all the segments in the gfa, and how the segments could connect with each other.
- query\_input\_path：A list file, containing the chromosome ID , the start and end location and line number of each fragment.
- query\_gaf\_path：The gaf file created with the query genome and the pan-genome graph file, including at least the reference and query genomes.
- output\_file\_path：Save the result of the Interval Conversionl.
- reference_name: The reference genome name for GFA File.

**Run command**

```bash
python "interval conversion.py" \
    --gaf_dir "your_gaf_dir" \
    --gfa_path "your_gfa_path" \
    --query_input_path "your_query_input_path" \
    --query_gaf_path "your_query_gaf_path" \
    --output_file_path "your_output_file_path" \
    --reference_name "reference_name"
```

### 2. Species Passing Information

Six required parameters need to be set when running Species Passing Information：

- gaf\_dir: Ditto
- gfa\_path: Ditto
- interval\_conversion：A Boolean value. True means that the input sequence will run  Interval Conversion first. The default value is False.
- query\_input\_path：A list file, containing the chromosome ID and the start and end location of each fragment.
- query\_gaf\_path：The gaf file created with the query genome and the pan-genome graph file, including at least the reference and query genomes. If interval\_conversion is false, this parameter is meaningless.
- output\_file\_path：Save the result of the Interval Conversionl.

**Run command**

```bash
python "species passing information.py" \
    --gaf_dir "your_gaf_dir" \
    --gfa_path "your_gfa_path" \
    --interval_conversion True \
    --query_input_path "your_query_input_path" \
    --query_gaf_path "your_query_gaf_path" \
    --output_file_path "your_output_file_path" \
```

### 3. SV Calling

Four required parameters need to be set when running SV Calling：

- gaf\_dir: Ditto
- gfa\_path: Ditto
- query\_bed\_path：This file is created with the pan-genome graph file, using for describing the bubbles in the graph. The parameter is: gfatools bubble graph.gfa > var.bed
- output\_dir：A directory path where the SV Calling results are saved.

**Run command**

```bash
python "SV calling.py" \
    --gaf_dir "your_gaf_dir" \
    --gfa_path "your_gfa_path" \
    --query_bed_path "your_bed_path" \
    --output_dir "your_output_dir" \
   
```

*last update:  4/19/2023*
