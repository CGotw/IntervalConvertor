# IntervalConvertor

Here we provide our bioinformatics tool, intervalConvertor. Using this tool, locations on any genome assemblies used to build the pangenome can be converted into standardized coordinates on the specified reference genome. 

## Installation

Some input files of the tools are generated using minigraph，therefore, minigraph needs to be installed first.

```bash
git clone https://github.com/lh3/minigraph
cd minigraph && make
```

Running environment required to use the **conda** installer.

```bash
#  The following commands install "Intervalconvertor" in a new conda environment called intervalconvertor
conda create -n intervalConvertor python=3.9
conda activate intervalConvertor
git clone https://github.com/CGotw/IntervalConvertor.git
cd IntervalConvertor
pip install -r requirements.txt
```

## Usage

Four required parameters need to be set when running IntervalConvertor：

- query\_input\_path：The loci chosen in the query genome.
- query\_gaf\_path：The alignment path between the query and reference genomes based on the pangenome graph.
- output\_file\_path：Save the result of the IntervalConvertor.
- reference_name: Reference genome name.

**Run command**

```bash
python "intervalconvertor.py" \
    --query_input_path "your_query_input_path" \
    --query_gaf_path "your_query_gaf_path" \
    --output_file_path "your_output_file_path" \
    --reference_name "reference_name"
```

*last update:  1/22/2025*
