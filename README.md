# sbx_demic

[Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension to perform co-assembly of reads from arbitrary groups of samples from a given project using [Megahit](https://github.com/voutcn/megahit) and run [DEMIC](https://sourceforge.net/projects/demic/files/)

*Important* -- This was tested with modified demic scripts version 1.0.2, please check the sourceforge.net site for any updates. One change I made was to save a couple of R data files for "manual" analysis in Rstudio if one so desires. See the git history for any other changes.

## Installation

0. Install sbx_coassembly from https://github.com/sunbeam-labs/sbx_coassembly
1. git clone https://github.com/sunbeam-labs/sbx_demic
2. cp sbx_demic $SUNBEAM_DIR/extensions/
3. cat sunbeam/extensions/sbx_demic/config.yml >> sunbeam_config.yml (the config.yml that your are using for your given project)
4. Edit sunbeam_config.yml to have desired parameters
    - Make sure that
    - *all.paired_end* is true if you have paired end reads
    - *sbx_demic.extras* has any parameters you want to pass to DEMIC.pl

## Configuration

If you'd like to coassemble all of your samples, no further input is needed!

If no grouping file is specified this extension by default co-assembles all samples. If you'd like to group specific samples, you need to provide a mapping file, then point to that mapping file in your config file. For example, if you have three samples from individual A (A_d1, A_d2, A_d3) and two from individual B (B_d1, B_d3), and you'd like to make one coassembly for each individual, you first need a mapping file like this (I'll call it `mapping.yml`):

    A: ['A_d1', 'A_d2', 'A_d3']
    B: ['B_d1', 'B_d3']

Here, the bracketed items in the list (i.e. `'A_d1'` or `'B_d3'` are your full sample names and must match the sample names Sunbeam knows. (Look in your `samples.csv` file if you're not sure of the names). The keys at the beginning of the lines are the group names--you can use any valid Python variable name here. After making this mapping file, make sure to edit your config file to point to this mapping file:

```
...
sbx_coassembly:
  threads: 4
  group_file: '/path/to/mapping.yml'
```

## Running

To generate coassembled samples, create a project, define your groupings, and use the `all_coassemble` target:

    sunbeam init --data_fp /path/to/reads/ /path/to/project/
    printf "A: ['A_d1', 'A_d2', 'A_d3']\nB: ['B_d1', 'B_d3']" > mapping.yml
    sunbeam config modify -i -f /path/to/project/sunbeam_config.yml -s 'sbx_coassembly: {{group_file: {/path/to/mapping.yml}}}'
    sunbeam run --profile /path/to/project/ all_coassemble

To run demic:

    sunbeam run --profile /path/to/project/ all_demic

See the wiki for a walkthrough of the full setup/installation and run processes.

### Trouble-shooting

- If you have trouble running with "--use-conda" it may be best just to install the needed packages into the sunbeam environment (e.g. `conda activate sunbeam && conda install --file sbx_demic_env.yml`) and then run sunbeam in two pieces (`sunbeam run all_coassemble --use-conda --configfile /path/to/config {rest of arguments for sunbeam}` and then `sunbeam run all_demic --configfile /path/to/config {rest of arguments for sunbeam}`). Please also create a new issue on GitHub detailing the error.

## References

- Sunbeam: https://github.com/sunbeam-labs/sunbeam

- DEMIC software: https://sourceforge.net/projects/demic/files/

- DEMIC publication: https://www.nature.com/articles/s41592-018-0182-0

## Legacy Installation

For sunbeam versions <3 or if `sunbeam extend` isn't working, you can use `git` directly to install an extension:

    git clone https://github.com/sunbeam-labs/sbx_coassembly.git extensions/sbx_coassembly

and then include it in the config for any given project with:

    cat extensions/sbx_coassembly/config.yml >> /path/to/project/sunbeam_config.yml
