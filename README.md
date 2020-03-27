# synta : SYNTAxic annotation pipeline for procaryotes


synta is a procaryotic genome annotation tool that identifies basic genomic features and writes standard file formats as output.

It was designed with speed in mind, and can annotate an _E. coli_ genome in less than 10 seconds using parallel processing.
It will NOT produce functional annotation, which is usually the time-consumming step of procaryotic annotation pipelines, and will write only the demanded file formats (gff3 by default).

# Install
## Using git and pip
Make a new directory, go into it and run :

`git clone https://github.com/axbazin/synta.git . ; pip install . `
or
`git clone https://github.com/axbazin/synta.git . ; python setup.py install`

You will need to have the dependencies in your PATH, they are listed in the 'requirements.txt' file.
If you run synta without one of those dependencies, it will tell you which is missing.

# requirements
It should run easily on macOS and linux.
It requires :
- python >= 3.6
- prodigal
- infernal
- aragorn

# Basic usage

The simplest command is :
`synta --fna MY_FASTA_FILE`

This command will use a single cpu and write a gff3 formated file as output.
Using one more cpu will make it faster. Providing with even more does not make it much faster.

`synta --fna MY_FASTA_FILE --cpu 2`

If you need other file formats, there are 5 of them that are currently provided : fna, faa, ffn, gff and gbff.
To obtain them, use the --format option by listing the formats that you want, separated by a comma, as such :

`synta --fna MY_FASTA_FILE --format gbff,fna,faa`

The example above will provide with gbff, fna and faa file formats. The names must be strictly separated by a comma and nothing else.

Should you want to customize some elements of the workflow, other options are accessible through the help.

`synta -h`

# The workflow

Synta is just a wrapper for other tools that parses their outputs and writes them in more practical and standard file formats. It predicts ORF with Prodigal, tRNA with Aragorn and rRNA with Infernal (using RFAM models). Then, it writes to the different output formats asked by the user. By default, it will write a gff3 file only.
