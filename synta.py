#!/usr/bin/env python3

from string import ascii_uppercase
from random import choice
import ast
import tempfile
from subprocess import Popen, PIPE
from multiprocessing import Pool
from io import TextIOWrapper
import argparse
import gzip
import time
import os
import logging
import sys

# minimum python version to use is 3.6, there are features that don't work with older versions in this code. (most notably, f-strings)
assert sys.version_info >= (3, 6)


class gene:
    def __init__(self, ID, contig, start, stop, strand, geneType, product=None, inference=None):
        self.ID = ID
        self.contig = contig
        self.start = int(start)
        self.stop = int(stop)
        self.type = geneType
        self.strand = strand
        self.inference = inference
        self.product = product

    def get_gff(self, version_numbers):
        """
            Takes the software's version numbers used for different types of genes (key:  gene type, value:  tuple with version string andsoftware name). Writes the gene's information in gff3 format
        """
        version, soft = version_numbers[self.type]

        Feature = f"{self.contig}\t{soft}: {version}\tgene\t{self.start}\t{self.stop}\t1\t{self.strand}\t0\tID=gene_{self.ID}\n"
        Feature += f'{self.contig}\t{soft}: {version}\t{self.type}\t{self.start}\t{self.stop}\t1\t{self.strand}\t0\tID={self.ID};Parent=gene_{self.ID};inference={self.inference}'
        Feature += f";product={self.product}\n" if self.product else "\n"
        return Feature

    def get_faa(self, dna_seq):
        """Takes a dna sequence string and returns the gene object's ID with the string  translated in faa format"""
        faa_seq = translate(dna_seq)
        seq = f">{self.ID}\n"
        j = 0
        seq += faa_seq[j: j+60] + "\n"
        j = 60
        while j < len(faa_seq):
            seq += faa_seq[j: j+60] + "\n"
            j += 60
        return seq

    def get_ffn(self, dna_seq):
        """Takes a dna sequence string and returns the gene object's ID with the string in ffn format"""
        seq = f">{self.ID}\n"
        j = 0
        seq += dna_seq[j: j+60] + "\n"
        j = 60
        while j < len(dna_seq):
            seq += dna_seq[j: j+60] + "\n"
            j += 60
        return seq


def reverse_complement(seq):
    """ reverse complement the given dna sequence """
    complement = {'A':  'T', 'C':  'G', 'G':  'C', 'T':  'A', 'N': 'N', 'R': 'Y', 'Y': 'R',
                  'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    # see https: //www.bioinformatics.org/sms/iupac.html for the code.
    # complement = {'A':  'T', 'C':  'G', 'G':  'C', 'T':  'A', 'N': 'N' } ## basic
    rcseq = ""
    for i in reversed(seq):
        rcseq += complement[i]
    return rcseq


def translate(seq):
    """ translates the given dna sequence with table code 11 of the ncbi (bacteria)"""
    start_table = {
        'ATA': 'M', 'ATC': 'M', 'ATT': 'M', 'ATG': 'M', 'ATN': 'M', 'ATR': 'M', 'ATY': 'M', 'ATS': 'M', 'ATW': 'M', 'ATK': 'M', 'ATM': 'M', 'ATB': 'M', 'ATD': 'M', 'ATH': 'M', 'ATV': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'ACN': 'T', 'ACR': 'T', 'ACY': 'T', 'ACS': 'T', 'ACW': 'T', 'ACK': 'T', 'ACM': 'T', 'ACB': 'T', 'ACD': 'T', 'ACH': 'T', 'ACV': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AAY': 'N', 'AAR': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R', 'AGY': 'S', 'AGR': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'M', 'CTT': 'L', 'CTB': 'L', 'CTY': 'L', 'CTW': 'L', 'CTM': 'L', 'CTH': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CCN': 'P', 'CCR': 'P', 'CCY': 'P', 'CCS': 'P', 'CCW': 'P', 'CCK': 'P', 'CCM': 'P', 'CCB': 'P', 'CCD': 'P', 'CCH': 'P', 'CCV': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CAY': 'H', 'CAR': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CGN': 'R', 'CGR': 'R', 'CGY': 'R', 'CGS': 'R', 'CGW': 'R', 'CGK': 'R', 'CGM': 'R', 'CGB': 'R', 'CGD': 'R', 'CGH': 'R', 'CGV': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'M', 'GTT': 'V', 'GTB': 'V', 'GTY': 'V', 'GTW': 'V', 'GTM': 'V', 'GTH': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GCN': 'A', 'GCR': 'A', 'GCY': 'A', 'GCS': 'A', 'GCW': 'A', 'GCK': 'A', 'GCM': 'A', 'GCB': 'A', 'GCD': 'A', 'GCH': 'A', 'GCV': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GAY': 'D', 'GAR': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GGN': 'G', 'GGR': 'G', 'GGY': 'G', 'GGS': 'G', 'GGW': 'G', 'GGK': 'G', 'GGM': 'G', 'GGB': 'G', 'GGD': 'G', 'GGH': 'G', 'GGV': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TCN': 'S', 'TCR': 'S', 'TCY': 'S', 'TCS': 'S', 'TCW': 'S', 'TCK': 'S', 'TCM': 'S', 'TCB': 'S', 'TCD': 'S', 'TCH': 'S', 'TCV': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'M', 'TTY': 'F',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ATY': 'I', 'ATH': 'I', 'ATW': 'I', 'ATM': 'I',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'ACN': 'T', 'ACR': 'T', 'ACY': 'T', 'ACS': 'T', 'ACW': 'T', 'ACK': 'T', 'ACM': 'T', 'ACB': 'T', 'ACD': 'T', 'ACH': 'T', 'ACV': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AAY': 'N', 'AAR': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R', 'AGY': 'S', 'AGR': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CTN': 'L', 'CTR': 'L', 'CTY': 'L', 'CTS': 'L', 'CTW': 'L', 'CTK': 'L', 'CTM': 'L', 'CTB': 'L', 'CTD': 'L', 'CTH': 'L', 'CTV': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CCN': 'P', 'CCR': 'P', 'CCY': 'P', 'CCS': 'P', 'CCW': 'P', 'CCK': 'P', 'CCM': 'P', 'CCB': 'P', 'CCD': 'P', 'CCH': 'P', 'CCV': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CAY': 'H', 'CAR': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CGN': 'R', 'CGR': 'R', 'CGY': 'R', 'CGS': 'R', 'CGW': 'R', 'CGK': 'R', 'CGM': 'R', 'CGB': 'R', 'CGD': 'R', 'CGH': 'R', 'CGV': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GTN': 'V', 'GTR': 'V', 'GTY': 'V', 'GTS': 'V', 'GTW': 'V', 'GTK': 'V', 'GTM': 'V', 'GTB': 'V', 'GTD': 'V', 'GTH': 'V', 'GTV': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GCN': 'A', 'GCR': 'A', 'GCY': 'A', 'GCS': 'A', 'GCW': 'A', 'GCK': 'A', 'GCM': 'A', 'GCB': 'A', 'GCD': 'A', 'GCH': 'A', 'GCV': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GAY': 'D', 'GAR': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GGN': 'G', 'GGR': 'G', 'GGY': 'G', 'GGS': 'G', 'GGW': 'G', 'GGK': 'G', 'GGM': 'G', 'GGB': 'G', 'GGD': 'G', 'GGH': 'G', 'GGV': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TCN': 'S', 'TCR': 'S', 'TCY': 'S', 'TCS': 'S', 'TCW': 'S', 'TCK': 'S', 'TCM': 'S', 'TCB': 'S', 'TCD': 'S', 'TCH': 'S', 'TCV': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L', 'TTY': 'F', 'TTR': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TAR': '*', 'TAY': 'Y',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W', 'TGY': 'C',
    }
    # code:  https: //www.bioinformatics.org/sms/iupac.html
    protein = ""
    if len(seq) % 3 == 0:
        protein = start_table[seq[0: 3]]
        for i in range(3, len(seq), 3):
            codon = seq[i: i + 3]
            try:
                protein += table[codon]
            except KeyError:  # codon was not planned for. Probably can't determine it.
                # print(codon)
                protein += 'X'  # X is for unknown
    else:
        print(len(seq))
        raise IndexError(
            "Given sequence length modulo 3 was different than 0, which is unexpected.")
    return protein


def is_compressed(file_or_file_path):
    """
        Checks is a file, or file path given is compressed or not
    """
    file = file_or_file_path
    if type(file) == str:
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except:
            return False
    if file.read(2).startswith(b'\x1f\x8b'):
        return True
    file.close()
    return False


def read_compressed_or_not(file_or_file_path):
    """
        Reads a file, compressed or not.
        Copied from http: //www.github.com/ggautreau/PPanGGOLiN.git's utils.py.
    """
    file = file_or_file_path
    if type(file) == str:
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except:
            return(file)
    if file.read(2).startswith(b'\x1f\x8b'):
        file.seek(0)
        logging.getLogger().info("Uncompressing the file: '" + file.name + "' ...")
        return(TextIOWrapper(gzip.open(filename=file, mode="r")))
    else:
        file.close()
        file = open(file.name, "r")
        return(file)


def write_compress_or_not(file_path, compress):
    """
        Returns a file-like object, compressed or not.
    """
    if compress:
        return gzip.open(file_path + ".gz", mode="wt")
    else:
        return open(file_path, "w")


def launch_aragorn(fnaFile, locustag):
    """ 
        launches Aragorn to annotate tRNAs. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    logging.getLogger().debug("Running Aragorn(tRNA).")
    cmd = ["aragorn", "-t", "-gcbact", "-l", "-w", fnaFile]
    logging.getLogger().debug(f"command for Aragorn:  '{' '.join(cmd)}'")
    p = Popen(cmd, stdout=PIPE)
    # loading the whole thing, reverting it to 'pop' in order.
    fileData = p.communicate()[0].decode().split("\n")[:: -1]
    geneObjs = []
    c = 0
    while len(fileData) != 0:
        line = fileData.pop()
        if line.startswith(">"):
            header = line.replace(">", "").split()[0]
            fileData.pop()  # then next line must be removed too.
        elif len(line) > 0:  # if the line isn't empty, there's data to get.
            lineData = line.split()
            start, stop = ast.literal_eval(lineData[2].replace("c", ""))
            c += 1
            geneObjs.append(gene(ID=locustag+'_tRNA_'+str(c).zfill(3),
                                 contig=header,
                                 start=start,
                                 stop=stop,
                                 strand="-" if lineData[2].startswith(
                                     "c") else "+",
                                 geneType="tRNA",
                                 inference="COORDINATES: profile: Aragorn",
                                 product=lineData[1] + lineData[4]))
    logging.getLogger().info(
        f"Done with Aragorn(tRNA). There are {len(geneObjs)} tRNAs.")
    return geneObjs


def launch_prodigal(fnaFile, locustag):
    """ 
        launches Prodigal to annotate CDS. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    logging.getLogger().debug("Running Prodigal(CDS).")
    cmd = ["prodigal", "-f", "sco", "-c", "-i", fnaFile, "-p", "single", "-q"]
    logging.getLogger().debug(f"command for Prodigal:  '{' '.join(cmd)}'")
    p = Popen(cmd, stdout=PIPE)

    geneObjs = []
    c = 0
    for line in p.communicate()[0].decode().split("\n"):
        if line.startswith("# Sequence Data: "):
            for data in line.split(";"):
                if data.startswith("seqhdr"):
                    header = data.split("=")[1].replace('"', "").split()[0]
                    # print(header)

        elif line.startswith(">"):
            c += 1
            lineData = line[1:].split("_")  # not considering the '>'
            geneObjs.append(gene(ID=locustag + "_CDS_" + str(c).zfill(4),
                                 contig=header,
                                 start=lineData[1],
                                 stop=lineData[2],
                                 strand=lineData[3],
                                 geneType="CDS",
                                 inference="ab initio prediction: Prodigal"))

    logging.getLogger().info(
        f"Done with Prodigal(CDS). There are {len(geneObjs)} CDSs.")
    return geneObjs


def launch_infernal(fnaFile, locustag):
    """ 
        launches Infernal in hmmer-only mode to annotate rRNAs. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    logging.getLogger().debug("Running Infernal(rRNA).")
    tmpFile = tempfile.NamedTemporaryFile(mode="r")
    cmd = ["cmscan", "--tblout", tmpFile.name, "--hmmonly", "--cpu",
           str(1), "--noali", f"{os.path.dirname(os.path.realpath(__file__))}/cmDB/rRNA_bact.cm", fnaFile]
    logging.getLogger().debug(f"command for Infernal:  '{' '.join(cmd)}'")
    p = Popen(cmd, stdout=open(os.devnull, "w"), stderr=PIPE)
    err = p.communicate()[1].decode().split()
    if err != []:
        if err[0] == 'Error: ':
            raise Exception(
                f"Infernal (cmscan) failed with error:  '{ ' '.join(err) }'. If you never used this script, you should press the .cm file using cmpress executable from Infernal. You should find the file in '{os.path.dirname(os.path.realpath(__file__))}/cmDB/'.")
        raise Exception(
            f"An error occurred with Infernal. Error is:  '{ ' '.join(err) }'.")
    # never managed to test what happens if the .cm files are compressed with a 'bad' version of infernal, so if that happens you are on your own.

    geneObjs = []
    c = 0
    for line in tmpFile:
        if not line.startswith("#"):
            c += 1
            lineData = line.split()
            strand = lineData[9]
            if strand == "-":
                start = lineData[8]
                stop = lineData[7]
            else:
                start = lineData[7]
                stop = lineData[8]
            geneObjs.append(gene(ID=locustag + "_rRNA_" + str(c).zfill(3),
                                 contig=lineData[2],
                                 start=start,
                                 stop=stop,
                                 strand=strand,
                                 geneType="rRNA",
                                 product=" ".join(lineData[17:])))

    logging.getLogger().info(
        f"Done with Infernal(rRNA). There are {len(geneObjs)} rRNAs.")
    return geneObjs


def write_output(outputbasename, contigs, genes, compress, format, cpu, versions):
    """
        Writes output file in the given formats.
        Takes in the basename for the output
        the contigs of the fna in a dictionnary
        the gene list to write sorted by contig size and start position
        whether to compress the output or not
        the formats to write as a string separated by ','
        the number of cpus to use.
        the version strings for the different softwares. (obtained with check_versions() )
    """
    with Pool(processes=cpu) as p:
        for forms in format.split(","):
            if forms == "gff":
                wgff = p.apply_async(func=write_gff, args=(
                    outputbasename, contigs, genes, compress, versions))
            elif forms == "fna":
                wfna = p.apply_async(func=write_fna, args=(
                    outputbasename, contigs, compress))
            elif forms == "ffn":
                wffn = p.apply_async(func=write_ffn, args=(
                    outputbasename, contigs, genes, compress))
            elif forms == "faa":
                wfaa = p.apply_async(func=write_faa, args=(
                    outputbasename, contigs, genes, compress))

        # there is a better way of writing this 'wait until all processes are done', isn't there ?
        for forms in format.split(","):
            if forms == "gff":
                wgff.get()
            elif forms == "fna":
                wfna.get()
            elif forms == "ffn":
                wffn.get()
            elif forms == "faa":
                wfaa.get()


def read_fasta(fnaFile):
    """
        Reads a fna file and stores it in a dictionnary with contigs as key and sequence as value.
    """
    logging.getLogger().debug("Reading fasta file and extracting sequence and contig names.")
    contigs = {}
    contig_name = ""
    contig_seq = ""
    for line in fnaFile:
        if line.startswith('>'):
            if contig_name != "" and contig_seq != "":
                contigs[contig_name] = contig_seq
            contig_seq = ""
            contig_name = line.split()[0][1:]
        else:
            contig_seq += line.strip()

    # processing the last contig
    if contig_name != "" and contig_seq != "":
        contigs[contig_name] = contig_seq
    logging.getLogger().debug(
        f"Done reading fasta. There are {len(contigs)} sequences.")
    return contigs


def write_gff(output, contigs, genes, compress, versions):
    """
        Writes a gff formated file.
    """
    logging.getLogger().debug("Writting GFF file ...")
    outfile = write_compress_or_not(output + ".gff", compress)
    outfile.write("##gff-version 3\n")
    precContig = ""
    for gene in genes:
        if gene.contig != precContig:
            outfile.write(
                f"##sequence-region {gene.contig} 1 {len(contigs[gene.contig])}\n")
            # GenBank origin is only for PanGBank. The circularity could be guessed with a gbff input.
            outfile.write(
                f"{gene.contig}\tGenBank\tregion\t1\t{len(contigs[gene.contig])}\t.\t+\t.\tID={gene.contig}\n")
        outfile.write(gene.get_gff(versions))
        precContig = gene.contig
    outfile.close()
    logging.getLogger().debug("Done writing GFF file.")


def write_ffn(output, contigs, genes, compress):
    """
        Writes a ffn formated file.
    """
    logging.getLogger().debug("Writting FFN file ...")
    outfile = write_compress_or_not(output + ".fnn", compress)
    for gene in genes:
        # ffn is only for coding sequences ! (even if some softs also include RNA sequences ... It was not designed as such)
        if gene.type == "CDS":
            if gene.strand == "+":
                outfile.write(gene.get_ffn(
                    contigs[gene.contig][gene.start-1: gene.stop]))
            elif gene.strand == "-":
                outfile.write(gene.get_ffn(reverse_complement(
                    contigs[gene.contig][gene.start-1: gene.stop])))
    outfile.close()
    logging.getLogger().debug("Done writing FFN file.")


def write_faa(output, contigs, genes, compress):
    """
        Writes a faa formated file.
    """
    logging.getLogger().debug("Writting FAA file ...")
    outfile = write_compress_or_not(output + ".faa", compress)
    for gene in genes:
        if gene.type == "CDS":  # faa is only for coding sequences !
            if gene.strand == "+":
                outfile.write(gene.get_faa(
                    contigs[gene.contig][gene.start-1: gene.stop]))
            elif gene.strand == "-":
                outfile.write(gene.get_faa(reverse_complement(
                    contigs[gene.contig][gene.start-1: gene.stop])))
    outfile.close()
    logging.getLogger().debug("Done writing FAA file.")


def write_fna(output, contigs, compress):
    """
        Writes a fna formated file.
    """
    logging.getLogger().debug("Writting FNA file ...")
    outfile = write_compress_or_not(output + ".fna", compress)
    for header in sorted(contigs.keys(), key=lambda x:  len(contigs[x]), reverse=True):
        outfile.write(f">{header}\n")
        j = 0
        while j < len(contigs[header]):
            outfile.write(f"{contigs[header][j: j+60]}\n")
            j += 60
    logging.getLogger().debug("Done writing FNA file.")
    outfile.close()


def write_tmp_fasta(contigs):
    """
        Writes a temporary fna formated file, and returns the file-like object.

        This is for the cases where the given file is compressed, then we write a temporary file for the annotation tools to read from. The file will be deleted when close() is called.
    """
    tmpFile = tempfile.NamedTemporaryFile(mode="w")
    for header in contigs.keys():
        tmpFile.write(f">{header}\n")
        logging.getLogger().debug(
            f"reading sequence '{header}' of length '{len(contigs[header])}'")
        j = 0
        while j < len(contigs[header]):
            tmpFile.write(contigs[header][j: j+60]+"\n")
            j += 60
    tmpFile.flush()  # force write what remains in the buffer.
    return tmpFile


def syntaxic_annotation(fastaFile, cpu, norna, locustag):
    """
        Runs the different softwares for the syntaxic annotation.

        Takes in the file-like object containing the uncompressed fasta sequences to annotate
        the number of cpus that we can use.
        whether to annotate rna or not
        the locustag to give gene IDs.
    """
    # launching tools for syntaxic annotation
    genes = []
    logging.getLogger().info(f"Launching syntaxic annotation.")

    with Pool(processes=cpu) as p:  # launching a process pool with all the accessible processes

        proGenes = p.apply_async(func=launch_prodigal,
                                 args=(fastaFile.name, locustag))
        logging.getLogger().debug("Started the process launching Prodigal")
        if not norna:
            araGenes = p.apply_async(
                func=launch_aragorn, args=(fastaFile.name, locustag))
            logging.getLogger().debug("Started the process launching Aragorn")
            infGenes = p.apply_async(
                func=launch_infernal, args=(fastaFile.name, locustag))
            logging.getLogger().debug("Started the process launching Infernal")
            genes.extend(araGenes.get())
            logging.getLogger().debug("Got the results from the process for Aragorn")
            genes.extend(infGenes.get())
            logging.getLogger().debug("Got the results from the process for Infernal")
        genes.extend(proGenes.get())
        logging.getLogger().debug("Got the results from the process for Prodigal")

    fastaFile.close()  # closing either tmp file or original fasta file.
    return genes


def overlap_filter(genes, contigs):
    """
        Removes the CDS that overlap with RNA genes.
    """
    tmpGenes = sorted(genes, key=lambda x:  (-len(contigs[x.contig]), x.start))
    rmGenes = set()

    for i, gene_i in enumerate(tmpGenes):
        if i+1 < len(tmpGenes):
            gene_j = tmpGenes[i+1]
            if gene_i.type != "CDS" and gene_j.type == "CDS" and gene_i.stop > gene_j.start and gene_i.contig == gene_j.contig:
                rmGenes.add(gene_j)
            elif gene_i.type == "CDS" and gene_j.type != "CDS" and gene_i.stop > gene_j.start and gene_i.contig == gene_j.contig:
                rmGenes.add(gene_i)

    for gene in rmGenes:
        tmpGenes.remove(gene)
    logging.getLogger().info(
        f"{len(rmGenes)} CDS were removed due to overlapping RNA features.")
    return tmpGenes


def mk_basename(output, afile):
    """
        Extracts the basename from a file. (name with the extension(s)))
    """
    outbasename = output
    if output[-1] != "/":
        outbasename += "/"
    if is_compressed(afile):
        outbasename += "".join(os.path.basename(afile).split(".")[: -2])
    else:
        outbasename += "".join(os.path.basename(afile).split(".")[: -1])
    return outbasename


def check_versions():
    """
        Checks if the programs we need exist in the PATH, and their version number. 
        Returns a dictionnary with the type of annotated genes as key, and a tuple with the version string and the software name as value.
    """
    # Prodigal
    logging.getLogger().info("Checking software versions")
    try:
        proVersion = Popen(
            ["prodigal", "-v"], stderr=PIPE).communicate()[1].decode().split()[1][1: -1].split(".")
    except FileNotFoundError:
        raise Exception(
            "Prodigal was not found in PATH. Please install it and/or add it to your PATH.")
    if not (int(proVersion[0]) == 2 and int(proVersion[1]) >= 6 and int(proVersion[2]) >= 2):
        raise Exception(
            f"Prodigal version is not good. It is read as {'.'.join(proVersion)} while version 2.6.2 at least is expected.")
    logging.getLogger().info(f"Prodigal's version is '{'.'.join(proVersion)}'")
    try:
        araVersion = Popen(["aragorn", "-h"], stdout=PIPE).communicate()[
            0].decode().split("\n")[1][9: 15].split(".")
    except FileNotFoundError:
        raise Exception(
            "Aragorn was not found in PATH. Please install it and/or add it to your PATH.")

    if not (int(araVersion[0]) == 1 and int(araVersion[1]) >= 2 and int(araVersion[2]) >= 38):
        raise Exception(
            f"Aragorn version is not good. It is read as {'.'.join(araVersion)} while version at least 1.2.38 is expected.")
    logging.getLogger().info(f"Aragorn's version is '{'.'.join(araVersion)}'")
    try:
        infVersion = Popen(["cmscan", "-h"], stdout=PIPE).communicate()[
            0].decode().split('\n')[1][11: 16].split(".")
    except FileNotFoundError:
        raise Exception(
            "Infernal (and notably its executable cmscan) was not found in PATH. Please install it and/or add it to your PATH.")
    if not (int(infVersion[0]) == 1 and int(infVersion[1]) >= 1 and int(infVersion[2]) >= 2):
        raise Exception(
            f"Aragorn version is not good. It is read as {'.'.join(infVersion)} while version at least 1.1.2 is expected.")
    logging.getLogger().info(f"Infernal's version is '{'.'.join(infVersion)}'")
    # indicating what soft is for what type of gene, and the soft's version number.
    # the version number, and the software name as it will be displayed in the output files.
    return {'CDS': ('.'.join(proVersion), "Prodigal"), 'tRNA': ('.'.join(araVersion), "ARAGORN"), 'rRNA': ('.'.join(infVersion), "INFERNAL")}


def read_gbff(gbffFile):
    """
        Reads a gbff from a file-like object.
        returns contigs with IDs as key and sequences as value.
    """
    contigs = {}
    # revert the order of the file, to read the first line first.
    lines = gbffFile.readlines()[::-1]
    while len(lines) != 0:
        line = lines.pop()
        # beginning of contig.## if multicontig, should happen right after '//'.
        if line.startswith('LOCUS'):
            while not line.startswith('FEATURES'):
                if line.startswith('VERSION'):
                    contigID = line[12:].strip()
                line = lines.pop()
        # start of the feature object.
        line = lines.pop()
        while not line.startswith("ORIGIN"):
            line = lines.pop()
        line = lines.pop()  # first sequence line.
        sequence = ""
        while not line.startswith('//'):
            sequence += line[10:].replace(' ', '').strip().upper()
            line = lines.pop()
        # ended the contig/sequence.
        contigs[contigID] = sequence
    return contigs


def cmdLine():
    """
        Functions that defines the command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--fna',  required=False, type=str,
                        help="fasta(.gz) file to annotate. Will be used if given.")
    parser.add_argument('--gbff', required=False, type=str,
                        help=".gbff (or .gbk) file to annotate."
                        " sequence will be used if no fna files are given.")
    parser.add_argument('--overlap', required=False, action='store_false',
                        default=True, help="Use to not remove genes overlapping with RNA features.")
    parser.add_argument('--compress', required=False, action='store_true',
                        default=False, help="Use to compress the output files.")
    parser.add_argument('--output', required=False, type=str, default="synta_outputdir"+time.strftime(
        "_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory path (optionnal)")
    parser.add_argument("--locustag", required=False, type=str,
                        help="Locustag to use for feature ID. If not provided, it will use generic random string IDs of length 12.")
    parser.add_argument("--basename", required=False, type=str,
                        help="Basename to use for output files. If not provided, it will be guessed from the input file name.")
    parser.add_argument("--norna", required=False, action="store_true",
                        default=False, help="Use to avoid annotating RNA features.")
    parser.add_argument("--cpu", required=False, type=int,
                        default=1, help="Number of cpus to use.")
    parser.add_argument("--format", required=False, type=str.lower, default="gff",
                        help="Different formats that you want as output, separated by a ','. Accepted strings are:  faa fna gff ffn.")
    parser.add_argument("--verbose", required=False, action="store_true",
                        default=False, help="Use to see the DEBUG log, which outputs uppon")
    args = parser.parse_args()

    # if any of them is not none, it's good.
    if args.fna is None and args.gbff is None:
        raise Exception(
            "You must provide at least a fna (with --fna) or a gbff (with --gbff) file to annotate from.")
    if args.locustag is None:
        args.locustag = ''.join(choice(ascii_uppercase) for i in range(12))

    return args


def main():
    start = time.time()
    args = cmdLine()
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(stream=sys.stdout, level=level,
                        format='\n%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    version_numbers = check_versions()

    if args.gbff:
        gbffFile = read_compressed_or_not(args.gbff)

    if args.fna is None:
        logging.getLogger().info("Reading the dna sequences from the gbff file.")
        contigs = read_gbff(gbffFile)
        fastaFile = write_tmp_fasta(contigs)
    else:
        logging.getLogger().info("Reading the dna sequences from the fna file.")
        fastaFile = read_compressed_or_not(args.fna)
        contigs = read_fasta(fastaFile)
        if is_compressed(args.fna):
            fastaFile = write_tmp_fasta(contigs)

    genes = syntaxic_annotation(
        fastaFile, args.cpu, args.norna, args.locustag)
    if args.overlap and not args.norna:
        # sorting and removing CDS that overlap.
        genes = overlap_filter(genes, contigs)
    else:
        # need to sort the genes by contig size, then by start postion.
        genes = sorted(
            genes, key=lambda x:  (-len(contigs[x.contig]), x.start))

    logging.getLogger().info("Writting output files...")
    if args.basename:
        outbasename = os.path.abspath(
            args.output + "/" + args.basename)
    elif args.fna:
        outbasename = mk_basename(args.output, args.fna)
    else:
        outbasename = mk_basename(args.output, args.gbff)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    write_output(outbasename, contigs, genes, args.compress,
                 args.format, args.cpu, version_numbers)

    logging.getLogger().info(
        f"There is a total of {len(genes)} annotated features.")
    logging.getLogger().info(f"Took {round(time.time() - start,2)} seconds.")
    # from here, all processes are done.


if __name__ == "__main__":
    main()
