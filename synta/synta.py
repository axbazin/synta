#!/usr/bin/env python3
#coding : utf-8
#PYTHON_ARGCOMPLETE_OK

# minimum python version to use is 3.6
import sys
assert sys.version_info >= (3, 6)

import pkg_resources
from string import ascii_uppercase
from random import choice
import ast
import tempfile
from subprocess import Popen, PIPE
from multiprocessing import Pool
import argparse
import time
import os
import logging
from collections import defaultdict

#flavor package, not required.
try:
    import argcomplete
except ImportError:
    pass

#own scripts
from synta.utilitaries.genetic_codes import genetic_codes
from synta.utilitaries.file_handlers import is_compressed, read_compressed_or_not, write_compress_or_not



class gene:
    def __init__(self, ID, contig, start, stop, strand, geneType, genetic_code = None, product=None, inference=None):
        self.ID = ID
        self.contig = contig
        self.start = int(start)
        self.stop = int(stop)
        self.type = geneType
        self.strand = strand
        self.inference = inference
        self.product = product
        self.genetic_code = genetic_code

        self.dbxref = None
        self.locus_tag = None
        self.protein_id = None

    def saveDBinfo(self, dbxref, locus, protein_id):
        self.dbxref = dbxref
        self.locus_tag = locus
        self.protein_id = protein_id

    def get_gff(self, version_numbers):
        """
            Takes the software's version numbers used for different types of genes (key:  gene type, value:  tuple with version string andsoftware name). Writes the gene's information in gff3 format
        """
        version, soft = version_numbers[self.type]

        Feature = f"{self.contig}\t{soft}:{version}\tgene\t{self.start}\t{self.stop}\t1\t{self.strand}\t0\tID=gene_{self.ID}\n"
        Feature += f'{self.contig}\t{soft}:{version}\t{self.type}\t{self.start}\t{self.stop}\t1\t{self.strand}\t0\tID={self.ID};Parent=gene_{self.ID};inference={self.inference}'

        if self.dbxref:
            Feature += f";db_ref={ ','.join(self.dbxref)}"
        if self.locus_tag is not None:
            Feature += f";locus_tag={self.locus_tag}"
        if self.protein_id is not None:
            Feature += f";protein_id={self.protein_id}"

        Feature += f";product={self.product}\n" if self.product else "\n"
        return Feature

    def get_faa(self, dna_seq, code):
        """Takes a dna sequence string and returns the gene object's ID with the string  translated in faa format"""
        faa_seq = translate(dna_seq, code)
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
    # see https://www.bioinformatics.org/sms/iupac.html for the code.
    # complement = {'A':  'T', 'C':  'G', 'G':  'C', 'T':  'A', 'N': 'N' } ## basic
    rcseq = ""
    for i in reversed(seq):
        rcseq += complement[i]
    return rcseq


def translate(seq, code):
    """ translates the given dna sequence with table code 11 of the ncbi (bacteria)"""
    # code:  https: //www.bioinformatics.org/sms/iupac.html
    start_table = code["start_table"]
    table = code["trans_table"]

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
                                 inference="COORDINATES:profile",
                                 product=lineData[1] + lineData[4]))
    logging.getLogger().info(
        f"Done with Aragorn(tRNA). There are {len(geneObjs)} tRNAs.")
    return geneObjs


def launch_prodigal(fnaFile, locustag, code):
    """
        launches Prodigal to annotate CDS. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects.
    """
    logging.getLogger().debug("Running Prodigal(CDS).")
    cmd = ["prodigal", "-f", "sco","-g",code, "-m", "-c", "-i", fnaFile, "-p", "single", "-q"]
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
                                 genetic_code=code,
                                 inference="ab initio prediction"))

    logging.getLogger().info(
        f"Done with Prodigal(CDS). There are {len(geneObjs)} CDSs.")
    return geneObjs


def launch_infernal(fnaFile, locustag, kingdom):
    """
        launches Infernal in hmmer-only mode to annotate rRNAs. Takes a fna file name and a locustag to give an ID to the found genes.
        returns the annotated genes in a list of gene objects
    """
    if kingdom == "bacteria":
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_bact.cm"
    elif kingdom == "archaea":
        modelfile = os.path.dirname(os.path.realpath(__file__)) + "/rRNA_DB/rRNA_arch.cm"

    logging.getLogger().debug("Running Infernal(rRNA).")
    tmpFile = tempfile.NamedTemporaryFile(mode="r", dir = "/dev/shm/")
    cmd = ["cmscan", "--tblout", tmpFile.name, "--hmmonly", "--cpu",
           str(1), "--noali", modelfile, fnaFile]
    logging.getLogger().debug(f"command for Infernal:  '{' '.join(cmd)}'")
    p = Popen(cmd, stdout=open(os.devnull, "w"), stderr=PIPE)
    err = p.communicate()[1].decode().split()
    if err != []:
        if err[0] == 'Error: ':
            raise Exception(
                f"Infernal (cmscan) failed with error:  '{ ' '.join(err) }'. If you never used this script, you should press the .cm file using cmpress executable from Infernal. You should find the file in '{os.path.dirname(os.path.realpath(__file__))}/rRNA_DB/'.")
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
                                 product=" ".join(lineData[17:]),
                                 inference="COORDINATES:profile"))

    logging.getLogger().info(
        f"Done with Infernal(rRNA). There are {len(geneObjs)} rRNAs.")
    return geneObjs


def write_output(outputbasename, contigs, genes, compress, formats, cpu, versions, code):
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
        for forms in formats.split(","):
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
                translation_table = genetic_codes(code)
                wfaa = p.apply_async(func=write_faa, args=(
                    outputbasename, contigs, genes, compress, translation_table))

        # there is a better way of writing this 'wait until all processes are done', isn't there ?
        for forms in formats.split(","):
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
    contig2genes = defaultdict(list)
    for gene in genes:
        contig2genes[gene.contig].append(gene)
    outfile = write_compress_or_not(output + ".gff", compress)
    outfile.write("##gff-version 3\n")
    for contig in sorted(contigs.keys(), key = lambda x : len(contigs[x])):#biggest contig first.
        outfile.write(
            f"##sequence-region {contig} 1 {len(contigs[contig])}\n")
        #The circularity could be presumed with a gbff input...
        outfile.write(
            f"{contig}\tsynta:{__version__}\tregion\t1\t{len(contigs[contig])}\t.\t+\t.\tID={contig}\n")
        for gene in contig2genes[contig]:
            outfile.write(gene.get_gff(versions))
    outfile.close()
    logging.getLogger().debug("Done writing GFF file.")


def write_ffn(output, contigs, genes, compress):
    """
        Writes a ffn formated file.
    """
    logging.getLogger().debug("Writting FFN file ...")
    outfile = write_compress_or_not(output + ".ffn", compress)
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


def write_faa(output, contigs, genes, compress, code):
    """
        Writes a faa formated file.
    """
    logging.getLogger().debug("Writting FAA file ...")
    outfile = write_compress_or_not(output + ".faa", compress)
    for gene in genes:
        if gene.type == "CDS":  # faa is only for coding sequences !
            if gene.strand == "+":
                outfile.write(gene.get_faa(
                    contigs[gene.contig][gene.start-1: gene.stop], code))
            elif gene.strand == "-":
                outfile.write(gene.get_faa(reverse_complement(
                    contigs[gene.contig][gene.start-1: gene.stop]), code))
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
    tmpFile = tempfile.NamedTemporaryFile(mode="w", dir = "/dev/shm/")
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


def syntaxic_annotation(fastaFile, cpu, norna, locustag, kingdom, code):
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
                                 args=(fastaFile.name, locustag, code))
        logging.getLogger().debug("Started the process launching Prodigal")
        if not norna:
            araGenes = p.apply_async(
                func=launch_aragorn, args=(fastaFile.name, locustag))
            logging.getLogger().debug("Started the process launching Aragorn")
            infGenes = p.apply_async(
                func=launch_infernal, args=(fastaFile.name, locustag, kingdom))
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
        Checks if the programs we need exist in the PATH, and their version number
        Returns a dictionnary with the type of annotated genes as key, and a tuple with the version string and the software name as value
    """
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
        Reads a file-like object containing a gbff file.
    """
    geneObjs = []
    contigs = {}
    frameshifts = 0
    pseudo = 0
    trans_except = 0
    incomplete = 0
    logging.getLogger().debug("Extracting genes informations from the given gbff")
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
        dbxref = set()
        protein_id = ""
        locus_tag = ""
        usefulInfo = False
        start = None
        end = None
        strand = None
        line = lines.pop()
        while not line.startswith("ORIGIN"):
            currType = line[5:21].strip()
            if currType != "":
                # anything needs a locus tag (in the format description)
                if usefulInfo:
                    newGene = gene(protein_id, contigID, start, end, strand, objType)
                    newGene.saveDBinfo(dbxref, locus_tag, protein_id)
                    geneObjs.append(newGene)
                usefulInfo = False
                objType = currType
                if objType in ['CDS']:  # only CDS for now
                    dbxref = set()
                    protein_id = ""
                    locus_tag = ""
                    try:
                        if not 'join' in line[21:]:
                            usefulInfo = True
                            if line[21:].startswith('complement('):
                                strand = "-"
                                start, end = line[32:].replace(
                                    ')', '').split("..")
                            else:
                                strand = "+"
                                start, end = line[21:].strip().split('..')
                            if '>' in start or '<' in start or '>' in end or '<' in end:
                                usefulInfo = False
                                incomplete += 1
                    except ValueError:
                        # print('Frameshift')
                        frameshifts += 1
                        # don't know what to do with that, ignoring for now.
                        # there is a protein with a frameshift mecanism.
            elif usefulInfo:  # current info goes to current objtype, if it's useful.
                if line[21:].startswith("/db_xref"):
                    dbxref.add(line.split("=")[1].replace('"', '').strip())
                elif line[21:].startswith('/locus_tag'):
                    locus_tag = line.split("=")[1].replace('"', '').strip()
                # remove pseudogenes.
                elif line[21:].startswith('/protein_id'):
                    protein_id = line.split('=')[1].replace('"', '').strip()
                # if it's a pseudogene, we're not keeping it.
                elif line[21:].startswith("/pseudo"):
                    pseudo += 1
                    usefulInfo = False
                # that's probably a codon stop into selenocystein.
                elif line[21:].startswith("/transl_except"):
                    trans_except += 1
                    usefulInfo = False
            line = lines.pop()
        # end of contig.
        line = lines.pop()  # first sequence line.
        # if the seq was to be gotten, it would be here.
        sequence = ""
        while not line.startswith('//'):
            sequence += line[10:].replace(" ", "").strip()
            line = lines.pop()
        # ended the contig/sequence.
        contigs[contigID] = sequence.upper()

    logging.getLogger().debug("Done extracting informations from the gbff")
    logging.getLogger().debug(
        f"There was {frameshifts} elements with multiple frames, {incomplete} incomplete genes, {pseudo} pseudogenes and {trans_except} translation exceptions.")
    return geneObjs, contigs


def compare_to_gbff(contigs, genes, gbffObjs):

    # sort the same way that 'genes'
    gbffObjs = sorted(
        gbffObjs, key=lambda x: (-len(contigs[x.contig]), x.start))

    gbffIndex = 0
    geneIndex = 0

    precContig = ""
    nomatchGbff = 0
    nomatchGenes = 0
    matchEnd = 0
    matchStart = 0

    perfect = 0
    while gbffIndex < len(gbffObjs) and geneIndex < len(genes):
        # print(gbffIndex, geneIndex)
        # print(gbffObjs[gbffIndex].contig, genes[geneIndex].contig )
        # print(gbffObjs[gbffIndex].strand, genes[geneIndex].strand)
        # print(gbffObjs[gbffIndex].start, genes[geneIndex].start)
        # print(gbffObjs[gbffIndex].stop, genes[geneIndex].stop)
        # print("#####")
        if gbffObjs[gbffIndex].contig == genes[geneIndex].contig and gbffObjs[gbffIndex].strand == genes[geneIndex].strand:  # then we can compare
            precContig = genes[geneIndex].contig
            if gbffObjs[gbffIndex].strand == "+":
                # start is different but stop is equal.
                if gbffObjs[gbffIndex].stop == genes[geneIndex].stop:
                    # if stop is equal, we consider the protein being the same !
                    genes[geneIndex].saveDBinfo(
                        dbxref=gbffObjs[gbffIndex].dbxref, locus=gbffObjs[gbffIndex].locus_tag, protein_id=gbffObjs[gbffIndex].protein_id)
                    matchEnd += 1
                    if gbffObjs[gbffIndex].start == genes[geneIndex].start:
                        matchStart += 1
                        perfect += 1
                    geneIndex += 1
                    gbffIndex += 1

                elif genes[geneIndex].start < gbffObjs[gbffIndex].start:
                    nomatchGenes += 1
                    geneIndex += 1

                elif gbffObjs[gbffIndex].start < genes[geneIndex].start:
                    nomatchGbff += 1
                    gbffIndex += 1
            else:  # 'start' is the actual stop.
                # start is different but stop is equal.
                if gbffObjs[gbffIndex].start == genes[geneIndex].start:
                    # if stop is equal, we consider the protein being the same !
                    genes[geneIndex].saveDBinfo(
                        dbxref=gbffObjs[gbffIndex].dbxref, locus=gbffObjs[gbffIndex].locus_tag, protein_id=gbffObjs[gbffIndex].protein_id)
                    matchEnd += 1

                    if gbffObjs[gbffIndex].stop == genes[geneIndex].stop:
                        matchStart += 1
                        perfect += 1
                    geneIndex += 1
                    gbffIndex += 1

                elif genes[geneIndex].stop < gbffObjs[gbffIndex].stop:
                    nomatchGenes += 1
                    geneIndex += 1

                elif gbffObjs[gbffIndex].stop < genes[geneIndex].stop:
                    nomatchGbff += 1
                    gbffIndex += 1
        else:
            if gbffObjs[gbffIndex].contig != genes[geneIndex].contig:
                if genes[geneIndex].contig != precContig:
                    nomatchGbff += 1
                    gbffIndex += 1
                elif gbffObjs[gbffIndex].contig != precContig:
                    nomatchGenes += 1
                    geneIndex += 1
            elif gbffObjs[gbffIndex].strand != genes[geneIndex].strand:
                if genes[geneIndex].start < gbffObjs[gbffIndex].start:
                    nomatchGenes += 1
                    geneIndex += 1

                elif gbffObjs[gbffIndex].start < genes[geneIndex].start:
                    nomatchGbff += 1
                    gbffIndex += 1
                else:  # they're equal on different strands ???
                    print("WHUT?")

    if gbffIndex < len(gbffObjs) - 1:
        nomatchGbff += len(gbffObjs)-1 - gbffIndex

    if geneIndex < len(genes) - 1:
        nomatchGenes += len(genes)-1 - geneIndex

    print(f"perfect : {perfect}, end match : {matchEnd}, start match : {matchStart}, noGeneMatch : {nomatchGenes}, no gbff match : {nomatchGbff}")
    print(f"total genes: {len(genes)}, total gbff: {len(gbffObjs)}")


def cmdLine():
    """
        Functions that defines the command line arguments.
    """
    parser = argparse.ArgumentParser(description = "Quickly annotates CDS, rRNA and tRNA genes in prokaryote genomes", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")

    required.add_argument('--fna',  required=False, type=str, help="fasta(.gz) file to annotate. Will be used in priority if given.")
    required.add_argument('--gbff', required=False, type=str, help=".gbff (or .gbk) file to annotate. The sequence will be used if no fna files are given.")

    expert = parser.add_argument_group(title = "Expert options")
    expert.add_argument('--overlap', required=False, action='store_false',default=True, help="Use to not remove genes overlapping with RNA features.")
    expert.add_argument("--norna", required=False, action="store_true", default=False, help="Use to avoid annotating RNA features.")
    expert.add_argument("--kingdom",required = False, type = str.lower, default = "bacteria", choices = ["bacteria","archaea"], help = "Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.")
    expert.add_argument("--compare", required=False, action="store_true", default=False, help="Use to link database references from the gbff to our annotations.")
    expert.add_argument("--translation_table",required=False, default="11", help = "Translation table to use for gene calling.")
    outlog = parser.add_argument_group(title = "Output and formating")
    outlog.add_argument('--compress', required=False, action='store_true',default=False, help="Compress the output files using gzip.")
    outlog.add_argument('--output', required=False, type=str, default="synta_outputdir"+time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S", time.localtime())+"_PID"+str(os.getpid()), help="Output directory path (optionnal)")
    outlog.add_argument("--format", required=False, type=str.lower, default="gff", help="Different formats that you want as output, separated by a ','. Accepted strings are:  faa fna gff ffn.")
    outlog.add_argument("--locustag", required=False, type=str, help="Locustag to use for feature ID. If not provided, it will use generic random string IDs of length 12.")
    outlog.add_argument("--basename", required=False, type=str, help="Basename to use for output files. If not provided, it will be guessed from the input file name.")
    outlog.add_argument("--log",required = False, action = "store_true", help = "Save log to file")
    misc = parser.add_argument_group(title = "Misc options")
    misc.add_argument("--cpu", required=False, type=int, default=1, help="Number of cpus to use.")
    misc.add_argument("--verbose", required=False, action="store_true", default=False, help="show the DEBUG log")
    misc.add_argument('--version', action='version', version='%(prog)s ' + pkg_resources.get_distribution("synta").version)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if "argcomplete" in sys.modules:
        argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # if any of them is not none, it's good.
    if args.fna is None and args.gbff is None:
        raise Exception(
            "You must provide at least a fna (with --fna) or a gbff (with --gbff) file to annotate from.")
    if args.locustag is None:
        args.locustag = ''.join(choice(ascii_uppercase) for i in range(12))
    if args.compare and args.gbff is None:
        raise Exception(
            "You asked to link database references from a given gbff but did not give any gbff.")
    if args.gbff is not None and not args.compare and args.fna is not None:
        logging.getLogger().warning("You provided a gbff but that was not useful since you provided also a fna and did not ask to realise a comparison. Ignoring the gbff file")
        args.gbff = None
    return args


def main():
    start = time.time()
    args = cmdLine()


    if not os.path.exists(args.output):
        os.makedirs(args.output)

    if args.basename:
        outbasename = os.path.abspath(args.output + "/" + args.basename)
    elif args.fna:
        outbasename = mk_basename(args.output, args.fna)
    else:
        outbasename = mk_basename(args.output, args.gbff)
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(stream=sys.stdout, level=level, format='%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    if args.log :
        fhandler = logging.FileHandler(filename = outbasename + ".log",mode = "w")
        fhandler.setFormatter(logging.Formatter(fmt = "%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s", datefmt='%Y-%m-%d %H:%M:%S'))
        fhandler.setLevel(level)
        logging.getLogger().addHandler(fhandler)
    logging.getLogger().info("version: synta " + __version__)

    version_numbers = check_versions()
    if args.gbff:
        gbffFile = read_compressed_or_not(args.gbff)

    if args.compare and args.fna is None:
        logging.getLogger().info(
            "Reading the dna sequences and the gene informations from the gbff file.")
        gbffObjs, contigs = read_gbff(gbffFile)
        fastaFile = write_tmp_fasta(contigs)
    elif args.compare:# args.fna is not none
        logging.getLogger().info("Reading the gene informations from the gbff file.")
        gbffObjs, _ = read_gbff(gbffFile)
    elif args.fna is None:# there will be no comparisons.
        logging.getLogger().info("Reading the dna sequences from the gbff file.")
        _, contigs = read_gbff(gbffFile)
        fastaFile = write_tmp_fasta(contigs)

    if args.fna:
        logging.getLogger().info("Reading the dna sequences from the fna file.")
        fastaFile = read_compressed_or_not(args.fna)
        contigs = read_fasta(fastaFile)
        if is_compressed(args.fna):
            fastaFile = write_tmp_fasta(contigs)

    genes = syntaxic_annotation(fastaFile, args.cpu, args.norna, args.locustag, args.kingdom, args.translation_table)
    if args.overlap and not args.norna:
        # sorting and removing CDS that overlap.
        genes = overlap_filter(genes, contigs)
    else:
        # need to sort the genes by contig size, then by start postion.
        genes = sorted(
            genes, key=lambda x:  (-len(contigs[x.contig]), x.start))

    if args.compare:
        logging.getLogger().info(
            "Adding databases informations to our annotations from the gbff file...")
        compare_to_gbff(contigs, genes, gbffObjs)

    logging.getLogger().info("Writting output files...")

    write_output(outbasename, contigs, genes, args.compress,
                 args.format, args.cpu, version_numbers, args.translation_table)

    logging.getLogger().info(
        f"There is a total of {len(genes)} annotated features.")
    logging.getLogger().info(f"Took {round(time.time() - start,2)} seconds.")
    # from here, all processes are done.


if __name__ == "__main__":
    main()
