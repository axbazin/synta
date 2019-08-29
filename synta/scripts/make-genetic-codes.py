#!/usr/bin/env python3
#coding: utf-8

from itertools import product

class genetic_code:
    def __init__(self):
        self.names = set()
        self.ID = None
        self.base_table = []
        self.start = []
        self.bases = [[],[],[]]

    def print(self):
        """
           Prints the basic original self.
        """
        print("name",self.names)
        print("ID",self.ID)
        print("AA table",self.base_table)
        print("start AA table",self.start)
        print("NT bases", self.bases)

    def assign_codons_start(self):
        """
        Assigns the proper amino acid to the start codon table.
        """
        for i, base in enumerate(self.base_table):
            if self.start[i] == "-":
                self.start[i] = base

    def make_trans_table(self):
        """
            Create the basic translation tables.
        """
        self.trans_table = {}
        self.start_table = {}
        for i in range(len(self.base_table)):
            self.trans_table[self.bases[0][i] + self.bases[1][i] + self.bases[2][i]] = self.base_table[i]
            self.start_table[self.bases[0][i] + self.bases[1][i] + self.bases[2][i]] = self.start[i]

    def compute(self):
        """
            Computes the translation tables of this genetic code.
        """
        self.assign_codons_start()
        self.make_trans_table()
        self.translation_table = extend_trans_tables(self.trans_table)
        self.start_translation_table = extend_trans_tables(self.start_table)

def extend_trans_tables(table):
        """
            Greedily testing every combination to see if all the string they mean give the same amino acid.
        """
        new_trans_table = {}
        iupac_extend = {"A":{"A"},"C":{"C"},"T":{"T"},"G":{"G"},"R":{"A","G"}, "Y":{"C","T"}, "S":{"G","C"}, "W":{"A","T"}, "K":{"G","T"}, "M":{"A","C"}, "B":{"C","G","T"}, "D":{"A","G","T"},"H":{"A","C","T"},"V":{"A","C","G"},"N":{"A","C","T","G"}}
        IUPAC = set(["A","T","G","C","N","R","Y","S","W","K","M","B","V","D","H"])
        #test if a combination give the same amino acid.
        for comb in product(IUPAC, repeat=3):
            currStrings = iupac_extend[comb[0]]#initialize with the corresponding nucleic letter meaning of the first letter.
            for letter in comb[1:]:
                newCurrString = set()
                values = iupac_extend[letter]
                for val in values:
                    for cs in currStrings:
                        newCurrString.add(cs + val)
                currStrings = newCurrString
            resulting_aa = set()
            for string in currStrings:
                resulting_aa.add(table[string])
            if len(resulting_aa) == 1:
                new_trans_table["".join(comb)] = resulting_aa.pop()#there is only one.

        return new_trans_table

def read_prt_file(filename):
    genetic_codes = {}
    currGenCode = genetic_code()
    with open(filename,"r") as f:
        for line in f:
            if line.startswith("--"):
                pass
            elif line.startswith(" {"):#starting something new
                if currGenCode.ID is not None:
                    currGenCode.compute()
                    genetic_codes[currGenCode.ID] = {"trans_table":currGenCode.translation_table, "start_table":currGenCode.start_translation_table}
                currGenCode = genetic_code()
            elif line.startswith("  name"):
                currGenCode.names.add(line[6:].replace('"','').replace(",","").strip())
            elif line.startswith("  id"):
                currGenCode.ID = line.strip().split()[1]
            elif line.startswith("  ncbieaa"):
                currGenCode.base_table = list(line.strip().split()[1].replace(",","").replace('"',""))
            elif line.startswith("  sncbieaa"):
                currGenCode.start = list(line.strip().split()[1].replace('"',""))
            elif line.startswith("  -- Base"):
                currGenCode.bases[int(line[9])-1] = list( line[10:].strip() )
    return genetic_codes

def main():
    filename = "gc.prt"
    print(read_prt_file(filename))


if __name__ == "__main__":
    main()