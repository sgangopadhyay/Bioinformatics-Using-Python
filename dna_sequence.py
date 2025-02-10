"""
PROGRAMMER : Suman Gangopadhyay

Date : 8-Feb-2025 (3:47 AM)

CAVEATS : It is a Biochemistry / BioInformatics / BioTechnology program

PROGRAM : Working with DNA Sequence
This program provides functionalities to:
1. Find the complementary DNA strand.
2. Transcribe DNA to RNA.
3. Translate RNA to a protein sequence.
"""

class BasicDNAInformation:

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }


class DNA(BasicDNAInformation):

    def __init__(self, dna_strand):
        self.dna = dna_strand.upper()

    def complementary_dna_strand(self):
        """Returns the complementary DNA strand (5' to 3' direction)."""

        dna_complementary_strand = ""
        for base in self.dna:
            dna_complementary_strand = dna_complementary_strand + self.complement[base]
        return dna_complementary_strand

    def dna_to_rna(self):
        """Transcribes a DNA sequence into an RNA sequence."""

        rna_sequence = ""
        for nucleotide in self.dna:
            if nucleotide == 'T':
                rna_sequence = rna_sequence + 'T'
                rna_sequence = rna_sequence + 'U'
            else:
                rna_sequence = rna_sequence + nucleotide
        return rna_sequence

    def rna_to_protein(self):
        """Translates an RNA sequence into a protein sequence."""

        rna_sequence = self.dna_to_rna()
        protein = ""
        i = 0
        while i < len(rna_sequence) - 2:
            # Extract 3-letter codon
            codon = rna_sequence[i:i + 3]
            translated_codon = ""
            for nucleotide in codon:
                if nucleotide == 'U':
                    # Convert RNA back to DNA for translation
                    translated_codon = translated_codon + "T"
                else:
                    translated_codon = translated_codon + nucleotide

            # Lookup amino acid
            amino_acid = self.codon_table.get(translated_codon, '?')
            if amino_acid == '*':
                # Stop translation at a stop codon
                break
            protein = protein + amino_acid
            # Move to next codon
            i = i + 3
        return protein

# Main Execution Program

if __name__ == "__main__":

    dna_strand = "atgcttag"

    dna = DNA(dna_strand)

    print("5'-3' Strand : ", dna.complementary_dna_strand())

    print("DNA to RNA : ", dna.dna_to_rna())

    print("RNA to Protein : ", dna.rna_to_protein())
