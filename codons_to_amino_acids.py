import random

# A dictionary mapping each amino acid to the codons that represent it.
AMINO_ACID_TO_BPS = {
    "alanine": ["CGA",  "CGG",  "CGT",  "CGC"],
    "arginine": ["GCA",  "GCG",  "GCT",  "GCC", "TCT",  "TCC"],
    "asparagine": ["TTA",  "TTG"],
    "aspartate": ["CTA",  "CTG"],
    "cysteine": ["ACA",  "ACG"],
    "glutamate": ["CTT",  "CTC"],
    "glutamine": ["GTT",  "GTC"],
    "glycine": ["CCA",  "CCG",  "CCT",  "CCC"],
    "histidine": ["GTA",  "GTG"],
    "isoleucine": ["TAA",  "TAG",  "TAT"],
    "leucine": ["AAT",  "AAC",  "GAA",  "GAG",  "GAT",  "GAC"],
    "lysine": ["TTT",  "TTC"],
    "methionine": ["TAC"],
    "phenylalanine": ["AAA",  "AAG"],
    "proline": ["GGA",  "GGG",  "GGT",  "GGC"],
    "serine": ["AGA",  "AGG",  "AGT",  "AGC",  "TCA",  "TCG"],
    "threonine": ["TGA",  "TGG",  "TGT",  "TGC"],
    "tryptophan": ["ACC"],
    "tyrosine": ["ATA",  "ATG"],
    "valine": ["CAA",  "CAG",  "CAT",  "CAC"],
    "stop": ["ATT",  "ATC",  "ACT"],
    "start": ["ATG"]
}

# Using the hand coded dictonary above, create the reverse lookup table so we can
# look up an amino acid from a 3 letter base pair sequence (a codon).
BP_TO_AMINO_ACID = {}
for codon, bps in AMINO_ACID_TO_BPS.items():
    for bp in bps:
        BP_TO_AMINO_ACID[bp] = codon

# It's very common to abbreviate the amino acid sequence, using single letters for each
# This dictionary allows us to go back and forth between the single letter and full name
AMINO_ACID_SINGLE_LETTER_MAP = {
    'A': 'alanine',
    'R': 'arginine',
    'N': 'asparagine',
    'D': 'aspartate',
    'C': 'cysteine',
    'E': 'glutamate',
    'Q': 'glutamine',
    'G': 'glycine',
    'H': 'histidine',
    'I': 'isoleucine',
    'L': 'leucine',
    'K': 'lysine',
    'M': 'methionine',
    'F': 'phenylalanine',
    'P': 'proline',
    'S': 'serine',
    'T': 'threonine',
    'W': 'tryptophan',
    'Y': 'tyrosine',
    'V': 'valine'
}


def main():
    '''
    This function demonstrates the use of the functions below. We create a random
    base pair sequence with 15 codons (45 base pairs), convert that sequence to
    the sequence of amino acids it represents (full names), then we compute the
    number of ways to generate the Trp Cage protein, and generate 2 sequences of
    base pairs that both represent Trp Cage.
    '''
    seq = generate_bp_codon_sequence(15)
    print(seq)

    amino_acid_seq = convert_to_amino_acids(seq)
    print(amino_acid_seq)

    # trp cage is the shortest known protein.
    trp_cage_amino_acid = 'NLYIQWLKDGGPSSGRPPPS'
    print(compute_number_of_encodings(trp_cage_amino_acid))
    print(random_valid_dna_sequence(trp_cage_amino_acid))
    print(random_valid_dna_sequence(trp_cage_amino_acid))


def generate_bp_codon_sequence(codon_length):
    '''
    Generate a string of base pairs representing a random sequence of codons
    starting with start codon and ending with a stop codon. The sequence will
    have codon_length number of codons, meaning 3*codon_length base pairs.
    '''
    # First codon is the start codon
    bps = [ AMINO_ACID_TO_BPS['start'][0] ]

    # We want a list of amino acids without start and stop to use in the loop
    amino_acids = [
        amino_acid for amino_acid in AMINO_ACID_TO_BPS.keys()
            if amino_acid not in ['stop', 'start']
    ]


    # Randomly choose the codon_length -2 amino_acids (-2 for the start and stop codon)
    for _ in range(codon_length - 2):
        codon = random.choice(amino_acids)
        bp = random.choice(AMINO_ACID_TO_BPS[codon])
        bps.append(bp)

    # Last codon is the stop codon.
    bps.append(random.choice(AMINO_ACID_TO_BPS['stop']))

    # Use join to make the list a single string
    return ''.join(bps)


def convert_to_amino_acids(bp_string):
    '''
    Assuming the input bp_string is only made up of codons, the first of which
    starts at the 0th character, return a sequence of amino. bp_string must have
    a length that is a multiple of 3.
    '''
    codons = []
    for n in range(0, len(bp_string), 3):
        codon_bps = bp_string[n : n + 3]
        codon = BP_TO_AMINO_ACID[codon_bps]
        codons.append(codon)

    return codons


def random_valid_dna_sequence(single_letter_amino_acids):
    '''
    Given a single letter amino acid sequence, return a valid DNA sequences that could
    make this amino acid sequence. When multiple codons could be used, pick one randomly
    '''
    dna_sequence = AMINO_ACID_TO_BPS['start'][0] # Start with the start codon.

    # Amino acid by amino acid, randomly select one of the matching codons
    for amino_acid_letter in single_letter_amino_acids:
        full_amino_acid = AMINO_ACID_SINGLE_LETTER_MAP[amino_acid_letter]
        codon = random.choice(AMINO_ACID_TO_BPS[full_amino_acid])
        dna_sequence += codon

    return dna_sequence + random.choice(AMINO_ACID_TO_BPS['stop'])


def compute_number_of_encodings(single_letter_amino_acids):
    '''
    Given a single letter amino acid code, compute how many different DNA sections
    could be used to represent the protein (not including start/stop codons)
    '''
    total = 1

    # For each amino acid, determine how many codons there are, and repeatedly multiply
    # total by that amount.
    for amino_acid_letter in single_letter_amino_acids:
        amino_acid_name = AMINO_ACID_SINGLE_LETTER_MAP[amino_acid_letter]
        total *= len(AMINO_ACID_TO_BPS[amino_acid_name])

    return total


def full_amino_acid_names(single_letter_amino_acids):
    '''
    Given a single letter amino acid sequence, return a list with the full names
    of the amino acids.
    '''
    return [ AMINO_ACID_SINGLE_LETTER_MAP[amino_acid_letter]
        for  amino_acid_letter in single_letter_amino_acids
    ]


if __name__ == '__main__':
    main()
