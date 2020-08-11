from itertools import product, combinations, groupby, permutations
from collections import Counter
import time


translation_table = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                     'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                     'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                     'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                     'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                     'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                     'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                     'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                     'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                     'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                     'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                     'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                     'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                     'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                     'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                     'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
                     }

# nomenclature for degenerate codons
expanded_code = {'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
                 'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'], 'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T'],
                 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
                 'N': ['A', 'C', 'G', 'T']
                 }
expanded_code_values = list(expanded_code.values())

# helpful for validating input
valid_nucleotides = {'A','C','G','T','W','S','M','K','R','Y','B','D','H','V','N'}
valid_aa = {'G','A','V','L','I','M','F','W','P','S','T','C','Y','N','Q','D','E','K','R','H','*'}

base_filter_cache = {}

def get_codon_for_amino_acids(amino_acids):
    # Validate that input is a set of amino acids
    validate_amino_acid_input(amino_acids)
    # Get list of possible codons
    possible_codons = get_possible_codons(amino_acids)
    # Get amino acid count for each possible codons
    codon_aa_count = get_amino_acid_count(possible_codons)
    # Get smallest number of amino acids
    smallest_aa_count = codon_aa_count[min(codon_aa_count, key=codon_aa_count.get)]
    # Get codons that have the smallest number of amino acids
    efficient_codons = set([key for (key, value) in codon_aa_count.items() if value == smallest_aa_count])
    # Calculate efficiency
    efficiency = len(amino_acids)/smallest_aa_count
    return(efficient_codons, efficiency)

def get_amino_acid_count(codons):
    codon_aa_count = {}
    # Iterate through codon list
    for codon in codons:
        # Get the bases the expanded code could code through
        base1 = expanded_code[codon[0]]
        base2 = expanded_code[codon[1]]
        base3 = expanded_code[codon[2]]
        # Itertools pairwise product forming
        possible_amino_acids = list(product(*[base1, base2, base3]))
        # Flatten the list of lists
        possible_amino_acids_flat = []
        for amino_acid in possible_amino_acids:
            possible_amino_acids_flat.append(''.join(amino_acid))
        # Translate base combinations into possible amino acids
        amino_acids = []
        for possible_amino_acids in possible_amino_acids_flat:
            amino_acids.append(translation_table[possible_amino_acids])
        # Remove duplicates
        amino_acids_dedup = list(set(amino_acids))
        # Get count of possible amino acids
        codon_aa_count[''.join(codon)] = len(amino_acids_dedup)
    return(codon_aa_count)

def validate_amino_acid_input(amino_acids):
    if type(amino_acids) is not set:
        print('Input must be a set, not a %s'%type(amino_acids))
        quit()
    if not amino_acids.issubset(valid_aa):
        print('%s is not a valid set of amino acids'%amino_acids)
        quit()
    if len(amino_acids) == 0:
        print('Input length must be greater than 0.')
        quit()

def get_possible_codons(amino_acids):
    codon1 = []
    codon2 = []
    codon3 = []
    # Get bases for each position for each amino acid
    for amino_acid in amino_acids:
        possible_codons = [key  for (key, value) in translation_table.items() if value == amino_acid]
        codon1.append(sorted(list(set([codon[0] for codon in possible_codons]))))
        codon2.append(sorted(list(set([codon[1] for codon in possible_codons]))))
        codon3.append(sorted(list(set([codon[2] for codon in possible_codons]))))
    # Remove duplicates
    codon1 = dedup_listoflists(codon1)
    codon2 = dedup_listoflists(codon2)
    codon3 = dedup_listoflists(codon3)
    # Get the bases for each position
    codon1_bases = get_bases(codon1)
    codon2_bases = get_bases(codon2)
    codon3_bases = get_bases(codon3)
    # Get all possible codons given the bases at each position
    all_possible_codons = get_all_possible_codons(codon1_bases, codon2_bases, codon3_bases)
    return(all_possible_codons)

def get_all_possible_codons(codon1_bases, codon2_bases, codon3_bases):
    # Get all possible codons
    all_possible_codons = list(product(*[codon1_bases, codon2_bases, codon3_bases]))
    # Flatten the list of lists
    all_possible_codons_flat = []
    for codon in all_possible_codons:
        all_possible_codons_flat.append(''.join(codon))
    # Remove duplicates
    all_possible_codons_dedup = list(set(all_possible_codons_flat))
    return(all_possible_codons)

def get_bases(codon_bases):
    # Calculate possible 3 base codon combinations
    intersects = list(product(*codon_bases))
    # Sort within each possible combination to read accurately from expanded_code dictionary
    intersects_sorted = list(map(lambda x: sorted(x), intersects)); intersects_sorted
    # Remove duplicates
    intersects_dedupe = dedup_listoflists(intersects_sorted)
    # Get expanded bases that include the required ones
    bases = base_filter(intersects_dedupe)
    # Get bases from expanded_bases
    expanded_bases = get_expanded_bases(bases)
    return(expanded_bases)

def get_expanded_bases(bases):
    expanded_codons = []
    # Get expanded_codon for each base
    for base in bases:
        expanded_codon = [key  for (key, value) in expanded_code.items() if value == base]
        expanded_codons.append(expanded_codon)
    # Flatten list of lists
    expanded_codons_flat = [codon for sublist in expanded_codons for codon in sublist]
    # Remove duplicates
    expanded_codons_dedup = list(set(expanded_codons_flat))
    return(expanded_codons_dedup)

def base_filter(intersects):
    bases = []
    for intersect in intersects:
        # Caching if available
        if tuple(intersect) in base_filter_cache:
            bases.extend(base_filter_cache[tuple(intersect)])
        else:
            intersect_list = []
            # Iterate through values from expanded_code dictionary
            for expanded_code_value in expanded_code_values:
                # Check if expanded_code value contains all bases from the intersect
                if all(item in expanded_code_value for item in intersect):
                    intersect_list.append(expanded_code_value)
            # Add to overall codon set
            bases.extend(intersect_list)
            # Cache for future reuse
            base_filter_cache[tuple(intersect)] = intersect_list
    return(bases)

def dedup_listoflists(listoflists):
    # Remove duplicates
    deduped = list(listoflists for listoflists,_ in groupby(sorted(listoflists)))
    return(deduped)

def truncate_list_of_amino_acids(amino_acids):
    # Validate that input is a set of amino acids
    validate_amino_acid_input(amino_acids)
    # Feed in amino acids to see the start line
    (efficient_codons, efficiency) = get_codon_for_amino_acids(amino_acids)
    # Check if 100%, and thus complete
    if efficiency == 1.0:
        return(efficient_codons)
    else :
        # Start testing subsets of the amino acids`
        truncated_amino_acids = test_truncated_list(amino_acids)
        return(truncated_amino_acids)

def test_truncated_list(amino_acids):
    # Validate if complete
    found_100 = 0
    # Set of qualifying (100%) amino acids
    amino_acid_100 = set()
    # Count down set size
    for size in range(len(amino_acids)-1, 0, -1):
        # Get possible permutations of amino acids given reduced set size
        amino_acid_permutations = list(permutations(amino_acids, size))
        # Remove duplicates
        amino_acid_permutations = list(set(frozenset(item) for item in amino_acid_permutations))
        # Iterate through permutations for this size
        for perm in amino_acid_permutations:
            # Run main program and proceed based on the results
            (efficient_codons, efficiency) = get_codon_for_amino_acids(set(perm))
            if efficiency == 1.0:
                found_100 = 1
                amino_acid_100.add(perm)
        if found_100 == 1:
            return(amino_acid_100)
