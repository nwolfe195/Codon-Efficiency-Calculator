# Codon Efficiency Calculator
 Given a set of amino acids, calculate the most efficient codons to encode those amino acids.

 truncate_list_of_amino_acids(amino_acids)
 def get_codon_for_amino_acids(amino_acids):

 This code contains two functions. The first, get_codon_for_amino_acids(amino_acids), takes a set of amino acids
 as an argument. Using an expanded set of bases beyond the usual ACGT representing bases that can be 2 or more
 bases, the script determines which codons utilizing the expanded bases could code for the provided argument set
 of amino acids, and selects the most efficient ones, ie the ones that code for the fewest amino acids beyond
 those given. It returns a set of the codons, as well as the efficiency rate they achieve.

 The second function, truncate_list_of_amino_acids(amino_acids), similarly takes a set of amino acids as an argument.
 This function uses the first function to test if the amino acids can be coded with 100% efficiency. If not, it
 iterates through decreasingly sized subsets of the original set, until such a set can be coded with 100% efficiency.

Usage
import codon_efficiency


amino_acids = {'A', 'I', 'V'}
codons, efficiency = codon_efficiency.get_codon_for_amino_acids(amino_acids)
print('Amino Acids: %s'%amino_acids)
print('Codons: %s'%codons)
print('Efficiency: %s'%efficiency)
truncated_amino_acids = codon_efficiency.truncate_list_of_amino_acids(amino_acids)
print('100% efficient truncated amino acids:')
print(truncated_amino_acids)
