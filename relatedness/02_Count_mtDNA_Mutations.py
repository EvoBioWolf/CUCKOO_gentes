import allel
import itertools
import numpy as np
import gzip
import argparse

# Function to read a VCF and return a GenotypeArray
def read_vcf(filename):
    callset = allel.read_vcf(filename)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    return genotypes, callset['samples']

# Function to compute pairwise comparisons
def pairwise_comparisons(genotypes, samples):
    results = []
    ploidy = genotypes.ploidy  # Get ploidy
    
    for id1, id2 in itertools.combinations(range(len(samples)), 2):
        g1 = genotypes[:, id1]
        g2 = genotypes[:, id2]

        # Filter out missing data
        non_missing = ~(g1.is_missing() | g2.is_missing())
        g1 = g1.compress(non_missing, axis=0)
        g2 = g2.compress(non_missing, axis=0)

        # Count the number of valid sites
        number_sites = non_missing.sum()

        # Count the number of SNPs correctly accounting for ploidy
        number_snps = np.count_nonzero(g1.to_n_alt() != g2.to_n_alt())
        
        results.append((samples[id1], samples[id2], number_sites, number_snps))
    return results

if __name__ == '__main__':
    # Argument parsing
    parser = argparse.ArgumentParser(description='Compute pairwise SNP differences from a VCF file.')
    parser.add_argument('--vcf', required=True, help='Input VCF file (gzipped)')
    parser.add_argument('--out', required=True, help='Output file for pairwise comparisons')
    args = parser.parse_args()

    # Read VCF
    genotypes, samples = read_vcf(args.vcf)

    # Perform pairwise comparisons
    comparisons = pairwise_comparisons(genotypes, samples)

    # Output results
    with open(args.out, 'w') as f:
        for comparison in comparisons:
            f.write(f'{comparison[0]} {comparison[1]} {comparison[2]} {comparison[3]}\n')

    print(f'Pairwise comparisons saved to {args.out}')

