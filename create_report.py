import argparse
import os
import openai
import pandas as pd

openai.api_key = os.getenv('OPENAI_API_KEY')

report_instruction_prompt = """You are provided with a detailed list of single nucleotide polymorphisms (SNPs) from an individual's genetic profile. Your role is to function as a distinguished lifestyle coach, earning an annual salary of $300,000 and holding a PhD in genetics. Deliver thorough, precise, and personalized lifestyle recommendations encompassing various domains such as nutrition, exercise, sleep hygiene, and other health-related practices. Ensure that each recommendation is deeply rooted in scientific evidence and meticulously tailored to align with the individual's unique genetic makeup."""

fetched_allele_df = pd.read_csv('data/snpedia_alleles_filtered.csv', index_col=0)
fetched_allele_df = fetched_allele_df[~fetched_allele_df.summary.isna()]


def generate_report(formatted_snps):
    client = openai.OpenAI()
    response = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": report_instruction_prompt},
            {"role": "user", "content": formatted_snps}
        ]
    )

    return response.choices[0].message.content


allele_normalization = {
    'GG': '(g;g)',
    'CC': '(c;c)',
    'TT': '(t;t)',
    'AA': '(a;a)',
    'AG': '(a;g)',
    'GA': '(a;g)',
    'CT': '(c;t)',
    'TC': '(c;t)',
    'GT': '(g;t)',
    'TG': '(g;t)',
    'CG': '(c;g)',
    'GC': '(c;g)',
    'AT': '(a;t)',
    'TA': '(t;a)',
    'AC': '(a;c)',
    'CA': '(a;c)',
}


def format_variant_information(found_variants_df, magnitude_threshold=2.):
    found_alleles_filtered = found_variants_df.query(f'magnitude>={magnitude_threshold}')
    retrieved_descriptions_formatted = '\n'.join((found_alleles_filtered.index + found_alleles_filtered.magnitude.apply(
        lambda x: f' Magnitude: {x} ') + 'Summary: ' + found_alleles_filtered.summary).tolist())
    return retrieved_descriptions_formatted


def find_variants(user_snp_df):
    formatted_alleles = (user_snp_df['rsid'] + user_snp_df.genotype.replace(allele_normalization)).values.tolist()
    found_alleles = list(set(fetched_allele_df.index).intersection(formatted_alleles))
    found_alleles_df = fetched_allele_df.loc[found_alleles].sort_values('magnitude', ascending=False)
    return found_alleles_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process SNP data.')
    parser.add_argument('snp_df_path', type=str, help='Path to user SNP data file')
    parser.add_argument('-t', '--threshold', type=float, default=2.0, help='Magnitude threshold')
    parser.add_argument('-o', '--output', type=str, default='report.txt', help='Path to output file')

    args = parser.parse_args()

    user_snp_df = pd.read_csv(args.snp_df_path, sep='\t', comment="#")
    user_snp_df.columns = ['rsid', 'chromosome', 'position', 'genotype']

    found_user_alleles = find_variants(user_snp_df)
    snp_descriptions_formatted = format_variant_information(found_user_alleles, magnitude_threshold=args.threshold)
    report = generate_report(snp_descriptions_formatted)

    with open(args.output, 'w') as f:
        f.write(report)
