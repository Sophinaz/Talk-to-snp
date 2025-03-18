
import pandas as pd
import json
import requests
import time
import re

from tqdm.auto import tqdm

SLEEP_TIME_SECONDS = 0.8
FETCH_BATCH_SIZE = 20
FETCH_ALLELES_BATCH_SIZE = 50


def fetch_snpedia_data(snps):
    base_url = "https://bots.snpedia.com/api.php"
    params = {
        'action': 'query',
        'titles': '|'.join(snps),
        'prop': 'revisions',
        'rvprop': 'content',
        'rvslots': 'main',  # Request the content of the main slot
        'format': 'json'
    }
    response = requests.get(base_url, params=params)
    if response.ok:
        data = response.json()
        # Extract content from the main slot of the revisions
        pages = data.get('query', {}).get('pages', {})
        for page_id, page in pages.items():
            if 'revisions' in page and page['revisions']:
                page_content = page['revisions'][0]['slots']['main']['*']
        return data
    else:
        return response.raise_for_status()

def fetch_data_in_batches(input_data, fetch_batch_size):
    fetched_data = []
    for i in tqdm(range(0, len(input_data), fetch_batch_size)):
        try:
            data = fetch_snpedia_data(input_data[i:i+fetch_batch_size])
            fetched_data.append(data)
            time.sleep(SLEEP_TIME_SECONDS)
        except KeyboardInterrupt:
            break
        except:
            print(f"Failed to fetch for batch {i}")
    return fetched_data

def structure_fetched_data(fetched_data):
    parsed_data = {}
    empty_count = 0
    for batch in fetched_data:
        for key in batch['query']['pages'].keys():
            page_content = batch['query']['pages'][key]
            if 'revisions' in page_content.keys():
                snp_text = page_content['revisions'][-1]['slots']['main']['*']
                snp_title = batch['query']['pages'][key]['title'].lower()
                parsed_data[snp_title] = snp_text
            else:
                empty_count += 1
    print(f'{empty_count} empty variants occured')
    result_df = pd.DataFrame.from_dict(parsed_data, orient='index')
    return result_df

def parse_wikitext_data(data_string):
    data_string = data_string.strip('{}')  # Remove the curly braces
    data_lines = data_string.split('\n')  # Split at newlines
    result = {}

    for line in data_lines:
        key_value = line.split('=')
        if len(key_value) == 2:
            key = key_value[0].strip('|')
            value = key_value[1]
            result[key] = value

    return pd.Series(result)


common_summary_blacklist = ['common in clinvar',
                            'common in complete genomics',
                            'common/normal',
                            'common genotype',
                            'normal',
                            'common on affy axiom data',
                            'common',
                            'Normal',
                            '',
                            'normal risk'
                            ]

if __name__ == '__main__':
    with open('data/snpedia_variants_list.txt', 'r') as f:
        snpedia_variants_list = [x.strip() for x in f.readlines()]

    fetched_data = fetch_data_in_batches(snpedia_variants_list[:6000], FETCH_BATCH_SIZE)
    snpedia_df = structure_fetched_data(fetched_data)

    snpedia_df.columns = ['snp_text']
    snpedia_df['variants'] = snpedia_df.apply({'snp_text':lambda x: list(set(re.findall('\([ACGT];[ACGT]\)',x)))})
    snpedia_df.to_csv('data/snpedia_fetched.csv')

    snpedia_variants_list_with_alleles = []
    for rsid, item in snpedia_df[snpedia_df.variants.apply(len)>0].iterrows():
        snpedia_variants_list_with_alleles.extend([rsid+variant for variant in item.variants])

    fetched_data_with_alleles = fetch_data_in_batches(snpedia_variants_list_with_alleles,
                                                      FETCH_ALLELES_BATCH_SIZE)
    snpedia_allele_df = structure_fetched_data(fetched_data_with_alleles)
    snpedia_allele_df.columns = ['wikitext_data']
    snpedia_allele_df = snpedia_allele_df['wikitext_data'].apply(parse_wikitext_data)
    snpedia_allele_df.to_csv('data/snpedia_alleles.csv')

    snpedia_allele_df_filtered=snpedia_allele_df[~snpedia_allele_df.summary.isin(common_summary_blacklist)]
    snpedia_allele_df_filtered = snpedia_allele_df_filtered[~snpedia_allele_df_filtered.summary.isna()]
    snpedia_allele_df_filtered['magnitude'] = pd.to_numeric(snpedia_allele_df_filtered.magnitude).fillna(1)
    snpedia_allele_df_filtered.to_csv('data/snpedia_alleles_filtered.csv')

