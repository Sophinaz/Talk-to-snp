import requests


def fetch_snps_from_category(base_url, category):
    """
    Fetch SNP names from a specific category using the MediaWiki API with pagination.
    Attempts to fetch up to 'limit' entries per request.
    """
    params = {
        'action': 'query',
        'list': 'categorymembers',
        'cmtitle': f'Category:{category}',
        'cmlimit': 500,  # Try for 5000 entries, but may be limited by server to 500
        'format': 'json'
    }
    snps = []

    while True:
        response = requests.get(base_url, params=params)
        if response.ok:
            data = response.json()
            snps.extend([item['title'] for item in data['query']['categorymembers']])

            # Check if there is more data to load
            if 'continue' in data:
                params['cmcontinue'] = data['continue']['cmcontinue']
            else:
                break
        else:
            print("Failed to fetch data")
            response.raise_for_status()

    return snps


if __name__ == '__main__':
    # URL to the SNPedia API endpoint
    base_url = "http://bots.snpedia.com/api.php"

    # Category of interest
    category = "Is_a_snp"

    # Fetch the SNP names using pagination with an attempt for a higher limit
    snp_names = fetch_snps_from_category(base_url, category)
    print("Retrieved SNP Names:", len(snp_names))
    print(snp_names[:50])  # Print first 50 names for brevity

    with open('data/snpedia_variants_list.txt', 'w') as f:
        f.write('\n'.join(snp_names))