"""
PubMed API integration for automated literature mining.

This module provides automated literature search and analysis
for aging and CRISPR research using the PubMed/NCBI APIs.
"""

import json
import logging
import re
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import requests

logger = logging.getLogger(__name__)


class PubMedMiner:
    """Automated literature mining from PubMed for aging research."""

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def __init__(self, api_key: Optional[str] = None, email: Optional[str] = None):
        """Initialize PubMed miner."""
        self.api_key = api_key
        self.email = email or "crispr-toolkit@example.com"
        self.rate_limit = 0.34 if api_key else 3.0  # Seconds between requests
        self.last_request = 0

        # Cache directory
        self.cache_dir = Path("data/literature")
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _rate_limit_request(self):
        """Enforce rate limiting for API requests."""
        elapsed = time.time() - self.last_request
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)
        self.last_request = time.time()

    def search_articles(
        self,
        query: str,
        max_results: int = 100,
        date_range_years: int = 5
    ) -> List[str]:
        """Search PubMed articles and return PMIDs."""
        logger.info(f"Searching PubMed for: {query}")

        self._rate_limit_request()

        # Construct search URL
        search_url = f"{self.BASE_URL}esearch.fcgi"

        # Date range
        end_date = datetime.now()
        start_date = end_date - timedelta(days=365 * date_range_years)

        params = {
            'db': 'pubmed',
            'term': query,
            'retmax': max_results,
            'datetype': 'pdat',
            'mindate': start_date.strftime('%Y/%m/%d'),
            'maxdate': end_date.strftime('%Y/%m/%d'),
            'retmode': 'json',
            'email': self.email
        }

        if self.api_key:
            params['api_key'] = self.api_key

        try:
            response = requests.get(search_url, params=params)
            response.raise_for_status()

            data = response.json()
            pmids = data.get('esearchresult', {}).get('idlist', [])

            logger.info(f"Found {len(pmids)} articles for query: {query}")
            return pmids

        except Exception as e:
            logger.error(f"Error searching PubMed: {e}")
            return []

    def fetch_article_details(self, pmids: List[str]) -> List[Dict[str, Any]]:
        """Fetch detailed article information for given PMIDs."""
        if not pmids:
            return []

        logger.info(f"Fetching details for {len(pmids)} articles")

        # Process in batches to avoid URL length limits
        batch_size = 200
        all_articles = []

        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i:i + batch_size]
            articles = self._fetch_batch_details(batch_pmids)
            all_articles.extend(articles)

            if i + batch_size < len(pmids):
                time.sleep(self.rate_limit)

        return all_articles

    def _fetch_batch_details(self, pmids: List[str]) -> List[Dict[str, Any]]:
        """Fetch details for a batch of PMIDs."""
        self._rate_limit_request()

        fetch_url = f"{self.BASE_URL}efetch.fcgi"

        params = {
            'db': 'pubmed',
            'id': ','.join(pmids),
            'retmode': 'xml',
            'email': self.email
        }

        if self.api_key:
            params['api_key'] = self.api_key

        try:
            response = requests.get(fetch_url, params=params)
            response.raise_for_status()

            return self._parse_pubmed_xml(response.text)

        except Exception as e:
            logger.error(f"Error fetching article details: {e}")
            return []

    def _parse_pubmed_xml(self, xml_content: str) -> List[Dict[str, Any]]:
        """Parse PubMed XML response into structured data."""
        articles = []

        try:
            root = ET.fromstring(xml_content)

            for article_elem in root.findall('.//PubmedArticle'):
                article = self._extract_article_info(article_elem)
                if article:
                    articles.append(article)

        except ET.ParseError as e:
            logger.error(f"Error parsing XML: {e}")

        return articles

    def _extract_article_info(self, article_elem) -> Optional[Dict[str, Any]]:
        """Extract information from a single article XML element."""
        try:
            # Basic article info
            medline_citation = article_elem.find('MedlineCitation')
            if medline_citation is None:
                return None

            pmid_elem = medline_citation.find('PMID')
            pmid = pmid_elem.text if pmid_elem is not None else None

            article_elem_inner = medline_citation.find('Article')
            if article_elem_inner is None:
                return None

            # Title
            title_elem = article_elem_inner.find('ArticleTitle')
            title = title_elem.text if title_elem is not None else ""

            # Abstract
            abstract_elem = article_elem_inner.find('Abstract/AbstractText')
            abstract = abstract_elem.text if abstract_elem is not None else ""

            # Authors
            authors = []
            author_list = article_elem_inner.find('AuthorList')
            if author_list is not None:
                for author in author_list.findall('Author'):
                    lastname = author.find('LastName')
                    firstname = author.find('ForeName')
                    if lastname is not None and firstname is not None:
                        authors.append(f"{firstname.text} {lastname.text}")

            # Journal
            journal_elem = article_elem_inner.find('Journal/Title')
            journal = journal_elem.text if journal_elem is not None else ""

            # Publication date
            pub_date = article_elem_inner.find('Journal/JournalIssue/PubDate')
            year = ""
            if pub_date is not None:
                year_elem = pub_date.find('Year')
                year = year_elem.text if year_elem is not None else ""

            # Keywords/MeSH terms
            keywords = []
            mesh_list = medline_citation.find('MeshHeadingList')
            if mesh_list is not None:
                for mesh_heading in mesh_list.findall('MeshHeading'):
                    descriptor = mesh_heading.find('DescriptorName')
                    if descriptor is not None:
                        keywords.append(descriptor.text)

            return {
                'pmid': pmid,
                'title': title,
                'abstract': abstract,
                'authors': authors,
                'journal': journal,
                'year': year,
                'keywords': keywords,
                'full_text': f"{title} {abstract}".lower()
            }

        except Exception as e:
            logger.error(f"Error extracting article info: {e}")
            return None

    def search_aging_crispr_literature(self) -> pd.DataFrame:
        """Search for aging and CRISPR-related literature."""
        logger.info("Searching for aging and CRISPR literature")

        # Define search queries for different aspects
        queries = {
            'aging_crispr': 'aging AND CRISPR AND (rejuvenation OR senescence)',
            'osk_reprogramming': 'OSK OR "Yamanaka factors" AND aging',
            'senolytic_crispr': 'senolytic AND CRISPR',
            'epigenetic_clock': 'epigenetic clock AND (CRISPR OR aging)',
            'cellular_rejuvenation': 'cellular rejuvenation AND genome editing',
            'senescence_intervention': 'cellular senescence AND intervention'
        }

        all_articles = []

        for category, query in queries.items():
            logger.info(f"Searching for {category}: {query}")

            # Check cache first
            cache_file = self.cache_dir / f"{category}_articles.json"
            if cache_file.exists():
                with open(cache_file, 'r') as f:
                    cached_articles = json.load(f)
                    for article in cached_articles:
                        article['search_category'] = category
                    all_articles.extend(cached_articles)
                continue

            # Search PubMed
            pmids = self.search_articles(query, max_results=50)
            articles = self.fetch_article_details(pmids)

            # Add category
            for article in articles:
                article['search_category'] = category

            # Cache results
            with open(cache_file, 'w') as f:
                json.dump(articles, f, indent=2)

            all_articles.extend(articles)

            # Rate limiting between categories
            time.sleep(1)

        df = pd.DataFrame(all_articles)

        if not df.empty:
            # Remove duplicates by PMID
            df = df.drop_duplicates(subset=['pmid'])

            # Save combined results
            output_path = self.cache_dir / "aging_crispr_literature.csv"
            df.to_csv(output_path, index=False)
            logger.info(f"Saved {len(df)} unique articles to {output_path}")

        return df

    def analyze_gene_literature(self, gene: str) -> Dict[str, Any]:
        """Analyze literature for a specific gene in aging context."""
        logger.info(f"Analyzing literature for gene: {gene}")

        # Search for gene + aging terms
        query = f'"{gene}" AND (aging OR senescence OR longevity OR rejuvenation)'
        pmids = self.search_articles(query, max_results=30)
        articles = self.fetch_article_details(pmids)

        if not articles:
            return {
                'gene': gene,
                'article_count': 0,
                'aging_relevance_score': 0.0,
                'top_keywords': [],
                'recent_articles': 0
            }

        # Analyze content
        all_text = ' '.join([
            article.get('full_text', '') for article in articles
        ])

        # Count aging-related terms
        aging_terms = [
            'senescence', 'aging', 'longevity', 'rejuvenation',
            'epigenetic clock', 'cellular aging', 'anti-aging',
            'lifespan', 'healthspan', 'age-related'
        ]

        term_counts = {}
        for term in aging_terms:
            count = len(re.findall(rf'\b{re.escape(term)}\b', all_text, re.IGNORECASE))
            term_counts[term] = count

        # Calculate relevance score
        total_terms = sum(term_counts.values())
        article_count = len(articles)
        relevance_score = min(1.0, (total_terms / max(1, article_count)) / 10)

        # Recent articles (last 2 years)
        current_year = datetime.now().year
        recent_articles = sum(
            1 for article in articles
            if article.get('year', '').isdigit() and
            int(article['year']) >= current_year - 2
        )

        # Top keywords
        all_keywords = []
        for article in articles:
            all_keywords.extend(article.get('keywords', []))

        keyword_counts = {}
        for keyword in all_keywords:
            keyword_counts[keyword] = keyword_counts.get(keyword, 0) + 1

        top_keywords = sorted(
            keyword_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:10]

        return {
            'gene': gene,
            'article_count': article_count,
            'aging_relevance_score': relevance_score,
            'term_frequencies': term_counts,
            'top_keywords': [kw[0] for kw in top_keywords],
            'recent_articles': recent_articles,
            'articles': articles
        }

    def create_literature_features(self, genes: List[str]) -> pd.DataFrame:
        """Create literature-based features for a list of genes."""
        logger.info(f"Creating literature features for {len(genes)} genes")

        features = []

        for gene in genes:
            try:
                analysis = self.analyze_gene_literature(gene)

                features.append({
                    'gene': gene,
                    'literature_count': analysis['article_count'],
                    'aging_relevance_score': analysis['aging_relevance_score'],
                    'recent_publications': analysis['recent_articles'],
                    'senescence_mentions': analysis['term_frequencies'].get('senescence', 0),
                    'longevity_mentions': analysis['term_frequencies'].get('longevity', 0),
                    'rejuvenation_mentions': analysis['term_frequencies'].get('rejuvenation', 0)
                })

                # Rate limiting
                time.sleep(self.rate_limit)

            except Exception as e:
                logger.error(f"Error analyzing {gene}: {e}")
                features.append({
                    'gene': gene,
                    'literature_count': 0,
                    'aging_relevance_score': 0.0,
                    'recent_publications': 0,
                    'senescence_mentions': 0,
                    'longevity_mentions': 0,
                    'rejuvenation_mentions': 0
                })

        df = pd.DataFrame(features)

        # Save features
        output_path = self.cache_dir / "gene_literature_features.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Saved literature features to {output_path}")

        return df


def search_aging_literature(
    genes: Optional[List[str]] = None,
    api_key: Optional[str] = None
) -> pd.DataFrame:
    """Convenience function to search aging literature."""
    miner = PubMedMiner(api_key=api_key)

    if genes:
        return miner.create_literature_features(genes)
    else:
        return miner.search_aging_crispr_literature()
