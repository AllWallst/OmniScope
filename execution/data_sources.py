import yfinance as yf
import feedparser
import requests
from bs4 import BeautifulSoup
import pandas as pd
from Bio import Entrez
from youtube_transcript_api import YouTubeTranscriptApi
from youtube_transcript_api._errors import TranscriptsDisabled, NoTranscriptFound
from textblob import TextBlob
import re
import urllib.parse

# Configure Entrez email for PubMed (Best practice)
Entrez.email = "info@allwallst.com" 

USER_AGENT = "OmniScope/1.0 (info@allwallst.com)"

def fetch_market_data(ticker_symbol, period="1mo"):
    """
    Fetches market data and company info using yfinance.
    """
    try:
        ticker = yf.Ticker(ticker_symbol)
        
        # 1. Info & History (Essential)
        try:
            info = ticker.info
            history = ticker.history(period=period)
        except Exception:
            # If info fails completely, return None as we need basic data
            return None

        # 2. Financials (Optional - prone to 404s)
        financials = pd.DataFrame()
        balance_sheet = pd.DataFrame()
        cashflow = pd.DataFrame()
        try:
            financials = ticker.financials
            balance_sheet = ticker.balance_sheet
            cashflow = ticker.cashflow
        except Exception:
            pass # Use empty DFs
            
        # 3. Recommendations & Holders (Optional)
        recommendations = pd.DataFrame()
        major_holders = pd.DataFrame()
        institutional_holders = pd.DataFrame()
        try:
           recommendations = ticker.recommendations
           major_holders = ticker.major_holders
           institutional_holders = ticker.institutional_holders
        except Exception:
           pass

        return {
            "info": info,
            "history": history,
            "financials": financials,
            "balance_sheet": balance_sheet,
            "cashflow": cashflow,
            "recommendations": recommendations,
            "major_holders": major_holders,
            "institutional_holders": institutional_holders
        }
    except Exception as e:
        print(f"Error fetching market data for {ticker_symbol}: {e}")
        return None

def fetch_news(query):
    """
    Fetches news from Google News RSS based on a query.
    """
    encoded_query = urllib.parse.quote(query)
    rss_url = f"https://news.google.com/rss/search?q={encoded_query}&hl=en-US&gl=US&ceid=US:en"
    feed = feedparser.parse(rss_url)
    
    news_items = []
    for entry in feed.entries[:10]:
        news_items.append({
            "title": entry.title,
            "link": entry.link,
            "published": entry.published,
            "source": entry.source.title if 'source' in entry else "Unknown"
        })
    return news_items

def fetch_sec_filings(ticker_symbol, limit=5):
    """
    Fetches recent 10-K and 10-Q filings from SEC EDGAR.
    Note: Real production use would require CIK lookup. 
    Here we attempt a direct search or use a known CIK mapping if available.
    For simplicity in this demo, we'll try to find the CIK via standard headers or search.
    However, browsing EDGAR programmatically often requires a CIK. 
    
    Strategy: Use the SEC's company tickers JSON to find CIK first.
    """
    headers = {"User-Agent": USER_AGENT}
    
    # 1. Get CIK (Mocking this part or doing a quick lookup is complex without a local DB)
    # We will try to rely on the SEC Company Tickers JSON
    try:
        tickers_json = requests.get("https://www.sec.gov/files/company_tickers.json", headers=headers).json()
        
        cik_str = None
        for key in tickers_json:
            entry = tickers_json[key]
            if entry['ticker'] == ticker_symbol.upper():
                cik_str = str(entry['cik_str']).zfill(10)
                break
        
        if not cik_str:
            return []

        # 2. Fetch submissions
        url = f"https://data.sec.gov/submissions/CIK{cik_str}.json"
        resp = requests.get(url, headers=headers)
        data = resp.json()
        
        filings = data.get('filings', {}).get('recent', {})
        results = []
        
        if filings:
            forms = filings.get('form', [])
            dates = filings.get('filingDate', [])
            accession_numbers = filings.get('accessionNumber', [])
            primary_documents = filings.get('primaryDocument', [])
            
            for i in range(min(limit, len(forms))):
                if forms[i] in ['10-K', '10-Q', '8-K', '13F-HR', '13F-HR/A']:
                    # Construct link
                    # https://www.sec.gov/Archives/edgar/data/{cik}/{accession}/{primaryDoc}
                    acc_no_dash = accession_numbers[i].replace('-', '')
                    link = f"https://www.sec.gov/Archives/edgar/data/{int(cik_str)}/{acc_no_dash}/{primary_documents[i]}"
                    
                    description = f"SEC Filing {forms[i]}"
                    if '13F' in forms[i]:
                        description += " (Institutional Holdings)"

                    results.append({
                        "form": forms[i],
                        "date": dates[i],
                        "link": link,
                        "description": description
                    })
                    
        return results
        
    except Exception as e:
        print(f"Error fetching SEC filings: {e}")
        return []

def fetch_pubmed_abstracts(term, max_results=5):
    """
    Searches PubMed for a term (biosciences/medical context).
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if not id_list:
            return []
            
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        # Parsing MEDLINE format is a bit raw, let's use XML instead for structured data
        handle.close()
        
        handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        records = Entrez.read(handle)
        
        abstracts = []
        for article in records['PubmedArticle']:
            try:
                medline = article['MedlineCitation']
                pmid = str(medline['PMID'])
                title = medline['Article']['ArticleTitle']
                abstract_list = medline['Article'].get('Abstract', {}).get('AbstractText', [])
                abstract_text = " ".join([str(x) for x in abstract_list])
                
                abstracts.append({
                    "title": title,
                    "abstract": abstract_text[:500] + "..." if len(abstract_text) > 500 else abstract_text,
                    "source": "PubMed",
                    "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
            except:
                continue
        return abstracts
        
    except Exception as e:
        print(f"Error fetching PubMed: {e}")
        return []

def fetch_congress_trading(ticker_symbol, api_key=None):
    """
    Fetches congressional trading data by scraping Quiver Quantitative.
    Source: https://www.quiverquant.com/congresstrading/stock/{TICKER}
    This is more direct than Capitol Trades and handles history better for specific stocks.
    """
    import requests
    import pandas as pd
    from io import StringIO
    
    url = f"https://www.quiverquant.com/congresstrading/stock/{ticker_symbol.strip().upper()}"
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"}
    
    debug_log = [f"Scraping {url}"]
    
    try:
        r = requests.get(url, headers=headers, timeout=10)
        if r.status_code != 200:
            debug_log.append(f"Status Code {r.status_code}")
            df_empty = pd.DataFrame()
            df_empty.attrs['debug_log'] = debug_log
            return df_empty
            
        try:
            dfs = pd.read_html(StringIO(r.text))
        except ImportError:
            # Fallback to bs4 if lxml is missing (common in some environments)
            # We know bs4 is installed since we use it for news
            dfs = pd.read_html(StringIO(r.text), flavor='bs4')
            
        if not dfs:
            debug_log.append("No tables found via read_html")
            df_empty = pd.DataFrame()
            df_empty.attrs['debug_log'] = debug_log
            return df_empty
            
        df = dfs[0]
        debug_log.append(f"Table found with {len(df)} rows. Columns: {df.columns.tolist()}")
        
        trades = []
        for idx, row in df.iterrows():
            try:
                # Parse Transaction (Type + Amount)
                raw_trans = str(row.get('Transaction', '')).strip()
                
                trans_type = "Unknown"
                amount = "Unknown"
                lower_trans = raw_trans.lower()
                
                # Logic to identify Sale/Purchase
                if "sale" in lower_trans:
                    trans_type = "Sale"
                    if '$' in raw_trans:
                        amount = raw_trans[raw_trans.find('$'):]
                    else:
                        parts = raw_trans.split(' ', 1)
                        if len(parts) > 1: amount = parts[1]
                elif "purchase" in lower_trans:
                    trans_type = "Purchase"
                    if '$' in raw_trans:
                        amount = raw_trans[raw_trans.find('$'):]
                    elif "exchange" in lower_trans:
                         trans_type = "Exchange"
                    else:
                         parts = raw_trans.split(' ', 1)
                         if len(parts) > 1: amount = parts[1]
                else:
                    parts = raw_trans.split(' ', 1)
                    if len(parts) > 0: trans_type = parts[0]
                    if len(parts) > 1: amount = parts[1]

                # Parse Politician
                raw_pol = str(row.get('Politician', ''))
                chamber = "Unknown"
                party = "Unknown"
                rep_name = raw_pol
                
                if "House" in raw_pol:
                    chamber = "House"
                    split_val = raw_pol.split("House")
                    if len(split_val) > 0: rep_name = split_val[0].strip()
                elif "Senate" in raw_pol:
                    chamber = "Senate"
                    split_val = raw_pol.split("Senate")
                    if len(split_val) > 0: rep_name = split_val[0].strip()
                
                if "/ D" in raw_pol: party = "Democrat"
                elif "/ R" in raw_pol: party = "Republican"
                    
                trades.append({
                    "Representative": rep_name,
                    "Chamber": chamber,
                    "Party": party,
                    "Transaction Date": row.get('Traded', 'N/A'),
                    "Disclosed": row.get('Filed', 'N/A'),
                    "Transaction": trans_type,
                    "Amount": amount,
                    "Price": "-" 
                })
            except Exception as e:
                debug_log.append(f"Error parsing row {idx}: {e}")
                continue
                
        if trades:
            debug_log.append(f"Successfully parsed {len(trades)} trades")
            res_df = pd.DataFrame(trades)
            res_df.attrs['debug_log'] = debug_log
            return res_df
            
        debug_log.append("Trades list empty after loop. Check parsing logic against row content.")
        if not df.empty:
             debug_log.append(f"First row raw: {df.iloc[0].to_dict()}")
             
        df_empty = pd.DataFrame()
        df_empty.attrs['debug_log'] = debug_log
        return df_empty

    except Exception as e:
        debug_log.append(f"Error scraping Quiver Quant: {e}")
        df_empty = pd.DataFrame()
        df_empty.attrs['debug_log'] = debug_log
        return df_empty
                
        if trades:
            return pd.DataFrame(trades)
            
        return pd.DataFrame()

    except Exception as e:
        print(f"Error scraping Quiver Quant: {e}")
        return pd.DataFrame()

def fetch_youtube_transcripts(ticker_symbol, limit=3):
    """
    Searches YouTube for '[Ticker] Earnings Call' and gets transcripts.
    This is experimental and relies on video search results + caption availability.
    """
    # 1. Search for video IDs (We need a way to search. YouTube Data API is best but requires Key.
    #    We can scrape search results lightly or ask user for key. 
    #    For this demo, we will use a workaround or skip if too complex without API key.
    #    Let's try a very simple scrape of youtube search results page to get IDs.)
    
    query = f"{ticker_symbol} earnings call"
    search_url = f"https://www.youtube.com/results?search_query={urllib.parse.quote(query)}"
    
    try:
        # Note: Scaping YouTube directly is brittle. 
        # A better approach without API key is tricky. 
        # For the sake of "Real Data" demo, we'll try to extract the first video ID found in the HTML.
        headers = {"User-Agent": "Mozilla/5.0"}
        resp = requests.get(search_url, headers=headers)
        
        # Regex to find video IDs
        video_ids = re.findall(r"watch\?v=(\S{11})", resp.text)
        unique_ids = []
        for vid in video_ids:
            if vid not in unique_ids:
                unique_ids.append(vid)
                
        results = []
        count = 0
        for vid in unique_ids:
            if count >= limit: break
            try:
                transcript = YouTubeTranscriptApi.get_transcript(vid)
                # Combine first few lines as preview
                full_text = " ".join([t['text'] for t in transcript[:20]]) 
                results.append({
                    "video_id": vid,
                    "link": f"https://www.youtube.com/watch?v={vid}",
                    "transcript_snippet": full_text + "..."
                })
                count += 1
            except (TranscriptsDisabled, NoTranscriptFound):
                continue
            except Exception as e:
                print(f"Error getting transcript for {vid}: {e}")
                continue
                
        return results

    except Exception as e:
        print(f"Error searching YouTube: {e}")
        return []

def analyze_sentiment(text_list):
    """
    Analyzes the sentiment of a list of text strings.
    Returns average polarity (-1 to 1) and subjective score.
    """
    if not text_list:
        return 0, 0
    
    total_polarity = 0
    total_subjectivity = 0
    count = 0
    
    for text in text_list:
        if text:
            blob = TextBlob(text)
            total_polarity += blob.sentiment.polarity
            total_subjectivity += blob.sentiment.subjectivity
            count += 1
            
    if count == 0:
        return 0, 0
        
    return total_polarity / count, total_subjectivity / count

def lookup_ticker(query):
    """
    Searches for a ticker symbol using Yahoo Finance's public autocomplete API.
    Returns the top match ticker (str) or None.
    """
    try:
        url = f"https://query2.finance.yahoo.com/v1/finance/search?q={urllib.parse.quote(query)}"
        headers = {'User-Agent': USER_AGENT}
        resp = requests.get(url, headers=headers)
        data = resp.json()
        
        if 'quotes' in data and len(data['quotes']) > 0:
            # Prefer equity results
            for quote in data['quotes']:
                 if quote.get('quoteType') == 'EQUITY' and quote.get('isYahooFinance', True):
                     return quote['symbol']
            # Fallback to first result
            return data['quotes'][0]['symbol']
            
        return None
    except Exception as e:
        print(f"Error looking up ticker: {e}")
        return None

def fetch_reddit_posts(query, limit=10):
    """
    Fetches recent Reddit posts mentioning the query from high-signal subreddits.
    Uses public Reddit JSON pages (No API key required, but strict rate limits).
    """
    # Search specific investing subreddits for higher signal-to-noise
    subreddits = "stocks+investing+wallstreetbets+stockmarket"
    url = f"https://www.reddit.com/r/{subreddits}/search.json?q={urllib.parse.quote(query)}&sort=new&restrict_sr=1&limit={limit}"
    
    headers = {
        "User-Agent": "OmniScope/1.0 (Contact: info@example.com)"
    }
    
    try:
        resp = requests.get(url, headers=headers)
        if resp.status_code != 200:
            print(f"Reddit error: {resp.status_code}")
            return []
            
        data = resp.json()
        posts = []
        
        for item in data.get('data', {}).get('children', []):
            post = item['data']
            posts.append({
                "title": post.get('title'),
                "url": f"https://www.reddit.com{post.get('permalink')}",
                "subreddit": post.get('subreddit'),
                "score": post.get('score'),
                "comments": post.get('num_comments'),
                "created_utc": post.get('created_utc'),
                "selftext": post.get('selftext', '')[:500] + "..." if len(post.get('selftext', '')) > 500 else post.get('selftext', '')
            })
            
        return posts
    except Exception as e:
        print(f"Error fetching Reddit: {e}")
        return []

def fetch_insider_transactions(ticker_symbol):
    """
    Fetches insider transactions using yfinance.
    """
    try:
        ticker = yf.Ticker(ticker_symbol)
        # insider_transactions is a DataFrame
        insiders = ticker.insider_transactions
        if insiders is not None and not insiders.empty:
            # Clean up for display if needed
            return insiders
        return pd.DataFrame()
    except Exception as e:
        print(f"Error fetching insider data: {e}")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error fetching insider data: {e}")
        return pd.DataFrame()

def fetch_esg_data(ticker_symbol):
    """
    Fetches sustainability/ESG scores.
    """
    try:
        ticker = yf.Ticker(ticker_symbol)
        esg = ticker.sustainability
        if esg is not None and not esg.empty:
            return esg
        return pd.DataFrame()
    except Exception as e:
        print(f"Error fetching ESG data: {e}")
        return pd.DataFrame()

def fetch_options_sentiment(ticker_symbol):
    """
    Calculates Put/Call Ratio for the nearest expiration date.
    Returns: Ratio (float), Expiration Date (str), Total Puts, Total Calls
    """
    try:
        ticker = yf.Ticker(ticker_symbol)
        expirations = ticker.options
        if not expirations:
            return None, None, 0, 0
            
        # Get nearest expiration
        nearest_date = expirations[0]
        chain = ticker.option_chain(nearest_date)
        
        calls = chain.calls
        puts = chain.puts
        
        total_call_vol = calls['volume'].sum() if not calls.empty else 0
        total_put_vol = puts['volume'].sum() if not puts.empty else 0
        
        if total_call_vol == 0:
            return 0.0, nearest_date, total_put_vol, total_call_vol
            
        pc_ratio = total_put_vol / total_call_vol
        return pc_ratio, nearest_date, total_put_vol, total_call_vol
        
    except Exception as e:
        print(f"Error fetching options data: {e}")
        return None, None, 0, 0
