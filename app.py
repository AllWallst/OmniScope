import streamlit as st
import pandas as pd
import sys
import os
import plotly.express as px

# Add project root to path to import local modules
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from projects.alpha_stream.execution import data_sources

st.set_page_config(layout="wide", page_title="OmniScope", page_icon="🔭")

# Load Custom CSS
def load_css():
    st.markdown("""
        <style>
        /* Modern Dark Theme with Glassmorphism */
        .stApp {
            background-color: #0e1117;
            font-family: 'Inter', sans-serif;
        }
        
        /* Sidebar styling */
        .css-1d391kg {
            background-color: #161b22;
        }
        
        /* Card styling for results */
        .doc-card {
            background-color: #ffffff;
            border-radius: 12px;
            padding: 20px;
            margin-bottom: 15px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            transition: transform 0.2s, box-shadow 0.2s;
            border-left: 5px solid #00D4FF; /* Cyan accent */
        }
        .doc-card:hover {
            transform: translateY(-2px);
            box-shadow: 0 8px 12px rgba(0,0,0,0.15);
        }
        .doc-meta {
            font-size: 0.8em;
            color: #666;
            margin-bottom: 5px;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }
        .doc-title {
            font-size: 1.1em;
            font-weight: 600;
            color: #1a73e8;
            text-decoration: none;
            display: block;
            margin-bottom: 8px;
        }
        .doc-title:hover {
            text-decoration: underline;
        }
        .doc-card p {
            color: #333333; /* Force dark text for readability on white card */
            font-size: 0.95em;
            line-height: 1.5;
        }
        </style>
    """, unsafe_allow_html=True)

load_css()

# Sidebar
with st.sidebar:
    st.header("OmniScope")
    st.markdown("---")
    
    st.subheader("Data Sources")
    use_sec = st.checkbox("Regulatory Filings (SEC)", value=True)
    use_news = st.checkbox("News & Events", value=True)
    use_transcripts = st.checkbox("Expert Call Transcripts", value=False)
    use_pubmed = st.checkbox("Life Sciences (PubMed)", value=False, help="Searches the National Library of Medicine (PubMed) for biomedical literature. Useful for tracking clinical trials, drug research, and medical studies related to the company.")
    use_reddit = st.checkbox("Social Listening (Reddit)", value=True, help="Scrapes r/stocks, r/investing, etc. for discussions.")
    
    st.markdown("---")
    st.subheader("Filters")
    date_range = st.date_input("Date Range", [])

# Main Layout
st.title("Market Intelligence Search")

# Search Bar
col_search, col_period = st.columns([4, 1])
with col_search:
    query = st.text_input("Search Ticker (e.g., AAPL) or Keyword", placeholder="Enter ticker or keyword...")
with col_period:
    period = st.selectbox("Graph Period", ["1d", "5d", "1mo", "6mo", "ytd", "1y", "5y", "10y", "max"], index=2)

# Results container
if query:
    st.subheader(f"Results for '{query}'")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # 1. Market Data (if ticker like)
        if query.isalpha() and len(query) < 6:
            with st.spinner("Fetching Market Data..."):
                market_data = data_sources.fetch_market_data(query, period=period)
                if market_data:
                    info = market_data['info']
                    current_price = info.get('currentPrice', 'N/A')
                    previous_close = info.get('previousClose', 0)
                    
                    # Calculate Delta
                    delta = None
                    if isinstance(current_price, (int, float)) and isinstance(previous_close, (int, float)):
                        delta = current_price - previous_close
                    
                    st.markdown(f"### {info.get('shortName', query.upper())} ({query.upper()})")
                    st.metric(
                        "Current Price", 
                        f"${current_price}", 
                        f"{delta:+.2f}" if delta is not None else None
                    )
                    st.caption(f"Recommendation: {info.get('recommendationKey', '').replace('_', ' ').title()}")
                    
                    # Enhanced Charting with Plotly
                    hist_data = market_data['history'].reset_index()
                    # Ensure Date column is formatted for Plotly
                    
                    fig = px.line(hist_data, x="Date", y="Close", title=f"{query.upper()} Price History")
                    fig.update_layout(
                        xaxis_title="", 
                        yaxis_title="Price (USD)",
                        margin=dict(l=0, r=0, t=30, b=0),
                        hovermode="x unified"
                    )
                    # Force y-axis to not start at 0 if the variation is small relative to price
                    # Plotly default is usually good, but we can enforce it:
                    fig.update_yaxes(autorange=True, fixedrange=False) 
                    
                    st.plotly_chart(fig, use_container_width=True)

        all_results = []
        
        # 2. News
        if use_news:
            with st.spinner("Fetching News..."):
                news = data_sources.fetch_news(query)
                for item in news:
                    item['type'] = 'News'
                    all_results.append(item)
                    
        # 3. SEC Filings
        if use_sec and query.isalpha() and len(query) < 6:
            with st.spinner("Fetching SEC Filings (including 13-F)..."):
                filings = data_sources.fetch_sec_filings(query)
                for item in filings:
                    item['title'] = item['description']
                    item['type'] = 'SEC Filing'
                    all_results.append(item)
        
        # 4. PubMed
        if use_pubmed:
            with st.spinner("Fetching Life Sciences Data..."):
                abstracts = data_sources.fetch_pubmed_abstracts(query)
                for item in abstracts:
                    item['type'] = 'PubMed'
                    all_results.append(item)

        # 5. Transcripts
        if use_transcripts and query.isalpha() and len(query) < 6:
            with st.spinner("Fetching Transcripts..."):
                transcripts = data_sources.fetch_youtube_transcripts(query)
                for item in transcripts:
                    item['type'] = 'Transcript'
                    item['title'] = f"{query.upper()} Earnings Call (Video)"
                    all_results.append(item)

        # 6. Reddit
        if use_reddit and len(query) > 1:
            with st.spinner("Listening to Social Media..."):
                reddit_posts = data_sources.fetch_reddit_posts(query)
                for item in reddit_posts:
                    item['type'] = 'Reddit'
                    item['source'] = f"r/{item['subreddit']}"
                    item['link'] = item['url']
                    item['description'] = item['selftext']
                    all_results.append(item)

        # Display Results in Tabs
        if not all_results and not market_data:
            st.info("No results found. Try a different query or enable more sources.")
        else:
            tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs(["Feed", "Financials", "Holders", "Analysis", "Social", "Insider", "Short Data", "Management"])
            
            with tab1:
                # Iterate and show documents
                for doc in all_results:
                    st.markdown(f"""
                    <div class="doc-card">
                        <div class="doc-meta">{doc.get('type')} | {doc.get('source', 'OmniScope')} | {doc.get('created_utc', doc.get('published', doc.get('date', '')))}</div>
                        <a href="{doc.get('link')}" target="_blank" class="doc-title">{doc.get('title')}</a>
                        <p>{doc.get('transcript_snippet', doc.get('description', doc.get('abstract', '')))}</p>
                    </div>
                    """, unsafe_allow_html=True)
            
            with tab2:
                if market_data:
                    st.subheader("Financial Statements")
                    st.markdown("#### Balance Sheet")
                    st.dataframe(market_data.get('balance_sheet'))
                    st.markdown("#### Income Statement")
                    st.dataframe(market_data.get('financials'))
                    st.markdown("#### Cash Flow")
                    st.dataframe(market_data.get('cashflow'))
                else:
                    st.info("No Financial data available.")

            with tab3:
                if market_data:
                    st.subheader("Major Holders")
                    st.dataframe(market_data.get('major_holders'))
                    st.subheader("Institutional Holders")
                    st.dataframe(market_data.get('institutional_holders'))
                else:
                    st.info("No Holders data available.")

            with tab4:
                col_a1, col_a2 = st.columns(2)
                with col_a1:
                    if market_data:
                        st.subheader("Analyst Recommendations")
                        st.dataframe(market_data.get('recommendations'))
                    else:
                        st.info("No Analysis data available.")
                
                with col_a2:
                    if market_data and query.isalpha():
                        st.subheader("ESG Risk Ratings")
                        with st.spinner("Fetching ESG Data..."):
                            esg_data = data_sources.fetch_esg_data(query)
                            if not esg_data.empty:
                                st.dataframe(esg_data)
                            else:
                                st.info("No ESG data found.")
                    
                    st.markdown("---")
                    st.subheader("Options Sentiment")
                    with st.spinner("Analyzing Options Chain..."):
                        pc_ratio, exp_date, puts, calls = data_sources.fetch_options_sentiment(query)
                        if pc_ratio is not None:
                            st.metric("Put/Call Ratio", f"{pc_ratio:.2f}", help=">1 Bearish, <1 Bullish")
                            st.caption(f"Expiration: {exp_date} | Calls: {calls:,} | Puts: {puts:,}")
                        else:
                            st.info("No Options data available.")

            with tab5:
                # Filter for Reddit results
                reddit_items = [x for x in all_results if x.get('type') == 'Reddit']
                if reddit_items:
                    st.subheader("Reddit Discussions")
                    for post in reddit_items:
                         st.markdown(f"""
                        <div class="doc-card">
                            <div class="doc-meta">r/{post.get('subreddit')} | Score: {post.get('score')} | Comments: {post.get('comments')}</div>
                            <a href="{post.get('link')}" target="_blank" class="doc-title">{post.get('title')}</a>
                            <p>{post.get('description')}</p>
                        </div>
                        """, unsafe_allow_html=True)
                else:
                    st.info("No recent social discussions found.")

            with tab6:
                if market_data and query.isalpha():
                    st.subheader("Insider Transactions")
                    with st.spinner("Fetching Insider Data..."):
                        insiders = data_sources.fetch_insider_transactions(query)
                        if not insiders.empty:
                            st.dataframe(insiders)
                        else:
                            st.info("No Insider transactions found.")
                else:
                    st.info("Insider data requires a valid ticker.")

            with tab7:
                if market_data:
                    st.subheader("Short Interest Analysis")
                    info = market_data.get('info', {})
                    
                    short_float = info.get('shortPercentOfFloat', 0)
                    short_ratio = info.get('shortRatio', 0)
                    shares_short = info.get('sharesShort', 0)
                    
                    col_s1, col_s2, col_s3 = st.columns(3)
                    
                    with col_s1:
                        val = f"{short_float * 100:.2f}%" if short_float else "N/A"
                        st.metric("Short % of Float", val)
                    with col_s2:
                        st.metric("Short Ratio (Days)", short_ratio)
                    with col_s3:
                        st.metric("Shares Short", f"{shares_short:,}" if shares_short else "N/A")
                        
                    st.markdown("---")
                    
                    # Squeeze Indicator
                    if short_float and short_float > 0.20:
                        st.error("⚠️ **High Short Squeeze Potential**: Short Interest > 20% of Float")
                    elif short_float and short_float > 0.10:
                        st.warning("⚠️ **Moderate Squeeze Potential**: Short Interest > 10% of Float")
                    else:
                        st.success("✅ **Low Squeeze Potential**: Short Interest is low.")
                        
                    st.caption("Data Source: Yahoo Finance")
                else:
                    st.info("Short data requires a valid ticker.")
            
            with tab8:
                if market_data:
                    st.subheader("Key Executives")
                    officers = market_data.get('info', {}).get('companyOfficers', [])
                    if officers:
                        # Normalize data for dataframe
                        clean_officers = []
                        for officer in officers:
                            clean_officers.append({
                                "Name": officer.get('name', 'N/A'),
                                "Title": officer.get('title', 'N/A'),
                                "Pay": f"${officer.get('totalPay', 0):,}" if officer.get('totalPay') else "N/A",
                                "Year Born": officer.get('yearBorn', 'N/A')
                            })
                        st.dataframe(pd.DataFrame(clean_officers))
                    else:
                        st.info("No Executive info available.")
                else:
                    st.info("Enter a ticker to see Management.")
                
    with col2:
        st.markdown("### Quick Insights")
        
        if query:
            st.write(f"Search analysis for **{query}**")
            
            # Real Sentiment Analysis
            sentences_to_analyze = []
            
            # Add news titles and reddit titles
            for item in all_results:
                if item.get('type') in ['News', 'Reddit']:
                    sentences_to_analyze.append(item.get('title', ''))
                    if item.get('type') == 'Reddit':
                        sentences_to_analyze.append(item.get('description', '')[:200]) # Add some body text
            
            # Add market analysis if available
            if market_data and not market_data['recommendations'].empty:
                 # Add some recommendation text to analysis if possible, 
                 # for now we stick to News headlines as primary real-time sentiment source
                 pass

            polarity, subjectivity = data_sources.analyze_sentiment(sentences_to_analyze)
            
            # Determine sentiment label
            if polarity > 0.1:
                sentiment_label = "Positive"
                sentiment_color = "green"
            elif polarity < -0.1:
                sentiment_label = "Negative"
                sentiment_color = "red"
            else:
                sentiment_label = "Neutral"
                sentiment_color = "gray"
                
            st.markdown(f"**Sentiment**: <span style='color:{sentiment_color};font-weight:bold'>{sentiment_label}</span> ({polarity:.2f})", unsafe_allow_html=True)
            st.write(f"**Subjectivity**: {subjectivity:.2f}")
            st.write(f"**Mentions (News & Social)**: {len(sentences_to_analyze)} data points analyzed")
        else:
            st.info("Select a document to view details.")

