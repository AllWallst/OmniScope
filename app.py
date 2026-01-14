import streamlit as st
import pandas as pd
import sys
import os
import plotly.express as px

# Add project root to path to import local modules
# sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

try:
    from execution import data_sources
    import importlib
    importlib.reload(data_sources)
except ImportError:
    # Fallback for local dev if run from different CWD
    try:
        from projects.alpha_stream.execution import data_sources
        import importlib
        importlib.reload(data_sources)
    except ImportError:
        pass

st.set_page_config(layout="wide", page_title="OmniScope", page_icon="🔭")

# Load Custom CSS
def load_css():
    st.markdown("""
        <style>
        /* Modern Apple-style Dark Theme */
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&display=swap');
        
        .stApp {
            background-color: #000000; /* Pure Black OLED friendly */
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        }
        
        /* Sidebar styling - Translucent and blurry */
        .css-1d391kg {
            background-color: rgba(28, 28, 30, 0.95);
            backdrop-filter: blur(20px);
            border-right: 1px solid rgba(255, 255, 255, 0.1);
        }
        
        /* Card styling - Apple Card Style */
        .doc-card {
            background-color: #1C1C1E; /* Dark Gray */
            color: #F5F5F7;
            border-radius: 18px; /* Apple rounded corners */
            padding: 24px;
            margin-bottom: 20px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.3);
            border: 1px solid rgba(255,255,255,0.05);
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }
        .doc-card:hover {
            transform: translateY(-2px);
            box-shadow: 0 8px 30px rgba(0,0,0,0.5);
            border: 1px solid rgba(255,255,255,0.1);
        }
        .doc-meta {
            font-size: 0.75em;
            color: #8E8E93; /* Apple System Gray */
            margin-bottom: 8px;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            font-weight: 500;
        }
        .doc-title {
            font-size: 1.25em;
            font-weight: 600;
            color: #2997FF; /* Apple System Blue */
            text-decoration: none;
            display: block;
            margin-bottom: 12px;
        }
        .doc-title:hover {
            color: #64D2FF; /* Lighter blue on hover */
            text-decoration: none;
        }
        .doc-card p {
            color: #D1D1D6; /* Light gray text */
            font-size: 0.95em;
            line-height: 1.6;
            font-weight: 400;
        }
        
        /* Segmented Control / Tabs Styling */
        .stTabs [data-baseweb="tab-list"] {
            gap: 12px;
            background-color: #1C1C1E;
            padding: 8px;
            border-radius: 16px;
            margin-bottom: 20px;
        }
        .stTabs [data-baseweb="tab"] {
            height: 40px;
            white-space: pre-wrap;
            background-color: transparent;
            border-radius: 10px;
            gap: 1px;
            padding: 8px 16px;
            color: #8E8E93;
            font-weight: 500;
            border: none;
            transition: all 0.2s ease;
        }
        .stTabs [aria-selected="true"] {
            background-color: #636366; /* System Gray 2 */
            color: #FFFFFF;
            font-weight: 600;
            box-shadow: 0 2px 8px rgba(0,0,0,0.2);
        }
        
        /* Input Fields */
        .stTextInput > div > div > input {
            background-color: #1C1C1E;
            color: white;
            border-radius: 12px;
            border: 1px solid rgba(255,255,255,0.1);
            padding: 10px 15px;
        }
        .stSelectbox > div > div > div {
             background-color: #1C1C1E;
             color: white;
             border-radius: 12px;
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
        # 1. Market Data (Smart Search)
        # Check if query looks like a ticker (short, 1 word) or needs lookup
        ticker = query.upper()
        
        # Heuristic: If it contains spaces OR is longer than 3 chars, it might be a company name.
        # "Ford" -> "F", "Apple" -> "AAPL", "Microsoft" -> "MSFT"
        if len(query) > 3 or " " in query:
             with st.spinner(f"Verifying ticker for '{query}'..."):
                found_ticker = data_sources.lookup_ticker(query)
                if found_ticker:
                    # Only switch and notify if it's different/better
                    if found_ticker != ticker:
                        ticker = found_ticker
                        st.toast(f"Resolved '{query}' to {ticker}")
                    else:
                        # It was already a correct ticker (e.g. "MSFT" -> "MSFT")
                        ticker = found_ticker
        
        if ticker.isalpha() and len(ticker) < 10: # Reasonable ticker length
            with st.spinner(f"Fetching Market Data for {ticker}..."):
                market_data = data_sources.fetch_market_data(ticker, period=period)
                if market_data:
                    info = market_data['info']
                    current_price = info.get('currentPrice', 'N/A')
                    previous_close = info.get('previousClose', 0)
                    
                    # Calculate Delta
                    delta = None
                    if isinstance(current_price, (int, float)) and isinstance(previous_close, (int, float)):
                        delta = current_price - previous_close
                    
                    st.markdown(f"### {info.get('shortName', ticker)} ({ticker})")
                    st.metric(
                        "Current Price", 
                        f"${current_price}", 
                        f"{delta:+.2f}" if delta is not None else None
                    )
                    st.caption(f"Recommendation: {info.get('recommendationKey', '').replace('_', ' ').title()}")
                    
                    # Enhanced Charting with Plotly
                    hist_data = market_data['history'].reset_index()
                    # Ensure Date column is formatted for Plotly
                    
                    fig = px.line(hist_data, x="Date", y="Close", title=f"{ticker} Price History")
                    fig.update_layout(
                        xaxis_title="", 
                        yaxis_title="Price (USD)",
                        margin=dict(l=0, r=0, t=30, b=0),
                        hovermode="x unified"
                    )
                    # Force y-axis to not start at 0 if the variation is small relative to price
                    # Plotly default is usually good, but we can enforce it:
                    fig.update_yaxes(autorange=True, fixedrange=False) 
                    
                    # Note: "use_container_width" is deprecated in favor of "width" in some versions of Streamlit/Plotly integration.
                    # We'll try to use the new "width" parameter if supported, or fall back explicitly.
                    # The warning suggests: use_container_width=True -> width="stretch" (or loop in the future).
                    try:
                       st.plotly_chart(fig, use_container_width=True) 
                       # If this still warns, it's safer to keep for now than breaking older versions with 'width'.
                       # User reported: "use_container_width will be removed... use width='stretch'".
                       # Let's try passing the config dict if possible, or just ignore for now as it's a warning.
                       # actually, let's use the explicit kwargs if we can.
                    except:
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
            tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs(["Feed", "Financials", "Holders", "Analysis", "Social", "Insider", "Short Data", "Management", "Congress"])
            
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
                    major_holders = market_data.get('major_holders')
                    if major_holders is not None and not major_holders.empty:
                        try:
                            # 1. Standardize Columns
                            # yfinance often returns Description as Index and 'Value' as text/float
                            major_holders.reset_index(inplace=True)
                            major_holders.columns = ["Description", "Value"]
                            
                            # 2. Process Rows
                            clean_rows = []
                            
                            # Mapping schema: Raw Key -> (Display Text, IsPercent)
                            # We check if the current description matches keys or partial text
                            formatting_map = {
                                "insidersPercentHeld": ("% of Shares Held by All Insider", True),
                                "institutionsPercentHeld": ("% of Shares Held by Institutions", True),
                                "institutionsFloatPercentHeld": ("% of Float Held by Institutions", True),
                                "institutionsCount": ("Number of Institutions Holding Shares", False)
                            }

                            for index, row in major_holders.iterrows():
                                desc = str(row.get('Description', ''))
                                val = row.get('Value', 0)
                                
                                # Default to current values
                                final_desc = desc
                                final_val = val
                                is_percent = True # Default behavior for this table unless specified
                                
                                # Check if we have a specific map for this row
                                if desc in formatting_map:
                                    final_desc, is_percent = formatting_map[desc]
                                elif "Count" in desc or "Number" in desc:
                                    is_percent = False
                                
                                # Apply Formatting
                                try:
                                    if is_percent:
                                        if isinstance(val, (float, int)):
                                            final_val = f"{val:.2%}"
                                    else:
                                        # Number format (comma separated)
                                        if isinstance(val, (float, int)):
                                            final_val = f"{int(val):,}"
                                except:
                                    pass # Keep original if format fails
                                    
                                clean_rows.append({"Value": final_val, "Description": final_desc})
                            
                            st.table(pd.DataFrame(clean_rows))
                        except Exception as e:
                            st.error(f"Error formatting holders: {e}")
                            st.dataframe(major_holders)
                    else:
                        st.info("No Major Holders data.")

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
                    
                    col_h1, col_h2 = st.columns([3, 1])
                    with col_h1:
                        st.subheader(f"${query.upper()} Insider Trading Activity")
                    with col_h2:
                        insider_period = st.selectbox("Graph Period", ["1mo", "6mo", "1y", "2y", "5y", "10y", "max"], index=4, key="insider_period")
                    
                    with st.spinner("Fetching Insider Data..."):
                        # Fetch insider data
                        insiders = data_sources.fetch_insider_transactions(query)
                        
                        # Fetch history specifically for this graph period to ensure alignment
                        insider_history = data_sources.fetch_market_data(query, period=insider_period)['history']
                        
                        if not insiders.empty and insider_history is not None:
                            # --- Data Processing for Graph ---
                            try:
                                # 1. Prepare Price History
                                history = insider_history.reset_index()
                                # Ensure timezone naive for comparison
                                if history['Date'].dt.tz is not None:
                                    history['Date'] = history['Date'].dt.tz_localize(None)
                                
                                # 2. Prepare Insider Data
                                df_insider = insiders.copy()
                                if 'Start Date' in df_insider.columns:
                                    # Ensure timezone naive
                                    df_insider['Date'] = pd.to_datetime(df_insider['Start Date'])
                                    if df_insider['Date'].dt.tz is not None:
                                        df_insider['Date'] = df_insider['Date'].dt.tz_localize(None)
                                    
                                    # Filter for relevant transaction types
                                    def get_trade_type(text):
                                        text = str(text).lower()
                                        if 'sale' in text or 'm - ' in text: 
                                            return 'Sale'
                                        elif 'purchase' in text or 'buy' in text:
                                            return 'Purchase'
                                        return 'Other'

                                    col_to_check = 'Text' if 'Text' in df_insider.columns else 'Transaction'
                                    df_insider['Type'] = df_insider[col_to_check].apply(get_trade_type)
                                    
                                    # Filter only Purchases and Sales
                                    trades = df_insider[df_insider['Type'].isin(['Purchase', 'Sale'])].copy()
                                    
                                    # Filter trades to be within the history range
                                    min_hist_date = history['Date'].min()
                                    max_hist_date = history['Date'].max()
                                    trades = trades[(trades['Date'] >= min_hist_date) & (trades['Date'] <= max_hist_date)]
                                    
                                    if not trades.empty:
                                        # Match trades to stock price
                                        history_lookup = history.set_index('Date')['Close'].sort_index()
                                        
                                        trade_prices = []
                                        for d in trades['Date']:
                                            try:
                                                # Use get_indexer to find nearest position
                                                idx = history_lookup.index.get_indexer([d], method='nearest')[0]
                                                price = history_lookup.iloc[idx]
                                                trade_prices.append(price)
                                            except:
                                                trade_prices.append(None)
                                        
                                        trades['Price'] = trade_prices
                                        trades = trades.dropna(subset=['Price'])

                                        # --- Plotting ---
                                        import plotly.graph_objects as go
                                        import numpy as np

                                        fig = go.Figure()

                                        # 1. Stock Price Line
                                        fig.add_trace(go.Scatter(
                                            x=history['Date'], 
                                            y=history['Close'],
                                            mode='lines',
                                            name='Price',
                                            line=dict(color='#8E8E93', width=1.5),
                                            hoverinfo='y+x'
                                        ))
                                        
                                        # 2. Insider Sales (Red)
                                        sales = trades[trades['Type'] == 'Sale']
                                        if not sales.empty:
                                            # Scale bubble size
                                            sizes = np.log1p(sales['Shares'].abs().astype(float))
                                            if sizes.max() != sizes.min():
                                                sizes = 5 + ((sizes - sizes.min()) / (sizes.max() - sizes.min())) * 15
                                            else:
                                                sizes = 10
                                            
                                            fig.add_trace(go.Scatter(
                                                x=sales['Date'],
                                                y=sales['Price'],
                                                mode='markers',
                                                name='Sales',
                                                marker=dict(
                                                    color='#FF453A',
                                                    size=sizes,
                                                    opacity=0.8,
                                                    line=dict(width=1, color='white')
                                                ),
                                                # Robust hovertemplate for nicer tooltips
                                                hovertemplate="<b>%{text}</b><br>Price: $%{y:.2f}<extra></extra>",
                                                text=sales.apply(lambda x: f"{x['Insider']}<br>Sold {x['Shares']:,}<br>{x['Date'].date()}", axis=1),
                                            ))

                                        # 3. Insider Purchases (Green)
                                        purchases = trades[trades['Type'] == 'Purchase']
                                        if not purchases.empty:
                                            sizes = np.log1p(purchases['Shares'].abs().astype(float))
                                            if sizes.max() != sizes.min():
                                                sizes = 5 + ((sizes - sizes.min()) / (sizes.max() - sizes.min())) * 15
                                            else:
                                                sizes = 10
                                                
                                            fig.add_trace(go.Scatter(
                                                x=purchases['Date'],
                                                y=purchases['Price'],
                                                mode='markers',
                                                name='Purchases',
                                                marker=dict(
                                                    color='#30D158',
                                                    size=sizes,
                                                    opacity=0.8,
                                                    line=dict(width=1, color='white')
                                                ),
                                                hovertemplate="<b>%{text}</b><br>Price: $%{y:.2f}<extra></extra>",
                                                text=purchases.apply(lambda x: f"{x['Insider']}<br>Bought {x['Shares']:,}<br>{x['Date'].date()}", axis=1),
                                            ))

                                        fig.update_layout(
                                            xaxis_title="",
                                            yaxis_title="",
                                            template="plotly_dark",
                                            paper_bgcolor='rgba(0,0,0,0)',
                                            plot_bgcolor='rgba(0,0,0,0)',
                                            margin=dict(l=0, r=0, t=10, b=0),
                                            legend=dict(
                                                yanchor="top",
                                                y=0.99,
                                                xanchor="left",
                                                x=0.01,
                                                bgcolor='rgba(0,0,0,0.5)'
                                            ),
                                            height=400,
                                            hovermode="x unified"
                                        )
                                        
                                        st.plotly_chart(fig, use_container_width=True)
                                        
                                        # --- Summary Text ---
                                        total_trades = len(trades)
                                        num_purchases = len(purchases)
                                        num_sales = len(sales)
                                        
                                        # Calculate time delta based on the selection, not just trades found
                                        if insider_period == "1mo": months_text = "1 month"
                                        elif insider_period == "6mo": months_text = "6 months"
                                        elif insider_period == "1y": months_text = "12 months"
                                        elif insider_period == "2y": months_text = "24 months"
                                        elif insider_period == "5y": months_text = "5 years"
                                        elif insider_period == "10y": months_text = "10 years"
                                        else: months_text = "recent history"

                                        summary_html = f"""
                                        <div style="background-color: #1C1C1E; padding: 16px; border-radius: 12px; border: 1px solid rgba(255,255,255,0.1); margin-top: 12px; margin-bottom: 20px;">
                                            <p style="margin: 0; font-size: 1.0em; line-height: 1.5; color: #E5E5EA;">
                                                <strong>{query.upper()} insiders</strong> have traded <strong>${query.upper()}</strong> stock on the open market 
                                                <strong>{total_trades} times</strong> in the past <strong>{months_text}</strong>. 
                                                Of those trades, <span style="color: #6EDB8E; font-weight: 600; background-color: rgba(48, 209, 88, 0.15); padding: 2px 6px; border-radius: 4px;">{num_purchases} have been purchases</span> 
                                                and <span style="color: #FF6961; font-weight: 600; background-color: rgba(255, 69, 58, 0.15); padding: 2px 6px; border-radius: 4px;">{num_sales} have been sales</span>.
                                            </p>
                                        </div>
                                        """
                                        st.markdown(summary_html, unsafe_allow_html=True)
                                        
                                    else:
                                        st.caption(f"No 'Purchase' or 'Sale' transactions found in the past {insider_period}.")
                            except Exception as e:
                                st.error(f"Error visualizing insider data: {e}")
                                
                            with st.expander("View Raw Insider Data"):
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
                                "Year Born": str(officer.get('yearBorn', 'N/A')) # Force string to avoid mixed-type Arrow errors
                            })
                        st.dataframe(pd.DataFrame(clean_officers))
                    else:
                        st.info("No Executive info available.")
                else:
                    st.info("Enter a ticker to see Management.")
            
            with tab9:
                if query.isalpha():
                    st.subheader(f"Congressional Trading on ${query.upper()}")
                    st.markdown("*Data Source: Quiver Quantitative*")
                    
                    with st.spinner("Fetching Congress Data..."):
                        # No API key needed now
                        congress_data = data_sources.fetch_congress_trading(query)
                        
                        if not congress_data.empty:
                            try:
                                st.dataframe(congress_data, use_container_width=True)
                            except:
                                st.dataframe(congress_data) # Fallback
                            
                            # Simple stats
                            dem_trades = len(congress_data[congress_data['Party'] == 'Democrat'])
                            rep_trades = len(congress_data[congress_data['Party'] == 'Republican'])
                            
                            col_c1, col_c2 = st.columns(2)
                            with col_c1:
                                st.metric("Democrat Trades", dem_trades)
                            with col_c2:
                                st.metric("Republican Trades", rep_trades)
                                
                            st.caption("Data sourced from Quiver Quantitative. Displays recent trading activity.")
                        else:
                            st.warning("No congressional trading activity found for this ticker on Quiver Quantitative.")
                else:
                    st.info("Enter a valid ticker symbol to view Congress trading data.")
                
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

