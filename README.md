# 🔭 OmniScope

## Overview
**OmniScope** is an advanced market intelligence and signal aggregation platform designed to provide institutional-grade insights to individual investors. By unifying data from regulatory filings, real-time market feeds, social media sentiment, and alternative data sources, OmniScope enables users to see the complete picture of any public company.

**"Thinking, deciding, and doing are separate concerns."** — OmniScope aggregates the "thinking" material so you can do the "deciding."

## 🚀 Key Features

### 1. Unified Search
*   **Intelligent Ticker Recognition**: Instantly pulls data for any US-listed company.
*   **Keyword Discovery**: Search across news and filings for specific topics (e.g., "AI", "Lithium").

### 2. Multi-Dimensional Data
*   **Real-Time Market Data**: Live prices, dynamic charts (1d to 5y+), and analyst recommendations.
*   **Regulatory Intelligence**: Direct access to SEC Filings (10-K, 10-Q, 8-K) and **13-F Institutional Holdings**.
*   **Fundamental Analysis**: Deep dive into Balance Sheets, Income Statements, and Cash Flows.

### 3. Alternative Signals (The "Alpha")
*   **Social Listening**: Real-time scraping of **Reddit** (`r/stocks`, `r/investing`, `r/wallstreetbets`) to gauge retail sentiment.
*   **Insider Trading**: Visualization of C-suite buy/sell activity (Form 4).
*   **Short Interest Analysis**: Squeeze potential indicators using Share Float, Short Ratio, and Short % of Float.
*   **Options Sentiment**: Put/Call Ratios to track "Smart Money" positioning.
*   **ESG Risk Scores**: Environmental, Social, and Governance risk ratings.

### 4. Advanced Analytics
*   **Sentiment Engine**: Powered by `TextBlob`, OmniScope analyzes thousands of headlines and social comments to generate real-time **Polarity** (Positive/Negative) and **Subjectivity** scores.
*   **Management Profiling**: Profiles of key executives, including compensation and tenure.

### 5. Congressional Trading
-   **Data Source**: Capitol Trades (Scraped).
-   **Features**:
    -   Tracks stock trades by US Senators and Representatives.
    -   Displays filtered transaction details (Purchases, Sales, Amounts).
    -   Data is scraped live from recent transactions (last ~500 trades).
    
    *Note*: Displays matches from the most recent global activity. Older trades for a specific ticker may not appear if they are not in the recent batch.

## 🛠️ Technology Stack
*   **Frontend**: Streamlit (Python)
*   **Data Visualization**: Plotly Express
*   **NLP/Sentiment**: TextBlob
*   **Data Providers**: 
    *   `yfinance` (Market Data, Options, ESG, Holders)
    *   `sec-edgar-api` / Direct HTTP (Regulatory Filings)
    *   `feedparser` (News)
    *   `requests` (Social Media & Web Scraping)

## 📦 Installation & Usage

1.  **Clone the Repository**
    ```bash
    git clone <repository_url>
    cd projects/alpha_stream
    ```

2.  **Install Dependencies**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Run the Application**
    ```bash
    streamlit run app.py
    ```

## ⚠️ Disclaimer
OmniScope is a research tool for informational purposes only. It does not constitute financial advice. All data is sourced from public APIs and may be subject to delays or inaccuracies. Always verify data before making investment decisions.

## 📄 License
**PolyForm Noncommercial License 1.0.0**

This software is free to use for personal, educational, and research purposes. 
**Commercial use is strictly prohibited** without a separate agreement. You cannot sell this software, sell services based on this software, or use it for commercial gain.
