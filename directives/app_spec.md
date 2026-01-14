# AlphaStream App Specification

## Overview
AlphaStream is a market intelligence platform mirroring AlphaSense capabilities using public data sources.

## Data Sources (Real Data)
1.  **Regulatory Filings (SEC)**
    -   Source: SEC EDGAR.
    -   Fetch method: HTTP requests to `sec.gov` with proper User-Agent.
    -   Content: 10-K, 10-Q links and summaries.
2.  **Company Data**
    -   Source: Yahoo Finance (`yfinance`).
    -   Content: Price history, Key statistics, Company info.
3.  **News & Events**
    -   Source: Google News RSS via `feedparser`.
    -   Content: Headlines, links, pub dates.
4.  **Life Sciences**
    -   Source: PubMed (NCBI).
    -   Content: Research paper abstracts, clinical trial mentions (if available via pubmed).
5.  **Transcripts**
    -   Source: YouTube (`youtube-transcript-api`).
    -   Strategy: Search for "[Ticker] Earnings Call" and retrieve auto-captions.

## UI/UX Flow
1.  **Global Search Bar**: Allows searching for Tickers (e.g., "AAPL") or Keywords (e.g., "Generative AI").
2.  **Sidebar Filters**:
    -   Date Range.
    -   Source Selection (Checkbox for each source type).
3.  **Results Feed**:
    -   Unified timeline of documents/items.
    -   Clicking an item opens it in the Document Viewer.
4.  **Document Viewer**:
    -   Right-hand column or expandable section showing full text/summary.
