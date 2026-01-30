# Gene Correlation Explorer (Web Edition)

A web-based application for analyzing gene correlations from DepMap CRISPR screen data. This is the JavaScript/web version of the original R Shiny "Gene Correlation Explorer" app.

## Features

- **Correlation Analysis**: Calculate pairwise correlations between genes based on CRISPR dependency data
- **Network Visualization**: Interactive network graph showing gene relationships with customizable styling
- **Cluster Detection**: Automatic clustering of correlated genes
- **Scatter Plot Inspection**: Detailed scatter plots for each gene pair with:
  - Hotspot mutation overlay (color-coded by mutation status: 0/1/2)
  - 3-panel stratification by mutation level
  - Statistical comparison table by cancer type
  - Customizable axis ranges and aspect ratio
  - Cell line highlighting and search
- **Data Export**: Download results as PNG, SVG, or CSV

## How to Run

### Option 1: Simple Local Server (Recommended)

The app requires a local web server due to browser security restrictions when loading data files.

**Using Python (most common):**
```bash
# Python 3
cd "correlation app green listed"
python -m http.server 8000

# Then open http://localhost:8000 in your browser
```

**Using Node.js:**
```bash
# Install serve globally (one time)
npm install -g serve

# Run the server
cd "correlation app green listed"
serve .

# Then open the URL shown in terminal
```

**Using PHP:**
```bash
cd "correlation app green listed"
php -S localhost:8000

# Then open http://localhost:8000 in your browser
```

### Option 2: VS Code Live Server

1. Install the "Live Server" extension in VS Code
2. Open the project folder in VS Code
3. Right-click on `index.html` and select "Open with Live Server"

### Option 3: Host on GitHub Pages

1. Push to GitHub
2. Go to repository Settings > Pages
3. Select "main" branch and save
4. Access via `https://yourusername.github.io/correlation-web-app/`

## Data Sources

The app uses pre-processed data from the DepMap (Cancer Dependency Map) project:

- **Gene Effect Data**: CRISPR knockout dependency scores (`CRISPRGeneEffect_25Q3.csv`)
- **Cell Line Metadata**: Cancer type and lineage information (`Model.csv`)
- **Hotspot Mutations**: Somatic mutation data (`OmicsSomaticMutationsMatrixHotspot_25Q3.csv`)

Data is stored in compressed format in the `web_data/` folder.

## File Structure

```
correlation-web-app/
├── index.html          # Main HTML file
├── app.js              # Application logic
├── README.md           # This file
└── web_data/           # Pre-processed data files
    ├── gene_effect_*.json.gz    # Gene effect data (chunked)
    ├── cell_lines.json          # Cell line metadata
    ├── mutations.json           # Hotspot mutation data
    └── gene_synonyms.json       # Gene name synonyms
```

## Usage

1. **Input Genes**: Enter gene symbols (one per line) or upload a CSV/TSV file with optional LFC/FDR statistics
2. **Set Parameters**:
   - Choose Analysis mode (within list) or Design mode (find correlated genes)
   - Set correlation cutoff (default: 0.5)
   - Set minimum cell lines (default: 10)
   - Set minimum slope (default: 0.1)
3. **Run Analysis**: Click "Run Analysis" to compute correlations
4. **Explore Results**:
   - **Network tab**: Interactive visualization of gene relationships
   - **Correlations tab**: Table of all gene pairs with correlation values
   - **Clusters tab**: Gene assignments to correlation clusters
   - **Summary tab**: Text summary of analysis results

## Inspect Modal Features

Click "Inspect" on any correlation to see detailed scatter plot:

- **Hotspot Mutation Modes**:
  - Color overlay: Points colored by mutation count (gray=WT, blue=1, red=2+)
  - 3-panel: Separate plots for each mutation level
  - Compare table: Statistical comparison by cancer type

- **Aspect Ratio**: Graphs are square by default (1:1 ratio), adjustable via controls

- **Export**: Download as PNG, SVG, or CSV with full mutation data

## Technical Details

- **Frontend**: Vanilla JavaScript, HTML5, CSS3
- **Visualization**: Plotly.js (scatter plots), vis-network (network graphs)
- **Data Processing**: Client-side JavaScript with pako for gzip decompression
- **Statistics**: Pearson correlation, linear regression, Fisher z-transformation

## Credits

- **Data Source**: [DepMap Portal](https://depmap.org/portal/) - Broad Institute
- **Original R Shiny App**: [Green Listed](https://greenlisted.cmm.se)
- **Lab**: [Wermeling Lab](https://wermelinglab.com), Karolinska Institutet

## Development

For AI-assisted development (Claude, GPT, etc.), see **[CLAUDE.md](CLAUDE.md)** which contains:
- Architecture overview and data flow
- Key data structures and methods
- Common modification patterns
- Testing instructions

## License

MIT License - Free to use and modify.
