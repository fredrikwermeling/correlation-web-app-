# Correlate - a CRISPR Screen Analysis Tool

A web-based application for analyzing gene correlations from DepMap CRISPR screen data. Part of the [Green Listed](https://greenlisted.cmm.se) tool family.

**Live App:** [https://fredrikwermeling.github.io/correlation-web-app-/](https://fredrikwermeling.github.io/correlation-web-app-/)

## Features

### Analysis Modes

1. **Analysis Mode** - Calculate pairwise correlations within your gene list
2. **Design Mode** - Find genes correlated with your input genes (expands the network)
3. **Mutation Analysis Mode** - Compare gene effects between wild-type and mutant cells for a selected hotspot mutation

### Core Functionality

- **Network Visualization**: Interactive network graph showing gene relationships
  - Customizable node size, font size, and edge width
  - Color nodes by gene effect (signed or absolute)
  - Color nodes by uploaded statistics (LFC, FDR)
  - Show/hide gene effect values in labels
  - Export as PNG or SVG with legend

- **Correlation Analysis**:
  - Pearson correlation with adjustable cutoff
  - Minimum slope filter
  - Minimum cell line threshold
  - Lineage/subtype filtering

- **Scatter Plot Inspection**: Detailed plots for each gene pair
  - Hotspot mutation overlay (color-coded by mutation status: 0/1/2)
  - 3-panel stratification by mutation level
  - Cancer type filtering
  - Cell line highlighting and search
  - Customizable axis ranges and aspect ratio

- **Mutation Analysis**:
  - Welch's t-test comparing gene effects between WT and mutant cells
  - Filter by p-value threshold
  - Gene effect distribution charts

- **Gene Synonym/Ortholog Support**:
  - Automatic lookup of human gene synonyms
  - Mouse-to-human ortholog mapping
  - Risk-tiered matching (low/mid confidence)

### Data Export

- Download network as PNG/SVG with legend
- Download correlation tables as CSV
- Download cluster assignments as CSV
- Download all data as ZIP archive

## Data Sources

The app uses data from the [DepMap (Cancer Dependency Map)](https://depmap.org/portal/) project, release 25Q3:

- **CRISPRGeneEffect.csv** - CRISPR knockout dependency scores
- **Model.csv** - Cell line metadata (cancer type, lineage)
- **OmicsSomaticMutationsMatrixHotspot.csv** - Somatic hotspot mutations

If you use this tool, please [acknowledge DepMap](https://depmap.org/portal/data_page/?tab=overview#how-to-cite) in your publications.

## How to Use

1. **Wait for data to load** - The app loads ~38MB of gene effect data on startup

2. **Input Genes**:
   - Paste gene symbols (one per line) in the text area, OR
   - Upload a CSV/TSV file with optional LFC/FDR statistics

3. **Set Parameters**:
   - Choose analysis mode
   - Set correlation cutoff (default: 0.5)
   - Set minimum cell lines (default: 50)
   - Set minimum slope (default: 0.1)
   - Optionally filter by lineage/subtype or hotspot mutation

4. **Run Analysis**: Click "Run Analysis"

5. **Explore Results**:
   - **Network tab**: Interactive visualization
   - **Correlations tab**: Table with Inspect button for each pair
   - **Clusters tab**: Gene assignments with mean effect statistics
   - **Mutation Analysis tab**: Differential gene effect results (when using Mutation Analysis mode)
   - **Summary tab**: Text summary

## File Structure

```
correlation-web-app/
├── index.html              # Main HTML file
├── app.js                  # Application logic (~4000 lines)
├── README.md               # This file
└── web_data/
    ├── geneEffects.bin.gz      # Gene effect matrix (binary, compressed)
    ├── metadata.json           # Gene names and cell line IDs
    ├── cellLineMetadata.json   # Lineage and subtype info
    ├── mutations.json          # Hotspot mutation data
    ├── synonyms.json           # Gene synonyms (low/mid risk)
    ├── orthologs.json          # Mouse-human orthologs
    └── correlate_logo.png      # App logo
```

## Technical Details

- **Frontend**: Vanilla JavaScript, HTML5, CSS3
- **Visualization**: Plotly.js (scatter plots), vis-network (network graphs)
- **Data Processing**: Client-side with pako for gzip decompression
- **Statistics**: Pearson correlation, linear regression, Welch's t-test, Fisher z-transformation

## Credits

- **Data**: [DepMap Portal](https://depmap.org/portal/) - Broad Institute
- **Tool Family**: [Green Listed](https://greenlisted.cmm.se)
- **Development**: [Wermeling Lab](https://wermelinglab.com), Karolinska Institutet

## Funding

This project was financially supported by the Swedish Cancer Society, the Swedish Research Council, and Karolinska Institutet.

## Related Tools

- [Green Listed](https://greenlisted.cmm.se) - CRISPR screen sgRNA lookup tool
- [YouTube Channel](https://www.youtube.com/@fredrikwermeling1330) - CRISPR-related tutorials

## License

MIT License - Free to use and modify.
