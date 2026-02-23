# Correlate - CRISPR Screen Correlation Analysis Tool

A web-based application for analyzing gene correlations from DepMap CRISPR screen data. Part of the [Green Listed](https://greenlisted.cmm.se) tool family.

**Live App:** [https://fredrikwermeling.github.io/correlation-web-app-/](https://fredrikwermeling.github.io/correlation-web-app-/)

## Features

### Analysis Modes

1. **Analysis Mode** - Calculate pairwise correlations within your gene list
2. **Design Mode** - Find genes correlated with your input genes (expands the network)
   - Optional expanded network: discover correlations between newly found genes
3. **Mutation Analysis Mode** - Compare gene effects between wild-type and mutant cells for a selected hotspot mutation (Welch's t-test)
4. **Synonym/Ortholog Lookup Mode** - Map gene symbols to DepMap identifiers using risk-tiered synonym and mouse-to-human ortholog matching

### Network Visualization

- Interactive network graph (vis-network) with draggable nodes
- Customizable node size, font size, and edge width
- Color nodes by gene effect (signed or absolute)
- Color nodes by uploaded statistics (LFC, FDR)
- Show gene effect values and SD in node labels
- Gene info tooltips on hover (via MyGene.info)
- Click nodes to hide them, with restore option
- Double-click nodes to open gene effect distribution
- Export as PNG or SVG with legend
- Physics-based or manual layout with auto-arrange

### Scatter Plot Inspection

- Detailed scatter plots for each gene pair via Inspect button
- Hotspot mutation overlay (color-coded by mutation level: 0/1/2)
- Three-panel stratification by mutation level
- By-tissue analysis with correlation comparison table (Fisher z-transformation)
- Cancer type filtering
- Cell line highlighting and search
- Customizable axis ranges and aspect ratio (including square mode)
- Export as PNG, SVG, or CSV

### Gene Effect Distribution

- Box plots of gene effect by cancer type (tissue)
- Box plots of gene effect by hotspot mutation status
- Tissue filter dropdown for focused analysis
- Sortable statistics tables
- Export as PNG, SVG, or CSV

### Mutation Analysis

- Welch's t-test comparing gene effects between WT and mutant cells
- Three mutation levels (0 = WT, 1 = one hotspot mutation, 2 = two or more)
- Filter by p-value threshold
- Gene effect distribution charts per mutation group
- **Compare by Cancer Type** - Δ GE for top significant genes across each tissue, with pinned "All" row; clickable values open inspect with tissue filter applied
- **Compare by Hotspot** - Cross-compare other hotspot mutations against the top significant genes, with pinned reference row
- Sortable compare table headers with pinned row support

### Gene Input

- Paste gene symbols directly (one per line, comma or space separated)
- Upload gene list from CSV/TSV file with optional LFC and FDR columns
- Test gene set for quick exploration
- Automatic gene synonym and ortholog resolution (low/mid risk tiers)

### Data Export

- Download network as PNG/SVG with legend
- Download scatter plots as PNG/SVG/CSV
- Download gene effect charts as PNG/SVG/CSV
- Download correlation and cluster tables as CSV
- Download all results as ZIP archive

## Data Sources

The app uses data from the [DepMap (Cancer Dependency Map)](https://depmap.org/portal/) project, release 25Q3:

- **CRISPRGeneEffect** - CRISPR knockout dependency scores
- **Model** - Cell line metadata (cancer type, lineage, subtype)
- **OmicsSomaticMutationsMatrixHotspot** - Somatic hotspot mutation levels (0/1/2)

If you use this tool, please [acknowledge DepMap](https://depmap.org/portal/data_page/?tab=overview#how-to-cite) in your publications.

## How to Use

1. **Wait for data to load** - The app loads ~38MB of gene effect data on startup

2. **Input Genes**:
   - Paste gene symbols (one per line) in the text area, OR
   - Upload a CSV/TSV file with optional LFC/FDR statistics columns

3. **Set Parameters**:
   - Choose analysis mode
   - Set correlation cutoff (default: 0.5)
   - Set minimum cell lines (default: 50)
   - Set minimum slope (default: 0.1)
   - Optionally filter by lineage/subtype or hotspot mutation

4. **Run Analysis**: Click "Run Analysis"

5. **Explore Results**:
   - **Network tab**: Interactive visualization with export options
   - **Correlations tab**: Sortable table with Inspect button for each pair
   - **Clusters tab**: Gene cluster assignments with mean effect and SD
   - **Mutation Analysis tab**: Differential gene effect results with compare tables (mutation mode)
   - **Synonyms/Orthologs tab**: Mapped gene symbols (synonym/ortholog mode)
   - **Summary tab**: Text summary of analysis parameters and results

## File Structure

```
correlation-web-app/
├── index.html                      # Main HTML file with UI and styling
├── app.js                          # Application logic
├── README.md                       # This file
└── web_data/
    ├── geneEffects.bin.gz          # Gene effect matrix (binary, gzip compressed)
    ├── metadata.json               # Gene names and cell line IDs
    ├── cellLineMetadata.json       # Lineage and subtype info
    ├── mutations.json              # Hotspot mutation data (levels 0/1/2)
    ├── synonyms.json               # Gene synonyms (low/mid risk tiers)
    ├── orthologs.json              # Mouse-to-human ortholog mapping
    └── correlate_logo.png          # App logo
```

## Technical Details

- **Frontend**: Vanilla JavaScript, HTML5, CSS3 (no build system, no backend)
- **Visualization**: [Plotly.js](https://plotly.com/javascript/) (scatter plots, box plots), [vis-network](https://visjs.github.io/vis-network/docs/network/) (network graphs)
- **Data Processing**: Client-side with [pako](https://github.com/nodeca/pako) for gzip decompression
- **Statistics**: Pearson correlation, linear regression, Welch's t-test, Fisher z-transformation

## Credits

- **Data**: [DepMap Portal](https://depmap.org/portal/) - Broad Institute
- **Tool Family**: [Green Listed](https://greenlisted.cmm.se)
- **Development**: [Wermeling Lab](https://wermelinglab.com), Karolinska Institutet

## Funding

This project was financially supported by the Swedish Cancer Society, the Swedish Research Council, and Karolinska Institutet.

## Related Tools

- [CoExpress](https://fredrikwermeling.github.io/CoExpress/) - Same analysis engine for gene expression (RNA-seq) data
- [Green Listed](https://greenlisted.cmm.se) - CRISPR screen sgRNA lookup tool
- [YouTube Channel](https://www.youtube.com/@fredrikwermeling1330) - CRISPR-related tutorials

## License

MIT License - Free to use and modify.
