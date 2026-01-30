# Claude/AI Development Guide

This document provides context for Claude or other AI assistants to understand and modify this project.

## Project Overview

**Gene Correlation Explorer (Web Edition)** is a client-side JavaScript web application that analyzes gene dependency correlations from DepMap CRISPR screening data. It's a web port of the original R Shiny app.

## Architecture

### Single-Page Application
- **No build system** - vanilla JavaScript, runs directly in browser
- **No backend** - all processing happens client-side
- **Data pre-processed** - large datasets are pre-chunked and compressed

### Key Files

| File | Purpose |
|------|---------|
| `index.html` | UI structure, CSS styling, modal definitions |
| `app.js` | All application logic in a single `GeneCorrelationApp` class |
| `web_data/` | Pre-processed data files |

## Data Flow

```
1. Page Load
   └── loadData() fetches compressed JSON files from web_data/
       ├── geneEffects.bin.gz (chunked gene dependency scores)
       ├── cellLineMetadata.json (cancer types, lineages)
       ├── mutations.json (hotspot mutation levels 0/1/2)
       └── orthologs.json (gene name synonyms)

2. User Input
   └── Genes entered via textarea or CSV upload
       └── validateGenes() checks against loaded gene list

3. Analysis
   └── runAnalysis() calculates pairwise correlations
       ├── pearsonWithSlope() for statistics
       ├── Clustering via correlation matrix
       └── Results stored in this.correlationResults, this.clusterResults

4. Visualization
   ├── buildNetwork() creates vis-network graph
   ├── renderSinglePanelPlot() / renderThreePanelPlot() for scatter plots
   └── Export functions for PNG/SVG/CSV
```

## Key Data Structures

### Gene Effect Data
```javascript
this.geneEffectData = {
    genes: ['TP53', 'BRCA1', ...],           // Gene symbols
    cellLines: ['ACH-000001', ...],           // Cell line IDs
    data: Float32Array                        // Flattened matrix [genes x cellLines]
}
// Access: data[geneIndex * numCellLines + cellLineIndex]
```

### Cell Line Metadata
```javascript
this.cellLineMetadata = {
    'ACH-000001': {
        name: 'NIHOVCAR3',
        lineage: 'Ovary/Fallopian Tube',
        subtype: 'Ovarian Epithelial Tumor'
    }
}
```

### Mutation Data
```javascript
this.mutationData = {
    genes: ['TP53', 'KRAS', ...],
    geneData: {
        'TP53': {
            mutations: { 'ACH-000001': 2, 'ACH-000002': 0, ... },
            counts: { 0: 800, 1: 150, 2: 50 }
        }
    }
}
// mutationLevel: 0 = WT, 1 = 1 hotspot mutation, 2 = 2+ mutations
```

### Correlation Results
```javascript
this.correlationResults = [
    { gene1: 'TP53', gene2: 'MDM2', correlation: 0.85, slope: 0.72, n: 1000, cluster: 1 },
    ...
]
```

## Important Methods

### Statistics
- `pearsonWithSlope(x, y)` - Returns `{correlation, slope, n}`
- `normalCDF(x)` - Standard normal CDF for p-values
- Fisher z-transformation used in `renderCompareTable()` for comparing correlations

### Visualization
- **Network**: vis-network library, nodes/edges defined in `buildNetwork()`
- **Scatter plots**: Plotly.js, three render modes:
  - `renderSinglePanelPlot()` - Color overlay by mutation
  - `renderThreePanelPlot()` - Side-by-side panels
  - `renderCompareTable()` - Statistical comparison HTML table

### Export
- `downloadNetworkPNG/SVG()` - Network with legend
- `downloadScatterPNG/SVG/CSV()` - Scatter plot exports
- `downloadAllData()` - ZIP file with all results

## UI Components

### Inspect Modal (`#inspectModal`)
Controls for scatter plot inspection:
- Axis ranges: `#scatterXmin/Xmax/Ymin/Ymax`
- Aspect ratio: `#aspectRatioX`, `#aspectRatioY`, `#forceSquare`
- Cancer filter: `#scatterCancerFilter`
- Cell line search: `#scatterCellSearch`
- Hotspot overlay: `#hotspotGene`, `#hotspotMode`

### Hotspot Modes
```javascript
// #hotspotMode options:
'none'          // No mutation overlay
'color'         // Single panel, points colored by mutation level
'three_panel'   // 3 separate panels (0/1/2 mutations)
'compare_table' // Statistical table by cancer type
```

## Common Modifications

### Adding a new hotspot mode
1. Add option to `#hotspotMode` select in `index.html`
2. Add case in `updateInspectPlot()` method
3. Create new render function if needed

### Modifying scatter plot layout
- Single panel: `renderSinglePanelPlot()` around line 2095
- Three panel: `renderThreePanelPlot()` around line 2283
- Plotly layout options control margins, titles, axes

### Adding new statistics
1. Calculate in appropriate render function
2. Add to annotations array for display on plot
3. Include in CSV export if relevant

### Network styling
- Node colors: `buildNetwork()` method, look for `color:` properties
- Edge colors: positive (blue #3182ce) / negative (red #e53e3e)
- Legend: `downloadNetworkPNG/SVG()` methods

## External Dependencies

```html
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>  <!-- Scatter plots -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/pako/2.1.0/pako.min.js"></script>  <!-- Gzip decompression -->
<script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>  <!-- Network graph -->
```

## Data Regeneration

If source data needs updating:

1. **Gene Effects**: Process `CRISPRGeneEffect_*.csv` into chunked binary format
2. **Cell Lines**: Extract from `Model.csv` into `cellLineMetadata.json`
3. **Mutations**: Process `OmicsSomaticMutationsMatrixHotspot_*.csv` into `mutations.json`
   - Count actual values (0/1/2) per cell line, not just binary

## Known Quirks

- Aspect ratio uses Plotly's `scaleanchor`/`scaleratio` - can behave unexpectedly with certain ranges
- Three-panel plot annotations positioned at `y: 1.02` to avoid title overlap
- Network legend is manually drawn on canvas for PNG export (not from vis-network)
- Cell line clicking stores in `this.clickedCells` Set, persists until cleared

## Testing the App

1. Run local server: `python -m http.server 8000`
2. Open http://localhost:8000
3. Click "Test Genes" for sample gene list
4. Click "Run Analysis"
5. Test scatter plot: Correlations tab > Inspect button
6. Test mutation overlay: Select gene from hotspot dropdown

## Original R Shiny App

The original `app.R` file is in the project folder (gitignored) for reference. The web version replicates most functionality but with some UI differences.
