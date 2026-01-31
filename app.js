/**
 * Gene Correlation Explorer - Web Edition
 * JavaScript implementation matching the R Shiny app functionality
 *
 * Author: Wermeling Lab / Karolinska Institutet
 * Based on: https://github.com/fredrikwermeling/correlation-app
 */

class CorrelationExplorer {
    constructor() {
        // Data storage
        this.metadata = null;
        this.cellLineMetadata = null;
        this.mutations = null;
        this.orthologs = null;
        this.geneEffects = null; // Float32Array [nGenes x nCellLines]
        this.nGenes = 0;
        this.nCellLines = 0;
        this.geneIndex = new Map(); // gene name -> row index
        this.geneNames = []; // gene names array

        // Analysis results
        this.results = null;
        this.network = null;

        // Current inspect state
        this.currentInspect = null;
        this.clickedCells = new Set();

        // Gene statistics (LFC, FDR)
        this.geneStats = null;
        this.statsFileData = null;

        // Synonyms/orthologs used
        this.synonymsUsed = [];

        this.init();
    }

    async init() {
        try {
            await this.loadData();
            this.setupUI();
            this.hideLoading();
        } catch (error) {
            console.error('Initialization error:', error);
            this.updateLoadingText('Error loading data: ' + error.message);
        }
    }

    updateLoadingText(text) {
        document.getElementById('loadingText').textContent = text;
    }

    hideLoading() {
        document.getElementById('loadingOverlay').classList.add('hidden');
    }

    async loadData() {
        this.updateLoadingText('Loading metadata...');

        // Load all JSON files in parallel
        const [metadataRes, cellLineRes, mutationsRes, orthologsRes] = await Promise.all([
            fetch('web_data/metadata.json'),
            fetch('web_data/cellLineMetadata.json'),
            fetch('web_data/mutations.json'),
            fetch('web_data/orthologs.json')
        ]);

        this.metadata = await metadataRes.json();
        this.cellLineMetadata = await cellLineRes.json();
        this.mutations = await mutationsRes.json();
        this.orthologs = await orthologsRes.json();

        this.nGenes = this.metadata.nGenes;
        this.nCellLines = this.metadata.nCellLines;
        this.geneNames = this.metadata.genes;

        // Build gene index
        this.metadata.genes.forEach((gene, idx) => {
            this.geneIndex.set(gene.toUpperCase(), idx);
        });

        // Load binary gene effects
        this.updateLoadingText('Loading gene effect matrix...');
        await this.loadGeneEffects();

        // Update reference status
        document.getElementById('referenceStatus').className = 'status-box status-success';
        document.getElementById('referenceStatus').innerHTML =
            `&#10003; ${this.nGenes.toLocaleString()} genes, ${this.nCellLines.toLocaleString()} cell lines loaded`;

        // Enable run button
        document.getElementById('runAnalysis').disabled = false;

        // Populate lineage filter if available
        this.populateLineageFilter();
    }

    async loadGeneEffects() {
        const response = await fetch('web_data/geneEffects.bin.gz');
        const compressedData = await response.arrayBuffer();

        this.updateLoadingText('Decompressing gene effect data...');

        // Decompress using pako
        const decompressed = pako.inflate(new Uint8Array(compressedData));

        // Convert to Int16Array
        const int16Data = new Int16Array(decompressed.buffer);

        // Convert to Float32Array and scale
        const scaleFactor = this.metadata.scaleFactor;
        const naValue = this.metadata.naValue;

        this.geneEffects = new Float32Array(int16Data.length);
        for (let i = 0; i < int16Data.length; i++) {
            if (int16Data[i] === naValue) {
                this.geneEffects[i] = NaN;
            } else {
                this.geneEffects[i] = int16Data[i] / scaleFactor;
            }
        }
    }

    populateLineageFilter() {
        const lineages = new Set();
        if (this.cellLineMetadata && this.cellLineMetadata.lineage) {
            Object.values(this.cellLineMetadata.lineage).forEach(l => {
                if (l) lineages.add(l);
            });
        }

        if (lineages.size > 0) {
            const select = document.getElementById('lineageFilter');
            select.innerHTML = '<option value="">All lineages</option>';
            Array.from(lineages).sort().forEach(lineage => {
                const option = document.createElement('option');
                option.value = lineage;
                option.textContent = lineage;
                select.appendChild(option);
            });
            document.getElementById('lineageFilterGroup').style.display = 'block';
        }
    }

    getCellLineName(cellLineId) {
        if (this.cellLineMetadata && this.cellLineMetadata.strippedCellLineName) {
            return this.cellLineMetadata.strippedCellLineName[cellLineId] ||
                   this.cellLineMetadata.cellLineName?.[cellLineId] ||
                   cellLineId;
        }
        return cellLineId;
    }

    getCellLineLineage(cellLineId) {
        if (this.cellLineMetadata && this.cellLineMetadata.lineage) {
            return this.cellLineMetadata.lineage[cellLineId] || '';
        }
        return '';
    }

    getGeneData(geneIndex) {
        const start = geneIndex * this.nCellLines;
        return this.geneEffects.subarray(start, start + this.nCellLines);
    }

    setupUI() {
        // Tab switching
        document.querySelectorAll('.nav-link').forEach(tab => {
            tab.addEventListener('click', () => {
                document.querySelectorAll('.nav-link').forEach(t => t.classList.remove('active'));
                document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
                tab.classList.add('active');
                document.getElementById('tab-' + tab.dataset.tab).classList.add('active');
            });
        });

        // Slider value displays
        document.getElementById('correlationCutoff').addEventListener('input', (e) => {
            document.getElementById('cutoffValue').textContent = parseFloat(e.target.value).toFixed(2);
        });

        document.getElementById('minSlope').addEventListener('input', (e) => {
            document.getElementById('slopeValue').textContent = parseFloat(e.target.value).toFixed(2);
        });

        // Gene textarea
        document.getElementById('geneTextarea').addEventListener('input', () => this.updateGeneCount());

        // Clear genes
        document.getElementById('clearGenes').addEventListener('click', () => {
            document.getElementById('geneTextarea').value = '';
            this.updateGeneCount();
        });

        // Load test genes
        document.getElementById('loadTestGenes').addEventListener('click', () => {
            const testGenes = ['TP53', 'BRCA1', 'BRCA2', 'MYC', 'KRAS', 'EGFR', 'PTEN',
                'RB1', 'APC', 'VHL', 'CDKN2A', 'NOTCH1', 'PIK3CA', 'BRAF',
                'ATM', 'ERBB2', 'CDK4', 'MDM2', 'NRAS', 'ARID1A'];
            document.getElementById('geneTextarea').value = testGenes.join('\n');
            this.updateGeneCount();
        });

        // Find synonyms button
        document.getElementById('findSynonyms').addEventListener('click', () => this.findSynonymsForMissingGenes());

        // Input method tabs
        document.querySelectorAll('.input-tab').forEach(tab => {
            tab.addEventListener('click', () => {
                document.querySelectorAll('.input-tab').forEach(t => t.classList.remove('active'));
                document.querySelectorAll('.input-panel').forEach(p => p.classList.remove('active'));
                tab.classList.add('active');
                document.getElementById('input-' + tab.dataset.input).classList.add('active');
            });
        });

        // Stats file upload
        document.getElementById('statsFileInput').addEventListener('change', (e) => this.handleStatsFileUpload(e));
        document.getElementById('loadStatsBtn').addEventListener('click', () => this.loadStatsFromFile());
        document.getElementById('loadTestStats').addEventListener('click', () => this.loadTestGenesWithStats());
        document.getElementById('downloadSampleStats').addEventListener('click', () => this.downloadSampleStatsFile());

        // Run analysis
        document.getElementById('runAnalysis').addEventListener('click', () => this.runAnalysis());

        // Network controls with slider bubble updates
        document.getElementById('netFontSize').addEventListener('input', (e) => {
            document.getElementById('fontSizeBubble').textContent = e.target.value;
            this.updateNetworkStyle();
        });
        document.getElementById('netNodeSize').addEventListener('input', (e) => {
            document.getElementById('nodeSizeBubble').textContent = e.target.value;
            this.updateNetworkStyle();
        });
        document.getElementById('netEdgeWidth').addEventListener('input', (e) => {
            document.getElementById('edgeWidthBubble').textContent = e.target.value;
            this.updateNetworkStyle();
        });
        document.getElementById('fitNetwork').addEventListener('click', () => {
            if (this.network) this.network.fit();
        });
        document.getElementById('showHiddenNodes').addEventListener('click', () => this.showHiddenNodes());
        document.getElementById('showGeneEffect').addEventListener('change', () => this.updateNetworkLabels());
        document.getElementById('downloadNetworkPNG').addEventListener('click', () => this.downloadNetworkPNG());
        document.getElementById('downloadNetworkSVG').addEventListener('click', () => this.downloadNetworkSVG());
        document.getElementById('downloadAllData').addEventListener('click', () => this.downloadAllData());

        // Color by stats controls
        document.getElementById('colorByStats').addEventListener('change', (e) => {
            document.getElementById('colorStatsOptions').style.display = e.target.checked ? 'block' : 'none';
            document.getElementById('legendNodeColor').style.display = e.target.checked ? 'block' : 'none';
            this.updateNetworkColors();
        });
        document.querySelectorAll('input[name="colorStatType"]').forEach(radio => {
            radio.addEventListener('change', () => this.updateNetworkColors());
        });
        document.querySelectorAll('input[name="colorScale"]').forEach(radio => {
            radio.addEventListener('change', () => this.updateNetworkColors());
        });
        // Stats label display (None/LFC/FDR)
        document.querySelectorAll('input[name="statsLabelDisplay"]').forEach(radio => {
            radio.addEventListener('change', () => this.updateNetworkLabelsWithStats());
        });

        // Table search
        document.getElementById('correlationsSearch').addEventListener('input', (e) => {
            this.filterTable('correlationsBody', e.target.value);
        });
        document.getElementById('clustersSearch').addEventListener('input', (e) => {
            this.filterTable('clustersBody', e.target.value);
        });

        // Download buttons
        document.getElementById('downloadCorrelations').addEventListener('click', () => this.downloadCSV('correlations'));
        document.getElementById('downloadClusters').addEventListener('click', () => this.downloadCSV('clusters'));
        document.getElementById('downloadSummary').addEventListener('click', () => this.downloadSummary());

        // Inspect modal
        document.getElementById('closeInspect').addEventListener('click', () => this.closeInspectModal());
        document.getElementById('closeInspectBtn').addEventListener('click', () => this.closeInspectModal());
        document.getElementById('inspectModal').addEventListener('click', (e) => {
            if (e.target.id === 'inspectModal') this.closeInspectModal();
        });

        // Inspect controls
        document.getElementById('resetAxes').addEventListener('click', () => this.resetInspectAxes());
        document.getElementById('clearHighlights').addEventListener('click', () => {
            document.getElementById('scatterCellSearch').value = '';
            this.clickedCells.clear();
            this.updateInspectPlot();
        });

        ['scatterXmin', 'scatterXmax', 'scatterYmin', 'scatterYmax'].forEach(id => {
            document.getElementById(id).addEventListener('change', () => this.updateInspectPlot());
        });

        document.getElementById('scatterCellSearch').addEventListener('input', () => this.updateInspectPlot());
        document.getElementById('scatterCancerFilter').addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('mutationFilterGene').addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('mutationFilterLevel').addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('hotspotGene').addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('hotspotMode').addEventListener('change', () => this.updateInspectPlot());

        document.getElementById('downloadScatterPNG').addEventListener('click', () => this.downloadScatterPNG());
        document.getElementById('downloadScatterSVG').addEventListener('click', () => this.downloadScatterSVG());
        document.getElementById('downloadScatterCSV').addEventListener('click', () => this.downloadScatterCSV());
        document.getElementById('downloadTissuePNG').addEventListener('click', () => this.downloadTissueChartPNG());
        document.getElementById('downloadTissueSVG').addEventListener('click', () => this.downloadTissueChartSVG());
        document.getElementById('downloadTissueCSV').addEventListener('click', () => this.downloadTissueTableCSV());
        document.getElementById('scatterFontSize')?.addEventListener('change', () => this.updateInspectPlot());

        // Aspect ratio control
        document.getElementById('aspectRatio')?.addEventListener('input', (e) => {
            document.getElementById('aspectRatioValue').textContent = parseFloat(e.target.value).toFixed(1);
            this.updateInspectPlot();
        });

        // Table header sorting
        document.querySelectorAll('.data-table th[data-sort]').forEach(th => {
            th.addEventListener('click', () => this.sortTable(th));
        });
    }

    updateGeneCount() {
        const text = document.getElementById('geneTextarea').value;
        const genes = text.split(/\s+/).filter(g => g.trim() !== '');
        const display = document.getElementById('geneCountDisplay');

        if (genes.length > 0) {
            display.innerHTML = `<strong>Genes entered:</strong> <span class="gene-count">${genes.length}</span>`;
            this.validateGenes(genes);
        } else {
            display.innerHTML = '';
            document.getElementById('geneValidationDisplay').innerHTML = '';
        }
    }

    validateGenes(genes) {
        const upperGenes = genes.map(g => g.toUpperCase().trim());
        const found = upperGenes.filter(g => this.geneIndex.has(g));
        const notFound = upperGenes.filter(g => !this.geneIndex.has(g));

        const display = document.getElementById('geneValidationDisplay');
        const synonymBtn = document.getElementById('findSynonyms');

        if (notFound.length === 0) {
            display.innerHTML = `<div class="status-box status-success">&#10003; All ${found.length} genes found in reference data</div>`;
            synonymBtn.style.display = 'none';
        } else {
            display.innerHTML = `<div class="status-box status-warning">
                <strong>${found.length} found</strong>, <strong>${notFound.length} not found</strong><br>
                <span>Not found: ${notFound.slice(0, 10).join(', ')}${notFound.length > 10 ? ` (+${notFound.length - 10} more)` : ''}</span>
            </div>`;
            synonymBtn.style.display = 'block';
            this.genesNotFound = notFound;
        }
    }

    findSynonymsForMissingGenes() {
        if (!this.genesNotFound || this.genesNotFound.length === 0) return;

        const replacements = [];
        const notFound = this.genesNotFound;

        notFound.forEach(gene => {
            // Check ortholog lookup
            const humanGene = this.orthologs?.mouseToHuman?.[gene];
            if (humanGene && this.geneIndex.has(humanGene.toUpperCase())) {
                replacements.push({ original: gene, replacement: humanGene.toUpperCase(), source: 'ortholog' });
            }
        });

        if (replacements.length > 0) {
            // Store synonyms used for summary
            this.synonymsUsed = replacements;

            // Update textarea
            const textarea = document.getElementById('geneTextarea');
            let text = textarea.value;

            replacements.forEach(r => {
                const regex = new RegExp(`\\b${r.original}\\b`, 'gi');
                text = text.replace(regex, r.replacement);
            });

            textarea.value = text;
            this.updateGeneCount();

            const msg = replacements.map(r => `${r.original} → ${r.replacement} [${r.source}]`).join(', ');
            alert(`Replaced ${replacements.length} gene(s):\n${msg}`);
        } else {
            alert('No synonyms or orthologs found for the missing genes');
        }
    }

    handleStatsFileUpload(event) {
        const file = event.target.files[0];
        if (!file) return;

        const reader = new FileReader();
        reader.onload = (e) => {
            const content = e.target.result;
            this.parseStatsFile(content, file.name);
        };
        reader.readAsText(file);
    }

    parseStatsFile(content, filename) {
        // Remove BOM (Byte Order Mark) if present (Excel UTF-8 exports)
        if (content.charCodeAt(0) === 0xFEFF) {
            content = content.slice(1);
        }

        // Normalize line endings (Windows \r\n, old Mac \r, Unix \n)
        content = content.replace(/\r\n/g, '\n').replace(/\r/g, '\n');

        const lines = content.trim().split('\n');
        if (lines.length < 2) {
            alert('File must have a header row and at least one data row');
            return;
        }

        const firstLine = lines[0];

        // Detect delimiter: tab, comma, or semicolon (European Excel)
        let delimiter = ',';
        if (firstLine.includes('\t')) {
            delimiter = '\t';
        } else if (firstLine.includes(';') && !firstLine.includes(',')) {
            delimiter = ';';
        }

        // Parse header
        const headers = firstLine.split(delimiter).map(h => h.trim().replace(/^"|"$/g, ''));
        this.statsFileData = {
            headers: headers,
            rows: []
        };

        // Parse data rows
        for (let i = 1; i < lines.length; i++) {
            const cols = lines[i].split(delimiter).map(c => c.trim().replace(/^"|"$/g, ''));
            if (cols.length === headers.length) {
                this.statsFileData.rows.push(cols);
            }
        }

        // Populate column selectors
        const geneColSelect = document.getElementById('statsGeneCol');
        const lfcColSelect = document.getElementById('statsLfcCol');
        const fdrColSelect = document.getElementById('statsFdrCol');

        geneColSelect.innerHTML = '';
        lfcColSelect.innerHTML = '<option value="">None</option>';
        fdrColSelect.innerHTML = '<option value="">None</option>';

        headers.forEach((h, idx) => {
            geneColSelect.innerHTML += `<option value="${idx}">${h}</option>`;
            lfcColSelect.innerHTML += `<option value="${idx}">${h}</option>`;
            fdrColSelect.innerHTML += `<option value="${idx}">${h}</option>`;
        });

        // Auto-select columns based on common names
        headers.forEach((h, idx) => {
            const hl = h.toLowerCase();
            if (hl.includes('gene') || hl === 'symbol' || hl === 'name') {
                geneColSelect.value = idx;
            }
            if (hl.includes('lfc') || hl.includes('log') || hl.includes('fold')) {
                lfcColSelect.value = idx;
            }
            if (hl.includes('fdr') || hl.includes('padj') || hl.includes('pval') || hl.includes('p.value')) {
                fdrColSelect.value = idx;
            }
        });

        document.getElementById('statsColumnSelect').style.display = 'block';
    }

    loadStatsFromFile() {
        if (!this.statsFileData) return;

        const geneColIdx = parseInt(document.getElementById('statsGeneCol').value);
        const lfcColIdx = document.getElementById('statsLfcCol').value;
        const fdrColIdx = document.getElementById('statsFdrCol').value;

        const genes = [];
        this.geneStats = new Map();

        this.statsFileData.rows.forEach(row => {
            const gene = row[geneColIdx]?.toUpperCase().trim();
            if (!gene) return;

            genes.push(gene);

            const stats = { gene };
            if (lfcColIdx !== '') {
                stats.lfc = parseFloat(row[parseInt(lfcColIdx)]) || null;
            }
            if (fdrColIdx !== '') {
                stats.fdr = parseFloat(row[parseInt(fdrColIdx)]) || null;
            }
            this.geneStats.set(gene, stats);
        });

        // Update textarea
        document.getElementById('geneTextarea').value = genes.join('\n');
        this.updateGeneCount();

        // Show stats controls if we have statistics
        const hasLfc = lfcColIdx !== '';
        const hasFdr = fdrColIdx !== '';
        if (hasLfc || hasFdr) {
            document.getElementById('statsControls').style.display = 'block';
        }

        // Show stats loaded message
        let msg = `Loaded ${genes.length} genes`;
        if (hasLfc || hasFdr) {
            msg += ` with statistics (${hasLfc ? 'LFC' : ''}${hasLfc && hasFdr ? ', ' : ''}${hasFdr ? 'FDR' : ''})`;
        }
        alert(msg);
    }

    loadTestGenesWithStats() {
        // Test data with LFC and FDR values
        const testData = [
            { gene: 'TP53', lfc: 1.8, fdr: 0.001 },
            { gene: 'MDM2', lfc: -2.3, fdr: 0.0001 },
            { gene: 'BRCA1', lfc: 0.3, fdr: 0.15 },
            { gene: 'MYC', lfc: -1.5, fdr: 0.02 },
            { gene: 'KRAS', lfc: 2.1, fdr: 0.005 },
            { gene: 'EGFR', lfc: -0.8, fdr: 0.08 },
            { gene: 'PTEN', lfc: 1.2, fdr: 0.03 },
            { gene: 'AKT1', lfc: -0.4, fdr: 0.25 },
            { gene: 'PIK3CA', lfc: 0.9, fdr: 0.04 },
            { gene: 'BRAF', lfc: -1.9, fdr: 0.002 },
            { gene: 'ATM', lfc: 0.7, fdr: 0.06 },
            { gene: 'CHEK2', lfc: 1.1, fdr: 0.01 },
            { gene: 'RB1', lfc: -0.6, fdr: 0.12 },
            { gene: 'CDKN2A', lfc: 1.4, fdr: 0.007 },
            { gene: 'SMAD4', lfc: -0.3, fdr: 0.35 },
            { gene: 'ARID1A', lfc: 0.5, fdr: 0.09 },
            { gene: 'NFE2L2', lfc: -1.2, fdr: 0.015 },
            { gene: 'KEAP1', lfc: 0.8, fdr: 0.05 },
            { gene: 'STK11', lfc: -0.9, fdr: 0.04 },
            { gene: 'CREBBP', lfc: 0.4, fdr: 0.18 }
        ];

        this.geneStats = new Map();
        const genes = [];

        testData.forEach(d => {
            genes.push(d.gene);
            this.geneStats.set(d.gene, d);
        });

        document.getElementById('geneTextarea').value = genes.join('\n');
        this.updateGeneCount();

        // Show stats controls
        document.getElementById('statsControls').style.display = 'block';

        alert(`Loaded ${genes.length} test genes with LFC and FDR statistics`);
    }

    downloadSampleStatsFile() {
        // Create sample CSV with gene stats format
        const sampleData = [
            ['Gene', 'LFC', 'FDR'],
            ['TP53', '1.8', '0.001'],
            ['MDM2', '-2.3', '0.0001'],
            ['BRCA1', '0.3', '0.15'],
            ['MYC', '-1.5', '0.02'],
            ['KRAS', '2.1', '0.005'],
            ['EGFR', '-0.8', '0.08'],
            ['PTEN', '1.2', '0.03'],
            ['AKT1', '-0.4', '0.25'],
            ['PIK3CA', '0.9', '0.04'],
            ['BRAF', '-1.9', '0.002'],
            ['ATM', '0.7', '0.06'],
            ['CHEK2', '1.1', '0.01'],
            ['RB1', '-0.6', '0.12'],
            ['CDKN2A', '1.4', '0.007'],
            ['SMAD4', '-0.3', '0.35'],
            ['ARID1A', '0.5', '0.09'],
            ['NFE2L2', '-1.2', '0.015'],
            ['KEAP1', '0.8', '0.05'],
            ['STK11', '-0.9', '0.04'],
            ['CREBBP', '0.4', '0.18']
        ];

        const csv = sampleData.map(row => row.join(',')).join('\n');
        this.downloadFile(csv, 'sample_genes_with_stats.csv', 'text/csv');
    }

    getGeneList() {
        const text = document.getElementById('geneTextarea').value;
        return text.split(/\s+/)
            .map(g => g.toUpperCase().trim())
            .filter(g => g !== '' && this.geneIndex.has(g));
    }

    getFilteredCellLineIndices() {
        const lineageFilter = document.getElementById('lineageFilter').value;
        if (!lineageFilter) {
            return Array.from({ length: this.nCellLines }, (_, i) => i);
        }

        const indices = [];
        this.metadata.cellLines.forEach((cellLine, idx) => {
            if (this.cellLineMetadata.lineage &&
                this.cellLineMetadata.lineage[cellLine] === lineageFilter) {
                indices.push(idx);
            }
        });
        return indices;
    }

    runAnalysis() {
        const geneList = this.getGeneList();
        const mode = document.querySelector('input[name="analysisMode"]:checked').value;
        const cutoff = parseFloat(document.getElementById('correlationCutoff').value);
        const minN = parseInt(document.getElementById('minCellLines').value);
        const minSlope = parseFloat(document.getElementById('minSlope').value);

        if (geneList.length === 0) {
            this.showStatus('error', 'Please enter at least one valid gene');
            return;
        }

        if (mode === 'analysis' && geneList.length < 2) {
            this.showStatus('error', 'Analysis mode requires at least 2 genes');
            return;
        }

        const cellLineIndices = this.getFilteredCellLineIndices();
        if (cellLineIndices.length < 10) {
            this.showStatus('error', 'Too few cell lines match the filter');
            return;
        }

        this.showStatus('info', 'Running correlation analysis...');

        // Use setTimeout to allow UI to update
        setTimeout(() => {
            try {
                this.results = this.calculateCorrelations(geneList, mode, cutoff, minN, minSlope, cellLineIndices);
                if (this.results.success) {
                    this.displayResults();
                    this.showStatus('success',
                        `&#10003; Analysis complete: ${this.results.correlations.length} correlations, ${this.results.clusters.length} genes in network`);
                } else {
                    this.showStatus('error', this.results.error);
                }
            } catch (error) {
                console.error('Analysis error:', error);
                this.showStatus('error', 'Analysis failed: ' + error.message);
            }
        }, 50);
    }

    calculateCorrelations(geneList, mode, cutoff, minN, minSlope, cellLineIndices) {
        const correlations = [];
        let targetGenes;

        if (mode === 'analysis') {
            // Analysis mode: correlate genes within the list
            targetGenes = geneList;
        } else {
            // Design mode: correlate against all genes
            targetGenes = Array.from(this.geneIndex.keys());
        }

        // Get gene data for input genes
        const inputData = new Map();
        geneList.forEach(gene => {
            const idx = this.geneIndex.get(gene);
            const fullData = this.getGeneData(idx);
            const filteredData = cellLineIndices.map(i => fullData[i]);
            inputData.set(gene, filteredData);
        });

        // Calculate correlations
        for (let i = 0; i < geneList.length; i++) {
            const gene1 = geneList[i];
            const data1 = inputData.get(gene1);

            const startJ = mode === 'analysis' ? i + 1 : 0;
            for (let j = startJ; j < targetGenes.length; j++) {
                const gene2 = targetGenes[j];
                if (gene1 === gene2) continue;
                // In analysis mode, startJ = i+1 already prevents duplicates

                let data2;
                if (inputData.has(gene2)) {
                    data2 = inputData.get(gene2);
                } else {
                    const idx = this.geneIndex.get(gene2);
                    const fullData = this.getGeneData(idx);
                    data2 = cellLineIndices.map(i => fullData[i]);
                }

                const result = this.pearsonWithSlope(data1, data2);
                if (result.n >= minN && Math.abs(result.correlation) >= cutoff && Math.abs(result.slope) >= minSlope) {
                    correlations.push({
                        gene1: gene1,
                        gene2: gene2,
                        correlation: Math.round(result.correlation * 1000) / 1000,
                        slope: Math.round(result.slope * 1000) / 1000,
                        n: result.n,
                        cluster: 0
                    });
                }
            }
        }

        if (correlations.length === 0) {
            return { success: false, error: `No correlations found above cutoff of ${cutoff}` };
        }

        // Assign clusters using simple connected components
        const clusters = this.findClusters(correlations);

        // Calculate mean effect for each gene
        const clusterData = clusters.map(gene => {
            const idx = this.geneIndex.get(gene);
            const data = this.getGeneData(idx);
            const validData = Array.from(data).filter(v => !isNaN(v));
            const mean = validData.reduce((a, b) => a + b, 0) / validData.length;
            const variance = validData.reduce((a, b) => a + (b - mean) ** 2, 0) / validData.length;
            const sd = Math.sqrt(variance);

            return {
                gene: gene,
                cluster: correlations.find(c => c.gene1 === gene || c.gene2 === gene)?.cluster || 0,
                meanEffect: Math.round(mean * 100) / 100,
                sdEffect: Math.round(sd * 100) / 100,
                inGeneList: geneList.includes(gene)
            };
        });

        return {
            success: true,
            correlations: correlations,
            clusters: clusterData,
            geneList: geneList,
            mode: mode,
            cutoff: cutoff,
            nCellLines: cellLineIndices.length
        };
    }

    pearsonWithSlope(x, y) {
        let sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0, n = 0;

        for (let i = 0; i < x.length; i++) {
            if (!isNaN(x[i]) && !isNaN(y[i])) {
                sumX += x[i];
                sumY += y[i];
                sumXY += x[i] * y[i];
                sumX2 += x[i] * x[i];
                sumY2 += y[i] * y[i];
                n++;
            }
        }

        if (n < 3) return { correlation: NaN, slope: NaN, n: 0 };

        const meanX = sumX / n;
        const meanY = sumY / n;
        const numerator = sumXY - n * meanX * meanY;
        const denomX = Math.sqrt(sumX2 - n * meanX * meanX);
        const denomY = Math.sqrt(sumY2 - n * meanY * meanY);

        const correlation = denomX * denomY === 0 ? NaN : numerator / (denomX * denomY);
        const slope = (sumX2 - n * meanX * meanX) === 0 ? NaN : numerator / (sumX2 - n * meanX * meanX);

        return { correlation, slope, n };
    }

    median(arr) {
        if (!arr || arr.length === 0) return NaN;
        const sorted = [...arr].filter(v => !isNaN(v)).sort((a, b) => a - b);
        const mid = Math.floor(sorted.length / 2);
        return sorted.length % 2 !== 0 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
    }

    findClusters(correlations) {
        const genes = new Set();
        correlations.forEach(c => {
            genes.add(c.gene1);
            genes.add(c.gene2);
        });

        const geneArray = Array.from(genes);
        const parent = {};
        geneArray.forEach(g => parent[g] = g);

        const find = (x) => {
            if (parent[x] !== x) parent[x] = find(parent[x]);
            return parent[x];
        };

        const union = (x, y) => {
            const px = find(x), py = find(y);
            if (px !== py) parent[px] = py;
        };

        correlations.forEach(c => union(c.gene1, c.gene2));

        // Assign cluster numbers
        const clusterMap = {};
        let clusterNum = 1;
        geneArray.forEach(g => {
            const root = find(g);
            if (!(root in clusterMap)) {
                clusterMap[root] = clusterNum++;
            }
        });

        correlations.forEach(c => {
            c.cluster = clusterMap[find(c.gene1)];
        });

        return geneArray;
    }

    showStatus(type, message) {
        const display = document.getElementById('analysisStatus');
        display.innerHTML = `<div class="status-box status-${type}">${message}</div>`;
    }

    displayResults() {
        this.displayNetwork();
        this.displayCorrelationsTable();
        this.displayClustersTable();
        this.displaySummary();
    }

    displayNetwork() {
        const container = document.getElementById('networkPlot');
        container.innerHTML = '';

        const nodes = [];
        const edges = [];
        const geneSet = new Set();

        const cutoff = this.results.cutoff;
        const edgeWidthBase = parseInt(document.getElementById('netEdgeWidth').value);
        const nodeSize = parseInt(document.getElementById('netNodeSize').value);
        const fontSize = parseInt(document.getElementById('netFontSize').value);

        // Create edges
        this.results.correlations.forEach((c, idx) => {
            geneSet.add(c.gene1);
            geneSet.add(c.gene2);

            const width = 1 + (Math.abs(c.correlation) - cutoff) / (1 - cutoff) * (edgeWidthBase * 3);
            edges.push({
                id: `edge_${idx}`,
                from: c.gene1,
                to: c.gene2,
                width: width,
                color: c.correlation > 0 ? '#3182ce' : '#e53e3e',
                title: `r = ${c.correlation.toFixed(3)}`,
                correlation: c.correlation
            });
        });

        // Create nodes
        geneSet.forEach(gene => {
            const cluster = this.results.clusters.find(c => c.gene === gene);
            const isInput = this.results.geneList.includes(gene);
            const geneStat = this.geneStats?.get(gene);

            // Build title with available information
            let titleLines = [gene];
            titleLines.push(`GE mean: ${cluster?.meanEffect || 'N/A'}`);
            titleLines.push(`GE SD: ${cluster?.sdEffect || 'N/A'}`);
            if (geneStat?.lfc !== undefined && geneStat?.lfc !== null) {
                titleLines.push(`LFC: ${geneStat.lfc.toFixed(3)}`);
            }
            if (geneStat?.fdr !== undefined && geneStat?.fdr !== null) {
                titleLines.push(`FDR: ${geneStat.fdr.toExponential(2)}`);
            }

            nodes.push({
                id: gene,
                label: gene,
                size: nodeSize,
                font: { size: fontSize, color: '#333' },
                color: {
                    background: this.results.mode === 'design' ?
                        (isInput ? '#16a34a' : '#86efac') : '#16a34a',
                    border: '#ffffff'
                },
                borderWidth: 2,
                title: titleLines.join('\n')
            });
        });

        const data = { nodes: new vis.DataSet(nodes), edges: new vis.DataSet(edges) };
        const options = {
            nodes: {
                shape: 'dot',
                scaling: {
                    min: 10,
                    max: 60
                },
                font: {
                    size: fontSize,
                    color: '#333'
                }
            },
            edges: {
                smooth: false
            },
            physics: {
                enabled: true,
                solver: 'forceAtlas2Based',
                forceAtlas2Based: {
                    gravitationalConstant: -50,
                    centralGravity: 0.01,
                    springLength: 100,
                    springConstant: 0.08
                },
                stabilization: { iterations: 150 }
            },
            interaction: {
                hover: true,
                tooltipDelay: 100,
                navigationButtons: true,
                keyboard: true
            }
        };

        this.network = new vis.Network(container, data, options);
        this.networkData = data;
        this.hiddenNodes = [];

        // Double-click to hide node
        this.network.on('doubleClick', (params) => {
            if (params.nodes.length > 0) {
                const nodeId = params.nodes[0];
                const node = this.networkData.nodes.get(nodeId);
                if (node) {
                    this.hiddenNodes.push(node);
                    this.networkData.nodes.remove(nodeId);
                }
            }
        });

        // Show legend
        document.getElementById('networkLegend').style.display = 'flex';
        const legendNodeType = document.getElementById('legendNodeType');
        if (this.results.mode === 'design') {
            legendNodeType.innerHTML = `
                <strong>Node Type:</strong>
                <span class="legend-item"><span class="legend-dot" style="background: #16a34a;"></span> Input</span>
                <span class="legend-item"><span class="legend-dot" style="background: #86efac;"></span> Correlated</span>
            `;
        } else {
            legendNodeType.innerHTML = '';
        }

        // Update edge thickness legend with actual data values
        this.updateEdgeLegend(edgeWidthBase, cutoff);
    }

    updateNetworkStyle() {
        if (!this.network || !this.networkData) return;

        const nodeSize = parseInt(document.getElementById('netNodeSize').value);
        const fontSize = parseInt(document.getElementById('netFontSize').value);
        const edgeWidthBase = parseInt(document.getElementById('netEdgeWidth').value);
        const cutoff = this.results?.cutoff || 0.5;

        // Update nodes
        const nodeUpdates = [];
        this.networkData.nodes.forEach(node => {
            nodeUpdates.push({
                id: node.id,
                size: nodeSize,
                font: { size: fontSize, color: '#333' }
            });
        });
        this.networkData.nodes.update(nodeUpdates);

        // Update edges
        const edgeUpdates = [];
        this.networkData.edges.forEach(edge => {
            const correlation = Math.abs(edge.correlation || 0.5);
            const width = 1 + (correlation - cutoff) / (1 - cutoff) * (edgeWidthBase * 3);
            edgeUpdates.push({
                id: edge.id,
                width: Math.max(1, width)
            });
        });
        this.networkData.edges.update(edgeUpdates);

        // Update legend for edge thickness
        this.updateEdgeLegend(edgeWidthBase, cutoff);
    }

    updateEdgeLegend(edgeWidthBase, cutoff) {
        const legendEl = document.getElementById('legendEdgeThickness');
        if (!legendEl) return;

        // Get actual correlation range from current network data
        const correlations = this.results.correlations.map(c => Math.abs(c.correlation));
        const rawMin = correlations.length > 0 ? Math.min(...correlations) : cutoff;
        const rawMax = correlations.length > 0 ? Math.max(...correlations) : 1.0;

        // Round min down and max up to ensure legend encompasses all data
        const minCorr = Math.floor(rawMin * 10) / 10;
        const maxCorr = Math.ceil(rawMax * 10) / 10;
        const midCorr = (minCorr + maxCorr) / 2;

        // Calculate widths for actual data range
        const widthMin = Math.max(1, 1 + (minCorr - cutoff) / (1 - cutoff) * (edgeWidthBase * 3));
        const widthMid = Math.max(1, 1 + (midCorr - cutoff) / (1 - cutoff) * (edgeWidthBase * 3));
        const widthMax = Math.max(1, 1 + (maxCorr - cutoff) / (1 - cutoff) * (edgeWidthBase * 3));

        legendEl.innerHTML = `
            <strong>Edge Thickness:</strong>
            <div class="legend-item"><span class="legend-line" style="background: #666; height: ${Math.round(widthMin)}px;"></span> r=${minCorr.toFixed(1)}</div>
            <div class="legend-item"><span class="legend-line" style="background: #666; height: ${Math.round(widthMid)}px;"></span> r=${midCorr.toFixed(1)}</div>
            <div class="legend-item"><span class="legend-line" style="background: #666; height: ${Math.round(widthMax)}px;"></span> r=${maxCorr.toFixed(1)}</div>
        `;

        // Store for use in PNG/SVG export
        this.edgeLegendValues = { minCorr, midCorr, maxCorr, widthMin, widthMid, widthMax };
    }

    displayCorrelationsTable() {
        const tbody = document.getElementById('correlationsBody');
        tbody.innerHTML = '';

        this.results.correlations
            .sort((a, b) => Math.abs(b.correlation) - Math.abs(a.correlation))
            .forEach((c) => {
                const tr = document.createElement('tr');
                const corrClass = c.correlation > 0 ? 'corr-positive' : 'corr-negative';
                tr.innerHTML = `
                    <td>${c.gene1}</td>
                    <td>${c.gene2}</td>
                    <td class="${corrClass}">${c.correlation.toFixed(3)}</td>
                    <td>${c.slope.toFixed(3)}</td>
                    <td>${c.n}</td>
                    <td>${c.cluster}</td>
                    <td>
                        <button class="btn btn-outline-primary btn-sm inspect-btn" data-gene1="${c.gene1}" data-gene2="${c.gene2}">Inspect</button>
                        <button class="btn btn-outline-info btn-sm tissue-btn" data-gene1="${c.gene1}" data-gene2="${c.gene2}" style="margin-left: 4px;">By tissue</button>
                    </td>
                `;
                // Add click handlers
                tr.querySelector('.inspect-btn').addEventListener('click', () => {
                    this.openInspectByGenes(c.gene1, c.gene2);
                });
                tr.querySelector('.tissue-btn').addEventListener('click', () => {
                    this.openByTissueByGenes(c.gene1, c.gene2);
                });
                tbody.appendChild(tr);
            });
    }

    displayClustersTable() {
        const tbody = document.getElementById('clustersBody');
        const thead = document.getElementById('clustersHead');
        tbody.innerHTML = '';

        // Check if we have stats
        const hasStats = this.geneStats && this.geneStats.size > 0;

        // Update header based on whether we have stats
        if (hasStats) {
            thead.innerHTML = `
                <tr>
                    <th data-sort="gene">Gene</th>
                    <th data-sort="cluster">Cluster</th>
                    <th data-sort="meanEffect">Mean Effect</th>
                    <th data-sort="sdEffect">SD</th>
                    <th data-sort="lfc">LFC</th>
                    <th data-sort="fdr">FDR</th>
                </tr>
            `;
        } else {
            thead.innerHTML = `
                <tr>
                    <th data-sort="gene">Gene</th>
                    <th data-sort="cluster">Cluster</th>
                    <th data-sort="meanEffect">Mean Effect</th>
                    <th data-sort="sdEffect">SD</th>
                </tr>
            `;
        }

        // Re-attach sorting event listeners
        thead.querySelectorAll('th[data-sort]').forEach(th => {
            th.addEventListener('click', () => this.sortTable(th));
        });

        this.results.clusters
            .sort((a, b) => a.cluster - b.cluster || a.gene.localeCompare(b.gene))
            .forEach(c => {
                const tr = document.createElement('tr');
                const geneStat = this.geneStats?.get(c.gene);

                if (hasStats) {
                    const lfc = geneStat?.lfc !== null && geneStat?.lfc !== undefined
                        ? geneStat.lfc.toFixed(2) : '-';
                    const fdr = geneStat?.fdr !== null && geneStat?.fdr !== undefined
                        ? geneStat.fdr.toExponential(2) : '-';

                    tr.innerHTML = `
                        <td>${c.gene}${c.inGeneList && this.results.mode === 'design' ? '*' : ''}</td>
                        <td>${c.cluster}</td>
                        <td>${c.meanEffect}</td>
                        <td>${c.sdEffect}</td>
                        <td>${lfc}</td>
                        <td>${fdr}</td>
                    `;
                } else {
                    tr.innerHTML = `
                        <td>${c.gene}${c.inGeneList && this.results.mode === 'design' ? '*' : ''}</td>
                        <td>${c.cluster}</td>
                        <td>${c.meanEffect}</td>
                        <td>${c.sdEffect}</td>
                    `;
                }
                tbody.appendChild(tr);
            });
    }

    displaySummary() {
        const text = document.getElementById('summaryText');
        const lineage = document.getElementById('lineageFilter').value || 'All lineages';

        // Build synonyms section if any were used
        let synonymsSection = '';
        if (this.synonymsUsed && this.synonymsUsed.length > 0) {
            synonymsSection = `
Synonyms/Orthologs Used:
${this.synonymsUsed.map(s => `  ${s.original} → ${s.replacement} (${s.source})`).join('\n')}
`;
        }

        // Build unrecognized genes section
        let unrecognizedSection = '';
        if (this.genesNotFound && this.genesNotFound.length > 0) {
            unrecognizedSection = `
Unrecognized Gene Names (${this.genesNotFound.length}):
${this.genesNotFound.join(', ')}
`;
        }

        // Calculate number of clusters
        const numClusters = this.results.correlations.length > 0
            ? Math.max(...this.results.correlations.map(c => c.cluster))
            : 0;

        text.textContent = `Gene Correlation Analysis Summary
================================

Analysis Mode: ${this.results.mode === 'analysis' ? 'Analysis (within gene list)' : 'Design (find correlated genes)'}
Correlation Cutoff: ${this.results.cutoff}
Minimum Cell Lines: ${document.getElementById('minCellLines').value}
Minimum Slope: ${document.getElementById('minSlope').value}
Lineage Filter: ${lineage}

Input Genes: ${this.results.geneList.length}
${this.results.geneList.join(', ')}
${synonymsSection}${unrecognizedSection}
Results:
- Total correlations found: ${this.results.correlations.length}
- Genes in network: ${this.results.clusters.length}
- Number of clusters: ${numClusters}
- Cell lines analyzed: ${this.results.nCellLines}
`;
    }

    filterTable(tableId, query) {
        const tbody = document.getElementById(tableId);
        const rows = tbody.querySelectorAll('tr');
        const lowerQuery = query.toLowerCase();

        rows.forEach(row => {
            const text = row.textContent.toLowerCase();
            row.style.display = text.includes(lowerQuery) ? '' : 'none';
        });
    }

    sortTable(th) {
        const table = th.closest('table');
        const tbody = table.querySelector('tbody');
        const rows = Array.from(tbody.querySelectorAll('tr'));
        const colIndex = Array.from(th.parentNode.children).indexOf(th);
        const sortKey = th.dataset.sort;
        const isNumeric = ['correlation', 'slope', 'n', 'cluster', 'meanEffect', 'sdEffect'].includes(sortKey);

        const currentDir = th.dataset.dir || 'asc';
        const newDir = currentDir === 'asc' ? 'desc' : 'asc';
        th.dataset.dir = newDir;

        rows.sort((a, b) => {
            const aVal = a.children[colIndex]?.textContent || '';
            const bVal = b.children[colIndex]?.textContent || '';

            if (isNumeric) {
                const aNum = parseFloat(aVal) || 0;
                const bNum = parseFloat(bVal) || 0;
                return newDir === 'asc' ? aNum - bNum : bNum - aNum;
            } else {
                return newDir === 'asc' ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
            }
        });

        rows.forEach(row => tbody.appendChild(row));
    }

    downloadCSV(type) {
        if (!this.results) return;

        let csv, filename;
        if (type === 'correlations') {
            csv = 'Gene1,Gene2,Correlation,Slope,N,Cluster\n';
            this.results.correlations.forEach(c => {
                csv += `${c.gene1},${c.gene2},${c.correlation},${c.slope},${c.n},${c.cluster}\n`;
            });
            filename = 'correlations.csv';
        } else {
            csv = 'Gene,Cluster,Mean_Effect,SD_Effect\n';
            this.results.clusters.forEach(c => {
                csv += `${c.gene},${c.cluster},${c.meanEffect},${c.sdEffect}\n`;
            });
            filename = 'clusters.csv';
        }

        this.downloadFile(csv, filename, 'text/csv');
    }

    downloadSummary() {
        const text = document.getElementById('summaryText').textContent;
        this.downloadFile(text, 'summary.txt', 'text/plain');
    }

    downloadFile(content, filename, mimeType) {
        const blob = new Blob([content], { type: mimeType });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        a.click();
        URL.revokeObjectURL(url);
    }

    downloadNetworkPNG() {
        if (!this.network) return;

        const networkCanvas = document.querySelector('#networkPlot canvas');
        if (!networkCanvas) {
            console.error('Network canvas not found');
            return;
        }

        // Create a new canvas that includes the legend - larger for publication quality
        const legendHeight = 160;
        const padding = 30;
        const networkWidth = networkCanvas.width;
        const networkHeight = networkCanvas.height;
        const totalWidth = networkWidth;
        const totalHeight = networkHeight + legendHeight + padding;

        const canvas = document.createElement('canvas');
        canvas.width = totalWidth;
        canvas.height = totalHeight;
        const ctx = canvas.getContext('2d');

        // Draw white background
        ctx.fillStyle = 'white';
        ctx.fillRect(0, 0, totalWidth, totalHeight);

        // Draw network
        ctx.drawImage(networkCanvas, 0, 0);

        // Draw legend background
        ctx.fillStyle = '#f9fafb';
        ctx.strokeStyle = '#e5e7eb';
        ctx.lineWidth = 1;
        ctx.fillRect(15, networkHeight + 10, totalWidth - 30, legendHeight - 10);
        ctx.strokeRect(15, networkHeight + 10, totalWidth - 30, legendHeight - 10);

        // Draw legend - LARGER fonts for publication
        const legendY = networkHeight + padding + 10;
        const titleFont = 'bold 16px Arial';
        const textFont = '14px Arial';
        const smallFont = '13px Arial';

        let legendX = 40;

        // Correlation legend
        ctx.font = titleFont;
        ctx.fillStyle = '#333';
        ctx.fillText('Correlation:', legendX, legendY);
        ctx.font = textFont;

        // Positive correlation line
        ctx.strokeStyle = '#3182ce';
        ctx.lineWidth = 4;
        ctx.beginPath();
        ctx.moveTo(legendX, legendY + 22);
        ctx.lineTo(legendX + 35, legendY + 22);
        ctx.stroke();
        ctx.fillStyle = '#333';
        ctx.fillText('Positive', legendX + 42, legendY + 27);

        // Negative correlation line
        ctx.strokeStyle = '#e53e3e';
        ctx.beginPath();
        ctx.moveTo(legendX, legendY + 48);
        ctx.lineTo(legendX + 35, legendY + 48);
        ctx.stroke();
        ctx.fillStyle = '#333';
        ctx.fillText('Negative', legendX + 42, legendY + 53);

        legendX += 160;

        // Edge thickness legend - use actual data values
        ctx.font = titleFont;
        ctx.fillText('Edge Thickness:', legendX, legendY);
        ctx.font = textFont;
        ctx.strokeStyle = '#666';

        const edgeWidthBase = parseInt(document.getElementById('netEdgeWidth').value) || 3;
        const legendVals = this.edgeLegendValues || { minCorr: 0.5, midCorr: 0.75, maxCorr: 1.0 };
        const cutoff = this.results?.cutoff || 0.5;

        // Min correlation
        ctx.lineWidth = Math.max(2, 2 + (legendVals.minCorr - cutoff) / (1 - cutoff) * (edgeWidthBase * 4));
        ctx.beginPath();
        ctx.moveTo(legendX, legendY + 22);
        ctx.lineTo(legendX + 35, legendY + 22);
        ctx.stroke();
        ctx.fillStyle = '#333';
        ctx.fillText(`r = ${legendVals.minCorr.toFixed(2)}`, legendX + 42, legendY + 27);

        // Mid correlation
        ctx.lineWidth = Math.max(2, 2 + (legendVals.midCorr - cutoff) / (1 - cutoff) * (edgeWidthBase * 4));
        ctx.beginPath();
        ctx.moveTo(legendX, legendY + 48);
        ctx.lineTo(legendX + 35, legendY + 48);
        ctx.stroke();
        ctx.fillText(`r = ${legendVals.midCorr.toFixed(2)}`, legendX + 42, legendY + 53);

        // Max correlation
        ctx.lineWidth = Math.max(2, 2 + (legendVals.maxCorr - cutoff) / (1 - cutoff) * (edgeWidthBase * 4));
        ctx.beginPath();
        ctx.moveTo(legendX, legendY + 74);
        ctx.lineTo(legendX + 35, legendY + 74);
        ctx.stroke();
        ctx.fillText(`r = ${legendVals.maxCorr.toFixed(2)}`, legendX + 42, legendY + 79);

        legendX += 160;

        // Node type legend (for design mode)
        if (this.results?.mode === 'design') {
            ctx.font = titleFont;
            ctx.fillStyle = '#333';
            ctx.fillText('Node Type:', legendX, legendY);
            ctx.font = textFont;

            // Input gene
            ctx.fillStyle = '#16a34a';
            ctx.beginPath();
            ctx.arc(legendX + 12, legendY + 25, 10, 0, Math.PI * 2);
            ctx.fill();
            ctx.fillStyle = '#333';
            ctx.fillText('Input', legendX + 28, legendY + 30);

            // Correlated gene
            ctx.fillStyle = '#86efac';
            ctx.beginPath();
            ctx.arc(legendX + 12, legendY + 52, 10, 0, Math.PI * 2);
            ctx.fill();
            ctx.fillStyle = '#333';
            ctx.fillText('Correlated', legendX + 28, legendY + 57);

            legendX += 140;
        }

        // Color by stats legend
        const colorByStats = document.getElementById('colorByStats').checked;
        if (colorByStats && this.geneStats && this.geneStats.size > 0) {
            const colorStatType = document.querySelector('input[name="colorStatType"]:checked')?.value || 'signed_lfc';
            const colorScale = document.querySelector('input[name="colorScale"]:checked')?.value || 'all';

            // Get stats based on scale option
            let stats;
            if (colorScale === 'network') {
                const networkGenes = [];
                this.networkData.nodes.forEach(node => {
                    const geneStat = this.geneStats.get(node.id);
                    if (geneStat) networkGenes.push(geneStat);
                });
                stats = networkGenes;
            } else {
                stats = Array.from(this.geneStats.values());
            }

            const scaleLabel = colorScale === 'network' ? ' (network)' : ' (all genes)';

            ctx.font = titleFont;
            ctx.fillStyle = '#333';
            ctx.fillText('Node Color:', legendX, legendY);
            ctx.font = textFont;

            // Draw gradient - LARGER
            const gradientWidth = 120;
            const gradientHeight = 18;
            const gradY = legendY + 18;

            if (colorStatType === 'signed_lfc') {
                const lfcValues = stats.map(s => s.lfc).filter(v => v !== null && !isNaN(v));
                const minLfc = Math.min(...lfcValues);
                const maxLfc = Math.max(...lfcValues);

                // Draw gradient
                const gradient = ctx.createLinearGradient(legendX, 0, legendX + gradientWidth, 0);
                gradient.addColorStop(0, '#2166ac');
                gradient.addColorStop(0.5, '#f7f7f7');
                gradient.addColorStop(1, '#b2182b');
                ctx.fillStyle = gradient;
                ctx.fillRect(legendX, gradY, gradientWidth, gradientHeight);
                ctx.strokeStyle = '#999';
                ctx.lineWidth = 1;
                ctx.strokeRect(legendX, gradY, gradientWidth, gradientHeight);

                // Labels
                ctx.fillStyle = '#333';
                ctx.font = smallFont;
                ctx.fillText(minLfc.toFixed(1), legendX, gradY + gradientHeight + 16);
                ctx.fillText('Signed LFC' + scaleLabel, legendX, gradY - 4);
                ctx.fillText(maxLfc.toFixed(1), legendX + gradientWidth - 20, gradY + gradientHeight + 16);
            } else if (colorStatType === 'abs_lfc') {
                const lfcValues = stats.map(s => Math.abs(s.lfc)).filter(v => v !== null && !isNaN(v));
                const maxLfc = Math.max(...lfcValues);

                const gradient = ctx.createLinearGradient(legendX, 0, legendX + gradientWidth, 0);
                gradient.addColorStop(0, '#f5f5f5');
                gradient.addColorStop(0.5, '#fdae61');
                gradient.addColorStop(1, '#d7191c');
                ctx.fillStyle = gradient;
                ctx.fillRect(legendX, gradY, gradientWidth, gradientHeight);
                ctx.strokeStyle = '#999';
                ctx.lineWidth = 1;
                ctx.strokeRect(legendX, gradY, gradientWidth, gradientHeight);

                ctx.fillStyle = '#333';
                ctx.font = smallFont;
                ctx.fillText('0', legendX, gradY + gradientHeight + 16);
                ctx.fillText('|LFC|' + scaleLabel, legendX, gradY - 4);
                ctx.fillText(maxLfc.toFixed(1), legendX + gradientWidth - 20, gradY + gradientHeight + 16);
            } else if (colorStatType === 'fdr') {
                const fdrValues = stats.map(s => s.fdr).filter(v => v !== null && !isNaN(v) && v > 0);
                const minFdr = Math.min(...fdrValues);

                const gradient = ctx.createLinearGradient(legendX, 0, legendX + gradientWidth, 0);
                gradient.addColorStop(0, '#d7191c');
                gradient.addColorStop(0.5, '#fdae61');
                gradient.addColorStop(1, '#f5f5f5');
                ctx.fillStyle = gradient;
                ctx.fillRect(legendX, gradY, gradientWidth, gradientHeight);
                ctx.strokeStyle = '#999';
                ctx.lineWidth = 1;
                ctx.strokeRect(legendX, gradY, gradientWidth, gradientHeight);

                ctx.fillStyle = '#333';
                ctx.font = smallFont;
                ctx.fillText(minFdr.toExponential(1), legendX - 5, gradY + gradientHeight + 16);
                ctx.fillText('FDR' + scaleLabel, legendX, gradY - 4);
                ctx.fillText('1', legendX + gradientWidth - 8, gradY + gradientHeight + 16);
            }

            // Missing indicator
            ctx.fillStyle = '#333';
            ctx.font = textFont;
            ctx.fillText('Missing:', legendX, gradY + gradientHeight + 42);
            ctx.fillStyle = '#cccccc';
            ctx.fillRect(legendX + 60, gradY + gradientHeight + 28, 18, 18);
            ctx.strokeStyle = '#999';
            ctx.strokeRect(legendX + 60, gradY + gradientHeight + 28, 18, 18);
        }

        // Create download link
        const dataURL = canvas.toDataURL('image/png');
        const a = document.createElement('a');
        a.href = dataURL;
        a.download = 'correlation_network.png';
        a.click();
    }

    downloadNetworkSVG() {
        if (!this.network || !this.networkData) return;

        // Build SVG from network data
        const container = document.getElementById('networkPlot');
        const width = container.clientWidth;
        const networkHeight = container.clientHeight;
        const legendHeight = 160;  // Larger for publication
        const totalHeight = networkHeight + legendHeight;

        // Get positions from vis.js
        const positions = this.network.getPositions();

        let svg = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${totalHeight}" viewBox="0 0 ${width} ${totalHeight}">
<defs>
    <linearGradient id="signedLfcGradient" x1="0%" y1="0%" x2="100%" y2="0%">
        <stop offset="0%" style="stop-color:#2166ac;stop-opacity:1" />
        <stop offset="50%" style="stop-color:#f7f7f7;stop-opacity:1" />
        <stop offset="100%" style="stop-color:#b2182b;stop-opacity:1" />
    </linearGradient>
    <linearGradient id="absLfcGradient" x1="0%" y1="0%" x2="100%" y2="0%">
        <stop offset="0%" style="stop-color:#f5f5f5;stop-opacity:1" />
        <stop offset="50%" style="stop-color:#fdae61;stop-opacity:1" />
        <stop offset="100%" style="stop-color:#d7191c;stop-opacity:1" />
    </linearGradient>
    <linearGradient id="fdrGradient" x1="0%" y1="0%" x2="100%" y2="0%">
        <stop offset="0%" style="stop-color:#d7191c;stop-opacity:1" />
        <stop offset="50%" style="stop-color:#fdae61;stop-opacity:1" />
        <stop offset="100%" style="stop-color:#f5f5f5;stop-opacity:1" />
    </linearGradient>
</defs>
<style>
  .node-label { font-family: Arial, sans-serif; font-size: 14px; fill: #333; }
  .legend-title { font-family: Arial, sans-serif; font-size: 16px; font-weight: bold; fill: #333; }
  .legend-text { font-family: Arial, sans-serif; font-size: 14px; fill: #333; }
  .legend-small { font-family: Arial, sans-serif; font-size: 13px; fill: #333; }
</style>
<rect width="100%" height="100%" fill="white"/>
`;

        // Draw edges
        this.networkData.edges.forEach(edge => {
            const from = positions[edge.from];
            const to = positions[edge.to];
            if (from && to) {
                // Convert from vis.js coordinates to SVG
                const x1 = from.x + width/2;
                const y1 = from.y + networkHeight/2;
                const x2 = to.x + width/2;
                const y2 = to.y + networkHeight/2;
                const color = edge.color || '#3182ce';
                const strokeWidth = edge.width || 2;
                svg += `  <line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${color}" stroke-width="${strokeWidth}" opacity="0.8"/>\n`;
            }
        });

        // Draw nodes
        const nodeSize = parseInt(document.getElementById('netNodeSize').value) || 25;
        const fontSize = parseInt(document.getElementById('netFontSize').value) || 16;
        this.networkData.nodes.forEach(node => {
            const pos = positions[node.id];
            if (pos) {
                const cx = pos.x + width/2;
                const cy = pos.y + networkHeight/2;
                const bgColor = node.color?.background || '#16a34a';
                svg += `  <circle cx="${cx}" cy="${cy}" r="${nodeSize/2}" fill="${bgColor}" stroke="white" stroke-width="2"/>\n`;

                // Handle multi-line labels
                const labelLines = (node.label || node.id).split('\n');
                labelLines.forEach((line, i) => {
                    const yOffset = cy + nodeSize/2 + 14 + (i * (fontSize - 2));
                    svg += `  <text x="${cx}" y="${yOffset}" text-anchor="middle" style="font-family: Arial; font-size: ${fontSize - 2}px; fill: #333;">${this.escapeXml(line)}</text>\n`;
                });
            }
        });

        // Draw legend - LARGER for publication
        const legendY = networkHeight + 35;
        let legendX = 40;

        // Legend background
        svg += `  <rect x="15" y="${networkHeight + 10}" width="${width - 30}" height="145" fill="#f9fafb" stroke="#e5e7eb" rx="4"/>\n`;

        // Correlation legend
        svg += `  <text x="${legendX}" y="${legendY}" class="legend-title">Correlation:</text>\n`;
        svg += `  <line x1="${legendX}" y1="${legendY + 22}" x2="${legendX + 35}" y2="${legendY + 22}" stroke="#3182ce" stroke-width="4"/>\n`;
        svg += `  <text x="${legendX + 42}" y="${legendY + 27}" class="legend-text">Positive</text>\n`;
        svg += `  <line x1="${legendX}" y1="${legendY + 48}" x2="${legendX + 35}" y2="${legendY + 48}" stroke="#e53e3e" stroke-width="4"/>\n`;
        svg += `  <text x="${legendX + 42}" y="${legendY + 53}" class="legend-text">Negative</text>\n`;

        legendX += 110;
        legendX += 160;

        // Edge thickness legend - use actual data values
        const cutoff = this.results?.cutoff || 0.5;
        const edgeWidthBase = parseInt(document.getElementById('netEdgeWidth').value) || 3;
        const legendVals = this.edgeLegendValues || { minCorr: 0.5, midCorr: 0.75, maxCorr: 1.0 };

        svg += `  <text x="${legendX}" y="${legendY}" class="legend-title">Edge Thickness:</text>\n`;

        const width1 = Math.max(2, 2 + (legendVals.minCorr - cutoff) / (1 - cutoff) * (edgeWidthBase * 4));
        const width2 = Math.max(2, 2 + (legendVals.midCorr - cutoff) / (1 - cutoff) * (edgeWidthBase * 4));
        const width3 = Math.max(2, 2 + (legendVals.maxCorr - cutoff) / (1 - cutoff) * (edgeWidthBase * 4));

        svg += `  <line x1="${legendX}" y1="${legendY + 22}" x2="${legendX + 35}" y2="${legendY + 22}" stroke="#666" stroke-width="${width1}"/>\n`;
        svg += `  <text x="${legendX + 42}" y="${legendY + 27}" class="legend-text">r = ${legendVals.minCorr.toFixed(2)}</text>\n`;
        svg += `  <line x1="${legendX}" y1="${legendY + 48}" x2="${legendX + 35}" y2="${legendY + 48}" stroke="#666" stroke-width="${width2}"/>\n`;
        svg += `  <text x="${legendX + 42}" y="${legendY + 53}" class="legend-text">r = ${legendVals.midCorr.toFixed(2)}</text>\n`;
        svg += `  <line x1="${legendX}" y1="${legendY + 74}" x2="${legendX + 35}" y2="${legendY + 74}" stroke="#666" stroke-width="${width3}"/>\n`;
        svg += `  <text x="${legendX + 42}" y="${legendY + 79}" class="legend-text">r = ${legendVals.maxCorr.toFixed(2)}</text>\n`;

        legendX += 160;

        // Node type legend (for design mode)
        if (this.results?.mode === 'design') {
            svg += `  <text x="${legendX}" y="${legendY}" class="legend-title">Node Type:</text>\n`;
            svg += `  <circle cx="${legendX + 12}" cy="${legendY + 25}" r="10" fill="#16a34a"/>\n`;
            svg += `  <text x="${legendX + 28}" y="${legendY + 30}" class="legend-text">Input</text>\n`;
            svg += `  <circle cx="${legendX + 12}" cy="${legendY + 52}" r="10" fill="#86efac"/>\n`;
            svg += `  <text x="${legendX + 28}" y="${legendY + 57}" class="legend-text">Correlated</text>\n`;

            legendX += 140;
        }

        // Color by stats legend
        const colorByStats = document.getElementById('colorByStats').checked;
        if (colorByStats && this.geneStats && this.geneStats.size > 0) {
            const colorStatType = document.querySelector('input[name="colorStatType"]:checked')?.value || 'signed_lfc';
            const colorScale = document.querySelector('input[name="colorScale"]:checked')?.value || 'all';

            // Get stats based on scale option
            let stats;
            if (colorScale === 'network') {
                const networkGenes = [];
                this.networkData.nodes.forEach(node => {
                    const geneStat = this.geneStats.get(node.id);
                    if (geneStat) networkGenes.push(geneStat);
                });
                stats = networkGenes;
            } else {
                stats = Array.from(this.geneStats.values());
            }

            const scaleLabel = colorScale === 'network' ? ' (network)' : ' (all genes)';

            svg += `  <text x="${legendX}" y="${legendY}" class="legend-title">Node Color:</text>\n`;

            const gradientWidth = 120;
            const gradientHeight = 18;
            const gradY = legendY + 18;

            if (colorStatType === 'signed_lfc') {
                const lfcValues = stats.map(s => s.lfc).filter(v => v !== null && !isNaN(v));
                const minLfc = Math.min(...lfcValues);
                const maxLfc = Math.max(...lfcValues);

                svg += `  <rect x="${legendX}" y="${gradY}" width="${gradientWidth}" height="${gradientHeight}" fill="url(#signedLfcGradient)" stroke="#999"/>\n`;
                svg += `  <text x="${legendX}" y="${gradY + gradientHeight + 16}" class="legend-small">${minLfc.toFixed(1)}</text>\n`;
                svg += `  <text x="${legendX}" y="${gradY - 4}" class="legend-small">Signed LFC${scaleLabel}</text>\n`;
                svg += `  <text x="${legendX + gradientWidth - 20}" y="${gradY + gradientHeight + 16}" class="legend-small">${maxLfc.toFixed(1)}</text>\n`;
            } else if (colorStatType === 'abs_lfc') {
                const lfcValues = stats.map(s => Math.abs(s.lfc)).filter(v => v !== null && !isNaN(v));
                const maxLfc = Math.max(...lfcValues);

                svg += `  <rect x="${legendX}" y="${gradY}" width="${gradientWidth}" height="${gradientHeight}" fill="url(#absLfcGradient)" stroke="#999"/>\n`;
                svg += `  <text x="${legendX}" y="${gradY + gradientHeight + 16}" class="legend-small">0</text>\n`;
                svg += `  <text x="${legendX}" y="${gradY - 4}" class="legend-small">|LFC|${scaleLabel}</text>\n`;
                svg += `  <text x="${legendX + gradientWidth - 20}" y="${gradY + gradientHeight + 16}" class="legend-small">${maxLfc.toFixed(1)}</text>\n`;
            } else if (colorStatType === 'fdr') {
                const fdrValues = stats.map(s => s.fdr).filter(v => v !== null && !isNaN(v) && v > 0);
                const minFdr = Math.min(...fdrValues);

                svg += `  <rect x="${legendX}" y="${gradY}" width="${gradientWidth}" height="${gradientHeight}" fill="url(#fdrGradient)" stroke="#999"/>\n`;
                svg += `  <text x="${legendX - 5}" y="${gradY + gradientHeight + 16}" class="legend-small">${minFdr.toExponential(1)}</text>\n`;
                svg += `  <text x="${legendX}" y="${gradY - 4}" class="legend-small">FDR${scaleLabel}</text>\n`;
                svg += `  <text x="${legendX + gradientWidth - 8}" y="${gradY + gradientHeight + 16}" class="legend-small">1</text>\n`;
            }

            // Missing indicator
            svg += `  <text x="${legendX}" y="${gradY + gradientHeight + 42}" class="legend-text">Missing:</text>\n`;
            svg += `  <rect x="${legendX + 60}" y="${gradY + gradientHeight + 28}" width="18" height="18" fill="#cccccc" stroke="#999" rx="2"/>\n`;
        }

        svg += '</svg>';

        // Download
        const blob = new Blob([svg], { type: 'image/svg+xml' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'correlation_network.svg';
        a.click();
        URL.revokeObjectURL(url);
    }

    escapeXml(str) {
        return str
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&apos;');
    }

    showHiddenNodes() {
        if (!this.network || !this.hiddenNodes || this.hiddenNodes.length === 0) {
            return;
        }

        // Get current style settings
        const nodeSize = parseInt(document.getElementById('netNodeSize').value);
        const fontSize = parseInt(document.getElementById('netFontSize').value);

        // Re-add hidden nodes with current style
        this.hiddenNodes.forEach(node => {
            node.size = nodeSize;
            node.font = { size: fontSize, color: '#333' };
            this.networkData.nodes.add(node);
        });
        this.hiddenNodes = [];
        this.network.fit();
    }

    updateNetworkLabels() {
        if (!this.network || !this.networkData) return;

        const showGE = document.getElementById('showGeneEffect').checked;
        const updates = [];

        this.networkData.nodes.forEach(node => {
            const cluster = this.results?.clusters?.find(c => c.gene === node.id);
            let label = node.id;

            if (showGE && cluster) {
                label = `${node.id}\n(GE:${cluster.meanEffect})`;
            }

            updates.push({
                id: node.id,
                label: label
            });
        });

        this.networkData.nodes.update(updates);
    }

    updateNetworkLabelsWithStats() {
        if (!this.network || !this.networkData) return;

        const statsDisplay = document.querySelector('input[name="statsLabelDisplay"]:checked')?.value || 'none';
        const showGE = document.getElementById('showGeneEffect').checked;
        const updates = [];

        this.networkData.nodes.forEach(node => {
            const cluster = this.results?.clusters?.find(c => c.gene === node.id);
            const geneStat = this.geneStats?.get(node.id);
            let label = node.id;

            // Add gene effect if checked
            if (showGE && cluster) {
                label = `${node.id}\n(GE:${cluster.meanEffect})`;
            }

            // Add stats label if selected
            if (statsDisplay === 'lfc' && geneStat?.lfc !== null && geneStat?.lfc !== undefined) {
                label += `\n[LFC:${geneStat.lfc.toFixed(1)}]`;
            } else if (statsDisplay === 'fdr' && geneStat?.fdr !== null && geneStat?.fdr !== undefined) {
                label += `\n[FDR:${geneStat.fdr.toExponential(1)}]`;
            }

            updates.push({
                id: node.id,
                label: label
            });
        });

        this.networkData.nodes.update(updates);
    }

    updateNetworkColors() {
        if (!this.network || !this.networkData) return;

        const colorByStats = document.getElementById('colorByStats').checked;
        const colorStatType = document.querySelector('input[name="colorStatType"]:checked')?.value || 'signed_lfc';
        const colorScale = document.querySelector('input[name="colorScale"]:checked')?.value || 'all';

        const updates = [];
        const colorLegend = document.getElementById('nodeColorLegendContent');
        const legendSection = document.getElementById('legendNodeColor');

        if (!colorByStats || !this.geneStats || this.geneStats.size === 0) {
            // Reset to default colors
            this.networkData.nodes.forEach(node => {
                const isInput = this.results?.geneList?.includes(node.id);
                updates.push({
                    id: node.id,
                    color: {
                        background: this.results?.mode === 'design' ?
                            (isInput ? '#16a34a' : '#86efac') : '#16a34a',
                        border: '#ffffff'
                    }
                });
            });
            if (colorLegend) colorLegend.innerHTML = '';
            if (legendSection) legendSection.style.display = 'none';
        } else {
            // Color by statistics
            if (legendSection) legendSection.style.display = 'block';

            // Get stats based on scale option
            let stats;
            if (colorScale === 'network') {
                // Only genes in the current network
                const networkGenes = [];
                this.networkData.nodes.forEach(node => {
                    const geneStat = this.geneStats.get(node.id);
                    if (geneStat) networkGenes.push(geneStat);
                });
                stats = networkGenes;
            } else {
                // All genes
                stats = Array.from(this.geneStats.values());
            }

            if (colorStatType === 'signed_lfc') {
                const lfcValues = stats.map(s => s.lfc).filter(v => v !== null && !isNaN(v));
                const minLfc = Math.min(...lfcValues);
                const maxLfc = Math.max(...lfcValues);
                const maxAbs = Math.max(Math.abs(minLfc), Math.abs(maxLfc));

                this.networkData.nodes.forEach(node => {
                    const geneStat = this.geneStats.get(node.id);
                    let bgColor = '#cccccc'; // Gray for missing

                    if (geneStat && geneStat.lfc !== null && !isNaN(geneStat.lfc)) {
                        // Blue (-) to White (0) to Red (+)
                        const normalized = (geneStat.lfc + maxAbs) / (2 * maxAbs);
                        bgColor = this.interpolateColor('#2166ac', '#f7f7f7', '#b2182b', normalized);
                    }

                    updates.push({
                        id: node.id,
                        color: { background: bgColor, border: '#ffffff' }
                    });
                });

                const scaleLabel = colorScale === 'network' ? ' (network)' : ' (all genes)';
                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">Signed LFC${scaleLabel}</div>
                    <div style="display: flex; align-items: center; gap: 4px;">
                        <span style="font-size: 10px;">${minLfc.toFixed(1)}</span>
                        <div style="width: 80px; height: 12px; background: linear-gradient(to right, #2166ac, #f7f7f7, #b2182b); border-radius: 2px;"></div>
                        <span style="font-size: 10px;">${maxLfc.toFixed(1)}</span>
                    </div>
                    <div class="legend-item" style="margin-top: 4px;">Missing: <span style="display: inline-block; width: 12px; height: 12px; background: #cccccc; border-radius: 2px; vertical-align: middle;"></span></div>
                `;
            } else if (colorStatType === 'abs_lfc') {
                const lfcValues = stats.map(s => Math.abs(s.lfc)).filter(v => v !== null && !isNaN(v));
                const maxLfc = Math.max(...lfcValues);

                this.networkData.nodes.forEach(node => {
                    const geneStat = this.geneStats.get(node.id);
                    let bgColor = '#cccccc';

                    if (geneStat && geneStat.lfc !== null && !isNaN(geneStat.lfc)) {
                        const normalized = Math.abs(geneStat.lfc) / maxLfc;
                        bgColor = this.interpolateColor('#f5f5f5', '#fdae61', '#d7191c', normalized);
                    }

                    updates.push({
                        id: node.id,
                        color: { background: bgColor, border: '#ffffff' }
                    });
                });

                const scaleLabel = colorScale === 'network' ? ' (network)' : ' (all genes)';
                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">|LFC|${scaleLabel}</div>
                    <div style="display: flex; align-items: center; gap: 4px;">
                        <span style="font-size: 10px;">0</span>
                        <div style="width: 80px; height: 12px; background: linear-gradient(to right, #f5f5f5, #fdae61, #d7191c); border-radius: 2px;"></div>
                        <span style="font-size: 10px;">${maxLfc.toFixed(1)}</span>
                    </div>
                    <div class="legend-item" style="margin-top: 4px;">Missing: <span style="display: inline-block; width: 12px; height: 12px; background: #cccccc; border-radius: 2px; vertical-align: middle;"></span></div>
                `;
            } else if (colorStatType === 'fdr') {
                const fdrValues = stats.map(s => s.fdr).filter(v => v !== null && !isNaN(v) && v > 0);
                const minFdr = Math.min(...fdrValues);

                this.networkData.nodes.forEach(node => {
                    const geneStat = this.geneStats.get(node.id);
                    let bgColor = '#cccccc';

                    if (geneStat && geneStat.fdr !== null && !isNaN(geneStat.fdr) && geneStat.fdr > 0) {
                        // Log scale: smaller FDR = more significant = redder
                        const logMin = Math.log10(minFdr);
                        const logVal = Math.log10(geneStat.fdr);
                        const normalized = 1 - Math.min(1, Math.max(0, (logVal - logMin) / (0 - logMin)));
                        bgColor = this.interpolateColor('#f5f5f5', '#fdae61', '#d7191c', normalized);
                    }

                    updates.push({
                        id: node.id,
                        color: { background: bgColor, border: '#ffffff' }
                    });
                });

                const scaleLabel = colorScale === 'network' ? ' (network)' : ' (all genes)';
                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">FDR${scaleLabel}</div>
                    <div style="display: flex; align-items: center; gap: 4px;">
                        <span style="font-size: 10px;">${minFdr.toExponential(1)}</span>
                        <div style="width: 80px; height: 12px; background: linear-gradient(to right, #d7191c, #fdae61, #f5f5f5); border-radius: 2px;"></div>
                        <span style="font-size: 10px;">1</span>
                    </div>
                    <div class="legend-item" style="margin-top: 4px;">Missing: <span style="display: inline-block; width: 12px; height: 12px; background: #cccccc; border-radius: 2px; vertical-align: middle;"></span></div>
                `;
            }
        }

        this.networkData.nodes.update(updates);
    }

    interpolateColor(color1, color2, color3, t) {
        // Parse hex colors
        const parseHex = (hex) => {
            const r = parseInt(hex.slice(1, 3), 16);
            const g = parseInt(hex.slice(3, 5), 16);
            const b = parseInt(hex.slice(5, 7), 16);
            return [r, g, b];
        };

        const [r1, g1, b1] = parseHex(color1);
        const [r2, g2, b2] = parseHex(color2);
        const [r3, g3, b3] = parseHex(color3);

        let r, g, b;
        if (t < 0.5) {
            const t2 = t * 2;
            r = Math.round(r1 + (r2 - r1) * t2);
            g = Math.round(g1 + (g2 - g1) * t2);
            b = Math.round(b1 + (b2 - b1) * t2);
        } else {
            const t2 = (t - 0.5) * 2;
            r = Math.round(r2 + (r3 - r2) * t2);
            g = Math.round(g2 + (g3 - g2) * t2);
            b = Math.round(b2 + (b3 - b2) * t2);
        }

        return `rgb(${r},${g},${b})`;
    }

    downloadAllData() {
        if (!this.results) return;

        // Create correlations CSV
        let correlationsCSV = 'Gene1,Gene2,Correlation,Slope,N,Cluster\n';
        this.results.correlations.forEach(c => {
            correlationsCSV += `${c.gene1},${c.gene2},${c.correlation},${c.slope},${c.n},${c.cluster}\n`;
        });

        // Create clusters CSV
        let clustersCSV = 'Gene,Cluster,Mean_Effect,SD_Effect\n';
        this.results.clusters.forEach(c => {
            clustersCSV += `${c.gene},${c.cluster},${c.meanEffect},${c.sdEffect}\n`;
        });

        // Create summary
        const lineage = document.getElementById('lineageFilter').value || 'All lineages';
        const summary = `Gene Correlation Analysis Summary
================================
Analysis Mode: ${this.results.mode === 'analysis' ? 'Analysis (within gene list)' : 'Design (find correlated genes)'}
Correlation Cutoff: ${this.results.cutoff}
Minimum Cell Lines: ${document.getElementById('minCellLines').value}
Minimum Slope: ${document.getElementById('minSlope').value}
Lineage Filter: ${lineage}

Input Genes: ${this.results.geneList.length}
${this.results.geneList.join(', ')}

Results:
- Total correlations found: ${this.results.correlations.length}
- Genes in network: ${this.results.clusters.length}
- Number of clusters: ${Math.max(...this.results.correlations.map(c => c.cluster))}
`;

        // For simplicity, download as separate files (ZIP would require a library)
        // Download correlations
        this.downloadFile(correlationsCSV, 'correlations.csv', 'text/csv');

        // Small delay then download clusters
        setTimeout(() => {
            this.downloadFile(clustersCSV, 'clusters.csv', 'text/csv');
        }, 500);

        // Small delay then download summary
        setTimeout(() => {
            this.downloadFile(summary, 'summary.txt', 'text/plain');
        }, 1000);
    }

    // Inspect Modal
    openInspectByGenes(gene1, gene2) {
        // Find the correlation entry by gene names
        const c = this.results.correlations.find(corr =>
            (corr.gene1 === gene1 && corr.gene2 === gene2) ||
            (corr.gene1 === gene2 && corr.gene2 === gene1)
        );
        if (!c) {
            console.error('Correlation not found for', gene1, gene2);
            return;
        }
        this.openInspect(c);
    }

    openByTissueByGenes(gene1, gene2) {
        // Find the correlation entry by gene names
        const c = this.results.correlations.find(corr =>
            (corr.gene1 === gene1 && corr.gene2 === gene2) ||
            (corr.gene1 === gene2 && corr.gene2 === gene1)
        );
        if (!c) {
            console.error('Correlation not found for', gene1, gene2);
            return;
        }

        // Set up currentInspect with the data needed for By tissue
        this.currentInspect = {
            gene1: c.gene1,
            gene2: c.gene2,
            correlation: c.correlation
        };

        // Get data for both genes
        const idx1 = this.geneIndex.get(c.gene1);
        const idx2 = this.geneIndex.get(c.gene2);
        const data1 = this.getGeneData(idx1);
        const data2 = this.getGeneData(idx2);

        // Prepare plot data
        const plotData = [];
        for (let i = 0; i < this.nCellLines; i++) {
            if (!isNaN(data1[i]) && !isNaN(data2[i])) {
                const cellLine = this.metadata.cellLines[i];
                plotData.push({
                    x: data1[i],
                    y: data2[i],
                    cellLineId: cellLine,
                    cellLineName: this.getCellLineName(cellLine),
                    lineage: this.getCellLineLineage(cellLine)
                });
            }
        }
        this.currentInspect.data = plotData;

        // Open the inspect modal and show By tissue view
        document.getElementById('inspectModal').classList.add('active');
        document.getElementById('inspectTitle').textContent =
            `${c.gene1} vs ${c.gene2} - By Tissue Breakdown`;

        // Hide the scatter plot controls (not needed for By tissue view)
        document.querySelector('.inspect-controls').style.display = 'none';

        // Hide scatter-specific download buttons
        document.getElementById('downloadScatterPNG').style.display = 'none';
        document.getElementById('downloadScatterSVG').style.display = 'none';
        document.getElementById('downloadScatterCSV').style.display = 'none';

        // Show tissue-specific download buttons
        document.getElementById('downloadTissuePNG').style.display = '';
        document.getElementById('downloadTissueSVG').style.display = '';
        document.getElementById('downloadTissueCSV').style.display = '';

        this.showByTissueModal();
    }

    openInspect(c) {
        // c is now the correlation object directly
        this.currentInspect = {
            gene1: c.gene1,
            gene2: c.gene2,
            correlation: c.correlation
        };
        this.clickedCells.clear();

        // Get data for both genes
        const idx1 = this.geneIndex.get(c.gene1);
        const idx2 = this.geneIndex.get(c.gene2);
        const data1 = this.getGeneData(idx1);
        const data2 = this.getGeneData(idx2);

        // Prepare plot data
        const plotData = [];
        for (let i = 0; i < this.nCellLines; i++) {
            if (!isNaN(data1[i]) && !isNaN(data2[i])) {
                const cellLine = this.metadata.cellLines[i];
                plotData.push({
                    x: data1[i],
                    y: data2[i],
                    cellLineId: cellLine,
                    cellLineName: this.getCellLineName(cellLine),
                    lineage: this.getCellLineLineage(cellLine)
                });
            }
        }
        this.currentInspect.data = plotData;

        // Set axis limits
        const xVals = plotData.map(d => d.x);
        const yVals = plotData.map(d => d.y);
        this.currentInspect.defaultXlim = [Math.min(...xVals), Math.max(...xVals)];
        this.currentInspect.defaultYlim = [Math.min(...yVals), Math.max(...yVals)];

        document.getElementById('scatterXmin').value = this.currentInspect.defaultXlim[0].toFixed(1);
        document.getElementById('scatterXmax').value = this.currentInspect.defaultXlim[1].toFixed(1);
        document.getElementById('scatterYmin').value = this.currentInspect.defaultYlim[0].toFixed(1);
        document.getElementById('scatterYmax').value = this.currentInspect.defaultYlim[1].toFixed(1);

        // Populate cancer filter
        const cancerFilter = document.getElementById('scatterCancerFilter');
        const cancerBox = document.getElementById('cancerFilterBox');
        const lineages = [...new Set(plotData.map(d => d.lineage).filter(l => l))].sort();
        if (lineages.length > 0) {
            cancerFilter.innerHTML = '<option value="">All cancer types</option>';
            lineages.forEach(l => {
                cancerFilter.innerHTML += `<option value="${l}">${l}</option>`;
            });
            cancerBox.style.display = 'block';
        } else {
            cancerBox.style.display = 'none';
        }

        // Populate hotspot genes (excluding HLA-A and HLA-B which have high variability)
        const hotspotSelect = document.getElementById('hotspotGene');
        const mutFilterGeneSelect = document.getElementById('mutationFilterGene');
        const excludedGenes = ['HLA-A', 'HLA-B'];
        if (this.mutations?.genes?.length > 0) {
            hotspotSelect.innerHTML = '<option value="">Select gene...</option>';
            mutFilterGeneSelect.innerHTML = '<option value="">No filter</option>';
            this.mutations.genes
                .filter(g => !excludedGenes.includes(g))
                .forEach(g => {
                    const count = this.mutations.geneCounts?.[g] || 0;
                    hotspotSelect.innerHTML += `<option value="${g}">${g} (${count} mut)</option>`;
                    mutFilterGeneSelect.innerHTML += `<option value="${g}">${g} (${count} mut)</option>`;
                });
            document.getElementById('mutationBox').style.display = 'block';
            document.getElementById('mutationFilterBox').style.display = 'block';
        } else {
            document.getElementById('mutationBox').style.display = 'none';
            document.getElementById('mutationFilterBox').style.display = 'none';
        }

        // Update title with slope
        const slope = this.pearsonWithSlope(plotData.map(d => d.x), plotData.map(d => d.y)).slope;
        document.getElementById('inspectTitle').textContent =
            `Correlation: ${c.gene1} vs ${c.gene2} | r=${c.correlation.toFixed(3)}, slope=${slope.toFixed(3)} (n=${plotData.length})`;

        // Show modal and render plot
        document.getElementById('inspectModal').classList.add('active');

        // Make sure controls are visible (may have been hidden by By tissue view)
        document.querySelector('.inspect-controls').style.display = '';
        document.getElementById('downloadScatterPNG').style.display = '';
        document.getElementById('downloadScatterSVG').style.display = '';
        document.getElementById('downloadScatterCSV').style.display = '';

        // Hide tissue-specific buttons
        document.getElementById('downloadTissuePNG').style.display = 'none';
        document.getElementById('downloadTissueSVG').style.display = 'none';
        document.getElementById('downloadTissueCSV').style.display = 'none';

        this.updateInspectPlot();
    }

    closeInspectModal() {
        document.getElementById('inspectModal').classList.remove('active');
        this.currentInspect = null;
    }

    resetInspectAxes() {
        if (!this.currentInspect) return;
        document.getElementById('scatterXmin').value = this.currentInspect.defaultXlim[0].toFixed(1);
        document.getElementById('scatterXmax').value = this.currentInspect.defaultXlim[1].toFixed(1);
        document.getElementById('scatterYmin').value = this.currentInspect.defaultYlim[0].toFixed(1);
        document.getElementById('scatterYmax').value = this.currentInspect.defaultYlim[1].toFixed(1);
        this.updateInspectPlot();
    }

    updateInspectPlot() {
        if (!this.currentInspect) return;

        const data = this.currentInspect.data;
        const gene1 = this.currentInspect.gene1;
        const gene2 = this.currentInspect.gene2;

        // Get filter settings
        const cancerFilter = document.getElementById('scatterCancerFilter').value;
        const mutFilterGene = document.getElementById('mutationFilterGene').value;
        const mutFilterLevel = document.getElementById('mutationFilterLevel').value;
        const searchTerms = document.getElementById('scatterCellSearch').value
            .split('\n').map(s => s.trim().toUpperCase()).filter(s => s);
        const fontSize = parseInt(document.getElementById('scatterFontSize')?.value) || 3;
        const hotspotGene = document.getElementById('hotspotGene').value;
        const hotspotMode = document.getElementById('hotspotMode').value;

        // Filter by cancer type
        let filteredData = cancerFilter ?
            data.filter(d => d.lineage === cancerFilter) : data;

        // Apply mutation filter (separate from overlay)
        if (mutFilterGene && this.mutations?.geneData?.[mutFilterGene] && mutFilterLevel !== 'all') {
            const filterMutations = this.mutations.geneData[mutFilterGene].mutations;
            filteredData = filteredData.filter(d => {
                const mutLevel = filterMutations[d.cellLineId] || 0;
                if (mutFilterLevel === '0') return mutLevel === 0;
                if (mutFilterLevel === '1+') return mutLevel >= 1;
                if (mutFilterLevel === '2+') return mutLevel >= 2;
                return true;
            });
        }

        // Get mutation info for overlay (separate gene)
        let mutationMap = new Map();
        if (hotspotGene && this.mutations?.geneData?.[hotspotGene]) {
            const geneData = this.mutations.geneData[hotspotGene];
            Object.entries(geneData.mutations).forEach(([cellLine, mutLevel]) => {
                mutationMap.set(cellLine, mutLevel);
            });
        }

        // Add mutation level to filtered data (for overlay coloring)
        filteredData = filteredData.map(d => ({
            ...d,
            mutationLevel: mutationMap.get(d.cellLineId) || 0
        }));

        // Show/hide plot and table based on mode
        const scatterPlot = document.getElementById('scatterPlot');
        const compareTable = document.getElementById('compareTable');

        if (hotspotMode === 'compare_table' && hotspotGene) {
            scatterPlot.style.display = 'none';
            compareTable.style.display = 'block';
            this.renderCompareTable(filteredData, gene1, gene2, hotspotGene);
            return;
        } else {
            scatterPlot.style.display = 'block';
            compareTable.style.display = 'none';
        }

        // Handle 3-panel mode
        if (hotspotMode === 'three_panel' && hotspotGene) {
            this.renderThreePanelPlot(filteredData, gene1, gene2, hotspotGene, searchTerms, fontSize);
            return;
        }

        // Single panel color mode
        this.renderSinglePanelPlot(filteredData, gene1, gene2, hotspotGene, hotspotMode, searchTerms, fontSize);
    }

    renderSinglePanelPlot(filteredData, gene1, gene2, hotspotGene, hotspotMode, searchTerms, fontSize) {
        // Calculate stats for each mutation group
        const wt = filteredData.filter(d => d.mutationLevel === 0);
        const mut1 = filteredData.filter(d => d.mutationLevel === 1);
        const mut2 = filteredData.filter(d => d.mutationLevel >= 2);

        const wtStats = this.pearsonWithSlope(wt.map(d => d.x), wt.map(d => d.y));
        const mut1Stats = this.pearsonWithSlope(mut1.map(d => d.x), mut1.map(d => d.y));
        const mut2Stats = this.pearsonWithSlope(mut2.map(d => d.x), mut2.map(d => d.y));
        const allStats = this.pearsonWithSlope(filteredData.map(d => d.x), filteredData.map(d => d.y));

        // Build traces
        const traces = [];
        const xRange = [parseFloat(document.getElementById('scatterXmin').value),
                       parseFloat(document.getElementById('scatterXmax').value)];
        const yRange = [parseFloat(document.getElementById('scatterYmin').value),
                       parseFloat(document.getElementById('scatterYmax').value)];

        if (hotspotMode === 'color' && hotspotGene) {
            // Color by mutation (0/1/2) mode with separate traces for legend
            const wtPct = (wt.length / filteredData.length * 100).toFixed(1);
            const mut1Pct = (mut1.length / filteredData.length * 100).toFixed(1);
            const mut2Pct = (mut2.length / filteredData.length * 100).toFixed(1);

            // WT trace (gray)
            traces.push({
                x: wt.map(d => d.x),
                y: wt.map(d => d.y),
                mode: 'markers',
                type: 'scatter',
                text: wt.map(d => `${d.cellLineName}<br>${d.lineage}<br>WT`),
                hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                marker: { color: '#9ca3af', size: 8, opacity: 0.6 },
                name: `WT (n=${wt.length}, ${wtPct}%)`
            });

            // 1 mut trace (blue)
            traces.push({
                x: mut1.map(d => d.x),
                y: mut1.map(d => d.y),
                mode: 'markers',
                type: 'scatter',
                text: mut1.map(d => `${d.cellLineName}<br>${d.lineage}<br>1 mutation`),
                hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                marker: { color: '#3b82f6', size: 10, opacity: 0.7 },
                name: `1 mut (n=${mut1.length}, ${mut1Pct}%)`
            });

            // 2 mut trace (red)
            traces.push({
                x: mut2.map(d => d.x),
                y: mut2.map(d => d.y),
                mode: 'markers',
                type: 'scatter',
                text: mut2.map(d => `${d.cellLineName}<br>${d.lineage}<br>2 mutations`),
                hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                marker: { color: '#dc2626', size: 11, opacity: 0.8 },
                name: `2 mut (n=${mut2.length}, ${mut2Pct}%)`
            });
        } else {
            // Default mode - all same color
            traces.push({
                x: filteredData.map(d => d.x),
                y: filteredData.map(d => d.y),
                mode: 'markers',
                type: 'scatter',
                text: filteredData.map(d => `${d.cellLineName}<br>${d.lineage}`),
                hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                marker: { color: '#3b82f6', size: 8, opacity: 0.7 },
                name: 'Cell lines',
                showlegend: false
            });
        }

        // Add highlights
        const highlightData = filteredData.filter(d =>
            searchTerms.some(term =>
                d.cellLineName.toUpperCase().includes(term) ||
                d.cellLineId.toUpperCase().includes(term)
            ) || this.clickedCells.has(d.cellLineName)
        );

        if (highlightData.length > 0) {
            traces.push({
                x: highlightData.map(d => d.x),
                y: highlightData.map(d => d.y),
                mode: 'markers+text',
                type: 'scatter',
                text: highlightData.map(d => `${d.cellLineName} (${d.lineage || 'Unknown'})`),
                textposition: 'top center',
                textfont: { size: fontSize * 3, color: '#000' },
                hovertemplate: '%{text}<extra></extra>',
                marker: {
                    color: '#f59e0b',
                    size: 12,
                    symbol: 'circle',
                    line: { color: '#000', width: 2 }
                },
                name: 'Highlighted',
                showlegend: false
            });
        }

        // Add regression line
        if (!isNaN(allStats.slope)) {
            const meanX = filteredData.reduce((a, d) => a + d.x, 0) / filteredData.length;
            const meanY = filteredData.reduce((a, d) => a + d.y, 0) / filteredData.length;
            const intercept = meanY - allStats.slope * meanX;

            traces.push({
                x: xRange,
                y: [allStats.slope * xRange[0] + intercept, allStats.slope * xRange[1] + intercept],
                mode: 'lines',
                type: 'scatter',
                line: { color: '#16a34a', width: 3 },
                fill: 'none',
                name: 'Regression',
                showlegend: false
            });
        }

        // Build title and annotations
        let titleText = `<b>${gene1} vs ${gene2}</b>`;
        let annotations = [];

        if (hotspotMode === 'color' && hotspotGene) {
            // Add mutation stats annotation - positioned above the plot
            const annotationLines = [
                `<b>Overall:</b> r=${allStats.correlation.toFixed(3)}, slope=${allStats.slope.toFixed(3)}, n=${filteredData.length}`,
                `<b>${hotspotGene}:</b> WT: n=${wt.length}, r=${wtStats.correlation.toFixed(3)}, slope=${wtStats.slope.toFixed(3)} | ` +
                `1mut: n=${mut1.length}, r=${mut1Stats.correlation.toFixed(3)}, slope=${mut1Stats.slope.toFixed(3)} | ` +
                `2mut: n=${mut2.length}, r=${mut2Stats.correlation.toFixed(3)}, slope=${mut2Stats.slope.toFixed(3)}`
            ];

            annotations.push({
                x: 0,
                y: 1.08,
                xref: 'paper',
                yref: 'paper',
                text: annotationLines.join('<br>'),
                showarrow: false,
                font: { size: 10, color: '#333' },
                align: 'left',
                xanchor: 'left'
            });
        }

        const layout = {
            title: {
                text: titleText,
                x: 0.5,
                y: 0.99,
                font: { size: 14 }
            },
            xaxis: {
                title: `${gene1} CRISPR Effect`,
                range: xRange,
                zeroline: true,
                zerolinecolor: '#ddd',
                constrain: 'domain'
            },
            yaxis: {
                title: `${gene2} CRISPR Effect`,
                range: yRange,
                zeroline: true,
                zerolinecolor: '#ddd',
                scaleanchor: 'x',
                scaleratio: parseFloat(document.getElementById('aspectRatio')?.value || 1),
                constrain: 'domain'
            },
            hovermode: 'closest',
            margin: { t: (hotspotMode === 'color' && hotspotGene) ? 110 : 50, r: 150, b: 60, l: 60 },
            showlegend: hotspotMode === 'color' && hotspotGene,
            legend: {
                x: 1.02,
                y: 0.5,
                xanchor: 'left',
                yanchor: 'middle',
                title: { text: hotspotGene, font: { size: 12 } }
            },
            annotations: annotations,
            plot_bgcolor: '#fafafa'
        };

        Plotly.newPlot('scatterPlot', traces, layout, { responsive: true });

        // Add click handler
        this.setupScatterClickHandler(filteredData);
    }

    renderThreePanelPlot(filteredData, gene1, gene2, hotspotGene, searchTerms, fontSize) {
        const wt = filteredData.filter(d => d.mutationLevel === 0);
        const mut1 = filteredData.filter(d => d.mutationLevel === 1);
        const mut2 = filteredData.filter(d => d.mutationLevel >= 2);

        const xRange = [parseFloat(document.getElementById('scatterXmin').value),
                       parseFloat(document.getElementById('scatterXmax').value)];
        const yRange = [parseFloat(document.getElementById('scatterYmin').value),
                       parseFloat(document.getElementById('scatterYmax').value)];

        const wtStats = this.pearsonWithSlope(wt.map(d => d.x), wt.map(d => d.y));
        const mut1Stats = this.pearsonWithSlope(mut1.map(d => d.x), mut1.map(d => d.y));
        const mut2Stats = this.pearsonWithSlope(mut2.map(d => d.x), mut2.map(d => d.y));

        // Calculate means and medians for each group
        const calcGroupStats = (data) => ({
            meanX: data.length > 0 ? data.reduce((a, d) => a + d.x, 0) / data.length : NaN,
            meanY: data.length > 0 ? data.reduce((a, d) => a + d.y, 0) / data.length : NaN,
            medianX: this.median(data.map(d => d.x)),
            medianY: this.median(data.map(d => d.y))
        });
        const wtExtra = calcGroupStats(wt);
        const mut1Extra = calcGroupStats(mut1);
        const mut2Extra = calcGroupStats(mut2);

        const traces = [];

        // Panel 1: WT (0 mut)
        traces.push({
            x: wt.map(d => d.x),
            y: wt.map(d => d.y),
            xaxis: 'x',
            yaxis: 'y',
            mode: 'markers',
            type: 'scatter',
            text: wt.map(d => `${d.cellLineName}<br>${d.lineage}`),
            hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
            marker: { color: '#9ca3af', size: 7, opacity: 0.6 },
            name: 'WT',
            showlegend: false
        });

        // Panel 2: 1 mutation
        traces.push({
            x: mut1.map(d => d.x),
            y: mut1.map(d => d.y),
            xaxis: 'x2',
            yaxis: 'y2',
            mode: 'markers',
            type: 'scatter',
            text: mut1.map(d => `${d.cellLineName}<br>${d.lineage}`),
            hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
            marker: { color: '#3b82f6', size: 8, opacity: 0.7 },
            name: '1 mut',
            showlegend: false
        });

        // Panel 3: 2 mutations
        traces.push({
            x: mut2.map(d => d.x),
            y: mut2.map(d => d.y),
            xaxis: 'x3',
            yaxis: 'y3',
            mode: 'markers',
            type: 'scatter',
            text: mut2.map(d => `${d.cellLineName}<br>${d.lineage}`),
            hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
            marker: { color: '#dc2626', size: 8, opacity: 0.7 },
            name: '2 mut',
            showlegend: false
        });

        // Add regression lines for each panel
        const addRegressionLine = (data, stats, xaxis, yaxis, color) => {
            if (data.length >= 3 && !isNaN(stats.slope)) {
                const meanX = data.reduce((a, d) => a + d.x, 0) / data.length;
                const meanY = data.reduce((a, d) => a + d.y, 0) / data.length;
                const intercept = meanY - stats.slope * meanX;
                traces.push({
                    x: xRange,
                    y: [stats.slope * xRange[0] + intercept, stats.slope * xRange[1] + intercept],
                    xaxis: xaxis,
                    yaxis: yaxis,
                    mode: 'lines',
                    type: 'scatter',
                    line: { color: color, width: 2 },
                    showlegend: false
                });
            }
        };

        addRegressionLine(wt, wtStats, 'x', 'y', '#16a34a');
        addRegressionLine(mut1, mut1Stats, 'x2', 'y2', '#16a34a');
        addRegressionLine(mut2, mut2Stats, 'x3', 'y3', '#16a34a');

        // Add highlighted cells for each panel
        const addHighlights = (data, xaxis, yaxis) => {
            const highlightData = data.filter(d =>
                searchTerms.some(term =>
                    d.cellLineName.toUpperCase().includes(term) ||
                    d.cellLineId.toUpperCase().includes(term)
                ) || this.clickedCells.has(d.cellLineName)
            );

            if (highlightData.length > 0) {
                traces.push({
                    x: highlightData.map(d => d.x),
                    y: highlightData.map(d => d.y),
                    xaxis: xaxis,
                    yaxis: yaxis,
                    mode: 'markers+text',
                    type: 'scatter',
                    text: highlightData.map(d => `${d.cellLineName} (${d.lineage || 'Unknown'})`),
                    textposition: 'top center',
                    textfont: { size: fontSize * 3, color: '#000' },
                    hovertemplate: '%{text}<extra></extra>',
                    marker: {
                        color: '#f59e0b',
                        size: 10,
                        symbol: 'circle',
                        line: { color: '#000', width: 2 }
                    },
                    showlegend: false
                });
            }
        };

        addHighlights(wt, 'x', 'y');
        addHighlights(mut1, 'x2', 'y2');
        addHighlights(mut2, 'x3', 'y3');

        const layout = {
            title: {
                text: `<b>${gene1} vs ${gene2} - ${hotspotGene} hotspot mutation stratification</b>`,
                x: 0.5,
                y: 0.98,
                font: { size: 14 }
            },
            grid: { rows: 1, columns: 3, pattern: 'independent' },
            xaxis: { title: `${gene1} Effect`, range: xRange, domain: [0, 0.28], constrain: 'domain' },
            yaxis: {
                title: `${gene2} Effect`, range: yRange,
                scaleanchor: 'x',
                scaleratio: parseFloat(document.getElementById('aspectRatio')?.value || 1),
                constrain: 'domain'
            },
            xaxis2: { title: `${gene1} Effect`, range: xRange, domain: [0.36, 0.64], constrain: 'domain' },
            yaxis2: {
                range: yRange, anchor: 'x2',
                scaleanchor: 'x2',
                scaleratio: parseFloat(document.getElementById('aspectRatio')?.value || 1),
                constrain: 'domain'
            },
            xaxis3: { title: `${gene1} Effect`, range: xRange, domain: [0.72, 1], constrain: 'domain' },
            yaxis3: {
                range: yRange, anchor: 'x3',
                scaleanchor: 'x3',
                scaleratio: parseFloat(document.getElementById('aspectRatio')?.value || 1),
                constrain: 'domain'
            },
            annotations: [
                { x: 0.14, y: 1.02, xref: 'paper', yref: 'paper',
                  text: `<b>WT (0 mut)</b> n=${wt.length}<br>r=${wtStats.correlation.toFixed(2)}, slope=${wtStats.slope.toFixed(2)}<br>mean: x=${wtExtra.meanX.toFixed(2)}, y=${wtExtra.meanY.toFixed(2)}<br>median: x=${wtExtra.medianX.toFixed(2)}, y=${wtExtra.medianY.toFixed(2)}`,
                  showarrow: false, font: { size: 9 } },
                { x: 0.5, y: 1.02, xref: 'paper', yref: 'paper',
                  text: `<b>1 mutation</b> n=${mut1.length}<br>r=${mut1Stats.correlation.toFixed(2)}, slope=${mut1Stats.slope.toFixed(2)}<br>mean: x=${mut1Extra.meanX.toFixed(2)}, y=${mut1Extra.meanY.toFixed(2)}<br>median: x=${mut1Extra.medianX.toFixed(2)}, y=${mut1Extra.medianY.toFixed(2)}`,
                  showarrow: false, font: { size: 9 } },
                { x: 0.86, y: 1.02, xref: 'paper', yref: 'paper',
                  text: `<b>2 mutations</b> n=${mut2.length}<br>r=${mut2Stats.correlation.toFixed(2)}, slope=${mut2Stats.slope.toFixed(2)}<br>mean: x=${mut2Extra.meanX.toFixed(2)}, y=${mut2Extra.meanY.toFixed(2)}<br>median: x=${mut2Extra.medianX.toFixed(2)}, y=${mut2Extra.medianY.toFixed(2)}`,
                  showarrow: false, font: { size: 9 } }
            ],
            margin: { t: 140, r: 30, b: 60, l: 60 },
            plot_bgcolor: '#fafafa'
        };

        Plotly.newPlot('scatterPlot', traces, layout, { responsive: true });
    }

    renderCompareTable(filteredData, gene1, gene2, hotspotGene) {
        // Group by cancer type (lineage) - comparing 0 vs 2 mutations only
        const lineageGroups = {};
        filteredData.forEach(d => {
            if (!d.lineage) return;
            if (!lineageGroups[d.lineage]) {
                lineageGroups[d.lineage] = { wt: [], mut: [] };
            }
            if (d.mutationLevel === 0) {
                lineageGroups[d.lineage].wt.push(d);
            } else if (d.mutationLevel >= 2) {
                lineageGroups[d.lineage].mut.push(d);
            }
            // Note: mutationLevel === 1 is excluded from comparison
        });

        // Calculate stats for each lineage
        const tableData = [];
        Object.entries(lineageGroups).forEach(([lineage, groups]) => {
            if (groups.wt.length >= 3 && groups.mut.length >= 3) {
                const wtStats = this.pearsonWithSlope(groups.wt.map(d => d.x), groups.wt.map(d => d.y));
                const mutStats = this.pearsonWithSlope(groups.mut.map(d => d.x), groups.mut.map(d => d.y));

                // Calculate delta and p-values using Fisher z-transformation for correlation difference
                const deltaR = mutStats.correlation - wtStats.correlation;
                const deltaSlope = mutStats.slope - wtStats.slope;

                // Fisher z-transformation for r difference p-value
                const z1 = 0.5 * Math.log((1 + wtStats.correlation) / (1 - wtStats.correlation));
                const z2 = 0.5 * Math.log((1 + mutStats.correlation) / (1 - mutStats.correlation));
                const se = Math.sqrt(1/(groups.wt.length - 3) + 1/(groups.mut.length - 3));
                const zDiff = (z2 - z1) / se;
                const pR = 2 * (1 - this.normalCDF(Math.abs(zDiff)));

                // Approximate p-value for slope difference (simplified)
                const pSlope = pR * 2; // Rough approximation

                tableData.push({
                    lineage,
                    nWT: groups.wt.length,
                    rWT: wtStats.correlation,
                    slopeWT: wtStats.slope,
                    nMut: groups.mut.length,
                    rMut: mutStats.correlation,
                    slopeMut: mutStats.slope,
                    deltaR,
                    pR,
                    deltaSlope,
                    pSlope
                });
            }
        });

        // Sort by deltaR
        tableData.sort((a, b) => a.pR - b.pR);

        // Build HTML table
        let html = `
            <h4 style="margin-bottom: 8px;">Mutation Effect on Correlation by Cancer Type</h4>
            <p style="font-size: 12px; color: #666; margin-bottom: 12px;">
                Comparing correlation between WT (0 mutations) vs Mutant (2 mutations) cells, stratified by cancer type.
                Note: Cells with exactly 1 mutation are excluded from this comparison.
            </p>
            <div style="overflow-x: auto;">
            <table class="data-table" style="width: 100%; font-size: 12px;">
                <thead>
                    <tr>
                        <th>Cancer Type</th>
                        <th>N (WT)</th>
                        <th>r (WT)</th>
                        <th>slope (WT)</th>
                        <th>N (Mut)</th>
                        <th>r (Mut)</th>
                        <th>slope (Mut)</th>
                        <th>Δr</th>
                        <th>p(Δr)</th>
                        <th>Δslope</th>
                        <th>p(Δslope)</th>
                    </tr>
                </thead>
                <tbody>
        `;

        tableData.forEach(row => {
            const deltaRColor = row.deltaR < 0 ? '#dc2626' : '#16a34a';
            const deltaSlopeColor = row.deltaSlope < 0 ? '#dc2626' : '#16a34a';

            html += `
                <tr>
                    <td>${row.lineage}</td>
                    <td>${row.nWT}</td>
                    <td>${row.rWT.toFixed(3)}</td>
                    <td>${row.slopeWT.toFixed(3)}</td>
                    <td>${row.nMut}</td>
                    <td>${row.rMut.toFixed(3)}</td>
                    <td>${row.slopeMut.toFixed(3)}</td>
                    <td style="color: ${deltaRColor}; font-weight: 600;">${row.deltaR.toFixed(3)}</td>
                    <td>${row.pR.toExponential(1)}</td>
                    <td style="color: ${deltaSlopeColor}; font-weight: 600;">${row.deltaSlope.toFixed(3)}</td>
                    <td>${row.pSlope.toExponential(1)}</td>
                </tr>
            `;
        });

        html += `
                </tbody>
            </table>
            </div>
            <div style="margin-top: 12px;">
                <button class="btn btn-success btn-sm" id="downloadCompareCSV">Download CSV</button>
            </div>
        `;

        document.getElementById('compareTable').innerHTML = html;

        // Add download handler
        document.getElementById('downloadCompareCSV')?.addEventListener('click', () => {
            let csv = 'Cancer Type,N (WT),r (WT),slope (WT),N (Mut),r (Mut),slope (Mut),Δr,p(Δr),Δslope,p(Δslope)\n';
            tableData.forEach(row => {
                csv += `"${row.lineage}",${row.nWT},${row.rWT.toFixed(4)},${row.slopeWT.toFixed(4)},${row.nMut},${row.rMut.toFixed(4)},${row.slopeMut.toFixed(4)},${row.deltaR.toFixed(4)},${row.pR.toExponential(2)},${row.deltaSlope.toFixed(4)},${row.pSlope.toExponential(2)}\n`;
            });
            this.downloadFile(csv, `${gene1}_vs_${gene2}_${hotspotGene}_mutation_comparison.csv`, 'text/csv');
        });
    }

    normalCDF(x) {
        // Approximation of standard normal CDF
        const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
        const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
        const sign = x < 0 ? -1 : 1;
        x = Math.abs(x) / Math.sqrt(2);
        const t = 1 / (1 + p * x);
        const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
        return 0.5 * (1 + sign * y);
    }

    setupScatterClickHandler(filteredData) {
        document.getElementById('scatterPlot').on('plotly_click', (eventData) => {
            if (eventData.points.length > 0) {
                const point = eventData.points[0];
                const matchingData = filteredData.find(d =>
                    Math.abs(d.x - point.x) < 0.001 && Math.abs(d.y - point.y) < 0.001
                );
                if (matchingData) {
                    if (this.clickedCells.has(matchingData.cellLineName)) {
                        this.clickedCells.delete(matchingData.cellLineName);
                    } else {
                        this.clickedCells.add(matchingData.cellLineName);
                    }
                    this.updateInspectPlot();
                }
            }
        });
    }

    downloadScatterPNG() {
        const hotspotGene = document.getElementById('hotspotGene').value;
        const hotspotMode = document.getElementById('hotspotMode').value;
        // Increase height for 3-panel to accommodate annotations without overlap
        const width = hotspotMode === 'three_panel' ? 1800 : 1000;
        const height = hotspotMode === 'three_panel' ? 800 : 1000;
        const suffix = hotspotGene && hotspotMode !== 'none' ? `_${hotspotGene}` : '';

        Plotly.downloadImage('scatterPlot', {
            format: 'png',
            width: width,
            height: height,
            filename: `scatter_${this.currentInspect.gene1}_vs_${this.currentInspect.gene2}${suffix}`
        });
    }

    downloadScatterSVG() {
        const hotspotGene = document.getElementById('hotspotGene').value;
        const hotspotMode = document.getElementById('hotspotMode').value;
        // Increase height for 3-panel to accommodate annotations without overlap
        const width = hotspotMode === 'three_panel' ? 1800 : 1000;
        const height = hotspotMode === 'three_panel' ? 800 : 1000;
        const suffix = hotspotGene && hotspotMode !== 'none' ? `_${hotspotGene}` : '';

        Plotly.downloadImage('scatterPlot', {
            format: 'svg',
            width: width,
            height: height,
            filename: `scatter_${this.currentInspect.gene1}_vs_${this.currentInspect.gene2}${suffix}`
        });
    }

    showByTissueModal() {
        if (!this.currentInspect) return;

        const { gene1, gene2, data } = this.currentInspect;

        // Group data by lineage
        const tissueGroups = {};
        data.forEach(d => {
            const lineage = d.lineage || 'Unknown';
            if (!tissueGroups[lineage]) {
                tissueGroups[lineage] = [];
            }
            tissueGroups[lineage].push(d);
        });

        // Calculate stats for each tissue
        const tissueStats = [];
        Object.entries(tissueGroups).forEach(([tissue, points]) => {
            if (points.length >= 3) {
                const xVals = points.map(d => d.x);
                const yVals = points.map(d => d.y);
                const stats = this.pearsonWithSlope(xVals, yVals);
                const meanX = xVals.reduce((a, b) => a + b, 0) / xVals.length;
                const meanY = yVals.reduce((a, b) => a + b, 0) / yVals.length;
                const sdX = Math.sqrt(xVals.reduce((a, b) => a + Math.pow(b - meanX, 2), 0) / xVals.length);
                const sdY = Math.sqrt(yVals.reduce((a, b) => a + Math.pow(b - meanY, 2), 0) / yVals.length);

                tissueStats.push({
                    tissue,
                    n: points.length,
                    correlation: stats.correlation,
                    meanX,
                    sdX,
                    meanY,
                    sdY
                });
            }
        });

        // Sort by correlation (highest first)
        tissueStats.sort((a, b) => b.correlation - a.correlation);

        // Store for CSV download
        this.currentTissueStats = tissueStats;

        if (tissueStats.length === 0) {
            alert('No tissue data available (need at least 3 samples per tissue)');
            return;
        }

        // Create horizontal bar chart
        const barColors = tissueStats.map(t => {
            if (t.correlation > 0) {
                // Green gradient for positive
                const intensity = Math.min(1, t.correlation);
                return `rgba(34, 197, 94, ${0.3 + intensity * 0.7})`;
            } else {
                // Red gradient for negative
                const intensity = Math.min(1, Math.abs(t.correlation));
                return `rgba(239, 68, 68, ${0.3 + intensity * 0.7})`;
            }
        });

        // Put n= text at the end of positive bars (right side), shown in hover for negative
        const trace = {
            type: 'bar',
            orientation: 'h',
            y: tissueStats.map(t => t.tissue),
            x: tissueStats.map(t => t.correlation),
            text: tissueStats.map(t => `n=${t.n}`),
            textposition: tissueStats.map(t => t.correlation >= 0 ? 'outside' : 'inside'),
            textfont: { size: 10 },
            insidetextanchor: 'start',
            marker: { color: barColors },
            hovertemplate: '%{y}<br>r=%{x:.2f}<br>n=%{text}<extra></extra>',
            cliponaxis: false
        };

        // Find longest tissue name for dynamic left margin
        const maxTissueLen = Math.max(...tissueStats.map(t => t.tissue.length));
        const leftMargin = Math.max(150, maxTissueLen * 7);

        const layout = {
            title: {
                text: `<b>${gene1} vs ${gene2}</b>`,
                font: { size: 14 }
            },
            xaxis: {
                title: 'Correlation',
                range: [-1.15, 1.15],
                zeroline: true,
                zerolinecolor: '#000',
                zerolinewidth: 1
            },
            yaxis: {
                automargin: true,
                tickfont: { size: 9 }
            },
            margin: { t: 50, r: 50, b: 50, l: leftMargin },
            plot_bgcolor: '#fafafa',
            showlegend: false,
            height: Math.max(400, tissueStats.length * 25 + 100)
        };

        // Build statistics table HTML
        let tableHtml = `
            <h4 style="margin: 0 0 10px 0;">Statistics by Lineage</h4>
            <div style="max-height: 500px; overflow-y: auto;">
            <table style="width: 100%; border-collapse: collapse; font-size: 11px;">
                <thead>
                    <tr style="background-color: #22c55e; color: white;">
                        <th style="padding: 6px; border: 1px solid #16a34a; text-align: left; position: sticky; top: 0; background-color: #22c55e;">Lineage</th>
                        <th style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e;">N</th>
                        <th style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e;">Corr</th>
                        <th style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e;">${gene1} (mean)</th>
                        <th style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e;">${gene1} (SD)</th>
                        <th style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e;">${gene2} (mean)</th>
                        <th style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e;">${gene2} (SD)</th>
                    </tr>
                </thead>
                <tbody>
        `;

        tissueStats.forEach(t => {
            const corrColor = t.correlation > 0 ?
                `rgba(34, 197, 94, ${Math.min(1, Math.abs(t.correlation))})` :
                `rgba(239, 68, 68, ${Math.min(1, Math.abs(t.correlation))})`;
            tableHtml += `
                <tr>
                    <td style="padding: 5px; border: 1px solid #ddd;">${t.tissue}</td>
                    <td style="padding: 5px; border: 1px solid #ddd; text-align: center;">${t.n}</td>
                    <td style="padding: 5px; border: 1px solid #ddd; text-align: center; background: ${corrColor}; color: ${Math.abs(t.correlation) > 0.5 ? 'white' : 'black'}">${t.correlation.toFixed(2)}</td>
                    <td style="padding: 5px; border: 1px solid #ddd; text-align: center;">${t.meanX.toFixed(2)}</td>
                    <td style="padding: 5px; border: 1px solid #ddd; text-align: center;">${t.sdX.toFixed(2)}</td>
                    <td style="padding: 5px; border: 1px solid #ddd; text-align: center;">${t.meanY.toFixed(2)}</td>
                    <td style="padding: 5px; border: 1px solid #ddd; text-align: center;">${t.sdY.toFixed(2)}</td>
                </tr>
            `;
        });

        tableHtml += '</tbody></table></div>';

        // Hide scatterPlot, use compareTable for side-by-side layout
        const scatterPlotEl = document.getElementById('scatterPlot');
        const compareTableEl = document.getElementById('compareTable');

        scatterPlotEl.style.display = 'none';

        // Create side-by-side layout in compareTable
        const chartHeight = Math.max(400, tissueStats.length * 22 + 80);
        compareTableEl.style.display = 'block';
        compareTableEl.style.maxHeight = 'none';
        compareTableEl.innerHTML = `
            <div style="display: flex; gap: 20px; align-items: flex-start;">
                <div style="flex: 0 0 45%;">
                    <div id="byTissueChart" style="height: ${chartHeight}px;"></div>
                </div>
                <div style="flex: 1; min-width: 0;">
                    ${tableHtml}
                </div>
            </div>
        `;

        // Render bar chart in the new container
        Plotly.newPlot('byTissueChart', [trace], { ...layout, height: chartHeight }, { responsive: true });
    }

    downloadScatterCSV() {
        if (!this.currentInspect) return;

        const hotspotGene = document.getElementById('hotspotGene').value;
        let header = 'CellLine,CellLineID,Lineage,Gene1_Effect,Gene2_Effect';
        if (hotspotGene && this.mutations?.geneData?.[hotspotGene]) {
            header += `,${hotspotGene}_mutation`;
        }
        header += '\n';

        let csv = header;
        const mutationData = hotspotGene && this.mutations?.geneData?.[hotspotGene]?.mutations;

        this.currentInspect.data.forEach(d => {
            let row = `"${d.cellLineName}","${d.cellLineId}","${d.lineage}",${d.x},${d.y}`;
            if (mutationData) {
                const mutLevel = mutationData[d.cellLineId] || 0;
                row += `,${mutLevel}`;
            }
            csv += row + '\n';
        });

        const suffix = hotspotGene ? `_${hotspotGene}` : '';
        this.downloadFile(csv,
            `scatter_${this.currentInspect.gene1}_vs_${this.currentInspect.gene2}${suffix}.csv`,
            'text/csv');
    }

    downloadTissueChartPNG() {
        if (!this.currentInspect) return;
        const chartEl = document.getElementById('byTissueChart');
        if (!chartEl) return;

        Plotly.downloadImage(chartEl, {
            format: 'png',
            width: 800,
            height: Math.max(400, (this.currentTissueStats?.length || 10) * 25 + 100),
            filename: `by_tissue_${this.currentInspect.gene1}_vs_${this.currentInspect.gene2}`
        });
    }

    downloadTissueChartSVG() {
        if (!this.currentInspect) return;
        const chartEl = document.getElementById('byTissueChart');
        if (!chartEl) return;

        Plotly.downloadImage(chartEl, {
            format: 'svg',
            width: 800,
            height: Math.max(400, (this.currentTissueStats?.length || 10) * 25 + 100),
            filename: `by_tissue_${this.currentInspect.gene1}_vs_${this.currentInspect.gene2}`
        });
    }

    downloadTissueTableCSV() {
        if (!this.currentInspect || !this.currentTissueStats) return;

        const { gene1, gene2 } = this.currentInspect;

        let csv = `Lineage,N,Correlation,${gene1}_Effect_mean,${gene1}_Effect_SD,${gene2}_Effect_mean,${gene2}_Effect_SD\n`;

        this.currentTissueStats.forEach(t => {
            csv += `"${t.tissue}",${t.n},${t.correlation.toFixed(4)},${t.meanX.toFixed(4)},${t.sdX.toFixed(4)},${t.meanY.toFixed(4)},${t.sdY.toFixed(4)}\n`;
        });

        this.downloadFile(csv,
            `by_tissue_${gene1}_vs_${gene2}.csv`,
            'text/csv');
    }
}

// Initialize app
const app = new CorrelationExplorer();
