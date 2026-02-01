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

        // Extended synonym lookup (low/mid risk)
        this.synonymLookup = null;

        // Network physics and layout state
        this.physicsEnabled = true;
        this.currentLayout = 0;

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
        const [metadataRes, cellLineRes, mutationsRes, orthologsRes, synonymsRes] = await Promise.all([
            fetch('web_data/metadata.json'),
            fetch('web_data/cellLineMetadata.json'),
            fetch('web_data/mutations.json'),
            fetch('web_data/orthologs.json'),
            fetch('web_data/synonyms.json')
        ]);

        this.metadata = await metadataRes.json();
        this.cellLineMetadata = await cellLineRes.json();
        this.mutations = await mutationsRes.json();
        this.orthologs = await orthologsRes.json();
        this.synonymLookup = await synonymsRes.json();

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
        // Count cell lines per lineage
        const lineageCounts = {};
        const subLineageCounts = {};

        if (this.cellLineMetadata && this.cellLineMetadata.lineage) {
            const cellLines = this.metadata.cellLines;
            cellLines.forEach(cellLine => {
                const lineage = this.cellLineMetadata.lineage[cellLine];
                const subLineage = this.cellLineMetadata.lineageSubtype?.[cellLine] || '';

                if (lineage) {
                    lineageCounts[lineage] = (lineageCounts[lineage] || 0) + 1;

                    if (subLineage) {
                        const key = `${lineage}|${subLineage}`;
                        subLineageCounts[key] = (subLineageCounts[key] || 0) + 1;
                    }
                }
            });
        }

        this.lineageCounts = lineageCounts;
        this.subLineageCounts = subLineageCounts;

        if (Object.keys(lineageCounts).length > 0) {
            const select = document.getElementById('lineageFilter');
            const total = this.metadata.cellLines.length;
            select.innerHTML = `<option value="">All lineages (n=${total})</option>`;

            Object.keys(lineageCounts).sort().forEach(lineage => {
                const option = document.createElement('option');
                option.value = lineage;
                option.textContent = `${lineage} (n=${lineageCounts[lineage]})`;
                select.appendChild(option);
            });
            document.getElementById('lineageFilterGroup').style.display = 'block';

            // Update sub-lineage when lineage changes
            select.addEventListener('change', () => this.updateSubLineageFilter());
        }

        // Also populate parameter hotspot filter
        this.populateParamHotspotFilter();
    }

    updateSubLineageFilter() {
        const lineage = document.getElementById('lineageFilter').value;
        const subSelect = document.getElementById('subLineageFilter');
        const isMutationMode = document.querySelector('input[name="analysisMode"]:checked')?.value === 'mutation';

        if (!lineage) {
            document.getElementById('subLineageFilterGroup').style.display = 'none';
            subSelect.innerHTML = '<option value="">All subtypes</option>';
            // Update hotspot counts for all lineages
            this.updateHotspotCountsForCurrentFilters();
            return;
        }

        // Find sub-lineages for this lineage
        const subLineages = {};
        Object.keys(this.subLineageCounts).forEach(key => {
            if (key.startsWith(lineage + '|')) {
                const subLineage = key.split('|')[1];
                subLineages[subLineage] = this.subLineageCounts[key];
            }
        });

        if (Object.keys(subLineages).length > 0) {
            const lineageCount = this.lineageCounts[lineage];
            subSelect.innerHTML = `<option value="">All subtypes (n=${lineageCount})</option>`;

            Object.keys(subLineages).sort().forEach(sub => {
                const option = document.createElement('option');
                option.value = sub;
                option.textContent = `${sub} (n=${subLineages[sub]})`;
                subSelect.appendChild(option);
            });
            document.getElementById('subLineageFilterGroup').style.display = 'block';

            // Add listener for sub-lineage changes (only add once)
            if (!subSelect.hasAttribute('data-listener-attached')) {
                subSelect.addEventListener('change', () => {
                    this.updateHotspotCountsForCurrentFilters();
                });
                subSelect.setAttribute('data-listener-attached', 'true');
            }
        } else {
            document.getElementById('subLineageFilterGroup').style.display = 'none';
        }

        // Update hotspot counts for selected lineage
        this.updateHotspotCountsForCurrentFilters();
    }

    updateHotspotCountsForCurrentFilters() {
        if (this.mutations?.geneData) {
            this.updateParamHotspotGeneCounts();

            // Also update mutation mode selector if in mutation mode
            const isMutationMode = document.querySelector('input[name="analysisMode"]:checked')?.value === 'mutation';
            if (isMutationMode) {
                this.populateMutationHotspotSelector();
            }
        }
    }

    populateParamHotspotFilter() {
        if (this.mutations && this.mutations.geneData) {
            const select = document.getElementById('paramHotspotGene');
            document.getElementById('paramHotspotFilterGroup').style.display = 'block';

            // Update level dropdown with counts when gene changes
            select.addEventListener('change', () => this.updateParamHotspotLevelCounts());

            // Initial population
            this.updateParamHotspotGeneCounts();
        }
    }

    updateParamHotspotGeneCounts() {
        const genes = Object.keys(this.mutations.geneData).sort();
        const select = document.getElementById('paramHotspotGene');
        const cellLines = this.metadata.cellLines;
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;
        const currentValue = select.value;

        select.innerHTML = '<option value="">No filter</option>';
        genes.forEach(gene => {
            // Count mutations for this gene (respecting lineage filter)
            const mutations = this.mutations.geneData[gene].mutations;
            let nMut = 0;
            cellLines.forEach(cl => {
                // Apply lineage filter
                if (lineageFilter && this.cellLineMetadata?.lineage?.[cl] !== lineageFilter) {
                    return;
                }
                // Apply sub-lineage filter
                if (subLineageFilter && this.cellLineMetadata?.lineageSubtype?.[cl] !== subLineageFilter) {
                    return;
                }
                if (mutations[cl] && mutations[cl] > 0) nMut++;
            });

            const option = document.createElement('option');
            option.value = gene;
            option.textContent = `${gene} (n=${nMut} mutated)`;
            select.appendChild(option);
        });

        // Restore selection if it was set
        if (currentValue) {
            select.value = currentValue;
        }

        // Also update level counts
        this.updateParamHotspotLevelCounts();
    }

    updateParamHotspotLevelCounts() {
        const gene = document.getElementById('paramHotspotGene').value;
        const levelSelect = document.getElementById('paramHotspotLevel');
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;

        if (!gene || !this.mutations?.geneData?.[gene]) {
            levelSelect.innerHTML = `
                <option value="all">All cells</option>
                <option value="0">Only WT (0 mutations)</option>
                <option value="1">Only 1 mutation</option>
                <option value="2">Only 2 mutations</option>
                <option value="1+2">Only mutated (1+2)</option>
            `;
            return;
        }

        // Count mutations for selected gene (respecting lineage filter)
        const mutations = this.mutations.geneData[gene].mutations;
        const cellLines = this.metadata.cellLines;
        let n0 = 0, n1 = 0, n2 = 0;

        cellLines.forEach(cellLine => {
            // Apply lineage filter
            if (lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== lineageFilter) {
                return;
            }
            // Apply sub-lineage filter
            if (subLineageFilter && this.cellLineMetadata?.lineageSubtype?.[cellLine] !== subLineageFilter) {
                return;
            }

            const level = mutations[cellLine] || 0;
            if (level === 0) n0++;
            else if (level === 1) n1++;
            else n2++;
        });

        const nMut = n1 + n2;
        const total = n0 + n1 + n2;

        levelSelect.innerHTML = `
            <option value="all">All cells (n=${total})</option>
            <option value="0">Only WT (n=${n0})</option>
            <option value="1">Only 1 mutation (n=${n1})</option>
            <option value="2">Only 2 mutations (n=${n2})</option>
            <option value="1+2">Only mutated 1+2 (n=${nMut})</option>
        `;
    }

    updateAnalysisModeUI() {
        const mode = document.querySelector('input[name="analysisMode"]:checked').value;
        const isMutationMode = mode === 'mutation';

        // Toggle visibility of correlation/slope params
        document.getElementById('correlationParams').style.display = isMutationMode ? 'none' : 'block';
        document.getElementById('slopeParams').style.display = isMutationMode ? 'none' : 'block';

        // Toggle visibility of mutation-specific params
        document.getElementById('mutationHotspotGroup').style.display = isMutationMode ? 'block' : 'none';
        document.getElementById('pValueThresholdGroup').style.display = isMutationMode ? 'block' : 'none';

        // Show/hide mutation tab
        document.getElementById('mutationTab').style.display = isMutationMode ? 'inline-block' : 'none';

        // Set default min cell lines based on mode
        const minCellLinesInput = document.getElementById('minCellLines');
        if (isMutationMode) {
            minCellLinesInput.value = '20';
            this.populateMutationHotspotSelector();
        } else {
            minCellLinesInput.value = '50';
        }
    }

    populateMutationHotspotSelector() {
        const select = document.getElementById('mutationHotspotSelect');
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;
        const currentValue = select.value;

        if (!this.mutations || !this.mutations.geneData) return;

        const genes = Object.keys(this.mutations.geneData).sort();
        const cellLines = this.metadata.cellLines;

        select.innerHTML = '<option value="">Select hotspot gene...</option>';
        genes.forEach(gene => {
            const mutations = this.mutations.geneData[gene].mutations;
            let nMut = 0;

            cellLines.forEach(cl => {
                // Apply lineage filter
                if (lineageFilter && this.cellLineMetadata?.lineage?.[cl] !== lineageFilter) {
                    return;
                }
                // Apply sub-lineage filter
                if (subLineageFilter && this.cellLineMetadata?.lineageSubtype?.[cl] !== subLineageFilter) {
                    return;
                }
                if (mutations[cl] && mutations[cl] > 0) nMut++;
            });

            const option = document.createElement('option');
            option.value = gene;
            option.textContent = `${gene} (${nMut} mutated cells)`;
            select.appendChild(option);
        });

        // Restore selection if it was set
        if (currentValue) {
            select.value = currentValue;
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
                'ATM', 'ERBB2', 'CDK4', 'MDM2', 'NRAS', 'ARID1A', 'TSC1', 'TSC2'];
            document.getElementById('geneTextarea').value = testGenes.join('\n');
            this.updateGeneCount();
        });

        // Find synonyms button
        document.getElementById('findSynonyms').addEventListener('click', async () => await this.findSynonymsForMissingGenes());

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

        // Analysis mode change
        document.querySelectorAll('input[name="analysisMode"]').forEach(radio => {
            radio.addEventListener('change', () => this.updateAnalysisModeUI());
        });

        // Mutation results search
        document.getElementById('mutationSearch').addEventListener('input', (e) => {
            this.filterMutationTable(e.target.value);
        });

        // Download mutation results
        document.getElementById('downloadMutationResults').addEventListener('click', () => {
            this.downloadMutationResults();
        });

        // Gene effect distribution modal
        document.getElementById('closeGeneEffect').addEventListener('click', () => {
            document.getElementById('geneEffectModal').style.display = 'none';
        });
        document.getElementById('geneEffectModal').addEventListener('click', (e) => {
            if (e.target.id === 'geneEffectModal') {
                document.getElementById('geneEffectModal').style.display = 'none';
            }
        });
        document.getElementById('downloadGeneEffectPNG').addEventListener('click', () => this.downloadGeneEffectPNG());
        document.getElementById('downloadGeneEffectSVG').addEventListener('click', () => this.downloadGeneEffectSVG());
        document.getElementById('downloadGeneEffectCSV').addEventListener('click', () => this.downloadGeneEffectCSV());

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
        document.getElementById('showGeneEffect').addEventListener('change', (e) => {
            document.getElementById('showGESDGroup').style.display = e.target.checked ? 'inline' : 'none';
            this.updateNetworkLabels();
        });
        document.getElementById('showGeneEffectSD').addEventListener('change', () => this.updateNetworkLabels());

        // Color by gene effect controls (mutually exclusive with stats)
        document.getElementById('colorByGeneEffect').addEventListener('change', (e) => {
            if (e.target.checked) {
                // Uncheck color by stats
                document.getElementById('colorByStats').checked = false;
                document.getElementById('colorStatsOptions').style.display = 'none';
            }
            document.getElementById('colorGEOptions').style.display = e.target.checked ? 'block' : 'none';
            this.updateNetworkColors();
        });
        document.querySelectorAll('input[name="colorGEType"]').forEach(radio => {
            radio.addEventListener('change', () => this.updateNetworkColors());
        });

        document.getElementById('downloadNetworkPNG').addEventListener('click', () => this.downloadNetworkPNG());
        document.getElementById('downloadNetworkSVG').addEventListener('click', () => this.downloadNetworkSVG());
        document.getElementById('downloadAllData').addEventListener('click', () => this.downloadAllData());

        // Color by stats controls (mutually exclusive with GE)
        document.getElementById('colorByStats').addEventListener('change', (e) => {
            if (e.target.checked) {
                // Uncheck color by gene effect
                document.getElementById('colorByGeneEffect').checked = false;
                document.getElementById('colorGEOptions').style.display = 'none';
            }
            document.getElementById('colorStatsOptions').style.display = e.target.checked ? 'block' : 'none';
            document.getElementById('legendNodeColor').style.display = e.target.checked ? 'block' : 'none';
            this.updateNetworkColors();
        });

        // Physics toggle and layout change buttons
        document.getElementById('togglePhysics').addEventListener('click', () => this.togglePhysics());
        document.getElementById('changeLayout').addEventListener('click', () => this.changeNetworkLayout());
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
        document.getElementById('compareAllMutationsBtn')?.addEventListener('click', () => this.showCompareAllMutations());

        // Aspect ratio control
        document.getElementById('aspectRatio')?.addEventListener('input', (e) => {
            document.getElementById('aspectRatioValue').textContent = parseFloat(e.target.value).toFixed(1);
            this.updateInspectPlot();
        });

        // Table header sorting
        document.querySelectorAll('.data-table th[data-sort]').forEach(th => {
            th.addEventListener('click', () => this.sortTable(th));
        });

        // Infographic modal
        document.getElementById('showInfoGraphic')?.addEventListener('click', () => {
            document.getElementById('infographicModal').style.display = 'flex';
        });
        document.getElementById('closeInfoGraphic')?.addEventListener('click', () => {
            document.getElementById('infographicModal').style.display = 'none';
        });
        document.getElementById('infographicModal')?.addEventListener('click', (e) => {
            if (e.target.id === 'infographicModal') {
                document.getElementById('infographicModal').style.display = 'none';
            }
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

    async findSynonymsForMissingGenes() {
        if (!this.genesNotFound || this.genesNotFound.length === 0) return;

        const replacements = [];
        const stillNotFound = [];
        const notFound = this.genesNotFound;

        // First pass: check local lookups (synonym table + ortholog table)
        notFound.forEach(gene => {
            const upperGene = gene.toUpperCase();

            // First check extended synonym lookup (low/mid risk from DepMap reference)
            if (this.synonymLookup) {
                const match = this.synonymLookup[upperGene];
                if (match && this.geneIndex.has(match.d.toUpperCase())) {
                    const sourceLabel = match.r === 'l' ? 'low-risk' : 'mid-risk';
                    replacements.push({
                        original: gene,
                        replacement: match.d.toUpperCase(),
                        source: sourceLabel
                    });
                    return;
                }
            }

            // Fallback: Check ortholog lookup (mouse to human)
            const humanGene = this.orthologs?.mouseToHuman?.[gene];
            if (humanGene && this.geneIndex.has(humanGene.toUpperCase())) {
                replacements.push({
                    original: gene,
                    replacement: humanGene.toUpperCase(),
                    source: 'ortholog'
                });
                return;
            }

            // Still not found - will try API
            stillNotFound.push(gene);
        });

        // Second pass: try MyGene.info API for remaining genes
        if (stillNotFound.length > 0) {
            const btn = document.getElementById('findSynonyms');
            const originalText = btn.textContent;
            btn.textContent = 'Searching API...';
            btn.disabled = true;

            try {
                const apiResults = await this.queryMyGeneAPI(stillNotFound);
                apiResults.forEach(r => {
                    if (r.replacement && this.geneIndex.has(r.replacement.toUpperCase())) {
                        replacements.push({
                            original: r.original,
                            replacement: r.replacement.toUpperCase(),
                            source: 'MyGene.info API'
                        });
                    }
                });
            } catch (error) {
                console.warn('MyGene.info API query failed:', error);
            } finally {
                btn.textContent = originalText;
                btn.disabled = false;
            }
        }

        if (replacements.length > 0) {
            // Store synonyms used for summary
            this.synonymsUsed = replacements;

            // Remove replaced genes from genesNotFound
            const replacedOriginals = new Set(replacements.map(r => r.original.toUpperCase()));
            this.genesNotFound = this.genesNotFound.filter(g => !replacedOriginals.has(g.toUpperCase()));

            // Update textarea
            const textarea = document.getElementById('geneTextarea');
            let text = textarea.value;

            replacements.forEach(r => {
                const regex = new RegExp(`\\b${r.original}\\b`, 'gi');
                text = text.replace(regex, r.replacement);
            });

            textarea.value = text;
            this.updateGeneCount();

            const msg = replacements.map(r => `${r.original} → ${r.replacement} [${r.source}]`).join('\n');
            alert(`Replaced ${replacements.length} gene(s):\n${msg}`);
        } else {
            alert('No synonyms or orthologs found for the missing genes');
        }
    }

    async queryMyGeneAPI(genes) {
        // Query MyGene.info API for gene synonyms
        // API docs: https://docs.mygene.info/en/latest/
        const results = [];

        // Query in batches of 10 to avoid overloading
        const batchSize = 10;
        for (let i = 0; i < genes.length; i += batchSize) {
            const batch = genes.slice(i, i + batchSize);

            await Promise.all(batch.map(async (gene) => {
                try {
                    const url = `https://mygene.info/v3/query?q=${encodeURIComponent(gene)}&scopes=symbol,alias&fields=symbol&species=human`;
                    const response = await fetch(url);

                    if (response.ok) {
                        const data = await response.json();
                        if (data.hits && data.hits.length > 0) {
                            // Take the first hit's symbol
                            const symbol = data.hits[0].symbol;
                            if (symbol && symbol.toUpperCase() !== gene.toUpperCase()) {
                                results.push({ original: gene, replacement: symbol });
                            }
                        }
                    }
                } catch (error) {
                    console.warn(`MyGene.info query failed for ${gene}:`, error);
                }
            }));

            // Small delay between batches to be nice to the API
            if (i + batchSize < genes.length) {
                await new Promise(resolve => setTimeout(resolve, 100));
            }
        }

        return results;
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
            { gene: 'CREBBP', lfc: 0.4, fdr: 0.18 },
            { gene: 'TSC1', lfc: 1.0, fdr: 0.02 },
            { gene: 'TSC2', lfc: 0.9, fdr: 0.03 }
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
        const hotspotGene = document.getElementById('paramHotspotGene').value;
        const hotspotLevel = document.getElementById('paramHotspotLevel').value;

        // Get mutation data for hotspot filter
        let mutationData = null;
        if (hotspotGene && this.mutations?.geneData?.[hotspotGene]) {
            mutationData = this.mutations.geneData[hotspotGene].mutations;
        }

        const indices = [];
        this.metadata.cellLines.forEach((cellLine, idx) => {
            // Check lineage filter
            if (lineageFilter) {
                if (!this.cellLineMetadata.lineage ||
                    this.cellLineMetadata.lineage[cellLine] !== lineageFilter) {
                    return;
                }
            }

            // Check hotspot mutation filter
            if (mutationData && hotspotLevel !== 'all') {
                const mutLevel = mutationData[cellLine] || 0;
                if (hotspotLevel === '0' && mutLevel !== 0) return;
                if (hotspotLevel === '1' && mutLevel !== 1) return;
                if (hotspotLevel === '2' && mutLevel < 2) return;
                if (hotspotLevel === '1+2' && mutLevel < 1) return;
            }

            indices.push(idx);
        });

        // Return all indices if no filters applied
        if (indices.length === 0 && !lineageFilter && (!hotspotGene || hotspotLevel === 'all')) {
            return Array.from({ length: this.nCellLines }, (_, i) => i);
        }

        return indices;
    }

    runAnalysis() {
        // Reset network settings to defaults when running new analysis
        this.resetNetworkSettings();

        const mode = document.querySelector('input[name="analysisMode"]:checked').value;

        // Handle mutation analysis mode separately
        if (mode === 'mutation') {
            this.runMutationAnalysis();
            return;
        }

        const geneList = this.getGeneList();
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

    runMutationAnalysis() {
        const hotspotGene = document.getElementById('mutationHotspotSelect').value;
        const minN = parseInt(document.getElementById('minCellLines').value);
        const pThreshold = parseFloat(document.getElementById('pValueThreshold').value);
        const lineageFilter = document.getElementById('lineageFilter').value;

        // Get additional hotspot filter (from parameter section)
        const additionalHotspot = document.getElementById('paramHotspotGene').value;
        const additionalHotspotLevel = document.getElementById('paramHotspotLevel').value;

        if (!hotspotGene) {
            this.showStatus('error', 'Please select a hotspot mutation');
            return;
        }

        this.showStatus('info', 'Running mutation analysis...');

        // Use setTimeout to allow UI to update
        setTimeout(() => {
            try {
                const analysisResult = this.calculateMutationAnalysis(hotspotGene, minN, lineageFilter, additionalHotspot, additionalHotspotLevel);

                // Filter by p-value threshold
                const significantResults = analysisResult.results.filter(r => r.p_mut < pThreshold || r.p_2 < pThreshold);

                // Sort by p-value (1+2 vs 0)
                significantResults.sort((a, b) => a.p_mut - b.p_mut);

                this.mutationResults = {
                    hotspotGene,
                    pThreshold,
                    minN,
                    lineageFilter,
                    additionalHotspot,
                    additionalHotspotLevel,
                    nWT: analysisResult.nWT,
                    nMut: analysisResult.nMut,
                    n2: analysisResult.n2,
                    allResults: analysisResult.results,
                    significantResults
                };

                this.displayMutationResults();

                // Switch to mutation tab
                document.querySelectorAll('.nav-link').forEach(link => link.classList.remove('active'));
                document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
                document.querySelector('[data-tab="mutation"]').classList.add('active');
                document.getElementById('tab-mutation').classList.add('active');

                this.showStatus('success',
                    `&#10003; Mutation analysis complete: ${significantResults.length} genes with p < ${pThreshold}`);
            } catch (error) {
                console.error('Mutation analysis error:', error);
                this.showStatus('error', 'Mutation analysis failed: ' + error.message);
            }
        }, 50);
    }

    calculateMutationAnalysis(hotspotGene, minN, lineageFilter, additionalHotspot, additionalHotspotLevel) {
        const mutationData = this.mutations.geneData[hotspotGene];
        if (!mutationData) {
            throw new Error(`No mutation data for ${hotspotGene}`);
        }

        // Get additional hotspot mutation data if specified
        const additionalMutData = additionalHotspot ? this.mutations.geneData[additionalHotspot] : null;

        const cellLines = this.metadata.cellLines;
        const results = [];

        // Categorize cell lines by mutation status and lineage filter
        const wtCellIndices = [];
        const mut1CellIndices = [];
        const mut2CellIndices = [];

        cellLines.forEach((cellLine, idx) => {
            // Check lineage filter
            if (lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== lineageFilter) {
                return;
            }

            // Check additional hotspot filter
            if (additionalMutData && additionalHotspotLevel !== 'all') {
                const addMutLevel = additionalMutData.mutations[cellLine] || 0;
                if (additionalHotspotLevel === '0' && addMutLevel !== 0) return;
                if (additionalHotspotLevel === '1' && addMutLevel !== 1) return;
                if (additionalHotspotLevel === '2' && addMutLevel < 2) return;
                if (additionalHotspotLevel === '1+2' && addMutLevel === 0) return;
            }

            const mutLevel = mutationData.mutations[cellLine] || 0;
            if (mutLevel === 0) {
                wtCellIndices.push(idx);
            } else if (mutLevel === 1) {
                mut1CellIndices.push(idx);
            } else {
                mut2CellIndices.push(idx);
            }
        });

        const mutAllCellIndices = [...mut1CellIndices, ...mut2CellIndices];

        // Check minimum cell count
        if (wtCellIndices.length < 3 || mutAllCellIndices.length < 3) {
            throw new Error(`Not enough cell lines: WT=${wtCellIndices.length}, Mutated=${mutAllCellIndices.length}`);
        }

        // Analyze each gene
        for (let geneIdx = 0; geneIdx < this.nGenes; geneIdx++) {
            const gene = this.geneNames[geneIdx];

            // Get gene effect values for each group
            const wtEffects = this.getGeneEffectsForCells(geneIdx, wtCellIndices);
            const mutAllEffects = this.getGeneEffectsForCells(geneIdx, mutAllCellIndices);
            const mut2Effects = this.getGeneEffectsForCells(geneIdx, mut2CellIndices);

            // Skip if not enough valid values
            if (wtEffects.length < minN || mutAllEffects.length < 3) continue;

            // Calculate statistics for WT vs 1+2
            const wtMean = this.mean(wtEffects);
            const mutMean = this.mean(mutAllEffects);
            const diff_mut = mutMean - wtMean;
            const tTest_mut = this.welchTTest(wtEffects, mutAllEffects);

            // Calculate statistics for WT vs 2 (if enough cells)
            let n_2 = mut2Effects.length;
            let mean_2 = NaN;
            let diff_2 = NaN;
            let p_2 = 1;

            if (mut2Effects.length >= 3) {
                mean_2 = this.mean(mut2Effects);
                diff_2 = mean_2 - wtMean;
                const tTest_2 = this.welchTTest(wtEffects, mut2Effects);
                p_2 = tTest_2.p;
            }

            results.push({
                gene,
                n_wt: wtEffects.length,
                mean_wt: wtMean,
                n_mut: mutAllEffects.length,
                mean_mut: mutMean,
                diff_mut,
                p_mut: tTest_mut.p,
                n_2,
                mean_2,
                diff_2,
                p_2
            });
        }

        return {
            results,
            nWT: wtCellIndices.length,
            nMut: mutAllCellIndices.length,
            n2: mut2CellIndices.length
        };
    }

    getGeneEffectsForCells(geneIdx, cellIndices) {
        const effects = [];
        for (const cellIdx of cellIndices) {
            const value = this.geneEffects[geneIdx * this.nCellLines + cellIdx];
            if (!isNaN(value)) {
                effects.push(value);
            }
        }
        return effects;
    }

    mean(arr) {
        if (arr.length === 0) return NaN;
        return arr.reduce((a, b) => a + b, 0) / arr.length;
    }

    variance(arr) {
        if (arr.length < 2) return 0;
        const m = this.mean(arr);
        return arr.reduce((acc, val) => acc + (val - m) ** 2, 0) / (arr.length - 1);
    }

    welchTTest(group1, group2) {
        // Welch's t-test for unequal variances
        const n1 = group1.length;
        const n2 = group2.length;
        const m1 = this.mean(group1);
        const m2 = this.mean(group2);
        const v1 = this.variance(group1);
        const v2 = this.variance(group2);

        if (n1 < 2 || n2 < 2) {
            return { t: NaN, df: NaN, p: 1 };
        }

        const se = Math.sqrt(v1 / n1 + v2 / n2);
        if (se === 0) {
            return { t: 0, df: n1 + n2 - 2, p: 1 };
        }

        const t = (m1 - m2) / se;

        // Welch-Satterthwaite degrees of freedom
        const df = Math.pow(v1 / n1 + v2 / n2, 2) /
            (Math.pow(v1 / n1, 2) / (n1 - 1) + Math.pow(v2 / n2, 2) / (n2 - 1));

        // Two-tailed p-value using t-distribution approximation
        const p = this.tDistributionPValue(Math.abs(t), df);

        return { t, df, p };
    }

    tDistributionPValue(t, df) {
        // Approximation of two-tailed p-value for t-distribution
        // Using normal approximation for large df, or beta approximation for small df
        if (df <= 0 || isNaN(t) || isNaN(df)) return 1;

        // For large df, approximate with normal distribution
        if (df > 100) {
            return 2 * (1 - this.normalCDF(t));
        }

        // Beta function approximation for t-distribution CDF
        const x = df / (df + t * t);
        const a = df / 2;
        const b = 0.5;

        // Incomplete beta function approximation
        const betaInc = this.incompleteBeta(x, a, b);
        return betaInc;
    }

    incompleteBeta(x, a, b) {
        // Simplified incomplete beta function for t-distribution p-value
        // This is an approximation suitable for statistical testing
        if (x <= 0) return 0;
        if (x >= 1) return 1;

        // Use continued fraction expansion (simplified)
        const bt = Math.exp(
            this.logGamma(a + b) - this.logGamma(a) - this.logGamma(b) +
            a * Math.log(x) + b * Math.log(1 - x)
        );

        if (x < (a + 1) / (a + b + 2)) {
            return bt * this.betaCF(x, a, b) / a;
        } else {
            return 1 - bt * this.betaCF(1 - x, b, a) / b;
        }
    }

    betaCF(x, a, b) {
        // Continued fraction for incomplete beta
        const maxIter = 100;
        const eps = 1e-10;

        let qab = a + b;
        let qap = a + 1;
        let qam = a - 1;
        let c = 1;
        let d = 1 - qab * x / qap;
        if (Math.abs(d) < 1e-30) d = 1e-30;
        d = 1 / d;
        let h = d;

        for (let m = 1; m <= maxIter; m++) {
            let m2 = 2 * m;
            let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1 + aa * d;
            if (Math.abs(d) < 1e-30) d = 1e-30;
            c = 1 + aa / c;
            if (Math.abs(c) < 1e-30) c = 1e-30;
            d = 1 / d;
            h *= d * c;

            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1 + aa * d;
            if (Math.abs(d) < 1e-30) d = 1e-30;
            c = 1 + aa / c;
            if (Math.abs(c) < 1e-30) c = 1e-30;
            d = 1 / d;
            let del = d * c;
            h *= del;

            if (Math.abs(del - 1) < eps) break;
        }

        return h;
    }

    logGamma(x) {
        // Lanczos approximation for log gamma function
        const g = 7;
        const c = [
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        ];

        if (x < 0.5) {
            return Math.log(Math.PI / Math.sin(Math.PI * x)) - this.logGamma(1 - x);
        }

        x -= 1;
        let a = c[0];
        for (let i = 1; i < g + 2; i++) {
            a += c[i] / (x + i);
        }

        const t = x + g + 0.5;
        return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
    }

    formatPValue(p) {
        // Format p-value with 1 decimal in exponent
        if (p >= 1 || isNaN(p)) return '-';
        if (p === 0 || p < 1e-300) return '1.0e-300';
        if (p < 0.001) {
            // Format as exponential with 1 decimal (e.g., 2.2e-10)
            const exp = Math.floor(Math.log10(p));
            const mantissa = p / Math.pow(10, exp);
            return `${mantissa.toFixed(1)}e${exp}`;
        }
        return p.toFixed(4);
    }

    displayMutationResults() {
        if (!this.mutationResults) return;

        const mr = this.mutationResults;
        const results = mr.significantResults;
        const tbody = document.getElementById('mutationTableBody');
        tbody.innerHTML = '';

        results.forEach(r => {
            const row = document.createElement('tr');
            row.innerHTML = `
                <td><a href="#" class="inspect-link" onclick="app.showGeneEffectDistribution('${r.gene}'); return false;">Inspect</a></td>
                <td>${r.gene}</td>
                <td>${r.n_wt}</td>
                <td>${r.mean_wt.toFixed(3)}</td>
                <td>${r.n_mut}</td>
                <td>${r.mean_mut.toFixed(3)}</td>
                <td class="${r.diff_mut < 0 ? 'negative' : 'positive'}">${r.diff_mut.toFixed(3)}</td>
                <td>${this.formatPValue(r.p_mut)}</td>
                <td>${r.n_2}</td>
                <td>${isNaN(r.mean_2) ? '-' : r.mean_2.toFixed(3)}</td>
                <td class="${r.diff_2 < 0 ? 'negative' : 'positive'}">${isNaN(r.diff_2) ? '-' : r.diff_2.toFixed(3)}</td>
                <td>${this.formatPValue(r.p_2)}</td>
            `;
            tbody.appendChild(row);
        });

        // Build settings summary
        let settingsText = `Hotspot: ${mr.hotspotGene} | `;
        settingsText += `WT: ${mr.nWT} cells | Mutated: ${mr.nMut} cells | `;
        settingsText += `Min cells: ${mr.minN} | p < ${mr.pThreshold}`;
        if (mr.lineageFilter) {
            settingsText += ` | Lineage: ${mr.lineageFilter}`;
        }
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            settingsText += ` | Filter: ${mr.additionalHotspot} ${mr.additionalHotspotLevel}`;
        }

        document.getElementById('mutationResultsCount').innerHTML =
            `<strong>${results.length} genes</strong> with p < ${mr.pThreshold}<br>
            <small style="color: #666;">${settingsText}</small>`;

        // Store for sorting
        this.mutationTableData = results;
    }

    filterMutationTable(query) {
        const tbody = document.getElementById('mutationTableBody');
        const rows = tbody.querySelectorAll('tr');
        const lowerQuery = query.toLowerCase();

        rows.forEach(row => {
            const gene = row.cells[0].textContent.toLowerCase();
            row.style.display = gene.includes(lowerQuery) ? '' : 'none';
        });
    }

    sortMutationTable(th) {
        if (!this.mutationTableData) return;

        const col = th.dataset.col;
        const currentDir = th.dataset.sortDir === 'asc' ? 'desc' : 'asc';
        th.dataset.sortDir = currentDir;

        // Update header indicators
        th.closest('tr').querySelectorAll('th').forEach(h => {
            h.textContent = h.textContent.replace(' ▲', '').replace(' ▼', '');
        });
        th.textContent += currentDir === 'asc' ? ' ▲' : ' ▼';

        // Sort data
        this.mutationTableData.sort((a, b) => {
            let valA = a[col];
            let valB = b[col];

            if (typeof valA === 'string') {
                return currentDir === 'asc' ?
                    valA.localeCompare(valB) : valB.localeCompare(valA);
            } else {
                if (isNaN(valA)) valA = currentDir === 'asc' ? Infinity : -Infinity;
                if (isNaN(valB)) valB = currentDir === 'asc' ? Infinity : -Infinity;
                return currentDir === 'asc' ? valA - valB : valB - valA;
            }
        });

        // Update significant results reference for export
        this.mutationResults.significantResults = this.mutationTableData;

        // Re-render table
        this.displayMutationResults();
    }

    downloadMutationResults() {
        if (!this.mutationResults) return;

        const mr = this.mutationResults;
        const results = mr.significantResults;

        // Build settings header
        let csv = '# Mutation Analysis Results\n';
        csv += `# Hotspot Mutation: ${mr.hotspotGene}\n`;
        csv += `# WT cells (0 mutations): ${mr.nWT}\n`;
        csv += `# Mutated cells (1+2 mutations): ${mr.nMut}\n`;
        csv += `# Cells with 2 mutations: ${mr.n2}\n`;
        csv += `# Min cell lines: ${mr.minN}\n`;
        csv += `# P-value threshold: ${mr.pThreshold}\n`;
        csv += `# Lineage filter: ${mr.lineageFilter || 'All lineages'}\n`;
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            csv += `# Additional hotspot filter: ${mr.additionalHotspot} = ${mr.additionalHotspotLevel}\n`;
        }
        csv += `# Date: ${new Date().toISOString().slice(0, 10)}\n`;
        csv += '#\n';

        const headers = ['Gene', 'N_WT', 'Mean_GE_WT', 'N_1+2', 'Mean_GE_1+2', 'Delta_GE', 'pValue_1+2_vs_0',
                        'N_2', 'Mean_GE_2', 'Delta_GE_2vs0', 'pValue_2_vs_0'];

        csv += headers.join(',') + '\n';
        results.forEach(r => {
            csv += [
                r.gene,
                r.n_wt,
                r.mean_wt.toFixed(4),
                r.n_mut,
                r.mean_mut.toFixed(4),
                r.diff_mut.toFixed(4),
                this.formatPValue(r.p_mut),
                r.n_2,
                isNaN(r.mean_2) ? '' : r.mean_2.toFixed(4),
                isNaN(r.diff_2) ? '' : r.diff_2.toFixed(4),
                this.formatPValue(r.p_2)
            ].join(',') + '\n';
        });

        const filename = `mutation_analysis_${mr.hotspotGene}_${new Date().toISOString().slice(0, 10)}.csv`;
        this.downloadFile(csv, filename, 'text/csv');
    }

    showGeneEffectDistribution(gene) {
        if (!this.mutationResults) return;

        const mr = this.mutationResults;
        const hotspotGene = mr.hotspotGene;
        const mutationData = this.mutations.geneData[hotspotGene];
        const geneIdx = this.geneIndex.get(gene.toUpperCase());

        if (geneIdx === undefined) {
            alert(`Gene ${gene} not found`);
            return;
        }

        // Collect data for each cell line
        const cellLines = this.metadata.cellLines;
        const data = { wt: [], mut1: [], mut2: [] };
        this.currentGeneEffectData = []; // Store for CSV export

        cellLines.forEach((cellLine, idx) => {
            // Apply same filters as mutation analysis
            if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) {
                return;
            }

            // Check additional hotspot filter
            if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
                const addMutData = this.mutations.geneData[mr.additionalHotspot];
                if (addMutData) {
                    const addMutLevel = addMutData.mutations[cellLine] || 0;
                    if (mr.additionalHotspotLevel === '0' && addMutLevel !== 0) return;
                    if (mr.additionalHotspotLevel === '1' && addMutLevel !== 1) return;
                    if (mr.additionalHotspotLevel === '2' && addMutLevel < 2) return;
                    if (mr.additionalHotspotLevel === '1+2' && addMutLevel === 0) return;
                }
            }

            const ge = this.geneEffects[geneIdx * this.nCellLines + idx];
            if (isNaN(ge)) return;

            const mutLevel = mutationData.mutations[cellLine] || 0;
            const cellName = this.getCellLineName(cellLine);
            const lineage = this.getCellLineLineage(cellLine);

            const point = { ge, cellLine, cellName, lineage, mutLevel };
            this.currentGeneEffectData.push(point);

            if (mutLevel === 0) {
                data.wt.push(point);
            } else if (mutLevel === 1) {
                data.mut1.push(point);
            } else {
                data.mut2.push(point);
            }
        });

        // Store current gene for downloads
        this.currentGeneEffectGene = gene;

        // Create jitter for y-axis
        const jitter = (base, spread = 0.15) => base + (Math.random() - 0.5) * spread;

        // Build Plotly traces
        const traces = [
            {
                x: data.wt.map(d => d.ge),
                y: data.wt.map(() => jitter(0)),
                mode: 'markers',
                type: 'scatter',
                name: `WT (n=${data.wt.length})`,
                marker: { color: '#888888', size: 8, opacity: 0.7 },
                text: data.wt.map(d => `${d.cellName}<br>${d.lineage}<br>GE: ${d.ge.toFixed(3)}`),
                hoverinfo: 'text'
            },
            {
                x: data.mut1.map(d => d.ge),
                y: data.mut1.map(() => jitter(1)),
                mode: 'markers',
                type: 'scatter',
                name: `1 mutation (n=${data.mut1.length})`,
                marker: { color: '#3b82f6', size: 8, opacity: 0.7 },
                text: data.mut1.map(d => `${d.cellName}<br>${d.lineage}<br>GE: ${d.ge.toFixed(3)}`),
                hoverinfo: 'text'
            },
            {
                x: data.mut2.map(d => d.ge),
                y: data.mut2.map(() => jitter(2)),
                mode: 'markers',
                type: 'scatter',
                name: `2 mutations (n=${data.mut2.length})`,
                marker: { color: '#dc2626', size: 8, opacity: 0.7 },
                text: data.mut2.map(d => `${d.cellName}<br>${d.lineage}<br>GE: ${d.ge.toFixed(3)}`),
                hoverinfo: 'text'
            }
        ];

        // Calculate means for each group
        const meanWT = data.wt.length > 0 ? this.mean(data.wt.map(d => d.ge)) : NaN;
        const meanMut1 = data.mut1.length > 0 ? this.mean(data.mut1.map(d => d.ge)) : NaN;
        const meanMut2 = data.mut2.length > 0 ? this.mean(data.mut2.map(d => d.ge)) : NaN;

        // Add mean lines
        const allGE = [...data.wt, ...data.mut1, ...data.mut2].map(d => d.ge);
        const xMin = Math.min(...allGE) - 0.1;
        const xMax = Math.max(...allGE) + 0.1;

        if (!isNaN(meanWT)) {
            traces.push({
                x: [meanWT, meanWT],
                y: [-0.3, 0.3],
                mode: 'lines',
                line: { color: '#888888', width: 3 },
                showlegend: false,
                hoverinfo: 'skip'
            });
        }
        if (!isNaN(meanMut1)) {
            traces.push({
                x: [meanMut1, meanMut1],
                y: [0.7, 1.3],
                mode: 'lines',
                line: { color: '#3b82f6', width: 3 },
                showlegend: false,
                hoverinfo: 'skip'
            });
        }
        if (!isNaN(meanMut2)) {
            traces.push({
                x: [meanMut2, meanMut2],
                y: [1.7, 2.3],
                mode: 'lines',
                line: { color: '#dc2626', width: 3 },
                showlegend: false,
                hoverinfo: 'skip'
            });
        }

        // Build subtitle with filter info
        let filterInfo = [];
        if (mr.lineageFilter) {
            filterInfo.push(`Lineage: ${mr.lineageFilter}`);
        }
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            filterInfo.push(`${mr.additionalHotspot}: ${mr.additionalHotspotLevel}`);
        }
        const subtitle = filterInfo.length > 0 ? filterInfo.join(' | ') : 'All lineages';

        const layout = {
            title: {
                text: `${gene} Gene Effect by ${hotspotGene} Mutation Status<br><sub style="font-size:12px;color:#666">${subtitle}</sub>`,
                font: { size: 16 }
            },
            xaxis: {
                title: `${gene} Gene Effect`,
                range: [xMin, xMax]
            },
            yaxis: {
                title: `${hotspotGene} Mutations`,
                tickmode: 'array',
                tickvals: [0, 1, 2],
                ticktext: [`0 WT (n=${data.wt.length})`, `1 (n=${data.mut1.length})`, `2 (n=${data.mut2.length})`],
                range: [-0.5, 2.5]
            },
            showlegend: false,
            margin: { t: 70, r: 30, b: 50, l: 120 }
        };

        // Show modal
        document.getElementById('geneEffectModal').style.display = 'flex';
        document.getElementById('geneEffectTitle').textContent = `${gene} Gene Effect by ${hotspotGene} Mutation`;

        Plotly.newPlot('geneEffectPlot', traces, layout, { responsive: true });
    }

    downloadGeneEffectPNG() {
        Plotly.downloadImage('geneEffectPlot', {
            format: 'png',
            width: 900,
            height: 550,
            filename: `gene_effect_${this.currentGeneEffectGene}_${this.mutationResults.hotspotGene}`
        });
    }

    downloadGeneEffectSVG() {
        Plotly.downloadImage('geneEffectPlot', {
            format: 'svg',
            width: 900,
            height: 550,
            filename: `gene_effect_${this.currentGeneEffectGene}_${this.mutationResults.hotspotGene}`
        });
    }

    downloadGeneEffectCSV() {
        if (!this.currentGeneEffectData || !this.currentGeneEffectGene) return;

        const mr = this.mutationResults;
        let csv = `# Gene Effect Distribution Data\n`;
        csv += `# Gene: ${this.currentGeneEffectGene}\n`;
        csv += `# Hotspot Mutation: ${mr.hotspotGene}\n`;
        csv += `# Lineage filter: ${mr.lineageFilter || 'All lineages'}\n`;
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            csv += `# Additional filter: ${mr.additionalHotspot} = ${mr.additionalHotspotLevel}\n`;
        }
        csv += `# Date: ${new Date().toISOString().slice(0, 10)}\n`;
        csv += '#\n';
        csv += 'CellLine,CellLineName,Lineage,GeneEffect,MutationLevel\n';

        this.currentGeneEffectData.forEach(d => {
            csv += `${d.cellLine},${d.cellName},${d.lineage},${d.ge.toFixed(4)},${d.mutLevel}\n`;
        });

        const filename = `gene_effect_${this.currentGeneEffectGene}_${mr.hotspotGene}_data.csv`;
        this.downloadFile(csv, filename, 'text/csv');
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

        // Calculate mean effect for each gene (both all cells and filtered cells)
        const clusterData = clusters.map(gene => {
            const idx = this.geneIndex.get(gene);
            const fullData = this.getGeneData(idx);

            // Stats for ALL cells
            const allValidData = Array.from(fullData).filter(v => !isNaN(v));
            const allMean = allValidData.length > 0 ? allValidData.reduce((a, b) => a + b, 0) / allValidData.length : NaN;
            const allVariance = allValidData.length > 0 ? allValidData.reduce((a, b) => a + (b - allMean) ** 2, 0) / allValidData.length : NaN;
            const allSd = Math.sqrt(allVariance);

            // Stats for FILTERED cells
            const filteredData = cellLineIndices.map(i => fullData[i]).filter(v => !isNaN(v));
            const filtMean = filteredData.length > 0 ? filteredData.reduce((a, b) => a + b, 0) / filteredData.length : NaN;
            const filtVariance = filteredData.length > 0 ? filteredData.reduce((a, b) => a + (b - filtMean) ** 2, 0) / filteredData.length : NaN;
            const filtSd = Math.sqrt(filtVariance);

            return {
                gene: gene,
                cluster: correlations.find(c => c.gene1 === gene || c.gene2 === gene)?.cluster || 0,
                meanEffect: Math.round(allMean * 100) / 100,
                sdEffect: Math.round(allSd * 100) / 100,
                meanEffectFiltered: Math.round(filtMean * 100) / 100,
                sdEffectFiltered: Math.round(filtSd * 100) / 100,
                nAll: allValidData.length,
                nFiltered: filteredData.length,
                inGeneList: geneList.includes(gene)
            };
        });

        // Check if filtering was applied
        const isFiltered = cellLineIndices.length < this.nCellLines;

        return {
            success: true,
            correlations: correlations,
            clusters: clusterData,
            geneList: geneList,
            mode: mode,
            cutoff: cutoff,
            nCellLines: cellLineIndices.length,
            isFiltered: isFiltered
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

        // Build reverse lookup: replacement gene -> original gene name
        const synonymLookup = new Map();
        if (this.synonymsUsed && this.synonymsUsed.length > 0) {
            this.synonymsUsed.forEach(s => {
                synonymLookup.set(s.replacement.toUpperCase(), s.original);
            });
        }

        // Create nodes
        geneSet.forEach(gene => {
            const cluster = this.results.clusters.find(c => c.gene === gene);
            const isInput = this.results.geneList.includes(gene);

            // Check if this gene is a synonym replacement
            const originalName = synonymLookup.get(gene.toUpperCase());
            const isSynonym = !!originalName;

            // Look up stats - try replacement name first, then original name
            let geneStat = this.geneStats?.get(gene);
            if (!geneStat && originalName) {
                geneStat = this.geneStats?.get(originalName);
            }

            // Build title with available information
            let titleLines = [gene];
            if (isSynonym) {
                titleLines.push(`(synonym of ${originalName})`);
            }
            titleLines.push(`GE mean: ${cluster?.meanEffect || 'N/A'}`);
            titleLines.push(`GE SD: ${cluster?.sdEffect || 'N/A'}`);
            if (geneStat?.lfc !== undefined && geneStat?.lfc !== null) {
                titleLines.push(`LFC: ${geneStat.lfc.toFixed(3)}`);
            }
            if (geneStat?.fdr !== undefined && geneStat?.fdr !== null) {
                titleLines.push(`FDR: ${geneStat.fdr.toExponential(2)}`);
            }

            // Add * to label if synonym
            const label = isSynonym ? `${gene}*` : gene;

            nodes.push({
                id: gene,
                label: label,
                size: nodeSize,
                font: { size: fontSize, color: '#333' },
                color: {
                    background: this.results.mode === 'design' ?
                        (isInput ? '#16a34a' : '#86efac') : '#16a34a',
                    border: '#ffffff'
                },
                borderWidth: 2,
                title: titleLines.join('\n'),
                isSynonym: isSynonym,
                originalName: originalName
            });
        });

        // Track if any synonyms are in the network for legend display
        this.hasSynonymsInNetwork = synonymLookup.size > 0 &&
            Array.from(geneSet).some(g => synonymLookup.has(g.toUpperCase()));

        const data = { nodes: new vis.DataSet(nodes), edges: new vis.DataSet(edges) };

        // Adjust stabilization iterations based on network size
        const nodeCount = nodes.length;
        const stabilizationIterations = nodeCount > 50 ? 300 : 150;

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
                    springConstant: 0.08,
                    damping: 0.4
                },
                stabilization: {
                    enabled: true,
                    iterations: stabilizationIterations,
                    updateInterval: 25
                }
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

        // Reset physics state for new network
        this.physicsEnabled = true;
        this.currentLayout = 0;
        const layoutBtn = document.getElementById('changeLayout');
        if (layoutBtn) layoutBtn.textContent = 'Change Layout';
        const physicsBtn = document.getElementById('togglePhysics');
        if (physicsBtn) {
            physicsBtn.textContent = 'Lock Nodes';
            physicsBtn.classList.remove('btn-active');
        }

        // For large networks, disable physics after stabilization
        if (nodeCount > 30) {
            this.network.once('stabilizationIterationsDone', () => {
                this.network.setOptions({ physics: { enabled: false } });
                this.physicsEnabled = false;
                if (physicsBtn) {
                    physicsBtn.textContent = 'Unlock Nodes';
                    physicsBtn.classList.add('btn-active');
                }
            });
        }

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

        // Click on edge to open inspect modal
        this.network.on('click', (params) => {
            if (params.edges.length > 0 && params.nodes.length === 0) {
                const edgeId = params.edges[0];
                const edge = this.networkData.edges.get(edgeId);
                if (edge) {
                    this.openInspectByGenes(edge.from, edge.to);
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
            legendNodeType.style.display = 'block';
        } else {
            legendNodeType.innerHTML = '';
            legendNodeType.style.display = 'none';
        }

        // Show synonym legend if synonyms are used in the network
        const legendSynonym = document.getElementById('legendSynonym');
        if (legendSynonym) {
            legendSynonym.style.display = this.hasSynonymsInNetwork ? 'block' : 'none';
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

        // Deduplicate correlations (A-B is same as B-A)
        const seenPairs = new Set();
        const uniqueCorrelations = this.results.correlations.filter(c => {
            const pairKey = [c.gene1, c.gene2].sort().join('|');
            if (seenPairs.has(pairKey)) {
                return false;
            }
            seenPairs.add(pairKey);
            return true;
        });

        uniqueCorrelations
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

        // Check if we have stats and if filtering was applied
        const hasStats = this.geneStats && this.geneStats.size > 0;
        const isFiltered = this.results.isFiltered;

        // Show filter info
        const filterInfoDiv = document.getElementById('clustersFilterInfo');
        const filterTextSpan = document.getElementById('clustersFilterText');
        if (filterInfoDiv && filterTextSpan) {
            const lineage = document.getElementById('lineageFilter').value;
            const subLineage = document.getElementById('subLineageFilter')?.value;

            if (isFiltered && (lineage || subLineage)) {
                let filterText = `Lineage: ${lineage || 'All'}`;
                if (subLineage) filterText += ` | Subtype: ${subLineage}`;
                filterText += ` | Filtered cells: ${this.results.nCellLines}`;
                filterTextSpan.textContent = filterText;
                filterInfoDiv.style.display = 'block';
            } else {
                filterInfoDiv.style.display = 'none';
            }
        }

        // Build header based on what data we have
        let headerCells = `
            <th data-sort="gene">Gene</th>
            <th data-sort="cluster">Cluster</th>
            <th data-sort="meanEffect">Mean (All)</th>
            <th data-sort="sdEffect">SD (All)</th>
        `;

        if (isFiltered) {
            headerCells += `
            <th data-sort="meanEffectFiltered">Mean (Filt)</th>
            <th data-sort="sdEffectFiltered">SD (Filt)</th>
            `;
        }

        if (hasStats) {
            headerCells += `
            <th data-sort="lfc">LFC</th>
            <th data-sort="fdr">FDR</th>
            `;
        }

        thead.innerHTML = `<tr>${headerCells}</tr>`;

        // Re-attach sorting event listeners
        thead.querySelectorAll('th[data-sort]').forEach(th => {
            th.addEventListener('click', () => this.sortTable(th));
        });

        this.results.clusters
            .sort((a, b) => a.cluster - b.cluster || a.gene.localeCompare(b.gene))
            .forEach(c => {
                const tr = document.createElement('tr');
                const geneStat = this.geneStats?.get(c.gene);

                let rowHtml = `
                    <td>${c.gene}${c.inGeneList && this.results.mode === 'design' ? '*' : ''}</td>
                    <td>${c.cluster}</td>
                    <td>${c.meanEffect}</td>
                    <td>${c.sdEffect}</td>
                `;

                if (isFiltered) {
                    rowHtml += `
                    <td>${c.meanEffectFiltered}</td>
                    <td>${c.sdEffectFiltered}</td>
                    `;
                }

                if (hasStats) {
                    const lfc = geneStat?.lfc !== null && geneStat?.lfc !== undefined
                        ? geneStat.lfc.toFixed(2) : '-';
                    const fdr = geneStat?.fdr !== null && geneStat?.fdr !== undefined
                        ? geneStat.fdr.toExponential(2) : '-';
                    rowHtml += `
                    <td>${lfc}</td>
                    <td>${fdr}</td>
                    `;
                }

                tr.innerHTML = rowHtml;
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
            const isFiltered = this.results.isFiltered;
            const lineage = document.getElementById('lineageFilter').value || 'All lineages';
            const subLineage = document.getElementById('subLineageFilter')?.value;

            // Add filter info as comments
            csv = `# Clusters Export\n`;
            csv += `# Lineage filter: ${lineage}\n`;
            if (subLineage) csv += `# Subtype filter: ${subLineage}\n`;
            csv += `# Filtered cell lines: ${this.results.nCellLines}\n`;
            csv += `# Date: ${new Date().toISOString().slice(0, 10)}\n`;
            csv += '#\n';

            if (isFiltered) {
                csv += 'Gene,Cluster,Mean_Effect_All,SD_Effect_All,Mean_Effect_Filtered,SD_Effect_Filtered\n';
                this.results.clusters.forEach(c => {
                    csv += `${c.gene},${c.cluster},${c.meanEffect},${c.sdEffect},${c.meanEffectFiltered},${c.sdEffectFiltered}\n`;
                });
            } else {
                csv += 'Gene,Cluster,Mean_Effect,SD_Effect\n';
                this.results.clusters.forEach(c => {
                    csv += `${c.gene},${c.cluster},${c.meanEffect},${c.sdEffect}\n`;
                });
            }
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

        // Color by Gene Effect legend
        const colorByGeneEffect = document.getElementById('colorByGeneEffect').checked;
        if (colorByGeneEffect && this.results?.clusters) {
            const colorGEType = document.querySelector('input[name="colorGEType"]:checked')?.value || 'signed';
            const effectValues = this.results.clusters.map(c => c.meanEffect).filter(v => !isNaN(v));

            ctx.font = titleFont;
            ctx.fillStyle = '#333';
            ctx.fillText('Node Color:', legendX, legendY);
            ctx.font = textFont;

            const gradientWidth = 120;
            const gradientHeight = 18;
            const gradY = legendY + 18;

            if (colorGEType === 'signed') {
                const minEffect = Math.min(...effectValues);
                const maxEffect = Math.max(...effectValues);

                // Red (negative) to White (0) to Blue (positive)
                const gradient = ctx.createLinearGradient(legendX, 0, legendX + gradientWidth, 0);
                gradient.addColorStop(0, '#b2182b');
                gradient.addColorStop(0.5, '#f7f7f7');
                gradient.addColorStop(1, '#2166ac');
                ctx.fillStyle = gradient;
                ctx.fillRect(legendX, gradY, gradientWidth, gradientHeight);
                ctx.strokeStyle = '#999';
                ctx.lineWidth = 1;
                ctx.strokeRect(legendX, gradY, gradientWidth, gradientHeight);

                ctx.fillStyle = '#333';
                ctx.font = smallFont;
                ctx.fillText(minEffect.toFixed(2), legendX, gradY + gradientHeight + 16);
                ctx.fillText('Gene Effect (+/−)', legendX, gradY - 4);
                ctx.fillText(maxEffect.toFixed(2), legendX + gradientWidth - 25, gradY + gradientHeight + 16);
            } else {
                const maxAbsEffect = Math.max(...effectValues.map(v => Math.abs(v)));

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
                ctx.fillText('|Gene Effect|', legendX, gradY - 4);
                ctx.fillText(maxAbsEffect.toFixed(2), legendX + gradientWidth - 25, gradY + gradientHeight + 16);
            }
            legendX += 170;
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

                // Draw gradient - Red (negative) to White (0) to Blue (positive)
                const gradient = ctx.createLinearGradient(legendX, 0, legendX + gradientWidth, 0);
                gradient.addColorStop(0, '#b2182b');
                gradient.addColorStop(0.5, '#f7f7f7');
                gradient.addColorStop(1, '#2166ac');
                ctx.fillStyle = gradient;
                ctx.fillRect(legendX, gradY, gradientWidth, gradientHeight);
                ctx.strokeStyle = '#999';
                ctx.lineWidth = 1;
                ctx.strokeRect(legendX, gradY, gradientWidth, gradientHeight);

                // Labels
                ctx.fillStyle = '#333';
                ctx.font = smallFont;
                ctx.fillText(minLfc.toFixed(1), legendX, gradY + gradientHeight + 16);
                ctx.fillText('LFC (+/−)' + scaleLabel, legendX, gradY - 4);
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
            legendX += 170;
        }

        // Synonym legend
        if (this.hasSynonymsInNetwork) {
            ctx.font = textFont;
            ctx.fillStyle = '#333';
            ctx.fillText('* = synonym/orthologue used', legendX, legendY + 25);
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
    <linearGradient id="signedGradient" x1="0%" y1="0%" x2="100%" y2="0%">
        <stop offset="0%" style="stop-color:#b2182b;stop-opacity:1" />
        <stop offset="50%" style="stop-color:#f7f7f7;stop-opacity:1" />
        <stop offset="100%" style="stop-color:#2166ac;stop-opacity:1" />
    </linearGradient>
    <linearGradient id="absGradient" x1="0%" y1="0%" x2="100%" y2="0%">
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

        // Color by Gene Effect legend
        const colorByGeneEffect = document.getElementById('colorByGeneEffect').checked;
        if (colorByGeneEffect && this.results?.clusters) {
            const colorGEType = document.querySelector('input[name="colorGEType"]:checked')?.value || 'signed';
            const effectValues = this.results.clusters.map(c => c.meanEffect).filter(v => !isNaN(v));

            svg += `  <text x="${legendX}" y="${legendY}" class="legend-title">Node Color:</text>\n`;

            const gradientWidth = 120;
            const gradientHeight = 18;
            const gradY = legendY + 18;

            if (colorGEType === 'signed') {
                const minEffect = Math.min(...effectValues);
                const maxEffect = Math.max(...effectValues);

                svg += `  <rect x="${legendX}" y="${gradY}" width="${gradientWidth}" height="${gradientHeight}" fill="url(#signedGradient)" stroke="#999"/>\n`;
                svg += `  <text x="${legendX}" y="${gradY + gradientHeight + 16}" class="legend-small">${minEffect.toFixed(2)}</text>\n`;
                svg += `  <text x="${legendX}" y="${gradY - 4}" class="legend-small">Gene Effect (+/−)</text>\n`;
                svg += `  <text x="${legendX + gradientWidth - 25}" y="${gradY + gradientHeight + 16}" class="legend-small">${maxEffect.toFixed(2)}</text>\n`;
            } else {
                const maxAbsEffect = Math.max(...effectValues.map(v => Math.abs(v)));

                svg += `  <rect x="${legendX}" y="${gradY}" width="${gradientWidth}" height="${gradientHeight}" fill="url(#absGradient)" stroke="#999"/>\n`;
                svg += `  <text x="${legendX}" y="${gradY + gradientHeight + 16}" class="legend-small">0</text>\n`;
                svg += `  <text x="${legendX}" y="${gradY - 4}" class="legend-small">|Gene Effect|</text>\n`;
                svg += `  <text x="${legendX + gradientWidth - 25}" y="${gradY + gradientHeight + 16}" class="legend-small">${maxAbsEffect.toFixed(2)}</text>\n`;
            }
            legendX += 170;
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

                svg += `  <rect x="${legendX}" y="${gradY}" width="${gradientWidth}" height="${gradientHeight}" fill="url(#signedGradient)" stroke="#999"/>\n`;
                svg += `  <text x="${legendX}" y="${gradY + gradientHeight + 16}" class="legend-small">${minLfc.toFixed(1)}</text>\n`;
                svg += `  <text x="${legendX}" y="${gradY - 4}" class="legend-small">LFC (+/−)${scaleLabel}</text>\n`;
                svg += `  <text x="${legendX + gradientWidth - 20}" y="${gradY + gradientHeight + 16}" class="legend-small">${maxLfc.toFixed(1)}</text>\n`;
            } else if (colorStatType === 'abs_lfc') {
                const lfcValues = stats.map(s => Math.abs(s.lfc)).filter(v => v !== null && !isNaN(v));
                const maxLfc = Math.max(...lfcValues);

                svg += `  <rect x="${legendX}" y="${gradY}" width="${gradientWidth}" height="${gradientHeight}" fill="url(#absGradient)" stroke="#999"/>\n`;
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
            legendX += 170;
        }

        // Synonym legend
        if (this.hasSynonymsInNetwork) {
            svg += `  <text x="${legendX}" y="${legendY + 25}" class="legend-text">* = synonym/orthologue used</text>\n`;
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
        const showSD = document.getElementById('showGeneEffectSD').checked;
        const updates = [];

        this.networkData.nodes.forEach(node => {
            const cluster = this.results?.clusters?.find(c => c.gene === node.id);
            // Add * suffix if this is a synonym
            const baseName = node.isSynonym ? `${node.id}*` : node.id;
            let label = baseName;

            if (showGE && cluster) {
                if (showSD && cluster.sdEffect) {
                    label = `${baseName}\n(GE:${cluster.meanEffect}±${cluster.sdEffect})`;
                } else {
                    label = `${baseName}\n(GE:${cluster.meanEffect})`;
                }
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
        const showSD = document.getElementById('showGeneEffectSD').checked;
        const updates = [];

        // Build reverse lookup for synonyms to find stats from original name
        const synonymLookup = new Map();
        if (this.synonymsUsed && this.synonymsUsed.length > 0) {
            this.synonymsUsed.forEach(s => {
                synonymLookup.set(s.replacement.toUpperCase(), s.original);
            });
        }

        this.networkData.nodes.forEach(node => {
            const cluster = this.results?.clusters?.find(c => c.gene === node.id);

            // Look up stats - try replacement name first, then original name
            let geneStat = this.geneStats?.get(node.id);
            const originalName = synonymLookup.get(node.id.toUpperCase());
            if (!geneStat && originalName) {
                geneStat = this.geneStats?.get(originalName);
            }

            // Add * suffix if this is a synonym
            const baseName = node.isSynonym ? `${node.id}*` : node.id;
            let label = baseName;

            // Add gene effect if checked
            if (showGE && cluster) {
                if (showSD && cluster.sdEffect) {
                    label = `${baseName}\n(GE:${cluster.meanEffect}±${cluster.sdEffect})`;
                } else {
                    label = `${baseName}\n(GE:${cluster.meanEffect})`;
                }
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
        const colorByGeneEffect = document.getElementById('colorByGeneEffect').checked;
        const colorStatType = document.querySelector('input[name="colorStatType"]:checked')?.value || 'signed_lfc';
        const colorGEType = document.querySelector('input[name="colorGEType"]:checked')?.value || 'signed';
        const colorScale = document.querySelector('input[name="colorScale"]:checked')?.value || 'all';

        const updates = [];
        const colorLegend = document.getElementById('nodeColorLegendContent');
        const legendSection = document.getElementById('legendNodeColor');

        // Color by gene effect (from DepMap data) - takes precedence
        if (colorByGeneEffect && this.results?.clusters) {
            if (legendSection) legendSection.style.display = 'block';

            // Build map of gene -> meanEffect
            const effectMap = new Map();
            this.results.clusters.forEach(c => effectMap.set(c.gene, c.meanEffect));

            const effectValues = this.results.clusters.map(c => c.meanEffect).filter(v => !isNaN(v));

            if (colorGEType === 'signed') {
                const minEffect = Math.min(...effectValues);
                const maxEffect = Math.max(...effectValues);
                const maxAbs = Math.max(Math.abs(minEffect), Math.abs(maxEffect));

                this.networkData.nodes.forEach(node => {
                    const effect = effectMap.get(node.id);
                    let bgColor = '#cccccc';

                    if (effect !== undefined && !isNaN(effect)) {
                        // Red (negative) to White (0) to Blue (positive)
                        const normalized = (effect + maxAbs) / (2 * maxAbs);
                        bgColor = this.interpolateColor('#b2182b', '#f7f7f7', '#2166ac', normalized);
                    }

                    updates.push({
                        id: node.id,
                        color: { background: bgColor, border: '#ffffff' }
                    });
                });

                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">Gene Effect (+/−)</div>
                    <div style="display: flex; align-items: center; gap: 4px;">
                        <span style="font-size: 10px;">${minEffect.toFixed(2)}</span>
                        <div style="width: 80px; height: 12px; background: linear-gradient(to right, #b2182b, #f7f7f7, #2166ac); border-radius: 2px;"></div>
                        <span style="font-size: 10px;">${maxEffect.toFixed(2)}</span>
                    </div>
                `;
            } else {
                // Absolute gene effect
                const maxAbsEffect = Math.max(...effectValues.map(v => Math.abs(v)));

                this.networkData.nodes.forEach(node => {
                    const effect = effectMap.get(node.id);
                    let bgColor = '#cccccc';

                    if (effect !== undefined && !isNaN(effect)) {
                        const normalized = Math.abs(effect) / maxAbsEffect;
                        bgColor = this.interpolateColor('#f5f5f5', '#fdae61', '#d7191c', normalized);
                    }

                    updates.push({
                        id: node.id,
                        color: { background: bgColor, border: '#ffffff' }
                    });
                });

                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">|Gene Effect|</div>
                    <div style="display: flex; align-items: center; gap: 4px;">
                        <span style="font-size: 10px;">0</span>
                        <div style="width: 80px; height: 12px; background: linear-gradient(to right, #f5f5f5, #fdae61, #d7191c); border-radius: 2px;"></div>
                        <span style="font-size: 10px;">${maxAbsEffect.toFixed(2)}</span>
                    </div>
                `;
            }

            this.networkData.nodes.update(updates);
            return;
        }

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
                        // Red (-) to White (0) to Blue (+)
                        const normalized = (geneStat.lfc + maxAbs) / (2 * maxAbs);
                        bgColor = this.interpolateColor('#b2182b', '#f7f7f7', '#2166ac', normalized);
                    }

                    updates.push({
                        id: node.id,
                        color: { background: bgColor, border: '#ffffff' }
                    });
                });

                const scaleLabel = colorScale === 'network' ? ' (network)' : ' (all genes)';
                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">LFC (+/−)${scaleLabel}</div>
                    <div style="display: flex; align-items: center; gap: 4px;">
                        <span style="font-size: 10px;">${minLfc.toFixed(1)}</span>
                        <div style="width: 80px; height: 12px; background: linear-gradient(to right, #b2182b, #f7f7f7, #2166ac); border-radius: 2px;"></div>
                        <span style="font-size: 10px;">${maxLfc.toFixed(1)}</span>
                    </div>
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

    togglePhysics() {
        if (!this.network) return;

        this.physicsEnabled = !this.physicsEnabled;
        this.network.setOptions({ physics: { enabled: this.physicsEnabled } });

        const btn = document.getElementById('togglePhysics');
        if (this.physicsEnabled) {
            btn.textContent = 'Lock Nodes';
            btn.classList.remove('btn-active');
        } else {
            btn.textContent = 'Unlock Nodes';
            btn.classList.add('btn-active');
        }
    }

    changeNetworkLayout() {
        if (!this.network) return;

        // Cycle through layout options
        const layouts = ['forceAtlas2Based', 'barnesHut', 'repulsion', 'hierarchical'];
        this.currentLayout = this.currentLayout || 0;
        this.currentLayout = (this.currentLayout + 1) % layouts.length;

        const layoutName = layouts[this.currentLayout];

        let options;
        if (layoutName === 'hierarchical') {
            options = {
                physics: { enabled: false },
                layout: {
                    hierarchical: {
                        enabled: true,
                        direction: 'UD',
                        sortMethod: 'hubsize',
                        nodeSpacing: 150,
                        levelSeparation: 150
                    }
                }
            };
        } else {
            options = {
                layout: { hierarchical: { enabled: false } },
                physics: {
                    enabled: true,
                    solver: layoutName,
                    stabilization: { iterations: 150 }
                }
            };
        }

        this.network.setOptions(options);
        this.physicsEnabled = layoutName !== 'hierarchical';

        // Update physics button state
        const btn = document.getElementById('togglePhysics');
        if (this.physicsEnabled) {
            btn.textContent = 'Lock Nodes';
            btn.classList.remove('btn-active');
        } else {
            btn.textContent = 'Unlock Nodes';
            btn.classList.add('btn-active');
        }

        // Show which layout is active
        const layoutNames = {
            'forceAtlas2Based': 'Force Atlas',
            'barnesHut': 'Barnes Hut',
            'repulsion': 'Repulsion',
            'hierarchical': 'Hierarchical'
        };
        document.getElementById('changeLayout').textContent = layoutNames[layoutName];
    }

    resetNetworkSettings() {
        // Reset checkboxes
        document.getElementById('showGeneEffect').checked = false;
        document.getElementById('showGeneEffectSD').checked = false;
        document.getElementById('colorByGeneEffect').checked = false;
        document.getElementById('colorByStats').checked = false;

        // Reset visibility of sub-options
        document.getElementById('showGESDGroup').style.display = 'none';
        document.getElementById('colorGEOptions').style.display = 'none';
        document.getElementById('colorStatsOptions').style.display = 'none';

        // Reset radio buttons
        document.querySelector('input[name="colorGEType"][value="signed"]').checked = true;
        document.querySelector('input[name="colorStatType"][value="signed_lfc"]').checked = true;
        document.querySelector('input[name="colorScale"][value="all"]').checked = true;
        document.querySelector('input[name="statsLabelDisplay"][value="none"]').checked = true;

        // Reset layout state
        this.currentLayout = 0;
        this.physicsEnabled = true;
        const layoutBtn = document.getElementById('changeLayout');
        if (layoutBtn) layoutBtn.textContent = 'Change Layout';
        const physicsBtn = document.getElementById('togglePhysics');
        if (physicsBtn) {
            physicsBtn.textContent = 'Lock Nodes';
            physicsBtn.classList.remove('btn-active');
        }

        // Reset tabs to show appropriate ones based on mode
        const mode = document.querySelector('input[name="analysisMode"]:checked').value;
        const mutationTab = document.getElementById('mutationTab');
        if (mode === 'mutation') {
            mutationTab.style.display = 'inline-block';
        } else {
            mutationTab.style.display = 'none';
            // If mutation tab was active, switch to network tab
            if (mutationTab.classList.contains('active')) {
                document.querySelectorAll('.nav-link').forEach(l => l.classList.remove('active'));
                document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
                document.querySelector('[data-tab="network"]').classList.add('active');
                document.getElementById('tab-network').classList.add('active');
            }
        }
    }

    downloadAllData() {
        if (!this.results) return;

        // Create correlations CSV
        let correlationsCSV = 'Gene1,Gene2,Correlation,Slope,N,Cluster\n';
        this.results.correlations.forEach(c => {
            correlationsCSV += `${c.gene1},${c.gene2},${c.correlation},${c.slope},${c.n},${c.cluster}\n`;
        });

        // Create clusters CSV
        let clustersCSV;
        if (this.results.isFiltered) {
            clustersCSV = 'Gene,Cluster,Mean_Effect_All,SD_Effect_All,Mean_Effect_Filtered,SD_Effect_Filtered\n';
            this.results.clusters.forEach(c => {
                clustersCSV += `${c.gene},${c.cluster},${c.meanEffect},${c.sdEffect},${c.meanEffectFiltered},${c.sdEffectFiltered}\n`;
            });
        } else {
            clustersCSV = 'Gene,Cluster,Mean_Effect,SD_Effect\n';
            this.results.clusters.forEach(c => {
                clustersCSV += `${c.gene},${c.cluster},${c.meanEffect},${c.sdEffect}\n`;
            });
        }

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

        // Create ZIP file using JSZip
        if (typeof JSZip === 'undefined') {
            // Fallback if JSZip not loaded
            console.warn('JSZip not loaded, downloading files separately');
            this.downloadFile(correlationsCSV, 'correlations.csv', 'text/csv');
            setTimeout(() => this.downloadFile(clustersCSV, 'clusters.csv', 'text/csv'), 300);
            setTimeout(() => this.downloadFile(summary, 'summary.txt', 'text/plain'), 600);
            return;
        }

        const zip = new JSZip();
        zip.file('correlations.csv', correlationsCSV);
        zip.file('clusters.csv', clustersCSV);
        zip.file('summary.txt', summary);

        // Add network images if network exists
        const addNetworkImages = async () => {
            if (this.network && this.networkData) {
                // Get PNG as base64
                const pngData = await this.getNetworkPNGData();
                if (pngData) {
                    // Remove data URL prefix to get just base64
                    const base64 = pngData.split(',')[1];
                    zip.file('network.png', base64, { base64: true });
                }

                // Get SVG
                const svgData = this.getNetworkSVGData();
                if (svgData) {
                    zip.file('network.svg', svgData);
                }
            }

            // Generate and download ZIP
            const content = await zip.generateAsync({ type: 'blob' });
            const url = URL.createObjectURL(content);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'correlation_analysis.zip';
            a.click();
            URL.revokeObjectURL(url);
        };

        addNetworkImages();
    }

    // Helper to get network PNG data for ZIP
    async getNetworkPNGData() {
        const networkCanvas = document.querySelector('#networkPlot canvas');
        if (!networkCanvas) return null;

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

        ctx.fillStyle = 'white';
        ctx.fillRect(0, 0, totalWidth, totalHeight);
        ctx.drawImage(networkCanvas, 0, 0);

        // Draw legend (simplified version)
        ctx.fillStyle = '#f9fafb';
        ctx.fillRect(15, networkHeight + 10, totalWidth - 30, legendHeight - 10);
        ctx.strokeStyle = '#e5e7eb';
        ctx.strokeRect(15, networkHeight + 10, totalWidth - 30, legendHeight - 10);

        const legendY = networkHeight + padding + 10;
        ctx.font = 'bold 14px Arial';
        ctx.fillStyle = '#333';
        ctx.fillText('Correlation: Blue=Positive, Red=Negative', 40, legendY);
        ctx.fillText(`Cutoff: ${this.results?.cutoff || 0.5}`, 40, legendY + 25);

        if (this.hasSynonymsInNetwork) {
            ctx.fillText('* = synonym/orthologue used', 40, legendY + 50);
        }

        return canvas.toDataURL('image/png');
    }

    // Helper to get network SVG data for ZIP
    getNetworkSVGData() {
        if (!this.network || !this.networkData) return null;

        const container = document.getElementById('networkPlot');
        const width = container.clientWidth;
        const networkHeight = container.clientHeight;
        const positions = this.network.getPositions();

        let svg = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${networkHeight}" viewBox="0 0 ${width} ${networkHeight}">
<rect width="100%" height="100%" fill="white"/>
`;

        // Draw edges
        this.networkData.edges.forEach(edge => {
            const from = positions[edge.from];
            const to = positions[edge.to];
            if (from && to) {
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
        const nodeSize = parseInt(document.getElementById('netNodeSize')?.value) || 25;
        this.networkData.nodes.forEach(node => {
            const pos = positions[node.id];
            if (pos) {
                const cx = pos.x + width/2;
                const cy = pos.y + networkHeight/2;
                const bgColor = node.color?.background || '#16a34a';
                svg += `  <circle cx="${cx}" cy="${cy}" r="${nodeSize/2}" fill="${bgColor}" stroke="white" stroke-width="2"/>\n`;
                svg += `  <text x="${cx}" y="${cy + nodeSize/2 + 14}" text-anchor="middle" style="font-family: Arial; font-size: 12px; fill: #333;">${this.escapeXml(node.label || node.id)}</text>\n`;
            }
        });

        svg += '</svg>';
        return svg;
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

        // Set axis limits with 10% padding on each side
        const xVals = plotData.map(d => d.x);
        const yVals = plotData.map(d => d.y);
        const xMin = Math.min(...xVals);
        const xMax = Math.max(...xVals);
        const yMin = Math.min(...yVals);
        const yMax = Math.max(...yVals);
        const xPadding = (xMax - xMin) * 0.1;
        const yPadding = (yMax - yMin) * 0.1;
        this.currentInspect.defaultXlim = [xMin - xPadding, xMax + xPadding];
        this.currentInspect.defaultYlim = [yMin - yPadding, yMax + yPadding];

        document.getElementById('scatterXmin').value = this.currentInspect.defaultXlim[0].toFixed(1);
        document.getElementById('scatterXmax').value = this.currentInspect.defaultXlim[1].toFixed(1);
        document.getElementById('scatterYmin').value = this.currentInspect.defaultYlim[0].toFixed(1);
        document.getElementById('scatterYmax').value = this.currentInspect.defaultYlim[1].toFixed(1);

        // Populate cancer filter with counts
        const cancerFilter = document.getElementById('scatterCancerFilter');
        const cancerBox = document.getElementById('cancerFilterBox');
        const lineageCounts = {};
        plotData.forEach(d => {
            if (d.lineage) {
                lineageCounts[d.lineage] = (lineageCounts[d.lineage] || 0) + 1;
            }
        });
        const lineages = Object.keys(lineageCounts).sort();
        if (lineages.length > 0) {
            cancerFilter.innerHTML = `<option value="">All cancer types (n=${plotData.length})</option>`;
            lineages.forEach(l => {
                cancerFilter.innerHTML += `<option value="${l}">${l} (n=${lineageCounts[l]})</option>`;
            });
            cancerBox.style.display = 'block';
        } else {
            cancerBox.style.display = 'none';
        }

        // Populate hotspot genes (excluding HLA-A and HLA-B which have high variability)
        // Count mutations based on cell lines with valid data (plotData)
        const hotspotSelect = document.getElementById('hotspotGene');
        const mutFilterGeneSelect = document.getElementById('mutationFilterGene');
        const excludedGenes = ['HLA-A', 'HLA-B'];
        const cellLinesInPlot = new Set(plotData.map(d => d.cellLineId));

        if (this.mutations?.genes?.length > 0) {
            hotspotSelect.innerHTML = '<option value="">Select gene...</option>';
            mutFilterGeneSelect.innerHTML = '<option value="">No filter</option>';
            this.mutations.genes
                .filter(g => !excludedGenes.includes(g))
                .forEach(g => {
                    // Count mutations only in cell lines with valid data
                    const mutData = this.mutations.geneData?.[g]?.mutations || {};
                    let count = 0;
                    cellLinesInPlot.forEach(cl => {
                        if (mutData[cl] && mutData[cl] > 0) count++;
                    });
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
                if (mutFilterLevel === '1') return mutLevel === 1;
                if (mutFilterLevel === '2') return mutLevel >= 2;
                if (mutFilterLevel === '1+2') return mutLevel >= 1;
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

        // Build filter description for title
        let filterParts = [];
        if (cancerFilter) {
            filterParts.push(`Cancer: ${cancerFilter}`);
        }
        if (mutFilterGene && mutFilterLevel !== 'all') {
            const levelText = mutFilterLevel === '0' ? 'WT' :
                              mutFilterLevel === '1' ? '1 mut' :
                              mutFilterLevel === '2' ? '2 mut' : '1+2 mut';
            filterParts.push(`${mutFilterGene}: ${levelText}`);
        }
        const filterDesc = filterParts.length > 0 ? filterParts.join(' | ') : '';

        // Show/hide plot and table based on mode
        const scatterPlot = document.getElementById('scatterPlot');
        const compareTable = document.getElementById('compareTable');

        if (hotspotMode === 'compare_table' && hotspotGene) {
            scatterPlot.style.display = 'none';
            compareTable.style.display = 'block';
            this.renderCompareTable(filteredData, gene1, gene2, hotspotGene, filterDesc);
            return;
        } else {
            scatterPlot.style.display = 'block';
            compareTable.style.display = 'none';
        }

        // Handle 3-panel mode
        if (hotspotMode === 'three_panel' && hotspotGene) {
            this.renderThreePanelPlot(filteredData, gene1, gene2, hotspotGene, searchTerms, fontSize, filterDesc);
            return;
        }

        // Single panel color mode
        this.renderSinglePanelPlot(filteredData, gene1, gene2, hotspotGene, hotspotMode, searchTerms, fontSize, filterDesc);
    }

    renderSinglePanelPlot(filteredData, gene1, gene2, hotspotGene, hotspotMode, searchTerms, fontSize, filterDesc = '') {
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

        // Calculate means and medians
        const meanX = filteredData.reduce((a, d) => a + d.x, 0) / filteredData.length;
        const meanY = filteredData.reduce((a, d) => a + d.y, 0) / filteredData.length;
        const medianX = this.median(filteredData.map(d => d.x));
        const medianY = this.median(filteredData.map(d => d.y));

        // Build title with all stats (like 3-panel style)
        let titleLines = [`<b>${gene1} vs ${gene2}</b>`];
        if (filterDesc) {
            titleLines.push(`<span style="font-size:11px;color:#666;">Filter: ${filterDesc}</span>`);
        }
        titleLines.push(`<span style="font-size:11px;">n=${filteredData.length}, r=${allStats.correlation.toFixed(3)}, slope=${allStats.slope.toFixed(3)}</span>`);
        titleLines.push(`<span style="font-size:10px;">mean: x=${meanX.toFixed(2)}, y=${meanY.toFixed(2)} | median: x=${medianX.toFixed(2)}, y=${medianY.toFixed(2)}</span>`);

        if (hotspotMode === 'color' && hotspotGene) {
            titleLines.push(`<span style="font-size:10px;"><b>${hotspotGene}:</b> WT: n=${wt.length}, r=${wtStats.correlation.toFixed(3)}, slope=${wtStats.slope.toFixed(3)} | ` +
                `1mut: n=${mut1.length}, r=${mut1Stats.correlation.toFixed(3)}, slope=${mut1Stats.slope.toFixed(3)} | ` +
                `2mut: n=${mut2.length}, r=${mut2Stats.correlation.toFixed(3)}, slope=${mut2Stats.slope.toFixed(3)}</span>`);
        }

        const titleText = titleLines.join('<br>');
        const annotations = [];

        // Calculate margin based on title lines
        const topMargin = 80 + (titleLines.length * 18);

        const layout = {
            title: {
                text: titleText,
                x: 0.5,
                y: 0.98,
                yanchor: 'top',
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
            margin: { t: topMargin, r: 150, b: 60, l: 60 },
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

    renderThreePanelPlot(filteredData, gene1, gene2, hotspotGene, searchTerms, fontSize, filterDesc = '') {
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

        // Build title with filter info
        let titleText = `<b>${gene1} vs ${gene2} - ${hotspotGene} hotspot mutation stratification</b>`;
        if (filterDesc) {
            titleText += `<br><span style="font-size: 11px; color: #666;">Filter: ${filterDesc}</span>`;
        }

        const layout = {
            title: {
                text: titleText,
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
                  text: `<b>WT (0 mut)</b> n=${wt.length}<br>r=${wtStats.correlation.toFixed(3)}, slope=${wtStats.slope.toFixed(3)}<br>mean: x=${wtExtra.meanX.toFixed(2)}, y=${wtExtra.meanY.toFixed(2)}<br>median: x=${wtExtra.medianX.toFixed(2)}, y=${wtExtra.medianY.toFixed(2)}`,
                  showarrow: false, font: { size: 9 } },
                { x: 0.5, y: 1.02, xref: 'paper', yref: 'paper',
                  text: `<b>1 mutation</b> n=${mut1.length}<br>r=${mut1Stats.correlation.toFixed(3)}, slope=${mut1Stats.slope.toFixed(3)}<br>mean: x=${mut1Extra.meanX.toFixed(2)}, y=${mut1Extra.meanY.toFixed(2)}<br>median: x=${mut1Extra.medianX.toFixed(2)}, y=${mut1Extra.medianY.toFixed(2)}`,
                  showarrow: false, font: { size: 9 } },
                { x: 0.86, y: 1.02, xref: 'paper', yref: 'paper',
                  text: `<b>2 mutations</b> n=${mut2.length}<br>r=${mut2Stats.correlation.toFixed(3)}, slope=${mut2Stats.slope.toFixed(3)}<br>mean: x=${mut2Extra.meanX.toFixed(2)}, y=${mut2Extra.meanY.toFixed(2)}<br>median: x=${mut2Extra.medianX.toFixed(2)}, y=${mut2Extra.medianY.toFixed(2)}`,
                  showarrow: false, font: { size: 9 } }
            ],
            margin: { t: filterDesc ? 160 : 140, r: 30, b: 60, l: 60 },
            plot_bgcolor: '#fafafa'
        };

        Plotly.newPlot('scatterPlot', traces, layout, { responsive: true });
    }

    renderCompareTable(filteredData, gene1, gene2, hotspotGene, filterDesc = '') {
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
        const filterInfo = filterDesc ? `<p style="font-size: 11px; color: #333; margin-bottom: 8px; background: #f0f9ff; padding: 4px 8px; border-radius: 4px;"><b>Filter:</b> ${filterDesc}</p>` : '';
        let html = `
            <h4 style="margin-bottom: 8px;">Effect of <span style="color: #0066cc;">${hotspotGene}</span> Mutation on ${gene1} vs ${gene2} Correlation</h4>
            ${filterInfo}
            <p style="font-size: 12px; color: #666; margin-bottom: 12px;">
                Comparing correlation between WT (0 ${hotspotGene} mutations) vs Mutant (2+ ${hotspotGene} mutations) cells, stratified by cancer type.
                Note: Cells with exactly 1 mutation are excluded from this comparison.
            </p>
            <div style="overflow-x: auto;">
            <table id="compareByCancerTable" class="data-table" style="width: 100%; font-size: 12px;">
                <thead>
                    <tr>
                        <th data-col="0" style="cursor: pointer;">Cancer Type ▼</th>
                        <th data-col="1" style="cursor: pointer;">N (WT)</th>
                        <th data-col="2" style="cursor: pointer;">r (WT)</th>
                        <th data-col="3" style="cursor: pointer;">slope (WT)</th>
                        <th data-col="4" style="cursor: pointer;">N (Mut)</th>
                        <th data-col="5" style="cursor: pointer;">r (Mut)</th>
                        <th data-col="6" style="cursor: pointer;">slope (Mut)</th>
                        <th data-col="7" style="cursor: pointer;">Δr</th>
                        <th data-col="8" style="cursor: pointer;">p(Δr)</th>
                        <th data-col="9" style="cursor: pointer;">Δslope</th>
                        <th data-col="10" style="cursor: pointer;">p(Δslope)</th>
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
                <button class="btn btn-primary btn-sm" id="backToGraphBtn2" style="margin-right: 8px;">← Back to Graph</button>
                <button class="btn btn-success btn-sm" id="downloadCompareCSV">Download CSV</button>
            </div>
        `;

        document.getElementById('compareTable').innerHTML = html;

        // Add back to graph handler
        document.getElementById('backToGraphBtn2')?.addEventListener('click', () => {
            document.getElementById('compareTable').style.display = 'none';
            document.getElementById('scatterPlot').style.display = 'block';
        });

        // Add download handler
        document.getElementById('downloadCompareCSV')?.addEventListener('click', () => {
            let csv = 'Cancer Type,N (WT),r (WT),slope (WT),N (Mut),r (Mut),slope (Mut),Δr,p(Δr),Δslope,p(Δslope)\n';
            tableData.forEach(row => {
                csv += `"${row.lineage}",${row.nWT},${row.rWT.toFixed(4)},${row.slopeWT.toFixed(4)},${row.nMut},${row.rMut.toFixed(4)},${row.slopeMut.toFixed(4)},${row.deltaR.toFixed(4)},${row.pR.toExponential(2)},${row.deltaSlope.toFixed(4)},${row.pSlope.toExponential(2)}\n`;
            });
            this.downloadFile(csv, `${gene1}_vs_${gene2}_${hotspotGene}_mutation_comparison.csv`, 'text/csv');
        });

        // Make table sortable
        this.setupSortableTable('compareByCancerTable');
    }

    showCompareAllMutations() {
        if (!this.currentInspect) return;

        const { gene1, gene2, data } = this.currentInspect;

        // Apply current filters
        const cancerFilter = document.getElementById('scatterCancerFilter').value;
        const mutFilterGene = document.getElementById('mutationFilterGene').value;
        const mutFilterLevel = document.getElementById('mutationFilterLevel').value;

        let filteredData = cancerFilter ?
            data.filter(d => d.lineage === cancerFilter) : data;

        // Apply mutation filter
        if (mutFilterGene && this.mutations?.geneData?.[mutFilterGene] && mutFilterLevel !== 'all') {
            const filterMutations = this.mutations.geneData[mutFilterGene].mutations;
            filteredData = filteredData.filter(d => {
                const mutLevel = filterMutations[d.cellLineId] || 0;
                if (mutFilterLevel === '0') return mutLevel === 0;
                if (mutFilterLevel === '1') return mutLevel === 1;
                if (mutFilterLevel === '2') return mutLevel >= 2;
                if (mutFilterLevel === '1+2') return mutLevel >= 1;
                return true;
            });
        }

        // Build filter description
        let filterParts = [];
        if (cancerFilter) {
            filterParts.push(`Cancer: ${cancerFilter}`);
        }
        if (mutFilterGene && mutFilterLevel !== 'all') {
            const levelText = mutFilterLevel === '0' ? 'WT' :
                              mutFilterLevel === '1' ? '1 mut' :
                              mutFilterLevel === '2' ? '2 mut' : '1+2 mut';
            filterParts.push(`${mutFilterGene}: ${levelText}`);
        }
        const filterDesc = filterParts.length > 0 ? filterParts.join(' | ') : '';

        // Show compare table
        document.getElementById('scatterPlot').style.display = 'none';
        document.getElementById('compareTable').style.display = 'block';

        this.renderMutationComparisonTable(filteredData, gene1, gene2, filterDesc);
    }

    renderMutationComparisonTable(filteredData, gene1, gene2, filterDesc = '') {
        // Compare how different hotspot mutations affect the correlation
        if (!this.mutations || !this.mutations.geneData) {
            document.getElementById('compareTable').innerHTML = '<p>No mutation data available.</p>';
            return;
        }

        const tableData = [];
        const mutationGenes = Object.keys(this.mutations.geneData).sort();

        mutationGenes.forEach(mutGene => {
            const mutations = this.mutations.geneData[mutGene].mutations;

            // Split data by mutation status for this gene
            const wt = filteredData.filter(d => (mutations[d.cellLineId] || 0) === 0);
            const mut2 = filteredData.filter(d => (mutations[d.cellLineId] || 0) >= 2);

            // Need at least 3 samples in each group
            if (wt.length >= 3 && mut2.length >= 3) {
                const wtStats = this.pearsonWithSlope(wt.map(d => d.x), wt.map(d => d.y));
                const mutStats = this.pearsonWithSlope(mut2.map(d => d.x), mut2.map(d => d.y));

                const deltaR = mutStats.correlation - wtStats.correlation;
                const deltaSlope = mutStats.slope - wtStats.slope;

                // Fisher z-transformation for correlation difference p-value
                const z1 = 0.5 * Math.log((1 + wtStats.correlation) / (1 - wtStats.correlation));
                const z2 = 0.5 * Math.log((1 + mutStats.correlation) / (1 - mutStats.correlation));
                const se = Math.sqrt(1/(wt.length - 3) + 1/(mut2.length - 3));
                const zDiff = (z2 - z1) / se;
                const pR = 2 * (1 - this.normalCDF(Math.abs(zDiff)));

                tableData.push({
                    mutGene,
                    nWT: wt.length,
                    rWT: wtStats.correlation,
                    slopeWT: wtStats.slope,
                    nMut: mut2.length,
                    rMut: mutStats.correlation,
                    slopeMut: mutStats.slope,
                    deltaR,
                    pR,
                    deltaSlope
                });
            }
        });

        // Sort by p-value (most significant first)
        tableData.sort((a, b) => a.pR - b.pR);

        // Build HTML table
        const filterInfo = filterDesc ? `<p style="font-size: 11px; color: #333; margin-bottom: 8px; background: #f0f9ff; padding: 4px 8px; border-radius: 4px;"><b>Filter:</b> ${filterDesc}</p>` : '';
        let html = `
            <h4 style="margin-bottom: 8px;">Effect of Different Hotspot Mutations on ${gene1} vs ${gene2} Correlation</h4>
            ${filterInfo}
            <p style="font-size: 12px; color: #666; margin-bottom: 12px;">
                Comparing correlation between WT (0 mutations) vs Mutant (2+ mutations) for each hotspot mutation gene.
                Sorted by p-value (most significant effect first).
            </p>
            <div style="overflow-x: auto;">
            <table id="compareMutationsTable" class="data-table" style="width: 100%; font-size: 12px;">
                <thead>
                    <tr>
                        <th data-col="0" style="cursor: pointer;">Mutation Gene ▼</th>
                        <th data-col="1" style="cursor: pointer;">N (WT)</th>
                        <th data-col="2" style="cursor: pointer;">r (WT)</th>
                        <th data-col="3" style="cursor: pointer;">slope (WT)</th>
                        <th data-col="4" style="cursor: pointer;">N (Mut)</th>
                        <th data-col="5" style="cursor: pointer;">r (Mut)</th>
                        <th data-col="6" style="cursor: pointer;">slope (Mut)</th>
                        <th data-col="7" style="cursor: pointer;">Δr</th>
                        <th data-col="8" style="cursor: pointer;">p(Δr)</th>
                        <th data-col="9" style="cursor: pointer;">Δslope</th>
                    </tr>
                </thead>
                <tbody>
        `;

        tableData.forEach(row => {
            const deltaRColor = row.deltaR < 0 ? '#dc2626' : '#16a34a';
            const deltaSlopeColor = row.deltaSlope < 0 ? '#dc2626' : '#16a34a';
            const pHighlight = row.pR < 0.05 ? 'background: #fef3c7;' : '';

            html += `
                <tr style="${pHighlight}">
                    <td><b>${row.mutGene}</b></td>
                    <td>${row.nWT}</td>
                    <td>${row.rWT.toFixed(3)}</td>
                    <td>${row.slopeWT.toFixed(3)}</td>
                    <td>${row.nMut}</td>
                    <td>${row.rMut.toFixed(3)}</td>
                    <td>${row.slopeMut.toFixed(3)}</td>
                    <td style="color: ${deltaRColor}; font-weight: 600;">${row.deltaR.toFixed(3)}</td>
                    <td>${row.pR.toExponential(1)}</td>
                    <td style="color: ${deltaSlopeColor}; font-weight: 600;">${row.deltaSlope.toFixed(3)}</td>
                </tr>
            `;
        });

        html += `
                </tbody>
            </table>
            </div>
            <p style="font-size: 11px; color: #666; margin-top: 8px;">
                Note: Rows highlighted in yellow have p < 0.05. This analysis may be biased as specific mutations select for cancer types.
            </p>
            <div style="margin-top: 12px;">
                <button class="btn btn-primary btn-sm" id="backToGraphBtn" style="margin-right: 8px;">← Back to Graph</button>
                <button class="btn btn-success btn-sm" id="downloadMutCompareCSV">Download CSV</button>
            </div>
        `;

        document.getElementById('compareTable').innerHTML = html;

        // Add back to graph handler
        document.getElementById('backToGraphBtn')?.addEventListener('click', () => {
            document.getElementById('compareTable').style.display = 'none';
            document.getElementById('scatterPlot').style.display = 'block';
        });

        // Add download handler
        document.getElementById('downloadMutCompareCSV')?.addEventListener('click', () => {
            let csv = 'Mutation Gene,N (WT),r (WT),slope (WT),N (Mut),r (Mut),slope (Mut),Δr,p(Δr),Δslope\n';
            tableData.forEach(row => {
                csv += `${row.mutGene},${row.nWT},${row.rWT.toFixed(4)},${row.slopeWT.toFixed(4)},${row.nMut},${row.rMut.toFixed(4)},${row.slopeMut.toFixed(4)},${row.deltaR.toFixed(4)},${row.pR.toExponential(2)},${row.deltaSlope.toFixed(4)}\n`;
            });
            this.downloadFile(csv, `${gene1}_vs_${gene2}_all_mutations_comparison.csv`, 'text/csv');
        });

        // Make table sortable
        this.setupSortableTable('compareMutationsTable');
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
                    slope: stats.slope,
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

        // Build statistics table HTML with sortable headers
        let tableHtml = `
            <h4 style="margin: 0 0 10px 0;">Statistics by Lineage</h4>
            <div style="max-height: 500px; overflow-y: auto;">
            <table id="byTissueTable" style="width: 100%; border-collapse: collapse; font-size: 11px;">
                <thead>
                    <tr style="background-color: #22c55e; color: white;">
                        <th data-col="0" style="padding: 6px; border: 1px solid #16a34a; text-align: left; position: sticky; top: 0; background-color: #22c55e; cursor: pointer;">Lineage ▼</th>
                        <th data-col="1" style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e; cursor: pointer;">N</th>
                        <th data-col="2" style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e; cursor: pointer;">Corr</th>
                        <th data-col="3" style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e; cursor: pointer;">Slope</th>
                        <th data-col="4" style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e; cursor: pointer;">${gene1} (mean)</th>
                        <th data-col="5" style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e; cursor: pointer;">${gene1} (SD)</th>
                        <th data-col="6" style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e; cursor: pointer;">${gene2} (mean)</th>
                        <th data-col="7" style="padding: 6px; border: 1px solid #16a34a; text-align: center; position: sticky; top: 0; background-color: #22c55e; cursor: pointer;">${gene2} (SD)</th>
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
                    <td style="padding: 5px; border: 1px solid #ddd; text-align: center;">${t.slope.toFixed(3)}</td>
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

        // Add sortable table headers
        this.setupSortableTable('byTissueTable');
    }

    setupSortableTable(tableId) {
        const table = document.getElementById(tableId);
        if (!table) return;

        const headers = table.querySelectorAll('th[data-col]');
        headers.forEach(th => {
            th.addEventListener('click', () => {
                const col = parseInt(th.dataset.col);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                const isAsc = th.dataset.dir !== 'asc';
                th.dataset.dir = isAsc ? 'asc' : 'desc';

                // Update header arrows
                headers.forEach(h => {
                    h.textContent = h.textContent.replace(/ [▲▼]$/, '');
                });
                th.textContent += isAsc ? ' ▲' : ' ▼';

                rows.sort((a, b) => {
                    const aVal = a.children[col]?.textContent || '';
                    const bVal = b.children[col]?.textContent || '';

                    // Try numeric sort first
                    const aNum = parseFloat(aVal);
                    const bNum = parseFloat(bVal);
                    if (!isNaN(aNum) && !isNaN(bNum)) {
                        return isAsc ? aNum - bNum : bNum - aNum;
                    }
                    return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                });

                rows.forEach(row => tbody.appendChild(row));
            });
        });
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
