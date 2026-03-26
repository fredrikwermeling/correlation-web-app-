/**
 * Gene Correlation Explorer - Web Edition
 * JavaScript implementation matching the R Shiny app functionality
 *
 * Author: Wermeling Lab / Karolinska Institutet
 * Based on: https://github.com/fredrikwermeling/correlation-app
 */

// Override Plotly default font: use Arial (universally available in Inkscape/Illustrator)
(function() {
    const _origNewPlot = Plotly.newPlot;
    Plotly.newPlot = function(div, data, layout, config) {
        layout = layout || {};
        if (!layout.font) layout.font = {};
        if (!layout.font.family) layout.font.family = 'Arial, Helvetica, sans-serif';
        return _origNewPlot.call(this, div, data, layout, config);
    };
})();

class CorrelationExplorer {
    static CATEGORY_COLORS = [
        '#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4',
        '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990',
        '#dcbeff', '#9A6324', '#800000', '#aaffc3', '#808000',
        '#000075', '#a9a9a9', '#e6beff', '#ffe119', '#ffd8b1'
    ];

    static PRIORITY_FUSION_GENES = new Set([
        'BCR','ABL1','ALK','EML4','EWSR1','FLI1','MYC','KMT2A','PML','RARA',
        'RET','ROS1','NTRK1','NTRK2','NTRK3','ETV6','RUNX1','BRAF','FGFR1',
        'FGFR2','FGFR3','TMPRSS2','ERG','SS18','PAX3','PAX7','NUTM1','RAF1',
        'MET','NRG1','PDGFRA','PDGFRB'
    ]);

    constructor() {
        // Cell Line Browser state
        this._clbSelectedCellLines = new Set();
        // CLB sort direction (true = ascending, false = descending)
        this._clbSortAsc = true;

        // Data storage
        this.metadata = null;
        this.cellLineMetadata = null;
        this.mutations = null;
        this.translocations = null;
        this.orthologs = null;
        this.geneEffects = null; // Float32Array [nGenes x nCellLines]
        this.nGenes = 0;
        this.nCellLines = 0;
        this.geneIndex = new Map(); // gene name -> row index
        this.geneNames = []; // gene names array

        // Analysis results
        this.results = null;
        this.network = null;
        this.savedNetworkView = null;

        // Current inspect state
        this.currentInspect = null;
        this.clickedCells = new Set();
        this._userLegendPosition = null;
        this._userTitlePosition = null;
        this._savedScatterTextSettings = null;
        this._userLabelPositions = new Map();

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

        // Select mode state
        this.selectMode = false;
        this.selectedNodes = new Set();

        this.init();
    }

    // Helper to read numeric value from input, handling locale comma as decimal separator
    getInputNum(id) {
        const val = document.getElementById(id).value.replace(',', '.');
        return parseFloat(val);
    }

    // Helper to format numbers with proper minus sign (− instead of -)
    formatNum(value, decimals = 3) {
        if (value === null || value === undefined || isNaN(value)) return '-';
        const formatted = value.toFixed(decimals);
        // Replace hyphen-minus with proper minus sign for negative numbers
        return formatted.replace(/^-/, '−');
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

        // Load essential JSON files in parallel (synonyms loaded lazily on demand)
        const [metadataRes, cellLineRes, mutationsRes, orthologsRes, translocationsRes] = await Promise.all([
            fetch('web_data/metadata.json'),
            fetch('web_data/cellLineMetadata.json'),
            fetch('web_data/mutations.json'),
            fetch('web_data/orthologs.json'),
            fetch('web_data/translocations.json').catch(() => null)
        ]);

        this.metadata = await metadataRes.json();
        this.cellLineMetadata = await cellLineRes.json();
        this.mutations = await mutationsRes.json();
        this.orthologs = await orthologsRes.json();
        if (translocationsRes && translocationsRes.ok) {
            this.translocations = await translocationsRes.json();
            // Pre-compute fusion gene counts for fast dropdown population
            if (this.translocations?.genes) {
                const PRIO = CorrelationExplorer.PRIORITY_FUSION_GENES;
                this._fusionGeneCounts = this.translocations.genes
                    .map(g => {
                        const td = this.translocations.geneData[g];
                        let nFused = 0;
                        if (td) { for (const v of Object.values(td.translocations)) { if (v >= 1) nFused++; } }
                        return { gene: g, nFused };
                    })
                    .filter(x => x.nFused > 0)
                    .sort((a, b) => {
                        const aPri = PRIO.has(a.gene) ? 1 : 0;
                        const bPri = PRIO.has(b.gene) ? 1 : 0;
                        if (aPri !== bPri) return bPri - aPri;
                        return b.nFused - a.nFused;
                    });
            }
        }
        // synonymLookup loaded lazily when "Find Synonyms" is clicked

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

        // Populate tissue exclude list
        this.excludedTissues = new Set();
        this.populateTissueExcludeList();
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
                const subLineage = this.cellLineMetadata.primaryDisease?.[cellLine] || '';

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

            Object.keys(lineageCounts).sort((a, b) => lineageCounts[b] - lineageCounts[a]).forEach(lineage => {
                const option = document.createElement('option');
                option.value = lineage;
                option.textContent = `${lineage} (n=${lineageCounts[lineage]})`;
                select.appendChild(option);
            });
            document.getElementById('lineageFilterGroup').style.display = 'block';

            // Update sub-lineage and cascade all dependent selectors when lineage changes
            select.addEventListener('change', () => {
                this.updateSubLineageFilter();
                this._refreshFilteredSelectors();
            });
        }

        // Also populate parameter hotspot filter
        this.populateParamHotspotFilter();

        // Populate parameter translocation filter
        this.populateParamTranslocationFilter();
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

            Object.keys(subLineages).sort((a, b) => subLineages[b] - subLineages[a]).forEach(sub => {
                const option = document.createElement('option');
                option.value = sub;
                option.textContent = `${sub} (n=${subLineages[sub]})`;
                subSelect.appendChild(option);
            });
            document.getElementById('subLineageFilterGroup').style.display = 'block';

            // Add listener for sub-lineage changes (only add once)
            if (!subSelect.hasAttribute('data-listener-attached')) {
                subSelect.addEventListener('change', () => {
                    this._refreshFilteredSelectors();
                });
                subSelect.setAttribute('data-listener-attached', 'true');
            }
        } else {
            document.getElementById('subLineageFilterGroup').style.display = 'none';
        }

        // Update hotspot counts for selected lineage
        this.updateHotspotCountsForCurrentFilters();
    }

    updateScatterSubtypeFilter() {
        const lineage = document.getElementById('scatterCancerFilter').value;
        const subSelect = document.getElementById('scatterSubtypeFilter');

        if (!lineage || !this.currentInspect?.data) {
            subSelect.style.display = 'none';
            subSelect.innerHTML = '<option value="">All subtypes</option>';
            return;
        }

        // Find subtypes for this lineage from current inspect data
        const subtypeCounts = {};
        this.currentInspect.data.forEach(d => {
            if (d.lineage === lineage) {
                const subtype = this.cellLineMetadata?.primaryDisease?.[d.cellLineId] || '';
                if (subtype) {
                    subtypeCounts[subtype] = (subtypeCounts[subtype] || 0) + 1;
                }
            }
        });

        const subtypes = Object.keys(subtypeCounts).sort();
        if (subtypes.length > 1) {
            const lineageCount = this.currentInspect.data.filter(d => d.lineage === lineage).length;
            subSelect.innerHTML = `<option value="">All subtypes (n=${lineageCount})</option>`;
            subtypes.forEach(sub => {
                subSelect.innerHTML += `<option value="${sub}">${sub} (n=${subtypeCounts[sub]})</option>`;
            });
            subSelect.style.display = 'block';
        } else {
            subSelect.style.display = 'none';
            subSelect.innerHTML = '<option value="">All subtypes</option>';
        }
    }

    _styleActiveFilters() {
        const ids = ['scatterCancerFilter', 'scatterSubtypeFilter', 'colorByCategory'];
        ids.forEach(id => {
            const el = document.getElementById(id);
            if (!el) return;
            if (el.value) {
                el.style.borderColor = '#3b82f6';
                el.style.background = '#eff6ff';
                el.style.fontWeight = '600';
            } else {
                el.style.borderColor = '';
                el.style.background = '';
                el.style.fontWeight = '';
            }
        });
    }

    updateScatterHotspotFilterCounts() {
        if (!this.currentInspect?.data) return;

        const cancerFilter = document.getElementById('scatterCancerFilter').value;
        const subtypeFilter = document.getElementById('scatterSubtypeFilter').value;

        // Get the cell lines visible in the scatter plot (respecting tissue/subtype filter)
        let filteredData = this.currentInspect.data;
        if (cancerFilter) {
            filteredData = filteredData.filter(d => d.lineage === cancerFilter);
        }
        if (subtypeFilter && this.cellLineMetadata?.primaryDisease) {
            filteredData = filteredData.filter(d =>
                this.cellLineMetadata.primaryDisease[d.cellLineId] === subtypeFilter);
        }
        const filteredCellLines = new Set(filteredData.map(d => d.cellLineId));

        // Update hotspot gene selector and mutation filter gene selector
        const hotspotSelect = document.getElementById('hotspotGene');
        const mutFilterGeneSelect = document.getElementById('mutationFilterGene');

        if (hotspotSelect && this.mutations?.genes?.length > 0) {
            const hotspotVal = hotspotSelect.value;
            const mutFilterVal = mutFilterGeneSelect?.value || '';
            hotspotSelect.innerHTML = '<option value="">Select gene...</option>';
            if (mutFilterGeneSelect) mutFilterGeneSelect.innerHTML = '<option value="">No filter</option>';

            this.mutations.genes.forEach(g => {
                const mutData = this.mutations.geneData?.[g]?.mutations || {};
                let count = 0;
                filteredCellLines.forEach(cl => { if (mutData[cl] > 0) count++; });
                hotspotSelect.innerHTML += `<option value="${g}"${g === hotspotVal ? ' selected' : ''}>${g} (${count} mut)</option>`;
                if (mutFilterGeneSelect) {
                    mutFilterGeneSelect.innerHTML += `<option value="${g}"${g === mutFilterVal ? ' selected' : ''}>${g} (${count} mut)</option>`;
                }
            });
        }

        // Update translocation/fusion datalists
        const transGeneDatalist = document.getElementById('translocationGeneList');
        const transFilterGeneDatalist = document.getElementById('translocationFilterGeneList');

        if (transGeneDatalist && this.translocations?.genes?.length > 0) {
            const geneCounts = [];
            for (const g of this.translocations.genes) {
                const transData = this.translocations.geneData?.[g]?.translocations || {};
                let count = 0;
                for (const cl of filteredCellLines) {
                    if (transData[cl] && transData[cl] > 0) count++;
                }
                if (count > 0) geneCounts.push({ gene: g, count });
            }
            geneCounts.sort((a, b) => {
                const aPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(a.gene) ? 1 : 0;
                const bPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(b.gene) ? 1 : 0;
                if (aPri !== bPri) return bPri - aPri;
                return b.count - a.count;
            });

            let transHtml = '';
            geneCounts.forEach(({ gene, count }) => {
                transHtml += `<option value="${gene}">${gene} (${count} fused)</option>`;
            });
            transGeneDatalist.innerHTML = transHtml;
            if (transFilterGeneDatalist) transFilterGeneDatalist.innerHTML = transHtml;
        }
    }

    adjustAxis(id, delta) {
        const el = document.getElementById(id);
        const current = parseFloat(el.value.replace(',', '.'));
        if (isNaN(current)) {
            el.value = delta.toFixed(1);
        } else {
            el.value = (current + delta).toFixed(1);
        }
        this.updateInspectPlot();
    }

    updateGeSubtypeFilter() {
        const lineage = document.getElementById('geTissueFilter')?.value;
        const subSelect = document.getElementById('geSubtypeFilter');
        if (!subSelect) return;

        if (!lineage || !this.cellLineMetadata?.primaryDisease) {
            subSelect.style.display = 'none';
            subSelect.innerHTML = '<option value="">All subtypes</option>';
            return;
        }

        // Gather primaryDisease counts from ALL cell lines in this lineage
        const cellLines = this.metadata?.cellLines || [];
        const subtypeCounts = {};
        let lineageTotal = 0;
        cellLines.forEach(cl => {
            if ((this.cellLineMetadata?.lineage?.[cl] || '') !== lineage) return;
            lineageTotal++;
            const subtype = this.cellLineMetadata.primaryDisease[cl] || '';
            if (subtype) {
                subtypeCounts[subtype] = (subtypeCounts[subtype] || 0) + 1;
            }
        });

        const subtypes = Object.keys(subtypeCounts).sort();
        if (subtypes.length > 1) {
            subSelect.innerHTML = `<option value="">All subtypes (n=${lineageTotal})</option>`;
            subtypes.forEach(sub => {
                subSelect.innerHTML += `<option value="${sub}">${sub} (n=${subtypeCounts[sub]})</option>`;
            });
            subSelect.style.display = '';
        } else {
            subSelect.style.display = 'none';
            subSelect.innerHTML = '<option value="">All subtypes</option>';
        }
    }

    updateHotspotCountsForCurrentFilters() {
        if (this.mutations?.geneData) {
            this.updateParamHotspotGeneCounts();

            // Also update mutation mode selector if in mutation mode
            const isMutationMode = document.querySelector('input[name="analysisMode"]:checked')?.value === 'mutation';
            if (isMutationMode) {
                this.populateMutationHotspotSelector();
                this.populateTranslocationHotspotSelector();
            }
        }
        if (this.translocations?.geneData) {
            this.updateParamTranslocationGeneCounts();
        }
    }

    _refreshFilteredSelectors() {
        // Called when lineage/sublineage changes — refresh all dependent selectors
        this.updateHotspotCountsForCurrentFilters();
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
        const select = document.getElementById('paramHotspotGene');
        const cellLines = this.metadata.cellLines;
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;
        const hasExcluded = this.excludedTissues && this.excludedTissues.size > 0;
        const currentValue = select.value;

        const filteredCLs = (lineageFilter || subLineageFilter || hasExcluded)
            ? cellLines.filter(cl => {
                if (lineageFilter && this.cellLineMetadata?.lineage?.[cl] !== lineageFilter) return false;
                if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cl] !== subLineageFilter) return false;
                if (hasExcluded) {
                    const lin = this.cellLineMetadata?.lineage?.[cl];
                    if (lin && this.excludedTissues.has(lin)) return false;
                }
                return true;
            }) : cellLines;

        const geneCounts = [];
        for (const gene of Object.keys(this.mutations.geneData)) {
            const mutations = this.mutations.geneData[gene].mutations;
            let nMut = 0;
            for (const cl of filteredCLs) {
                if (mutations[cl] && mutations[cl] > 0) nMut++;
            }
            if (nMut > 0) geneCounts.push({ gene, count: nMut });
        }
        geneCounts.sort((a, b) => b.count - a.count);

        select.innerHTML = '<option value="">No filter</option>';
        for (const { gene, count } of geneCounts) {
            const option = document.createElement('option');
            option.value = gene;
            option.textContent = `${gene} (n=${count} mutated)`;
            select.appendChild(option);
        }

        if (currentValue) select.value = currentValue;

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
            if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== subLineageFilter) {
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

    populateParamTranslocationFilter() {
        if (this.translocations && this.translocations.geneData) {
            const input = document.getElementById('paramTranslocationGene');
            document.getElementById('paramTranslocationFilterGroup').style.display = 'block';

            input.addEventListener('input', () => {
                // Only act on valid gene names or empty (clear filter)
                const val = input.value.trim();
                if (val === '' || this.translocations.geneData[val]) {
                    this.updateParamTranslocationLevelCounts();
                }
            });

            this.updateParamTranslocationGeneCounts();
        }
    }

    updateParamTranslocationGeneCounts() {
        if (!this.translocations?.geneData) return;
        const genes = Object.keys(this.translocations.geneData);
        const input = document.getElementById('paramTranslocationGene');
        const datalist = document.getElementById('paramTranslocationGeneList');
        const cellLines = this.metadata.cellLines;
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;
        const hasExcluded = this.excludedTissues && this.excludedTissues.size > 0;
        const currentValue = input.value;

        // Build filtered cell line list once
        const filteredCLs = (lineageFilter || subLineageFilter || hasExcluded)
            ? cellLines.filter(cl => {
                if (lineageFilter && this.cellLineMetadata?.lineage?.[cl] !== lineageFilter) return false;
                if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cl] !== subLineageFilter) return false;
                if (hasExcluded) {
                    const lin = this.cellLineMetadata?.lineage?.[cl];
                    if (lin && this.excludedTissues.has(lin)) return false;
                }
                return true;
            }) : cellLines;

        // Count fusions per gene, skip genes with 0 fusions
        const geneCounts = [];
        for (const gene of genes) {
            const translocations = this.translocations.geneData[gene].translocations;
            let nFused = 0;
            for (const cl of filteredCLs) {
                if (translocations[cl] && translocations[cl] > 0) nFused++;
            }
            if (nFused > 0) geneCounts.push({ gene, nFused });
        }
        geneCounts.sort((a, b) => {
            const aPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(a.gene) ? 1 : 0;
            const bPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(b.gene) ? 1 : 0;
            if (aPri !== bPri) return bPri - aPri;
            return b.nFused - a.nFused;
        });

        let html = '';
        geneCounts.forEach(({ gene, nFused }) => {
            html += `<option value="${gene}">${gene} (n=${nFused} fused)</option>`;
        });
        datalist.innerHTML = html;

        if (currentValue) input.value = currentValue;

        this.updateParamTranslocationLevelCounts();
    }

    updateParamTranslocationLevelCounts() {
        const gene = document.getElementById('paramTranslocationGene').value;
        const levelSelect = document.getElementById('paramTranslocationLevel');
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;

        if (!gene || !this.translocations?.geneData?.[gene]) {
            levelSelect.innerHTML = `
                <option value="all">All cells</option>
                <option value="0">Only WT (no fusion)</option>
                <option value="1">Only 1 fusion partner</option>
                <option value="2">Only 2+ fusion partners</option>
                <option value="1+2">Only fused (1+2)</option>
            `;
            return;
        }

        const translocations = this.translocations.geneData[gene].translocations;
        const cellLines = this.metadata.cellLines;
        let n0 = 0, n1 = 0, n2 = 0;

        cellLines.forEach(cellLine => {
            if (lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== lineageFilter) return;
            if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== subLineageFilter) return;

            const level = translocations[cellLine] || 0;
            if (level === 0) n0++;
            else if (level === 1) n1++;
            else n2++;
        });

        const nFused = n1 + n2;
        const total = n0 + n1 + n2;

        levelSelect.innerHTML = `
            <option value="all">All cells (n=${total})</option>
            <option value="0">Only WT (n=${n0})</option>
            <option value="1">Only 1 partner (n=${n1})</option>
            <option value="2">Only 2+ partners (n=${n2})</option>
            <option value="1+2">Only fused 1+2 (n=${nFused})</option>
        `;
    }

    updateAnalysisModeUI() {
        const mode = document.querySelector('input[name="analysisMode"]:checked').value;
        const isMutationMode = mode === 'mutation';
        const isDesignMode = mode === 'design';
        const isSynonymMode = mode === 'synonym';

        // Toggle visibility of correlation/slope params (hide for mutation and synonym modes)
        const hideParams = isMutationMode || isSynonymMode;
        document.getElementById('correlationParams').style.display = hideParams ? 'none' : 'block';
        document.getElementById('slopeParams').style.display = hideParams ? 'none' : 'block';

        // Hide min cell lines and filters for synonym mode
        document.getElementById('minCellLinesGroup').style.display = isSynonymMode ? 'none' : 'block';
        if (isSynonymMode) {
            document.getElementById('lineageFilterGroup').style.display = 'none';
            document.getElementById('subLineageFilterGroup').style.display = 'none';
            document.getElementById('paramHotspotFilterGroup').style.display = 'none';
            document.getElementById('paramTranslocationFilterGroup').style.display = 'none';
        } else {
            // Restore filters if data is loaded
            if (this.cellLineMetadata && this.cellLineMetadata.lineage) {
                document.getElementById('lineageFilterGroup').style.display = 'block';
                // Ensure subtype filter is shown if a lineage is already selected
                this.updateSubLineageFilter();
            }
            if (this.mutations && this.mutations.geneData) {
                document.getElementById('paramHotspotFilterGroup').style.display = 'block';
            }
            if (this.translocations && this.translocations.geneData) {
                document.getElementById('paramTranslocationFilterGroup').style.display = 'block';
            }
        }

        // Toggle visibility of mutation-specific params
        document.getElementById('mutationHotspotGroup').style.display = isMutationMode ? 'block' : 'none';
        document.getElementById('pValueThresholdGroup').style.display = isMutationMode ? 'block' : 'none';

        // Show/hide design mode hint
        document.getElementById('designModeHint').style.display = isDesignMode ? 'block' : 'none';

        // Show/hide design expand option
        document.getElementById('designExpandOption').style.display = isDesignMode ? 'block' : 'none';

        // Show/hide mutation tab
        document.getElementById('mutationTab').style.display = isMutationMode ? 'inline-block' : 'none';

        // Show/hide synonyms tab
        document.getElementById('synonymsTab').style.display = isSynonymMode ? 'inline-block' : 'none';

        // Disable/enable gene input elements for mutation mode (not needed, but keep Run button active)
        const geneInputElements = document.querySelectorAll('#geneTextarea, #manualStatsTextarea, #statsFileInput, .input-tab, .stats-sub-tab, #loadTestGenes, #clearGenes, #loadManualStatsBtn, #loadTestStats, #downloadSampleStats');
        geneInputElements.forEach(el => {
            if (el) {
                el.disabled = isMutationMode;
                el.style.opacity = isMutationMode ? '0.5' : '1';
                el.style.pointerEvents = isMutationMode ? 'none' : 'auto';
            }
        });
        // Grey out the input panels but not the Run button
        const inputPanels = document.querySelectorAll('.input-panel, .input-tabs');
        inputPanels.forEach(el => {
            if (el) {
                el.style.opacity = isMutationMode ? '0.5' : '1';
            }
        });

        // Set default min cell lines based on mode
        const minCellLinesInput = document.getElementById('minCellLines');
        if (isMutationMode) {
            minCellLinesInput.value = '5';
            this.populateMutationHotspotSelector();
            this.populateTranslocationHotspotSelector();
            this.updateMutAnalysisTypeUI();
        } else {
            minCellLinesInput.value = '25';
        }
    }

    updateMutAnalysisTypeUI() {
        const subType = document.querySelector('input[name="mutAnalysisType"]:checked')?.value || 'hotspot';
        const isTranslocation = subType === 'translocation';
        document.getElementById('hotspotAnalysisControls').style.display = isTranslocation ? 'none' : '';
        document.getElementById('translocationAnalysisControls').style.display = isTranslocation ? '' : 'none';
        // Hide translocation radio if no data
        if (!this.translocations?.geneData) {
            document.getElementById('mutAnalysisTypeSelector').style.display = 'none';
        }
    }

    populateMutationHotspotSelector() {
        const select = document.getElementById('mutationHotspotSelect');
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;
        const hasExcluded = this.excludedTissues && this.excludedTissues.size > 0;
        const currentValue = select.value;

        if (!this.mutations || !this.mutations.geneData) return;

        const genes = Object.keys(this.mutations.geneData);
        const cellLines = this.metadata.cellLines;

        // Build filtered cell line list
        const filteredCLs = (lineageFilter || subLineageFilter || hasExcluded)
            ? cellLines.filter(cl => {
                if (lineageFilter && this.cellLineMetadata?.lineage?.[cl] !== lineageFilter) return false;
                if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cl] !== subLineageFilter) return false;
                if (hasExcluded) {
                    const lin = this.cellLineMetadata?.lineage?.[cl];
                    if (lin && this.excludedTissues.has(lin)) return false;
                }
                return true;
            }) : cellLines;

        // Count mutations per gene, skip genes with 0
        const geneCounts = [];
        for (const gene of genes) {
            const mutations = this.mutations.geneData[gene].mutations;
            let nMut = 0;
            for (const cl of filteredCLs) {
                if (mutations[cl] && mutations[cl] > 0) nMut++;
            }
            if (nMut > 0) geneCounts.push({ gene, count: nMut });
        }
        geneCounts.sort((a, b) => b.count - a.count);

        select.innerHTML = '<option value="">Select hotspot gene...</option>';
        geneCounts.forEach(({ gene, count }) => {
            const option = document.createElement('option');
            option.value = gene;
            option.textContent = `${gene} (${count} mutated)`;
            select.appendChild(option);
        });

        if (currentValue) select.value = currentValue;
    }

    populateTranslocationHotspotSelector() {
        const input = document.getElementById('translocationHotspotSelect');
        if (!input) return;
        if (!this.translocations || !this.translocations.geneData) return;

        const datalist = document.getElementById('translocationHotspotList');
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;
        const hasExcluded = this.excludedTissues && this.excludedTissues.size > 0;
        const currentValue = input.value;

        const genes = Object.keys(this.translocations.geneData);
        const cellLines = this.metadata.cellLines;

        // Build filtered cell line list once
        const filteredCLs = (lineageFilter || subLineageFilter || hasExcluded)
            ? cellLines.filter(cl => {
                if (lineageFilter && this.cellLineMetadata?.lineage?.[cl] !== lineageFilter) return false;
                if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cl] !== subLineageFilter) return false;
                if (hasExcluded) {
                    const lin = this.cellLineMetadata?.lineage?.[cl];
                    if (lin && this.excludedTissues.has(lin)) return false;
                }
                return true;
            }) : cellLines;

        // Count fusions per gene, skip genes with 0 fusions
        const geneCounts = [];
        for (const gene of genes) {
            const translocations = this.translocations.geneData[gene].translocations;
            let nFused = 0;
            for (const cl of filteredCLs) {
                if (translocations[cl] && translocations[cl] > 0) nFused++;
            }
            if (nFused > 0) geneCounts.push({ gene, nFused });
        }
        geneCounts.sort((a, b) => {
            const aPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(a.gene) ? 1 : 0;
            const bPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(b.gene) ? 1 : 0;
            if (aPri !== bPri) return bPri - aPri;
            return b.nFused - a.nFused;
        });

        let html = '';
        geneCounts.forEach(({ gene, nFused }) => {
            html += `<option value="${gene}">${gene} (${nFused} fused cells)</option>`;
        });
        datalist.innerHTML = html;

        if (currentValue) input.value = currentValue;
    }

    _showUpsetPlot() {
        const data = this._oncoprintData;
        if (!data || !data.topGenes) return;

        const activeFilters = Object.entries(this._oncoprintFilters || {}).filter(([, v]) => v !== 'none');
        const hotspotGene = document.getElementById('mutationHotspotSelect')?.value;
        let upsetGenes, upsetLabel;
        if (activeFilters.length >= 2) {
            upsetGenes = activeFilters.map(([gene]) => data.topGenes.find(g => g.gene === gene)).filter(Boolean);
            if (hotspotGene && !upsetGenes.find(g => g.gene === hotspotGene)) {
                const hg = data.topGenes.find(g => g.gene === hotspotGene);
                if (hg) upsetGenes.unshift(hg);
            }
            upsetLabel = upsetGenes.map(g => g.gene).join(', ');
        } else {
            const used = new Set();
            upsetGenes = [];
            if (hotspotGene) {
                const hg = data.topGenes.find(g => g.gene === hotspotGene);
                if (hg) { upsetGenes.push(hg); used.add(hotspotGene); }
            }
            for (const g of data.topGenes) {
                if (upsetGenes.length >= 5) break;
                if (!used.has(g.gene)) { upsetGenes.push(g); used.add(g.gene); }
            }
            upsetLabel = hotspotGene ? `${hotspotGene} + Top ${upsetGenes.length - 1}` : `Top ${upsetGenes.length} most mutated`;
        }
        if (upsetGenes.length < 2) {
            alert('Select at least 2 genes in the oncoprint (include or exclude) to generate an UpSet plot, or use top 5 by default.');
            return;
        }

        const cls = data.allFilteredCLs || data.sortedCLs;
        const intersections = new Map();
        for (const cl of cls) {
            let key = '';
            for (const g of upsetGenes) {
                key += (g.muts[cl] > 0) ? '1' : '0';
            }
            if (!intersections.has(key)) intersections.set(key, { count: 0, key });
            intersections.get(key).count++;
        }

        const sorted = [...intersections.values()].sort((a, b) => b.count - a.count);
        const nGenes = upsetGenes.length;
        const nBars = sorted.length;
        const barW = 30;
        const matrixH = nGenes * 20;
        const plotW = Math.max(400, nBars * barW + 100);
        const plotH = 250 + matrixH;

        document.getElementById('upsetPopup')?.remove();

        const popup = document.createElement('div');
        popup.id = 'upsetPopup';
        popup.style.cssText = `position:fixed; z-index:10001; background:white; border:1px solid #d1d5db; border-radius:8px; box-shadow:0 8px 24px rgba(0,0,0,0.15); display:flex; flex-direction:column; max-width:90vw; max-height:85vh;`;
        popup.style.left = '50px';
        popup.style.top = '50px';

        const lineageF = document.getElementById('lineageFilter')?.value || '';
        const subLineageF = document.getElementById('subLineageFilter')?.value || '';
        const hasExcl = this.excludedTissues && this.excludedTissues.size > 0;
        let filterCtx = '';
        if (lineageF) {
            filterCtx = lineageF;
            if (subLineageF) filterCtx += ` · ${subLineageF}`;
        } else if (hasExcl) {
            const allLin = this.cellLineMetadata?.lineage ? [...new Set(Object.values(this.cellLineMetadata.lineage))] : [];
            const incl = allLin.filter(t => !this.excludedTissues.has(t));
            filterCtx = incl.length <= 4 ? incl.join(', ') : `${incl.length} tissues`;
        }

        let html = `<div id="upsetDragHandle" style="display:flex; justify-content:space-between; align-items:center; padding:6px 10px; background:#f0fdf4; border-radius:8px 8px 0 0; cursor:move; user-select:none;">`;
        html += `<span style="font-weight:600; font-size:12px;">UpSet — ${upsetLabel}</span>`;
        html += `<span style="font-size:10px; color:#6b7280;">${cls.length} cell lines${filterCtx ? ' · ' + filterCtx : ''}</span>`;
        html += `<button onclick="document.getElementById('upsetPopup').remove()" style="background:none;border:none;font-size:16px;cursor:pointer;color:#999;">&times;</button>`;
        html += `</div>`;
        html += `<div style="padding:10px; overflow:auto; flex:1;">`;
        html += `<div id="upsetPlotDiv" style="width:${plotW}px; height:${plotH}px;"></div>`;
        html += `</div>`;
        html += `<div style="display:flex; gap:6px; padding:6px 10px; border-top:1px solid #e5e7eb; align-items:center; flex-wrap:wrap;">`;
        html += `<label style="font-size:10px;cursor:pointer;color:#374151;"><input type="checkbox" id="upsetShowCounts" ${this._upsetShowCounts ? 'checked' : ''} onchange="app._upsetToggle('counts')" style="margin:0 3px 0 0;vertical-align:middle;"> Counts</label>`;
        html += `<label style="font-size:10px;cursor:pointer;color:#374151;"><input type="checkbox" id="upsetShowPct" ${this._upsetShowPct ? 'checked' : ''} onchange="app._upsetToggle('pct')" style="margin:0 3px 0 0;vertical-align:middle;"> %</label>`;
        html += `<label style="font-size:10px;cursor:pointer;color:#374151;"><input type="checkbox" id="upsetShowNames" ${this._upsetShowNames ? 'checked' : ''} onchange="app._upsetToggle('names')" style="margin:0 3px 0 0;vertical-align:middle;"> Names</label>`;
        html += `<span style="border-left:1px solid #d1d5db;height:14px;"></span>`;
        html += `<button onclick="app._upsetExport('svg')" style="font-size:10px;padding:2px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:#f9fafb;">SVG</button>`;
        html += `<button onclick="app._upsetExport('png')" style="font-size:10px;padding:2px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:#f9fafb;">PNG</button>`;
        html += `</div>`;
        popup.innerHTML = html;
        document.body.appendChild(popup);

        const dh = document.getElementById('upsetDragHandle');
        dh.addEventListener('mousedown', (e) => {
            const dx = e.clientX - popup.offsetLeft, dy = e.clientY - popup.offsetTop;
            const onMove = (e2) => { popup.style.left = (e2.clientX - dx) + 'px'; popup.style.top = (e2.clientY - dy) + 'px'; };
            const onUp = () => { document.removeEventListener('mousemove', onMove); document.removeEventListener('mouseup', onUp); };
            document.addEventListener('mousemove', onMove);
            document.addEventListener('mouseup', onUp);
        });

        this._upsetShowNames = this._upsetShowNames || false;
        this._upsetShowCounts = this._upsetShowCounts !== false;

        const barX = sorted.map((_, i) => i);
        const barY = sorted.map(s => s.count);
        const barLabels = sorted.map(s => {
            const bits = s.key.split('');
            const names = bits.map((b, i) => b === '1' ? upsetGenes[i].gene : null).filter(Boolean);
            return names.length === 0 ? 'None' : names.join(' \u2229 ');
        });

        const matrixTraces = [];
        for (let col = 0; col < nBars; col++) {
            const bits = sorted[col].key.split('');
            const activeRows = [];
            bits.forEach((b, row) => { if (b === '1') activeRows.push(row); });
            if (activeRows.length >= 2) {
                matrixTraces.push({
                    x: activeRows.map(() => col), y: activeRows.map(r => -(r + 1)),
                    mode: 'lines', line: { color: '#374151', width: 2 },
                    showlegend: false, hoverinfo: 'skip', xaxis: 'x', yaxis: 'y2'
                });
            }
        }
        const dotX = [], dotY = [], dotColor = [], dotText = [];
        for (let col = 0; col < nBars; col++) {
            const bits = sorted[col].key.split('');
            bits.forEach((b, row) => {
                dotX.push(col); dotY.push(-(row + 1));
                dotColor.push(b === '1' ? '#374151' : '#d1d5db');
                dotText.push(`${upsetGenes[row].gene}: ${b === '1' ? 'Mutated' : 'WT'}`);
            });
        }
        matrixTraces.push({
            x: dotX, y: dotY, mode: 'markers',
            marker: { color: dotColor, size: 10 },
            text: dotText, hoverinfo: 'text', showlegend: false, xaxis: 'x', yaxis: 'y2'
        });

        // Determine which bar matches the active oncoprint filter selection
        const oncoFilters = Object.entries(this._oncoprintFilters || {}).filter(([, v]) => v !== 'none');
        const barColors = sorted.map(s => {
            if (oncoFilters.length === 0) return '#3b82f6';
            const bits = s.key.split('');
            let matches = true;
            for (const [gene, state] of oncoFilters) {
                const idx = upsetGenes.findIndex(g => g.gene === gene);
                if (idx < 0) continue;
                const isMut = bits[idx] === '1';
                if (state === 'mut' && !isMut) { matches = false; break; }
                if (state === 'wt' && isMut) { matches = false; break; }
            }
            return matches ? '#16a34a' : '#cbd5e1';
        });

        const traces = [{
            x: barX, y: barY, type: 'bar',
            marker: { color: barColors, line: { color: barColors.map(c => c === '#16a34a' ? '#15803d' : 'transparent'), width: barColors.map(c => c === '#16a34a' ? 2 : 0) } },
            text: this._upsetShowNames ? barLabels : barLabels.map(() => ''),
            textposition: 'inside', textangle: -90, textfont: { size: 9, color: 'white' },
            customdata: barLabels,
            hovertemplate: '%{customdata}<br>%{y} cell lines<extra></extra>',
            showlegend: false
        }, ...matrixTraces];

        const totalCLs = cls.length;
        const showCounts = this._upsetShowCounts;
        const showPct = this._upsetShowPct;
        const countAnnotations = (showCounts || showPct) ? barX.map((x, i) => {
            const parts = [];
            if (showCounts) parts.push(String(barY[i]));
            if (showPct) parts.push(`${(barY[i] / totalCLs * 100).toFixed(1)}%`);
            return { x, y: barY[i], text: parts.join(' \u00b7 '), showarrow: false, font: { size: 9, color: '#374151' }, yanchor: 'bottom', yshift: 2 };
        }) : [];

        const layout = {
            font: { family: 'Arial, Helvetica, sans-serif' },
            annotations: countAnnotations,
            width: plotW, height: plotH,
            margin: { t: 20, b: 10, l: 60, r: 20 },
            xaxis: { showticklabels: false, showgrid: false, zeroline: false, range: [-0.5, nBars - 0.5] },
            yaxis: { title: 'Cell lines', domain: [0.4, 1], anchor: 'x' },
            yaxis2: {
                domain: [0, 0.35], anchor: 'x',
                tickvals: upsetGenes.map((_, i) => -(i + 1)),
                ticktext: upsetGenes.map(g => g.gene),
                showgrid: false, zeroline: false, range: [-(nGenes + 0.5), -0.5]
            },
            bargap: 0.3, plot_bgcolor: '#fafafa', paper_bgcolor: 'white'
        };

        Plotly.newPlot('upsetPlotDiv', traces, layout, { responsive: false, displayModeBar: false });

        const esc = (e) => { if (e.key === 'Escape') { e.stopImmediatePropagation(); popup.remove(); document.removeEventListener('keydown', esc); } };
        document.addEventListener('keydown', esc);
    }

    _upsetToggle(what) {
        if (what === 'counts') this._upsetShowCounts = !this._upsetShowCounts;
        if (what === 'pct') this._upsetShowPct = !this._upsetShowPct;
        if (what === 'names') this._upsetShowNames = !this._upsetShowNames;
        this._showUpsetPlot();
    }

    async _upsetExport(format) {
        const plotEl = document.getElementById('upsetPlotDiv');
        if (!plotEl?.data) return;
        const w = plotEl.layout?.width || 500;
        const h = plotEl.layout?.height || 400;
        const filename = 'upset_plot';

        if (format === 'svg') {
            const svgDataUrl = await Plotly.toImage(plotEl, { format: 'svg', width: w, height: h });
            let svgStr = svgDataUrl.indexOf('base64,') > -1 ? atob(svgDataUrl.split('base64,')[1]) : decodeURIComponent(svgDataUrl.split(',').slice(1).join(','));
            if (this._finalizeSvgForExport) svgStr = this._finalizeSvgForExport(svgStr);
            const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = `${filename}.svg`;
            a.click();
            URL.revokeObjectURL(a.href);
        } else {
            const url = await Plotly.toImage(plotEl, { format: 'png', width: w * 2, height: h * 2, scale: 2 });
            const a = document.createElement('a');
            a.href = url;
            a.download = `${filename}.png`;
            a.click();
        }
    }

    updateTBSelectionCount() {
        const popup = document.getElementById('tissueBreakdownPopup');
        if (!popup) return;
        const tissueCount = popup.querySelectorAll('.tb-check:checked').length;
        const subCount = popup.querySelectorAll('.tb-sub-check:checked').length;
        const label = document.getElementById('tbSelectionCount');
        const parts = [];
        if (tissueCount > 0) parts.push(`${tissueCount} tissue${tissueCount > 1 ? 's' : ''}`);
        if (subCount > 0) parts.push(`${subCount} subtype${subCount > 1 ? 's' : ''}`);
        if (label) label.textContent = parts.length === 0 ? '0 selected' : parts.join(', ');
    }

    showOncoprint(context) {
        document.getElementById('oncoprintPopup')?.remove();

        if (!this.mutations?.geneData) return;

        let filteredCLs;
        let filterLabel = '';
        if (context === 'clb' && this._clbVisibleCellLines) {
            filteredCLs = [...this._clbVisibleCellLines];
            const clbTissue = document.getElementById('clbTissueFilter')?.value;
            if (clbTissue) filterLabel = clbTissue;
        } else {
            const lineageFilter = document.getElementById('lineageFilter')?.value || '';
            const subLineageFilter = document.getElementById('subLineageFilter')?.value || '';
            const paramHotspot = document.getElementById('paramHotspotGene')?.value || '';
            const paramHotspotLevel = document.getElementById('paramHotspotLevel')?.value || 'all';
            const paramTrans = document.getElementById('paramTranslocationGene')?.value || '';
            const paramTransLevel = document.getElementById('paramTranslocationLevel')?.value || 'all';
            const paramHotspotMuts = paramHotspot && paramHotspotLevel !== 'all' ? this.mutations?.geneData?.[paramHotspot]?.mutations : null;
            const paramTransMuts = paramTrans && paramTransLevel !== 'all' ? this.translocations?.geneData?.[paramTrans]?.translocations : null;
            filteredCLs = this.metadata.cellLines.filter(cl => {
                if (lineageFilter && this.cellLineMetadata?.lineage?.[cl] !== lineageFilter) return false;
                if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cl] !== subLineageFilter) return false;
                if (this.excludedTissues && this.excludedTissues.size > 0) {
                    const lin = this.cellLineMetadata?.lineage?.[cl];
                    if (lin && this.excludedTissues.has(lin)) return false;
                }
                if (paramHotspotMuts) {
                    const level = paramHotspotMuts[cl] || 0;
                    if (paramHotspotLevel === '0' && level !== 0) return false;
                    if (paramHotspotLevel === '1+2' && level < 1) return false;
                    if (paramHotspotLevel === '1' && level !== 1) return false;
                    if (paramHotspotLevel === '2' && level < 2) return false;
                }
                if (paramTransMuts) {
                    const level = paramTransMuts[cl] || 0;
                    if (paramTransLevel === '0' && level !== 0) return false;
                    if (paramTransLevel === '1+2' && level < 1) return false;
                }
                return true;
            });
            const filterParts = [];
            if (lineageFilter) filterParts.push(lineageFilter);
            if (paramHotspot && paramHotspotLevel !== 'all') filterParts.push(`${paramHotspot} ${paramHotspotLevel === '0' ? 'WT' : 'Mut'}`);
            filterLabel = filterParts.join(' · ');
            if (!filterLabel && this.excludedTissues && this.excludedTissues.size > 0) filterLabel = 'filtered tissues';
        }

        const maxCLs = 200;
        const clsToShow = filteredCLs.slice(0, maxCLs);

        const maxGenes = 25;
        const geneCounts = [];
        for (const gene of Object.keys(this.mutations.geneData)) {
            const muts = this.mutations.geneData[gene].mutations;
            let n = 0;
            for (const cl of clsToShow) {
                if (muts[cl] && muts[cl] > 0) n++;
            }
            if (n > 0) geneCounts.push({ gene, n, muts });
        }
        geneCounts.sort((a, b) => b.n - a.n);
        const topGenes = geneCounts.slice(0, maxGenes);
        if (topGenes.length === 0) return;

        if (!this._oncoprintFilters) this._oncoprintFilters = {};

        const clMutCounts = new Map();
        for (const cl of clsToShow) {
            let count = 0;
            for (const { muts } of topGenes) {
                if (muts[cl] > 0) count++;
            }
            clMutCounts.set(cl, count);
        }
        const sortedCLs = [...clsToShow].sort((a, b) => clMutCounts.get(b) - clMutCounts.get(a));

        const cellW = Math.max(3, Math.min(8, Math.floor(500 / sortedCLs.length)));
        const cellH = 14;
        const labelW = 72;
        const boxW = 10;
        const boxGap = 1;
        const boxAreaW = boxW * 2 + boxGap + 4;
        const gridW = sortedCLs.length * cellW;
        const gridH = topGenes.length * cellH;
        const countW = 30;
        const totalW = boxAreaW + labelW + gridW + countW;
        const footerH = 50;
        const totalH = gridH + footerH;

        const popup = document.createElement('div');
        popup.id = 'oncoprintPopup';
        popup.style.cssText = `position:fixed; z-index:10000; background:white; border:1px solid #d1d5db; border-radius:8px; box-shadow:0 8px 24px rgba(0,0,0,0.15); display:flex; flex-direction:column; max-width:90vw; max-height:85vh;`;
        popup.style.right = '20px';
        popup.style.top = '20px';

        const currentHotspot = document.getElementById('mutationHotspotSelect').value;

        let html = `<div id="oncoprintDragHandle" style="display:flex; justify-content:space-between; align-items:center; padding:6px 10px; background:#f0fdf4; border-radius:8px 8px 0 0; cursor:move; user-select:none;">`;
        html += `<span style="font-weight:600; font-size:12px;">Oncoprint — Top ${topGenes.length} hotspot genes</span>`;
        html += `<span style="font-size:10px; color:#6b7280;">${sortedCLs.length} cell lines${filterLabel ? ' \u00b7 ' + filterLabel : ''}${filteredCLs.length > maxCLs ? ` (showing ${maxCLs} of ${filteredCLs.length})` : ''}</span>`;
        html += `<button onclick="document.getElementById('oncoprintPopup').remove()" style="background:none;border:none;font-size:16px;cursor:pointer;color:#999;">&times;</button>`;
        html += `</div>`;
        html += `<div style="padding:6px 10px; overflow:auto; flex:1;">`;
        html += `<div style="font-size:9px; color:#9ca3af; margin-bottom:4px;"><span style="color:#16a34a;">&#9632;</span> include \u00b7 <span style="color:#dc2626;">&#9632;</span> exclude \u00b7 <span style="display:inline-block;width:8px;height:8px;background:#3b82f6;vertical-align:middle;"></span> 1 mut \u00b7 <span style="display:inline-block;width:8px;height:8px;background:#1e40af;vertical-align:middle;"></span> 2 mut \u00b7 <span style="display:inline-block;width:8px;height:8px;background:#f3f4f6;border:1px solid #d1d5db;vertical-align:middle;"></span> WT</div>`;
        html += `<canvas id="oncoprintCanvas" width="${totalW}" height="${totalH}" style="cursor:pointer;"></canvas>`;
        html += `<div id="oncoprintStatus" style="font-size:10px; margin-top:4px; display:flex; gap:6px; align-items:center; flex-wrap:wrap;"></div>`;
        html += `</div>`;
        html += `<div style="display:flex; gap:4px; padding:6px 10px; border-top:1px solid #e5e7eb;">`;
        html += `<button onclick="app._oncoprintExport('svg')" style="font-size:10px;padding:2px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:#f9fafb;">SVG</button>`;
        html += `<button onclick="app._oncoprintExport('png')" style="font-size:10px;padding:2px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:#f9fafb;">PNG</button>`;
        html += `<button onclick="app._oncoprintExport('csv')" style="font-size:10px;padding:2px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:#f9fafb;">CSV</button>`;
        html += `<span style="border-left:1px solid #d1d5db;height:16px;margin:0 2px;"></span>`;
        html += `<button onclick="app._showUpsetPlot()" style="font-size:10px;padding:2px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:#f0fdf4;color:#16a34a;font-weight:500;">UpSet</button>`;
        html += `</div>`;
        popup.innerHTML = html;
        document.body.appendChild(popup);

        const dragHandle = document.getElementById('oncoprintDragHandle');
        let dragX, dragY;
        dragHandle.addEventListener('mousedown', (e) => {
            const rect = popup.getBoundingClientRect();
            popup.style.left = rect.left + 'px';
            popup.style.top = rect.top + 'px';
            popup.style.right = 'auto';
            dragX = e.clientX - rect.left;
            dragY = e.clientY - rect.top;
            const onMove = (e2) => { popup.style.left = (e2.clientX - dragX) + 'px'; popup.style.top = (e2.clientY - dragY) + 'px'; };
            const onUp = () => { document.removeEventListener('mousemove', onMove); document.removeEventListener('mouseup', onUp); };
            document.addEventListener('mousemove', onMove);
            document.addEventListener('mouseup', onUp);
        });

        this._oncoprintData = { topGenes, sortedCLs, allFilteredCLs: filteredCLs, cellW, cellH, boxAreaW, labelW, boxW, boxGap };
        this._oncoprintContext = context;

        const self = this;
        const canvas = document.getElementById('oncoprintCanvas');
        const ctx = canvas.getContext('2d');

        const drawOncoprint = () => {
            ctx.clearRect(0, 0, totalW, totalH);
            ctx.fillStyle = '#f9fafb';
            ctx.fillRect(0, 0, totalW, totalH);

            topGenes.forEach((g, rowIdx) => {
                const y = rowIdx * cellH;
                const isSelected = g.gene === currentHotspot;
                const filterState = self._oncoprintFilters[g.gene] || 'none';
                const bx1 = 2;
                const bx2 = 2 + boxW + boxGap;
                const by = y + 2;
                const bh = cellH - 4;

                ctx.fillStyle = filterState === 'mut' ? '#16a34a' : '#e5e7eb';
                ctx.fillRect(bx1, by, boxW, bh);
                ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 0.5;
                ctx.strokeRect(bx1, by, boxW, bh);
                if (filterState === 'mut') {
                    ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5;
                    ctx.beginPath();
                    ctx.moveTo(bx1 + 2, by + bh / 2); ctx.lineTo(bx1 + boxW / 2 - 1, by + bh - 3); ctx.lineTo(bx1 + boxW - 2, by + 2);
                    ctx.stroke();
                }

                ctx.fillStyle = filterState === 'wt' ? '#dc2626' : '#e5e7eb';
                ctx.fillRect(bx2, by, boxW, bh);
                ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 0.5;
                ctx.strokeRect(bx2, by, boxW, bh);
                if (filterState === 'wt') {
                    ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5;
                    ctx.beginPath();
                    ctx.moveTo(bx2 + 2, by + 2); ctx.lineTo(bx2 + boxW - 2, by + bh - 2);
                    ctx.moveTo(bx2 + boxW - 2, by + 2); ctx.lineTo(bx2 + 2, by + bh - 2);
                    ctx.stroke();
                }

                const hasFilter = filterState !== 'none';
                ctx.fillStyle = isSelected ? '#5a9f4a' : filterState === 'mut' ? '#16a34a' : filterState === 'wt' ? '#dc2626' : '#374151';
                ctx.font = (isSelected || hasFilter) ? 'bold 10px Arial' : '10px Arial';
                ctx.textAlign = 'right';
                ctx.textBaseline = 'middle';
                ctx.fillText(g.gene, boxAreaW + labelW - 4, y + cellH / 2);

                ctx.fillStyle = '#9ca3af';
                ctx.font = '8px Arial';
                ctx.textAlign = 'left';
                ctx.fillText(`${g.n}`, boxAreaW + labelW + gridW + 2, y + cellH / 2);

                sortedCLs.forEach((cl, colIdx) => {
                    const x = boxAreaW + labelW + colIdx * cellW;
                    const mutLevel = g.muts[cl] || 0;
                    if (mutLevel > 0) {
                        ctx.fillStyle = mutLevel >= 2 ? '#1e40af' : '#3b82f6';
                    } else {
                        ctx.fillStyle = '#f3f4f6';
                    }
                    ctx.fillRect(x, y + 1, cellW - 1, cellH - 2);
                });
            });

            const activeFilters = Object.entries(self._oncoprintFilters).filter(([, v]) => v !== 'none');
            const statusEl = document.getElementById('oncoprintStatus');
            if (activeFilters.length === 0) {
                statusEl.innerHTML = '<span style="color:#9ca3af;">Click <span style="color:#16a34a;">&#9632;</span> to include or <span style="color:#dc2626;">&#9632;</span> to exclude mutated cells.</span>';
            } else {
                let matchCount = 0;
                for (const cl of filteredCLs) {
                    let passes = true;
                    for (const [gene, state] of activeFilters) {
                        const muts = self.mutations.geneData[gene]?.mutations;
                        const isMut = muts && muts[cl] > 0;
                        if (state === 'mut' && !isMut) { passes = false; break; }
                        if (state === 'wt' && isMut) { passes = false; break; }
                    }
                    if (passes) matchCount++;
                }
                const tags = activeFilters.map(([gene, state]) =>
                    `<span style="display:inline-flex;align-items:center;gap:2px;padding:1px 6px;border-radius:10px;font-size:10px;background:${state === 'mut' ? '#dcfce7' : '#fef2f2'};color:${state === 'mut' ? '#16a34a' : '#dc2626'};border:1px solid ${state === 'mut' ? '#86efac' : '#fecaca'};">${gene} ${state === 'mut' ? '\u2713' : '\u2717'}<button onclick="app._oncoprintClearGene('${gene}')" style="background:none;border:none;cursor:pointer;font-size:10px;color:#999;padding:0 0 0 2px;">\u00d7</button></span>`
                ).join('');
                statusEl.innerHTML = `${tags} <span style="color:#6b7280;">${matchCount}/${filteredCLs.length} CLs</span> <button onclick="app._oncoprintApplyFilters()" style="padding:2px 8px;font-size:10px;background:#5a9f4a;color:white;border:none;border-radius:4px;cursor:pointer;">Apply</button> <button onclick="app._oncoprintClearAll()" style="padding:2px 8px;font-size:10px;background:#f3f4f6;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;">Clear</button>`;
            }
        };

        drawOncoprint();

        canvas.addEventListener('click', (e) => {
            const rect = canvas.getBoundingClientRect();
            const y = e.clientY - rect.top;
            const x = e.clientX - rect.left;
            const rowIdx = Math.floor(y / cellH);
            if (rowIdx < 0 || rowIdx >= topGenes.length) return;

            const gene = topGenes[rowIdx].gene;
            const bx1 = 2, bx2 = 2 + boxW + boxGap;

            if (x >= bx1 && x <= bx1 + boxW) {
                if (this._oncoprintFilters[gene] === 'mut') delete this._oncoprintFilters[gene];
                else this._oncoprintFilters[gene] = 'mut';
                this._oncoprintSyncFilters();
                drawOncoprint();
            } else if (x >= bx2 && x <= bx2 + boxW) {
                if (this._oncoprintFilters[gene] === 'wt') delete this._oncoprintFilters[gene];
                else this._oncoprintFilters[gene] = 'wt';
                this._oncoprintSyncFilters();
                drawOncoprint();
            } else if (x < boxAreaW + labelW) {
                document.getElementById('mutationHotspotSelect').value = gene;
                document.getElementById('tissueBreakdownBtn').style.display = 'inline-block';
                this.showOncoprint(this._oncoprintContext);
            }
        });

        canvas.addEventListener('mousemove', (e) => {
            const canvasRect = canvas.getBoundingClientRect();
            const y = e.clientY - canvasRect.top;
            const x = e.clientX - canvasRect.left;
            const rowIdx = Math.floor(y / cellH);
            if (rowIdx >= 0 && rowIdx < topGenes.length) {
                canvas.style.cursor = 'pointer';
                const bx1 = 2, bx2 = 2 + boxW + boxGap;
                const colIdx = Math.floor((x - boxAreaW - labelW) / cellW);
                if (x >= bx1 && x <= bx1 + boxW) {
                    const isMut = this._oncoprintFilters[topGenes[rowIdx].gene] === 'mut';
                    canvas.title = `${topGenes[rowIdx].gene} \u2014 ${isMut ? 'remove include filter' : 'require mutated'}`;
                } else if (x >= bx2 && x <= bx2 + boxW) {
                    const isWt = this._oncoprintFilters[topGenes[rowIdx].gene] === 'wt';
                    canvas.title = `${topGenes[rowIdx].gene} \u2014 ${isWt ? 'remove exclude filter' : 'exclude mutated'}`;
                } else if (x < boxAreaW + labelW) {
                    canvas.title = `${topGenes[rowIdx].gene} (${topGenes[rowIdx].n} mut) \u2014 click to set as hotspot`;
                } else if (colIdx >= 0 && colIdx < sortedCLs.length) {
                    canvas.title = `${topGenes[rowIdx].gene} \u00b7 ${this.getCellLineName(sortedCLs[colIdx])} \u00b7 ${topGenes[rowIdx].muts[sortedCLs[colIdx]] > 0 ? 'Mutated' : 'WT'}`;
                } else {
                    canvas.title = '';
                }
            } else {
                canvas.style.cursor = 'default';
                canvas.title = '';
            }
        });

        const escHandler = (e) => { if (e.key === 'Escape') { popup.remove(); document.removeEventListener('keydown', escHandler); } };
        document.addEventListener('keydown', escHandler);

        setTimeout(() => {
            const outsideHandler = (e) => {
                const upsetPopup = document.getElementById('upsetPopup');
                if (!popup.contains(e.target) && e.target.id !== 'oncoprintBtn' && (!upsetPopup || !upsetPopup.contains(e.target))) {
                    popup.remove();
                    document.removeEventListener('mousedown', outsideHandler);
                }
            };
            document.addEventListener('mousedown', outsideHandler);
        }, 100);
    }

    _oncoprintSyncFilters() {
        const filters = Object.entries(this._oncoprintFilters || {}).filter(([, v]) => v !== 'none');
        this._activeOncoprintFilters = filters.length > 0 ? filters.map(([gene, state]) => ({ gene, state })) : null;

        const controls = document.getElementById('paramHotspotControls');
        const label = document.getElementById('oncoprintFilterLabel');
        if (controls && label) {
            if (filters.length > 0) {
                controls.style.opacity = '0.3';
                controls.style.pointerEvents = 'none';
                const parts = filters.map(([gene, state]) => `${gene} ${state === 'mut' ? 'Mut' : 'WT'}`);
                label.innerHTML = `Oncoprint filter: ${parts.join(', ')} <span style="font-size:9px;color:#9ca3af;">(click to clear)</span>`;
                label.style.display = '';
                label.onclick = () => { this._oncoprintClearAll(); };
            } else {
                controls.style.opacity = '';
                controls.style.pointerEvents = '';
                label.style.display = 'none';
            }
        }

        if (this._oncoprintContext === 'clb') {
            this.renderCellLineList();
            this.updateClbFilterCounts();
        } else {
            this._updateLineageFilterCounts();
        }
    }

    _updateLineageFilterCounts() {
        const select = document.getElementById('lineageFilter');
        if (!select || !this.cellLineMetadata?.lineage) return;
        const currentVal = select.value;
        const cellLines = this.metadata.cellLines;
        const counts = {};
        let total = 0;
        for (const cl of cellLines) {
            if (!this._cellLinePassesOncoprintFilters(cl)) continue;
            if (this.excludedTissues && this.excludedTissues.size > 0) {
                const lin = this.cellLineMetadata.lineage[cl];
                if (lin && this.excludedTissues.has(lin)) continue;
            }
            const lin = this.cellLineMetadata.lineage[cl];
            if (lin) counts[lin] = (counts[lin] || 0) + 1;
            total++;
        }
        select.innerHTML = `<option value="">All lineages (n=${total})</option>`;
        Object.keys(counts).sort((a, b) => counts[b] - counts[a]).forEach(lin => {
            const opt = document.createElement('option');
            opt.value = lin;
            opt.textContent = `${lin} (n=${counts[lin]})`;
            select.appendChild(opt);
        });
        if (currentVal) select.value = currentVal;
    }

    _oncoprintClearGene(gene) {
        delete this._oncoprintFilters[gene];
        this._oncoprintSyncFilters();
        this.showOncoprint(this._oncoprintContext);
    }

    _oncoprintClearAll() {
        this._oncoprintFilters = {};
        this._oncoprintSyncFilters();
        this.showOncoprint(this._oncoprintContext);
    }

    _oncoprintApplyFilters() {
        const filters = Object.entries(this._oncoprintFilters || {}).filter(([, v]) => v !== 'none');
        if (filters.length === 0) return;

        this._activeOncoprintFilters = filters.map(([gene, state]) => ({ gene, state }));

        const mutGenes = filters.filter(([, v]) => v === 'mut');
        const wtGenes = filters.filter(([, v]) => v === 'wt');

        if (mutGenes.length > 0) {
            document.getElementById('mutationHotspotSelect').value = mutGenes[0][0];
            document.getElementById('tissueBreakdownBtn').style.display = 'inline-block';
        }

        if (wtGenes.length > 0) {
            document.getElementById('paramHotspotGene').value = wtGenes[0][0];
            document.getElementById('paramHotspotLevel').value = '0';
        } else if (mutGenes.length > 1) {
            document.getElementById('paramHotspotGene').value = mutGenes[1][0];
            document.getElementById('paramHotspotLevel').value = '1+2';
        }

        const mode = document.querySelector('input[name="analysisMode"]:checked')?.value;
        if (mode === 'mutation') {
            this.runAnalysis();
        }
        document.getElementById('oncoprintPopup')?.remove();
    }

    _cellLinePassesOncoprintFilters(cellLine) {
        if (!this._activeOncoprintFilters || this._activeOncoprintFilters.length === 0) return true;
        for (const { gene, state } of this._activeOncoprintFilters) {
            const muts = this.mutations?.geneData?.[gene]?.mutations;
            const isMut = muts && muts[cellLine] > 0;
            if (state === 'mut' && !isMut) return false;
            if (state === 'wt' && isMut) return false;
        }
        return true;
    }

    _oncoprintExport(format) {
        const canvas = document.getElementById('oncoprintCanvas');
        if (!canvas) return;

        if (format === 'csv') {
            const data = this._oncoprintData;
            if (!data) return;
            let csv = 'Gene,' + data.sortedCLs.map(cl => this.getCellLineName(cl)).join(',') + '\n';
            data.topGenes.forEach(g => {
                csv += g.gene + ',' + data.sortedCLs.map(cl => g.muts[cl] || 0).join(',') + '\n';
            });
            const blob = new Blob([csv], { type: 'text/csv' });
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = 'oncoprint.csv';
            a.click();
            URL.revokeObjectURL(a.href);
        } else if (format === 'png') {
            const scale = 4;
            const exportCanvas = document.createElement('canvas');
            exportCanvas.width = canvas.width * scale;
            exportCanvas.height = canvas.height * scale;
            const ectx = exportCanvas.getContext('2d');
            ectx.scale(scale, scale);
            ectx.fillStyle = 'white';
            ectx.fillRect(0, 0, canvas.width, canvas.height);
            ectx.drawImage(canvas, 0, 0);
            const a = document.createElement('a');
            a.href = exportCanvas.toDataURL('image/png');
            a.download = 'oncoprint.png';
            a.click();
        } else {
            const w = canvas.width, h = canvas.height;
            let svg = `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">\n`;
            svg += `<rect width="${w}" height="${h}" fill="white"/>\n`;
            const data = this._oncoprintData;
            if (data) {
                const { topGenes: tg, sortedCLs: cls, cellW, cellH: cH, boxAreaW: baw, labelW: lw, boxW: bw, boxGap: bg } = data;
                tg.forEach((g, rowIdx) => {
                    const y = rowIdx * cH;
                    const fs = this._oncoprintFilters[g.gene] || 'none';
                    svg += `<rect x="2" y="${y + 2}" width="${bw}" height="${cH - 4}" fill="${fs === 'mut' ? '#16a34a' : '#e5e7eb'}" stroke="#9ca3af" stroke-width="0.5"/>\n`;
                    svg += `<rect x="${2 + bw + bg}" y="${y + 2}" width="${bw}" height="${cH - 4}" fill="${fs === 'wt' ? '#dc2626' : '#e5e7eb'}" stroke="#9ca3af" stroke-width="0.5"/>\n`;
                    const labelColor = fs === 'mut' ? '#16a34a' : fs === 'wt' ? '#dc2626' : '#374151';
                    svg += `<text x="${baw + lw - 4}" y="${y + cH / 2}" text-anchor="end" dominant-baseline="central" font-family="Arial" font-size="10" fill="${labelColor}" ${fs !== 'none' ? 'font-weight="bold"' : ''}>${g.gene}</text>\n`;
                    svg += `<text x="${baw + lw + cls.length * cellW + 2}" y="${y + cH / 2}" font-family="Arial" font-size="8" fill="#9ca3af" dominant-baseline="central">${g.n}</text>\n`;
                    cls.forEach((cl, colIdx) => {
                        const x = baw + lw + colIdx * cellW;
                        const mutLevel = g.muts[cl] || 0;
                        svg += `<rect x="${x}" y="${y + 1}" width="${cellW - 1}" height="${cH - 2}" fill="${mutLevel >= 2 ? '#1e40af' : mutLevel > 0 ? '#3b82f6' : '#f3f4f6'}"/>\n`;
                    });
                });
            }
            svg += '</svg>';
            const blob = new Blob([svg], { type: 'image/svg+xml' });
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = 'oncoprint.svg';
            a.click();
            URL.revokeObjectURL(a.href);
        }
    }

    getTissueBreakdownForHotspot(gene) {
        if (!this.mutations?.geneData?.[gene] || !this.cellLineMetadata?.lineage) return [];
        const mutations = this.mutations.geneData[gene].mutations;
        const cellLines = this.metadata.cellLines;
        const tissueMap = {};

        cellLines.forEach(cl => {
            const lineage = this.cellLineMetadata.lineage[cl];
            if (!lineage) return;
            if (!tissueMap[lineage]) tissueMap[lineage] = { lineage, nMut: 0, nWT: 0 };
            if (mutations[cl] && mutations[cl] > 0) {
                tissueMap[lineage].nMut++;
            } else {
                tissueMap[lineage].nWT++;
            }
        });

        return Object.values(tissueMap).sort((a, b) => b.nMut - a.nMut);
    }

    getTissueBreakdownForTranslocation(gene) {
        if (!this.translocations?.geneData?.[gene] || !this.cellLineMetadata?.lineage) return [];
        const translocations = this.translocations.geneData[gene].translocations;
        const cellLines = this.metadata.cellLines;
        const tissueMap = {};

        cellLines.forEach(cl => {
            const lineage = this.cellLineMetadata.lineage[cl];
            if (!lineage) return;
            if (!tissueMap[lineage]) tissueMap[lineage] = { lineage, nMut: 0, nWT: 0 };
            if (translocations[cl] && translocations[cl] > 0) {
                tissueMap[lineage].nMut++;
            } else {
                tissueMap[lineage].nWT++;
            }
        });

        return Object.values(tissueMap).sort((a, b) => b.nMut - a.nMut);
    }

    showTissueBreakdownPopup(type) {
        this.hideTissueBreakdownPopup();
        const isTransloc = type === 'translocation';
        const gene = isTransloc
            ? document.getElementById('translocationHotspotSelect').value
            : document.getElementById('mutationHotspotSelect').value;
        if (!gene) return;

        const breakdown = isTransloc
            ? this.getTissueBreakdownForTranslocation(gene)
            : this.getTissueBreakdownForHotspot(gene);
        if (breakdown.length === 0) return;

        const currentLineage = document.getElementById('lineageFilter').value;
        const mutLabel = isTransloc ? 'Fused' : 'Mut';

        // Pre-compute sub-tissue breakdown for each tissue
        const subBreakdowns = {};
        if (this.cellLineMetadata?.primaryDisease) {
            const cellLines = this.metadata.cellLines;
            const mutSource = isTransloc ? this.translocations?.geneData?.[gene]?.translocations
                : this.mutations?.geneData?.[gene]?.mutations;
            if (mutSource) {
                cellLines.forEach(cl => {
                    const lin = this.cellLineMetadata.lineage?.[cl];
                    const sub = this.cellLineMetadata.primaryDisease?.[cl];
                    if (!lin || !sub) return;
                    if (!subBreakdowns[lin]) subBreakdowns[lin] = {};
                    if (!subBreakdowns[lin][sub]) subBreakdowns[lin][sub] = { nMut: 0, nWT: 0 };
                    if (mutSource[cl] > 0) subBreakdowns[lin][sub].nMut++;
                    else subBreakdowns[lin][sub].nWT++;
                });
            }
        }

        const popup = document.createElement('div');
        popup.id = 'tissueBreakdownPopup';
        popup.style.cssText = 'position: fixed; z-index: 10000; background: white; border: 1px solid #d1d5db; border-radius: 8px; box-shadow: 0 8px 24px rgba(0,0,0,0.15); padding: 0; min-width: 340px; max-width: 420px; display: flex; flex-direction: column;';

        // Position near the button, clamped to viewport
        const btn = document.getElementById(isTransloc ? 'translocationTissueBreakdownBtn' : 'tissueBreakdownBtn');
        const rect = btn.getBoundingClientRect();
        const vw = window.innerWidth;
        const vh = window.innerHeight;
        const margin = 10;
        const popupWidth = 400;

        // Try below the button first; if not enough space, go above
        let top = rect.bottom + 6;
        let maxH = vh - top - margin;
        if (maxH < 200) {
            // Not enough room below — place above
            maxH = rect.top - 6 - margin;
            top = rect.top - 6 - Math.min(480, maxH);
            maxH = Math.min(480, maxH);
        } else {
            maxH = Math.min(480, maxH);
        }

        popup.style.top = Math.max(margin, top) + 'px';
        popup.style.left = Math.max(margin, Math.min(rect.left - 100, vw - popupWidth - margin)) + 'px';
        popup.style.maxHeight = maxH + 'px';

        const maxMut = Math.max(...breakdown.map(t => t.nMut));

        // Header
        let html = `<div style="padding: 10px 14px 8px; border-bottom: 1px solid #e5e7eb; display: flex; justify-content: space-between; align-items: center;">
            <span style="font-weight: 600; font-size: 13px; color: #1f2937;">${gene} — Tissue Breakdown</span>
            <button id="tbCloseBtn" style="background: none; border: none; cursor: pointer; font-size: 18px; color: #6b7280; line-height: 1; padding: 0 2px;">&times;</button>
        </div>`;

        // Table
        html += `<div style="overflow-y: auto; flex: 1; padding: 4px 0;">
            <table style="width: 100%; border-collapse: collapse; font-size: 12px;">
            <thead><tr style="position: sticky; top: 0; background: #f9fafb; z-index: 1;">
                <th style="padding: 4px 8px; text-align: left; width: 28px;"><input type="checkbox" id="tbSelectAll" title="Select all"></th>
                <th style="padding: 4px 4px; text-align: left;">Tissue</th>
                <th style="padding: 4px 6px; text-align: right; color: #dc2626;">${mutLabel}</th>
                <th style="padding: 4px 6px; text-align: right; color: #6b7280;">WT</th>
                <th style="padding: 4px 8px; text-align: left; width: 90px;"></th>
                <th style="padding: 4px 8px; text-align: right; width: 44px;">%</th>
            </tr></thead><tbody>`;

        breakdown.forEach(t => {
            const total = t.nMut + t.nWT;
            const pct = total > 0 ? (t.nMut / total * 100).toFixed(1) : '0.0';
            const barWidth = maxMut > 0 ? (t.nMut / maxMut * 100) : 0;
            const isCurrentFilter = (currentLineage && t.lineage === currentLineage);
            const rowBg = isCurrentFilter ? 'background: #ecfdf5;' : '';
            const checked = isCurrentFilter ? ' checked' : '';
            const hasSubs = subBreakdowns[t.lineage] && Object.keys(subBreakdowns[t.lineage]).length > 1;

            html += `<tr class="tb-row" data-tissue="${t.lineage}" style="cursor: pointer; ${rowBg}" onmouseenter="this.style.background=this.style.background||'#f3f4f6'" onmouseleave="this.style.background='${isCurrentFilter ? '#ecfdf5' : ''}'">
                <td style="padding: 3px 8px;"><input type="checkbox" class="tb-check" value="${t.lineage}"${checked}></td>
                <td style="padding: 3px 4px; font-weight: ${t.nMut > 0 ? '500' : '400'}; color: ${t.nMut > 0 ? '#1f2937' : '#9ca3af'};">${hasSubs ? '<span class="tb-expand" style="font-size:9px;color:#9ca3af;margin-right:2px;">&#9654;</span>' : ''}${t.lineage}</td>
                <td style="padding: 3px 6px; text-align: right; color: #dc2626; font-weight: 600;">${t.nMut}</td>
                <td style="padding: 3px 6px; text-align: right; color: #6b7280;">${t.nWT}</td>
                <td style="padding: 3px 8px;"><div style="background: #fee2e2; border-radius: 2px; height: 10px; width: 100%;"><div style="background: #ef4444; border-radius: 2px; height: 10px; width: ${barWidth}%;"></div></div></td>
                <td style="padding: 3px 8px; text-align: right; color: #6b7280; font-size: 11px;">${pct}%</td>
            </tr>`;

            // Sub-tissue rows (hidden by default)
            if (hasSubs) {
                const subs = Object.entries(subBreakdowns[t.lineage]).sort((a, b) => b[1].nMut - a[1].nMut);
                subs.forEach(([sub, counts]) => {
                    const subTotal = counts.nMut + counts.nWT;
                    const subPct = subTotal > 0 ? (counts.nMut / subTotal * 100).toFixed(1) : '0.0';
                    const subBarW = maxMut > 0 ? (counts.nMut / maxMut * 100) : 0;
                    html += `<tr class="tb-sub-row" data-parent="${t.lineage}" data-subtype="${sub}" style="display:none; cursor:pointer; background:#fafafa;" onmouseenter="this.style.background='#f0f0f0'" onmouseleave="this.style.background='#fafafa'">
                        <td style="padding: 2px 8px 2px 16px;"><input type="checkbox" class="tb-sub-check" value="${sub}" data-parent="${t.lineage}"></td>
                        <td style="padding: 2px 4px 2px 8px; font-size:11px; color:#6b7280;">${sub}</td>
                        <td style="padding: 2px 6px; text-align: right; color: #dc2626; font-size:11px;">${counts.nMut}</td>
                        <td style="padding: 2px 6px; text-align: right; color: #6b7280; font-size:11px;">${counts.nWT}</td>
                        <td style="padding: 2px 8px;"><div style="background: #fee2e2; border-radius: 2px; height: 8px; width: 100%;"><div style="background: #f87171; border-radius: 2px; height: 8px; width: ${subBarW}%;"></div></div></td>
                        <td style="padding: 2px 8px; text-align: right; color: #9ca3af; font-size: 10px;">${subPct}%</td>
                    </tr>`;
                });
            }
        });

        html += `</tbody></table></div>`;

        // Footer
        html += `<div style="padding: 8px 14px; border-top: 1px solid #e5e7eb; display: flex; justify-content: space-between; align-items: center;">
            <span id="tbSelectionCount" style="font-size: 11px; color: #6b7280;">0 selected</span>
            <div style="display: flex; gap: 6px;">
                <button id="tbClearBtn" style="padding: 4px 12px; font-size: 12px; background: #f3f4f6; color: #374151; border: 1px solid #d1d5db; border-radius: 4px; cursor: pointer;">Clear</button>
                <button id="tbApplyBtn" style="padding: 4px 12px; font-size: 12px; background: #5a9f4a; color: white; border: none; border-radius: 4px; cursor: pointer; font-weight: 500;">Apply Filter</button>
            </div>
        </div>`;

        popup.innerHTML = html;
        document.body.appendChild(popup);

        // Row click: toggle checkbox, or expand sub-tissues if clicking the tissue name
        popup.querySelectorAll('.tb-row').forEach(row => {
            row.addEventListener('click', (e) => {
                if (e.target.type === 'checkbox') return;
                const tissue = row.dataset.tissue;
                const expandArrow = row.querySelector('.tb-expand');

                // If clicking tissue name area and has sub-tissues, toggle expansion
                if (expandArrow && (e.target.closest('td') === row.cells[1])) {
                    const subRows = popup.querySelectorAll(`.tb-sub-row[data-parent="${tissue}"]`);
                    const isExpanded = subRows[0]?.style.display !== 'none';
                    subRows.forEach(sr => sr.style.display = isExpanded ? 'none' : '');
                    expandArrow.innerHTML = isExpanded ? '&#9654;' : '&#9660;';
                    return;
                }

                const cb = row.querySelector('.tb-check');
                cb.checked = !cb.checked;
                this.updateTBSelectionCount();
            });
        });

        // Checkbox change (tissue)
        popup.querySelectorAll('.tb-check').forEach(cb => {
            cb.addEventListener('change', () => this.updateTBSelectionCount());
        });

        // Sub-tissue row click
        popup.querySelectorAll('.tb-sub-row').forEach(row => {
            row.addEventListener('click', (e) => {
                if (e.target.type === 'checkbox') return;
                const cb = row.querySelector('.tb-sub-check');
                cb.checked = !cb.checked;
                this.updateTBSelectionCount();
            });
        });
        popup.querySelectorAll('.tb-sub-check').forEach(cb => {
            cb.addEventListener('change', () => this.updateTBSelectionCount());
        });

        // Select all
        document.getElementById('tbSelectAll').addEventListener('change', (e) => {
            popup.querySelectorAll('.tb-check').forEach(cb => { cb.checked = e.target.checked; });
            popup.querySelectorAll('.tb-sub-check').forEach(cb => { cb.checked = false; });
            this.updateTBSelectionCount();
        });

        // Close button
        document.getElementById('tbCloseBtn').addEventListener('click', () => this.hideTissueBreakdownPopup());

        // Clear button
        document.getElementById('tbClearBtn').addEventListener('click', () => {
            popup.querySelectorAll('.tb-check').forEach(cb => { cb.checked = false; });
            document.getElementById('tbSelectAll').checked = false;
            this.updateTBSelectionCount();
        });

        // Apply button
        document.getElementById('tbApplyBtn').addEventListener('click', () => {
            const selectedTissues = [...popup.querySelectorAll('.tb-check:checked')].map(cb => cb.value);
            const selectedSubtypes = [...popup.querySelectorAll('.tb-sub-check:checked')].map(cb => ({
                subtype: cb.value,
                parent: cb.dataset.parent
            }));
            this.applyTissueBreakdownSelection(selectedTissues, selectedSubtypes);
            this.hideTissueBreakdownPopup();
        });

        // Escape key
        this._tbEscHandler = (e) => { if (e.key === 'Escape') this.hideTissueBreakdownPopup(); };
        document.addEventListener('keydown', this._tbEscHandler);

        // Click outside
        setTimeout(() => {
            this._tbOutsideHandler = (e) => {
                const p = document.getElementById('tissueBreakdownPopup');
                if (p && !p.contains(e.target) && e.target.id !== 'tissueBreakdownBtn' && e.target.id !== 'translocationTissueBreakdownBtn') {
                    this.hideTissueBreakdownPopup();
                }
            };
            document.addEventListener('mousedown', this._tbOutsideHandler);
        }, 0);

        this.updateTBSelectionCount();
    }

    hideTissueBreakdownPopup() {
        const popup = document.getElementById('tissueBreakdownPopup');
        if (popup) popup.remove();
        if (this._tbEscHandler) {
            document.removeEventListener('keydown', this._tbEscHandler);
            this._tbEscHandler = null;
        }
        if (this._tbOutsideHandler) {
            document.removeEventListener('mousedown', this._tbOutsideHandler);
            this._tbOutsideHandler = null;
        }
    }

    applyTissueBreakdownSelection(selectedTissues, selectedSubtypes = []) {
        const lineageSelect = document.getElementById('lineageFilter');
        const lineageGroup = document.getElementById('lineageFilterGroup');
        const selectedSet = new Set(selectedTissues);

        if (selectedTissues.length === 0) {
            // Clear all tissue filters
            lineageSelect.value = '';
            lineageSelect.disabled = false;
            lineageSelect.style.opacity = '';
            this.excludedTissues = new Set();
            document.querySelectorAll('#tissueExcludeList input[type="checkbox"]').forEach(cb => {
                cb.checked = false;
            });
            // Remove override label
            const overrideLabel = lineageGroup?.querySelector('.tb-override-label');
            if (overrideLabel) overrideLabel.remove();
        } else if (selectedTissues.length === 1) {
            // Single tissue: use lineage dropdown directly
            lineageSelect.value = selectedTissues[0];
            lineageSelect.disabled = false;
            lineageSelect.style.opacity = '';
            this.excludedTissues = new Set();
            document.querySelectorAll('#tissueExcludeList input[type="checkbox"]').forEach(cb => {
                cb.checked = false;
            });
            const overrideLabel = lineageGroup?.querySelector('.tb-override-label');
            if (overrideLabel) overrideLabel.remove();
        } else {
            // Multiple tissues: disable lineage dropdown, use excludedTissues
            lineageSelect.value = '';
            lineageSelect.disabled = true;
            lineageSelect.style.opacity = '0.5';
            const allLineages = this.cellLineMetadata?.lineage
                ? new Set(Object.values(this.cellLineMetadata.lineage))
                : new Set();
            const excludeSet = new Set();
            allLineages.forEach(t => {
                if (!selectedSet.has(t)) excludeSet.add(t);
            });
            this.excludedTissues = excludeSet;
            document.querySelectorAll('#tissueExcludeList input[type="checkbox"]').forEach(cb => {
                cb.checked = excludeSet.has(cb.value);
            });
            // Show override label
            let overrideLabel = lineageGroup?.querySelector('.tb-override-label');
            if (!overrideLabel && lineageGroup) {
                overrideLabel = document.createElement('div');
                overrideLabel.className = 'tb-override-label';
                overrideLabel.style.cssText = 'font-size: 11px; color: #5a9f4a; margin-top: 2px; cursor: pointer;';
                overrideLabel.title = 'Click to clear tissue selection';
                overrideLabel.addEventListener('click', () => {
                    this.applyTissueBreakdownSelection([]);
                });
                lineageGroup.appendChild(overrideLabel);
            }
            if (overrideLabel) {
                overrideLabel.textContent = `Tissue filter: ${selectedTissues.join(', ')} (click to clear)`;
            }
        }

        // Handle sub-tissue selection
        if (selectedSubtypes.length > 0) {
            const parents = [...new Set(selectedSubtypes.map(s => s.parent))];
            if (parents.length === 1 && (selectedTissues.length === 0 || (selectedTissues.length === 1 && selectedTissues[0] === parents[0]))) {
                lineageSelect.value = parents[0];
                lineageSelect.disabled = false;
                lineageSelect.style.opacity = '';
                this.excludedTissues = new Set();
                this._pendingSubtypeSelection = selectedSubtypes.map(s => s.subtype);
            }
            const subtypeNames = selectedSubtypes.map(s => s.subtype);
            let overrideLabel = lineageGroup?.querySelector('.tb-override-label');
            if (!overrideLabel && lineageGroup) {
                overrideLabel = document.createElement('div');
                overrideLabel.className = 'tb-override-label';
                overrideLabel.style.cssText = 'font-size: 11px; color: #5a9f4a; margin-top: 2px; cursor: pointer;';
                overrideLabel.title = 'Click to clear tissue selection';
                overrideLabel.addEventListener('click', () => { this.applyTissueBreakdownSelection([]); });
                lineageGroup.appendChild(overrideLabel);
            }
            if (overrideLabel) {
                overrideLabel.textContent = `Subtype: ${subtypeNames.join(', ')} (click to clear)`;
            }
        }

        // Trigger UI updates
        this.updateSubLineageFilter();

        // Apply pending subtype selection after dropdown is populated
        if (this._pendingSubtypeSelection) {
            const subSelect = document.getElementById('subLineageFilter');
            if (subSelect && this._pendingSubtypeSelection.length === 1) {
                subSelect.value = this._pendingSubtypeSelection[0];
            }
            this._pendingSubtypeSelection = null;
        }

        this.updateHotspotCountsForCurrentFilters();
    }

    showGeneralTissueBreakdownPopup() {
        this.hideTissueBreakdownPopup();
        if (!this.cellLineMetadata?.lineage) return;

        // Count cell lines per tissue
        const tissueCounts = {};
        const cellLines = this.metadata?.cellLines || [];
        cellLines.forEach(cellLine => {
            const lineage = this.cellLineMetadata.lineage[cellLine] || 'Unknown';
            if (!tissueCounts[lineage]) tissueCounts[lineage] = 0;
            tissueCounts[lineage]++;
        });

        const breakdown = Object.entries(tissueCounts)
            .map(([lineage, count]) => ({ lineage, count }))
            .sort((a, b) => b.count - a.count);

        if (breakdown.length === 0) return;

        const currentLineage = document.getElementById('lineageFilter').value;
        const maxCount = Math.max(...breakdown.map(t => t.count));

        const popup = document.createElement('div');
        popup.id = 'tissueBreakdownPopup';
        popup.style.cssText = 'position: fixed; z-index: 10000; background: white; border: 1px solid #d1d5db; border-radius: 8px; box-shadow: 0 8px 24px rgba(0,0,0,0.15); padding: 0; min-width: 300px; max-width: 380px; display: flex; flex-direction: column;';

        const btn = document.getElementById('generalTissueBreakdownBtn');
        const rect = btn.getBoundingClientRect();
        const vw = window.innerWidth;
        const vh = window.innerHeight;
        const margin = 10;
        let top = rect.bottom + 6;
        let maxH = vh - top - margin;
        if (maxH < 200) { maxH = rect.top - 6 - margin; top = rect.top - 6 - Math.min(480, maxH); maxH = Math.min(480, maxH); } else { maxH = Math.min(480, maxH); }
        popup.style.top = Math.max(margin, top) + 'px';
        popup.style.left = Math.max(margin, Math.min(rect.left - 50, vw - 380 - margin)) + 'px';
        popup.style.maxHeight = maxH + 'px';

        let html = `<div style="padding: 10px 14px 8px; border-bottom: 1px solid #e5e7eb; display: flex; justify-content: space-between; align-items: center;">
            <span style="font-weight: 600; font-size: 13px; color: #1f2937;">Cell Lines per Tissue (${cellLines.length} total)</span>
            <button id="tbCloseBtn" style="background: none; border: none; cursor: pointer; font-size: 18px; color: #6b7280; line-height: 1; padding: 0 2px;">&times;</button>
        </div>`;
        html += `<div style="overflow-y: auto; flex: 1; padding: 4px 0;">
            <table style="width: 100%; border-collapse: collapse; font-size: 12px;">
            <thead><tr style="position: sticky; top: 0; background: #f9fafb; z-index: 1;">
                <th style="padding: 4px 8px; text-align: left;">Tissue</th>
                <th style="padding: 4px 6px; text-align: right;">N</th>
                <th style="padding: 4px 8px; text-align: left; width: 90px;"></th>
            </tr></thead><tbody>`;

        breakdown.forEach(t => {
            const barWidth = maxCount > 0 ? (t.count / maxCount * 100) : 0;
            const isCurrentFilter = (currentLineage && t.lineage === currentLineage);
            const rowBg = isCurrentFilter ? 'background: #ecfdf5;' : '';
            html += `<tr style="cursor: pointer; ${rowBg}" data-tissue="${t.lineage}" onmouseenter="this.style.background=this.style.background||'#f3f4f6'" onmouseleave="this.style.background='${isCurrentFilter ? '#ecfdf5' : ''}'">
                <td style="padding: 3px 8px; font-weight: 500;">${t.lineage}</td>
                <td style="padding: 3px 6px; text-align: right; font-weight: 600;">${t.count}</td>
                <td style="padding: 3px 8px;"><div style="background: #e2f4de; border-radius: 2px; height: 10px; width: 100%;"><div style="background: #5a9f4a; border-radius: 2px; height: 10px; width: ${barWidth}%;"></div></div></td>
            </tr>`;
        });

        html += `</tbody></table></div>`;
        popup.innerHTML = html;
        document.body.appendChild(popup);

        // Click row to set lineage filter
        popup.querySelectorAll('tr[data-tissue]').forEach(row => {
            row.addEventListener('click', () => {
                const tissue = row.dataset.tissue;
                const lineageSelect = document.getElementById('lineageFilter');
                lineageSelect.value = lineageSelect.value === tissue ? '' : tissue;
                this.updateSubLineageFilter();
                this.hideTissueBreakdownPopup();
            });
        });

        document.getElementById('tbCloseBtn').addEventListener('click', () => this.hideTissueBreakdownPopup());
        this._tbEscHandler = (e) => { if (e.key === 'Escape') this.hideTissueBreakdownPopup(); };
        document.addEventListener('keydown', this._tbEscHandler);
        setTimeout(() => {
            this._tbOutsideHandler = (e) => {
                const p = document.getElementById('tissueBreakdownPopup');
                if (p && !p.contains(e.target) && e.target.id !== 'generalTissueBreakdownBtn') this.hideTissueBreakdownPopup();
            };
            document.addEventListener('mousedown', this._tbOutsideHandler);
        }, 0);
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
        this.setupCellLineBrowserEvents();

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

        // Clear stats genes
        document.getElementById('clearStatsGenes')?.addEventListener('click', () => {
            document.getElementById('manualStatsTextarea').value = 'Gene\tLFC\tFDR\n';
            this.geneStats = null;
            this.updateGeneCount();
        });

        // Load test genes (20 genes)
        document.getElementById('loadTestGenes').addEventListener('click', () => {
            const mode = document.querySelector('input[name="analysisMode"]:checked')?.value;
            let testGenes;
            if (mode === 'synonym') {
                // Mix of mouse orthologs, old aliases, and valid human genes to demo synonym lookup
                testGenes = ['Trp53', 'Brca1', 'ERBB2', 'p21', 'Rb1', 'PTEN',
                    'Myc', 'Kras', 'Braf', 'Akt1', 'Bcl2', 'mTOR',
                    'CD8a', 'Foxp3', 'PD-1', 'PD-L1', 'CTLA4'];
            } else {
                testGenes = ['TP53', 'BRCA1', 'BRCA2', 'MYC', 'KRAS', 'EGFR', 'PTEN',
                    'RB1', 'APC', 'CDKN2A', 'NOTCH1', 'PIK3CA', 'BRAF',
                    'ATM', 'ERBB2', 'CDK4', 'MDM2', 'NRAS', 'TSC1', 'TSC2',
                    'BCR', 'ABL1'];
            }
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

        // Stats sub-tabs (file upload vs manual entry)
        document.querySelectorAll('.stats-sub-tab').forEach(tab => {
            tab.addEventListener('click', () => {
                document.querySelectorAll('.stats-sub-tab').forEach(t => {
                    t.classList.remove('active');
                    t.style.background = 'transparent';
                    t.style.color = '#6b7280';
                });
                document.querySelectorAll('.stats-sub-panel').forEach(p => {
                    p.classList.remove('active');
                    p.style.display = 'none';
                });
                tab.classList.add('active');
                tab.style.background = '#5a9f4a';
                tab.style.color = 'white';
                const panel = document.getElementById('stats-' + tab.dataset.statsInput);
                panel.classList.add('active');
                panel.style.display = 'block';
            });
        });

        // Manual stats entry - now auto-loaded on Run, no button needed

        // Prevent keyboard events from propagating to network when typing in textareas/inputs
        // Use event delegation to handle all inputs including dynamically added ones
        document.addEventListener('keydown', (e) => {
            const el = e.target;
            if (el.tagName === 'TEXTAREA' || el.tagName === 'INPUT') {
                // Stop propagation to prevent vis-network from capturing arrow keys etc.
                e.stopPropagation();

                // Handle Tab key in textareas - insert tab character instead of changing focus
                if (e.key === 'Tab' && el.tagName === 'TEXTAREA') {
                    e.preventDefault();
                    const start = el.selectionStart;
                    const end = el.selectionEnd;
                    el.value = el.value.substring(0, start) + '\t' + el.value.substring(end);
                    el.selectionStart = el.selectionEnd = start + 1;
                }
            }
        }, true); // Use capture phase to intercept before vis-network

        // Run analysis
        document.getElementById('runAnalysis').addEventListener('click', () => this.runAnalysis());

        // Reset App button
        document.getElementById('resetAppBtn')?.addEventListener('click', () => location.reload());

        // Reset excluded tissues button
        document.getElementById('resetExcludedTissuesBtn')?.addEventListener('click', () => {
            document.querySelectorAll('#tissueExcludeList input[type="checkbox"]').forEach(cb => cb.checked = false);
            this.excludedTissues = new Set();
        });

        // Analysis mode change
        document.querySelectorAll('input[name="analysisMode"]').forEach(radio => {
            radio.addEventListener('change', () => this.updateAnalysisModeUI());
        });

        // Tissue breakdown button (mutations)
        document.getElementById('tissueBreakdownBtn').addEventListener('click', () => this.showTissueBreakdownPopup('mutation'));
        document.getElementById('mutationHotspotSelect').addEventListener('change', () => {
            const btn = document.getElementById('tissueBreakdownBtn');
            btn.style.display = document.getElementById('mutationHotspotSelect').value ? 'inline-block' : 'none';
        });

        // Tissue breakdown button (translocations)
        document.getElementById('translocationTissueBreakdownBtn')?.addEventListener('click', () => this.showTissueBreakdownPopup('translocation'));
        document.getElementById('translocationHotspotSelect')?.addEventListener('input', () => {
            const btn = document.getElementById('translocationTissueBreakdownBtn');
            const val = document.getElementById('translocationHotspotSelect').value.trim();
            if (btn) btn.style.display = (val && this.translocations?.geneData?.[val]) ? 'inline-block' : 'none';
        });

        // General tissue breakdown button (for lineage filter area)
        document.getElementById('generalTissueBreakdownBtn')?.addEventListener('click', () => this.showGeneralTissueBreakdownPopup());

        // Oncoprint button
        document.getElementById('oncoprintBtn')?.addEventListener('click', () => this.showOncoprint());

        // Mutation analysis sub-type selector
        document.querySelectorAll('input[name="mutAnalysisType"]').forEach(radio => {
            radio.addEventListener('change', () => this.updateMutAnalysisTypeUI());
        });

        // Parameter translocation filter
        document.getElementById('paramTranslocationGene')?.addEventListener('input', () => {
            const val = document.getElementById('paramTranslocationGene').value.trim();
            if (val === '' || this.translocations?.geneData?.[val]) {
                this.updateParamTranslocationLevelCounts();
            }
        });
        document.getElementById('paramTranslocationLevel')?.addEventListener('change', () => {
            // Counts stay the same; filter applies on Run
        });

        // Mutation results search
        document.getElementById('mutationSearch').addEventListener('input', (e) => {
            this.filterMutationTable(e.target.value);
        });

        // Download mutation results
        document.getElementById('downloadMutationResults').addEventListener('click', () => {
            this.downloadMutationResults();
        });
        document.getElementById('exportMutTablePNG')?.addEventListener('click', () => this._exportMutationTable('png'));
        document.getElementById('exportMutTableSVG')?.addEventListener('click', () => this._exportMutationTable('svg'));

        // Synonyms search
        document.getElementById('synonymsSearch').addEventListener('input', (e) => {
            this.filterSynonymsTable(e.target.value);
        });

        // Download synonyms results
        document.getElementById('downloadSynonyms').addEventListener('click', () => {
            this.downloadSynonymsCSV();
        });

        // Gene effect distribution modal
        document.getElementById('closeGeneEffect').addEventListener('click', () => {
            document.getElementById('geneEffectModal').style.display = 'none';
            this._geHighlightCellLine = null;
        });
        document.getElementById('geneEffectModal').addEventListener('click', (e) => {
            if (e.target.id === 'geneEffectModal') {
                document.getElementById('geneEffectModal').style.display = 'none';
                this._geHighlightCellLine = null;
            }
        });
        document.getElementById('downloadGeneEffectPNG').addEventListener('click', () => this.downloadGeneEffectPNG());
        document.getElementById('downloadGeneEffectSVG').addEventListener('click', () => this.downloadGeneEffectSVG());

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
        document.getElementById('restoreAllNodes').addEventListener('click', () => this.showHiddenNodes());
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

        // Physics toggle, layout change, and remove mode buttons
        document.getElementById('togglePhysics').addEventListener('click', () => this.togglePhysics());
        document.getElementById('changeLayout').addEventListener('click', () => this.changeNetworkLayout());
        document.getElementById('toggleRemoveMode').addEventListener('click', () => this.toggleRemoveMode());
        document.getElementById('toggleSelectMode').addEventListener('click', () => this.toggleSelectMode());
        document.getElementById('clearSelectedNodes').addEventListener('click', () => this.clearSelectedNodes());
        document.getElementById('showUncorrelatedGenes').addEventListener('change', () => {
            if (this.results) {
                this.displayNetwork();
                // Re-apply colors after network rebuild if color by GE is active
                if (document.getElementById('colorByGeneEffect').checked) {
                    setTimeout(() => this.updateNetworkColors(), 100);
                }
            }
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

        // Global keyboard handler for closing modals with Enter or Escape
        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape' || e.key === 'Enter') {
                // Close gene effect modal first (it sits on top of CLB at z-index 1200)
                if (e.key === 'Escape') {
                    const geneEffectModal = document.getElementById('geneEffectModal');
                    if (geneEffectModal && geneEffectModal.style.display !== 'none') {
                        geneEffectModal.style.display = 'none';
                        this._geHighlightCellLine = null;
                        return;
                    }
                }
                // Close cell line browser if open (Escape only)
                if (e.key === 'Escape') {
                    const clbModal = document.getElementById('cellLineBrowserModal');
                    if (clbModal && clbModal.style.display !== 'none') {
                        this.closeCellLineBrowser();
                        return;
                    }
                }
                // Close inspect modal if open
                const inspectModal = document.getElementById('inspectModal');
                if (inspectModal && inspectModal.classList.contains('active')) {
                    this.closeInspectModal();
                    return;
                }
                // Close infographic modal if open
                const infographicModal = document.getElementById('infographicModal');
                if (infographicModal && infographicModal.style.display !== 'none') {
                    infographicModal.style.display = 'none';
                    return;
                }
            }
        });

        // Inspect controls
        document.getElementById('resetAxes').addEventListener('click', () => this.resetInspectAxes());
        document.getElementById('clearHighlights').addEventListener('click', () => {
            document.getElementById('scatterCellSearch').value = '';
            this.clickedCells.clear();
            this._userLabelPositions.clear();
            document.getElementById('colorByCategory').value = '';
            this._styleActiveFilters();
            this.updateInspectPlot();
        });

        document.getElementById('resetAllFilters').addEventListener('click', () => this.resetAllInspectFilters());

        ['scatterXmin', 'scatterXmax', 'scatterYmin', 'scatterYmax'].forEach(id => {
            const el = document.getElementById(id);
            el.addEventListener('change', () => {
                const val = el.value.trim();
                if (val === '' || val === '-' || val === '.' || val === '-.') return;
                this.updateInspectPlot();
            });
            el.addEventListener('keydown', (e) => { if (e.key === 'Enter') this.updateInspectPlot(); });
        });

        document.getElementById('showCorrelationLine').addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('showZeroLines')?.addEventListener('change', () => this.updateInspectPlot());

        document.getElementById('scatterCellSearch').addEventListener('input', () => this.updateInspectPlot());
        document.getElementById('colorByCategory').addEventListener('change', () => {
            this._styleActiveFilters();
            this.updateInspectPlot();
        });
        document.getElementById('scatterCancerFilter').addEventListener('change', () => {
            this.updateScatterSubtypeFilter();
            // Show/hide "Color by subtype" option based on whether a tissue is selected
            const subtypeOpt = document.getElementById('colorBySubtypeOption');
            if (subtypeOpt) {
                subtypeOpt.style.display = document.getElementById('scatterCancerFilter').value ? '' : 'none';
                // Reset color-by if subtype was selected but tissue filter cleared
                if (!document.getElementById('scatterCancerFilter').value && document.getElementById('colorByCategory').value === 'subtype') {
                    document.getElementById('colorByCategory').value = '';
                }
            }
            this._styleActiveFilters();
            this.updateScatterHotspotFilterCounts();
            this.updateInspectPlot();
        });
        document.getElementById('scatterSubtypeFilter').addEventListener('change', () => {
            this._styleActiveFilters();
            this.updateScatterHotspotFilterCounts();
            this.updateInspectPlot();
        });
        document.getElementById('mutationFilterGene').addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('mutationFilterLevel').addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('hotspotGene').addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('hotspotMode').addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('translocationGene')?.addEventListener('input', () => {
            const val = document.getElementById('translocationGene').value.trim();
            if (val === '' || this.translocations?.geneData?.[val]) this.updateInspectPlot();
        });
        document.getElementById('translocationMode')?.addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('translocationFilterGene')?.addEventListener('input', () => {
            const val = document.getElementById('translocationFilterGene').value.trim();
            if (val === '' || this.translocations?.geneData?.[val]) this.updateInspectPlot();
        });
        document.getElementById('translocationFilterLevel')?.addEventListener('change', () => this.updateInspectPlot());

        document.getElementById('downloadScatterPNG').addEventListener('click', () => this.downloadScatterPNG());
        document.getElementById('downloadScatterSVG').addEventListener('click', () => this.downloadScatterSVG());
        document.getElementById('downloadScatterCSV').addEventListener('click', () => this.downloadScatterCSV());
        document.getElementById('scatterTextSettingsBtn')?.addEventListener('click', () => this.openTextSettings('scatterPlot'));
        document.getElementById('geTextSettingsBtn')?.addEventListener('click', () => this.openTextSettings('geneEffectPlot'));
        document.getElementById('networkAaBtn')?.addEventListener('click', () => this.openNetworkTextSettings());
        document.getElementById('networkNodeBorder')?.addEventListener('change', (e) => this.toggleNetworkBorder(e.target.checked));
        this._initTextSettingsDrag();
        document.getElementById('restoreFromSvgInput')?.addEventListener('change', (e) => {
            if (e.target.files[0]) this.restoreFromExport(e.target.files[0]);
            e.target.value = '';
        });
        document.getElementById('restoreFromExportInput')?.addEventListener('change', (e) => {
            if (e.target.files[0]) this.restoreFromExport(e.target.files[0]);
            e.target.value = '';
        });
        document.getElementById('downloadTissuePNG').addEventListener('click', () => this.downloadTissueChartPNG());
        document.getElementById('downloadTissueSVG').addEventListener('click', () => this.downloadTissueChartSVG());
        document.getElementById('downloadTissueCSV').addEventListener('click', () => this.downloadTissueTableCSV());
        document.getElementById('scatterFontSize')?.addEventListener('change', () => this.updateInspectPlot());
        document.getElementById('compareAllMutationsBtn')?.addEventListener('click', () => this.showCompareAllMutations());
        document.getElementById('compareAllTranslocationsBtn')?.addEventListener('click', () => this.showCompareAllTranslocations());
        document.getElementById('compareAllCancerTypesBtn')?.addEventListener('click', () => this.showCompareAllCancerTypes());
        document.getElementById('updateInspectGenes')?.addEventListener('click', () => this.updateInspectGenes());

        // Plot size controls
        ['plotWidth', 'plotHeight'].forEach(id => {
            const el = document.getElementById(id);
            el.addEventListener('change', () => this.updateInspectPlot());
            el.addEventListener('keydown', (e) => { if (e.key === 'Enter') this.updateInspectPlot(); });
        });

        // Table header sorting (Ctrl/Cmd+click to copy column)
        document.querySelectorAll('.data-table th[data-sort]').forEach(th => {
            th.addEventListener('click', (e) => {
                if (e.ctrlKey || e.metaKey) {
                    const colIndex = Array.from(th.parentNode.children).indexOf(th);
                    this.copyColumnToClipboard(th.closest('table'), colIndex);
                } else {
                    this.sortTable(th);
                }
            });
            th.title = 'Click to sort. Ctrl/Cmd+click to copy column.';
        });

        // Copy genes buttons
        document.getElementById('copyCorrelationGenes')?.addEventListener('click', () => this.copyGeneColumn('correlationsTable'));
        document.getElementById('copyClustersGenes')?.addEventListener('click', () => this.copyGeneColumn('clustersTable'));
        document.getElementById('copyMutationGenes')?.addEventListener('click', () => this.copyGeneColumn('mutationTable'));


        // Table filter toggles
        document.getElementById('filterCorrelationsToggle')?.addEventListener('click', () => this.toggleTableFilters('correlationsTable'));
        document.getElementById('filterClustersToggle')?.addEventListener('click', () => this.toggleTableFilters('clustersTable'));
        document.getElementById('filterMutationToggle')?.addEventListener('click', () => this.toggleTableFilters('mutationTable'));

        // Enrichr buttons
        document.querySelectorAll('.enrichrBtn').forEach(btn => {
            btn.addEventListener('click', () => this.openEnrichr(btn.dataset.source));
        });
        document.getElementById('mutEnrichrBtn')?.addEventListener('click', () => this.openEnrichr('mutations'));
        document.getElementById('enrichrCloseBtn')?.addEventListener('click', () => {
            document.getElementById('enrichrModal').style.display = 'none';
        });
        document.getElementById('enrichrModal')?.addEventListener('click', (e) => {
            if (e.target.id === 'enrichrModal') {
                document.getElementById('enrichrModal').style.display = 'none';
            }
        });
        document.getElementById('enrichrDownloadBtn')?.addEventListener('click', () => this.downloadEnrichrCSV());

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

        // Gene Effect modal controls
        document.getElementById('geneEffectSearchBtn')?.addEventListener('click', () => {
            const gene = document.getElementById('geneEffectSearch').value.trim().toUpperCase();
            if (!gene) return;
            if (this.geneEffectViewMode === 'mutation' && this.mutationResults) {
                this.showGeneEffectDistribution(gene);
            } else {
                this.showGeneEffectAnalysis(gene, this.currentGEView || 'tissue');
            }
        });
        document.getElementById('geneEffectSearch')?.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') {
                const gene = e.target.value.trim().toUpperCase();
                if (!gene) return;
                if (this.geneEffectViewMode === 'mutation' && this.mutationResults) {
                    this.showGeneEffectDistribution(gene);
                } else {
                    this.showGeneEffectAnalysis(gene, this.currentGEView || 'tissue');
                }
            }
        });
        document.getElementById('geViewTissue')?.addEventListener('click', () => {
            if (this.geneEffectViewMode === 'mutation') {
                // Switch from mutation inspect to full gene effect analysis
                const gene = document.getElementById('geneEffectSearch').value.trim().toUpperCase() || this.currentGeneEffectGene;
                if (gene) this.openGeneEffectModal(gene, 'tissue');
            } else {
                this.switchGeneEffectView('tissue');
            }
        });
        document.getElementById('geViewHotspot')?.addEventListener('click', () => {
            if (this.geneEffectViewMode === 'mutation') {
                const gene = document.getElementById('geneEffectSearch').value.trim().toUpperCase() || this.currentGeneEffectGene;
                if (gene) this.openGeneEffectModal(gene, 'hotspot');
            } else {
                this.switchGeneEffectView('hotspot');
            }
        });
        document.getElementById('geTissueFilter')?.addEventListener('change', () => {
            this.updateGeSubtypeFilter();
            if (this.geneEffectViewMode === 'mutation' && this.mutationResults && this.currentGeneEffectGene) {
                this.showGeneEffectDistribution(this.currentGeneEffectGene);
            } else {
                this.switchGeneEffectView(this.currentGEView || 'tissue');
            }
        });
        // Inspect-level subtype filter
        document.getElementById('geSubtypeFilter')?.addEventListener('change', () => {
            if (this.geneEffectViewMode === 'mutation' && this.currentGeneEffectGene) {
                this.showGeneEffectDistribution(this.currentGeneEffectGene);
            }
        });
        // Inspect-level hotspot filter
        document.getElementById('geHotspotFilter')?.addEventListener('change', () => {
            const warn = document.getElementById('geHotspotBiasWarning');
            if (warn) warn.style.display = document.getElementById('geHotspotFilter').value ? 'inline' : 'none';
            if (this.geneEffectViewMode === 'mutation' && this.currentGeneEffectGene) {
                this.showGeneEffectDistribution(this.currentGeneEffectGene);
            } else if (this.geneEffectViewMode === 'geneEffect') {
                this.switchGeneEffectView(this.currentGEView || 'tissue');
            }
        });
        // Inspect-level fusion filter
        document.getElementById('geFusionFilter')?.addEventListener('change', () => {
            if (this.geneEffectViewMode === 'mutation' && this.currentGeneEffectGene) {
                this.showGeneEffectDistribution(this.currentGeneEffectGene);
            } else if (this.geneEffectViewMode === 'geneEffect') {
                this.switchGeneEffectView(this.currentGEView || 'tissue');
            }
        });
        // Hotspot gene selector (Y axis mutation gene)
        document.getElementById('geHotspotGeneSelect')?.addEventListener('change', () => {
            if (this.geneEffectViewMode === 'mutation' && this.currentGeneEffectGene && this.mutationResults) {
                const newHotspot = document.getElementById('geHotspotGeneSelect').value;
                const isTransloc = this.mutationResults.isTranslocation;
                const hasData = isTransloc
                    ? this.translocations?.geneData?.[newHotspot]
                    : this.mutations?.geneData?.[newHotspot];
                if (newHotspot && hasData) {
                    this.mutationResults.hotspotGene = newHotspot;
                    this.showGeneEffectDistribution(this.currentGeneEffectGene);
                }
            }
        });
        // Cell line search in gene effect modal
        let geCellLineSearchTimer;
        document.getElementById('geCellLineSearch')?.addEventListener('input', () => {
            clearTimeout(geCellLineSearchTimer);
            geCellLineSearchTimer = setTimeout(() => {
                if (this.currentGeneEffect) {
                    this.switchGeneEffectView(this.currentGEView || 'tissue');
                }
            }, 200);
        });
        // Inline compare buttons
        document.getElementById('geCompareByTissueBtn')?.addEventListener('click', () => this.showInlineCompareByTissue());
        document.getElementById('geCompareByHotspotBtn')?.addEventListener('click', () => this.showInlineCompareByHotspot());
        document.getElementById('geCompareByTranslocationBtn')?.addEventListener('click', () => this.showInlineCompareByTranslocation());
        document.getElementById('geInlineCompareClose')?.addEventListener('click', () => {
            document.getElementById('geInlineCompareTable').style.display = 'none';
        });
        document.getElementById('geResetFiltersBtn')?.addEventListener('click', () => {
            document.getElementById('geTissueFilter').value = '';
            document.getElementById('geSubtypeFilter').value = '';
            const geSubEl = document.getElementById('geSubtypeFilter');
            if (geSubEl) geSubEl.style.display = 'none';
            document.getElementById('geHotspotFilter').value = '';
            document.getElementById('geFusionFilter').value = '';
            // Clear analysis-level filters so inspect truly shows ALL cell lines
            if (this.mutationResults) {
                this.mutationResults.lineageFilter = '';
                this.mutationResults.subLineageFilter = '';
                this.mutationResults.excludedTissues = new Set();
                this.mutationResults.additionalHotspot = '';
                this.mutationResults.additionalHotspotLevel = 'all';
            }
            if (this.geneEffectViewMode === 'mutation' && this.currentGeneEffectGene) {
                this.showGeneEffectDistribution(this.currentGeneEffectGene);
            }
        });
        // Full-screen compare modal buttons
        document.getElementById('mutCompareByTissueBtn')?.addEventListener('click', () => this.showMutationCompareByTissue());
        document.getElementById('mutCompareByHotspotBtn')?.addEventListener('click', () => this.showMutationCompareByHotspot());
        document.getElementById('mutCompareByFusionBtn')?.addEventListener('click', () => this.showMutationCompareByFusion());
        document.getElementById('mutCompareModalClose')?.addEventListener('click', () => {
            document.getElementById('mutCompareModal').style.display = 'none';
        });
        document.getElementById('mutCompareMinN')?.addEventListener('change', () => {
            if (this._compareModalMode && document.getElementById('mutCompareModal').style.display !== 'none') {
                this.renderCompareModal();
            }
        });
        document.getElementById('mutCompareExportCSV')?.addEventListener('click', () => this.downloadCompareCSV());
        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape' && document.getElementById('mutCompareModal')?.style.display !== 'none') {
                document.getElementById('mutCompareModal').style.display = 'none';
            }
        });
        document.getElementById('geShowAllBtn')?.addEventListener('click', () => {
            this.showAllGeneEffect();
        });
        document.getElementById('geTableSearch')?.addEventListener('input', (e) => {
            this.filterGETable(e.target.value);
        });
        document.getElementById('downloadGeneEffectPNG')?.addEventListener('click', () => this.downloadGeneEffectChartPNG());
        document.getElementById('downloadGeneEffectSVG')?.addEventListener('click', () => this.downloadGeneEffectChartSVG());
        document.getElementById('downloadGETableCSV')?.addEventListener('click', () => this.downloadGETableCSV());
        document.getElementById('downloadGECellLineCSV')?.addEventListener('click', () => this.downloadGECellLineCSV());

        // Chart width control — works for both inspect and detailed views
        document.getElementById('geAspectRatio')?.addEventListener('input', (e) => {
            const ratio = parseFloat(e.target.value);
            document.getElementById('geAspectRatioValue').textContent = ratio.toFixed(1);
            this.geChartWidthRatio = ratio;
            const container = document.getElementById('geChartContainer');
            if (container) {
                container.style.flex = `0 0 ${Math.round(ratio * 55)}%`;
            }
            // Re-render the detailed view with new width if in that mode
            if (this.geDetailedView) {
                this.showGEDetailedView(this.geDetailedView.group, this.geDetailedView.mode);
            } else {
                // Trigger Plotly resize for the inspect plot
                const plotEl = document.getElementById('geneEffectPlot');
                if (plotEl?.data?.length) Plotly.Plots.resize(plotEl);
            }
        });

        // Chart height control for inspect plot
        document.getElementById('geHeightRatio')?.addEventListener('input', (e) => {
            const ratio = parseFloat(e.target.value);
            document.getElementById('geHeightRatioValue').textContent = ratio.toFixed(1);
            this.geChartHeightRatio = ratio;
            const plotEl = document.getElementById('geneEffectPlot');
            if (plotEl?.data?.length) {
                const baseHeight = 400;
                Plotly.relayout(plotEl, { height: Math.round(baseHeight * ratio) });
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
            // Find fuzzy suggestions for not-found genes
            const suggestions = this._findGeneSuggestions(notFound);
            let notFoundHtml = notFound.slice(0, 10).map(g => {
                const sugg = suggestions.get(g);
                if (sugg && sugg.length > 0) {
                    const links = sugg.slice(0, 3).map(s =>
                        `<a href="#" style="color: #0066cc;" onclick="app.replaceGeneInTextarea('${g}', '${s}'); return false;">${s}</a>`
                    ).join(', ');
                    return `${g} → ${links}?`;
                }
                return g;
            }).join(', ');
            if (notFound.length > 10) notFoundHtml += ` (+${notFound.length - 10} more)`;

            display.innerHTML = `<div class="status-box status-warning">
                <strong>${found.length} found</strong>, <strong>${notFound.length} not found</strong><br>
                <span>Not found: ${notFoundHtml}</span>
            </div>`;
            synonymBtn.style.display = 'block';
            this.genesNotFound = notFound;
        }
    }

    _findGeneSuggestions(notFoundGenes) {
        const suggestions = new Map();
        if (!this.geneIndex || this.geneIndex.size === 0) return suggestions;

        const allGenes = Array.from(this.geneIndex.keys());

        for (const query of notFoundGenes) {
            const upper = query.toUpperCase();
            const matches = [];

            // 1. Check prefix matches
            for (const gene of allGenes) {
                if (gene.startsWith(upper) || upper.startsWith(gene)) {
                    matches.push(gene);
                    if (matches.length >= 3) break;
                }
            }

            // 2. Check edit distance ≤ 2 (only if few/no prefix matches found)
            if (matches.length < 3) {
                for (const gene of allGenes) {
                    if (matches.includes(gene)) continue;
                    if (Math.abs(gene.length - upper.length) > 2) continue;
                    const dist = this._editDistance(upper, gene);
                    if (dist <= 2) {
                        matches.push(gene);
                        if (matches.length >= 3) break;
                    }
                }
            }

            if (matches.length > 0) {
                suggestions.set(query, matches);
            }
        }
        return suggestions;
    }

    _editDistance(a, b) {
        if (a.length === 0) return b.length;
        if (b.length === 0) return a.length;
        // Quick bail for very different lengths
        if (Math.abs(a.length - b.length) > 2) return 3;

        const m = a.length, n = b.length;
        let prev = new Array(n + 1);
        let curr = new Array(n + 1);
        for (let j = 0; j <= n; j++) prev[j] = j;

        for (let i = 1; i <= m; i++) {
            curr[0] = i;
            for (let j = 1; j <= n; j++) {
                curr[j] = a[i - 1] === b[j - 1]
                    ? prev[j - 1]
                    : 1 + Math.min(prev[j - 1], prev[j], curr[j - 1]);
            }
            [prev, curr] = [curr, prev];
        }
        return prev[n];
    }

    replaceGeneInTextarea(oldGene, newGene) {
        const textarea = document.getElementById('geneTextarea');
        const lines = textarea.value.split('\n');
        const updatedLines = lines.map(line => {
            if (line.trim().toUpperCase() === oldGene.toUpperCase()) {
                return newGene;
            }
            return line;
        });
        textarea.value = updatedLines.join('\n');
        this.updateGeneCount();
    }

    async findSynonymsForMissingGenes() {
        if (!this.genesNotFound || this.genesNotFound.length === 0) return;

        // Lazy load synonyms if not already loaded
        if (!this.synonymLookup) {
            const btn = document.getElementById('findSynonyms');
            const originalText = btn.textContent;
            btn.textContent = 'Loading synonyms...';
            btn.disabled = true;
            try {
                const synonymsRes = await fetch('web_data/synonyms.json');
                this.synonymLookup = await synonymsRes.json();
            } catch (e) {
                console.error('Failed to load synonyms:', e);
                this.synonymLookup = {};
            }
            btn.textContent = originalText;
            btn.disabled = false;
        }

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

            // Update geneTextarea
            const textarea = document.getElementById('geneTextarea');
            let text = textarea.value;

            replacements.forEach(r => {
                const regex = new RegExp(`\\b${r.original}\\b`, 'gi');
                text = text.replace(regex, r.replacement);
            });

            textarea.value = text;

            // Also update manualStatsTextarea if it has content
            const manualTextarea = document.getElementById('manualStatsTextarea');
            if (manualTextarea && manualTextarea.value.trim()) {
                let manualText = manualTextarea.value;
                replacements.forEach(r => {
                    const regex = new RegExp(`\\b${r.original}\\b`, 'gi');
                    manualText = manualText.replace(regex, r.replacement);
                });
                manualTextarea.value = manualText;
            }

            // Also update geneStats if loaded
            if (this.geneStats && this.geneStats.size > 0) {
                replacements.forEach(r => {
                    const oldKey = r.original.toUpperCase();
                    const newKey = r.replacement.toUpperCase();
                    if (this.geneStats.has(oldKey)) {
                        const stats = this.geneStats.get(oldKey);
                        stats.gene = newKey;
                        this.geneStats.delete(oldKey);
                        this.geneStats.set(newKey, stats);
                    }
                });
            }

            this.updateGeneCount();

            const msg = replacements.map(r => `${r.original} → ${r.replacement} [${r.source}]`).join('\n');
            alert(`Replaced ${replacements.length} gene(s):\n${msg}`);
        } else {
            alert('No synonyms or orthologs found for the missing genes');
        }
    }

    // Resolve a single gene to its canonical name using synonym/ortholog lookup
    async resolveGeneSynonym(gene) {
        const upperGene = gene.toUpperCase();

        // If already in dataset, return as-is
        if (this.geneIndex.has(upperGene)) {
            return { gene: upperGene, source: 'direct' };
        }

        // Load synonyms if not yet loaded
        if (!this.synonymLookup) {
            try {
                const synonymsRes = await fetch('web_data/synonyms.json');
                this.synonymLookup = await synonymsRes.json();
            } catch (e) {
                console.warn('Failed to load synonyms:', e);
                this.synonymLookup = {};
            }
        }

        // Check local synonym lookup (loaded from synonyms.json)
        if (this.synonymLookup) {
            const match = this.synonymLookup[upperGene];
            if (match && this.geneIndex.has(match.toUpperCase())) {
                return { gene: match.toUpperCase(), source: 'synonym' };
            }
        }

        // Check orthologs
        if (this.orthologs && this.orthologs[upperGene]) {
            const humanGene = this.orthologs[upperGene];
            if (this.geneIndex.has(humanGene.toUpperCase())) {
                return { gene: humanGene.toUpperCase(), source: 'ortholog' };
            }
        }

        // Check previously used synonyms from the session
        if (this.synonymsUsed && this.synonymsUsed.length > 0) {
            const found = this.synonymsUsed.find(s => s.original.toUpperCase() === upperGene);
            if (found && this.geneIndex.has(found.replacement.toUpperCase())) {
                return { gene: found.replacement.toUpperCase(), source: found.source };
            }
        }

        return null;
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

    loadManualStats() {
        const text = document.getElementById('manualStatsTextarea').value.trim();
        if (!text) {
            alert('Please enter gene data in the text area.');
            return;
        }

        const lines = text.split('\n').filter(l => l.trim());
        if (lines.length === 0) return;

        // Check if first line is a header
        const firstLine = lines[0].toLowerCase();
        const headerKeywords = ['gene', 'symbol', 'lfc', 'logfc', 'log2fc', 'fdr', 'padj', 'pvalue', 'p-value'];
        let hasHeader = headerKeywords.some(kw => firstLine.includes(kw));

        // Detect separator from data lines (not header) - check which separator gives consistent column count
        const dataLines = hasHeader ? lines.slice(1) : lines;
        const testLine = dataLines.find(l => l.trim()) || lines[0];

        let separator = '\t';
        if (testLine.includes('\t') && testLine.split('\t').length >= 2) {
            separator = '\t';
        } else if (testLine.includes(',') && testLine.split(',').length >= 2) {
            separator = ',';
        } else if (testLine.includes(';') && testLine.split(';').length >= 2) {
            separator = ';';
        }

        // Parse header to detect columns (use tab for header if it has tabs, otherwise same as data)
        const headerSep = lines[0].includes('\t') ? '\t' : separator;
        const firstParts = lines[0].split(headerSep).map(p => p.trim().toLowerCase());

        // Determine column indices
        let geneCol = 0, lfcCol = -1, fdrCol = -1;
        if (hasHeader) {
            firstParts.forEach((col, i) => {
                if (['gene', 'symbol', 'geneid', 'gene_symbol', 'genename'].includes(col)) geneCol = i;
                else if (['lfc', 'logfc', 'log2fc', 'log2foldchange', 'log_fc'].includes(col)) lfcCol = i;
                else if (['fdr', 'padj', 'adj.p.val', 'adjpval', 'q', 'qvalue', 'q.value'].includes(col)) fdrCol = i;
            });
        } else {
            // Assume: Gene, LFC, FDR order if 3 columns; Gene only if 1 column
            const testParts = testLine.split(separator);
            if (testParts.length >= 2) lfcCol = 1;
            if (testParts.length >= 3) fdrCol = 2;
        }

        const genes = [];
        this.geneStats = new Map();
        const startLine = hasHeader ? 1 : 0;

        for (let i = startLine; i < lines.length; i++) {
            const parts = lines[i].split(separator).map(p => p.trim());
            const gene = parts[geneCol]?.toUpperCase();
            if (!gene) continue;

            genes.push(gene);

            const stats = { gene };
            if (lfcCol >= 0 && parts[lfcCol]) {
                stats.lfc = parseFloat(parts[lfcCol]) || null;
            }
            if (fdrCol >= 0 && parts[fdrCol]) {
                stats.fdr = parseFloat(parts[fdrCol]) || null;
            }
            this.geneStats.set(gene, stats);
        }

        // Update textarea in paste panel
        document.getElementById('geneTextarea').value = genes.join('\n');
        this.updateGeneCount();

        // Show stats controls if we have statistics
        const hasLfc = lfcCol >= 0;
        const hasFdr = fdrCol >= 0;
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
        // Test data with LFC and FDR values (22 genes — matches simple test gene list)
        const testData = [
            { gene: 'TP53', lfc: 1.8, fdr: 0.001 },
            { gene: 'BRCA1', lfc: 0.3, fdr: 0.15 },
            { gene: 'BRCA2', lfc: -0.5, fdr: 0.11 },
            { gene: 'MYC', lfc: -1.5, fdr: 0.02 },
            { gene: 'KRAS', lfc: 2.1, fdr: 0.005 },
            { gene: 'EGFR', lfc: -0.8, fdr: 0.08 },
            { gene: 'PTEN', lfc: 1.2, fdr: 0.03 },
            { gene: 'RB1', lfc: -0.6, fdr: 0.12 },
            { gene: 'APC', lfc: 0.7, fdr: 0.09 },
            { gene: 'CDKN2A', lfc: 1.4, fdr: 0.007 },
            { gene: 'NOTCH1', lfc: -1.1, fdr: 0.03 },
            { gene: 'PIK3CA', lfc: 0.9, fdr: 0.04 },
            { gene: 'BRAF', lfc: -1.9, fdr: 0.002 },
            { gene: 'ATM', lfc: 0.7, fdr: 0.06 },
            { gene: 'ERBB2', lfc: 1.6, fdr: 0.008 },
            { gene: 'CDK4', lfc: -0.3, fdr: 0.22 },
            { gene: 'MDM2', lfc: -2.3, fdr: 0.0001 },
            { gene: 'NRAS', lfc: 1.3, fdr: 0.015 },
            { gene: 'TSC1', lfc: 1.0, fdr: 0.02 },
            { gene: 'TSC2', lfc: 0.9, fdr: 0.03 },
            { gene: 'BCR', lfc: -0.7, fdr: 0.14 },
            { gene: 'ABL1', lfc: 0.6, fdr: 0.07 }
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
            ['BRCA1', '0.3', '0.15'],
            ['BRCA2', '-0.5', '0.11'],
            ['MYC', '-1.5', '0.02'],
            ['KRAS', '2.1', '0.005'],
            ['EGFR', '-0.8', '0.08'],
            ['PTEN', '1.2', '0.03'],
            ['RB1', '-0.6', '0.12'],
            ['APC', '0.7', '0.09'],
            ['CDKN2A', '1.4', '0.007'],
            ['NOTCH1', '-1.1', '0.03'],
            ['PIK3CA', '0.9', '0.04'],
            ['BRAF', '-1.9', '0.002'],
            ['ATM', '0.7', '0.06'],
            ['ERBB2', '1.6', '0.008'],
            ['CDK4', '-0.3', '0.22'],
            ['MDM2', '-2.3', '0.0001'],
            ['NRAS', '1.3', '0.015'],
            ['TSC1', '1.0', '0.02'],
            ['TSC2', '0.9', '0.03'],
            ['BCR', '-0.7', '0.14'],
            ['ABL1', '0.6', '0.07']
        ];

        const csv = sampleData.map(row => row.join(',')).join('\n');
        this.downloadFile(csv, 'sample_genes_with_stats.csv', 'text/csv');
    }

    getGeneList() {
        const text = document.getElementById('geneTextarea').value;
        const genes = text.split(/\s+/)
            .map(g => g.toUpperCase().trim())
            .filter(g => g !== '' && this.geneIndex.has(g));
        // Remove duplicates
        return [...new Set(genes)];
    }

    getFilteredCellLineIndices() {
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;
        const hotspotGene = document.getElementById('paramHotspotGene').value;
        const hotspotLevel = document.getElementById('paramHotspotLevel').value;
        const transGene = document.getElementById('paramTranslocationGene').value;
        const transLevel = document.getElementById('paramTranslocationLevel').value;

        // Get mutation data for hotspot filter
        let mutationData = null;
        if (hotspotGene && this.mutations?.geneData?.[hotspotGene]) {
            mutationData = this.mutations.geneData[hotspotGene].mutations;
        }

        // Get translocation data for fusion filter
        let transData = null;
        if (transGene && this.translocations?.geneData?.[transGene]) {
            transData = this.translocations.geneData[transGene].translocations;
        }

        const indices = [];
        this.metadata.cellLines.forEach((cellLine, idx) => {
            // Check excluded tissues
            if (this.excludedTissues && this.excludedTissues.size > 0) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                if (lineage && this.excludedTissues.has(lineage)) {
                    return;
                }
            }

            // Check lineage filter
            if (lineageFilter) {
                if (!this.cellLineMetadata.lineage ||
                    this.cellLineMetadata.lineage[cellLine] !== lineageFilter) {
                    return;
                }
            }

            // Check sublineage filter
            if (subLineageFilter) {
                if (!this.cellLineMetadata.primaryDisease ||
                    this.cellLineMetadata.primaryDisease[cellLine] !== subLineageFilter) {
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

            // Check translocation filter
            if (transData && transLevel !== 'all') {
                const tLevel = transData[cellLine] || 0;
                if (transLevel === '0' && tLevel !== 0) return;
                if (transLevel === '1' && tLevel !== 1) return;
                if (transLevel === '2' && tLevel < 2) return;
                if (transLevel === '1+2' && tLevel < 1) return;
            }

            // Check oncoprint multi-gene filters
            if (!this._cellLinePassesOncoprintFilters(cellLine)) return;

            indices.push(idx);
        });

        // Return all indices if no filters applied
        if (indices.length === 0 && !lineageFilter && !subLineageFilter &&
            (!hotspotGene || hotspotLevel === 'all') &&
            (!transGene || transLevel === 'all')) {
            return Array.from({ length: this.nCellLines }, (_, i) => i);
        }

        return indices;
    }

    runAnalysis() {
        // Reset network settings to defaults when running new analysis
        this.resetNetworkSettings();

        // Auto-load manual stats if the manual textarea has content (beyond default header)
        const manualTextarea = document.getElementById('manualStatsTextarea');
        const manualContent = manualTextarea?.value.trim();
        // Check if there's actual data beyond the header line
        if (manualContent && manualContent.split('\n').filter(l => l.trim()).length > 1) {
            this.loadManualStats();
        }

        const mode = document.querySelector('input[name="analysisMode"]:checked').value;

        // Handle mutation analysis mode separately
        if (mode === 'mutation') {
            this.runMutationAnalysis();
            return;
        }

        // Handle synonym/ortholog lookup mode
        if (mode === 'synonym') {
            this.runSynonymLookup();
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
        if (cellLineIndices.length < minN) {
            this.showStatus('error', `Too few cell lines for analysis (${cellLineIndices.length} available, ${minN} required). Adjust filters or reduce "Min Cell Lines" setting.`);
            return;
        }
        if (cellLineIndices.length < 10) {
            this.showStatus('error', `Too few cell lines match the filter (${cellLineIndices.length} found). Adjust filter settings.`);
            return;
        }

        const expandNetwork = mode === 'design' && document.getElementById('designExpandNetwork')?.checked;
        this.showStatus('info', expandNetwork ? 'Running correlation analysis (expanded network)...' : 'Running correlation analysis...');

        // Use setTimeout to allow UI to update
        setTimeout(() => {
            try {
                this.results = this.calculateCorrelations(geneList, mode, cutoff, minN, minSlope, cellLineIndices, expandNetwork);
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
        // Reset inspect-level filters and close modals
        const geTissueEl = document.getElementById('geTissueFilter');
        if (geTissueEl) geTissueEl.value = '';
        const geHotspotEl = document.getElementById('geHotspotFilter');
        if (geHotspotEl) geHotspotEl.value = '';
        document.getElementById('geneEffectModal').style.display = 'none';
        document.getElementById('geInlineCompareTable').style.display = 'none';
        this.currentGeneEffectGene = null;
        this.geneEffectViewMode = null;

        // Check mutation analysis sub-type
        const mutAnalysisType = document.querySelector('input[name="mutAnalysisType"]:checked')?.value || 'hotspot';
        const isTranslocation = mutAnalysisType === 'translocation';

        const hotspotGene = isTranslocation
            ? document.getElementById('translocationHotspotSelect').value
            : document.getElementById('mutationHotspotSelect').value;
        const minN = parseInt(document.getElementById('minCellLines').value);
        const pThreshold = this.getInputNum('pValueThreshold');
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;

        // Get additional hotspot filter (from parameter section)
        const additionalHotspot = document.getElementById('paramHotspotGene').value;
        const additionalHotspotLevel = document.getElementById('paramHotspotLevel').value;

        // Get additional translocation filter (from parameter section)
        const additionalTransGene = document.getElementById('paramTranslocationGene').value;
        const additionalTransLevel = document.getElementById('paramTranslocationLevel').value;

        if (!hotspotGene) {
            this.showStatus('error', isTranslocation ? 'Please select a translocation/fusion gene' : 'Please select a hotspot mutation');
            return;
        }
        if (isTranslocation && !this.translocations?.geneData?.[hotspotGene]) {
            this.showStatus('error', `"${hotspotGene}" is not a valid fusion gene. Please select from the list.`);
            return;
        }

        this.showStatus('info', isTranslocation ? 'Running fusion analysis...' : 'Running mutation analysis...');

        // Use setTimeout to allow UI to update
        setTimeout(() => {
            try {
                const analysisResult = isTranslocation
                    ? this.calculateTranslocationAnalysis(hotspotGene, minN, lineageFilter, subLineageFilter, additionalHotspot, additionalHotspotLevel, additionalTransGene, additionalTransLevel)
                    : this.calculateMutationAnalysis(hotspotGene, minN, lineageFilter, subLineageFilter, additionalHotspot, additionalHotspotLevel, additionalTransGene, additionalTransLevel);

                // Filter by p-value threshold
                const significantResults = analysisResult.results.filter(r => r.p_mut < pThreshold || r.p_2 < pThreshold || r.p_2v1 < pThreshold || (r.p_fused !== undefined && r.p_fused < pThreshold));

                // Sort by Δ GE (1+2 vs 0) ascending — most negative first
                significantResults.sort((a, b) => (a.diff_mut || 0) - (b.diff_mut || 0));

                // Close compare modal on new analysis
                if (document.getElementById('mutCompareModal')) {
                    document.getElementById('mutCompareModal').style.display = 'none';
                }

                this.mutationResults = {
                    hotspotGene,
                    pThreshold,
                    minN,
                    lineageFilter,
                    subLineageFilter,
                    additionalHotspot,
                    additionalHotspotLevel,
                    additionalTransGene,
                    additionalTransLevel,
                    isTranslocation,
                    excludedTissues: new Set(this.excludedTissues),
                    nWT: analysisResult.nWT,
                    nMut: analysisResult.nMut,
                    n2: analysisResult.n2,
                    hasFusionData: analysisResult.hasFusionData,
                    nFused: analysisResult.nFused,
                    nWTFusion: analysisResult.nWTFusion,
                    allResults: analysisResult.results,
                    significantResults
                };

                this.displayMutationResults(true);

                // Switch to mutation tab
                document.querySelectorAll('.nav-link').forEach(link => link.classList.remove('active'));
                document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
                document.querySelector('[data-tab="mutation"]').classList.add('active');
                document.getElementById('tab-mutation').classList.add('active');

                const analysisLabel = isTranslocation ? 'Fusion' : 'Mutation';
                const nSkipped = analysisResult.nSkippedMinN || 0;
                let statusMsg = `&#10003; ${analysisLabel} analysis complete: ${significantResults.length} genes with p < ${pThreshold}`;
                if (significantResults.length === 0 && nSkipped > 0) {
                    statusMsg = `&#9888; ${analysisLabel} analysis: no significant genes found. ${nSkipped} genes were skipped because WT group had fewer than ${minN} cell lines. Try lowering "Min Cell Lines" (currently ${minN}).`;
                    this.showStatus('warning', statusMsg);
                } else if (nSkipped > 0 && nSkipped > analysisResult.results.length) {
                    statusMsg += ` (${nSkipped} genes skipped — WT < ${minN} cells)`;
                    this.showStatus('success', statusMsg);
                } else {
                    this.showStatus('success', statusMsg);
                }
            } catch (error) {
                console.error('Mutation analysis error:', error);
                this.showStatus('error', 'Mutation analysis failed: ' + error.message);
            }
        }, 50);
    }

    async runSynonymLookup() {
        const textarea = document.getElementById('geneTextarea');
        const text = textarea.value.trim();

        if (!text) {
            this.showStatus('error', 'Please enter at least one gene');
            return;
        }

        const inputGenes = text.split(/[\n,\s]+/).map(g => g.trim()).filter(g => g.length > 0);

        if (inputGenes.length === 0) {
            this.showStatus('error', 'Please enter at least one gene');
            return;
        }

        this.showStatus('info', 'Looking up synonyms and orthologs...');

        // Load synonyms if not loaded
        if (!this.synonymLookup) {
            try {
                const synonymsRes = await fetch('web_data/synonyms.json');
                this.synonymLookup = await synonymsRes.json();
            } catch (e) {
                console.error('Failed to load synonyms:', e);
                this.synonymLookup = {};
            }
        }

        // Process each input gene
        const results = [];
        const genesForAPI = [];

        for (const gene of inputGenes) {
            const upperGene = gene.toUpperCase();
            const result = {
                input: gene,
                status: 'Not Found',
                official: '-',
                lowRisk: [],
                midRisk: [],
                orthologs: []
            };

            // Check if it's already a valid gene in the dataset
            if (this.geneIndex.has(upperGene)) {
                result.status = 'Valid';
                result.official = upperGene;
            }

            // Check synonym lookup (this maps aliases to official symbols)
            if (this.synonymLookup) {
                const match = this.synonymLookup[upperGene];
                if (match) {
                    const risk = match.r === 'l' ? 'lowRisk' : 'midRisk';
                    result[risk].push(match.d);
                    if (result.status === 'Not Found' && this.geneIndex.has(match.d.toUpperCase())) {
                        result.status = 'Synonym Found';
                        result.official = match.d;
                    }
                }

                // Also find all synonyms that point TO this gene (reverse lookup)
                Object.entries(this.synonymLookup).forEach(([alias, data]) => {
                    if (data.d.toUpperCase() === upperGene) {
                        const risk = data.r === 'l' ? 'lowRisk' : 'midRisk';
                        if (!result[risk].includes(alias)) {
                            result[risk].push(alias);
                        }
                    }
                });
            }

            // Check ortholog lookup (mouse to human)
            if (this.orthologs?.mouseToHuman) {
                const humanGene = this.orthologs.mouseToHuman[gene] || this.orthologs.mouseToHuman[upperGene];
                if (humanGene) {
                    result.orthologs.push(`${gene} → ${humanGene}`);
                    if (result.status === 'Not Found' && this.geneIndex.has(humanGene.toUpperCase())) {
                        result.status = 'Ortholog Found';
                        result.official = humanGene;
                    }
                }

                // Reverse lookup: find mouse genes that map to this human gene
                Object.entries(this.orthologs.mouseToHuman).forEach(([mouse, human]) => {
                    if (human.toUpperCase() === upperGene) {
                        const orthStr = `${mouse} → ${human}`;
                        if (!result.orthologs.includes(orthStr)) {
                            result.orthologs.push(orthStr);
                        }
                    }
                });
            }

            // If still not found, queue for API lookup
            if (result.status === 'Not Found') {
                genesForAPI.push({ gene, result });
            }

            results.push(result);
        }

        // Query MyGene.info API for remaining genes
        if (genesForAPI.length > 0) {
            try {
                const apiResults = await this.queryMyGeneAPI(genesForAPI.map(g => g.gene));
                apiResults.forEach((apiResult, idx) => {
                    if (apiResult.replacement) {
                        const result = genesForAPI[idx].result;
                        result.lowRisk.push(`${apiResult.replacement} (API)`);
                        if (this.geneIndex.has(apiResult.replacement.toUpperCase())) {
                            result.status = 'API Match';
                            result.official = apiResult.replacement;
                        }
                    }
                });
            } catch (error) {
                console.warn('MyGene.info API query failed:', error);
            }
        }

        // Store results and display
        this.synonymResults = results;
        this.displaySynonymResults();

        // Switch to synonyms tab
        document.querySelectorAll('.nav-link').forEach(link => link.classList.remove('active'));
        document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
        document.querySelector('[data-tab="synonyms"]').classList.add('active');
        document.getElementById('tab-synonyms').classList.add('active');

        this.showStatus('success', `&#10003; Synonym lookup complete: ${results.length} genes processed`);
    }

    displaySynonymResults() {
        const tbody = document.getElementById('synonymsTableBody');
        tbody.innerHTML = '';

        if (!this.synonymResults || this.synonymResults.length === 0) {
            tbody.innerHTML = '<tr><td colspan="6" style="padding: 40px; color: var(--gray-500);">No results</td></tr>';
            return;
        }

        // Check if any result has mid-risk data
        const hasMidRisk = this.synonymResults.some(r => r.midRisk.length > 0);

        // Show/hide the mid-risk column header and help text
        const table = document.getElementById('synonymsTable');
        const midRiskHeader = table.querySelector('thead th:nth-child(5)');
        if (midRiskHeader) midRiskHeader.style.display = hasMidRisk ? '' : 'none';
        const helpText = document.getElementById('synonymsHelpText');
        if (helpText) {
            helpText.innerHTML = hasMidRisk
                ? '<strong>Low Risk:</strong> Exact synonyms from HGNC. <strong>Mid Risk:</strong> Less certain synonyms. <strong>Orthologs:</strong> Mouse-to-human gene mappings.'
                : '<strong>Low Risk:</strong> Exact synonyms from HGNC. <strong>Orthologs:</strong> Mouse-to-human gene mappings.';
        }

        this.synonymResults.forEach(r => {
            const statusColor = r.status === 'Valid' ? '#16a34a' :
                               r.status === 'Not Found' ? '#dc2626' : '#f59e0b';
            const tr = document.createElement('tr');
            tr.innerHTML = `
                <td>${r.input}</td>
                <td style="color: ${statusColor}; font-weight: 500;">${r.status}</td>
                <td>${r.official}</td>
                <td style="border-left: 2px solid #16a34a;">${r.lowRisk.length > 0 ? r.lowRisk.join(', ') : '-'}</td>
                <td style="border-left: 2px solid #f59e0b; display: ${hasMidRisk ? '' : 'none'};">${r.midRisk.length > 0 ? r.midRisk.join(', ') : '-'}</td>
                <td style="border-left: 2px solid #6366f1;">${r.orthologs.length > 0 ? r.orthologs.join(', ') : '-'}</td>
            `;
            tbody.appendChild(tr);
        });
    }

    filterSynonymsTable(searchTerm) {
        const tbody = document.getElementById('synonymsTableBody');
        if (!tbody) return;

        const term = searchTerm.toLowerCase().trim();
        const rows = tbody.querySelectorAll('tr');

        rows.forEach(row => {
            const text = row.textContent.toLowerCase();
            row.style.display = term === '' || text.includes(term) ? '' : 'none';
        });
    }

    downloadSynonymsCSV() {
        if (!this.synonymResults || this.synonymResults.length === 0) return;

        const hasMidRisk = this.synonymResults.some(r => r.midRisk.length > 0);

        let csv = hasMidRisk
            ? 'Input_Gene,Status,Official_Symbol,Low_Risk_Synonyms,Mid_Risk_Synonyms,Orthologs\n'
            : 'Input_Gene,Status,Official_Symbol,Low_Risk_Synonyms,Orthologs\n';

        this.synonymResults.forEach(r => {
            if (hasMidRisk) {
                csv += `"${r.input}","${r.status}","${r.official}","${r.lowRisk.join('; ')}","${r.midRisk.join('; ')}","${r.orthologs.join('; ')}"\n`;
            } else {
                csv += `"${r.input}","${r.status}","${r.official}","${r.lowRisk.join('; ')}","${r.orthologs.join('; ')}"\n`;
            }
        });

        this.downloadFile(csv, 'synonym_ortholog_lookup.csv', 'text/csv');
    }

    calculateMutationAnalysis(hotspotGene, minN, lineageFilter, subLineageFilter, additionalHotspot, additionalHotspotLevel, additionalTransGene, additionalTransLevel) {
        const mutationData = this.mutations.geneData[hotspotGene];
        if (!mutationData) {
            throw new Error(`No mutation data for ${hotspotGene}`);
        }

        // Get additional hotspot mutation data if specified
        const additionalMutData = additionalHotspot ? this.mutations.geneData[additionalHotspot] : null;

        // Get additional translocation filter data if specified
        const additionalTransData = (additionalTransGene && this.translocations?.geneData?.[additionalTransGene])
            ? this.translocations.geneData[additionalTransGene].translocations : null;

        const cellLines = this.metadata.cellLines;
        const results = [];

        // Categorize cell lines by mutation status and lineage filter
        const wtCellIndices = [];
        const mut1CellIndices = [];
        const mut2CellIndices = [];

        cellLines.forEach((cellLine, idx) => {
            // Check excluded tissues
            if (this.excludedTissues && this.excludedTissues.size > 0) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                if (lineage && this.excludedTissues.has(lineage)) {
                    return;
                }
            }

            // Check lineage filter
            if (lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== lineageFilter) {
                return;
            }

            // Check sublineage filter
            if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== subLineageFilter) {
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

            // Check additional translocation filter
            if (additionalTransData && additionalTransLevel !== 'all') {
                const tLevel = additionalTransData[cellLine] || 0;
                if (additionalTransLevel === '0' && tLevel !== 0) return;
                if (additionalTransLevel === '1' && tLevel !== 1) return;
                if (additionalTransLevel === '2' && tLevel < 2) return;
                if (additionalTransLevel === '1+2' && tLevel < 1) return;
            }

            // Check oncoprint multi-gene filters
            if (!this._cellLinePassesOncoprintFilters(cellLine)) return;

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

        // Also categorize by fusion status for the same hotspot gene (if translocation data exists)
        const hotspotTransData = this.translocations?.geneData?.[hotspotGene]?.translocations;
        let fusedCellIndices = null;
        let wtFusionCellIndices = null;
        if (hotspotTransData) {
            fusedCellIndices = [];
            wtFusionCellIndices = [];
            // Use the same filtered cell set (union of wt + mut indices)
            const allFilteredIndices = [...wtCellIndices, ...mutAllCellIndices];
            allFilteredIndices.forEach(idx => {
                const cellLine = cellLines[idx];
                const tLevel = hotspotTransData[cellLine] || 0;
                if (tLevel >= 1) {
                    fusedCellIndices.push(idx);
                } else {
                    wtFusionCellIndices.push(idx);
                }
            });
        }
        const hasFusionData = fusedCellIndices && fusedCellIndices.length >= 3 && wtFusionCellIndices.length >= 3;

        // Analyze each gene
        let nSkippedMinN = 0;
        for (let geneIdx = 0; geneIdx < this.nGenes; geneIdx++) {
            const gene = this.geneNames[geneIdx];

            // Get gene effect values for each group
            const wtEffects = this.getGeneEffectsForCells(geneIdx, wtCellIndices);
            const mut1Effects = this.getGeneEffectsForCells(geneIdx, mut1CellIndices);
            const mutAllEffects = this.getGeneEffectsForCells(geneIdx, mutAllCellIndices);
            const mut2Effects = this.getGeneEffectsForCells(geneIdx, mut2CellIndices);

            // Skip if not enough valid values
            if (wtEffects.length < minN || mutAllEffects.length < 3) { nSkippedMinN++; continue; }

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

            // Calculate statistics for 2 vs 1 (dose-response)
            let n_1 = mut1Effects.length;
            let mean_1 = NaN;
            let diff_2v1 = NaN;
            let p_2v1 = 1;

            if (mut1Effects.length >= 3 && mut2Effects.length >= 3) {
                mean_1 = this.mean(mut1Effects);
                diff_2v1 = this.mean(mut2Effects) - mean_1;
                const tTest_2v1 = this.welchTTest(mut1Effects, mut2Effects);
                p_2v1 = tTest_2v1.p;
            }

            // Calculate fusion statistics (if hotspot gene has translocation data)
            let n_fused = 0, mean_fused = NaN, diff_fused = NaN, p_fused = 1, n_wt_fusion = 0;
            if (hasFusionData) {
                const fusedEffects = this.getGeneEffectsForCells(geneIdx, fusedCellIndices);
                const wtFusionEffects = this.getGeneEffectsForCells(geneIdx, wtFusionCellIndices);
                n_fused = fusedEffects.length;
                n_wt_fusion = wtFusionEffects.length;
                if (fusedEffects.length >= 3 && wtFusionEffects.length >= 3) {
                    mean_fused = this.mean(fusedEffects);
                    const wtFusionMean = this.mean(wtFusionEffects);
                    diff_fused = mean_fused - wtFusionMean;
                    const tTest_fused = this.welchTTest(wtFusionEffects, fusedEffects);
                    p_fused = tTest_fused.p;
                }
            }

            results.push({
                gene,
                n_wt: wtEffects.length,
                mean_wt: wtMean,
                n_1,
                mean_1,
                n_mut: mutAllEffects.length,
                mean_mut: mutMean,
                diff_mut,
                p_mut: tTest_mut.p,
                n_2,
                mean_2,
                diff_2,
                p_2,
                diff_2v1,
                p_2v1,
                n_fused,
                mean_fused,
                diff_fused,
                p_fused
            });
        }

        return {
            results,
            nWT: wtCellIndices.length,
            nMut: mutAllCellIndices.length,
            n2: mut2CellIndices.length,
            hasFusionData,
            nFused: fusedCellIndices?.length || 0,
            nWTFusion: wtFusionCellIndices?.length || 0,
            nSkippedMinN
        };
    }

    calculateTranslocationAnalysis(hotspotGene, minN, lineageFilter, subLineageFilter, additionalHotspot, additionalHotspotLevel, additionalTransGene, additionalTransLevel) {
        const transData = this.translocations.geneData[hotspotGene];
        if (!transData) {
            throw new Error(`No translocation data for ${hotspotGene}`);
        }

        const additionalMutData = additionalHotspot ? this.mutations?.geneData?.[additionalHotspot] : null;
        const additionalTransFilterData = (additionalTransGene && additionalTransGene !== hotspotGene && this.translocations?.geneData?.[additionalTransGene])
            ? this.translocations.geneData[additionalTransGene].translocations : null;

        const cellLines = this.metadata.cellLines;
        const results = [];

        const wtCellIndices = [];
        const mut1CellIndices = [];
        const mut2CellIndices = [];

        cellLines.forEach((cellLine, idx) => {
            if (this.excludedTissues && this.excludedTissues.size > 0) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                if (lineage && this.excludedTissues.has(lineage)) return;
            }
            if (lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== lineageFilter) return;
            if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== subLineageFilter) return;

            if (additionalMutData && additionalHotspotLevel !== 'all') {
                const addMutLevel = additionalMutData.mutations[cellLine] || 0;
                if (additionalHotspotLevel === '0' && addMutLevel !== 0) return;
                if (additionalHotspotLevel === '1' && addMutLevel !== 1) return;
                if (additionalHotspotLevel === '2' && addMutLevel < 2) return;
                if (additionalHotspotLevel === '1+2' && addMutLevel === 0) return;
            }

            if (additionalTransFilterData && additionalTransLevel !== 'all') {
                const tLevel = additionalTransFilterData[cellLine] || 0;
                if (additionalTransLevel === '0' && tLevel !== 0) return;
                if (additionalTransLevel === '1' && tLevel !== 1) return;
                if (additionalTransLevel === '2' && tLevel < 2) return;
                if (additionalTransLevel === '1+2' && tLevel < 1) return;
            }

            // Check oncoprint multi-gene filters
            if (!this._cellLinePassesOncoprintFilters(cellLine)) return;

            const transLevel = transData.translocations[cellLine] || 0;
            if (transLevel === 0) {
                wtCellIndices.push(idx);
            } else if (transLevel === 1) {
                mut1CellIndices.push(idx);
            } else {
                mut2CellIndices.push(idx);
            }
        });

        const mutAllCellIndices = [...mut1CellIndices, ...mut2CellIndices];

        if (wtCellIndices.length < 3 || mutAllCellIndices.length < 3) {
            throw new Error(`Not enough cell lines: WT=${wtCellIndices.length}, Fused=${mutAllCellIndices.length}`);
        }

        let nSkippedMinN = 0;
        for (let geneIdx = 0; geneIdx < this.nGenes; geneIdx++) {
            const gene = this.geneNames[geneIdx];

            const wtEffects = this.getGeneEffectsForCells(geneIdx, wtCellIndices);
            const mutAllEffects = this.getGeneEffectsForCells(geneIdx, mutAllCellIndices);
            const mut2Effects = this.getGeneEffectsForCells(geneIdx, mut2CellIndices);

            if (wtEffects.length < minN || mutAllEffects.length < 3) { nSkippedMinN++; continue; }

            const wtMean = this.mean(wtEffects);
            const mutMean = this.mean(mutAllEffects);
            const diff_mut = mutMean - wtMean;
            const tTest_mut = this.welchTTest(wtEffects, mutAllEffects);

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
            n2: mut2CellIndices.length,
            hasFusionData: false,
            nSkippedMinN
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

    displayMutationResults(resetSortIndicator = false) {
        if (!this.mutationResults) return;

        const mr = this.mutationResults;
        const results = mr.significantResults;
        const hasFusion = mr.hasFusionData && mr.isTranslocation;
        const tbody = document.getElementById('mutationTableBody');
        tbody.innerHTML = '';

        // Build dynamic header with hotspot gene name
        if (resetSortIndicator) {
            this._mutTableSortCol = 'diff_mut';
            this._mutTableSortDir = 'asc';
        }
        const sortCol = this._mutTableSortCol || 'diff_mut';
        const sortDir = this._mutTableSortDir || 'asc';
        const hg = mr.hotspotGene || 'Hotspot';
        const isT = mr.isTranslocation;
        const wtLabel = isT ? `No ${hg} Fusion` : `${hg} WT`;
        const mutLbl = isT ? `${hg} Fused` : `${hg} Mut`;
        const thead = document.querySelector('#mutationTable thead');
        const thStyle = 'cursor: pointer;';
        const sortClick = 'onclick="app.sortMutationTable(this, event)"';
        const tip = 'title="Click to sort. Ctrl/Cmd+click to copy column."';
        const arrow = (col) => col === sortCol ? (sortDir === 'asc' ? ' ▲' : ' ▼') : '';
        const sortAttr = (col) => col === sortCol ? ` data-sort-dir="${sortDir}"` : '';
        const cols = [
            { col: 'gene', label: 'Gene', style: '' },
            { col: 'n_wt', label: `N (${wtLabel})`, style: 'border-left: 2px solid #2563eb;' },
            { col: 'mean_wt', label: `Mean GE (${wtLabel})`, style: '' },
            { col: 'n_mut', label: `N (${mutLbl} 1+2)`, style: 'border-left: 2px solid #f97316;' },
            { col: 'mean_mut', label: `Mean GE (${mutLbl} 1+2)`, style: '' },
            { col: 'diff_mut', label: 'Δ GE (1+2v0)', style: 'border-left: 2px solid #d1d5db;' },
            { col: 'p_mut', label: 'p-value', style: '' },
            { col: 'n_2', label: `N (${mutLbl} 2${isT ? '+' : ''})`, style: 'border-left: 2px solid #dc2626;', cls: 'mut2-col' },
            { col: 'mean_2', label: `Mean GE (${mutLbl} 2${isT ? '+' : ''})`, style: '', cls: 'mut2-col' },
            { col: 'diff_2', label: 'Δ GE (2v0)', style: 'border-left: 2px solid #d1d5db;', cls: 'mut2-col' },
            { col: 'p_2', label: 'p-value (2v0)', style: '', cls: 'mut2-col' },
            { col: 'diff_2v1', label: 'Δ GE (2v1)', style: 'border-left: 2px solid #d1d5db;', cls: 'mut2-col' },
            { col: 'p_2v1', label: 'p (2v1)', style: '', cls: 'mut2-col' }
        ];
        let headerHTML = '<tr><th></th>';
        cols.forEach(c => {
            headerHTML += `<th ${sortClick} data-col="${c.col}" ${tip} class="${c.cls || ''}" style="${thStyle} ${c.style}"${sortAttr(c.col)}>${c.label}${arrow(c.col)}</th>`;
        });
        if (hasFusion) {
            const fusionCols = [
                { col: 'n_fused', label: 'N (Fused)', style: 'border-left: 2px solid #8b5cf6;' },
                { col: 'mean_fused', label: 'Mean GE (F)', style: '' },
                { col: 'diff_fused', label: 'Δ GE (F)', style: 'border-left: 2px solid #d1d5db;' },
                { col: 'p_fused', label: 'p (F)', style: '' }
            ];
            fusionCols.forEach(c => {
                headerHTML += `<th ${sortClick} data-col="${c.col}" class="fusion-col" ${tip} style="${thStyle} ${c.style}"${sortAttr(c.col)}>${c.label}${arrow(c.col)}</th>`;
            });
        }
        headerHTML += '</tr>';
        thead.innerHTML = headerHTML;

        // Show/hide Compare by Fusion button
        const hasTranslocations = this.translocations?.genes?.length > 0;
        document.getElementById('mutCompareByFusionBtn').style.display = hasTranslocations ? '' : 'none';

        results.forEach(r => {
            const row = document.createElement('tr');
            let html = `
                <td><a href="#" class="inspect-link" onclick="app.showGeneEffectDistribution('${r.gene}'); return false;">Inspect</a></td>
                <td class="gene-hover" data-gene="${r.gene}"><a href="#" style="color: var(--green-700); text-decoration: none; cursor: pointer;" onclick="app.openGeneEffectModal('${r.gene}', 'tissue'); return false;" onmouseover="this.style.textDecoration='underline'" onmouseout="this.style.textDecoration='none'">${r.gene}</a></td>
                <td style="border-left: 2px solid #2563eb;">${r.n_wt}</td>
                <td>${r.mean_wt.toFixed(2)}</td>
                <td style="border-left: 2px solid #f97316;">${r.n_mut}</td>
                <td>${r.mean_mut.toFixed(2)}</td>
                <td class="${r.diff_mut < 0 ? 'negative' : 'positive'}">${r.diff_mut.toFixed(2)}</td>
                <td>${this.formatPValue(r.p_mut)}</td>
                <td class="mut2-col" style="border-left: 2px solid #dc2626;">${r.n_2}</td>
                <td class="mut2-col">${isNaN(r.mean_2) ? '-' : r.mean_2.toFixed(2)}</td>
                <td class="mut2-col ${r.diff_2 < 0 ? 'negative' : 'positive'}">${isNaN(r.diff_2) ? '-' : r.diff_2.toFixed(2)}</td>
                <td class="mut2-col">${this.formatPValue(r.p_2)}</td>
                <td class="mut2-col ${r.diff_2v1 < 0 ? 'negative' : 'positive'}" style="border-left: 2px solid #d1d5db;">${isNaN(r.diff_2v1) ? '-' : r.diff_2v1.toFixed(2)}</td>
                <td class="mut2-col">${isNaN(r.diff_2v1) ? '-' : this.formatPValue(r.p_2v1)}</td>
            `;
            if (hasFusion) {
                html += `
                    <td class="fusion-col" style="border-left: 2px solid #8b5cf6;">${r.n_fused || 0}</td>
                    <td class="fusion-col">${isNaN(r.mean_fused) ? '-' : r.mean_fused.toFixed(2)}</td>
                    <td class="fusion-col ${r.diff_fused < 0 ? 'negative' : 'positive'}">${isNaN(r.diff_fused) ? '-' : r.diff_fused.toFixed(2)}</td>
                    <td class="fusion-col">${isNaN(r.diff_fused) ? '-' : this.formatPValue(r.p_fused)}</td>
                `;
            }
            row.innerHTML = html;
            tbody.appendChild(row);
        });

        // Build settings summary
        const typeLabel = mr.isTranslocation ? 'Fusion Gene' : 'Hotspot';
        const mutLabel = mr.isTranslocation ? 'Fused' : 'Mutated';
        let settingsText = `${typeLabel}: ${mr.hotspotGene} | `;
        settingsText += `WT: ${mr.nWT} cells | ${mutLabel}: ${mr.nMut} cells`;
        if (hasFusion) {
            settingsText += ` | Fused: ${mr.nFused} cells`;
        }
        settingsText += ` | `;
        settingsText += `Min cells: ${mr.minN} | p < ${mr.pThreshold}`;
        if (mr.lineageFilter) {
            let lineageText = mr.lineageFilter;
            if (mr.subLineageFilter) {
                lineageText += ` (${mr.subLineageFilter})`;
            }
            settingsText += ` | Lineage: ${lineageText}`;
        }
        // Collect all mutation filters into one label
        const mutFilterParts = [];
        const shownGenes = new Set();
        if (mr.hotspotGene) shownGenes.add(mr.hotspotGene);
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            const ll = { '0': 'WT', '1': 'Mut', '2': 'Mut', '1+2': 'Mut' };
            mutFilterParts.push(`${mr.additionalHotspot} ${ll[mr.additionalHotspotLevel] || mr.additionalHotspotLevel}`);
            shownGenes.add(mr.additionalHotspot);
        }
        if (mr.additionalTransGene && mr.additionalTransLevel !== 'all') {
            const ll = { '0': 'WT', '1': 'Fused', '2': 'Fused', '1+2': 'Fused' };
            mutFilterParts.push(`${mr.additionalTransGene} ${ll[mr.additionalTransLevel] || mr.additionalTransLevel}`);
            shownGenes.add(mr.additionalTransGene);
        }
        if (this._activeOncoprintFilters && this._activeOncoprintFilters.length > 0) {
            for (const f of this._activeOncoprintFilters) {
                if (!shownGenes.has(f.gene)) {
                    mutFilterParts.push(`${f.gene} ${f.state === 'mut' ? 'Mut' : 'WT'}`);
                    shownGenes.add(f.gene);
                }
            }
        }
        if (mutFilterParts.length > 0) {
            settingsText += ` | Mutation filter: ${mutFilterParts.join(', ')}`;
        }
        if (mr.excludedTissues && mr.excludedTissues.size > 0) {
            // Show included tissues (inverse of excluded) for clarity
            const allLineages = this.cellLineMetadata?.lineage
                ? [...new Set(Object.values(this.cellLineMetadata.lineage))].sort()
                : [];
            const included = allLineages.filter(t => !mr.excludedTissues.has(t));
            if (included.length <= 5) {
                settingsText += ` | Tissues: ${included.join(', ')}`;
            } else {
                settingsText += ` | ${mr.excludedTissues.size} tissues excluded`;
            }
        }

        document.getElementById('mutationResultsCount').innerHTML =
            `<strong>${results.length} genes</strong> with p &lt; ${mr.pThreshold}<br>
            <small style="color: #666;">${settingsText}</small>`;

        // Store for sorting
        this.mutationTableData = results;

        // Show/hide mut2 column toggle
        const mut2Toggle = document.getElementById('mut2ColToggle');
        if (mut2Toggle) {
            mut2Toggle.style.display = '';
            const showMut2 = document.getElementById('showMut2Cols');
            if (showMut2 && !showMut2.checked) {
                document.querySelectorAll('#mutationTable .mut2-col').forEach(el => el.style.display = 'none');
            }
        }

        // Attach gene tooltips
        this.attachGeneTooltips(tbody);
    }

    toggleMut2Columns() {
        const show = document.getElementById('showMut2Cols')?.checked;
        document.querySelectorAll('#mutationTable .mut2-col').forEach(el => {
            el.style.display = show ? '' : 'none';
        });
    }

    filterMutationTable(query) {
        const tbody = document.getElementById('mutationTableBody');
        const rows = tbody.querySelectorAll('tr');
        const lowerQuery = query.toLowerCase();

        rows.forEach(row => {
            // Gene name is in cells[1] (cells[0] is the Inspect link)
            const gene = row.cells[1]?.textContent.toLowerCase() || '';
            row.style.display = gene.includes(lowerQuery) ? '' : 'none';
        });
    }

    sortMutationTable(th, event) {
        if (event && (event.ctrlKey || event.metaKey)) {
            const colIndex = Array.from(th.parentNode.children).indexOf(th);
            this.copyColumnToClipboard(th.closest('table'), colIndex);
            return;
        }
        if (!this.mutationTableData) return;

        const col = th.dataset.col;
        const currentDir = th.dataset.sortDir === 'asc' ? 'desc' : 'asc';
        this._mutTableSortCol = col;
        this._mutTableSortDir = currentDir;

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

    async _exportMutationTable(format) {
        const table = document.getElementById('mutationTable');
        const settingsEl = document.getElementById('mutationResultsCount');
        if (!table) return;

        const filename = `mutation_analysis_${this.mutationResults?.hotspotGene || 'table'}`;

        const headers = [];
        table.querySelectorAll('thead th').forEach(th => {
            if (th.offsetParent !== null || th.offsetWidth > 0) headers.push(th.textContent.trim());
        });
        const rows = [];
        table.querySelectorAll('tbody tr').forEach(tr => {
            if (tr.style.display === 'none') return;
            const cells = [];
            tr.querySelectorAll('td').forEach(td => {
                if (td.offsetParent !== null || td.offsetWidth > 0) cells.push(td.textContent.trim());
            });
            if (cells.length > 0) rows.push(cells);
        });
        const settingsText = settingsEl ? settingsEl.textContent.trim() : '';

        if (format === 'png') {
            const tableContainer = table.closest('.table-container') || table.parentElement;
            try {
                const origMaxH = tableContainer.style.maxHeight;
                const origOverflow = tableContainer.style.overflow;
                tableContainer.style.maxHeight = 'none';
                tableContainer.style.overflow = 'visible';
                const canvas = await html2canvas(tableContainer, { scale: 2, backgroundColor: '#ffffff', scrollY: -window.scrollY });
                tableContainer.style.maxHeight = origMaxH;
                tableContainer.style.overflow = origOverflow;
                const a = document.createElement('a');
                a.href = canvas.toDataURL('image/png');
                a.download = `${filename}.png`;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
            } catch (e) {
                console.error('Table PNG export failed:', e);
                alert('PNG export failed. Try SVG instead.');
            }
        } else {
            const cellW = 70, cellH = 20, headerH = 28, pad = 10, titleH = settingsText ? 20 : 0;
            const nCols = headers.length;
            const nRows = rows.length;
            const firstColW = 90;
            const totalW = firstColW + (nCols - 1) * cellW + pad * 2;
            const totalH = titleH + headerH + nRows * cellH + pad * 2;

            let svg = `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${totalW}" height="${totalH}" viewBox="0 0 ${totalW} ${totalH}">\n`;
            svg += `<rect width="${totalW}" height="${totalH}" fill="white"/>\n`;
            if (settingsText) {
                svg += `<text x="${pad}" y="${pad + 12}" font-family="Arial" font-size="10" fill="#374151">${this.escapeXml(settingsText)}</text>\n`;
            }
            const startY = pad + titleH;
            svg += `<rect x="${pad}" y="${startY}" width="${totalW - pad * 2}" height="${headerH}" fill="#f0fdf4"/>\n`;
            headers.forEach((h, i) => {
                const x = pad + (i === 0 ? firstColW / 2 : firstColW + (i - 1) * cellW + cellW / 2);
                svg += `<text x="${x}" y="${startY + headerH / 2 + 4}" font-family="Arial" font-size="9" font-weight="bold" fill="#1a4a1a" text-anchor="middle">${this.escapeXml(h)}</text>\n`;
            });
            rows.forEach((row, ri) => {
                const y = startY + headerH + ri * cellH;
                if (ri % 2 === 1) svg += `<rect x="${pad}" y="${y}" width="${totalW - pad * 2}" height="${cellH}" fill="#f9fafb"/>\n`;
                row.forEach((cell, ci) => {
                    const x = pad + (ci === 0 ? firstColW / 2 : firstColW + (ci - 1) * cellW + cellW / 2);
                    const color = cell.match(/^-?\d/) && parseFloat(cell) < 0 ? '#dc2626' : parseFloat(cell) > 0 && ci > 0 ? '#16a34a' : '#374151';
                    svg += `<text x="${x}" y="${y + cellH / 2 + 4}" font-family="Arial" font-size="9" fill="${color}" text-anchor="middle">${this.escapeXml(cell)}</text>\n`;
                });
                svg += `<line x1="${pad}" y1="${y + cellH}" x2="${totalW - pad}" y2="${y + cellH}" stroke="#e5e7eb" stroke-width="0.5"/>\n`;
            });
            svg += '</svg>';
            const blob = new Blob([svg], { type: 'image/svg+xml' });
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = `${filename}.svg`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(a.href);
        }
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
        if (mr.excludedTissues && mr.excludedTissues.size > 0) {
            const allLineages = this.cellLineMetadata?.lineage
                ? [...new Set(Object.values(this.cellLineMetadata.lineage))].sort()
                : [];
            const included = allLineages.filter(t => !mr.excludedTissues.has(t));
            csv += `# Tissues: ${included.join(', ')}\n`;
        } else {
            csv += `# Lineage filter: ${mr.lineageFilter || 'All lineages'}\n`;
        }
        if (mr.subLineageFilter) {
            csv += `# Subtype filter: ${mr.subLineageFilter}\n`;
        }
        // Collect all mutation filters
        const allMutFilters = [];
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            const ll = { '0': 'WT', '1': 'Mut', '2': 'Mut', '1+2': 'Mut' };
            allMutFilters.push(`${mr.additionalHotspot} ${ll[mr.additionalHotspotLevel] || mr.additionalHotspotLevel}`);
        }
        if (this._activeOncoprintFilters) {
            const shown = new Set([mr.hotspotGene, mr.additionalHotspot].filter(Boolean));
            for (const f of this._activeOncoprintFilters) {
                if (!shown.has(f.gene)) {
                    allMutFilters.push(`${f.gene} ${f.state === 'mut' ? 'Mut' : 'WT'}`);
                    shown.add(f.gene);
                }
            }
        }
        if (allMutFilters.length > 0) {
            csv += `# Mutation filter: ${allMutFilters.join(', ')}\n`;
        }
        csv += `# Date: ${new Date().toISOString().slice(0, 10)}\n`;
        csv += '#\n';

        const hasFusion = mr.hasFusionData && mr.isTranslocation;
        let headers = ['Gene', 'N_WT', 'Mean_GE_WT', 'N_1+2', 'Mean_GE_1+2', 'Delta_GE', 'pValue_1+2_vs_0',
                        'N_2', 'Mean_GE_2', 'Delta_GE_2vs0', 'pValue_2_vs_0',
                        'Delta_GE_2vs1', 'pValue_2_vs_1'];
        if (hasFusion) {
            headers.push('N_Fused', 'Mean_GE_Fused', 'Delta_GE_Fused', 'pValue_Fused');
        }

        csv += headers.join(',') + '\n';
        results.forEach(r => {
            const row = [
                r.gene,
                r.n_wt,
                r.mean_wt.toFixed(2),
                r.n_mut,
                r.mean_mut.toFixed(2),
                r.diff_mut.toFixed(2),
                this.formatPValue(r.p_mut),
                r.n_2,
                isNaN(r.mean_2) ? '' : r.mean_2.toFixed(2),
                isNaN(r.diff_2) ? '' : r.diff_2.toFixed(2),
                this.formatPValue(r.p_2),
                isNaN(r.diff_2v1) ? '' : r.diff_2v1.toFixed(2),
                isNaN(r.diff_2v1) ? '' : this.formatPValue(r.p_2v1)
            ];
            if (hasFusion) {
                row.push(
                    r.n_fused || 0,
                    isNaN(r.mean_fused) ? '' : r.mean_fused.toFixed(2),
                    isNaN(r.diff_fused) ? '' : r.diff_fused.toFixed(2),
                    isNaN(r.diff_fused) ? '' : this.formatPValue(r.p_fused)
                );
            }
            csv += row.join(',') + '\n';
        });

        const filename = `mutation_analysis_${mr.hotspotGene}_${new Date().toISOString().slice(0, 10)}.csv`;
        this.downloadFile(csv, filename, 'text/csv');
    }

    showGeneEffectDistribution(gene, tissueOverride, inspectHotspotOverride) {
        if (!this.mutationResults) return;

        const mr = this.mutationResults;
        const hotspotGene = mr.hotspotGene;
        const isTranslocation = mr.isTranslocation;
        const mutationData = isTranslocation
            ? this.translocations.geneData[hotspotGene]
            : this.mutations.geneData[hotspotGene];
        const geneIdx = this.geneIndex.get(gene.toUpperCase());

        if (geneIdx === undefined) {
            alert(`Gene ${gene} not found`);
            return;
        }

        // Get tissue filter - use override if provided, otherwise read from dropdown
        const inspectTissueFilter = tissueOverride !== undefined ? tissueOverride : (document.getElementById('geTissueFilter')?.value || '');

        // Compute inspect-level additional hotspot filter (from dropdown in inspect modal)
        const inspectHotspot = inspectHotspotOverride !== undefined
            ? inspectHotspotOverride
            : (document.getElementById('geHotspotFilter')?.value || '');

        // Compute inspect-level fusion filter
        const inspectFusion = document.getElementById('geFusionFilter')?.value || '';
        const inspFusionData = inspectFusion ? this.translocations?.geneData?.[inspectFusion]?.translocations : null;

        // Compute inspect-level subtype filter
        const inspectSubtype = document.getElementById('geSubtypeFilter')?.value || '';

        // Collect data for each cell line
        const cellLines = this.metadata.cellLines;
        const data = { wt: [], mut1: [], mut2: [] };
        this.currentGeneEffectData = []; // Store for CSV export

        cellLines.forEach((cellLine, idx) => {
            // If a tissue filter is active (from compare table or dropdown), it overrides the lineage filter
            if (inspectTissueFilter) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine] || '';
                if (lineage !== inspectTissueFilter) return;
            } else {
                // Apply same filters as mutation analysis
                if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) {
                    return;
                }
                // Apply excluded tissues from analysis snapshot
                if (mr.excludedTissues && mr.excludedTissues.size > 0) {
                    const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                    if (lineage && mr.excludedTissues.has(lineage)) return;
                }
            }

            // Check sublineage filter (analysis-level)
            if (!inspectTissueFilter && mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) {
                return;
            }

            // Check inspect-level subtype filter
            if (inspectSubtype && this.cellLineMetadata?.primaryDisease) {
                if (this.cellLineMetadata.primaryDisease[cellLine] !== inspectSubtype) return;
            }

            // Check additional hotspot filter
            if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
                const addMutData = this.mutations?.geneData?.[mr.additionalHotspot];
                if (addMutData) {
                    const addMutLevel = addMutData.mutations[cellLine] || 0;
                    if (mr.additionalHotspotLevel === '0' && addMutLevel !== 0) return;
                    if (mr.additionalHotspotLevel === '1' && addMutLevel !== 1) return;
                    if (mr.additionalHotspotLevel === '2' && addMutLevel < 2) return;
                    if (mr.additionalHotspotLevel === '1+2' && addMutLevel === 0) return;
                }
            }

            // Check additional translocation filter
            if (mr.additionalTransGene && mr.additionalTransLevel !== 'all') {
                const addTransData = this.translocations?.geneData?.[mr.additionalTransGene]?.translocations;
                if (addTransData) {
                    const tLevel = addTransData[cellLine] || 0;
                    if (mr.additionalTransLevel === '0' && tLevel !== 0) return;
                    if (mr.additionalTransLevel === '1' && tLevel !== 1) return;
                    if (mr.additionalTransLevel === '2' && tLevel < 2) return;
                    if (mr.additionalTransLevel === '1+2' && tLevel < 1) return;
                }
            }

            // Check oncoprint multi-gene filters
            if (!this._cellLinePassesOncoprintFilters(cellLine)) return;

            // Check inspect-level additional hotspot filter
            if (inspectHotspot) {
                const inspHotData = this.mutations?.geneData?.[inspectHotspot];
                if (inspHotData) {
                    const inspMutLevel = inspHotData.mutations[cellLine] || 0;
                    if (inspMutLevel === 0) return;
                }
            }

            // Check inspect-level fusion filter
            if (inspectFusion && inspFusionData) {
                const fusionLevel = inspFusionData[cellLine] || 0;
                if (fusionLevel < 1) return;
            }

            const ge = this.geneEffects[geneIdx * this.nCellLines + idx];
            if (isNaN(ge)) return;

            const mutLevel = isTranslocation
                ? (mutationData.translocations[cellLine] || 0)
                : (mutationData.mutations[cellLine] || 0);
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
        this.currentInspectHotspot = inspectHotspot;

        // Create jitter for y-axis
        const jitter = (base, spread = 0.15) => base + (Math.random() - 0.5) * spread;

        // Labels and colors depend on translocation mode
        const mut1Label = isTranslocation ? '1 fusion partner' : '1 mutation';
        const mut2Label = isTranslocation ? '2+ fusion partners' : '2 mutations';
        const color1 = '#3b82f6';
        const color2 = '#dc2626';

        // Build hover text with fusion partner info for translocation mode
        const makeHoverText = (d) => {
            let text = `${d.cellName}<br>${d.lineage}<br>GE: ${d.ge.toFixed(2)}`;
            if (isTranslocation && mutationData.partners) {
                const partners = mutationData.partners[d.cellLine];
                if (partners && partners.length > 0) {
                    text += `<br>Partners: ${partners.join(', ')}`;
                }
            }
            return text;
        };

        // Build Plotly traces
        const traces = [
            {
                x: data.wt.map(d => d.ge),
                y: data.wt.map(() => jitter(0)),
                mode: 'markers',
                type: 'scatter',
                name: `WT (n=${data.wt.length})`,
                marker: { color: '#888888', size: 8, opacity: 0.7 },
                text: data.wt.map(d => makeHoverText(d)),
                hoverinfo: 'text'
            },
            {
                x: data.mut1.map(d => d.ge),
                y: data.mut1.map(() => jitter(1)),
                mode: 'markers',
                type: 'scatter',
                name: `${mut1Label} (n=${data.mut1.length})`,
                marker: { color: color1, size: 8, opacity: 0.7 },
                text: data.mut1.map(d => makeHoverText(d)),
                hoverinfo: 'text'
            },
            {
                x: data.mut2.map(d => d.ge),
                y: data.mut2.map(() => jitter(2)),
                mode: 'markers',
                type: 'scatter',
                name: `${mut2Label} (n=${data.mut2.length})`,
                marker: { color: color2, size: 8, opacity: 0.7 },
                text: data.mut2.map(d => makeHoverText(d)),
                hoverinfo: 'text'
            }
        ];

        // Calculate stats for each group
        const calcStats = (arr) => {
            if (arr.length === 0) return { n: 0, mean: NaN, median: NaN };
            const values = arr.map(d => d.ge).sort((a, b) => a - b);
            const n = values.length;
            const mean = values.reduce((a, b) => a + b, 0) / n;
            const median = n % 2 === 0 ? (values[n/2 - 1] + values[n/2]) / 2 : values[Math.floor(n/2)];
            return { n, mean, median };
        };

        const wtStats = calcStats(data.wt);
        const mut1Stats = calcStats(data.mut1);
        const mut2Stats = calcStats(data.mut2);
        const mutAllStats = calcStats([...data.mut1, ...data.mut2]);

        const meanWT = wtStats.mean;
        const meanMut1 = mut1Stats.mean;
        const meanMut2 = mut2Stats.mean;

        // Calculate p-values
        let pWTvsMut = NaN, pWTvs2 = NaN;
        if (wtStats.n >= 3 && mutAllStats.n >= 3) {
            const tTest = this.welchTTest(data.wt.map(d => d.ge), [...data.mut1, ...data.mut2].map(d => d.ge));
            pWTvsMut = tTest.p;
        }
        if (wtStats.n >= 3 && mut2Stats.n >= 3) {
            const tTest = this.welchTTest(data.wt.map(d => d.ge), data.mut2.map(d => d.ge));
            pWTvs2 = tTest.p;
        }

        // Add mean lines
        const allGE = [...data.wt, ...data.mut1, ...data.mut2].map(d => d.ge);
        const xMin = Math.min(0, Math.min(...allGE)) - 0.1;
        const xMax = Math.max(0, Math.max(...allGE)) + 0.1;

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
                line: { color: color1, width: 3 },
                showlegend: false,
                hoverinfo: 'skip'
            });
        }
        if (!isNaN(meanMut2)) {
            traces.push({
                x: [meanMut2, meanMut2],
                y: [1.7, 2.3],
                mode: 'lines',
                line: { color: color2, width: 3 },
                showlegend: false,
                hoverinfo: 'skip'
            });
        }

        // Build subtitle with filter info
        let filterInfo = [];
        if (mr.lineageFilter) {
            let lineageText = mr.lineageFilter;
            if (mr.subLineageFilter) {
                lineageText += ` (${mr.subLineageFilter})`;
            }
            filterInfo.push(`Lineage: ${lineageText}`);
        }
        if (inspectTissueFilter) {
            let tissueText = `Tissue: ${inspectTissueFilter}`;
            if (inspectSubtype) tissueText += ` (${inspectSubtype})`;
            filterInfo.push(tissueText);
        }
        if (inspectSubtype && !inspectTissueFilter) {
            filterInfo.push(`Subtype: ${inspectSubtype}`);
        }
        if (mr.excludedTissues && mr.excludedTissues.size > 0 && !inspectTissueFilter) {
            const allLineages = this.cellLineMetadata?.lineage
                ? [...new Set(Object.values(this.cellLineMetadata.lineage))].sort()
                : [];
            const included = allLineages.filter(t => !mr.excludedTissues.has(t));
            filterInfo.push(`Tissues: ${included.join(', ')}`);
        }
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            const ll = { '0': 'WT', '1': 'Mut', '2': 'Mut', '1+2': 'Mut' };
            filterInfo.push(`${mr.additionalHotspot} ${ll[mr.additionalHotspotLevel] || mr.additionalHotspotLevel}`);
        }
        if (mr.additionalTransGene && mr.additionalTransLevel !== 'all') {
            const ll = { '0': 'WT', '1': 'Fused', '2': 'Fused', '1+2': 'Fused' };
            filterInfo.push(`${mr.additionalTransGene} ${ll[mr.additionalTransLevel] || mr.additionalTransLevel}`);
        }
        // Include oncoprint multi-gene filters
        if (this._activeOncoprintFilters && this._activeOncoprintFilters.length > 0) {
            const shown = new Set([mr.hotspotGene, mr.additionalHotspot, mr.additionalTransGene].filter(Boolean));
            for (const f of this._activeOncoprintFilters) {
                if (!shown.has(f.gene)) {
                    filterInfo.push(`${f.gene} ${f.state === 'mut' ? 'Mut' : 'WT'}`);
                }
            }
        }
        if (inspectHotspot) {
            filterInfo.push(`Also ${inspectHotspot}-mutated`);
        }
        if (inspectFusion) {
            filterInfo.push(`Also ${inspectFusion}-fused`);
        }
        const lineageText = filterInfo.length > 0 ? filterInfo.join(' | ') : 'All lineages';

        // Build stats text for subtitle - condensed to one line
        const formatP = (p) => isNaN(p) ? '-' : (p < 0.001 ? p.toExponential(1) : p.toFixed(3));
        const fusedLabel = isTranslocation ? 'Fused' : 'Mut';
        let statsLine = `WT: n=${wtStats.n}, mean=${wtStats.mean.toFixed(2)}, med=${wtStats.median.toFixed(2)}  ·  ${fusedLabel}: n=${mutAllStats.n}, mean=${mutAllStats.mean.toFixed(2)}, med=${mutAllStats.median.toFixed(2)}`;
        statsLine += `  ·  p(WT vs ${fusedLabel}): ${formatP(pWTvsMut)}`;
        if (mut2Stats.n >= 3) {
            statsLine += `  ·  p(WT vs 2${isTranslocation ? '+' : ''}): ${formatP(pWTvs2)}`;
        }

        // Combine lineage info and stats in subtitle
        const subtitle = `${lineageText}<br>${statsLine}`;

        const statusLabel = isTranslocation ? 'Fusion Status' : 'Mutation Status';
        const yAxisTitle = isTranslocation ? `${hotspotGene} Fusions` : `${hotspotGene} Mutations`;
        const tick0Label = isTranslocation ? '0 No fusion' : '0 WT';
        const tick1Label = isTranslocation ? '1 partner' : '1';
        const tick2Label = isTranslocation ? '2+ partners' : '2';

        const titleText = `${gene} Gene Effect by ${hotspotGene} ${statusLabel}`;
        const subtitleText = subtitle;

        const layout = {
            annotations: [{
                text: `<b>${titleText}</b><br><span style="font-size:10px;color:#666">${subtitleText}</span>`,
                xref: 'paper',
                yref: 'paper',
                x: 0.5,
                y: 1.35,
                xanchor: 'center',
                yanchor: 'top',
                showarrow: false,
                font: { size: 13 }
            }],
            xaxis: {
                title: `${gene} Gene Effect`,
                range: [xMin, xMax]
            },
            yaxis: {
                title: yAxisTitle,
                tickmode: 'array',
                tickvals: [0, 1, 2],
                ticktext: [`${tick0Label} (n=${data.wt.length})`, `${tick1Label} (n=${data.mut1.length})`, `${tick2Label} (n=${data.mut2.length})`],
                range: [-0.5, 2.5]
            },
            showlegend: false,
            margin: { t: 160, r: 30, b: 55, l: 160 },
            height: Math.round(400 * (this.geChartHeightRatio || 1))
        };

        // Show modal
        document.getElementById('geneEffectModal').style.display = 'flex';
        document.getElementById('geneEffectTitle').textContent = `${gene} Gene Effect by ${hotspotGene} ${isTranslocation ? 'Fusion' : 'Mutation'}`;

        // Populate tissue filter dropdown with ALL lineages (inspect can override analysis filters)
        const tissueFilterEl = document.getElementById('geTissueFilter');
        if (tissueFilterEl) {
            const allLineages = [...new Set(cellLines.map(cl => this.cellLineMetadata?.lineage?.[cl]).filter(Boolean))].sort();
            let tHtml = '<option value="">All tissues</option>';
            for (const l of allLineages) {
                const sel = l === inspectTissueFilter ? ' selected' : '';
                tHtml += `<option value="${l}"${sel}>${l}</option>`;
            }
            tissueFilterEl.innerHTML = tHtml;

            // Pre-select lineage filter from analysis params if no inspect override
            if (!inspectTissueFilter && mr.lineageFilter) {
                tissueFilterEl.value = mr.lineageFilter;
                this.updateGeSubtypeFilter();
                if (mr.subLineageFilter) {
                    const geSubEl = document.getElementById('geSubtypeFilter');
                    if (geSubEl) geSubEl.value = mr.subLineageFilter;
                }
            } else {
                this.updateGeSubtypeFilter();
                // Restore subtype selection after repopulating dropdown
                if (inspectSubtype) {
                    const geSubEl = document.getElementById('geSubtypeFilter');
                    if (geSubEl) geSubEl.value = inspectSubtype;
                }
            }
        }

        // Populate inspect-level hotspot filter dropdown
        const hotspotFilterEl = document.getElementById('geHotspotFilter');
        if (hotspotFilterEl && this.mutations?.genes) {
            let hHtml = '<option value="">No hotspot filter</option>';
            for (const g of this.mutations.genes) {
                if (g === hotspotGene) continue;
                const sel = g === inspectHotspot ? ' selected' : '';
                hHtml += `<option value="${g}"${sel}>${g}</option>`;
            }
            hotspotFilterEl.innerHTML = hHtml;
        }

        // Populate inspect-level fusion filter dropdown (using pre-computed cache)
        const fusionFilterEl = document.getElementById('geFusionFilter');
        if (fusionFilterEl && this._fusionGeneCounts?.length > 0) {
            const currentFusion = fusionFilterEl.value;
            let html = '<option value="">No fusion filter</option>';
            for (const { gene, nFused } of this._fusionGeneCounts) {
                if (gene === hotspotGene) continue;
                const sel = gene === currentFusion ? ' selected' : '';
                html += `<option value="${gene}"${sel}>${gene} (${nFused} fused)</option>`;
            }
            fusionFilterEl.innerHTML = html;
        }

        // Populate and show hotspot gene selector (Y axis mutation/fusion)
        const hotspotGeneSelectEl = document.getElementById('geHotspotGeneSelect');
        if (hotspotGeneSelectEl) {
            const geneList = isTranslocation
                ? (this.translocations?.genes || [])
                : (this.mutations?.genes || []);
            let gHtml = '';
            for (const g of geneList) {
                const sel = g === hotspotGene ? ' selected' : '';
                gHtml += `<option value="${g}"${sel}>${g}</option>`;
            }
            hotspotGeneSelectEl.innerHTML = gHtml;
        }
        document.getElementById('geHotspotGeneGroup').style.display = '';

        // Show gene search bar so user can change the gene (#12)
        document.getElementById('geSearchBar').style.display = '';
        document.getElementById('geneEffectSearch').value = gene.toUpperCase();
        document.getElementById('geneEffectSummary').style.display = 'none';
        document.getElementById('geTableContainer').style.display = 'none';
        document.getElementById('geByTissueView').style.display = 'block';
        document.getElementById('geByHotspotView').style.display = 'none';

        // Make chart container full width
        document.getElementById('geChartContainer').style.flex = '1';

        // Mark this as mutation analysis view
        this.geneEffectViewMode = 'mutation';

        // Show mutation inspect controls, hide non-mutation view buttons
        document.getElementById('geHotspotFilter').style.display = '';
        document.getElementById('geFusionFilter').style.display =
            this.translocations?.genes?.length > 0 ? '' : 'none';
        document.getElementById('geCompareButtons').style.display = '';
        document.getElementById('geResetFiltersBtn').style.display = '';
        document.getElementById('geCompareByTranslocationBtn').style.display =
            this.translocations?.genes?.length > 0 ? '' : 'none';
        document.getElementById('geViewTissue').style.display = 'none';
        document.getElementById('geViewHotspot').style.display = 'none';
        // Hide the "View:" label too (previous sibling span)
        const viewLabel = document.getElementById('geViewTissue').previousElementSibling;
        if (viewLabel && viewLabel.textContent.trim() === 'View:') viewLabel.style.display = 'none';
        if (!this._keepInlineCompare) {
            document.getElementById('geInlineCompareTable').style.display = 'none';
        }
        this._keepInlineCompare = false;

        // Apply current width ratio to container
        const container = document.getElementById('geChartContainer');
        const ratio = this.geChartWidthRatio || 1;
        if (container) {
            container.style.flex = `0 0 ${Math.round(ratio * 55)}%`;
        }

        Plotly.newPlot('geneEffectPlot', traces, layout, {
            responsive: true,
            edits: { annotationPosition: true, legendPosition: true }
        });

        // Show chart width and Y range controls
        this.updateShowAllButton();
    }

    _exportMutationInspectChart(format) {
        if (this.geneEffectViewMode !== 'mutation') return;
        const plotEl = document.getElementById('geneEffectPlot');
        if (!plotEl || !plotEl.data) return;

        const filename = `gene_effect_${this.currentGeneEffectGene}_${this.mutationResults.hotspotGene}`;

        // Export at on-screen size, then post-process SVG to expand viewBox to fit all content
        const exportWidth = plotEl.offsetWidth;
        const exportHeight = plotEl.offsetHeight;

        Plotly.toImage(plotEl, {
            format: 'svg',
            width: exportWidth,
            height: exportHeight
        }).then(async svgDataUrl => {
            // Decode SVG
            let svgString;
            if (svgDataUrl.indexOf('base64,') > -1) {
                svgString = atob(svgDataUrl.split('base64,')[1]);
            } else {
                svgString = decodeURIComponent(svgDataUrl.split(',').slice(1).join(','));
            }

            const parser = new DOMParser();
            const svgDoc = parser.parseFromString(svgString, 'image/svg+xml');
            const svgEl = svgDoc.documentElement;

            // Expand viewBox to fit all content:
            // 1. Insert into DOM so getBBox works
            // 2. Remove all clipPaths that might restrict measurement
            // 3. Measure, then restore clipPaths for final output
            const measurer = document.createElement('div');
            measurer.style.cssText = 'position:absolute; left:-99999px; top:-99999px;';
            document.body.appendChild(measurer);
            const measureSvg = svgEl.cloneNode(true);
            measureSvg.style.overflow = 'visible';
            // Temporarily strip all clip-path attributes so getBBox sees full extent
            measureSvg.querySelectorAll('[clip-path]').forEach(el => {
                el.removeAttribute('clip-path');
            });
            measurer.appendChild(measureSvg);

            try {
                const bbox = measureSvg.getBBox();
                const pad = 10;
                const vbX = Math.min(0, Math.floor(bbox.x - pad));
                const vbY = Math.min(0, Math.floor(bbox.y - pad));
                const vbW = Math.max(exportWidth, Math.ceil(bbox.x + bbox.width + pad)) - vbX;
                const vbH = Math.max(exportHeight, Math.ceil(bbox.y + bbox.height + pad)) - vbY;

                svgEl.setAttribute('viewBox', `${vbX} ${vbY} ${vbW} ${vbH}`);
                svgEl.setAttribute('width', String(vbW));
                svgEl.setAttribute('height', String(vbH));
            } catch (e) {
                console.warn('getBBox failed, keeping original dimensions', e);
            }
            document.body.removeChild(measurer);

            svgString = new XMLSerializer().serializeToString(svgEl);

            const meta = this._buildExportMetadata('mutation_inspect', {
                gene: this.currentGeneEffectGene,
                hotspotGene: this.mutationResults?.hotspotGene,
                isTranslocation: this.mutationResults?.isTranslocation || false,
                lineageFilter: this.mutationResults?.lineageFilter || '',
                textSettings: this._capturePlotTextSettings('geneEffectPlot'),
                geChartWidthRatio: this.geChartWidthRatio || 1.0,
                oncoprintFilters: this._activeOncoprintFilters || null
            });
            const metaJson = JSON.stringify(meta);

            const a = document.createElement('a');
            if (format === 'svg') {
                svgString = svgString.replace('</svg>', `<metadata><correlate-meta>${metaJson}</correlate-meta></metadata></svg>`);
                svgString = await this._finalizeSvgForExport(svgString);
                const blob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
                a.href = URL.createObjectURL(blob);
                a.download = `${filename}.svg`;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
            } else {
                // Render the fixed SVG to canvas at 4x for publication quality PNG
                const svgBlob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
                const svgUrl = URL.createObjectURL(svgBlob);
                const img = new Image();
                img.onload = async () => {
                    const scale = 4;
                    const canvas = document.createElement('canvas');
                    canvas.width = img.naturalWidth * scale;
                    canvas.height = img.naturalHeight * scale;
                    const ctx = canvas.getContext('2d');
                    ctx.scale(scale, scale);
                    ctx.fillStyle = 'white';
                    ctx.fillRect(0, 0, img.naturalWidth, img.naturalHeight);
                    ctx.drawImage(img, 0, 0);
                    URL.revokeObjectURL(svgUrl);
                    const pngDataUrl = canvas.toDataURL('image/png');
                    const pngResp = await fetch(pngDataUrl);
                    const pngBuf = await pngResp.arrayBuffer();
                    const pngWithMeta = this._addPngTextChunk(pngBuf, 'correlate-meta', metaJson);
                    const blob = new Blob([pngWithMeta], { type: 'image/png' });
                    a.href = URL.createObjectURL(blob);
                    a.download = `${filename}.png`;
                    document.body.appendChild(a);
                    a.click();
                    document.body.removeChild(a);
                    URL.revokeObjectURL(a.href);
                };
                img.src = svgUrl;
            }
        });
    }

    downloadGeneEffectPNG() {
        this._exportMutationInspectChart('png');
    }

    downloadGeneEffectSVG() {
        this._exportMutationInspectChart('svg');
    }

    downloadGeneEffectCSV() {
        if (!this.currentGeneEffectData || !this.currentGeneEffectGene) return;

        const mr = this.mutationResults;
        let csv = `# Gene Effect Distribution Data\n`;
        csv += `# Gene: ${this.currentGeneEffectGene}\n`;
        csv += `# Hotspot Mutation: ${mr.hotspotGene}\n`;
        if (mr.excludedTissues && mr.excludedTissues.size > 0) {
            const allLineages = this.cellLineMetadata?.lineage
                ? [...new Set(Object.values(this.cellLineMetadata.lineage))].sort()
                : [];
            const included = allLineages.filter(t => !mr.excludedTissues.has(t));
            csv += `# Tissues: ${included.join(', ')}\n`;
        } else {
            csv += `# Lineage filter: ${mr.lineageFilter || 'All lineages'}\n`;
        }
        if (mr.subLineageFilter) {
            csv += `# Subtype filter: ${mr.subLineageFilter}\n`;
        }
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            csv += `# Additional filter: ${mr.additionalHotspot} = ${mr.additionalHotspotLevel}\n`;
        }
        csv += `# Date: ${new Date().toISOString().slice(0, 10)}\n`;
        csv += '#\n';
        csv += 'CellLine,CellLineName,Lineage,GeneEffect,MutationLevel\n';

        this.currentGeneEffectData.forEach(d => {
            csv += `${d.cellLine},${d.cellName},${d.lineage},${d.ge.toFixed(2)},${d.mutLevel}\n`;
        });

        const filename = `gene_effect_${this.currentGeneEffectGene}_${mr.hotspotGene}_data.csv`;
        this.downloadFile(csv, filename, 'text/csv');
    }

    calculateCorrelations(geneList, mode, cutoff, minN, minSlope, cellLineIndices, expandNetwork = false) {
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

        // Calculate correlations (first pass: input genes vs target genes)
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

        // Second pass for expanded network: find correlations between discovered genes
        if (mode === 'design' && expandNetwork && correlations.length > 0) {
            // Collect all discovered genes (not in original input)
            const discoveredGenes = new Set();
            correlations.forEach(c => {
                if (!geneList.includes(c.gene2)) discoveredGenes.add(c.gene2);
            });

            const discoveredArray = Array.from(discoveredGenes);
            if (discoveredArray.length > 1) {
                // Cache gene data for discovered genes
                const discoveredData = new Map();
                discoveredArray.forEach(gene => {
                    const idx = this.geneIndex.get(gene);
                    const fullData = this.getGeneData(idx);
                    discoveredData.set(gene, cellLineIndices.map(i => fullData[i]));
                });

                // Find correlations between discovered genes (pairwise)
                for (let i = 0; i < discoveredArray.length; i++) {
                    const gene1 = discoveredArray[i];
                    const data1 = discoveredData.get(gene1);

                    for (let j = i + 1; j < discoveredArray.length; j++) {
                        const gene2 = discoveredArray[j];
                        const data2 = discoveredData.get(gene2);

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
            }
        }

        if (correlations.length === 0) {
            return { success: false, error: `No correlations found (cutoff: ${cutoff}, min slope: ${minSlope}, min cells: ${minN}). Try lowering thresholds or adjusting filters.` };
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

        // Add uncorrelated input genes so they appear in the clusters table
        const correlatedGenes = new Set(clusterData.map(c => c.gene));
        geneList.forEach(gene => {
            if (!correlatedGenes.has(gene) && this.geneIndex.has(gene)) {
                const idx = this.geneIndex.get(gene);
                const fullData = this.getGeneData(idx);
                const allValidData = Array.from(fullData).filter(v => !isNaN(v));
                const allMean = allValidData.length > 0 ? allValidData.reduce((a, b) => a + b, 0) / allValidData.length : NaN;
                const allVariance = allValidData.length > 0 ? allValidData.reduce((a, b) => a + (b - allMean) ** 2, 0) / allValidData.length : NaN;
                const allSd = Math.sqrt(allVariance);
                const filteredData = cellLineIndices.map(i => fullData[i]).filter(v => !isNaN(v));
                const filtMean = filteredData.length > 0 ? filteredData.reduce((a, b) => a + b, 0) / filteredData.length : NaN;
                const filtVariance = filteredData.length > 0 ? filteredData.reduce((a, b) => a + (b - filtMean) ** 2, 0) / filteredData.length : NaN;
                const filtSd = Math.sqrt(filtVariance);

                clusterData.push({
                    gene: gene,
                    cluster: '-',
                    meanEffect: Math.round(allMean * 100) / 100,
                    sdEffect: Math.round(allSd * 100) / 100,
                    meanEffectFiltered: Math.round(filtMean * 100) / 100,
                    sdEffectFiltered: Math.round(filtSd * 100) / 100,
                    nAll: allValidData.length,
                    nFiltered: filteredData.length,
                    inGeneList: true,
                    hasCorrelation: false
                });
            }
        });
        // Mark correlated genes
        clusterData.forEach(c => {
            if (c.hasCorrelation === undefined) c.hasCorrelation = true;
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
        // Reset network settings to defaults
        this.resetNetworkSettings();

        // Switch to network tab FIRST so vis.js can calculate layout in visible container
        document.querySelectorAll('.nav-link').forEach(link => link.classList.remove('active'));
        document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
        document.querySelector('[data-tab="network"]').classList.add('active');
        document.getElementById('tab-network').classList.add('active');

        // Display all results
        this.displayNetwork();
        this.displayCorrelationsTable();
        this.displayClustersTable();
        this.displaySummary();
    }

    resetNetworkSettings() {
        // Reset sliders to default values
        document.getElementById('netFontSize').value = 16;
        document.getElementById('fontSizeBubble').textContent = '16';
        document.getElementById('netNodeSize').value = 25;
        document.getElementById('nodeSizeBubble').textContent = '25';
        document.getElementById('netEdgeWidth').value = 3;
        document.getElementById('edgeWidthBubble').textContent = '3';

        // Reset checkboxes
        document.getElementById('showGeneEffect').checked = false;
        document.getElementById('showGeneEffectSD').checked = false;
        document.getElementById('colorByGeneEffect').checked = false;
        document.getElementById('colorAbsoluteGE').checked = false;
        document.getElementById('colorByLFC').checked = false;
        document.getElementById('colorByFDR').checked = false;

        // Hide dependent options
        document.getElementById('showGESDGroup').style.display = 'none';
        document.getElementById('colorAbsoluteGroup').style.display = 'none';
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
                        (isInput ? '#5a9f4a' : '#a8d89a') : '#5a9f4a',
                    border: '#000000'
                },
                borderWidth: document.getElementById('networkNodeBorder')?.checked === false ? 0 : 2,
                title: titleLines.join('\n'),
                isSynonym: isSynonym,
                originalName: originalName
            });
        });

        // Add uncorrelated input genes as isolated nodes if checkbox is checked
        const showUncorrelated = document.getElementById('showUncorrelatedGenes')?.checked;
        if (showUncorrelated && this.results.geneList) {
            this.results.geneList.forEach(gene => {
                if (!geneSet.has(gene) && this.geneIndex.has(gene)) {
                    const idx = this.geneIndex.get(gene);
                    const fullData = this.getGeneData(idx);
                    const validData = Array.from(fullData).filter(v => !isNaN(v));
                    const meanEffect = validData.length > 0 ? (validData.reduce((a, b) => a + b, 0) / validData.length) : NaN;
                    const sd = validData.length > 0 ? Math.sqrt(validData.reduce((a, b) => a + (b - meanEffect) ** 2, 0) / validData.length) : NaN;

                    const originalName = synonymLookup.get(gene.toUpperCase());
                    const isSynonym = !!originalName;
                    const label = isSynonym ? `${gene}*` : gene;

                    let titleLines = [gene];
                    if (isSynonym) titleLines.push(`(synonym of ${originalName})`);
                    titleLines.push(`GE mean: ${isNaN(meanEffect) ? 'N/A' : meanEffect.toFixed(2)}`);
                    titleLines.push(`GE SD: ${isNaN(sd) ? 'N/A' : sd.toFixed(2)}`);
                    titleLines.push('(no correlations found)');

                    nodes.push({
                        id: gene,
                        label: label,
                        size: nodeSize,
                        font: { size: fontSize, color: '#999' },
                        color: { background: '#d1d5db', border: '#9ca3af' },
                        borderWidth: document.getElementById('networkNodeBorder')?.checked === false ? 0 : 2,
                        borderWidthSelected: 3,
                        title: titleLines.join('\n'),
                        isSynonym: isSynonym,
                        originalName: originalName || null
                    });
                    geneSet.add(gene);
                }
            });
        }

        // Track if any synonyms are in the network for legend display
        this.hasSynonymsInNetwork = synonymLookup.size > 0 &&
            Array.from(geneSet).some(g => synonymLookup.has(g.toUpperCase()));

        const data = { nodes: new vis.DataSet(nodes), edges: new vis.DataSet(edges) };

        // Adjust stabilization iterations based on network size
        const nodeCount = nodes.length;
        const stabilizationIterations = nodeCount > 50 ? 300 : 150;

        const options = {
            autoResize: false,  // Prevent auto-fit when container resizes
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

        // Add filter banner after vis.Network has set up its canvas
        const filterText = this._getNetworkFilterText();
        if (filterText) {
            const banner = document.createElement('div');
            banner.className = 'network-filter-banner';
            banner.textContent = filterText;
            container.appendChild(banner);

            // Make banner draggable
            let bDragX = 0, bDragY = 0, bDragging = false;
            banner.addEventListener('mousedown', (e) => {
                bDragging = true;
                const rect = banner.getBoundingClientRect();
                bDragX = e.clientX - rect.left;
                bDragY = e.clientY - rect.top;
                banner.style.transform = 'none';
                banner.style.left = rect.left - container.getBoundingClientRect().left + 'px';
                banner.style.top = rect.top - container.getBoundingClientRect().top + 'px';
                e.preventDefault();
            });
            document.addEventListener('mousemove', (e) => {
                if (!bDragging) return;
                const cr = container.getBoundingClientRect();
                const newLeft = e.clientX - cr.left - bDragX;
                const newTop = e.clientY - cr.top - bDragY;
                banner.style.left = newLeft + 'px';
                banner.style.top = newTop + 'px';
                this._netBannerPos = { x: newLeft, y: newTop };
            });
            document.addEventListener('mouseup', () => { bDragging = false; });
        }

        this.hiddenNodes = [];
        this.removeMode = false;
        this.selectMode = false;
        this.selectedNodes.clear();

        // Reset physics state for new network
        this.physicsEnabled = true;
        this.currentLayout = 0;
        const layoutBtn = document.getElementById('changeLayout');
        if (layoutBtn) layoutBtn.textContent = 'Layout';
        const physicsBtn = document.getElementById('togglePhysics');
        if (physicsBtn) {
            physicsBtn.textContent = 'Lock';
            physicsBtn.classList.remove('btn-active');
        }
        const removeModeBtn = document.getElementById('toggleRemoveMode');
        if (removeModeBtn) {
            removeModeBtn.classList.remove('btn-active');
            removeModeBtn.style.backgroundColor = '';
            removeModeBtn.style.borderColor = '';
            removeModeBtn.style.color = '';
        }
        const selectModeBtn = document.getElementById('toggleSelectMode');
        if (selectModeBtn) {
            selectModeBtn.classList.remove('btn-active');
            selectModeBtn.style.backgroundColor = '';
            selectModeBtn.style.borderColor = '';
            selectModeBtn.style.color = '';
        }
        document.getElementById('selectedNodesList').style.display = 'none';

        // After stabilization: resolve edge crossings, then lock large networks
        this.network.once('stabilizationIterationsDone', () => {
            this.resolveEdgeCrossings();
            if (nodeCount > 30) {
                this.network.setOptions({ physics: { enabled: false } });
                this.physicsEnabled = false;
                if (physicsBtn) {
                    physicsBtn.textContent = 'Unlock Nodes';
                    physicsBtn.classList.add('btn-active');
                }
            }
        });

        // Track state for click vs double-click vs drag
        let clickTimeout = null;
        let isDragging = false;
        let dragStartPos = null;

        // Track drag start
        this.network.on('dragStart', (params) => {
            isDragging = true;
            dragStartPos = params.pointer.canvas;
        });

        // Track drag end
        this.network.on('dragEnd', (params) => {
            // Check if actually moved (more than 5 pixels)
            if (dragStartPos) {
                const dx = params.pointer.canvas.x - dragStartPos.x;
                const dy = params.pointer.canvas.y - dragStartPos.y;
                const distance = Math.sqrt(dx * dx + dy * dy);
                // Only count as drag if moved more than 5 pixels
                if (distance > 5) {
                    isDragging = true;
                } else {
                    isDragging = false;
                }
            }
            // Reset drag state after a short delay to allow click to check it
            setTimeout(() => {
                isDragging = false;
                dragStartPos = null;
            }, 100);
        });

        // Hover over node to show gene info tooltip
        this.network.on('hoverNode', (params) => {
            const nodeId = params.node;
            const domEvent = params.event && params.event.center ?
                { clientX: params.event.center.x, clientY: params.event.center.y } :
                { clientX: params.pointer?.DOM?.x || 0, clientY: params.pointer?.DOM?.y || 0 };
            // Add offset for the network container position
            const container = document.getElementById('networkPlot');
            if (container) {
                const rect = container.getBoundingClientRect();
                domEvent.clientX += rect.left;
                domEvent.clientY += rect.top;
            }
            clearTimeout(this._networkTooltipTimer);
            this._networkTooltipTimer = setTimeout(() => {
                this.showGeneTooltip(domEvent, nodeId);
            }, 400);
        });
        this.network.on('blurNode', () => {
            clearTimeout(this._networkTooltipTimer);
            this.hideGeneTooltip();
        });

        // Double-click to open Gene Effect (node) or Inspect (edge)
        this.network.on('doubleClick', (params) => {
            clearTimeout(this._networkTooltipTimer);
            this.hideGeneTooltip();
            if (params.nodes.length > 0) {
                // Node double-clicked - open Gene Effect analysis
                const nodeId = params.nodes[0];
                this.openGeneEffectFromNetwork(nodeId);
            } else if (params.edges.length > 0) {
                // Edge double-clicked - open correlation inspect
                const edgeId = params.edges[0];
                const edge = this.networkData.edges.get(edgeId);
                if (edge) {
                    this.openInspectByGenes(edge.from, edge.to);
                }
            }
        });

        // Single click - used in Remove Mode and Select Mode
        this.network.on('click', (params) => {
            // Skip if we just finished dragging
            if (isDragging) {
                return;
            }

            if (params.nodes.length > 0) {
                const nodeId = params.nodes[0];

                // Select Mode: toggle gene selection and update input gene box
                if (this.selectMode) {
                    if (this.selectedNodes.has(nodeId)) {
                        this.selectedNodes.delete(nodeId);
                        // Restore original node color
                        const isInput = this.results.geneList.includes(nodeId);
                        this.networkData.nodes.update({
                            id: nodeId,
                            color: {
                                background: this.results.mode === 'design' ?
                                    (isInput ? '#5a9f4a' : '#a8d89a') : '#5a9f4a',
                                border: '#000000'
                            }
                        });
                    } else {
                        this.selectedNodes.add(nodeId);
                        // Highlight selected node
                        this.networkData.nodes.update({
                            id: nodeId,
                            color: { background: '#2563eb', border: '#1d4ed8' }
                        });
                    }
                    this.updateSelectedNodesList();
                    // Update gene textarea with only selected genes
                    document.getElementById('geneTextarea').value = Array.from(this.selectedNodes).join('\n');
                    this.updateGeneCount();
                    return;
                }

                // Remove Mode: remove node from network AND from gene input list
                if (this.removeMode) {
                    const node = this.networkData.nodes.get(nodeId);
                    if (node) {
                        this.hiddenNodes.push(node);
                        this.networkData.nodes.remove(nodeId);
                        this.updateRemovedNodesList();
                        // Also remove gene from input list
                        const textarea = document.getElementById('geneTextarea');
                        const genes = textarea.value.split(/[\n\r]+/).map(g => g.trim()).filter(g => g);
                        const remaining = genes.filter(g => g.toUpperCase() !== nodeId.toUpperCase());
                        textarea.value = remaining.join('\n');
                        this.updateGeneCount();
                    }
                    return;
                }
            }
        });

        // Show legend
        document.getElementById('networkLegend').style.display = 'flex';
        const legendNodeType = document.getElementById('legendNodeType');
        if (this.results.mode === 'design') {
            legendNodeType.innerHTML = `
                <strong>Node Type:</strong>
                <span class="legend-item"><span class="legend-dot" style="background: #5a9f4a;"></span> Input</span>
                <span class="legend-item"><span class="legend-dot" style="background: #a8d89a;"></span> Correlated</span>
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

    resolveEdgeCrossings() {
        if (!this.network || !this.networkData) return;

        const positions = this.network.getPositions();
        const edges = [];
        this.networkData.edges.forEach(e => edges.push(e));

        // Line segment intersection test (excludes shared endpoints)
        const cross = (a, b, c, d) => {
            const det = (c.x - a.x) * (d.y - a.y) - (d.x - a.x) * (c.y - a.y);
            const det2 = (c.x - b.x) * (d.y - b.y) - (d.x - b.x) * (c.y - b.y);
            if (det * det2 >= 0) return false;
            const det3 = (a.x - c.x) * (b.y - c.y) - (b.x - c.x) * (a.y - c.y);
            const det4 = (a.x - d.x) * (b.y - d.y) - (b.x - d.x) * (a.y - d.y);
            return det3 * det4 < 0;
        };

        const countCrossings = () => {
            let n = 0;
            for (let i = 0; i < edges.length; i++) {
                for (let j = i + 1; j < edges.length; j++) {
                    const e1 = edges[i], e2 = edges[j];
                    if (e1.from === e2.from || e1.from === e2.to || e1.to === e2.from || e1.to === e2.to) continue;
                    const p1 = positions[e1.from], p2 = positions[e1.to], p3 = positions[e2.from], p4 = positions[e2.to];
                    if (p1 && p2 && p3 && p4 && cross(p1, p2, p3, p4)) n++;
                }
            }
            return n;
        };

        let crossings = countCrossings();
        if (crossings === 0) return;

        // Try swapping pairs of nodes to reduce crossings
        const nodeIds = Object.keys(positions);
        let improved = true;
        while (improved && crossings > 0) {
            improved = false;
            for (let i = 0; i < nodeIds.length && !improved; i++) {
                for (let j = i + 1; j < nodeIds.length && !improved; j++) {
                    const a = nodeIds[i], b = nodeIds[j];
                    // Swap positions
                    const tmp = { x: positions[a].x, y: positions[a].y };
                    positions[a] = { x: positions[b].x, y: positions[b].y };
                    positions[b] = tmp;

                    const newCrossings = countCrossings();
                    if (newCrossings < crossings) {
                        crossings = newCrossings;
                        improved = true;
                    } else {
                        // Undo swap
                        const tmp2 = { x: positions[a].x, y: positions[a].y };
                        positions[a] = { x: positions[b].x, y: positions[b].y };
                        positions[b] = tmp2;
                    }
                }
            }
        }

        // Apply resolved positions
        for (const nodeId of nodeIds) {
            this.network.moveNode(nodeId, positions[nodeId].x, positions[nodeId].y);
        }
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

    toggleNetworkBorder(showBorder) {
        if (!this.network || !this.networkData) return;
        const borderWidth = showBorder ? 2 : 0;
        const nodeUpdates = [];
        this.networkData.nodes.forEach(node => {
            nodeUpdates.push({ id: node.id, borderWidth: borderWidth });
        });
        this.networkData.nodes.update(nodeUpdates);
        this.network.setOptions({ nodes: { borderWidth: borderWidth } });
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
                    <td class="gene-hover" data-gene="${c.gene1}" style="cursor: help;">${c.gene1}</td>
                    <td class="gene-hover" data-gene="${c.gene2}" style="cursor: help;">${c.gene2}</td>
                    <td class="${corrClass}">${c.correlation.toFixed(3)}</td>
                    <td>${c.slope.toFixed(3)}</td>
                    <td>${c.n}</td>
                    <td>${c.cluster}</td>
                    <td style="white-space: nowrap;">
                        <button class="btn btn-sm inspect-btn" style="padding: 2px 6px; font-size: 10px; background: #5a9f4a; color: white;" data-gene1="${c.gene1}" data-gene2="${c.gene2}">Correlate</button>
                        <button class="btn btn-sm tissue-btn" style="padding: 2px 6px; font-size: 10px; margin-left: 4px; background: #6b7280; color: white;" data-gene1="${c.gene1}" data-gene2="${c.gene2}">By Tissue</button>
                        <button class="btn btn-sm hotspot-btn" style="padding: 2px 6px; font-size: 10px; margin-left: 4px; background: #6b7280; color: white;" data-gene1="${c.gene1}" data-gene2="${c.gene2}">By Hotspot</button>
                    </td>
                `;
                // Add click handlers
                tr.querySelector('.inspect-btn').addEventListener('click', () => {
                    this.openInspectByGenes(c.gene1, c.gene2);
                });
                tr.querySelector('.tissue-btn').addEventListener('click', () => {
                    this.openByTissueByGenes(c.gene1, c.gene2);
                });
                tr.querySelector('.hotspot-btn').addEventListener('click', () => {
                    this.openInspectWithHotspot(c.gene1, c.gene2);
                });
                tbody.appendChild(tr);
            });

        // Attach gene tooltip handlers
        this.attachGeneTooltips(tbody);
    }

    openInspectWithHotspot(gene1, gene2) {
        // Open the inspect modal and pre-select hotspot overlay
        this.openInspectByGenes(gene1, gene2);
        // After opening, switch to hotspot overlay mode
        setTimeout(() => {
            const hotspotSelect = document.getElementById('hotspotGene');
            const hotspotMode = document.getElementById('hotspotMode');
            if (hotspotSelect && hotspotSelect.options.length > 1) {
                hotspotSelect.selectedIndex = 1; // Select first available hotspot gene
                if (hotspotMode) hotspotMode.value = 'color';
                hotspotSelect.dispatchEvent(new Event('change'));
            }
        }, 300);
    }

    attachGeneTooltips(container) {
        container.querySelectorAll('.gene-hover').forEach(el => {
            el.addEventListener('mouseenter', (e) => {
                this._tooltipTimer = setTimeout(() => {
                    this.showGeneTooltip(e, el.dataset.gene);
                }, 400);
            });
            el.addEventListener('mouseleave', () => {
                clearTimeout(this._tooltipTimer);
                this.hideGeneTooltip();
            });
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

        // Check if there are any uncorrelated genes
        const hasUncorrelated = this.results.clusters.some(c => c.hasCorrelation === false);

        // Build header based on what data we have
        let headerCells = `
            <th data-sort="gene">Gene</th>
            <th data-sort="cluster">Cluster</th>
        `;
        if (hasUncorrelated) {
            headerCells += `<th data-sort="hasCorrelation">Corr</th>`;
        }
        headerCells += `
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

        // Add actions column
        headerCells += `<th style="text-align: center;">Analyze</th>`;

        thead.innerHTML = `<tr>${headerCells}</tr>`;

        // Re-attach sorting event listeners (Ctrl/Cmd+click to copy column)
        thead.querySelectorAll('th[data-sort]').forEach(th => {
            th.addEventListener('click', (e) => {
                if (e.ctrlKey || e.metaKey) {
                    const colIndex = Array.from(th.parentNode.children).indexOf(th);
                    this.copyColumnToClipboard(th.closest('table'), colIndex);
                } else {
                    this.sortTable(th);
                }
            });
            th.title = 'Click to sort. Ctrl/Cmd+click to copy column.';
        });

        this.results.clusters
            .sort((a, b) => {
                // Sort correlated genes first (by cluster), then uncorrelated
                const aCorr = a.hasCorrelation === false ? 1 : 0;
                const bCorr = b.hasCorrelation === false ? 1 : 0;
                if (aCorr !== bCorr) return aCorr - bCorr;
                const aCluster = a.cluster === '-' ? Infinity : a.cluster;
                const bCluster = b.cluster === '-' ? Infinity : b.cluster;
                return aCluster - bCluster || a.gene.localeCompare(b.gene);
            })
            .forEach(c => {
                const tr = document.createElement('tr');
                if (c.hasCorrelation === false) {
                    tr.style.opacity = '0.7';
                }
                const geneStat = this.geneStats?.get(c.gene);

                let rowHtml = `
                    <td class="gene-hover" data-gene="${c.gene}" style="cursor: help;">${c.gene}${c.inGeneList && this.results.mode === 'design' ? '*' : ''}</td>
                    <td>${c.cluster}</td>
                `;
                if (hasUncorrelated) {
                    rowHtml += `<td style="text-align: center; color: ${c.hasCorrelation === false ? '#dc2626' : '#16a34a'}; font-weight: 600;">${c.hasCorrelation === false ? 'No' : 'Yes'}</td>`;
                }
                rowHtml += `
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

                // Add analyze buttons
                rowHtml += `
                    <td style="text-align: center; white-space: nowrap;">
                        <button class="btn btn-sm gene-effect-btn" style="padding: 2px 6px; font-size: 10px; background: #5a9f4a; color: white;" data-gene="${c.gene}">Gene Effect</button>
                        <button class="btn btn-sm tissue-btn" style="padding: 2px 6px; font-size: 10px; margin-left: 4px; background: #6b7280; color: white;" data-gene="${c.gene}">By Tissue</button>
                        <button class="btn btn-sm hotspot-btn" style="padding: 2px 6px; font-size: 10px; margin-left: 4px; background: #6b7280; color: white;" data-gene="${c.gene}">By Hotspot</button>
                    </td>
                `;

                tr.innerHTML = rowHtml;
                tbody.appendChild(tr);
            });

        // Add event listeners to buttons
        tbody.querySelectorAll('.gene-effect-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                this.openGeneEffectModal(btn.dataset.gene, 'tissue');
                this._applyParamFiltersToGEModal();
            });
        });
        tbody.querySelectorAll('.tissue-btn').forEach(btn => {
            btn.addEventListener('click', () => this.openGeneEffectModal(btn.dataset.gene, 'tissue'));
        });
        tbody.querySelectorAll('.hotspot-btn').forEach(btn => {
            btn.addEventListener('click', () => this.openGeneEffectModal(btn.dataset.gene, 'hotspot'));
        });

        // Attach gene tooltips
        this.attachGeneTooltips(tbody);
    }

    displaySummary() {
        const text = document.getElementById('summaryText');
        const lineage = document.getElementById('lineageFilter').value || 'All lineages';
        const subtype = document.getElementById('subLineageFilter')?.value;

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

        // Build lineage filter text
        let lineageText = lineage;
        if (subtype) {
            lineageText += ` (${subtype})`;
        }

        // Format date and time
        const now = new Date();
        const dateTimeStr = now.toLocaleString('en-US', {
            year: 'numeric',
            month: 'short',
            day: 'numeric',
            hour: '2-digit',
            minute: '2-digit',
            hour12: false
        });

        // Build excluded tissues section
        let excludedTissuesText = 'None';
        if (this.excludedTissues && this.excludedTissues.size > 0) {
            excludedTissuesText = [...this.excludedTissues].sort().join(', ');
        }

        // Build hotspot/translocation filter section
        let hotspotFilterText = '';
        const paramHotspotGene = document.getElementById('paramHotspotGene')?.value;
        const paramHotspotLevel = document.getElementById('paramHotspotLevel')?.value;
        if (paramHotspotGene) {
            hotspotFilterText = `\nHotspot Mutation Filter: ${paramHotspotGene} (${paramHotspotLevel || 'all'})`;
        }
        const paramTranslocGene = document.getElementById('paramTranslocationGene')?.value;
        const paramTranslocLevel = document.getElementById('paramTranslocationLevel')?.value;
        if (paramTranslocGene) {
            hotspotFilterText += `\nTranslocation Filter: ${paramTranslocGene} (${paramTranslocLevel || 'all'})`;
        }

        // P-value threshold (for mutation mode)
        let pValueText = '';
        if (this.results.mode === 'mutation') {
            pValueText = `\nP-value Threshold: ${document.getElementById('pValueThreshold')?.value || '0.001'}`;
        }

        text.textContent = `Gene Correlation Analysis Summary
================================
Run: ${dateTimeStr}

Analysis Mode: ${this.results.mode === 'analysis' ? 'Analysis (within gene list)' : this.results.mode === 'design' ? 'Design (find correlated genes)' : this.results.mode === 'mutation' ? 'Mutation Analysis' : this.results.mode}
Correlation Cutoff: ${this.results.cutoff}
Minimum Cell Lines: ${document.getElementById('minCellLines').value}
Minimum Slope: ${document.getElementById('minSlope').value}
Lineage Filter: ${lineageText}
Excluded Tissues: ${excludedTissuesText}${hotspotFilterText}${pValueText}

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
        const numericColumns = [
            'correlation', 'slope', 'n', 'cluster',
            'meanEffect', 'sdEffect', 'meanEffectFiltered', 'sdEffectFiltered',
            'lfc', 'fdr'
        ];
        const isNumeric = numericColumns.includes(sortKey);

        const currentDir = th.dataset.dir || 'asc';
        const newDir = currentDir === 'asc' ? 'desc' : 'asc';
        th.dataset.dir = newDir;

        rows.sort((a, b) => {
            let aVal = a.children[colIndex]?.textContent?.trim() || '';
            let bVal = b.children[colIndex]?.textContent?.trim() || '';

            // Remove any asterisks or other markers from gene names
            aVal = aVal.replace(/\*$/, '');
            bVal = bVal.replace(/\*$/, '');

            if (isNumeric) {
                const aNum = parseFloat(aVal);
                const bNum = parseFloat(bVal);
                // Handle NaN values - put them at the end
                if (isNaN(aNum) && isNaN(bNum)) return 0;
                if (isNaN(aNum)) return 1;
                if (isNaN(bNum)) return -1;
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
            csv += `# Analysis mode: ${this.results.mode === 'design' ? 'Design (find correlated genes)' : 'Analysis (within gene list)'}\n`;
            csv += `# Lineage filter: ${lineage}\n`;
            if (subLineage) csv += `# Subtype filter: ${subLineage}\n`;
            csv += `# Filtered cell lines: ${this.results.nCellLines}\n`;
            csv += `# Date: ${new Date().toISOString().slice(0, 10)}\n`;
            if (this.results.mode === 'design') {
                csv += `# Gene_Type: Input = user-provided gene, Correlated = found by correlation analysis\n`;
            }
            csv += '#\n';

            const isDesignMode = this.results.mode === 'design';
            const geneList = this.results.geneList || [];
            const hasStats = this.geneStats && this.geneStats.size > 0;

            // Build header matching screen columns
            let header = isDesignMode ? 'Gene,Gene_Type,' : 'Gene,';
            header += 'Cluster,Has_Correlation,Mean_Effect_All,SD_Effect_All';
            if (isFiltered) header += ',Mean_Effect_Filtered,SD_Effect_Filtered';
            if (hasStats) header += ',LFC,FDR';
            csv += header + '\n';

            this.results.clusters.forEach(c => {
                const geneType = geneList.includes(c.gene) ? 'Input' : 'Correlated';
                const hasCorr = c.hasCorrelation === false ? 'No' : 'Yes';
                let row = isDesignMode ? `${c.gene},${geneType},` : `${c.gene},`;
                row += `${c.cluster},${hasCorr},${c.meanEffect},${c.sdEffect}`;
                if (isFiltered) row += `,${c.meanEffectFiltered},${c.sdEffectFiltered}`;
                if (hasStats) {
                    const geneStat = this.geneStats.get(c.gene);
                    const lfc = geneStat?.lfc !== null && geneStat?.lfc !== undefined ? geneStat.lfc.toFixed(4) : '';
                    const fdr = geneStat?.fdr !== null && geneStat?.fdr !== undefined ? geneStat.fdr.toExponential(2) : '';
                    row += `,${lfc},${fdr}`;
                }
                csv += row + '\n';
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

    _getNetworkFilterText() {
        const parts = [];
        const lineage = document.getElementById('lineageFilter')?.value;
        const subLineage = document.getElementById('subLineageFilter')?.value;
        if (lineage) {
            parts.push(subLineage ? `${lineage} / ${subLineage}` : lineage);
        }
        if (this.excludedTissues && this.excludedTissues.size > 0) {
            parts.push(`${this.excludedTissues.size} tissue${this.excludedTissues.size > 1 ? 's' : ''} excluded`);
        }
        const hotspotGene = document.getElementById('paramHotspotGene')?.value;
        const hotspotLevel = document.getElementById('paramHotspotLevel')?.value;
        if (hotspotGene) {
            const levelLabel = hotspotLevel === '1+2' ? 'mut' : hotspotLevel === '0' ? 'WT' : `level ${hotspotLevel}`;
            parts.push(`${hotspotGene} ${levelLabel}`);
        }
        const translocGene = document.getElementById('paramTranslocationGene')?.value;
        const translocLevel = document.getElementById('paramTranslocationLevel')?.value;
        if (translocGene) {
            const levelLabel = translocLevel === '1+2' ? 'fused' : translocLevel === '0' ? 'not fused' : `level ${translocLevel}`;
            parts.push(`${translocGene} ${levelLabel}`);
        }
        // Include oncoprint multi-gene filters
        if (this._activeOncoprintFilters && this._activeOncoprintFilters.length > 0) {
            const shown = new Set();
            if (hotspotGene) shown.add(hotspotGene);
            if (translocGene) shown.add(translocGene);
            for (const f of this._activeOncoprintFilters) {
                if (!shown.has(f.gene)) {
                    parts.push(`${f.gene} ${f.state === 'mut' ? 'Mut' : 'WT'}`);
                    shown.add(f.gene);
                }
            }
        }
        if (parts.length === 0) return null;
        return `Filters: ${parts.join('  \u00b7  ')}  \u00b7  n=${this.results.nCellLines}`;
    }

    downloadNetworkPNG() {
        if (!this.network) return;

        const networkCanvas = document.querySelector('#networkPlot canvas');
        if (!networkCanvas) {
            console.error('Network canvas not found');
            return;
        }

        // Create a high-resolution canvas for publication quality (300 DPI)
        const pngScale = 2;
        const container = document.getElementById('networkPlot');
        const cssWidth = container.clientWidth;
        const cssHeight = container.clientHeight;
        const legendHeight = 160;
        const padding = 30;
        const filterText = this._getNetworkFilterText();
        const svgBannerFs = this._netBannerFontSize || 12; const bannerFs = svgBannerFs; const filterBannerHeight = filterText ? bannerFs + 14 : 0;
        const totalWidth = cssWidth;
        const totalHeight = filterBannerHeight + cssHeight + legendHeight + padding;

        const transparentBg = document.getElementById('exportNetworkTransparentBg')?.checked;

        const canvas = document.createElement('canvas');
        canvas.width = totalWidth * pngScale;
        canvas.height = totalHeight * pngScale;
        const ctx = canvas.getContext('2d');
        ctx.scale(pngScale, pngScale);

        // Draw background (skip for transparent)
        if (!transparentBg) {
            ctx.fillStyle = 'white';
            ctx.fillRect(0, 0, totalWidth, totalHeight);
        }

        // Draw filter banner at top if active
        if (filterText) {
            ctx.font = `${bannerFs}px Arial`;
            ctx.fillStyle = '#374151';
            if (this._netBannerPos) {
                ctx.textAlign = 'left';
                ctx.fillText(filterText, this._netBannerPos.x, this._netBannerPos.y + bannerFs);
            } else {
                ctx.textAlign = 'center';
                ctx.fillText(filterText, totalWidth / 2, bannerFs + 4);
            }
            ctx.textAlign = 'left';
        }

        // Draw network scaled from canvas pixels to CSS dimensions
        ctx.drawImage(networkCanvas, 0, filterBannerHeight, cssWidth, cssHeight);

        // Draw legend background
        const legendTop = filterBannerHeight + cssHeight;
        if (!transparentBg) {
            ctx.fillStyle = '#f9fafb';
            ctx.strokeStyle = '#e5e7eb';
            ctx.lineWidth = 1;
            ctx.fillRect(15, legendTop + 10, totalWidth - 30, legendHeight - 10);
            ctx.strokeRect(15, legendTop + 10, totalWidth - 30, legendHeight - 10);
        }

        // Draw legend
        const legendY = legendTop + padding + 10;
        const titleFont = 'bold 16px Arial';
        const textFont = '14px Arial';
        const smallFont = '13px Arial';

        // Calculate total legend width to center it
        let totalLegendWidth = 160 + 160; // Correlation + Edge Thickness
        if (this.results?.mode === 'design') totalLegendWidth += 140;
        if (document.getElementById('colorByGeneEffect').checked && this.results?.clusters) totalLegendWidth += 170;
        if (document.getElementById('colorByStats').checked && this.geneStats && this.geneStats.size > 0) totalLegendWidth += 200;

        let legendX = Math.max(40, (totalWidth - totalLegendWidth) / 2);

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
            ctx.fillStyle = '#5a9f4a';
            ctx.beginPath();
            ctx.arc(legendX + 12, legendY + 25, 10, 0, Math.PI * 2);
            ctx.fill();
            ctx.fillStyle = '#333';
            ctx.fillText('Input', legendX + 28, legendY + 30);

            // Correlated gene
            ctx.fillStyle = '#a8d89a';
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
                ctx.fillText('Gene Effect', legendX, gradY - 4);
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

        // Embed metadata and create download link
        const meta = this._buildExportMetadata('network', {
            geneList: this.getGeneList(),
            mode: this.results?.mode,
            cutoff: this.results?.cutoff,
            nCellLines: this.results?.nCellLines,
            networkSettings: this._captureNetworkSettings(),
            oncoprintFilters: this._activeOncoprintFilters || null
        });
        const dataURL = canvas.toDataURL('image/png');
        fetch(dataURL).then(r => r.arrayBuffer()).then(buf => {
            const pngWithMeta = this._addPngTextChunk(buf, 'correlate-meta', JSON.stringify(meta));
            const blob = new Blob([pngWithMeta], { type: 'image/png' });
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = 'correlation_network.png';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(a.href);
        });
    }

    async downloadNetworkSVG() {
        if (!this.network || !this.networkData) return;

        const transparentBg = document.getElementById('exportNetworkTransparentBg')?.checked;

        // Build SVG from network data
        const container = document.getElementById('networkPlot');
        const width = container.clientWidth;
        const networkHeight = container.clientHeight;
        const legendHeight = 160;  // Larger for publication
        const filterText = this._getNetworkFilterText();
        const svgBannerFs = this._netBannerFontSize || 12; const bannerFs = svgBannerFs; const filterBannerHeight = filterText ? bannerFs + 14 : 0;
        const totalHeight = filterBannerHeight + networkHeight + legendHeight;

        // Get positions from vis.js and convert to DOM coordinates
        const positions = this.network.getPositions();
        const domPositions = {};
        for (const nodeId in positions) {
            const canvasPos = positions[nodeId];
            const domPos = this.network.canvasToDOM({ x: canvasPos.x, y: canvasPos.y });
            domPositions[nodeId] = { x: domPos.x, y: domPos.y + filterBannerHeight };
        }

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
${transparentBg ? '' : '<rect width="100%" height="100%" fill="white"/>'}
${filterText ? `<text x="${this._netBannerPos ? this._netBannerPos.x : width / 2}" y="${this._netBannerPos ? this._netBannerPos.y + svgBannerFs : svgBannerFs + 4}" text-anchor="${this._netBannerPos ? 'start' : 'middle'}" style="font-family: Arial, sans-serif; font-size: ${svgBannerFs}px; fill: #374151;">${this.escapeXml(filterText)}</text>` : ''}
`;

        // Get current scale for sizing elements
        const scale = this.network.getScale();

        // Draw edges — read stored width so SVG matches vis.js rendering
        this.networkData.edges.forEach(edge => {
            const from = domPositions[edge.from];
            const to = domPositions[edge.to];
            if (from && to) {
                const color = typeof edge.color === 'object' ? (edge.color?.color || '#3182ce') : (edge.color || '#3182ce');
                const strokeWidth = (edge.width || 1) * scale;
                svg += `  <line x1="${from.x}" y1="${from.y}" x2="${to.x}" y2="${to.y}" stroke="${color}" stroke-width="${strokeWidth}" opacity="0.8"/>\n`;
            }
        });

        // Draw nodes — read stored size/font so SVG matches vis.js rendering
        this.networkData.nodes.forEach(node => {
            const pos = domPositions[node.id];
            if (pos) {
                const bgColor = node.color?.background || '#5a9f4a';
                const nodeRadius = (node.size || 25) * scale;
                const nodeFontSize = (node.font?.size || 16) * scale;
                svg += `  <circle cx="${pos.x}" cy="${pos.y}" r="${nodeRadius}" fill="${bgColor}" stroke="white" stroke-width="${2 * scale}"/>\n`;

                // Handle multi-line labels
                const labelLines = (node.label || node.id).split('\n');
                labelLines.forEach((line, i) => {
                    const yOffset = pos.y + nodeRadius + 14 * scale + (i * nodeFontSize);
                    svg += `  <text x="${pos.x}" y="${yOffset}" text-anchor="middle" style="font-family: Arial; font-size: ${nodeFontSize}px; fill: ${node.font?.color || '#333'};">${this.escapeXml(line)}</text>\n`;
                });
            }
        });

        // Draw legend - LARGER for publication
        const legendTop = filterBannerHeight + networkHeight;
        const legendY = legendTop + 35;

        // Calculate total legend width to center it
        let totalLegendWidth = 160 + 160; // Correlation + Edge Thickness
        if (this.results?.mode === 'design') totalLegendWidth += 140;
        if (document.getElementById('colorByGeneEffect').checked && this.results?.clusters) totalLegendWidth += 170;
        if (document.getElementById('colorByStats').checked && this.geneStats && this.geneStats.size > 0) totalLegendWidth += 200;

        let legendX = Math.max(40, (width - totalLegendWidth) / 2);

        // Legend background
        if (!transparentBg) {
            svg += `  <rect x="15" y="${legendTop + 10}" width="${width - 30}" height="145" fill="#f9fafb" stroke="#e5e7eb" rx="4"/>\n`;
        }

        // Correlation legend
        svg += `  <text x="${legendX}" y="${legendY}" class="legend-title">Correlation:</text>\n`;
        svg += `  <line x1="${legendX}" y1="${legendY + 22}" x2="${legendX + 35}" y2="${legendY + 22}" stroke="#3182ce" stroke-width="4"/>\n`;
        svg += `  <text x="${legendX + 42}" y="${legendY + 27}" class="legend-text">Positive</text>\n`;
        svg += `  <line x1="${legendX}" y1="${legendY + 48}" x2="${legendX + 35}" y2="${legendY + 48}" stroke="#e53e3e" stroke-width="4"/>\n`;
        svg += `  <text x="${legendX + 42}" y="${legendY + 53}" class="legend-text">Negative</text>\n`;

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
            svg += `  <circle cx="${legendX + 12}" cy="${legendY + 25}" r="10" fill="#5a9f4a"/>\n`;
            svg += `  <text x="${legendX + 28}" y="${legendY + 30}" class="legend-text">Input</text>\n`;
            svg += `  <circle cx="${legendX + 12}" cy="${legendY + 52}" r="10" fill="#a8d89a"/>\n`;
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

        // Embed metadata
        const meta = this._buildExportMetadata('network', {
            geneList: this.getGeneList(),
            mode: this.results?.mode,
            cutoff: this.results?.cutoff,
            nCellLines: this.results?.nCellLines,
            networkSettings: this._captureNetworkSettings(),
            oncoprintFilters: this._activeOncoprintFilters || null
        });
        svg += `<metadata><correlate-meta>${JSON.stringify(meta)}</correlate-meta></metadata>`;
        svg += '</svg>';

        // Sanitize for Illustrator compatibility + optional text outlining
        svg = await this._finalizeSvgForExport(svg);

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
            this.updateRemovedNodesList();
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
        this.updateRemovedNodesList();
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

            // Use visible network genes for scale (not all clusters)
            const visibleEffects = [];
            this.networkData.nodes.forEach(node => {
                const effect = effectMap.get(node.id);
                if (effect !== undefined && !isNaN(effect)) visibleEffects.push(effect);
            });
            const effectValues = visibleEffects.length > 0 ? visibleEffects : this.results.clusters.map(c => c.meanEffect).filter(v => !isNaN(v));

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
                        color: { background: bgColor, border: '#000000' }
                    });
                });

                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">Gene Effect</div>
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
                        color: { background: bgColor, border: '#000000' }
                    });
                });

                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">Absolute Gene Effect</div>
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
                            (isInput ? '#5a9f4a' : '#a8d89a') : '#5a9f4a',
                        border: '#000000'
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
                        color: { background: bgColor, border: '#000000' }
                    });
                });

                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">LFC</div>
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
                        color: { background: bgColor, border: '#000000' }
                    });
                });

                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">Absolute LFC</div>
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
                        color: { background: bgColor, border: '#000000' }
                    });
                });

                if (colorLegend) colorLegend.innerHTML = `
                    <div class="legend-item">FDR</div>
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
            btn.textContent = 'Lock';
            btn.classList.remove('btn-active');
        } else {
            btn.textContent = 'Unlock';
            btn.classList.add('btn-active');
        }
    }

    toggleRemoveMode() {
        // Deactivate select mode if active (without re-entering toggleRemoveMode)
        if (this.selectMode) {
            this.selectMode = false;
            this.restoreSelectedNodeColors();
            this.selectedNodes.clear();
            this.updateSelectedNodesList();
            const selBtn = document.getElementById('toggleSelectMode');
            if (selBtn) { selBtn.classList.remove('btn-active'); selBtn.style.backgroundColor = ''; selBtn.style.borderColor = ''; selBtn.style.color = ''; }
        }

        this.removeMode = !this.removeMode;
        const btn = document.getElementById('toggleRemoveMode');
        const helpText = document.getElementById('networkHelpText');
        if (this.removeMode) {
            btn.classList.add('btn-active');
            btn.style.backgroundColor = '#dc2626';
            btn.style.borderColor = '#dc2626';
            btn.style.color = 'white';
            if (helpText) helpText.textContent = 'Click node to remove from network and gene list';
        } else {
            btn.classList.remove('btn-active');
            btn.style.backgroundColor = '';
            btn.style.borderColor = '';
            btn.style.color = '';
            if (helpText) helpText.textContent = 'Double-click node for Gene Effect, edge for Correlation';
        }
    }

    toggleSelectMode() {
        // Deactivate remove mode if active (without re-entering toggleSelectMode)
        if (this.removeMode) {
            this.removeMode = false;
            const rmBtn = document.getElementById('toggleRemoveMode');
            if (rmBtn) { rmBtn.classList.remove('btn-active'); rmBtn.style.backgroundColor = ''; rmBtn.style.borderColor = ''; rmBtn.style.color = ''; }
        }

        this.selectMode = !this.selectMode;
        const btn = document.getElementById('toggleSelectMode');
        const helpText = document.getElementById('networkHelpText');
        if (this.selectMode) {
            this.selectedNodes.clear();
            btn.classList.add('btn-active');
            btn.style.backgroundColor = '#2563eb';
            btn.style.borderColor = '#2563eb';
            btn.style.color = 'white';
            if (helpText) helpText.textContent = 'Click nodes to select them for the gene list';
        } else {
            // Restore original node colors for any selected nodes
            this.restoreSelectedNodeColors();
            this.selectedNodes.clear();
            this.updateSelectedNodesList();
            btn.classList.remove('btn-active');
            btn.style.backgroundColor = '';
            btn.style.borderColor = '';
            btn.style.color = '';
            if (helpText) helpText.textContent = 'Double-click node for Gene Effect, edge for Correlation';
        }
    }

    restoreSelectedNodeColors() {
        if (!this.networkData) return;
        this.selectedNodes.forEach(nodeId => {
            const isInput = this.results?.geneList?.includes(nodeId);
            this.networkData.nodes.update({
                id: nodeId,
                color: {
                    background: this.results?.mode === 'design' ?
                        (isInput ? '#5a9f4a' : '#a8d89a') : '#5a9f4a',
                    border: '#000000'
                }
            });
        });
    }

    updateSelectedNodesList() {
        const listEl = document.getElementById('selectedNodesList');
        const textEl = document.getElementById('selectedNodesText');
        if (!listEl || !textEl) return;

        if (this.selectedNodes.size > 0) {
            textEl.textContent = Array.from(this.selectedNodes).join(', ');
            listEl.style.display = 'block';
        } else {
            listEl.style.display = 'none';
        }
    }

    clearSelectedNodes() {
        this.restoreSelectedNodeColors();
        this.selectedNodes.clear();
        this.updateSelectedNodesList();
        // Restore the original gene list in textarea
        if (this.results?.geneList) {
            document.getElementById('geneTextarea').value = this.results.geneList.join('\n');
            this.updateGeneCount();
        }
    }

    updateRemovedNodesList() {
        const listEl = document.getElementById('removedNodesList');
        const textEl = document.getElementById('removedNodesText');
        if (!listEl || !textEl) return;

        if (this.hiddenNodes && this.hiddenNodes.length > 0) {
            // Create clickable list of removed nodes
            const nodeLinks = this.hiddenNodes.map((node, idx) =>
                `<span class="restore-node" data-idx="${idx}" style="cursor: pointer; text-decoration: underline; color: #2563eb;">${node.id}</span>`
            ).join(', ');
            textEl.innerHTML = nodeLinks;
            listEl.style.display = 'block';

            // Add click handlers to restore individual nodes
            listEl.querySelectorAll('.restore-node').forEach(el => {
                el.addEventListener('click', (e) => {
                    const idx = parseInt(e.target.dataset.idx);
                    this.restoreNode(idx);
                });
            });
        } else {
            listEl.style.display = 'none';
        }
    }

    restoreNode(idx) {
        if (!this.hiddenNodes || idx < 0 || idx >= this.hiddenNodes.length) return;
        const node = this.hiddenNodes.splice(idx, 1)[0];
        if (node && this.networkData) {
            this.networkData.nodes.add(node);
        }
        this.updateRemovedNodesList();
    }

    changeNetworkLayout() {
        if (!this.network) return;

        // Cycle through layout options: Default (Barnes Hut), Force Atlas, Hierarchical
        const layouts = ['default', 'forceAtlas2Based', 'hierarchical'];
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
        } else if (layoutName === 'default') {
            // Restore original Barnes Hut with original stabilization
            const nodeCount = this.networkData?.nodes?.length || 50;
            const stabilizationIterations = nodeCount > 50 ? 300 : 150;
            options = {
                layout: { hierarchical: { enabled: false } },
                physics: {
                    enabled: true,
                    solver: 'barnesHut',
                    stabilization: { iterations: stabilizationIterations }
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
            btn.textContent = 'Lock';
            btn.classList.remove('btn-active');
        } else {
            btn.textContent = 'Unlock';
            btn.classList.add('btn-active');
        }

        // Show which layout is active
        const layoutBtn = document.getElementById('changeLayout');
        const layoutNames = { 'default': 'Default', 'forceAtlas2Based': 'Force Atlas', 'hierarchical': 'Hierarchical' };
        layoutBtn.textContent = layoutNames[layoutName];

        // Re-center after layout change
        setTimeout(() => this.network?.fit(), 500);
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
        if (layoutBtn) layoutBtn.textContent = 'Layout';
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

        // Create clusters CSV (matching screen columns)
        const hasStats = this.geneStats && this.geneStats.size > 0;
        const isDesignMode = this.results.mode === 'design';
        let clustersHeader = isDesignMode ? 'Gene,Gene_Type,' : 'Gene,';
        clustersHeader += 'Cluster,Has_Correlation,Mean_Effect_All,SD_Effect_All';
        if (this.results.isFiltered) clustersHeader += ',Mean_Effect_Filtered,SD_Effect_Filtered';
        if (hasStats) clustersHeader += ',LFC,FDR';
        let clustersCSV = clustersHeader + '\n';
        const geneList = this.results.geneList || [];
        this.results.clusters.forEach(c => {
            const geneType = geneList.includes(c.gene) ? 'Input' : 'Correlated';
            const hasCorr = c.hasCorrelation === false ? 'No' : 'Yes';
            let row = isDesignMode ? `${c.gene},${geneType},` : `${c.gene},`;
            row += `${c.cluster},${hasCorr},${c.meanEffect},${c.sdEffect}`;
            if (this.results.isFiltered) row += `,${c.meanEffectFiltered},${c.sdEffectFiltered}`;
            if (hasStats) {
                const geneStat = this.geneStats.get(c.gene);
                const lfc = geneStat?.lfc !== null && geneStat?.lfc !== undefined ? geneStat.lfc.toFixed(4) : '';
                const fdr = geneStat?.fdr !== null && geneStat?.fdr !== undefined ? geneStat.fdr.toExponential(2) : '';
                row += `,${lfc},${fdr}`;
            }
            clustersCSV += row + '\n';
        });

        // Create summary
        const lineage = document.getElementById('lineageFilter').value || 'All lineages';
        const subtype = document.getElementById('subLineageFilter')?.value;
        let lineageText = lineage;
        if (subtype) {
            lineageText += ` (${subtype})`;
        }
        const summary = `Gene Correlation Analysis Summary
================================
Analysis Mode: ${this.results.mode === 'analysis' ? 'Analysis (within gene list)' : 'Design (find correlated genes)'}
Correlation Cutoff: ${this.results.cutoff}
Minimum Cell Lines: ${document.getElementById('minCellLines').value}
Minimum Slope: ${document.getElementById('minSlope').value}
Lineage Filter: ${lineageText}

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

        // Calculate total legend width to center it
        let totalLegendWidth = 160 + 160; // Correlation + Edge Thickness
        if (this.results?.mode === 'design') totalLegendWidth += 140;
        if (document.getElementById('colorByGeneEffect').checked && this.results?.clusters) totalLegendWidth += 170;
        if (document.getElementById('colorByStats').checked && this.geneStats && this.geneStats.size > 0) totalLegendWidth += 200;

        let legendX = Math.max(40, (totalWidth - totalLegendWidth) / 2);

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
            ctx.fillStyle = '#5a9f4a';
            ctx.beginPath();
            ctx.arc(legendX + 12, legendY + 25, 10, 0, Math.PI * 2);
            ctx.fill();
            ctx.fillStyle = '#333';
            ctx.fillText('Input', legendX + 28, legendY + 30);

            // Correlated gene
            ctx.fillStyle = '#a8d89a';
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
                ctx.fillText('Gene Effect', legendX, gradY - 4);
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

            const gradientWidth = 120;
            const gradientHeight = 18;
            const gradY = legendY + 18;

            if (colorStatType === 'signed_lfc') {
                const lfcValues = stats.map(s => s.lfc).filter(v => v !== null && !isNaN(v));
                const minLfc = Math.min(...lfcValues);
                const maxLfc = Math.max(...lfcValues);

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
                ctx.fillText(minLfc.toFixed(1), legendX, gradY + gradientHeight + 16);
                ctx.fillText(`LFC (+/−)${scaleLabel}`, legendX, gradY - 4);
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
                ctx.fillText(`|LFC|${scaleLabel}`, legendX, gradY - 4);
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
                ctx.fillText(`FDR${scaleLabel}`, legendX, gradY - 4);
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

        return canvas.toDataURL('image/png');
    }

    // Helper to get network SVG data for ZIP
    getNetworkSVGData() {
        if (!this.network || !this.networkData) return null;

        const container = document.getElementById('networkPlot');
        const width = container.clientWidth;
        const networkHeight = container.clientHeight;
        const legendHeight = 160;  // Larger for publication
        const filterText = this._getNetworkFilterText();
        const svgBannerFs = this._netBannerFontSize || 12; const bannerFs = svgBannerFs; const filterBannerHeight = filterText ? bannerFs + 14 : 0;
        const totalHeight = filterBannerHeight + networkHeight + legendHeight;

        // Get positions from vis.js and convert to DOM coordinates
        const positions = this.network.getPositions();
        const domPositions = {};
        for (const nodeId in positions) {
            const canvasPos = positions[nodeId];
            const domPos = this.network.canvasToDOM({ x: canvasPos.x, y: canvasPos.y });
            domPositions[nodeId] = { x: domPos.x, y: domPos.y + filterBannerHeight };
        }

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
${filterText ? `<text x="${this._netBannerPos ? this._netBannerPos.x : width / 2}" y="${this._netBannerPos ? this._netBannerPos.y + svgBannerFs : svgBannerFs + 4}" text-anchor="${this._netBannerPos ? 'start' : 'middle'}" style="font-family: Arial, sans-serif; font-size: ${svgBannerFs}px; fill: #374151;">${this.escapeXml(filterText)}</text>` : ''}
`;

        // Get current scale for sizing elements
        const scale = this.network.getScale();

        // Draw edges — read stored width so SVG matches vis.js rendering
        this.networkData.edges.forEach(edge => {
            const from = domPositions[edge.from];
            const to = domPositions[edge.to];
            if (from && to) {
                const color = typeof edge.color === 'object' ? (edge.color?.color || '#3182ce') : (edge.color || '#3182ce');
                const strokeWidth = (edge.width || 1) * scale;
                svg += `  <line x1="${from.x}" y1="${from.y}" x2="${to.x}" y2="${to.y}" stroke="${color}" stroke-width="${strokeWidth}" opacity="0.8"/>\n`;
            }
        });

        // Draw nodes — read stored size/font so SVG matches vis.js rendering
        this.networkData.nodes.forEach(node => {
            const pos = domPositions[node.id];
            if (pos) {
                const bgColor = node.color?.background || '#5a9f4a';
                const nodeRadius = (node.size || 25) * scale;
                const nodeFontSize = (node.font?.size || 16) * scale;
                svg += `  <circle cx="${pos.x}" cy="${pos.y}" r="${nodeRadius}" fill="${bgColor}" stroke="white" stroke-width="${2 * scale}"/>\n`;

                // Handle multi-line labels
                const labelLines = (node.label || node.id).split('\n');
                labelLines.forEach((line, i) => {
                    const yOffset = pos.y + nodeRadius + 14 * scale + (i * nodeFontSize);
                    svg += `  <text x="${pos.x}" y="${yOffset}" text-anchor="middle" style="font-family: Arial; font-size: ${nodeFontSize}px; fill: ${node.font?.color || '#333'};">${this.escapeXml(line)}</text>\n`;
                });
            }
        });

        // Draw legend - LARGER for publication
        const legendTopZip = filterBannerHeight + networkHeight;
        const legendY = legendTopZip + 35;

        // Calculate total legend width to center it
        let totalLegendWidth = 160 + 160; // Correlation + Edge Thickness
        if (this.results?.mode === 'design') totalLegendWidth += 140;
        if (document.getElementById('colorByGeneEffect').checked && this.results?.clusters) totalLegendWidth += 170;
        if (document.getElementById('colorByStats').checked && this.geneStats && this.geneStats.size > 0) totalLegendWidth += 200;

        let legendX = Math.max(40, (width - totalLegendWidth) / 2);

        // Legend background
        svg += `  <rect x="15" y="${legendTopZip + 10}" width="${width - 30}" height="145" fill="#f9fafb" stroke="#e5e7eb" rx="4"/>\n`;

        // Correlation legend
        svg += `  <text x="${legendX}" y="${legendY}" class="legend-title">Correlation:</text>\n`;
        svg += `  <line x1="${legendX}" y1="${legendY + 22}" x2="${legendX + 35}" y2="${legendY + 22}" stroke="#3182ce" stroke-width="4"/>\n`;
        svg += `  <text x="${legendX + 42}" y="${legendY + 27}" class="legend-text">Positive</text>\n`;
        svg += `  <line x1="${legendX}" y1="${legendY + 48}" x2="${legendX + 35}" y2="${legendY + 48}" stroke="#e53e3e" stroke-width="4"/>\n`;
        svg += `  <text x="${legendX + 42}" y="${legendY + 53}" class="legend-text">Negative</text>\n`;

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
            svg += `  <circle cx="${legendX + 12}" cy="${legendY + 25}" r="10" fill="#5a9f4a"/>\n`;
            svg += `  <text x="${legendX + 28}" y="${legendY + 30}" class="legend-text">Input</text>\n`;
            svg += `  <circle cx="${legendX + 12}" cy="${legendY + 52}" r="10" fill="#a8d89a"/>\n`;
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
        this._userLegendPosition = null;
        this._userTitlePosition = null;
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

        // Calculate axis limits (needed if user clicks tissue to go to Inspect)
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
        this._userLegendPosition = null;
        this._userTitlePosition = null;
        this._userLabelPositions = new Map();
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

        // Set gene inputs
        document.getElementById('inspectGeneX').value = c.gene1;
        document.getElementById('inspectGeneY').value = c.gene2;

        // Get filters from parameters section to carry over
        const paramLineageFilter = document.getElementById('lineageFilter').value;
        const paramHotspotGene = document.getElementById('paramHotspotGene').value;

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
            cancerFilter.innerHTML = `<option value="">All tissues (n=${plotData.length})</option>`;
            lineages.forEach(l => {
                cancerFilter.innerHTML += `<option value="${l}">${l} (n=${lineageCounts[l]})</option>`;
            });
            // Pre-select the lineage filter from parameters if it exists
            if (paramLineageFilter && lineages.includes(paramLineageFilter)) {
                cancerFilter.value = paramLineageFilter;
                this.updateScatterSubtypeFilter();
                // Also pre-select subtype from parameters
                const paramSubtype = document.getElementById('subLineageFilter')?.value;
                if (paramSubtype) {
                    document.getElementById('scatterSubtypeFilter').value = paramSubtype;
                }
            } else {
                // Reset subtype filter
                document.getElementById('scatterSubtypeFilter').innerHTML = '<option value="">All subtypes</option>';
                document.getElementById('scatterSubtypeFilter').style.display = 'none';
            }
            cancerBox.style.display = 'block';
        } else {
            cancerBox.style.display = 'none';
            document.getElementById('scatterSubtypeFilter').style.display = 'none';
        }

        // Reset color-by dropdown
        document.getElementById('colorByCategory').value = '';
        // Show/hide subtype option based on initial cancer filter
        const subtypeOpt = document.getElementById('colorBySubtypeOption');
        if (subtypeOpt) subtypeOpt.style.display = (paramLineageFilter && lineages.includes(paramLineageFilter)) ? '' : 'none';
        this._styleActiveFilters();

        // Populate hotspot genes
        // Count mutations based on cell lines with valid data (plotData)
        const hotspotSelect = document.getElementById('hotspotGene');
        const mutFilterGeneSelect = document.getElementById('mutationFilterGene');
        const cellLinesInPlot = new Set(plotData.map(d => d.cellLineId));

        if (this.mutations?.genes?.length > 0) {
            hotspotSelect.innerHTML = '<option value="">Select gene...</option>';
            mutFilterGeneSelect.innerHTML = '<option value="">No filter</option>';
            this.mutations.genes
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
            // Pre-select the hotspot gene from parameters if it exists
            if (paramHotspotGene && this.mutations.genes.includes(paramHotspotGene)) {
                hotspotSelect.value = paramHotspotGene;
            }
            document.getElementById('mutationBox').style.display = 'block';
            document.getElementById('mutationFilterBox').style.display = 'block';
        } else {
            document.getElementById('mutationBox').style.display = 'none';
            document.getElementById('mutationFilterBox').style.display = 'none';
        }

        // Populate translocation/fusion selectors (datalists, sorted by count desc)
        const transGeneInput = document.getElementById('translocationGene');
        const transFilterGeneInput = document.getElementById('translocationFilterGene');
        const transGeneDatalist = document.getElementById('translocationGeneList');
        const transFilterGeneDatalist = document.getElementById('translocationFilterGeneList');

        if (this.translocations?.genes?.length > 0) {
            transGeneInput.value = '';
            transFilterGeneInput.value = '';
            const geneCounts = [];
            for (const g of this.translocations.genes) {
                const transData = this.translocations.geneData?.[g]?.translocations || {};
                let count = 0;
                for (const cl of cellLinesInPlot) {
                    if (transData[cl] && transData[cl] > 0) count++;
                }
                if (count > 0) geneCounts.push({ gene: g, count });
            }
            geneCounts.sort((a, b) => {
                const aPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(a.gene) ? 1 : 0;
                const bPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(b.gene) ? 1 : 0;
                if (aPri !== bPri) return bPri - aPri;
                return b.count - a.count;
            });

            let transHtml = '';
            geneCounts.forEach(({ gene, count }) => {
                transHtml += `<option value="${gene}">${gene} (${count} fused)</option>`;
            });
            transGeneDatalist.innerHTML = transHtml;
            transFilterGeneDatalist.innerHTML = transHtml;
            document.getElementById('translocationBox').style.display = 'block';
            document.getElementById('translocationFilterBox').style.display = 'block';
            document.getElementById('compareAllTranslocationsBtn').style.display = '';
        } else {
            document.getElementById('translocationBox').style.display = 'none';
            document.getElementById('translocationFilterBox').style.display = 'none';
            document.getElementById('compareAllTranslocationsBtn').style.display = 'none';
        }

        // Calculate stats for ALL cells (unfiltered) for the title
        const allCellsStats = this.pearsonWithSlope(plotData.map(d => d.x), plotData.map(d => d.y));
        document.getElementById('inspectTitle').textContent =
            `${c.gene1} vs ${c.gene2} | r=${this.formatNum(allCellsStats.correlation)}, slope=${this.formatNum(allCellsStats.slope)}, n=${plotData.length} (all cells)`;

        // Show modal and render plot
        document.getElementById('inspectModal').classList.add('active');

        // Make sure controls are visible (may have been hidden by By tissue view)
        document.querySelector('.inspect-controls').style.display = '';
        document.querySelector('.inspect-layout').style.display = '';
        document.getElementById('byTissueContainer').style.display = 'none';
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

    resetAllInspectFilters() {
        // Reset cancer type and subtype filters
        document.getElementById('scatterCancerFilter').value = '';
        document.getElementById('scatterSubtypeFilter').value = '';
        document.getElementById('scatterSubtypeFilter').style.display = 'none';

        // Reset mutation filters
        document.getElementById('mutationFilterGene').value = '';
        document.getElementById('mutationFilterLevel').value = 'all';

        // Reset hotspot overlay
        document.getElementById('hotspotGene').value = '';
        document.getElementById('hotspotMode').value = 'color';
        document.getElementById('mutationCautionScatter').style.display = 'none';

        // Reset translocation filters
        document.getElementById('translocationFilterGene').value = '';
        document.getElementById('translocationFilterLevel').value = 'all';
        document.getElementById('translocationGene').value = '';
        document.getElementById('translocationMode').value = 'color';

        // Reset color-by
        document.getElementById('colorByCategory').value = '';

        this._styleActiveFilters();
        // Update the plot
        this.updateInspectPlot();
    }

    updateInspectPlot() {
        if (!this.currentInspect) return;

        const data = this.currentInspect.data;
        const gene1 = this.currentInspect.gene1;
        const gene2 = this.currentInspect.gene2;

        // Get filter settings
        const cancerFilter = document.getElementById('scatterCancerFilter').value;
        const subtypeFilter = document.getElementById('scatterSubtypeFilter').value;
        const mutFilterGene = document.getElementById('mutationFilterGene').value;
        const mutFilterLevel = document.getElementById('mutationFilterLevel').value;
        const searchTerms = document.getElementById('scatterCellSearch').value
            .split('\n').map(s => s.trim().toUpperCase()).filter(s => s);
        const fontSize = parseInt(document.getElementById('scatterFontSize')?.value) || 3;
        const hotspotGene = document.getElementById('hotspotGene').value;
        const hotspotMode = document.getElementById('hotspotMode').value;
        const transOverlayGene = document.getElementById('translocationGene').value;
        const transOverlayMode = document.getElementById('translocationMode').value;
        const transFilterGene = document.getElementById('translocationFilterGene').value;
        const transFilterLevel = document.getElementById('translocationFilterLevel').value;
        const colorByCategory = document.getElementById('colorByCategory').value;

        // Show/hide mutation caution message
        const cautionEl = document.getElementById('mutationCautionScatter');
        if (cautionEl) {
            cautionEl.style.display = (hotspotGene && hotspotMode !== 'none') ? 'block' : 'none';
        }

        // Filter by cancer type
        let filteredData = cancerFilter ?
            data.filter(d => d.lineage === cancerFilter) : data;

        // Filter by subtype
        if (subtypeFilter && this.cellLineMetadata?.primaryDisease) {
            filteredData = filteredData.filter(d =>
                this.cellLineMetadata.primaryDisease[d.cellLineId] === subtypeFilter
            );
        }

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

        // Apply translocation filter
        if (transFilterGene && this.translocations?.geneData?.[transFilterGene] && transFilterLevel !== 'all') {
            const filterTrans = this.translocations.geneData[transFilterGene].translocations;
            filteredData = filteredData.filter(d => {
                const tLevel = filterTrans[d.cellLineId] || 0;
                if (transFilterLevel === '0') return tLevel === 0;
                if (transFilterLevel === '1') return tLevel === 1;
                if (transFilterLevel === '2') return tLevel >= 2;
                if (transFilterLevel === '1+2') return tLevel >= 1;
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

        // Get translocation info for overlay
        let translocationMap = new Map();
        let translocationPartnersMap = new Map();
        if (transOverlayGene && this.translocations?.geneData?.[transOverlayGene]) {
            const tData = this.translocations.geneData[transOverlayGene];
            Object.entries(tData.translocations).forEach(([cellLine, level]) => {
                translocationMap.set(cellLine, level);
            });
            if (tData.partners) {
                Object.entries(tData.partners).forEach(([cellLine, partners]) => {
                    translocationPartnersMap.set(cellLine, partners);
                });
            }
        }

        // Add mutation level to filtered data (for overlay coloring)
        filteredData = filteredData.map(d => ({
            ...d,
            mutationLevel: mutationMap.get(d.cellLineId) || 0,
            translocationLevel: translocationMap.get(d.cellLineId) || 0,
            translocationPartners: translocationPartnersMap.get(d.cellLineId) || []
        }));

        // Build filter description for title
        let filterParts = [];
        if (cancerFilter) {
            let cancerText = cancerFilter;
            if (subtypeFilter) cancerText += ` / ${subtypeFilter}`;
            filterParts.push(`Cancer: ${cancerText}`);
        }
        if (mutFilterGene && mutFilterLevel !== 'all') {
            const levelText = mutFilterLevel === '0' ? 'WT' :
                              mutFilterLevel === '1' ? '1 mut' :
                              mutFilterLevel === '2' ? '2 mut' : '1+2 mut';
            filterParts.push(`${mutFilterGene}: ${levelText}`);
        }
        if (transFilterGene && transFilterLevel !== 'all') {
            const levelText = transFilterLevel === '0' ? 'No fusion' :
                              transFilterLevel === '1' ? '1 partner' :
                              transFilterLevel === '2' ? '2+ partners' : '1+2 fused';
            filterParts.push(`${transFilterGene}: ${levelText}`);
        }
        const filterDesc = filterParts.length > 0 ? filterParts.join(' | ') : '';

        // Show/hide plot and table based on mode
        const scatterPlot = document.getElementById('scatterPlot');
        const compareTable = document.getElementById('compareTable');

        if (hotspotMode === 'compare_table' && hotspotGene) {
            scatterPlot.style.display = 'none';
            compareTable.style.display = 'block';
            const legend = document.getElementById('colorByLegend');
            if (legend) legend.style.display = 'none';
            this.renderCompareTable(filteredData, gene1, gene2, hotspotGene, filterDesc);
            return;
        } else if (transOverlayMode === 'compare_table' && transOverlayGene) {
            scatterPlot.style.display = 'none';
            compareTable.style.display = 'block';
            const legend = document.getElementById('colorByLegend');
            if (legend) legend.style.display = 'none';
            this.renderCompareTable(filteredData, gene1, gene2, transOverlayGene, filterDesc, true);
            return;
        } else {
            scatterPlot.style.display = 'block';
            compareTable.style.display = 'none';
        }

        // Handle 3-panel mode
        if (hotspotMode === 'three_panel' && hotspotGene) {
            this.renderThreePanelPlot(filteredData, gene1, gene2, hotspotGene, searchTerms, fontSize, filterDesc, false, colorByCategory);
            return;
        }
        if (transOverlayMode === 'three_panel' && transOverlayGene) {
            this.renderThreePanelPlot(filteredData, gene1, gene2, transOverlayGene, searchTerms, fontSize, filterDesc, true, colorByCategory);
            return;
        }

        // Single panel color mode
        this.renderSinglePanelPlot(filteredData, gene1, gene2, hotspotGene, hotspotMode, searchTerms, fontSize, filterDesc, transOverlayGene, transOverlayMode, colorByCategory);
    }

    renderSinglePanelPlot(filteredData, gene1, gene2, hotspotGene, hotspotMode, searchTerms, fontSize, filterDesc = '', transOverlayGene = '', transOverlayMode = 'none', colorByCategory = '') {
        // Calculate stats for each mutation group
        const wt = filteredData.filter(d => d.mutationLevel === 0);
        const mut1 = filteredData.filter(d => d.mutationLevel === 1);
        const mut2 = filteredData.filter(d => d.mutationLevel >= 2);

        const wtStats = this.pearsonWithSlope(wt.map(d => d.x), wt.map(d => d.y));
        const mut1Stats = this.pearsonWithSlope(mut1.map(d => d.x), mut1.map(d => d.y));
        const mut2Stats = this.pearsonWithSlope(mut2.map(d => d.x), mut2.map(d => d.y));
        const allStats = this.pearsonWithSlope(filteredData.map(d => d.x), filteredData.map(d => d.y));

        // Color-by tracking
        let colorByCategories = null;
        let colorByColors = null;
        const colorByTraceIndices = {};

        // Build traces
        const traces = [];
        const xRange = [this.getInputNum('scatterXmin'),
                       this.getInputNum('scatterXmax')];
        const yRange = [this.getInputNum('scatterYmin'),
                       this.getInputNum('scatterYmax')];

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
        } else if (transOverlayMode === 'color' && transOverlayGene) {
            // Color by translocation/fusion level (purple tones)
            const tWT = filteredData.filter(d => d.translocationLevel === 0);
            const t1 = filteredData.filter(d => d.translocationLevel === 1);
            const t2 = filteredData.filter(d => d.translocationLevel >= 2);
            const tWTPct = (tWT.length / filteredData.length * 100).toFixed(1);
            const t1Pct = (t1.length / filteredData.length * 100).toFixed(1);
            const t2Pct = (t2.length / filteredData.length * 100).toFixed(1);

            const makeTransHover = (d, label) => {
                let text = `${d.cellLineName}<br>${d.lineage}<br>${label}`;
                if (d.translocationPartners && d.translocationPartners.length > 0) {
                    text += `<br>Partners: ${d.translocationPartners.join(', ')}`;
                }
                return text;
            };

            traces.push({
                x: tWT.map(d => d.x), y: tWT.map(d => d.y),
                mode: 'markers', type: 'scatter',
                text: tWT.map(d => makeTransHover(d, 'No fusion')),
                hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                marker: { color: '#9ca3af', size: 8, opacity: 0.6 },
                name: `No fusion (n=${tWT.length}, ${tWTPct}%)`
            });
            traces.push({
                x: t1.map(d => d.x), y: t1.map(d => d.y),
                mode: 'markers', type: 'scatter',
                text: t1.map(d => makeTransHover(d, '1 fusion partner')),
                hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                marker: { color: '#3b82f6', size: 10, opacity: 0.7 },
                name: `1 partner (n=${t1.length}, ${t1Pct}%)`
            });
            traces.push({
                x: t2.map(d => d.x), y: t2.map(d => d.y),
                mode: 'markers', type: 'scatter',
                text: t2.map(d => makeTransHover(d, '2+ fusion partners')),
                hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                marker: { color: '#dc2626', size: 11, opacity: 0.8 },
                name: `2+ partners (n=${t2.length}, ${t2Pct}%)`
            });
        } else if (colorByCategory === 'tissue' || colorByCategory === 'subtype') {
            // Color by tissue or subtype
            const categoryMap = {};
            filteredData.forEach(d => {
                let cat;
                if (colorByCategory === 'subtype') {
                    cat = this.cellLineMetadata?.primaryDisease?.[d.cellLineId] || d.lineage || 'Unknown';
                } else {
                    cat = d.lineage || 'Unknown';
                }
                if (!categoryMap[cat]) categoryMap[cat] = [];
                categoryMap[cat].push(d);
            });
            // Sort categories by count descending
            colorByCategories = Object.keys(categoryMap).sort((a, b) => categoryMap[b].length - categoryMap[a].length);
            colorByColors = CorrelationExplorer.CATEGORY_COLORS;
            colorByCategories.forEach((cat, i) => {
                const catData = categoryMap[cat];
                const color = i < colorByColors.length ? colorByColors[i] : '#999';
                colorByTraceIndices[cat] = [traces.length];
                traces.push({
                    x: catData.map(d => d.x),
                    y: catData.map(d => d.y),
                    mode: 'markers',
                    type: 'scatter',
                    text: catData.map(d => `${d.cellLineName}<br>${cat}`),
                    hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                    marker: { color: color, size: 8, opacity: 0.8 },
                    name: `${cat} (${catData.length})`,
                    legendgroup: cat,
                    showlegend: true
                });
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
                mode: 'markers',
                type: 'scatter',
                hovertemplate: highlightData.map(d => `${d.cellLineName} (${d.lineage || 'Unknown'})<extra></extra>`),
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

        // Build annotations for highlighted cell labels (draggable)
        const highlightAnnotations = [];
        highlightData.forEach(d => {
            const saved = this._userLabelPositions.get(d.cellLineName);
            highlightAnnotations.push({
                x: d.x, y: d.y,
                xref: 'x', yref: 'y',
                text: `${d.cellLineName} (${d.lineage || 'Unknown'})`,
                showarrow: true,
                arrowhead: 0,
                arrowcolor: '#999',
                ax: saved ? saved.ax : 0,
                ay: saved ? saved.ay : -25,
                font: { size: fontSize * 3, color: '#000' },
                bgcolor: 'rgba(255,255,255,0.7)',
                borderpad: 2
            });
        });

        // Add regression line
        if (!isNaN(allStats.slope) && document.getElementById('showCorrelationLine')?.checked !== false) {
            const meanX = filteredData.reduce((a, d) => a + d.x, 0) / filteredData.length;
            const meanY = filteredData.reduce((a, d) => a + d.y, 0) / filteredData.length;
            const intercept = meanY - allStats.slope * meanX;

            traces.push({
                x: xRange,
                y: [allStats.slope * xRange[0] + intercept, allStats.slope * xRange[1] + intercept],
                mode: 'lines',
                type: 'scatter',
                line: { color: '#5a9f4a', width: 3 },
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

        // Build title - condensed to avoid overlapping with data
        let titleLines = [`<b>${gene1} vs ${gene2}</b>`];
        if (filterDesc) {
            titleLines.push(`<span style="font-size:10px;color:#666;">${filterDesc}</span>`);
        }
        titleLines.push(`<span style="font-size:10px;">n=${filteredData.length}, r=${allStats.correlation.toFixed(3)}, slope=${allStats.slope.toFixed(3)} | mean(${meanX.toFixed(2)}, ${meanY.toFixed(2)}) med(${medianX.toFixed(2)}, ${medianY.toFixed(2)})</span>`);

        if (hotspotMode === 'color' && hotspotGene) {
            titleLines.push(`<span style="font-size:10px;"><b>${hotspotGene}:</b> WT n=${wt.length} r=${wtStats.correlation.toFixed(3)} | 1mut n=${mut1.length} r=${mut1Stats.correlation.toFixed(3)} | 2mut n=${mut2.length} r=${mut2Stats.correlation.toFixed(3)}</span>`);
        } else if (transOverlayMode === 'color' && transOverlayGene) {
            const tWT = filteredData.filter(d => d.translocationLevel === 0);
            const tFused = filteredData.filter(d => d.translocationLevel >= 1);
            const tWTStats = this.pearsonWithSlope(tWT.map(d => d.x), tWT.map(d => d.y));
            const tFusedStats = this.pearsonWithSlope(tFused.map(d => d.x), tFused.map(d => d.y));
            titleLines.push(`<span style="font-size:10px;">${transOverlayGene}: No fusion n=${tWT.length} r=${tWTStats.correlation.toFixed(3)} | Fused n=${tFused.length} r=${tFusedStats.correlation.toFixed(3)}</span>`);
        }

        const titleText = titleLines.join('<br>');

        // Title as draggable annotation
        const titleAnnotation = {
            x: this._userTitlePosition ? this._userTitlePosition.x : 0.5,
            y: this._userTitlePosition ? this._userTitlePosition.y : 1.0,
            xref: 'paper', yref: 'paper',
            xanchor: this._userTitlePosition ? 'auto' : 'center',
            yanchor: this._userTitlePosition ? 'auto' : 'bottom',
            text: titleText,
            showarrow: false,
            font: { size: 14 }
        };
        const annotations = [titleAnnotation, ...highlightAnnotations];

        // Calculate margin based on title lines
        const topMargin = 80 + (titleLines.length * 18);

        const showZero = document.getElementById('showZeroLines')?.checked !== false;

        const layout = {
            xaxis: {
                title: `${gene1} Gene Effect`,
                range: xRange,
                zeroline: showZero,
                zerolinecolor: showZero ? '#000' : '#ddd',
                zerolinewidth: showZero ? 2 : 0
            },
            yaxis: {
                title: `${gene2} Gene Effect`,
                range: yRange,
                zeroline: showZero,
                zerolinecolor: showZero ? '#000' : '#ddd',
                zerolinewidth: showZero ? 2 : 0
            },
            hovermode: 'closest',
            margin: { t: topMargin, r: ((hotspotMode === 'color' && hotspotGene) || (transOverlayMode === 'color' && transOverlayGene)) ? 240 : 30, b: colorByCategory ? 100 : 60, l: 60, autoexpand: false },
            showlegend: (hotspotMode === 'color' && hotspotGene) || (transOverlayMode === 'color' && transOverlayGene) || !!colorByCategory,
            legend: colorByCategory ? {
                orientation: 'h',
                x: 0.5,
                y: -0.15,
                xanchor: 'center',
                yanchor: 'top',
                bgcolor: 'rgba(255,255,255,0.85)',
                bordercolor: '#ddd',
                borderwidth: 1,
                font: { size: 10 },
                tracegroupgap: 0,
                entrywidth: 120,
                entrywidthmode: 'pixels'
            } : {
                x: 1.02,
                y: 1,
                xanchor: 'left',
                yanchor: 'top',
                bgcolor: 'rgba(255,255,255,0.85)',
                bordercolor: '#ddd',
                borderwidth: 1,
                title: { text: (transOverlayMode === 'color' && transOverlayGene) ? `${transOverlayGene} (fusion)` : hotspotGene, font: { size: 11 } },
                font: { size: 11 }
            },
            annotations: annotations,
            plot_bgcolor: '#fafafa'
        };

        // Apply user-dragged legend position if available
        if (this._userLegendPosition) {
            layout.legend.x = this._userLegendPosition.x;
            layout.legend.y = this._userLegendPosition.y;
            layout.legend.xanchor = 'auto';
            layout.legend.yanchor = 'auto';
        }

        // Apply plot dimensions (square by default)
        // Values refer to the plot area; margins are added on top
        const plotContainer = document.getElementById('scatterPlot');
        const widthEl = document.getElementById('plotWidth');
        const heightEl = document.getElementById('plotHeight');
        const m = layout.margin;
        let plotAreaW = parseInt(widthEl?.value);
        let plotAreaH = parseInt(heightEl?.value);
        if (isNaN(plotAreaW) || isNaN(plotAreaH)) {
            // Default: square, fit within viewport
            const defaultSize = 400;
            if (isNaN(plotAreaW)) plotAreaW = defaultSize;
            if (isNaN(plotAreaH)) plotAreaH = defaultSize;
        }
        // Show current values in the inputs
        if (widthEl) widthEl.value = plotAreaW;
        if (heightEl) heightEl.value = plotAreaH;
        layout.width = plotAreaW + m.l + m.r;
        layout.height = plotAreaH + m.t + m.b;
        plotContainer.style.width = layout.width + 'px';
        plotContainer.style.height = layout.height + 'px';

        Plotly.newPlot('scatterPlot', traces, layout, {
            responsive: false,
            edits: { annotationPosition: true, annotationTail: true, legendPosition: true }
        }).then(plotEl => {
            // Listen for legend and title drag events
            let isProgrammaticRelayout = false;
            plotEl.on('plotly_relayout', (relayoutData) => {
                if (isProgrammaticRelayout) return;
                if (relayoutData['legend.x'] !== undefined && relayoutData['legend.y'] !== undefined) {
                    this._userLegendPosition = {
                        x: relayoutData['legend.x'],
                        y: relayoutData['legend.y']
                    };
                }
                if (relayoutData['annotations[0].x'] !== undefined && relayoutData['annotations[0].y'] !== undefined) {
                    this._userTitlePosition = {
                        x: relayoutData['annotations[0].x'],
                        y: relayoutData['annotations[0].y']
                    };
                }
                // Capture dragged cell label positions (annotations[1+])
                for (let i = 0; i < highlightData.length; i++) {
                    const idx = i + 1; // offset by 1 for title annotation
                    const axKey = `annotations[${idx}].ax`;
                    const ayKey = `annotations[${idx}].ay`;
                    if (relayoutData[axKey] !== undefined || relayoutData[ayKey] !== undefined) {
                        const current = this._userLabelPositions.get(highlightData[i].cellLineName) || { ax: 0, ay: -25 };
                        this._userLabelPositions.set(highlightData[i].cellLineName, {
                            ax: relayoutData[axKey] !== undefined ? relayoutData[axKey] : current.ax,
                            ay: relayoutData[ayKey] !== undefined ? relayoutData[ayKey] : current.ay
                        });
                    }
                }
            });

            // Auto-reposition vertical legend outside the actual plot domain
            if (!this._userLegendPosition && layout.showlegend && layout.legend?.orientation !== 'h') {
                const fl = plotEl._fullLayout;
                if (fl && fl.xaxis && fl.yaxis) {
                    const xDomain = fl.xaxis.domain;
                    const yDomain = fl.yaxis.domain;
                    // Place legend just outside the plot to the right
                    isProgrammaticRelayout = true;
                    Plotly.relayout('scatterPlot', {
                        'legend.x': xDomain[1] + 0.02,
                        'legend.y': yDomain[1],
                        'legend.xanchor': 'left',
                        'legend.yanchor': 'top'
                    }).then(() => { isProgrammaticRelayout = false; });
                }
            }

            // Add click handler
            this.setupScatterClickHandler(filteredData);

            // Build custom HTML legend for color-by mode
            if (colorByCategories) {
                this._buildColorByLegend(colorByCategories, colorByColors, colorByTraceIndices);
            } else {
                const legend = document.getElementById('colorByLegend');
                if (legend) legend.style.display = 'none';
            }
        });
    }

    renderThreePanelPlot(filteredData, gene1, gene2, hotspotGene, searchTerms, fontSize, filterDesc = '', isFusion = false, colorByCategory = '') {
        const levelField = isFusion ? 'translocationLevel' : 'mutationLevel';
        const wt = filteredData.filter(d => d[levelField] === 0);
        const mut1 = filteredData.filter(d => d[levelField] === 1);
        const mut2 = filteredData.filter(d => d[levelField] >= 2);

        const xRange = [this.getInputNum('scatterXmin'),
                       this.getInputNum('scatterXmax')];
        const yRange = [this.getInputNum('scatterYmin'),
                       this.getInputNum('scatterYmax')];

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

        // Panel labels
        const panelLabels = isFusion
            ? ['No fusion', '1 partner', '2+ partners']
            : ['WT', '1 mut', '2 mut'];

        // Build category map for color-by mode (shared across panels)
        let categoryOrder = null;
        if (colorByCategory === 'tissue' || colorByCategory === 'subtype') {
            const catCounts = {};
            filteredData.forEach(d => {
                const cat = colorByCategory === 'subtype'
                    ? (this.cellLineMetadata?.primaryDisease?.[d.cellLineId] || d.lineage || 'Unknown')
                    : (d.lineage || 'Unknown');
                catCounts[cat] = (catCounts[cat] || 0) + 1;
            });
            categoryOrder = Object.keys(catCounts).sort((a, b) => catCounts[b] - catCounts[a]);
        }

        const getCategory = (d) => {
            if (colorByCategory === 'subtype') return this.cellLineMetadata?.primaryDisease?.[d.cellLineId] || d.lineage || 'Unknown';
            return d.lineage || 'Unknown';
        };

        // Build panel traces with optional category coloring
        const panelConfigs = [
            { data: wt, xaxis: 'x', yaxis: 'y', color: '#9ca3af', size: 7, opacity: 0.6 },
            { data: mut1, xaxis: 'x2', yaxis: 'y2', color: '#3b82f6', size: 8, opacity: 0.7 },
            { data: mut2, xaxis: 'x3', yaxis: 'y3', color: '#dc2626', size: 8, opacity: 0.7 }
        ];

        const colors = CorrelationExplorer.CATEGORY_COLORS;
        const colorByTraceIndices = {};

        panelConfigs.forEach((cfg, i) => {
            if (categoryOrder) {
                // Color by category within each panel
                categoryOrder.forEach((cat, ci) => {
                    const catData = cfg.data.filter(d => getCategory(d) === cat);
                    if (catData.length === 0) return;
                    if (!colorByTraceIndices[cat]) colorByTraceIndices[cat] = [];
                    colorByTraceIndices[cat].push(traces.length);
                    traces.push({
                        x: catData.map(d => d.x), y: catData.map(d => d.y),
                        xaxis: cfg.xaxis, yaxis: cfg.yaxis,
                        mode: 'markers', type: 'scatter',
                        text: catData.map(d => `${d.cellLineName}<br>${cat}`),
                        hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                        marker: { color: ci < colors.length ? colors[ci] : '#999', size: cfg.size, opacity: 0.8 },
                        name: `${cat} (${catData.length})`,
                        showlegend: i === 0,
                        legendgroup: cat
                    });
                });
            } else {
                traces.push({
                    x: cfg.data.map(d => d.x), y: cfg.data.map(d => d.y),
                    xaxis: cfg.xaxis, yaxis: cfg.yaxis,
                    mode: 'markers', type: 'scatter',
                    text: cfg.data.map(d => `${d.cellLineName}<br>${d.lineage}`),
                    hovertemplate: '%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<extra></extra>',
                    marker: { color: cfg.color, size: cfg.size, opacity: cfg.opacity },
                    name: panelLabels[i], showlegend: false
                });
            }
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

        if (document.getElementById('showCorrelationLine')?.checked !== false) {
            addRegressionLine(wt, wtStats, 'x', 'y', '#5a9f4a');
            addRegressionLine(mut1, mut1Stats, 'x2', 'y2', '#5a9f4a');
            addRegressionLine(mut2, mut2Stats, 'x3', 'y3', '#5a9f4a');
        }

        // Add highlighted cells for each panel
        const threePanelHighlightAnnotations = [];
        const addHighlights = (data, xaxis, yaxis) => {
            const xrefAxis = xaxis === 'x' ? 'x' : xaxis;
            const yrefAxis = yaxis === 'y' ? 'y' : yaxis;
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
                    mode: 'markers',
                    type: 'scatter',
                    hovertemplate: highlightData.map(d => `${d.cellLineName} (${d.lineage || 'Unknown'})<extra></extra>`),
                    marker: {
                        color: '#f59e0b',
                        size: 10,
                        symbol: 'circle',
                        line: { color: '#000', width: 2 }
                    },
                    showlegend: false
                });
                highlightData.forEach(d => {
                    const saved = this._userLabelPositions.get(d.cellLineName);
                    threePanelHighlightAnnotations.push({
                        x: d.x, y: d.y,
                        xref: xrefAxis, yref: yrefAxis,
                        text: `${d.cellLineName} (${d.lineage || 'Unknown'})`,
                        showarrow: true,
                        arrowhead: 0,
                        arrowcolor: '#999',
                        ax: saved ? saved.ax : 0,
                        ay: saved ? saved.ay : -25,
                        font: { size: fontSize * 3, color: '#000' },
                        bgcolor: 'rgba(255,255,255,0.7)',
                        borderpad: 2
                    });
                });
            }
        };

        addHighlights(wt, 'x', 'y');
        addHighlights(mut1, 'x2', 'y2');
        addHighlights(mut2, 'x3', 'y3');

        // Build title with filter info
        const stratLabel = isFusion ? 'fusion stratification' : 'hotspot mutation stratification';
        let titleText = `<b>${gene1} vs ${gene2} - ${hotspotGene} ${stratLabel}</b>`;
        if (filterDesc) {
            titleText += `<br><span style="font-size: 11px; color: #666;">Filter: ${filterDesc}</span>`;
        }

        // Annotation labels for panels
        const annotLabels = isFusion
            ? [`<b>No fusion</b>`, `<b>1 partner</b>`, `<b>2+ partners</b>`]
            : [`<b>WT (0 mut)</b>`, `<b>1 mutation</b>`, `<b>2 mutations</b>`];

        // Title annotation (draggable)
        const titleAnnotation = {
            x: this._userTitlePosition ? this._userTitlePosition.x : 0.5,
            y: this._userTitlePosition ? this._userTitlePosition.y : 1.12,
            xref: 'paper', yref: 'paper',
            xanchor: this._userTitlePosition ? 'auto' : 'center',
            yanchor: this._userTitlePosition ? 'auto' : 'bottom',
            text: titleText,
            showarrow: false,
            font: { size: 14 }
        };

        const layout = {
            grid: { rows: 1, columns: 3, pattern: 'independent' },
            xaxis: { title: `${gene1} Effect`, range: xRange, domain: [0, 0.28] },
            yaxis: { title: `${gene2} Effect`, range: yRange },
            xaxis2: { title: `${gene1} Effect`, range: xRange, domain: [0.36, 0.64] },
            yaxis2: { range: yRange, anchor: 'x2' },
            xaxis3: { title: `${gene1} Effect`, range: xRange, domain: [0.72, 1] },
            yaxis3: { range: yRange, anchor: 'x3' },
            annotations: [
                titleAnnotation,
                { x: 0.14, y: 1.02, xref: 'paper', yref: 'paper',
                  text: `${annotLabels[0]} n=${wt.length}<br>r=${wtStats.correlation.toFixed(3)}, slope=${wtStats.slope.toFixed(3)}<br>mean: x=${wtExtra.meanX.toFixed(2)}, y=${wtExtra.meanY.toFixed(2)}<br>median: x=${wtExtra.medianX.toFixed(2)}, y=${wtExtra.medianY.toFixed(2)}`,
                  showarrow: false, font: { size: 9 } },
                { x: 0.5, y: 1.02, xref: 'paper', yref: 'paper',
                  text: `${annotLabels[1]} n=${mut1.length}<br>r=${mut1Stats.correlation.toFixed(3)}, slope=${mut1Stats.slope.toFixed(3)}<br>mean: x=${mut1Extra.meanX.toFixed(2)}, y=${mut1Extra.meanY.toFixed(2)}<br>median: x=${mut1Extra.medianX.toFixed(2)}, y=${mut1Extra.medianY.toFixed(2)}`,
                  showarrow: false, font: { size: 9 } },
                { x: 0.86, y: 1.02, xref: 'paper', yref: 'paper',
                  text: `${annotLabels[2]} n=${mut2.length}<br>r=${mut2Stats.correlation.toFixed(3)}, slope=${mut2Stats.slope.toFixed(3)}<br>mean: x=${mut2Extra.meanX.toFixed(2)}, y=${mut2Extra.meanY.toFixed(2)}<br>median: x=${mut2Extra.medianX.toFixed(2)}, y=${mut2Extra.medianY.toFixed(2)}`,
                  showarrow: false, font: { size: 9 } },
                ...threePanelHighlightAnnotations
            ],
            margin: { t: filterDesc ? 160 : 140, r: 30, b: categoryOrder ? 100 : 60, l: 60 },
            showlegend: !!categoryOrder,
            legend: {
                orientation: 'h',
                x: 0.5,
                y: -0.15,
                xanchor: 'center',
                yanchor: 'top',
                bgcolor: 'rgba(255,255,255,0.85)',
                bordercolor: '#ddd',
                borderwidth: 1,
                font: { size: 10 },
                tracegroupgap: 0,
                entrywidth: 120,
                entrywidthmode: 'pixels'
            },
            plot_bgcolor: '#fafafa'
        };

        // Apply user-dragged legend position if available
        if (this._userLegendPosition && layout.showlegend) {
            layout.legend.x = this._userLegendPosition.x;
            layout.legend.y = this._userLegendPosition.y;
            layout.legend.xanchor = 'auto';
            layout.legend.yanchor = 'auto';
        }

        // Collect all highlight data across panels for relayout tracking
        const allThreePanelHighlightData = [wt, mut1, mut2].flatMap(data =>
            data.filter(d =>
                searchTerms.some(term =>
                    d.cellLineName.toUpperCase().includes(term) ||
                    d.cellLineId.toUpperCase().includes(term)
                ) || this.clickedCells.has(d.cellLineName)
            )
        );

        // Apply plot dimensions (wide default for three-panel)
        const plotContainer3 = document.getElementById('scatterPlot');
        const m3 = layout.margin;
        let plotAreaW3 = parseInt(document.getElementById('plotWidth')?.value);
        let plotAreaH3 = parseInt(document.getElementById('plotHeight')?.value);
        if (isNaN(plotAreaW3)) plotAreaW3 = Math.max(600, Math.min(1000, window.innerWidth - 500));
        if (isNaN(plotAreaH3)) plotAreaH3 = Math.max(300, Math.min(500, window.innerHeight - 300));
        layout.width = plotAreaW3 + m3.l + m3.r;
        layout.height = plotAreaH3 + m3.t + m3.b;
        plotContainer3.style.width = layout.width + 'px';
        plotContainer3.style.height = layout.height + 'px';

        Plotly.newPlot('scatterPlot', traces, layout, {
            responsive: false,
            edits: { annotationPosition: true, annotationTail: true, legendPosition: true }
        }).then(plotEl => {
            plotEl.on('plotly_relayout', (relayoutData) => {
                if (relayoutData['legend.x'] !== undefined && relayoutData['legend.y'] !== undefined) {
                    this._userLegendPosition = {
                        x: relayoutData['legend.x'],
                        y: relayoutData['legend.y']
                    };
                }
                if (relayoutData['annotations[0].x'] !== undefined && relayoutData['annotations[0].y'] !== undefined) {
                    this._userTitlePosition = {
                        x: relayoutData['annotations[0].x'],
                        y: relayoutData['annotations[0].y']
                    };
                }
                // Capture dragged cell label positions (annotations[4+] in three-panel: title + 3 panel stats)
                const labelOffset = 4;
                for (let i = 0; i < allThreePanelHighlightData.length; i++) {
                    const idx = i + labelOffset;
                    const axKey = `annotations[${idx}].ax`;
                    const ayKey = `annotations[${idx}].ay`;
                    if (relayoutData[axKey] !== undefined || relayoutData[ayKey] !== undefined) {
                        const current = this._userLabelPositions.get(allThreePanelHighlightData[i].cellLineName) || { ax: 0, ay: -25 };
                        this._userLabelPositions.set(allThreePanelHighlightData[i].cellLineName, {
                            ax: relayoutData[axKey] !== undefined ? relayoutData[axKey] : current.ax,
                            ay: relayoutData[ayKey] !== undefined ? relayoutData[ayKey] : current.ay
                        });
                    }
                }
            });
            this.setupScatterClickHandler(filteredData);

            // Build custom HTML legend for color-by mode
            if (categoryOrder) {
                this._buildColorByLegend(categoryOrder, colors, colorByTraceIndices);
            } else {
                const legend = document.getElementById('colorByLegend');
                if (legend) legend.style.display = 'none';
            }
        });
    }

    _buildColorByLegend(categories, colors, categoryTraceIndices) {
        const container = document.getElementById('colorByLegend');
        if (!container) return;

        if (!categories || categories.length === 0) {
            container.style.display = 'none';
            return;
        }

        container.style.display = 'block';
        let html = '<div style="display: flex; gap: 2px; margin-bottom: 4px;">';
        html += '<button id="colorByShowAll" class="btn btn-secondary btn-sm" style="font-size: 9px; padding: 1px 5px;">Show all</button>';
        html += '<button id="colorByHideAll" class="btn btn-secondary btn-sm" style="font-size: 9px; padding: 1px 5px;">Hide all</button>';
        html += '</div>';
        html += '<div style="display: flex; flex-wrap: wrap; gap: 2px;">';
        categories.forEach((cat, i) => {
            const color = i < colors.length ? colors[i] : '#999';
            html += `<span class="color-by-chip" data-cat-index="${i}" style="display: inline-flex; align-items: center; gap: 3px; padding: 1px 5px; border-radius: 3px; font-size: 10px; cursor: pointer; border: 1px solid ${color}; background: rgba(255,255,255,0.9); white-space: nowrap; user-select: none;">`;
            html += `<span style="width: 8px; height: 8px; border-radius: 50%; background: ${color}; display: inline-block; flex-shrink: 0;"></span>`;
            html += `${cat}`;
            html += `</span>`;
        });
        html += '</div>';
        container.innerHTML = html;

        // Track hidden categories
        const hidden = new Set();

        // Click chip to toggle
        container.querySelectorAll('.color-by-chip').forEach(chip => {
            chip.addEventListener('click', () => {
                const idx = parseInt(chip.dataset.catIndex);
                const cat = categories[idx];
                const indices = categoryTraceIndices[cat];
                if (!indices) return;

                if (hidden.has(cat)) {
                    hidden.delete(cat);
                    chip.style.opacity = '1';
                    chip.style.textDecoration = 'none';
                    Plotly.restyle('scatterPlot', { visible: true }, indices);
                } else {
                    hidden.add(cat);
                    chip.style.opacity = '0.3';
                    chip.style.textDecoration = 'line-through';
                    Plotly.restyle('scatterPlot', { visible: false }, indices);
                }
            });
        });

        // Show all
        document.getElementById('colorByShowAll')?.addEventListener('click', () => {
            hidden.clear();
            container.querySelectorAll('.color-by-chip').forEach(c => {
                c.style.opacity = '1';
                c.style.textDecoration = 'none';
            });
            const allIndices = Object.values(categoryTraceIndices).flat();
            if (allIndices.length > 0) Plotly.restyle('scatterPlot', { visible: true }, allIndices);
        });

        // Hide all
        document.getElementById('colorByHideAll')?.addEventListener('click', () => {
            categories.forEach(cat => hidden.add(cat));
            container.querySelectorAll('.color-by-chip').forEach(c => {
                c.style.opacity = '0.3';
                c.style.textDecoration = 'line-through';
            });
            const allIndices = Object.values(categoryTraceIndices).flat();
            if (allIndices.length > 0) Plotly.restyle('scatterPlot', { visible: false }, allIndices);
        });
    }

    renderCompareTable(filteredData, gene1, gene2, hotspotGene, filterDesc = '', isFusion = false) {
        // Group by cancer type (lineage) - comparing 0 vs 2+ level only
        const levelField = isFusion ? 'translocationLevel' : 'mutationLevel';
        const lineageGroups = {};
        filteredData.forEach(d => {
            if (!d.lineage) return;
            if (!lineageGroups[d.lineage]) {
                lineageGroups[d.lineage] = { wt: [], mut: [] };
            }
            if (d[levelField] === 0) {
                lineageGroups[d.lineage].wt.push(d);
            } else if (d[levelField] >= (isFusion ? 1 : 2)) {
                lineageGroups[d.lineage].mut.push(d);
            }
            // Note: for hotspot, mutationLevel === 1 is excluded from comparison
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
        const typeLabel = isFusion ? 'Fusion' : 'Mutation';
        const wtLabel = isFusion ? 'No fusion' : 'WT';
        const mutLabel = isFusion ? 'Fused (1+)' : 'Mut';
        const filterInfo = filterDesc ? `<p style="font-size: 11px; color: #333; margin-bottom: 8px; background: #f0f9ff; padding: 4px 8px; border-radius: 4px;"><b>Filter:</b> ${filterDesc}</p>` : '';
        const exclusionNote = isFusion ? '' : ' Note: Cells with exactly 1 mutation are excluded from this comparison.';
        const wtDesc = isFusion ? `no ${hotspotGene} fusions` : `0 ${hotspotGene} mutations`;
        const mutDesc = isFusion ? `${hotspotGene} fused (1+)` : `2+ ${hotspotGene} mutations`;
        let html = `
            <h4 style="margin-bottom: 8px;">Effect of <span style="color: #0066cc;">${hotspotGene}</span> ${typeLabel} on ${gene1} vs ${gene2} Correlation</h4>
            ${filterInfo}
            <p style="font-size: 11px; color: #666; margin-bottom: 8px;">
                Comparing correlation between ${wtLabel} (${wtDesc}) vs ${mutLabel} (${mutDesc}) cells, stratified by cancer type.${exclusionNote}
                <strong>Click a cancer type</strong> to view its scatter plot with the ${hotspotGene} ${typeLabel.toLowerCase()} overlay.
            </p>
            <p style="font-size: 10px; color: #0c4a6e; background: #f0f9ff; padding: 4px 8px; border-radius: 4px; margin-bottom: 12px;">
                <b>Statistics:</b> p(Δr) uses Fisher z-transformation to compare correlations. p(Δslope) is an approximation based on correlation difference.
            </p>
            <div style="overflow-x: auto;">
            <table id="compareByCancerTable" class="data-table" style="width: 100%; font-size: 12px;">
                <thead>
                    <tr>
                        <th data-col="0" style="cursor: pointer;">Cancer Type ▼</th>
                        <th data-col="1" style="cursor: pointer; border-left: 2px solid #2563eb;">N (${wtLabel})</th>
                        <th data-col="2" style="cursor: pointer;">r (${wtLabel})</th>
                        <th data-col="3" style="cursor: pointer;">slope (${wtLabel})</th>
                        <th data-col="4" style="cursor: pointer; border-left: 2px solid #dc2626;">N (${mutLabel})</th>
                        <th data-col="5" style="cursor: pointer;">r (${mutLabel})</th>
                        <th data-col="6" style="cursor: pointer;">slope (${mutLabel})</th>
                        <th data-col="7" style="cursor: pointer; border-left: 2px solid #6b7280;">Δr</th>
                        <th data-col="8" style="cursor: pointer;">p(Δr)</th>
                        <th data-col="9" style="cursor: pointer;">Δslope</th>
                        <th data-col="10" style="cursor: pointer;">p(Δslope)</th>
                    </tr>
                </thead>
                <tbody>
        `;

        tableData.forEach(row => {
            const deltaRColor = row.deltaR < 0 ? '#dc2626' : '#5a9f4a';
            const deltaSlopeColor = row.deltaSlope < 0 ? '#dc2626' : '#5a9f4a';

            html += `
                <tr class="clickable-row" data-lineage="${row.lineage}" style="cursor: pointer;" title="Click to view ${row.lineage} scatter plot with ${hotspotGene} overlay">
                    <td style="color: var(--green-700); font-weight: 500;">${row.lineage}</td>
                    <td style="border-left: 2px solid #2563eb;">${row.nWT}</td>
                    <td>${row.rWT.toFixed(3)}</td>
                    <td>${row.slopeWT.toFixed(3)}</td>
                    <td style="border-left: 2px solid #dc2626;">${row.nMut}</td>
                    <td>${row.rMut.toFixed(3)}</td>
                    <td>${row.slopeMut.toFixed(3)}</td>
                    <td style="border-left: 2px solid #6b7280; color: ${deltaRColor}; font-weight: 600;">${row.deltaR.toFixed(3)}</td>
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
            let csv = `# Correlation: ${gene1} vs ${gene2}\n`;
            csv += `# ${isFusion ? 'Fusion' : 'Hotspot'} filter: ${hotspotGene}\n`;
            csv += `# Comparing ${wtLabel} vs ${mutLabel} by cancer type\n`;
            csv += `Cancer Type,N (${wtLabel}),r (${wtLabel}),slope (${wtLabel}),N (${mutLabel}),r (${mutLabel}),slope (${mutLabel}),Δr,p(Δr),Δslope,p(Δslope)\n`;
            tableData.forEach(row => {
                csv += `"${row.lineage}",${row.nWT},${row.rWT.toFixed(4)},${row.slopeWT.toFixed(4)},${row.nMut},${row.rMut.toFixed(4)},${row.slopeMut.toFixed(4)},${row.deltaR.toFixed(4)},${row.pR.toExponential(2)},${row.deltaSlope.toFixed(4)},${row.pSlope.toExponential(2)}\n`;
            });
            this.downloadFile(csv, `correlation_${gene1}_vs_${gene2}_by_${hotspotGene}_${isFusion ? 'fusion' : 'mutation'}.csv`, 'text/csv');
        });

        // Make table sortable
        this.setupSortableTable('compareByCancerTable');

        // Make rows clickable - clicking a cancer type shows scatter with that cancer + overlay filter
        document.querySelectorAll('#compareByCancerTable .clickable-row').forEach(row => {
            row.addEventListener('click', () => {
                const lineage = row.dataset.lineage;
                // Set cancer type filter
                document.getElementById('scatterCancerFilter').value = lineage;
                this.updateScatterSubtypeFilter();
                this._styleActiveFilters();
                // Keep the gene as color overlay
                if (isFusion) {
                    document.getElementById('translocationGene').value = hotspotGene;
                    document.getElementById('translocationMode').value = 'color';
                } else {
                    document.getElementById('hotspotGene').value = hotspotGene;
                    document.getElementById('hotspotMode').value = 'color';
                }
                // Switch back to scatter plot
                document.getElementById('compareTable').style.display = 'none';
                document.getElementById('scatterPlot').style.display = 'block';
                this.updateInspectPlot();
            });
        });
    }

    showCompareAllMutations() {
        if (!this.currentInspect) return;

        const { gene1, gene2, data } = this.currentInspect;

        // Apply current filters
        const cancerFilter = document.getElementById('scatterCancerFilter').value;
        const subtypeFilter = document.getElementById('scatterSubtypeFilter').value;
        const mutFilterGene = document.getElementById('mutationFilterGene').value;
        const mutFilterLevel = document.getElementById('mutationFilterLevel').value;

        let filteredData = cancerFilter ?
            data.filter(d => d.lineage === cancerFilter) : data;

        // Apply subtype filter
        if (subtypeFilter && this.cellLineMetadata?.primaryDisease) {
            filteredData = filteredData.filter(d =>
                this.cellLineMetadata.primaryDisease[d.cellLineId] === subtypeFilter
            );
        }

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
            let cancerText = `Cancer: ${cancerFilter}`;
            if (subtypeFilter) {
                cancerText += ` (${subtypeFilter})`;
            }
            filterParts.push(cancerText);
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
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                <h4 style="margin: 0;">Mutations affecting ${gene1} vs ${gene2}</h4>
                <div>
                    <button class="btn btn-primary btn-sm" id="backToGraphBtn" style="margin-right: 8px;">← Back to Graph</button>
                    <button class="btn btn-success btn-sm" id="downloadMutCompareCSV">Download CSV</button>
                </div>
            </div>
            ${filterInfo}
            <p style="font-size: 11px; color: #666; margin-bottom: 6px;">
                Comparing WT (0 mutations) vs Mutant (2+ mutations). Sorted by p-value.
            </p>
            <p style="font-size: 10px; color: #0c4a6e; background: #f0f9ff; padding: 4px 8px; border-radius: 4px; margin-bottom: 8px;">
                <b>Statistics:</b> p(Δr) uses Fisher z-transformation to test if correlations differ significantly between WT and mutant cells.
            </p>
            <p style="font-size: 11px; color: #059669; margin-bottom: 8px;">Click any row to color scatter plot by that mutation</p>
            <div class="table-container" style="max-height: 380px; overflow-y: auto;">
            <table id="compareMutationsTable" class="data-table" style="width: 100%; font-size: 11px;">
                <thead>
                    <tr>
                        <th data-sort="mutGene" data-type="string" style="cursor: pointer;">Mutation ↕</th>
                        <th data-sort="nWT" data-type="number" style="cursor: pointer; border-left: 2px solid #2563eb;">N(WT) ↕</th>
                        <th data-sort="rWT" data-type="number" style="cursor: pointer;">r(WT) ↕</th>
                        <th data-sort="nMut" data-type="number" style="cursor: pointer; border-left: 2px solid #dc2626;">N(Mut) ↕</th>
                        <th data-sort="rMut" data-type="number" style="cursor: pointer;">r(Mut) ↕</th>
                        <th data-sort="deltaR" data-type="number" style="cursor: pointer; border-left: 2px solid #6b7280;">Δr ↕</th>
                        <th data-sort="pR" data-type="number" style="cursor: pointer;">p(Δr) ↕</th>
                    </tr>
                </thead>
                <tbody>
        `;

        tableData.forEach(row => {
            const deltaRColor = row.deltaR < 0 ? '#dc2626' : '#5a9f4a';
            const pHighlight = row.pR < 0.05 ? 'background: #fef3c7;' : '';

            html += `
                <tr class="clickable-mut-row" data-mut-gene="${row.mutGene}" style="${pHighlight} cursor: pointer;">
                    <td><b>${row.mutGene}</b></td>
                    <td style="text-align: center; border-left: 2px solid #2563eb;">${row.nWT}</td>
                    <td style="text-align: center;">${row.rWT.toFixed(3)}</td>
                    <td style="text-align: center; border-left: 2px solid #dc2626;">${row.nMut}</td>
                    <td style="text-align: center;">${row.rMut.toFixed(3)}</td>
                    <td style="text-align: center; border-left: 2px solid #6b7280; color: ${deltaRColor}; font-weight: 600;">${row.deltaR.toFixed(3)}</td>
                    <td style="text-align: center;">${row.pR.toExponential(1)}</td>
                </tr>
            `;
        });

        html += `
                </tbody>
            </table>
            </div>
            <p style="font-size: 11px; color: #666; margin-top: 8px;">
                Yellow = p < 0.05. This analysis may be biased as mutations select for cancer types.
            </p>
        `;

        document.getElementById('compareTable').innerHTML = html;

        // Add back to graph handler
        document.getElementById('backToGraphBtn')?.addEventListener('click', () => {
            document.getElementById('compareTable').style.display = 'none';
            document.getElementById('scatterPlot').style.display = 'block';
        });

        // Add download handler
        document.getElementById('downloadMutCompareCSV')?.addEventListener('click', () => {
            let csv = 'Mutation_Gene,N_WT,r_WT,slope_WT,N_Mut,r_Mut,slope_Mut,Delta_r,p_Delta_r,Delta_slope\n';
            tableData.forEach(row => {
                csv += `${row.mutGene},${row.nWT},${row.rWT.toFixed(4)},${row.slopeWT.toFixed(4)},${row.nMut},${row.rMut.toFixed(4)},${row.slopeMut.toFixed(4)},${row.deltaR.toFixed(4)},${row.pR.toExponential(2)},${row.deltaSlope.toFixed(4)}\n`;
            });
            this.downloadFile(csv, `${gene1}_vs_${gene2}_all_mutations_comparison.csv`, 'text/csv');
        });

        // Add row click handlers to filter by mutation
        document.querySelectorAll('#compareMutationsTable .clickable-mut-row').forEach(row => {
            row.addEventListener('click', () => {
                const mutGene = row.dataset.mutGene;
                // Set the hotspot overlay and go back to scatter plot
                document.getElementById('hotspotGene').value = mutGene;
                document.getElementById('hotspotMode').value = 'color';
                document.getElementById('compareTable').style.display = 'none';
                document.getElementById('scatterPlot').style.display = 'block';
                this.updateInspectPlot();
            });
        });

        // Make table sortable
        this.setupSortableTable('compareMutationsTable');
    }

    showCompareAllTranslocations() {
        if (!this.currentInspect) return;

        const { gene1, gene2, data } = this.currentInspect;

        // Apply current filters
        const cancerFilter = document.getElementById('scatterCancerFilter').value;
        const subtypeFilter = document.getElementById('scatterSubtypeFilter').value;
        const mutFilterGene = document.getElementById('mutationFilterGene').value;
        const mutFilterLevel = document.getElementById('mutationFilterLevel').value;
        const transFilterGene = document.getElementById('translocationFilterGene').value;
        const transFilterLevel = document.getElementById('translocationFilterLevel').value;

        let filteredData = cancerFilter ?
            data.filter(d => d.lineage === cancerFilter) : data;

        if (subtypeFilter && this.cellLineMetadata?.primaryDisease) {
            filteredData = filteredData.filter(d =>
                this.cellLineMetadata.primaryDisease[d.cellLineId] === subtypeFilter
            );
        }

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

        if (transFilterGene && this.translocations?.geneData?.[transFilterGene] && transFilterLevel !== 'all') {
            const filterTrans = this.translocations.geneData[transFilterGene].translocations;
            filteredData = filteredData.filter(d => {
                const tLevel = filterTrans[d.cellLineId] || 0;
                if (transFilterLevel === '0') return tLevel === 0;
                if (transFilterLevel === '1') return tLevel === 1;
                if (transFilterLevel === '2') return tLevel >= 2;
                if (transFilterLevel === '1+2') return tLevel >= 1;
                return true;
            });
        }

        let filterParts = [];
        if (cancerFilter) {
            let cancerText = `Cancer: ${cancerFilter}`;
            if (subtypeFilter) cancerText += ` (${subtypeFilter})`;
            filterParts.push(cancerText);
        }
        if (mutFilterGene && mutFilterLevel !== 'all') {
            const levelText = mutFilterLevel === '0' ? 'WT' :
                              mutFilterLevel === '1' ? '1 mut' :
                              mutFilterLevel === '2' ? '2 mut' : '1+2 mut';
            filterParts.push(`${mutFilterGene}: ${levelText}`);
        }
        if (transFilterGene && transFilterLevel !== 'all') {
            const levelText = transFilterLevel === '0' ? 'WT' :
                              transFilterLevel === '1' ? '1 fusion' :
                              transFilterLevel === '2' ? '2 fusions' : '1+2 fusions';
            filterParts.push(`${transFilterGene}: ${levelText}`);
        }
        const filterDesc = filterParts.length > 0 ? filterParts.join(' | ') : '';

        document.getElementById('scatterPlot').style.display = 'none';
        document.getElementById('compareTable').style.display = 'block';

        this.renderTranslocationComparisonTable(filteredData, gene1, gene2, filterDesc);
    }

    renderTranslocationComparisonTable(filteredData, gene1, gene2, filterDesc = '') {
        if (!this.translocations || !this.translocations.geneData) {
            document.getElementById('compareTable').innerHTML = '<p>No fusion data available.</p>';
            return;
        }

        const tableData = [];
        const transGenes = Object.keys(this.translocations.geneData);

        // Pre-filter: only process genes that have ≥3 fused cells in filteredData
        const cellIdSet = new Set(filteredData.map(d => d.cellLineId));
        const candidateGenes = transGenes.filter(tGene => {
            const transData = this.translocations.geneData[tGene].translocations;
            let fusedCount = 0;
            for (const cl of cellIdSet) {
                if (transData[cl] && transData[cl] > 0) fusedCount++;
                if (fusedCount >= 3) return true;
            }
            return false;
        });

        candidateGenes.forEach(tGene => {
            const transData = this.translocations.geneData[tGene].translocations;

            const wt = filteredData.filter(d => (transData[d.cellLineId] || 0) === 0);
            const fused = filteredData.filter(d => (transData[d.cellLineId] || 0) >= 1);

            if (wt.length >= 3 && fused.length >= 3) {
                const wtStats = this.pearsonWithSlope(wt.map(d => d.x), wt.map(d => d.y));
                const fusedStats = this.pearsonWithSlope(fused.map(d => d.x), fused.map(d => d.y));

                const deltaR = fusedStats.correlation - wtStats.correlation;
                const deltaSlope = fusedStats.slope - wtStats.slope;

                const z1 = 0.5 * Math.log((1 + wtStats.correlation) / (1 - wtStats.correlation));
                const z2 = 0.5 * Math.log((1 + fusedStats.correlation) / (1 - fusedStats.correlation));
                const se = Math.sqrt(1/(wt.length - 3) + 1/(fused.length - 3));
                const zDiff = (z2 - z1) / se;
                const pR = 2 * (1 - this.normalCDF(Math.abs(zDiff)));

                tableData.push({
                    tGene,
                    nWT: wt.length,
                    rWT: wtStats.correlation,
                    slopeWT: wtStats.slope,
                    nFused: fused.length,
                    rFused: fusedStats.correlation,
                    slopeFused: fusedStats.slope,
                    deltaR,
                    pR,
                    deltaSlope
                });
            }
        });

        tableData.sort((a, b) => a.pR - b.pR);

        const filterInfo = filterDesc ? `<p style="font-size: 11px; color: #333; margin-bottom: 8px; background: #f0f9ff; padding: 4px 8px; border-radius: 4px;"><b>Filter:</b> ${filterDesc}</p>` : '';
        let html = `
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                <h4 style="margin: 0;">Fusions affecting ${gene1} vs ${gene2}</h4>
                <div>
                    <button class="btn btn-primary btn-sm" id="backToGraphBtnTrans" style="margin-right: 8px;">← Back to Graph</button>
                    <button class="btn btn-success btn-sm" id="downloadTransCompareCSV">Download CSV</button>
                </div>
            </div>
            ${filterInfo}
            <p style="font-size: 11px; color: #666; margin-bottom: 6px;">
                Comparing WT (0 fusions) vs Fused (1+ fusions). Sorted by p-value.
            </p>
            <p style="font-size: 10px; color: #0c4a6e; background: #f0f9ff; padding: 4px 8px; border-radius: 4px; margin-bottom: 8px;">
                <b>Statistics:</b> p(Δr) uses Fisher z-transformation to test if correlations differ significantly between WT and fused cells.
            </p>
            <p style="font-size: 11px; color: #059669; margin-bottom: 8px;">Click any row to color scatter plot by that fusion</p>
            <div class="table-container" style="max-height: 380px; overflow-y: auto;">
            <table id="compareTranslocationsTable" class="data-table" style="width: 100%; font-size: 11px;">
                <thead>
                    <tr>
                        <th data-sort="tGene" data-type="string" style="cursor: pointer;">Fusion Gene ↕</th>
                        <th data-sort="nWT" data-type="number" style="cursor: pointer; border-left: 2px solid #2563eb;">N(WT) ↕</th>
                        <th data-sort="rWT" data-type="number" style="cursor: pointer;">r(WT) ↕</th>
                        <th data-sort="nFused" data-type="number" style="cursor: pointer; border-left: 2px solid #dc2626;">N(Fused) ↕</th>
                        <th data-sort="rFused" data-type="number" style="cursor: pointer;">r(Fused) ↕</th>
                        <th data-sort="deltaR" data-type="number" style="cursor: pointer; border-left: 2px solid #6b7280;">Δr ↕</th>
                        <th data-sort="pR" data-type="number" style="cursor: pointer;">p(Δr) ↕</th>
                    </tr>
                </thead>
                <tbody>
        `;

        tableData.forEach(row => {
            const deltaRColor = row.deltaR < 0 ? '#dc2626' : '#5a9f4a';
            const pHighlight = row.pR < 0.05 ? 'background: #fef3c7;' : '';

            html += `
                <tr class="clickable-trans-row" data-trans-gene="${row.tGene}" style="${pHighlight} cursor: pointer;">
                    <td><b>${row.tGene}</b></td>
                    <td style="text-align: center; border-left: 2px solid #2563eb;">${row.nWT}</td>
                    <td style="text-align: center;">${row.rWT.toFixed(3)}</td>
                    <td style="text-align: center; border-left: 2px solid #dc2626;">${row.nFused}</td>
                    <td style="text-align: center;">${row.rFused.toFixed(3)}</td>
                    <td style="text-align: center; border-left: 2px solid #6b7280; color: ${deltaRColor}; font-weight: 600;">${row.deltaR.toFixed(3)}</td>
                    <td style="text-align: center;">${row.pR.toExponential(1)}</td>
                </tr>
            `;
        });

        html += `
                </tbody>
            </table>
            </div>
            <p style="font-size: 11px; color: #666; margin-top: 8px;">
                Yellow = p < 0.05. This analysis may be biased as fusions select for cancer types.
            </p>
        `;

        document.getElementById('compareTable').innerHTML = html;

        document.getElementById('backToGraphBtnTrans')?.addEventListener('click', () => {
            document.getElementById('compareTable').style.display = 'none';
            document.getElementById('scatterPlot').style.display = 'block';
        });

        document.getElementById('downloadTransCompareCSV')?.addEventListener('click', () => {
            let csv = 'Fusion_Gene,N_WT,r_WT,slope_WT,N_Fused,r_Fused,slope_Fused,Delta_r,p_Delta_r,Delta_slope\n';
            tableData.forEach(row => {
                csv += `${row.tGene},${row.nWT},${row.rWT.toFixed(4)},${row.slopeWT.toFixed(4)},${row.nFused},${row.rFused.toFixed(4)},${row.slopeFused.toFixed(4)},${row.deltaR.toFixed(4)},${row.pR.toExponential(2)},${row.deltaSlope.toFixed(4)}\n`;
            });
            this.downloadFile(csv, `${gene1}_vs_${gene2}_fusion_comparison.csv`, 'text/csv');
        });

        document.querySelectorAll('#compareTranslocationsTable .clickable-trans-row').forEach(row => {
            row.addEventListener('click', () => {
                const tGene = row.dataset.transGene;
                document.getElementById('translocationGene').value = tGene;
                document.getElementById('translocationMode').value = 'color';
                document.getElementById('compareTable').style.display = 'none';
                document.getElementById('scatterPlot').style.display = 'block';
                this.updateInspectPlot();
            });
        });

        this.setupSortableTable('compareTranslocationsTable');
    }

    updateInspectGenes() {
        const gene1 = document.getElementById('inspectGeneX').value.trim().toUpperCase();
        const gene2 = document.getElementById('inspectGeneY').value.trim().toUpperCase();

        if (!gene1 || !gene2) {
            alert('Please enter both X and Y genes.');
            return;
        }

        if (!this.geneIndex.has(gene1)) {
            alert(`Gene "${gene1}" not found in the dataset.`);
            return;
        }
        if (!this.geneIndex.has(gene2)) {
            alert(`Gene "${gene2}" not found in the dataset.`);
            return;
        }

        // Re-open inspect with the new genes
        this.openInspect({ gene1, gene2, correlation: null });

        // Update title
        document.getElementById('inspectTitle').textContent = `Correlation: ${gene1} vs ${gene2}`;
    }

    showCompareAllCancerTypes() {
        if (!this.currentInspect) return;

        const { gene1, gene2, data } = this.currentInspect;

        // Apply current filters (mutation filter only, since we're stratifying by cancer)
        const mutFilterGene = document.getElementById('mutationFilterGene').value;
        const mutFilterLevel = document.getElementById('mutationFilterLevel').value;

        let filteredData = [...data];

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
        let filterDesc = '';
        if (mutFilterGene && mutFilterLevel !== 'all') {
            const levelText = mutFilterLevel === '0' ? 'WT' :
                              mutFilterLevel === '1' ? '1 mut' :
                              mutFilterLevel === '2' ? '2 mut' : '1+2 mut';
            filterDesc = `${mutFilterGene}: ${levelText}`;
        }

        // Show compare table
        document.getElementById('scatterPlot').style.display = 'none';
        document.getElementById('compareTable').style.display = 'block';

        this.renderCancerTypeComparisonTable(filteredData, gene1, gene2, filterDesc);
    }

    renderCancerTypeComparisonTable(filteredData, gene1, gene2, filterDesc = '') {
        // Group data by cancer type
        const cancerGroups = {};
        filteredData.forEach(d => {
            const lineage = d.lineage || 'Unknown';
            if (!cancerGroups[lineage]) cancerGroups[lineage] = [];
            cancerGroups[lineage].push(d);
        });

        const tableData = [];

        Object.entries(cancerGroups).forEach(([cancerType, data]) => {
            if (data.length >= 5) { // Need at least 5 samples for meaningful correlation
                const xVals = data.map(d => d.x);
                const yVals = data.map(d => d.y);
                const stats = this.pearsonWithSlope(xVals, yVals);

                // Calculate mean gene effects
                const meanX = xVals.reduce((a, b) => a + b, 0) / xVals.length;
                const meanY = yVals.reduce((a, b) => a + b, 0) / yVals.length;

                tableData.push({
                    cancerType,
                    n: data.length,
                    correlation: stats.correlation,
                    slope: stats.slope,
                    meanX,
                    meanY
                });
            }
        });

        // Sort by absolute correlation by default
        tableData.sort((a, b) => Math.abs(b.correlation) - Math.abs(a.correlation));

        // Build HTML table
        let html = `
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                <h4 style="margin: 0;">${gene1} vs ${gene2} - By Cancer Type</h4>
                <div>
                    <button class="btn btn-primary btn-sm" id="backToGraphBtn2" style="margin-right: 8px;">← Back to Graph</button>
                    <button class="btn btn-success btn-sm" id="downloadCancerCompareCSV">Download CSV</button>
                </div>
            </div>
            ${filterDesc ? `<p style="font-size: 11px; color: #666; margin-bottom: 10px;">Filter: ${filterDesc}</p>` : ''}
            <p style="font-size: 11px; color: #059669; margin-bottom: 8px;">Click any row to view that cancer type in the scatter plot</p>
            <div class="table-container" style="max-height: 400px; overflow-y: auto;">
            <table class="data-table" id="compareCancerTypesTable" style="font-size: 12px;">
                <thead>
                    <tr>
                        <th data-sort="cancerType" data-type="string" style="cursor: pointer;">Cancer Type ↕</th>
                        <th data-sort="n" data-type="number" style="cursor: pointer;">N ↕</th>
                        <th data-sort="correlation" data-type="number" style="cursor: pointer;">Correlation ↕</th>
                        <th data-sort="slope" data-type="number" style="cursor: pointer;">Slope ↕</th>
                        <th data-sort="meanX" data-type="number" style="cursor: pointer;">${gene1} Mean ↕</th>
                        <th data-sort="meanY" data-type="number" style="cursor: pointer;">${gene2} Mean ↕</th>
                    </tr>
                </thead>
                <tbody>
        `;

        tableData.forEach(row => {
            const rColor = Math.abs(row.correlation) >= 0.5 ? (row.correlation > 0 ? '#16a34a' : '#dc2626') : '#374151';
            const xColor = row.meanX < -0.5 ? '#dc2626' : row.meanX > 0 ? '#16a34a' : '#374151';
            const yColor = row.meanY < -0.5 ? '#dc2626' : row.meanY > 0 ? '#16a34a' : '#374151';
            html += `
                <tr class="clickable-row" data-cancer-type="${row.cancerType}" style="cursor: pointer;">
                    <td>${row.cancerType}</td>
                    <td style="text-align: center;">${row.n}</td>
                    <td style="text-align: center; font-weight: 600; color: ${rColor};">${row.correlation.toFixed(3)}</td>
                    <td style="text-align: center;">${row.slope.toFixed(3)}</td>
                    <td style="text-align: center; color: ${xColor};">${row.meanX.toFixed(2)}</td>
                    <td style="text-align: center; color: ${yColor};">${row.meanY.toFixed(2)}</td>
                </tr>
            `;
        });

        html += `
                </tbody>
            </table>
            </div>
            <p style="font-size: 11px; color: #666; margin-top: 8px;">
                Showing cancer types with N ≥ 5. Correlations ≥ |0.5| are highlighted.
            </p>
        `;

        document.getElementById('compareTable').innerHTML = html;

        // Add back to graph handler
        document.getElementById('backToGraphBtn2')?.addEventListener('click', () => {
            document.getElementById('compareTable').style.display = 'none';
            document.getElementById('scatterPlot').style.display = 'block';
        });

        // Add download handler
        document.getElementById('downloadCancerCompareCSV')?.addEventListener('click', () => {
            let csv = `Cancer_Type,N,Correlation,Slope,${gene1}_Mean,${gene2}_Mean\n`;
            tableData.forEach(row => {
                csv += `"${row.cancerType}",${row.n},${row.correlation.toFixed(4)},${row.slope.toFixed(4)},${row.meanX.toFixed(2)},${row.meanY.toFixed(2)}\n`;
            });
            this.downloadFile(csv, `${gene1}_vs_${gene2}_by_cancer_type.csv`, 'text/csv');
        });

        // Add row click handlers to filter by cancer type
        document.querySelectorAll('#compareCancerTypesTable .clickable-row').forEach(row => {
            row.addEventListener('click', () => {
                const cancerType = row.dataset.cancerType;
                // Set the cancer filter and go back to scatter plot
                document.getElementById('scatterCancerFilter').value = cancerType;
                document.getElementById('compareTable').style.display = 'none';
                document.getElementById('scatterPlot').style.display = 'block';
                this.updateScatterSubtypeFilter();
                this._styleActiveFilters();
                this.updateInspectPlot();
            });
        });

        // Make table sortable
        this.setupSortableTable('compareCancerTypesTable');
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

    async _exportScatterChart(format) {
        const plotEl = document.getElementById('scatterPlot');
        if (!plotEl || !plotEl.data) return;

        const hotspotGene = document.getElementById('hotspotGene').value;
        const hotspotMode = document.getElementById('hotspotMode').value;
        const transGene = document.getElementById('translocationGene').value;
        const transMode = document.getElementById('translocationMode').value;

        let suffix = '';
        if (hotspotGene && hotspotMode !== 'none') suffix = `_${hotspotGene}`;
        else if (transGene && transMode !== 'none') suffix = `_${transGene}`;
        const filename = `scatter_${this.currentInspect.gene1}_vs_${this.currentInspect.gene2}${suffix}`;

        // Read export settings
        const transparentBg = document.getElementById('exportTransparentBg')?.checked;
        const whitePlotBg = document.getElementById('exportWhitePlotBg')?.checked;

        // Determine export backgrounds
        const origPaperBg = plotEl.layout.paper_bgcolor || '#fff';
        const origPlotBg = plotEl.layout.plot_bgcolor || '#fafafa';
        const exportPaperBg = transparentBg ? 'rgba(0,0,0,0)' : origPaperBg;
        const exportPlotBg = transparentBg ? 'rgba(0,0,0,0)' : (whitePlotBg ? '#fff' : origPlotBg);

        // Temporarily set export backgrounds on the live plot so Plotly.toImage renders them
        await Plotly.relayout(plotEl, { paper_bgcolor: exportPaperBg, plot_bgcolor: exportPlotBg });

        // Always export as SVG first so we can post-process (remove legend clipPath),
        // then convert to PNG via canvas if needed.
        const svgDataUrl = await Plotly.toImage(plotEl, {
            format: 'svg',
            width: plotEl.layout.width,
            height: plotEl.layout.height
        });

        // Restore original backgrounds
        await Plotly.relayout(plotEl, { paper_bgcolor: origPaperBg, plot_bgcolor: origPlotBg });

        // Decode SVG
        let svgString;
        if (svgDataUrl.indexOf('base64,') > -1) {
            svgString = atob(svgDataUrl.split('base64,')[1]);
        } else {
            svgString = decodeURIComponent(svgDataUrl.split(',').slice(1).join(','));
        }

        // Fix legend: remove clipPath, measure text with canvas, widen rect to fit
        svgString = svgString.replace(/clip-path="url\(#legend[^"]*\)"/g, '');
        const parser = new DOMParser();
        const svgDoc = parser.parseFromString(svgString, 'image/svg+xml');
        const legendBg = svgDoc.querySelector('.legend rect.bg');
        const legendTexts = svgDoc.querySelectorAll('.legend .legendtext, .legend .legendtitletext');
        if (legendBg && legendTexts.length) {
            // Use canvas.measureText with Arial for accurate width
            const measureCanvas = document.createElement('canvas');
            const mctx = measureCanvas.getContext('2d');
            let maxRight = 0;
            legendTexts.forEach(t => {
                const fs = t.style.fontSize || '11px';
                const x = parseFloat(t.getAttribute('x')) || 0;
                mctx.font = `${fs} "Arial", "Helvetica", sans-serif`;
                const right = x + mctx.measureText(t.textContent).width;
                if (right > maxRight) maxRight = right;
            });
            if (maxRight > 0) {
                // Scale up by 1.3x — canvas.measureText uses a fallback font that's
                // narrower than the actual font rendered in the SVG
                const newWidth = Math.ceil(maxRight * 1.3);
                const oldWidth = parseFloat(legendBg.getAttribute('width'));
                if (newWidth > oldWidth) {
                    legendBg.setAttribute('width', String(newWidth));
                    // Expand SVG width if legend extends beyond edge
                    const svgEl = svgDoc.documentElement;
                    const legendTransform = svgDoc.querySelector('.legend')?.getAttribute('transform');
                    const tMatch = legendTransform?.match(/translate\(([\d.]+)/);
                    if (tMatch) {
                        const legendX = parseFloat(tMatch[1]);
                        const needed = legendX + newWidth + 5;
                        const svgW = parseFloat(svgEl.getAttribute('width'));
                        if (needed > svgW) {
                            const w = Math.ceil(needed);
                            svgEl.setAttribute('width', String(w));
                            svgEl.setAttribute('viewBox', `0 0 ${w} ${svgEl.getAttribute('height')}`);
                            // Expand paper background rect
                            const paperRect = svgEl.querySelector('rect');
                            if (paperRect) paperRect.setAttribute('width', String(w));
                        }
                    }
                }
            }
            // Remove clipPath definitions too
            svgDoc.querySelectorAll('clipPath[id^="legend"]').forEach(cp => cp.remove());
        }
        svgString = new XMLSerializer().serializeToString(svgDoc.documentElement);

        const state = this.captureAppState();
        const stateJson = JSON.stringify(state);
        const meta = this._buildExportMetadata('scatter', { gene1: this.currentInspect?.gene1, gene2: this.currentInspect?.gene2 });
        const metaJson = JSON.stringify(meta);

        const a = document.createElement('a');
        if (format === 'svg') {
            // Embed app state + metadata for re-opening
            svgString = svgString.replace('</svg>', `<metadata><correlate-state>${stateJson}</correlate-state><correlate-meta>${metaJson}</correlate-meta></metadata></svg>`);
            svgString = await this._finalizeSvgForExport(svgString);
            const blob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
            a.href = URL.createObjectURL(blob);
            a.download = filename + '.svg';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
        } else {
            // Render the fixed SVG to canvas at 4x for publication quality PNG
            const svgBlob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
            const svgUrl = URL.createObjectURL(svgBlob);
            const img = new Image();
            img.onload = async () => {
                const scale = 4;
                const canvas = document.createElement('canvas');
                canvas.width = img.naturalWidth * scale;
                canvas.height = img.naturalHeight * scale;
                const ctx = canvas.getContext('2d');
                ctx.scale(scale, scale);
                if (!transparentBg) {
                    ctx.fillStyle = 'white';
                    ctx.fillRect(0, 0, img.naturalWidth, img.naturalHeight);
                }
                ctx.drawImage(img, 0, 0);
                URL.revokeObjectURL(svgUrl);
                const pngDataUrl = canvas.toDataURL('image/png');
                const pngResp = await fetch(pngDataUrl);
                const pngBuf = await pngResp.arrayBuffer();
                const pngWithMeta = this._addPngTextChunk(pngBuf, 'correlate-state', stateJson);
                const blob = new Blob([pngWithMeta], { type: 'image/png' });
                a.href = URL.createObjectURL(blob);
                a.download = filename + '.png';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                URL.revokeObjectURL(a.href);
            };
            img.src = svgUrl;
        }
    }

    downloadScatterPNG() {
        this._exportScatterChart('png');
    }

    downloadScatterSVG() {
        this._exportScatterChart('svg');
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
                    <tr style="background-color: #5a9f4a; color: white;">
                        <th data-col="0" style="padding: 6px; border: 1px solid #5a9f4a; text-align: left; position: sticky; top: 0; background-color: #5a9f4a; cursor: pointer;">Lineage ▼</th>
                        <th data-col="1" style="padding: 6px; border: 1px solid #5a9f4a; text-align: center; position: sticky; top: 0; background-color: #5a9f4a; cursor: pointer;">N</th>
                        <th data-col="2" style="padding: 6px; border: 1px solid #5a9f4a; text-align: center; position: sticky; top: 0; background-color: #5a9f4a; cursor: pointer;">Corr</th>
                        <th data-col="3" style="padding: 6px; border: 1px solid #5a9f4a; text-align: center; position: sticky; top: 0; background-color: #5a9f4a; cursor: pointer;">Slope</th>
                        <th data-col="4" style="padding: 6px; border: 1px solid #5a9f4a; text-align: center; position: sticky; top: 0; background-color: #5a9f4a; cursor: pointer;">${gene1} (mean)</th>
                        <th data-col="5" style="padding: 6px; border: 1px solid #5a9f4a; text-align: center; position: sticky; top: 0; background-color: #5a9f4a; cursor: pointer;">${gene1} (SD)</th>
                        <th data-col="6" style="padding: 6px; border: 1px solid #5a9f4a; text-align: center; position: sticky; top: 0; background-color: #5a9f4a; cursor: pointer;">${gene2} (mean)</th>
                        <th data-col="7" style="padding: 6px; border: 1px solid #5a9f4a; text-align: center; position: sticky; top: 0; background-color: #5a9f4a; cursor: pointer;">${gene2} (SD)</th>
                    </tr>
                </thead>
                <tbody>
        `;

        tissueStats.forEach(t => {
            const corrColor = t.correlation > 0 ?
                `rgba(34, 197, 94, ${Math.min(1, Math.abs(t.correlation))})` :
                `rgba(239, 68, 68, ${Math.min(1, Math.abs(t.correlation))})`;
            const escapedTissue = t.tissue.replace(/'/g, "\\'");
            tableHtml += `
                <tr style="cursor: pointer;" onclick="app.switchToInspectWithTissue('${escapedTissue}')" title="Click to view scatter plot filtered by ${t.tissue}" onmouseover="this.style.backgroundColor='#f0f9ff'" onmouseout="this.style.backgroundColor=''">
                    <td style="padding: 5px; border: 1px solid #ddd; color: #0066cc; text-decoration: underline;">${t.tissue}</td>
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
        tableHtml += '<p style="font-size: 11px; color: #666; margin-top: 8px; font-style: italic;">Click a lineage to view its scatter plot</p>';

        // Hide inspect-layout, use byTissueContainer for full-width layout
        document.querySelector('.inspect-layout').style.display = 'none';
        const tissueContainer = document.getElementById('byTissueContainer');

        const chartHeight = Math.max(400, tissueStats.length * 22 + 80);
        tissueContainer.style.display = 'block';
        tissueContainer.innerHTML = `
            <div style="display: flex; gap: 20px; align-items: flex-start;">
                <div style="flex: 1 1 50%; min-width: 400px;">
                    <div id="byTissueChart" style="height: ${chartHeight}px;"></div>
                </div>
                <div style="flex: 1 1 50%; min-width: 0; overflow-x: auto;">
                    ${tableHtml}
                </div>
            </div>
        `;

        // Render bar chart in the new container
        Plotly.newPlot('byTissueChart', [trace], { ...layout, height: chartHeight }, { responsive: true });

        // Add sortable table headers
        this.setupSortableTable('byTissueTable');
    }

    switchToInspectWithTissue(tissue) {
        // Switch from By Tissue view back to Inspect scatter plot with the selected tissue preset
        if (!this.currentInspect || !this.currentInspect.data) return;

        const { gene1, gene2, data } = this.currentInspect;

        // Restore inspect layout and hide tissue container
        document.querySelector('.inspect-layout').style.display = '';
        document.getElementById('byTissueContainer').style.display = 'none';
        document.getElementById('scatterPlot').style.display = '';

        // Show the scatter plot controls again
        document.querySelector('.inspect-controls').style.display = '';

        // Show scatter-specific download buttons
        document.getElementById('downloadScatterPNG').style.display = '';
        document.getElementById('downloadScatterSVG').style.display = '';
        document.getElementById('downloadScatterCSV').style.display = '';

        // Hide tissue-specific download buttons
        document.getElementById('downloadTissuePNG').style.display = 'none';
        document.getElementById('downloadTissueSVG').style.display = 'none';
        document.getElementById('downloadTissueCSV').style.display = 'none';

        // Set axis range inputs from defaults
        if (this.currentInspect.defaultXlim && this.currentInspect.defaultYlim) {
            document.getElementById('scatterXmin').value = this.currentInspect.defaultXlim[0].toFixed(1);
            document.getElementById('scatterXmax').value = this.currentInspect.defaultXlim[1].toFixed(1);
            document.getElementById('scatterYmin').value = this.currentInspect.defaultYlim[0].toFixed(1);
            document.getElementById('scatterYmax').value = this.currentInspect.defaultYlim[1].toFixed(1);
        }

        // Populate cancer filter with counts (same as openInspect)
        const cancerFilter = document.getElementById('scatterCancerFilter');
        const cancerBox = document.getElementById('cancerFilterBox');
        const lineageCounts = {};
        data.forEach(d => {
            if (d.lineage) {
                lineageCounts[d.lineage] = (lineageCounts[d.lineage] || 0) + 1;
            }
        });
        const lineages = Object.keys(lineageCounts).sort();
        if (lineages.length > 0) {
            cancerFilter.innerHTML = `<option value="">All tissues (n=${data.length})</option>`;
            lineages.forEach(l => {
                cancerFilter.innerHTML += `<option value="${l}">${l} (n=${lineageCounts[l]})</option>`;
            });
            // Set the filter to the selected tissue
            cancerFilter.value = tissue;
            cancerBox.style.display = 'block';
            // Update subtype filter for the selected tissue
            this.updateScatterSubtypeFilter();
            this._styleActiveFilters();
        } else {
            cancerBox.style.display = 'none';
            document.getElementById('scatterSubtypeFilter').style.display = 'none';
        }

        // Populate hotspot genes (same as openInspect)
        const hotspotSelect = document.getElementById('hotspotGene');
        const mutFilterGeneSelect = document.getElementById('mutationFilterGene');
        const cellLinesInPlot = new Set(data.map(d => d.cellLineId));

        if (this.mutations?.genes?.length > 0) {
            hotspotSelect.innerHTML = '<option value="">Select gene...</option>';
            mutFilterGeneSelect.innerHTML = '<option value="">No filter</option>';
            this.mutations.genes.forEach(g => {
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

        // Populate translocation/fusion selectors (datalists, sorted by count desc)
        const transGeneInput2 = document.getElementById('translocationGene');
        const transFilterGeneInput2 = document.getElementById('translocationFilterGene');
        const transGeneDatalist2 = document.getElementById('translocationGeneList');
        const transFilterGeneDatalist2 = document.getElementById('translocationFilterGeneList');

        if (this.translocations?.genes?.length > 0) {
            transGeneInput2.value = '';
            transFilterGeneInput2.value = '';
            const geneCounts2 = [];
            for (const g of this.translocations.genes) {
                const transData = this.translocations.geneData?.[g]?.translocations || {};
                let count = 0;
                for (const cl of cellLinesInPlot) {
                    if (transData[cl] && transData[cl] > 0) count++;
                }
                if (count > 0) geneCounts2.push({ gene: g, count });
            }
            geneCounts2.sort((a, b) => {
                const aPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(a.gene) ? 1 : 0;
                const bPri = CorrelationExplorer.PRIORITY_FUSION_GENES.has(b.gene) ? 1 : 0;
                if (aPri !== bPri) return bPri - aPri;
                return b.count - a.count;
            });

            let transHtml2 = '';
            geneCounts2.forEach(({ gene, count }) => {
                transHtml2 += `<option value="${gene}">${gene} (${count} fused)</option>`;
            });
            transGeneDatalist2.innerHTML = transHtml2;
            transFilterGeneDatalist2.innerHTML = transHtml2;
            document.getElementById('translocationBox').style.display = 'block';
            document.getElementById('translocationFilterBox').style.display = 'block';
            document.getElementById('compareAllTranslocationsBtn').style.display = '';
        } else {
            document.getElementById('translocationBox').style.display = 'none';
            document.getElementById('translocationFilterBox').style.display = 'none';
            document.getElementById('compareAllTranslocationsBtn').style.display = 'none';
        }

        // Reset hotspot mode to default
        document.getElementById('hotspotMode').value = 'none';
        document.getElementById('mutationFilterLevel').value = 'all';

        // Calculate stats for ALL cells (unfiltered) for the title
        const allCellsStats = this.pearsonWithSlope(data.map(d => d.x), data.map(d => d.y));
        document.getElementById('inspectTitle').textContent =
            `${gene1} vs ${gene2} | r=${this.formatNum(allCellsStats.correlation)}, slope=${this.formatNum(allCellsStats.slope)}, n=${data.length} (all cells)`;

        // Show the scatter plot and hide compareTable
        document.getElementById('scatterPlot').style.display = 'block';
        document.getElementById('compareTable').style.display = 'none';

        // Re-render the scatter plot with the filter
        this.updateInspectPlot();
    }

    setupSortableTable(tableId) {
        const table = document.getElementById(tableId);
        if (!table) return;

        // Support both data-col and data-sort formats
        const headers = table.querySelectorAll('th[data-col], th[data-sort]');
        const allHeaders = table.querySelectorAll('th');

        headers.forEach(th => {
            th.addEventListener('click', () => {
                // Get column index - either from data-col or from position
                let col;
                if (th.dataset.col !== undefined) {
                    col = parseInt(th.dataset.col);
                } else {
                    col = Array.from(allHeaders).indexOf(th);
                }

                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                const isAsc = th.dataset.dir !== 'asc';
                th.dataset.dir = isAsc ? 'asc' : 'desc';

                // Update header arrows
                allHeaders.forEach(h => {
                    h.textContent = h.textContent.replace(/ [▲▼↕]$/, '');
                    if (h !== th && (h.dataset.col !== undefined || h.dataset.sort !== undefined)) {
                        h.textContent += ' ↕';
                    }
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
        let header = 'CellLine,CellLineID,Lineage,Subtype,Gene1_Effect,Gene2_Effect';
        if (hotspotGene && this.mutations?.geneData?.[hotspotGene]) {
            header += `,${hotspotGene}_mutation`;
        }
        header += '\n';

        let csv = header;
        const mutationData = hotspotGene && this.mutations?.geneData?.[hotspotGene]?.mutations;

        this.currentInspect.data.forEach(d => {
            const subtype = this.cellLineMetadata?.primaryDisease?.[d.cellLineId] || '';
            let row = `"${d.cellLineName}","${d.cellLineId}","${d.lineage}","${subtype}",${d.x},${d.y}`;
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
        this._exportTissueChart('png');
    }

    downloadTissueChartSVG() {
        this._exportTissueChart('svg');
    }

    async _exportTissueChart(format) {
        if (!this.currentInspect) return;
        const chartEl = document.getElementById('byTissueChart');
        if (!chartEl) return;

        const filename = `by_tissue_${this.currentInspect.gene1}_vs_${this.currentInspect.gene2}`;
        const w = 800;
        const h = Math.max(400, (this.currentTissueStats?.length || 10) * 25 + 100);
        const meta = this._buildExportMetadata('tissue_chart', {
            gene1: this.currentInspect.gene1,
            gene2: this.currentInspect.gene2
        });
        const metaJson = JSON.stringify(meta);

        const svgDataUrl = await Plotly.toImage(chartEl, { format: 'svg', width: w, height: h });
        let svgStr;
        if (svgDataUrl.indexOf('base64,') > -1) svgStr = atob(svgDataUrl.split('base64,')[1]);
        else svgStr = decodeURIComponent(svgDataUrl.split(',').slice(1).join(','));

        const a = document.createElement('a');
        if (format === 'svg') {
            svgStr = svgStr.replace('</svg>', `<metadata><correlate-meta>${metaJson}</correlate-meta></metadata></svg>`);
            svgStr = await this._finalizeSvgForExport(svgStr);
            const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
            a.href = URL.createObjectURL(blob);
            a.download = `${filename}.svg`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
        } else {
            const svgBlob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
            const svgUrl = URL.createObjectURL(svgBlob);
            const img = new Image();
            img.onload = async () => {
                const scale = 4;
                const canvas = document.createElement('canvas');
                canvas.width = img.naturalWidth * scale;
                canvas.height = img.naturalHeight * scale;
                const ctx = canvas.getContext('2d');
                ctx.scale(scale, scale);
                ctx.fillStyle = 'white';
                ctx.fillRect(0, 0, img.naturalWidth, img.naturalHeight);
                ctx.drawImage(img, 0, 0);
                URL.revokeObjectURL(svgUrl);
                const pngDataUrl = canvas.toDataURL('image/png');
                const pngResp = await fetch(pngDataUrl);
                const pngBuf = await pngResp.arrayBuffer();
                const pngWithMeta = this._addPngTextChunk(pngBuf, 'correlate-meta', metaJson);
                const blob = new Blob([pngWithMeta], { type: 'image/png' });
                a.href = URL.createObjectURL(blob);
                a.download = `${filename}.png`;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                URL.revokeObjectURL(a.href);
            };
            img.src = svgUrl;
        }
    }

    downloadTissueTableCSV() {
        if (!this.currentInspect || !this.currentTissueStats) return;

        const { gene1, gene2 } = this.currentInspect;

        let csv = `Lineage,N,Correlation,${gene1}_Effect_mean,${gene1}_Effect_SD,${gene2}_Effect_mean,${gene2}_Effect_SD\n`;

        this.currentTissueStats.forEach(t => {
            csv += `"${t.tissue}",${t.n},${t.correlation.toFixed(4)},${t.meanX.toFixed(2)},${t.sdX.toFixed(2)},${t.meanY.toFixed(2)},${t.sdY.toFixed(2)}\n`;
        });

        this.downloadFile(csv,
            `by_tissue_${gene1}_vs_${gene2}.csv`,
            'text/csv');
    }

    // ============================================================
    // Gene Effect Modal Methods
    // ============================================================

    openGeneEffectModal(gene, view = 'tissue') {
        const geneUpper = gene.toUpperCase();
        if (!this.geneIndex.has(geneUpper)) {
            alert(`Gene "${gene}" not found in the dataset.`);
            return;
        }

        // Store current gene and view
        this.currentGeneEffect = {
            gene: geneUpper,
            data: []
        };
        this.currentGEView = view;

        // Get gene effect data for all cell lines
        const geneIdx = this.geneIndex.get(geneUpper);
        const geneData = this.getGeneData(geneIdx);

        for (let i = 0; i < this.nCellLines; i++) {
            if (!isNaN(geneData[i])) {
                const cellLine = this.metadata.cellLines[i];
                this.currentGeneEffect.data.push({
                    cellLineId: cellLine,
                    cellLineName: this.getCellLineName(cellLine),
                    lineage: this.getCellLineLineage(cellLine),
                    geneEffect: geneData[i]
                });
            }
        }

        // Calculate overall stats
        const allEffects = this.currentGeneEffect.data.map(d => d.geneEffect);
        const mean = allEffects.reduce((a, b) => a + b, 0) / allEffects.length;
        const sd = Math.sqrt(allEffects.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / allEffects.length);

        // Update UI
        document.getElementById('geneEffectTitle').textContent = `${geneUpper} - Gene Effect Analysis`;
        document.getElementById('geneEffectSearch').value = geneUpper;
        document.getElementById('geneEffectCurrentGene').textContent = '';
        const geCellLineSearch = document.getElementById('geCellLineSearch');
        if (geCellLineSearch) geCellLineSearch.value = '';
        document.getElementById('geSummaryGene').textContent = geneUpper;
        document.getElementById('geSummaryMean').textContent = mean.toFixed(2);
        document.getElementById('geSummarySD').textContent = sd.toFixed(2);
        document.getElementById('geSummaryN').textContent = allEffects.length;
        document.getElementById('geneEffectSummary').style.display = 'block';

        // Populate tissue filter dropdown
        const tissueFilter = document.getElementById('geTissueFilter');
        if (tissueFilter) {
            const currentValue = tissueFilter.value;
            const lineages = [...new Set(this.currentGeneEffect.data.map(d => d.lineage).filter(Boolean))].sort();
            tissueFilter.innerHTML = '<option value="">All tissues</option>';
            lineages.forEach(l => {
                const opt = document.createElement('option');
                opt.value = l;
                opt.textContent = l;
                tissueFilter.appendChild(opt);
            });
            tissueFilter.value = currentValue || '';
        }

        // Restore visibility of UI elements (in case hidden by mutation analysis inspect)
        document.getElementById('geSearchBar').style.display = '';
        document.getElementById('geTableContainer').style.display = '';
        document.getElementById('geChartContainer').style.flex = '0 0 55%';

        // Mark this as regular gene effect view
        this.geneEffectViewMode = 'geneEffect';

        // Show hotspot/fusion filters in standalone mode (hide mutation-inspect-only controls)
        document.getElementById('geCompareButtons').style.display = 'none';
        document.getElementById('geResetFiltersBtn').style.display = 'none';
        document.getElementById('geInlineCompareTable').style.display = 'none';
        document.getElementById('geHotspotGeneGroup').style.display = 'none';

        // Populate and show hotspot filter
        const hotspotFilterEl = document.getElementById('geHotspotFilter');
        if (hotspotFilterEl && this.mutations?.genes) {
            let hHtml = '<option value="">No hotspot filter</option>';
            for (const g of this.mutations.genes) {
                hHtml += `<option value="${g}">${g}</option>`;
            }
            hotspotFilterEl.innerHTML = hHtml;
            hotspotFilterEl.value = '';
            hotspotFilterEl.style.display = '';
        }
        const warn = document.getElementById('geHotspotBiasWarning');
        if (warn) warn.style.display = 'none';

        // Populate and show fusion filter
        const fusionFilterEl = document.getElementById('geFusionFilter');
        if (fusionFilterEl && this._fusionGeneCounts?.length > 0) {
            let fHtml = '<option value="">No fusion filter</option>';
            for (const { gene, nFused } of this._fusionGeneCounts) {
                fHtml += `<option value="${gene}">${gene} (${nFused} fused)</option>`;
            }
            fusionFilterEl.innerHTML = fHtml;
            fusionFilterEl.value = '';
            fusionFilterEl.style.display = '';
        }

        // Restore view buttons (may have been hidden by mutation inspect)
        document.getElementById('geViewTissue').style.display = '';
        document.getElementById('geViewHotspot').style.display = '';
        const viewLabel = document.getElementById('geViewTissue').previousElementSibling;
        if (viewLabel && viewLabel.textContent.trim() === 'View:') viewLabel.style.display = '';

        // Show modal
        document.getElementById('geneEffectModal').style.display = 'flex';
        this._updateGeOncoprintLabel();

        // Render the selected view
        this.switchGeneEffectView(view);
    }

    _updateGeOncoprintLabel() {
        const el = document.getElementById('geOncoprintLabel');
        const hotspotSelect = document.getElementById('geHotspotFilter');
        const fusionSelect = document.getElementById('geFusionFilter');
        if (!el) return;
        if (this._activeOncoprintFilters && this._activeOncoprintFilters.length > 0) {
            const filteredN = this.getGETissueFilteredData().length;
            const totalN = this.currentGeneEffect?.data?.length || 0;
            const tags = this._activeOncoprintFilters.map(f => {
                const bg = f.state === 'mut' ? '#dcfce7' : '#fef2f2';
                const color = f.state === 'mut' ? '#16a34a' : '#dc2626';
                return `<span style="background:${bg};color:${color};padding:1px 6px;border-radius:10px;font-size:10px;">${f.gene} ${f.state === 'mut' ? 'Mut' : 'WT'}</span>`;
            }).join(' ');
            el.innerHTML = tags + ` <span style="font-size:10px;color:#6b7280;">(${filteredN}/${totalN} cell lines)</span>`;
            el.style.display = 'inline-flex';
            if (hotspotSelect) { hotspotSelect.style.opacity = '0.3'; hotspotSelect.style.pointerEvents = 'none'; }
            if (fusionSelect) { fusionSelect.style.opacity = '0.3'; fusionSelect.style.pointerEvents = 'none'; }
        } else {
            el.style.display = 'none';
            if (hotspotSelect) { hotspotSelect.style.opacity = ''; hotspotSelect.style.pointerEvents = ''; }
            if (fusionSelect) { fusionSelect.style.opacity = ''; fusionSelect.style.pointerEvents = ''; }
        }
    }

    switchGeneEffectView(view) {
        this.currentGEView = view;

        // Reset detailed view state and search
        this.geDetailedView = null;
        this.updateShowAllButton();
        const searchInput = document.getElementById('geTableSearch');
        if (searchInput) searchInput.value = '';

        // Update button styles
        const tissueBtn = document.getElementById('geViewTissue');
        const hotspotBtn = document.getElementById('geViewHotspot');

        // Update statistics explanation text
        const statsExplanation = document.getElementById('geStatsExplanationText');

        if (view === 'tissue') {
            tissueBtn.style.background = '#5a9f4a';
            tissueBtn.style.color = 'white';
            tissueBtn.classList.remove('btn-secondary');
            hotspotBtn.style.background = '';
            hotspotBtn.style.color = '';
            hotspotBtn.classList.add('btn-secondary');
            document.getElementById('geByTissueView').style.display = 'block';
            document.getElementById('geByHotspotView').style.display = 'none';
            if (statsExplanation) statsExplanation.textContent = "p-values: Welch's t-test comparing each cancer type vs all other cell lines.";
            this.renderGeneEffectByTissue();
        } else {
            hotspotBtn.style.background = '#5a9f4a';
            hotspotBtn.style.color = 'white';
            hotspotBtn.classList.remove('btn-secondary');
            tissueBtn.style.background = '';
            tissueBtn.style.color = '';
            tissueBtn.classList.add('btn-secondary');
            document.getElementById('geByTissueView').style.display = 'none';
            document.getElementById('geByHotspotView').style.display = 'block';
            if (statsExplanation) statsExplanation.textContent = "Shows 3 mutation levels: 0 (WT, blue), 1 (orange), 2 (red). p-value: Welch's t-test comparing 1+2 combined vs WT.";
            this.renderGeneEffectByHotspot();
        }
    }

    filterGETable(searchTerm) {
        const tbody = document.getElementById('geneEffectTableBody');
        if (!tbody) return;

        const term = searchTerm.toLowerCase().trim();
        const rows = tbody.querySelectorAll('tr');

        rows.forEach(row => {
            const groupCell = row.querySelector('td:first-child');
            if (groupCell) {
                const text = groupCell.textContent.toLowerCase();
                row.style.display = term === '' || text.includes(term) ? '' : 'none';
            }
        });
    }

    showGeneEffectAnalysis(gene, view = 'tissue') {
        this.openGeneEffectModal(gene, view);
    }

    _applyParamFiltersToGEModal() {
        // Carry active filters (from parameter section OR cell line browser) into the GE modal
        const lineage = document.getElementById('lineageFilter')?.value || document.getElementById('clbTissueFilter')?.value;
        if (lineage) {
            const geTissue = document.getElementById('geTissueFilter');
            if (geTissue) {
                geTissue.value = lineage;
                this.updateGeSubtypeFilter?.();
                const subLineage = document.getElementById('subLineageFilter')?.value || document.getElementById('clbSubtypeFilter')?.value;
                if (subLineage) {
                    const geSub = document.getElementById('geSubtypeFilter');
                    if (geSub) geSub.value = subLineage;
                }
            }
        }
        // Hotspot/translocation — skip if oncoprint handles multi-gene filtering
        if (!this._activeOncoprintFilters || this._activeOncoprintFilters.length === 0) {
            const paramHotspot = document.getElementById('paramHotspotGene')?.value;
            const clbHotspot = document.getElementById('clbHotspotFilter')?.value;
            const hotspot = paramHotspot || clbHotspot;
            if (hotspot) {
                const geHotspot = document.getElementById('geHotspotFilter');
                if (geHotspot) geHotspot.value = hotspot;
            }
            const paramTrans = document.getElementById('paramTranslocationGene')?.value;
            const clbTrans = document.getElementById('clbTranslocationFilter')?.value;
            const trans = paramTrans || clbTrans;
            if (trans) {
                const geFusion = document.getElementById('geFusionFilter');
                if (geFusion) geFusion.value = trans;
            }
        }
        // Show oncoprint filter label in GE modal
        this._updateGeOncoprintLabel();
        // Re-render with the applied filters
        this.switchGeneEffectView(this.currentGEView || 'tissue');
        // Update cell line count to reflect filters
        const filteredData = this.getGETissueFilteredData();
        const totalData = this.currentGeneEffect?.data?.length || 0;
        const nEl = document.getElementById('geSummaryN');
        if (nEl && filteredData.length < totalData) {
            nEl.textContent = `${filteredData.length} of ${totalData}`;
        }
    }

    getGETissueFilteredData() {
        if (!this.currentGeneEffect) return [];
        let data = this.currentGeneEffect.data;
        const tissueFilter = document.getElementById('geTissueFilter')?.value;
        if (tissueFilter) {
            data = data.filter(d => d.lineage === tissueFilter);
        }
        const hotspotGene = document.getElementById('geHotspotFilter')?.value;
        if (hotspotGene && this.mutations?.geneData?.[hotspotGene]) {
            const mutData = this.mutations.geneData[hotspotGene].mutations || {};
            data = data.filter(d => (mutData[d.cellLineId] || 0) >= 1);
        }
        const fusionGene = document.getElementById('geFusionFilter')?.value;
        if (fusionGene && this.translocations?.geneData?.[fusionGene]) {
            const transData = this.translocations.geneData[fusionGene].translocations || {};
            data = data.filter(d => (transData[d.cellLineId] || 0) >= 1);
        }
        // Oncoprint multi-gene filters
        if (this._activeOncoprintFilters && this._activeOncoprintFilters.length > 0) {
            data = data.filter(d => this._cellLinePassesOncoprintFilters(d.cellLineId));
        }
        return data;
    }

    renderGeneEffectByTissue() {
        if (!this.currentGeneEffect) return;

        const data = this.getGETissueFilteredData();
        const gene = this.currentGeneEffect.gene;

        // Determine if we should group by subtype (when a lineage filter is active)
        const tissueFilter = document.getElementById('geTissueFilter')?.value;
        const groupBySubtype = !!tissueFilter && !!this.cellLineMetadata?.primaryDisease;

        // Get all gene effects for comparison
        const allEffects = data.map(d => d.geneEffect);

        // Group data by cancer type or subtype
        const groupedData = {};
        data.forEach(d => {
            const groupKey = groupBySubtype
                ? (this.cellLineMetadata.primaryDisease[d.cellLineId] || 'Unknown')
                : (d.lineage || 'Unknown');
            if (!groupedData[groupKey]) groupedData[groupKey] = [];
            groupedData[groupKey].push({
                geneEffect: d.geneEffect,
                cellLineName: d.cellLineName || d.cellLineId,
                cellLineId: d.cellLineId
            });
        });

        // Lower min group size when mutation/fusion filters are active (small groups expected)
        const hasActiveFilter = document.getElementById('geHotspotFilter')?.value || document.getElementById('geFusionFilter')?.value;
        const minGroupSize = parseInt(document.getElementById('geMinGroupSize')?.value) || 1;

        // Calculate stats for each group including p-value vs all cells
        const stats = [];
        Object.entries(groupedData).forEach(([groupName, cellData]) => {
            if (cellData.length >= minGroupSize) {
                const effects = cellData.map(c => c.geneEffect);
                const mean = effects.reduce((a, b) => a + b, 0) / effects.length;
                const sd = Math.sqrt(effects.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / effects.length);

                // Calculate p-value comparing this group to all other cells
                const otherEffects = data.filter(d => {
                    const key = groupBySubtype
                        ? (this.cellLineMetadata.primaryDisease[d.cellLineId] || 'Unknown')
                        : (d.lineage || 'Unknown');
                    return key !== groupName;
                }).map(d => d.geneEffect);
                const tTest = this.welchTTest(effects, otherEffects);

                stats.push({
                    group: groupName,
                    n: cellData.length,
                    mean,
                    sd,
                    pValue: tTest.p,
                    cellData: cellData
                });
            }
        });

        // Sort by median gene effect for chart display
        stats.sort((a, b) => a.mean - b.mean);
        this.currentGEStats = stats;

        // Create box plot traces for each group
        const traces = stats.map((s, idx) => ({
            type: 'box',
            name: `${s.group} (n=${s.n})`,
            x: s.cellData.map(c => c.geneEffect),
            text: s.cellData.map(c => c.cellLineName),
            boxpoints: 'all',
            jitter: 0.3,
            pointpos: 0,
            marker: {
                color: s.mean < -0.5 ? 'rgba(220, 38, 38, 0.6)' : s.mean > 0 ? 'rgba(34, 197, 94, 0.5)' : 'rgba(107, 114, 128, 0.5)',
                size: 4
            },
            line: { color: s.pValue < 0.05 ? '#1f2937' : '#9ca3af' },
            fillcolor: s.mean < -0.5 ? 'rgba(220, 38, 38, 0.2)' : s.mean > 0 ? 'rgba(34, 197, 94, 0.15)' : 'rgba(156, 163, 175, 0.2)',
            hovertemplate: '<b>%{text}</b><br>Gene Effect: %{x:.3f}<extra></extra>'
        }));

        // Calculate dynamic sizing
        const numEntries = stats.length;
        const tickFontSize = numEntries > 25 ? 7 : numEntries > 15 ? 8 : 9;
        const boxHeight = numEntries > 25 ? 18 : numEntries > 15 ? 22 : 28;
        const chartHeight = Math.max(350, numEntries * boxHeight + 80);

        const layout = {
            title: { text: `${gene} by ${groupBySubtype ? 'Disease Subtype' : 'Cancer Type'}`, font: { size: 13 } },
            xaxis: {
                title: `${gene} Gene Effect`,
                zeroline: true,
                zerolinecolor: '#374151',
                zerolinewidth: 2
            },
            yaxis: {
                automargin: true,
                tickfont: { size: tickFontSize }
            },
            margin: { t: 40, b: 50, l: 10, r: 30 },
            height: chartHeight,
            showlegend: false,
            paper_bgcolor: 'white',
            plot_bgcolor: 'white'
        };

        Plotly.newPlot('geneEffectPlot', traces, layout, { responsive: true, edits: { annotationPosition: true } });

        // Highlight cell line if requested (from CLB gene link or cell line search)
        const highlightCl = this._geHighlightCellLine;
        const searchCl = document.getElementById('geCellLineSearch')?.value?.trim().toUpperCase();
        if (highlightCl || searchCl) {
            const annotations = [];
            for (let i = 0; i < stats.length; i++) {
                const s = stats[i];
                for (const c of s.cellData) {
                    const nameUpper = (c.cellLineName || '').toUpperCase();
                    const idUpper = (c.cellLineId || '').toUpperCase();
                    const isHighlight = (highlightCl && (c.cellLineId === highlightCl));
                    const isSearch = (searchCl && searchCl.length >= 2 && (nameUpper.includes(searchCl) || idUpper.includes(searchCl)));
                    if (isHighlight || isSearch) {
                        annotations.push({
                            x: c.geneEffect, y: `${s.group} (n=${s.n})`,
                            text: c.cellLineName, showarrow: true, arrowhead: 2, arrowsize: 1,
                            arrowcolor: '#dc2626', ax: 40, ay: -25,
                            font: { size: 10, color: '#dc2626', weight: 'bold' },
                            bgcolor: 'rgba(255,255,255,0.85)', borderpad: 2
                        });
                    }
                }
            }
            if (annotations.length > 0) Plotly.relayout('geneEffectPlot', { annotations });
        }

        // Render table (remove cellData from stats for table)
        const tableStats = stats.map(s => ({
            group: s.group,
            n: s.n,
            mean: s.mean,
            sd: s.sd,
            pValue: s.pValue
        }));
        tableStats.sort((a, b) => a.pValue - b.pValue); // Sort by p-value (low to high)
        this.currentGEStats = tableStats;
        this._geGroupBySubtype = groupBySubtype;
        this.renderGETable(tableStats, 'tissue');
    }

    renderGeneEffectByHotspot() {
        if (!this.currentGeneEffect) {
            document.getElementById('geneEffectHotspotPlot').innerHTML = '<div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #6b7280;">No gene selected</div>';
            return;
        }

        if (!this.mutations?.genes || this.mutations.genes.length === 0) {
            document.getElementById('geneEffectHotspotPlot').innerHTML = '<div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #6b7280;">No hotspot mutation data available</div>';
            document.getElementById('geneEffectTableBody').innerHTML = '<tr><td colspan="9" style="text-align: center; padding: 20px; color: #6b7280;">No mutation data</td></tr>';
            return;
        }

        const data = this.getGETissueFilteredData();
        const gene = this.currentGeneEffect.gene;

        // Calculate stats for each hotspot gene, keeping cell-level data for box plots
        // Now showing 3 levels: 0 (WT), 1, and 2 mutations
        const hotspotStats = [];

        this.mutations.genes.forEach(hotspotGene => {
            const mutData = this.mutations.geneData?.[hotspotGene]?.mutations || {};

            const cellData0 = []; // WT (0 mutations)
            const cellData1 = []; // 1 mutation
            const cellData2 = []; // 2 mutations

            data.forEach(d => {
                const mutLevel = mutData[d.cellLineId] || 0;
                const cellInfo = {
                    geneEffect: d.geneEffect,
                    cellLineName: d.cellLineName || d.cellLineId,
                    cellLineId: d.cellLineId
                };
                if (mutLevel === 0) {
                    cellData0.push(cellInfo);
                } else if (mutLevel === 1) {
                    cellData1.push(cellInfo);
                } else {
                    cellData2.push(cellInfo);
                }
            });

            // Require at least 3 in WT and at least 1 mutant cell (combine 1+2 for significance)
            const mutCellData = [...cellData1, ...cellData2];
            if (mutCellData.length >= 1 && cellData0.length >= 3) {
                const effects0 = cellData0.map(c => c.geneEffect);
                const effects1 = cellData1.map(c => c.geneEffect);
                const effects2 = cellData2.map(c => c.geneEffect);
                const effectsMut = mutCellData.map(c => c.geneEffect);

                const mean0 = effects0.reduce((a, b) => a + b, 0) / effects0.length;
                const sd0 = Math.sqrt(effects0.reduce((a, b) => a + Math.pow(b - mean0, 2), 0) / effects0.length);

                const mean1 = effects1.length > 0 ? effects1.reduce((a, b) => a + b, 0) / effects1.length : NaN;
                const sd1 = effects1.length > 1 ? Math.sqrt(effects1.reduce((a, b) => a + Math.pow(b - mean1, 2), 0) / effects1.length) : NaN;

                const mean2 = effects2.length > 0 ? effects2.reduce((a, b) => a + b, 0) / effects2.length : NaN;
                const sd2 = effects2.length > 1 ? Math.sqrt(effects2.reduce((a, b) => a + Math.pow(b - mean2, 2), 0) / effects2.length) : NaN;

                // p-value comparing WT vs 1+2 combined (for ranking)
                const tTest = effectsMut.length >= 3 ? this.welchTTest(effects0, effectsMut) : { p: 1 };
                const diff = effectsMut.length > 0 ? effectsMut.reduce((a, b) => a + b, 0) / effectsMut.length - mean0 : 0;

                hotspotStats.push({
                    group: hotspotGene,
                    n0: cellData0.length,
                    n1: cellData1.length,
                    n2: cellData2.length,
                    nMut: mutCellData.length,
                    mean0,
                    mean1,
                    mean2,
                    sd0,
                    sd1,
                    sd2,
                    diff,
                    pValue: tTest.p,
                    cellData0,
                    cellData1,
                    cellData2
                });
            }
        });

        if (hotspotStats.length === 0) {
            document.getElementById('geneEffectHotspotPlot').innerHTML = '<div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #6b7280;">No hotspots with sufficient data (need ≥3 mutant and ≥3 WT cells)</div>';
            document.getElementById('geneEffectTableBody').innerHTML = '<tr><td colspan="9" style="text-align: center; padding: 20px; color: #6b7280;">Insufficient data</td></tr>';
            return;
        }

        // Sort by p-value (most significant first) for better visualization
        hotspotStats.sort((a, b) => a.pValue - b.pValue);

        // Take top 20 most significant for box plot display
        const topStats = hotspotStats.slice(0, 20);

        // Create box plots for 3 mutation levels (0, 1, 2) for each hotspot
        const traces = [];
        const yCategories = [];

        // Track which legend items need to be shown
        let show0Legend = true, show1Legend = true, show2Legend = true;

        topStats.forEach((s, idx) => {
            const yLabel = `${s.group}`;
            yCategories.push(yLabel);

            // Add traces in reverse order so WT appears first (top) in each group
            // 2 mutations trace (red) - added first, appears at bottom
            if (s.cellData2.length > 0) {
                traces.push({
                    type: 'box',
                    name: '2 mutations',
                    legendgroup: '2',
                    showlegend: show2Legend,
                    y: Array(s.cellData2.length).fill(yLabel),
                    x: s.cellData2.map(c => c.geneEffect),
                    orientation: 'h',
                    boxpoints: 'outliers',
                    marker: { color: '#dc2626', size: 4, outliercolor: '#991b1b' },
                    line: { color: '#991b1b', width: 1.5 },
                    fillcolor: 'rgba(220, 38, 38, 0.6)',
                    hoverinfo: 'x',
                    offsetgroup: '2'
                });
                show2Legend = false;
            }

            // 1 mutation trace (orange)
            if (s.cellData1.length > 0) {
                traces.push({
                    type: 'box',
                    name: '1 mutation',
                    legendgroup: '1',
                    showlegend: show1Legend,
                    y: Array(s.cellData1.length).fill(yLabel),
                    x: s.cellData1.map(c => c.geneEffect),
                    orientation: 'h',
                    boxpoints: 'outliers',
                    marker: { color: '#f97316', size: 4, outliercolor: '#c2410c' },
                    line: { color: '#c2410c', width: 1.5 },
                    fillcolor: 'rgba(249, 115, 22, 0.6)',
                    hoverinfo: 'x',
                    offsetgroup: '1'
                });
                show1Legend = false;
            }

            // WT trace (blue) - 0 mutations - added last, appears at top
            traces.push({
                type: 'box',
                name: '0 (WT)',
                legendgroup: '0',
                showlegend: show0Legend,
                y: Array(s.cellData0.length).fill(yLabel),
                x: s.cellData0.map(c => c.geneEffect),
                orientation: 'h',
                boxpoints: 'outliers',
                marker: { color: '#2563eb', size: 4, outliercolor: '#1e40af' },
                line: { color: '#1e40af', width: 1.5 },
                fillcolor: 'rgba(37, 99, 235, 0.6)',
                hoverinfo: 'x',
                offsetgroup: '0'
            });
            show0Legend = false;
        });

        // Calculate dynamic sizing
        const numEntries = topStats.length;
        const tickFontSize = numEntries > 15 ? 8 : 9;
        const boxHeight = numEntries > 15 ? 35 : 45;
        const chartHeight = Math.max(400, numEntries * boxHeight + 100);

        const layout = {
            title: { text: `${gene} Gene Effect by Hotspot Mutation`, font: { size: 13 } },
            xaxis: {
                title: `${gene} Gene Effect`,
                zeroline: true,
                zerolinecolor: '#374151',
                zerolinewidth: 2
            },
            yaxis: {
                automargin: true,
                tickfont: { size: tickFontSize },
                categoryorder: 'array',
                categoryarray: yCategories.slice().reverse()
            },
            boxmode: 'group',
            boxgap: 0.1,
            boxgroupgap: 0.05,
            margin: { t: 50, b: 50, l: 10, r: 30 },
            height: chartHeight,
            showlegend: true,
            legend: { x: 0.5, y: 1.0, xanchor: 'center', yanchor: 'bottom', orientation: 'h', font: { size: 10 }, bgcolor: 'rgba(255,255,255,0.8)', traceorder: 'reversed' },
            paper_bgcolor: 'white',
            plot_bgcolor: 'white'
        };

        Plotly.newPlot('geneEffectHotspotPlot', traces, layout, { responsive: true, edits: { annotationPosition: true } });

        // Store stats for table (without cell data)
        const tableStats = hotspotStats.map(s => ({
            group: s.group,
            n0: s.n0,
            n1: s.n1,
            n2: s.n2,
            nMut: s.nMut,
            mean0: s.mean0,
            mean1: s.mean1,
            mean2: s.mean2,
            sd0: s.sd0,
            sd1: s.sd1,
            sd2: s.sd2,
            diff: s.diff,
            pValue: s.pValue
        }));

        // Sort by p-value (low to high) for table
        tableStats.sort((a, b) => a.pValue - b.pValue);
        this.currentGEStats = tableStats;

        // Render table
        this.renderGETable(tableStats, 'hotspot');
    }

    openGeneEffectFromNetwork(gene) {
        this.openGeneEffectModal(gene, 'tissue');
        this._applyParamFiltersToGEModal();
    }

    renderGETable(stats, mode) {
        const tbody = document.getElementById('geneEffectTableBody');
        const thead = document.getElementById('geTableHead');
        this.currentGETableMode = mode;

        const headerStyle = 'cursor: pointer; user-select: none;';
        const sortIcon = ' ↕';

        if (mode === 'tissue') {
            const groupLabel = this._geGroupBySubtype ? 'Disease Subtype' : 'Cancer Type';
            thead.innerHTML = `<tr>
                <th style="${headerStyle}" data-sort="group" data-type="string">${groupLabel}${sortIcon}</th>
                <th style="${headerStyle}" data-sort="n" data-type="number">N${sortIcon}</th>
                <th style="${headerStyle}" data-sort="mean" data-type="number">Mean GE${sortIcon}</th>
                <th style="${headerStyle}" data-sort="sd" data-type="number">SD${sortIcon}</th>
                <th style="${headerStyle}" data-sort="pValue" data-type="number">p-value${sortIcon}</th>
            </tr>`;
            this.renderGETableBody(stats, mode);
        } else {
            const hg = this.mutationResults?.hotspotGene || 'Hotspot';
            const isT = this.mutationResults?.isTranslocation;
            const wtLbl = isT ? 'No Fusion' : `${hg} WT`;
            const mutLbl = isT ? 'Fused' : `${hg} Mut`;
            thead.innerHTML = `<tr>
                <th style="${headerStyle}" data-sort="group" data-type="string">Gene${sortIcon}</th>
                <th style="${headerStyle}; border-left: 2px solid #2563eb;" data-sort="n0" data-type="number">N (${wtLbl})${sortIcon}</th>
                <th style="${headerStyle}" data-sort="mean0" data-type="number">GE (${wtLbl})${sortIcon}</th>
                <th style="${headerStyle}; border-left: 2px solid #f97316;" data-sort="n1" data-type="number">N (${mutLbl} 1)${sortIcon}</th>
                <th style="${headerStyle}" data-sort="mean1" data-type="number">GE (${mutLbl} 1)${sortIcon}</th>
                <th style="${headerStyle}; border-left: 2px solid #dc2626;" data-sort="n2" data-type="number">N (${mutLbl} 2${isT ? '+' : ''})${sortIcon}</th>
                <th style="${headerStyle}" data-sort="mean2" data-type="number">GE (${mutLbl} 2${isT ? '+' : ''})${sortIcon}</th>
                <th style="${headerStyle}; border-left: 2px solid #6b7280;" data-sort="diff" data-type="number">Δ GE${sortIcon}</th>
                <th style="${headerStyle}" data-sort="pValue" data-type="number">p-value${sortIcon}</th>
            </tr>`;
            this.renderGETableBody(stats, mode);
        }

        // Add click handlers for sorting (Ctrl/Cmd+click to copy column)
        thead.querySelectorAll('th').forEach(th => {
            th.addEventListener('click', (e) => {
                if (e.ctrlKey || e.metaKey) {
                    const colIndex = Array.from(th.parentNode.children).indexOf(th);
                    this.copyColumnToClipboard(th.closest('table'), colIndex);
                } else {
                    this.sortGETable(th.dataset.sort, th.dataset.type);
                }
            });
            th.title = 'Click to sort. Ctrl/Cmd+click to copy column.';
        });
    }

    renderGETableBody(stats, mode) {
        const tbody = document.getElementById('geneEffectTableBody');
        tbody.innerHTML = '';

        if (mode === 'tissue') {
            stats.forEach(s => {
                const color = s.pValue < 0.05 ? (s.mean < -0.5 ? '#dc2626' : s.mean > 0.5 ? '#16a34a' : '#374151') : '#6b7280';
                const pStr = s.pValue < 0.001 ? '<0.001' : s.pValue.toFixed(3);
                tbody.innerHTML += `<tr class="clickable-row" data-group="${s.group}" style="cursor: pointer;">
                    <td>${s.group}</td>
                    <td style="text-align: center;">${s.n}</td>
                    <td style="text-align: center; font-weight: 500; color: ${color};">${s.mean.toFixed(2)}</td>
                    <td style="text-align: center;">${s.sd.toFixed(2)}</td>
                    <td style="text-align: center; ${s.pValue < 0.05 ? 'font-weight: 600;' : ''}">${pStr}</td>
                </tr>`;
            });
        } else {
            stats.forEach(s => {
                const pStr = s.pValue < 0.001 ? '<0.001' : s.pValue.toFixed(3);
                const mean1Str = isNaN(s.mean1) ? '-' : s.mean1.toFixed(2);
                const mean2Str = isNaN(s.mean2) ? '-' : s.mean2.toFixed(2);
                const diffStr = s.diff !== undefined ? s.diff.toFixed(2) : '-';
                const diffColor = s.diff > 0 ? '#16a34a' : s.diff < 0 ? '#dc2626' : '#374151';
                tbody.innerHTML += `<tr class="clickable-row" data-group="${s.group}" style="cursor: pointer;">
                    <td>${s.group}</td>
                    <td style="text-align: center; color: #2563eb; border-left: 2px solid #2563eb;">${s.n0}</td>
                    <td style="text-align: center; color: #2563eb;">${s.mean0.toFixed(2)}</td>
                    <td style="text-align: center; color: #f97316; border-left: 2px solid #f97316;">${s.n1 || '-'}</td>
                    <td style="text-align: center; color: #f97316;">${mean1Str}</td>
                    <td style="text-align: center; color: #dc2626; border-left: 2px solid #dc2626;">${s.n2 || '-'}</td>
                    <td style="text-align: center; color: #dc2626;">${mean2Str}</td>
                    <td style="text-align: center; color: ${diffColor}; font-weight: 500; border-left: 2px solid #6b7280;">${diffStr}</td>
                    <td style="text-align: center; ${s.pValue < 0.05 ? 'font-weight: 600;' : ''}">${pStr}</td>
                </tr>`;
            });
        }

        // Add click handlers for detailed view
        tbody.querySelectorAll('.clickable-row').forEach(row => {
            row.addEventListener('click', () => {
                const group = row.dataset.group;
                this.showGEDetailedView(group, mode);
            });
        });
    }

    showGEDetailedView(group, mode) {
        if (!this.currentGeneEffect) return;

        const data = this.currentGeneEffect.data;
        const gene = this.currentGeneEffect.gene;

        // Helper to calculate stats
        const calcStats = (values) => {
            const n = values.length;
            if (n === 0) return { n: 0, mean: NaN, sd: NaN };
            const mean = values.reduce((a, b) => a + b, 0) / n;
            const sd = n > 1 ? Math.sqrt(values.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / (n - 1)) : 0;
            return { n, mean, sd };
        };

        // Helper to format stats annotation
        const formatStats = (stats, label) => {
            if (stats.n === 0) return `${label}: n=0`;
            return `${label}: n=${stats.n}, GE=${stats.mean.toFixed(2)}, SD=${stats.sd.toFixed(2)}`;
        };

        if (mode === 'hotspot') {
            // Filter by hotspot mutation - show 3 levels (0, 1, 2)
            const mutData = this.mutations?.geneData?.[group]?.mutations || {};
            const data0 = data.filter(d => (mutData[d.cellLineId] || 0) === 0);
            const data1 = data.filter(d => (mutData[d.cellLineId] || 0) === 1);
            const data2 = data.filter(d => (mutData[d.cellLineId] || 0) >= 2);

            const effects0 = data0.map(d => d.geneEffect);
            const effects1 = data1.map(d => d.geneEffect);
            const effects2 = data2.map(d => d.geneEffect);
            const stats0 = calcStats(effects0);
            const stats1 = calcStats(effects1);
            const stats2 = calcStats(effects2);

            // Calculate p-value (WT vs 1+2 combined)
            const mutAllEffects = [...effects1, ...effects2];
            let pValue = NaN;
            if (stats0.n >= 3 && mutAllEffects.length >= 3) {
                const tTest = this.welchTTest(effects0, mutAllEffects);
                pValue = tTest.p;
            }

            // Create three traces for 0, 1, 2 mutations
            const traces = [
                {
                    type: 'box',
                    name: `0 (WT), n=${stats0.n}`,
                    y: effects0,
                    text: data0.map(d => d.cellLineName),
                    customdata: data0.map(d => d.lineage || 'Unknown'),
                    boxpoints: 'all',
                    jitter: 0.3,
                    pointpos: 0,
                    marker: { color: '#2563eb', size: 5 },
                    line: { color: '#1e40af', width: 2 },
                    fillcolor: 'rgba(37, 99, 235, 0.4)',
                    hovertemplate: '<b>%{text}</b><br>%{customdata}<br>Gene Effect: %{y:.3f}<extra>0 (WT)</extra>'
                }
            ];

            // Add 1 mutation trace if data exists
            if (data1.length > 0) {
                traces.push({
                    type: 'box',
                    name: `1 mut, n=${stats1.n}`,
                    y: effects1,
                    text: data1.map(d => d.cellLineName),
                    customdata: data1.map(d => d.lineage || 'Unknown'),
                    boxpoints: 'all',
                    jitter: 0.3,
                    pointpos: 0,
                    marker: { color: '#f97316', size: 5 },
                    line: { color: '#c2410c', width: 2 },
                    fillcolor: 'rgba(249, 115, 22, 0.4)',
                    hovertemplate: '<b>%{text}</b><br>%{customdata}<br>Gene Effect: %{y:.3f}<extra>1 mutation</extra>'
                });
            }

            // Add 2 mutations trace if data exists
            if (data2.length > 0) {
                traces.push({
                    type: 'box',
                    name: `2 mut, n=${stats2.n}`,
                    y: effects2,
                    text: data2.map(d => d.cellLineName),
                    customdata: data2.map(d => d.lineage || 'Unknown'),
                    boxpoints: 'all',
                    jitter: 0.3,
                    pointpos: 0,
                    marker: { color: '#dc2626', size: 5 },
                    line: { color: '#991b1b', width: 2 },
                    fillcolor: 'rgba(220, 38, 38, 0.4)',
                    hovertemplate: '<b>%{text}</b><br>%{customdata}<br>Gene Effect: %{y:.3f}<extra>2 mutations</extra>'
                });
            }

            // Build stats as separate annotations for each row to ensure they display properly
            const statsAnnotations = [];
            const pStr = !isNaN(pValue) ? (pValue < 0.001 ? pValue.toExponential(2) : pValue.toFixed(4)) : null;

            // Row 1: 0 (WT) stats
            statsAnnotations.push(`0 (WT): n=${stats0.n}, GE=${stats0.mean.toFixed(2)}, SD=${stats0.sd.toFixed(2)}`);
            // Row 2: 1 mut stats
            if (stats1.n > 0) {
                statsAnnotations.push(`1 mut: n=${stats1.n}, GE=${stats1.mean.toFixed(2)}${stats1.n > 1 ? `, SD=${stats1.sd.toFixed(2)}` : ''}`);
            }
            // Row 3: 2 mut stats
            if (stats2.n > 0) {
                statsAnnotations.push(`2 mut: n=${stats2.n}, GE=${stats2.mean.toFixed(2)}${stats2.n > 1 ? `, SD=${stats2.sd.toFixed(2)}` : ''}`);
            }
            // Row 4: p-value
            if (pStr) {
                statsAnnotations.push(`p-value (0 vs 1+2): ${pStr}`);
            }
            const statsText = statsAnnotations.join('<br>');

            // Apply width ratio from slider (default 1.0 = full width)
            const widthRatio = this.geChartWidthRatio || 1.0;
            const plotId = this.currentGEView === 'tissue' ? 'geneEffectPlot' : 'geneEffectHotspotPlot';
            const container = document.getElementById(plotId);
            const containerWidth = container ? container.offsetWidth || 500 : 500;
            const chartWidth = Math.round(containerWidth * widthRatio);

            const layout = {
                title: { text: `${gene} gene effect by ${group} mutation status`, font: { size: 14 } },
                yaxis: { title: 'Gene Effect', zeroline: true, zerolinecolor: '#374151' },
                showlegend: false,
                height: 520,
                width: chartWidth,
                margin: { t: 50, b: 170, l: 60, r: 30 },
                paper_bgcolor: 'white',
                plot_bgcolor: 'white',
                annotations: [{
                    x: 0,
                    y: -0.38,
                    xref: 'paper',
                    yref: 'paper',
                    xanchor: 'left',
                    text: statsText,
                    showarrow: false,
                    font: { size: 10, family: 'monospace' },
                    align: 'left'
                }]
            };

            // Mark that we're in detailed view
            this.geDetailedView = { mode: 'hotspot', group };

            Plotly.newPlot(plotId, traces, layout, { responsive: true });
            this.updateShowAllButton();
            return;
        }

        // For tissue mode, show comparison with all cells
        const filteredData = data.filter(d => (d.lineage || 'Unknown') === group);
        const allEffects = data.map(d => d.geneEffect);
        const groupEffects = filteredData.map(d => d.geneEffect);
        const allStats = calcStats(allEffects);
        const groupStats = calcStats(groupEffects);

        // Calculate p-value comparing group to all cells
        let pValue = NaN;
        if (groupStats.n >= 3 && allStats.n >= 3) {
            const tTest = this.welchTTest(groupEffects, allEffects);
            pValue = tTest.p;
        }

        // Create traces for All Cells and the selected tissue (include cancer type in hover)
        const traces = [
            {
                type: 'box',
                name: `All cells, n=${allStats.n}`,
                y: allEffects,
                text: data.map(d => d.cellLineName),
                customdata: data.map(d => d.lineage || 'Unknown'),
                boxpoints: 'all',
                jitter: 0.3,
                pointpos: 0,
                marker: { color: '#6b7280', size: 4 },
                line: { color: '#374151', width: 2 },
                fillcolor: 'rgba(107, 114, 128, 0.3)',
                hovertemplate: '<b>%{text}</b><br>%{customdata}<br>Gene Effect: %{y:.3f}<extra>All cells</extra>'
            },
            {
                type: 'box',
                name: `${group}, n=${groupStats.n}`,
                y: groupEffects,
                text: filteredData.map(d => d.cellLineName),
                boxpoints: 'all',
                jitter: 0.3,
                pointpos: 0,
                marker: {
                    color: groupEffects.map(v => v < -0.5 ? '#dc2626' : v > 0 ? '#16a34a' : '#6b7280'),
                    size: 6
                },
                line: { color: '#5a9f4a', width: 2 },
                fillcolor: 'rgba(90, 159, 74, 0.4)',
                hovertemplate: '<b>%{text}</b><br>Gene Effect: %{y:.3f}<extra>' + group + '</extra>'
            }
        ];

        // Build stats as separate rows using <br> for proper line breaks
        const statsAnnotations = [];
        statsAnnotations.push(`All cells: n=${allStats.n}, GE=${allStats.mean.toFixed(2)}, SD=${allStats.sd.toFixed(2)}`);
        statsAnnotations.push(`${group}: n=${groupStats.n}, GE=${groupStats.mean.toFixed(2)}, SD=${groupStats.sd.toFixed(2)}`);
        if (!isNaN(pValue)) {
            statsAnnotations.push(`p-value: ${pValue < 0.001 ? pValue.toExponential(2) : pValue.toFixed(4)}`);
        }
        const statsText = statsAnnotations.join('<br>');

        // Apply width ratio from slider (default 1.0 = full width)
        const widthRatio = this.geChartWidthRatio || 1.0;
        const plotId = this.currentGEView === 'tissue' ? 'geneEffectPlot' : 'geneEffectHotspotPlot';
        const container = document.getElementById(plotId);
        const containerWidth = container ? container.offsetWidth || 500 : 500;
        const chartWidth = Math.round(containerWidth * widthRatio);

        const layout = {
            title: { text: `${gene} gene effect in ${group}`, font: { size: 14 } },
            yaxis: { title: 'Gene Effect', zeroline: true, zerolinecolor: '#374151', zerolinewidth: 2 },
            showlegend: false,
            height: 500,
            width: chartWidth,
            margin: { t: 50, b: 140, l: 60, r: 30 },
            paper_bgcolor: 'white',
            plot_bgcolor: 'white',
            annotations: [{
                x: 0,
                y: -0.30,
                xref: 'paper',
                yref: 'paper',
                xanchor: 'left',
                text: statsText,
                showarrow: false,
                font: { size: 10, family: 'monospace' },
                align: 'left'
            }]
        };

        // Mark that we're in detailed view
        this.geDetailedView = { mode: 'tissue', group };

        Plotly.newPlot(plotId, traces, layout, { responsive: true });
        this.updateShowAllButton();
    }

    updateShowAllButton() {
        const btn = document.getElementById('geShowAllBtn');
        const ratioControl = document.getElementById('geAspectRatioControl');
        const isInspect = this.geneEffectViewMode === 'mutation';
        if (btn) {
            btn.style.display = this.geDetailedView ? 'inline-block' : 'none';
        }
        if (ratioControl) {
            // Show controls for detailed views AND mutation inspect views
            ratioControl.style.display = (this.geDetailedView || isInspect) ? 'flex' : 'none';
            // Height control only relevant for inspect plot (not tissue/hotspot overview)
            const heightControl = document.getElementById('geHeightControl');
            if (heightControl) {
                heightControl.style.display = isInspect ? 'inline-flex' : 'none';
            }
        }
    }

    showAllGeneEffect() {
        this.geDetailedView = null;
        this.updateShowAllButton();
        if (this.currentGEView === 'tissue') {
            this.renderGeneEffectByTissue();
        } else {
            this.renderGeneEffectByHotspot();
        }
    }

    sortGETable(sortKey, sortType) {
        if (!this.currentGEStats) return;

        // Toggle sort direction
        if (this.geTableSortKey === sortKey) {
            this.geTableSortDir = this.geTableSortDir === 'asc' ? 'desc' : 'asc';
        } else {
            this.geTableSortKey = sortKey;
            this.geTableSortDir = 'desc'; // Default to descending for new column
        }

        const dir = this.geTableSortDir === 'asc' ? 1 : -1;

        this.currentGEStats.sort((a, b) => {
            let aVal = a[sortKey];
            let bVal = b[sortKey];

            if (sortType === 'string') {
                return dir * aVal.localeCompare(bVal);
            } else {
                return dir * (aVal - bVal);
            }
        });

        this.renderGETableBody(this.currentGEStats, this.currentGETableMode);

        // Update header to show sort direction
        const thead = document.getElementById('geTableHead');
        thead.querySelectorAll('th').forEach(th => {
            const baseText = th.textContent.replace(/ [↑↓↕]$/, '');
            if (th.dataset.sort === sortKey) {
                th.textContent = baseText + (this.geTableSortDir === 'asc' ? ' ↑' : ' ↓');
            } else {
                th.textContent = baseText + ' ↕';
            }
        });
    }

    async downloadGeneEffectChartPNG() {
        // Mutation inspect mode handled by downloadGeneEffectPNG
        if (this.geneEffectViewMode === 'mutation') return;
        if (!this.currentGeneEffect) return;
        const plotId = this.currentGEView === 'tissue' ? 'geneEffectPlot' : 'geneEffectHotspotPlot';
        const plotEl = document.getElementById(plotId);
        const chartWidth = plotEl?._fullLayout?.width || 800;
        const chartHeight = Math.max(plotEl?.scrollHeight || 0, plotEl?._fullLayout?.height || 0, this.geDetailedView ? 650 : 550) + 40;
        const filename = `gene_effect_${this.currentGeneEffect.gene}_by_${this.currentGEView}`;
        const meta = this._buildExportMetadata('gene_effect', {
            gene: this.currentGeneEffect.gene, view: this.currentGEView,
            textSettings: this._capturePlotTextSettings(this.currentGEView === 'tissue' ? 'geneEffectPlot' : 'geneEffectHotspotPlot')
        });
        const metaJson = JSON.stringify(meta);

        const svgDataUrl = await Plotly.toImage(plotEl, { format: 'svg', width: chartWidth, height: chartHeight });
        let svgStr;
        if (svgDataUrl.indexOf('base64,') > -1) svgStr = atob(svgDataUrl.split('base64,')[1]);
        else svgStr = decodeURIComponent(svgDataUrl.split(',').slice(1).join(','));

        const svgBlob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
        const svgUrl = URL.createObjectURL(svgBlob);
        const img = new Image();
        img.onload = async () => {
            const scale = 4;
            const canvas = document.createElement('canvas');
            canvas.width = img.naturalWidth * scale;
            canvas.height = img.naturalHeight * scale;
            const ctx = canvas.getContext('2d');
            ctx.scale(scale, scale);
            ctx.fillStyle = 'white';
            ctx.fillRect(0, 0, img.naturalWidth, img.naturalHeight);
            ctx.drawImage(img, 0, 0);
            URL.revokeObjectURL(svgUrl);
            const pngDataUrl = canvas.toDataURL('image/png');
            const pngResp = await fetch(pngDataUrl);
            const pngBuf = await pngResp.arrayBuffer();
            const pngWithMeta = this._addPngTextChunk(pngBuf, 'correlate-meta', metaJson);
            const blob = new Blob([pngWithMeta], { type: 'image/png' });
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = `${filename}.png`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(a.href);
        };
        img.src = svgUrl;
    }

    downloadGeneEffectChartSVG() {
        // Mutation inspect mode handled by downloadGeneEffectSVG
        if (this.geneEffectViewMode === 'mutation') return;
        if (!this.currentGeneEffect) return;
        const plotId = this.currentGEView === 'tissue' ? 'geneEffectPlot' : 'geneEffectHotspotPlot';
        const plotEl = document.getElementById(plotId);
        const chartWidth = plotEl?._fullLayout?.width || 800;
        const chartHeight = Math.max(plotEl?.scrollHeight || 0, plotEl?._fullLayout?.height || 0, this.geDetailedView ? 650 : 550) + 40;
        const filename = `gene_effect_${this.currentGeneEffect.gene}_by_${this.currentGEView}`;

        // Use toImage + post-process to expand viewBox and remove clipPaths that crop content
        Plotly.toImage(plotEl, {
            format: 'svg',
            width: chartWidth,
            height: chartHeight
        }).then(async svgDataUrl => {
            let svgString;
            if (svgDataUrl.indexOf('base64,') > -1) {
                svgString = atob(svgDataUrl.split('base64,')[1]);
            } else {
                svgString = decodeURIComponent(svgDataUrl.split(',').slice(1).join(','));
            }

            const parser = new DOMParser();
            const svgDoc = parser.parseFromString(svgString, 'image/svg+xml');
            const svgEl = svgDoc.documentElement;

            // Measure true bounding box by temporarily removing clipPaths
            const measurer = document.createElement('div');
            measurer.style.cssText = 'position:absolute;left:-99999px;top:-99999px;';
            document.body.appendChild(measurer);
            const measureSvg = svgEl.cloneNode(true);
            measureSvg.style.overflow = 'visible';
            measureSvg.querySelectorAll('[clip-path]').forEach(el => el.removeAttribute('clip-path'));
            measurer.appendChild(measureSvg);

            try {
                const bbox = measureSvg.getBBox();
                const pad = 10;
                const newX = Math.min(0, bbox.x - pad);
                const newY = Math.min(0, bbox.y - pad);
                const newW = Math.max(chartWidth, bbox.x + bbox.width + pad) - newX;
                const newH = Math.max(chartHeight, bbox.y + bbox.height + pad) - newY;
                svgEl.setAttribute('viewBox', `${newX} ${newY} ${newW} ${newH}`);
                svgEl.setAttribute('width', newW);
                svgEl.setAttribute('height', newH);
            } catch (e) {
                // fallback: just use original dimensions
            }
            document.body.removeChild(measurer);

            // Remove clipPaths on the plot area that crop y-axis labels
            svgEl.querySelectorAll('clipPath').forEach(cp => {
                const rect = cp.querySelector('rect');
                if (rect && parseFloat(rect.getAttribute('x') || 0) === 0 && parseFloat(rect.getAttribute('y') || 0) === 0) {
                    rect.setAttribute('x', -500);
                    rect.setAttribute('width', parseFloat(rect.getAttribute('width') || 0) + 600);
                    rect.setAttribute('y', -50);
                    rect.setAttribute('height', parseFloat(rect.getAttribute('height') || 0) + 100);
                }
            });

            let finalSvg = new XMLSerializer().serializeToString(svgEl);
            const meta = this._buildExportMetadata('gene_effect', {
                gene: this.currentGeneEffect?.gene, view: this.currentGEView,
                textSettings: this._capturePlotTextSettings(this.currentGEView === 'tissue' ? 'geneEffectPlot' : 'geneEffectHotspotPlot')
            });
            finalSvg = finalSvg.replace('</svg>', `<metadata><correlate-meta>${JSON.stringify(meta)}</correlate-meta></metadata></svg>`);
            finalSvg = await this._finalizeSvgForExport(finalSvg);
            const blob = new Blob([finalSvg], { type: 'image/svg+xml' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename + '.svg';
            a.click();
            URL.revokeObjectURL(url);
        });
    }

    downloadGETableCSV() {
        // Mutation inspect mode: export by group (WT / 1 mutation / 2 mutations)
        if (this.geneEffectViewMode === 'mutation' && this.currentGeneEffectData && this.currentGeneEffectGene) {
            return this.downloadGeneEffectCSV();
        }

        if (!this.currentGeneEffect || !this.currentGEStats) return;

        const gene = this.currentGeneEffect.gene;
        let csv = '';

        if (this.currentGEView === 'tissue') {
            csv = 'Cancer_Type,N,Mean_GE,SD,p_value\n';
            this.currentGEStats.forEach(s => {
                csv += `"${s.group}",${s.n},${s.mean.toFixed(2)},${s.sd.toFixed(2)},${s.pValue.toFixed(6)}\n`;
            });
        } else {
            csv = 'Hotspot,N_Mutant,N_WT,Mean_WT,Mean_Mutant,Delta_GE,SD_Mutant,p_value\n';
            this.currentGEStats.forEach(s => {
                csv += `"${s.group}",${s.nMut},${s.nWT},${s.wtMean.toFixed(2)},${s.mutMean.toFixed(2)},${s.diff.toFixed(2)},${s.mutSD.toFixed(2)},${s.pValue.toFixed(6)}\n`;
            });
        }

        this.downloadFile(csv, `gene_effect_${gene}_by_${this.currentGEView}.csv`, 'text/csv');
    }

    downloadGECellLineCSV() {
        // Mutation inspect mode: export per cell line
        if (this.geneEffectViewMode === 'mutation' && this.currentGeneEffectData && this.currentGeneEffectGene) {
            const mr = this.mutationResults;
            const gene = this.currentGeneEffectGene;
            let csv = 'Cell_Line_ID,Cell_Line_Name,Cancer_Type,Cancer_Subtype,Gene_Effect,Mutation_Level\n';
            this.currentGeneEffectData.forEach(d => {
                const subtype = this.getCellLineSublineage?.(d.cellLine) || '';
                csv += `"${d.cellLine}","${d.cellName}","${d.lineage}","${subtype}",${d.ge.toFixed(2)},${d.mutLevel}\n`;
            });
            this.downloadFile(csv, `gene_effect_${gene}_${mr.hotspotGene}_by_cell_line.csv`, 'text/csv');
            return;
        }

        if (!this.currentGeneEffect) return;

        const gene = this.currentGeneEffect.gene;
        const data = this.currentGeneEffect.data;

        // Get sublineage info
        let csv = 'Cell_Line_ID,Cell_Line_Name,Cancer_Type,Cancer_Subtype,Gene_Effect\n';

        data.forEach(d => {
            const subtype = this.getCellLineSublineage?.(d.cellLineId) || '';
            csv += `"${d.cellLineId}","${d.cellLineName}","${d.lineage}","${subtype}",${d.geneEffect.toFixed(2)}\n`;
        });

        this.downloadFile(csv, `gene_effect_${gene}_by_cell_line.csv`, 'text/csv');
    }

    getCellLineSublineage(cellLineId) {
        if (!this.cellLineMetadata) return '';
        return this.cellLineMetadata.primaryDisease?.[cellLineId] || '';
    }

    // ===== Copy Column / Copy Genes functionality (#9/#10) =====

    copyColumnToClipboard(tableEl, colIndex) {
        const tbody = tableEl.querySelector('tbody');
        const rows = Array.from(tbody.querySelectorAll('tr')).filter(r => r.style.display !== 'none');
        const values = rows.map(r => {
            const cell = r.children[colIndex];
            return cell ? cell.textContent.trim() : '';
        }).filter(v => v !== '' && v !== '-');

        const text = values.join('\n');
        navigator.clipboard.writeText(text).then(() => {
            this.showCopyNotification(`Copied ${values.length} values`);
        });
    }

    copyGeneColumn(tableId) {
        const table = document.getElementById(tableId);
        if (!table) return;
        const tbody = table.querySelector('tbody');
        const rows = Array.from(tbody.querySelectorAll('tr')).filter(r => r.style.display !== 'none');

        // Find gene column indices
        const headers = Array.from(table.querySelectorAll('thead th'));
        const geneColIndices = [];
        headers.forEach((th, idx) => {
            const sort = th.dataset.sort || th.dataset.col || '';
            if (sort === 'gene' || sort === 'gene1' || sort === 'gene2') {
                geneColIndices.push(idx);
            }
        });
        // Fallback: mutation table has gene in col 1 (col 0 is inspect link)
        if (geneColIndices.length === 0) {
            if (tableId === 'mutationTable') geneColIndices.push(1);
            else geneColIndices.push(0);
        }

        // Collect genes from all gene columns, deduplicate
        const geneSet = new Set();
        rows.forEach(r => {
            geneColIndices.forEach(colIdx => {
                const cell = r.children[colIdx];
                if (cell) {
                    const gene = cell.textContent.trim().replace(/\*$/, '');
                    if (gene) geneSet.add(gene);
                }
            });
        });

        const genes = [...geneSet];
        navigator.clipboard.writeText(genes.join('\n')).then(() => {
            this.showCopyNotification(`Copied ${genes.length} genes`);
        });
    }

    showCopyNotification(message) {
        const notification = document.createElement('div');
        notification.textContent = message;
        notification.style.cssText = 'position: fixed; top: 20px; right: 20px; background: #5a9f4a; color: white; padding: 8px 16px; border-radius: 6px; font-size: 13px; z-index: 10000; transition: opacity 0.3s; box-shadow: 0 2px 8px rgba(0,0,0,0.15);';
        document.body.appendChild(notification);
        setTimeout(() => {
            notification.style.opacity = '0';
            setTimeout(() => notification.remove(), 300);
        }, 1500);
    }

    // ===== Table Column Filters (#14) =====

    toggleTableFilters(tableId) {
        const table = document.getElementById(tableId);
        if (!table) return;

        let filterRow = table.querySelector('.filter-row');
        if (filterRow) {
            const isHidden = filterRow.style.display === 'none';
            filterRow.style.display = isHidden ? '' : 'none';
            if (!isHidden) {
                // Clear all filters when hiding
                filterRow.querySelectorAll('input').forEach(input => { input.value = ''; });
                this.applyColumnFilters(tableId);
            }
            return;
        }

        // Create filter row
        const thead = table.querySelector('thead');
        const headerRow = thead.querySelector('tr');
        const headers = headerRow.querySelectorAll('th');

        filterRow = document.createElement('tr');
        filterRow.className = 'filter-row';
        filterRow.style.background = '#f9fafb';

        const numericSortKeys = [
            'correlation', 'slope', 'n', 'cluster',
            'meanEffect', 'sdEffect', 'meanEffectFiltered', 'sdEffectFiltered',
            'lfc', 'fdr', 'hasCorrelation',
            'r', 'p'
        ];
        const numericColKeys = [
            'n_wt', 'mean_wt', 'n_mut', 'mean_mut', 'diff_mut', 'p_mut',
            'n_2', 'mean_2', 'diff_2', 'p_2', 'diff_2v1', 'p_2v1',
            'n_fused', 'mean_fused', 'diff_fused', 'p_fused'
        ];

        headers.forEach((th, idx) => {
            const td = document.createElement('td');
            td.style.padding = '2px 4px';

            const sortKey = th.dataset.sort || th.dataset.col || '';
            const isNumeric = numericSortKeys.includes(sortKey) || numericColKeys.includes(sortKey);

            if (sortKey) {
                const input = document.createElement('input');
                input.type = 'text';
                input.style.cssText = 'width: 100%; font-size: 10px; padding: 2px 4px; border: 1px solid #d1d5db; border-radius: 3px; box-sizing: border-box;';
                input.placeholder = isNumeric ? '>0.5 or <-1' : 'filter';
                input.dataset.colIndex = idx;
                input.dataset.isNumeric = isNumeric;
                input.addEventListener('input', () => this.applyColumnFilters(tableId));
                td.appendChild(input);
            }

            filterRow.appendChild(td);
        });

        thead.appendChild(filterRow);
    }

    applyColumnFilters(tableId) {
        const table = document.getElementById(tableId);
        if (!table) return;

        const filterRow = table.querySelector('.filter-row');
        if (!filterRow) return;

        const filters = [];
        filterRow.querySelectorAll('input').forEach(input => {
            const val = input.value.trim();
            if (val) {
                filters.push({
                    colIndex: parseInt(input.dataset.colIndex),
                    isNumeric: input.dataset.isNumeric === 'true',
                    condition: val
                });
            }
        });

        const tbody = table.querySelector('tbody');
        const rows = tbody.querySelectorAll('tr');

        rows.forEach(row => {
            let show = true;

            filters.forEach(filter => {
                if (!show) return;
                const cell = row.children[filter.colIndex];
                if (!cell) return;
                const cellText = cell.textContent.trim();

                if (filter.isNumeric) {
                    // Replace proper minus sign with hyphen for parsing
                    const cellNum = parseFloat(cellText.replace('−', '-').replace(/^[<>]/g, ''));
                    if (isNaN(cellNum)) { show = false; return; }

                    const match = filter.condition.match(/^([><=!]+)?\s*([\d.eE\-+]+)$/);
                    if (match) {
                        const op = match[1] || '>=';
                        const threshold = parseFloat(match[2]);
                        if (isNaN(threshold)) return;

                        if (op === '>') show = cellNum > threshold;
                        else if (op === '<') show = cellNum < threshold;
                        else if (op === '>=') show = cellNum >= threshold;
                        else if (op === '<=') show = cellNum <= threshold;
                        else if (op === '=' || op === '==') show = Math.abs(cellNum - threshold) < 1e-10;
                        else if (op === '!=') show = Math.abs(cellNum - threshold) >= 1e-10;
                    }
                } else {
                    show = cellText.toLowerCase().includes(filter.condition.toLowerCase());
                }
            });

            row.style.display = show ? '' : 'none';
        });
    }

    // ===== Gene Info Tooltips (#15) =====

    async fetchGeneInfo(gene) {
        // Cache results
        if (!this.geneInfoCache) this.geneInfoCache = {};
        if (this.geneInfoCache[gene]) return this.geneInfoCache[gene];

        try {
            const res = await fetch(`https://mygene.info/v3/query?q=symbol:${gene}&species=human&fields=symbol,name,summary&size=1`);
            const data = await res.json();
            if (data.hits && data.hits.length > 0) {
                const hit = data.hits[0];
                const info = {
                    symbol: hit.symbol || gene,
                    name: hit.name || '',
                    summary: hit.summary || ''
                };
                this.geneInfoCache[gene] = info;
                return info;
            }
        } catch (e) {
            console.warn('Gene info fetch failed:', e);
        }
        return null;
    }

    showGeneTooltip(event, gene) {
        // Remove existing tooltip
        this.hideGeneTooltip();

        const tooltip = document.createElement('div');
        tooltip.id = 'geneTooltip';
        tooltip.style.cssText = 'position: fixed; z-index: 10001; background: white; border: 1px solid #d1d5db; border-radius: 8px; padding: 10px 14px; max-width: 350px; box-shadow: 0 4px 12px rgba(0,0,0,0.15); font-size: 11px; line-height: 1.4;';
        tooltip.innerHTML = `<div style="color: #6b7280;">Loading ${gene} info...</div>`;

        // Position near the mouse
        const x = Math.min(event.clientX + 10, window.innerWidth - 370);
        const y = Math.min(event.clientY + 10, window.innerHeight - 200);
        tooltip.style.left = x + 'px';
        tooltip.style.top = y + 'px';

        document.body.appendChild(tooltip);

        this.fetchGeneInfo(gene).then(info => {
            const el = document.getElementById('geneTooltip');
            if (!el) return;

            if (!info) {
                el.innerHTML = `<b>${gene}</b><br><span style="color: #999;">No info available</span>`;
                return;
            }

            let html = `<div style="margin-bottom: 4px;"><b style="color: #5a9f4a; font-size: 13px;">${info.symbol}</b> <span style="color: #374151;">${info.name}</span></div>`;

            if (info.summary) {
                const shortSummary = info.summary.length > 200 ? info.summary.substring(0, 200) + '...' : info.summary;
                html += `<div style="color: #4b5563;">${shortSummary}</div>`;
            }

            el.innerHTML = html;

            // Reposition if needed
            const rect = el.getBoundingClientRect();
            if (rect.bottom > window.innerHeight) {
                el.style.top = (window.innerHeight - rect.height - 10) + 'px';
            }
        });
    }

    hideGeneTooltip() {
        const existing = document.getElementById('geneTooltip');
        if (existing) existing.remove();
    }

    // ===== Inline Compare by Tissue/Hotspot (in inspect modal) =====

    showInlineCompareByTissue() {
        if (!this.mutationResults || !this.currentGeneEffectGene) return;
        this._inlineSortCol = null;
        this._inlineSortAsc = true;
        const mr = this.mutationResults;
        const gene = this.currentGeneEffectGene;
        const hotspotGene = mr.hotspotGene;
        const isTranslocation = mr.isTranslocation;
        const mutationData = isTranslocation
            ? this.translocations?.geneData?.[hotspotGene]
            : this.mutations?.geneData?.[hotspotGene];
        if (!mutationData) return;

        const geneIdx = this.geneIndex.get(gene.toUpperCase());
        if (geneIdx === undefined) return;

        const cellLines = this.metadata.cellLines;
        const inspectHotspot = document.getElementById('geHotspotFilter')?.value || '';
        const inspectFusion = document.getElementById('geFusionFilter')?.value || '';
        const inspFusionData = inspectFusion ? this.translocations?.geneData?.[inspectFusion]?.translocations : null;
        const inspectSubtype = document.getElementById('geSubtypeFilter')?.value || '';

        // Determine if we should group by subtype (when a lineage filter is active)
        const activeTissueFilter = document.getElementById('geTissueFilter')?.value || '';
        const activeLineage = activeTissueFilter || mr.lineageFilter;
        const groupBySubtype = !!activeLineage && !!this.cellLineMetadata?.primaryDisease;

        // Gather cell lines respecting all filters
        const groupMap = {};
        const allWT = [], allMut = [];
        cellLines.forEach((cellLine, idx) => {
            // Apply tissue/lineage filters
            if (activeTissueFilter) {
                if ((this.cellLineMetadata?.lineage?.[cellLine] || '') !== activeTissueFilter) return;
            } else {
                if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) return;
                if (mr.excludedTissues && mr.excludedTissues.size > 0) {
                    const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                    if (lineage && mr.excludedTissues.has(lineage)) return;
                }
            }
            if (!activeTissueFilter && mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) return;
            if (inspectSubtype && this.cellLineMetadata?.primaryDisease?.[cellLine] !== inspectSubtype) return;
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
            if (inspectHotspot) {
                const inspHotData = this.mutations.geneData[inspectHotspot];
                if (inspHotData) {
                    const inspMutLevel = inspHotData.mutations[cellLine] || 0;
                    if (inspMutLevel === 0) return;
                }
            }
            if (inspectFusion && inspFusionData) {
                if ((inspFusionData[cellLine] || 0) < 1) return;
            }

            const ge = this.geneEffects[geneIdx * this.nCellLines + idx];
            if (isNaN(ge)) return;

            const groupKey = groupBySubtype
                ? (this.cellLineMetadata.primaryDisease[cellLine] || 'Unknown')
                : (this.cellLineMetadata?.lineage?.[cellLine] || 'Unknown');
            const mutLevel = isTranslocation
                ? (mutationData.translocations[cellLine] || 0)
                : (mutationData.mutations[cellLine] || 0);
            if (!groupMap[groupKey]) groupMap[groupKey] = { wt: [], mut: [] };
            if (mutLevel === 0) { groupMap[groupKey].wt.push(ge); allWT.push(ge); }
            else { groupMap[groupKey].mut.push(ge); allMut.push(ge); }
        });

        // Build rows
        const rows = [];
        const mutLabel = isTranslocation ? 'Fused' : 'Mut';
        if (allWT.length > 0 && allMut.length > 0) {
            const meanWT = allWT.reduce((a, b) => a + b, 0) / allWT.length;
            const meanMut = allMut.reduce((a, b) => a + b, 0) / allMut.length;
            rows.push({ label: 'All', nWT: allWT.length, meanWT, nMut: allMut.length, meanMut, delta: meanMut - meanWT, isAll: true, tissue: '' });
        }
        Object.entries(groupMap).forEach(([group, groups]) => {
            if (groups.wt.length > 0 && groups.mut.length > 0) {
                const meanWT = groups.wt.reduce((a, b) => a + b, 0) / groups.wt.length;
                const meanMut = groups.mut.reduce((a, b) => a + b, 0) / groups.mut.length;
                rows.push({ label: group, nWT: groups.wt.length, meanWT, nMut: groups.mut.length, meanMut, delta: meanMut - meanWT, tissue: group });
            }
        });

        const allRow = rows.filter(r => r.isAll);
        const otherRows = rows.filter(r => !r.isAll).sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta));

        const groupLabel = groupBySubtype ? 'Subtype' : 'Tissue';
        this._inlineCompareData = { title: `${gene} — Δ GE by ${groupLabel} (${hotspotGene} WT vs ${mutLabel})`, headers: [groupLabel, 'N(WT)', 'GE(WT)', `N(${mutLabel})`, `GE(${mutLabel})`, 'Δ GE'], refRows: allRow, sortableRows: otherRows, mode: 'tissue' };
        this._renderInlineCompareTable();
    }

    showInlineCompareByHotspot() {
        if (!this.mutationResults || !this.currentGeneEffectGene) return;
        this._inlineSortCol = null;
        this._inlineSortAsc = true;
        const mr = this.mutationResults;
        const gene = this.currentGeneEffectGene;
        const mainHotspot = mr.hotspotGene;
        const isTranslocation = mr.isTranslocation;
        const mainMutData = isTranslocation
            ? this.translocations?.geneData?.[mainHotspot]
            : this.mutations?.geneData?.[mainHotspot];
        if (!mainMutData) return;

        const geneIdx = this.geneIndex.get(gene.toUpperCase());
        if (geneIdx === undefined) return;

        const cellLines = this.metadata.cellLines;
        const tissueFilter = document.getElementById('geTissueFilter')?.value || '';
        const inspectSubtype = document.getElementById('geSubtypeFilter')?.value || '';

        const baseCells = [];
        cellLines.forEach((cellLine, idx) => {
            if (tissueFilter) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine] || '';
                if (lineage !== tissueFilter) return;
            } else {
                if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) return;
                if (mr.excludedTissues && mr.excludedTissues.size > 0) {
                    const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                    if (lineage && mr.excludedTissues.has(lineage)) return;
                }
            }
            if (!tissueFilter && mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) return;
            if (inspectSubtype && this.cellLineMetadata?.primaryDisease?.[cellLine] !== inspectSubtype) return;
            if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
                const addMutData = this.mutations?.geneData?.[mr.additionalHotspot];
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
            const mainLevel = isTranslocation
                ? (mainMutData.translocations[cellLine] || 0)
                : (mainMutData.mutations[cellLine] || 0);
            baseCells.push({ cellLine, idx, ge, mainMut: mainLevel });
        });

        const rows = [];
        const mutLabel = isTranslocation ? 'Fused' : 'Mut';
        const noneWT = baseCells.filter(c => c.mainMut === 0).map(c => c.ge);
        const noneMut = baseCells.filter(c => c.mainMut >= 1).map(c => c.ge);
        if (noneWT.length > 0 && noneMut.length > 0) {
            const meanWT = noneWT.reduce((a, b) => a + b, 0) / noneWT.length;
            const meanMut = noneMut.reduce((a, b) => a + b, 0) / noneMut.length;
            rows.push({ label: 'None', nWT: noneWT.length, meanWT, nMut: noneMut.length, meanMut, delta: meanMut - meanWT, isRef: true, hotspot: '' });
        }

        // Iterate over hotspot mutation genes only (fusions handled by showInlineCompareByTranslocation)
        this.mutations?.genes?.forEach(hGene => {
            if (hGene === mainHotspot && !isTranslocation) return;
            const hMutData = this.mutations.geneData[hGene];
            if (!hMutData) return;
            const filtered = baseCells.filter(c => (hMutData.mutations[c.cellLine] || 0) >= 1);
            const wtGE = filtered.filter(c => c.mainMut === 0).map(c => c.ge);
            const mutGE = filtered.filter(c => c.mainMut >= 1).map(c => c.ge);
            if (wtGE.length > 0 && mutGE.length > 0) {
                const meanWT = wtGE.reduce((a, b) => a + b, 0) / wtGE.length;
                const meanMut = mutGE.reduce((a, b) => a + b, 0) / mutGE.length;
                rows.push({ label: hGene, nWT: wtGE.length, meanWT, nMut: mutGE.length, meanMut, delta: meanMut - meanWT, hotspot: hGene });
            }
        });

        const refRow = rows.filter(r => r.isRef);
        const otherRows = rows.filter(r => !r.isRef).sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta));

        const typeLabel = isTranslocation ? 'Fusion' : 'Hotspot';
        this._inlineCompareData = { title: `${gene} — Δ GE by Additional ${typeLabel} (${mainHotspot} WT vs ${mutLabel})`, headers: [`${typeLabel} Filter`, 'N(WT)', 'GE(WT)', `N(${mutLabel})`, `GE(${mutLabel})`, 'Δ GE'], refRows: refRow, sortableRows: otherRows, mode: 'hotspot' };
        this._renderInlineCompareTable();
    }

    showInlineCompareByTranslocation() {
        if (!this.mutationResults || !this.currentGeneEffectGene) return;
        if (!this.translocations?.genes?.length) return;
        this._inlineSortCol = null;
        this._inlineSortAsc = true;
        const mr = this.mutationResults;
        const gene = this.currentGeneEffectGene;
        const mainHotspot = mr.hotspotGene;
        const isTranslocation = mr.isTranslocation;
        const mainMutData = isTranslocation
            ? this.translocations?.geneData?.[mainHotspot]
            : this.mutations?.geneData?.[mainHotspot];
        if (!mainMutData) return;

        const geneIdx = this.geneIndex.get(gene.toUpperCase());
        if (geneIdx === undefined) return;

        const cellLines = this.metadata.cellLines;
        const tissueFilter = document.getElementById('geTissueFilter')?.value || '';
        const inspectSubtype = document.getElementById('geSubtypeFilter')?.value || '';

        const baseCells = [];
        cellLines.forEach((cellLine, idx) => {
            if (tissueFilter) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine] || '';
                if (lineage !== tissueFilter) return;
            } else {
                if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) return;
                if (mr.excludedTissues && mr.excludedTissues.size > 0) {
                    const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                    if (lineage && mr.excludedTissues.has(lineage)) return;
                }
            }
            if (!tissueFilter && mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) return;
            if (inspectSubtype && this.cellLineMetadata?.primaryDisease?.[cellLine] !== inspectSubtype) return;
            if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
                const addMutData = this.mutations?.geneData?.[mr.additionalHotspot];
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
            const mainLevel = isTranslocation
                ? (mainMutData.translocations[cellLine] || 0)
                : (mainMutData.mutations[cellLine] || 0);
            baseCells.push({ cellLine, idx, ge, mainMut: mainLevel });
        });

        const rows = [];
        const mutLabel = isTranslocation ? 'Fused' : 'Mut';
        const noneWT = baseCells.filter(c => c.mainMut === 0).map(c => c.ge);
        const noneMut = baseCells.filter(c => c.mainMut >= 1).map(c => c.ge);
        if (noneWT.length > 0 && noneMut.length > 0) {
            const meanWT = noneWT.reduce((a, b) => a + b, 0) / noneWT.length;
            const meanMut = noneMut.reduce((a, b) => a + b, 0) / noneMut.length;
            rows.push({ label: 'None', nWT: noneWT.length, meanWT, nMut: noneMut.length, meanMut, delta: meanMut - meanWT, isRef: true, fusion: '' });
        }

        // Iterate ONLY over translocation genes — pre-filter to genes with fusions in baseCells
        const baseCellIds = new Set(baseCells.map(c => c.cellLine));
        for (const tGene of this.translocations.genes) {
            if (tGene === mainHotspot && isTranslocation) continue;
            const tData = this.translocations.geneData[tGene];
            if (!tData) continue;
            // Quick check: does this gene have any fused cells in our base set?
            let hasFused = false;
            for (const cl of baseCellIds) {
                if (tData.translocations[cl] && tData.translocations[cl] > 0) { hasFused = true; break; }
            }
            if (!hasFused) continue;
            const filtered = baseCells.filter(c => (tData.translocations[c.cellLine] || 0) >= 1);
            const wtGE = filtered.filter(c => c.mainMut === 0).map(c => c.ge);
            const mutGE = filtered.filter(c => c.mainMut >= 1).map(c => c.ge);
            if (wtGE.length > 0 && mutGE.length > 0) {
                const meanWT = wtGE.reduce((a, b) => a + b, 0) / wtGE.length;
                const meanMut = mutGE.reduce((a, b) => a + b, 0) / mutGE.length;
                rows.push({ label: tGene, nWT: wtGE.length, meanWT, nMut: mutGE.length, meanMut, delta: meanMut - meanWT, fusion: tGene });
            }
        }

        const refRow = rows.filter(r => r.isRef);
        const otherRows = rows.filter(r => !r.isRef).sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta));

        this._inlineCompareData = { title: `${gene} — Δ GE by Fusion (${mainHotspot} WT vs ${mutLabel})`, headers: ['Fusion Gene', 'N(WT)', 'GE(WT)', `N(${mutLabel})`, `GE(${mutLabel})`, 'Δ GE'], refRows: refRow, sortableRows: otherRows, mode: 'fusion' };
        this._renderInlineCompareTable();
    }

    sortInlineCompare(colIndex) {
        if (!this._inlineCompareData) return;
        if (this._inlineSortCol === colIndex) {
            this._inlineSortAsc = !this._inlineSortAsc;
        } else {
            this._inlineSortCol = colIndex;
            this._inlineSortAsc = true;
        }
        const d = this._inlineCompareData;
        const keyMap = ['label', 'nWT', 'meanWT', 'nMut', 'meanMut', 'delta'];
        const key = keyMap[colIndex] || 'delta';
        const asc = this._inlineSortAsc;
        d.sortableRows.sort((a, b) => {
            const va = key === 'label' ? a[key].toLowerCase() : a[key];
            const vb = key === 'label' ? b[key].toLowerCase() : b[key];
            return asc ? (va < vb ? -1 : va > vb ? 1 : 0) : (va > vb ? -1 : va < vb ? 1 : 0);
        });
        this._renderInlineCompareTable();
    }

    _renderInlineCompareTable() {
        if (!this._inlineCompareData) return;
        const { title, headers, refRows, sortableRows, mode } = this._inlineCompareData;
        const allRows = [...refRows, ...sortableRows];

        const container = document.getElementById('geInlineCompareTable');
        const titleEl = document.getElementById('geInlineCompareTitle');
        const bodyEl = document.getElementById('geInlineCompareBody');

        titleEl.textContent = title;

        let maxAbs = 0;
        allRows.forEach(r => { const abs = Math.abs(r.delta); if (abs > maxAbs) maxAbs = abs; });
        if (maxAbs === 0) maxAbs = 1;

        let html = '<table style="border-collapse:collapse; font-size:11px;">';
        html += '<thead><tr>';
        headers.forEach((h, i) => {
            let arrow = '';
            if (this._inlineSortCol === i) arrow = this._inlineSortAsc ? ' ▲' : ' ▼';
            html += `<th onclick="app.sortInlineCompare(${i})" style="text-align:left; padding:2px 6px; border-bottom:2px solid #6366f1; font-size:10px; cursor:pointer; white-space:nowrap;">${h}${arrow}</th>`;
        });
        html += '</tr></thead><tbody>';

        allRows.forEach(row => {
            const delta = row.delta;
            const intensity = Math.min(Math.abs(delta) / maxAbs, 1);
            let r, g, b;
            if (delta < 0) {
                r = 255; g = Math.round(255 - 140 * intensity); b = Math.round(255 - 140 * intensity);
            } else {
                r = Math.round(255 - 140 * intensity); g = Math.round(255 - 50 * intensity); b = Math.round(255 - 140 * intensity);
            }
            const bgColor = `rgb(${r},${g},${b})`;
            const bold = row.isAll || row.isRef ? 'font-weight:600;' : '';
            const rowMode = row.clickMode || mode;
            const clickVal = rowMode === 'tissue' ? row.tissue : rowMode === 'fusion' ? (row.fusion ?? '') : row.hotspot;
            const clickFn = rowMode === 'tissue'
                ? `app.onInlineCompareTissueClick('${(clickVal || '').replace(/'/g, "\\'")}')`
                : rowMode === 'fusion'
                ? `app.onInlineCompareFusionClick('${(clickVal || '').replace(/'/g, "\\'")}')`
                : `app.onInlineCompareHotspotClick('${(clickVal || '').replace(/'/g, "\\'")}')`;

            html += `<tr onclick="${clickFn}" style="cursor:pointer; ${bold}">`;
            html += `<td style="padding:2px 6px; border-bottom:1px solid #e5e7eb;">${row.label}</td>`;
            html += `<td style="padding:2px 6px; border-bottom:1px solid #e5e7eb; text-align:right;">${row.nWT}</td>`;
            html += `<td style="padding:2px 6px; border-bottom:1px solid #e5e7eb; text-align:right; color:#6b7280;">${row.meanWT !== undefined ? row.meanWT.toFixed(2) : ''}</td>`;
            html += `<td style="padding:2px 6px; border-bottom:1px solid #e5e7eb; text-align:right;">${row.nMut}</td>`;
            html += `<td style="padding:2px 6px; border-bottom:1px solid #e5e7eb; text-align:right; color:#6b7280;">${row.meanMut !== undefined ? row.meanMut.toFixed(2) : ''}</td>`;
            html += `<td style="padding:2px 6px; border-bottom:1px solid #e5e7eb; text-align:right; background:${bgColor};">${delta.toFixed(2)}</td>`;
            html += '</tr>';
        });

        html += '</tbody></table>';
        bodyEl.innerHTML = html;
        container.style.display = '';
    }

    onInlineCompareTissueClick(tissue) {
        // If we're in subtype grouping mode, set the subtype filter instead
        const activeTissueFilter = document.getElementById('geTissueFilter')?.value || '';
        const activeLineage = activeTissueFilter || this.mutationResults?.lineageFilter;
        if (activeLineage && this.cellLineMetadata?.primaryDisease) {
            // Check if the clicked value is a subtype (not a lineage)
            const lineages = this.cellLineMetadata?.lineage ? new Set(Object.values(this.cellLineMetadata.lineage)) : new Set();
            if (!lineages.has(tissue)) {
                document.getElementById('geSubtypeFilter').value = tissue;
                this._keepInlineCompare = true;
                this.showGeneEffectDistribution(this.currentGeneEffectGene);
                return;
            }
        }
        document.getElementById('geTissueFilter').value = tissue;
        this.updateGeSubtypeFilter();
        this._keepInlineCompare = true;
        this.showGeneEffectDistribution(this.currentGeneEffectGene);
    }

    onInlineCompareHotspotClick(hotspot) {
        document.getElementById('geHotspotFilter').value = hotspot;
        this._keepInlineCompare = true;
        this.showGeneEffectDistribution(this.currentGeneEffectGene);
    }

    onInlineCompareFusionClick(fusion) {
        document.getElementById('geFusionFilter').value = fusion;
        this._keepInlineCompare = true;
        this.showGeneEffectDistribution(this.currentGeneEffectGene);
    }

    updateExcludedTissues() {
        const checkboxes = document.querySelectorAll('#tissueExcludeList input[type="checkbox"]');
        this.excludedTissues = new Set();
        checkboxes.forEach(cb => {
            if (cb.checked) {
                this.excludedTissues.add(cb.value);
            }
        });
    }

    populateTissueExcludeList() {
        const container = document.getElementById('tissueExcludeList');
        if (!container || !this.cellLineMetadata?.lineage) return;

        // Get unique lineages
        const lineages = [...new Set(Object.values(this.cellLineMetadata.lineage))].sort();
        container.innerHTML = '';

        lineages.forEach(lineage => {
            const label = document.createElement('label');
            label.style.cssText = 'display: flex; align-items: center; gap: 4px; font-size: 11px; padding: 1px 0; cursor: pointer;';
            const cb = document.createElement('input');
            cb.type = 'checkbox';
            cb.value = lineage;
            cb.addEventListener('change', () => this.updateExcludedTissues());
            label.appendChild(cb);
            label.appendChild(document.createTextNode(lineage));
            container.appendChild(label);
        });
    }
    // ===== Full-Screen Compare Modal =====

    showMutationCompareByTissue() {
        if (!this.mutationResults) return;
        const mr = this.mutationResults;
        const hotspotGene = mr.hotspotGene;
        const isTranslocation = mr.isTranslocation;
        const mutationData = isTranslocation
            ? this.translocations?.geneData?.[hotspotGene]
            : this.mutations?.geneData?.[hotspotGene];
        if (!mutationData) return;

        const cellLines = this.metadata.cellLines;

        // Determine if we should group by subtype (when a lineage filter is active)
        const groupBySubtype = !!(mr.lineageFilter) && !!this.cellLineMetadata?.primaryDisease;

        // Build group -> { cellLine indices by mutation status } map
        const groupMap = {};
        cellLines.forEach((cellLine, idx) => {
            // Apply lineage/sublineage filters from analysis params
            if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) return;
            if (mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) return;
            if (mr.excludedTissues && mr.excludedTissues.size > 0) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                if (lineage && mr.excludedTissues.has(lineage)) return;
            }
            if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
                const addMutData = this.mutations?.geneData?.[mr.additionalHotspot];
                if (addMutData) {
                    const addMutLevel = addMutData.mutations[cellLine] || 0;
                    if (mr.additionalHotspotLevel === '0' && addMutLevel !== 0) return;
                    if (mr.additionalHotspotLevel === '1' && addMutLevel !== 1) return;
                    if (mr.additionalHotspotLevel === '2' && addMutLevel < 2) return;
                    if (mr.additionalHotspotLevel === '1+2' && addMutLevel === 0) return;
                }
            }
            const groupKey = groupBySubtype
                ? (this.cellLineMetadata.primaryDisease[cellLine] || 'Unknown')
                : (this.cellLineMetadata?.lineage?.[cellLine] || 'Unknown');
            const mutLevel = isTranslocation
                ? (mutationData.translocations[cellLine] || 0)
                : (mutationData.mutations[cellLine] || 0);
            if (!groupMap[groupKey]) groupMap[groupKey] = { wt: [], mut: [], total: 0 };
            groupMap[groupKey].total++;
            if (mutLevel === 0) groupMap[groupKey].wt.push(idx);
            else groupMap[groupKey].mut.push(idx);
        });

        // Also build "All" column
        const allWT = [], allMut = [];
        Object.values(groupMap).forEach(t => { allWT.push(...t.wt); allMut.push(...t.mut); });

        // Build columns: "All" + each group
        const cols = [{ label: 'All', wtIdx: allWT, mutIdx: allMut, totalCells: allWT.length + allMut.length, nWT: allWT.length, nMut: allMut.length, isRef: true }];
        Object.entries(groupMap).forEach(([group, data]) => {
            cols.push({ label: group, wtIdx: data.wt, mutIdx: data.mut, totalCells: data.total, nWT: data.wt.length, nMut: data.mut.length, tissue: group });
        });

        const typeLabel = isTranslocation ? 'Translocation/Fusion' : 'Hotspot Mutational';
        const groupLabel = groupBySubtype ? 'Subtype' : 'Tissue';
        // Store data for rendering
        this._compareModalData = {
            cols,
            genes: mr.significantResults.map(r => r.gene),
            hotspotGene,
            mode: 'tissue',
            title: `Compare by ${groupLabel} — ${hotspotGene} ${typeLabel} Analysis`,
            isTranslocation
        };
        this._compareModalMode = 'tissue';
        this._compareSortCol = null;
        this._compareSortAsc = true;

        this.renderCompareModal();
        document.getElementById('mutCompareModal').style.display = '';
    }

    showMutationCompareByHotspot() {
        if (!this.mutationResults) return;
        const mr = this.mutationResults;
        const mainHotspot = mr.hotspotGene;
        const isTranslocation = mr.isTranslocation;
        const mainMutData = isTranslocation
            ? this.translocations?.geneData?.[mainHotspot]
            : this.mutations?.geneData?.[mainHotspot];
        if (!mainMutData) return;

        const cellLines = this.metadata.cellLines;

        // Build base filtered cell indices
        const baseCells = [];
        cellLines.forEach((cellLine, idx) => {
            if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) return;
            if (mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) return;
            if (mr.excludedTissues && mr.excludedTissues.size > 0) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                if (lineage && mr.excludedTissues.has(lineage)) return;
            }
            if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
                const addMutData = this.mutations?.geneData?.[mr.additionalHotspot];
                if (addMutData) {
                    const addMutLevel = addMutData.mutations[cellLine] || 0;
                    if (mr.additionalHotspotLevel === '0' && addMutLevel !== 0) return;
                    if (mr.additionalHotspotLevel === '1' && addMutLevel !== 1) return;
                    if (mr.additionalHotspotLevel === '2' && addMutLevel < 2) return;
                    if (mr.additionalHotspotLevel === '1+2' && addMutLevel === 0) return;
                }
            }
            const mainMut = isTranslocation
                ? (mainMutData.translocations[cellLine] || 0)
                : (mainMutData.mutations[cellLine] || 0);
            baseCells.push({ cellLine, idx, mainMut });
        });

        // "All" column — no hotspot filter
        const allWT = baseCells.filter(c => c.mainMut === 0).map(c => c.idx);
        const allMut = baseCells.filter(c => c.mainMut >= 1).map(c => c.idx);
        const cols = [{ label: 'All', wtIdx: allWT, mutIdx: allMut, totalCells: baseCells.length, nWT: allWT.length, nMut: allMut.length, isRef: true }];

        // For each other hotspot gene
        this.mutations?.genes?.forEach(hGene => {
            if (hGene === mainHotspot && !isTranslocation) return;
            const hMutData = this.mutations.geneData[hGene];
            if (!hMutData) return;
            const filtered = baseCells.filter(c => (hMutData.mutations[c.cellLine] || 0) >= 1);
            const wt = filtered.filter(c => c.mainMut === 0).map(c => c.idx);
            const mut = filtered.filter(c => c.mainMut >= 1).map(c => c.idx);
            if (wt.length > 0 || mut.length > 0) {
                cols.push({ label: hGene, wtIdx: wt, mutIdx: mut, totalCells: filtered.length, nWT: wt.length, nMut: mut.length, hotspot: hGene });
            }
        });

        const typeLabel = isTranslocation ? 'Translocation/Fusion' : 'Hotspot Mutational';
        this._compareModalData = {
            cols,
            genes: mr.significantResults.map(r => r.gene),
            hotspotGene: mainHotspot,
            mode: 'hotspot',
            title: `Compare by Hotspot — ${mainHotspot} ${typeLabel} Analysis`,
            isTranslocation
        };
        this._compareModalMode = 'hotspot';
        this._compareSortCol = null;
        this._compareSortAsc = true;

        this.renderCompareModal();
        document.getElementById('mutCompareModal').style.display = '';
    }

    showMutationCompareByFusion() {
        if (!this.mutationResults) return;
        if (!this.translocations?.genes?.length) return;
        const mr = this.mutationResults;
        const mainHotspot = mr.hotspotGene;
        const isTranslocation = mr.isTranslocation;
        const mainMutData = isTranslocation
            ? this.translocations?.geneData?.[mainHotspot]
            : this.mutations?.geneData?.[mainHotspot];
        if (!mainMutData) return;

        const cellLines = this.metadata.cellLines;

        // Build base filtered cell indices
        const baseCells = [];
        cellLines.forEach((cellLine, idx) => {
            if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) return;
            if (mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) return;
            if (mr.excludedTissues && mr.excludedTissues.size > 0) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                if (lineage && mr.excludedTissues.has(lineage)) return;
            }
            if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
                const addMutData = this.mutations?.geneData?.[mr.additionalHotspot];
                if (addMutData) {
                    const addMutLevel = addMutData.mutations[cellLine] || 0;
                    if (mr.additionalHotspotLevel === '0' && addMutLevel !== 0) return;
                    if (mr.additionalHotspotLevel === '1' && addMutLevel !== 1) return;
                    if (mr.additionalHotspotLevel === '2' && addMutLevel < 2) return;
                    if (mr.additionalHotspotLevel === '1+2' && addMutLevel === 0) return;
                }
            }
            if (mr.additionalTransGene && mr.additionalTransLevel !== 'all') {
                const addTransData = this.translocations?.geneData?.[mr.additionalTransGene]?.translocations;
                if (addTransData) {
                    const tLevel = addTransData[cellLine] || 0;
                    if (mr.additionalTransLevel === '0' && tLevel !== 0) return;
                    if (mr.additionalTransLevel === '1' && tLevel !== 1) return;
                    if (mr.additionalTransLevel === '2' && tLevel < 2) return;
                    if (mr.additionalTransLevel === '1+2' && tLevel < 1) return;
                }
            }
            // Check oncoprint multi-gene filters
            if (!this._cellLinePassesOncoprintFilters(cellLine)) return;
            const mainMut = isTranslocation
                ? (mainMutData.translocations[cellLine] || 0)
                : (mainMutData.mutations[cellLine] || 0);
            baseCells.push({ cellLine, idx, mainMut });
        });

        // "All" column — no fusion filter
        const allWT = baseCells.filter(c => c.mainMut === 0).map(c => c.idx);
        const allMut = baseCells.filter(c => c.mainMut >= 1).map(c => c.idx);
        const cols = [{ label: 'All', wtIdx: allWT, mutIdx: allMut, totalCells: baseCells.length, nWT: allWT.length, nMut: allMut.length, isRef: true }];

        // Pre-filter: only fusion genes with fused cells in baseCells
        const baseCellIds = new Set(baseCells.map(c => c.cellLine));
        for (const tGene of this.translocations.genes) {
            if (tGene === mainHotspot && isTranslocation) continue;
            const tData = this.translocations.geneData[tGene];
            if (!tData) continue;
            let hasFused = false;
            for (const cl of baseCellIds) {
                if (tData.translocations[cl] && tData.translocations[cl] > 0) { hasFused = true; break; }
            }
            if (!hasFused) continue;
            const filtered = baseCells.filter(c => (tData.translocations[c.cellLine] || 0) >= 1);
            const wt = filtered.filter(c => c.mainMut === 0).map(c => c.idx);
            const mut = filtered.filter(c => c.mainMut >= 1).map(c => c.idx);
            if (wt.length > 0 || mut.length > 0) {
                cols.push({ label: tGene, wtIdx: wt, mutIdx: mut, totalCells: filtered.length, nWT: wt.length, nMut: mut.length, hotspot: tGene });
            }
        }

        const typeLabel = isTranslocation ? 'Translocation/Fusion' : 'Hotspot Mutational';
        this._compareModalData = {
            cols,
            genes: mr.significantResults.map(r => r.gene),
            hotspotGene: mainHotspot,
            mode: 'hotspot',
            title: `Compare by Fusion — ${mainHotspot} ${typeLabel} Analysis`,
            isTranslocation
        };
        this._compareModalMode = 'hotspot';
        this._compareSortCol = null;
        this._compareSortAsc = true;

        this.renderCompareModal();
        document.getElementById('mutCompareModal').style.display = '';
    }

    renderCompareModal() {
        if (!this._compareModalData) return;
        const d = this._compareModalData;
        const minN = parseInt(document.getElementById('mutCompareMinN')?.value) || 5;

        // Filter columns by minN (only requires min N mutated cells)
        const filteredCols = d.cols.filter(c => c.mutIdx.length >= minN);
        this._compareModalCols = filteredCols;

        // Compute delta matrix: genes × cols
        const deltaMatrix = [];
        d.genes.forEach(gene => {
            const geneIdx = this.geneIndex.get(gene.toUpperCase());
            if (geneIdx === undefined) { deltaMatrix.push(null); return; }
            const row = {};
            filteredCols.forEach(col => {
                const wtVals = col.wtIdx.map(i => this.geneEffects[geneIdx * this.nCellLines + i]).filter(v => !isNaN(v));
                const mutVals = col.mutIdx.map(i => this.geneEffects[geneIdx * this.nCellLines + i]).filter(v => !isNaN(v));
                if (wtVals.length >= 3 && mutVals.length >= minN) {
                    const meanWT = wtVals.reduce((a, b) => a + b, 0) / wtVals.length;
                    const meanMut = mutVals.reduce((a, b) => a + b, 0) / mutVals.length;
                    row[col.label] = meanMut - meanWT;
                } else {
                    row[col.label] = null;
                }
            });
            deltaMatrix.push(row);
        });

        // Build gene list with indices
        let geneRows = d.genes.map((gene, i) => ({ gene, deltas: deltaMatrix[i], idx: i })).filter(r => r.deltas !== null);

        // Sort if active
        if (this._compareSortCol !== null) {
            const colLabel = this._compareSortCol;
            const asc = this._compareSortAsc;
            geneRows.sort((a, b) => {
                const va = a.deltas[colLabel];
                const vb = b.deltas[colLabel];
                if (va === null && vb === null) return 0;
                if (va === null) return 1;
                if (vb === null) return -1;
                return asc ? va - vb : vb - va;
            });
        }

        // Find global max |delta| for color scaling
        let maxAbs = 0;
        geneRows.forEach(r => {
            Object.values(r.deltas).forEach(v => { if (v !== null && Math.abs(v) > maxAbs) maxAbs = Math.abs(v); });
        });
        if (maxAbs === 0) maxAbs = 1;

        // Title and info
        document.getElementById('mutCompareModalTitle').textContent = d.title;
        const isTrans = d.isTranslocation;
        const mutLabel = isTrans ? 'fused' : 'mutated';
        const mutNoun = isTrans ? 'fusion' : 'mutation';
        const mutAbbr = isTrans ? 'Fus' : 'Mut';
        const modeLabel = d.mode === 'tissue' ? 'tissue/cancer type' : (isTrans ? 'fusion gene' : 'hotspot mutation');
        document.getElementById('mutCompareModalInfo').innerHTML =
            `<b>Δ GE = Mean GE(${mutLabel}) − Mean GE(WT)</b> for ${d.hotspotGene} ${mutNoun}, stratified by ${modeLabel}. ` +
            `<span style="color:#dc2626;">Red = more essential when ${mutLabel}</span>, <span style="color:#16a34a;">Green = less essential</span>. ` +
            `${geneRows.length} genes × ${filteredCols.length} ${d.mode === 'tissue' ? 'tissues' : 'hotspots'} | Min ${mutLabel} cells: ${minN} | ` +
            `Click cell to inspect, hover column header for N(WT)/N(${mutAbbr})`;

        // Build table HTML
        let html = '<table style="border-collapse:collapse; font-size:11px; width:auto; max-width:100%; margin:0 auto;">';
        html += '<thead style="position:sticky; top:0; z-index:1;"><tr><th style="padding:4px 8px; background:#f0fdf4; border-bottom:2px solid #5a9f4a; position:sticky; left:0; z-index:2; text-align:left;">Gene</th>';
        filteredCols.forEach((col, ci) => {
            let arrow = '';
            if (this._compareSortCol === col.label) arrow = this._compareSortAsc ? ' ▲' : ' ▼';
            const isRef = col.isRef ? 'font-weight:700;' : '';
            html += `<th onclick="app.sortCompareModal('${col.label.replace(/'/g, "\\'")}')" onmouseenter="app.showColumnTooltip(event, ${ci})" onmouseleave="app.hideColumnTooltip()" style="padding:4px 6px; background:#f0fdf4; border-bottom:2px solid #5a9f4a; cursor:pointer; white-space:nowrap; font-size:10px; ${isRef}">${col.label}${arrow}<br><span style="font-weight:400; font-size:9px; color:#6b7280;">${col.nWT}/${col.nMut}</span></th>`;
        });
        html += '</tr></thead><tbody>';

        geneRows.forEach(r => {
            html += '<tr>';
            html += `<td onmouseenter="app.showGeneTooltip(event, '${r.gene}')" onmouseleave="app.hideGeneTooltip()" style="padding:3px 8px; border-bottom:1px solid #e5e7eb; position:sticky; left:0; background:white; font-weight:500; cursor:pointer; white-space:nowrap; color:#5a9f4a;" onclick="app.openCompareInspect('${r.gene}', '', '', '${d.mode}')">${r.gene}</td>`;
            filteredCols.forEach(col => {
                const v = r.deltas[col.label];
                if (v === null) {
                    html += '<td style="padding:3px 6px; border-bottom:1px solid #e5e7eb; text-align:center; color:#ccc;">-</td>';
                } else {
                    const intensity = Math.min(Math.abs(v) / maxAbs, 1);
                    let red, green, blue;
                    if (v < 0) {
                        // Negative (more essential) → white to red
                        red = 255;
                        green = Math.round(255 - 140 * intensity);
                        blue = Math.round(255 - 140 * intensity);
                    } else {
                        // Positive (less essential) → white to green
                        red = Math.round(255 - 140 * intensity);
                        green = Math.round(255 - 50 * intensity);
                        blue = Math.round(255 - 140 * intensity);
                    }
                    const bgColor = `rgb(${red},${green},${blue})`;
                    const clickArg = d.mode === 'tissue' ? `'${r.gene}', '${(col.tissue || '').replace(/'/g, "\\'")}', '', 'tissue'` : `'${r.gene}', '', '${(col.hotspot || '').replace(/'/g, "\\'")}', 'hotspot'`;
                    html += `<td onclick="app.openCompareInspect(${clickArg})" style="padding:3px 6px; border-bottom:1px solid #e5e7eb; text-align:center; background:${bgColor}; cursor:pointer; font-size:10px;" title="${v.toFixed(3)}">${v.toFixed(2)}</td>`;
                }
            });
            html += '</tr>';
        });

        html += '</tbody></table>';
        document.getElementById('mutCompareModalBody').innerHTML = html;
    }

    downloadCompareCSV() {
        if (!this._compareModalData || !this._compareModalCols) return;
        const d = this._compareModalData;
        const filteredCols = this._compareModalCols;
        const minN = parseInt(document.getElementById('mutCompareMinN')?.value) || 5;

        // Recompute delta matrix (same logic as renderCompareModal)
        const geneRows = [];
        d.genes.forEach(gene => {
            const geneIdx = this.geneIndex.get(gene.toUpperCase());
            if (geneIdx === undefined) return;
            const row = { gene };
            filteredCols.forEach(col => {
                const wtVals = col.wtIdx.map(i => this.geneEffects[geneIdx * this.nCellLines + i]).filter(v => !isNaN(v));
                const mutVals = col.mutIdx.map(i => this.geneEffects[geneIdx * this.nCellLines + i]).filter(v => !isNaN(v));
                if (wtVals.length >= 3 && mutVals.length >= minN) {
                    const meanWT = wtVals.reduce((a, b) => a + b, 0) / wtVals.length;
                    const meanMut = mutVals.reduce((a, b) => a + b, 0) / mutVals.length;
                    row[col.label] = (meanMut - meanWT).toFixed(4);
                } else {
                    row[col.label] = '';
                }
            });
            geneRows.push(row);
        });

        // Build CSV
        const colLabels = filteredCols.map(c => c.label);
        let csv = `# Compare by ${d.mode === 'tissue' ? 'Tissue' : 'Hotspot'} — ${d.hotspotGene} Mutation\n`;
        csv += `# Min N: ${minN}\n`;
        csv += `# Date: ${new Date().toISOString().slice(0, 10)}\n`;
        csv += '#\n';
        csv += 'Gene,' + colLabels.map(l => `"${l}"`).join(',') + '\n';
        geneRows.forEach(r => {
            csv += r.gene + ',' + colLabels.map(l => r[l]).join(',') + '\n';
        });

        const filename = `compare_${d.mode}_${d.hotspotGene}_${new Date().toISOString().slice(0, 10)}.csv`;
        this.downloadFile(csv, filename, 'text/csv');
    }

    sortCompareModal(colLabel) {
        if (this._compareSortCol === colLabel) {
            this._compareSortAsc = !this._compareSortAsc;
        } else {
            this._compareSortCol = colLabel;
            this._compareSortAsc = true;
        }
        this.renderCompareModal();
    }

    showColumnTooltip(event, colIdx) {
        this.hideColumnTooltip();
        if (!this._compareModalCols || !this._compareModalCols[colIdx]) return;
        const col = this._compareModalCols[colIdx];
        const tooltip = document.createElement('div');
        tooltip.id = 'columnTooltip';
        tooltip.style.cssText = 'position:fixed; z-index:10001; background:white; border:1px solid #d1d5db; border-radius:8px; padding:8px 12px; max-width:250px; box-shadow:0 4px 12px rgba(0,0,0,0.15); font-size:11px; line-height:1.5;';
        const mutAbbr = this._compareModalData?.isTranslocation ? 'Fus' : 'Mut';
        tooltip.innerHTML = `<b>${col.label}</b><br>Total cell lines: ${col.totalCells}<br>N(WT): ${col.nWT}<br>N(${mutAbbr}): ${col.nMut}`;
        const x = Math.min(event.clientX + 10, window.innerWidth - 270);
        const y = Math.min(event.clientY + 10, window.innerHeight - 100);
        tooltip.style.left = x + 'px';
        tooltip.style.top = y + 'px';
        document.body.appendChild(tooltip);
    }

    hideColumnTooltip() {
        const existing = document.getElementById('columnTooltip');
        if (existing) existing.remove();
    }

    // ===== Enrichr Pathway Analysis =====

    getGenesFromTable(source) {
        const geneSet = new Set();
        const tableMap = {
            'correlations':   { bodyId: 'correlationsBody', geneCols: [0, 1] },
            'clusters':       { bodyId: 'clustersBody', geneCols: [0] },
            'mutations':      { bodyId: 'mutationTableBody', geneCols: [1] }
        };
        const config = tableMap[source];
        if (!config) return [];

        // For mutations, apply the enrichr filter dropdown
        if (source === 'mutations' && this.mutationResults) {
            return this._getFilteredMutationGenes();
        }

        const tbody = document.getElementById(config.bodyId);
        if (!tbody) return [];
        tbody.querySelectorAll('tr').forEach(row => {
            if (row.style.display === 'none') return;
            config.geneCols.forEach(colIdx => {
                const cell = row.cells[colIdx];
                if (cell) {
                    const gene = cell.textContent.trim().replace(/\*$/, '');
                    if (gene) geneSet.add(gene);
                }
            });
        });
        return [...geneSet];
    }

    _getFilteredMutationGenes() {
        const filter = document.getElementById('mutEnrichrFilter')?.value || 'all';
        const mr = this.mutationResults;
        if (!mr) return [];
        const pThreshold = mr.pThreshold || 0.05;

        // "filtered" — use visible rows from table (respects column filters)
        if (filter === 'filtered') {
            const geneSet = new Set();
            const tbody = document.getElementById('mutationTableBody');
            if (tbody) {
                tbody.querySelectorAll('tr').forEach(row => {
                    if (row.style.display === 'none') return;
                    const cell = row.cells[1];
                    if (cell) geneSet.add(cell.textContent.trim());
                });
            }
            return [...geneSet];
        }

        // Get sorted results (respects current table sort order)
        const sortCol = this._mutTableSortCol || 'diff_mut';
        const sortDir = this._mutTableSortDir || 'asc';
        let results = [...(mr.significantResults || [])];
        results.sort((a, b) => {
            const va = a[sortCol] ?? 0;
            const vb = b[sortCol] ?? 0;
            if (sortCol === 'gene') return sortDir === 'asc' ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
            return sortDir === 'asc' ? va - vb : vb - va;
        });

        // Top-N based on current sort
        const topMatch = filter.match(/^top(\d+)$/);
        if (topMatch) {
            return results.slice(0, parseInt(topMatch[1])).map(r => r.gene);
        }

        // Filter by specific comparison
        return results.filter(r => {
            switch (filter) {
                case 'p_mut':    return r.p_mut < pThreshold;
                case 'p_2':      return r.p_2 < pThreshold;
                case 'p_2v1':    return r.p_2v1 < pThreshold;
                case 'diff_neg': return r.p_mut < pThreshold && r.diff_mut < 0;
                case 'diff_pos': return r.p_mut < pThreshold && r.diff_mut > 0;
                default:         return true;
            }
        }).map(r => r.gene);
    }

    async openEnrichr(source) {
        const genes = this.getGenesFromTable(source);
        if (genes.length < 2) {
            this.showCopyNotification('Need at least 2 genes for Enrichr analysis');
            return;
        }

        const modal = document.getElementById('enrichrModal');
        const content = document.getElementById('enrichrContent');
        const title = document.getElementById('enrichrTitle');
        title.textContent = `Enrichr — ${genes.length} genes`;
        content.innerHTML = '<div style="text-align:center; padding:60px; color:#aaa;"><div style="font-size:24px; margin-bottom:12px;">⏳</div>Submitting to Enrichr...</div>';
        modal.style.display = 'block';

        try {
            await this.submitToEnrichr(genes);
        } catch (err) {
            content.innerHTML = `<div style="text-align:center; padding:60px; color:#ef4444;">Failed to connect to Enrichr. Check internet connection.<br><small style="color:#888;">${err.message}</small></div>`;
        }
    }

    async submitToEnrichr(genes) {
        const formData = new FormData();
        formData.append('list', genes.join('\n'));
        formData.append('description', 'Correlate');

        const addRes = await fetch('https://maayanlab.cloud/Enrichr/addList', {
            method: 'POST',
            body: formData
        });
        if (!addRes.ok) throw new Error(`Enrichr addList failed: ${addRes.status}`);
        const { userListId } = await addRes.json();

        const libraries = [
            { key: 'GO_Biological_Process_2023', label: 'GO:BP' },
            { key: 'GO_Molecular_Function_2023', label: 'GO:MF' },
            { key: 'GO_Cellular_Component_2023', label: 'GO:CC' },
            { key: 'KEGG_2021_Human', label: 'KEGG' },
            { key: 'Reactome_2022', label: 'Reactome' }
        ];

        const content = document.getElementById('enrichrContent');
        content.innerHTML = '<div style="text-align:center; padding:60px; color:#aaa;"><div style="font-size:24px; margin-bottom:12px;">⏳</div>Fetching enrichment results...</div>';

        const results = {};
        await Promise.all(libraries.map(async (lib) => {
            const res = await fetch(`https://maayanlab.cloud/Enrichr/enrich?userListId=${userListId}&backgroundType=${lib.key}`);
            if (!res.ok) throw new Error(`Enrichr enrich failed for ${lib.key}: ${res.status}`);
            const data = await res.json();
            results[lib.key] = data[lib.key] || [];
        }));

        this._enrichrData = { genes, results, libraries };
        this._enrichrSortState = {};
        this.renderEnrichrResults(libraries[0].key);
    }

    renderEnrichrResults(activeLibrary) {
        const { results, libraries } = this._enrichrData;
        const tabsEl = document.getElementById('enrichrTabs');
        const contentEl = document.getElementById('enrichrContent');

        // Render tabs
        tabsEl.innerHTML = libraries.map(lib => {
            const active = lib.key === activeLibrary;
            const total = results[lib.key]?.length || 0;
            const sig = (results[lib.key] || []).filter(r => r[6] < 0.05).length;
            return `<button data-lib="${lib.key}" style="padding:5px 12px; font-size:12px; border:1px solid ${active ? '#5a9f4a' : '#555'}; background:${active ? '#5a9f4a' : 'transparent'}; color:${active ? '#fff' : '#ccc'}; border-radius:4px; cursor:pointer;">${lib.label} (${sig}/${total})</button>`;
        }).join('');

        tabsEl.querySelectorAll('button').forEach(btn => {
            btn.addEventListener('click', () => this.renderEnrichrResults(btn.dataset.lib));
        });

        const rows = results[activeLibrary] || [];
        if (rows.length === 0) {
            contentEl.innerHTML = '<div style="text-align:center; padding:60px; color:#aaa;">No pathways found for this library.</div>';
            return;
        }

        // Parse rows: [rank, term, pvalue, zscore, combined_score, overlapping_genes, adj_pvalue, old_p, old_adj_p]
        let parsed = rows.map(r => ({
            rank: r[0],
            term: r[1],
            pValue: r[2],
            zScore: r[3],
            combinedScore: r[4],
            genes: r[5],
            adjPValue: r[6]
        })).filter(r => r.adjPValue < 0.05);

        if (parsed.length === 0) {
            contentEl.innerHTML = '<div style="text-align:center; padding:60px; color:#aaa;">No significant pathways (adj p &lt; 0.05)</div>';
            return;
        }

        // Sort
        const sortState = this._enrichrSortState[activeLibrary] || { col: 'adjPValue', asc: true };
        this._enrichrSortState[activeLibrary] = sortState;
        this._enrichrActiveLibrary = activeLibrary;

        parsed.sort((a, b) => {
            let va = a[sortState.col], vb = b[sortState.col];
            if (typeof va === 'string') return sortState.asc ? va.localeCompare(vb) : vb.localeCompare(va);
            return sortState.asc ? va - vb : vb - va;
        });

        const columns = [
            { key: 'rank', label: '#' },
            { key: 'term', label: 'Term' },
            { key: 'pValue', label: 'P-value' },
            { key: 'adjPValue', label: 'Adj P-value' },
            { key: 'zScore', label: 'Z-score' },
            { key: 'combinedScore', label: 'Combined' },
            { key: 'genesDisplay', label: 'Genes' },
            { key: 'geneCount', label: '# Genes' }
        ];

        const arrow = (col) => {
            if (sortState.col !== col) return '';
            return sortState.asc ? ' ▲' : ' ▼';
        };

        let html = '<table style="width:100%; border-collapse:collapse; font-size:12px;">';
        html += '<thead><tr>';
        columns.forEach(c => {
            const sortable = c.key !== 'genesDisplay';
            const sortKey = c.key === 'geneCount' ? 'geneCount' : c.key;
            html += `<th data-enrichr-sort="${sortable ? sortKey : ''}" style="padding:6px 8px; text-align:left; border-bottom:2px solid #555; color:#aaa; cursor:${sortable ? 'pointer' : 'default'}; white-space:nowrap; font-size:11px;">${c.label}${arrow(sortKey)}</th>`;
        });
        html += '</tr></thead><tbody>';

        parsed.forEach((row, i) => {
            const geneList = Array.isArray(row.genes) ? row.genes.join(', ') : String(row.genes);
            const geneCount = Array.isArray(row.genes) ? row.genes.length : 0;
            const isLong = geneList.length > 60;
            const truncatedGenes = isLong ? geneList.substring(0, 60) + '...' : geneList;

            html += `<tr style="border-bottom:1px solid #333;">`;
            html += `<td style="padding:5px 8px; color:#888;">${row.rank}</td>`;
            html += `<td style="padding:5px 8px; max-width:350px; overflow:hidden; text-overflow:ellipsis;" title="${row.term}">${row.term}</td>`;
            html += `<td style="padding:5px 8px; font-family:monospace; font-size:11px;">${row.pValue.toExponential(2)}</td>`;
            html += `<td style="padding:5px 8px; font-family:monospace; font-size:11px; color:#5a9f4a; font-weight:bold;">${row.adjPValue.toExponential(2)}</td>`;
            html += `<td style="padding:5px 8px;">${row.zScore.toFixed(2)}</td>`;
            html += `<td style="padding:5px 8px;">${row.combinedScore.toFixed(1)}</td>`;
            if (isLong) {
                html += `<td class="enrichr-gene-cell" data-row-idx="${i}" style="padding:5px 8px; max-width:200px; font-size:11px; cursor:pointer;" title="Click to expand"><span class="enrichr-genes-short">${truncatedGenes} <span style="color:#7cabcf;">[${geneCount}]</span></span><span class="enrichr-genes-full" style="display:none; white-space:normal; word-break:break-word;">${geneList} <span style="color:#7cabcf;">[collapse]</span></span></td>`;
            } else {
                html += `<td style="padding:5px 8px; max-width:200px; overflow:hidden; text-overflow:ellipsis; font-size:11px;">${geneList}</td>`;
            }
            html += `<td style="padding:5px 8px; text-align:center;">${geneCount}</td>`;
            html += '</tr>';

            // Store geneCount for sorting
            row.geneCount = geneCount;
        });

        html += '</tbody></table>';
        contentEl.innerHTML = html;

        // Wire sort clicks
        contentEl.querySelectorAll('th[data-enrichr-sort]').forEach(th => {
            const col = th.dataset.enrichrSort;
            if (!col) return;
            th.addEventListener('click', () => {
                const st = this._enrichrSortState[activeLibrary];
                if (st.col === col) {
                    st.asc = !st.asc;
                } else {
                    st.col = col;
                    st.asc = true;
                }
                this.renderEnrichrResults(activeLibrary);
            });
        });

        // Wire gene list expand/collapse clicks
        contentEl.querySelectorAll('.enrichr-gene-cell').forEach(td => {
            td.addEventListener('click', () => {
                const short = td.querySelector('.enrichr-genes-short');
                const full = td.querySelector('.enrichr-genes-full');
                if (!short || !full) return;
                const isExpanded = full.style.display !== 'none';
                short.style.display = isExpanded ? '' : 'none';
                full.style.display = isExpanded ? 'none' : '';
                td.title = isExpanded ? 'Click to expand' : 'Click to collapse';
            });
        });
    }

    downloadEnrichrCSV() {
        if (!this._enrichrData || !this._enrichrActiveLibrary) return;
        const lib = this._enrichrActiveLibrary;
        const rows = this._enrichrData.results[lib] || [];
        if (rows.length === 0) return;

        let csv = 'Rank,Term,P-value,Adjusted P-value,Z-score,Combined Score,Genes,Num Genes\n';
        rows.forEach(r => {
            const genes = Array.isArray(r[5]) ? r[5].join(';') : String(r[5]);
            const geneCount = Array.isArray(r[5]) ? r[5].length : 0;
            const term = String(r[1]).replace(/"/g, '""');
            csv += `${r[0]},"${term}",${r[2]},${r[6]},${r[3]},${r[4]},"${genes}",${geneCount}\n`;
        });

        const libLabel = this._enrichrData.libraries.find(l => l.key === lib)?.label || lib;
        this.downloadFile(csv, `enrichr_${libLabel}.csv`, 'text/csv');
    }

    openCompareInspect(gene, tissue, hotspot, mode) {
        // Close compare modal
        document.getElementById('mutCompareModal').style.display = 'none';
        // Open inspect with appropriate filters
        if (mode === 'tissue' && tissue) {
            this.showGeneEffectDistribution(gene, tissue, '');
        } else if (mode === 'hotspot' && hotspot) {
            this.showGeneEffectDistribution(gene, '', hotspot);
        } else {
            this.showGeneEffectDistribution(gene);
        }
    }

    // ===== Cell Line Browser =====

    setupCellLineBrowserEvents() {
        document.getElementById('showCellLineBrowser').addEventListener('click', (e) => {
            e.preventDefault();
            this.openCellLineBrowser();
        });
        document.getElementById('clbCloseBtn').addEventListener('click', () => this.closeCellLineBrowser());
        document.getElementById('cellLineBrowserModal').addEventListener('click', (e) => {
            if (e.target.id === 'cellLineBrowserModal') this.closeCellLineBrowser();
        });

        let clbSearchTimer;
        document.getElementById('clbSearch').addEventListener('input', () => {
            clearTimeout(clbSearchTimer);
            clbSearchTimer = setTimeout(() => this.renderCellLineList(), 150);
        });
        document.getElementById('clbTissueFilter').addEventListener('change', () => {
            this.renderCellLineList();
        });
        document.getElementById('clbSubtypeFilter').addEventListener('change', () => this.renderCellLineList());
        document.getElementById('clbHotspotFilter').addEventListener('input', () => {
            const val = document.getElementById('clbHotspotFilter').value.trim();
            if (val === '' || this.mutations?.geneData?.[val]) this.renderCellLineList();
        });
        document.getElementById('clbTranslocationFilter').addEventListener('input', () => {
            const val = document.getElementById('clbTranslocationFilter').value.trim();
            if (val === '' || this.translocations?.geneData?.[val]) this.renderCellLineList();
        });
        document.getElementById('clbOncoprintBtn')?.addEventListener('click', () => this.showOncoprint('clb'));

        let clbGeneTimer;
        document.getElementById('clbSortGene').addEventListener('input', () => {
            clearTimeout(clbGeneTimer);
            clbGeneTimer = setTimeout(() => {
                const hasGene = document.getElementById('clbSortGene').value.trim() !== '';
                document.getElementById('clbSortDir').style.display = hasGene ? '' : 'none';
                this.renderCellLineList();
            }, 200);
        });

        document.getElementById('clbSortDir').addEventListener('click', () => {
            this._clbSortAsc = !this._clbSortAsc;
            document.getElementById('clbSortDir').innerHTML = this._clbSortAsc ? '&#x25B2;' : '&#x25BC;';
            this.renderCellLineList();
        });

        document.getElementById('clbResetFilters').addEventListener('click', () => {
            document.getElementById('clbSearch').value = '';
            document.getElementById('clbTissueFilter').value = '';
            document.getElementById('clbSubtypeFilter').value = '';
            document.getElementById('clbHotspotFilter').value = '';
            document.getElementById('clbTranslocationFilter').value = '';
            document.getElementById('clbSortGene').value = '';
            document.getElementById('clbSortDir').style.display = 'none';
            this._oncoprintFilters = {};
            this._activeOncoprintFilters = null;
            this._oncoprintSyncFilters?.();
            this.renderCellLineList();
        });

        document.getElementById('clbTopN').addEventListener('change', () => {
            if (this._clbInspectedCellLine) this.showCellLineDetail(this._clbInspectedCellLine);
        });

        document.getElementById('clbListContainer').addEventListener('click', (e) => {
            const entry = e.target.closest('.clb-entry');
            if (!entry) return;
            const clId = entry.dataset.clid;
            if (e.target.type === 'checkbox') {
                if (e.target.checked) this._clbSelectedCellLines.add(clId);
                else this._clbSelectedCellLines.delete(clId);
                entry.classList.toggle('clb-selected', e.target.checked);
                this.updateClbSelectionCount();
            } else {
                this.showCellLineDetail(clId);
            }
        });

        document.getElementById('clbSelectVisible').addEventListener('click', () => {
            this._clbVisibleCellLines.forEach(cl => this._clbSelectedCellLines.add(cl));
            this.renderCellLineList();
            this.updateClbSelectionCount();
        });
        document.getElementById('clbDeselectAll').addEventListener('click', () => {
            this._clbSelectedCellLines.clear();
            this.renderCellLineList();
            this.updateClbSelectionCount();
        });
        document.getElementById('clbDetailGeneLists').addEventListener('click', (e) => {
            const link = e.target.closest('.clb-gene-link');
            if (link) {
                e.preventDefault();
                this._geHighlightCellLine = this._clbInspectedCellLine;
                this.openGeneEffectModal(link.dataset.gene, 'tissue');
                this._applyParamFiltersToGEModal();
                return;
            }
            const enrichrBtn = e.target.closest('.clb-enrichr-btn');
            if (enrichrBtn) {
                const genes = this._clbGeneLists?.[enrichrBtn.dataset.list];
                if (genes && genes.length >= 2) {
                    const modal = document.getElementById('enrichrModal');
                    const content = document.getElementById('enrichrContent');
                    document.getElementById('enrichrTitle').textContent = `Enrichr — ${genes.length} genes`;
                    content.innerHTML = '<div style="text-align:center; padding:60px; color:#aaa;"><div style="font-size:24px; margin-bottom:12px;">⏳</div>Submitting to Enrichr...</div>';
                    modal.style.display = 'block';
                    this.submitToEnrichr(genes).catch(err => {
                        content.innerHTML = `<div style="text-align:center; padding:60px; color:#ef4444;">Failed to connect to Enrichr.<br><small style="color:#888;">${err.message}</small></div>`;
                    });
                } else {
                    this.showCopyNotification('Need at least 2 genes for Enrichr analysis');
                }
                return;
            }
        });

        document.getElementById('clbExportMinimal').addEventListener('click', () => this.exportCellLineBrowserCSV('minimal'));
        document.getElementById('clbExportFull').addEventListener('click', () => this.exportCellLineBrowserCSV('full'));

        // Gene tooltips on gene links in detail panel
        const geneLists = document.getElementById('clbDetailGeneLists');
        geneLists.addEventListener('mouseenter', (e) => {
            const link = e.target.closest('.clb-gene-link');
            if (!link) return;
            this._clbGeneTooltipTimer = setTimeout(() => {
                this.showGeneTooltip(e, link.dataset.gene);
            }, 400);
        }, true);
        geneLists.addEventListener('mouseleave', (e) => {
            const link = e.target.closest('.clb-gene-link');
            if (!link) return;
            clearTimeout(this._clbGeneTooltipTimer);
            this.hideGeneTooltip();
        }, true);

        // Click handler for gene links in top section (mutations/fusions)
        const detailTop = document.getElementById('clbDetailTop');
        detailTop.addEventListener('click', (e) => {
            const link = e.target.closest('.clb-gene-link');
            if (!link) return;
            e.preventDefault();
            this._geHighlightCellLine = this._clbInspectedCellLine;
            this.openGeneEffectModal(link.dataset.gene, 'tissue');
            this._applyParamFiltersToGEModal();
        });

        // Gene tooltips on mutation/fusion gene names in top section
        detailTop.addEventListener('mouseenter', (e) => {
            const link = e.target.closest('.gene-hover');
            if (!link) return;
            this._clbGeneTooltipTimer = setTimeout(() => {
                this.showGeneTooltip(e, link.dataset.gene);
            }, 400);
        }, true);
        detailTop.addEventListener('mouseleave', (e) => {
            const link = e.target.closest('.gene-hover');
            if (!link) return;
            clearTimeout(this._clbGeneTooltipTimer);
            this.hideGeneTooltip();
        }, true);
    }

    openCellLineBrowser() {
        this._oncoprintFilters = {};
        this._activeOncoprintFilters = null;
        document.getElementById('clbSearch').value = '';
        document.getElementById('clbSortGene').value = '';
        document.getElementById('clbTissueFilter').value = '';
        document.getElementById('clbSubtypeFilter').value = '';
        document.getElementById('clbSubtypeFilter').style.display = 'none';
        document.getElementById('clbHotspotFilter').value = '';
        document.getElementById('clbTranslocationFilter').value = '';
        document.getElementById('clbDetailPanel').classList.remove('active');
        document.getElementById('clbDetailContent').style.display = 'none';
        document.getElementById('clbDetailPlaceholder').style.display = '';
        this._clbInspectedCellLine = null;

        this.renderCellLineList();
        this.updateClbSelectionCount();
        this.updateClbFilterCounts(this.metadata.cellLines);
        document.getElementById('cellLineBrowserModal').style.display = 'flex';
    }

    closeCellLineBrowser() {
        document.getElementById('cellLineBrowserModal').style.display = 'none';
        this._oncoprintFilters = {};
        this._activeOncoprintFilters = null;
        this._oncoprintSyncFilters?.();
        document.getElementById('oncoprintPopup')?.remove();
        document.querySelectorAll('.nav-link').forEach(t => t.classList.remove('active'));
        document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
        document.querySelector('[data-tab="network"]').classList.add('active');
        document.getElementById('tab-network').classList.add('active');
    }

    renderCellLineList() {
        const search = document.getElementById('clbSearch').value.trim().toLowerCase();
        const tissue = document.getElementById('clbTissueFilter').value;
        const subtype = document.getElementById('clbSubtypeFilter').value;
        const hotspotGene = document.getElementById('clbHotspotFilter').value;
        const transGene = document.getElementById('clbTranslocationFilter').value;

        const hotspotMuts = hotspotGene && this.mutations?.geneData?.[hotspotGene]?.mutations;
        const transMuts = transGene && this.translocations?.geneData?.[transGene]?.translocations;

        let filtered = this.metadata.cellLines.filter(cl => {
            if (tissue && this.getCellLineLineage(cl) !== tissue) return false;
            if (subtype && this.getCellLineSublineage(cl) !== subtype) return false;
            if (hotspotMuts && !(hotspotMuts[cl] >= 1)) return false;
            if (transMuts && !(transMuts[cl] >= 1)) return false;
            if (!this._cellLinePassesOncoprintFilters(cl)) return false;
            if (search) {
                const name = this.getCellLineName(cl).toLowerCase();
                const lin = this.getCellLineLineage(cl).toLowerCase();
                if (!name.includes(search) && !lin.includes(search) && !cl.toLowerCase().includes(search)) return false;
            }
            return true;
        });

        const sortGene = document.getElementById('clbSortGene').value.trim().toUpperCase();
        const sortGeneIdx = sortGene ? this.geneIndex.get(sortGene) : undefined;
        let geMap = null;
        if (sortGeneIdx !== undefined) {
            geMap = new Map();
            for (const cl of filtered) {
                const clIdx = this.metadata.cellLines.indexOf(cl);
                if (clIdx >= 0) {
                    const val = this.geneEffects[sortGeneIdx * this.nCellLines + clIdx];
                    geMap.set(cl, (!isNaN(val) && val !== -999) ? val : NaN);
                }
            }
            const dir = this._clbSortAsc ? 1 : -1;
            filtered.sort((a, b) => {
                const va = geMap.get(a), vb = geMap.get(b);
                if (isNaN(va) && isNaN(vb)) return 0;
                if (isNaN(va)) return 1;
                if (isNaN(vb)) return -1;
                return (va - vb) * dir;
            });
        } else {
            filtered.sort((a, b) => this.getCellLineName(a).localeCompare(this.getCellLineName(b)));
        }
        this._clbVisibleCellLines = filtered;
        this._updateClbActiveFilterLabel();

        const container = document.getElementById('clbList');
        if (filtered.length === 0) {
            container.innerHTML = '<div class="clb-no-results">No cell lines match your filters</div>';
            this.updateClbFilterCounts();
            return;
        }

        const html = filtered.map(cl => {
            const name = this.getCellLineName(cl);
            const lin = this.getCellLineLineage(cl);
            const selected = this._clbSelectedCellLines.has(cl);
            const inspected = this._clbInspectedCellLine === cl;
            const cls = ['clb-entry'];
            if (selected) cls.push('clb-selected');
            if (inspected) cls.push('clb-inspected');
            const geVal = geMap ? geMap.get(cl) : null;
            const geStr = geVal !== null && !isNaN(geVal) ? `<span style="font-size:10px; color:#666; margin-left:auto; flex-shrink:0;">${geVal.toFixed(2)}</span>` : '';
            return `<div class="${cls.join(' ')}" data-clid="${cl}">` +
                `<input type="checkbox"${selected ? ' checked' : ''}>` +
                `<span class="clb-entry-name" title="${name}">${name}</span>` +
                `<span class="clb-entry-tissue">${lin}</span>${geStr}</div>`;
        }).join('');
        container.innerHTML = html;
        this.updateClbFilterCounts();
    }

    updateClbSubtypeFilter() {
        const tissue = document.getElementById('clbTissueFilter').value;
        const subSelect = document.getElementById('clbSubtypeFilter');

        if (!tissue) {
            subSelect.style.display = 'none';
            subSelect.innerHTML = '<option value="">All subtypes</option>';
            return;
        }

        const prefix = `${tissue}|`;
        const subtypes = Object.keys(this.subLineageCounts)
            .filter(k => k.startsWith(prefix))
            .map(k => ({ name: k.slice(prefix.length), count: this.subLineageCounts[k] }))
            .sort((a, b) => a.name.localeCompare(b.name));

        if (subtypes.length === 0) {
            subSelect.style.display = 'none';
            subSelect.innerHTML = '<option value="">All subtypes</option>';
            return;
        }

        subSelect.innerHTML = '<option value="">All subtypes</option>';
        subtypes.forEach(s => {
            const opt = document.createElement('option');
            opt.value = s.name;
            opt.textContent = `${s.name} (n=${s.count})`;
            subSelect.appendChild(opt);
        });
        subSelect.style.display = '';
    }

    _updateClbActiveFilterLabel() {
        const el = document.getElementById('clbActiveFilters');
        if (!el) return;
        const parts = [];
        const tissue = document.getElementById('clbTissueFilter')?.value;
        const subtype = document.getElementById('clbSubtypeFilter')?.value;
        const hotspot = document.getElementById('clbHotspotFilter')?.value;
        const trans = document.getElementById('clbTranslocationFilter')?.value;
        if (tissue) parts.push(`<span style="background:#dbeafe;color:#1e40af;padding:1px 6px;border-radius:10px;">${tissue}${subtype ? ' · ' + subtype : ''}</span>`);
        if (hotspot) parts.push(`<span style="background:#dcfce7;color:#16a34a;padding:1px 6px;border-radius:10px;">${hotspot} mutated</span>`);
        if (trans) parts.push(`<span style="background:#fae8ff;color:#7c3aed;padding:1px 6px;border-radius:10px;">${trans} fused</span>`);
        if (this._activeOncoprintFilters) {
            const shown = new Set([hotspot, trans].filter(Boolean));
            for (const f of this._activeOncoprintFilters) {
                if (!shown.has(f.gene)) {
                    const bg = f.state === 'mut' ? '#dcfce7' : '#fef2f2';
                    const color = f.state === 'mut' ? '#16a34a' : '#dc2626';
                    parts.push(`<span style="background:${bg};color:${color};padding:1px 6px;border-radius:10px;">${f.gene} ${f.state === 'mut' ? 'Mut' : 'WT'}</span>`);
                }
            }
        }
        if (parts.length > 0) {
            const n = this._clbVisibleCellLines?.length || 0;
            const total = this.metadata?.cellLines?.length || 0;
            el.innerHTML = `<span style="color:#6b7280;">Showing ${n} of ${total}:</span> ${parts.join(' ')}`;
            el.style.display = 'flex';
        } else {
            el.style.display = 'none';
        }
    }

    updateClbFilterCounts() {
        const tissue = document.getElementById('clbTissueFilter').value;
        const subtype = document.getElementById('clbSubtypeFilter').value;
        const hotspotGene = document.getElementById('clbHotspotFilter').value;
        const transGene = document.getElementById('clbTranslocationFilter').value;

        const hotspotMuts = hotspotGene && this.mutations?.geneData?.[hotspotGene]?.mutations;
        const transMuts = transGene && this.translocations?.geneData?.[transGene]?.translocations;
        const allCls = this.metadata.cellLines;

        const getBaseSet = (excludeFilter) => {
            return allCls.filter(cl => {
                if (excludeFilter !== 'tissue' && tissue && this.getCellLineLineage(cl) !== tissue) return false;
                if (excludeFilter !== 'subtype' && subtype && this.getCellLineSublineage(cl) !== subtype) return false;
                if (excludeFilter !== 'hotspot' && hotspotMuts && !(hotspotMuts[cl] >= 1)) return false;
                if (excludeFilter !== 'translocation' && transMuts && !(transMuts[cl] >= 1)) return false;
                if (excludeFilter !== 'oncoprint' && !this._cellLinePassesOncoprintFilters(cl)) return false;
                return true;
            });
        };

        // Update tissue filter counts
        const tissueBase = getBaseSet('tissue');
        const tissueSelect = document.getElementById('clbTissueFilter');
        const tissueVal = tissueSelect.value;
        tissueSelect.innerHTML = `<option value="">All tissues (n=${tissueBase.length})</option>`;
        const tissueCounts = {};
        tissueBase.forEach(cl => {
            const lin = this.getCellLineLineage(cl);
            tissueCounts[lin] = (tissueCounts[lin] || 0) + 1;
        });
        Object.keys(tissueCounts).sort((a, b) => tissueCounts[b] - tissueCounts[a]).forEach(lin => {
            const opt = document.createElement('option');
            opt.value = lin;
            opt.textContent = `${lin} (n=${tissueCounts[lin]})`;
            if (lin === tissueVal) opt.selected = true;
            tissueSelect.appendChild(opt);
        });

        // Update subtype filter counts
        const subtypeBase = getBaseSet('subtype');
        const subSelect = document.getElementById('clbSubtypeFilter');
        const subVal = subSelect.value;
        if (tissue) {
            const subCounts = {};
            subtypeBase.forEach(cl => {
                if (this.getCellLineLineage(cl) === tissue) {
                    const sub = this.getCellLineSublineage(cl);
                    if (sub) subCounts[sub] = (subCounts[sub] || 0) + 1;
                }
            });
            subSelect.innerHTML = '<option value="">All subtypes</option>';
            Object.keys(subCounts).sort((a, b) => subCounts[b] - subCounts[a]).forEach(sub => {
                const opt = document.createElement('option');
                opt.value = sub;
                opt.textContent = `${sub} (n=${subCounts[sub]})`;
                if (sub === subVal) opt.selected = true;
                subSelect.appendChild(opt);
            });
            subSelect.style.display = Object.keys(subCounts).length > 0 ? '' : 'none';
        } else {
            subSelect.style.display = 'none';
            subSelect.innerHTML = '<option value="">All subtypes</option>';
        }

        // Update hotspot filter counts (datalist)
        if (this.mutations?.geneData) {
            const hotspotBase = getBaseSet('hotspot');
            const hotspotBaseSet = new Set(hotspotBase);
            const hotspotDatalist = document.getElementById('clbHotspotList');
            if (hotspotDatalist) {
                const geneCounts = [];
                Object.keys(this.mutations.geneData).forEach(gene => {
                    const muts = this.mutations.geneData[gene].mutations;
                    let n = 0;
                    for (const [cl, v] of Object.entries(muts)) {
                        if (v > 0 && hotspotBaseSet.has(cl)) n++;
                    }
                    if (n > 0) geneCounts.push({ gene, n });
                });
                geneCounts.sort((a, b) => b.n - a.n);
                let hHtml = '';
                geneCounts.forEach(({ gene, n }) => {
                    hHtml += `<option value="${gene}">${gene} (n=${n})</option>`;
                });
                hotspotDatalist.innerHTML = hHtml;
            }
        }

        // Update translocation filter counts (datalist)
        if (this._fusionGeneCounts && this.translocations?.geneData) {
            const transBase = getBaseSet('translocation');
            const transBaseSet = new Set(transBase);
            const transDatalist = document.getElementById('clbTranslocationList');
            if (transDatalist) {
                const geneCounts = [];
                this._fusionGeneCounts.forEach(({ gene }) => {
                    const td = this.translocations.geneData[gene]?.translocations;
                    if (!td) return;
                    let n = 0;
                    for (const [cl, v] of Object.entries(td)) {
                        if (v >= 1 && transBaseSet.has(cl)) n++;
                    }
                    if (n > 0) geneCounts.push({ gene, n });
                });
                let tHtml = '';
                geneCounts.forEach(({ gene, n }) => {
                    tHtml += `<option value="${gene}">${gene} (n=${n})</option>`;
                });
                transDatalist.innerHTML = tHtml;
            }
        }
    }

    showCellLineDetail(cellLineId) {
        this._clbInspectedCellLine = cellLineId;
        document.querySelectorAll('#clbList .clb-entry').forEach(el => {
            el.classList.toggle('clb-inspected', el.dataset.clid === cellLineId);
        });

        const panel = document.getElementById('clbDetailPanel');
        const content = document.getElementById('clbDetailContent');
        const placeholder = document.getElementById('clbDetailPlaceholder');
        panel.classList.add('active');
        placeholder.style.display = 'none';
        content.style.display = '';

        const name = this.getCellLineName(cellLineId);
        const lineage = this.getCellLineLineage(cellLineId);
        const sublineage = this.getCellLineSublineage(cellLineId);
        const N = parseInt(document.getElementById('clbTopN').value) || 10;

        const mutGenes = [];
        if (this.mutations?.geneData) {
            for (const gene of Object.keys(this.mutations.geneData)) {
                if (this.mutations.geneData[gene].mutations?.[cellLineId] >= 1) mutGenes.push(gene);
            }
        }

        const fusionGenes = [];
        if (this.translocations?.geneData) {
            for (const gene of Object.keys(this.translocations.geneData)) {
                if (this.translocations.geneData[gene].translocations?.[cellLineId] >= 1) fusionGenes.push(gene);
            }
        }

        const clIdx = this.metadata.cellLines.indexOf(cellLineId);
        let sum = 0, count = 0, min = Infinity, max = -Infinity;
        const geneVals = [];
        if (clIdx >= 0) {
            for (let g = 0; g < this.nGenes; g++) {
                const val = this.geneEffects[g * this.nCellLines + clIdx];
                if (!isNaN(val) && val !== -999) {
                    sum += val; count++;
                    if (val < min) min = val;
                    if (val > max) max = val;
                    geneVals.push({ gene: this.geneNames[g], val });
                }
            }
        }
        const mean = count > 0 ? sum / count : NaN;
        geneVals.sort((a, b) => a.val - b.val);

        let top = `<h4>${name}</h4>`;
        top += `<div class="clb-detail-id">${cellLineId}</div>`;
        top += `<div class="clb-detail-section">`;
        top += `<div class="clb-stat-row"><span class="clb-stat-label">Tissue</span><span class="clb-stat-value">${lineage || '-'}</span></div>`;
        top += `<div class="clb-stat-row"><span class="clb-stat-label">Subtype</span><span class="clb-stat-value">${sublineage || '-'}</span></div>`;
        top += `</div>`;

        top += `<div class="clb-detail-section"><strong>Hotspot Mutations (${mutGenes.length})</strong>`;
        top += `<div style="color:var(--gray-500); font-size:11px;">${mutGenes.length > 0 ? mutGenes.map(g => `<span class="gene-hover clb-gene-link" data-gene="${g}" style="cursor:help;">${g}</span>`).join(', ') : 'None'}</div></div>`;

        top += `<div class="clb-detail-section"><strong>Fusions (${fusionGenes.length})</strong>`;
        top += `<div style="color:var(--gray-500); font-size:11px;">${fusionGenes.length > 0 ? fusionGenes.map(g => `<span class="gene-hover clb-gene-link" data-gene="${g}" style="cursor:help;">${g}</span>`).join(', ') : 'None'}</div></div>`;

        top += `<div class="clb-detail-section"><strong>Gene Effect Stats</strong>`;
        top += `<div class="clb-stat-row"><span class="clb-stat-label">Genes</span><span class="clb-stat-value">${count.toLocaleString()}</span></div>`;
        top += `<div class="clb-stat-row"><span class="clb-stat-label">Mean GE</span><span class="clb-stat-value">${this.formatNum(mean)}</span></div>`;
        top += `<div class="clb-stat-row"><span class="clb-stat-label">Range</span><span class="clb-stat-value">${count > 0 ? this.formatNum(min) + ' to ' + this.formatNum(max) : '-'}</span></div>`;
        top += `</div>`;
        document.getElementById('clbDetailTop').innerHTML = top;

        const bottomN = geneVals.slice(0, N);
        const topN = geneVals.slice(-N).reverse();

        this._clbGeneLists = {
            bottom: bottomN.map(g => g.gene),
            top: topN.map(g => g.gene)
        };

        let gl = `<div class="clb-detail-section"><strong>Most Depleted (Bottom ${N})</strong> <button class="clb-enrichr-btn" data-list="bottom">Enrichr</button>`;
        gl += `<div style="font-size:11px;">`;
        bottomN.forEach(({ gene, val }) => {
            gl += `<div class="clb-stat-row"><span class="clb-stat-label"><a class="clb-gene-link" data-gene="${gene}" href="#">${gene}</a></span><span class="clb-stat-value">${this.formatNum(val)}</span></div>`;
        });
        gl += `</div></div>`;

        gl += `<div class="clb-detail-section"><strong>Most Enriched (Top ${N})</strong> <button class="clb-enrichr-btn" data-list="top">Enrichr</button>`;
        gl += `<div style="font-size:11px;">`;
        topN.forEach(({ gene, val }) => {
            gl += `<div class="clb-stat-row"><span class="clb-stat-label"><a class="clb-gene-link" data-gene="${gene}" href="#">${gene}</a></span><span class="clb-stat-value">${this.formatNum(val)}</span></div>`;
        });
        gl += `</div></div>`;

        const filteredCls = this._clbVisibleCellLines;
        const filteredIndices = filteredCls.map(cl => this.metadata.cellLines.indexOf(cl)).filter(i => i >= 0);
        if (filteredIndices.length >= 3 && clIdx >= 0) {
            const zScores = [];
            for (let g = 0; g < this.nGenes; g++) {
                const offset = g * this.nCellLines;
                let s = 0, s2 = 0, n = 0;
                for (const ci of filteredIndices) {
                    const v = this.geneEffects[offset + ci];
                    if (!isNaN(v) && v !== -999) { s += v; s2 += v * v; n++; }
                }
                if (n < 3) continue;
                const mu = s / n;
                const sd = Math.sqrt(s2 / n - mu * mu);
                if (sd < 1e-6) continue;
                const val = this.geneEffects[offset + clIdx];
                if (isNaN(val) || val === -999) continue;
                zScores.push({ gene: this.geneNames[g], z: (val - mu) / sd, val });
            }
            zScores.sort((a, b) => a.z - b.z);

            const filterParts = [];
            const tissueVal = document.getElementById('clbTissueFilter').value;
            const subtypeVal = document.getElementById('clbSubtypeFilter').value;
            const hotspotVal = document.getElementById('clbHotspotFilter').value;
            const transVal = document.getElementById('clbTranslocationFilter').value;
            if (tissueVal) filterParts.push(tissueVal);
            if (subtypeVal) filterParts.push(subtypeVal);
            if (hotspotVal) filterParts.push(hotspotVal + ' Mut');
            if (transVal) filterParts.push(transVal + ' Fused');
            if (this._activeOncoprintFilters) {
                const shown = new Set([hotspotVal, transVal].filter(Boolean));
                for (const f of this._activeOncoprintFilters) {
                    if (!shown.has(f.gene)) filterParts.push(`${f.gene} ${f.state === 'mut' ? 'Mut' : 'WT'}`);
                }
            }
            const filterLabel = `filtered (n=${filteredIndices.length})${filterParts.length > 0 ? ': ' + filterParts.join(', ') : ''}`;

            const extremeLow = zScores.slice(0, N);
            const extremeHigh = zScores.slice(-N).reverse();

            this._clbGeneLists.uniqueLow = extremeLow.map(g => g.gene);
            this._clbGeneLists.uniqueHigh = extremeHigh.map(g => g.gene);

            gl += `<div class="clb-detail-section"><strong>Uniquely Depleted vs ${filterLabel}</strong> <button class="clb-enrichr-btn" data-list="uniqueLow">Enrichr</button>`;
            gl += `<div style="font-size:10px; color:var(--gray-500); margin-bottom:3px;">Lowest z-score vs visible cell lines (n=${filteredIndices.length})</div>`;
            gl += `<div style="font-size:11px;">`;
            extremeLow.forEach(({ gene, z, val }) => {
                gl += `<div class="clb-stat-row"><span class="clb-stat-label"><a class="clb-gene-link" data-gene="${gene}" href="#">${gene}</a></span><span class="clb-stat-value">${this.formatNum(val)} <span style="color:#888;">(z=${this.formatNum(z, 1)})</span></span></div>`;
            });
            gl += `</div></div>`;

            gl += `<div class="clb-detail-section"><strong>Uniquely Enriched vs ${filterLabel}</strong> <button class="clb-enrichr-btn" data-list="uniqueHigh">Enrichr</button>`;
            gl += `<div style="font-size:10px; color:var(--gray-500); margin-bottom:3px;">Highest z-score vs visible cell lines (n=${filteredIndices.length})</div>`;
            gl += `<div style="font-size:11px;">`;
            extremeHigh.forEach(({ gene, z, val }) => {
                gl += `<div class="clb-stat-row"><span class="clb-stat-label"><a class="clb-gene-link" data-gene="${gene}" href="#">${gene}</a></span><span class="clb-stat-value">${this.formatNum(val)} <span style="color:#888;">(z=${this.formatNum(z, 1)})</span></span></div>`;
            });
            gl += `</div></div>`;
        }

        document.getElementById('clbDetailGeneLists').innerHTML = gl;
    }

    updateClbSelectionCount() {
        document.getElementById('clbSelectionCount').textContent = `${this._clbSelectedCellLines.size} selected`;
    }

    exportCellLineBrowserCSV(mode) {
        if (this._clbSelectedCellLines.size === 0) {
            alert('Select at least one cell line to export.');
            return;
        }

        const selectedIds = [...this._clbSelectedCellLines];
        const clIndices = selectedIds.map(cl => this.metadata.cellLines.indexOf(cl)).filter(i => i >= 0);
        const clNames = clIndices.map(i => this.getCellLineName(this.metadata.cellLines[i]));
        const namePart = clNames.length <= 4
            ? '_' + clNames.map(n => n.replace(/[^A-Za-z0-9]/g, '')).join('_')
            : `_${clNames.length}cl`;

        if (mode === 'minimal') {
            const headerParts = ['Gene_Effect'];
            clNames.forEach(n => headerParts.push(n.replace(/,/g, '')));
            const lines = [`# Gene Effect (CRISPR DepMap) matrix for ${clNames.length} selected cell lines`, `# Source: DepMap 25Q3 CRISPRGeneEffect`, `# Date: ${new Date().toISOString().slice(0, 10)}`, headerParts.join(',')];
            for (let g = 0; g < this.nGenes; g++) {
                const row = [this.geneNames[g]];
                for (let c = 0; c < clIndices.length; c++) {
                    const val = this.geneEffects[g * this.nCellLines + clIndices[c]];
                    row.push((!isNaN(val) && val !== -999) ? val.toFixed(4) : '');
                }
                lines.push(row.join(','));
            }
            this.downloadFile(lines.join('\n'), `correlate${namePart}.csv`, 'text/csv');
        } else {
            const N = parseInt(document.getElementById('clbTopN').value) || 10;
            const filteredCls = this._clbVisibleCellLines;
            const filteredIndices = filteredCls.map(cl => this.metadata.cellLines.indexOf(cl)).filter(i => i >= 0);

            const filterParts = [];
            const tissueVal = document.getElementById('clbTissueFilter').value;
            const subtypeVal = document.getElementById('clbSubtypeFilter').value;
            const hotspotVal = document.getElementById('clbHotspotFilter').value;
            const transVal = document.getElementById('clbTranslocationFilter').value;
            if (tissueVal) filterParts.push(tissueVal);
            if (subtypeVal) filterParts.push(subtypeVal);
            if (hotspotVal) filterParts.push(hotspotVal + ' mut');
            if (transVal) filterParts.push(transVal + ' fus');
            const filterLabel = filterParts.length > 0 ? filterParts.join(' / ') : `all (n=${filteredIndices.length})`;

            const lines = ['CellLine,CellLineID,Tissue,Subtype,Category,Rank,Gene,GE,ZScore'];

            for (let c = 0; c < clIndices.length; c++) {
                const clId = this.metadata.cellLines[clIndices[c]];
                const clName = clNames[c].replace(/,/g, '');
                const tissue = this.getCellLineLineage(clId);
                const subtype = this.getCellLineSublineage(clId);
                const ci = clIndices[c];

                if (this.mutations?.geneData) {
                    for (const gene of Object.keys(this.mutations.geneData)) {
                        if (this.mutations.geneData[gene].mutations?.[clId] >= 1) {
                            lines.push(`${clName},${clId},${tissue},${subtype},Hotspot Mutation,,${gene},,`);
                        }
                    }
                }

                if (this.translocations?.geneData) {
                    for (const gene of Object.keys(this.translocations.geneData)) {
                        if (this.translocations.geneData[gene].translocations?.[clId] >= 1) {
                            lines.push(`${clName},${clId},${tissue},${subtype},Fusion,,${gene},,`);
                        }
                    }
                }

                const geneVals = [];
                for (let g = 0; g < this.nGenes; g++) {
                    const val = this.geneEffects[g * this.nCellLines + ci];
                    if (!isNaN(val) && val !== -999) geneVals.push({ gene: this.geneNames[g], val, idx: g });
                }
                geneVals.sort((a, b) => a.val - b.val);

                geneVals.slice(0, N).forEach(({ gene, val }, i) => {
                    lines.push(`${clName},${clId},${tissue},${subtype},Most Depleted,${i + 1},${gene},${val.toFixed(4)},`);
                });

                geneVals.slice(-N).reverse().forEach(({ gene, val }, i) => {
                    lines.push(`${clName},${clId},${tissue},${subtype},Most Enriched,${i + 1},${gene},${val.toFixed(4)},`);
                });

                if (filteredIndices.length >= 3) {
                    const zScores = [];
                    for (let g = 0; g < this.nGenes; g++) {
                        const offset = g * this.nCellLines;
                        let s = 0, s2 = 0, n = 0;
                        for (const fi of filteredIndices) {
                            const v = this.geneEffects[offset + fi];
                            if (!isNaN(v) && v !== -999) { s += v; s2 += v * v; n++; }
                        }
                        if (n < 3) continue;
                        const mu = s / n;
                        const sd = Math.sqrt(s2 / n - mu * mu);
                        if (sd < 1e-6) continue;
                        const val = this.geneEffects[offset + ci];
                        if (isNaN(val) || val === -999) continue;
                        zScores.push({ gene: this.geneNames[g], z: (val - mu) / sd, val });
                    }
                    zScores.sort((a, b) => a.z - b.z);

                    zScores.slice(0, N).forEach(({ gene, val, z }, i) => {
                        lines.push(`${clName},${clId},${tissue},${subtype},Uniquely Depleted vs ${filterLabel},${i + 1},${gene},${val.toFixed(4)},${z.toFixed(2)}`);
                    });

                    zScores.slice(-N).reverse().forEach(({ gene, val, z }, i) => {
                        lines.push(`${clName},${clId},${tissue},${subtype},Uniquely Enriched vs ${filterLabel},${i + 1},${gene},${val.toFixed(4)},${z.toFixed(2)}`);
                    });
                }
            }

            this.downloadFile(lines.join('\n'), `correlate_report${namePart}.csv`, 'text/csv');
        }
    }

    // ============================================================
    // SVG Export Sanitization & Text Outlining
    // ============================================================

    _flattenTspan(tspan) {
        // Recursively flatten nested tspans: collapse style attributes upward
        // If a tspan has exactly one child tspan (and no text siblings), merge them
        while (true) {
            const children = [...tspan.childNodes];
            const childTspans = children.filter(c => c.nodeType === 1 && c.tagName === 'tspan');
            const textNodes = children.filter(c => c.nodeType === 3 && c.textContent.trim());
            // Only flatten if there's exactly one child tspan and no direct text
            if (childTspans.length === 1 && textNodes.length === 0) {
                const child = childTspans[0];
                // Merge child's style into parent
                const childStyle = child.getAttribute('style') || '';
                if (childStyle) {
                    const parentStyle = tspan.getAttribute('style') || '';
                    tspan.setAttribute('style', parentStyle ? parentStyle.replace(/;?\s*$/, '; ') + childStyle : childStyle);
                }
                // Copy child's attributes (except style which we merged)
                for (const attr of child.attributes) {
                    if (attr.name !== 'style' && !tspan.hasAttribute(attr.name)) {
                        tspan.setAttribute(attr.name, attr.value);
                    }
                }
                // Replace child tspan with its contents
                while (child.firstChild) {
                    tspan.insertBefore(child.firstChild, child);
                }
                child.remove();
                // Continue loop to flatten further nested levels
            } else {
                break;
            }
        }
    }

    sanitizeSvgForIllustrator(svgStr) {
        const parser = new DOMParser();
        const doc = parser.parseFromString(svgStr, 'image/svg+xml');
        const svgEl = doc.documentElement;

        // Remove <style> blocks — inline all styles instead
        doc.querySelectorAll('style').forEach(styleEl => {
            const rules = [];
            try {
                const sheet = new CSSStyleSheet();
                sheet.replaceSync(styleEl.textContent);
                for (const rule of sheet.cssRules) {
                    if (rule.selectorText) rules.push(rule);
                }
            } catch (e) {
                // Fallback: just remove the style block
            }
            for (const rule of rules) {
                try {
                    const els = doc.querySelectorAll(rule.selectorText);
                    els.forEach(el => {
                        for (let i = 0; i < rule.style.length; i++) {
                            const prop = rule.style[i];
                            if (!el.style.getPropertyValue(prop)) {
                                el.style.setProperty(prop, rule.style.getPropertyValue(prop));
                            }
                        }
                    });
                } catch (e) { /* skip invalid selectors */ }
            }
            styleEl.remove();
        });

        // Clean all elements: convert colors, fix fonts, strip Plotly artifacts
        const allEls = doc.querySelectorAll('*');
        const rgbToHex = (r, g, b) => '#' + [r, g, b].map(x => parseInt(x).toString(16).padStart(2, '0')).join('');
        const rgbRe = /rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)/g;

        allEls.forEach(el => {
            // Convert rgb in style attributes to hex
            let style = el.getAttribute('style');
            if (style) {
                style = style.replace(rgbRe, (_, r, g, b) => rgbToHex(r, g, b));
                style = style.replace(/white-space:\s*pre;?\s*/gi, '');
                style = style.replace(/cursor:\s*[^;]+;?\s*/gi, '');
                style = style.replace(/pointer-events:\s*[^;]+;?\s*/gi, '');
                style = style.replace(/;\s*;/g, ';').replace(/^\s*;\s*/, '').replace(/;\s*$/, '');
                if (style.trim()) {
                    el.setAttribute('style', style);
                } else {
                    el.removeAttribute('style');
                }
            }
            // Convert rgb in fill/stroke attributes
            ['fill', 'stroke', 'color', 'stop-color'].forEach(attr => {
                const val = el.getAttribute(attr);
                if (val) {
                    const replaced = val.replace(rgbRe, (_, r, g, b) => rgbToHex(r, g, b));
                    if (replaced !== val) el.setAttribute(attr, replaced);
                }
            });

            // font-size: remove "px" suffix for Illustrator compatibility
            const fs = el.style?.fontSize;
            if (fs && fs.endsWith('px')) {
                el.style.fontSize = parseFloat(fs).toString();
            }
            const fsAttr = el.getAttribute('font-size');
            if (fsAttr && fsAttr.endsWith('px')) {
                el.setAttribute('font-size', parseFloat(fsAttr).toString());
            }

            // Strip Plotly-specific class and data- attributes
            el.removeAttribute('class');
            [...el.attributes].forEach(attr => {
                if (attr.name.startsWith('data-')) el.removeAttribute(attr.name);
            });
        });

        // Remove legend/scrollbox clip-paths only
        doc.querySelectorAll('[clip-path]').forEach(el => {
            const cp = el.getAttribute('clip-path');
            if (cp && (cp.includes('legend') || cp.includes('scrollbox'))) {
                el.removeAttribute('clip-path');
            }
        });
        doc.querySelectorAll('clipPath').forEach(cp => {
            const id = cp.getAttribute('id') || '';
            if (id.includes('legend') || id.includes('scrollbox')) cp.remove();
        });

        // Remove empty <defs> and empty <g> groups
        doc.querySelectorAll('defs').forEach(d => { if (!d.children.length) d.remove(); });
        let removed = true;
        while (removed) {
            removed = false;
            doc.querySelectorAll('g').forEach(g => {
                if (!g.children.length && !g.textContent.trim()) { g.remove(); removed = true; }
            });
        }

        // Flatten deeply nested Plotly tspans
        doc.querySelectorAll('text').forEach(textEl => {
            const topTspans = textEl.querySelectorAll(':scope > tspan');
            topTspans.forEach(topTs => {
                this._flattenTspan(topTs);
            });
        });

        // Fix SVG root element for Inkscape/Illustrator compatibility
        const vb = svgEl.getAttribute('viewBox');
        if (vb) {
            const parts = vb.trim().split(/\s+/).map(Number);
            if (parts.length === 4) {
                const vbW = parts[2];
                const vbH = parts[3];
                const maxMm = 280;
                const ptW = vbW;
                const ptH = vbH;
                const mmW = ptW / 72 * 25.4;
                if (mmW > maxMm) {
                    const scale = maxMm / mmW;
                    svgEl.setAttribute('width', (ptW * scale).toFixed(1) + 'pt');
                    svgEl.setAttribute('height', (ptH * scale).toFixed(1) + 'pt');
                } else {
                    svgEl.setAttribute('width', ptW + 'pt');
                    svgEl.setAttribute('height', ptH + 'pt');
                }
            }
        }
        svgEl.removeAttribute('class');
        const rootStyle = svgEl.getAttribute('style');
        if (!rootStyle || rootStyle.trim() === '') svgEl.removeAttribute('style');

        let result = new XMLSerializer().serializeToString(svgEl);
        result = '<?xml version="1.0" encoding="UTF-8"?>\n' + result;
        return result;
    }

    _finalizeSvgForExport(svgStr) {
        return this.sanitizeSvgForIllustrator(svgStr);
    }

    // ============================================================
    // Export Metadata & State Restore
    // ============================================================

    captureAppState() {
        const plotEl = document.getElementById('scatterPlot');
        const lay = plotEl?.layout;
        const anns = lay?.annotations || [];
        const titleAnn = anns.find(a => a._tsRole === 'title');
        const xAnn = anns.find(a => a._tsRole === 'xlabel');
        const yAnn = anns.find(a => a._tsRole === 'ylabel');

        return {
            version: 1,
            gene1: this.currentInspect?.gene1,
            gene2: this.currentInspect?.gene2,
            xRange: [document.getElementById('scatterXmin')?.value, document.getElementById('scatterXmax')?.value],
            yRange: [document.getElementById('scatterYmin')?.value, document.getElementById('scatterYmax')?.value],
            plotWidth: document.getElementById('plotWidth')?.value,
            plotHeight: document.getElementById('plotHeight')?.value,
            hotspotGene: document.getElementById('hotspotGene')?.value,
            hotspotMode: document.getElementById('hotspotMode')?.value,
            transGene: document.getElementById('translocationGene')?.value,
            transMode: document.getElementById('translocationMode')?.value,
            cancerFilter: document.getElementById('scatterCancerFilter')?.value,
            searchTerms: document.getElementById('cellLineSearch')?.value,
            mutFilterGene: document.getElementById('mutationFilterGene')?.value,
            mutFilterLevel: document.getElementById('mutationFilterLevel')?.value,
            colorBy: document.getElementById('colorBySelect')?.value,
            textSettings: (() => {
                const annText = titleAnn?.text || '';
                const titleSizeMatch = annText.match(/^<span style="font-size:(\d+)px">/);
                const subSizeMatch = annText.match(/<br><span style="font-size:(\d+)px/);
                return {
                    titleFontSize: titleSizeMatch ? parseInt(titleSizeMatch[1]) : (titleAnn?.font?.size || 25),
                    subtitleSize: subSizeMatch ? parseInt(subSizeMatch[1]) : 15,
                    xLabelFontSize: xAnn?.font?.size,
                    yLabelFontSize: yAnn?.font?.size,
                    xTickSize: lay?.xaxis?.tickfont?.size,
                    yTickSize: lay?.yaxis?.tickfont?.size,
                    legendSize: lay?.legend?.font?.size,
                    markerSize: plotEl?.data?.[0]?.marker?.size,
                    fontFamily: lay?.font?.family
                };
            })(),
            showZeroLines: document.getElementById('showZeroLines')?.checked,
            showCorrelationLine: document.getElementById('showCorrelationLine')?.checked,
            oncoprintFilters: this._activeOncoprintFilters || null
        };
    }

    async restoreFromExport(file) {
        if (file.name.toLowerCase().endsWith('.png')) {
            return this.restoreFromPng(file);
        }
        return this.restoreFromSvg(file);
    }

    async restoreFromPng(file) {
        const buf = await file.arrayBuffer();
        const stateJson = this._readPngTextChunk(buf, 'correlate-state');
        if (stateJson) {
            let state;
            try { state = JSON.parse(stateJson); } catch (e) { alert('Invalid state data in PNG.'); return; }
            return this._restoreFromState(state);
        }
        const metaJson = this._readPngTextChunk(buf, 'correlate-meta');
        if (metaJson) {
            let meta;
            try { meta = JSON.parse(metaJson); } catch (e) { alert('Invalid metadata in PNG.'); return; }
            return this._handleExportMeta(meta);
        }
        alert('No Correlate data found in this PNG file.');
    }

    async restoreFromSvg(file) {
        const text = await file.text();
        const stateMatch = text.match(/<correlate-state>([\s\S]*?)<\/correlate-state>/);
        if (stateMatch) {
            let state;
            try { state = JSON.parse(stateMatch[1]); } catch (e) { alert('Invalid state data in SVG.'); return; }
            return this._restoreFromState(state);
        }
        const metaMatch = text.match(/<correlate-meta>([\s\S]*?)<\/correlate-meta>/);
        if (metaMatch) {
            let meta;
            try { meta = JSON.parse(metaMatch[1]); } catch (e) { alert('Invalid metadata in SVG.'); return; }
            return this._handleExportMeta(meta);
        }
        alert('No Correlate data found in this file.');
    }

    _resetForRestore() {
        const resetEl = (id, val) => { const el = document.getElementById(id); if (el) el.value = val; };
        resetEl('lineageFilter', '');
        resetEl('subLineageFilter', '');
        resetEl('paramHotspotGene', '');
        resetEl('paramHotspotLevel', 'all');
        resetEl('paramTranslocationGene', '');
        resetEl('paramTranslocationLevel', 'all');
        resetEl('scatterCancerFilter', '');
        resetEl('cellLineSearch', '');
        resetEl('hotspotGene', '');
        resetEl('hotspotMode', 'none');
        resetEl('translocationGene', '');
        resetEl('translocationMode', 'none');
        resetEl('mutationFilterGene', '');
        resetEl('mutationFilterLevel', 'all');
        resetEl('colorBySelect', '');
        const zeroCb = document.getElementById('showZeroLines');
        if (zeroCb) zeroCb.checked = true;
        const corrCb = document.getElementById('showCorrelationLine');
        if (corrCb) corrCb.checked = true;
    }

    _handleExportMeta(meta) {
        this._resetForRestore();

        // Scatter-like exports with gene pair
        if (meta.gene1 && meta.gene2) {
            return this._restoreFromState(meta);
        }

        // Network exports
        if (meta.graphType === 'network' && meta.geneList && meta.geneList.length > 0) {
            document.getElementById('geneTextarea').value = meta.geneList.join('\n');
            if (meta.mode) {
                const radio = document.querySelector(`input[name="analysisMode"][value="${meta.mode}"]`);
                if (radio) radio.checked = true;
            }
            if (meta.cutoff != null) {
                document.getElementById('correlationCutoff').value = meta.cutoff;
                const label = document.getElementById('correlationCutoff')?.nextElementSibling;
                if (label) label.textContent = meta.cutoff;
            }
            // Restore network visual settings
            const ns = meta.networkSettings;
            if (ns) {
                if (ns.fontSize) document.getElementById('netFontSize').value = ns.fontSize;
                if (ns.nodeSize) document.getElementById('netNodeSize').value = ns.nodeSize;
                if (ns.edgeWidth) document.getElementById('netEdgeWidth').value = ns.edgeWidth;
                if (ns.minSlope) document.getElementById('minSlope').value = ns.minSlope;
                if (ns.minCellLines) document.getElementById('minCellLines').value = ns.minCellLines;
                if (ns.lineageFilter) document.getElementById('lineageFilter').value = ns.lineageFilter;
                if (ns.subLineageFilter) {
                    const subEl = document.getElementById('subLineageFilter');
                    if (subEl) subEl.value = ns.subLineageFilter;
                }
                const cbGE = document.getElementById('colorByGeneEffect');
                if (cbGE) cbGE.checked = ns.colorByGeneEffect || false;
                if (ns.colorGEType) {
                    const r = document.querySelector(`input[name="colorGEType"][value="${ns.colorGEType}"]`);
                    if (r) r.checked = true;
                }
                const cbStats = document.getElementById('colorByStats');
                if (cbStats) cbStats.checked = ns.colorByStats || false;
                if (ns.colorStatType) {
                    const r = document.querySelector(`input[name="colorStatType"][value="${ns.colorStatType}"]`);
                    if (r) r.checked = true;
                }
                if (ns.colorScale) {
                    const r = document.querySelector(`input[name="colorScale"][value="${ns.colorScale}"]`);
                    if (r) r.checked = true;
                }
                const tbg = document.getElementById('exportNetworkTransparentBg');
                if (tbg) tbg.checked = ns.transparentBg || false;
            }
            // Restore oncoprint filters
            if (meta.oncoprintFilters && meta.oncoprintFilters.length > 0) {
                this._activeOncoprintFilters = meta.oncoprintFilters;
                this._oncoprintFilters = {};
                for (const f of meta.oncoprintFilters) this._oncoprintFilters[f.gene] = f.state;
                this._oncoprintSyncFilters?.();
            }
            this.runAnalysis();
            return;
        }

        // Mutation inspect exports
        if (meta.graphType === 'mutation_inspect' && meta.gene && meta.hotspotGene) {
            const mutRadio = document.querySelector('input[name="analysisMode"][value="mutation"]');
            if (mutRadio) mutRadio.checked = true;
            // Switch to mutation tab so UI is visible
            document.querySelectorAll('.nav-link').forEach(link => link.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
            const mutTab = document.querySelector('[data-tab="mutation"]');
            if (mutTab) mutTab.classList.add('active');
            const mutContent = document.getElementById('tab-mutation');
            if (mutContent) mutContent.classList.add('active');
            if (meta.isTranslocation) {
                const radio = document.querySelector('input[name="mutAnalysisType"][value="translocation"]');
                if (radio) { radio.checked = true; this.updateMutationAnalysisType?.(); }
                const sel = document.getElementById('translocationHotspotSelect');
                if (sel) sel.value = meta.hotspotGene;
            } else {
                const radio = document.querySelector('input[name="mutAnalysisType"][value="hotspot"]');
                if (radio) { radio.checked = true; this.updateMutationAnalysisType?.(); }
                this.populateMutationHotspotSelector?.();
                const sel = document.getElementById('mutationHotspotSelect');
                if (sel) sel.value = meta.hotspotGene;
            }
            if (meta.lineageFilter) {
                document.getElementById('lineageFilter').value = meta.lineageFilter;
            }
            this.runAnalysis();
            const waitForResults = () => {
                if (this.mutationResults && this.mutationResults.hotspotGene === meta.hotspotGene) {
                    this.showGeneEffectDistribution(meta.gene);
                    // Apply text settings after plot renders
                    if (meta.textSettings) {
                        setTimeout(() => {
                            this._savedScatterTextSettings = meta.textSettings;
                            if (meta.geChartWidthRatio) this.geChartWidthRatio = meta.geChartWidthRatio;
                        }, 200);
                    }
                } else {
                    setTimeout(waitForResults, 200);
                }
            };
            setTimeout(waitForResults, 300);
            return;
        }

        // Gene effect exports
        if (meta.graphType === 'gene_effect' && meta.gene) {
            this.openGeneEffectModal(meta.gene, meta.view || 'tissue');
            return;
        }

        // Fallback
        const graphLabels = {
            'network': 'Correlation Network',
            'tissue_chart': 'By-Tissue Chart',
            'gene_effect': 'Gene Effect Chart',
            'mutation_inspect': 'Mutation Inspect',
            'scatter': 'Scatter Plot'
        };
        const label = graphLabels[meta.graphType] || meta.graphType || 'Unknown';
        const details = [];
        if (meta.gene) details.push(`Gene: ${meta.gene}`);
        if (meta.hotspotGene) details.push(`Hotspot: ${meta.hotspotGene}`);
        if (meta.view) details.push(`View: ${meta.view}`);
        if (meta.geneList) details.push(`Genes: ${meta.geneList.join(', ')}`);
        if (meta.cutoff) details.push(`Cutoff: ${meta.cutoff}`);
        if (meta.date) details.push(`Exported: ${meta.date.slice(0, 10)}`);
        const info = details.length ? '\n' + details.join('\n') : '';
        alert(`Correlate export: ${label}${info}\n\nRestore not yet supported for this graph type.`);
    }

    async _restoreFromState(state) {
        if (!state.gene1 || !state.gene2) { alert('Missing gene information in state.'); return; }
        this._resetForRestore();

        // Restore oncoprint multi-gene filters if saved
        if (state.oncoprintFilters && state.oncoprintFilters.length > 0) {
            this._activeOncoprintFilters = state.oncoprintFilters;
            this._oncoprintFilters = {};
            for (const f of state.oncoprintFilters) {
                this._oncoprintFilters[f.gene] = f.state;
            }
            this._oncoprintSyncFilters?.();
        }

        // Open inspect with the genes
        this.openInspect({ gene1: state.gene1, gene2: state.gene2, correlation: null });

        // Restore settings after a short delay to let the plot render
        setTimeout(() => {
            if (state.xRange) {
                document.getElementById('scatterXmin').value = state.xRange[0];
                document.getElementById('scatterXmax').value = state.xRange[1];
            }
            if (state.yRange) {
                document.getElementById('scatterYmin').value = state.yRange[0];
                document.getElementById('scatterYmax').value = state.yRange[1];
            }
            if (state.plotWidth) document.getElementById('plotWidth').value = state.plotWidth;
            if (state.plotHeight) document.getElementById('plotHeight').value = state.plotHeight;
            if (state.hotspotGene) document.getElementById('hotspotGene').value = state.hotspotGene;
            if (state.hotspotMode) document.getElementById('hotspotMode').value = state.hotspotMode;
            if (state.transGene) document.getElementById('translocationGene').value = state.transGene;
            if (state.transMode) document.getElementById('translocationMode').value = state.transMode;
            if (state.cancerFilter) document.getElementById('scatterCancerFilter').value = state.cancerFilter;
            if (state.searchTerms) document.getElementById('cellLineSearch').value = state.searchTerms;
            if (state.mutFilterGene) document.getElementById('mutationFilterGene').value = state.mutFilterGene;
            if (state.mutFilterLevel) document.getElementById('mutationFilterLevel').value = state.mutFilterLevel;
            if (state.colorBy) {
                const colorByEl = document.getElementById('colorBySelect');
                if (colorByEl) colorByEl.value = state.colorBy;
            }
            if (state.showZeroLines != null) document.getElementById('showZeroLines').checked = state.showZeroLines;
            if (state.showCorrelationLine != null) document.getElementById('showCorrelationLine').checked = state.showCorrelationLine;

            if (state.textSettings) {
                if (state.textSettings.titleFontSize && state.textSettings.titleFontSize < 20 && !state.textSettings.subtitleSize) {
                    state.textSettings.titleFontSize = 25;
                    state.textSettings.subtitleSize = 15;
                }
                this._savedScatterTextSettings = state.textSettings;
            }

            this.updateInspectPlot();
        }, 300);
    }

    // --- PNG metadata: read/write tEXt chunks ---

    _readPngTextChunk(arrayBuffer, keyword) {
        const view = new DataView(arrayBuffer);
        let offset = 8;
        while (offset < view.byteLength) {
            const length = view.getUint32(offset);
            const typeBytes = new Uint8Array(arrayBuffer, offset + 4, 4);
            const type = String.fromCharCode(...typeBytes);
            if (type === 'tEXt') {
                const data = new Uint8Array(arrayBuffer, offset + 8, length);
                const nullIdx = data.indexOf(0);
                if (nullIdx > 0) {
                    const key = new TextDecoder().decode(data.slice(0, nullIdx));
                    if (key === keyword) {
                        return new TextDecoder().decode(data.slice(nullIdx + 1));
                    }
                }
            }
            if (type === 'IEND') break;
            offset += 12 + length;
        }
        return null;
    }

    _addPngTextChunk(pngArrayBuffer, keyword, text) {
        const keyBytes = new TextEncoder().encode(keyword);
        const textBytes = new TextEncoder().encode(text);
        const chunkDataLen = keyBytes.length + 1 + textBytes.length;
        const chunkData = new Uint8Array(chunkDataLen);
        chunkData.set(keyBytes, 0);
        chunkData[keyBytes.length] = 0;
        chunkData.set(textBytes, keyBytes.length + 1);

        const typeBytes = new Uint8Array([0x74, 0x45, 0x58, 0x74]); // tEXt
        const crcInput = new Uint8Array(4 + chunkDataLen);
        crcInput.set(typeBytes, 0);
        crcInput.set(chunkData, 4);
        const crc = this._crc32(crcInput);

        const chunk = new Uint8Array(12 + chunkDataLen);
        new DataView(chunk.buffer).setUint32(0, chunkDataLen);
        chunk.set(typeBytes, 4);
        chunk.set(chunkData, 8);
        new DataView(chunk.buffer).setUint32(8 + chunkDataLen, crc);

        const src = new Uint8Array(pngArrayBuffer);
        const result = new Uint8Array(src.length + chunk.length);
        result.set(src.slice(0, src.length - 12), 0);
        result.set(chunk, src.length - 12);
        result.set(src.slice(src.length - 12), src.length - 12 + chunk.length);
        return result.buffer;
    }

    _crc32(data) {
        if (!this._crc32Table) {
            const table = new Uint32Array(256);
            for (let n = 0; n < 256; n++) {
                let c = n;
                for (let k = 0; k < 8; k++) {
                    c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
                }
                table[n] = c;
            }
            this._crc32Table = table;
        }
        let crc = 0xFFFFFFFF;
        for (let i = 0; i < data.length; i++) {
            crc = this._crc32Table[(crc ^ data[i]) & 0xFF] ^ (crc >>> 8);
        }
        return (crc ^ 0xFFFFFFFF) >>> 0;
    }

    _captureNetworkSettings() {
        return {
            fontSize: parseInt(document.getElementById('netFontSize')?.value) || 16,
            nodeSize: parseInt(document.getElementById('netNodeSize')?.value) || 25,
            edgeWidth: parseInt(document.getElementById('netEdgeWidth')?.value) || 3,
            colorByGeneEffect: document.getElementById('colorByGeneEffect')?.checked || false,
            colorGEType: document.querySelector('input[name="colorGEType"]:checked')?.value || 'signed',
            colorByStats: document.getElementById('colorByStats')?.checked || false,
            colorStatType: document.querySelector('input[name="colorStatType"]:checked')?.value || 'signed_lfc',
            colorScale: document.querySelector('input[name="colorScale"]:checked')?.value || 'all',
            lineageFilter: document.getElementById('lineageFilter')?.value || '',
            subLineageFilter: document.getElementById('subLineageFilter')?.value || '',
            minSlope: document.getElementById('minSlope')?.value || '0.1',
            minCellLines: document.getElementById('minCellLines')?.value || '50',
            transparentBg: document.getElementById('exportNetworkTransparentBg')?.checked || false
        };
    }

    _capturePlotTextSettings(plotDivId) {
        const plotEl = document.getElementById(plotDivId);
        if (!plotEl?.layout) return null;
        const lay = plotEl.layout;
        const anns = lay.annotations || [];
        const titleAnn = anns.find(a => a._tsRole === 'title');
        const xAnn = anns.find(a => a._tsRole === 'xlabel');
        const yAnn = anns.find(a => a._tsRole === 'ylabel');
        const annText = titleAnn?.text || '';
        const titleSizeMatch = annText.match(/font-size:(\d+)px/);
        const subSizeMatch = annText.match(/<br><span style="font-size:(\d+)px/);
        return {
            titleFontSize: titleSizeMatch ? parseInt(titleSizeMatch[1]) : (titleAnn?.font?.size || 25),
            subtitleSize: subSizeMatch ? parseInt(subSizeMatch[1]) : 15,
            xLabelFontSize: xAnn?.font?.size,
            yLabelFontSize: yAnn?.font?.size,
            xTickSize: lay.xaxis?.tickfont?.size,
            yTickSize: lay.yaxis?.tickfont?.size,
            legendSize: lay.legend?.font?.size,
            markerSize: plotEl.data?.[0]?.marker?.size,
            fontFamily: lay.font?.family,
            plotWidth: lay.width,
            plotHeight: lay.height
        };
    }

    _buildExportMetadata(graphType, extra = {}) {
        return {
            app: 'Correlate',
            version: document.getElementById('versionBadge')?.textContent || '',
            graphType,
            date: new Date().toISOString(),
            ...extra
        };
    }

    // ===== Text Settings (Aa) Panel =====

    openTextSettings(plotDivId) {
        const panel = document.getElementById('textSettingsPanel');
        const body = document.getElementById('textSettingsBody');
        this._textSettingsPlotId = plotDivId;

        const plotEl = document.getElementById(plotDivId);
        if (!plotEl?.layout) { body.innerHTML = '<div style="color:#6b7280;">No plot to configure.</div>'; panel.style.display = 'block'; return; }

        const layout = plotEl.layout;

        // Detect if title/labels are annotation-based (scatter) or layout-based (gene effect)
        const anns = layout.annotations || [];
        const titleAnn = anns.find(a => a._tsRole === 'title');
        const xLabelAnn = anns.find(a => a._tsRole === 'xlabel');
        const yLabelAnn = anns.find(a => a._tsRole === 'ylabel');
        const ann0 = titleAnn || anns[0];
        const usesAnnotationTitle = !!titleAnn || (ann0 && !ann0._gateAnnotation && ann0.xref === 'paper');
        // For annotation-based titles, the real title size is in the inline span, not annotation font.size
        let titleSize;
        if (usesAnnotationTitle) {
            const titleSizeMatch = (ann0?.text || '').match(/^<span style="font-size:(\d+)px">/);
            titleSize = titleSizeMatch ? parseInt(titleSizeMatch[1]) : (ann0?.font?.size || 14);
        } else {
            titleSize = layout.title?.font?.size || 14;
        }

        // Extract plain text from HTML annotation text
        const stripHtml = (html) => html ? html.replace(/<[^>]*>/g, '').replace(/\s+/g, ' ').trim() : '';

        // Get title text — for annotations it's HTML with <b> and <br>
        let titleText = '';
        let subtitleText = '';
        let subtitleSize = 15;
        if (usesAnnotationTitle) {
            const raw = ann0.text || '';
            // Split on <br> — first part is title (strip <b>), rest is subtitle
            const parts = raw.split(/<br\s*\/?>/i);
            titleText = stripHtml(parts[0]);
            subtitleText = parts.slice(1).map(stripHtml).join('\n');
            // Extract subtitle font-size from subtitle parts (skip title inline span)
            const subtitleHtml = parts.slice(1).join('<br>');
            const sizeMatch = subtitleHtml.match(/font-size:\s*(\d+)px/);
            if (sizeMatch) subtitleSize = parseInt(sizeMatch[1]);
        } else {
            titleText = typeof layout.title === 'string' ? layout.title : (layout.title?.text || '');
        }

        const xLabel = xLabelAnn ? (xLabelAnn.text || '') : (typeof layout.xaxis?.title === 'string' ? layout.xaxis.title : (layout.xaxis?.title?.text || ''));
        const yLabel = yLabelAnn ? (yLabelAnn.text || '') : (typeof layout.yaxis?.title === 'string' ? layout.yaxis.title : (layout.yaxis?.title?.text || ''));
        const xLabelSize = xLabelAnn ? (xLabelAnn.font?.size || 20) : (layout.xaxis?.title?.font?.size || 12);
        const yLabelSize = yLabelAnn ? (yLabelAnn.font?.size || 20) : (layout.yaxis?.title?.font?.size || 12);
        const xTickSize = layout.xaxis?.tickfont?.size || 17;
        const yTickSize = layout.yaxis?.tickfont?.size || 17;
        const legendSize = layout.legend?.font?.size || 17;
        const markerSize = plotEl.data?.[0]?.marker?.size || 10;
        const hasLegend = layout.showlegend !== false && plotEl.data?.some(t => t.showlegend !== false && t.name);

        // Track visibility states
        this._tsVisible = {
            title: ann0?.visible !== false && (usesAnnotationTitle || !!titleText),
            xLabel: !!xLabel,
            yLabel: !!yLabel,
            legend: hasLegend
        };
        // Store original texts for restore
        this._tsOriginal = { titleText, subtitleText, subtitleSize, xLabel, yLabel, usesAnnotationTitle };

        const sizeRow = (label, id, val, min, max, toggleId, visible) => `
            <div style="display:flex; align-items:center; margin-bottom:5px; gap:4px;">
                ${toggleId ? `<input type="checkbox" id="${toggleId}" ${visible ? 'checked' : ''} onchange="app._tsToggle('${toggleId}')" style="margin:0;">` : '<span style="width:15px;"></span>'}
                <span style="color:#374151;flex:1;min-width:55px;font-size:11px;">${label}</span>
                <div style="display:flex; align-items:center;">
                    <button onclick="app._tsStep('${id}',-1)" style="width:20px;height:20px;border:1px solid #d1d5db;background:#f9fafb;border-radius:4px 0 0 4px;cursor:pointer;font-size:12px;line-height:1;">−</button>
                    <input type="number" id="${id}" value="${val}" min="${min}" max="${max}" style="width:36px;text-align:center;border:1px solid #d1d5db;border-left:none;border-right:none;font-size:10px;padding:1px;-moz-appearance:textfield;appearance:textfield;" oninput="app._tsApply()">
                    <button onclick="app._tsStep('${id}',1)" style="width:20px;height:20px;border:1px solid #d1d5db;background:#f9fafb;border-radius:0 4px 4px 0;cursor:pointer;font-size:12px;line-height:1;">+</button>
                </div>
            </div>`;

        const textRow = (label, id, val, useTextarea = false) => {
            const biButtons = `<span style="display:inline-flex;gap:1px;margin-left:4px;vertical-align:middle;">` +
                `<button id="${id}_bold" onclick="app._tsToggleBold('${id}')" style="width:20px;height:20px;border:1px solid #d1d5db;border-radius:3px;cursor:pointer;font-weight:bold;font-size:11px;background:#f9fafb;line-height:1;" title="Bold">B</button>` +
                `<button id="${id}_italic" onclick="app._tsToggleItalic('${id}')" style="width:20px;height:20px;border:1px solid #d1d5db;border-radius:3px;cursor:pointer;font-style:italic;font-size:11px;background:#f9fafb;line-height:1;" title="Italic">I</button>` +
                `</span>`;
            return useTextarea ? `
            <div style="margin-bottom:5px;">
                <label style="font-size:10px;color:#6b7280;">${label}${biButtons}</label>
                <textarea id="${id}" rows="3" style="width:100%;border:1px solid #d1d5db;border-radius:4px;padding:3px 6px;font-size:11px;margin-top:1px;box-sizing:border-box;resize:vertical;" oninput="app._tsApplyText()">${this._escapeAttr(val)}</textarea>
            </div>` : `
            <div style="margin-bottom:5px;">
                <label style="font-size:10px;color:#6b7280;">${label}${biButtons}</label>
                <input type="text" id="${id}" value="${this._escapeAttr(val)}" style="width:100%;border:1px solid #d1d5db;border-radius:4px;padding:3px 6px;font-size:11px;margin-top:1px;box-sizing:border-box;" oninput="app._tsApplyText()">
            </div>`;
        };

        // Detect current font from layout
        const currentFont = layout.font?.family || layout.title?.font?.family || 'Arial, Helvetica, sans-serif';

        body.innerHTML = `
            <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;">
                <span style="font-weight:600;color:#1f2937;font-size:11px;">Scale All</span>
                <div style="display:flex;align-items:center;gap:2px;">
                    <button onclick="app._tsScaleAll(-1)" style="width:24px;height:22px;border:1px solid #d1d5db;background:#f0fdf4;border-radius:4px 0 0 4px;cursor:pointer;font-size:13px;font-weight:bold;">−</button>
                    <button onclick="app._tsScaleAll(1)" style="width:24px;height:22px;border:1px solid #d1d5db;background:#f0fdf4;border-radius:0 4px 4px 0;cursor:pointer;font-size:13px;font-weight:bold;">+</button>
                </div>
            </div>
            <div style="display:flex;align-items:center;gap:6px;margin-bottom:6px;">
                <span style="font-weight:600;color:#1f2937;font-size:11px;">Font</span>
                <select id="ts_fontFamily" onchange="app._tsApplyFont()" style="flex:1;font-size:11px;padding:2px 4px;border:1px solid #d1d5db;border-radius:4px;">
                    <option value="Arial, Helvetica, sans-serif" ${currentFont.includes('Arial') ? 'selected' : ''} style="font-family:Arial,sans-serif;">Arial</option>
                    <option value="Open Sans, sans-serif" ${currentFont.includes('Open Sans') ? 'selected' : ''} style="font-family:'Open Sans',sans-serif;">Open Sans</option>
                    <option value="Helvetica, Arial, sans-serif" ${currentFont.includes('Helvetica') && !currentFont.includes('Arial,') ? 'selected' : ''} style="font-family:Helvetica,sans-serif;">Helvetica</option>
                    <option value="Georgia, serif" ${currentFont.includes('Georgia') ? 'selected' : ''} style="font-family:Georgia,serif;">Georgia</option>
                    <option value="Times New Roman, serif" ${currentFont.includes('Times') ? 'selected' : ''} style="font-family:'Times New Roman',serif;">Times New Roman</option>
                    <option value="Roboto Mono, monospace" ${currentFont.includes('Roboto Mono') ? 'selected' : ''} style="font-family:'Roboto Mono',monospace;">Roboto Mono</option>
                    <option value="Courier New, monospace" ${currentFont.includes('Courier') ? 'selected' : ''} style="font-family:'Courier New',monospace;">Courier New</option>
                </select>
            </div>
            <div style="border-top:1px solid #e5e7eb;margin:6px 0;"></div>
            <div style="font-weight:600;margin-bottom:4px;color:#1f2937;font-size:11px;">Text Content</div>
            ${textRow('Title', 'ts_titleText', titleText)}
            ${subtitleText ? textRow('Subtitle', 'ts_subtitleText', subtitleText, true) : ''}
            ${textRow('X Axis', 'ts_xLabelText', xLabel)}
            ${textRow('Y Axis', 'ts_yLabelText', yLabel)}
            <div style="border-top:1px solid #e5e7eb;margin:6px 0;"></div>
            <div style="font-weight:600;margin-bottom:4px;color:#1f2937;font-size:11px;">Font Sizes &amp; Visibility</div>
            ${sizeRow('Title', 'ts_title', titleSize, 6, 48, 'ts_titleVis', this._tsVisible.title)}
            ${subtitleText ? sizeRow('Subtitle', 'ts_subtitle', subtitleSize, 6, 30, null, true) : ''}
            ${sizeRow('X Label', 'ts_xlabel', xLabelSize, 6, 36, 'ts_xLabelVis', this._tsVisible.xLabel)}
            ${sizeRow('Y Label', 'ts_ylabel', yLabelSize, 6, 36, 'ts_yLabelVis', this._tsVisible.yLabel)}
            ${sizeRow('X Tick', 'ts_xtick', xTickSize, 6, 30, null, true)}
            ${sizeRow('Y Tick', 'ts_ytick', yTickSize, 6, 30, null, true)}
            ${sizeRow('Legend', 'ts_legend', legendSize, 6, 30, 'ts_legendVis', this._tsVisible.legend)}
            <div style="border-top:1px solid #e5e7eb;margin:6px 0;"></div>
            <div style="font-weight:600;margin-bottom:4px;color:#1f2937;font-size:11px;">Markers</div>
            ${sizeRow('Size', 'ts_marker', markerSize, 1, 40, null, true)}
            ${this._tsHasBoxTraces(plotEl) ? `
            <div style="border-top:1px solid #e5e7eb;margin:6px 0;"></div>
            <div style="font-weight:600;margin-bottom:4px;color:#1f2937;font-size:11px;">Box Color Scheme</div>
            <div style="display:flex;flex-wrap:wrap;gap:4px;">
                <button onclick="app._tsColorScheme('essential')" class="ts-color-btn" title="Red/Gray/Green by mean GE" style="font-size:10px;padding:3px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:linear-gradient(90deg,#dc2626 33%,#9ca3af 33%,#9ca3af 66%,#22c55e 66%);">Essential</button>
                <button onclick="app._tsColorScheme('bw')" class="ts-color-btn" title="Black & white" style="font-size:10px;padding:3px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:#e5e7eb;">B&W</button>
                <button onclick="app._tsColorScheme('blue')" class="ts-color-btn" title="Blue gradient by mean GE" style="font-size:10px;padding:3px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:linear-gradient(90deg,#1e40af,#93c5fd);">Blue</button>
                <button onclick="app._tsColorScheme('redblue')" class="ts-color-btn" title="Red-Blue diverging" style="font-size:10px;padding:3px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:linear-gradient(90deg,#dc2626,#f5f5f5 50%,#2563eb);">Red-Blue</button>
                <button onclick="app._tsColorScheme('viridis')" class="ts-color-btn" title="Viridis continuous" style="font-size:10px;padding:3px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:linear-gradient(90deg,#440154,#31688e,#35b779,#fde725);">Viridis</button>
                <button onclick="app._tsColorScheme('steelblue')" class="ts-color-btn" title="Uniform steelblue" style="font-size:10px;padding:3px 8px;border:1px solid #d1d5db;border-radius:4px;cursor:pointer;background:#4682b4;color:white;">Uniform</button>
            </div>
            ` : ''}
            <div style="border-top:1px solid #e5e7eb;margin:6px 0;"></div>
            <div style="font-size:10px;color:#9ca3af;">Drag title, axis labels, and annotations on plot to reposition.<br>Click an annotation, then use arrow keys to nudge (Shift = larger steps).</div>
        `;
        panel.style.display = 'block';

        // Initialize B/I button states from current formatting
        this._tsInitBoldItalic(plotEl, usesAnnotationTitle, ann0, xLabelAnn, yLabelAnn);

        // Setup arrow key listener for selected annotation
        this._tsSetupArrowKeys(plotEl);
    }

    openNetworkTextSettings() {
        const panel = document.getElementById('textSettingsPanel');
        const body = document.getElementById('textSettingsBody');
        this._textSettingsPlotId = '__network__';

        const sizeRow = (label, id, val, min, max) => `
            <div style="display:flex; align-items:center; margin-bottom:5px; gap:4px;">
                <span style="width:15px;"></span>
                <span style="color:#374151;flex:1;min-width:55px;font-size:11px;">${label}</span>
                <div style="display:flex; align-items:center;">
                    <button onclick="app._netTsStep('${id}',-1)" style="width:20px;height:20px;border:1px solid #d1d5db;background:#f9fafb;border-radius:4px 0 0 4px;cursor:pointer;font-size:12px;line-height:1;">−</button>
                    <input type="number" id="${id}" value="${val}" min="${min}" max="${max}" style="width:36px;text-align:center;border:1px solid #d1d5db;border-left:none;border-right:none;font-size:10px;padding:1px;-moz-appearance:textfield;appearance:textfield;" oninput="app._netTsApply()">
                    <button onclick="app._netTsStep('${id}',1)" style="width:20px;height:20px;border:1px solid #d1d5db;background:#f9fafb;border-radius:0 4px 4px 0;cursor:pointer;font-size:12px;line-height:1;">+</button>
                </div>
            </div>`;

        const colorRow = (label, id, val) => `
            <div style="display:flex; align-items:center; margin-bottom:5px; gap:4px;">
                <span style="width:15px;"></span>
                <span style="color:#374151;flex:1;min-width:55px;font-size:11px;">${label}</span>
                <input type="color" id="${id}" value="${val}" style="width:28px;height:22px;border:1px solid #d1d5db;border-radius:4px;padding:0;cursor:pointer;" oninput="app._netTsApply()">
            </div>`;

        // Current values
        const legendEl = document.getElementById('networkLegend');
        const currentFont = this._netFontFamily || 'Arial, sans-serif';
        const legendFontSize = this._netLegendFontSize || 11;
        const legendColor = this._netLegendColor || '#374151';
        const bannerFontSize = this._netBannerFontSize || 10;
        const bannerColor = this._netBannerColor || '#374151';
        const labelColor = this._netLabelColor || '#333333';
        const nodeColor = this._netNodeColor || '#5a9f4a';

        // Read current legend text from DOM
        const legendSections = legendEl?.querySelectorAll('.legend-section') || [];
        let legendTexts = [];
        legendSections.forEach(sec => {
            const strong = sec.querySelector('strong');
            if (strong) legendTexts.push(strong.textContent.replace(/:$/, ''));
        });

        // Read filter banner text
        const banner = document.querySelector('.network-filter-banner');
        const bannerText = banner?.textContent || '';

        body.innerHTML = `
            <div style="display:flex;align-items:center;gap:6px;margin-bottom:6px;">
                <span style="font-weight:600;color:#1f2937;font-size:11px;">Font</span>
                <select id="net_ts_font" onchange="app._netTsApply()" style="flex:1;font-size:11px;padding:2px 4px;border:1px solid #d1d5db;border-radius:4px;">
                    <option value="Arial, sans-serif" ${currentFont.includes('Arial') ? 'selected' : ''}>Arial</option>
                    <option value="Helvetica, Arial, sans-serif" ${currentFont.includes('Helvetica') && !currentFont.startsWith('Arial') ? 'selected' : ''}>Helvetica</option>
                    <option value="Open Sans, sans-serif" ${currentFont.includes('Open Sans') ? 'selected' : ''}>Open Sans</option>
                    <option value="Georgia, serif" ${currentFont.includes('Georgia') ? 'selected' : ''}>Georgia</option>
                    <option value="Times New Roman, serif" ${currentFont.includes('Times') ? 'selected' : ''}>Times New Roman</option>
                    <option value="Courier New, monospace" ${currentFont.includes('Courier') ? 'selected' : ''}>Courier New</option>
                </select>
            </div>
            <div style="border-top:1px solid #e5e7eb;margin:6px 0;"></div>
            <div style="font-weight:600;margin-bottom:4px;color:#1f2937;font-size:11px;">Legend</div>
            ${sizeRow('Font Size', 'net_ts_legendSize', legendFontSize, 6, 30)}
            ${colorRow('Text Color', 'net_ts_legendColor', legendColor)}
            <div style="border-top:1px solid #e5e7eb;margin:6px 0;"></div>
            <div style="font-weight:600;margin-bottom:4px;color:#1f2937;font-size:11px;">Filter Banner</div>
            ${sizeRow('Font Size', 'net_ts_bannerSize', bannerFontSize, 6, 24)}
            ${colorRow('Text Color', 'net_ts_bannerColor', bannerColor)}
            ${bannerText ? `<div style="margin-bottom:5px;">
                <label style="font-size:10px;color:#6b7280;">Banner Text</label>
                <input type="text" id="net_ts_bannerText" value="${this._escapeAttr(bannerText)}" style="width:100%;border:1px solid #d1d5db;border-radius:4px;padding:3px 6px;font-size:11px;margin-top:1px;box-sizing:border-box;" oninput="app._netTsApply()">
            </div>` : ''}
            <div style="border-top:1px solid #e5e7eb;margin:6px 0;"></div>
            <div style="font-weight:600;margin-bottom:4px;color:#1f2937;font-size:11px;">Colors</div>
            ${colorRow('Node Label', 'net_ts_labelColor', labelColor)}
            ${colorRow('Node Fill', 'net_ts_nodeColor', nodeColor)}
        `;

        panel.style.display = 'block';
    }

    _netTsStep(id, dir) {
        const el = document.getElementById(id);
        if (!el) return;
        el.value = Math.max(parseInt(el.min) || 1, Math.min(parseInt(el.max) || 99, (parseInt(el.value) || 0) + dir));
        this._netTsApply();
    }

    _netTsApply() {
        if (!this.network || !this.networkData) return;

        // Font family
        const fontFamily = document.getElementById('net_ts_font')?.value || 'Arial, sans-serif';
        this._netFontFamily = fontFamily;

        // Node label color
        const labelColor = document.getElementById('net_ts_labelColor')?.value || '#333333';
        const nodeColor = document.getElementById('net_ts_nodeColor')?.value || '#5a9f4a';
        this._netLabelColor = labelColor;
        this._netNodeColor = nodeColor;

        // Sync node font from slider (slider still controls node label size)
        const fontSize = parseInt(document.getElementById('netFontSize').value) || 16;

        // Update all nodes
        const nodeUpdates = [];
        this.networkData.nodes.forEach(node => {
            const update = {
                id: node.id,
                font: { size: fontSize, color: labelColor, face: fontFamily }
            };
            if (nodeColor !== '#5a9f4a' && node.color?.background !== '#d1d5db') {
                update.color = { ...node.color, background: nodeColor };
            }
            nodeUpdates.push(update);
        });
        this.networkData.nodes.update(nodeUpdates);

        // Legend styling
        const legendFontSize = parseInt(document.getElementById('net_ts_legendSize')?.value) || 11;
        const legendColor = document.getElementById('net_ts_legendColor')?.value || '#374151';
        this._netLegendFontSize = legendFontSize;
        this._netLegendColor = legendColor;

        const legendEl = document.getElementById('networkLegend');
        if (legendEl) {
            legendEl.style.fontSize = legendFontSize + 'px';
            legendEl.style.color = legendColor;
            legendEl.style.fontFamily = fontFamily;
            // Apply to all text within legend
            legendEl.querySelectorAll('strong, div, span').forEach(el => {
                el.style.fontSize = legendFontSize + 'px';
                el.style.color = legendColor;
                el.style.fontFamily = fontFamily;
            });
        }

        // Filter banner styling
        const bannerFontSize = parseInt(document.getElementById('net_ts_bannerSize')?.value) || 10;
        const bannerColor = document.getElementById('net_ts_bannerColor')?.value || '#374151';
        this._netBannerFontSize = bannerFontSize;
        this._netBannerColor = bannerColor;

        const banner = document.querySelector('.network-filter-banner');
        if (banner) {
            banner.style.fontSize = bannerFontSize + 'px';
            banner.style.color = bannerColor;
            banner.style.fontFamily = fontFamily;
            const bannerText = document.getElementById('net_ts_bannerText')?.value;
            if (bannerText !== undefined) banner.textContent = bannerText;
        }
    }

    _escapeAttr(s) {
        return s.replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
    }

    _tsInitBoldItalic(plotEl, usesAnn, titleAnn, xAnn, yAnn) {
        const setActive = (btnId, active) => {
            const btn = document.getElementById(btnId);
            if (btn) {
                btn.style.background = active ? '#dbeafe' : '#f9fafb';
                btn.style.borderColor = active ? '#3b82f6' : '#d1d5db';
            }
        };
        // Title bold/italic
        if (usesAnn && titleAnn) {
            const text = titleAnn.text || '';
            setActive('ts_titleText_bold', /<b>/i.test(text));
            setActive('ts_titleText_italic', /<i>/i.test(text));
        }
        // X/Y label bold/italic (check annotation font or layout font)
        const checkAnnStyle = (ann, prefix) => {
            if (ann) {
                const t = ann.text || '';
                setActive(`${prefix}_bold`, /<b>/i.test(t) || ann.font?.weight === 'bold');
                setActive(`${prefix}_italic`, /<i>/i.test(t) || ann.font?.style === 'italic');
            }
        };
        checkAnnStyle(xAnn, 'ts_xLabelText');
        checkAnnStyle(yAnn, 'ts_yLabelText');
    }

    _tsToggleBold(inputId) {
        const btn = document.getElementById(inputId + '_bold');
        const isActive = btn.style.borderColor === 'rgb(59, 130, 246)';
        btn.style.background = isActive ? '#f9fafb' : '#dbeafe';
        btn.style.borderColor = isActive ? '#d1d5db' : '#3b82f6';
        this._tsApplyText();
    }

    _tsToggleItalic(inputId) {
        const btn = document.getElementById(inputId + '_italic');
        const isActive = btn.style.borderColor === 'rgb(59, 130, 246)';
        btn.style.background = isActive ? '#f9fafb' : '#dbeafe';
        btn.style.borderColor = isActive ? '#d1d5db' : '#3b82f6';
        this._tsApplyText();
    }

    _tsBtnActive(btnId) {
        const btn = document.getElementById(btnId);
        return btn && btn.style.borderColor === 'rgb(59, 130, 246)';
    }

    _tsHasBoxTraces(plotEl) {
        return plotEl.data?.some(t => t.type === 'box');
    }

    _tsColorScheme(scheme) {
        const plotEl = document.getElementById(this._textSettingsPlotId);
        if (!plotEl?.data) return;

        // Get stats for mean-based schemes
        const stats = this.currentGEStats || [];
        const means = plotEl.data.map((t, i) => stats[i]?.mean ?? 0);
        const minMean = Math.min(...means);
        const maxMean = Math.max(...means);
        const range = Math.max(Math.abs(minMean), Math.abs(maxMean)) || 1;

        const markerColors = [];
        const fillColors = [];
        const lineColors = [];

        for (let i = 0; i < plotEl.data.length; i++) {
            if (plotEl.data[i].type !== 'box') continue;
            const mean = means[i];
            const t = (mean - minMean) / (maxMean - minMean || 1); // 0..1 normalized

            let mc, fc, lc;
            switch (scheme) {
                case 'essential':
                    mc = mean < -0.5 ? 'rgba(220,38,38,0.6)' : mean > 0 ? 'rgba(34,197,94,0.5)' : 'rgba(107,114,128,0.5)';
                    fc = mean < -0.5 ? 'rgba(220,38,38,0.2)' : mean > 0 ? 'rgba(34,197,94,0.15)' : 'rgba(156,163,175,0.2)';
                    lc = (stats[i]?.pValue ?? 1) < 0.05 ? '#1f2937' : '#9ca3af';
                    break;
                case 'bw':
                    mc = 'rgba(80,80,80,0.5)';
                    fc = 'rgba(200,200,200,0.3)';
                    lc = '#374151';
                    break;
                case 'blue': {
                    // Darker blue = more negative (essential)
                    const bInt = Math.round(255 - (1 - t) * 200);
                    mc = `rgba(30,64,${bInt},0.6)`;
                    fc = `rgba(30,64,${bInt},0.2)`;
                    lc = '#1e3a5f';
                    break;
                }
                case 'redblue': {
                    // Red for negative, blue for positive, white at zero
                    const frac = (mean + range) / (2 * range); // 0=most negative, 1=most positive
                    const r = Math.round(220 * (1 - frac));
                    const b = Math.round(220 * frac);
                    const g = Math.round(80 * (1 - Math.abs(frac - 0.5) * 2));
                    mc = `rgba(${r},${g},${b},0.6)`;
                    fc = `rgba(${r},${g},${b},0.2)`;
                    lc = '#374151';
                    break;
                }
                case 'viridis': {
                    // Approximate viridis: purple → teal → yellow
                    const vColors = [
                        [68,1,84], [59,82,139], [33,145,140], [94,201,98], [253,231,37]
                    ];
                    const idx = Math.min(Math.floor(t * 4), 3);
                    const frac = (t * 4) - idx;
                    const c0 = vColors[idx], c1 = vColors[idx + 1];
                    const r = Math.round(c0[0] + (c1[0] - c0[0]) * frac);
                    const g = Math.round(c0[1] + (c1[1] - c0[1]) * frac);
                    const b = Math.round(c0[2] + (c1[2] - c0[2]) * frac);
                    mc = `rgba(${r},${g},${b},0.7)`;
                    fc = `rgba(${r},${g},${b},0.25)`;
                    lc = `rgb(${r},${g},${b})`;
                    break;
                }
                case 'steelblue':
                    mc = 'rgba(70,130,180,0.6)';
                    fc = 'rgba(70,130,180,0.2)';
                    lc = '#2c5f8a';
                    break;
            }
            markerColors.push(mc);
            fillColors.push(fc);
            lineColors.push(lc);
        }

        // Apply via restyle
        const boxIndices = [];
        for (let i = 0; i < plotEl.data.length; i++) {
            if (plotEl.data[i].type === 'box') boxIndices.push(i);
        }
        for (let j = 0; j < boxIndices.length; j++) {
            Plotly.restyle(plotEl, {
                'marker.color': markerColors[j],
                fillcolor: fillColors[j],
                'line.color': lineColors[j]
            }, [boxIndices[j]]);
        }
    }

    _tsStep(inputId, direction) {
        const inp = document.getElementById(inputId);
        if (!inp) return;
        const cur = parseInt(inp.value) || 10;
        inp.value = Math.max(parseInt(inp.min) || 1, Math.min(parseInt(inp.max) || 48, cur + direction));
        this._tsApply();
    }

    _tsScaleAll(direction) {
        const ids = ['ts_title', 'ts_subtitle', 'ts_xlabel', 'ts_ylabel', 'ts_xtick', 'ts_ytick', 'ts_legend', 'ts_marker'];
        ids.forEach(id => {
            const inp = document.getElementById(id);
            if (!inp) return;
            const cur = parseInt(inp.value) || 10;
            inp.value = Math.max(parseInt(inp.min) || 1, Math.min(parseInt(inp.max) || 48, cur + direction));
        });
        this._tsApply();
    }

    _tsApplyFont() {
        const plotEl = document.getElementById(this._textSettingsPlotId);
        if (!plotEl?.layout) return;
        const family = document.getElementById('ts_fontFamily')?.value || 'Arial, Helvetica, sans-serif';

        const updates = {
            'font.family': family,
            'title.font.family': family,
            'xaxis.title.font.family': family,
            'yaxis.title.font.family': family,
            'xaxis.tickfont.family': family,
            'yaxis.tickfont.family': family,
            'legend.font.family': family
        };

        // Also apply to annotations
        const anns = plotEl.layout.annotations || [];
        for (let i = 0; i < anns.length; i++) {
            updates[`annotations[${i}].font.family`] = family;
        }

        Plotly.relayout(plotEl, updates);
    }

    _tsFindAnn(plotEl, role) {
        const anns = plotEl.layout?.annotations || [];
        return anns.findIndex(a => a._tsRole === role);
    }

    _tsToggle(checkboxId) {
        const plotEl = document.getElementById(this._textSettingsPlotId);
        if (!plotEl?.layout) return;
        const checked = document.getElementById(checkboxId)?.checked;

        if (checkboxId === 'ts_titleVis') {
            const idx = this._tsFindAnn(plotEl, 'title');
            if (idx >= 0) {
                Plotly.relayout(plotEl, { [`annotations[${idx}].visible`]: checked });
            } else if (this._tsOriginal.usesAnnotationTitle && plotEl.layout.annotations?.length > 0) {
                Plotly.relayout(plotEl, { 'annotations[0].visible': checked });
            } else {
                Plotly.relayout(plotEl, { 'title.text': checked ? (this._tsOriginal.titleText || ' ') : '' });
            }
        } else if (checkboxId === 'ts_xLabelVis') {
            const idx = this._tsFindAnn(plotEl, 'xlabel');
            if (idx >= 0) {
                Plotly.relayout(plotEl, { [`annotations[${idx}].visible`]: checked });
            } else {
                Plotly.relayout(plotEl, { 'xaxis.title.text': checked ? this._tsOriginal.xLabel : '' });
            }
        } else if (checkboxId === 'ts_yLabelVis') {
            const idx = this._tsFindAnn(plotEl, 'ylabel');
            if (idx >= 0) {
                Plotly.relayout(plotEl, { [`annotations[${idx}].visible`]: checked });
            } else {
                Plotly.relayout(plotEl, { 'yaxis.title.text': checked ? this._tsOriginal.yLabel : '' });
            }
        } else if (checkboxId === 'ts_legendVis') {
            Plotly.relayout(plotEl, { showlegend: checked });
        }
    }

    _tsApplyText() {
        const plotEl = document.getElementById(this._textSettingsPlotId);
        if (!plotEl?.layout) return;

        const titleText = document.getElementById('ts_titleText')?.value || '';
        const subtitleEl = document.getElementById('ts_subtitleText');
        const xLabel = document.getElementById('ts_xLabelText')?.value || '';
        const yLabel = document.getElementById('ts_yLabelText')?.value || '';

        const titleBold = this._tsBtnActive('ts_titleText_bold');
        const titleItalic = this._tsBtnActive('ts_titleText_italic');
        const subBold = this._tsBtnActive('ts_subtitleText_bold');
        const subItalic = this._tsBtnActive('ts_subtitleText_italic');
        const xBold = this._tsBtnActive('ts_xLabelText_bold');
        const xItalic = this._tsBtnActive('ts_xLabelText_italic');
        const yBold = this._tsBtnActive('ts_yLabelText_bold');
        const yItalic = this._tsBtnActive('ts_yLabelText_italic');

        const wrapBI = (text, bold, italic) => {
            let s = text;
            if (bold) s = `<b>${s}</b>`;
            if (italic) s = `<i>${s}</i>`;
            return s;
        };

        const updates = {};

        // Title — wrap in inline font-size span; annotation font.size controls line spacing
        const titleIdx = this._tsFindAnn(plotEl, 'title');
        const titleSizeVal = parseInt(document.getElementById('ts_title')?.value) || 14;
        if (titleIdx >= 0) {
            const subSize = parseInt(document.getElementById('ts_subtitle')?.value) || this._tsOriginal.subtitleSize || 10;
            let html = `<span style="font-size:${titleSizeVal}px">${wrapBI(titleText, titleBold, titleItalic)}</span>`;
            if (subtitleEl) {
                const lines = subtitleEl.value.split('\n').filter(l => l.trim());
                lines.forEach(line => { html += `<br><span style="font-size:${subSize}px;color:#666">${wrapBI(line, subBold, subItalic)}</span>`; });
            }
            updates[`annotations[${titleIdx}].text`] = html;
            updates[`annotations[${titleIdx}].font.size`] = Math.round(subSize * 0.85);
        } else if (this._tsOriginal.usesAnnotationTitle && plotEl.layout.annotations?.length > 0) {
            const subSize = parseInt(document.getElementById('ts_subtitle')?.value) || this._tsOriginal.subtitleSize || 10;
            let html = `<span style="font-size:${titleSizeVal}px">${wrapBI(titleText, titleBold, titleItalic)}</span>`;
            if (subtitleEl) {
                const lines = subtitleEl.value.split('\n').filter(l => l.trim());
                lines.forEach(line => { html += `<br><span style="font-size:${subSize}px;color:#666">${wrapBI(line, subBold, subItalic)}</span>`; });
            }
            updates['annotations[0].text'] = html;
            updates['annotations[0].font.size'] = Math.round(subSize * 0.85);
        } else {
            updates['title.text'] = wrapBI(titleText, titleBold, titleItalic);
        }

        // X/Y axis labels — annotation-based or native
        const xIdx = this._tsFindAnn(plotEl, 'xlabel');
        if (xIdx >= 0) {
            updates[`annotations[${xIdx}].text`] = wrapBI(xLabel, xBold, xItalic);
        } else {
            updates['xaxis.title.text'] = wrapBI(xLabel, xBold, xItalic);
        }
        const yIdx = this._tsFindAnn(plotEl, 'ylabel');
        if (yIdx >= 0) {
            updates[`annotations[${yIdx}].text`] = wrapBI(yLabel, yBold, yItalic);
        } else {
            updates['yaxis.title.text'] = wrapBI(yLabel, yBold, yItalic);
        }

        Plotly.relayout(plotEl, updates);
    }

    _tsApply() {
        const plotEl = document.getElementById(this._textSettingsPlotId);
        if (!plotEl?.layout) return;

        const getVal = (id) => { const v = parseInt(document.getElementById(id)?.value); return isNaN(v) ? null : v; };

        const xStandoff = getVal('ts_xStandoff');
        const yStandoff = getVal('ts_yStandoff');
        const updates = {
            'xaxis.tickfont.size': getVal('ts_xtick'),
            'yaxis.tickfont.size': getVal('ts_ytick'),
            'legend.font.size': getVal('ts_legend')
        };

        // Find annotation indices by role
        const anns = plotEl.layout.annotations || [];
        const titleIdx = anns.findIndex(a => a._tsRole === 'title');
        const xLabelIdx = anns.findIndex(a => a._tsRole === 'xlabel');
        const yLabelIdx = anns.findIndex(a => a._tsRole === 'ylabel');

        const titleSize = getVal('ts_title');
        const subtitleSize = getVal('ts_subtitle');
        if (titleIdx >= 0) {
            // Annotation font.size controls <br> line spacing — set to subtitle-based
            if (subtitleSize) updates[`annotations[${titleIdx}].font.size`] = Math.round(subtitleSize * 0.85);
            // Update inline font sizes: first match = title, rest = subtitle
            const raw = anns[titleIdx].text || '';
            let isFirst = true;
            const updatedText = raw.replace(/font-size:\s*\d+px/g, (match) => {
                if (isFirst) { isFirst = false; return `font-size:${titleSize}px`; }
                return subtitleSize ? `font-size:${subtitleSize}px` : match;
            });
            if (updatedText !== raw) updates[`annotations[${titleIdx}].text`] = updatedText;
        } else if (this._tsOriginal?.usesAnnotationTitle && anns.length > 0) {
            if (subtitleSize) updates['annotations[0].font.size'] = Math.round(subtitleSize * 0.85);
            const raw = anns[0].text || '';
            let isFirst = true;
            const updatedText = raw.replace(/font-size:\s*\d+px/g, (match) => {
                if (isFirst) { isFirst = false; return `font-size:${titleSize}px`; }
                return subtitleSize ? `font-size:${subtitleSize}px` : match;
            });
            if (updatedText !== raw) updates['annotations[0].text'] = updatedText;
        } else {
            updates['title.font.size'] = titleSize;
        }

        // Axis labels — annotation-based or native
        const xLabelSize = getVal('ts_xlabel');
        const yLabelSize = getVal('ts_ylabel');
        if (xLabelIdx >= 0) {
            if (xLabelSize) updates[`annotations[${xLabelIdx}].font.size`] = xLabelSize;
        } else {
            if (xLabelSize) updates['xaxis.title.font.size'] = xLabelSize;
        }
        if (yLabelIdx >= 0) {
            if (yLabelSize) updates[`annotations[${yLabelIdx}].font.size`] = yLabelSize;
        } else {
            if (yLabelSize) updates['yaxis.title.font.size'] = yLabelSize;
        }

        Plotly.relayout(plotEl, updates);

        const markerSize = getVal('ts_marker');
        if (markerSize && plotEl.data) {
            const indices = [];
            for (let i = 0; i < plotEl.data.length; i++) {
                if (plotEl.data[i]?.marker) indices.push(i);
            }
            if (indices.length > 0) Plotly.restyle(plotEl, { 'marker.size': markerSize }, indices);
        }

        // Persist settings so they survive plot re-renders
        if (this._textSettingsPlotId === 'scatterPlot') {
            this._savedScatterTextSettings = {
                titleFontSize: titleSize,
                subtitleSize: subtitleSize || this._tsOriginal?.subtitleSize || 15,
                xLabelFontSize: getVal('ts_xlabel') || 20,
                yLabelFontSize: getVal('ts_ylabel') || 20,
                xTickSize: getVal('ts_xtick') || 17,
                yTickSize: getVal('ts_ytick') || 17,
                legendSize: getVal('ts_legend') || 17,
                markerSize: markerSize || 10,
                fontFamily: plotEl.layout?.font?.family || null
            };
        }
    }

    _tsSetupArrowKeys(plotEl) {
        // Allow arrow keys to nudge annotations when plot is focused
        if (this._tsArrowHandler) {
            document.removeEventListener('keydown', this._tsArrowHandler);
        }
        this._tsSelectedAnnotation = null;

        // Listen for annotation clicks to select one
        plotEl.removeAllListeners?.('plotly_clickannotation');
        plotEl.on('plotly_clickannotation', (ev) => {
            this._tsSelectedAnnotation = ev.index;
            // Brief visual feedback
            const el = plotEl.querySelector('.annotation[data-index="' + ev.index + '"]');
            if (el) { el.style.outline = '2px solid #3b82f6'; setTimeout(() => el.style.outline = '', 1500); }
        });

        this._tsArrowHandler = (e) => {
            if (this._tsSelectedAnnotation == null) return;
            const idx = this._tsSelectedAnnotation;
            const ann = plotEl.layout.annotations?.[idx];
            if (!ann || ann.xref !== 'paper') return;

            const step = e.shiftKey ? 0.01 : 0.005;
            let dx = 0, dy = 0;
            if (e.key === 'ArrowLeft') dx = -step;
            else if (e.key === 'ArrowRight') dx = step;
            else if (e.key === 'ArrowUp') dy = step;
            else if (e.key === 'ArrowDown') dy = -step;
            else return;

            e.preventDefault();
            const upd = {};
            upd[`annotations[${idx}].x`] = (ann.x || 0.5) + dx;
            upd[`annotations[${idx}].y`] = (ann.y || 1.0) + dy;
            Plotly.relayout(plotEl, upd);
        };
        document.addEventListener('keydown', this._tsArrowHandler);
    }

    _initTextSettingsDrag() {
        const panel = document.getElementById('textSettingsPanel');
        const handle = document.getElementById('textSettingsDragHandle');
        if (!panel || !handle) return;
        let dragging = false, ox, oy;
        handle.addEventListener('mousedown', (e) => {
            dragging = true;
            ox = e.clientX - panel.getBoundingClientRect().left;
            oy = e.clientY - panel.getBoundingClientRect().top;
            e.preventDefault();
        });
        document.addEventListener('mousemove', (e) => {
            if (!dragging) return;
            panel.style.left = (e.clientX - ox) + 'px';
            panel.style.top = (e.clientY - oy) + 'px';
            panel.style.right = 'auto';
        });
        document.addEventListener('mouseup', () => { dragging = false; });
    }
}

// Initialize app
const app = new CorrelationExplorer();
