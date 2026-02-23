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
        this.savedNetworkView = null;

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
                if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cl] !== subLineageFilter) {
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
                if (subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cl] !== subLineageFilter) {
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

    showTissueBreakdownPopup() {
        this.hideTissueBreakdownPopup();
        const gene = document.getElementById('mutationHotspotSelect').value;
        if (!gene) return;

        const breakdown = this.getTissueBreakdownForHotspot(gene);
        if (breakdown.length === 0) return;

        const currentLineage = document.getElementById('lineageFilter').value;

        const popup = document.createElement('div');
        popup.id = 'tissueBreakdownPopup';
        popup.style.cssText = 'position: fixed; z-index: 10000; background: white; border: 1px solid #d1d5db; border-radius: 8px; box-shadow: 0 8px 24px rgba(0,0,0,0.15); padding: 0; min-width: 340px; max-width: 420px; display: flex; flex-direction: column;';

        // Position near the button, clamped to viewport
        const btn = document.getElementById('tissueBreakdownBtn');
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
                <th style="padding: 4px 6px; text-align: right; color: #dc2626;">Mut</th>
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

            html += `<tr class="tb-row" data-tissue="${t.lineage}" style="cursor: pointer; ${rowBg}" onmouseenter="this.style.background=this.style.background||'#f3f4f6'" onmouseleave="this.style.background='${isCurrentFilter ? '#ecfdf5' : ''}'">
                <td style="padding: 3px 8px;"><input type="checkbox" class="tb-check" value="${t.lineage}"${checked}></td>
                <td style="padding: 3px 4px; font-weight: ${t.nMut > 0 ? '500' : '400'}; color: ${t.nMut > 0 ? '#1f2937' : '#9ca3af'};">${t.lineage}</td>
                <td style="padding: 3px 6px; text-align: right; color: #dc2626; font-weight: 600;">${t.nMut}</td>
                <td style="padding: 3px 6px; text-align: right; color: #6b7280;">${t.nWT}</td>
                <td style="padding: 3px 8px;"><div style="background: #fee2e2; border-radius: 2px; height: 10px; width: 100%;"><div style="background: #ef4444; border-radius: 2px; height: 10px; width: ${barWidth}%;"></div></div></td>
                <td style="padding: 3px 8px; text-align: right; color: #6b7280; font-size: 11px;">${pct}%</td>
            </tr>`;
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

        // Row click toggles checkbox
        popup.querySelectorAll('.tb-row').forEach(row => {
            row.addEventListener('click', (e) => {
                if (e.target.type === 'checkbox') return;
                const cb = row.querySelector('.tb-check');
                cb.checked = !cb.checked;
                this.updateTBSelectionCount();
            });
        });

        // Checkbox change
        popup.querySelectorAll('.tb-check').forEach(cb => {
            cb.addEventListener('change', () => this.updateTBSelectionCount());
        });

        // Select all
        document.getElementById('tbSelectAll').addEventListener('change', (e) => {
            popup.querySelectorAll('.tb-check').forEach(cb => { cb.checked = e.target.checked; });
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
            const selected = [...popup.querySelectorAll('.tb-check:checked')].map(cb => cb.value);
            this.applyTissueBreakdownSelection(selected);
            this.hideTissueBreakdownPopup();
        });

        // Escape key
        this._tbEscHandler = (e) => { if (e.key === 'Escape') this.hideTissueBreakdownPopup(); };
        document.addEventListener('keydown', this._tbEscHandler);

        // Click outside
        setTimeout(() => {
            this._tbOutsideHandler = (e) => {
                const p = document.getElementById('tissueBreakdownPopup');
                if (p && !p.contains(e.target) && e.target.id !== 'tissueBreakdownBtn') {
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

    updateTBSelectionCount() {
        const popup = document.getElementById('tissueBreakdownPopup');
        if (!popup) return;
        const count = popup.querySelectorAll('.tb-check:checked').length;
        const label = document.getElementById('tbSelectionCount');
        if (label) label.textContent = count === 0 ? '0 selected' : `${count} selected`;
    }

    applyTissueBreakdownSelection(selectedTissues) {
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

        // Trigger UI updates
        this.updateSubLineageFilter();
        this.updateHotspotCountsForCurrentFilters();
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
        // Tab switching - preserve network view state
        document.querySelectorAll('.nav-link').forEach(tab => {
            tab.addEventListener('click', () => {
                // Save network view state before switching away
                if (this.network) {
                    this.savedNetworkView = {
                        scale: this.network.getScale(),
                        position: this.network.getViewPosition()
                    };
                }

                document.querySelectorAll('.nav-link').forEach(t => t.classList.remove('active'));
                document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
                tab.classList.add('active');
                document.getElementById('tab-' + tab.dataset.tab).classList.add('active');

                // Restore network view state when switching back to network tab
                if (tab.dataset.tab === 'network' && this.network && this.savedNetworkView) {
                    const container = document.getElementById('networkPlot');
                    const savedView = this.savedNetworkView;

                    const restoreView = () => {
                        if (savedView && this.network) {
                            this.network.moveTo({
                                position: savedView.position,
                                scale: savedView.scale,
                                animation: false
                            });
                        }
                    };

                    // Wait for container to have valid dimensions
                    const waitForDimensions = () => {
                        const width = container.clientWidth;
                        const height = container.clientHeight;
                        if (width > 0 && height > 0) {
                            // Tell vis.js to use the correct size
                            this.network.setSize(width + 'px', height + 'px');
                            // Restore view multiple times to ensure it sticks
                            restoreView();
                            setTimeout(restoreView, 100);
                            setTimeout(restoreView, 300);
                        } else {
                            // Container not ready yet, try again
                            requestAnimationFrame(waitForDimensions);
                        }
                    };

                    requestAnimationFrame(waitForDimensions);
                }
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
            const testGenes = ['TP53', 'BRCA1', 'BRCA2', 'MYC', 'KRAS', 'EGFR', 'PTEN',
                'RB1', 'APC', 'CDKN2A', 'NOTCH1', 'PIK3CA', 'BRAF',
                'ATM', 'ERBB2', 'CDK4', 'MDM2', 'NRAS', 'TSC1', 'TSC2'];
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

        // Analysis mode change
        document.querySelectorAll('input[name="analysisMode"]').forEach(radio => {
            radio.addEventListener('change', () => this.updateAnalysisModeUI());
        });

        // Tissue breakdown button
        document.getElementById('tissueBreakdownBtn').addEventListener('click', () => this.showTissueBreakdownPopup());
        document.getElementById('mutationHotspotSelect').addEventListener('change', () => {
            const btn = document.getElementById('tissueBreakdownBtn');
            btn.style.display = document.getElementById('mutationHotspotSelect').value ? 'inline-block' : 'none';
        });

        // Mutation results search
        document.getElementById('mutationSearch').addEventListener('input', (e) => {
            this.filterMutationTable(e.target.value);
        });

        // Download mutation results
        document.getElementById('downloadMutationResults').addEventListener('click', () => {
            this.downloadMutationResults();
        });

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
        });
        document.getElementById('geneEffectModal').addEventListener('click', (e) => {
            if (e.target.id === 'geneEffectModal') {
                document.getElementById('geneEffectModal').style.display = 'none';
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
            if (this.results) this.displayNetwork();
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
                // Close inspect modal if open
                const inspectModal = document.getElementById('inspectModal');
                if (inspectModal && inspectModal.classList.contains('active')) {
                    this.closeInspectModal();
                    return;
                }
                // Close gene effect modal if open
                const geneEffectModal = document.getElementById('geneEffectModal');
                if (geneEffectModal && geneEffectModal.style.display !== 'none') {
                    geneEffectModal.style.display = 'none';
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
            this.updateInspectPlot();
        });

        document.getElementById('resetAllFilters').addEventListener('click', () => this.resetAllInspectFilters());

        ['scatterXmin', 'scatterXmax', 'scatterYmin', 'scatterYmax'].forEach(id => {
            const el = document.getElementById(id);
            el.addEventListener('change', () => this.updateInspectPlot());
            el.addEventListener('keydown', (e) => { if (e.key === 'Enter') this.updateInspectPlot(); });
        });

        document.getElementById('scatterCellSearch').addEventListener('input', () => this.updateInspectPlot());
        document.getElementById('scatterCancerFilter').addEventListener('change', () => {
            this.updateScatterSubtypeFilter();
            this.updateInspectPlot();
        });
        document.getElementById('scatterSubtypeFilter').addEventListener('change', () => this.updateInspectPlot());
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
        document.getElementById('compareAllCancerTypesBtn')?.addEventListener('click', () => this.showCompareAllCancerTypes());
        document.getElementById('updateInspectGenes')?.addEventListener('click', () => this.updateInspectGenes());

        // Aspect ratio control
        document.getElementById('aspectRatio')?.addEventListener('input', (e) => {
            document.getElementById('aspectRatioValue').textContent = parseFloat(e.target.value).toFixed(1);
            this.updateInspectPlot();
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

        // Mutation analysis compare buttons (#16)
        document.getElementById('mutCompareByTissueBtn')?.addEventListener('click', () => this.showMutationCompareByTissue());
        document.getElementById('mutCompareByHotspotBtn')?.addEventListener('click', () => this.showMutationCompareByHotspot());
        document.getElementById('mutCompareCloseBtn')?.addEventListener('click', () => {
            document.getElementById('mutationCompareTable').style.display = 'none';
            document.getElementById('mutCompareCloseBtn').style.display = 'none';
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
            if (this.geneEffectViewMode === 'mutation' && this.mutationResults && this.currentGeneEffectGene) {
                this.showGeneEffectDistribution(this.currentGeneEffectGene);
            } else {
                this.switchGeneEffectView(this.currentGEView || 'tissue');
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

        // Aspect ratio control for detailed views
        document.getElementById('geAspectRatio')?.addEventListener('input', (e) => {
            const ratio = parseFloat(e.target.value);
            document.getElementById('geAspectRatioValue').textContent = ratio.toFixed(1);
            this.geChartWidthRatio = ratio;
            // Re-render the detailed view with new width
            if (this.geDetailedView) {
                this.showGEDetailedView(this.geDetailedView.group, this.geDetailedView.mode);
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
        // Test data with LFC and FDR values (20 genes)
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

        // Get mutation data for hotspot filter
        let mutationData = null;
        if (hotspotGene && this.mutations?.geneData?.[hotspotGene]) {
            mutationData = this.mutations.geneData[hotspotGene].mutations;
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

            indices.push(idx);
        });

        // Return all indices if no filters applied
        if (indices.length === 0 && !lineageFilter && !subLineageFilter && (!hotspotGene || hotspotLevel === 'all')) {
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
        const hotspotGene = document.getElementById('mutationHotspotSelect').value;
        const minN = parseInt(document.getElementById('minCellLines').value);
        const pThreshold = this.getInputNum('pValueThreshold');
        const lineageFilter = document.getElementById('lineageFilter').value;
        const subLineageFilter = document.getElementById('subLineageFilter')?.value;

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
                const analysisResult = this.calculateMutationAnalysis(hotspotGene, minN, lineageFilter, subLineageFilter, additionalHotspot, additionalHotspotLevel);

                // Filter by p-value threshold
                const significantResults = analysisResult.results.filter(r => r.p_mut < pThreshold || r.p_2 < pThreshold);

                // Sort by p-value (1+2 vs 0)
                significantResults.sort((a, b) => a.p_mut - b.p_mut);

                this.mutationResults = {
                    hotspotGene,
                    pThreshold,
                    minN,
                    lineageFilter,
                    subLineageFilter,
                    additionalHotspot,
                    additionalHotspotLevel,
                    excludedTissues: new Set(this.excludedTissues),
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

    calculateMutationAnalysis(hotspotGene, minN, lineageFilter, subLineageFilter, additionalHotspot, additionalHotspotLevel) {
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
                <td class="gene-hover" data-gene="${r.gene}"><a href="#" style="color: var(--green-700); text-decoration: none; cursor: pointer;" onclick="app.openGeneEffectModal('${r.gene}', 'tissue'); return false;" onmouseover="this.style.textDecoration='underline'" onmouseout="this.style.textDecoration='none'">${r.gene}</a></td>
                <td style="border-left: 2px solid #2563eb;">${r.n_wt}</td>
                <td>${r.mean_wt.toFixed(3)}</td>
                <td style="border-left: 2px solid #f97316;">${r.n_mut}</td>
                <td>${r.mean_mut.toFixed(3)}</td>
                <td class="${r.diff_mut < 0 ? 'negative' : 'positive'}">${r.diff_mut.toFixed(3)}</td>
                <td>${this.formatPValue(r.p_mut)}</td>
                <td style="border-left: 2px solid #dc2626;">${r.n_2}</td>
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
            let lineageText = mr.lineageFilter;
            if (mr.subLineageFilter) {
                lineageText += ` (${mr.subLineageFilter})`;
            }
            settingsText += ` | Lineage: ${lineageText}`;
        }
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            settingsText += ` | Filter: ${mr.additionalHotspot} ${mr.additionalHotspotLevel}`;
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

        // Attach gene tooltips
        this.attachGeneTooltips(tbody);
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

    showGeneEffectDistribution(gene, tissueOverride) {
        if (!this.mutationResults) return;

        const mr = this.mutationResults;
        const hotspotGene = mr.hotspotGene;
        const mutationData = this.mutations.geneData[hotspotGene];
        const geneIdx = this.geneIndex.get(gene.toUpperCase());

        if (geneIdx === undefined) {
            alert(`Gene ${gene} not found`);
            return;
        }

        // Get tissue filter - use override if provided, otherwise read from dropdown
        const inspectTissueFilter = tissueOverride !== undefined ? tissueOverride : (document.getElementById('geTissueFilter')?.value || '');

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

            // Check sublineage filter
            if (!inspectTissueFilter && mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) {
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
            let lineageText = mr.lineageFilter;
            if (mr.subLineageFilter) {
                lineageText += ` (${mr.subLineageFilter})`;
            }
            filterInfo.push(`Lineage: ${lineageText}`);
        }
        if (inspectTissueFilter) {
            filterInfo.push(`Tissue: ${inspectTissueFilter}`);
        }
        if (mr.excludedTissues && mr.excludedTissues.size > 0 && !inspectTissueFilter) {
            const allLineages = this.cellLineMetadata?.lineage
                ? [...new Set(Object.values(this.cellLineMetadata.lineage))].sort()
                : [];
            const included = allLineages.filter(t => !mr.excludedTissues.has(t));
            filterInfo.push(`Tissues: ${included.join(', ')}`);
        }
        if (mr.additionalHotspot && mr.additionalHotspotLevel !== 'all') {
            filterInfo.push(`${mr.additionalHotspot}: ${mr.additionalHotspotLevel}`);
        }
        const lineageText = filterInfo.length > 0 ? filterInfo.join(' | ') : 'All lineages';

        // Build stats text for subtitle
        const formatP = (p) => isNaN(p) ? '-' : (p < 0.001 ? p.toExponential(1) : p.toFixed(3));
        let statsLine1 = `WT: n=${wtStats.n}, mean=${wtStats.mean.toFixed(2)}, med=${wtStats.median.toFixed(2)}`;
        statsLine1 += `  ·  Mut: n=${mutAllStats.n}, mean=${mutAllStats.mean.toFixed(2)}, med=${mutAllStats.median.toFixed(2)}`;
        let statsLine2 = `p(WT vs Mut): ${formatP(pWTvsMut)}`;
        if (mut2Stats.n >= 3) {
            statsLine2 += `  ·  p(WT vs 2): ${formatP(pWTvs2)}`;
        }

        // Combine lineage info and stats in subtitle
        const subtitle = `${lineageText}<br>${statsLine1}<br>${statsLine2}`;

        const layout = {
            title: {
                text: `${gene} Gene Effect by ${hotspotGene} Mutation Status<br><sub style="font-size:11px;color:#666">${subtitle}</sub>`,
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
                range: [-0.5, 2.5],
                automargin: true
            },
            showlegend: false,
            margin: { t: 100, r: 30, b: 50, l: 30 }
        };

        // Show modal
        document.getElementById('geneEffectModal').style.display = 'flex';
        document.getElementById('geneEffectTitle').textContent = `${gene} Gene Effect by ${hotspotGene} Mutation`;

        // Populate tissue filter dropdown with available lineages (respecting excluded tissues)
        const tissueFilterEl = document.getElementById('geTissueFilter');
        if (tissueFilterEl) {
            const allLineages = [...new Set(cellLines.map(cl => this.cellLineMetadata?.lineage?.[cl]).filter(Boolean))].sort();
            const visibleLineages = mr.excludedTissues && mr.excludedTissues.size > 0
                ? allLineages.filter(l => !mr.excludedTissues.has(l))
                : allLineages;
            const defaultLabel = mr.excludedTissues && mr.excludedTissues.size > 0
                ? `Filtered tissues (${visibleLineages.length})`
                : 'All tissues';
            tissueFilterEl.innerHTML = `<option value="">${defaultLabel}</option>`;
            visibleLineages.forEach(l => {
                const opt = document.createElement('option');
                opt.value = l;
                opt.textContent = l;
                tissueFilterEl.appendChild(opt);
            });
            tissueFilterEl.value = inspectTissueFilter;
        }

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

        Plotly.newPlot('geneEffectPlot', traces, layout, { responsive: true });
    }

    _exportMutationInspectChart(format) {
        if (this.geneEffectViewMode !== 'mutation') return;
        const plotEl = document.getElementById('geneEffectPlot');
        if (!plotEl || !plotEl.data) return;
        const exportWidth = 1200;
        const exportHeight = 600;
        const filename = `gene_effect_${this.currentGeneEffectGene}_${this.mutationResults.hotspotGene}`;

        // Deep-copy data and layout from live chart
        const data = JSON.parse(JSON.stringify(plotEl.data));
        const layout = JSON.parse(JSON.stringify(plotEl.layout));

        // Override layout for export
        layout.margin = Object.assign({}, layout.margin, { l: 200 });
        layout.width = exportWidth;
        layout.height = exportHeight;
        // Add standoff to push y-axis title away from tick labels
        if (layout.yaxis) {
            const titleText = typeof layout.yaxis.title === 'string' ? layout.yaxis.title : layout.yaxis.title?.text || '';
            layout.yaxis.title = { text: titleText, standoff: 40 };
            layout.yaxis.automargin = true;
        }

        // Render into a temporary off-screen div for a clean export
        const tempDiv = document.createElement('div');
        tempDiv.style.position = 'absolute';
        tempDiv.style.left = '-10000px';
        tempDiv.style.top = '0';
        tempDiv.style.width = exportWidth + 'px';
        tempDiv.style.height = exportHeight + 'px';
        document.body.appendChild(tempDiv);

        Plotly.newPlot(tempDiv, data, layout, { staticPlot: true }).then(() => {
            return Plotly.downloadImage(tempDiv, {
                format,
                width: exportWidth,
                height: exportHeight,
                filename
            });
        }).then(() => {
            Plotly.purge(tempDiv);
            document.body.removeChild(tempDiv);
        }).catch(() => {
            try { Plotly.purge(tempDiv); document.body.removeChild(tempDiv); } catch(e) {}
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
            csv += `${d.cellLine},${d.cellName},${d.lineage},${d.ge.toFixed(4)},${d.mutLevel}\n`;
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

        // Display all results
        this.displayNetwork();
        this.displayCorrelationsTable();
        this.displayClustersTable();
        this.displaySummary();

        // Switch to network tab
        document.querySelectorAll('.nav-link').forEach(link => link.classList.remove('active'));
        document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
        document.querySelector('[data-tab="network"]').classList.add('active');
        document.getElementById('tab-network').classList.add('active');
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
                    border: '#ffffff'
                },
                borderWidth: 2,
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
                        borderWidth: 2,
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
                                border: '#ffffff'
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
                    <td class="gene-hover" data-gene="${c.gene1}" style="cursor: help;">${c.gene1}</td>
                    <td class="gene-hover" data-gene="${c.gene2}" style="cursor: help;">${c.gene2}</td>
                    <td class="${corrClass}">${c.correlation.toFixed(3)}</td>
                    <td>${c.slope.toFixed(3)}</td>
                    <td>${c.n}</td>
                    <td>${c.cluster}</td>
                    <td style="white-space: nowrap;">
                        <button class="btn btn-sm inspect-btn" style="padding: 2px 6px; font-size: 10px; background: #5a9f4a; color: white;" data-gene1="${c.gene1}" data-gene2="${c.gene2}">Correlate</button>
                        <button class="btn btn-sm tissue-btn" style="padding: 2px 6px; font-size: 10px; margin-left: 4px; background: #6b7280; color: white;" data-gene1="${c.gene1}" data-gene2="${c.gene2}">By Tissue</button>
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

        // Attach gene tooltip handlers
        this.attachGeneTooltips(tbody);
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
                        <button class="btn btn-sm tissue-btn" style="padding: 2px 6px; font-size: 10px; background: #5a9f4a; color: white;" data-gene="${c.gene}">By Tissue</button>
                        <button class="btn btn-sm hotspot-btn" style="padding: 2px 6px; font-size: 10px; margin-left: 4px; background: #6b7280; color: white;" data-gene="${c.gene}">By Hotspot</button>
                    </td>
                `;

                tr.innerHTML = rowHtml;
                tbody.appendChild(tr);
            });

        // Add event listeners to buttons
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

        text.textContent = `Gene Correlation Analysis Summary
================================
Run: ${dateTimeStr}

Analysis Mode: ${this.results.mode === 'analysis' ? 'Analysis (within gene list)' : 'Design (find correlated genes)'}
Correlation Cutoff: ${this.results.cutoff}
Minimum Cell Lines: ${document.getElementById('minCellLines').value}
Minimum Slope: ${document.getElementById('minSlope').value}
Lineage Filter: ${lineageText}

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

        // Get positions from vis.js and convert to DOM coordinates
        const positions = this.network.getPositions();
        const domPositions = {};
        for (const nodeId in positions) {
            const canvasPos = positions[nodeId];
            const domPos = this.network.canvasToDOM({ x: canvasPos.x, y: canvasPos.y });
            domPositions[nodeId] = domPos;
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
`;

        // Get current scale for sizing elements
        const scale = this.network.getScale();

        // Draw edges
        this.networkData.edges.forEach(edge => {
            const from = domPositions[edge.from];
            const to = domPositions[edge.to];
            if (from && to) {
                const color = edge.color || '#3182ce';
                const strokeWidth = (edge.width || 2) * scale;
                svg += `  <line x1="${from.x}" y1="${from.y}" x2="${to.x}" y2="${to.y}" stroke="${color}" stroke-width="${strokeWidth}" opacity="0.8"/>\n`;
            }
        });

        // Draw nodes
        const nodeSize = (parseInt(document.getElementById('netNodeSize').value) || 25) * scale;
        const fontSize = (parseInt(document.getElementById('netFontSize').value) || 16) * scale;
        this.networkData.nodes.forEach(node => {
            const pos = domPositions[node.id];
            if (pos) {
                const bgColor = node.color?.background || '#5a9f4a';
                svg += `  <circle cx="${pos.x}" cy="${pos.y}" r="${nodeSize/2}" fill="${bgColor}" stroke="white" stroke-width="${2 * scale}"/>\n`;

                // Handle multi-line labels
                const labelLines = (node.label || node.id).split('\n');
                labelLines.forEach((line, i) => {
                    const yOffset = pos.y + nodeSize/2 + 14 * scale + (i * (fontSize - 2 * scale));
                    svg += `  <text x="${pos.x}" y="${yOffset}" text-anchor="middle" style="font-family: Arial; font-size: ${fontSize - 2 * scale}px; fill: #333;">${this.escapeXml(line)}</text>\n`;
                });
            }
        });

        // Draw legend - LARGER for publication
        const legendY = networkHeight + 35;

        // Calculate total legend width to center it
        let totalLegendWidth = 160 + 160; // Correlation + Edge Thickness (note: 110+160 is used but let's use similar calc)
        if (this.results?.mode === 'design') totalLegendWidth += 140;
        if (document.getElementById('colorByGeneEffect').checked && this.results?.clusters) totalLegendWidth += 170;
        if (document.getElementById('colorByStats').checked && this.geneStats && this.geneStats.size > 0) totalLegendWidth += 200;

        let legendX = Math.max(40, (width - totalLegendWidth) / 2);

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
                        color: { background: bgColor, border: '#ffffff' }
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
                        color: { background: bgColor, border: '#ffffff' }
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
                        color: { background: bgColor, border: '#ffffff' }
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
                    border: '#ffffff'
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

        // Cycle through layout options: Barnes Hut (default), Force Atlas, Hierarchical
        const layouts = ['barnesHut', 'forceAtlas2Based', 'hierarchical'];
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
            btn.textContent = 'Lock';
            btn.classList.remove('btn-active');
        } else {
            btn.textContent = 'Unlock';
            btn.classList.add('btn-active');
        }

        // Show which layout is active
        const layoutBtn = document.getElementById('changeLayout');
        const layoutNames = { 'barnesHut': 'Barnes Hut', 'forceAtlas2Based': 'Force Atlas', 'hierarchical': 'Hierarchical' };
        layoutBtn.textContent = layoutNames[layoutName];
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
        const totalHeight = networkHeight + legendHeight;

        // Get positions from vis.js and convert to DOM coordinates
        const positions = this.network.getPositions();
        const domPositions = {};
        for (const nodeId in positions) {
            const canvasPos = positions[nodeId];
            const domPos = this.network.canvasToDOM({ x: canvasPos.x, y: canvasPos.y });
            domPositions[nodeId] = domPos;
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
`;

        // Get current scale for sizing elements
        const scale = this.network.getScale();

        // Draw edges
        this.networkData.edges.forEach(edge => {
            const from = domPositions[edge.from];
            const to = domPositions[edge.to];
            if (from && to) {
                const color = edge.color || '#3182ce';
                const strokeWidth = (edge.width || 2) * scale;
                svg += `  <line x1="${from.x}" y1="${from.y}" x2="${to.x}" y2="${to.y}" stroke="${color}" stroke-width="${strokeWidth}" opacity="0.8"/>\n`;
            }
        });

        // Draw nodes
        const nodeSize = (parseInt(document.getElementById('netNodeSize').value) || 25) * scale;
        const fontSize = (parseInt(document.getElementById('netFontSize').value) || 16) * scale;
        this.networkData.nodes.forEach(node => {
            const pos = domPositions[node.id];
            if (pos) {
                const bgColor = node.color?.background || '#5a9f4a';
                svg += `  <circle cx="${pos.x}" cy="${pos.y}" r="${nodeSize/2}" fill="${bgColor}" stroke="white" stroke-width="${2 * scale}"/>\n`;

                // Handle multi-line labels
                const labelLines = (node.label || node.id).split('\n');
                labelLines.forEach((line, i) => {
                    const yOffset = pos.y + nodeSize/2 + 14 * scale + (i * (fontSize - 2 * scale));
                    svg += `  <text x="${pos.x}" y="${yOffset}" text-anchor="middle" style="font-family: Arial; font-size: ${fontSize - 2 * scale}px; fill: #333;">${this.escapeXml(line)}</text>\n`;
                });
            }
        });

        // Draw legend - LARGER for publication
        const legendY = networkHeight + 35;

        // Calculate total legend width to center it
        let totalLegendWidth = 160 + 160; // Correlation + Edge Thickness
        if (this.results?.mode === 'design') totalLegendWidth += 140;
        if (document.getElementById('colorByGeneEffect').checked && this.results?.clusters) totalLegendWidth += 170;
        if (document.getElementById('colorByStats').checked && this.geneStats && this.geneStats.size > 0) totalLegendWidth += 200;

        let legendX = Math.max(40, (width - totalLegendWidth) / 2);

        // Legend background
        svg += `  <rect x="15" y="${networkHeight + 10}" width="${width - 30}" height="145" fill="#f9fafb" stroke="#e5e7eb" rx="4"/>\n`;

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
            cancerFilter.innerHTML = `<option value="">All cancer types (n=${plotData.length})</option>`;
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

        // Calculate stats for ALL cells (unfiltered) for the title
        const allCellsStats = this.pearsonWithSlope(plotData.map(d => d.x), plotData.map(d => d.y));
        document.getElementById('inspectTitle').textContent =
            `${c.gene1} vs ${c.gene2} | r=${this.formatNum(allCellsStats.correlation)}, slope=${this.formatNum(allCellsStats.slope)}, n=${plotData.length} (all cells)`;

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

        addRegressionLine(wt, wtStats, 'x', 'y', '#5a9f4a');
        addRegressionLine(mut1, mut1Stats, 'x2', 'y2', '#5a9f4a');
        addRegressionLine(mut2, mut2Stats, 'x3', 'y3', '#5a9f4a');

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
            <p style="font-size: 11px; color: #666; margin-bottom: 8px;">
                Comparing correlation between WT (0 ${hotspotGene} mutations) vs Mutant (2+ ${hotspotGene} mutations) cells, stratified by cancer type.
                Note: Cells with exactly 1 mutation are excluded from this comparison.
                <strong>Click a cancer type</strong> to view its scatter plot with the ${hotspotGene} mutation overlay.
            </p>
            <p style="font-size: 10px; color: #0c4a6e; background: #f0f9ff; padding: 4px 8px; border-radius: 4px; margin-bottom: 12px;">
                <b>Statistics:</b> p(Δr) uses Fisher z-transformation to compare correlations. p(Δslope) is an approximation based on correlation difference.
            </p>
            <div style="overflow-x: auto;">
            <table id="compareByCancerTable" class="data-table" style="width: 100%; font-size: 12px;">
                <thead>
                    <tr>
                        <th data-col="0" style="cursor: pointer;">Cancer Type ▼</th>
                        <th data-col="1" style="cursor: pointer; border-left: 2px solid #2563eb;">N (WT)</th>
                        <th data-col="2" style="cursor: pointer;">r (WT)</th>
                        <th data-col="3" style="cursor: pointer;">slope (WT)</th>
                        <th data-col="4" style="cursor: pointer; border-left: 2px solid #dc2626;">N (Mut)</th>
                        <th data-col="5" style="cursor: pointer;">r (Mut)</th>
                        <th data-col="6" style="cursor: pointer;">slope (Mut)</th>
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
            csv += `# Hotspot filter: ${hotspotGene}\n`;
            csv += `# Comparing WT (0 mutations) vs Mutant (2 mutations) by cancer type\n`;
            csv += 'Cancer Type,N (WT),r (WT),slope (WT),N (Mut),r (Mut),slope (Mut),Δr,p(Δr),Δslope,p(Δslope)\n';
            tableData.forEach(row => {
                csv += `"${row.lineage}",${row.nWT},${row.rWT.toFixed(4)},${row.slopeWT.toFixed(4)},${row.nMut},${row.rMut.toFixed(4)},${row.slopeMut.toFixed(4)},${row.deltaR.toFixed(4)},${row.pR.toExponential(2)},${row.deltaSlope.toFixed(4)},${row.pSlope.toExponential(2)}\n`;
            });
            this.downloadFile(csv, `correlation_${gene1}_vs_${gene2}_by_${hotspotGene}_mutation.csv`, 'text/csv');
        });

        // Make table sortable
        this.setupSortableTable('compareByCancerTable');

        // Make rows clickable - clicking a cancer type shows scatter with that cancer + hotspot filter
        document.querySelectorAll('#compareByCancerTable .clickable-row').forEach(row => {
            row.addEventListener('click', () => {
                const lineage = row.dataset.lineage;
                // Set cancer type filter
                document.getElementById('scatterCancerFilter').value = lineage;
                this.updateScatterSubtypeFilter();
                // Keep the hotspot gene as color overlay
                document.getElementById('hotspotGene').value = hotspotGene;
                document.getElementById('hotspotMode').value = 'color';
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
                    <td style="text-align: center; color: ${xColor};">${row.meanX.toFixed(3)}</td>
                    <td style="text-align: center; color: ${yColor};">${row.meanY.toFixed(3)}</td>
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
                csv += `"${row.cancerType}",${row.n},${row.correlation.toFixed(4)},${row.slope.toFixed(4)},${row.meanX.toFixed(4)},${row.meanY.toFixed(4)}\n`;
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

    switchToInspectWithTissue(tissue) {
        // Switch from By Tissue view back to Inspect scatter plot with the selected tissue preset
        if (!this.currentInspect || !this.currentInspect.data) return;

        const { gene1, gene2, data } = this.currentInspect;

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
            cancerFilter.innerHTML = `<option value="">All cancer types (n=${data.length})</option>`;
            lineages.forEach(l => {
                cancerFilter.innerHTML += `<option value="${l}">${l} (n=${lineageCounts[l]})</option>`;
            });
            // Set the filter to the selected tissue
            cancerFilter.value = tissue;
            cancerBox.style.display = 'block';
            // Update subtype filter for the selected tissue
            this.updateScatterSubtypeFilter();
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
        document.getElementById('geSummaryGene').textContent = geneUpper;
        document.getElementById('geSummaryMean').textContent = mean.toFixed(3);
        document.getElementById('geSummarySD').textContent = sd.toFixed(3);
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

        // Show modal
        document.getElementById('geneEffectModal').style.display = 'flex';

        // Render the selected view
        this.switchGeneEffectView(view);
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

    getGETissueFilteredData() {
        if (!this.currentGeneEffect) return [];
        const tissueFilter = document.getElementById('geTissueFilter')?.value;
        if (tissueFilter) {
            return this.currentGeneEffect.data.filter(d => d.lineage === tissueFilter);
        }
        return this.currentGeneEffect.data;
    }

    renderGeneEffectByTissue() {
        if (!this.currentGeneEffect) return;

        const data = this.getGETissueFilteredData();
        const gene = this.currentGeneEffect.gene;

        // Get all gene effects for comparison
        const allEffects = data.map(d => d.geneEffect);

        // Group data by cancer type (keep full data for box plots)
        const groupedData = {};
        data.forEach(d => {
            const lineage = d.lineage || 'Unknown';
            if (!groupedData[lineage]) groupedData[lineage] = [];
            groupedData[lineage].push({
                geneEffect: d.geneEffect,
                cellLineName: d.cellLineName || d.cellLineId,
                cellLineId: d.cellLineId
            });
        });

        // Calculate stats for each group including p-value vs all cells
        const stats = [];
        Object.entries(groupedData).forEach(([lineage, cellData]) => {
            if (cellData.length >= 3) {
                const effects = cellData.map(c => c.geneEffect);
                const mean = effects.reduce((a, b) => a + b, 0) / effects.length;
                const sd = Math.sqrt(effects.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / effects.length);

                // Calculate p-value comparing this group to all other cells
                const otherEffects = data.filter(d => (d.lineage || 'Unknown') !== lineage).map(d => d.geneEffect);
                const tTest = this.welchTTest(effects, otherEffects);

                stats.push({
                    group: lineage,
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

        // Create box plot traces for each cancer type
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
            title: { text: `${gene} by Cancer Type`, font: { size: 13 } },
            xaxis: {
                title: 'Gene Effect',
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

        Plotly.newPlot('geneEffectPlot', traces, layout, { responsive: true });

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
                title: 'Gene Effect',
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

        Plotly.newPlot('geneEffectHotspotPlot', traces, layout, { responsive: true });

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
    }

    renderGETable(stats, mode) {
        const tbody = document.getElementById('geneEffectTableBody');
        const thead = document.getElementById('geTableHead');
        this.currentGETableMode = mode;

        const headerStyle = 'cursor: pointer; user-select: none;';
        const sortIcon = ' ↕';

        if (mode === 'tissue') {
            thead.innerHTML = `<tr>
                <th style="${headerStyle}" data-sort="group" data-type="string">Cancer Type${sortIcon}</th>
                <th style="${headerStyle}" data-sort="n" data-type="number">N${sortIcon}</th>
                <th style="${headerStyle}" data-sort="mean" data-type="number">Mean GE${sortIcon}</th>
                <th style="${headerStyle}" data-sort="sd" data-type="number">SD${sortIcon}</th>
                <th style="${headerStyle}" data-sort="pValue" data-type="number">p-value${sortIcon}</th>
            </tr>`;
            this.renderGETableBody(stats, mode);
        } else {
            thead.innerHTML = `<tr>
                <th style="${headerStyle}" data-sort="group" data-type="string">Hotspot${sortIcon}</th>
                <th style="${headerStyle}; border-left: 2px solid #2563eb;" data-sort="n0" data-type="number">n(0)${sortIcon}</th>
                <th style="${headerStyle}" data-sort="mean0" data-type="number">GE(0)${sortIcon}</th>
                <th style="${headerStyle}; border-left: 2px solid #f97316;" data-sort="n1" data-type="number">n(1)${sortIcon}</th>
                <th style="${headerStyle}" data-sort="mean1" data-type="number">GE(1)${sortIcon}</th>
                <th style="${headerStyle}; border-left: 2px solid #dc2626;" data-sort="n2" data-type="number">n(2)${sortIcon}</th>
                <th style="${headerStyle}" data-sort="mean2" data-type="number">GE(2)${sortIcon}</th>
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
                    <td style="text-align: center; font-weight: 500; color: ${color};">${s.mean.toFixed(3)}</td>
                    <td style="text-align: center;">${s.sd.toFixed(3)}</td>
                    <td style="text-align: center; ${s.pValue < 0.05 ? 'font-weight: 600;' : ''}">${pStr}</td>
                </tr>`;
            });
        } else {
            stats.forEach(s => {
                const pStr = s.pValue < 0.001 ? '<0.001' : s.pValue.toFixed(3);
                const mean1Str = isNaN(s.mean1) ? '-' : s.mean1.toFixed(3);
                const mean2Str = isNaN(s.mean2) ? '-' : s.mean2.toFixed(3);
                const diffStr = s.diff !== undefined ? s.diff.toFixed(3) : '-';
                const diffColor = s.diff > 0 ? '#16a34a' : s.diff < 0 ? '#dc2626' : '#374151';
                tbody.innerHTML += `<tr class="clickable-row" data-group="${s.group}" style="cursor: pointer;">
                    <td>${s.group}</td>
                    <td style="text-align: center; color: #2563eb; border-left: 2px solid #2563eb;">${s.n0}</td>
                    <td style="text-align: center; color: #2563eb;">${s.mean0.toFixed(3)}</td>
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
            return `${label}: n=${stats.n}, GE=${stats.mean.toFixed(3)}, SD=${stats.sd.toFixed(3)}`;
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
            statsAnnotations.push(`0 (WT): n=${stats0.n}, GE=${stats0.mean.toFixed(3)}, SD=${stats0.sd.toFixed(3)}`);
            // Row 2: 1 mut stats
            if (stats1.n > 0) {
                statsAnnotations.push(`1 mut: n=${stats1.n}, GE=${stats1.mean.toFixed(3)}${stats1.n > 1 ? `, SD=${stats1.sd.toFixed(3)}` : ''}`);
            }
            // Row 3: 2 mut stats
            if (stats2.n > 0) {
                statsAnnotations.push(`2 mut: n=${stats2.n}, GE=${stats2.mean.toFixed(3)}${stats2.n > 1 ? `, SD=${stats2.sd.toFixed(3)}` : ''}`);
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
        statsAnnotations.push(`All cells: n=${allStats.n}, GE=${allStats.mean.toFixed(3)}, SD=${allStats.sd.toFixed(3)}`);
        statsAnnotations.push(`${group}: n=${groupStats.n}, GE=${groupStats.mean.toFixed(3)}, SD=${groupStats.sd.toFixed(3)}`);
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
        if (btn) {
            btn.style.display = this.geDetailedView ? 'inline-block' : 'none';
        }
        if (ratioControl) {
            ratioControl.style.display = this.geDetailedView ? 'flex' : 'none';
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

    downloadGeneEffectChartPNG() {
        // Mutation inspect mode handled by downloadGeneEffectPNG
        if (this.geneEffectViewMode === 'mutation') return;
        if (!this.currentGeneEffect) return;
        const plotId = this.currentGEView === 'tissue' ? 'geneEffectPlot' : 'geneEffectHotspotPlot';
        const plotEl = document.getElementById(plotId);
        const chartWidth = plotEl?._fullLayout?.width || 800;
        const chartHeight = plotEl?._fullLayout?.height || (this.geDetailedView ? 650 : 550);
        Plotly.downloadImage(plotId, {
            format: 'png',
            width: chartWidth,
            height: chartHeight,
            filename: `gene_effect_${this.currentGeneEffect.gene}_by_${this.currentGEView}`
        });
    }

    downloadGeneEffectChartSVG() {
        // Mutation inspect mode handled by downloadGeneEffectSVG
        if (this.geneEffectViewMode === 'mutation') return;
        if (!this.currentGeneEffect) return;
        const plotId = this.currentGEView === 'tissue' ? 'geneEffectPlot' : 'geneEffectHotspotPlot';
        const plotEl = document.getElementById(plotId);
        const chartWidth = plotEl?._fullLayout?.width || 800;
        const chartHeight = plotEl?._fullLayout?.height || (this.geDetailedView ? 650 : 550);
        Plotly.downloadImage(plotId, {
            format: 'svg',
            width: chartWidth,
            height: chartHeight,
            filename: `gene_effect_${this.currentGeneEffect.gene}_by_${this.currentGEView}`
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
                csv += `"${s.group}",${s.n},${s.mean.toFixed(4)},${s.sd.toFixed(4)},${s.pValue.toFixed(6)}\n`;
            });
        } else {
            csv = 'Hotspot,N_Mutant,N_WT,Mean_WT,Mean_Mutant,Delta_GE,SD_Mutant,p_value\n';
            this.currentGEStats.forEach(s => {
                csv += `"${s.group}",${s.nMut},${s.nWT},${s.wtMean.toFixed(4)},${s.mutMean.toFixed(4)},${s.diff.toFixed(4)},${s.mutSD.toFixed(4)},${s.pValue.toFixed(6)}\n`;
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
                csv += `"${d.cellLine}","${d.cellName}","${d.lineage}","${subtype}",${d.ge.toFixed(4)},${d.mutLevel}\n`;
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
            csv += `"${d.cellLineId}","${d.cellLineName}","${d.lineage}","${subtype}",${d.geneEffect.toFixed(4)}\n`;
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
            'lfc', 'fdr', 'hasCorrelation'
        ];
        const numericColKeys = [
            'n_wt', 'mean_wt', 'n_mut', 'mean_mut', 'diff_mut', 'p_mut',
            'n_2', 'mean_2', 'diff_2', 'p_2'
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
                input.placeholder = isNumeric ? '>0.5' : 'filter';
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

    // ===== Mutation Analysis Compare by Tissue/Hotspot (#16) =====

    showMutationCompareByTissue() {
        if (!this.mutationResults) return;
        const mr = this.mutationResults;
        const hotspotGene = mr.hotspotGene;
        const mutationData = this.mutations.geneData[hotspotGene];
        if (!mutationData) return;

        // Get all unique lineages and track all WT/mut indices
        const cellLines = this.metadata.cellLines;
        const lineageMap = {};
        const allWT = [], allMut = [];
        cellLines.forEach((cellLine, idx) => {
            if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) return;
            if (mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) return;
            if (this.excludedTissues && this.excludedTissues.size > 0) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                if (lineage && this.excludedTissues.has(lineage)) return;
            }

            const lineage = this.cellLineMetadata?.lineage?.[cellLine] || 'Unknown';
            if (!lineageMap[lineage]) lineageMap[lineage] = { wt: [], mut: [] };
            const mutLevel = mutationData.mutations[cellLine] || 0;
            if (mutLevel === 0) { lineageMap[lineage].wt.push(idx); allWT.push(idx); }
            else { lineageMap[lineage].mut.push(idx); allMut.push(idx); }
        });

        const topGenes = mr.significantResults.slice(0, 10).map(r => r.gene);
        if (topGenes.length === 0) { alert('No significant genes to compare.'); return; }

        // Helper to calculate gene diffs for a set of WT/mut indices
        const calcGeneDiffs = (wtIndices, mutIndices) => {
            const geneDiffs = {};
            topGenes.forEach(gene => {
                const geneIdx = this.geneIndex.get(gene);
                if (geneIdx === undefined) return;
                const geneData = this.getGeneData(geneIdx);
                const wtEffects = wtIndices.map(i => geneData[i]).filter(v => !isNaN(v));
                const mutEffects = mutIndices.map(i => geneData[i]).filter(v => !isNaN(v));
                if (wtEffects.length >= 3 && mutEffects.length >= 1) {
                    const meanWT = wtEffects.reduce((a, b) => a + b, 0) / wtEffects.length;
                    const meanMut = mutEffects.reduce((a, b) => a + b, 0) / mutEffects.length;
                    geneDiffs[gene] = meanMut - meanWT;
                }
            });
            return geneDiffs;
        };

        // Build "All" row from pooled data
        const allGeneDiffs = calcGeneDiffs(allWT, allMut);
        const allRow = { lineage: 'All', nWT: allWT.length, nMut: allMut.length, geneDiffs: allGeneDiffs, isAll: true };

        // Build per-tissue rows
        const tableData = [];
        Object.entries(lineageMap).forEach(([lineage, groups]) => {
            if (groups.wt.length < 3 || groups.mut.length < 1) return;
            const geneDiffs = calcGeneDiffs(groups.wt, groups.mut);
            if (Object.keys(geneDiffs).length === 0) return;
            tableData.push({ lineage, nWT: groups.wt.length, nMut: groups.mut.length, geneDiffs });
        });

        // Sort by absolute average delta across genes
        const avgAbsDelta = (row) => { const vals = Object.values(row.geneDiffs); return vals.length > 0 ? vals.reduce((a, b) => a + Math.abs(b), 0) / vals.length : 0; };
        tableData.sort((a, b) => avgAbsDelta(b) - avgAbsDelta(a));

        const allRows = [allRow, ...tableData];

        // Build HTML
        const buildRow = (row) => {
            const rowStyle = row.isAll ? 'background: #f9fafb; border-bottom: 2px solid #d1d5db;' : '';
            const pinned = row.isAll ? ' data-pinned="1"' : '';
            const tissue = row.isAll ? '' : row.lineage;
            return `<tr style="${rowStyle}"${pinned} data-nwt="${row.nWT}" data-nmut="${row.nMut}">
                <td><b>${row.lineage}</b></td>
                <td style="text-align: center;">${row.nWT}</td>
                <td style="text-align: center;">${row.nMut}</td>
                ${topGenes.map(g => {
                    const d = row.geneDiffs[g];
                    if (d === undefined) return '<td style="text-align: center; color: #999;">-</td>';
                    const c = d < 0 ? '#dc2626' : '#16a34a';
                    return `<td style="text-align: center; color: ${c}; cursor: pointer; text-decoration: underline;" data-val="${d}" onclick="app.openCompareInspect('${g}', '${tissue}')">${d.toFixed(3)}</td>`;
                }).join('')}
            </tr>`;
        };

        let html = `
            <div style="margin-bottom: 10px;">
                <h4 style="margin: 0;">${hotspotGene} mutation effect by Cancer Type</h4>
            </div>
            <div style="background: #f0f9ff; border: 1px solid #0ea5e9; border-radius: 4px; padding: 8px 12px; margin-bottom: 10px; font-size: 11px; color: #0c4a6e;">
                <b>What this shows:</b> For each cancer type, we split cell lines into WT (no ${hotspotGene} mutation) and mutated groups, then calculate the difference in gene effect (Δ GE = mean mutated − mean WT) for the top ${topGenes.length} most significant genes. <b>Red</b> = lower gene effect in mutated cells (stronger dependency), <b>green</b> = higher. Click a value to inspect that gene filtered by tissue.
            </div>
            <table class="data-table" id="mutCompareTable" style="font-size: 11px; width: 100%;">
                <thead><tr>
                    <th style="cursor: pointer;" data-col="0">Cancer Type ▼</th>
                    <th style="cursor: pointer;" data-col="1">N (WT)</th>
                    <th style="cursor: pointer;" data-col="2">N (Mut)</th>
                    ${topGenes.map((g, i) => `<th style="font-size: 10px; cursor: pointer;" data-col="${3 + i}">${g}</th>`).join('')}
                </tr></thead>
                <tbody>
                    ${allRows.map(row => buildRow(row)).join('')}
                </tbody>
            </table>
        `;

        const container = document.getElementById('mutationCompareTable');
        container.innerHTML = html;
        container.style.display = 'block';
        document.getElementById('mutCompareCloseBtn').style.display = 'inline-block';
        this.attachCompareTableSort('mutCompareTable');
    }

    showMutationCompareByHotspot() {
        if (!this.mutationResults) return;
        const mr = this.mutationResults;
        const mainHotspot = mr.hotspotGene;

        // Get top significant genes
        const topGenes = mr.significantResults.slice(0, 10).map(r => r.gene);
        if (topGenes.length === 0) { alert('No significant genes to compare.'); return; }

        // Get cell line indices that match the mutation analysis filters
        const cellLines = this.metadata.cellLines;
        const filteredIndices = [];
        cellLines.forEach((cellLine, idx) => {
            if (mr.lineageFilter && this.cellLineMetadata?.lineage?.[cellLine] !== mr.lineageFilter) return;
            if (mr.subLineageFilter && this.cellLineMetadata?.primaryDisease?.[cellLine] !== mr.subLineageFilter) return;
            if (this.excludedTissues && this.excludedTissues.size > 0) {
                const lineage = this.cellLineMetadata?.lineage?.[cellLine];
                if (lineage && this.excludedTissues.has(lineage)) return;
            }
            filteredIndices.push(idx);
        });

        // Helper to calculate gene diffs for given WT/mut indices
        const calcGeneDiffs = (wtIdx, mutIdx) => {
            const geneDiffs = {};
            topGenes.forEach(gene => {
                const geneIdx = this.geneIndex.get(gene);
                if (geneIdx === undefined) return;
                const geneData = this.getGeneData(geneIdx);
                const wtEffects = wtIdx.map(i => geneData[i]).filter(v => !isNaN(v));
                const mutEffects = mutIdx.map(i => geneData[i]).filter(v => !isNaN(v));
                if (wtEffects.length >= 3 && mutEffects.length >= 3) {
                    const meanWT = wtEffects.reduce((a, b) => a + b, 0) / wtEffects.length;
                    const meanMut = mutEffects.reduce((a, b) => a + b, 0) / mutEffects.length;
                    geneDiffs[gene] = meanMut - meanWT;
                }
            });
            return geneDiffs;
        };

        // Build reference row for the main hotspot
        const mainMutData = this.mutations.geneData[mainHotspot]?.mutations || {};
        const mainWT = filteredIndices.filter(i => (mainMutData[cellLines[i]] || 0) === 0);
        const mainMut = filteredIndices.filter(i => (mainMutData[cellLines[i]] || 0) >= 1);
        const mainRow = {
            hotspotGene: `${mainHotspot} (reference)`, nWT: mainWT.length, nMut: mainMut.length,
            geneDiffs: calcGeneDiffs(mainWT, mainMut), isRef: true
        };

        // For each other hotspot mutation gene, calculate diffs for top genes
        const tableData = [];
        this.mutations.genes.forEach(hotspotGene => {
            if (hotspotGene === mainHotspot) return;
            const mutData = this.mutations.geneData[hotspotGene]?.mutations || {};
            const wtIndices = filteredIndices.filter(i => (mutData[cellLines[i]] || 0) === 0);
            const mutIndices = filteredIndices.filter(i => (mutData[cellLines[i]] || 0) >= 1);
            if (wtIndices.length < 3 || mutIndices.length < 3) return;
            const geneDiffs = calcGeneDiffs(wtIndices, mutIndices);
            if (Object.keys(geneDiffs).length === 0) return;
            tableData.push({ hotspotGene, nWT: wtIndices.length, nMut: mutIndices.length, geneDiffs });
        });

        const avgAbsDelta = (row) => { const vals = Object.values(row.geneDiffs); return vals.length > 0 ? vals.reduce((a, b) => a + Math.abs(b), 0) / vals.length : 0; };
        tableData.sort((a, b) => avgAbsDelta(b) - avgAbsDelta(a));

        const allRows = [mainRow, ...tableData.slice(0, 30)];

        // Build HTML
        const buildRow = (row) => {
            const rowStyle = row.isRef ? 'background: #f9fafb; border-bottom: 2px solid #d1d5db;' : '';
            const pinned = row.isRef ? ' data-pinned="1"' : '';
            const label = row.isRef ? `${mainHotspot} (reference)` : row.hotspotGene;
            return `<tr style="${rowStyle}"${pinned} data-nwt="${row.nWT}" data-nmut="${row.nMut}">
                <td><b>${label}</b></td>
                <td style="text-align: center;">${row.nWT}</td>
                <td style="text-align: center;">${row.nMut}</td>
                ${topGenes.map(g => {
                    const d = row.geneDiffs[g];
                    if (d === undefined) return '<td style="text-align: center; color: #999;">-</td>';
                    const c = d < 0 ? '#dc2626' : '#16a34a';
                    return `<td style="text-align: center; color: ${c}; cursor: pointer; text-decoration: underline;" data-val="${d}" onclick="app.openCompareInspect('${g}', '')">${d.toFixed(3)}</td>`;
                }).join('')}
            </tr>`;
        };

        let html = `
            <div style="margin-bottom: 10px;">
                <h4 style="margin: 0;">Co-occurring mutations affecting ${mainHotspot}-sensitive genes</h4>
            </div>
            <div style="background: #f0f9ff; border: 1px solid #0ea5e9; border-radius: 4px; padding: 8px 12px; margin-bottom: 10px; font-size: 11px; color: #0c4a6e;">
                <b>What this shows:</b> For each other hotspot mutation, we check whether cells carrying that mutation also show altered gene effects for the top ${topGenes.length} genes identified in the ${mainHotspot} analysis. This helps identify co-occurring mutations that may amplify or counteract the effect. <b>Red</b> = lower gene effect in mutated cells, <b>green</b> = higher. Click a value to inspect that gene.
            </div>
            <table class="data-table" id="mutCompareTable" style="font-size: 11px; width: 100%;">
                <thead><tr>
                    <th style="cursor: pointer;" data-col="0">Hotspot ▼</th>
                    <th style="cursor: pointer;" data-col="1">N (WT)</th>
                    <th style="cursor: pointer;" data-col="2">N (Mut)</th>
                    ${topGenes.map((g, i) => `<th style="font-size: 10px; cursor: pointer;" data-col="${3 + i}">${g}</th>`).join('')}
                </tr></thead>
                <tbody>
                    ${allRows.map(row => buildRow(row)).join('')}
                </tbody>
            </table>
        `;

        const container = document.getElementById('mutationCompareTable');
        container.innerHTML = html;
        container.style.display = 'block';
        document.getElementById('mutCompareCloseBtn').style.display = 'inline-block';
        this.attachCompareTableSort('mutCompareTable');
    }

    // ===== Tissue exclusion from analysis (#17) =====

    attachCompareTableSort(tableId) {
        const table = document.getElementById(tableId);
        if (!table) return;
        const headers = table.querySelectorAll('thead th');
        headers.forEach(th => {
            th.addEventListener('click', () => {
                const colIdx = parseInt(th.dataset.col);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                // Separate pinned rows (All / reference) from the rest
                const pinnedRows = rows.filter(r => r.dataset.pinned === '1');
                const sortableRows = rows.filter(r => r.dataset.pinned !== '1');

                const dir = th._sortDir === 'asc' ? 'desc' : 'asc';
                th._sortDir = dir;
                // Clear sort indicators from other headers
                headers.forEach(h => { h.textContent = h.textContent.replace(/ [▲▼]$/, ''); });
                th.textContent += dir === 'asc' ? ' ▲' : ' ▼';

                sortableRows.sort((a, b) => {
                    const cellA = a.children[colIdx];
                    const cellB = b.children[colIdx];
                    let valA, valB;
                    if (colIdx === 0) {
                        valA = cellA.textContent.trim().toLowerCase();
                        valB = cellB.textContent.trim().toLowerCase();
                        return dir === 'asc' ? valA.localeCompare(valB) : valB.localeCompare(valA);
                    }
                    valA = parseFloat(cellA.dataset.val ?? cellA.textContent.replace(/[^0-9.\-eE+]/g, '')) || 0;
                    valB = parseFloat(cellB.dataset.val ?? cellB.textContent.replace(/[^0-9.\-eE+]/g, '')) || 0;
                    return dir === 'asc' ? valA - valB : valB - valA;
                });

                // Re-append: pinned first, then sorted
                pinnedRows.forEach(r => tbody.appendChild(r));
                sortableRows.forEach(r => tbody.appendChild(r));
            });
        });
    }

    openCompareInspect(gene, tissue) {
        this.showGeneEffectDistribution(gene, tissue || '');
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
}

// Initialize app
const app = new CorrelationExplorer();
