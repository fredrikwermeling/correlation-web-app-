// Tour of the app: one scrolling page that presents Correlate the way a
// product page would, with real screenshots of the app (tour_img/*.webp,
// captured from the live views) beside short feature text, and a button
// under each that opens the same view with the same example loaded.
// Reads window.app only. Counts come from the loaded data where available.
(function () {
    'use strict';

    const A = () => window.app || null;
    const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, c =>
        ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
    const num = (n) => Number(n).toLocaleString('en-US');
    const gi = (g) => `<i>${esc(g)}</i>`;
    const IMG = 'tour_img/';
    // The pictures are cached by URL; the script's own ?v= tag renews them.
    const IMG_V = (() => { const m = (document.currentScript?.src || [...document.scripts].map(x => x.src).find(x => /tour\.js/.test(x)) || '').match(/[?&]v=([^&]+)/); return m ? '?v=' + m[1] : ''; })();

    // ------------------------------------------------------------ data
    const D = {
        nLines() { try { return A().metadata.cellLines.length; } catch (e) { return 1208; } },
        nGenes() { try { return A().nGenes || A().metadata.genes.length; } catch (e) { return 18531; } },
        nExprLines() { try { return A().expressionMetadata?.nCellLines || 1718; } catch (e) { return 1718; } },
        lineId(name) {
            try {
                const map = A()._buildCellLineNameToIdMap();
                return map.get(String(name).toUpperCase()) || null;
            } catch (e) { return null; }
        },
    };

    // ------------------------------------------------------------ actions
    // Every action reproduces the state shown in the screenshot beside it.
    const HERO_SET = ['TP53', 'MDM2', 'MDM4', 'PPM1D', 'USP7', 'CDKN1A', 'TP53BP1', 'CHEK2', 'ATM', 'TSC1', 'TSC2'];
    const setVal = (id, v) => {
        const el = document.getElementById(id); if (!el) return;
        el.value = v;
        el.dispatchEvent(new Event('input', { bubbles: true }));
        el.dispatchEvent(new Event('change', { bubbles: true }));
    };
    function runExampleSet(o) {
        o = o || {};
        const a = A();
        document.getElementById('modeGeneSetBtn')?.click();
        // The plain list tab, and none of the screen-results options a
        // previous example may have left on.
        [...document.querySelectorAll('.input-tab')].find(b => /paste/i.test(b.textContent))?.click();
        for (const id of ['showUncorrelatedGenes', 'colorByStats']) {
            const cb = document.getElementById(id);
            if (cb?.checked) cb.click();
        }
        const basis = document.querySelector(`input[name="analysisBasis"][value="${o.basis || 'ge'}"]`);
        if (basis && !basis.checked) { basis.checked = true; basis.dispatchEvent(new Event('change', { bubbles: true })); }
        if (o.cutoff) setVal('correlationCutoff', o.cutoff);
        // Written after the tab switch has settled: leaving the screen-results
        // tab refills the gene box from its table.
        const before = a.results;
        setTimeout(() => {
            const ta = document.getElementById('geneTextarea');
            if (ta) { ta.value = HERO_SET.filter(g => a.geneIndex.has(g)).join('\n'); ta.dispatchEvent(new Event('input', { bubbles: true })); }
            setTimeout(() => document.getElementById('runAnalysis')?.click(), 150);
        }, 300);
        if (!o.colorGE) return;
        let tries = 0;
        const arm = () => {
            if (a.results && a.results !== before && a.results.success) {
                const cb = document.getElementById('colorByGeneEffect');
                if (cb && !cb.checked) { cb.checked = true; cb.dispatchEvent(new Event('change', { bubbles: true })); }
                return;
            }
            if (tries++ < 80) setTimeout(arm, 250);
        };
        setTimeout(arm, 600);
    }
    const wait = (ms) => new Promise(r => setTimeout(r, ms));
    const tick = (id) => { const cb = document.getElementById(id); if (cb && !cb.checked) cb.click(); };
    const until = async (test, max) => { for (let i = 0; i < (max || 80); i++) { if (test()) return true; await wait(250); } return false; };
    // The "Test Genes" set with pasted screen statistics, nodes colored by
    // LFC, genes without a correlation shown too.
    async function runStatsExample(a, withStats) {
        document.getElementById('modeGeneSetBtn')?.click();
        [...document.querySelectorAll('.input-tab')].find(b => /with stats/i.test(b.textContent))?.click();
        await wait(300);
        a.loadTestGenesWithStats();
        setVal('correlationCutoff', 0.4);
        const before = a.results;
        await wait(200);
        document.getElementById('runAnalysis')?.click();
        await until(() => a.results && a.results !== before && a.results.success);
        if (withStats) { tick('showUncorrelatedGenes'); await wait(800); tick('colorByStats'); }
    }
    async function openMelanoma(a, hotspotGene) {
        a.openCellLineBrowser(); await wait(500);
        setVal('clbTissueFilter', 'Skin'); await wait(500);
        setVal('clbSubtypeFilter', 'Melanoma'); await wait(500);
        setVal('clbHotspotFilter', hotspotGene || ''); await wait(700);
    }
    async function compareMelanoma(a) {
        await openMelanoma(a, 'BRAF');
        a._resetGEInspectScope();
        a._geInspectScope = 'all';
        a._runBusy({ button: 'clbInspectGEBtn', label: 'Comparing' }, () => a.inspectSelectionGE());
    }
    async function pcaMelanoma(a) {
        await openMelanoma(a, 'BRAF');
        [...document.querySelectorAll('#clbModalCard button')].find(b => /PCA \/ UMAP/.test(b.textContent))?.click();
        await wait(800);
        setVal('clbUmapMethod', 'pca'); setVal('clbUmapDataType', 'expr');
        tick('clbUmapShowLoadings');
        const before = a._clbUmapData;
        document.getElementById('clbUmapRunBtn')?.click();
        await until(() => a._clbUmapData && a._clbUmapData !== before, 120);
        await wait(500);
        const cg = document.getElementById('clbUmapColorGene');
        if (cg) { cg.value = 'MITF'; cg.dispatchEvent(new Event('input', { bubbles: true })); }
        setVal('clbUmapColorGeneType', 'expr');
        document.getElementById('clbUmapColorGeneBtn')?.click();
    }
    async function heatmapMelanoma(a) {
        await openMelanoma(a, 'BRAF');
        document.getElementById('clbHeatmapBtn')?.click();
        await wait(800);
        setVal('hmLineage', 'Skin'); await wait(400);
        setVal('hmSubtype', 'Melanoma'); await wait(400);
        setVal('hmHotspotFilter', 'BRAF'); await wait(400);
        // Lineage as the one block, TP53 wild-type first inside it, then
        // clustered: the layout of the figure.
        a._hmAnnRows = [
            { mode: 'lineage', gene: null, sortDir: 'desc', sortKey: 'mag' },
            { mode: 'hotspot', gene: 'TP53', sortDir: 'asc', sortKey: 'cat-asc' },
            { mode: 'cluster', gene: null, sortDir: 'asc', sortKey: null },
        ];
        a._hmRenderAnnRowsBlock?.();
        setVal('hmPreset', 'tp53_members');
    }
    async function runMutationExample(a) {
        document.getElementById('modeMutationBtn')?.click();
        await wait(700);
        setVal('mutationHotspotSelect', 'BRAF');
        setVal('lineageFilter', 'Skin');
        await wait(500);
        const before = a.mutationResults;
        document.getElementById('runMutationAnalysisBtn')?.click();
        await until(() => a.mutationResults && a.mutationResults !== before, 120);
        await wait(400);
        a.showGeneEffectDistribution('BRAF');
    }
    // A one-time card beside the network listing what hovering and clicking
    // it does; kept off the drawing itself, whose layout varies per run.
    // Closed with its own button or by the next run.
    async function networkHint(a) {
        document.getElementById('tourNetHint')?.remove();
        const before = a.results;
        if (!await until(() => a.results && a.results !== before && a.results.success, 120)) return;
        await wait(600);
        const wrap = document.getElementById('networkWrap');
        if (!wrap || !wrap.offsetParent) return;
        const card = document.createElement('div');
        card.id = 'tourNetHint';
        card.innerHTML = `<button type="button" aria-label="Close">&times;</button>
            <b>Explore the network</b>
            <ul>
                <li><span>Hover</span> a gene or an edge to see what it is.</li>
                <li><span>Double-click a gene</span> for its gene effect across the cell lines.</li>
                <li><span>Double-click an edge</span> for the correlation plot behind it.</li>
                <li><span>Shift-click a gene</span> to open it in the Gene Browser.</li>
            </ul>`;
        card.querySelector('button').addEventListener('click', () => card.remove());
        document.body.appendChild(card);
        // Beside the network's left edge, pointing at it, in page coordinates
        // so it scrolls with the network. No room on the left (narrow window):
        // inside the network's top-left corner instead.
        const placeCard = () => {
            if (!card.isConnected) { removeEventListener('resize', placeCard); return; }
            const r = wrap.getBoundingClientRect();
            const w = card.offsetWidth, gap = 14;
            const beside = r.left - w - gap >= 8;
            card.classList.toggle('tour-hint-beside', beside);
            card.style.left = `${Math.round((beside ? r.left - w - gap : r.left + 10) + scrollX)}px`;
            card.style.top = `${Math.round(r.top + scrollY + (beside ? 40 : 10))}px`;
        };
        placeCard();
        addEventListener('resize', placeCard);
        document.getElementById('runAnalysis')?.addEventListener('click', () => card.remove(), { once: true });
        wrap.scrollIntoView({ behavior: 'smooth', block: 'center' });
    }
    // The same kind of card for the correlation plot: it sits beside the
    // popout's lower left corner and goes when the popout closes.
    async function plotHint() {
        document.getElementById('tourPlotHint')?.remove();
        const modal = document.getElementById('inspectModal');
        const vis = (m) => getComputedStyle(m).display !== 'none' && m.getClientRects().length > 0;
        const shown = () => modal && vis(modal);
        if (!await until(shown, 40)) return;
        await wait(900);
        const card = document.createElement('div');
        card.id = 'tourPlotHint';
        card.innerHTML = `<button type="button" aria-label="Close">&times;</button>
            <b>Explore the plot</b>
            <ul>
                <li><span>Hover a dot</span> to see which cell line it is.</li>
                <li><span>Click a dot</span> to label it, <span>Shift-click</span> to open the cell line's page.</li>
                <li><span>Drag</span> over an area to zoom in, <span>double-click</span> to zoom out.</li>
                <li>Use the panel on the left to filter, color the dots, highlight cell lines or compare tissues.</li>
                <li><span>Settings</span> under the plot changes sizes and fonts.</li>
            </ul>`;
        card.querySelector('button').addEventListener('click', () => card.remove());
        document.body.appendChild(card);
        // Follow the popout: gone when it closes, hidden under anything opened on top.
        const iv = setInterval(() => {
            if (!card.isConnected) { clearInterval(iv); return; }
            if (!shown()) { card.remove(); clearInterval(iv); return; }
            const top = [...document.querySelectorAll('.modal-overlay')].some(m => m !== modal && vis(m)
                && (+getComputedStyle(m).zIndex || 0) >= 1300);
            card.style.visibility = top ? 'hidden' : '';
        }, 400);
    }
    const ACTIONS = {
        hero: (a) => { runExampleSet({ cutoff: 0.5 }); return networkHint(a); },
        heroStats: (a) => { networkHint(a); return runStatsExample(a, true); },
        heroMrna: (a) => { runExampleSet({ cutoff: 0.3, basis: 'expr' }); return networkHint(a); },
        pair: (a) => { a.openInspectByGenes('TP53', 'MDM2'); return plotHint(); },
        eda2r: (a) => { a._inspectPreset = { colorByGene: 'EDA2R', colorByGeneType: 'expr' }; a.openInspectByGenes('TP53', 'MDM2'); return plotHint(); },
        mutation: (a) => runMutationExample(a),
        compare: (a) => compareMelanoma(a),
        pca: (a) => pcaMelanoma(a),
        heatmap: (a) => heatmapMelanoma(a),
        myc: (a) => a.openGeneBrowser('MYC'),
        tp53seq: async (a) => { a.openGeneBrowser('TP53'); await wait(2500); document.getElementById('gbStructCard')?.scrollIntoView({ behavior: 'smooth', block: 'start' }); },
        tsc1: (a) => a.openGeneBrowser('TSC1'),
        clb: (a) => openMelanoma(a, 'BRAF'),
        oncoprint: async (a) => {
            await openMelanoma(a, '');
            a.showOncoprint('clb', 'hotspot');
            await wait(700);
            // The picks shown in the figure, made through the grid's own click
            // handler: BRAF mutated kept, TP53 and CDKN2A mutated left out.
            const d = a._oncoprintData;
            const cv = [...document.querySelectorAll('#oncoprintPopup canvas')].find(c => c.getBoundingClientRect().width > 0 && c.getBoundingClientRect().width < 260);
            if (!d || !cv) return;
            const r = cv.getBoundingClientRect();
            for (const [gene, state] of [['BRAF', 'mut'], ['TP53', 'wt'], ['CDKN2A', 'wt']]) {
                const row = d.topGenes.findIndex(g => g.gene === gene);
                if (row < 0) continue;
                const x = 2 + (state === 'mut' ? d.boxW / 2 : d.boxW + d.boxGap + d.boxW / 2);
                cv.dispatchEvent(new MouseEvent('click', { bubbles: true, clientX: r.left + x, clientY: r.top + row * d.cellH + d.cellH / 2 }));
            }
        },
        wiki: (a) => a.openCellLineWiki(D.lineId('A375') || 'ACH-000219'),
    };

    // ------------------------------------------------------------ page
    // The V1 build (scripts/build_v1.py) leaves some tools out, and the tour
    // shows only what this build has; each test is that tool's markup.
    const HAS = {
        mrnaNet: () => !!document.querySelector('input[name="analysisBasis"]'),
        pca: () => !!document.getElementById('clbUmapSection'),
        heatmap: () => !!document.getElementById('heatmapModal'),
        cutoff: () => !!document.getElementById('cutoffModal'),
    };
    const shot = (file, alt, act, o) => {
        o = o || {};
        // A figure is a picture; the one way to open its view is the button
        // under it, so there is never a hidden click target and a button for
        // the same thing.
        return `<figure class="tour-shot${o.cls ? ' ' + o.cls : ''}">
            <div class="tour-shot-bar"><span></span><span></span><span></span></div>
            <img src="${IMG}${file}.webp${IMG_V}" alt="${esc(alt)}" loading="${o.eager ? 'eager' : 'lazy'}" decoding="async">
            ${o.cap ? `<figcaption>${o.cap}</figcaption>` : ''}
            ${act && o.btn ? `<div class="tour-shot-act">${tryBtn(act, o.btn)}</div>` : ''}
        </figure>`;
    };
    const tryBtn = (act, label) => `<button type="button" class="tour-try" data-act="${act}">${label}<span aria-hidden="true"> &rarr;</span></button>`;
    const checks = (items) => `<ul class="tour-checks">${items.map(t => `<li>${t}</li>`).join('')}</ul>`;

    // fx: an optional transition played while the page glides to the section
    // ('cube' or 'zoom'). Only two sections use one, so it stays a surprise.
    const nextBtn = (id, label, fx) => `<button type="button" class="tour-next" data-next="#${id}"${fx ? ` data-fx="${fx}"` : ''}>Next: ${label} <span aria-hidden="true">&darr;</span></button>`;

    function pageHTML() {
        const nL = num(D.nLines()), nE = num(D.nExprLines()), nG = num(D.nGenes()), nSet = HERO_SET.length;
        const mrnaNet = HAS.mrnaNet(), pca = HAS.pca(), heatmap = HAS.heatmap();
        const ways = ['against all other cell lines', pca && 'as a PCA', heatmap && 'as a heatmap of a gene set'].filter(Boolean);
        const otherTools = ['mutation analysis', 'gene effect or expression of one gene', 'the correlation plot',
            heatmap && 'the gene set heatmap', HAS.cutoff() && 'the threshold scan', 'Paste from AI', 'opening a saved figure'].filter(Boolean);
        return `
<header class="tour-top">
    <div class="tour-wrap tour-top-in">
        <div class="tour-brand"><b>Correlate</b><span>Tour of the app</span></div>
        <nav class="tour-nav" aria-label="Sections">
            <a href="#tour-network">Gene set analysis</a>
            <a href="#tour-corr">Correlation plot</a>
            <a href="#tour-cells">Cell Line Browser</a>
            <a href="#tour-gene">Gene Browser</a>
            <a href="#tour-more">Other features</a>
        </nav>
        <button type="button" class="tour-close" id="tourClose" title="Close the tour (Esc)">Close tour <span aria-hidden="true">&times;</span></button>
    </div>
</header>

<section class="tour-hero">
    <div class="tour-wrap tour-hero-grid">
        <div class="tour-hero-text">
        <h1>Correlate brings together CRISPR knockout screens, mRNA expression, mutations, copy number and drug response from DepMap, adds gene and cell line information from many other databases, and connects them.</h1>
        <p class="tour-sub">Use it to:</p>
        <ul class="tour-hero-list">
            <li class="tour-hero-fn"><button type="button" class="tour-modebtn" data-next="#tour-network">Gene set<br>analysis</button><span>Find connections within your gene list, for example hits from a CRISPR screen or a transcriptomics experiment, from unbiased correlations of CRISPR gene effect${mrnaNet ? ' or mRNA' : ''} across the cell lines.</span></li>
            <li class="tour-hero-fn"><button type="button" class="tour-modebtn" data-next="#tour-cells">Cell Line<br>Browser</button><span>Find the right cell lines for your project, by mutations, mutation patterns and cancer type, and compare omics data between groups of cell lines.</span></li>
            <li class="tour-hero-fn"><button type="button" class="tour-modebtn" data-next="#tour-gene">Gene<br>Browser</button><span>Get extensive information about any gene: dependency, expression, mutations, protein complexes, drugs and literature, gathered from DepMap, UniProt, Cellosaurus and other sources.</span></li>
            <li>Export any view for AI, to discuss the data with your favorite language model.</li>
        </ul>
        </div>
        <div class="tour-hero-shot">
            ${shot('net_p53', 'A network built from eleven p53 pathway and TSC genes', 'hero', { eager: true, cap: `A network from a list of ${nSet} genes: the p53 pathway and the TSC complex. Blue edges join genes that cells need together, red edges join opposites: cells that need ${gi('MDM2')} tolerate losing ${gi('TP53')}. In the app, hover over genes and edges for details and double-click them to open their plots.`, btn: 'Try it: build and explore this network' })}
        </div>
    </div>
    <div class="tour-wrap tour-next-row">${nextBtn('tour-network', 'Gene set analysis', 'cube')}</div>
</section>

<section class="tour-stats" id="tour-stats">
    <div class="tour-wrap tour-stats-in">
        <div><b>${nL}</b><span>cell lines with DepMap CRISPR gene effect</span></div>
        <div><b>${nE}</b><span>cell lines with DepMap mRNA expression</span></div>
        <div><b>${nG}</b><span>genes</span></div>
        <div><b>10+</b><span>other databases integrated, such as Cellosaurus, UniProt, Ensembl, Reactome, EBI Complex Portal and COSMIC</span></div>
    </div>
</section>

<section class="tour-sec" id="tour-network">
    <div class="tour-wrap">
        <div class="tour-text tour-text-wide">
            <p class="tour-eyebrow">1 &middot; Gene set analysis</p>
            <h2>Paste a gene set, get a network</h2>
            <p>Correlate calculates the correlation (Pearson r) between every pair of genes in your set across the cell lines, on CRISPR gene effect (${nL} cell lines)${mrnaNet ? ' or on mRNA expression' : ''}. Pairs whose correlation passes the cutoff you set are joined by a line in the network: blue for a positive correlation, red for a negative one, thicker for a stronger one. Genes that end up joined often work together, in the same complex or pathway.</p>
            ${checks([
                'Hover over a gene or an edge for details. Double-click a gene for its gene effect across the cell lines, or an edge for the correlation plot behind it.',
                'Paste your own screen results (LFC and FDR) next to the gene list, and Correlate shows the set as a network colored by your results, or as a correlation matrix.',
                'Expand the set to find new partners that are not in your list, and send genes to Enrichr to see which known gene sets (Reactome, GO and others) they overlap with.',
            ])}
        </div>
        <div class="tour-two tour-two-net${mrnaNet ? '' : ' tour-two-one'}">
            ${shot('net_stats', 'A list of 22 cancer genes with screen statistics pasted alongside, nodes colored by LFC', 'heroStats', { cap: 'A network of 22 genes with statistics from your own experiment added: each gene is colored by its value, here the log fold change (LFC) from a screen.', btn: 'Build it with the screen results' })}
            ${mrnaNet ? shot('net_mrna', 'The same eleven genes built on mRNA expression', 'heroMrna', { cap: `The ${nSet} genes from the top of the page, as a network built on mRNA expression data.`, btn: 'Build this set on mRNA' }) : ''}
        </div>
    </div>
    <div class="tour-wrap tour-next-row">${nextBtn('tour-corr', 'Correlation plot')}</div>
</section>

<section class="tour-sec tour-alt" id="tour-corr">
    <div class="tour-wrap">
        <div class="tour-text tour-text-wide">
            <p class="tour-eyebrow">2 &middot; Correlation plot</p>
            <h2>Double-click an edge to see what the connection is based on</h2>
            <p>An edge in the network connects two genes. Double-click it and the correlation plot opens: each dot is a cell line, placed by its value for one gene against the other. The correlation (Pearson r) says how closely the two genes rise and fall together across the cell lines: 1 means they always move together, 0 means no relation, and −1 means one goes up when the other goes down. The p-value says how likely an r this strong would be by chance. The plot can also be opened on its own, from <b>Other</b> and then <b>Correlation</b>, to explore the correlation between any two genes without building a network.</p>
        </div>
        <div class="tour-wide-shot tour-edge-shot">
            ${shot('edge_to_plot', 'Double-clicking the edge between TP53 and MDM2 opens their correlation plot', 'pair', { cap: `The edge between ${gi('TP53')} and ${gi('MDM2')} is red: a negative correlation, r −0.72 across ${nL} cell lines. MDM2 holds p53 in check, so cell lines with working p53 depend on ${gi('MDM2')} (negative gene effect), while knocking out ${gi('TP53')} lets them grow faster (positive gene effect).`, btn: 'Open TP53 vs MDM2, then hover, click and change the settings' })}
        </div>
        <div class="tour-feature tour-flip">
            <div class="tour-text">
                <h3>Then ask why</h3>
                ${checks([
                    'Filter to a cancer type or a mutation. r and p follow the filter.',
                    `Color the dots by tissue, by a mutation, or by a third gene. <b>Suggest</b> finds the gene that explains the pattern best.`,
                    'Highlight cell lines by typing their names or clicking their dots, and compare tissues or mutation groups side by side.',
                    '<b>Export for AI</b> sends the plot with the data behind it, so you can ask an AI what could explain the pattern.',
                ])}
            </div>
            ${shot('scatter_eda2r', 'TP53 against MDM2 gene effect, colored by EDA2R mRNA', 'eda2r', { cap: `The same plot, each dot a cell line, colored by its ${gi('EDA2R')} mRNA (a p53 target gene): the cell lines that depend on ${gi('MDM2')} are the ones with active p53.`, btn: 'Open it colored by EDA2R' })}
        </div>
    </div>
    <div class="tour-wrap tour-next-row">${nextBtn('tour-cells', 'Cell Line Browser', 'zoom')}</div>
</section>

<section class="tour-sec" id="tour-cells">
    <div class="tour-wrap">
        <div class="tour-feature">
            <div class="tour-text">
                <p class="tour-eyebrow">3 &middot; Cell Line Browser</p>
                <h2>Find the right cell lines for your project</h2>
                <p>Narrow the ${nL} cell lines down by cancer type, mutation, fusion, copy number, sex, drug response, or how much they express or need any gene. Every cell line has its own Wiki, a comprehensive summary.</p>
                ${checks([
                    `Example: melanoma cell lines with mutated ${gi('BRAF')}. The count updates as you go.`,
                    `From the list: compare the selection with the rest${pca ? ', a PCA or UMAP' : ''}${heatmap ? ', a heatmap of any gene set' : ''}.`,
                    'Download a cell line’s full tables, or the list on screen, as CSV.',
                ])}
            </div>
            ${shot('clb_mel', 'The Cell Line Browser filtered to melanoma cell lines with mutated BRAF, A375 selected', 'clb', { cap: `Melanoma cell lines with mutated ${gi('BRAF')}, A375 opened in the side panel.`, btn: 'Open the browser with these filters' })}
        </div>
        <div class="tour-feature tour-flip tour-flip-wide">
            <div class="tour-text">
                <h3>Select cell lines by their mutations</h3>
                ${checks([
                    'The mutation grid (oncoprint) shows the most often mutated genes across the cell lines on screen, one column per cell line.',
                    'Tick green to keep the cell lines that carry a mutation and red to leave them out, combine several genes, then apply. Copy number and fusions have their own grids.',
                ])}
            </div>
            ${shot('oncoprint_mel', 'The mutation grid for the 67 melanoma cell lines, BRAF mutated kept, TP53 and CDKN2A mutated left out', 'oncoprint', { cap: `The 67 melanoma cell lines: ${gi('BRAF')}-mutated kept (green), ${gi('TP53')}- and ${gi('CDKN2A')}-mutated left out (red), which leaves 38.`, btn: 'Open the mutation grid with these picks' })}
        </div>
        <div class="tour-text tour-text-wide" style="margin-top:56px;">
            <h3>${ways.length > 1 ? 'Compare the cell lines in different ways' : 'Compare them with the other cell lines'}</h3>
            <p>The same 50 <i>BRAF</i>-mutated melanoma cell lines${ways.length > 1 ? `, ${ways.length === 3 ? 'three' : 'two'} ways: ${ways.slice(0, -1).join(', ')}${ways.length === 3 ? ',' : ''} and ${ways[ways.length - 1]}` : ` ${ways[0]}`}.</p>
        </div>
        <div class="tour-two tour-two-tight ${pca ? 'tour-two-pca' : 'tour-two-one'}">
            ${shot('volcano_mel', 'Volcano plot comparing the 50 BRAF-mutant melanoma cell lines with every other cell line', 'compare', { cap: `Compare: the 50 selected cell lines against the other 1,158. ${gi('BRAF')}, ${gi('SOX10')} and ${gi('MAPK1')} are their dependencies.`, btn: 'Compare these cell lines with the rest' })}
            ${pca ? shot('pca_mitf', 'PCA of the melanoma cell lines on mRNA, colored by MITF expression, with loadings', 'pca', { cap: `PCA of the same 50 cell lines on mRNA, colored by ${gi('MITF')} expression. The loadings name the genes behind the axes.`, btn: 'Open this PCA' }) : ''}
        </div>
        ${heatmap ? `<div class="tour-wide-shot">
            ${shot('heatmap_tp53', 'Heatmap of the TP53 pathway members across the 50 BRAF-mutated melanoma cell lines, with a TP53 mutation row', 'heatmap', { cap: `Heatmap: each column is one of the same 50 cell lines, each row a gene of the ${gi('TP53')} pathway. Cell lines with wild-type ${gi('TP53')} depend on ${gi('MDM2')} and ${gi('USP7')}; the ${gi('TP53')}-mutated ones, on the right, do not.`, btn: 'Open this heatmap' })}
        </div>` : ''}
        <div class="tour-feature">
            <div class="tour-text">
                <h3>The Wiki: a comprehensive summary of every cell line</h3>
                ${checks([
                    'Many sections, gathered from DepMap, Cellosaurus and other sources: origin and history, classification, mutations, copy number, dependencies, drug response and more.',
                    'Which pathways are altered, pointing to the ones involved in transforming the cell line.',
                    'Suggested antigens for FACS: surface proteins the cell line expresses strongly and fairly uniquely.',
                    'Open it from the Cell Line Browser, or Shift-click a cell line in any plot.',
                ])}
            </div>
            ${shot('wiki_a375', 'The Wiki for A375: overview, history and cancer classification', 'wiki', { cap: 'The start of the A375 Wiki: an overview, its history and derived cell lines, and how it is classified.', btn: 'Open the A375 Wiki' })}
        </div>
    </div>
    <div class="tour-wrap tour-next-row">${nextBtn('tour-gene', 'Gene Browser')}</div>
</section>

<section class="tour-sec tour-alt" id="tour-gene">
    <div class="tour-wrap">
        <div class="tour-feature tour-flip">
            <div class="tour-text">
                <p class="tour-eyebrow">4 &middot; Gene Browser</p>
                <h2>A comprehensive summary of every gene</h2>
                <p>Type a gene and get a summary over many sections: how much cells need it and which cancers need it most, how it is expressed, which genes follow it, which complexes and pathways it belongs to, its mutations, structure and sequence. Save it as a PDF.</p>
                ${checks([
                    'A verdict for dependency and for expression, ranked against every gene.',
                    'Co-dependencies and co-expression, with complex partners ringed.',
                    'Mutations, fusions, copy number, drug response, structure and sequence.',
                ])}
            </div>
            ${shot('gb_myc', 'The Gene Browser page for MYC: description, chromosome position and records', 'myc', { cap: `${gi('MYC')} in the Gene Browser: description, chromosome position and the records it links to.`, btn: 'Open MYC in the Gene Browser' })}
        </div>
        <div class="tour-feature">
            <div class="tour-text">
                <h3>Which genes work with it?</h3>
                ${checks([
                    'Co-dependencies: the genes whose knockout affects the same cell lines, with known complex and pathway partners marked.',
                    'Protein complexes from the EBI Complex Portal, showing how closely each member follows the gene.',
                ])}
            </div>
            ${shot('gb_codep', 'Co-dependencies of TSC1 as a table and a fan figure', 'tsc1', { cap: `Co-dependencies of ${gi('TSC1')}, complex partners ringed.`, btn: 'Open TSC1 in the Gene Browser' })}
        </div>
        <div class="tour-text tour-text-wide" style="margin-top:56px;">
            <h3>Exons and sequence</h3>
            ${checks([
                'The exon structure of the reference transcript (MANE Select) or any other transcript, with mutation hotspots marked.',
                'A sequence viewer for the gene, mRNA, coding sequence or protein, with search and the exon junctions marked.',
                'Download GenBank or FASTA files that open in SnapGene and other sequence editors.',
            ])}
        </div>
        <div class="tour-wide-shot" style="margin-top:0;">
            ${shot('gb_seq', 'The TP53 exon structure and its genomic sequence in the sequence viewer', 'tp53seq', { cap: `${gi('TP53')}: the exons of its MANE Select transcript, and the gene in the sequence viewer, exons in green between the introns and the start codon marked.`, btn: 'Open TP53 in the Gene Browser' })}
        </div>
    </div>
    <div class="tour-wrap tour-next-row">${nextBtn('tour-more', 'Other features')}</div>
</section>

<section class="tour-sec" id="tour-more">
    <div class="tour-wrap">
        <p class="tour-eyebrow">Other features</p>
        <p class="tour-lead-note">Found under <b>Other</b> in the app: ${otherTools.slice(0, -1).join(', ')} and ${otherTools[otherTools.length - 1]}.</p>
        <div class="tour-feature">
            <div class="tour-text">
                <h3>Mutation analysis</h3>
                <p>Pick an alteration from curated lists of cancer-relevant hotspot mutations, functional losses and amplifications, or any fusion. The cell lines that carry it are compared with those that do not, and every gene is ranked by the difference in gene effect or mRNA. Any gene in the table opens as a strip plot.</p>
            </div>
            ${shot('strip_braf', 'BRAF gene effect in skin cell lines by BRAF mutation status', 'mutation', { cap: `${gi('BRAF')} gene effect in skin cell lines, split by ${gi('BRAF')} mutation: the mutated cell lines depend on it, the wild-type ones do not.`, btn: 'Open BRAF mutation analysis in skin' })}
        </div>
        <div class="tour-feature tour-flip">
            <div class="tour-text">
                <h3>Export images, tables and methods</h3>
                <p>Under every view: Export image (PNG, SVG, PDF, PowerPoint at print size) and CSV. An exported PNG or SVG reopens in Correlate with its analysis.</p>
                <p><b>Methods</b> writes the methods for the view you are looking at, from its actual settings: the DepMap release, which cell lines were included, the filters and the statistics used. You get one paragraph ready to paste into a materials and methods section, and a longer plain-language version that explains what the figure shows. Copy either or download them as a text file.</p>
                <p><b>Export for AI</b> writes the data behind the view into one file made for a language model. Upload it to ChatGPT, Claude or Gemini and discuss the results: which of these dependencies has a drug, which cell line fits my experiment, what links these genes?</p>
            </div>
            ${shot('export_dialog', 'The Export image dialog', null, { cls: 'tour-shot-dialog', cap: 'Export image: what to export, format, size in cm, resolution, background.' })}
        </div>
        <button type="button" class="tour-close-end" id="tourCloseEnd">Close the tour and start using Correlate</button>
        <p class="tour-foot">Data: DepMap 26Q1 (Broad Institute), with EBI Complex Portal, Reactome, PRISM and the other sources listed under Reference Data. Correlate is developed in the Wermeling lab.</p>
    </div>
</section>`;
    }

    // ------------------------------------------------------------ style
    const CSS = `
#tourPlotHint { position: fixed; bottom: 80px; left: 40px; z-index: 1350; width: 300px; }
#tourNetHint { position: absolute; z-index: 1200; width: 250px; }
#tourNetHint.tour-hint-beside::after { content: ''; position: absolute; right: -9px; top: 22px; width: 14px; height: 14px; background: #fff; border-top: 1px solid #b9dd9a; border-right: 1px solid #b9dd9a; transform: rotate(45deg); }
#tourNetHint, #tourPlotHint { background: rgba(255,255,255,.96); border: 1px solid #b9dd9a; border-left: 4px solid #5d9239; border-radius: 8px; box-shadow: 0 6px 20px rgba(17,24,39,.14); padding: 10px 30px 8px 12px; font-size: 13px; line-height: 1.45; color: #374151; }
#tourNetHint b, #tourPlotHint b { display: block; font-size: 13.5px; color: #2f4d1a; margin-bottom: 4px; }
#tourNetHint ul, #tourPlotHint ul { margin: 0; padding: 0 0 0 16px; }
#tourNetHint li, #tourPlotHint li { margin: 0 0 3px; }
#tourNetHint li span, #tourPlotHint li span { font-weight: 600; color: #111827; }
#tourNetHint button, #tourPlotHint button { position: absolute; top: 4px; right: 6px; border: 0; background: none; font-size: 18px; line-height: 1; color: #6b7280; cursor: pointer; padding: 2px 4px; }
#tourNetHint button:hover, #tourPlotHint button:hover { color: #111827; }
html.popout-open #tourNetHint { display: none; }
#tourPage { position: fixed; inset: 0; z-index: 1500; overflow-y: auto; overflow-x: hidden; background: #fff; color: #111827; display: none; -webkit-overflow-scrolling: touch; }
#tourPage.active { display: block; }
#tourPage * { box-sizing: border-box; }
#tourPage .tour-wrap { max-width: 1120px; margin: 0 auto; padding: 0 24px; }
#tourPage .tour-top { position: sticky; top: 0; z-index: 5; background: rgba(255,255,255,.94); backdrop-filter: blur(8px); border-bottom: 1px solid #e5e7eb; }
#tourPage .tour-top-in { display: flex; align-items: center; gap: 18px; height: 54px; }
#tourPage .tour-brand { display: flex; align-items: baseline; gap: 10px; white-space: nowrap; }
#tourPage .tour-brand b { font-size: 20px; color: #3b5e23; font-weight: 800; letter-spacing: -.01em; }
#tourPage .tour-brand span { font-size: 13px; color: #6b7280; }
#tourPage .tour-nav { display: flex; gap: 2px; margin-left: auto; overflow-x: auto; scrollbar-width: none; }
#tourPage .tour-nav::-webkit-scrollbar { display: none; }
#tourPage .tour-nav a { font-size: 13.5px; color: #374151; text-decoration: none; padding: 7px 10px; border-radius: 6px; white-space: nowrap; }
#tourPage .tour-nav a:hover { background: #f3f4f6; }
#tourPage .tour-nav a.on { color: #3b5e23; background: #eef6e6; font-weight: 600; }
#tourPage .tour-close { flex: none; display: inline-flex; align-items: center; gap: 6px; height: 36px; padding: 0 14px; border: 1px solid #5d9239; background: #5d9239; border-radius: 999px; font: inherit; font-size: 14px; font-weight: 700; color: #fff; cursor: pointer; }
#tourPage .tour-close span { font-size: 20px; line-height: 1; font-weight: 400; }
#tourPage .tour-close:hover { background: #4f7f30; border-color: #4f7f30; }
#tourPage .tour-close-end { display: block; margin: 48px auto 22px; font: inherit; font-size: 16px; font-weight: 700; color: #fff; background: #374151; border: 0; border-radius: 10px; padding: 13px 28px; cursor: pointer; }
#tourPage .tour-close-end:hover { background: #111827; }
#tourPage .tour-eyebrow { font-size: 12px; font-weight: 700; letter-spacing: .12em; text-transform: uppercase; color: #5d9239; margin: 0 0 10px; }
#tourPage .tour-hero { background: linear-gradient(160deg, #1f3a12 0%, #3b5e23 55%, #4f7f30 100%); color: #fff; padding: 48px 0 26px; position: relative; display: flow-root; }
#tourPage .tour-hero .tour-wrap { max-width: 1280px; }
#tourPage .tour-hero-grid { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 500px); gap: 56px; align-items: center; }
#tourPage .tour-hero .tour-eyebrow { color: #b9dd9a; }
#tourPage .tour-hero h1 { font-size: 31px; line-height: 1.25; font-weight: 700; letter-spacing: -.01em; margin: 0 0 16px; }
#tourPage .tour-sub { font-size: 19px; line-height: 1.6; color: #e3efd8; margin: 0 0 4px; }
#tourPage .tour-hero-list { list-style: none; margin: 0; padding: 0; }
#tourPage .tour-hero-list li { font-size: 18px; line-height: 1.5; color: #fff; padding: 4px 0 4px 26px; position: relative; }
#tourPage .tour-hero-list li::before { content: ''; position: absolute; left: 4px; top: 13px; width: 9px; height: 9px; border-radius: 50%; background: #b9dd9a; }
/* The three tool lines carry the app's own buttons, so the text columns line
   up on one edge and the dot is dropped. */
#tourPage .tour-hero-list li.tour-hero-fn { display: grid; grid-template-columns: 124px minmax(0, 1fr); gap: 14px; align-items: center; padding: 6px 0 6px 0; }
#tourPage .tour-hero-list li.tour-hero-fn::before { display: none; }
#tourPage .tour-modebtn { font: inherit; font-size: 13px; font-weight: 600; line-height: 1.2; text-align: center; color: #3b5e23; background: #fff; border: 1px solid #b9dd9a; border-radius: 6px; padding: 7px 4px; width: 100%; cursor: pointer; box-shadow: 0 1px 2px rgba(0,0,0,.06); }
#tourPage .tour-modebtn:hover { background: #eef6e6; }
#tourPage .tour-cta { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 18px; }
#tourPage .tour-try { font: inherit; font-size: 15px; font-weight: 700; color: #fff; background: #5d9239; border: 1px solid #5d9239; border-radius: 8px; padding: 11px 18px; cursor: pointer; box-shadow: 0 1px 2px rgba(0,0,0,.08); }
#tourPage .tour-try:hover { background: #4f7f30; border-color: #4f7f30; }
#tourPage .tour-hero .tour-try { background: #fff; color: #2f4d1a; border-color: #fff; }
#tourPage .tour-hero .tour-try:hover { background: #eef6e6; }
#tourPage .tour-stats { background: #fff; padding: 40px 0 10px; }
#tourPage .tour-stats-in { display: grid; grid-template-columns: repeat(4, 1fr); gap: 18px; border-top: 1px solid #e5e7eb; border-bottom: 1px solid #e5e7eb; padding: 22px 0; }
#tourPage .tour-stats-in div { text-align: center; }
#tourPage .tour-stats-in b { display: block; font-size: 32px; font-weight: 800; color: #3b5e23; letter-spacing: -.02em; line-height: 1.1; }
#tourPage .tour-stats-in span { font-size: 13px; color: #6b7280; }
#tourPage .tour-sec { padding: 72px 0; }
#tourPage .tour-alt { background: #f6faf2; }
#tourPage .tour-feature { display: grid; grid-template-columns: minmax(0, 5fr) minmax(0, 7fr); gap: 48px; align-items: center; }
#tourPage .tour-feature + .tour-feature, #tourPage .tour-two + .tour-feature, #tourPage .tour-feature + .tour-two, #tourPage .tour-two + .tour-wide-shot, #tourPage .tour-wide-shot + .tour-feature { margin-top: 56px; }
#tourPage .tour-flip .tour-text { order: 2; }
#tourPage .tour-flip-wide { grid-template-columns: minmax(0, 7fr) minmax(0, 5fr); }
#tourPage .tour-text h2 { font-size: 32px; line-height: 1.15; font-weight: 800; letter-spacing: -.02em; margin: 0 0 14px; color: #111827; }
#tourPage .tour-text h3 { font-size: 21px; font-weight: 800; margin: 0 0 12px; color: #111827; letter-spacing: -.01em; }
#tourPage .tour-text p { font-size: 16.5px; line-height: 1.6; color: #374151; margin: 0 0 14px; }
#tourPage .tour-lead-note { font-size: 16px; line-height: 1.6; color: #374151; max-width: 820px; margin: 0 0 36px; }
#tourPage .tour-text-wide { max-width: 820px; margin-bottom: 30px; }
#tourPage .tour-checks { list-style: none; margin: 0; padding: 0; }
#tourPage .tour-checks li { position: relative; padding-left: 28px; font-size: 15px; line-height: 1.55; color: #374151; margin: 0 0 10px; }
#tourPage .tour-checks li::before { content: ''; position: absolute; left: 0; top: 4px; width: 18px; height: 18px; border-radius: 50%; background: #eef6e6; border: 1px solid #b9dd9a; }
#tourPage .tour-checks li::after { content: ''; position: absolute; left: 6px; top: 8px; width: 5px; height: 8px; border: solid #3b5e23; border-width: 0 2px 2px 0; transform: rotate(45deg); }
#tourPage .tour-shot { margin: 0; background: #fff; border: 1px solid #d9dee5; border-radius: 12px; overflow: hidden; box-shadow: 0 14px 40px rgba(17,24,39,.16), 0 2px 6px rgba(17,24,39,.06); }
#tourPage .tour-shot-bar { height: 22px; background: #f3f4f6; border-bottom: 1px solid #e5e7eb; display: flex; align-items: center; gap: 5px; padding: 0 10px; }
#tourPage .tour-shot-bar span { width: 8px; height: 8px; border-radius: 50%; background: #d1d5db; }
#tourPage .tour-shot img { display: block; width: 100%; height: auto; }
#tourPage .tour-shot figcaption { font-size: 13px; line-height: 1.45; color: #6b7280; padding: 9px 12px 10px; border-top: 1px solid #eef0f3; background: #fff; }
#tourPage .tour-shot-act { padding: 2px 12px 12px; background: #fff; }
#tourPage .tour-shot-act .tour-try { width: 100%; }
#tourPage .tour-hero-shot .tour-shot { box-shadow: 0 30px 70px rgba(0,0,0,.35); }
#tourPage .tour-hero .tour-shot-act .tour-try { background: #5d9239; color: #fff; border-color: #5d9239; }
#tourPage .tour-hero .tour-shot-act .tour-try:hover { background: #4f7f30; border-color: #4f7f30; }
#tourPage .tour-next-row { display: flex; justify-content: center; }
#tourPage .tour-next { font: inherit; font-size: 15px; font-weight: 700; color: #3b5e23; background: #fff; border: 1px solid #b9dd9a; border-radius: 999px; padding: 9px 22px; margin-top: 48px; cursor: pointer; box-shadow: 0 1px 2px rgba(0,0,0,.06); }
#tourPage .tour-next:hover { background: #eef6e6; }
#tourPage .tour-hero .tour-next { margin-top: 26px; color: #fff; background: rgba(255,255,255,.12); border-color: rgba(255,255,255,.45); }
#tourPage .tour-hero .tour-next:hover { background: rgba(255,255,255,.22); }
#tourPage .tour-stack { display: flex; flex-direction: column; gap: 18px; }
#tourPage .tour-shot-dialog { max-width: 420px; margin: 0 auto; }
#tourPage .tour-two { display: grid; grid-template-columns: 1fr 1fr; gap: 36px; align-items: start; }
#tourPage .tour-two-tight { gap: 24px; }
/* Paired figures share one frame height and keep their buttons on one
   line, whatever the pictures' own proportions. */
#tourPage .tour-two .tour-shot { display: flex; flex-direction: column; }
#tourPage .tour-two .tour-shot img { aspect-ratio: 5 / 6; object-fit: contain; padding: 6px; }
#tourPage .tour-two .tour-shot figcaption { flex: 1; }
#tourPage .tour-two-net { align-items: stretch; }
#tourPage .tour-two-pca { grid-template-columns: minmax(0, 99fr) minmax(0, 133fr); align-items: stretch; }
#tourPage .tour-two-pca .tour-shot img { aspect-ratio: auto; }
#tourPage .tour-two-net .tour-shot img { aspect-ratio: 1 / 1; }
#tourPage .tour-two-one { grid-template-columns: minmax(0, 560px); justify-content: center; }
#tourPage .tour-stack .tour-shot img { aspect-ratio: 4 / 3; object-fit: contain; padding: 6px; }
#tourPage .tour-card-fig h3 { font-size: 19px; font-weight: 800; margin: 18px 0 8px; color: #111827; }
#tourPage .tour-card-fig p { font-size: 15px; line-height: 1.55; color: #374151; margin: 0 0 14px; }
#tourPage .tour-wide-shot { margin-top: 40px; }
#tourPage .tour-edge-shot { max-width: 920px; margin-left: auto; margin-right: auto; }
#tourPage .tour-cards { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; margin-top: 6px; }
#tourPage .tour-card { background: #fff; border: 1px solid #e5e7eb; border-radius: 10px; padding: 14px 16px; }
#tourPage .tour-card b { display: block; font-size: 15px; color: #3b5e23; margin-bottom: 4px; }
#tourPage .tour-card span { font-size: 13.5px; line-height: 1.5; color: #4b5563; }
#tourPage .tour-more { padding: 56px 0 40px; border-top: 1px solid #e5e7eb; }
#tourPage .tour-foot { text-align: center; font-size: 12.5px; color: #9ca3af; margin: 10px 0 0; }
@media (max-width: 900px) {
  #tourPage .tour-feature, #tourPage .tour-two, #tourPage .tour-cards, #tourPage .tour-hero-grid { grid-template-columns: 1fr; }
  #tourPage .tour-feature { gap: 26px; }
  #tourPage .tour-flip .tour-text { order: 0; }
  #tourPage .tour-stats-in { grid-template-columns: 1fr 1fr; }
  #tourPage .tour-hero h1 { font-size: 26px; }
  #tourPage .tour-text h2 { font-size: 26px; }
  #tourPage .tour-sec { padding: 48px 0; }
}
@media (max-width: 640px) {
  #tourPage .tour-wrap { padding: 0 16px; }
  #tourPage .tour-brand span, #tourPage .tour-nav { display: none; }
  #tourPage .tour-close { margin-left: auto; }
  #tourPage .tour-hero { padding-top: 40px; }
  #tourPage .tour-hero h1 { font-size: 22px; }
  #tourPage .tour-sub { font-size: 16px; }
  #tourPage .tour-hero-list li.tour-hero-fn { grid-template-columns: minmax(0, 1fr); gap: 6px; }
  #tourPage .tour-modebtn { width: auto; justify-self: start; padding: 6px 12px; }
  #tourPage .tour-stats-in b { font-size: 24px; }
  #tourPage .tour-try { width: 100%; text-align: center; min-height: 44px; }
}`;

    // ------------------------------------------------------------ shell
    let built = false, io = null;
    function build() {
        if (built) return;
        built = true;
        const style = document.createElement('style');
        style.textContent = CSS;
        document.head.appendChild(style);
        const page = document.createElement('div');
        page.id = 'tourPage';
        page.setAttribute('role', 'dialog');
        page.setAttribute('aria-label', 'Tour of the app');
        document.body.appendChild(page);
        // A deliberate scroll of the reader's own wins over a running glide.
        // A trackpad keeps sending small wheel events after a gesture and
        // while a finger rests on it, so only a real push counts, and not one
        // in the same moment as the click: otherwise Next looks dead.
        const stopGlide = (e) => {
            if (e.type === 'wheel' && Math.abs(e.deltaY || 0) < 40) return;
            if (performance.now() - glideT0 < 400) return;
            cancelAnimationFrame(scrollRaf);
        };
        page.addEventListener('wheel', stopGlide, { passive: true });
        page.addEventListener('touchmove', stopGlide, { passive: true });
        page.addEventListener('click', (e) => {
            const act = e.target.closest('[data-act]');
            if (act) { e.preventDefault(); run(act.getAttribute('data-act')); return; }
            const link = e.target.closest('.tour-nav a, [data-next]');
            if (link) {
                e.preventDefault();
                const t = page.querySelector(link.getAttribute('href') || link.getAttribute('data-next'));
                // Land with the section's heading at the top, just under the
                // sticky bar, not with the section's empty top padding.
                const head = t?.querySelector('.tour-eyebrow, h2, h3') || t;
                const bar = page.querySelector('.tour-top')?.offsetHeight || 54;
                if (!head) return;
                const top = page.scrollTop + head.getBoundingClientRect().top - page.getBoundingClientRect().top - bar - 16;
                const dur = Math.min(1500, Math.max(600, Math.abs(top - page.scrollTop) * 0.6));
                sectionFx(link.closest('section'), t, link.getAttribute('data-fx'), dur);
                glideTo(page, top, dur);
            }
        });
        document.addEventListener('keydown', (e) => { if (e.key === 'Escape' && isOpen()) close(); });
    }
    // The browser's own smooth scroll crosses a long page in a few hundred ms,
    // which reads as a jump, so the glide is run by hand over a duration that
    // grows with the distance and stops as soon as the reader scrolls.
    let scrollRaf = 0, glideT0 = 0;
    function glideTo(page, top, dur) {
        cancelAnimationFrame(scrollRaf);
        const from = page.scrollTop, d = top - from, t0 = performance.now();
        glideT0 = t0;
        if (!d) return;
        if (reduced()) { page.scrollTop = top; return; }
        const step = (now) => {
            const p = Math.min(1, (now - t0) / dur);
            const e = p < 0.5 ? 4 * p * p * p : 1 - Math.pow(-2 * p + 2, 3) / 2;
            page.scrollTop = from + d * e;
            if (p < 1) scrollRaf = requestAnimationFrame(step);
        };
        scrollRaf = requestAnimationFrame(step);
    }
    const reduced = () => !!window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

    // Two of the Next buttons turn the page rather than just scrolling it.
    // The animations carry no fill, so the sections keep their plain styling
    // once they are over and nothing is left behind to clean up.
    function sectionFx(from, to, kind, dur) {
        if (!kind || !from || !to || reduced() || !from.animate) return;
        const a = [...from.querySelectorAll(':scope > .tour-wrap')];
        const b = [...to.querySelectorAll(':scope > .tour-wrap')];
        if (!a.length || !b.length) return;
        const play = (els, frames, opt, origin) => els.map(el => {
            if (origin) el.style.transformOrigin = origin;
            return el.animate(frames, opt);
        });
        if (kind === 'cube') {
            from.style.perspective = to.style.perspective = '1500px';
            const done = () => { from.style.perspective = to.style.perspective = ''; };
            play(a, [{ transform: 'rotateX(0deg)', opacity: 1 }, { transform: 'rotateX(-72deg)', opacity: .1 }],
                 { duration: Math.round(dur * .75), easing: 'cubic-bezier(.55,0,.75,.35)' }, '50% 100%');
            play(b, [{ transform: 'rotateX(72deg)', opacity: .1 }, { transform: 'rotateX(0deg)', opacity: 1 }],
                 { duration: dur, easing: 'cubic-bezier(.2,.7,.25,1)' }, '50% 0%').forEach(an => { an.onfinish = done; an.oncancel = done; });
        } else if (kind === 'zoom') {
            play(a, [{ transform: 'scale(1)', opacity: 1 }, { transform: 'scale(1.07)', opacity: 0 }],
                 { duration: Math.round(dur * .7), easing: 'cubic-bezier(.5,0,.8,.4)' });
            play(b, [{ transform: 'scale(.92)', opacity: 0, filter: 'blur(6px)' }, { transform: 'scale(1)', opacity: 1, filter: 'blur(0px)' }],
                 { duration: dur, easing: 'cubic-bezier(.2,.7,.25,1)' });
        }
    }

    const isOpen = () => !!document.getElementById('tourPage')?.classList.contains('active');

    function render() {
        const page = document.getElementById('tourPage');
        page.innerHTML = pageHTML();
        page.querySelector('#tourClose').addEventListener('click', close);
        page.querySelector('#tourCloseEnd')?.addEventListener('click', close);
        // Section links follow the scroll.
        if (io) io.disconnect();
        const links = [...page.querySelectorAll('.tour-nav a')];
        io = new IntersectionObserver((entries) => {
            entries.forEach(en => {
                if (!en.isIntersecting) return;
                links.forEach(l => l.classList.toggle('on', l.getAttribute('href') === '#' + en.target.id));
            });
        }, { root: page, rootMargin: '-40% 0px -55% 0px', threshold: 0 });
        links.forEach(l => { const t = page.querySelector(l.getAttribute('href')); if (t) io.observe(t); });
    }

    function run(name) {
        const fn = ACTIONS[name];
        if (!fn || !A()) return;
        close();
        try { const r = fn(A()); if (r && r.catch) r.catch(e => console.warn('Tour action failed:', e)); }
        catch (e) { console.warn('Tour action failed:', e); }
    }

    function open() {
        document.getElementById('tourNetHint')?.remove();
        document.getElementById('tourPlotHint')?.remove();
        build();
        render();
        const page = document.getElementById('tourPage');
        page.classList.add('active');
        page.scrollTop = 0;
    }
    function close() {
        document.getElementById('tourPage')?.classList.remove('active');
    }

    window.CorrelateTour = { open, close, isOpen };
})();
