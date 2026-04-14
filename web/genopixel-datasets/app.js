const FACETS = [
  { key: 'tissue', label: 'Tissue', field: 'tissues', type: 'array' },
  { key: 'disease', label: 'Disease', field: 'diseases', type: 'array' },
  { key: 'organism', label: 'Organism', field: 'organisms', type: 'array' },
  { key: 'tissue_type', label: 'Tissue type', field: 'tissue_types', type: 'array' },
  { key: 'journal', label: 'Journal', field: 'journal', type: 'scalar' },
];
const ACTIVE_DATASET_SYNC_INTERVAL_MS = 30_000;

const state = {
  catalog: null,
  datasets: [],
  filters: Object.fromEntries(FACETS.map((facet) => [facet.key, new Set()])),
  facetSearch: { tissue: '', disease: '', journal: '' },
  ranges: { yearMin: '', yearMax: '', cellMin: '', cellMax: '' },
  expandedFacets: new Set(),
  search: '',
  sort: 'smart',
  page: 1,
  pageSize: 24,
  selectedRow: null,
  selectedVariantRow: null,
  detailCache: new Map(),
  detailLoading: false,
  detailError: '',
  detailActionLoading: false,
  detailActionMessage: '',
  detailActionTone: 'success',
  activeDataset: null,
  activeDatasetLastSyncAt: '',
  activeDatasetSyncError: '',
  activeDatasetSyncInFlight: false,
};
let activeDatasetSyncTimer = null;

const elements = {
  searchInput: document.querySelector('#searchInput'),
  sortSelect: document.querySelector('#sortSelect'),
  clearAllButton: document.querySelector('#clearAllButton'),
  activeDatasetBanner: document.querySelector('#activeDatasetBanner'),
  resultSummary: document.querySelector('#resultSummary'),
  activeChips: document.querySelector('#activeChips'),
  facetSectionsTop: document.querySelector('#facetSectionsTop'),
  facetSectionsBottom: document.querySelector('#facetSectionsBottom'),
  resultsPanel: document.querySelector('.results-panel'),
  resultsContainer: document.querySelector('#resultsContainer'),
  pagination: document.querySelector('#pagination'),
  statusBanner: document.querySelector('#statusBanner'),
  detailTitle: document.querySelector('#detailTitle'),
  detailBody: document.querySelector('#detailBody'),
  detailDrawer: document.querySelector('#detailDrawer'),
  detailBackdrop: document.querySelector('#detailBackdrop'),
  closeDetailButton: document.querySelector('#closeDetailButton'),
  openFiltersButton: document.querySelector('#openFiltersButton'),
  closeFiltersButton: document.querySelector('#closeFiltersButton'),
  yearMin: document.querySelector('#yearMin'),
  yearMax: document.querySelector('#yearMax'),
};

function formatNumber(value) {
  if (value === null || value === undefined || value === '') {
    return '—';
  }
  return new Intl.NumberFormat('en-US').format(Number(value));
}

function escapeHtml(value) {
  return String(value ?? '')
    .replaceAll('&', '&amp;')
    .replaceAll('<', '&lt;')
    .replaceAll('>', '&gt;')
    .replaceAll('"', '&quot;')
    .replaceAll("'", '&#39;');
}

function listPreview(values, limit = 2) {
  if (!Array.isArray(values) || values.length === 0) {
    return '—';
  }
  if (values.length <= limit) {
    return values.join(', ');
  }
  return `${values.slice(0, limit).join(', ')} +${values.length - limit}`;
}

function tokenize(query) {
  return query
    .trim()
    .toLowerCase()
    .split(/\s+/)
    .filter(Boolean);
}

function getActiveDatasetRow() {
  const active = state.activeDataset;
  if (!active || !active.loaded) {
    return null;
  }
  const row = Number(active.all_excel_row);
  return Number.isFinite(row) ? row : null;
}

function getActiveVariantRow() {
  const active = state.activeDataset;
  if (!active || !active.loaded) {
    return null;
  }
  if (active.multiple_excel_row === null || active.multiple_excel_row === undefined || active.multiple_excel_row === '') {
    return null;
  }
  const row = Number(active.multiple_excel_row);
  return Number.isFinite(row) ? row : null;
}

function isRuntimeActiveDatasetRow(allExcelRow) {
  return getActiveDatasetRow() === allExcelRow;
}

function getDatasetByRow(allExcelRow) {
  return state.datasets.find((dataset) => dataset.all_excel_row === allExcelRow) || null;
}

function clearSearchFacetAndYearFilters() {
  state.search = '';
  FACETS.forEach((facet) => state.filters[facet.key].clear());
  state.facetSearch.tissue = '';
  state.facetSearch.disease = '';
  state.facetSearch.journal = '';
  state.ranges.yearMin = '';
  state.ranges.yearMax = '';
  elements.searchInput.value = '';
  elements.yearMin.value = '';
  elements.yearMax.value = '';
}

function getActiveDatasetVisibility(filteredDatasets) {
  const activeRow = getActiveDatasetRow();
  if (activeRow === null) {
    return {
      hasActiveDataset: false,
      activeRow: null,
      inCatalog: false,
      inFiltered: false,
      inCurrentPage: false,
      pageForActive: null,
    };
  }

  const inCatalog = state.datasets.some((dataset) => dataset.all_excel_row === activeRow);
  const filteredIndex = filteredDatasets.findIndex((dataset) => dataset.all_excel_row === activeRow);
  const inFiltered = filteredIndex >= 0;
  const pageForActive = inFiltered ? Math.floor(filteredIndex / state.pageSize) + 1 : null;
  const pageStart = (state.page - 1) * state.pageSize;
  const inCurrentPage = inFiltered && filteredIndex >= pageStart && filteredIndex < pageStart + state.pageSize;

  return {
    hasActiveDataset: true,
    activeRow,
    inCatalog,
    inFiltered,
    inCurrentPage,
    pageForActive,
  };
}

async function fetchActiveDatasetState() {
  const response = await fetch('/api/genopixel-runtime/active-dataset');
  if (!response.ok) {
    throw new Error(`Active dataset request failed with status ${response.status}.`);
  }
  return response.json();
}

async function refreshActiveDataset({ silent = true, render: shouldRender = true } = {}) {
  if (state.activeDatasetSyncInFlight) {
    return;
  }
  state.activeDatasetSyncInFlight = true;
  try {
    const payload = await fetchActiveDatasetState();
    if (payload.loaded) {
      state.activeDataset = payload;
      state.activeDatasetLastSyncAt = new Date().toISOString();
      state.activeDatasetSyncError = '';
    } else {
      state.activeDataset = null;
      state.activeDatasetLastSyncAt = new Date().toISOString();
      state.activeDatasetSyncError = '';
    }
  } catch (error) {
    if (state.activeDataset) {
      state.activeDatasetSyncError = 'Live active-dataset sync is unavailable. Showing last known active dataset.';
    } else if (!silent) {
      setStatus(error.message || 'Failed to refresh active dataset state.', 'error');
    }
  } finally {
    state.activeDatasetSyncInFlight = false;
    if (shouldRender) {
      render();
    }
  }
}

function startActiveDatasetSync() {
  if (activeDatasetSyncTimer !== null) {
    window.clearInterval(activeDatasetSyncTimer);
  }
  activeDatasetSyncTimer = window.setInterval(() => {
    if (document.visibilityState !== 'visible') {
      return;
    }
    refreshActiveDataset({ silent: true });
  }, ACTIVE_DATASET_SYNC_INTERVAL_MS);
}

function buildDatasetIndex(dataset) {
  const titleSearch = (dataset.title || '').toLowerCase();
  const authorJournalSearch = [dataset.author, dataset.journal].filter(Boolean).join(' ').toLowerCase();
  const searchableChunks = [
    dataset.title,
    dataset.author,
    dataset.journal,
    dataset.doi,
    dataset.cellxgene_doi,
    dataset.project,
    ...(dataset.tissues || []),
    ...(dataset.diseases || []),
    ...(dataset.organisms || []),
    ...(dataset.tissue_types || []),
  ].filter(Boolean);

  return {
    ...dataset,
    _titleSearch: titleSearch,
    _authorJournalSearch: authorJournalSearch,
    _searchBlob: searchableChunks.join(' ').toLowerCase(),
    _facetValues: {
      project: dataset.project ? [dataset.project.toLowerCase()] : [],
      organism: (dataset.organisms || []).map((value) => value.toLowerCase()),
      tissue: (dataset.tissues || []).map((value) => value.toLowerCase()),
      tissue_type: (dataset.tissue_types || []).map((value) => value.toLowerCase()),
      disease: (dataset.diseases || []).map((value) => value.toLowerCase()),
      journal: dataset.journal ? [dataset.journal.toLowerCase()] : [],
      merged: dataset.merged ? [dataset.merged.toLowerCase()] : [],
    },
  };
}

function computeSearchTier(dataset, tokens) {
  if (!tokens.length) {
    return 0;
  }
  if (!tokens.every((token) => dataset._searchBlob.includes(token))) {
    return -1;
  }
  if (tokens.every((token) => dataset._titleSearch.includes(token))) {
    return 3;
  }
  if (tokens.every((token) => dataset._authorJournalSearch.includes(token))) {
    return 2;
  }
  return 1;
}

function normalizeSelection(values) {
  return [...values].map((value) => value.toLowerCase());
}

function datasetPassesFacet(dataset, facetKey) {
  const selected = state.filters[facetKey];
  if (!selected || selected.size === 0) {
    return true;
  }
  const selectedValues = normalizeSelection(selected);
  const datasetValues = dataset._facetValues[facetKey] || [];
  return selectedValues.some((value) => datasetValues.includes(value));
}

function datasetPassesRanges(dataset) {
  const yearValue = Number(dataset.year || 0);
  const { yearMin, yearMax } = state.ranges;

  if (yearMin !== '' && yearValue < Number(yearMin)) {
    return false;
  }
  if (yearMax !== '' && yearValue > Number(yearMax)) {
    return false;
  }
  return true;
}

function datasetPassesSearch(dataset, tokens) {
  return computeSearchTier(dataset, tokens) >= 0;
}

function datasetPassesAllFilters(dataset, tokens, ignoredFacet = '') {
  if (!datasetPassesRanges(dataset)) {
    return false;
  }
  if (!datasetPassesSearch(dataset, tokens)) {
    return false;
  }
  return FACETS.every((facet) => {
    if (facet.key === ignoredFacet) {
      return true;
    }
    return datasetPassesFacet(dataset, facet.key);
  });
}

function compareDatasets(left, right, tokens) {
  const queryActive = tokens.length > 0;
  const sort = state.sort;

  if (sort === 'title_asc') {
    return (left.title || '').localeCompare(right.title || '');
  }
  if (sort === 'cells_desc') {
    return (right.cell_counts || 0) - (left.cell_counts || 0) || (right.year || 0) - (left.year || 0);
  }
  if (sort === 'year_desc') {
    return (right.year || 0) - (left.year || 0) || (right.cell_counts || 0) - (left.cell_counts || 0) || (left.title || '').localeCompare(right.title || '');
  }

  if (queryActive) {
    const leftTier = computeSearchTier(left, tokens);
    const rightTier = computeSearchTier(right, tokens);
    return rightTier - leftTier || (right.year || 0) - (left.year || 0) || (right.cell_counts || 0) - (left.cell_counts || 0) || (left.title || '').localeCompare(right.title || '');
  }

  return (right.year || 0) - (left.year || 0) || (right.cell_counts || 0) - (left.cell_counts || 0) || (left.title || '').localeCompare(right.title || '');
}

function getFilteredDatasets() {
  const tokens = tokenize(state.search);
  return state.datasets
    .filter((dataset) => datasetPassesAllFilters(dataset, tokens))
    .sort((left, right) => compareDatasets(left, right, tokens));
}

function getFacetCounts(facetKey) {
  const tokens = tokenize(state.search);
  const counts = new Map();

  state.datasets.forEach((dataset) => {
    if (!datasetPassesAllFilters(dataset, tokens, facetKey)) {
      return;
    }
    const values = dataset._facetValues[facetKey] || [];
    [...new Set(values)].forEach((value) => {
      counts.set(value, (counts.get(value) || 0) + 1);
    });
  });

  const labelLookup = new Map();
  state.datasets.forEach((dataset) => {
    const values = Array.isArray(dataset[FACETS.find((facet) => facet.key === facetKey)?.field])
      ? dataset[FACETS.find((facet) => facet.key === facetKey).field]
      : [dataset[FACETS.find((facet) => facet.key === facetKey).field]].filter(Boolean);
    values.forEach((value) => {
      labelLookup.set(value.toLowerCase(), value);
    });
  });

  return [...counts.entries()]
    .map(([value, count]) => ({ value, label: labelLookup.get(value) || value, count }))
    .sort((left, right) => right.count - left.count || left.label.localeCompare(right.label));
}

function setStatus(message = '', tone = 'default') {
  if (!message) {
    elements.statusBanner.hidden = true;
    elements.statusBanner.textContent = '';
    return;
  }
  elements.statusBanner.hidden = false;
  elements.statusBanner.dataset.tone = tone;
  elements.statusBanner.textContent = message;
}

function renderFacets() {
  const renderFacetGroup = (facetGroup) => facetGroup.map((facet) => {
    const counts = getFacetCounts(facet.key);
    const facetQuery = String(state.facetSearch[facet.key] || '').trim().toLowerCase();
    const filteredCounts = facetQuery
      ? counts.filter((entry) => entry.label.toLowerCase().includes(facetQuery))
      : counts;
    const expanded = state.expandedFacets.has(facet.key);
    const visible = facetQuery ? filteredCounts : (expanded ? filteredCounts : filteredCounts.slice(0, 10));
    const selected = state.filters[facet.key];
    return `
      <section class="facet-card">
        <div class="range-header">
          <div>
            <h3>${escapeHtml(facet.label)}</h3>
            <p>${formatNumber(counts.length)} values</p>
          </div>
        </div>
        ${['tissue', 'disease', 'journal'].includes(facet.key) ? `
          <div class="facet-search-shell">
            <input
              class="facet-search-input"
              type="search"
              data-facet-search="${facet.key}"
              value="${escapeHtml(state.facetSearch[facet.key] || '')}"
              placeholder="Type to search ${escapeHtml(facet.label.toLowerCase())}"
            />
          </div>
        ` : ''}
        <div class="facet-list">
          ${visible.length ? visible.map((entry) => `
            <div class="facet-option">
              <label>
                <input
                  type="checkbox"
                  data-facet-key="${facet.key}"
                  data-facet-value="${escapeHtml(entry.label)}"
                  ${selected.has(entry.label) ? 'checked' : ''}
                />
                <span>${escapeHtml(entry.label)}</span>
              </label>
              <span class="facet-count">${formatNumber(entry.count)}</span>
            </div>
          `).join('') : `<p class="facet-empty">No matching ${escapeHtml(facet.label.toLowerCase())} values.</p>`}
        </div>
        ${counts.length > 10 && !facetQuery ? `<button class="toggle-button" data-toggle-facet="${facet.key}" type="button">${expanded ? 'Show fewer' : 'Show all'}</button>` : ''}
      </section>
    `;
  }).join('');

  elements.facetSectionsTop.innerHTML = renderFacetGroup(FACETS.filter((facet) => facet.key !== 'journal'));
  elements.facetSectionsBottom.innerHTML = renderFacetGroup(FACETS.filter((facet) => facet.key === 'journal'));
}

function buildActiveChips() {
  const chips = [];
  if (state.search.trim()) {
    chips.push({ kind: 'search', label: `Search: ${state.search.trim()}` });
  }
  FACETS.forEach((facet) => {
    state.filters[facet.key].forEach((value) => {
      chips.push({ kind: facet.key, label: `${facet.label}: ${value}` });
    });
  });
  if (state.ranges.yearMin !== '' || state.ranges.yearMax !== '') {
    chips.push({ kind: 'range', label: `Year: ${state.ranges.yearMin || 'min'} to ${state.ranges.yearMax || 'max'}` });
  }
  return chips;
}

function renderActiveFilters(filteredDatasets) {
  elements.resultSummary.hidden = true;
  elements.resultSummary.textContent = '';

  const chips = buildActiveChips();
  elements.activeChips.innerHTML = chips.length
    ? chips.map((chip) => `<button class="chip" data-chip-kind="${chip.kind}" data-chip-label="${escapeHtml(chip.label)}" type="button">${escapeHtml(chip.label)} ×</button>`).join('')
    : '<span class="result-summary">No active filters.</span>';
}

function renderActiveDatasetBanner(filteredDatasets) {
  const activeRow = getActiveDatasetRow();
  if (activeRow === null) {
    elements.activeDatasetBanner.hidden = true;
    elements.activeDatasetBanner.innerHTML = '';
    return;
  }

  const active = state.activeDataset || {};
  const visibility = getActiveDatasetVisibility(filteredDatasets);
  const activeDataset = getDatasetByRow(activeRow);
  const activeTitle = active.title || activeDataset?.title || `Dataset ${activeRow}`;

  let tone = 'info';
  let helperText = '';
  if (!visibility.inCatalog) {
    tone = 'warn';
    helperText = 'Active dataset not found in current catalog snapshot.';
  } else if (!visibility.inFiltered) {
    tone = 'warn';
    helperText = 'Active dataset is hidden by current filters.';
  } else if (!visibility.inCurrentPage) {
    helperText = `Active dataset is on page ${visibility.pageForActive}.`;
  }

  const syncWarning = state.activeDatasetSyncError
    ? `<p class="active-dataset-warning">${escapeHtml(state.activeDatasetSyncError)}</p>`
    : '';

  elements.activeDatasetBanner.hidden = false;
  elements.activeDatasetBanner.innerHTML = `
    <section class="active-dataset-banner" data-tone="${escapeHtml(tone)}">
      <div class="active-dataset-copy">
        <p class="active-dataset-label">Runtime active dataset</p>
        <p class="active-dataset-title">${escapeHtml(activeTitle)}</p>
        ${helperText ? `<p class="active-dataset-warning">${escapeHtml(helperText)}</p>` : ''}
        ${syncWarning}
      </div>
      <button class="active-dataset-jump" data-active-dataset-action="show" type="button">Show active dataset</button>
    </section>
  `;
}

function renderResults(filteredDatasets) {
  if (!state.catalog) {
    elements.resultsContainer.innerHTML = '<div class="loading-state">Loading catalog…</div>';
    return;
  }

  if (filteredDatasets.length === 0) {
    elements.resultsContainer.innerHTML = `
      <section class="empty-state">
        <h3>No datasets match the current search and filters.</h3>
        <p>Try broadening the query, clearing a facet, or widening the year and cell-count ranges.</p>
      </section>
    `;
    return;
  }

  const start = (state.page - 1) * state.pageSize;
  const visibleRows = filteredDatasets.slice(start, start + state.pageSize);

  const desktopRows = visibleRows.map((dataset) => {
    const isSelected = state.selectedRow === dataset.all_excel_row;
    const isRuntimeActive = isRuntimeActiveDatasetRow(dataset.all_excel_row);
    const runtimeVariantRow = isRuntimeActive ? getActiveVariantRow() : null;
    const rowClasses = ['result-row'];
    if (isSelected) {
      rowClasses.push('is-active');
    }
    if (isRuntimeActive) {
      rowClasses.push('is-runtime-active');
    }
    if (isSelected && isRuntimeActive) {
      rowClasses.push('is-active-runtime');
    }

    return `
      <button class="${rowClasses.join(' ')}" data-row="${dataset.all_excel_row}" type="button">
        <div class="result-primary">
          <h3>${escapeHtml(dataset.title || 'Untitled dataset')}</h3>
          <div class="badge-row">
            ${isRuntimeActive ? '<span class="badge runtime-active-badge">Active dataset</span>' : ''}
            <span class="badge">${escapeHtml(dataset.author || 'Unknown author')}</span>
            <span class="badge">${escapeHtml(dataset.journal || 'No journal')}</span>
          </div>
        </div>
        <div class="result-cell">${escapeHtml(dataset.author || '—')}</div>
        <div class="result-cell">${escapeHtml(dataset.year || '—')}</div>
        <div class="result-cell">${escapeHtml(listPreview(dataset.tissues))}</div>
        <div class="result-cell">${escapeHtml(listPreview(dataset.diseases))}</div>
        <div class="result-cell">${escapeHtml(listPreview(dataset.organisms))}</div>
        <div class="result-cell">${formatNumber(dataset.cell_counts)}</div>
      </button>
    `;
  }).join('');

  const mobileCards = visibleRows.map((dataset) => {
    const isSelected = state.selectedRow === dataset.all_excel_row;
    const isRuntimeActive = isRuntimeActiveDatasetRow(dataset.all_excel_row);
    const runtimeVariantRow = isRuntimeActive ? getActiveVariantRow() : null;
    const rowClasses = ['card-row'];
    if (isSelected) {
      rowClasses.push('is-active');
    }
    if (isRuntimeActive) {
      rowClasses.push('is-runtime-active');
    }
    if (isSelected && isRuntimeActive) {
      rowClasses.push('is-active-runtime');
    }

    return `
      <button class="${rowClasses.join(' ')}" data-row="${dataset.all_excel_row}" type="button">
        <div class="result-primary">
          <h3>${escapeHtml(dataset.title || 'Untitled dataset')}</h3>
          ${isRuntimeActive ? `<div class="badge-row">
            <span class="badge runtime-active-badge">Active dataset</span>
          </div>` : ''}
        </div>
        <div class="metric"><span>Author</span><strong>${escapeHtml(dataset.author || '—')}</strong></div>
        <div class="metric"><span>Year</span><strong>${escapeHtml(dataset.year || '—')}</strong></div>
        <div class="metric"><span>Tissue</span><strong>${escapeHtml(listPreview(dataset.tissues))}</strong></div>
        <div class="metric"><span>Disease</span><strong>${escapeHtml(listPreview(dataset.diseases))}</strong></div>
        <div class="metric"><span>Organism</span><strong>${escapeHtml(listPreview(dataset.organisms))}</strong></div>
        <div class="metric"><span>Cells</span><strong>${formatNumber(dataset.cell_counts)}</strong></div>
      </button>
    `;
  }).join('');

  elements.resultsContainer.innerHTML = `
    <div class="results-head">
      <span>Dataset</span>
      <span>Author</span>
      <span>Year</span>
      <span>Tissue</span>
      <span>Disease</span>
      <span>Organism</span>
      <span>Cells</span>
    </div>
    ${desktopRows}
    <div class="card-list">${mobileCards}</div>
  `;
}

function renderPagination(filteredDatasets) {
  const pageCount = Math.max(1, Math.ceil(filteredDatasets.length / state.pageSize));
  if (state.page > pageCount) {
    state.page = pageCount;
  }

  if (pageCount <= 1) {
    elements.pagination.innerHTML = '';
    return;
  }

  const buttons = [];
  const start = Math.max(1, state.page - 2);
  const end = Math.min(pageCount, start + 4);
  buttons.push(`<button class="page-button" type="button" data-page="${Math.max(1, state.page - 1)}">Previous</button>`);
  for (let page = start; page <= end; page += 1) {
    buttons.push(`<button class="page-button ${page === state.page ? 'is-current' : ''}" type="button" data-page="${page}">${page}</button>`);
  }
  buttons.push(`<button class="page-button" type="button" data-page="${Math.min(pageCount, state.page + 1)}">Next</button>`);
  elements.pagination.innerHTML = buttons.join('');
}

async function loadDetail(allExcelRow) {
  state.selectedRow = allExcelRow;
  const activeVariantRow = getActiveVariantRow();
  state.selectedVariantRow = (isRuntimeActiveDatasetRow(allExcelRow) && activeVariantRow !== null)
    ? activeVariantRow
    : null;
  state.detailLoading = true;
  state.detailError = '';
  state.detailActionLoading = false;
  state.detailActionMessage = '';
  state.detailActionTone = 'success';
  document.body.classList.add('detail-open');
  renderDetail();

  if (state.detailCache.has(allExcelRow)) {
    state.detailLoading = false;
    renderDetail();
    return;
  }

  try {
    const response = await fetch(`/api/genopixel-catalog/datasets/${allExcelRow}`);
    if (!response.ok) {
      throw new Error(`Failed to load dataset ${allExcelRow}.`);
    }
    const payload = await response.json();
    state.detailCache.set(allExcelRow, payload.dataset);
  } catch (error) {
    state.detailError = error.message;
  } finally {
    state.detailLoading = false;
    renderDetail();
  }
}

function detailFieldRow(label, value) {
  return `
    <dt>${escapeHtml(label)}</dt>
    <dd>${escapeHtml(value)}</dd>
  `;
}

function renderDetail() {
  if (!state.selectedRow) {
    elements.detailTitle.textContent = 'Select a dataset';
    elements.detailBody.innerHTML = '<p class="detail-placeholder">Choose a row to inspect the workbook metadata and any variant files.</p>';
    return;
  }

  const detail = state.detailCache.get(state.selectedRow);

  if (state.detailLoading && !detail) {
    elements.detailTitle.textContent = `Dataset ${state.selectedRow}`;
    elements.detailBody.innerHTML = '<div class="loading-state">Loading dataset detail…</div>';
    return;
  }

  if (state.detailError) {
    elements.detailTitle.textContent = `Dataset ${state.selectedRow}`;
    elements.detailBody.innerHTML = `<section class="status-banner">${escapeHtml(state.detailError)}</section>`;
    return;
  }

  if (!detail) {
    return;
  }

  const activeParentRow = getActiveDatasetRow();
  const activeVariantRow = getActiveVariantRow();
  const isActiveParent = activeParentRow === state.selectedRow;
  const variantList = detail.variants || [];
  if (isActiveParent && state.selectedVariantRow === null && activeVariantRow !== null) {
    state.selectedVariantRow = activeVariantRow;
  }
  const activeVariantMissing = isActiveParent
    && activeVariantRow !== null
    && variantList.length > 0
    && !variantList.some((variant) => variant.multiple_excel_row === activeVariantRow);

  elements.detailTitle.textContent = detail.title || `Dataset ${state.selectedRow}`;
  const detailActionStatus = state.detailActionMessage
    ? `<section class="status-banner detail-action-status" data-tone="${escapeHtml(state.detailActionTone)}">${escapeHtml(state.detailActionMessage)}</section>`
    : '';
  const runtimeActiveDetailBadge = isActiveParent
    ? '<div class="detail-runtime-badge"><span class="inline-badge runtime-active-badge">Runtime active dataset</span></div>'
    : '';
  const variantRows = variantList.map((variant) => {
    const isRuntimeActiveVariant = isActiveParent && activeVariantRow === variant.multiple_excel_row;
    const isSelectedVariant = state.selectedVariantRow === variant.multiple_excel_row;
    return `
    <div class="variant-row ${isRuntimeActiveVariant ? 'is-runtime-active' : ''} ${isSelectedVariant ? 'is-selected' : ''}">
      <div class="variant-select-cell">
        <input
          type="radio"
          name="subdatasetSelection"
          aria-label="Select sub dataset row ${variant.multiple_excel_row}"
          data-variant-row="${variant.multiple_excel_row}"
          ${state.selectedVariantRow === variant.multiple_excel_row ? 'checked' : ''}
        />
      </div>
      <div>
        ${escapeHtml(variant.description || '—')}
        ${isSelectedVariant && !isRuntimeActiveVariant ? '<div class="variant-active-indicator"><span class="inline-badge selected-file-badge">Selected file</span></div>' : ''}
        ${isRuntimeActiveVariant ? '<div class="variant-active-indicator"><span class="inline-badge runtime-active-badge">Active file</span></div>' : ''}
      </div>
      <div>${formatNumber(variant.cell_counts)}</div>
    </div>
  `;
  }).join('');

  elements.detailBody.innerHTML = `
    <section class="detail-section">
      <dl class="field-list">
        ${detailFieldRow('Author', detail.author || '—')}
        ${detailFieldRow('Year', detail.year || '—')}
        ${detailFieldRow('Journal', detail.journal || '—')}
        ${detailFieldRow('DOI', detail.doi || '—')}
      </dl>
      <div class="detail-link-row">
        ${detail.doi ? `<a class="link-button detail-link-button" href="${escapeHtml(detail.doi)}" target="_blank" rel="noreferrer">Open DOI</a>` : ''}
      </div>
      <dl class="field-list">
        ${detailFieldRow('Tissue', (detail.tissues || []).join(', ') || '—')}
        ${detailFieldRow('Tissue type', (detail.tissue_types || []).join(', ') || '—')}
        ${detailFieldRow('Disease', (detail.diseases || []).join(', ') || '—')}
        ${detailFieldRow('Organism', (detail.organisms || []).join(', ') || '—')}
        ${detailFieldRow('Cell counts', formatNumber(detail.cell_counts))}
      </dl>
      ${runtimeActiveDetailBadge}
    </section>

    ${variantList.length ? `
      <section class="variant-table-wrap">
        <div class="subdataset-header">
          <div class="subdataset-copy">
            <h3>Choose a sub dataset</h3>
            <p class="table-caption">This project includes too many cells. For your convenience, it is splitted into smaller sub datasets.</p>
          </div>
          <button class="subdataset-action" data-detail-action="analyze-subdataset" type="button" ${state.detailActionLoading ? 'disabled' : ''}>
            ${state.detailActionLoading ? 'Loading dataset...' : 'Analyze'}
          </button>
        </div>
        ${detailActionStatus}
        ${activeVariantMissing ? `<p class="active-dataset-warning detail-active-variant-warning">The active sub-dataset was not found in the current detail payload.</p>` : ''}
        <div class="variant-header">
          <span>Select</span>
          <span>Description</span>
          <span>Cells</span>
        </div>
        ${variantRows}
      </section>
    ` : `
      <section class="detail-section detail-action-panel">
        <button class="subdataset-action" data-detail-action="analyze-dataset" type="button" ${state.detailActionLoading ? 'disabled' : ''}>
          ${state.detailActionLoading ? 'Loading dataset...' : 'Analyze this dataset'}
        </button>
        ${detailActionStatus}
      </section>
    `}
  `;
}

function render() {
  renderFacets();
  const filteredDatasets = getFilteredDatasets();
  const pageCount = Math.max(1, Math.ceil(filteredDatasets.length / state.pageSize));
  state.page = Math.min(Math.max(1, state.page), pageCount);
  renderActiveDatasetBanner(filteredDatasets);
  renderActiveFilters(filteredDatasets);
  renderResults(filteredDatasets);
  renderPagination(filteredDatasets);
  renderDetail();
}

function clearAll() {
  clearSearchFacetAndYearFilters();
  state.sort = 'smart';
  state.page = 1;
  elements.sortSelect.value = 'smart';
  state.ranges.cellMin = '';
  state.ranges.cellMax = '';
  render();
}

function focusDatasetRowInResults(allExcelRow, { forceShow = false, smooth = true } = {}) {
  let filteredDatasets = getFilteredDatasets();
  let activeIndex = filteredDatasets.findIndex((dataset) => dataset.all_excel_row === allExcelRow);

  if (activeIndex < 0 && forceShow) {
    clearSearchFacetAndYearFilters();
    state.page = 1;
    filteredDatasets = getFilteredDatasets();
    activeIndex = filteredDatasets.findIndex((dataset) => dataset.all_excel_row === allExcelRow);
  }

  if (activeIndex < 0) {
    return false;
  }

  state.page = Math.floor(activeIndex / state.pageSize) + 1;
  setStatus('');
  render();

  if (elements.resultsPanel instanceof HTMLElement) {
    elements.resultsPanel.scrollTo({ top: 0, behavior: smooth ? 'smooth' : 'auto' });
  }
  const rowElement = elements.resultsContainer.querySelector(`[data-row="${allExcelRow}"]`);
  if (rowElement instanceof HTMLElement) {
    rowElement.scrollIntoView({ behavior: smooth ? 'smooth' : 'auto', block: 'center' });
  }
  return true;
}

async function showActiveDataset() {
  const activeRow = getActiveDatasetRow();
  if (activeRow === null) {
    return;
  }

  if (!getDatasetByRow(activeRow)) {
    setStatus('Active dataset not found in current catalog snapshot.', 'error');
    render();
    return;
  }
  const focused = focusDatasetRowInResults(activeRow, { forceShow: true });
  if (!focused) {
    setStatus('Active dataset not found in current catalog snapshot.', 'error');
    render();
    return;
  }
  if (state.selectedRow !== activeRow) {
    await loadDetail(activeRow);
  } else {
    renderDetail();
  }
  closeFilters();
}

function closeDetail() {
  document.body.classList.remove('detail-open');
}

function closeFilters() {
  document.body.classList.remove('filters-open');
}

function handleChipClick(kind, label) {
  if (kind === 'search') {
    state.search = '';
    elements.searchInput.value = '';
  } else if (kind === 'range') {
    if (label.startsWith('Year:')) {
      state.ranges.yearMin = '';
      state.ranges.yearMax = '';
      elements.yearMin.value = '';
      elements.yearMax.value = '';
    }
  } else {
    const [, value] = label.split(': ');
    state.filters[kind].delete(value);
  }
  state.page = 1;
  render();
}

function attachEvents() {
  elements.searchInput.addEventListener('input', (event) => {
    state.search = event.target.value;
    state.page = 1;
    render();
  });

  elements.sortSelect.addEventListener('change', (event) => {
    state.sort = event.target.value;
    state.page = 1;
    render();
  });

  elements.clearAllButton.addEventListener('click', clearAll);
  elements.closeDetailButton.addEventListener('click', closeDetail);
  elements.detailBackdrop.addEventListener('click', () => {
    closeDetail();
    closeFilters();
  });
  elements.openFiltersButton.addEventListener('click', () => document.body.classList.add('filters-open'));
  elements.closeFiltersButton.addEventListener('click', closeFilters);

  [
    ['yearMin', elements.yearMin],
    ['yearMax', elements.yearMax],
  ].forEach(([key, input]) => {
    input.addEventListener('change', (event) => {
      state.ranges[key] = event.target.value.trim();
      state.page = 1;
      render();
    });
  });

  const handleFacetChange = (event) => {
    const target = event.target;
    if (!(target instanceof HTMLInputElement) || target.type !== 'checkbox') {
      return;
    }
    const facetKey = target.dataset.facetKey;
    const facetValue = target.dataset.facetValue;
    if (!facetKey || !facetValue) {
      return;
    }
    if (target.checked) {
      state.filters[facetKey].add(facetValue);
    } else {
      state.filters[facetKey].delete(facetValue);
    }
    state.page = 1;
    render();
  };

  const handleFacetToggle = (event) => {
    const target = event.target;
    if (!(target instanceof HTMLElement)) {
      return;
    }
    const facetKey = target.dataset.toggleFacet;
    if (!facetKey) {
      return;
    }
    if (state.expandedFacets.has(facetKey)) {
      state.expandedFacets.delete(facetKey);
    } else {
      state.expandedFacets.add(facetKey);
    }
    render();
  };

  const handleFacetSearch = (event) => {
    const target = event.target;
    if (!(target instanceof HTMLInputElement)) {
      return;
    }
    const facetKey = target.dataset.facetSearch;
    if (!facetKey) {
      return;
    }
    const selectionStart = target.selectionStart ?? target.value.length;
    const selectionEnd = target.selectionEnd ?? selectionStart;
    state.facetSearch[facetKey] = target.value;
    renderFacets();
    const nextInput = document.querySelector(`[data-facet-search="${facetKey}"]`);
    if (nextInput instanceof HTMLInputElement) {
      nextInput.focus({ preventScroll: true });
      nextInput.setSelectionRange(selectionStart, selectionEnd);
    }
  };

  [elements.facetSectionsTop, elements.facetSectionsBottom].forEach((container) => {
    container.addEventListener('change', handleFacetChange);
    container.addEventListener('click', handleFacetToggle);
    container.addEventListener('input', handleFacetSearch);
  });

  elements.activeChips.addEventListener('click', (event) => {
    const target = event.target;
    if (!(target instanceof HTMLElement)) {
      return;
    }
    const chipButton = target.closest('.chip');
    if (!(chipButton instanceof HTMLButtonElement)) {
      return;
    }
    handleChipClick(chipButton.dataset.chipKind, chipButton.dataset.chipLabel);
  });

  elements.activeDatasetBanner.addEventListener('click', (event) => {
    const target = event.target;
    if (!(target instanceof HTMLElement)) {
      return;
    }
    const actionButton = target.closest('[data-active-dataset-action]');
    if (!(actionButton instanceof HTMLButtonElement)) {
      return;
    }
    if (actionButton.dataset.activeDatasetAction === 'show') {
      showActiveDataset();
    }
  });

  elements.resultsContainer.addEventListener('click', (event) => {
    const target = event.target;
    if (!(target instanceof HTMLElement)) {
      return;
    }
    const rowButton = target.closest('[data-row]');
    if (!(rowButton instanceof HTMLButtonElement)) {
      return;
    }
    const allExcelRow = Number(rowButton.dataset.row);
    if (!Number.isFinite(allExcelRow)) {
      return;
    }
    loadDetail(allExcelRow);
    closeFilters();
  });

  elements.pagination.addEventListener('click', (event) => {
    const target = event.target;
    if (!(target instanceof HTMLElement)) {
      return;
    }
    const pageButton = target.closest('[data-page]');
    if (!(pageButton instanceof HTMLButtonElement)) {
      return;
    }
    state.page = Number(pageButton.dataset.page);
    render();
    window.scrollTo({ top: 0, behavior: 'smooth' });
  });

  elements.detailBody.addEventListener('change', (event) => {
    const target = event.target;
    if (!(target instanceof HTMLInputElement) || target.type !== 'radio') {
      return;
    }
    const variantRow = Number(target.dataset.variantRow);
    if (!Number.isFinite(variantRow)) {
      return;
    }
    state.selectedVariantRow = variantRow;
    state.detailActionMessage = '';
    state.detailActionTone = 'success';
    renderDetail();
  });

  elements.detailBody.addEventListener('click', async (event) => {
    const target = event.target;
    if (!(target instanceof HTMLElement)) {
      return;
    }
    const actionButton = target.closest('[data-detail-action]');
    if (!(actionButton instanceof HTMLButtonElement) || !state.selectedRow || state.detailActionLoading) {
      return;
    }

    try {
      const detail = state.detailCache.get(state.selectedRow);
      let h5adPath = detail?.resolved_h5ad_path;
      let multipleExcelRow = null;

      if (actionButton.dataset.detailAction === 'analyze-subdataset') {
        if (!state.selectedVariantRow) {
          throw new Error('Please make a choice first.');
        }
        const selectedVariant = (detail?.variants || []).find((variant) => variant.multiple_excel_row === state.selectedVariantRow);
        if (!selectedVariant) {
          throw new Error('Please make a choice first.');
        }
        h5adPath = selectedVariant.resolved_h5ad_path;
        multipleExcelRow = selectedVariant.multiple_excel_row;
      }

      if (!h5adPath) {
        throw new Error('Resolved h5ad path is not available for this dataset.');
      }

      state.detailActionLoading = true;
      state.detailActionMessage = '';
      state.detailActionTone = 'success';
      renderDetail();

      const response = await fetch(`/api/genopixel-catalog/datasets/${state.selectedRow}/analyze`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ h5ad_path: h5adPath, multiple_excel_row: multipleExcelRow }),
      });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(payload.detail || `Analyze request failed with status ${response.status}.`);
      }
      state.detailActionMessage = payload.message || 'data is loaded, happy analysis';
      state.detailActionTone = 'success';

      if (payload && Number.isFinite(Number(payload.all_excel_row))) {
        state.activeDataset = {
          loaded: true,
          all_excel_row: Number(payload.all_excel_row),
          multiple_excel_row: payload.multiple_excel_row ?? null,
          title: payload.active_dataset?.title || detail?.title || '',
          h5ad_path: payload.h5ad_path || h5adPath,
          loaded_at: payload.active_dataset?.loaded_at || new Date().toISOString(),
          backed: payload.backed ?? payload.active_dataset?.backed ?? null,
          total_cells: state.activeDataset?.total_cells ?? null,
        };
        state.activeDatasetLastSyncAt = new Date().toISOString();
        state.activeDatasetSyncError = '';
      }
      await refreshActiveDataset({ silent: true, render: false });
      const activeRow = getActiveDatasetRow();
      if (activeRow !== null) {
        const focused = focusDatasetRowInResults(activeRow, { forceShow: true, smooth: false });
        if (!focused) {
          setStatus('Active dataset not found in current catalog snapshot.', 'error');
        }
      }
    } catch (error) {
      state.detailActionMessage = error.message || 'Failed to load the dataset.';
      state.detailActionTone = 'error';
    } finally {
      state.detailActionLoading = false;
      render();
    }
  });

  window.addEventListener('focus', () => {
    refreshActiveDataset({ silent: true });
  });
  document.addEventListener('visibilitychange', () => {
    if (document.visibilityState === 'visible') {
      refreshActiveDataset({ silent: true });
    }
  });
}

async function loadCatalog() {
  setStatus('Loading catalog…');
  try {
    const response = await fetch('/api/genopixel-catalog/catalog');
    if (!response.ok) {
      throw new Error(`Catalog request failed with status ${response.status}.`);
    }
    const payload = await response.json();
    state.catalog = payload;
    state.datasets = payload.datasets.map(buildDatasetIndex);

    const url = new URL(window.location.href);
    const query = url.searchParams.get('q');
    if (query) {
      state.search = query;
      elements.searchInput.value = query;
    }

    if (state.datasets.length > 0) {
      const years = state.datasets.map((dataset) => dataset.year).filter(Boolean);
      elements.yearMin.placeholder = String(Math.min(...years));
      elements.yearMax.placeholder = String(Math.max(...years));
    }

    await refreshActiveDataset({ silent: true, render: false });
    setStatus('');
    render();
  } catch (error) {
    setStatus(error.message || 'Failed to load the catalog.', 'error');
    elements.resultsContainer.innerHTML = `
      <section class="empty-state">
        <h3>The dataset browser could not load the catalog.</h3>
        <p>${escapeHtml(error.message || 'Unknown error.')}</p>
      </section>
    `;
  }
}

attachEvents();
startActiveDatasetSync();
loadCatalog();
