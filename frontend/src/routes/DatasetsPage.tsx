"use client"

import { useEffect, useState, useCallback, useMemo, useRef } from "react"
import { useAuth } from "react-oidc-context"
import { AppSidebar } from "@/components/layout/AppSidebar"
import "./DatasetsPage.css"
import {
  listDatasets,
  getDataset,
  analyzeDataset,
  getActiveDataset,
  type Dataset,
  type ActiveDatasetResponse,
  type DatasetVariant,
} from "@/services/catalogService"

// ---- Constants ----

const FACETS = [
  { key: "tissue", label: "Tissue", field: "tissues", type: "array" },
  { key: "disease", label: "Disease", field: "diseases", type: "array" },
  { key: "organism", label: "Organism", field: "organisms", type: "array" },
  { key: "tissue_type", label: "Tissue type", field: "tissue_types", type: "array" },
  { key: "journal", label: "Journal", field: "journal", type: "scalar" },
] as const

type FacetKey = "tissue" | "disease" | "organism" | "tissue_type" | "journal"
type Filters = Record<FacetKey, Set<string>>

const PAGE_SIZE = 24
const SYNC_INTERVAL_MS = 30_000

// ---- Helpers ----

function formatNumber(value?: number | null): string {
  if (value === null || value === undefined) return "—"
  return new Intl.NumberFormat("en-US").format(Number(value))
}

function listPreview(values?: string[], limit = 2): string {
  if (!Array.isArray(values) || values.length === 0) return "—"
  if (values.length <= limit) return values.join(", ")
  return `${values.slice(0, limit).join(", ")} +${values.length - limit}`
}

function tokenize(query: string): string[] {
  return query.trim().toLowerCase().split(/\s+/).filter(Boolean)
}

function emptyFilters(): Filters {
  return { tissue: new Set(), disease: new Set(), organism: new Set(), tissue_type: new Set(), journal: new Set() }
}

// ---- Search-indexed dataset ----

interface IndexedDataset extends Dataset {
  _titleSearch: string
  _authorJournalSearch: string
  _searchBlob: string
  _facetValues: Record<FacetKey, string[]>
}

function buildIndex(dataset: Dataset): IndexedDataset {
  const titleSearch = (dataset.title || "").toLowerCase()
  const authorJournalSearch = [dataset.author, dataset.journal].filter(Boolean).join(" ").toLowerCase()
  const chunks = [
    dataset.title, dataset.author, dataset.journal,
    dataset.doi, dataset.cellxgene_doi, dataset.project,
    ...(dataset.tissues || []),
    ...(dataset.diseases || []),
    ...(dataset.organisms || []),
    ...(dataset.tissue_types || []),
  ].filter(Boolean)
  return {
    ...dataset,
    _titleSearch: titleSearch,
    _authorJournalSearch: authorJournalSearch,
    _searchBlob: chunks.join(" ").toLowerCase(),
    _facetValues: {
      tissue: (dataset.tissues || []).map(v => v.toLowerCase()),
      disease: (dataset.diseases || []).map(v => v.toLowerCase()),
      organism: (dataset.organisms || []).map(v => v.toLowerCase()),
      tissue_type: (dataset.tissue_types || []).map(v => v.toLowerCase()),
      journal: dataset.journal ? [dataset.journal.toLowerCase()] : [],
    },
  }
}

function searchTier(dataset: IndexedDataset, tokens: string[]): number {
  if (!tokens.length) return 0
  if (!tokens.every(t => dataset._searchBlob.includes(t))) return -1
  if (tokens.every(t => dataset._titleSearch.includes(t))) return 3
  if (tokens.every(t => dataset._authorJournalSearch.includes(t))) return 2
  return 1
}

// ---- FacetCard ----

interface FacetEntry { value: string; label: string; count: number }

interface FacetCardProps {
  facetKey: FacetKey
  label: string
  counts: FacetEntry[]
  selected: Set<string>
  searchValue: string
  expanded: boolean
  onToggle: (key: FacetKey, label: string, checked: boolean) => void
  onSearchChange: (key: FacetKey, value: string) => void
  onToggleExpand: (key: FacetKey) => void
}

function FacetCard({ facetKey, label, counts, selected, searchValue, expanded, onToggle, onSearchChange, onToggleExpand }: FacetCardProps) {
  const query = (searchValue || "").trim().toLowerCase()
  const filtered = query ? counts.filter(e => e.label.toLowerCase().includes(query)) : counts
  const visible = query ? filtered : (expanded ? filtered : filtered.slice(0, 10))
  const hasSearch = ["tissue", "disease", "journal"].includes(facetKey)
  return (
    <section className="facet-card">
      <div className="range-header">
        <div>
          <h3>{label}</h3>
          <p>{counts.length.toLocaleString()} values</p>
        </div>
      </div>
      {hasSearch && (
        <div className="facet-search-shell">
          <input
            className="facet-search-input"
            type="search"
            placeholder={`Type to search ${label.toLowerCase()}`}
            value={searchValue || ""}
            onChange={e => onSearchChange(facetKey, e.target.value)}
          />
        </div>
      )}
      <div className="facet-list">
        {visible.length > 0 ? visible.map(entry => (
          <div className="facet-option" key={entry.value}>
            <label>
              <input
                type="checkbox"
                checked={selected.has(entry.label)}
                onChange={e => onToggle(facetKey, entry.label, e.target.checked)}
              />
              <span>{entry.label}</span>
            </label>
            <span className="facet-count">{entry.count.toLocaleString()}</span>
          </div>
        )) : (
          <p className="facet-empty">No matching {label.toLowerCase()} values.</p>
        )}
      </div>
      {counts.length > 10 && !query && (
        <button className="toggle-button" onClick={() => onToggleExpand(facetKey)}>
          {expanded ? "Show fewer" : "Show all"}
        </button>
      )}
    </section>
  )
}

// ---- ActiveDatasetBanner ----

interface BannerProps {
  active: ActiveDatasetResponse
  syncError: string
  filteredDatasets: IndexedDataset[]
  allDatasets: IndexedDataset[]
  page: number
  onShow: () => void
}

function ActiveDatasetBanner({ active, syncError, filteredDatasets, allDatasets, page, onShow }: BannerProps) {
  const activeRow = active.loaded ? Number(active.all_excel_row) : null
  if (!activeRow) return null

  const inCatalog = allDatasets.some(d => d.all_excel_row === activeRow)
  const filteredIdx = filteredDatasets.findIndex(d => d.all_excel_row === activeRow)
  const inFiltered = filteredIdx >= 0
  const pageForActive = inFiltered ? Math.floor(filteredIdx / PAGE_SIZE) + 1 : null
  const pageStart = (page - 1) * PAGE_SIZE
  const inCurrentPage = inFiltered && filteredIdx >= pageStart && filteredIdx < pageStart + PAGE_SIZE

  let tone = "info"
  let helperText = ""
  if (!inCatalog) {
    tone = "warn"
    helperText = "Active dataset not found in current catalog snapshot."
  } else if (!inFiltered) {
    tone = "warn"
    helperText = "Active dataset is hidden by current filters."
  } else if (!inCurrentPage) {
    helperText = `Active dataset is on page ${pageForActive}.`
  }

  const title = active.title || `Dataset ${activeRow}`

  return (
    <div className="active-dataset-banner-shell">
      <section className="active-dataset-banner" data-tone={tone}>
        <div className="active-dataset-copy">
          <p className="active-dataset-label">Runtime active dataset</p>
          <p className="active-dataset-title">{title}</p>
          {helperText && <p className="active-dataset-warning">{helperText}</p>}
          {syncError && <p className="active-dataset-warning">{syncError}</p>}
        </div>
        <button className="active-dataset-jump" onClick={onShow}>Show active dataset</button>
      </section>
    </div>
  )
}

// ---- Detail drawer body ----

interface DrawerProps {
  selectedRow: number | null
  currentDetail: Dataset | null
  detailLoading: boolean
  detailError: string
  selectedVariantRow: number | null
  activeParentRow: number | null
  activeVariantRow: number | null
  detailActionLoading: boolean
  detailActionMessage: string
  detailActionTone: "success" | "error"
  onClose: () => void
  onSelectVariant: (row: number) => void
  onAnalyze: (row: number, variantRow?: number) => void
}

function DetailDrawer({
  selectedRow, currentDetail, detailLoading, detailError,
  selectedVariantRow, activeParentRow, activeVariantRow,
  detailActionLoading, detailActionMessage, detailActionTone,
  onClose, onSelectVariant, onAnalyze,
}: DrawerProps) {
  const isActiveParent = activeParentRow === selectedRow
  const variantList: DatasetVariant[] = currentDetail?.variants || []
  const activeVariantMissing = isActiveParent
    && activeVariantRow !== null
    && variantList.length > 0
    && !variantList.some(v => v.multiple_excel_row === activeVariantRow)

  return (
    <>
      <div className="detail-header">
        <div>
          <p className="panel-kicker">Dataset detail</p>
          <h2>{currentDetail?.title || (selectedRow ? `Dataset ${selectedRow}` : "Select a dataset")}</h2>
        </div>
        <button className="ghost-button" onClick={onClose}>Close</button>
      </div>
      <div className="detail-body">
        {!selectedRow && (
          <p className="detail-placeholder">Choose a row to inspect the workbook metadata and any variant files.</p>
        )}
        {selectedRow && detailLoading && !currentDetail && (
          <div className="loading-state">Loading dataset detail…</div>
        )}
        {selectedRow && detailError && (
          <section className="status-banner">{detailError}</section>
        )}
        {currentDetail && (
          <>
            <section className="detail-section">
              <dl className="field-list">
                <dt>Author</dt><dd>{currentDetail.author || "—"}</dd>
                <dt>Year</dt><dd>{currentDetail.year || "—"}</dd>
                <dt>Journal</dt><dd>{currentDetail.journal || "—"}</dd>
                <dt>DOI</dt><dd>{currentDetail.doi || "—"}</dd>
              </dl>
              {currentDetail.doi && (
                <div className="detail-link-row">
                  <a className="link-button detail-link-button" href={`https://doi.org/${currentDetail.doi}`} target="_blank" rel="noreferrer">Open DOI</a>
                </div>
              )}
              <dl className="field-list">
                <dt>Tissue</dt><dd>{(currentDetail.tissues || []).join(", ") || "—"}</dd>
                <dt>Tissue type</dt><dd>{(currentDetail.tissue_types || []).join(", ") || "—"}</dd>
                <dt>Disease</dt><dd>{(currentDetail.diseases || []).join(", ") || "—"}</dd>
                <dt>Organism</dt><dd>{(currentDetail.organisms || []).join(", ") || "—"}</dd>
                <dt>Cell counts</dt><dd>{formatNumber(currentDetail.cell_counts)}</dd>
              </dl>
              {isActiveParent && (
                <div className="detail-runtime-badge">
                  <span className="inline-badge runtime-active-badge">Runtime active dataset</span>
                </div>
              )}
            </section>

            {variantList.length > 0 ? (
              <section className="variant-table-wrap">
                <div className="subdataset-header">
                  <div className="subdataset-copy">
                    <h3>Choose a sub dataset</h3>
                    <p className="table-caption">This project includes too many cells. For your convenience, it is splitted into smaller sub datasets.</p>
                  </div>
                  <button
                    className="subdataset-action"
                    onClick={() => selectedRow && onAnalyze(selectedRow, selectedVariantRow ?? undefined)}
                    disabled={detailActionLoading || !selectedVariantRow}
                  >
                    {detailActionLoading ? "Loading dataset..." : "Analyze"}
                  </button>
                </div>
                {detailActionMessage && (
                  <section className="status-banner detail-action-status" data-tone={detailActionTone}>
                    {detailActionMessage}
                  </section>
                )}
                {activeVariantMissing && (
                  <p className="active-dataset-warning detail-active-variant-warning">
                    The active sub-dataset was not found in the current detail payload.
                  </p>
                )}
                <div className="variant-header">
                  <span>Select</span>
                  <span>Description</span>
                  <span>Cells</span>
                </div>
                {variantList.map(variant => {
                  const isRuntimeActive = isActiveParent && activeVariantRow === variant.multiple_excel_row
                  const isSelected = selectedVariantRow === variant.multiple_excel_row
                  return (
                    <div
                      key={variant.multiple_excel_row}
                      className={`variant-row ${isRuntimeActive ? "is-runtime-active" : ""} ${isSelected ? "is-selected" : ""}`}
                    >
                      <div className="variant-select-cell">
                        <input
                          type="radio"
                          name="subdatasetSelection"
                          checked={isSelected}
                          onChange={() => onSelectVariant(variant.multiple_excel_row)}
                        />
                      </div>
                      <div>
                        {variant.description || "—"}
                        {isSelected && !isRuntimeActive && (
                          <div className="variant-active-indicator">
                            <span className="inline-badge selected-file-badge">Selected file</span>
                          </div>
                        )}
                        {isRuntimeActive && (
                          <div className="variant-active-indicator">
                            <span className="inline-badge runtime-active-badge">Active file</span>
                          </div>
                        )}
                      </div>
                      <div>{formatNumber(variant.cell_counts ?? variant.cell_count)}</div>
                    </div>
                  )
                })}
              </section>
            ) : (
              <section className="detail-section detail-action-panel">
                <button
                  className="subdataset-action"
                  onClick={() => selectedRow && onAnalyze(selectedRow)}
                  disabled={detailActionLoading}
                >
                  {detailActionLoading ? "Loading dataset..." : "Analyze this dataset"}
                </button>
                {detailActionMessage && (
                  <section className="status-banner detail-action-status" data-tone={detailActionTone}>
                    {detailActionMessage}
                  </section>
                )}
              </section>
            )}
          </>
        )}
      </div>
    </>
  )
}

// ---- Main page ----

export default function DatasetsPage() {
  const auth = useAuth()
  const idToken = auth.user?.id_token ?? ""

  // Catalog
  const [allDatasets, setAllDatasets] = useState<IndexedDataset[]>([])
  const [catalogLoading, setCatalogLoading] = useState(true)
  const [catalogError, setCatalogError] = useState<string | null>(null)
  const [yearMinPlaceholder, setYearMinPlaceholder] = useState("")
  const [yearMaxPlaceholder, setYearMaxPlaceholder] = useState("")

  // Filters
  const [filters, setFilters] = useState<Filters>(emptyFilters)
  const [facetSearch, setFacetSearch] = useState<Record<string, string>>({ tissue: "", disease: "", journal: "" })
  const [yearMin, setYearMin] = useState("")
  const [yearMax, setYearMax] = useState("")
  const [expandedFacets, setExpandedFacets] = useState<Set<string>>(new Set())
  const [search, setSearch] = useState("")
  const [sort, setSort] = useState("smart")
  const [page, setPage] = useState(1)

  // Detail
  const detailCacheRef = useRef<Map<number, Dataset>>(new Map())
  const [currentDetail, setCurrentDetail] = useState<Dataset | null>(null)
  const [selectedRow, setSelectedRow] = useState<number | null>(null)
  const [selectedVariantRow, setSelectedVariantRow] = useState<number | null>(null)
  const [detailLoading, setDetailLoading] = useState(false)
  const [detailError, setDetailError] = useState("")
  const [detailActionLoading, setDetailActionLoading] = useState(false)
  const [detailActionMessage, setDetailActionMessage] = useState("")
  const [detailActionTone, setDetailActionTone] = useState<"success" | "error">("success")
  const [detailOpen, setDetailOpen] = useState(false)
  const [filtersOpen, setFiltersOpen] = useState(false)

  // Active dataset
  const [activeDataset, setActiveDataset] = useState<ActiveDatasetResponse | null>(null)
  const [syncError, setSyncError] = useState("")
  const syncInFlight = useRef(false)

  // Status
  const [statusMessage, setStatusMessage] = useState<{ text: string; tone: string } | null>(null)

  // Load catalog
  useEffect(() => {
    if (!idToken) return
    setCatalogLoading(true)
    listDatasets(idToken, { page_size: 9999 })
      .then(resp => {
        const indexed = resp.datasets.map(buildIndex)
        setAllDatasets(indexed)
        const years = indexed.map(d => d.year).filter((y): y is number => !!y)
        if (years.length) {
          setYearMinPlaceholder(String(Math.min(...years)))
          setYearMaxPlaceholder(String(Math.max(...years)))
        }
        setCatalogLoading(false)
      })
      .catch(err => {
        setCatalogError(err instanceof Error ? err.message : "Failed to load catalog")
        setCatalogLoading(false)
      })
  }, [idToken])

  // Active dataset polling
  const refreshActive = useCallback(async (silent = true) => {
    if (syncInFlight.current || !idToken) return
    syncInFlight.current = true
    try {
      const payload = await getActiveDataset(idToken)
      setActiveDataset(payload.loaded ? payload : null)
      setSyncError("")
    } catch (err) {
      if (!silent) setSyncError(err instanceof Error ? err.message : "Sync failed")
    } finally {
      syncInFlight.current = false
    }
  }, [idToken])

  useEffect(() => {
    if (!idToken) return
    refreshActive(true)
    const timer = setInterval(() => {
      if (document.visibilityState === "visible") refreshActive(true)
    }, SYNC_INTERVAL_MS)
    const onFocus = () => refreshActive(true)
    const onVisibility = () => { if (document.visibilityState === "visible") refreshActive(true) }
    window.addEventListener("focus", onFocus)
    document.addEventListener("visibilitychange", onVisibility)
    return () => {
      clearInterval(timer)
      window.removeEventListener("focus", onFocus)
      document.removeEventListener("visibilitychange", onVisibility)
    }
  }, [idToken, refreshActive])

  // ---- Derived: filtered + sorted datasets ----

  const tokens = useMemo(() => tokenize(search), [search])

  const passesRanges = useCallback((d: IndexedDataset) => {
    const y = Number(d.year || 0)
    if (yearMin !== "" && y < Number(yearMin)) return false
    if (yearMax !== "" && y > Number(yearMax)) return false
    return true
  }, [yearMin, yearMax])

  const passesFacet = useCallback((d: IndexedDataset, key: FacetKey, ignore = "") => {
    if (key === ignore) return true
    const sel = filters[key]
    if (!sel || sel.size === 0) return true
    const selVals = [...sel].map(v => v.toLowerCase())
    return selVals.some(v => (d._facetValues[key] || []).includes(v))
  }, [filters])

  const passesAll = useCallback((d: IndexedDataset, toks: string[], ignore = "") => {
    if (!passesRanges(d)) return false
    if (toks.length > 0 && searchTier(d, toks) < 0) return false
    return FACETS.every(f => passesFacet(d, f.key, ignore))
  }, [passesRanges, passesFacet])

  const compareDs = useCallback((a: IndexedDataset, b: IndexedDataset, toks: string[]) => {
    if (sort === "title_asc") return (a.title || "").localeCompare(b.title || "")
    if (sort === "cells_desc") return (b.cell_counts || 0) - (a.cell_counts || 0) || (b.year || 0) - (a.year || 0)
    if (sort === "year_desc") return (b.year || 0) - (a.year || 0) || (b.cell_counts || 0) - (a.cell_counts || 0) || (a.title || "").localeCompare(b.title || "")
    if (toks.length > 0) {
      const at = searchTier(a, toks), bt = searchTier(b, toks)
      return bt - at || (b.year || 0) - (a.year || 0) || (b.cell_counts || 0) - (a.cell_counts || 0) || (a.title || "").localeCompare(b.title || "")
    }
    return (b.year || 0) - (a.year || 0) || (b.cell_counts || 0) - (a.cell_counts || 0) || (a.title || "").localeCompare(b.title || "")
  }, [sort])

  const filteredDatasets = useMemo(() =>
    allDatasets.filter(d => passesAll(d, tokens)).sort((a, b) => compareDs(a, b, tokens)),
    [allDatasets, tokens, passesAll, compareDs]
  )

  // Facet counts (per-facet, excluding that facet's own filter)
  const getFacetCounts = useCallback((key: FacetKey): FacetEntry[] => {
    const counts = new Map<string, number>()
    const labels = new Map<string, string>()
    allDatasets.forEach(d => {
      if (!passesAll(d, tokens, key)) return
      ;[...new Set(d._facetValues[key] || [])].forEach(v => {
        counts.set(v, (counts.get(v) || 0) + 1)
      })
    })
    const facetDef = FACETS.find(f => f.key === key)!
    allDatasets.forEach(d => {
      const fieldValue = (d as unknown as Record<string, unknown>)[facetDef.field]
      const raws: string[] = facetDef.type === "array"
        ? Array.isArray(fieldValue)
          ? fieldValue.filter((value): value is string => typeof value === "string")
          : []
        : typeof fieldValue === "string" && fieldValue
          ? [fieldValue]
          : []
      raws.forEach(v => labels.set(v.toLowerCase(), v))
    })
    return [...counts.entries()]
      .map(([v, c]) => ({ value: v, label: labels.get(v) || v, count: c }))
      .sort((a, b) => b.count - a.count || a.label.localeCompare(b.label))
  }, [allDatasets, tokens, passesAll])

  // Precompute facet counts
  const facetCounts = useMemo(() => ({
    tissue: getFacetCounts("tissue"),
    disease: getFacetCounts("disease"),
    organism: getFacetCounts("organism"),
    tissue_type: getFacetCounts("tissue_type"),
    journal: getFacetCounts("journal"),
  }), [getFacetCounts])

  // Page bounds
  const totalPages = Math.max(1, Math.ceil(filteredDatasets.length / PAGE_SIZE))
  const safePage = Math.min(Math.max(1, page), totalPages)
  const pageStart = (safePage - 1) * PAGE_SIZE
  const visibleDatasets = filteredDatasets.slice(pageStart, pageStart + PAGE_SIZE)

  // Active dataset helpers
  const activeRow = activeDataset?.loaded ? Number(activeDataset.all_excel_row) : null
  const activeVariantRow = activeDataset?.loaded && activeDataset.multiple_excel_row != null
    ? Number(activeDataset.multiple_excel_row) : null

  // ---- Actions ----

  const clearAll = useCallback(() => {
    setSearch("")
    setFilters(emptyFilters())
    setFacetSearch({ tissue: "", disease: "", journal: "" })
    setYearMin("")
    setYearMax("")
    setSort("smart")
    setPage(1)
  }, [])

  const toggleFilter = useCallback((key: FacetKey, label: string, checked: boolean) => {
    setFilters(prev => {
      const next = { ...prev, [key]: new Set(prev[key]) }
      checked ? next[key].add(label) : next[key].delete(label)
      return next
    })
    setPage(1)
  }, [])

  const toggleFacetExpand = useCallback((key: FacetKey) => {
    setExpandedFacets(prev => {
      const next = new Set(prev)
      next.has(key) ? next.delete(key) : next.add(key)
      return next
    })
  }, [])

  const loadDetail = useCallback(async (row: number) => {
    setSelectedRow(row)
    const curActiveVariantRow = activeDataset?.loaded && activeDataset.multiple_excel_row != null
      ? Number(activeDataset.multiple_excel_row) : null
    setSelectedVariantRow(activeRow === row && curActiveVariantRow !== null ? curActiveVariantRow : null)
    setDetailActionLoading(false)
    setDetailActionMessage("")
    setDetailOpen(true)

    const cached = detailCacheRef.current.get(row)
    if (cached) {
      setCurrentDetail(cached)
      return
    }

    setDetailLoading(true)
    setDetailError("")
    setCurrentDetail(null)
    try {
      const detail = await getDataset(idToken, row)
      detailCacheRef.current.set(row, detail)
      setCurrentDetail(detail)
    } catch (err) {
      setDetailError(err instanceof Error ? err.message : "Failed to load dataset detail.")
    } finally {
      setDetailLoading(false)
    }
  }, [idToken, activeDataset, activeRow])

  const closeDetail = useCallback(() => setDetailOpen(false), [])

  const handleAnalyze = useCallback(async (row: number, variantRow?: number) => {
    setDetailActionLoading(true)
    setDetailActionMessage("")
    try {
      const resp = await analyzeDataset(idToken, row, variantRow)
      setDetailActionMessage(resp.message || "Data is loaded, happy analysis")
      setDetailActionTone("success")
      await refreshActive(true)
    } catch (err) {
      setDetailActionMessage(err instanceof Error ? err.message : "Failed to load the dataset.")
      setDetailActionTone("error")
    } finally {
      setDetailActionLoading(false)
    }
  }, [idToken, refreshActive])

  const showActiveDataset = useCallback(async () => {
    if (activeRow === null) return
    const inCatalog = allDatasets.some(d => d.all_excel_row === activeRow)
    if (!inCatalog) {
      setStatusMessage({ text: "Active dataset not found in current catalog snapshot.", tone: "error" })
      return
    }
    // Clear all filters, then find the active dataset's page in the clean sorted list
    const cleanFilters = emptyFilters()
    const cleanTokens: string[] = []
    const cleanFiltered = allDatasets
      .filter(d => {
        // passes all filters when all empty
        const y = Number(d.year || 0)
        if ("" !== "" && y < Number("")) return false
        if ("" !== "" && y > Number("")) return false
        return true
      })
      .sort((a, b) => compareDs(a, b, cleanTokens))
    const idx = cleanFiltered.findIndex(d => d.all_excel_row === activeRow)
    setSearch("")
    setFilters(cleanFilters)
    setFacetSearch({ tissue: "", disease: "", journal: "" })
    setYearMin("")
    setYearMax("")
    if (idx >= 0) setPage(Math.floor(idx / PAGE_SIZE) + 1)
    setStatusMessage(null)
    setFiltersOpen(false)
    await loadDetail(activeRow)
  }, [activeRow, allDatasets, compareDs, loadDetail])

  // Active chips
  const activeChips = useMemo(() => {
    const chips: { kind: string; label: string }[] = []
    if (search.trim()) chips.push({ kind: "search", label: `Search: ${search.trim()}` })
    FACETS.forEach(f => {
      filters[f.key].forEach(v => chips.push({ kind: f.key, label: `${f.label}: ${v}` }))
    })
    if (yearMin !== "" || yearMax !== "") {
      chips.push({ kind: "range", label: `Year: ${yearMin || "min"} to ${yearMax || "max"}` })
    }
    return chips
  }, [search, filters, yearMin, yearMax])

  const removeChip = useCallback((kind: string, label: string) => {
    if (kind === "search") {
      setSearch("")
    } else if (kind === "range") {
      setYearMin("")
      setYearMax("")
    } else {
      const [, value] = label.split(": ")
      setFilters(prev => {
        const next = { ...prev, [kind]: new Set(prev[kind as FacetKey]) }
        next[kind as FacetKey].delete(value)
        return next
      })
    }
    setPage(1)
  }, [])

  // ---- Render ----

  const facetsTop = FACETS.filter(f => f.key !== "journal")
  const facetsBottom = FACETS.filter(f => f.key === "journal")

  return (
    <div className="flex h-screen bg-background">
      <AppSidebar activeTab="datasets" />
      <div className="flex-1 min-w-0 overflow-auto">
        <div className={`gp-datasets ${detailOpen ? "detail-open" : ""} ${filtersOpen ? "filters-open" : ""}`}>
          <div className="page-shell">
            <header className="hero" style={{ textAlign: "center" }}>
              <div className="hero-copy">
                <div className="hero-text">
                  <h1 style={{ fontSize: "clamp(1.5rem, 2.6vw, 2.4rem)" }}>Dataset browser for GenoPixel catalog</h1>
                  <p>Search, facet, and inspect the GenoPixel catalog before analysis.</p>
                </div>
              </div>
            </header>

            <div className="app-grid">
          {/* Filters sidebar */}
          <aside className="filters-panel">
            <div className="panel-header">
              <div>
                <p className="panel-kicker">Filters</p>
                <h2>Refine the catalog</h2>
              </div>
              <button className="ghost-button close-mobile" onClick={() => setFiltersOpen(false)}>Close</button>
            </div>

            <div className="facet-sections">
              {facetsTop.map(f => (
                <FacetCard
                  key={f.key}
                  facetKey={f.key}
                  label={f.label}
                  counts={facetCounts[f.key]}
                  selected={filters[f.key]}
                  searchValue={facetSearch[f.key] || ""}
                  expanded={expandedFacets.has(f.key)}
                  onToggle={toggleFilter}
                  onSearchChange={(key, val) => setFacetSearch(prev => ({ ...prev, [key]: val }))}
                  onToggleExpand={toggleFacetExpand}
                />
              ))}
            </div>

            <section className="range-card">
              <div className="range-header">
                <h3>Year</h3>
                <p>Inclusive range</p>
              </div>
              <div className="range-grid">
                <label>
                  <span>From</span>
                  <input
                    type="number"
                    inputMode="numeric"
                    value={yearMin}
                    placeholder={yearMinPlaceholder}
                    onChange={e => { setYearMin(e.target.value.trim()); setPage(1) }}
                  />
                </label>
                <label>
                  <span>To</span>
                  <input
                    type="number"
                    inputMode="numeric"
                    value={yearMax}
                    placeholder={yearMaxPlaceholder}
                    onChange={e => { setYearMax(e.target.value.trim()); setPage(1) }}
                  />
                </label>
              </div>
            </section>

            <div className="facet-sections">
              {facetsBottom.map(f => (
                <FacetCard
                  key={f.key}
                  facetKey={f.key}
                  label={f.label}
                  counts={facetCounts[f.key]}
                  selected={filters[f.key]}
                  searchValue={facetSearch[f.key] || ""}
                  expanded={expandedFacets.has(f.key)}
                  onToggle={toggleFilter}
                  onSearchChange={(key, val) => setFacetSearch(prev => ({ ...prev, [key]: val }))}
                  onToggleExpand={toggleFacetExpand}
                />
              ))}
            </div>
          </aside>

          {/* Results panel */}
          <main className="results-panel">
            <section className="toolbar">
              <div className="search-shell">
                <label className="search-label" htmlFor="gp-search-input">
                  Search author, year, journal, or title
                </label>
                <div className="search-row">
                  <input
                    id="gp-search-input"
                    type="search"
                    placeholder="Start with a author name, year of publication, journal name or title"
                    value={search}
                    onChange={e => { setSearch(e.target.value); setPage(1) }}
                  />
                  <button className="ghost-button mobile-only" onClick={() => setFiltersOpen(true)}>Filters</button>
                </div>
              </div>
              <div className="toolbar-controls">
                <label className="sort-control">
                  <span>Sort</span>
                  <select value={sort} onChange={e => { setSort(e.target.value); setPage(1) }}>
                    <option value="smart">Smart sort</option>
                    <option value="year_desc">Newest year</option>
                    <option value="cells_desc">Most cells</option>
                    <option value="title_asc">Title A-Z</option>
                  </select>
                </label>
                <button className="clear-button" onClick={clearAll}>Clear all</button>
              </div>
            </section>

            <section className="active-filters">
              {activeDataset?.loaded && (
                <ActiveDatasetBanner
                  active={activeDataset}
                  syncError={syncError}
                  filteredDatasets={filteredDatasets}
                  allDatasets={allDatasets}
                  page={safePage}
                  onShow={showActiveDataset}
                />
              )}
              <div className="chip-row">
                {activeChips.length > 0 ? (
                  activeChips.map(chip => (
                    <button
                      key={chip.label}
                      className="chip"
                      data-kind={chip.kind}
                      onClick={() => removeChip(chip.kind, chip.label)}
                    >
                      {chip.label} ×
                    </button>
                  ))
                ) : (
                  <span className="result-summary">No active filters.</span>
                )}
              </div>
            </section>

            {statusMessage && (
              <section className="status-banner" data-tone={statusMessage.tone}>
                {statusMessage.text}
              </section>
            )}

            <section className="results-container" aria-live="polite">
              {catalogLoading && (
                <div className="loading-state">Loading catalog…</div>
              )}
              {!catalogLoading && catalogError && (
                <section className="empty-state">
                  <h3>The dataset browser could not load the catalog.</h3>
                  <p>{catalogError}</p>
                </section>
              )}
              {!catalogLoading && !catalogError && filteredDatasets.length === 0 && (
                <section className="empty-state">
                  <h3>No datasets match the current search and filters.</h3>
                  <p>Try broadening the query, clearing a facet, or widening the year and cell-count ranges.</p>
                </section>
              )}
              {!catalogLoading && !catalogError && filteredDatasets.length > 0 && (
                <>
                  <div className="results-head">
                    <span>Dataset</span>
                    <span>Author</span>
                    <span>Year</span>
                    <span>Tissue</span>
                    <span>Disease</span>
                    <span>Organism</span>
                    <span>Cells</span>
                  </div>
                  {visibleDatasets.map(d => {
                    const isSel = selectedRow === d.all_excel_row
                    const isActive = activeRow === d.all_excel_row
                    const cls = ["result-row", isSel && "is-active", isActive && "is-runtime-active", isSel && isActive && "is-active-runtime"].filter(Boolean).join(" ")
                    return (
                      <button key={d.all_excel_row} className={cls} onClick={() => loadDetail(d.all_excel_row)}>
                        <div className="result-primary">
                          <h3>{d.title || "Untitled dataset"}</h3>
                          <div className="badge-row">
                            {isActive && <span className="badge runtime-active-badge">Active dataset</span>}
                            <span className="badge">{d.author || "Unknown author"}</span>
                            <span className="badge">{d.journal || "No journal"}</span>
                          </div>
                        </div>
                        <div className="result-cell">{d.author || "—"}</div>
                        <div className="result-cell">{d.year || "—"}</div>
                        <div className="result-cell">{listPreview(d.tissues)}</div>
                        <div className="result-cell">{listPreview(d.diseases)}</div>
                        <div className="result-cell">{listPreview(d.organisms)}</div>
                        <div className="result-cell">{formatNumber(d.cell_counts)}</div>
                      </button>
                    )
                  })}
                  <div className="card-list">
                    {visibleDatasets.map(d => {
                      const isSel = selectedRow === d.all_excel_row
                      const isActive = activeRow === d.all_excel_row
                      const cls = ["card-row", isSel && "is-active", isActive && "is-runtime-active", isSel && isActive && "is-active-runtime"].filter(Boolean).join(" ")
                      return (
                        <button key={d.all_excel_row} className={cls} onClick={() => loadDetail(d.all_excel_row)}>
                          <div className="result-primary">
                            <h3>{d.title || "Untitled dataset"}</h3>
                            {isActive && <div className="badge-row"><span className="badge runtime-active-badge">Active dataset</span></div>}
                          </div>
                          <div className="metric"><span>Author</span><strong>{d.author || "—"}</strong></div>
                          <div className="metric"><span>Year</span><strong>{d.year || "—"}</strong></div>
                          <div className="metric"><span>Tissue</span><strong>{listPreview(d.tissues)}</strong></div>
                          <div className="metric"><span>Disease</span><strong>{listPreview(d.diseases)}</strong></div>
                          <div className="metric"><span>Organism</span><strong>{listPreview(d.organisms)}</strong></div>
                          <div className="metric"><span>Cells</span><strong>{formatNumber(d.cell_counts)}</strong></div>
                        </button>
                      )
                    })}
                  </div>
                </>
              )}
            </section>

            {/* Pagination */}
            {totalPages > 1 && (
              <nav className="pagination">
                <button
                  className="page-button"
                  onClick={() => { setPage(p => Math.max(1, p - 1)); window.scrollTo({ top: 0, behavior: "smooth" }) }}
                >
                  Previous
                </button>
                {Array.from({ length: Math.min(5, totalPages) }, (_, i) => {
                  const start = Math.max(1, safePage - 2)
                  const p = start + i
                  if (p > totalPages) return null
                  return (
                    <button
                      key={p}
                      className={`page-button ${p === safePage ? "is-current" : ""}`}
                      onClick={() => { setPage(p); window.scrollTo({ top: 0, behavior: "smooth" }) }}
                    >
                      {p}
                    </button>
                  )
                })}
                <button
                  className="page-button"
                  onClick={() => { setPage(p => Math.min(totalPages, p + 1)); window.scrollTo({ top: 0, behavior: "smooth" }) }}
                >
                  Next
                </button>
              </nav>
            )}
          </main>
            </div>
          </div>

          {/* Detail overlay */}
          <div
            className="detail-backdrop"
            onClick={() => { closeDetail(); setFiltersOpen(false) }}
          />
          <aside className="detail-drawer" aria-live="polite">
            <DetailDrawer
              selectedRow={selectedRow}
              currentDetail={currentDetail}
              detailLoading={detailLoading}
              detailError={detailError}
              selectedVariantRow={selectedVariantRow}
              activeParentRow={activeRow}
              activeVariantRow={activeVariantRow}
              detailActionLoading={detailActionLoading}
              detailActionMessage={detailActionMessage}
              detailActionTone={detailActionTone}
              onClose={closeDetail}
              onSelectVariant={row => { setSelectedVariantRow(row); setDetailActionMessage("") }}
              onAnalyze={handleAnalyze}
            />
          </aside>
        </div>
      </div>
    </div>
  )
}
