/**
 * GenoPixel Catalog Service
 * Wraps the catalog API: dataset listing, search/filter, active-dataset management.
 */

let CATALOG_API_URL = ""

async function loadCatalogApiUrl(): Promise<string> {
  if (CATALOG_API_URL) return CATALOG_API_URL
  const response = await fetch("/aws-exports.json")
  const config = await response.json()
  if (!config.catalogApiUrl) throw new Error("catalogApiUrl not found in aws-exports.json")
  // Ensure no trailing slash
  CATALOG_API_URL = config.catalogApiUrl.replace(/\/$/, "")
  return CATALOG_API_URL
}

function authHeaders(idToken: string): HeadersInit {
  return {
    "Content-Type": "application/json",
    Authorization: `Bearer ${idToken}`,
  }
}

// ---- Types ----

export interface DatasetVariant {
  multiple_excel_row: number
  file?: string
  tissue?: string
  disease?: string
  description?: string
  cell_count?: number
  cell_counts?: number
}

export interface Dataset {
  all_excel_row: number
  project: string
  doi?: string
  cellxgene_doi?: string
  title: string
  author?: string
  year?: number
  journal?: string
  tissues: string[]
  tissue_types: string[]
  diseases: string[]
  organisms: string[]
  cell_counts?: number
  merged: boolean
  primary_file?: string
  variant_count?: number
  variants?: DatasetVariant[]
}

export interface FacetCounts {
  project?: Record<string, number>
  organism?: Record<string, number>
  tissue?: Record<string, number>
  tissue_type?: Record<string, number>
  disease?: Record<string, number>
  journal?: Record<string, number>
  merged?: Record<string, number>
}

export interface CatalogResponse {
  datasets: Dataset[]
  total: number
  page: number
  page_size: number
  facets: FacetCounts
}

export interface ActiveDatasetResponse {
  loaded: boolean
  all_excel_row?: number
  multiple_excel_row?: number
  title?: string
  total_cells?: number
  primary_file?: string
  selected_at?: string
}

// ---- Catalog cache ----
// The full dataset list changes only when the metadata Excel is re-uploaded to S3.
// Cache it in memory for CATALOG_CACHE_TTL_MS so navigating away and back is instant.

const CATALOG_CACHE_TTL_MS = 5 * 60 * 1000

interface CatalogCacheEntry {
  response: CatalogResponse
  cachedAt: number
  idToken: string
}

let _catalogCache: CatalogCacheEntry | null = null

function getCachedCatalog(idToken: string): CatalogResponse | null {
  if (!_catalogCache) return null
  if (_catalogCache.idToken !== idToken) return null
  if (Date.now() - _catalogCache.cachedAt > CATALOG_CACHE_TTL_MS) return null
  return _catalogCache.response
}

// ---- API calls ----

export async function listDatasets(
  idToken: string,
  params: {
    search?: string
    project?: string
    organism?: string
    tissue?: string
    tissue_type?: string
    disease?: string
    merged?: boolean
    page?: number
    page_size?: number
  } = {}
): Promise<CatalogResponse> {
  const base = await loadCatalogApiUrl()
  const qs = new URLSearchParams()
  if (params.search) qs.set("search", params.search)
  if (params.project) qs.set("project", params.project)
  if (params.organism) qs.set("organism", params.organism)
  if (params.tissue) qs.set("tissue", params.tissue)
  if (params.tissue_type) qs.set("tissue_type", params.tissue_type)
  if (params.disease) qs.set("disease", params.disease)
  if (params.merged !== undefined) qs.set("merged", String(params.merged))
  if (params.page) qs.set("page", String(params.page))
  if (params.page_size) qs.set("page_size", String(params.page_size))

  const isFullFetch = (params.page_size ?? 0) >= 9999 && !params.search &&
    !params.organism && !params.tissue && !params.tissue_type &&
    !params.disease && !params.project && params.merged === undefined

  if (isFullFetch) {
    const cached = getCachedCatalog(idToken)
    if (cached) return cached
  }

  const url = `${base}/api/catalog${qs.toString() ? "?" + qs.toString() : ""}`
  const resp = await fetch(url, { headers: authHeaders(idToken) })
  if (!resp.ok) throw new Error(`Catalog API error: ${resp.status}`)
  const data: CatalogResponse = await resp.json()

  if (isFullFetch) {
    _catalogCache = { response: data, cachedAt: Date.now(), idToken }
  }

  return data
}

export async function getDataset(
  idToken: string,
  allExcelRow: number
): Promise<Dataset> {
  const base = await loadCatalogApiUrl()
  const resp = await fetch(`${base}/api/catalog/${allExcelRow}`, {
    headers: authHeaders(idToken),
  })
  if (!resp.ok) throw new Error(`Catalog API error: ${resp.status}`)
  return resp.json()
}

export async function analyzeDataset(
  idToken: string,
  allExcelRow: number,
  multipleExcelRow?: number
): Promise<{ message: string; all_excel_row: number }> {
  const base = await loadCatalogApiUrl()
  const body: Record<string, unknown> = {}
  if (multipleExcelRow !== undefined) body.multiple_excel_row = multipleExcelRow

  const resp = await fetch(`${base}/api/catalog/${allExcelRow}/analyze`, {
    method: "POST",
    headers: authHeaders(idToken),
    body: JSON.stringify(body),
  })
  if (!resp.ok) throw new Error(`Catalog API error: ${resp.status}`)
  return resp.json()
}

export async function getActiveDataset(
  idToken: string
): Promise<ActiveDatasetResponse> {
  const base = await loadCatalogApiUrl()
  const resp = await fetch(`${base}/api/catalog/active-dataset`, {
    headers: authHeaders(idToken),
  })
  if (!resp.ok) throw new Error(`Catalog API error: ${resp.status}`)
  return resp.json()
}
