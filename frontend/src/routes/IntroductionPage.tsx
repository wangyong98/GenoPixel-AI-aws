"use client"

import { AppSidebar } from "@/components/layout/AppSidebar"

const EXAMPLE_PROMPTS = [
  "Show me the data.",
  "Load the active dataset and summarize the key metadata.",
  "Plot UMAP colored by author_cell_type.",
  "Show top cell counts by author_cell_type.",
  "Generate a dot plot for CD3E, CD4, CD8A grouped by author_cell_type.",
  "Compare cell type proportions by disease.",
]

const cardStyle: React.CSSProperties = {
  background: "#ffffff",
  border: "1px solid #dce4ea",
  borderRadius: "12px",
  padding: "1.75rem 2rem",
  boxShadow: "0 2px 12px rgba(26, 53, 80, 0.06)",
}

const kicerStyle: React.CSSProperties = {
  fontSize: "0.7rem",
  fontWeight: 700,
  textTransform: "uppercase",
  letterSpacing: "0.18em",
  color: "#0d6f68",
  marginBottom: "0.6rem",
}

export default function IntroductionPage() {
  return (
    <div
      className="flex h-screen"
      style={{ fontFamily: '"Avenir Next", "Avenir", "Segoe UI", "Helvetica Neue", Arial, sans-serif' }}
    >
      <AppSidebar activeTab="introduction" />
      <main className="flex-1 min-w-0 overflow-auto">
        <div className="max-w-5xl mx-auto px-8 py-8">
        <section style={cardStyle}>
          <p style={kicerStyle}>Overview</p>
          <h1
            className="text-2xl font-bold"
            style={{ color: "#1a3550", letterSpacing: "-0.01em", marginBottom: "0.85rem" }}
          >
            GenoPixel AI — Single-Cell Analysis Assistant
          </h1>
          <p className="leading-7 text-[0.97rem]" style={{ color: "#3a5068" }}>
            GenoPixel AI is an AI-powered assistant for exploring curated single-cell RNA-seq
            datasets in AnnData (<code
              className="px-1.5 py-0.5 rounded text-[0.88em]"
              style={{ background: "#e4f2f1", color: "#0d6f68", fontFamily: "'IBM Plex Mono', monospace" }}
            >.h5ad</code>) format. It integrates with Scanpy to generate publication-quality
            visualizations directly from natural language requests.
          </p>
          <p className="mt-3 leading-7 text-[0.97rem]" style={{ color: "#3a5068" }}>
            Designed for scientists and clinicians working with single-cell genomics data in
            biotech, biopharma, and research settings — no programming required.
          </p>
        </section>

        <section className="mt-4" style={cardStyle}>
          <p style={kicerStyle}>Workflow</p>
          <ol className="mt-1 space-y-3">
            {[
              { step: "1", label: "Select a dataset", desc: "Browse the catalog in Datasets and click Analyze to set your active dataset." },
              { step: "2", label: "Open Chat", desc: "The active dataset is automatically pre-loaded when you start a conversation." },
              { step: "3", label: "Ask questions", desc: "Request UMAP plots, cell type summaries, gene expression comparisons, and more." },
            ].map(({ step, label, desc }) => (
              <li key={step} className="flex gap-3 items-start">
                <span
                  className="shrink-0 w-6 h-6 rounded-full flex items-center justify-center text-xs font-bold mt-0.5"
                  style={{ background: "#0d6f68", color: "#fff" }}
                >
                  {step}
                </span>
                <div>
                  <span className="font-semibold text-[0.95rem]" style={{ color: "#1a3550" }}>
                    {label}
                  </span>
                  <span className="text-[0.93rem]" style={{ color: "#5b7389" }}>
                    {" — "}{desc}
                  </span>
                </div>
              </li>
            ))}
          </ol>
        </section>

        <section className="mt-4" style={cardStyle}>
          <p style={kicerStyle}>Example Prompts</p>
          <ul className="mt-1 space-y-2">
            {EXAMPLE_PROMPTS.map(prompt => (
              <li key={prompt} className="flex items-center gap-2.5">
                <span
                  className="shrink-0 w-1.5 h-1.5 rounded-full"
                  style={{ background: "#0d6f68" }}
                />
                <code
                  className="text-[0.88em] px-2.5 py-1 rounded-md leading-relaxed"
                  style={{
                    background: "#f1f5f8",
                    color: "#1e2d3d",
                    border: "1px solid #dce4ea",
                    fontFamily: "'IBM Plex Mono', 'SFMono-Regular', Consolas, monospace",
                  }}
                >
                  {prompt}
                </code>
              </li>
            ))}
          </ul>
        </section>
        </div>
      </main>
    </div>
  )
}
