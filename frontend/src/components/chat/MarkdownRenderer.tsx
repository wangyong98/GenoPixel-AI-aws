"use client"

import { useState } from "react"
import ReactMarkdown, { defaultUrlTransform } from "react-markdown"
import remarkGfm from "remark-gfm"
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter"
import { oneLight } from "react-syntax-highlighter/dist/esm/styles/prism"
import { Copy, Check } from "lucide-react"

function completePartialMarkdown(text: string): string {
  const fenceCount = (text.match(/^```/gm) || []).length
  if (fenceCount % 2 !== 0) return text + "\n```"
  return text
}

function normalizeDataImageMarkdown(text: string): string {
  // Some model responses split image markdown across lines:
  // ![alt]
  // (data:image/png;base64,...)
  // Re-join it so markdown can parse and render the image.
  const rejoined = text.replace(
    /!\[([^\]]*)\]\s*\n+\s*\((data:image\/[a-zA-Z0-9.+-]+;base64,[^)]+)\)/gi,
    "![$1]($2)"
  )

  // If the base64 payload itself was split with whitespace/newlines, collapse it.
  return rejoined.replace(
    /\((data:image\/[a-zA-Z0-9.+-]+;base64,[a-zA-Z0-9+/=\s]+)\)/gi,
    (_, dataUri: string) => `(${dataUri.replace(/\s+/g, "")})`
  )
}

function isSafeDataImageUrl(url: string): boolean {
  return /^data:image\/[a-zA-Z0-9.+-]+;base64,[a-zA-Z0-9+/=]+$/i.test(url)
}

function markdownUrlTransform(url: string): string {
  const trimmed = url.trim()
  if (isSafeDataImageUrl(trimmed)) return trimmed
  return defaultUrlTransform(trimmed)
}

function CopyButton({ text }: { text: string }) {
  const [copied, setCopied] = useState(false)
  const handleCopy = () => {
    navigator.clipboard.writeText(text)
    setCopied(true)
    setTimeout(() => setCopied(false), 2000)
  }
  return (
    <button
      onClick={handleCopy}
      className="p-1 text-gray-400 hover:text-gray-600 transition-colors"
      aria-label="Copy code"
    >
      {copied ? <Check size={14} /> : <Copy size={14} />}
    </button>
  )
}

// react-markdown v10 + React 19 has overly strict component types for element-specific refs.
// Using Record<string, ...> to avoid the type mismatch on pre, p, th, td, etc.
// eslint-disable-next-line @typescript-eslint/no-explicit-any
const components: Record<string, any> = {
  code({ className, children }: { className?: string; children?: React.ReactNode }) {
    const match = /language-(\w+)/.exec(className || "")
    const codeString = String(children).replace(/\n$/, "")
    if (match) {
      return (
        <div className="my-2 rounded-md overflow-hidden border border-gray-300 bg-white">
          <div className="flex items-center justify-between px-3 py-1 bg-gray-100 border-b border-gray-300">
            <span className="text-xs text-gray-500">{match[1]}</span>
            <CopyButton text={codeString} />
          </div>
          <SyntaxHighlighter
            style={oneLight}
            language={match[1]}
            PreTag="div"
            customStyle={{
              margin: 0,
              padding: "0.75rem",
              fontSize: "0.8rem",
              background: "white",
            }}
          >
            {codeString}
          </SyntaxHighlighter>
        </div>
      )
    }
    return (
      <code className="px-1 py-0.5 bg-gray-200/60 rounded text-[0.85em] font-mono">{children}</code>
    )
  },
  pre({ children }: { children?: React.ReactNode }) {
    return <>{children}</>
  },
  img({ src, alt }: { src?: string; alt?: string }) {
    if (!src) return null
    return (
      <img
        src={src}
        alt={alt || "Generated figure"}
        className="my-2 max-w-full h-auto rounded-md border border-gray-200"
        loading="lazy"
      />
    )
  },
}

export function MarkdownRenderer({ content }: { content: string }) {
  if (!content) return null
  const normalized = completePartialMarkdown(normalizeDataImageMarkdown(content))
  return (
    <div className="markdown-body leading-relaxed [&_p]:my-1.5 [&_ul]:my-1.5 [&_ul]:pl-5 [&_ul]:list-disc [&_ol]:my-1.5 [&_ol]:pl-5 [&_ol]:list-decimal [&_li]:my-0.5 [&_h1]:text-lg [&_h1]:font-semibold [&_h1]:mt-3 [&_h1]:mb-1.5 [&_h2]:text-base [&_h2]:font-semibold [&_h2]:mt-2.5 [&_h2]:mb-1 [&_h3]:text-sm [&_h3]:font-semibold [&_h3]:mt-2 [&_h3]:mb-1 [&_blockquote]:border-l-2 [&_blockquote]:border-gray-300 [&_blockquote]:pl-3 [&_blockquote]:my-1.5 [&_blockquote]:text-gray-600 [&_table]:my-2 [&_table]:min-w-full [&_table]:border-collapse [&_table]:text-xs [&_th]:px-2 [&_th]:py-1 [&_th]:bg-gray-100 [&_th]:border [&_th]:border-gray-300 [&_th]:text-left [&_th]:font-medium [&_td]:px-2 [&_td]:py-1 [&_td]:border [&_td]:border-gray-300 [&>*:first-child]:mt-0 [&>*:last-child]:mb-0">
      <ReactMarkdown remarkPlugins={[remarkGfm]} components={components} urlTransform={markdownUrlTransform}>
        {normalized}
      </ReactMarkdown>
    </div>
  )
}
