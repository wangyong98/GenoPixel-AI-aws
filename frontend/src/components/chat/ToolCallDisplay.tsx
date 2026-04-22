"use client"

import { useState } from "react"
import { Wrench, Loader2, CheckCircle2, ChevronRight, ChevronDown } from "lucide-react"
import type { ToolRenderProps } from "@/hooks/useToolRenderer"

function extractDataImageUri(text?: string): string | null {
  if (!text) return null

  const normalized = text.replace(/\\n/g, "\n").replace(/\\"/g, '"')
  const match = normalized.match(/data:image\/[a-zA-Z0-9.+-]+;base64,[a-zA-Z0-9+/=\s]+/i)
  if (!match) return null

  return match[0].replace(/\s+/g, "")
}

export function ToolCallDisplay({ name, args, status, result }: ToolRenderProps) {
  const [expanded, setExpanded] = useState(false)
  const figureSrc = extractDataImageUri(result)

  return (
    <div className="my-1 text-sm">
      <button
        onClick={() => setExpanded(!expanded)}
        className="flex items-center gap-1.5 px-2 py-1 rounded hover:bg-gray-200/50 transition-colors w-full text-left"
      >
        {expanded ? (
          <ChevronDown size={12} className="text-gray-400" />
        ) : (
          <ChevronRight size={12} className="text-gray-400" />
        )}
        <Wrench size={12} className="text-gray-400" />
        <span className="text-gray-600">{name}</span>
        {status === "streaming" && (
          <Loader2 size={12} className="animate-spin text-blue-500 ml-auto" />
        )}
        {status === "executing" && (
          <Loader2 size={12} className="animate-spin text-amber-500 ml-auto" />
        )}
        {status === "complete" && <CheckCircle2 size={12} className="text-green-500 ml-auto" />}
      </button>

      {status === "complete" && figureSrc && (
        <div className="ml-6 mt-2">
          <img
            src={figureSrc}
            alt={`${name} figure`}
            className="max-w-full h-auto rounded-md border border-gray-200 bg-white"
            loading="lazy"
          />
        </div>
      )}

      {expanded && (
        <div className="ml-6 mt-1 border-l-2 border-gray-200 pl-3 space-y-2">
          {args && (
            <div>
              <div className="text-xs text-gray-400">Input</div>
              <pre className="text-xs text-gray-600 whitespace-pre-wrap break-words mt-0.5">
                {args}
              </pre>
            </div>
          )}
          {result && (
            <div>
              <div className="text-xs text-gray-400">Result</div>
              <pre className="text-xs text-gray-600 whitespace-pre-wrap break-words mt-0.5">
                {result}
              </pre>
            </div>
          )}
        </div>
      )}
    </div>
  )
}
