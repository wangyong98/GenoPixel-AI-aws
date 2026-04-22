"use client"

import { useState } from "react"
import { ThumbsUp, ThumbsDown } from "lucide-react"
import { Message } from "./types"
import { FeedbackDialog } from "./FeedbackDialog"
import { getToolRenderer } from "@/hooks/useToolRenderer"
import { MarkdownRenderer } from "./MarkdownRenderer"

interface ChatMessageProps {
  message: Message
  sessionId: string
  onFeedbackSubmit: (feedbackType: "positive" | "negative", comment: string) => Promise<void>
}

export function ChatMessage({
  message,
  sessionId: _sessionId,
  onFeedbackSubmit,
}: ChatMessageProps) {
  const [isDialogOpen, setIsDialogOpen] = useState(false)
  const [selectedFeedbackType, setSelectedFeedbackType] = useState<"positive" | "negative">(
    "positive"
  )
  const [feedbackSubmitted, setFeedbackSubmitted] = useState(false)

  const formatTime = (timestamp: string) => {
    return new Date(timestamp).toLocaleTimeString([], {
      hour: "2-digit",
      minute: "2-digit",
    })
  }

  const handleFeedbackClick = (type: "positive" | "negative") => {
    setSelectedFeedbackType(type)
    setIsDialogOpen(true)
  }

  const handleFeedbackSubmit = async (comment: string) => {
    await onFeedbackSubmit(selectedFeedbackType, comment)
    setFeedbackSubmitted(true)
  }

  const renderAssistantContent = () => {
    // If segments exist, render them in order (interleaved text + tools)
    if (message.segments && message.segments.length > 0) {
      return message.segments.map((seg, i) => {
        if (seg.type === "text") {
          return <MarkdownRenderer key={i} content={seg.content} />
        }
        const render = getToolRenderer(seg.toolCall.name)
        if (!render) return null
        return (
          <div key={seg.toolCall.toolUseId} className="my-1">
            {render({
              name: seg.toolCall.name,
              args: seg.toolCall.input,
              status: seg.toolCall.status,
              result: seg.toolCall.result,
            })}
          </div>
        )
      })
    }
    // Fallback: just render content as markdown
    return <MarkdownRenderer content={message.content} />
  }

  return (
    <div className={`flex flex-col ${message.role === "user" ? "items-end" : "items-start"}`}>
      <div
        className={`max-w-[80%] break-words ${
          message.role === "user"
            ? "p-3 rounded-lg rounded-br-none whitespace-pre-wrap"
            : ""
        }`}
        style={
          message.role === "user"
            ? { background: "#1a3550", color: "#ffffff" }
            : { color: "#1e2d3d" }
        }
      >
        {message.role === "assistant" ? renderAssistantContent() : message.content}
      </div>

      {/* Timestamp and Feedback buttons for assistant messages */}
      <div className="flex items-center gap-2 mt-1 px-1">
        <div className="text-xs" style={{ color: "#5f6d75" }}>
          {formatTime(message.timestamp)}
        </div>

        {/* Show feedback buttons only for assistant messages with content */}
        {message.role === "assistant" && message.content && (
          <div className="flex items-center gap-1 ml-2">
            <button
              onClick={() => handleFeedbackClick("positive")}
              disabled={feedbackSubmitted}
              className="p-1 rounded transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
              style={{ color: "#5f6d75" }}
              onMouseEnter={e => (e.currentTarget.style.color = "#0d6f68")}
              onMouseLeave={e => (e.currentTarget.style.color = "#5f6d75")}
              aria-label="Positive feedback"
              title="Good response"
            >
              <ThumbsUp size={14} />
            </button>
            <button
              onClick={() => handleFeedbackClick("negative")}
              disabled={feedbackSubmitted}
              className="p-1 rounded transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
              style={{ color: "#5f6d75" }}
              onMouseEnter={e => (e.currentTarget.style.color = "#b67b2e")}
              onMouseLeave={e => (e.currentTarget.style.color = "#5f6d75")}
              aria-label="Negative feedback"
              title="Bad response"
            >
              <ThumbsDown size={14} />
            </button>
            {feedbackSubmitted && (
              <span className="text-xs ml-1" style={{ color: "#5f6d75" }}>
                Thanks for your feedback!
              </span>
            )}
          </div>
        )}
      </div>

      {/* Feedback Dialog */}
      <FeedbackDialog
        isOpen={isDialogOpen}
        onClose={() => setIsDialogOpen(false)}
        onSubmit={handleFeedbackSubmit}
        feedbackType={selectedFeedbackType}
      />
    </div>
  )
}
