"use client"

import { useCallback, useEffect, useRef, useState } from "react"
import { ChatInput } from "./ChatInput"
import { ChatMessages } from "./ChatMessages"
import { Message, MessageSegment, ToolCall } from "./types"

import { useGlobal } from "@/app/context/GlobalContext"
import { AgentCoreClient } from "@/lib/agentcore-client"
import type { AgentPattern } from "@/lib/agentcore-client"
import { submitFeedback } from "@/services/feedbackService"
import { useAuth } from "react-oidc-context"
import { useDefaultTool } from "@/hooks/useToolRenderer"
import { ToolCallDisplay } from "./ToolCallDisplay"

const RECENT_CHATS_STORAGE_KEY = "genopixel_recent_chats_v1"
const MAX_RECENT_CHATS = 20

type StoredRecentChat = {
  id: string
  name: string
  messages: Message[]
  updatedAt: string
}

export type RecentChatSummary = {
  id: string
  name: string
  updatedAt: string
}

export type ChatControls = {
  onNewChat: () => void
  canStartNewChat: boolean
  recentChats: RecentChatSummary[]
  currentSessionId: string
  onSelectChat: (chatId: string) => void
}

interface ChatInterfaceProps {
  onChatControlsChange?: (controls: ChatControls) => void
}

export default function ChatInterface({ onChatControlsChange }: ChatInterfaceProps) {
  const [messages, setMessages] = useState<Message[]>([])
  const [input, setInput] = useState("")
  const [error, setError] = useState<string | null>(null)
  const [client, setClient] = useState<AgentCoreClient | null>(null)
  const [sessionId, setSessionId] = useState<string>(() => crypto.randomUUID())
  const [recentChats, setRecentChats] = useState<StoredRecentChat[]>([])

  const { isLoading, setIsLoading } = useGlobal()
  const auth = useAuth()

  // Ref for message container to enable auto-scrolling
  const messagesEndRef = useRef<HTMLDivElement>(null)

  // Register default tool renderer (wildcard "*")
  useDefaultTool(({ name, args, status, result }) => (
    <ToolCallDisplay name={name} args={args} status={status} result={result} />
  ))

  const sessionTitleFromMessages = useCallback((history: Message[]) => {
    const firstUserText =
      history.find(m => m.role === "user" && m.content.trim().length > 0)?.content?.trim() || "New chat"
    return firstUserText.length > 50 ? `${firstUserText.slice(0, 50)}…` : firstUserText
  }, [])

  const persistRecentChats = useCallback((value: StoredRecentChat[]) => {
    try {
      window.localStorage.setItem(RECENT_CHATS_STORAGE_KEY, JSON.stringify(value))
    } catch (err) {
      console.error("Failed to persist recent chats:", err)
    }
  }, [])

  const upsertRecentChat = useCallback((chatId: string, history: Message[]) => {
    if (!history.length) return
    const name = sessionTitleFromMessages(history)
    const updatedAt = history[history.length - 1]?.timestamp || new Date().toISOString()

    setRecentChats(prev => {
      const next = [
        { id: chatId, name, messages: history, updatedAt },
        ...prev.filter(c => c.id !== chatId),
      ]
        .sort((a, b) => b.updatedAt.localeCompare(a.updatedAt))
        .slice(0, MAX_RECENT_CHATS)

      persistRecentChats(next)
      return next
    })
  }, [persistRecentChats, sessionTitleFromMessages])

  const onSelectRecentChat = useCallback((chatId: string) => {
    const target = recentChats.find(c => c.id === chatId)
    if (!target) return
    setError(null)
    setInput("")
    setSessionId(target.id)
    setMessages(target.messages || [])
  }, [recentChats])

  // Load agent configuration and create client on mount
  useEffect(() => {
    async function loadConfig() {
      try {
        const response = await fetch("/aws-exports.json")
        if (!response.ok) {
          throw new Error("Failed to load configuration")
        }
        const config = await response.json()

        if (!config.agentRuntimeArn) {
          throw new Error("Agent Runtime ARN not found in configuration")
        }

        const agentClient = new AgentCoreClient({
          runtimeArn: config.agentRuntimeArn,
          region: config.awsRegion || "us-east-1",
          pattern: (config.agentPattern || "strands-single-agent") as AgentPattern,
        })

        setClient(agentClient)
      } catch (err) {
        const errorMessage = err instanceof Error ? err.message : "Unknown error"
        setError(`Configuration error: ${errorMessage}`)
        console.error("Failed to load agent configuration:", err)
      }
    }

    loadConfig()
  }, [])

  useEffect(() => {
    try {
      const raw = window.localStorage.getItem(RECENT_CHATS_STORAGE_KEY)
      if (!raw) return
      const parsed = JSON.parse(raw) as StoredRecentChat[]
      if (!Array.isArray(parsed)) return
      const cleaned = parsed
        .filter(c => c && typeof c.id === "string" && Array.isArray(c.messages))
        .sort((a, b) => (b.updatedAt || "").localeCompare(a.updatedAt || ""))
        .slice(0, MAX_RECENT_CHATS)
      setRecentChats(cleaned)
    } catch (err) {
      console.error("Failed to load recent chats:", err)
    }
  }, [])

  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: "smooth" })
  }, [messages])

  const sendMessage = async (userMessage: string) => {
    if (!userMessage.trim() || !client) return

    // Clear any previous errors
    setError(null)

    // Add user message to chat
    const newUserMessage: Message = {
      role: "user",
      content: userMessage,
      timestamp: new Date().toISOString(),
    }

    setMessages(prev => [...prev, newUserMessage])
    setInput("")
    setIsLoading(true)

    // Create placeholder for assistant response
    const assistantResponse: Message = {
      role: "assistant",
      content: "",
      timestamp: new Date().toISOString(),
    }

    setMessages(prev => [...prev, assistantResponse])

    try {
      // Get auth token from react-oidc-context
      const accessToken = auth.user?.access_token

      if (!accessToken) {
        throw new Error("Authentication required. Please log in again.")
      }

      const segments: MessageSegment[] = []
      const toolCallMap = new Map<string, ToolCall>()

      const updateMessage = () => {
        // Build content from text segments for backward compat
        const content = segments
          .filter((s): s is Extract<MessageSegment, { type: "text" }> => s.type === "text")
          .map(s => s.content)
          .join("")

        setMessages(prev => {
          const updated = [...prev]
          updated[updated.length - 1] = {
            ...updated[updated.length - 1],
            content,
            segments: [...segments],
          }
          return updated
        })
      }

      // User identity is extracted server-side from the validated JWT token,
      // not passed as a parameter — prevents impersonation via prompt injection.
      await client.invoke(userMessage, sessionId, accessToken, event => {
        switch (event.type) {
          case "text": {
            // If text arrives after a tool segment, mark all pending tools as complete
            const prev = segments[segments.length - 1]
            if (prev && prev.type === "tool") {
              for (const tc of toolCallMap.values()) {
                if (tc.status === "streaming" || tc.status === "executing") {
                  tc.status = "complete"
                }
              }
            }
            // Append to last text segment, or create new one
            const last = segments[segments.length - 1]
            if (last && last.type === "text") {
              last.content += event.content
            } else {
              segments.push({ type: "text", content: event.content })
            }
            updateMessage()
            break
          }
          case "tool_use_start": {
            const tc: ToolCall = {
              toolUseId: event.toolUseId,
              name: event.name,
              input: "",
              status: "streaming",
            }
            toolCallMap.set(event.toolUseId, tc)
            segments.push({ type: "tool", toolCall: tc })
            updateMessage()
            break
          }
          case "tool_use_delta": {
            const tc = toolCallMap.get(event.toolUseId)
            if (tc) {
              tc.input += event.input
            }
            updateMessage()
            break
          }
          case "tool_result": {
            const tc = toolCallMap.get(event.toolUseId)
            if (tc) {
              tc.result = event.result
              tc.status = "complete"
            }
            updateMessage()
            break
          }
          case "message": {
            if (event.role === "assistant") {
              for (const tc of toolCallMap.values()) {
                if (tc.status === "streaming") tc.status = "executing"
              }
              updateMessage()
            }
            break
          }
        }
      })
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : "Unknown error"
      setError(`Failed to get response: ${errorMessage}`)
      console.error("Error invoking AgentCore:", err)

      // Update the assistant message with error
      setMessages(prev => {
        const updated = [...prev]
        updated[updated.length - 1] = {
          ...updated[updated.length - 1],
          content:
            "I apologize, but I encountered an error processing your request. Please try again.",
        }
        return updated
      })
    } finally {
      setIsLoading(false)
    }
  }

  // Handle form submission
  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault()

    sendMessage(input)
  }

  // Handle feedback submission
  const handleFeedbackSubmit = async (
    messageContent: string,
    feedbackType: "positive" | "negative",
    comment: string
  ) => {
    try {
      // Use ID token for API Gateway Cognito authorizer (not access token)
      const idToken = auth.user?.id_token

      if (!idToken) {
        throw new Error("Authentication required. Please log in again.")
      }

      await submitFeedback(
        {
          sessionId,
          message: messageContent,
          feedbackType,
          comment: comment || undefined,
        },
        idToken
      )

      console.log("Feedback submitted successfully")
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : "Unknown error"
      console.error("Error submitting feedback:", err)
      setError(`Failed to submit feedback: ${errorMessage}`)
    }
  }

  // Start a new chat by clearing messages and generating a fresh session ID.
  // A new UUID is required so the backend treats this as a distinct conversation context.
  const startNewChat = useCallback(() => {
    setMessages([])
    setInput("")
    setError(null)
    setSessionId(crypto.randomUUID())
  }, [])

  // Check if this is the initial state (no messages)
  const isInitialState = messages.length === 0

  // Check if there are any assistant messages
  const hasAssistantMessages = messages.some(message => message.role === "assistant")

  useEffect(() => {
    if (!messages.length) return
    upsertRecentChat(sessionId, messages)
  }, [sessionId, messages, upsertRecentChat])

  useEffect(() => {
    onChatControlsChange?.({
      onNewChat: startNewChat,
      canStartNewChat: hasAssistantMessages,
      recentChats: recentChats.map(c => ({ id: c.id, name: c.name, updatedAt: c.updatedAt })),
      currentSessionId: sessionId,
      onSelectChat: onSelectRecentChat,
    })
  }, [
    onChatControlsChange,
    startNewChat,
    hasAssistantMessages,
    recentChats,
    sessionId,
    onSelectRecentChat,
  ])

  return (
    <div className="flex flex-col h-full w-full">
      {error && (
        <div className="bg-red-50 border-l-4 border-red-500 p-4 mx-4 mt-2 flex-none">
          <p className="text-sm text-red-700">{error}</p>
        </div>
      )}

      {/* Conditional layout based on whether there are messages */}
      {isInitialState ? (
        // Initial state - input in the middle
        <>
          {/* Empty space above */}
          <div className="grow" />

          {/* Centered welcome message */}
          <div className="text-center mb-6">
            <h2 className="text-2xl font-bold text-gray-800">Welcome to GenoPixel Chat</h2>
            <p className="text-gray-600 mt-2">Ask me "show me the data" to get started</p>
          </div>

          {/* Centered input */}
          <div className="px-4 mb-16 max-w-4xl mx-auto w-full">
            <ChatInput
              input={input}
              setInput={setInput}
              handleSubmit={handleSubmit}
              isLoading={isLoading}
            />
          </div>

          {/* Empty space below */}
          <div className="grow" />
        </>
      ) : (
        // Chat in progress - normal layout
        <>
          {/* Scrollable message area */}
          <div className="grow overflow-hidden">
            <div className="max-w-4xl mx-auto w-full h-full">
              <ChatMessages
                messages={messages}
                messagesEndRef={messagesEndRef}
                sessionId={sessionId}
                onFeedbackSubmit={handleFeedbackSubmit}
              />
            </div>
          </div>

          {/* Fixed input area at bottom */}
          <div className="flex-none">
            <div className="max-w-4xl mx-auto w-full">
              <ChatInput
                input={input}
                setInput={setInput}
                handleSubmit={handleSubmit}
                isLoading={isLoading}
              />
            </div>
          </div>
        </>
      )}
    </div>
  )
}
