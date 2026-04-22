"use client"
// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

import ChatInterface from "@/components/chat/ChatInterface"
import type { ChatControls } from "@/components/chat/ChatInterface"
import { AppSidebar } from "@/components/layout/AppSidebar"
import { Button } from "@/components/ui/button"
import { useAuth } from "@/hooks/useAuth"
import { GlobalContextProvider } from "@/app/context/GlobalContext"
import { useState } from "react"

export default function ChatPage() {
  const { isAuthenticated, signIn } = useAuth()
  const [chatControls, setChatControls] = useState<ChatControls>({
    onNewChat: () => {},
    canStartNewChat: false,
    recentChats: [],
    currentSessionId: "",
    onSelectChat: () => {},
  })

  if (!isAuthenticated) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen gap-4">
        <p className="text-4xl">Please sign in</p>
        <Button onClick={() => signIn()}>Sign In</Button>
      </div>
    )
  }

  return (
    <GlobalContextProvider>
      <div className="flex h-screen">
        <AppSidebar
          activeTab="chat"
          showChatActions
          onNewChat={chatControls.onNewChat}
          canStartNewChat={chatControls.canStartNewChat}
          recentChats={chatControls.recentChats}
          currentSessionId={chatControls.currentSessionId}
          onSelectChat={chatControls.onSelectChat}
        />
        <div className="relative flex-1 min-w-0 overflow-hidden">
          <ChatInterface onChatControlsChange={setChatControls} />
        </div>
      </div>
    </GlobalContextProvider>
  )
}
