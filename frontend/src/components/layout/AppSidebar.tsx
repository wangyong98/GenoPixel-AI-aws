import { Button } from "@/components/ui/button"
import { Plus, MessageSquare, Database, BookOpenText, LogOut } from "lucide-react"
import { useAuth } from "@/hooks/useAuth"
import { useNavigate } from "react-router-dom"
import {
  AlertDialog,
  AlertDialogAction,
  AlertDialogCancel,
  AlertDialogContent,
  AlertDialogDescription,
  AlertDialogFooter,
  AlertDialogHeader,
  AlertDialogTitle,
  AlertDialogTrigger,
} from "@/components/ui/alert-dialog"

export type AppTab = "introduction" | "datasets" | "chat"

type RecentChatSummary = {
  id: string
  name: string
  updatedAt: string
}

type AppSidebarProps = {
  activeTab: AppTab
  onNewChat?: () => void
  canStartNewChat?: boolean
  showChatActions?: boolean
  recentChats?: RecentChatSummary[]
  currentSessionId?: string
  onSelectChat?: (chatId: string) => void
}

export function AppSidebar({
  activeTab,
  onNewChat,
  canStartNewChat = false,
  showChatActions = false,
  recentChats = [],
  currentSessionId,
  onSelectChat,
}: AppSidebarProps) {
  const { isAuthenticated, signOut } = useAuth()
  const navigate = useNavigate()

  return (
    <aside
      className="w-16 md:w-64 shrink-0 flex flex-col"
      style={{
        background: "#ffffff",
        borderRight: "1px solid #e4eaef",
        fontFamily: '"Avenir Next", "Avenir", "Segoe UI", "Helvetica Neue", Arial, sans-serif',
      }}
    >
      <div className="px-4 py-5" style={{ borderBottom: "1px solid #e4eaef" }}>
        <p
          className="text-xs font-bold uppercase tracking-[0.20em]"
          style={{ color: "#0d6f68" }}
        >
          <span className="md:hidden">GP</span>
          <span className="hidden md:inline">GenoPixel</span>
        </p>
      </div>

      <nav className="p-3 space-y-1">
        <Button
          variant={activeTab === "introduction" ? "default" : "ghost"}
          className="w-full justify-start"
          onClick={() => navigate("/")}
        >
          <BookOpenText className="h-4 w-4 mr-2" />
          <span className="hidden md:inline">Introduction</span>
        </Button>
        <Button
          variant={activeTab === "datasets" ? "default" : "ghost"}
          className="w-full justify-start"
          onClick={() => navigate("/datasets")}
        >
          <Database className="h-4 w-4 mr-2" />
          <span className="hidden md:inline">Datasets</span>
        </Button>
        <Button
          variant={activeTab === "chat" ? "default" : "ghost"}
          className="w-full justify-start"
          onClick={() => navigate("/chat")}
        >
          <MessageSquare className="h-4 w-4 mr-2" />
          <span className="hidden md:inline">Chat</span>
        </Button>
      </nav>

      {showChatActions && (
        <div className="px-3 pb-2 space-y-3">
          <Button
            onClick={onNewChat}
            variant="outline"
            className="w-full justify-start"
            disabled={!canStartNewChat}
          >
            <Plus className="h-4 w-4 mr-2" />
            <span className="hidden md:inline">New Chat</span>
          </Button>

          <div className="hidden md:block">
            <p
              className="px-2 text-[11px] font-semibold uppercase tracking-[0.16em]"
              style={{ color: "#5f6d75" }}
            >
              Recents
            </p>
            <div className="mt-1 space-y-1">
              {recentChats.length === 0 && (
                <p className="px-2 py-1 text-xs" style={{ color: "#5f6d75" }}>
                  No recent chats yet.
                </p>
              )}
              {recentChats.map(chat => (
                <Button
                  key={chat.id}
                  variant={chat.id === currentSessionId ? "secondary" : "ghost"}
                  className="w-full justify-start text-left h-auto py-2"
                  onClick={() => onSelectChat?.(chat.id)}
                  title={chat.name}
                >
                  <span className="truncate text-sm">{chat.name}</span>
                </Button>
              ))}
            </div>
          </div>
        </div>
      )}

      <div className="mt-auto p-3" style={{ borderTop: "1px solid #e4eaef" }}>
        {isAuthenticated && (
          <AlertDialog>
            <AlertDialogTrigger asChild>
              <Button variant="outline" className="w-full justify-start">
                <LogOut className="h-4 w-4 mr-2" />
                <span className="hidden md:inline">Logout</span>
              </Button>
            </AlertDialogTrigger>
            <AlertDialogContent>
              <AlertDialogHeader>
                <AlertDialogTitle>Confirm Logout</AlertDialogTitle>
                <AlertDialogDescription>
                  Are you sure you want to log out? You will need to sign in again to access your
                  account.
                </AlertDialogDescription>
              </AlertDialogHeader>
              <AlertDialogFooter>
                <AlertDialogCancel>Cancel</AlertDialogCancel>
                <AlertDialogAction onClick={() => signOut()}>Confirm</AlertDialogAction>
              </AlertDialogFooter>
            </AlertDialogContent>
          </AlertDialog>
        )}
      </div>
    </aside>
  )
}
