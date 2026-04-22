import { Button } from "@/components/ui/button"
import { Plus } from "lucide-react"
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

type AppTab = "introduction" | "datasets" | "chat"

type ChatHeaderProps = {
  activeTab: AppTab
  onNewChat?: () => void
  canStartNewChat?: boolean
  showChatActions?: boolean
}

export function ChatHeader({
  activeTab,
  onNewChat,
  canStartNewChat = false,
  showChatActions = false,
}: ChatHeaderProps) {
  const { isAuthenticated, signOut } = useAuth()
  const navigate = useNavigate()

  return (
    <header className="flex items-center justify-between px-3 py-2 border-b w-full bg-background shrink-0">
      <div className="flex items-center gap-1">
        <Button
          variant={activeTab === "introduction" ? "default" : "ghost"}
          size="sm"
          className="rounded-b-none border-b-0 mb-[-1px]"
          onClick={() => navigate("/")}
        >
          Introduction
        </Button>
        <Button
          variant={activeTab === "datasets" ? "default" : "ghost"}
          size="sm"
          className="rounded-b-none border-b-0 mb-[-1px]"
          onClick={() => navigate("/datasets")}
        >
          Datasets
        </Button>
        <Button
          variant={activeTab === "chat" ? "default" : "ghost"}
          size="sm"
          className="rounded-b-none border-b-0 mb-[-1px]"
          onClick={() => navigate("/chat")}
        >
          Chat
        </Button>
      </div>
      {showChatActions && (
        <div className="flex items-center gap-2">
          <Button
            onClick={onNewChat}
            variant="outline"
            className="gap-2"
            disabled={!canStartNewChat}
          >
            <Plus className="h-4 w-4" />
            New Chat
          </Button>
          {isAuthenticated && (
            <AlertDialog>
              <AlertDialogTrigger asChild>
                <Button variant="outline">Logout</Button>
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
      )}
    </header>
  )
}
