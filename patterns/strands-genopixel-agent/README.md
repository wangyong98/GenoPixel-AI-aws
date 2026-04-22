# Strands GenoPixel Agent Pattern

This pattern uses the [Strands Agents](https://github.com/strands-agents/strands-agents) framework to build a single agent with Gateway tool access, Code Interpreter, and AgentCore Memory for conversation history.

## Features

- **Token-Level Streaming**: True token-by-token streaming via `agent.stream_async()`
- **AgentCore Memory**: Conversation history persisted across requests via `AgentCoreMemorySessionManager`
- **Code Interpreter**: Secure Python execution via `StrandsCodeInterpreterTools`
- **Gateway Integration**: Access Lambda-based tools through AgentCore Gateway (MCP protocol with OAuth2 auth)
- **Secure Identity**: User identity extracted from validated JWT token (`RequestContext`), not from payload

## Architecture

```
User Request
    |
BedrockAgentCoreApp (basic_agent.py)
    |
Strands Agent (BedrockModel)
    |
    +-- AgentCore Memory (conversation history)
    |     AgentCoreMemorySessionManager
    |
    +-- Code Interpreter
    |     StrandsCodeInterpreterTools (execute_python_securely)
    |
    +-- Gateway MCP Client (streamable HTTP)
          Lambda-based tools via AgentCore Gateway
```

## File Structure

```
patterns/strands-genopixel-agent/
├── basic_agent.py                # Main entrypoint (BedrockAgentCoreApp)
├── tools/                        # GenoPixel and gateway tool integrations
│   ├── gp_tools.py
│   ├── code_interpreter.py
│   └── gateway.py
├── requirements.txt              # Pinned dependencies
└── Dockerfile                    # Container build (Python 3.13)
```

## Available Tools

| Tool | Source | Description |
|------|--------|-------------|
| `execute_python_securely` | Code Interpreter | Execute Python code in a secure sandbox |
| Gateway tools | AgentCore Gateway | Lambda-based tools discovered via MCP |

## Model

- **Configured in**: `infra-cdk/config.yaml` → `backend.model_id`
- **Runtime env var**: `BEDROCK_MODEL_ID`
- **Current default/pin**: `us.anthropic.claude-haiku-4-5-20251001-v1:0`

## Dataset Preload Safeguards

To keep chat responsive when EFS is unavailable and h5ad files fall back to S3:

- `ALLOW_S3_PRELOAD` (default: `false`)
  - `false`: auto-preload on chat invocation only happens for EFS paths.
  - `true`: allow auto-preload from S3 fallback paths.
- `PRELOAD_MAX_BYTES` (default: `1000000000`)
  - Skip auto-preload when file size exceeds this limit.
- `S3_PRELOAD_BACKED` (default: `true`)
  - When S3 auto-preload is enabled, load in AnnData backed mode.
- `S3_FALLBACK_BACKED` (default: `true`)
  - For explicit `load_dataset` calls on S3 fallback files, use backed mode.

## Streaming Events

The agent yields SSE `data: {json}` lines via `agent.stream_async()`. The frontend parser at `frontend/src/lib/agentcore-client/parsers/strands.ts` handles these event types:

| Event | Format | Description |
|-------|--------|-------------|
| Text | `{"data": "text"}` | Token-level text content |
| Tool use start | `{"current_tool_use": {...}, "delta": {"toolUse": {"input": ""}}}` | Tool invocation begins |
| Tool use delta | `{"current_tool_use": {...}, "delta": {"toolUse": {"input": "..."}}}` | Streaming tool input |
| Tool result | `{"message": {"role": "user", "content": [{"toolResult": {...}}]}}` | Tool execution result |
| Result | `{"result": {"stop_reason": "end_turn"}}` | Agent finished |
| Lifecycle | `{"init_event_loop": true}` / `{"start_event_loop": true}` | Agent lifecycle events |

## Memory Integration

This pattern uses **AgentCore Memory** for conversation persistence:

1. `MEMORY_ID` environment variable provides the memory resource ID
2. `AgentCoreMemoryConfig` is initialized with `memory_id`, `session_id`, and `actor_id` (user ID)
3. `AgentCoreMemorySessionManager` handles storing/retrieving conversation history
4. Memory is tied to the `runtimeSessionId` from the client

## Security

- **User identity**: Extracted from the validated JWT token via `RequestContext`, not from the payload body
- **STACK_NAME validation**: Validated for alphanumeric format before use in SSM parameter paths
- **Payload validation**: Required fields (`prompt`, `runtimeSessionId`) validated before processing
- **Gateway auth**: OAuth2 client credentials flow via Cognito for machine-to-machine authentication

## Deployment

```bash
cd infra-cdk
# Set pattern in config.yaml:
#   backend:
#     pattern: strands-genopixel-agent
#     model_id: us.anthropic.claude-haiku-4-5-20251001-v1:0
#     deployment_type: docker
cdk deploy
```

## Dependencies

```
strands-agents==1.24.0
mcp==1.26.0
bedrock-agentcore[strands-agents]==1.2.0
PyJWT[crypto]>=2.10.1
```
