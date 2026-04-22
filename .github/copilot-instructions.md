# AI Agent Instructions – GenoPixel AI AWS

## Purpose and architecture

GenoPixel AI is a single-cell genomics analysis assistant deployed entirely on AWS. Users browse datasets through a React frontend (Amplify), select one, and chat with an AI agent that loads and analyses the h5ad file using Scanpy.

There is no local dev stack. All runtime components are AWS services.

## Key files

- Agent entrypoint: `patterns/strands-genopixel-agent/basic_agent.py`
- Agent tools: `patterns/strands-genopixel-agent/tools/gp_tools.py`
- Agent skills (system prompt fragments): `patterns/strands-genopixel-agent/skills/`
- Shared Python modules (copied into container at build): `Docker/genopixel/`
- Catalog Lambda: `infra-cdk/lambdas/genopixel-catalog/index.py`
- CDK infrastructure: `infra-cdk/lib/`
- Deployment config: `infra-cdk/config.yaml`
- Frontend: `frontend/`
- Metadata Excel: `data/cellxgene_HCA_final_webUI.xlsx`

## AWS services

- **Cognito** — user authentication
- **Amplify** — React frontend hosting
- **AgentCore Runtime** — agent Docker container (Strands + Claude Haiku via Bedrock)
- **AgentCore Gateway** — MCP tool routing into the agent
- **AgentCore Memory** — short-term conversation history per session
- **Lambda + API Gateway** — catalog API (dataset list, active-dataset selection)
- **DynamoDB** — stores each user's active dataset selection (24-hour TTL)
- **EFS** — h5ad files mounted into the agent container at `/mnt/genopixel/h5ad`
- **S3** — h5ad fallback storage + metadata Excel at `metadata/metadata.xlsx`

## Deploy

```bash
cd infra-cdk
npm install
npx cdk deploy --all
```

After deploy, upload the metadata Excel:

```bash
aws s3 cp data/cellxgene_HCA_final_webUI.xlsx s3://<H5AD_S3_BUCKET>/metadata/metadata.xlsx
```

## Agent behaviour

The agent is a Strands `Agent` wrapping Claude Haiku on Bedrock. On each invocation it:

1. Reads the user's active dataset selection from DynamoDB (`_try_preload_active_dataset`)
2. Stores the selection in `RUNTIME_STATE` as a pending selection (even if the h5ad is not yet in memory)
3. Creates the agent with the full system prompt (base + skills + dataset context)
4. Streams the response back via `BedrockAgentCoreApp`

Skills in `patterns/strands-genopixel-agent/skills/` are loaded at module import time and appended to the system prompt. Editing a skill requires a container redeploy.

## Dataset flow

1. User selects a dataset in the Datasets tab → frontend calls `POST /api/catalog/{row}/analyze`
2. Catalog Lambda writes the selection to DynamoDB (key: Cognito `sub`)
3. User sends a chat message → agent reads DynamoDB, pre-loads h5ad from EFS (or S3 fallback)
4. `get_active_dataset_info()` returns the loaded dataset or a pending selection with load instructions
5. Agent calls `load_dataset()` if the h5ad is not yet in memory

## Metadata Excel

Two sheets:

- **all** — parent datasets: `title`, `author`, `file` (h5ad filename), `tissue`, `disease`, `organism`, `project`, `journal`, `cell_counts`, `merged`, `year`, `doi`, `cellxgene_doi`
- **multiple** — variant rows linked to parents via `cellxgene_doi`

The catalog Lambda caches the parsed result in memory and invalidates on S3 ETag change. Uploading a new Excel to S3 is sufficient to update the catalog — no Lambda redeploy needed.

## Troubleshooting

- **Catalog 503**: metadata Excel not uploaded to S3, or `H5AD_S3_BUCKET` env var missing on the Lambda
- **Agent "no active dataset"**: user hasn't clicked "Analyze this data", or DynamoDB TTL expired
- **h5ad not found**: filename in Excel `file` column doesn't match any file on EFS or S3
- **CloudFormation errors**: check the Events tab for the failing stack in the AWS Console
- **Container cold start slow (~20s)**: normal — Scanpy and scientific Python imports are heavy; container stays warm between requests
