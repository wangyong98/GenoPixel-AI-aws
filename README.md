# GenoPixel AI — AWS

Single-cell genomics analysis assistant deployed on AWS Bedrock AgentCore. Users browse and select datasets through a React frontend, then chat with an AI agent that loads and analyses the selected h5ad file using Scanpy.

## Architecture

```text
User browser
  │
  ├─ React frontend (AWS Amplify)
  │     ├─ Datasets tab  →  Catalog API (Lambda + API Gateway)
  │     └─ Chat tab      →  AgentCore Runtime (Docker container)
  │
AWS
  ├─ Cognito              — user authentication
  ├─ Amplify              — frontend hosting + custom domain
  ├─ AgentCore Runtime    — agent container (Strands + Claude Haiku)
  ├─ AgentCore Gateway    — MCP tool routing
  ├─ AgentCore Memory     — short-term conversation history
  ├─ Lambda               — catalog API handler
  ├─ API Gateway          — catalog API endpoint
  ├─ DynamoDB             — active dataset selection per user
  ├─ EFS                  — h5ad file storage (mounted into agent container)
  └─ S3                   — h5ad files (fallback) + metadata Excel
```

## Repository Structure

```text
infra-cdk/                    — CDK infrastructure (deploy this)
  config.yaml                 — deployment configuration
  lib/
    backend-stack.ts          — AgentCore, Lambda, DynamoDB, EFS, S3
    amplify-hosting-stack.ts  — Amplify frontend hosting
  lambdas/
    genopixel-catalog/
      index.py                — catalog API handler (reads metadata Excel from S3)

patterns/
  strands-genopixel-agent/
    basic_agent.py            — agent entrypoint (Strands + BedrockAgentCoreApp)
    tools/
      gp_tools.py             — Scanpy plotting and dataset tools
    skills/                   — system prompt skill fragments loaded at startup
      genopixel-tool-usage/
      genopixel-plot-formatting/
      scanpy-single-cell-analysis/

Docker/genopixel/             — shared Python modules copied into the agent container
  gp_h5ad_loader.py
  gp_models.py
  gp_runtime_state.py
  gp_plot_response_formatter.py

frontend/                     — React chat + dataset browser UI
data/                         — source metadata Excel (upload to S3 before deploying)
  cellxgene_HCA_final_webUI.xlsx
gateway/                      — AgentCore Gateway shared utilities
scripts/                      — deployment helper scripts
tests/                        — integration tests
```

## Prerequisites

- AWS CLI configured for your account
- Node.js 18+ and npm
- AWS CDK CLI: `npm install -g aws-cdk`
- Python 3.11+ (for scripts and tests)
- Docker (for building the agent container image)

## Deployment

All commands run from `infra-cdk/`.

### 1. Configure

Edit `infra-cdk/config.yaml`:

```yaml
stack_name_base: genopixel
admin_user_email: you@example.com
custom_domain: www.yoursite.com   # omit to use raw Amplify URL

backend:
  pattern: strands-genopixel-agent
  model_id: us.anthropic.claude-haiku-4-5-20251001-v1:0
  deployment_type: docker
  network_mode: VPC               # required for EFS access
  vpc:
    vpc_id: <your-vpc-id>
    subnet_ids:
      - <subnet-a>
      - <subnet-b>
```

### 2. Install and bootstrap

```bash
cd infra-cdk
npm install
npx cdk bootstrap   # first time only
```

### 3. Deploy

```bash
npx cdk deploy --all
```

This creates three stacks in order: Cognito → Backend → Amplify.

### 4. Upload metadata Excel

```bash
aws s3 cp data/cellxgene_HCA_final_webUI.xlsx s3://<H5AD_S3_BUCKET>/metadata/metadata.xlsx
```

The S3 bucket name is output by the CDK deploy. The catalog Lambda reads the Excel from this fixed key on every request and re-parses when the ETag changes.

### 5. Upload h5ad files

```bash
aws s3 cp <file>.h5ad s3://<H5AD_S3_BUCKET>/
# or sync a directory
aws s3 sync <local-h5ad-dir>/ s3://<H5AD_S3_BUCKET>/
```

Files are served to the agent via EFS (mounted at `/mnt/genopixel/h5ad` inside the container). S3 is the fallback if a file is not on EFS.

## Dataset Metadata Excel

The catalog is driven by `data/cellxgene_HCA_final_webUI.xlsx` with two sheets:

- **all** — one row per parent dataset. Required columns: `title`, `author`, `file` (h5ad filename), `tissue`, `disease`, `organism`, `project`, `journal`, `cell_counts`, `merged`, `year`, `doi`, `cellxgene_doi`
- **multiple** — variant rows for merged datasets, linked to parents via `cellxgene_doi`. Columns: `publication`, `file`, `cell_counts`, `tissue`, `disease`, `organism`, `description`

To update the catalog, upload a new Excel to S3 — no redeploy needed.

## Agent Skills

Skills are system prompt fragments in `patterns/strands-genopixel-agent/skills/`. They are loaded at agent startup and appended to the system prompt in this order:

1. `genopixel-tool-usage` — tool selection rules and follow-up suggestion logic
2. `genopixel-plot-formatting` — how to render plot responses
3. `scanpy-single-cell-analysis` — biological interpretation hints

Edit the `SKILL.md` files to change agent behaviour, then redeploy the container.

## Updating

### Agent code or skills

```bash
cd infra-cdk
npx cdk deploy BackendStack
```

### Frontend only

```bash
python3 scripts/deploy-frontend.py
```

### Infrastructure only (no container rebuild)

```bash
npx cdk deploy --all --no-rollback
```

## Useful CDK Commands

```bash
npx cdk diff          # preview changes before deploying
npx cdk synth         # emit CloudFormation templates
npx cdk destroy --all # tear down all resources
```

## Troubleshooting

**Catalog returns 503**: the metadata Excel has not been uploaded to S3 yet, or `H5AD_S3_BUCKET` env var is not set on the Lambda.

**Agent says "no active dataset"**: the user has not clicked "Analyze this data" on a dataset in the Datasets tab, or the DynamoDB selection has expired (24-hour TTL).

**h5ad file not found**: the filename in the Excel `file` column does not match any file on EFS or in S3. Upload the file and verify the name matches exactly.

**Cold start slow**: the agent container imports Scanpy and related scientific libraries on startup (~20s). The container stays warm between requests; subsequent invocations are fast.

**CloudFormation deployment errors**: check the Events tab for the failing stack in the AWS Console.
