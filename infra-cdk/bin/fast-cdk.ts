#!/usr/bin/env node
import * as cdk from "aws-cdk-lib"
import { FastMainStack } from "../lib/fast-main-stack"
import { ConfigManager } from "../lib/utils/config-manager"

// Load configuration using ConfigManager
const configManager = new ConfigManager("config.yaml")

// Initial props consist of configuration parameters
const props = configManager.getProps()

const app = new cdk.App()

// Deploy the new Amplify-based stack that solves the circular dependency
// GenoPixel targets us-east-1 (best AgentCore/Bedrock availability).
// Use AWS_DEFAULT_REGION env var to override if needed.
const amplifyStack = new FastMainStack(app, props.stack_name_base, {
  config: props,
  env: {
    account: process.env.CDK_DEFAULT_ACCOUNT,
    // Always deploy to us-east-1 (best AgentCore availability).
    // Override by setting AWS_DEFAULT_REGION in the shell before running cdk commands.
    region: process.env.AWS_DEFAULT_REGION ?? "us-east-1",
  },
})

app.synth()
