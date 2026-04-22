// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

import { Routes, Route } from "react-router-dom"
import ChatPage from "./ChatPage"
import DatasetsPage from "./DatasetsPage"
import IntroductionPage from "./IntroductionPage"

export default function AppRoutes() {
  return (
    <Routes>
      <Route path="/" element={<IntroductionPage />} />
      <Route path="/introduction" element={<IntroductionPage />} />
      <Route path="/datasets" element={<DatasetsPage />} />
      <Route path="/chat" element={<ChatPage />} />
    </Routes>
  )
}
