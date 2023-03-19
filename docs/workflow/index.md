---
layout: default-cylinter
title: Workflow
nav_order: 3
has_children: true
---

# Workflow

| Input Directory Structure <br /> (see [Input Files](input#input-directory-structure) for details) | Pipeline Diagram | Output Directory Structure <br /> (see [Output Files](output#output-directory-structure) for details) |
| :-- | :-: | :-- |
| <code>INPUT_DIR<br>├── config.yml<br>├── csv/<br>├── markers.csv<br>├── mask/<br>├── seg/<br>└── tif/<br></code> | <img src="{{ site.baseurl }}/assets/images/1-overview.png" alt="CyLinter" width="850"/> | <code>OUTPUT_DIR<br>├── area/<br>├── checkpoints/<br>├── clustering/<br>├── contrast/<br>├── cycles/<br>├── intensity/<br>├── metaQC/<br>├── PCA/<br>├── pruning/<br>└──  ROIs/<br></code>
