---
layout: default
title: Overview
nav_order: 1
description: ""
permalink: /
last_modified_date: 2021-01-15
---

# ACRONYM

Quality Control Tools for Quantitative Immunofluorescence Microscopy Data.
{: .fs-6 .fw-300 }

ACRONYM is an interactive quality control pipeline for vetting and filtering low-quality, single-cell data from multiplexed whole tissue and tissue microarray images. It comprises multiple modules designed to improve dataset quality by identifying and censoring cells impacted by optical and image preprocessing artifacts. Modules include those for region of interest (ROI) selection, assignment of signal intensity and segmentation area cutoffs, correlation of cells across imaging cycles, and those for single-cell clustering and visualization. The pipeline is instantiated as a configurable [Python](https://www.python.org) Class object with the results of each module being cached to allow for dynamic restarts and iterative quality control strategies.

ACRONYM development is led by [Greg Baker](https://github.com/gjbaker) at [Laboratory of Systems Pharmacology](https://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/about/), Harvard Medical School.

## Funding

This work is supported by:

* NIH grant U54CA225088: Systems Pharmacology of Therapeutic and Adverse Responses to Immune Checkpoint and Small Molecule Drugs
* Ludwig Center at Harvard Medical School and the Ludwig Cancer Research Foundation
