---
layout: default
title: Overview
nav_order: 1
description: ""
permalink: /
last_modified_date: 2021-01-15
---

# Cytolinter

Cytolinter vets and filters low-quality, single-cell data from multiplexed whole tissue and tissue microarray images. It comprises multiple modules designed to improve dataset quality by cross-referencing tabular single-cell data with the microscopy images from which they were derived to identify and censor cells impacted by optical and image preprocessing artifacts. The pipeline is instantiated as a configurable [Python](https://www.python.org) Class object whose partial results are cached after completion of each module to allow for dynamic restarts and iterative QC strategies.

Cytolinter development is led by [Greg Baker](https://github.com/gjbaker) at [Laboratory of Systems Pharmacology](https://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/about/), Harvard Medical School.

## Funding

This work is supported by:

* NIH grant U54CA225088: Systems Pharmacology of Therapeutic and Adverse Responses to Immune Checkpoint and Small Molecule Drugs
* Ludwig Center at Harvard Medical School and the Ludwig Cancer Research Foundation
