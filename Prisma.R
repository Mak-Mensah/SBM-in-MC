
# ============================================================
# PRISMA FLOWCHART FOR META-ANALYSIS
# UPDATED VERSION: 67 STUDIES INCLUDED
# ============================================================

# 1. Install and load required packages
packages <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")

new_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


# 2. Define PRISMA flowchart
dot_code <- '
digraph PRISMA {

  graph [
    rankdir = TB,
    bgcolor = "white",
    nodesep = 0.65,
    ranksep = 0.55,
    splines = polyline,
    pad = 0.05,
    outputorder = edgesfirst
  ]

  node [
    shape = box,
    style = filled,
    fontname = "Arial",
    fontsize = 11,
    color = "#555555",
    penwidth = 0.8,
    fontcolor = "black",
    margin = "0.18,0.12"
  ]

  edge [
    color = "black",
    penwidth = 1.0,
    arrowsize = 0.8
  ]


  # ----------------------------------------------------------
  # NODES
  # ----------------------------------------------------------

  records [
    label = "Records Identified (n = 2000):\\n- Web of Science (150)\\n- Google Scholar (250)\\n- ScienceDirect (1100)\\n- Grey Literature (500)",
    fillcolor = "#CDECE8",
    width = 3.0
  ]

  duplicates [
    label = "Records Removed Before Screening\\n(Duplicates)\\n(n = 300)",
    fillcolor = "#F2C5CB",
    width = 3.0
  ]

  screened [
    label = "Records Screened by Title/Abstract\\n(n = 1700)",
    fillcolor = "#F3E8C9",
    width = 3.0
  ]

  excluded_screening [
    label = "Records Excluded - Not Meeting\\nInclusion Criteria (n = 1250)",
    fillcolor = "#F3E8C9",
    width = 3.0
  ]

  sought [
    label = "Full-Text Articles Sought for Retrieval\\n(n = 450)",
    fillcolor = "#F3E8C9",
    width = 3.4
  ]

  not_retrievable [
    label = "Full-Text Not Retrievable\\n(n = 20)",
    fillcolor = "#F3E8C9",
    width = 3.0
  ]

  assessed [
    label = "Full-Text Articles Assessed\\nfor Eligibility (n = 430)",
    fillcolor = "#F3E8C9",
    width = 3.0
  ]

  excluded_fulltext [
    label = "Full-Text Articles Excluded:\\nNo Extractable Statistical Data (n = 363)",
    fillcolor = "#F3E8C9",
    width = 3.6
  ]

  included [
    label = "Studies Included in Meta-Analysis\\n(n = 67)",
    fillcolor = "#C9E0FA",
    color = "#1857B7",
    penwidth = 1.8,
    width = 3.5
  ]


  # ----------------------------------------------------------
  # EDGES / FLOW
  # ----------------------------------------------------------

  records -> duplicates;
  duplicates -> screened;

  screened -> excluded_screening;
  screened -> sought;

  sought -> not_retrievable;
  sought -> assessed;

  assessed -> excluded_fulltext;
  assessed -> included;


  # ----------------------------------------------------------
  # ALIGNMENT
  # ----------------------------------------------------------

  { rank = same; excluded_screening; sought; }
  { rank = same; not_retrievable; assessed; }
  { rank = same; excluded_fulltext; included; }

  excluded_screening -> sought [style = invis, weight = 10];
  not_retrievable -> assessed [style = invis, weight = 10];
  excluded_fulltext -> included [style = invis, weight = 10];

}


# 3. Generate flowchart
'

prisma_plot <- grViz(dot_code)

# Display flowchart
prisma_plot


# 4. Export SVG (vector format)
svg_output <- export_svg(prisma_plot)

writeLines(
  svg_output,
  "PRISMA_Flowchart_MetaAnalysis_67.svg"
)


# 5. Export high-resolution PNG
rsvg_png(
  charToRaw(svg_output),
  file = "PRISMA_Flowchart_MetaAnalysis_67.png",
  width = 1600,
  height = 1300
)


# 6. Completion message
cat(
  "\nUpdated PRISMA flowchart successfully generated!\n",
  "Final studies included: 67\n",
  "Full-text articles excluded: 363\n",
  "\nFiles saved in: ", getwd(), "\n",
  "- PRISMA_Flowchart_MetaAnalysis_67.png\n",
  "- PRISMA_Flowchart_MetaAnalysis_67.svg\n"
)
