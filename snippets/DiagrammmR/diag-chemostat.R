library("DiagrammeR")
library("DiagrammeRsvg")
library("rsvg")
library("magrittr")

diag <- grViz("digraph chemostat {
         graph [rankdir = 'LR']
           node [shape = box, penwidth=2, fontname = 'Helvetica']
             X P
           node [shape = octagon, penwidth=0.5, style='rounded', fixedsize=25, fontsize=8]
             Source
           node [shape = none, fontsize=10]
             import P_export X_export growth
           node [fontsize=20]
             D
           edge [penwidth=1.5]
             P -> growth
             growth -> X
             X -> X_export
             Source -> import -> P
           edge [tailport='n']
             P -> P_export
           edge [penwidth=0.7, tailport = 's', headport = 's', constraint = false, color=tomato]
             X -> growth
           edge [penwidth=0.7, tailport = 'e', headport = 'n', constraint = true, color=blue]
             D -> import
             D -> P_export
             D -> X_export
}")

diag

diag %>% export_svg %>% charToRaw %>% rsvg_svg("chemostat.svg")

