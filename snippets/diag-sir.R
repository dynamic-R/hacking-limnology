library(DiagrammeR)

grViz("digraph feedback {
         graph [rankdir = 'LR']
           node [shape = box, penwidth=2, fontname = 'Helvetica', style='filled', color='dodgerblue']
             S, I, R
           node [shape = none, fontsize=10, style='']
             infection, recovery, mortality
            node [shape = octagon, penwidth=0.5, style='rounded', fixedsize=25, fontsize=8]
             X
           edge [penwidth=1.5]
             S -> infection -> I
             I -> recovery -> R
             I -> mortality
             mortality -> X
}")

