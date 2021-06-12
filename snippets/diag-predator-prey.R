grViz("digraph lotka {
         graph [rankdir = 'LR']
           node [shape = box, penwidth=2, fontname = 'Helvetica']
             Predator, Prey
           node [shape = octagon, penwidth=0.5, style='rounded', fixedsize=25, fontsize=8]
             Source Sink
           node [shape = none, fontsize=10]
             growth death grazing
           edge [penwidth=1.5]
             Source -> growth -> Prey -> grazing -> Predator -> death -> Sink
           edge [penwidth=0.7, tailport = 'n', headport = 'n', constraint = false, color=tomato]
             Prey -> growth
             Predator -> grazing
           edge [penwidth=0.7, tailport = 's', headport = 's', constraint = false, color=tomato]
             Prey -> grazing
             Predator -> death

}")
