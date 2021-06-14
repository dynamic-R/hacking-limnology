library("deSolve")
library("gifski")
library("gganimate")
library("dplyr")
library("tidyr")

mytheme <- theme(axis.text=element_text(size=rel(1.5)),
                 axis.title=element_text(size=rel(1.5)),
                 legend.text=element_text(size=rel(1.2)),
                 legend.title=element_text(size=rel(1.2)))

lv <- function(t, y, p) {
  dx <-  0.2 * y[1] - 0.2*y[1]*y[2]
  dy <- -0.2 * y[2] + 0.2*y[1]*y[2]
  list(c(dx, dy))
}

out <- ode(y=c(Prey=1, Predator=2), 0:100, lv, NULL)

p <- out %>%
  as.data.frame() %>%
  pivot_longer(cols=-1, names_to = "Population",
               values_to = "Abundance") %>%
  ggplot(aes(x=time, y=Abundance, color=Population)) +
  geom_line(size=1.5) + geom_point(size=4) + ylim(0, 2.5) + mytheme
p

anim <- p + transition_reveal(time)
anim_save("predprey_ts.gif", anim, width=600, height=400)


p <- out %>%
  as.data.frame() %>%
  ggplot(aes(x=Prey, y=Predator)) + geom_path(size=1, color="navy") +
  geom_point(size=5, color="tomato") + ylim(0, 2.5) + mytheme
p

anim <- p + transition_reveal(time)
anim_save("predprey_state.gif", anim, width=400, height=400)
