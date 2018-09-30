# island mouse multistate projection
library(tidyverse)
library(cowplot)
library(statmod)

# functions ----

beta_shape_parms = function(mean, sd){
  avg = mean
  var = sd*sd
  
  a = ifelse(var < avg*(1-avg), avg*(((avg*(1-avg))/var)-1), 100*avg)
  b = ifelse(var < avg*(1-avg), (1-avg)*(((avg*(1-avg))/var)-1), 100*(1-avg))
  
  return(list(a = a, b = b))
}

get_beta_vals = function(n, mean, sd){
  x = beta_shape_parms(mean, sd)
  rbeta(n, x$a, x$b)
}

# inputs ----

reps = 1000
nyrs = 15

env = F

# setup simulation
pops = 11
ecotype = c(rep("Coastal", 4), rep("Mountain", 3), rep("Palms", 4))

all.means = matrix(c(0.6, 0.3, 0.05,
                   0.4, 0.5, 0.2,
                   0, 0.2, 0.75), nrow = 3, ncol = 3, byrow = T)
all.sd = matrix(c(0.05, 0.1, 0.02,
              0.08, 0.05, 0.04,
              0, 0.09, 0.05), nrow = 3, ncol = 3, byrow = T)

# initial states
state = array(NA, dim = c(reps, pops, nyrs))
for(i in 1:reps){
  state[i,,1] = c(3, 3, 2, 3, 2, 2, 1, 2, 3, 2, 2)
}

# replicate-level means and sd
t.m = array(NA, dim = c(reps, 3, 3))
t.sd = array(NA, dim = c(reps, 3, 3))

for(i in 1:3){
  for(t in 1:3){
    t.m[1:reps, i, t] = get_beta_vals(n = reps, 
                                    mean = all.means[i,t], 
                                    sd = all.sd[i,t])
    t.sd[1:reps, i, t] = rinvgauss(n = reps, mean = all.sd[i,t])
  }
}

# draw values for each year
tp = array(NA, dim = c(reps, nyrs, 3, 3))
for(r in 1:reps){
  for(i in 1:3){
    for(t in 1:3){
      tp[r,,i,t] = ifelse(env, get_beta_vals(n = nyrs,
                                      mean = t.m[r,i,t],
                                      sd = t.sd[r,i,t]),
                          t.m[r,i,t])
    }
  }
}

# remove NaNs
tp[is.na(tp)] = 0

# make sure colSums == 1
for(r in 1:reps){
  for(t in 1:nyrs){
    for(i in 1:3){
      tp[r,t,,i] = tp[r,t,,i]/sum(tp[r,t,,i])
    }
  }
}

# projection ----

for(r in 1:reps){
  for(i in 1:pops){
    for(t in 2:nyrs){
      state[r,i,t] = which(rmultinom(1, 1, tp[r,t,,state[r,i,t-1]]) == 1)
    }
  }
}

# make dataframe for plotting ----
res = expand.grid(rep = c(1:reps),
                  pop = c(1:pops),
                  year = c(1:nyrs))

res <- res %>% 
  mutate(ecotype = ecotype[pop],
         state = c(state))

# plot proportion in each state over time
res %>% 
  group_by(rep, year, ecotype) %>% 
  count(state) %>% 
  mutate(prop = n/sum(n)) %>% 
  ungroup() %>% 
  group_by(year, ecotype, state) %>% 
  summarize(mprop = median(prop),
            lcl = quantile(prop, probs = 0.025),
            ucl = quantile(prop, probs = 0.975)) %>% 
  ungroup() %>% 
  complete(year, ecotype, state, fill = list(mprop = 0)) %>% 
  mutate(state = factor(c("Extirpated", "Low", "High")[state],
                        levels = c("Extirpated", "Low", "High"))) %>% 
  ggplot(aes(x = year, col = state, fill = state)) +
  geom_linerange(aes(ymin = lcl, ymax = ucl), alpha= 0.9) +
  geom_line(aes(y = mprop), lwd = 1.5) +
  facet_grid(state~ecotype, scales= "free") +
  ylim(0,1) +
  scale_fill_viridis_d(end = 0.7, name = "State") +
  scale_color_viridis_d(end = 0.7, name = "State") +
  scale_x_continuous(breaks = seq(0, nyrs, 2)) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  xlab("Year") +
  ylab("Proportion of populations in each state") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)



res %>% 
  group_by(year, pop) %>% 
  count(state) %>% 
  mutate(prop = n/sum(n)) %>% 
  ungroup() %>% 
  group_by(year, pop, state) %>% 
  summarize(mprop = median(prop),
            lcl = quantile(prop, probs = 0.025),
            ucl = quantile(prop, probs = 0.975)) %>% 
  ungroup() %>% 
  complete(year, pop, state, fill = list(mprop = 0)) %>% 
  mutate(state = factor(c("Extirpated", "Low", "High")[state],
                        levels = c("Extirpated", "Low", "High"))) %>% 
  ggplot(aes(x = year, col = state, fill = state)) +
  geom_linerange(aes(ymin = lcl, ymax = ucl), alpha= 0.9) +
  geom_line(aes(y = mprop), lwd = 1.5) +
  facet_grid(state~ecotype, scales= "free") +
  ylim(0,1) +
  scale_fill_viridis_d(end = 0.7, name = "State") +
  scale_color_viridis_d(end = 0.7, name = "State") +
  scale_x_continuous(breaks = seq(0, nyrs, 2)) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  xlab("Year") +
  ylab("Proportion of populations in each state") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
