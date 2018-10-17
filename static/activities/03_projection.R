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

# t[2,3], t[1,3], t[1,2] +
b_noise = 1.5

# t[3,3], t[3,2]  -
b_temp = -1.8

# input future noise and temp
noise = 2
temp = 0.5

# setup simulation
pops = 11
ecotype = c(rep("Coastal", 4), rep("Mountain", 3), rep("Palms", 4))

all.means = matrix(c(0.9, 0.3, 0.05,
                   0.1, 0.5, 0.2,
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
      
      # effect of noise
      if(i == 1 && t == 2 | i == 1 && t == 3 | i == 2 && t == 3){
        tp[r,,i,t] = plogis(qlogis(tp[r,,i,t]) + b_noise*noise)
      }
      
      # effect of temp range
      if(i == 3 && t == 2 | i == 3 && t == 3){
        tp[r,,i,t] = plogis(qlogis(tp[r,,i,t]) + b_temp*temp)
      }
      
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

# plot number in each state over time

res %>% 
  group_by(rep, year) %>% 
  count(state) %>% 
  ungroup() %>% 
  group_by(year, state) %>% 
  summarize(mprop = median(n),
            lcl = quantile(n, probs = 0.025),
            ucl = quantile(n, probs = 0.975)) %>% 
  ungroup() %>% 
  complete(year, state, fill = list(mprop = 0)) %>% 
  mutate(state = factor(c("Extirpated", "Low", "High")[state],
                        levels = c("Extirpated", "Low", "High"))) %>% 
  ggplot(aes(x = year, col = state, fill = state)) +
  #geom_linerange(aes(ymin = lcl, ymax = ucl), alpha= 0.9) +
  geom_line(aes(y = lcl), lty = 2, lwd = 1) +
  geom_line(aes(y = ucl), lty = 2, lwd = 1) +
  geom_line(aes(y = mprop), lwd = 1.5) +
  facet_grid(~state, scales= "free") +
  scale_fill_viridis_d(end = 0.7, name = "State") +
  scale_color_viridis_d(end = 0.7, name = "State") +
  scale_x_continuous(breaks = seq(0, nyrs, 2)) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  xlab("Year") +
  ylab("Number of populations in each state (11 total)") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)


# number in final state by ecotype

res %>% 
  group_by(rep, year, ecotype) %>% 
  count(state) %>% 
  ungroup() %>% 
  group_by(year, ecotype, state) %>% 
  summarize(mprop = median(n),
            lcl = quantile(n, probs = 0.025),
            ucl = quantile(n, probs = 0.975)) %>% 
  ungroup() %>% 
  complete(year, ecotype, state, fill = list(mprop = 0)) %>% 
  filter(year == nyrs) %>% 
  mutate(state = factor(c("Extirpated", "Low", "High")[state],
                        levels = c("Extirpated", "Low", "High")),
         ecotype = c("Coastal" = "Coastal (4 total)",
                     "Mountain" = "Mountain (3 total)",
                     "Palms" = "Paradise Palms (4 total)")[ecotype]) %>% 
  ggplot(aes(x = state, y = mprop, col= state, fill = state)) +
  geom_bar(stat = "identity", alpha = 0.5, width = 0.75) +
  geom_linerange(aes(ymin = lcl, ymax = ucl), lwd = 1.5) +
  scale_fill_viridis_d(end = 0.7, name = "State") +
  scale_color_viridis_d(end = 0.7, name = "State") +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  xlab("Population state in final year") +
  ylab("Number of populations") +
  facet_wrap(~ecotype)


# prop at least one high
res %>% 
  filter(year == nyrs) %>% 
  group_by(rep, ecotype) %>% 
  summarize(h = ifelse(3 %in% state,1, 0)) %>% 
  ungroup() %>% 
  group_by(rep) %>% 
  summarize(hs = ifelse(sum(h) == 3, 1, 0)) %>% 
  pull(hs) %>% 
  mean()


# extinction probs
res %>% 
  filter(year == nyrs) %>% 
  group_by(rep, ecotype) %>% 
  summarize(h = ifelse(all(state == 1), 1, 0)) %>%
  ungroup() %>% 
  group_by(ecotype) %>% 
  summarize(pext = mean(h))

res %>% 
  filter(year == nyrs) %>% 
  group_by(rep) %>% 
  summarize(h = ifelse(all(state == 1), 1, 0)) %>%
  ungroup() %>% 
  summarize(pext = mean(h))


#### compare scenarios
all_res <- res1 %>% 
  bind_rows(res2) %>% 
  bind_rows(res3)


all_res %>% 
  group_by(rep, year, scenario) %>% 
  count(state) %>% 
  ungroup() %>% 
  group_by(year, state, scenario) %>% 
  summarize(mprop = median(n),
            lcl = quantile(n, probs = 0.025),
            ucl = quantile(n, probs = 0.975)) %>% 
  ungroup() %>% 
  complete(year, state, scenario, fill = list(mprop = 0)) %>% 
  mutate(state = factor(c("Extirpated", "Low", "High")[state],
                        levels = c("Extirpated", "Low", "High")),
         scenario = factor(c("Status quo", "Noise reduction", "Noise increase")[scenario],
                           levels = c("Status quo", "Noise reduction", "Noise increase"))) %>% 
  ggplot(aes(x = year, col = scenario, fill = scenario)) +
  #geom_linerange(aes(ymin = lcl, ymax = ucl), alpha= 0.9) +
  geom_line(aes(y = lcl), lty = 2, lwd = 1) +
  geom_line(aes(y = ucl), lty = 2, lwd = 1) +
  geom_line(aes(y = mprop), lwd = 1.5) +
  facet_grid(~state, scales= "free") +
  scale_fill_viridis_d(end = 0.8, name = "Scenario", direction = -1) +
  scale_color_viridis_d(end = 0.8, name = "Scenario", direction = -1) +
  scale_x_continuous(breaks = seq(0, nyrs, 2)) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,1,1), size = 14)) +
  xlab("Year") +
  ylab("Number of populations in each state (11 total)") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)



all_res %>% 
  group_by(rep, year, ecotype, scenario) %>% 
  count(state) %>% 
  ungroup() %>% 
  group_by(year, ecotype, state, scenario) %>% 
  summarize(mprop = median(n),
            lcl = quantile(n, probs = 0.025),
            ucl = quantile(n, probs = 0.975)) %>% 
  ungroup() %>% 
  complete(year, ecotype, state, scenario, fill = list(mprop = 0)) %>% 
  filter(year == nyrs) %>% 
  mutate(state = factor(c("Extirpated", "Low", "High")[state],
                        levels = c("Extirpated", "Low", "High")),
         ecotype = c("Coastal" = "Coastal (4 total)",
                     "Mountain" = "Mountain (3 total)",
                     "Palms" = "Paradise Palms (4 total)")[ecotype],
         scenario = factor(c("Status quo", "Noise reduction", "Noise increase")[scenario],
                           levels = c("Status quo", "Noise reduction", "Noise increase"))) %>% 
  ggplot(aes(x = state, y = mprop, col= scenario, fill = scenario)) +
  geom_bar(stat = "identity", alpha = 0.5, width = 0.75, 
           position = position_dodge(width = 0.75)) +
  geom_linerange(aes(ymin = lcl, ymax = ucl), lwd = 1.5,
                 position = position_dodge(width = 0.75)) +
  scale_fill_viridis_d(end = 0.8, name = "Scenario", direction = -1) +
  scale_color_viridis_d(end = 0.8, name = "Scenario", direction = -1) +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "top",
        strip.text = element_text(margin = margin(1,1,1,1), size = 14)) +
  xlab("Population state in final year") +
  ylab("Number of populations") +
  facet_wrap(~ecotype)
