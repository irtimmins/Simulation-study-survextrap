###################################

store_directory <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/"
store_dir_cetux <- paste0(store_directory , "simsurvextrap_", "cetux_sigma/")
store_dir_niv <- paste0(store_directory , "simsurvextrap_", "niv_sigma/")

setwd(store_dir_cetux)

test <- readRDS("scen5/model1.rds")
test2 <- readRDS("scen6/model1.rds")

test_draws <- readRDS("scen5_sigma_draws.rds")
test2_draws <- readRDS("scen6_sigma_draws.rds")


test_draws
test2_draws
test$km %>% head()
test2$km %>% head()

test$mspline
test2$mspline
plot(test)
plot(test2)


test$fit_method
test2$fit_method

test
test2

sigma_post_all <- readRDS("sigma_post_density_all.rds")
sigma_post_all <- sigma_post_all %>%
  filter(scenario_id %in% c(5,6,7))
sigma_post_all %>%
  filter(scenario_id == 5)

sigma_post_all %>%
  filter(scenario_id == 6)

sigma_post_all %>%
  filter(scenario_id == 7)

test_draws <- as_draws_df(test$stanfit) %>%
  as_tibble() %>%
  select(starts_with("hsd")) %>%
  rename(sigma = `hsd[1]`)

test2_draws <- as_draws_df(test2$stanfit) %>%
  as_tibble() %>%
  select(starts_with("hsd")) %>%
  rename(sigma = `hsd[1]`)

test_draws
test2_draws
scenarios <- read_csv("scenarios_cetux_OS.csv")
scenarios %>% slice(5:6) %>% View()
post_dens1 <- density(test_draws$sigma, bw = 0.15)
post_dens2 <- density(test2_draws$sigma, bw = 0.15)
post_df1 <- tibble(x = post_dens1$x,
                  dist = post_dens1$y) %>%
  mutate(test = "test1")
post_df2 <- tibble(x = post_dens2$x,
                   dist = post_dens2$y) %>%
  mutate(test = "test2")

rbind(post_df1, post_df2) %>%
  ggplot(aes(x = x, y = dist, colour = test))+
  geom_line()





