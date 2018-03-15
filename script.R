library(tidyverse)
subject <- rep(paste0("Subj", 1:6), 3)
trt <- rep(c("A", "B"), 9)
week <- rep(1:3, each = 6)
value <- sample.int(10, size = 18, replace = TRUE)
df <- data_frame(subject, trt, week, value)

ggplot(df, aes(x = week, y = subject, fill = value)) +
   facet_grid(trt ~ ., scale = "free") +
   geom_tile()
