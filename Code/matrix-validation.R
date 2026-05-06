library(dplyr)
library(tidyr)
library(ggplot2)
library(haven)
library(janitor)

# ---------------------------------------------------------------------------
# Load LSOG data
# ---------------------------------------------------------------------------

tmp = tempfile(fileext = ".rda")
download.file(
  url  = "https://www.dropbox.com/scl/fi/z8l5gjj0jh5347wa4fkto/22100-0001-Data.rda?rlkey=xdak54wcl4uy0mpdsb20l83y2&st=hnduu3zx&dl=1",
  destfile = tmp,
  mode = "wb")
load(tmp)

lsog = da22100.0001
glimpse(lsog)
