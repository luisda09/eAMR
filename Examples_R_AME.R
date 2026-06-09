
# Install from CRAN
install.packages("AMR")

# Then load it
library(AMR)

# Now load the example data
data("example_isolates")
head(example_isolates)


# Load example soil data from phyloseq
library(phyloseq)
data("soilrep")
soilrep  # 24 soil communities, 16,825 OTUs [web:63]

# ---- 2. LOAD DATA ----------------------------------------------------------
data("example_isolates")
df <- example_isolates

# ---- 3. GENERAL DATA EXPLORATION -------------------------------------------
# View structure
str(df)
head(df)

# Summary statistics
summary(df)

# Count of isolates
nrow(df)  # total isolates
ncol(df)  # number of variables

# Microorganism distribution
df %>%
  count(mo, sort = TRUE) %>%
  mutate(mo_name = mo_name(mo)) %>%
  head(10)

# ---- 4. RESISTANCE PATTERNS BY ANTIBIOTIC ----------------------------------
# Get all antibiotic columns (SIR class)
ab_cols <- df %>% select(where(is.sir)) %>% names()

# Calculate resistance % for each antibiotic
resistance_by_ab <- df %>%
  summarise(across(all_of(ab_cols), ~ resistance(.x) * 100)) %>%
  pivot_longer(everything(), names_to = "antibiotic", values_to = "resistance_pct") %>%
  mutate(antibiotic_name = ab_name(antibiotic)) %>%
  arrange(desc(resistance_pct))

# View top 10 most resistant antibiotics
resistance_by_ab %>% head(10)

# ---- 5. VISUALIZATION: RESISTANCE BY ANTIBIOTIC ----------------------------
resistance_by_ab %>%
  head(15) %>%
  mutate(antibiotic_name = factor(antibiotic_name, levels = rev(antibiotic_name))) %>%
  ggplot(aes(x = antibiotic_name, y = resistance_pct, fill = resistance_pct)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Resistance Percentage by Antibiotic",
    x = "Antibiotic",
    y = "Resistance (%)",
    fill = "Resistance %"
  ) +
  scale_fill_gradient(low = "green", high = "red") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

# ---- 6. MULTI-DRUG RESISTANCE (MDR) ANALYSIS -------------------------------
# Calculate MDR status
df <- df %>%
  mutate(
    mdr_count = rowSums(select(., all_of(ab_cols)) == "R", na.rm = TRUE),
    is_mdr = mdr_count >= 3  # MDR = resistant to ≥3 antibiotic classes
  )

# MDR prevalence
mdr_summary <- df %>%
  summarise(
    total_isolates = n(),
    mdr_isolates = sum(is_mdr),
    mdr_percentage = mean(is_mdr) * 100
  )

mdr_summary

# MDR by microorganism
df %>%
  group_by(mo) %>%
  summarise(
    total = n(),
    mdr = sum(is_mdr),
    mdr_pct = mean(is_mdr) * 100
  ) %>%
  arrange(desc(mdr_pct)) %>%
  top_n(10, mdr_pct) %>%
  mutate(mo_name = mo_name(mo))

# ---- 7. RESISTANCE BY WARD / LOCATION --------------------------------------
resistance_by_ward <- df %>%
  group_by(ward) %>%
  summarise(
    total = n(),
    resistance_pct = mean(select(., all_of(ab_cols)) == "R", na.rm = TRUE) * 100
  )

resistance_by_ward

# ---- 8. RESISTANCE BY MICROORGANISM ----------------------------------------
resistance_by_mo <- df %>%
  group_by(mo) %>%
  summarise(
    total = n(),
    resistance_pct = mean(select(., all_of(ab_cols)) == "R", na.rm = TRUE) * 100
  ) %>%
  arrange(desc(resistance_pct)) %>%
  mutate(mo_name = mo_name(mo))

resistance_by_mo %>% head(10)

# ---- 9. SUSCEPTIBILITY BY ANTIBIOTIC ---------------------------------------
susceptibility_by_ab <- df %>%
  summarise(across(all_of(ab_cols), ~ susceptibility(.x) * 100)) %>%
  pivot_longer(everything(), names_to = "antibiotic", values_to = "susceptibility_pct") %>%
  mutate(antibiotic_name = ab_name(antibiotic)) %>%
  arrange(desc(susceptibility_pct))

# ---- 10. CO-RESISTANCE PATTERNS --------------------------------------------
# Example: resistance to PEN vs RIF
corresistance <- df %>%
  mutate(
    PEN_R = `PEN:RIF` == "R",
    RIF_R = RIF == "R",
    both_R = PEN_R & RIF_R
  ) %>%
  summarise(
    PEN_R_only = sum(PEN_R & !RIF_R),
    RIF_R_only = sum(RIF_R & !PEN_R),
    both_resistant = sum(both_R),
    neither = sum(!PEN_R & !RIF_R)
  )

corresistance

# ---- 11. TRENDS OVER TIME (if date available) -------------------------------
if ("date" %in% names(df)) {
  df %>%
    mutate(year = year(date), month = month(date)) %>%
    group_by(year) %>%
    summarise(
      total = n(),
      resistance_pct = mean(select(., all_of(ab_cols)) == "R", na.rm = TRUE) * 100
    ) %>%
    ggplot(aes(x = year, y = resistance_pct)) +
    geom_line() +
    geom_point() +
    labs(
      title = "Resistance Trend Over Time",
      x = "Year",
      y = "Resistance (%)"
    ) +
    theme_bw()
}


########      to conect this script to github   #############
#2.5 Connect RStudio and GitHub
#GitHub needs to validate who you are before you can connect it and RStudio. We can do this by 
#generating and saving a Personal Access Token (PAT). You need to do this once for every RStudio project.

usethis::create_github_token()

#token generate: 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

#the run:

gitcreds::gitcreds_set()


# another extra test for more

##this is another proof1gh

#a new change more   

## The other to see the branch issue


