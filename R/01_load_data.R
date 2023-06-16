combined_data <- readRDS(gzcon(url("https://github.com/DidDrog11/SL_lassa_ELISA/raw/main/output/ELISA_output.rds")))

trap_data <- readRDS(gzcon(url("https://github.com/DidDrog11/chapter_3/raw/observed_data/data/data_for_export/chapter_4_extract.rds")))

write_rds(combined_data, here("input", "ELISA_data.rds"))
write_rds(trap_data, here("input", "trap_data.rds"))
