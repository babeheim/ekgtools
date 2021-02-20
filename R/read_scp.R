
hexvec <- function(x){
  output <- rep("", nchar(x) / 2)
  for (i in 1:length(output)) output[i] <- substr(x, (2 * i) - 1, 2 * i)
  output
}

hex2dec <- function(hexvec) {strtoi(paste(rev(hexvec), collapse = ""), 16L)}
hex2char <- function(hexvec) {rawToChar(as.raw(strtoi(hexvec, 16L)))}

# scp/45fe70c3745f069.scp HR is 61
# scp/2f5fef25b6c2918.scp is 78

read_scp <- function(filepath) {
  system_size <- file.info(filepath)$size
  x <- readBin(filepath, "raw", n = system_size)
  out <- list()
  out$source_filename <- filepath
  # first two bytes are the CRC for the remainder of the file
  crc <- strtoi(paste(rev(x[1:2]), collapse = ""), 16L)
  # next four bytes are the file size
  file_size <- hex2dec(x[3:6])

  if (file_size != system_size) stop("internal file_size incorrect, possibly corrupted file!")

  # section zero
  start <- 6
  section_zero_start <- start
  section_zero_crc <- hex2dec(x[start + 1:2])
  section_zero_number <- hex2dec(x[start + 3:4])
  expect_true(section_zero_number == 0)
  section_zero_length <- hex2dec(x[start + 5:8])
  section_zero_section_version <- hex2dec(x[start + 9])
  section_zero_protocol_version <- hex2dec(x[start + 10])
  section_zero_reserved <- hex2char(x[start + 11:16])

  section_zero_body <- x[start + 17:section_zero_length]
  n_sections <- length(section_zero_body)/10
  sec <- list()
  for (i in 1:n_sections) {
    ob <- list()
    ob$section_number <- hex2dec(section_zero_body[10*(i - 1) + 1:2])
    ob$section_length <- hex2dec(section_zero_body[10*(i - 1) + 3:6])
    ob$section_start <- hex2dec(section_zero_body[10*(i - 1) + 7:10])
    sec[[i]] <- ob
  }
  # 2, 4, 4
  toc <- bind_rows(sec) %>% as.data.frame()

  o <- order(toc$section_start)
  toc <- toc[o,]

  drop <- which(toc$section_start == 0)
  toc <- toc[-drop,]

  expect_true(all(c(7, toc$section_start + toc$section_length) == c(toc$section_start, (length(x) + 1))))

  section_one_start <- toc$section_start[which(toc$section_number == 1)]
  section_eight_start <- toc$section_start[which(toc$section_number == 8)]
  section_133_start <- toc$section_start[which(toc$section_number == 133)]

  # section 1 - patient metadata
  start <- section_one_start - 1
  section_one_crc <- hex2dec(x[start + 1:2])
  section_one_number <- hex2dec(x[start + 3:4])
  section_one_length <- hex2dec(x[start + 5:8])
  expect_true(section_one_length == toc$section_length[which(toc$section_number == 1)])
  section_one_section_version <- hex2dec(x[start + 9])
  section_one_protocol_version <- hex2dec(x[start + 10])
  section_one_reserved <- hex2char(x[start + 11:16])

  section_one_body <- x[start + 17:section_one_length]

  # 1 byte tag, 2 bytes length of value, value bytes
  i <- 1
  tag_start <- 0
  tag_stop <- 3
  tags <- list()
  dec_tags <- c(8, 10, 15, 29, 211, 213, 230)

  while (tag_stop < length(section_one_body)) {
    entry <- list()
    entry$tag <- hex2dec(section_one_body[tag_start + 1])
    entry$length <- hex2dec(section_one_body[tag_start + 2:3])
    if (entry$tag %in% dec_tags) {
      entry$value <- as.character(hex2dec(section_one_body[tag_start + 3 + 1:entry$length]))
    } else if (entry$tag %in% c(4, 6, 7)) {
      data_point <- as.character(hex2dec(section_one_body[tag_start + 3 + 1:2]))
      entry$value <- data_point
    } else if (entry$tag %in% c(5, 25)) {
      year <- as.character(hex2dec(section_one_body[tag_start + 3 + 1:2]))
      if (nchar(year) == 1) year <- paste0("0", year)
      month <- as.character(hex2dec(section_one_body[tag_start + 3 + 3]))
      if (nchar(month) == 1) month <- paste0("0", month)
      day <- as.character(hex2dec(section_one_body[tag_start + 3 + 4]))
      if (nchar(day) == 1) day <- paste0("0", day)
      entry$value <- paste(year, month, day, sep = "-")
    } else if (entry$tag %in% c(26)) {
      hour <- as.character(hex2dec(section_one_body[tag_start + 3 + 1]))
      if (nchar(hour) == 1) hour <- paste0("0", hour)
      min <- as.character(hex2dec(section_one_body[tag_start + 3 + 2]))
      if (nchar(min) == 1) min <- paste0("0", min)
      sec <- as.character(hex2dec(section_one_body[tag_start + 3 + 3]))
      if (nchar(sec) == 1) sec <- paste0("0", sec)
      entry$value <- paste(hour, min, sec, sep = ":")
    } else if (entry$tag %in% c(239)) {
      hr <- as.character(hex2dec(section_one_body[tag_start + 3 + 1]))
      entry$value <- hr
    } else if (entry$tag %in% c(14)) {
      inst <- as.character(hex2dec(section_one_body[tag_start + 3 + 1:2]))
      dep <- as.character(hex2dec(section_one_body[tag_start + 3 + 3:4]))
      device <- as.character(hex2dec(section_one_body[tag_start + 3 + 5:6]))
      type <- as.character(hex2dec(section_one_body[tag_start + 3 + 7]))
      manu <- as.character(hex2dec(section_one_body[tag_start + 3 + 8]))
      model <- hex2char(section_one_body[tag_start + 3 + 9:14])
      proto <- as.character(hex2dec(section_one_body[tag_start + 3 + 15]))
      entry$value <- paste(inst, dep, device, type, manu, model, proto, sep = "-")
    } else {  
      entry$value <- hex2char(section_one_body[tag_start + 3 + 1:entry$length])
    }

    tags[[i]] <- entry
    tag_stop <- tag_start + 3 + entry$length
    tag_start <- tag_stop
    i <- i + 1
  }

  tags <- bind_rows(tags) %>% as.data.frame()

  # weird
  drop <- which(tags$length == 0)
  if (length(drop) > 0) tags <- tags[-drop,] 

  for (i in 1:nrow(tags)) {
    if (tags$tag[i] == 0) out$last_name_1 <- tags$value[i]
    if (tags$tag[i] == 1) out$first_name <- tags$value[i]
    if (tags$tag[i] == 2) out$patient_id <- tags$value[i]
    if (tags$tag[i] == 3) out$last_name_2 <- tags$value[i]
    if (tags$tag[i] == 4) out$age <- round(as.numeric(tags$value[i])/365)
    if (tags$tag[i] == 5) out$date_of_birth <- tags$value[i]
    if (tags$tag[i] == 6) out$height <- tags$value[i]
    if (tags$tag[i] == 7) out$weight <- tags$value[i]
    if (tags$tag[i] == 8) out$sex <- tags$value[i]
    if (tags$tag[i] == 14) out$device_id <- tags$value[i]
    if (tags$tag[i] == 25) out$scan_date <- tags$value[i]
    if (tags$tag[i] == 26) out$scan_time <- tags$value[i]
    if (tags$tag[i] == 214) out$raw_filename <- tags$value[i]
    if (tags$tag[i] == 215) out$adv <- tags$value[i]
    if (tags$tag[i] == 216) out$serial_number <- tags$value[i]
    if (tags$tag[i] == 200) out$community <- tags$value[i]
    if (tags$tag[i] == 226) out$version <- tags$value[i]
    if (tags$tag[i] == 239) out$heart_rate <- tags$value[i]
  }

  # section 133 - qr, prs, and heart rate?

  start <- section_133_start - 1
  section_133_crc <- hex2dec(x[start + 1:2])
  section_133_number <- hex2dec(x[start + 3:4])
  section_133_length <- hex2dec(x[start + 5:8])
  expect_true(section_133_length == toc$section_length[which(toc$section_number == 133)])
  section_133_section_version <- hex2dec(x[start + 9])
  section_133_protocol_version <- hex2dec(x[start + 10])
  section_133_reserved <- hex2char(x[start + 11:16])

  if (section_133_length != 0) {
    section_133_body <- x[start + 17:section_133_length]
    pr_string <- section_133_body[275]
    out$pr <- hex2dec(pr_string)
    qrs_string <- section_133_body[285]
    out$qrs <- hex2dec(qrs_string)
  }

  # section 8 - diagnoses
  start <- section_eight_start - 1
  section_eight_crc <- hex2dec(x[start + 1:2])
  section_eight_number <- hex2dec(x[start + 3:4])
  section_eight_length <- hex2dec(x[start + 5:8])
  expect_true(section_eight_length == toc$section_length[which(toc$section_number == 8)])
  section_eight_section_version <- hex2dec(x[start + 9])
  section_eight_protocol_version <- hex2dec(x[start + 10])
  section_eight_reserved <- hex2char(x[start + 11:16])

  if (section_eight_length != 0) {
    section_eight_body <- x[start + 17:section_eight_length]
    statement_year <- hex2dec(section_eight_body[2:3])
    statement_month <- hex2dec(section_eight_body[4])
    if (nchar(statement_month) == 1) statement_month <- paste0("0", statement_month)
    statement_day <- hex2dec(section_eight_body[5])
    if (nchar(statement_day) == 1) statement_day <- paste0("0", statement_day)
    statement_hour <- hex2dec(section_eight_body[6])
    if (nchar(statement_hour) == 1) statement_hour <- paste0("0", statement_hour)
    statement_min <- hex2dec(section_eight_body[7])
    if (nchar(statement_min) == 1) statement_min <- paste0("0", statement_min)
    statement_sec <- hex2dec(section_eight_body[8])
    if (nchar(statement_sec) == 1) statement_sec <- paste0("0", statement_sec)
    out$statement_timestamp <- paste(
      paste(statement_year, statement_month, statement_day, sep = "-"),
      paste(statement_hour, statement_min, statement_sec, sep = ":")
    )
    n_statements <- hex2dec(section_eight_body[9])

    if (n_statements > 0) {
      statements <- list()
      stat_start <- 9
      stat_stop <- 12
      for (i in 1:n_statements) {
        stat <- list()
        stat$n <- hex2dec(section_eight_body[stat_start + 1])
        stat$length <- hex2dec(section_eight_body[stat_start + 2:3])
        stat$text <- hex2char(section_eight_body[stat_start + 4:(stat$length + 1)])
        statements[[i]] <- stat
        stat_stop <- stat_start + 3 + stat$length
        stat_start <- stat_stop
      }

      statements <- bind_rows(statements) %>% as.data.frame()

      out$clinical_statements <- paste(statements$text, collapse = "; ")
    } 
  }

  return(out)
}
