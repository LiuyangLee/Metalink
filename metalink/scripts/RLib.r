###
#' @Date: 2022-07-15 23:40:20
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2024-06-03 15:11:20
#' @FilePath: /liliuyang/R/RLib/R/RLib.r
#' @Description:

###
# Minimal helper library bundled with Metalink
# Only functions used by the Metalink .cmd.R scripts are kept here.
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

split_lineage_fill <- function(x, n = 9, fill = "") {
  res <- tstrsplit(x, ";", fixed = TRUE, fill = fill)
  if (length(res) < n) {
    res <- c(res, rep(list(fill), n - length(res)))
  }
  res <- res[1:n]
  lapply(res, function(v) replace(v, is.na(v) | is.null(v), fill))
}


gtdb_tax = function(Ref_metadata){
  Ref_metadata$Domain = sapply(as.character(Ref_metadata$Taxonomy),
                               function(x) unlist(strsplit(x, ";"))[1])
  Ref_metadata$Phylum = sapply(as.character(Ref_metadata$Taxonomy),
                               function(x) unlist(strsplit(x, ";"))[2])
  Ref_metadata$Class = sapply(as.character(Ref_metadata$Taxonomy),
                              function(x) unlist(strsplit(x, ";"))[3])
  Ref_metadata$Order = sapply(as.character(Ref_metadata$Taxonomy),
                              function(x) unlist(strsplit(x, ";"))[4])
  Ref_metadata$Family = sapply(as.character(Ref_metadata$Taxonomy),
                               function(x) unlist(strsplit(x, ";"))[5])
  Ref_metadata$Genus = sapply(as.character(Ref_metadata$Taxonomy),
                              function(x) unlist(strsplit(x, ";"))[6])
  Ref_metadata$Species = sapply(as.character(Ref_metadata$Taxonomy),
                                function(x) unlist(strsplit(x, ";"))[7])
  return(Ref_metadata)
}


rename_GTDB_taxonomy_dataframe = function(data) {
  # 0-row input (e.g. a vMAG/mMAG absent from this sample's coverm output):
  # nothing to rename. Keep the expected schema (NewTax column) so downstream
  # select()/mutate() calls still work; case_when() on empty columns would
  # otherwise fail with a character/logical type conflict.
  if (nrow(data) == 0) {
    if (!"NewTax" %in% colnames(data)) data$NewTax <- character(0)
    return(data)
  }
  data[is.na(data)] <- ""
  prescale_data = data %>%
    mutate(., Species = (ifelse((is.na(.$Species)| (grepl("^s__.*unclassified", .$Species)& !grepl("s__root species", .$Species))  |.$Species == "no support" | .$Species == "" ), "s__",.$Species))) %>%
    mutate(., Genus = (ifelse((is.na(.$Genus)| (grepl("^g__.*unclassified", .$Genus)& !grepl("g__root genus", .$Genus))  |.$Genus == "no support" | .$Genus == "" ), "g__", .$Genus))) %>%
    mutate(., Family = (ifelse((is.na(.$Family)| (grepl("^f__.*unclassified", .$Family)& !grepl("f__root family", .$Family))  | .$Family == "no support" | .$Family == "" ), "f__", .$Family))) %>%
    mutate(., Order = (ifelse((is.na(.$Order)| (grepl("^o__.*unclassified", .$Order)& !grepl("o__root order", .$Order))  | .$Order == "no support" | .$Order == "" ), "o__", .$Order))) %>%
    mutate(., Class = (ifelse((is.na(.$Class)| (grepl("^c__.*unclassified", .$Class)& !grepl("c__root class", .$Class))  |.$Class == "no support" | .$Class == "" ), "c__", .$Class))) %>%
    mutate(., Phylum = (ifelse((is.na(.$Phylum)| (grepl("^p__.*unclassified", .$Phylum)& !grepl("p__root phylum", .$Phylum))  |.$Phylum == "no support" | .$Phylum == "" ), "p__", .$Phylum))) %>%
    mutate(., Domain = (ifelse((is.na(.$Domain)| (grepl(".*unclassified", .$Domain) & !grepl("d__root domain", .$Domain)) |.$Domain == "no support" | .$Domain == "" ), "d__unclassified", .$Domain)))
  scale_data = prescale_data %>%
    mutate(NewTax = ifelse(Species == "s__", "New species", "Known species")) %>%
    mutate(NewTax = ifelse(Genus == "g__" & Species == "s__", "New genus", NewTax)) %>%
    mutate(NewTax = ifelse(Family == "f__" & Genus == "g__" & Species == "s__" ,"New family", NewTax)) %>%
    mutate(NewTax = ifelse(Order == "o__" & Family == "f__" & Genus == "g__" & Species == "s__", "New order", NewTax)) %>%
    mutate(NewTax = ifelse(Class == "c__" & Order == "o__" & Family == "f__" & Genus == "g__" & Species == "s__", "New class", NewTax)) %>%
    mutate(NewTax = ifelse(Phylum == "p__" & Class == "c__" & Order == "o__" & Family == "f__" & Genus == "g__" & Species == "s__", "New phylum", NewTax)) %>%
    mutate(NewTax = case_when(
      Domain == "d__root domain" ~ "Root",
      Domain == "d__unclassified"        ~ "Unclassified",
      Domain == "Unclassified"           ~ "Unclassified",
      Domain == "Unclassified_Archaea"   ~ "Unclassified_Archaea",
      Domain == "Unclassified_Bacteria"  ~ "Unclassified_Bacteria",
      Domain == "Unclassified_Eukaryota" ~ "Unclassified_Eukaryota",
      Domain == "Unclassified_Virus"     ~ "Unclassified_Virus",
      TRUE                              ~ NewTax
    )) %>%
    dplyr::mutate(
      Phylum = dplyr::case_when(
        Phylum == "p__" & NewTax == "Unclassified" ~ paste0("p__", Domain, " phylum"),
        Phylum == "p__" & NewTax == "Unclassified_Archaea" ~ paste0("p__", Domain, " phylum"),
        Phylum == "p__" & NewTax == "Unclassified_Bacteria" ~ paste0("p__", Domain, " phylum"),
        Phylum == "p__" & NewTax == "Unclassified_Eukaryota" ~ paste0("p__", Domain, " phylum"),
        Phylum == "p__" & NewTax == "Unclassified_Virus" ~ paste0("p__", Domain, " phylum"),
        Phylum == "p__" & NewTax == "New phylum" ~ paste0("p__unclassified ", Domain, " phylum"),
        TRUE                                     ~ Phylum
      ),
      Class = dplyr::case_when(
        Class == "c__" & NewTax == "Unclassified" ~ paste0("c__", Domain, " class"),
        Class == "c__" & NewTax == "Unclassified_Archaea" ~ paste0("c__", Domain, " class"),
        Class == "c__" & NewTax == "Unclassified_Bacteria" ~ paste0("c__", Domain, " class"),
        Class == "c__" & NewTax == "Unclassified_Eukaryota" ~ paste0("c__", Domain, " class"),
        Class == "c__" & NewTax == "Unclassified_Virus" ~ paste0("c__", Domain, " class"),
        Class == "c__" & NewTax == "New phylum"  ~ paste0("c__unclassified ", Domain, " class"),
        Class == "c__" & NewTax == "New class"   ~ paste0("c__unclassified ", Phylum, " class"),
        TRUE                                     ~ Class
      ),
      Order = dplyr::case_when(
        Order == "o__" & NewTax == "Unclassified" ~ paste0("o__", Domain, " order"),
        Order == "o__" & NewTax == "Unclassified_Archaea" ~ paste0("o__", Domain, " order"),
        Order == "o__" & NewTax == "Unclassified_Bacteria" ~ paste0("o__", Domain, " order"),
        Order == "o__" & NewTax == "Unclassified_Eukaryota" ~ paste0("o__", Domain, " order"),
        Order == "o__" & NewTax == "Unclassified_Virus" ~ paste0("o__", Domain, " order"),
        Order == "o__" & NewTax == "New phylum"  ~ paste0("o__unclassified ", Domain, " order"),
        Order == "o__" & NewTax == "New class"   ~ paste0("o__unclassified ", Phylum, " order"),
        Order == "o__" & NewTax == "New order"   ~ paste0("o__unclassified ", Class, " order"),
        TRUE                                     ~ Order
      ),
      Family = dplyr::case_when(
        Family == "f__" & NewTax == "Unclassified" ~ paste0("f__", Domain, " family"),
        Family == "f__" & NewTax == "Unclassified_Archaea" ~ paste0("f__", Domain, " family"),
        Family == "f__" & NewTax == "Unclassified_Bacteria" ~ paste0("f__", Domain, " family"),
        Family == "f__" & NewTax == "Unclassified_Eukaryota" ~ paste0("f__", Domain, " family"),
        Family == "f__" & NewTax == "Unclassified_Virus" ~ paste0("f__", Domain, " family"),
        Family == "f__" & NewTax == "New phylum" ~ paste0("f__unclassified ", Domain, " family"),
        Family == "f__" & NewTax == "New class"  ~ paste0("f__unclassified ", Phylum, " family"),
        Family == "f__" & NewTax == "New order"  ~ paste0("f__unclassified ", Class, " family"),
        Family == "f__" & NewTax == "New family" ~ paste0("f__unclassified ", Order, " family"),
        TRUE                                     ~ Family
      ),
      Genus = dplyr::case_when(
        Genus == "g__" & NewTax == "Unclassified" ~ paste0("g__", Domain, " genus"),
        Genus == "g__" & NewTax == "Unclassified_Archaea" ~ paste0("g__", Domain, " genus"),
        Genus == "g__" & NewTax == "Unclassified_Bacteria" ~ paste0("g__", Domain, " genus"),
        Genus == "g__" & NewTax == "Unclassified_Eukaryota" ~ paste0("g__", Domain, " genus"),
        Genus == "g__" & NewTax == "Unclassified_Virus" ~ paste0("g__", Domain, " genus"),
        Genus == "g__" & NewTax == "New phylum"  ~ paste0("g__unclassified ", Domain, " genus"),
        Genus == "g__" & NewTax == "New class"   ~ paste0("g__unclassified ", Phylum, " genus"),
        Genus == "g__" & NewTax == "New order"   ~ paste0("g__unclassified ", Class, " genus"),
        Genus == "g__" & NewTax == "New family"  ~ paste0("g__unclassified ", Order, " genus"),
        Genus == "g__" & NewTax == "New genus"  ~ paste0("g__unclassified ", Family, " genus"),
        TRUE                                     ~ Genus
      ),
      Species = dplyr::case_when(
        Species == "s__" & NewTax == "Unclassified" ~ paste0("s__", Domain, " species"),
        Species == "s__" & NewTax == "Unclassified_Archaea" ~ paste0("s__", Domain, " species"),
        Species == "s__" & NewTax == "Unclassified_Bacteria" ~ paste0("s__", Domain, " species"),
        Species == "s__" & NewTax == "Unclassified_Eukaryota" ~ paste0("s__", Domain, " species"),
        Species == "s__" & NewTax == "Unclassified_Virus" ~ paste0("s__", Domain, " species"),
        Species == "s__" & NewTax == "New phylum"  ~ paste0("s__unclassified ", Domain, " species"),
        Species == "s__" & NewTax == "New class"   ~ paste0("s__unclassified ", Phylum, " species"),
        Species == "s__" & NewTax == "New order"   ~ paste0("s__unclassified ", Class, " species"),
        Species == "s__" & NewTax == "New family"  ~ paste0("s__unclassified ", Order, " species"),
        Species == "s__" & NewTax == "New genus"   ~ paste0("s__unclassified ", Family, " species"),
        Species == "s__" & NewTax == "New species" ~ paste0("s__unclassified ", Genus, " species"),
        TRUE                                       ~ Species
      )
  ) %>%
    # 填充病毒和真核中间层级为空的
    dplyr::mutate(
      ## 2. Phylum
      Phylum = dplyr::case_when(
        Phylum != "p__" ~ Phylum,
        Phylum == "p__" & NewTax %in% c("New class","New order","New family","New genus","New species","Known species") &
          !grepl("^[a-z]__.*unclassified", Domain)                          ~ paste0("p__unclassified ", Domain, " phylum"),
        Phylum == "p__" & NewTax %in% c("New class","New order","New family","New genus",
                                        "New species","Known species") &
          grepl("^[a-z]__.*unclassified", Domain)                           ~ paste0("p__unclassified ", Domain, " phylum"),
        TRUE                                                                ~ Phylum
      )
    ) %>%
    dplyr::mutate(
      ## 1. Class
      Class = dplyr::case_when(
        Class != "c__" ~ Class,
        Class == "c__" & NewTax %in% c("New order","New family","New genus","New species","Known species") &
          !grepl("^[a-z]__.*unclassified", Phylum)                     ~ paste0("c__unclassified ", Phylum, " class"),
        Class == "c__" & NewTax %in% c("New order","New family","New genus","New species","Known species") &
          grepl("^[a-z]__.*unclassified", Phylum) &
          !grepl("^[a-z]__.*unclassified", Domain)                     ~ paste0("c__unclassified ", Domain, " class"),
        Class == "c__" & NewTax %in% c("New order","New family","New genus","New species","Known species") &
          grepl("^[a-z]__.*unclassified", Phylum) &
          grepl("^[a-z]__.*unclassified", Domain)                      ~ paste0("c__unclassified ", Domain, " class"),
        TRUE                                                          ~ Class
      )
    ) %>%
    dplyr::mutate(
      ## 2. Order
      Order = dplyr::case_when(
        Order != "o__" ~ Order,
        Order == "o__" & NewTax %in% c("New family","New genus","New species","Known species") &
          !grepl("^[a-z]__.*unclassified", Class)                      ~ paste0("o__unclassified ", Class, " order"),
        Order == "o__" & NewTax %in% c("New family","New genus","New species","Known species") &
          grepl("^[a-z]__.*unclassified", Class) &
          !grepl("^[a-z]__.*unclassified", Phylum)                     ~ paste0("o__unclassified ", Phylum, " order"),
        Order == "o__" & NewTax %in% c("New family","New genus","New species","Known species") &
          grepl("^[a-z]__.*unclassified", Class) &
          grepl("^[a-z]__.*unclassified", Phylum) &
          !grepl("^[a-z]__.*unclassified", Domain)                     ~ paste0("o__unclassified ", Domain, " order"),
        Order == "o__" & NewTax %in% c("New family","New genus","New species","Known species") &
          grepl("^[a-z]__.*unclassified", Class) &
          grepl("^[a-z]__.*unclassified", Phylum) &
          grepl("^[a-z]__.*unclassified", Domain)                      ~ paste0("o__unclassified ", Domain, " order"),
        TRUE                                                          ~ Order
      )
    ) %>%
    dplyr::mutate(
      ## 3. Family
      Family = dplyr::case_when(
        Family != "f__" ~ Family,
        Family == "f__" & NewTax %in% c("New genus","New species","Known species") &
          !grepl("^[a-z]__.*unclassified", Order)                      ~ paste0("f__unclassified ", Order, " family"),
        Family == "f__" & NewTax %in% c("New genus","New species","Known species") &
          grepl("^[a-z]__.*unclassified", Order) &
          !grepl("^[a-z]__.*unclassified", Class)                      ~ paste0("f__unclassified ", Class, " family"),
        Family == "f__" & NewTax %in% c("New genus","New species","Known species") &
          grepl("^[a-z]__.*unclassified", Order) &
          grepl("^[a-z]__.*unclassified", Class) &
          !grepl("^[a-z]__.*unclassified", Phylum)                     ~ paste0("f__unclassified ", Phylum, " family"),
        Family == "f__" & NewTax %in% c("New genus","New species","Known species") &
          grepl("^[a-z]__.*unclassified", Order) &
          grepl("^[a-z]__.*unclassified", Class) &
          grepl("^[a-z]__.*unclassified", Phylum) &
          !grepl("^[a-z]__.*unclassified", Domain)                     ~ paste0("f__unclassified ", Domain, " family"),
        Family == "f__" & NewTax %in% c("New genus","New species","Known species") &
          grepl("^[a-z]__.*unclassified", Order) &
          grepl("^[a-z]__.*unclassified", Class) &
          grepl("^[a-z]__.*unclassified", Phylum) &
          grepl("^[a-z]__.*unclassified", Domain)                      ~ paste0("f__unclassified ", Domain, " family"),
        TRUE                                                          ~ Family
      )
    ) %>%
    dplyr::mutate(
      ## 4. Genus
      Genus = dplyr::case_when(
        Genus != "g__" ~ Genus,
        Genus == "g__" & NewTax %in% c("New species","Known species") &
          !grepl("^[a-z]__.*unclassified", Family)                      ~ paste0("g__unclassified ", Family, " genus"),
        Genus == "g__" & NewTax %in% c("New species","Known species") &
          grepl("^[a-z]__.*unclassified", Family) &
          !grepl("^[a-z]__.*unclassified", Order)                       ~ paste0("g__unclassified ", Order, " genus"),
        Genus == "g__" & NewTax %in% c("New species","Known species") &
          grepl("^[a-z]__.*unclassified", Family) &
          grepl("^[a-z]__.*unclassified", Order) &
          !grepl("^[a-z]__.*unclassified", Class)                       ~ paste0("g__unclassified ", Class, " genus"),
        Genus == "g__" & NewTax %in% c("New species","Known species") &
          grepl("^[a-z]__.*unclassified", Family) &
          grepl("^[a-z]__.*unclassified", Order) &
          grepl("^[a-z]__.*unclassified", Class) &
          !grepl("^[a-z]__.*unclassified", Phylum)                      ~ paste0("g__unclassified ", Phylum, " genus"),
        Genus == "g__" & NewTax %in% c("New species","Known species") &
          grepl("^[a-z]__.*unclassified", Family) &
          grepl("^[a-z]__.*unclassified", Order) &
          grepl("^[a-z]__.*unclassified", Class) &
          grepl("^[a-z]__.*unclassified", Phylum) &
          !grepl("^[a-z]__.*unclassified", Domain)                      ~ paste0("g__unclassified ", Domain, " genus"),
        Genus == "g__" & NewTax %in% c("New species","Known species") &
          grepl("^[a-z]__.*unclassified", Family) &
          grepl("^[a-z]__.*unclassified", Order) &
          grepl("^[a-z]__.*unclassified", Class) &
          grepl("^[a-z]__.*unclassified", Phylum) &
          grepl("^[a-z]__.*unclassified", Domain)                       ~ paste0("g__unclassified ", Domain, " genus"),
        TRUE                                                           ~ Genus
      )
    ) %>%
    dplyr::mutate(across(c(Domain, Phylum, Class, Order, Family, Genus, Species), as.character))
  return(scale_data)
}



rename_strain_dataframe = function(data) {
  scale_data = data %>%
    mutate(Strain = case_when(
      Domain == "d__root domain"  ~ "t__root strain",
      is.na(Strain) ~ "t__",
      grepl("^t__.*unclassified", Strain)   ~ "t__",
      Strain == "no support"   ~ "t__",
      Strain == ""    ~ "t__",
      TRUE                              ~ Strain
    )) %>%
    dplyr::mutate(
      Strain = dplyr::case_when(
        Strain == "t__" & NewTax == "Unclassified" ~           paste0("t__",              Domain, " strain"),
        Strain == "t__" & NewTax == "Unclassified_Archaea" ~   paste0("t__",              Domain, " strain"),
        Strain == "t__" & NewTax == "Unclassified_Bacteria" ~  paste0("t__",              Domain, " strain"),
        Strain == "t__" & NewTax == "Unclassified_Eukaryota" ~ paste0("t__",              Domain, " strain"),
        Strain == "t__" & NewTax == "Unclassified_Virus" ~     paste0("t__",              Domain, " strain"),
        Strain == "t__" & NewTax == "New phylum"  ~            paste0("t__unclassified ", Domain, " strain"),
        Strain == "t__" & NewTax == "New class"   ~            paste0("t__unclassified ", Phylum, " strain"),
        Strain == "t__" & NewTax == "New order"   ~            paste0("t__unclassified ", Class,  " strain"),
        Strain == "t__" & NewTax == "New family"  ~            paste0("t__unclassified ", Order,  " strain"),
        Strain == "t__" & NewTax == "New genus"   ~            paste0("t__unclassified ", Family, " strain"),
        Strain == "t__" & NewTax == "New species" ~            paste0("t__unclassified ", Genus,  " strain"),
        Strain == "t__" & NewTax == "Known species" ~          paste0("t__unclassified ", Species," strain"),
        TRUE                                       ~ Strain
        )
    )
  return(scale_data)
}


