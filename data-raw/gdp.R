gdp <-
    readr::read_csv("~/Desktop/Home/Study/Columbia/Research/Richard/Fokianos/simulation study/TED.csv") |>
    dplyr::select(-c(REGION, MEASURE)) |>
    tidyr::pivot_longer(
        -c(ISO, COUNTRY, INDICATOR),
        names_to = "YEAR",
        values_to = "values"
    ) |>
    dplyr::mutate("YEAR" = gsub("X", "", YEAR) |> as.numeric()) |>
    dplyr::mutate("values" = gsub(",", "", values) |> as.numeric()) |>
    dplyr::filter(
        COUNTRY %in% c(
            "Austria", "Belgium", "Denmark",
            "Finland", "France", "Germany",
            "Greece", "Iceland", "Ireland",
            "Italy", "Luxembourg", "Netherlands",
            "Norway", "Portugal", "Spain",
            "Sweden", "Switzerland", "United Kingdom",
            "Canada", "United States", "Australia",
            "New Zealand", "Japan"
        )
    ) |>
    dplyr::rename_with(tolower) |>
    dplyr::filter(indicator == "Real GDP", year >= 1980, year <= 2019) |>
    dplyr::select(-iso, -indicator) |>
    dplyr::rename(gdp = values) |>
    tsibble::as_tsibble(index = year, key = country)

usethis::use_data(gdp, overwrite = TRUE)
