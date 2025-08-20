# ===== Setup: מתקין/טוען חבילות רק אם חסרות =====
ensure_pkg <- function(pkgs) {
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss)) {
    install.packages(miss, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  invisible(lapply(pkgs, require, character.only = TRUE))
}

ensure_pkg(c("tidyverse", "forcats"))

# סגנון בסיס לגרפים שתומך בעברית (בחר פונט קיים אצלך, למשל Arial)
ggplot2::theme_set(ggplot2::theme_minimal(base_family = "Arial"))

# ===== קריאת הדאטה מ-GitHub =====
gh_raw <- "https://raw.githubusercontent.com/Roybin12/data-ISOC-IL/refs/heads/main/data_processed.csv"

# מומלץ: readr שומר שמות עמודות כמו שהם, כולל עברית ורווחים
data <- readr::read_csv(gh_raw,
                        locale = readr::locale(encoding = "UTF-8"),
                        show_col_types = FALSE)

#------------------------------------------------------------------------------
# --- שלב הכנה לניתוח: ניקוי, טיפוסי משתנים, וקבוצות גיל ---

# עמודות הדמוגרפיה ושאלות הסקר כמצופה בדאטה
demographic_cols <- c("מגדר","גיל","לאום","הגדרה דתית")
survey_cols      <- c("סוג חיבור האינטרנט","ספק האינטרנט",
                      "ספק טלויזיה ישראלי","ספק טלוויזיה בינלאומי",
                      "צריכת שירותי מוזיקה")

# 0) הגיינה לשמות עמודות, למקרה שהגיעו עם נקודות במקום רווחים
names(data) <- trimws(gsub("\\s+", " ", gsub("\\.", " ", names(data))))

# 1) ולידציה: אם חסר משהו נעצור עם הודעה ברורה
missing_cols <- setdiff(c(demographic_cols, survey_cols), names(data))
if (length(missing_cols)) {
  stop("חסרות העמודות: ", paste(missing_cols, collapse = ", "))
}

# 2) ניקוי ערכי מחרוזת: רווחים מיותרים -> NA
data <- data |>
  dplyr::mutate(dplyr::across(where(is.character), ~ dplyr::na_if(trimws(.), "")))

# 3) המרת עמודות סקר לפקטורים עם תווית חסר בעברית
data <- data |>
  dplyr::mutate(dplyr::across(dplyr::all_of(survey_cols),
                              ~ forcats::fct_explicit_na(as.factor(.), na_level = "— חסר —")))

# 4) דמוגרפיה: מגדר/לאום/הגדרה דתית כפקטורים; גיל מספרי + יצירת קבוצות גיל
data <- data |>
  dplyr::mutate(
    dplyr::across(c("מגדר","לאום","הגדרה דתית"),
                  ~ forcats::fct_explicit_na(as.factor(.), na_level = "— חסר —")),
    גיל = suppressWarnings(as.numeric(גיל)),
    קבוצת_גיל = cut(
      גיל,
      breaks = c(-Inf, 24, 34, 44, 54, Inf),
      labels = c("עד 24","25–34","35–44","45–54","55+")
    ),
    קבוצת_גיל = forcats::fct_explicit_na(קבוצת_גיל, na_level = "— חסר —")
  )

# 5) סדר רמות לפי שכיחות עבור משתני הסקר, שימושי לגרפים
data <- data |>
  dplyr::mutate(
    `סוג חיבור האינטרנט`   = forcats::fct_infreq(`סוג חיבור האינטרנט`),
    `ספק האינטרנט`         = forcats::fct_infreq(`ספק האינטרנט`),
    `ספק טלויזיה ישראלי`   = forcats::fct_infreq(`ספק טלויזיה ישראלי`),
    `ספק טלויזיה בינלאומי` = forcats::fct_infreq(`ספק טלוויזיה בינלאומי`),
    `צריכת שירותי מוזיקה`  = forcats::fct_infreq(`צריכת שירותי מוזיקה`)
  )

# 6) פונקציית עזר אופציונלית לגרפים: צמצום ל-Top-N + עטיפת תוויות (לגרף נקי כשיש הרבה קטגוריות)
collapse_top_n <- function(x, N = 8, other_label = "אחר", wrap = 18) {
  x_fac <- forcats::fct_explicit_na(as.factor(x), na_level = "— חסר —")
  top_levels <- names(sort(table(x_fac[x_fac != "— חסר —"]), decreasing = TRUE))[seq_len(min(N, length(unique(x_fac))))]
  x_top <- forcats::fct_other(x_fac, keep = as.character(top_levels), other_level = other_label)
  list(
    factor = x_top,
    label  = stringr::str_wrap(as.character(x_top), width = wrap)
  )
}

message("הדאטה מוכן לניתוח: פקטורים הוגדרו, קבוצות גיל נוצרו, ושמות העמודות נקיים ✅")
#------------------------------------------------------------------------------
# עזר לחישוב p-value ו-Cramer's V
chi_cramers <- function(x, y) {
  tab <- table(x, y, useNA = "no")
  chi1 <- suppressWarnings(chisq.test(tab))
  chi  <- if (any(chi1$expected < 5, na.rm = TRUE))
    suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 10000))
  else chi1
  n <- sum(tab); k <- min(nrow(tab)-1, ncol(tab)-1)
  V <- sqrt(as.numeric(chi$statistic) / (n * k))
  list(chi = chi, p = as.numeric(chi$p.value), V = V)
}

# עזר לצמצום קטגוריות ל-Top-N + תוויות עטופות (לגרפים כשיש הרבה ערכים)
collapse_top_n <- function(x, N = 8, other_label = "אחר", wrap = 18) {
  x_fac <- forcats::fct_explicit_na(as.factor(x), na_level = "— חסר —")
  top_levels <- names(sort(table(x_fac[x_fac != "— חסר —"]), decreasing = TRUE))[seq_len(min(N, length(unique(x_fac))))]
  x_top <- forcats::fct_other(x_fac, keep = as.character(top_levels), other_level = other_label)
  list(factor = x_top, label = stringr::str_wrap(as.character(x_top), width = wrap))
}

#------------------------------------------------------------------------------
#שלב1
var <- "סוג חיבור האינטרנט"
demographic_cols <- c("מגדר","קבוצת_גיל","לאום","הגדרה דתית")

# טבלת שכיחויות כללית
freq_tbl_conn <- data |>
  dplyr::count(.data[[var]], name = "n") |>
  dplyr::mutate(pct = 100 * n / sum(n)) |>
  dplyr::arrange(dplyr::desc(n))
print(freq_tbl_conn)

# גרף כללי
p_conn_all <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(.data[[var]]))) +
  ggplot2::geom_bar(width = 0.75) +
  ggplot2::labs(title = "התפלגות סוג חיבור האינטרנט", x = "סוג חיבור", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::coord_flip()
print(p_conn_all)

# גרף לפי מגדר (אפשר לשכפל גם לשאר הדמוגרפיות)
p_conn_gender <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(.data[[var]]))) +
  ggplot2::geom_bar() +
  ggplot2::facet_wrap(~ `מגדר`, ncol = 2, scales = "free_y") +
  ggplot2::labs(title = "סוג חיבור האינטרנט לפי מגדר", x = "סוג חיבור", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::coord_flip()
print(p_conn_gender)

# סיכום משמעותיות (תמציתי)
summ1 <- do.call(rbind, lapply(demographic_cols, function(g) {
  out <- chi_cramers(data[[g]], data[[var]])
  c(group = g, p_value = out$p, cramers_v = out$V)
}))
summ1_df <- data.frame(
  השוואה = paste(var, "×", summ1[, "group"]),
  `משתנה דמוגרפי` = summ1[, "group"],
  `p-ערך` = ifelse(as.numeric(summ1[, "p_value"]) < 0.001, "<0.001", sprintf("%.3f", as.numeric(summ1[, "p_value"]))),
  `Cramer's V` = sprintf("%.3f", as.numeric(summ1[, "cramers_v"])),
  מובהקות = ifelse(as.numeric(summ1[, "p_value"]) < 0.05, "מובהק", "לא מובהק"),
  check.names = FALSE
)
print(summ1_df, row.names = FALSE)

# לפי קבוצת גיל
p_conn_age <- ggplot(data, aes(x = forcats::fct_infreq(`סוג חיבור האינטרנט`))) +
  geom_bar() +
  facet_wrap(~ `קבוצת_גיל`, ncol = 2, scales = "free_y") +
  labs(title = "סוג חיבור האינטרנט לפי קבוצת גיל", x = "סוג חיבור", y = "תדירות") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()
print(p_conn_age)

# לפי לאום
p_conn_eth <- ggplot(data, aes(x = forcats::fct_infreq(`סוג חיבור האינטרנט`))) +
  geom_bar() +
  facet_wrap(~ `לאום`, ncol = 2, scales = "free_y") +
  labs(title = "סוג חיבור האינטרנט לפי לאום", x = "סוג חיבור", y = "תדירות") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()
print(p_conn_eth)

# לפי הגדרה דתית
p_conn_relig <- ggplot(data, aes(x = forcats::fct_infreq(`סוג חיבור האינטרנט`))) +
  geom_bar() +
  facet_wrap(~ `הגדרה דתית`, ncol = 2, scales = "free_y") +
  labs(title = "סוג חיבור האינטרנט לפי הגדרה דתית", x = "סוג חיבור", y = "תדירות") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()
print(p_conn_relig)


#------------------------------------------------------------------------------
#שלב2

var <- "ספק האינטרנט"
demographic_cols <- c("מגדר","קבוצת_גיל","לאום","הגדרה דתית")

freq_tbl_isp <- data |>
  dplyr::count(.data[[var]], name = "n") |>
  dplyr::mutate(pct = 100 * n / sum(n)) |>
  dplyr::arrange(dplyr::desc(n))
print(freq_tbl_isp)

p_isp_all <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(.data[[var]]))) +
  ggplot2::geom_bar(width = 0.75) +
  ggplot2::labs(title = "התפלגות ספק האינטרנט", x = "ספק", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::coord_flip()
print(p_isp_all)

p_isp_gender <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(.data[[var]]))) +
  ggplot2::geom_bar() +
  ggplot2::facet_wrap(~ `מגדר`, ncol = 2, scales = "free_y") +
  ggplot2::labs(title = "ספק האינטרנט לפי מגדר", x = "ספק", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::coord_flip()
print(p_isp_gender)

summ2 <- do.call(rbind, lapply(demographic_cols, function(g) {
  out <- chi_cramers(data[[g]], data[[var]])
  c(group = g, p_value = out$p, cramers_v = out$V)
}))
summ2_df <- data.frame(
  השוואה = paste(var, "×", summ2[, "group"]),
  `משתנה דמוגרפי` = summ2[, "group"],
  `p-ערך` = ifelse(as.numeric(summ2[, "p_value"]) < 0.001, "<0.001", sprintf("%.3f", as.numeric(summ2[, "p_value"]))),
  `Cramer's V` = sprintf("%.3f", as.numeric(summ2[, "cramers_v"])),
  מובהקות = ifelse(as.numeric(summ2[, "p_value"]) < 0.05, "מובהק", "לא מובהק"),
  check.names = FALSE
)
print(summ2_df, row.names = FALSE)

# לפי קבוצת גיל
p_isp_age <- ggplot(data, aes(x = forcats::fct_infreq(`ספק האינטרנט`))) +
  geom_bar() +
  facet_wrap(~ `קבוצת_גיל`, ncol = 2, scales = "free_y") +
  labs(title = "ספק האינטרנט לפי קבוצת גיל", x = "ספק", y = "תדירות") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()
print(p_isp_age)

# לפי לאום
p_isp_eth <- ggplot(data, aes(x = forcats::fct_infreq(`ספק האינטרנט`))) +
  geom_bar() +
  facet_wrap(~ `לאום`, ncol = 2, scales = "free_y") +
  labs(title = "ספק האינטרנט לפי לאום", x = "ספק", y = "תדירות") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()
print(p_isp_eth)

# לפי הגדרה דתית
p_isp_relig <- ggplot(data, aes(x = forcats::fct_infreq(`ספק האינטרנט`))) +
  geom_bar() +
  facet_wrap(~ `הגדרה דתית`, ncol = 2, scales = "free_y") +
  labs(title = "ספק האינטרנט לפי הגדרה דתית", x = "ספק", y = "תדירות") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()
print(p_isp_relig)


#------------------------------------------------------------------------------
#שלב 3
var <- "ספק טלויזיה ישראלי"
demographic_cols <- c("מגדר","קבוצת_גיל","לאום","הגדרה דתית")

# כל הערכים – גרף כללי
p_tv_il_all_full <- ggplot2::ggplot(
  data |> dplyr::mutate(lbl = stringr::str_wrap(as.character(forcats::fct_explicit_na(as.factor(.data[[var]]), "— חסר —")), width = 18)),
  ggplot2::aes(x = forcats::fct_infreq(lbl))
) +
  ggplot2::geom_bar(width = 0.7) +
  ggplot2::labs(title = "התפלגות ספק הטלויזיה הישראלי (כל הערכים)", x = "ספק", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = 7)) +
  ggplot2::coord_flip()
print(p_tv_il_all_full)

# גרסת Top-N + "אחר" (מומלץ לפרסום)
tmp_il <- collapse_top_n(data[[var]], N = 8)
p_tv_il_all_top <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(tmp_il$label))) +
  ggplot2::geom_bar(width = 0.75) +
  ggplot2::labs(title = "התפלגות ספק הטלויזיה הישראלי (Top-8 + אחר)", x = "ספק", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::coord_flip()
print(p_tv_il_all_top)

# מפורק לפי מגדר (Top-N)
p_tv_il_gender <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(tmp_il$label))) +
  ggplot2::geom_bar() +
  ggplot2::facet_wrap(~ `מגדר`, ncol = 2, scales = "free_y") +
  ggplot2::labs(title = "ספק הטלוויזיה הישראלי לפי מגדר (Top-8 + אחר)", x = "ספק", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::coord_flip()
print(p_tv_il_gender)

# טבלת שכיחויות תמציתית (Top-N + אחר)
freq_tbl_tv_il <- data |>
  dplyr::count(ספק = tmp_il$factor, name = "n") |>
  dplyr::mutate(pct = round(100 * n / sum(n), 1)) |>
  dplyr::arrange(dplyr::desc(n))
print(freq_tbl_tv_il)

# נכין את אובייקט ה-Top-N אם לא קיים
if (!exists("tmp_il")) {
  tmp_il <- collapse_top_n(data$`ספק טלויזיה ישראלי`, N = 8)  # שנה N לפי הצורך
}

# לפי קבוצת גיל
p_tv_il_age <- ggplot(data, aes(x = forcats::fct_infreq(tmp_il$label))) +
  geom_bar() +
  facet_wrap(~ `קבוצת_גיל`, ncol = 2, scales = "free_y") +
  labs(title = "ספק הטלוויזיה הישראלי לפי קבוצת גיל (Top-N + אחר)", x = "ספק", y = "תדירות") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 8)) +
  coord_flip()
print(p_tv_il_age)

# לפי לאום
p_tv_il_eth <- ggplot(data, aes(x = forcats::fct_infreq(tmp_il$label))) +
  geom_bar() +
  facet_wrap(~ `לאום`, ncol = 2, scales = "free_y") +
  labs(title = "ספק הטלוויזיה הישראלי לפי לאום (Top-N + אחר)", x = "ספק", y = "תדירות") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 8)) +
  coord_flip()
print(p_tv_il_eth)

# לפי הגדרה דתית
p_tv_il_relig <- ggplot(data, aes(x = forcats::fct_infreq(tmp_il$label))) +
  geom_bar() +
  facet_wrap(~ `הגדרה דתית`, ncol = 2, scales = "free_y") +
  labs(title = "ספק הטלוויזיה הישראלי לפי הגדרה דתית (Top-N + אחר)", x = "ספק", y = "תדירות") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 8)) +
  coord_flip()
print(p_tv_il_relig)


#------------------------------------------------------------------------------
#שלב 4

var <- "ספק טלויזיה בינלאומי"
demographic_cols <- c("מגדר","קבוצת_גיל","לאום","הגדרה דתית")

# כל הערכים – גרף כללי
p_tv_int_all_full <- ggplot2::ggplot(
  data |> dplyr::mutate(lbl = stringr::str_wrap(as.character(forcats::fct_explicit_na(as.factor(.data[[var]]), "— חסר —")), width = 18)),
  ggplot2::aes(x = forcats::fct_infreq(lbl))
) +
  ggplot2::geom_bar(width = 0.7) +
  ggplot2::labs(title = "התפלגות ספק הטלוויזיה הבינלאומי (כל הערכים)", x = "ספק", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = 7)) +
  ggplot2::coord_flip()
print(p_tv_int_all_full)

# גרסת Top-N + "אחר"
tmp_int <- collapse_top_n(data[[var]], N = 8)
p_tv_int_all_top <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(tmp_int$label))) +
  ggplot2::geom_bar(width = 0.75) +
  ggplot2::labs(title = "ספק הטלוויזיה הבינלאומי (Top-8 + אחר)", x = "ספק", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::coord_flip()
print(p_tv_int_all_top)

# מפורק לפי קבוצת גיל (Top-N)
p_tv_int_age <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(tmp_int$label))) +
  ggplot2::geom_bar() +
  ggplot2::facet_wrap(~ `קבוצת_גיל`, ncol = 2, scales = "free_y") +
  ggplot2::labs(title = "ספק הטלוויזיה הבינלאומי לפי קבוצת גיל (Top-8 + אחר)", x = "ספק", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::coord_flip()
print(p_tv_int_age)

# טבלת שכיחויות תמציתית
freq_tbl_tv_int <- data |>
  dplyr::count(ספק = tmp_int$factor, name = "n") |>
  dplyr::mutate(pct = round(100 * n / sum(n), 1)) |>
  dplyr::arrange(dplyr::desc(n))
print(freq_tbl_tv_int)

# --- שלב 4: ספק טלויזיה בינלאומי לפי מגדר ---

var <- "ספק טלויזיה בינלאומי"
N   <- 8  # שנה לפי הצורך

# אם טרם הוגדר אובייקט ה-Top-N לשלב 4:
if (!exists("tmp_int")) {
  tmp_int <- collapse_top_n(data[[var]], N = N)
}

# 4A) גרף Top-N + "אחר" לפי מגדר (מומלץ למסמך)
p_tv_int_gender <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(tmp_int$label))) +
  ggplot2::geom_bar() +
  ggplot2::facet_wrap(~ `מגדר`, ncol = 2, scales = "free_y") +
  ggplot2::labs(title = sprintf("ספק הטלוויזיה הבינלאומי לפי מגדר (Top-%d + אחר)", N),
                x = "ספק", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                 axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::coord_flip()
print(p_tv_int_gender)

#------------------------------------------------------------------------------
#שלב5

var <- "צריכת שירותי מוזיקה"
demographic_cols <- c("מגדר","קבוצת_גיל","לאום","הגדרה דתית")

# הצגה מלאה (כל הערכים)
p_music_all_full <- ggplot2::ggplot(
  data |> dplyr::mutate(lbl = stringr::str_wrap(as.character(forcats::fct_explicit_na(as.factor(.data[[var]]), "— חסר —")), width = 18)),
  ggplot2::aes(x = forcats::fct_infreq(lbl))
) +
  ggplot2::geom_bar(width = 0.7) +
  ggplot2::labs(title = "צריכת שירותי מוזיקה (כל הערכים)", x = "שירות", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = 7)) +
  ggplot2::coord_flip()
print(p_music_all_full)

# גרסת Top-N + "אחר"
tmp_music <- collapse_top_n(data[[var]], N = 8)
p_music_all_top <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(tmp_music$label))) +
  ggplot2::geom_bar(width = 0.75) +
  ggplot2::labs(title = "צריכת שירותי מוזיקה (Top-8 + אחר)", x = "שירות", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::coord_flip()
print(p_music_all_top)

# מפורק לפי מגדר (Top-N)
p_music_gender <- ggplot2::ggplot(data, ggplot2::aes(x = forcats::fct_infreq(tmp_music$label))) +
  ggplot2::geom_bar() +
  ggplot2::facet_wrap(~ `מגדר`, ncol = 2, scales = "free_y") +
  ggplot2::labs(title = "צריכת שירותי מוזיקה לפי מגדר (Top-8 + אחר)", x = "שירות", y = "תדירות") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::coord_flip()
print(p_music_gender)

# טבלת שכיחויות תמציתית
freq_tbl_music <- data |>
  dplyr::count(שירות = tmp_music$factor, name = "n") |>
  dplyr::mutate(pct = round(100 * n / sum(n), 1)) |>
  dplyr::arrange(dplyr::desc(n))
print(freq_tbl_music)
