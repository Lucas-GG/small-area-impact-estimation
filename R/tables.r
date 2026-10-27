# Manuscript tables: one format for every exhibit table (results/tables/).
#   - 9 pt, fits a 6.5 in portrait text width (checked when LaTeX is available);
#   - booktabs style, no vertical rules; stub (row labels) left, numbers right;
#   - true minus signs; header labels may break onto two lines ("\n" in the label);
#   - the caption (results/CAPTIONS.md) is a brief title; every clarification (sample,
#     definitions, Monte Carlo SEs, reference values) is a source note under the table;
#   - written as HTML, LaTeX (.tex) and Word (.docx); a full-precision CSV is written by the caller.
# Usage: build <- function(fmt) {...gt object, labels via mt_label(x, fmt)...};
#        save_manuscript_table(build, "results/tables/table_name")

#' True minus signs in formatted numbers.
mt_minus <- function(x) gsub("-", "−", x, fixed = TRUE)

#' Column label that may break onto two lines: "\n" in `x` becomes a line break in HTML and
#' LaTeX, a space in Word (which wraps headers itself).
mt_label <- function(x, fmt) {
  if (!grepl("\n", x, fixed = TRUE)) return(x)
  switch(fmt,
    html = gt::html(gsub("\n", "<br>", x, fixed = TRUE)),
    latex = gsub("\n", "⏎", x, fixed = TRUE),   # placeholder, expanded after gtsave
    docx = gsub("\n", " ", x, fixed = TRUE))
}

#' Common table style: 9 pt in LaTeX (12 px in HTML), tight padding, no stub rule,
#' stub left and other columns right-aligned.
mt_style <- function(g, stub = NULL) {
  g <- gt::tab_options(g, table.font.size = gt::px(12), column_labels.font.size = gt::px(12),
    source_notes.font.size = gt::px(10.5), data_row.padding = gt::px(1), column_labels.padding = gt::px(2),
    row_group.padding = gt::px(1), source_notes.padding = gt::px(1), table.width = gt::pct(100),
    stub.border.style = "none", stub.border.width = gt::px(0))
  g <- gt::cols_align(g, "right", columns = gt::everything())
  if (!is.null(stub)) g <- gt::cols_align(g, "left", columns = tidyselect::all_of(stub))
  g
}

#' LaTeX post-processing of a gtsave()d table: column padding, no vertical rules, two-line
#' header cells, plain row-group labels (gt's \linewidth-wide group box plus padding overflows).
mt_fix_latex <- function(path, tabcolsep = "3.5pt") {
  tex <- readLines(path, encoding = "UTF-8")
  tex <- sub("\\fontsize{9.0pt}{11.0pt}\\selectfont",
    paste0("\\fontsize{9.0pt}{11.0pt}\\selectfont\\setlength\\tabcolsep{", tabcolsep, "}"), tex, fixed = TRUE)
  spec <- grep("\\begin{tabular", tex, fixed = TRUE)
  tex[spec] <- gsub("|", "", tex[spec], fixed = TRUE)
  hdr <- grepl("⏎", tex, fixed = TRUE)
  tex[hdr] <- vapply(tex[hdr], function(line) {
    cells <- strsplit(line, " & ", fixed = TRUE)[[1]]
    cells <- vapply(cells, function(cl) {
      if (!grepl("⏎", cl, fixed = TRUE)) return(cl)
      tail <- regmatches(cl, regexpr("\\s*\\\\\\\\\\s*$", cl))
      core <- sub("\\s*\\\\\\\\\\s*$", "", cl)
      paste0("\\begin{tabular}[b]{@{}r@{}}", gsub("⏎", "\\\\", core, fixed = TRUE), "\\end{tabular}",
        if (length(tail)) tail else "")
    }, "")
    paste(cells, collapse = " & ")
  }, "")
  grp <- grepl("\\parbox{\\linewidth}{", tex, fixed = TRUE)
  tex[grp] <- sub("}} \\\\", "} \\\\", sub("{>{\\raggedright\\arraybackslash}m{1\\linewidth}}{\\parbox{\\linewidth}{",
    "{l}{", tex[grp], fixed = TRUE), fixed = TRUE)
  writeLines(tex, path, useBytes = TRUE)
  invisible(path)
}

#' Reorder an OOXML element's children into schema order (Word rejects or rebuilds out-of-order ones).
mt_order_children <- function(node, order) {
  kids <- xml2::xml_children(node)
  ord <- order(match(xml2::xml_name(kids), order, nomatch = length(order) + 1L))
  for (n in kids) xml2::xml_remove(n)
  for (n in kids[ord]) xml2::xml_add_child(node, n)
  invisible(node)
}

#' Schema order of paragraph and cell properties and cell margins in every table of a parsed document.xml
#' (gt writes w:keepNext after w:spacing, and w:bottom before w:left in w:tblCellMar).
mt_order_xml <- function(doc) {
  ns <- c(w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main")
  ppr_order <- c("pStyle", "keepNext", "keepLines", "pageBreakBefore", "framePr", "widowControl", "numPr",
    "suppressLineNumbers", "pBdr", "shd", "tabs", "suppressAutoHyphens", "kinsoku", "wordWrap", "overflowPunct",
    "topLinePunct", "autoSpaceDE", "autoSpaceDN", "bidi", "adjustRightInd", "snapToGrid", "spacing", "ind",
    "contextualSpacing", "mirrorIndents", "suppressOverlap", "jc", "textDirection", "textAlignment",
    "textboxTightWrap", "outlineLvl", "divId", "cnfStyle", "rPr", "sectPr", "pPrChange")
  mar_order <- c("top", "start", "left", "bottom", "end", "right")
  tcpr_order <- c("cnfStyle", "tcW", "gridSpan", "hMerge", "vMerge", "tcBorders", "shd", "noWrap", "tcMar",
    "textDirection", "tcFitText", "vAlign", "hideMark")
  for (n in xml2::xml_find_all(doc, "//w:tbl//w:pPr", ns)) mt_order_children(n, ppr_order)
  for (n in xml2::xml_find_all(doc, "//w:tcPr", ns)) mt_order_children(n, tcpr_order)
  for (n in xml2::xml_find_all(doc, "//w:tblCellMar|//w:tcMar", ns)) mt_order_children(n, mar_order)
  invisible(doc)
}

#' Apply mt_order_xml() to a saved .docx in place; changes nothing else (widths, fonts).
mt_fix_order <- function(path) {
  dir <- tempfile("mtdocx"); dir.create(dir)
  utils::unzip(path, exdir = dir)
  doc_path <- file.path(dir, "word", "document.xml")
  doc <- mt_order_xml(xml2::read_xml(doc_path))
  xml2::write_xml(doc, doc_path, options = "as_xml")
  files <- list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  out <- normalizePath(path, mustWork = FALSE); unlink(out)
  zip::zip(out, files, root = dir, mode = "mirror")
  unlink(dir, recursive = TRUE)
  invisible(path)
}

#' Word post-processing of a gtsave()d table. gt writes no column grid or cell widths and puts
#' tcPr children out of schema order, so Word rebuilds the layout and the table falls apart on
#' any edit. This adds a w:tblGrid and w:tcW widths filling `textwidth_in` (columns sized by their
#' longest single-column entry), orders tcPr children, drops vertical rules, uses left/right in
#' place of start/end, lets the table break across pages (keep-with-next only on header rows), and
#' sets the font size to 9 pt. Idempotent.
mt_fix_docx <- function(path, textwidth_in = 6.5, font_pt = 9) {
  dir <- tempfile("mtdocx"); dir.create(dir)
  utils::unzip(path, exdir = dir)
  doc_path <- file.path(dir, "word", "document.xml")
  doc <- xml2::read_xml(doc_path)
  ns <- c(w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main")
  wa <- function(node, attr) xml2::xml_attr(node, paste0("w:", attr), ns)
  tcpr_order <- c("cnfStyle", "tcW", "gridSpan", "hMerge", "vMerge", "tcBorders", "shd", "noWrap", "tcMar",
    "textDirection", "tcFitText", "vAlign", "hideMark")
  total <- round(textwidth_in * 1440)
  for (tbl in xml2::xml_find_all(doc, "//w:tbl", ns)) {
    rows <- xml2::xml_find_all(tbl, "w:tr", ns)
    span_of <- function(tc) { g <- xml2::xml_find_first(tc, "w:tcPr/w:gridSpan", ns)
      if (inherits(g, "xml_missing")) 1L else as.integer(wa(g, "val")) }
    ncol <- max(vapply(rows, function(r) sum(vapply(xml2::xml_find_all(r, "w:tc", ns), span_of, 1L)), 1L))
    len <- rep(8, ncol)   # floor: a column always fits a short label (e.g. "ETWFE") on one line
    for (r in rows) {
      pos <- 1L
      for (tc in xml2::xml_find_all(r, "w:tc", ns)) {
        k <- span_of(tc)
        if (k == 1L) len[pos] <- max(len[pos], nchar(xml2::xml_text(tc)))
        pos <- pos + k
      }
    }
    w <- floor(total * (len + 2) / sum(len + 2)); w[ncol] <- total - sum(w[-ncol])
    tblpr <- xml2::xml_find_first(tbl, "w:tblPr", ns)
    for (n in xml2::xml_find_all(tbl, "w:tblGrid", ns)) xml2::xml_remove(n)
    for (n in xml2::xml_find_all(tblpr, "w:tblW|w:tblLayout", ns)) xml2::xml_remove(n)
    xml2::xml_add_child(tblpr, "w:tblW", "w:w" = total, "w:type" = "dxa", .where = 0)
    xml2::xml_add_child(tblpr, "w:tblLayout", "w:type" = "fixed")
    grid <- xml2::xml_add_sibling(tblpr, "w:tblGrid", .where = "after")
    for (x in w) xml2::xml_add_child(grid, "w:gridCol", "w:w" = x)
    for (r in rows) {
      pos <- 1L
      for (tc in xml2::xml_find_all(r, "w:tc", ns)) {
        k <- span_of(tc)
        pr <- xml2::xml_find_first(tc, "w:tcPr", ns)
        if (inherits(pr, "xml_missing")) pr <- xml2::xml_add_child(tc, "w:tcPr", .where = 0)
        for (n in xml2::xml_find_all(pr, "w:tcW", ns)) xml2::xml_remove(n)
        xml2::xml_add_child(pr, "w:tcW", "w:w" = sum(w[pos:(pos + k - 1L)]), "w:type" = "dxa")
        # No vertical rules.
        for (n in xml2::xml_find_all(pr, "w:tcBorders/w:start|w:tcBorders/w:end|w:tcBorders/w:left|w:tcBorders/w:right", ns))
          xml2::xml_remove(n)
        mt_order_children(pr, tcpr_order)
        pos <- pos + k
      }
    }
    for (n in xml2::xml_find_all(tblpr, "w:tblCellMar/w:start|w:tblCellMar/w:end", ns))
      xml2::xml_set_name(n, if (xml2::xml_name(n) == "start") "left" else "right")   # keeps the w namespace
    mt_order_children(tblpr, c("tblStyle", "tblpPr", "tblOverlap", "bidiVisual", "tblStyleRowBandSize",
      "tblStyleColBandSize", "tblW", "jc", "tblCellSpacing", "tblInd", "tblBorders", "shd", "tblLayout",
      "tblCellMar", "tblLook"))
    # gt sets keep-with-next on every paragraph, which makes the whole table one unbreakable block:
    # Word moves it to the next page and leaves the caption behind. Keep only the header rows
    # (repeated on each page) with the first body row; rows themselves still do not split.
    for (r in rows) if (inherits(xml2::xml_find_first(r, "w:trPr/w:tblHeader", ns), "xml_missing"))
      for (n in xml2::xml_find_all(r, ".//w:pPr/w:keepNext", ns)) xml2::xml_remove(n)
    for (jc in xml2::xml_find_all(tbl, ".//w:pPr/w:jc", ns)) {
      v <- wa(jc, "val")
      if (identical(v, "start")) xml2::xml_set_attr(jc, "w:val", "left")
      if (identical(v, "end")) xml2::xml_set_attr(jc, "w:val", "right")
    }
    for (sz in xml2::xml_find_all(tbl, ".//w:rPr/w:sz", ns)) xml2::xml_set_attr(sz, "w:val", 2 * font_pt)
  }
  mt_order_xml(doc)
  xml2::write_xml(doc, doc_path, options = "as_xml")
  xml2::read_xml(doc_path)   # fails here, before the .docx is replaced, if the XML is not well formed
  files <- list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  out <- normalizePath(path, mustWork = FALSE); unlink(out)
  zip::zip(out, files, root = dir, mode = "mirror")
  unlink(dir, recursive = TRUE)
  invisible(path)
}

#' Natural width (inches) of a LaTeX table at its own font size, and the number of overfull
#' boxes when set in a `textwidth_in` text block. NA when xelatex is unavailable.
mt_latex_width <- function(path, textwidth_in = 6.5) {
  xelatex <- Sys.which("xelatex")
  if (!nzchar(xelatex)) return(list(width_in = NA_real_, overfull = NA_integer_))
  t <- paste(readLines(path, encoding = "UTF-8"), collapse = "\n")
  dir <- tempfile("mtwidth"); dir.create(dir)
  body <- substr(t, regexpr("\\begin{tabular*}", t, fixed = TRUE), nchar(t))
  body <- substr(body, 1, regexpr("\\end{tabular*}", body, fixed = TRUE) + nchar("\\end{tabular*}") - 1)
  body <- sub("\\begin{tabular*}{1\\linewidth}{@{\\extracolsep{\\fill}}", "\\begin{tabular}{", body, fixed = TRUE)
  body <- sub("\\end{tabular*}", "\\end{tabular}", body, fixed = TRUE)
  pre <- regmatches(t, regexpr("\\\\fontsize\\{[^\n]*", t))
  pkgs <- "\\usepackage{fontspec}\\usepackage{booktabs}\\usepackage{array}\\usepackage{caption}"
  writeLines(c(paste0("\\documentclass{article}", pkgs, "\\begin{document}\\newsavebox\\tb\\sbox\\tb{", pre, body, "}"),
    "\\typeout{TABLEWIDTH=\\the\\wd\\tb}\\end{document}"), file.path(dir, "m.tex"), useBytes = TRUE)
  writeLines(c(paste0("\\documentclass{article}\\usepackage[textwidth=", textwidth_in, "in]{geometry}", pkgs,
    "\\begin{document}"), t, "\\end{document}"), file.path(dir, "f.tex"), useBytes = TRUE)
  run <- function(f) system2(xelatex, c("-interaction=nonstopmode", paste0("-output-directory=", dir), file.path(dir, f)),
    stdout = FALSE, stderr = FALSE)
  run("m.tex"); run("f.tex")
  mlog <- readLines(file.path(dir, "m.log"), warn = FALSE); flog <- readLines(file.path(dir, "f.log"), warn = FALSE)
  w <- as.numeric(sub(".*TABLEWIDTH=([0-9.]+)pt.*", "\\1", grep("TABLEWIDTH=", mlog, value = TRUE)[1]))
  list(width_in = w / 72.27, overfull = sum(grepl("^Overfull \\\\hbox", flog)))
}

#' Write a manuscript table as HTML, LaTeX and Word from a builder function(fmt) that returns
#' a gt object; returns the LaTeX width check and warns when the table exceeds the text width.
save_manuscript_table <- function(build, stem, textwidth_in = 6.5, tabcolsep = "3.5pt") {
  gt::gtsave(build("html"), paste0(stem, ".html"))
  gt::gtsave(build("docx"), paste0(stem, ".docx"))
  mt_fix_docx(paste0(stem, ".docx"), textwidth_in)
  tex <- paste0(stem, ".tex")
  gt::gtsave(build("latex"), tex)
  mt_fix_latex(tex, tabcolsep)
  chk <- mt_latex_width(tex, textwidth_in)
  if (isTRUE(chk$width_in > textwidth_in) || isTRUE(chk$overfull > 0))
    warning(sprintf("%s: %.2f in wide (%d overfull) for a %.1f in text width", basename(stem),
      chk$width_in, chk$overfull, textwidth_in), call. = FALSE)
  invisible(chk)
}
