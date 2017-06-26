# Steps to take in the Git repo:
# 1. Make new orphan branch (git checkout --orphan code)
# 2. Keep only R files and purl('vignette/maindocument.Rmd')
# 3. git commit the files
# 4. git archive --format=zip --output=code-archive.zip code
# 5. Then run the below code
send_to_figshare <- function() {
    rfigshare::fs_new_article(
        "Analysis code for poster on triacylglycerol DNL on components of the Metabolic Syndrome in PROMISE",
        "Analysis code for a poster for the American Diabetes Association Scientific Sessions in 2017",
        type = "fileset",
        # there was a problem with the authors...
        # authors = c("First Last"),
        tags = c(
            "insulin sensitivity",
            "beta-cell function",
            "fatty acids",
            "triacylglycerols",
            "cohort",
            "longitudinal"
        ),
        categories = c("Diseases", "Pathogenesis", "Epidemiology"),
        links = "https://github.com/lwjohnst86/tagDiabetes",
        files = "code-archive.zip",
        visibility = "draft"
    )
}
