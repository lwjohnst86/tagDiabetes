# Steps to take in the Git repo:
# 1. Make new orphan branch (git checkout --orphan code)
# 2. Keep only R files and purl('vignette/maindocument.Rmd')
# 3. git commit the files
# 4. git archive --format=zip --output=code-archive.zip code
# 5. Then run the below code
send_to_figshare <- function() {
    rfigshare::fs_new_article(
        "Code for the analysis of TAGFA on diabetes pathogenesis",
        "Analysis code for a manuscript on the role of triacylglyceride fatty acid composition on the pathogenesis of diabetes, based on a longitudinal cohort (PROMISE)",
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
