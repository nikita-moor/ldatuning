packages <- c("ldatuning", "tm")
lapply(packages, require, character.only = TRUE)
rm(packages)

data("crude", package = "tm")

dtm <- tm::DocumentTermMatrix(crude)

LDA_models <- FindTopicsNumber(
  dtm,
  topics = 8:11,
  metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
  return_models = TRUE
)

LDA_models
