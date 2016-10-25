# # test the package!
#
# needs::needs(dplyr, purrr, boot)
#
# # load source
# source("R/bootr.r")
# source("R/bootr_output.r")
# data = MatchIt::lalonde
# f = re78 ~ treat + educ + black + age + nodegree
# fm = ~ treat + educ + black + age+ nodegree
# lm(f, data) %>% summary
#
# out = bootr(numruns = 1000, formula = f, formula_match = fm, data = data)
