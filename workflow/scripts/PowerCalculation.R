library(genpwr)
pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=1000, Case.Rate=c(0.005,0.015,0.025), k=NULL,
                  MAF=seq(0.05, 0.5, 0.01), OR=c(6,10,15),Alpha=5E-8,
                  True.Model='Additive', 
                  Test.Model='Additive')
pw$OR <- factor(pw$OR)

ggplot2::ggplot(pw) + aes(x=MAF,y=`Power_at_Alpha_5e-08`,color = factor(Case.Rate)) + geom_line() + facet_grid(~OR,labeller = label_both) +
  ylab('Power') + labs(color = 'Mtb MAF') + xlab('Human MAF')
