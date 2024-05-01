

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(survival)


# Load and Treatment Data -------------------------------------------------

data <- read.csv(file = "recuperaInadimplencia.csv")
data <- data[,2:4]



# Descptive Analysis ------------------------------------------------------

# histograma da var resposta Tempo até o Pagamento da dívida
hist(data$tempo_pag, xlab="Tempo Pagamento", 
     main = "Histograma do Tempo de Pagamento")

# curva de Kaplan-Meier estimada
ekm <- survfit(Surv(time = tempo_pag, event = status) ~ 1, data = data)
ekm_df <- data.frame(time = ekm$time,
                     surv = ekm$surv)

ggplot(data = ekm_df, mapping = aes(time, surv)) + 
    geom_line() + labs(x = "Tempos", y = "S(t)") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    geom_hline(yintercept = 0.76, linetype="dashed", color = "red") +
    geom_hline(yintercept = 0.21, linetype="dashed", color = "red") +
    geom_text(label = "1-p0", x = 0, y = 0.81, size = 8/.pt) +
    geom_text(label = "p1", x = 0, y = 0.19, size = 8/.pt) +
    theme_bw()

# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo_pag, event = status) ~ consulta, data = data)
ekm_df <- data.frame(time = ekm$time,
                     surv = ekm$surv,
                     consulta = as.factor(c(rep(1,44),rep(2,61))))

ggplot(data = ekm_df, mapping = aes(time, surv)) + 
    geom_step(aes(color = consulta)) + labs(x = "Tempos", y = "S(t)") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_bw() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))

# gráfico alisado
ggplot(data = ekm_df, mapping = aes(time, surv)) + 
    geom_line(aes(color = consulta)) + labs(x = "Tempos", y = "S(t)") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_bw() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))