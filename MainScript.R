

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(survival)


# Load and Treatment Data -------------------------------------------------

data <- read.csv(file = "recuperaInadimplencia.csv")
data <- data[,2:4]



# Descriptive Analysis ------------------------------------------------------

# response var histogram Time until debt payment
# histograma da var resposta Tempo até o Pagamento da dívida
hist(data$tempo_pag, xlab="Tempo Pagamento", 
     main = "Histograma do Tempo de Pagamento")

# estimated Kaplan-Meier curve
# curva de Kaplan-Meier estimada
ekm <- survfit(Surv(time = tempo_pag, event = status) ~ 1, data = data)
ekm_df <- data.frame(time = ekm$time,
                     surv = ekm$surv)

ggplot(data = ekm_df, mapping = aes(time, surv)) + 
    geom_line() + labs(x = "Tempos", y = "S(t)") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    geom_hline(yintercept = 0.76, linetype="dashed", color = "red") +
    geom_hline(yintercept = 0.215, linetype="dashed", color = "red") +
    geom_text(label = "1-p0", x = 0, y = 0.81, size = 8/.pt) +
    geom_text(label = "p1", x = 0, y = 0.19, size = 8/.pt) +
    theme_bw()

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo_pag, event = status) ~ consulta, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     consulta = c(rep(0,44),rep(1,61)) )

ggplot(data = ekm_df, mapping = aes(tempo, km)) + 
    geom_step(aes(color = as.factor(consulta))) + 
    labs(x = "Tempos", y = "S(t)", color = "Consulta") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_bw() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))

# smoothed graphic
# gráfico alisado
ggplot(data = ekm_df, mapping = aes(tempo, km)) + 
    geom_line(aes(color = as.factor(consulta))) + 
    labs(x = "Tempos", y = "S(t)", color = "Consulta") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_bw() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))


# Modeling ----------------------------------------------------------------

# verisimilitude of GTDL Gamma
# verossimilhança da GTDL Gamma
veroGTDL <- function(x, par) {
    
    # parameters
    alph <- par[1]
    lambd <- exp(par[2])
    thet <- exp(par[3])
    bet <- par[4:length(par)]
    # covariates
    X <- as.matrix(x[x$tempo > 0, 3:ncol(x)])
    X0 <- as.matrix(x[x$tempo == 0, 3:ncol(x)])
    cens <- x[x$tempo > 0,1]
    tempo <- x[x$tempo > 0,2]
    # model
    aux1 <- sum( (X0%*%bet) - log( 1 + exp(X0%*%bet) ) )
    aux2 <- log(lambd)*sum(cens) - 
        sum( log( 1 + exp(X%*%bet) ) ) +
        sum( cens*(alph*tempo + X%*%bet)) -
        sum( cens*log( 1 + exp(alph*tempo + X%*%bet) ) ) - 
        sum( (cens + (1/thet))*log( 1 + ((thet*lambd/alph)*log( (1 + exp(alph*tempo + X%*%bet))/(1 + exp(X%*%bet)) )) ) )
    
    aux <- aux1 + aux2
    return(-aux)
}

# GTDL Zero Gamma Inflated Survival Model
# modelo de sobrevivência da GTDL Gama zero inflacionada
sobrevGTDL <- function(x, par) {
    
    alph <- par[1]
    lambd <- par[2]
    thet <- par[3]
    bet <- par[4:length(par)]
    tempo <- x$tempo
    X <- as.matrix(x[, 2:ncol(x)])
    st <- (1/(1 + exp(X%*%bet) )) * (1 + ( (lambd*thet/alph) * log( (1 + exp(alph*tempo + X%*%bet) )/(1 + exp(X%*%bet))) ) )^(-1/thet)
    
    return(st)
}



# Aplication --------------------------------------------------------------

# parameter estimation
# estimação dos parâmetros
emv <- optim(par=c(0.1,0.3,0.5,2), veroGTDL, x=data)
epar <- c(emv$par[1], exp(emv$par[2]), exp(emv$par[3]), emv$par[4])

ekm_df$gtdl <- sobrevGTDL(x = ekm_df[,c(1,3)], par = epar)

# Graph of survival curves by Kaplan Meyer (KM) and GTDL
# Gráfico das curvas de sobrevivência por Kaplan Meyer (KM) e por GTDL
ggplot(data = ekm_df, mapping = aes(tempo, gtdl)) +
    geom_line(aes(linetype = as.factor(consulta))) +
    geom_line(aes(x = tempo, y = km, color = as.factor(consulta))) +
    labs(x = "Tempos", y = "S(t)", color = "Consulta por KM", linetype = "Consulta por GTDL") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_bw() + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(color = "white"))

# Testing
# X <- data.frame(tempo = data$tempo_pag,
#                 consulta = as.factor(data$consulta),
#                 gtdl = sobrevGTDL(x = data[,2:3], par = epar))
