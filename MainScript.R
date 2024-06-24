rm(list = ls())

# Libraries ---------------------------------------------------------------

library(ggplot2)
library(survival)
library(msm) # pacote pra usar o método delta

# Load and Treatment Data -------------------------------------------------

dados <- read.table(file = "dados_censura_versao4.txt", header = TRUE)


# Descriptive Analysis ------------------------------------------------------

# barplot of time to fail 
aux <- data.frame(tempo = dados$tempo[dados$tempo < 25])
ggplot(aux, aes(x = tempo)) + 
    geom_bar() + ylim(0,2400)


# Kaplan-Meier estimate graphs

# without covar
resultado <- survfit(Surv(time = tempo, event = status)~1, data = dados)

plot(resultado, xlim=c(0, 60), ylim=c(0, 1),
     xlab = " Tempo (meses)", ylab = "Estimativa de S(t)", bty = "n")

# with covar type of debt
resultado <- survfit(Surv(time = tempo, event = status)~factor(x4), data = dados)

plot(resultado, xlim=c(0, 60), ylim=c(0, 1),
     xlab = " Tempo (meses)", ylab = "Estimativa de S(t)",
     lty = c(2, 2), col = c("black", "blue"),  bty="n")
legend(7, 1,
       title="Tipo de dívida", c("0 - Banco", "1 - outros segmentos"),
       lty = c(2, 2), col = c("black", "blue"), bty = "n")


# estimated Kaplan-Meier curve
# curva de Kaplan-Meier estimada
ekm <- survfit(Surv(time = tempo, event = status) ~ 1, data = dados)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv)

# graph with limits
ggplot(data = ekm_df, mapping = aes(tempo, km)) + 
    geom_line() + labs(x = "Tempos", y = "S(t)") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    geom_hline(yintercept = 0.76, linetype="dashed", color = "red") +
    geom_hline(yintercept = 0.215, linetype="dashed", color = "red") +
    geom_text(label = "1-p0", x = 0, y = 0.81, size = 8/.pt) +
    geom_text(label = "p1", x = 0, y = 0.19, size = 8/.pt) +
    theme_classic()

# graph of Kaplan-Meyer adjust
ggplot(data = ekm_df, mapping = aes(tempo, km)) + 
    geom_step() + 
    labs(x = "Tempo (meses)", y = "Estimativa de S(t)") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))


# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com COVARIÁVEL DE CONSULTA
ekm <- survfit(Surv(time = tempo, event = status) ~ x1, data = dados)
# ekm_df <- data.frame(tempo = c(0,0, ekm$time), km = c(1,1, ekm$surv), 
#                      covar = c(0,1, rep(0,ekm$strata[1]),rep(1,ekm$strata[2])))
ekm_df <- data.frame(tempo = ekm$time, km = ekm$surv, 
                     covar = c(rep(0,ekm$strata[1]),rep(1,ekm$strata[2])))


ggplot(data = ekm_df, mapping = aes(tempo, km)) + 
    geom_step(aes(color = as.factor(covar))) + 
    labs(x = "Tempo (meses)", y = "Estimativa de S(t)", 
         color = "Informação de Consulta") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                       values = c("red","blue")) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_light() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))

# smoothed graphic
# gráfico alisado
ggplot(data = ekm_df, mapping = aes(tempo, km)) + 
    geom_line(aes(color = as.factor(covar))) + 
    labs(x = "Tempo (meses)", y = "Estimativa de S(t)", 
         color = "Informação de Consulta") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                       values = c("red","blue")) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_light() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))

LogRankX1 = survdiff(formula = Surv(tempo, status)~factor(x1), data = dados, rho=0)
LogRankX1

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com COVARIÁVEL DE TIPO DE DIVIDA
ekm <- survfit(Surv(time = tempo, event = status) ~ x4, data = dados)
# ekm_df <- data.frame(tempo = c(0,0, ekm$time), km = c(1,1, ekm$surv), 
#                      covar = c(0,1, rep(0,ekm$strata[1]),rep(1,ekm$strata[2])))
ekm_df <- data.frame(tempo = ekm$time, km = ekm$surv, 
                     covar = c(rep(0,ekm$strata[1]),rep(1,ekm$strata[2])))

ggplot(data = ekm_df, mapping = aes(tempo, km)) + 
    geom_step(aes(color = as.factor(covar))) + 
    labs(x = "Tempo (meses)", y = "Estimativa de S(t)", 
         color = "Tipo de Dívida") +
    scale_color_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                       values = c("red","blue")) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_light() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))

# smoothed graphic
# gráfico alisado
ggplot(data = ekm_df, mapping = aes(tempo, km)) + 
    geom_line(aes(color = as.factor(covar))) + 
    labs(x = "Tempo (meses)", y = "Estimativa de S(t)",
         color = "Tipo de Dívida") +
    scale_color_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                       values = c("red","blue")) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))

LogRankX4 = survdiff(formula = Surv(tempo, status)~factor(x4), data = dados, rho=0)
LogRankX4

#ggpubr::ggarrange(g1,g2)


# Model 1 GTDL Gamma ------------------------------------------------------

# likelihoof of GTDL Gamma
veroGTDLg <- function(x, par) {
    
    # parameters
    alph <- par[1]
    lambd <- exp(par[2])
    thet <- exp(par[3])
    bet <- par[4:length(par)]
    
    # covariates
    X <- as.matrix(x[, 3:ncol(x)])
    cens <- x[,1]
    tempo <- x[,2]
    
    # model
    aux2 <- log(lambd)*sum(cens) +
        sum( cens*(X%*%bet + alph*tempo)) -
        sum( cens*log( 1 + exp(X%*%bet + alph*tempo)) ) - 
        sum( (cens + (1/thet))*log( 1 + ((thet*lambd/alph)*log( (1 + exp(X%*%bet + alph*tempo))/(1 + exp(X%*%bet)) )) ) )
    
    return(-aux2)
}

# GTDL Gamma Survival Model
sobrevGTDLg <- function(x, par) {
    
    alph <- par[1]
    lambd <- par[2]
    thet <- par[3]
    bet <- par[4:length(par)]
    tempo <- x[,1]
    X <- as.matrix(x[, 2:ncol(x)])
    st <- (1 + ( (lambd*thet/alph)*log( (1+ exp(X%*%bet + alph*tempo))/(1 + exp(X%*%bet))) ) )^(-1/thet)
    return(st)
}



# Aplication Model 1 ------------------------------------------------------


# QUERY VARIABLE
# VARIÁVEL CONSULTA
data <- dados[dados$tempo > 0,1:3]
#data$tempo[data$tempo == 0] <- 10^-3

# parameter estimation of GTDL Gamma
emvg <- optim(par=c(0.02,0.5,0.5,-1), veroGTDLg, x=data, hessian = TRUE);emvg
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo, event = status) ~ x1, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     covar = c(rep(0,ekm$strata[1]),rep(1,ekm$strata[2])) )

# GTDL survival curve estimation
# estimação da curva de sobrevivência GTDL
gtdl_df <- data.frame(tempo = rep(seq(1,60,1),2),
                      covar = c(rep(0,60),rep(1,60)))
gtdl_df$surv <- sobrevGTDLg(x = gtdl_df, par = egpar)


# Graph of survival curves by Kaplan Meyer (KM) and GTDL
# Gráfico das curvas de sobrevivência por Kaplan Meyer (KM) e por GTDL 
ggplot(data = ekm_df, mapping = aes(x = tempo, y = km)) +
    geom_step(aes(color = as.factor(covar))) +
    geom_line(data = gtdl_df, aes(x = tempo, y = surv, linetype = as.factor(covar))) +
    labs(x = "Tempo (meses)", y = "Estimativa de S(t)",
         color = "Consulta por Kaplan-Meier", linetype = "Consulta por GTDL Gama") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                       values = c("red","blue")) +
    scale_linetype_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_light() + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(color = "white"))


# estimation of the confidence interval of the parameters
# estimação do intervalo de confiança dos parâmetros
sd <- sqrt(diag(solve(emvg$hessian)))

ic_alpha <- numeric()
ic_alpha[1] <- egpar[1] -qnorm(0.975)*sd[1]
ic_alpha[2] <- egpar[1] +qnorm(0.975)*sd[1]
round(ic_alpha,3);egpar[1];sd[1]

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = solve(emvg$hessian)[2,2])
ic_lambda <- numeric()
ic_lambda[1] <- egpar[2] -qnorm(0.975)*sd
ic_lambda[2] <- egpar[2] +qnorm(0.975)*sd
round(ic_lambda,3);egpar[2];sd

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_theta <- numeric()
ic_theta[1] <- egpar[3] -qnorm(0.975)*sd
ic_theta[2] <- egpar[3] +qnorm(0.975)*sd
round(ic_theta,3);egpar[3];sd

sd <- sqrt(diag(solve(emvg$hessian)))
ic_beta1g <- numeric()
ic_beta1g[1] <- egpar[4] -qnorm(0.975)*sd[4]
ic_beta1g[2] <- egpar[4] +qnorm(0.975)*sd[4]
round(ic_beta1g,3);egpar[4];sd[4]


# VARIABLE TYPE OF DEBT
# VARIÁVEL TIPO DE DÍVIDA

data <- dados[dados$tempo > 0,c(1,2,6)]
#data$tempo[data$tempo == 0] <- 10^-3

emvg <- optim(par=c(0.02, 0.5, 0.5, -1), veroGTDLg, x=data, hessian = TRUE);emvg
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])

ekm <- survfit(Surv(time = tempo, event = status) ~ x4, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     covar = c(rep(0,ekm$strata[1]),rep(1,ekm$strata[2])) )

gtdl_df <- data.frame(tempo = rep(seq(1,60,1),2),
                      covar = c(rep(0,60),rep(1,60)))
gtdl_df$surv <- sobrevGTDLg(x = gtdl_df, par = egpar)


ggplot(data = ekm_df, mapping = aes(x = tempo, y = km)) +
    geom_step(aes(color = as.factor(covar))) +
    geom_line(data = gtdl_df, aes(x = tempo, y = surv, linetype = as.factor(covar))) +
    labs(x="Tempo (meses)", y="Estimativa de S(t)", 
         color="Tipo de Dívida por Kaplan-Meier",
         linetype="Tipo de Dívida por GTDL Gama") +
    scale_color_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                       values = c("red","blue")) +
    scale_linetype_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_light() + 
    theme(legend.position = c(0.75,0.7),
          legend.background = element_rect(color = "white"))

ggpubr::ggarrange(g1,g2)


# estimation of the confidence interval of the parameters
# estimação do intervalo de confiança dos parâmetros
sd <- sqrt(diag(solve(emvg$hessian)))

ic_alpha <- numeric()
ic_alpha[1] <- egpar[1] -qnorm(0.975)*sd[1]
ic_alpha[2] <- egpar[1] +qnorm(0.975)*sd[1]
round(ic_alpha,3);egpar[1];sd[1]

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = solve(emvg$hessian)[2,2])
ic_lambda <- numeric()
ic_lambda[1] <- egpar[2] -qnorm(0.975)*sd
ic_lambda[2] <- egpar[2] +qnorm(0.975)*sd
round(ic_lambda,3);egpar[2];sd

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_theta <- numeric()
ic_theta[1] <- egpar[3] -qnorm(0.975)*sd
ic_theta[2] <- egpar[3] +qnorm(0.975)*sd
round(ic_theta,3);egpar[3];sd

sd <- sqrt(diag(solve(emvg$hessian)))
ic_beta1g <- numeric()
ic_beta1g[1] <- egpar[4] -qnorm(0.975)*sd[4]
ic_beta1g[2] <- egpar[4] +qnorm(0.975)*sd[4]
round(ic_beta1g,3);egpar[4];sd[4]



# Model 2 GTDL p0 ---------------------------------------------------------


# likelihoof of GTDL Gamma p0
veroGTDLp0 <- function(x, par) {
    
    # parameters
    p0 <- par[1]
    alph <- par[2]
    lambd <- exp(par[3])
    thet <- exp(par[4])
    bet <- par[5:length(par)]
    n0 <- nrow(x[x$tempo == 0,])
    n1 <- nrow(x[x$tempo > 0,])
    
    # covariates
    X <- as.matrix(x[x$tempo > 0, 3:ncol(x)])
    cens <- x[x$tempo > 0,1]
    tempo <- x[x$tempo > 0,2]
    
    # model
    aux2 <- n0*log(p0) + n1*log(1-p0) + log(lambd)*sum(cens) +
        sum( cens*(alph*tempo + X%*%bet)) -
        sum( cens*log( 1 + exp(alph*tempo + X%*%bet) ) ) - 
        sum( (cens + (1/thet))*log( 1 + ((thet*lambd/alph)*log( (1 + exp(alph*tempo + X%*%bet))/(1 + exp(X%*%bet)) )) ) )
    
    return(-aux2)
}

# GTDL Gamma p0 Survival Model
sobrevGTDLp0 <- function(x, par) {
    
    p0 <- par[1]
    alph <- par[2]
    lambd <- par[3]
    thet <- par[4]
    bet <- par[5:length(par)]
    tempo <- x$tempo
    X <- as.matrix(x[, 2:ncol(x)])
    st <- (1-p0) * (1 + ( (lambd*thet/alph) * log( (1 + exp(alph*tempo + X%*%bet) )/(1 + exp(X%*%bet))) ) )^(-1/thet)
    return(st)
}



# Aplication Model 2 ------------------------------------------------------

# QUERY VARIABLE
# VARIÁVEL CONSULTA

data <- dados[,1:3]

# parameter estimation of GTDL Gamma p0
emvp0 <- optim(par=c(0.5,0.1,-0.5,-0.5,0.1), veroGTDLp0, x=data, hessian = TRUE);emvp0
ep0par <- c(emvp0$par[1], emvp0$par[2], exp(emvp0$par[3]), exp(emvp0$par[4]), emvp0$par[5])

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo, event = status) ~ x1, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     consulta = c(rep(0,44),rep(1,61)) )

# GTDL survival curve estimation
# estimação da curva de sobrevivência GTDL
ekm_df$gtdlp0 <- sobrevGTDLp0(x = ekm_df[,c(1,3)], par = ep0par)

# Graph of survival curves by Kaplan Meyer (KM) and GTDL
# Gráfico das curvas de sobrevivência por Kaplan Meyer (KM) e por GTDL
g1 <- ggplot(data = ekm_df, mapping = aes(tempo, gtdlp0[,1])) +
    geom_line(aes(linetype = as.factor(consulta))) +
    geom_line(aes(x = tempo, y = km, color = as.factor(consulta))) +
    labs(x = "Tempo (meses)", y = "Estimativa de S(t)",
         color = "Consulta por Kaplan-Meier",
         linetype = "Consulta por GTDL Gama p0") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                       values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(color = "white"))


# estimation of the confidence interval of the parameters
# estimação do intervalo de confiança dos parâmetros
sd <- sqrt(diag(solve(emvp0$hessian)))

ic_p0 <- numeric()
ic_p0[1] <- ep0par[1] -qnorm(0.975)*sd[1]
ic_p0[2] <- ep0par[1] +qnorm(0.975)*sd[1]
round(ic_p0,3);round(ep0par[1],3);round(sd[1],3)

ic_alpha <- numeric()
ic_alpha[1] <- ep0par[2] -qnorm(0.975)*sd[2]
ic_alpha[2] <- ep0par[2] +qnorm(0.975)*sd[2]
round(ic_alpha,3);round(ep0par[2],3);round(sd[2],3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvp0$par[2], cov = solve(emvp0$hessian)[2,2])
ic_lambda <- numeric()
ic_lambda[1] <- ep0par[3] -qnorm(0.975)*sd
ic_lambda[2] <- ep0par[3] +qnorm(0.975)*sd
round(ic_lambda,3);round(ep0par[3],3);round(sd,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvp0$par[3], cov = solve(emvp0$hessian)[3,3])
ic_theta <- numeric()
ic_theta[1] <- ep0par[4] -qnorm(0.975)*sd
ic_theta[2] <- ep0par[4] +qnorm(0.975)*sd
round(ic_theta,3);round(ep0par[4],3);round(sd,3)

sd <- sqrt(diag(solve(emvp0$hessian)))
ic_beta1g <- numeric()
ic_beta1g[1] <- ep0par[5] -qnorm(0.975)*sd[5]
ic_beta1g[2] <- ep0par[5] +qnorm(0.975)*sd[5]
round(ic_beta1g,3);round(ep0par[5],3);round(sd[5],3)



# VARIABLE TYPE OF DEBT
# VARIÁVEL TIPO DE DÍVIDA

data <- dados[,c(1,2,6)]

emvp0 <- optim(par=c(0.25,0.1,-0.5,-0.5,0.1), veroGTDLp0, x=data, hessian = TRUE);emvp0
ep0par <- c(emvp0$par[1], emvp0$par[2], exp(emvp0$par[3]), exp(emvp0$par[4]), emvp0$par[5])

ekm <- survfit(Surv(tempo, status) ~ x4, data = data)
ekm_df <- data.frame(tempo = ekm$time, km = ekm$surv, 
                     covar = c(rep(0,ekm$strata[[1]]),rep(1,ekm$strata[[2]])))

ekm_df$gtdlp0 <- sobrevGTDLp0(x = ekm_df[,c(1,3)], par = ep0par)

g2 <- ggplot(data = ekm_df, mapping = aes(x = tempo, y = gtdlp0)) +
    geom_line(aes(linetype = as.factor(covar))) +
    geom_line(aes(x=tempo, y=km, color=factor(covar))) +
    labs(x="Tempo (meses)", y="Estimativa de S(t)", 
         color="Tipo de Dívida por Kaplan-Meier",
         linetype="Tipo de Dívida por GTDL Gama p0") +
    scale_color_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                       values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_classic() + 
    theme(legend.position = c(0.75,0.7),
          legend.background = element_rect(color = "white"))

ggpubr::ggarrange(g1,g2)


# estimation of the confidence interval of the parameters
# estimação do intervalo de confiança dos parâmetros
sd <- sqrt(diag(solve(emvp0$hessian)))

ic_p0 <- numeric()
ic_p0[1] <- ep0par[1] -qnorm(0.975)*sd[1]
ic_p0[2] <- ep0par[1] +qnorm(0.975)*sd[1]
round(ic_p0,3);round(ep0par[1],3);round(sd[1],3)

ic_alpha <- numeric()
ic_alpha[1] <- ep0par[2] -qnorm(0.975)*sd[2]
ic_alpha[2] <- ep0par[2] +qnorm(0.975)*sd[2]
round(ic_alpha,3);round(ep0par[2],3);round(sd[2],3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvp0$par[2], cov = solve(emvp0$hessian)[2,2])
ic_lambda <- numeric()
ic_lambda[1] <- ep0par[3] -qnorm(0.975)*sd
ic_lambda[2] <- ep0par[3] +qnorm(0.975)*sd
round(ic_lambda,3);round(ep0par[3],3);round(sd,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvp0$par[3], cov = solve(emvp0$hessian)[3,3])
ic_theta <- numeric()
ic_theta[1] <- ep0par[4] -qnorm(0.975)*sd
ic_theta[2] <- ep0par[4] +qnorm(0.975)*sd
round(ic_theta,3);round(ep0par[4],3);round(sd,3)

sd <- sqrt(diag(solve(emvp0$hessian)))
ic_beta1g <- numeric()
ic_beta1g[1] <- ep0par[5] -qnorm(0.975)*sd[5]
ic_beta1g[2] <- ep0par[5] +qnorm(0.975)*sd[5]
round(ic_beta1g,3);round(ep0par[5],3);round(sd[5],3)



# Model 3 GTDL Zero v1 ----------------------------------------------------

# a versão 1 do GTDL Zero Inflacionada utiliza duas verossimilhanças
# uma é para a inflação de zeros e a outra para o GTDL, além disso a v1
# obtém apenas um valor para o parâmetro alpha
# que é o parâmetro responsável pelo ajuste da fração de cura


# likelihood of zero inflation
# verossimilhança da inflação de zeros
veroZerov1 <- function(x, beta) {
    
    # covariates
    X0 <- as.matrix(x[x$tempo == 0, 3:ncol(x)])
    X <- as.matrix(x[x$tempo > 0, 3:ncol(x)])
    
    # model
    aux1 <- sum((X0%*%beta)) - sum(log( 1 + exp(X0%*%beta))) - sum(log(1 + exp(X%*%beta)))
    return(-aux1)
}

# likelihood of GTDL Gamma
# verossimilhança da GTDL Gamma
veroGTDLv1 <- function(x, par) {
    
    # parameters
    alph <- par[1]
    lambd <- exp(par[2])
    thet <- exp(par[3])
    bet <- par[4:length(par)]
    
    # covariates
    X <- as.matrix(x[x$tempo > 0, 3:ncol(x)])
    cens <- x[x$tempo > 0,1]
    tempo <- x[x$tempo > 0,2]
    
    # model
    aux2 <- log(lambd)*sum(cens) +
        sum( cens*(alph*tempo + X%*%bet)) -
        sum( cens*log( 1 + exp(alph*tempo + X%*%bet) ) ) - 
        sum( (cens + (1/thet))*log( 1 + ((thet*lambd/alph)*log( (1 + exp(alph*tempo + X%*%bet))/(1 + exp(X%*%bet)) )) ) )
    
    return(-aux2)
}

# GTDL Zero Gamma Inflated Survival Model
# modelo de sobrevivência da GTDL Gama zero inflacionada
sobrevGTDLv1 <- function(x, par) {
    
    bet0 <- par[1:2]
    alph <- par[3]
    lambd <- par[4]
    thet <- par[5]
    bet <- par[6:length(par)]
    tempo <- x$tempo
    X <- as.matrix(x[, 2:ncol(x)])
    X0 <- matrix(data = c(rep(1,nrow(x)),X), nrow = nrow(x), ncol = ncol(x))
    st <- (1/(1 + exp(X0%*%bet0) )) * (1 + ( (lambd*thet/alph) * log( (1 + exp(alph*tempo + X%*%bet) )/(1 + exp(X%*%bet))) ) )^(-1/thet)
    
    return(st)
}



# Aplication Model 3 --------------------------------------------------------------


# QUERY VARIABLE
# VARIÁVEL CONSULTA

data <- dados[,1:3]
data0$covar <- data[,3]

# parameter estimation of GTDL Gamma
# estimação dos parâmetros do GTDL Gamma
emvg <- optim(par=c(0.1, -0.5, -0.5, 0.1), veroGTDLv1, x=data, hessian = TRUE);emvg
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])

# parameter estimation of zero inflation
# estimação dos parâmetros da inflação de zeros
emv0 <- optim(par=c(-0.5, 0.5), veroZerov1, x=data0, hessian = TRUE);emv0
e0par <- c(emv0$par[1], emv0$par[2])

# combinando os parâmetros estimados
epar <- c(e0par, egpar); epar

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo, event = status) ~ x1, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     consulta = c(rep(0,44),rep(1,61)) )

# GTDL survival curve estimation
# estimação da curva de sobrevivência GTDL
ekm_df$gtdl <- sobrevGTDLv1(x = ekm_df[,c(1,3)], par = epar)

# Graph of survival curves by Kaplan Meyer (KM) and GTDL
# Gráfico das curvas de sobrevivência por Kaplan Meyer (KM) e por GTDL
g1 <- ggplot(data = ekm_df, mapping = aes(tempo, gtdl)) +
    geom_line(aes(linetype = as.factor(consulta))) +
    geom_line(aes(x = tempo, y = km, color = as.factor(consulta))) +
    labs(x = "Tempo (meses)", y = "Estimativa de S(t)", 
         color = "Consulta por Kaplan-Meier", 
         linetype = "Consulta por GTDL \nZero Inflacionado v1") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), 
                       values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), 
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(color = "white"))


# estimation of the confidence interval of the parameters
# estimação do intervalo de confiança dos parâmetros
sd0 <- sqrt(diag(solve(emv0$hessian)))
sdg <- sqrt(diag(solve(emvg$hessian)))
sd <- c(sd0,sdg); sd

ic_beta00 <- numeric()
ic_beta00[1] <- epar[1] -qnorm(0.975)*sd[1]
ic_beta00[2] <- epar[1] +qnorm(0.975)*sd[1]
round(ic_beta00,3);round(epar[1],3);round(sd[1],3)

ic_beta10 <- numeric()
ic_beta10[1] <- epar[2] -qnorm(0.975)*sd[2]
ic_beta10[2] <- epar[2] +qnorm(0.975)*sd[2]
round(ic_beta10,3);round(epar[2],3);round(sd[2],3)

ic_alpha <- numeric()
ic_alpha[1] <- epar[3] -qnorm(0.975)*sd[3]
ic_alpha[2] <- epar[3] +qnorm(0.975)*sd[3]
round(ic_alpha,3);round(epar[3],3);round(sd[3],3)


# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = solve(emvg$hessian)[2,2])
ic_lambda <- numeric()
ic_lambda[1] <- epar[4] -qnorm(0.975)*sd
ic_lambda[2] <- epar[4] +qnorm(0.975)*sd
round(ic_lambda,3);round(epar[4],3);round(sd,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_theta <- numeric()
ic_theta[1] <- epar[5] -qnorm(0.975)*sd
ic_theta[2] <- epar[5] +qnorm(0.975)*sd
round(ic_theta,3);round(epar[5],3);round(sd,3)

sd <- c(sd0,sdg)
ic_beta1g <- numeric()
ic_beta1g[1] <- epar[6] -qnorm(0.975)*sd[6]
ic_beta1g[2] <- epar[6] +qnorm(0.975)*sd[6]
round(ic_beta1g,3);round(epar[6],3);round(sd[6],3)



# VARIABLE TYPE OF DEBT
# VARIÁVEL TIPO DE DÍVIDA

data <- dados[,c(1,2,6)]
data0$covar <- data[,3]

emvg <- optim(par=c(0.1,-0.5,-0.5,0.1), veroGTDLv1, x=data, hessian = TRUE);emvg
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])

emv0 <- optim(par=c(-0.5,0.5), veroZerov1, x=data0, hessian = TRUE);emv0
e0par <- c(emv0$par[1], emv0$par[2])

epar <- c(e0par, egpar); epar

ekm <- survfit(Surv(tempo, status) ~ x4, data = data)
ekm_df <- data.frame(tempo = ekm$time, km = ekm$surv, 
                     covar = c(rep(0,ekm$strata[[1]]),rep(1,ekm$strata[[2]])))

ekm_df$gtdl0 <- sobrevGTDLv1(ekm_df[,c(1,3)],epar)

g2 <- ggplot(data = ekm_df, mapping = aes(x = tempo, y = gtdl0)) +
    geom_line(aes(linetype = as.factor(covar))) +
    geom_line(aes(x=tempo, y=km, color=factor(covar))) +
    labs(x="Tempo (meses)", y="Estimativa de S(t)", 
         color="Tipo de Dívida por Kaplan-Meier",
         linetype="Tipo de Dívida por GTDL \nZero Inflacionado v1") +
    scale_color_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"), 
                       values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"), 
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_classic() + 
    theme(legend.position = c(0.75,0.7),
          legend.background = element_rect(color = "white"))

ggpubr::ggarrange(g1,g2)

sd0 <- sqrt(diag(solve(emv0$hessian)))
sdg <- sqrt(diag(solve(emvg$hessian)))
sd <- c(sd0,sdg); sd

ic_beta00 <- numeric()
ic_beta00[1] <- epar[1] -qnorm(0.975)*sd[1]
ic_beta00[2] <- epar[1] +qnorm(0.975)*sd[1]
round(ic_beta00,3);round(epar[1],3);round(sd[1],3)

ic_beta10 <- numeric()
ic_beta10[1] <- epar[2] -qnorm(0.975)*sd[2]
ic_beta10[2] <- epar[2] +qnorm(0.975)*sd[2]
round(ic_beta10,3);round(epar[2],3);round(sd[2],3)

ic_alpha <- numeric()
ic_alpha[1] <- epar[3] -qnorm(0.975)*sd[3]
ic_alpha[2] <- epar[3] +qnorm(0.975)*sd[3]
round(ic_alpha,3);round(epar[3],3);round(sd[3],3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = solve(emvg$hessian)[2,2])
ic_lambda <- numeric()
ic_lambda[1] <- epar[4] -qnorm(0.975)*sd
ic_lambda[2] <- epar[4] +qnorm(0.975)*sd
round(ic_lambda,3);round(epar[4],3);round(sd,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_theta <- numeric()
ic_theta[1] <- epar[5] -qnorm(0.975)*sd
ic_theta[2] <- epar[5] +qnorm(0.975)*sd
round(ic_theta,3);round(epar[5],3);round(sd,3)

sd <- c(sd0,sdg)
ic_beta1g <- numeric()
ic_beta1g[1] <- epar[6] -qnorm(0.975)*sd[6]
ic_beta1g[2] <- epar[6] +qnorm(0.975)*sd[6]
round(ic_beta1g,3);round(epar[6],3);round(sd[6],3)



# Model 4 GTDL Zero v2 ----------------------------------------------------------------

# a versão 2 do GTDL Zero Inflacionada faz regressão no parâmetro alpha
# que é o parâmetro responsável pelo ajuste da fração de cura

# likelihood of zero inflation
# verossimilhança da inflação de zeros
veroZerov2 <- function(x, beta) {
    
    # covariates
    X0 <- as.matrix(x[x$tempo == 0, 3:ncol(x)])
    X <- as.matrix(x[x$tempo > 0, 3:ncol(x)])
     
    # model
    aux1 <- sum((X0%*%beta)) - sum(log( 1 + exp(X0%*%beta))) - sum(log(1 + exp(X%*%beta)))
    return(-aux1)
}

# likelihood of GTDL Gamma
# verossimilhança da GTDL Gamma
veroGTDLv2 <- function(x, par) {
    
    # parameters
    alph <- par[1:2]
    lambd <- exp(par[3])
    thet <- exp(par[4])
    bet <- par[5:length(par)]
    
    # covariates
    X <- as.matrix(x[x$tempo > 0, 3:ncol(x)])
    xalpha <- matrix(c(rep(1,nrow(X)),X), ncol = ncol(X)+1)
    cens <- x[x$tempo > 0,1]
    tempo <- x[x$tempo > 0,2]
    
    # model
    aux2 <- log(lambd)*sum(cens) +
        sum( cens*((xalpha%*%alph)*tempo + X%*%bet)) -
        sum( cens*log( 1 + exp((xalpha%*%alph)*tempo + X%*%bet) ) ) - 
        sum( (cens + (1/thet))*log( 1 + ((thet*lambd/(xalpha%*%alph))*log( (1 + exp((xalpha%*%alph)*tempo + X%*%bet))/(1 + exp(X%*%bet)) )) ) )

    return(-aux2)
}

# GTDL Zero Gamma Inflated Survival Model
# modelo de sobrevivência da GTDL Gama zero inflacionada
sobrevGTDLv2 <- function(x, par) {
    
    bet0 <- par[1:2]
    alph <- par[3:4]
    lambd <- par[5]
    thet <- par[6]
    bet <- par[7:length(par)]
    tempo <- x$tempo
    X <- as.matrix(x[, 2:ncol(x)])
    X0 <- matrix(data = c(rep(1,nrow(x)),X), nrow = nrow(x), ncol = ncol(x))
    xalpha <- matrix(c(rep(1,nrow(X)),X), ncol = ncol(X)+1)
    st <- (1/(1 + exp(X0%*%bet0) )) * (1 + ( (lambd*thet/(xalpha%*%alph)) * log( (1 + exp((xalpha%*%alph)*tempo + X%*%bet) )/(1 + exp(X%*%bet))) ) )^(-1/thet)
    
    return(st)
}



# Aplication Model 4 --------------------------------------------------------------


# QUERY VARIABLE
# VARIÁVEL CONSULTA

data <- dados[,1:3]
data0$covar <- data[,3]

# parameter estimation of GTDL Gamma
# estimação dos parâmetros do GTDL Gamma
emvg <- optim(par=c(0.1, 0.5, -0.5, -0.5, 0.1), veroGTDLv2, x=data, hessian = TRUE);emvg
egpar <- c(emvg$par[1], emvg$par[2], exp(emvg$par[3]), exp(emvg$par[4]), emvg$par[5])

# parameter estimation of zero inflation
# estimação dos parâmetros da inflação de zeros
emv0 <- optim(par=c(2, 1), veroZerov2, x=data0, hessian = TRUE);emv0
e0par <- c(emv0$par[1], emv0$par[2])

# combinando os parâmetros estimados
epar <- c(e0par, egpar); epar

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo, event = status) ~ x1, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     consulta = c(rep(0,44),rep(1,61)) )

# GTDL survival curve estimation
# estimação da curva de sobrevivência GTDL
ekm_df$gtdl <- sobrevGTDLv2(x = ekm_df[,c(1,3)], par = epar)

# Graph of survival curves by Kaplan Meyer (KM) and GTDL
# Gráfico das curvas de sobrevivência por Kaplan Meyer (KM) e por GTDL
g1 <- ggplot(data = ekm_df, mapping = aes(tempo, gtdl)) +
    geom_line(aes(linetype = as.factor(consulta))) +
    geom_line(aes(x = tempo, y = km, color = as.factor(consulta))) +
    labs(x = "Tempo (meses)", y = "Estimativa de S(t)",
         color = "Consulta por Kaplan-Meier",
         linetype = "Consulta por GTDL \nZero Inflacionado v2") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                       values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(color = "white"))


# estimation of the confidence interval of the parameters
# estimação do intervalo de confiança dos parâmetros
sd0 <- sqrt(diag(solve(emv0$hessian)))
sdg <- sqrt(diag(solve(emvg$hessian)))
sd <- c(sd0,sdg); sd

ic_beta00 <- numeric()
ic_beta00[1] <- epar[1] -qnorm(0.975)*sd[1]
ic_beta00[2] <- epar[1] +qnorm(0.975)*sd[1]
round(ic_beta00,3);round(epar[1],3);round(sd[1],3)

ic_beta10 <- numeric()
ic_beta10[1] <- epar[2] -qnorm(0.975)*sd[2]
ic_beta10[2] <- epar[2] +qnorm(0.975)*sd[2]
round(ic_beta10,3);round(epar[2],3);round(sd[2],3)

ic_alpha0 <- numeric()
ic_alpha0[1] <- epar[3] -qnorm(0.975)*sd[3]
ic_alpha0[2] <- epar[3] +qnorm(0.975)*sd[3]
round(ic_alpha0,3);round(epar[3],3);round(sd[3],3)

ic_alpha1 <- numeric()
ic_alpha1[1] <- epar[4] -qnorm(0.975)*sd[4]
ic_alpha1[2] <- epar[4] +qnorm(0.975)*sd[4]
round(ic_alpha1,3);round(epar[4],3);round(sd[4],3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_lambda <- numeric()
ic_lambda[1] <- epar[5] -qnorm(0.975)*sd
ic_lambda[2] <- epar[5] +qnorm(0.975)*sd
round(ic_lambda,3);round(epar[5],3);round(sd,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[4], cov = solve(emvg$hessian)[4,4])
ic_theta <- numeric()
ic_theta[1] <- epar[6] -qnorm(0.975)*sd
ic_theta[2] <- epar[6] +qnorm(0.975)*sd
round(ic_theta,3);round(epar[6],3);round(sd,3)

sd <- c(sd0,sdg)
ic_beta1g <- numeric()
ic_beta1g[1] <- epar[7] -qnorm(0.975)*sd[7]
ic_beta1g[2] <- epar[7] +qnorm(0.975)*sd[7]
round(ic_beta1g,3);round(epar[7],3);round(sd[7],3)



# VARIABLE TYPE OF DEBT
# VARIÁVEL TIPO DE DÍVIDA

data <- dados[,c(1,2,6)]
data0$covar <- data[,3]

emvg <- optim(par=c(-0.1, -0.1, -0.5, -0.5, 0.1), veroGTDLv2, x=data, hessian = TRUE, method = "BFGS");emvg
#emvg <- optim(par=c(0.1, 0.5, -0.5, -0.5, 0.1), veroGTDLv2, x=data, hessian = TRUE);emvg
egpar <- c(emvg$par[1], emvg$par[2], exp(emvg$par[3]), exp(emvg$par[4]), emvg$par[5])

emv0 <- optim(par=c(2,1), veroZerov2, x=data0, hessian = TRUE);emv0
e0par <- c(emv0$par[1], emv0$par[2])

epar <- c(e0par, egpar); epar

ekm <- survfit(Surv(tempo, status) ~ x4, data = data)
ekm_df <- data.frame(tempo = ekm$time, km = ekm$surv, 
                     covar = c(rep(0,ekm$strata[[1]]),rep(1,ekm$strata[[2]])))

ekm_df$gtdl0 <- sobrevGTDLv2(ekm_df[,c(1,3)],epar)

g2 <- ggplot(data = ekm_df, mapping = aes(x = tempo, y = gtdl0)) +
    geom_line(aes(linetype = as.factor(covar))) +
    geom_line(aes(x=tempo, y=km, color=factor(covar))) +
    labs(x="Tempo (meses)", y="Estimativa de S(t)", 
         color="Tipo de Dívida por Kaplan-Meier",
         linetype="Tipo de Dívida por GTDL \nZero Inflacionado v2") +
    scale_color_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                       values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_classic() + 
    theme(legend.position = c(0.75,0.7),
          legend.background = element_rect(color = "white"))

ggpubr::ggarrange(g1,g2)

sd0 <- sqrt(diag(solve(emv0$hessian)))
sdg <- sqrt(diag(solve(emvg$hessian)))
sd <- c(sd0,sdg); sd

ic_beta00 <- numeric()
ic_beta00[1] <- epar[1] -qnorm(0.975)*sd[1]
ic_beta00[2] <- epar[1] +qnorm(0.975)*sd[1]
round(ic_beta00,3);round(epar[1],3);round(sd[1],3)

ic_beta10 <- numeric()
ic_beta10[1] <- epar[2] -qnorm(0.975)*sd[2]
ic_beta10[2] <- epar[2] +qnorm(0.975)*sd[2]
round(ic_beta10,3);round(epar[2],3);round(sd[2],3)

ic_alpha0 <- numeric()
ic_alpha0[1] <- epar[3] -qnorm(0.975)*sd[3]
ic_alpha0[2] <- epar[3] +qnorm(0.975)*sd[3]
round(ic_alpha0,3);round(epar[3],3);round(sd[3],3)

ic_alpha1 <- numeric()
ic_alpha1[1] <- epar[4] -qnorm(0.975)*sd[4]
ic_alpha1[2] <- epar[4] +qnorm(0.975)*sd[4]
round(ic_alpha1,3);round(epar[4],3);round(sd[4],3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_lambda <- numeric()
ic_lambda[1] <- epar[5] -qnorm(0.975)*sd
ic_lambda[2] <- epar[5] +qnorm(0.975)*sd
round(ic_lambda,3);round(epar[5],3);round(sd,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[4], cov = solve(emvg$hessian)[4,4])
ic_theta <- numeric()
ic_theta[1] <- epar[6] -qnorm(0.975)*sd
ic_theta[2] <- epar[6] +qnorm(0.975)*sd
round(ic_theta,3);round(epar[6],3);round(sd,3)

sd <- c(sd0,sdg)
ic_beta1g <- numeric()
ic_beta1g[1] <- epar[7] -qnorm(0.975)*sd[7]
ic_beta1g[2] <- epar[7] +qnorm(0.975)*sd[7]
round(ic_beta1g,3);round(epar[7],3);round(sd[7],3)





# Model 5 GTDL Final -----------------------------------------------------------------


# model final com covariáveis significativas


# likelihood of zero inflation
# verossimilhança da inflação de zeros
veroZeroFinal <- function(x, gamm) {
    
    # covariates
    X0 <- as.matrix(x[x$tempo == 0, 2:ncol(x)])
    X <- as.matrix(x[x$tempo > 0, 2:ncol(x)])
    
    # model
    aux1 <- sum((X0%*%gamm)) - sum(log( 1 + exp(X0%*%gamm))) - sum(log(1 + exp(X%*%gamm)))
    return(-aux1)
}

# likelihood of GTDL Gamma
# verossimilhança da GTDL Gamma
veroGTDLFinal <- function(x, par) {
    
    # parameters
    alph <- par[1:2]
    lambd <- exp(par[3])
    thet <- exp(par[4])
    bet <- par[5:length(par)]
    
    # covariates
    X <- as.matrix(x[x$tempo > 0, 3])
    xalpha <- matrix(c(rep(1,nrow(X)),X), ncol = 2)
    cens <- x[x$tempo > 0,1]
    tempo <- x[x$tempo > 0,2]
    
    # model
    aux2 <- log(lambd)*sum(cens) +
        sum( cens*((xalpha%*%alph)*tempo + X%*%bet)) -
        sum( cens*log( 1 + exp((xalpha%*%alph)*tempo + X%*%bet) ) ) - 
        sum( (cens + (1/thet))*log( 1 + ((thet*lambd/(xalpha%*%alph))*log( (1 + exp((xalpha%*%alph)*tempo + X%*%bet))/(1 + exp(X%*%bet)) )) ) )
    
    return(-aux2)
}

# GTDL Zero Gamma Inflated Survival Model
# modelo de sobrevivência da GTDL Gama zero inflacionada
sobrevGTDL0Final <- function(x, par) {
    
    gamm0 <- par[1:3]
    alph <- par[4:5]
    lambd <- par[6]
    thet <- par[7]
    bet <- par[8:length(par)]
    tempo <- x$tempo
    X <- as.matrix(x[, 3])
    auX <- as.matrix(x[, 2:3])
    X0 <- matrix(data = c(rep(1,nrow(x)),auX), nrow = nrow(x), ncol = ncol(x))
    xalpha <- matrix(c(rep(1,nrow(X)),X), ncol = 2)
    st <- (1/(1 + exp(X0%*%gamm0) )) * (1 + ( (lambd*thet/(xalpha%*%alph)) * log( (1 + exp((xalpha%*%alph)*tempo + X%*%bet) )/(1 + exp(X%*%bet))) ) )^(-1/thet)
    
    return(st)
}



# Aplication Model 5 --------------------------------------------------------------


data <- dados[,c(1,2,6)]
data0 <- data.frame(tempo = dados$tempo,interc = rep(1,nrow(dados)),
                    consulta = dados$x1, tipoDiv = dados$x4)

# parameter estimation of zero inflation
# estimação dos parâmetros da inflação de zeros
emv0 <- optim(par=c(2, 1, 2), veroZeroFinal, x=data0, hessian = TRUE);emv0
e0par <- c(emv0$par[1], emv0$par[2], emv0$par[3])

# parameter estimation of GTDL Gamma
# estimação dos parâmetros do GTDL Gamma
emvg <- optim(par=c(0.1, 0.5, -0.5, -0.5, 0.1), veroGTDLFinal, x=data, hessian = TRUE);emvg
egpar <- c(emvg$par[1], emvg$par[2], exp(emvg$par[3]), exp(emvg$par[4]), emvg$par[5])

# combinando os parâmetros estimados
epar <- c(e0par, egpar); epar

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo, event = status) ~ x1 + x4, data = dados)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     consulta = c( rep(0,32), rep(0,34), rep(1,61), rep(1,61) ),
                     tipoDiv = c( rep(0,32), rep(1,34), rep(0,61), rep(1,61) ) )

# GTDL survival curve estimation
# estimação da curva de sobrevivência GTDL
ekm_df$gtdl <- sobrevGTDL0Final(x = ekm_df[,c(1,3,4)], par = epar)


# estimation of the confidence interval of the parameters
# estimação do intervalo de confiança dos parâmetros
sd0 <- sqrt(diag(solve(emv0$hessian)))
sdg <- sqrt(diag(solve(emvg$hessian)))
sd <- c(sd0,sdg); sd

ic_gama00 <- numeric()
ic_gama00[1] <- epar[1] -qnorm(0.975)*sd[1]
ic_gama00[2] <- epar[1] +qnorm(0.975)*sd[1]
round(ic_gama00,3);round(epar[1],3);round(sd[1],3)

ic_gamma1 <- numeric()
ic_gamma1[1] <- epar[2] -qnorm(0.975)*sd[2]
ic_gamma1[2] <- epar[2] +qnorm(0.975)*sd[2]
round(ic_gamma1,3);round(epar[2],3);round(sd[2],3)

ic_gamma2 <- numeric()
ic_gamma2[1] <- epar[3] -qnorm(0.975)*sd[3]
ic_gamma2[2] <- epar[3] +qnorm(0.975)*sd[3]
round(ic_gamma2,3);round(epar[3],3);round(sd[3],3)

ic_alpha0 <- numeric()
ic_alpha0[1] <- epar[4] -qnorm(0.975)*sd[4]
ic_alpha0[2] <- epar[4] +qnorm(0.975)*sd[4]
round(ic_alpha0,3);round(epar[4],3);round(sd[4],3)

ic_alpha1 <- numeric()
ic_alpha1[1] <- epar[5] -qnorm(0.975)*sd[5]
ic_alpha1[2] <- epar[5] +qnorm(0.975)*sd[5]
round(ic_alpha1,3);round(epar[5],3);round(sd[5],3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_lambda <- numeric()
ic_lambda[1] <- epar[6] -qnorm(0.975)*sd
ic_lambda[2] <- epar[6] +qnorm(0.975)*sd
round(ic_lambda,3);round(epar[6],3);round(sd,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[4], cov = solve(emvg$hessian)[4,4])
ic_theta <- numeric()
ic_theta[1] <- epar[7] -qnorm(0.975)*sd
ic_theta[2] <- epar[7] +qnorm(0.975)*sd
round(ic_theta,3);round(epar[7],3);round(sd,3)

sd <- c(sd0,sdg)
ic_beta1 <- numeric()
ic_beta1[1] <- epar[8] -qnorm(0.975)*sd[8]
ic_beta1[2] <- epar[8] +qnorm(0.975)*sd[8]
round(ic_beta1,3);round(epar[8],3);round(sd[8],3)





# Adequação do Modelo -----------------------------------------------------


# Resíduo de Cox-Snell para o Modelo GTDL Gamma Zero Inflacionado
residuo <- function(x, par) {
    
    # parâmetros estimados
    gamm0 <- par[1:3]
    alph <- par[4:5]
    lambd <- par[6]
    thet <- par[7]
    bet <- par[8:length(par)]
    # covariáveis
    tempo <- x$tempo
    X <- as.matrix(x[, 3])
    auX <- as.matrix(x[, 2:3])
    X0 <- matrix(data = c(rep(1,nrow(x)),auX), nrow = nrow(x), ncol = ncol(x))
    xalpha <- matrix(c(rep(1,nrow(X)),X), ncol = 2)
    # resíduo estimado
    ei = log(1 + exp(X0%*%gamm0)) + (1/thet)*log(1 + ( (lambd*thet/(xalpha%*%alph)) * log( (1 + exp((xalpha%*%alph)*tempo + X%*%bet) )/(1 + exp(X%*%bet))) ))
    
    return(ei)
}

data <- dados[,c(2,3,6)]; colnames(data) <- c("tempo","consulta","tipoDiv")
ei <- residuo(data, epar)

df_residuo <- data.frame(tempo = data$tempo,
                         residuo = ei,
                         cens = dados$status,
                         res_exp = exp(-ei))
ekm <- survfit(Surv(time = residuo, event = cens) ~ 1, data = df_residuo)
ekm_df <- data.frame(km = ekm$surv, ei = ekm$time,
                          sobre_exp_ei = exp(-ekm$time))

# gráfico de adequação 1
ggplot(data = ekm_df, aes(x = ei, y = km)) +
    geom_line(aes(linetype = "solid")) + 
    geom_line(aes(x = residuo, y = res_exp,  linetype = "dashed"), data = df_residuo) +
    labs(x = "Resíduos", y = "S(t) estimada", linetype = "Curvas de Sobrevivência") +
    scale_linetype_manual(labels = c("Kaplan-Meier","Modelo Exponencial Padrão"),
                          values = c("solid","dashed")) +
    ylim(0:1) + theme_classic() + 
    theme(legend.position = c(0.75,0.7),
          legend.background = element_rect(color = "white"))


# gráfico de adequação 2
ggplot(data = ekm_df, aes(x = km, y = sobre_exp_ei)) +
    geom_point() +
    geom_abline(aes(intercept = 0, slope = 1)) +
    xlim(0,1) + ylim(0,1) +
    theme_classic() + 
    theme(legend.position = c(0.75,0.7),
          legend.background = element_rect(color = "white"))




 # AIC e BIC --------------------------------------------




# critérios para seleção do modelo
lv <- -veroGTDLp0(data, ep0par)
maic <- function(lv,k) { -2*lv + 2*k }
mbic <- function(lv,k,n) { -2*lv + k*log(n) }
maic(lv,5)
mbic(lv,5,9645)



