rm(list = ls())

# Libraries ---------------------------------------------------------------

library(ggplot2)
library(survival)
library(msm) # pacote pra usar o método delta

# Load and Treatment Data -------------------------------------------------

dados <- read.table(file = "dados_censura_versao4.txt", header = TRUE)
data0 <- data.frame(status = dados$status,
                    tempo = dados$tempo,
                    interc = rep(1,nrow(dados)))


# Descriptive Analysis ------------------------------------------------------

# response var histogram Time until debt payment
# histograma da var resposta Tempo até o Pagamento da dívida
hist(data$tempo, xlab="Tempo Pagamento", 
     main = "Histograma do Tempo de Pagamento")

# estimated Kaplan-Meier curve
# curva de Kaplan-Meier estimada
ekm <- survfit(Surv(time = tempo, event = status) ~ 1, data = data)
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
    labs(x = "Tempos", y = "S(t)") +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo, event = status) ~ consulta, data = data)
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
    labs(x = "Tempos", y = "S(t)", color = "Inform. Consulta") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(2,4)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))



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
    tempo <- x$tempo
    X <- as.matrix(x[, 2:ncol(x)])
    st <- (1 + ( (lambd*thet/alph)*log( (1+ exp(X%*%bet + alph*tempo))/(1 + exp(X%*%bet))) ) )^(-1/thet)
    return(st)
}



# Aplication Model 1 ------------------------------------------------------


# soma 10^-3 nos tempos 0
# espero do gráfico que as curvas se ultrapassem

data$tempo[data$tempo == 0] <- 10^-3

# parameter estimation of GTDL Gamma
emvg <- optim(par=c(0.02,0.5,0.5,-1), veroGTDLg, x=data, hessian = TRUE)
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo, event = status) ~ consulta, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     consulta = c(rep(0,44),rep(1,61)) )

# GTDL survival curve estimation
# estimação da curva de sobrevivência GTDL
ekm_df$gtdlg <- sobrevGTDLg(x = ekm_df[,c(1,3)], par = egpar)

# Graph of survival curves by Kaplan Meyer (KM) and GTDL
# Gráfico das curvas de sobrevivência por Kaplan Meyer (KM) e por GTDL
ggplot(data = ekm_df, mapping = aes(tempo, gtdlg)) +
    geom_line(aes(linetype = as.factor(consulta))) +
    geom_line(aes(x = tempo, y = km, color = as.factor(consulta))) +
    labs(x = "Tempos", y = "S(t)", color = "Consulta por KM", linetype = "Consulta por GTDL Gama") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(color = "white"))

# estimation of the confidence interval of the parameters
# estimação do intervalo de confiança dos parâmetros
sd <- sqrt(diag(solve(emvg$hessian)))

ic_alpha <- numeric()
ic_alpha[1] <- egpar[1] -qnorm(0.975)*sd[1]
ic_alpha[2] <- egpar[1] +qnorm(0.975)*sd[1]
round(ic_alpha,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = solve(emvg$hessian)[2,2])
ic_lambda <- numeric()
ic_lambda[1] <- egpar[2] -qnorm(0.975)*sd
ic_lambda[2] <- egpar[2] +qnorm(0.975)*sd
round(ic_lambda,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_theta <- numeric()
ic_theta[1] <- egpar[3] -qnorm(0.975)*sd
ic_theta[2] <- egpar[3] +qnorm(0.975)*sd
round(ic_theta,3)

sd <- sqrt(diag(solve(emvg$hessian)))
ic_beta1g <- numeric()
ic_beta1g[1] <- egpar[4] -qnorm(0.975)*sd[4]
ic_beta1g[2] <- egpar[4] +qnorm(0.975)*sd[4]
round(ic_beta1g,3)



# Model 2 GTDL p0 ---------------------------------------------------------

# REVER A CONSTRUÇÃO DO GRÁFICO, TALVEZ NÃO SEJA O MODELO

# likelihoof of GTDL Gamma p0
veroGTDLp0 <- function(x, par) {
    
    # parameters
    alph <- par[1]
    lambd <- exp(par[2])
    thet <- exp(par[3])
    bet <- par[4:length(par)]
    p0 <- par[5]
    n0 <- nrow(x[x$tempo == 0,])
    n1 <- nrow(x[x$tempo > 0,])
    
    # covariates
    X <- as.matrix(x[, 3:ncol(x)])
    cens <- x[,1]
    tempo <- x[,2]
    
    # model
    aux2 <- n0*log(p0) + n1*log(1-p0) + log(lambd)*sum(cens) +
        sum( cens*(X%*%bet + alph*tempo)) -
        sum( cens*log( 1 + exp(X%*%bet + alph*tempo)) ) - 
        sum( (cens + (1/thet))*log( 1 + ((thet*lambd/alph)*log( (1 + exp(X%*%bet + alph*tempo))/(1 + exp(X%*%bet)) )) ))
    
    return(-aux2)
}

# GTDL Gamma p0 Survival Model
sobrevGTDLp0 <- function(x, par) {
    
    alph <- par[1]
    lambd <- par[2]
    thet <- par[3]
    bet <- par[4:length(par)]
    p0 <- par[5]
    tempo <- x$tempo
    X <- as.matrix(x[, 2:ncol(x)])
    st <- (1-p0)*(1 + ( (lambd*thet/alph)*log( (1+ exp(X%*%bet + alph*tempo))/(1 + exp(X%*%bet))) ) )^(-1/thet)
    st
}



# Aplication Model 2 ------------------------------------------------------


# parameter estimation of GTDL Gamma p0
emvp0 <- optim(par=c(0.02,0.5,0.5,0.1,0.25), veroGTDLp0, x=data, hessian = TRUE)
ep0par <- c(emvp0$par[1], exp(emvp0$par[2]), exp(emvp0$par[3]), emvp0$par[4], emvp0$par[5])

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo, event = status) ~ consulta, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     consulta = c(rep(0,44),rep(1,61)) )

# GTDL survival curve estimation
# estimação da curva de sobrevivência GTDL
ekm_df$gtdlp0 <- sobrevGTDLp0(x = ekm_df[,c(1,3)], par = ep0par)

# Graph of survival curves by Kaplan Meyer (KM) and GTDL
# Gráfico das curvas de sobrevivência por Kaplan Meyer (KM) e por GTDL
ggplot(data = ekm_df, mapping = aes(tempo, gtdlp0[,1])) +
    geom_line(aes(linetype = as.factor(consulta))) +
    geom_line(aes(x = tempo, y = km, color = as.factor(consulta))) +
    labs(x = "Tempos", y = "S(t)", color = "Consulta por KM", linetype = "Consulta por GTDL Gama p0") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(color = "white"))

# critérios para seleção do modelo
lv <- -veroGTDLp0(data, ep0par)
maic <- function(lv,k) { -2*lv + 2*k }
mbic <- function(lv,k,n) { -2*lv + k*log(n) }
maic(lv,5)
mbic(lv,5,9645)


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

# VARIÁVEL CONSULTA

data <- dados[,1:3]
data0$covar <- data[,3]

# parameter estimation of GTDL Gamma
# estimação dos parâmetros do GTDL Gamma
emvg <- optim(par=c(0.1,-0.5,-0.5,0.1), veroGTDLv1, x=data, hessian = TRUE)
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])

# parameter estimation of zero inflation
# estimação dos parâmetros da inflação de zeros
emv0 <- optim(par=c(-0.5,0.5), veroZerov1, x=data0, hessian = TRUE)
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
ggplot(data = ekm_df, mapping = aes(tempo, gtdl)) +
    geom_line(aes(linetype = as.factor(consulta))) +
    geom_line(aes(x = tempo, y = km, color = as.factor(consulta))) +
    labs(x = "Tempos", y = "S(t)", color = "Consulta por KM", linetype = "Consulta por GTDL Zero Inflacionado v1") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(1,2)) +
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
round(ic_beta00,3)

ic_beta10 <- numeric()
ic_beta10[1] <- epar[2] -qnorm(0.975)*sd[2]
ic_beta10[2] <- epar[2] +qnorm(0.975)*sd[2]
round(ic_beta10,3)

ic_alpha <- numeric()
ic_alpha[1] <- epar[3] -qnorm(0.975)*sd[3]
ic_alpha[2] <- epar[3] +qnorm(0.975)*sd[3]
round(ic_alpha,3)


# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = solve(emvg$hessian)[2,2])
ic_lambda <- numeric()
ic_lambda[1] <- epar[4] -qnorm(0.975)*sd
ic_lambda[2] <- epar[4] +qnorm(0.975)*sd
round(ic_lambda,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_theta <- numeric()
ic_theta[1] <- epar[5] -qnorm(0.975)*sd
ic_theta[2] <- epar[5] +qnorm(0.975)*sd
round(ic_theta,3)

sd <- c(sd0,sdg)
ic_beta1g <- numeric()
ic_beta1g[1] <- epar[6] -qnorm(0.975)*sd[6]
ic_beta1g[2] <- epar[6] +qnorm(0.975)*sd[6]
round(ic_beta1g,3)




# VARIÁVEL TIPO DE DÍVIDA

data <- dados[,c(1,2,6)]
data0$covar <- data[,3]


emvg <- optim(par=c(0.1,-0.5,-0.5,0.1), veroGTDLv1, x=data, hessian = TRUE)
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])

emv0 <- optim(par=c(-0.5,0.5), veroZerov1, x=data0, hessian = TRUE)
e0par <- c(emv0$par[1], emv0$par[2])

epar <- c(e0par, egpar); epar


ekm <- survfit(Surv(tempo, status) ~ x4, data = data)
ekm_df <- data.frame(tempo = ekm$time, km = ekm$surv, 
                     covar = c(rep(0,ekm$strata[[1]]),rep(1,ekm$strata[[2]])))

ekm_df$gtdl0 <- sobrevGTDLv1(ekm_df[,c(1,3)],epar)

ggplot(data = ekm_df, mapping = aes(x = tempo, y = gtdl0)) +
    geom_line(aes(linetype = as.factor(covar))) +
    geom_line(aes(x=tempo, y=km, color=factor(covar))) +
    labs(x="Tempo em meses", y="S(t)", color="Tipo de Dívida por KM",
         linetype="Tipo de Dívida por GTDL Zero Inflacionado") +
    scale_color_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"), values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"), values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_classic() + 
    theme(legend.position = c(0.75,0.7),
          legend.background = element_rect(color = "white"))


# estimation of the confidence interval of the parameters
# estimação do intervalo de confiança dos parâmetros
sd0 <- sqrt(diag(solve(emv0$hessian)))
sdg <- sqrt(diag(solve(emvg$hessian)))
sd <- c(sd0,sdg); sd

ic_beta00 <- numeric()
ic_beta00[1] <- epar[1] -qnorm(0.975)*sd[1]
ic_beta00[2] <- epar[1] +qnorm(0.975)*sd[1]
round(ic_beta00,3)

ic_beta10 <- numeric()
ic_beta10[1] <- epar[2] -qnorm(0.975)*sd[2]
ic_beta10[2] <- epar[2] +qnorm(0.975)*sd[2]
round(ic_beta10,3)

ic_alpha <- numeric()
ic_alpha[1] <- epar[3] -qnorm(0.975)*sd[3]
ic_alpha[2] <- epar[3] +qnorm(0.975)*sd[3]
round(ic_alpha,3)


# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = solve(emvg$hessian)[2,2])
ic_lambda <- numeric()
ic_lambda[1] <- epar[4] -qnorm(0.975)*sd
ic_lambda[2] <- epar[4] +qnorm(0.975)*sd
round(ic_lambda,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_theta <- numeric()
ic_theta[1] <- epar[5] -qnorm(0.975)*sd
ic_theta[2] <- epar[5] +qnorm(0.975)*sd
round(ic_theta,3)

sd <- c(sd0,sdg)
ic_beta1g <- numeric()
ic_beta1g[1] <- epar[6] -qnorm(0.975)*sd[6]
ic_beta1g[2] <- epar[6] +qnorm(0.975)*sd[6]
round(ic_beta1g,3)




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

# parameter estimation of GTDL Gamma
# estimação dos parâmetros do GTDL Gamma
emvg <- optim(par=c(0.1,0.5,-0.5,-0.5,0.1), veroGTDLv2, x=data, hessian = TRUE)
egpar <- c(emvg$par[1], emvg$par[2], exp(emvg$par[3]), exp(emvg$par[4]), emvg$par[5])

# parameter estimation of zero inflation
# estimação dos parâmetros da inflação de zeros
emv0 <- optim(par=c(2,1), veroZerov2, x=data0, hessian = TRUE)
e0par <- c(emv0$par[1], emv0$par[2])

# combinando os parâmetros estimados
epar <- c(e0par, egpar); epar

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo, event = status) ~ consulta, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     consulta = c(rep(0,44),rep(1,61)) )

# GTDL survival curve estimation
# estimação da curva de sobrevivência GTDL
ekm_df$gtdl <- sobrevGTDLv2(x = ekm_df[,c(1,3)], par = epar)

# Graph of survival curves by Kaplan Meyer (KM) and GTDL
# Gráfico das curvas de sobrevivência por Kaplan Meyer (KM) e por GTDL
ggplot(data = ekm_df, mapping = aes(tempo, gtdl)) +
    geom_line(aes(linetype = as.factor(consulta))) +
    geom_line(aes(x = tempo, y = km, color = as.factor(consulta))) +
    labs(x = "Tempos", y = "S(t)", color = "Consulta por KM", linetype = "Consulta por GTDL Zero Inflacionado v2") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(2,4)) +
    scale_linetype_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(1,2)) +
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
round(ic_beta00,3)

ic_beta10 <- numeric()
ic_beta10[1] <- epar[2] -qnorm(0.975)*sd[2]
ic_beta10[2] <- epar[2] +qnorm(0.975)*sd[2]
round(ic_beta10,3)

ic_alpha0 <- numeric()
ic_alpha0[1] <- epar[3] -qnorm(0.975)*sd[3]
ic_alpha0[2] <- epar[3] +qnorm(0.975)*sd[3]
round(ic_alpha0,3)

ic_alpha1 <- numeric()
ic_alpha1[1] <- epar[4] -qnorm(0.975)*sd[4]
ic_alpha1[2] <- epar[4] +qnorm(0.975)*sd[4]
round(ic_alpha1,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
ic_lambda <- numeric()
ic_lambda[1] <- epar[5] -qnorm(0.975)*sd
ic_lambda[2] <- epar[5] +qnorm(0.975)*sd
round(ic_lambda,3)

# desvio padrão pelo método delta
sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[4], cov = solve(emvg$hessian)[4,4])
ic_theta <- numeric()
ic_theta[1] <- epar[6] -qnorm(0.975)*sd
ic_theta[2] <- epar[6] +qnorm(0.975)*sd
round(ic_theta,3)

sd <- c(sd0,sdg)
ic_beta1g <- numeric()
ic_beta1g[1] <- epar[7] -qnorm(0.975)*sd[7]
ic_beta1g[2] <- epar[7] +qnorm(0.975)*sd[7]
round(ic_beta1g,3)



# COVARIÁVEL DE TIPO DE DÍVIDA --------------------------------------------




resultado <- survfit(Surv(time = tempo, event = status)~factor(x4), data = dados_censura_versao4)

resultado=survfit(Surv(time, delta)~factor(x4), conf.type="plain")
plot(resultado, xlim=c(0, 60), ylim=c(0, 1), xlab = " Tempo (meses)", ylab = "Estimativa de S(t)", lty = c(2, 2), col = c("black", "blue"),  bty="n")
legend(7, 1, title="Tipo de dívida", c("0 - Banco", "1 - outros segmentos"), lty = c(2, 2), col = c("black", "blue"), bty = "n")
TesteC4=survdiff(Surv(time, delta)~factor(x4), rho=0)
TesteC4
