rm(list = ls())

# Libraries ---------------------------------------------------------------

library(ggplot2)
library(survival)


# Load and Treatment Data -------------------------------------------------

data <- read.csv(file = "recuperaInadimplencia.csv")
data <- data[,2:4]
data0 <- data.frame(status = data$status,
                    tempo_pag = data$tempo_pag,
                    interc = rep(1,nrow(data)),
                    consulta = data$consulta)


# Descriptive Analysis ------------------------------------------------------

# response var histogram Time until debt payment
# histograma da var resposta Tempo até o Pagamento da dívida
hist(data$tempo_pag, xlab="Tempo Pagamento", 
     main = "Histograma do Tempo de Pagamento")

# estimated Kaplan-Meier curve
# curva de Kaplan-Meier estimada
ekm <- survfit(Surv(time = tempo_pag, event = status) ~ 1, data = data)
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
    labs(x = "Tempos", y = "S(t)", color = "Inform. Consulta") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"), values = c(2,4)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,60,5)) +
    theme_classic() + 
    theme(legend.position = c(0.8,0.8),
          legend.background = element_rect(color = "white"))



# Model 1 GTDL Gamma ------------------------------------------------------

# likelihoof of GRTDL Gamma
veroGTDLg <- function(x, par) {
    
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


# parameter estimation of GTDL Gamma
emvg <- optim(par=c(0.1,0.5,0.5,0.1), veroGTDLg, x=data, hessian = TRUE)
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo_pag, event = status) ~ consulta, data = data)
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



# Modeling ----------------------------------------------------------------


# likelihood of zero inflation
# verossimilhança da inflação de zeros
veroZero <- function(x, beta) {
    
    # covariates
    X0 <- as.matrix(x[x$tempo == 0, 3:ncol(x)])
    X <- as.matrix(x[x$tempo > 0, 3:ncol(x)])
     
    # model
    aux1 <- sum((X0%*%beta)) - sum(log( 1 + exp(X0%*%beta))) - sum(log(1 + exp(X%*%beta)))
    return(-aux1)
}

# likelihood of GTDL Gamma
# verossimilhança da GTDL Gamma
veroGTDL <- function(x, par) {
    
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
sobrevGTDL <- function(x, par) {
    
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



# Aplication --------------------------------------------------------------

# parameter estimation of GTDL Gamma
# estimação dos parâmetros do GTDL Gamma
emvg <- optim(par=c(0.1,0.5,-0.5,-0.5,0.1), veroGTDL, x=data, hessian = TRUE)
egpar <- c(emvg$par[1], emvg$par[2], exp(emvg$par[3]), exp(emvg$par[4]), emvg$par[5])

# parameter estimation of zero inflation
# estimação dos parâmetros da inflação de zeros
emv0 <- optim(par=c(2,1), veroZero, x=data0, hessian = TRUE)
e0par <- c(emv0$par[1], emv0$par[2])

# combinando os parâmetros estimados
epar <- c(e0par, egpar)

# estimated Kaplan-Meier curve with query covariate
# curva de Kaplan-Meier estimada com covariável de consulta
ekm <- survfit(Surv(time = tempo_pag, event = status) ~ consulta, data = data)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     consulta = c(rep(0,44),rep(1,61)) )

# GTDL survival curve estimation
# estimação da curva de sobrevivência GTDL
ekm_df$gtdl <- sobrevGTDL(x = ekm_df[,c(1,3)], par = epar)

# Graph of survival curves by Kaplan Meyer (KM) and GTDL
# Gráfico das curvas de sobrevivência por Kaplan Meyer (KM) e por GTDL
ggplot(data = ekm_df, mapping = aes(tempo, gtdl)) +
    geom_line(aes(linetype = as.factor(consulta))) +
    geom_line(aes(x = tempo, y = km, color = as.factor(consulta))) +
    labs(x = "Tempos", y = "S(t)", color = "Consulta por KM", linetype = "Consulta por GTDL") +
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
sdg;sd0

ic_alpha <- numeric()
ic_alpha[1] <- emv$par[1] -qnorm(0.975)*sd[1]
ic_alpha[2] <- emv$par[1] +qnorm(0.975)*sd[1]
ic_alpha

ic_lambda <- numeric()
ic_lambda[1] <- emv$par[2] -qnorm(0.975)*sd[2]
ic_lambda[2] <- emv$par[2] +qnorm(0.975)*sd[2]
ic_lambda

ic_theta <- numeric()
ic_theta[1] <- emv$par[3] -qnorm(0.975)*sd[3]
ic_theta[2] <- emv$par[3] +qnorm(0.975)*sd[3]
ic_theta

ic_beta <- numeric()
ic_beta[1] <- emv$par[4] -qnorm(0.975)*sd[4]
ic_beta[2] <- emv$par[4] +qnorm(0.975)*sd[4]
ic_beta
