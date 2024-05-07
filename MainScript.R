

# Libraries ---------------------------------------------------------------

library(ggplot2)
library(survival)


# Load and Treatment Data -------------------------------------------------

data <- read.csv(file = "recuperaInadimplencia.csv")
data <- data.frame(status = data$status,
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
    aux1 <- sum( (X0%*%bet) ) - sum( log( 1 + exp(X0%*%bet) ) )
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
emv <- optim(par=c(0.1,0.3,0.5,2), veroGTDL, x=data, hessian = TRUE)
epar <- c(emv$par[1], exp(emv$par[2]), exp(emv$par[3]), emv$par[4])

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
emv <- optim(par=c(0.1,0.3,0.5,2), veroGTDL, x=data, hessian = TRUE)
emv$par
sd <- sqrt(diag(solve(emv$hessian)))
sd

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
