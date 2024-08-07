
rm(list = ls())


# libraries ---------------------------------------------------------------

library(survival)
library(ggplot2)
library(msm)


# final version -----------------------------------------------------------


# verossimilhança da inflação de zeros
veroZero <- function(w, gamma) {
    
    W0 <- as.matrix(w[w$tempo == 0, 2:ncol(w)])
    W <- as.matrix(w[w$tempo > 0, 2:ncol(w)])
    aux1 <- sum((W0%*%gamma)) - sum(log( 1 + exp(W0%*%gamma))) - sum(log(1 + exp(W%*%gamma)))
    return(-aux1)
}

# verossimilhança do GTDL com fragilidade Gama
veroGTDL <- function(x, par) {
    
    alpha <- par[1]
    lambda <- exp(par[2])
    theta <- exp(par[3])
    beta <- par[4]
    
    cens <- x[x$tempo > 0,1]
    tempo <- x[x$tempo > 0,2]
    X <- as.matrix(x[x$tempo > 0,3])
    
    aux2 <- log(lambda)*sum(cens) +
        sum( cens*(alpha*tempo + X%*%beta)) -
        sum( cens*log( 1 + exp(alpha*tempo + X%*%beta) ) ) - 
        sum( (cens + (1/theta))*log( 1 + ((theta*lambda/alpha)*log( (1 + exp(alpha*tempo + X%*%beta))/(1 + exp(X%*%beta)) )) ) )
    return(-aux2)
}


# função de sobrevivência do modelo com regressão apenas na inflação de zero
sobrevGTDLv1 <- function(t, x, par) {
    
    gamma <- par[1:2]
    alpha <- par[3]
    lambda <- par[4]
    theta <- par[5]
    beta <- par[6]
    
    tempo <- t
    X <- as.matrix(x)
    X0 <- matrix(data = c(rep(1,nrow(X)),X), nrow = nrow(X), ncol = 2)
    
    st <- (1/(1 + exp(X0%*%gamma) )) * (1 + ( (lambda*theta/alpha) * log( (1 + exp(alpha*tempo + X%*%beta) )/(1 + exp(X%*%beta))) ) )^(-1/theta)
    return(st)
}


# função que calcula a proporção de zeros inflacionados dados os níveis da covar
p0f <- function(w, gamma) {
    
    w <- as.matrix(w)
    W <- matrix(data = c(rep(1,nrow(w)),w), nrow = nrow(w), ncol = 2)
    
    p <- exp(W%*%gamma)/(1+exp(W%*%gamma))
    return(p)
}


# função que calcula a fração de cura dados os níveis da covar
p1f <- function(x, par) {
    gamma <- par[1:2]
    alpha <- par[3]
    lambda <- par[4]
    theta <- par[5]
    beta1 <- par[6]
    
    X <- as.matrix(x)
    W <- matrix(data = c(rep(1,nrow(X)),X), nrow = nrow(X), ncol = 2)
    
    p = (1/(1+exp(W%*%gamma))) * (1 + ((lambda*theta)/alpha)*log(1/(1+exp(X%*%beta1))) )^(-1/theta)
    return(p)
}


# função criada para usar no comando uniroot e encontrar a raiz
froot <- function(t,x,par,u) sobrevGTDLv1(t,x,par) - 1 + u


# parâmetros fixados para simulação
set.seed(2024)
n = 100
X = rbinom(n, size = 1, prob = 0.5)
p0 = p0f(X, gamma = c(-1.029,-0.302))
p1 = p1f(X, par = c(-1.029,-0.302,-0.140,0.492,0.399,-0.932))


t <- numeric()
delta <- numeric()
epar = sd <- matrix(nrow = 1000, ncol = 6)
ic_gamma0 = ic_gamma1 = ic_alpha = ic_beta1 = ic_lambda = ic_theta = numeric()
contador = numeric(6)


# iteração pra simulação e estimação dos parâmetros
for (j in 1:1000) {
    
    # processo de geração dos tempos
    for (i in 1:n) {
        
        u1 <- runif(1,0,1)
        u2 <- runif(1, min = p0[i], max = (1-p1[i]))
        
        tf <- ifelse(test = u1 < p0[i],
                     yes = 0,
                     no = ifelse(test = (1 - p1[i]) < u1,
                                 yes = Inf,
                                 no = uniroot(f = froot, x=X[i], u=u2,
                                              par = c(-1.029,-0.302,-0.140,
                                                      0.492,0.399,-0.932),
                                              interval = c(0,1080))$root
                     )
        )
        
        tf <- as.numeric(tf)
        tc <- runif(1, 0, 65)
        t[i] <- min(tf,tc)
        
        delta[i] <- ifelse(t[i] == tf, 1, 0)
    }
    
    # processo de estimação dos parâmetros
    dadoSimulado <- data.frame(cens = delta, tempo = t, x1 = X)
    dado0 <- data.frame(tempo = t, w0 = rep(1,n), w1 = X)
    
    emvg <- optim(par=c(-0.1,-0.5,-0.5,1), veroGTDL, x=dadoSimulado, hessian = TRUE, method = "BFGS")
    egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])
    emv0 <- optim(par=c(1,1), veroZero, w=dado0, hessian = TRUE)
    e0par <- c(emv0$par[1], emv0$par[2])
    epar[j,] <- c(e0par, egpar)
    sd0 <- sqrt(diag(solve(emv0$hessian)))
    sdg <- sqrt(diag(solve(emvg$hessian)))
    #sd[j,] = c(sd0,sdg)
    
    
    # intervalos de confiança
    
    ic_gamma0[1] <- round(e0par[1] -qnorm(0.975)*sd0[1],3)
    ic_gamma0[2] <- round(e0par[1] +qnorm(0.975)*sd0[1],3)
    g0p <- ifelse(test = ic_gamma0[1] <= -1.029 & -1.029 <= ic_gamma0[2], yes = 1, no = 0)
    contador[1] <- contador[1] + g0p
    
    ic_gamma1[1] <- round(e0par[2] -qnorm(0.975)*sd0[2],3)
    ic_gamma1[2] <- round(e0par[2] +qnorm(0.975)*sd0[2],3)
    g1p <- ifelse(test = ic_gamma1[1] <= -0.302 & -0.302 <= ic_gamma1[2], yes = 1, no = 0)
    contador[2] <- contador[2] + g1p
    
    ic_alpha[1] <- round(egpar[1] -qnorm(0.975)*sdg[1],3)
    ic_alpha[2] <- round(egpar[1] +qnorm(0.975)*sdg[1],3)
    ap <- ifelse(test = ic_alpha[1] <= -0.140 & -0.140 <= ic_alpha[2], yes = 1, no = 0)
    contador[3] <- contador[3] + ap
    
    ic_beta1[1] <- round(egpar[4] -qnorm(0.975)*sdg[4],3)
    ic_beta1[2] <- round(egpar[4] +qnorm(0.975)*sdg[4],3)
    b1p <- ifelse(test = ic_beta1[1] <= -0.932 & -0.932 <= ic_beta1[2], yes = 1, no = 0)
    contador[6] <- contador[6] + b1p
    
    # desvio padrão pelo método delta
    sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = solve(emvg$hessian)[2,2])
    ic_lambda[1] <- round(egpar[2] -qnorm(0.975)*sd,3)
    ic_lambda[2] <- round(egpar[2] +qnorm(0.975)*sd,3)
    lp <- ifelse(test = ic_lambda[1] <= 0.492 & 0.492 <= ic_lambda[2], yes = 1, no = 0)
    contador[4] <- contador[4] + lp
    
    # desvio padrão pelo método delta
    sd <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
    ic_theta[1] <- round(egpar[3] -qnorm(0.975)*sd,3)
    ic_theta[2] <- round(egpar[3] +qnorm(0.975)*sd,3)
    tp <- ifelse(test = ic_theta[1] <= 0.399 & 0.399 <= ic_theta[2], yes = 1, no = 0)
    contador[5] <- contador[5] + tp
    
    
}


# análise gráfica
eppar <- c(mean(epar[,1]),mean(epar[,2]),mean(epar[,3]),mean(epar[,4]),
           mean(epar[,5]),mean(epar[,6]))

ekm <- survfit(Surv(time = tempo, event = cens) ~ X, data = dadoSimulado)
ekm_df <- data.frame(tempo = ekm$time, km = ekm$surv,
                     covar = c(rep(0,ekm$strata[1]),rep(1,ekm$strata[2])) )

dadoSimulado$gtdl_surv <- sobrevGTDLv1(t = dadoSimulado$tempo, x = dadoSimulado$x1, par = eppar)


ggplot(data = ekm_df, mapping = aes(x = tempo, y = km)) +
    geom_step(aes(color = as.factor(covar))) +
    geom_line(data = dadoSimulado, aes(x = tempo, y = gtdl_surv, linetype = as.factor(x1))) +
    labs(x="Tempo (meses)", y="Estimativa de S(t)", 
         color="Tipo de Dívida por Kaplan-Meier",
         linetype="Tipo de Dívida por GTDL \nZero Inflacionado") +
    scale_color_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                       values = c("red","blue")) +
    scale_linetype_manual(labels = c("0 - Bancos", "1 - Outros Segmentos"),
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_light() + 
    theme(legend.position = c(0.7,0.7),
          legend.background = element_rect(color = "white"))


# estimativas pontuais

Vies <- matrix(nrow = 6, ncol = 3)
Variancia <- matrix(nrow = 6, ncol = 3)
EQM <- matrix(nrow = 6, ncol = 3)
PC <- matrix(nrow = 6, ncol = 3)
real <- c(-1.029,-0.302,-0.140,0.492,0.399,-0.932)

vies <- function(xpars, parametro) mean(xpars) - parametro
eqm <- function(xpars, parametro) sum((xpars-parametro)^2)/1000


# construção da tabela de verificação dos estimadores 
for (k in 1:6) {
    
    Vies[k,1] <- round(vies(xpars = epar[,k], parametro = real[k]),3)
    EQM[k,1] <- round(eqm(xpars = epar[,k], parametro = real[k]),3)
    Variancia[k,1] <- round(var(epar[,k]),3)
    PC[k,1] <- contador[k]/1000
}


