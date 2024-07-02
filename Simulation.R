
veroZero <- function(w, gamma) {
    
    W0 <- as.matrix(w[w$tempo == 0, 2:ncol(w)])
    W <- as.matrix(w[w$tempo > 0, 2:ncol(w)])
    aux1 <- sum((W0%*%gamma)) - sum(log( 1 + exp(W0%*%gamma))) - sum(log(1 + exp(W%*%gamma)))
    return(-aux1)
}

veroGTDL <- function(x, par) {
    
    alpha <- par[1]
    lambda <- exp(par[2])
    theta <- exp(par[3])
    beta <- par[4:length(par)]
    
    X <- as.matrix(x[x$tempo > 0, 3:ncol(x)])
    cens <- x[x$tempo > 0,1]
    tempo <- x[x$tempo > 0,2]
    
    aux2 <- log(lambda)*sum(cens) +
        sum( cens*(alpha*tempo + X%*%beta)) -
        sum( cens*log( 1 + exp(alpha*tempo + X%*%beta) ) ) - 
        sum( (cens + (1/theta))*log( 1 + ((theta*lambda/alpha)*log( (1 + exp(alpha*tempo + X%*%beta))/(1 + exp(X%*%beta)) )) ) )
    return(-aux2)
}


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


p0f <- function(w, gamma) {
    
    w <- as.matrix(w)
    W <- matrix(data = c(rep(1,nrow(w)),w), nrow = nrow(w), ncol = 2)
    
    p <- exp(W%*%gamma)/(1+exp(W%*%gamma))
    return(p)
}

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

froot <- function(t,x,par,u) sobrevGTDLv1(t,x,par) - 1 + u


set.seed(2024)
n = 1000
X = rbinom(n, size = 1, prob = 0.5)
epar = c(-1.0293802, -0.3018615, -0.1403741,  0.4924747,  0.3988917, -0.9315740)
e0par = c(-1.0293802, -0.3018615)
alpha = rep(c(-0.1403741),n)
p0 = p0f(X, gamma = e0par)
p1 = p1f(X, par = epar)


t <- numeric()
delta <- numeric()

for (i in 1:n) {
    
    u1 <- runif(1,0,1)
    u2 <- runif(1, min = p0[i], max = (1-p1[i]))
    u3 <- runif(1, min = p0[i], max = 1)
    #cat("u1= ", u1, "u2= ", u2, "u3=", u3, )
    tf <- ifelse(test = alpha[i] < 0,
                 yes = ifelse(test = u1 < p0[i],
                              yes = 0,
                              no = ifelse(test = (1 - p1[i]) < u1,
                                          yes = Inf,
                                          no = uniroot(f = froot, x=X[i], u=u2,
                                                       par = epar, maxiter = 10000,
                                                       interval = c(0,1080), 
                                                       extendInt = "yes")$root
                              )
                 ),
                 no = ifelse(test = u1 < p0[i] , 
                             yes = 0, 
                             no = uniroot(f = froot, x=X[i], u = u3,
                                          par = c(-1.0293802, -0.3018615, 
                                                  0.1403741,  0.4924747,  
                                                  0.3988917, -0.9315740), 
                                          maxiter = 10000, interval = c(0,1080),
                                          extendInt = "yes")$root)
    )
    tf <- as.numeric(tf)
    tc <- runif(1, 0, 65)
    t[i] <- min(tf,tc)
    
    delta[i] <- ifelse(t[i] == tf, 1, 0)
    
}


library(survival)
library(ggplot2)

# processo de estimação dos parâmetros
dadoSimul <- data.frame(cens = delta, tempo = t, X)
data0 <- data.frame(tempo = dadoSimul$t, x0 = rep(1,nrow(dadoSimul)),
                    X = dadoSimul$X)

emvg <- optim(par=c(-0.2,1,0.5,-1.5), veroGTDL, x=dadoSimul, hessian = TRUE);emvg
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])
emv0 <- optim(par=c(1,1), veroZero, w=data0, hessian = TRUE);emv0
e0par <- c(emv0$par[1], emv0$par[2])
epar <- c(e0par, egpar); epar

ekm <- survfit(Surv(time = tempo, event = cens) ~ X, data = dadoSimul)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     covar = c(rep(0,ekm$strata[1]),rep(1,ekm$strata[2])) )

gtdl_df <- data.frame(tempo = t,covar = X)
#gtdl_df <- data.frame(tempo = rep(seq(0,60,1),2),covar = c(rep(0,61),rep(1,61)))
gtdl_df$surv <- sobrevGTDLv1(t = gtdl_df$tempo ,x = gtdl_df$covar, par = epar)


ggplot(data = ekm_df, mapping = aes(x = tempo, y = km)) +
    geom_step(aes(color = as.factor(covar))) +
    geom_line(data = gtdl_df, aes(x = tempo, y = surv, linetype = as.factor(covar))) +
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
