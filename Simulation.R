

# Model p0 ---------------------------------------------------------------


# GTDL Gamma p0 Survival Model
sobrevGTDLp0 <- function(t, x, par) {
    
    p0 <- par[1]
    alph <- par[2]
    lambd <- par[3]
    thet <- par[4]
    bet <- par[5:length(par)]
    tempo <- t
    X <- as.matrix(x)
    st <- (1-p0) * (1 + ( (lambd*thet/alph) * log( (1 + exp(alph*tempo + X%*%bet) )/(1 + exp(X%*%bet))) ) )^(-1/thet)
    return(st)
}


froot <- function(t,x,par,u) sobrevGTDLp0(t,x,par) - 1 + u



alpha = -1
p0 = 0.13
p1 = 0.18

ifelse(test = alpha < 0,
       yes = ifelse(test = runif(1,0,1) < p0,
                    yes = 0,
                    no = ifelse(test = 1 - p1 < runif(1,0,1),
                                yes = Inf,
                                no = uniroot(f = froot, x=1, par=c(0.2,1,0.5,0.5,1), 
                                             u=runif(1, min = p0, max = 1-p1), 
                                             interval = c(0,60))$root
                                )
                    ),
       no = ifelse(test = runif(1,0,1) < p0 , 
                   yes = 0, 
                   no = uniroot(f = froot, x=1, par=c(0.2,1,0.5,0.5,1), 
                                u=runif(n = 1, min = p0, max = 1), 
                                interval = c(0,1080))$root)
       )



# Model 0 v2 -------------------------------------------------------------




# GTDL Zero Gamma Inflated Survival Model
# modelo de sobrevivência da GTDL Gama zero inflacionada
sobrevGTDLv2 <- function(t, x, par) {
    
    gamma <- par[1:2]
    alpha <- par[3:4]
    lambda <- par[5]
    theta <- par[6]
    beta1 <- par[7]
    tempo <- t
    x <- as.matrix(x)
    X <- matrix(data = c(rep(1,nrow(x)),x), nrow = nrow(x), ncol = 2)
    #xalpha <- matrix(c(rep(1,nrow(X)),X), ncol = 2)
    st <- (1/(1 + exp(X%*%gamma) )) * (1 + ( (lambda*theta/(X%*%alpha)) * log( (1 + exp((X%*%alpha)*tempo + x%*%beta1) )/(1 + exp(x%*%beta1))) ) )^(-1/theta)
    
    return(st)
}


p0f <- function(x, beta) {
    
    x <- as.matrix(x)
    X <- matrix(data = c(rep(1,nrow(x)),x), nrow = nrow(x), ncol = 2)
    
    p <- exp(X%*%beta)/(1+exp(X%*%beta))
    return(p)
}

p1f <- function(x, par) {
    gamma <- par[1:2]
    alpha <- par[3:4]
    lambda <- par[5]
    theta <- par[6]
    beta1 <- par[7]
    
    x <- as.matrix(x)
    X <- matrix(data = c(rep(1,nrow(x)),x), nrow = nrow(x), ncol = 2)
    
    p = (1/(1+exp(X%*%gamma))) * (1 + ((lambda*theta)/(X%*%alpha))*log(1/(1+exp(x%*%beta1))) )^(-1/theta)
    return(p)
}

alfa <- function(x,par) x%*%par
froot <- function(t,x,par,u) sobrevGTDLv2(t,x,par) - 1 + u

# chutes iniciais
set.seed(2024)
n = 20000
#X = c(0,1)
X = rbinom(n, size = 1, prob = 0.5)
#(1,0,1,1) e c(-0.1327, -0.0058)
#alpha = alfa(x = matrix(c(1,0,1,1),ncol=2,byrow=T), par = c(-0.1327, -0.0058))
alpha = alfa(x = matrix(c(rep(1,2000),X), ncol = 2), par = c(-0.13269312, -0.00575169))
p0 = p0f(X, beta = c(-0.8534, -0.3231))
p1 = p1f(X, par = c(-0.8534, -0.3231, -0.13269312, -0.00575169, 0.4197, 0.6441, -0.0979))

# processo de geração dos tempos

t <- numeric()
delta <- numeric()

for (i in 1:20000) {
    
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
                                                           par = c(-0.85341178, -0.32305341, 
                                                                   -0.13269312, -0.00575169,
                                                                   0.41971047,  0.64408736, 
                                                                   -0.09790160),
                                                           interval = c(0,1080), 
                                                           extendInt = "yes")$root
                                  )
                     ),
                     no = ifelse(test = u1 < p0[i] , 
                                 yes = 0, 
                                 no = uniroot(f = froot, x=X[i], u = u3,
                                              par = c(-0.85341178, -0.32305341, 
                                                      -0.13269312, -0.00575169,
                                                      0.41971047,  0.64408736, -0.09790160), 
                                              interval = c(0,1080))$root)
        )
        tf <- as.numeric(tf)
        tc <- runif(1, 0, 20)
        t[i] <- min(tf,tc)
        
        delta[i] <- ifelse(t[i] == tf, 1, 0)

}


# processo de estimação dos parâmetros
X = rbinom(2000, size = 1, prob = 0.4709)
dadoSimul <- data.frame(cens = delta, tempo = t, X)
data0 <- data.frame(tempo = dadoSimul$t, x0 = rep(1,nrow(dadoSimul)),
                    X = dadoSimul$X)

emvg <- optim(c(-0.13269312, -0.00575169,
                0.41971047,  0.64408736, -0.09790160), veroGTDLv2, x=dadoSimul, hessian = TRUE, method = "BFGS");emvg
egpar <- c(emvg$par[1], emvg$par[2], exp(emvg$par[3]), exp(emvg$par[4]), emvg$par[5])
emv0 <- optim(par=c(2, 1), veroZerov2, x=data0, hessian = TRUE);emv0
e0par <- c(emv0$par[1], emv0$par[2])

epar <- c(e0par, egpar); epar

ekm <- survfit(Surv(time = tempo, event = delta) ~ X, data = dadoSimul)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     covar = c(rep(0,ekm$strata[[1]]),rep(1,ekm$strata[[2]])) )

ekm_df$gtdl0 <- sobrevGTDLv2(x = ekm_df[,c(1,3)], par = epar)

ggplot(data = ekm_df, mapping = aes(x = tempo, y = gtdl0)) +
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




# Model 0 v1 --------------------------------------------------------------



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


sobrevGTDLv1 <- function(t, x, par) {
    
    gamm <- par[1:2]
    alph <- par[3]
    lambd <- par[4]
    thet <- par[5]
    bet <- par[6]
    tempo <- t
    X <- as.matrix(x)
    X0 <- matrix(data = c(rep(1,nrow(X)),X), nrow = nrow(X), ncol = 2)
    st <- (1/(1 + exp(X0%*%gamm) )) * (1 + ( (lambd*thet/alph) * log( (1 + exp(alph*tempo + X%*%bet) )/(1 + exp(X%*%bet))) ) )^(-1/thet)
    
    return(st)
}


p0f <- function(x, beta) {
    
    x <- as.matrix(x)
    X <- matrix(data = c(rep(1,nrow(x)),x), nrow = nrow(x), ncol = 2)
    
    p <- exp(X%*%beta)/(1+exp(X%*%beta))
    return(p)
}

p1f <- function(x, par) {
    gamma <- par[1:2]
    alpha <- par[3]
    lambda <- par[4]
    theta <- par[5]
    beta1 <- par[6]
    
    x <- as.matrix(x)
    X <- matrix(data = c(rep(1,nrow(x)),x), nrow = nrow(x), ncol = 2)
    
    p = (1/(1+exp(X%*%gamma))) * (1 + ((lambda*theta)/alpha)*log(1/(1+exp(x%*%beta1))) )^(-1/theta)
    return(p)
}

froot <- function(t,x,par,u) sobrevGTDLv1(t,x,par) - 1 + u

# chutes iniciais
set.seed(2024)
n = 2000
#X = c(0,1)
X = rbinom(n, size = 1, prob = 0.5)
epar = c(-1.0293802, -0.3018615, -0.1403741,  0.4924747,  0.3988917, -0.9315740)
e0par = c(-1.0293802, -0.3018615)
#(1,0,1,1) e c(-0.1327, -0.0058)
#alpha = alfa(x = matrix(c(1,0,1,1),ncol=2,byrow=T), par = c(-0.1327, -0.0058))
alpha = rep(c(-0.1403741, 0.1403741),1000)
p0 = p0f(X, beta = e0par)
p1 = p1f(X, par = epar)

# processo de geração dos tempos

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
                                          par = epar, maxiter = 10000,
                                          interval = c(0,1080),
                                          extendInt = "yes")$root)
    )
    tf <- as.numeric(tf)
    tc <- runif(1, 0, 20)
    t[i] <- min(tf,tc)
    
    delta[i] <- ifelse(t[i] == tf, 1, 0)
    
}


# processo de estimação dos parâmetros
dadoSimul <- data.frame(cens = delta, tempo = t, X)
data0 <- data.frame(tempo = dadoSimul$t, x0 = rep(1,nrow(dadoSimul)),
                    X = dadoSimul$X)

emvg <- optim(par=c(0.1, -0.5, -0.5, 0.1), veroGTDLv1, x=dadoSimul, hessian = TRUE);emvg
egpar <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])
emv0 <- optim(par=c(-0.5, 0.5), veroZero, x=data0, hessian = TRUE);emv0
e0par <- c(emv0$par[1], emv0$par[2])

epar <- c(e0par, egpar); epar

ekm <- survfit(Surv(time = tempo, event = cens) ~ X, data = dadoSimul)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv,
                     covar = c(rep(0,ekm$strata[1]),rep(1,ekm$strata[2])) )

gtdl_df <- data.frame(tempo = rep(seq(0,60,1),2),
                      covar = c(rep(0,61),rep(1,61)))
gtdl_df$surv <- sobrevGTDLv1(t = gtdl_df$tempo ,x = gtdl_df$covar, par = epar)


ggplot(data = ekm_df, mapping = aes(x = tempo, y = km)) +
    geom_step(aes(color = as.factor(consulta))) +
    geom_line(data = gtdl_df, aes(x = tempo, y = surv, linetype = as.factor(covar))) +
    labs(x="Tempo (meses)", y="Estimativa de S(t)", 
         color="Inform. Consulta por Kaplan-Meier",
         linetype="Inform. Consulta por GTDL \nZero Inflacionado") +
    scale_color_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                       values = c("red","blue")) +
    scale_linetype_manual(labels = c("0 - Com consulta", "1 - Sem consulta"),
                          values = c(1,2)) +
    ylim(0:1) + scale_x_continuous(breaks = seq(0,90,5)) +
    theme_light() + 
    theme(legend.position = c(0.7,0.7),
          legend.background = element_rect(color = "white"))
