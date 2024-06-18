

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



# Model Final -------------------------------------------------------------




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
X = c(0,1)
alpha = alfa(x = matrix(c(1,0,1,1),ncol=2,byrow=T), par = c(-0.1327, -0.0058))
p0 = p0f(x = X, beta = c(-0.8534, -0.3231))
p1 = p1f(x = X, par = c(-0.8534, -0.3231, -0.1327, -0.0058, 0.4197, 0.6441, -0.0979))

t <- numeric()
delta <- numeric()

for (i in 1:2000) {
    
    for (j in 1:2) {
        u1 <- runif(1,0,1)
        u2 <- runif(1, min = p0[j], max = 1-p1[j])
        u3 <- runif(n = 1, min = p0[j], max = 1)
        
        tf <- ifelse(test = alpha[j] < 0,
                     yes = ifelse(test = u1 < p0[j],
                                  yes = 0,
                                  no = ifelse(test = 1 - p1[j] < u1,
                                              yes = Inf,
                                              no = uniroot(f = froot, x=X[j], u=u2,
                                                           par = c(-0.85341178, -0.32305341, 
                                                                   -0.13269312, -0.00575169,
                                                                   0.41971047,  0.64408736, 
                                                                   -0.09790160),
                                                           interval = c(0,1080), 
                                                           extendInt = "yes")$root
                                  )
                     ),
                     no = ifelse(test = u1 < p0[j] , 
                                 yes = 0, 
                                 no = uniroot(f = froot, x=X[j], u = u3,
                                              par = c(-0.85341178, -0.32305341, 
                                                      -0.13269312, -0.00575169,
                                                      0.41971047,  0.64408736, -0.09790160), 
                                              interval = c(0,1080), extendInt = "yes")$root)
        )
        tf <- as.numeric(tf)
        tc <- runif(1, 0, 1080)
        t[i] <- min(tf,tc)
        
        delta[i] <- ifelse(t[i] == tf, 1, 0)
    }
    
}


