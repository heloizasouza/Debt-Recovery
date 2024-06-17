

# buscar comando de carregar as funções do script Main


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










