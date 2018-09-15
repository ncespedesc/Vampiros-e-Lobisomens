###################################################################################
#                    modelagem matematica lobisomem vampiro humanos               #
###################################################################################


# Formulacao com fracoes (proporcoes) 

# Carregando pacotes ----
library(deSolve)
library(tidyverse)
library(reshape2)

# Parametros ----

# Taxas por ano 

natalidade <- 0 #taxa matalidade mundial 2015
Motalidade.Hum <- 0.0001 #taxa de mortalidade
beta.vamp <- 0.001 #taxa de infecao vampiros
beta.lobi <- 0.001 #taxa de infecao lobisomn
prev.vamp <- 0.0001 # taxa de prevencao que um humano vire vamp porque se suicida ou le disparam 
prev.lobi <- 0.0001 #taxa de prevencao que um humano vire lobi porque se suicida ou le disparam 
gamma.vamp  <- 0#365/21 # periodo de latencia para virar vampiro 3 sem
gamma.lobi <- 0 #365/21 # periodo de latencia lobi 3 sem 
letha.homen.mata.vampiro <- 0 # taxa de humonaos que matam vampiros 
letha.lobi.mata.vampiro <- 0 #taxa de lobi que matam vampiros 
letha.homenm.mata.lobi <- 0 # taxa de humonaos que matam lobi
letha.vampi.mata.lobi <- 0 # taxa de vampiro que matam lobi
mortalidade.lobi <- 0 # taxa de morte natural dos lobi

# o desolve precisa um conjunto de parametrros pra souber o nome da equacao 

par.SVW <- c(
  natalidade = natalidade,
  Motalidade.Hum = Motalidade.Hum,
  beta.vamp = beta.vamp,
  beta.lobi = beta.lobi,
  prev.vamp = prev.vamp,
  prev.lobi = prev.lobi
  # gamma.vamp =  gamma.vamp,
  # gamma.lobi = gamma.lobi,
  
  # gamma.vamp =  gamma.lobi,
  # gamma.lobi = gamma.vamp,
  
  # letha.homen.mata.vampiro = letha.homen.mata.vampiro,
  # letha.lobi.mata.vampiro = letha.lobi.mata.vampiro,
  # letha.homenm.mata.lobi = letha.homenm.mata.lobi,
  # letha.vampi.mata.lobi = letha.vampi.mata.lobi,
  # mortalidade.lobi = mortalidade.lobi
)

# Calculando R0 ----
#
R0 <- beta.vamp/(gamma.vamp)
R0 

# Variaveis e condicao inicial ----

S  <- 1000
Iv <- 0
Iw <- 0
V <- 1
W <- 1



# state.SIR <- c(s=0.9999,i=0.0001,r=0)

state.SVW <- c(S = S, Iv = Iv,Iw = Iw)



# Tempo de simulacao ----
# Varie o tempo de simulacao: 20, 100
tsim <- 1000
Dt <- 1

# Funcao para o modelo SIR ----
# Termo de transmissao frequencia dependente 
SIRS <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # # rate of change
    # ds <- mu - beta*s*i - mu*s + + omegaR*r
    # di <- beta*s*i - (mu+gama)*i
    # dr <- gama*i - mu*r - omegaR*r
    
    ds <- natalidade*S - Motalidade.Hum*S - beta.vamp *S*V -beta.lobi*S*W
    
    dIv <- beta.vamp*S*V - prev.vamp* Iv 
    # dV <- gamma.vamp*Iv - letha.homen.mata.vampiro*V - letha.lobi.mata.vampiro*V
    
    dIw <- beta.lobi*S*W - prev.lobi*Iw 
    # dW <- gamma.lobi*Iw - letha.homenm.mata.lobi*W - letha.vampi.mata.lobi*W - mortalidade.lobi*W 
    
    
    # return the output of the model
    return(list(c(ds, dIv, dIw)))
    
  })
}

tempos <- seq(from=0,to=tsim,by=Dt)

# modSIRS <- ode(y = state.SIR, times = tempos, func = SIRS, parms = par.SIRS, method = "ode45")
modSIRS <- ode(y = state.SVW, times = tempos, func = SIRS, parms = par.SVW, method = "ode45")


modSIRS <- as.data.frame(modSIRS)

modSIRS %>%
  gather(key = 'compartimento', value = 'valor', -time)%>%
  ggplot(aes(x= time, y = valor))+
  geom_line()+
  facet_wrap('compartimento', scales = 'fixed')
