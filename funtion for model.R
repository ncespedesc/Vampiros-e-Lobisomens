###################################################################################
#                    modelagem matematica lobisomem vampiro humanos               #
##################################################################################

# Parametros ----

natalidade <- 0.001 #taxa matalidade mundial 2015
Motalidade.Hum <- 0.0001 #taxa de mortalidade
beta.vamp <- 0.001 #taxa de infecao vampiros
beta.lobi <- 0.001 #taxa de infecao lobisomn
prev.vamp <- 0.001 # taxa de prevencao que um humano vire vamp porque se suicida ou le disparam 
prev.lobi <- 0.001 #taxa de prevencao que um humano vire lobi porque se suicida ou le disparam 
# gamma.vamp  <- 0#365/21 # periodo de latencia para virar vampiro 3 sem
# gamma.lobi <- 0 #365/21 # periodo de latencia lobi 3 sem 
letha.homen.mata.vampiro <- 0.1 # taxa de humonaos que matam vampiros 
letha.lobi.mata.vampiro <- 0.1 #taxa de lobi que matam vampiros 
letha.homenm.mata.lobi <- 0.1 # taxa de humonaos que matam lobi
letha.vampi.mata.lobi <- 0.1 # taxa de vampiro que matam lobi
mortalidade.lobi <- 0.1

#funtion form model ----

#paramoteros do modelo (correr funcao primeiro )

HVW.model(natalidade ,
          Motalidade.Hum ,
          beta.vamp ,
          beta.lobi ,
          prev.vamp ,
          prev.lobi , 
          # gamma.vamp =  gamma.vamp,
          # gamma.lobi = gamma.lobi,
          letha.homen.mata.vampiro ,
          letha.lobi.mata.vampiro ,
          letha.homenm.mata.lobi ,
          letha.vampi.mata.lobi ,
          mortalidade.lobi )


# funtion para o modelo ----

HVW.model  <- function(natalidade = natalidade,
                       Motalidade.Hum = Motalidade.Hum,
                       beta.vamp = beta.vamp,
                       beta.lobi = beta.lobi,
                       prev.vamp = prev.vamp,
                       prev.lobi = prev.lobi, 
                       # gamma.vamp =  gamma.vamp,
                       # gamma.lobi = gamma.lobi,
                       letha.homen.mata.vampiro = letha.homen.mata.vampiro,
                       letha.lobi.mata.vampiro = letha.lobi.mata.vampiro,
                       letha.homenm.mata.lobi = letha.homenm.mata.lobi,
                       letha.vampi.mata.lobi = letha.vampi.mata.lobi,
                       mortalidade.lobi = mortalidade.lobi) {

# Formulacao com fracoes (proporcoes) 

# Carregando pacotes ----
  
  if(!(require(deSolve))){install.packages("deSolve")}; library(deSolve)
  if(!(require(tidyverse))){install.packages("tidyverse")}; library(tidyverse)
  if(!(require(reshape2))){install.packages("reshape2")}; library(reshape2)
  if(!(require(ggthemes))){install.packages("ggthemes")}; library(ggthemes)




# o desolve precisa um conjunto de parametrros pra souber o nome da equacao 

par.SVW <- c(
  natalidade = natalidade,
  Motalidade.Hum = Motalidade.Hum,
  beta.vamp = beta.vamp,
  beta.lobi = beta.lobi,
  prev.vamp = prev.vamp,
  prev.lobi = prev.lobi, 
  # gamma.vamp =  gamma.vamp,
  # gamma.lobi = gamma.lobi,
  letha.homen.mata.vampiro = letha.homen.mata.vampiro,
  letha.lobi.mata.vampiro = letha.lobi.mata.vampiro,
  letha.homenm.mata.lobi = letha.homenm.mata.lobi,
  letha.vampi.mata.lobi = letha.vampi.mata.lobi,
  mortalidade.lobi = mortalidade.lobi
)

# Calculando R0 ----
#
# R0 <- beta.vamp/(gamma.vamp)
# R0 

# Variaveis e condicao inicial ----

S  <- 1000
Iv <- 0
Iw <- 0
V <- 2
W <- 2



# state.SIR <- c(s=0.9999,i=0.0001,r=0)

state.SVW <- c(S = S, Iv = Iv,Iw = Iw)



# Tempo de simulacao ----

tsim <- 1000
Dt <- 1

# Funcao para o modelo SIR ----
# Termo de transmissao frequencia dependente 
SIRS <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # # rate of change
    
    ds <- natalidade*S - Motalidade.Hum*S - beta.vamp *S*V -beta.lobi*S*W
    
    dIv <- beta.vamp*S*V - prev.vamp* Iv - letha.lobi.mata.vampiro*V - letha.homen.mata.vampiro*V 
    
    
    dIw <- beta.lobi*S*W - prev.lobi*Iw - letha.vampi.mata.lobi*W - letha.homenm.mata.lobi*W - mortalidade.lobi*W 
    
    # return the output of the model
    return(list(c(ds, dIv, dIw)))
    
  })
}

tempos <- seq(from=0,to=tsim,by=Dt)

# modSIRS <- ode(y = state.SIR, times = tempos, func = SIRS, parms = par.SIRS, method = "ode45")
modSIRS <- ode(y = state.SVW, times = tempos, func = SIRS, parms = par.SVW, method = "ode45")


modSIRS <- as.data.frame(modSIRS)
names(modSIRS) <- c("Time", "Humans", "Vampires", "Wolfman")



plot.vamp.lobi.hum <- modSIRS %>%
                        gather(key = 'Population', value = 'valor', -Time)  %>% 
                        ggplot( )+
                        geom_line(aes(x= Time, y = valor, colour = Population), size = 1)+
                        ylim(0,1000)+
                        # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
                        # ggthemes::theme_igray()+
                        xlab("Population ")+ ylab("Time (years)")+ ggtitle("Interration of Populations")+
                        theme(axis.text.x = element_text(angle = 45, hjust = 1) ,text = element_text(size = 17, face = "bold") )
                      
                      



return(plot.vamp.lobi.hum)

}

