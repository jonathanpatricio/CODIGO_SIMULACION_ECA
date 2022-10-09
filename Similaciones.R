#Funciones----------------------------------------------------------------------
Pob_ECA <- function(n, t, b0, b1, b2, b3, b4, b5, var_v0i, var_v1i, cov_v0iv1i, var_eij) {
  # Semilla para fijar los números aleatorios
  set.seed(16)
  
  ## Librerías
  library(MASS)
  library(dplyr)
  
  ## Parámetros del modelo
  
  # Tratamientos
  treat <- as.numeric(as.numeric(cbind(rep(0, n), rep(1, n))))
  
  treat_2 <- as.numeric(cbind(rep(0, n), rep(1, n), rep(0, n)))
  treat_3 <- as.numeric(cbind(rep(0, n), rep(0, n), rep(1, n)))
  
  
  # Tamaño total
  n_2 <- n*2
  n_3 <- n*3
  
  # t_ij
  t_j <- ((1:t)-1)/(t-1)
  tij <- matrix(rep(t_j), nrow=n_2, ncol=t, byrow = TRUE)
  tij_ <- matrix(rep(t_j), nrow=n_3, ncol=t, byrow = TRUE)
  
  # Errores eij
  cov_eij <- 0 # Covarianza de los errores
  mu_eij <- rep(0, t) # promedio de los errores
  sigma_eij <- matrix(rep(cov_eij, t*t), nrow=t, ncol=t); diag(sigma_eij) <- rep(var_eij, t) # Matriz de varianzas y covarianzas
  
  eij <- mvrnorm(n_2, mu_eij, sigma_eij) #Errores de acuerdo con los parámetros establecidos
  eij_ <- mvrnorm(n_3, mu_eij, sigma_eij) #Errores de acuerdo con los parámetros establecidos
  #plot(eij, pch='.')
  
  # Efectos aleatorios
  mu_v0iv1i  <- rep(0, 2) # promedio de los efectos aleatorios
  sigma_v0iv1i  <- matrix(rep(cov_v0iv1i, 4), nrow=2, ncol=2); diag(sigma_v0iv1i) <- c(var_v0i, var_v1i) # Matriz de varianzas y covarianzas
  
  v_i <- mvrnorm(n_2, mu_v0iv1i, sigma_v0iv1i) #efectos aleatorios de acuerdo con los parámetros establecidos
  v_i_ <- mvrnorm(n_3, mu_v0iv1i, sigma_v0iv1i) #efectos aleatorios de acuerdo con los parámetros establecidos
  
  # Efecto aleatorio del intercepto v_0i
  v0i <- v_i[,1]
  v0i_ <- v_i_[,1]
  
  # Efecto aleatorio de la pendiente v_1i
  v1i <- v_i[,2]
  v1i_ <- v_i_[,2]
  
  # Beta_0
  b0_ <- b0 
  b0 <- b0 
  
  # Beta_1
  b1_ <- b1 
  b1 <- b1 
  
  # Beta_2
  b2_ <- b2 * treat_2
  b2 <- b2 * treat
  
  # Beta_3
  b3_ <- b3 * treat_2
  b3 <- b3 * treat
  
  # Beta_4
  b4_ <- b4 * treat_3
  
  # Beta_5
  b5_ <- b5 * treat_3
  
  # estimación y_ij a lo ancho
  yij_2_treat <- as.data.frame(b0   + b2  +       v0i  + b1  * tij   + b3  * tij   +              v1i  * tij  + eij  )
  yij_3_treat <- as.data.frame(b0_  + b2_ + b4_ + v0i_ + b1_ * tij_  + b3_ * tij_  + b5_ * tij_ + v1i_ * tij_ + eij_ )
  
  grupos <- as.numeric(as.numeric(cbind(rep(1, n), rep(2, n))))
  yij_2_treat <- mutate(yij_2_treat, treat = as.factor(grupos), ID = 1:n_2)
  
  grupos_ <- as.factor(cbind(rep(1, n), rep(2, n), rep(3, n)))
  yij_3_treat <- mutate(yij_3_treat, treat = grupos_, ID = 1:n_3)
  
  # estimación y_ij a lo largo
  yij_2_treat_long <- reshape(data = yij_2_treat,varying = 1:t, v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
  yij_2_treat_long <- arrange(yij_2_treat_long,ID,tiempo)
  yij_2_treat_long$tiempo <- as.numeric(yij_2_treat_long$tiempo)
  yij_2_treat_long$tiempo <- (yij_2_treat_long$tiempo-1)/(t-1)
  
  
  yij_3_treat_long <- reshape(data = yij_3_treat,varying = 1:t, v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
  yij_3_treat_long <- arrange(yij_3_treat_long,ID,tiempo)
  yij_3_treat_long$tiempo <- as.numeric(yij_3_treat_long$tiempo)
  yij_3_treat_long$tiempo <- (yij_3_treat_long$tiempo-1)/(t-1)
  
  # Matriz de correlación 
  cor_2_treat <- cor(yij_2_treat[,1:t])
  cor_3_treat <- cor(yij_3_treat[,1:t])
  
  # Listado de resultados
  resul_2_treat <- list(yij_2 = yij_2_treat, yij_long_2 = yij_2_treat_long, v_i_2 = v_i, eij_2 = eij, cor_2 = cor_2_treat)
  resul_3_treat <- list(yij_3 = yij_3_treat, yij_long_3 = yij_3_treat_long, v_i_3 = v_i_, eij_3 = eij_, cor_3 = cor_3_treat)
  
  return(list(Treat_2 = resul_2_treat, Treat_3 = resul_3_treat ))
  
}
var_cov_v0iv1i <- function(sd_v0i, sd_v1i, cor_v0iv1i) {
  var_v0i <- sd_v0i * sd_v0i  
  var_v1i <- sd_v1i * sd_v1i
  var_cov_v0iv1i <- cor_v0iv1i * ( sd_v0i * sd_v1i)
  return(list(var_v0i = var_v0i,  var_v1i = var_v1i, var_cov_v0iv1i = var_cov_v0iv1i))
}
Comp_poder_2_treat <- function(yij_2_treat, sample_min, sample_max, repeticiones, t, k) {
  
  #Semilla
  set.seed(16)
  
  #Librerías
  library(dplyr)
  library(gee)
  library(ggplot2)
  library(lmerTest)
  library(lme4)
  
  #Matrices
  Gee_intercanbiable <- matrix(0,(sample_max-sample_min),repeticiones)
  Gee_AR1            <- matrix(0,(sample_max-sample_min),repeticiones)
  Gee_unstructured   <- matrix(0,(sample_max-sample_min),repeticiones)
  Mixto_intercepto   <- matrix(0,(sample_max-sample_min),repeticiones)
  Mixto_pen_inter    <- matrix(0,(sample_max-sample_min),repeticiones)
  
  #Bucle
  for (i in (sample_min:sample_max)*k) {for(j in 1:repeticiones){
    
    sample_treat_1 <- yij_2_treat[sample(x = 1:(nrow(yij_2_treat)/2), size = i, replace = FALSE),]
    sample_treat_2 <- yij_2_treat[sample(x = (nrow(yij_2_treat)/2+1):nrow(yij_2_treat), size = i, replace = FALSE),]
    sample <- bind_rows(sample_treat_1,sample_treat_2)
    
    sample_long <- reshape(data = sample,varying = 1:t, v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
    sample_long <- arrange(sample_long,ID,tiempo)
    sample_long$tiempo <- as.numeric(sample_long$tiempo)
    sample_long$tiempo <- (sample_long$tiempo-1)/(t-1)
    
    #Modelos
    intercanbiable <- gee(yij ~ treat + tiempo + treat * tiempo, id = ID, data = sample_long, family = gaussian, corstr = "exchangeable")
    AR1            <- gee(yij ~ treat + tiempo + treat * tiempo, id = ID, data = sample_long, family = gaussian, corstr = "AR-M", Mv = 1)
    unstructured   <- gee(yij ~ treat + tiempo + treat * tiempo, id = ID, data = sample_long, family = gaussian, corstr = "unstructured")
    intercepto     <- lmer(yij ~ treat + tiempo + treat * tiempo + (1|ID), data = sample_long, REML = FALSE)
    pen_intercepto <- lmer(yij ~ treat + tiempo + treat * tiempo + (tiempo|ID), data = sample_long, REML = FALSE)
    
    #Completando las matrices con la decisión de la hipótesis
    Gee_intercanbiable[(i/k)-(sample_min),j] <-if(pnorm(as.matrix(intercanbiable$coefficients)[4,]/sqrt(intercanbiable$robust.variance[4,4]), 0, 1)        < 0.05) 1 else 0
    Gee_AR1[(i/k)-(sample_min),j]            <-if(pnorm(as.matrix(AR1$coefficients)[4,]/sqrt(AR1$robust.variance[4,4]), 0, 1)                              < 0.05) 1 else 0
    Gee_unstructured[(i/k)-(sample_min),j]   <-if(pnorm(as.matrix(unstructured$coefficients)[4,]/sqrt(unstructured$robust.variance[4,4]), 0, 1)            < 0.05) 1 else 0
    Mixto_intercepto[(i/k)-(sample_min),j]   <-if(pt(as.matrix(intercepto@beta)[4,]/sqrt(intercepto@vcov_beta[4,4]), df = as.matrix(intercepto@Gp)[2,]-1)     < 0.05) 1 else 0
    Mixto_pen_inter[(i/k)-(sample_min),j]    <-if(pt(as.matrix(pen_intercepto@beta)[4,]/sqrt(pen_intercepto@vcov_beta[4,4]), df = as.matrix(pen_intercepto@Gp)[2,]-1)     < 0.05) 1 else 0
    
  }}
  
  #Base de datos
  Gee_inter  <- as.matrix(apply(X = Gee_intercanbiable, MARGIN = 1, FUN = mean))
  Gee_AR     <- as.matrix(apply(X = Gee_AR1,            MARGIN = 1, FUN = mean))
  Gee_unst   <- as.matrix(apply(X = Gee_unstructured,   MARGIN = 1, FUN = mean))
  Mixto_inte <- as.matrix(apply(X = Mixto_intercepto,   MARGIN = 1, FUN = mean))
  Mixto_pen_ <- as.matrix(apply(X = Mixto_pen_inter,    MARGIN = 1, FUN = mean))
  
  Base <- as.data.frame(cbind(ID = (sample_min:(sample_max-1)*k),Gee_inter))
  Base <- mutate(Base,Gee_AR)
  Base <- mutate(Base,Gee_unst)
  Base <- mutate(Base,Mixto_inte)
  Base <- mutate(Base,Mixto_pen_)
  
  #Gráfica
  Base_largo <- reshape(data = Base, varying = 2:6, v.names = "Poder", timevar= "Modelo", idvar = "ID", direction = "long")
  colnames(Base_largo) <- c("n","modelo","Acepta_HO")
  Base_largo<-arrange(Base_largo,n,modelo)
  Base_largo$modelo <- factor(Base_largo$modelo, labels = c("Gee exchangeable", "Gee AR(1)", "Gee unstructured","Mixto intercepto aleatorio", "Mixto intercepto y pendiente aleatoria"))
  
  Grafico <- ggplot(data = Base_largo, aes(x = n, y = Acepta_HO, color = modelo)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = loess, se = FALSE) +
    theme_classic() +
    labs(title="P", y="Poder", x="Tamaño de muestra", caption="Fuente: Simulación", size = 2 )
  
  return(list(Base = Base, Gee_intercanbiable = Gee_intercanbiable, Gee_AR1 = Gee_AR1,
              Gee_unstructured = Gee_unstructured, Mixto_intercepto = Mixto_intercepto,
              Mixto_pen_inter = Mixto_pen_inter, Base_largo = Base_largo, 
              Grafico = Grafico ))
  
}

var_cov_v0iv1i(sd_v0i = 10.1, 
               sd_v1i = 8, 
               cor_v0iv1i = -0.75)


# Escenario 7
# Población
población <- Pob_ECA(n = 1000000, 
                     t = 4, 
                     b0 = 44.9, 
                     b1 = 0.11, 
                     b2 = -2.3, 
                     b3 = 0, 
                     b4 = -2.3, 
                     b5 = 0, 
                     var_v0i = 102.01, 
                     var_v1i = 64, 
                     cov_v0iv1i = -60.6, 
                     var_eij = 22.09)

yij_2_treat <-población$Treat_2$yij_2


#Escenarios --------------------------------------------------------------------

# Escenario_7.1
Escenario_7.1 <- Comp_poder_2_treat(yij_2_treat = yij_2_treat, 
                                    sample_min = 10, 
                                    sample_max = 30, 
                                    repeticiones = 2500, 
                                    t = 4,
                                    k = 1)

library(xlsx)
write.table(Escenario_7.1$Base, "Datos/Escenario 7/Base.txt")
write.table(Escenario_7.1$Gee_intercanbiable, "Datos/Escenario 7/Gee_intercanbiable.txt")
write.table(Escenario_7.1$Gee_AR1, "Datos/Escenario 7/Gee_AR1.txt")
write.table(Escenario_7.1$Gee_unstructured, "Datos/Escenario 7/Gee_unstructured.txt")
write.table(Escenario_7.1$Mixto_intercepto, "Datos/Escenario 7/Mixto_intercepto.txt")
write.table(Escenario_7.1$Mixto_pen_inter, "Datos/Escenario 7/Mixto_pen_inter.txt")


# Escenario_7.2
Escenario_7.2 <- Comp_poder_2_treat(yij_2_treat = yij_2_treat, 
                                    sample_min = 3, 
                                    sample_max = 6, 
                                    repeticiones = 2500, 
                                    t = 4,
                                    k = 10)


write.table(Escenario_7.2$Base, "Datos/Escenario 7/Base.txt")
write.table(Escenario_7.2$Gee_intercanbiable, "Datos/Escenario 7/Gee_intercanbiable.txt")
write.table(Escenario_7.2$Gee_AR1, "Datos/Escenario 7/Gee_AR1.txt")
write.table(Escenario_7.2$Gee_unstructured, "Datos/Escenario 7/Gee_unstructured.txt")
write.table(Escenario_7.2$Mixto_intercepto, "Datos/Escenario 7/Mixto_intercepto.txt")
write.table(Escenario_7.2$Mixto_pen_inter, "Datos/Escenario 7/Mixto_pen_inter.txt")


# Escenario_7.3
Escenario_7.3 <- Comp_poder_2_treat(yij_2_treat = yij_2_treat, 
                                    sample_min = 3.5, 
                                    sample_max = 5.5, 
                                    repeticiones = 2500, 
                                    t = 4,
                                    k = 10)


write.table(Escenario_7.3$Base, "Datos/Escenario 7/Base.txt")
write.table(Escenario_7.3$Gee_intercanbiable, "Datos/Escenario 7/Gee_intercanbiable.txt")
write.table(Escenario_7.3$Gee_AR1, "Datos/Escenario 7/Gee_AR1.txt")
write.table(Escenario_7.3$Gee_unstructured, "Datos/Escenario 7/Gee_unstructured.txt")
write.table(Escenario_7.3$Mixto_intercepto, "Datos/Escenario 7/Mixto_intercepto.txt")
write.table(Escenario_7.3$Mixto_pen_inter, "Datos/Escenario 7/Mixto_pen_inter.txt")


# Escenario_7.4
Escenario_7.4 <- Comp_poder_2_treat(yij_2_treat = yij_2_treat, 
                                    sample_min = 5, 
                                    sample_max = 10, 
                                    repeticiones = 2500, 
                                    t = 4,
                                    k = 10)


write.table(Escenario_7.4$Base, "Datos/Escenario 7/Base.txt")
write.table(Escenario_7.4$Gee_intercanbiable, "Datos/Escenario 7/Gee_intercanbiable.txt")
write.table(Escenario_7.4$Gee_AR1, "Datos/Escenario 7/Gee_AR1.txt")
write.table(Escenario_7.4$Gee_unstructured, "Datos/Escenario 7/Gee_unstructured.txt")
write.table(Escenario_7.4$Mixto_intercepto, "Datos/Escenario 7/Mixto_intercepto.txt")
write.table(Escenario_7.4$Mixto_pen_inter, "Datos/Escenario 7/Mixto_pen_inter.txt")


# Escenario_7.5
Escenario_7.5 <- Comp_poder_2_treat(yij_2_treat = yij_2_treat, 
                                    sample_min = 1, 
                                    sample_max = 11, 
                                    repeticiones = 2500, 
                                    t = 4,
                                    k = 100)


write.table(Escenario_7.5$Base, "Datos/Escenario 7/Base.txt")
write.table(Escenario_7.5$Gee_intercanbiable, "Datos/Escenario 7/Gee_intercanbiable.txt")
write.table(Escenario_7.5$Gee_AR1, "Datos/Escenario 7/Gee_AR1.txt")
write.table(Escenario_7.5$Gee_unstructured, "Datos/Escenario 7/Gee_unstructured.txt")
write.table(Escenario_7.5$Mixto_intercepto, "Datos/Escenario 7/Mixto_intercepto.txt")
write.table(Escenario_7.5$Mixto_pen_inter, "Datos/Escenario 7/Mixto_pen_inter.txt")

# Escenario_7.6
Escenario_7.6 <- Comp_poder_2_treat(yij_2_treat = yij_2_treat, 
                                    sample_min = 1, 
                                    sample_max = 11, 
                                    repeticiones = 2500, 
                                    t = 4,
                                    k = 100)


write.table(Escenario_7.6$Base, "Datos/Escenario 7/Base.txt")
write.table(Escenario_7.6$Gee_intercanbiable, "Datos/Escenario 7/Gee_intercanbiable.txt")
write.table(Escenario_7.6$Gee_AR1, "Datos/Escenario 7/Gee_AR1.txt")
write.table(Escenario_7.6$Gee_unstructured, "Datos/Escenario 7/Gee_unstructured.txt")
write.table(Escenario_7.6$Mixto_intercepto, "Datos/Escenario 7/Mixto_intercepto.txt")
write.table(Escenario_7.6$Mixto_pen_inter, "Datos/Escenario 7/Mixto_pen_inter.txt")




#graficos de comparación--------------------------------------------------------
plot(Datos.unificados$Gee_int ~ Datos.unificados$ID, type = "b",          pch = "1", ylim = c(0.04, 0.12))
lines(Escenario_7.5$Base$Mixto_inte ~ Escenario_7.5$Base$ID, type = "b", pch = "2", col = 2)
lines(Escenario_7.5$Base$Gee_AR ~ Escenario_7.5$Base$ID, type = "b",     pch = "3", col = 3)
lines(Escenario_7.5$Base$Gee_unst ~ Escenario_7.5$Base$ID, type = "b",   pch = "4", col = 4)
lines(Escenario_7.5$Base$Mixto_pen_ ~ Escenario_7.5$Base$ID, type = "b", pch = "5", col = 5)

library(ggplot2)
ggplot(data = Datos.unificados, aes(x = log(ID), y = Gee_int )) +
  geom_line()





Escenario_5.2$Grafico +
  geom_hline(yintercept = 0.05)
