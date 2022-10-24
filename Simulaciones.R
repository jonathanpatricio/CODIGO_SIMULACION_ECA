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
Comp_2_treat <- function(yij_2_treat, sample_min, sample_max, repeticiones, t, k) {
  
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
  
  Base <- as.data.frame(cbind(ID = (sample_min:(sample_max-1)*k)))
  Base <- mutate(Base,Gee_inter)
  Base <- mutate(Base,Gee_AR)
  Base <- mutate(Base,Gee_unst)
  Base <- mutate(Base,Mixto_inte)
  Base <- mutate(Base,Mixto_pen_)
  
  #Gráfica
  #Base_largo <- reshape(data = Base, varying = 2:6, v.names = "Poder", timevar= "Modelo", idvar = "ID", direction = "long")
  #colnames(Base_largo) <- c("n","modelo","Acepta_HO")
  #Base_largo<-arrange(Base_largo,n,modelo)
  #Base_largo$modelo <- factor(Base_largo$modelo, labels = c("Gee exchangeable", "Gee AR(1)", "Gee unstructured","Mixto intercepto aleatorio", "Mixto intercepto y pendiente aleatoria"))
  
  #Grafico <- ggplot(data = Base_largo, aes(x = n, y = Acepta_HO, color = modelo)) +
   # geom_point(alpha = 0.3, size = 1) +
    #geom_smooth(method = loess, se = FALSE) +
    #theme_classic() +
    #labs(title="P", y="Poder", x="Tamaño de muestra", caption="Fuente: Simulación", size = 2 )
  
  return(list(Base = Base, Gee_intercanbiable = Gee_intercanbiable, Gee_AR1 = Gee_AR1,
              Gee_unstructured = Gee_unstructured, Mixto_intercepto = Mixto_intercepto,
              Mixto_pen_inter = Mixto_pen_inter))#Base_largo = Base_largo, 
              #Grafico = Grafico ))
  
}
Comp_3_treat <- function(yij_3_treat, sample_min, sample_max, repeticiones, t, k) {
  
  #Semilla
  set.seed(16)
  
  #Librerías
  library(dplyr)
  library(gee)
  library(ggplot2)
  library(lmerTest)
  library(lme4)
  
  #Matrices
  Gee_intercanbiable_treat2 <- matrix(0,(sample_max-sample_min),repeticiones)
  Gee_AR1_treat2            <- matrix(0,(sample_max-sample_min),repeticiones)
  Gee_unstructured_treat2   <- matrix(0,(sample_max-sample_min),repeticiones)
  Mixto_intercepto_treat2   <- matrix(0,(sample_max-sample_min),repeticiones)
  Mixto_pen_inter_treat2    <- matrix(0,(sample_max-sample_min),repeticiones)
  
  Gee_intercanbiable_treat3 <- matrix(0,(sample_max-sample_min),repeticiones)
  Gee_AR1_treat3            <- matrix(0,(sample_max-sample_min),repeticiones)
  Gee_unstructured_treat3   <- matrix(0,(sample_max-sample_min),repeticiones)
  Mixto_intercepto_treat3   <- matrix(0,(sample_max-sample_min),repeticiones)
  Mixto_pen_inter_treat3    <- matrix(0,(sample_max-sample_min),repeticiones)
  
  #Bucle
  for (i in (sample_min:sample_max)*k) {for(j in 1:repeticiones){
    
    sample_treat_1 <- yij_3_treat[sample(x = 1:(nrow(yij_3_treat)/3), size = i, replace = FALSE),]
    sample_treat_2 <- yij_3_treat[sample(x = (nrow(yij_3_treat)/3+1):((nrow(yij_3_treat)/3)*2), size = i, replace = FALSE),]
    sample_treat_3 <- yij_3_treat[sample(x = ((nrow(yij_3_treat)/3)*2+1):nrow(yij_3_treat), size = i, replace = FALSE),]
    
    sample <- bind_rows(sample_treat_1,sample_treat_2,sample_treat_3)
    
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
    Gee_intercanbiable_treat2[(i/k)-(sample_min),j] <-if(pnorm(as.matrix(intercanbiable$coefficients)[5,]/sqrt(intercanbiable$robust.variance[5,5]), 0, 1)        < 0.05) 1 else 0
    Gee_AR1_treat2[(i/k)-(sample_min),j]            <-if(pnorm(as.matrix(AR1$coefficients)[5,]/sqrt(AR1$robust.variance[5,5]), 0, 1)                              < 0.05) 1 else 0
    Gee_unstructured_treat2[(i/k)-(sample_min),j]   <-if(pnorm(as.matrix(unstructured$coefficients)[5,]/sqrt(unstructured$robust.variance[5,5]), 0, 1)            < 0.05) 1 else 0
    Mixto_intercepto_treat2[(i/k)-(sample_min),j]   <-if(pt(as.matrix(intercepto@beta)[5,]/sqrt(intercepto@vcov_beta[5,5]), df = as.matrix(intercepto@Gp)[2,]-1)     < 0.05) 1 else 0
    Mixto_pen_inter_treat2[(i/k)-(sample_min),j]    <-if(pt(as.matrix(pen_intercepto@beta)[5,]/sqrt(pen_intercepto@vcov_beta[5,5]), df = as.matrix(pen_intercepto@Gp)[2,]-1)     < 0.05) 1 else 0
    
    Gee_intercanbiable_treat3[(i/k)-(sample_min),j] <-if(pnorm(as.matrix(intercanbiable$coefficients)[6,]/sqrt(intercanbiable$robust.variance[6,6]), 0, 1)        < 0.05) 1 else 0
    Gee_AR1_treat3[(i/k)-(sample_min),j]            <-if(pnorm(as.matrix(AR1$coefficients)[6]/sqrt(AR1$robust.variance[6,6]), 0, 1)                              < 0.05) 1 else 0
    Gee_unstructured_treat3[(i/k)-(sample_min),j]   <-if(pnorm(as.matrix(unstructured$coefficients)[6,]/sqrt(unstructured$robust.variance[6,6]), 0, 1)            < 0.05) 1 else 0
    Mixto_intercepto_treat3[(i/k)-(sample_min),j]   <-if(pt(as.matrix(intercepto@beta)[6,]/sqrt(intercepto@vcov_beta[6,6]), df = as.matrix(intercepto@Gp)[2,]-1)     < 0.05) 1 else 0
    Mixto_pen_inter_treat3[(i/k)-(sample_min),j]    <-if(pt(as.matrix(pen_intercepto@beta)[6,]/sqrt(pen_intercepto@vcov_beta[6,6]), df = as.matrix(pen_intercepto@Gp)[2,]-1)     < 0.05) 1 else 0
    
  }}
  
  #Base de datos
  Gee_inter_treat2  <- as.matrix(apply(X = Gee_intercanbiable_treat2, MARGIN = 1, FUN = mean))
  Gee_AR_treat2     <- as.matrix(apply(X = Gee_AR1_treat2,            MARGIN = 1, FUN = mean))
  Gee_unst_treat2   <- as.matrix(apply(X = Gee_unstructured_treat2,   MARGIN = 1, FUN = mean))
  Mixto_inte_treat2 <- as.matrix(apply(X = Mixto_intercepto_treat2,   MARGIN = 1, FUN = mean))
  Mixto_pen__treat2 <- as.matrix(apply(X = Mixto_pen_inter_treat2,    MARGIN = 1, FUN = mean))
  
  Gee_inter_treat3  <- as.matrix(apply(X = Gee_intercanbiable_treat3, MARGIN = 1, FUN = mean))
  Gee_AR_treat3     <- as.matrix(apply(X = Gee_AR1_treat3,            MARGIN = 1, FUN = mean))
  Gee_unst_treat3   <- as.matrix(apply(X = Gee_unstructured_treat3,   MARGIN = 1, FUN = mean))
  Mixto_inte_treat3 <- as.matrix(apply(X = Mixto_intercepto_treat3,   MARGIN = 1, FUN = mean))
  Mixto_pen__treat3 <- as.matrix(apply(X = Mixto_pen_inter_treat3,    MARGIN = 1, FUN = mean))
  
  Base <- as.data.frame(cbind(ID = (sample_min:(sample_max-1)*k)))
  Base <- mutate(Base,Gee_inter_treat2)
  Base <- mutate(Base,Gee_AR_treat2)
  Base <- mutate(Base,Gee_unst_treat2)
  Base <- mutate(Base,Mixto_inte_treat2)
  Base <- mutate(Base,Mixto_pen__treat2)
  
  Base <- mutate(Base,Gee_inter_treat3)
  Base <- mutate(Base,Gee_AR_treat3)
  Base <- mutate(Base,Gee_unst_treat3)
  Base <- mutate(Base,Mixto_inte_treat3)
  Base <- mutate(Base,Mixto_pen__treat3)
  
  #Gráfica
  #Base_largo <- reshape(data = Base, varying = 2:11, v.names = "Poder", timevar= "Modelo", idvar = "ID", direction = "long")
  #colnames(Base_largo) <- c("n","modelo","Acepta_HO")
  #Base_largo<-arrange(Base_largo,n,modelo)
  #Base_largo$modelo <- factor(Base_largo$modelo, labels = c("Gee exchangeable: treat 2", "Gee AR(1): treat 2", "Gee unstructured: treat 2",
   #                                                         "Mixto intercepto aleatorio: treat 2", "Mixto intercepto y pendiente aleatoria: treat 2",
    #                                                        "Gee exchangeable: treat 3", "Gee AR(1): treat 3", "Gee unstructured: treat 3",
     #                                                       "Mixto intercepto aleatorio: treat 3", "Mixto intercepto y pendiente aleatoria: treat 3"))
  
  #Grafico <- ggplot(data = Base_largo, aes(x = n, y = Acepta_HO, color = modelo)) +
   # geom_point(alpha = 0.3, size = 1) +
    #geom_smooth(method = loess, se = FALSE) +
    #theme_classic() +
    #labs(title="P", y="Poder", x="Tamaño de muestra", caption="Fuente: Simulación", size = 2 )
  
  return(list(Base = Base))# Grafico = Grafico, Base_largo = Base_largo))
  
}
var_cov_v0iv1i <- function(sd_v0i, sd_v1i, cor_v0iv1i) {
  var_v0i <- sd_v0i * sd_v0i  
  var_v1i <- sd_v1i * sd_v1i
  var_cov_v0iv1i <- cor_v0iv1i * ( sd_v0i * sd_v1i)
  return(list(var_v0i = var_v0i,  var_v1i = var_v1i, var_cov_v0iv1i = var_cov_v0iv1i))
}

Comp_2_treat_missing <- function(yij_2_treat, sample_min, sample_max, repeticiones, t, k, m) {
  
  #Semilla
  set.seed(16)
  
  #Librerías
  library(dplyr)
  library(gee)
  library(ggplot2)
  library(lmerTest)
  library(lme4)
  
  library(survey)
  library(sampling)

  #Matrices
  Gee_intercanbiable <- matrix(0,(sample_max-sample_min),repeticiones)
  #Gee_AR1            <- matrix(0,(sample_max-sample_min),repeticiones)
  Gee_unstructured   <- matrix(0,(sample_max-sample_min),repeticiones)
  Mixto_intercepto   <- matrix(0,(sample_max-sample_min),repeticiones)
  Mixto_pen_inter    <- matrix(0,(sample_max-sample_min),repeticiones)
  
  missing_treat_1    <- matrix(0,(sample_max-sample_min),repeticiones)
  missing_treat_2    <- matrix(0,(sample_max-sample_min),repeticiones)
  
  #Bucle
  for (i in (sample_min:sample_max)*k) {for(j in 1:repeticiones){
    
    sample_treat_1 <- yij_2_treat[sample(x = 1:(nrow(yij_2_treat)/2), size = i, replace = FALSE),]
    
    # Eliminando las observaciones paara simular las pérdidas de seguimiento para el tratamiento 1 
    
    t1 <- as.data.frame(cbind(ID = sample_treat_1$ID, t1 = sample_treat_1$V1) )
    
    t2 <- cbind(ID = sample_treat_1$ID, t2 = sample_treat_1$V2, umbral = cut(x = sample_treat_1$V2, breaks = c(-Inf, quantile(x = sample_treat_1$V2, probs = 0.10), quantile(x = sample_treat_1$V2, probs = 0.90),Inf))) 
    t2_borrar <- sampling::strata(data = t2, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_1)*m/(t-1))*0.4), ceiling((nrow(sample_treat_1)*m/(t-1))*0.2), ceiling((nrow(sample_treat_1)*m/(t-1))*0.4)))[,2]
    t2 <- as.data.frame(t2[-t2_borrar,] )
    
    t3 <- cbind(ID = sample_treat_1$ID, t3 = sample_treat_1$V3, umbral = cut(x = sample_treat_1$V3, breaks = c(-Inf, quantile(x = sample_treat_1$V3, probs = 0.10), quantile(x = sample_treat_1$V3, probs = 0.90),Inf))) 
    t3_borrar <- sampling::strata(data = t3, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_1)*m/(t-1))*0.4), ceiling((nrow(sample_treat_1)*m/(t-1))*0.2), ceiling((nrow(sample_treat_1)*m/(t-1))*0.4)))[,2]
    t3_borrar <- c(t2_borrar,t3_borrar)
    t3 <- as.data.frame(t3[-t3_borrar,])
    
    t4 <- cbind(ID = sample_treat_1$ID, t4 = sample_treat_1$V4, umbral = cut(x = sample_treat_1$V4, breaks = c(-Inf, quantile(x = sample_treat_1$V4, probs = 0.10), quantile(x = sample_treat_1$V4, probs = 0.90),Inf))) 
    t4_borrar <- sampling::strata(data = t4, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_1)*m/(t-1))*0.4), ceiling((nrow(sample_treat_1)*m/(t-1))*0.2), ceiling((nrow(sample_treat_1)*m/(t-1))*0.4)))[,2]
    t4_borrar <- c(t3_borrar,t4_borrar)
    t4 <- as.data.frame(t4[-t4_borrar,])  
    
    t5 <- cbind(ID = sample_treat_1$ID, t5 = sample_treat_1$V5, umbral = cut(x = sample_treat_1$V5, breaks = c(-Inf, quantile(x = sample_treat_1$V5, probs = 0.10), quantile(x = sample_treat_1$V5, probs = 0.90),Inf))) 
    t5_borrar <- sampling::strata(data = t5, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_1)*m/(t-1))*0.4), ceiling((nrow(sample_treat_1)*m/(t-1))*0.2), ceiling((nrow(sample_treat_1)*m/(t-1))*0.4)))[,2]
    t5_borrar <- c(t4_borrar,t5_borrar)
    t5 <- as.data.frame(t5[-t5_borrar,])  
    
    t6 <- cbind(ID = sample_treat_1$ID, t6 = sample_treat_1$V6, umbral = cut(x = sample_treat_1$V6, breaks = c(-Inf, quantile(x = sample_treat_1$V6, probs = 0.10), quantile(x = sample_treat_1$V6, probs = 0.90),Inf))) 
    t6_borrar <- sampling::strata(data = t6, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_1)*m/(t-1))*0.4), ceiling((nrow(sample_treat_1)*m/(t-1))*0.2), ceiling((nrow(sample_treat_1)*m/(t-1))*0.4)))[,2]
    t6_borrar <- c(t5_borrar,t6_borrar)
    t6 <- as.data.frame(t6[-t6_borrar,]) 
    
    t7 <- cbind(ID = sample_treat_1$ID, t7 = sample_treat_1$V7, umbral = cut(x = sample_treat_1$V7, breaks = c(-Inf, quantile(x = sample_treat_1$V7, probs = 0.10), quantile(x = sample_treat_1$V7, probs = 0.90),Inf))) 
    t7_borrar <- sampling::strata(data = t7, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_1)*m/(t-1))*0.4), ceiling((nrow(sample_treat_1)*m/(t-1))*0.2), ceiling((nrow(sample_treat_1)*m/(t-1))*0.4)))[,2]
    t7_borrar <- c(t6_borrar,t7_borrar)
    t7 <- as.data.frame(t7[-t7_borrar,]) 
    
    t8 <- cbind(ID = sample_treat_1$ID, t8 = sample_treat_1$V8, umbral = cut(x = sample_treat_1$V8, breaks = c(-Inf, quantile(x = sample_treat_1$V8, probs = 0.10), quantile(x = sample_treat_1$V8, probs = 0.90),Inf))) 
    t8_borrar <- sampling::strata(data = t8, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_1)*m/(t-1))*0.4), ceiling((nrow(sample_treat_1)*m/(t-1))*0.2), ceiling((nrow(sample_treat_1)*m/(t-1))*0.4)))[,2]
    t8_borrar <- c(t7_borrar,t8_borrar)
    t8 <- as.data.frame(t8[-t8_borrar,]) 
    
    sample_treat_1_perdidas <- full_join(x = t1,    y = t2, by = "ID")
    sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t3, by = "ID")
    sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t4, by = "ID")
    sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t5, by = "ID")
    sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t6, by = "ID")
    sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t7, by = "ID")
    sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t8, by = "ID")
    sample_treat_1_perdidas <- mutate(sample_treat_1_perdidas, treat = rep(1, nrow(sample_treat_1)))
    
    sample_treat_1_perdidas <- sample_treat_1_perdidas[,-c(4,6,8,10,12,14,16)]
    
    missing_treat_1[(i/k)-(sample_min),j]    <- length(t8_borrar)/nrow(sample_treat_1_perdidas)
    
    # Eliminando las observaciones paara simular las pérdidas de seguimiento para el tratamiento 2
    
    sample_treat_2 <- yij_2_treat[sample(x = (nrow(yij_2_treat)/2+1):nrow(yij_2_treat), size = i, replace = FALSE),]
    
    t1 <- as.data.frame(cbind(ID = sample_treat_2$ID, t1 = sample_treat_2$V1) )
    
    t2 <- cbind(ID = sample_treat_2$ID, t2 = sample_treat_2$V2, umbral = cut(x = sample_treat_2$V2, breaks = c(-Inf, quantile(x = sample_treat_2$V2, probs = 0.10), quantile(x = sample_treat_2$V2, probs = 0.90),Inf))) 
    t2_borrar <- sampling::strata(data = t2, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_2)*m/(t-1))*0.4), ceiling((nrow(sample_treat_2)*m/(t-1))*0.2), ceiling((nrow(sample_treat_2)*m/(t-1))*0.4)))[,2]
    t2 <- as.data.frame(t2[-t2_borrar,] )
    
    t3 <- cbind(ID = sample_treat_2$ID, t3 = sample_treat_2$V3, umbral = cut(x = sample_treat_2$V3, breaks = c(-Inf, quantile(x = sample_treat_2$V3, probs = 0.10), quantile(x = sample_treat_2$V3, probs = 0.90),Inf))) 
    t3_borrar <- sampling::strata(data = t3, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_2)*m/(t-1))*0.4), ceiling((nrow(sample_treat_2)*m/(t-1))*0.2), ceiling((nrow(sample_treat_2)*m/(t-1))*0.4)))[,2]
    t3_borrar <- c(t2_borrar,t3_borrar)
    t3 <- as.data.frame(t3[-t3_borrar,])
    
    t4 <- cbind(ID = sample_treat_2$ID, t4 = sample_treat_2$V4, umbral = cut(x = sample_treat_2$V4, breaks = c(-Inf, quantile(x = sample_treat_2$V4, probs = 0.10), quantile(x = sample_treat_2$V4, probs = 0.90),Inf))) 
    t4_borrar <- sampling::strata(data = t4, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_2)*m/(t-1))*0.4), ceiling((nrow(sample_treat_2)*m/(t-1))*0.2), ceiling((nrow(sample_treat_2)*m/(t-1))*0.4)))[,2]
    t4_borrar <- c(t3_borrar,t4_borrar)
    t4 <- as.data.frame(t4[-t4_borrar,])  
    
    t5 <- cbind(ID = sample_treat_2$ID, t5 = sample_treat_2$V5, umbral = cut(x = sample_treat_2$V5, breaks = c(-Inf, quantile(x = sample_treat_2$V5, probs = 0.10), quantile(x = sample_treat_2$V5, probs = 0.90),Inf))) 
    t5_borrar <- sampling::strata(data = t5, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_2)*m/(t-1))*0.4), ceiling((nrow(sample_treat_2)*m/(t-1))*0.2), ceiling((nrow(sample_treat_2)*m/(t-1))*0.4)))[,2]
    t5_borrar <- c(t4_borrar,t5_borrar)
    t5 <- as.data.frame(t5[-t5_borrar,])  
    
    t6 <- cbind(ID = sample_treat_2$ID, t6 = sample_treat_2$V6, umbral = cut(x = sample_treat_2$V6, breaks = c(-Inf, quantile(x = sample_treat_2$V6, probs = 0.10), quantile(x = sample_treat_2$V6, probs = 0.90),Inf))) 
    t6_borrar <- sampling::strata(data = t6, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_2)*m/(t-1))*0.4), ceiling((nrow(sample_treat_2)*m/(t-1))*0.2), ceiling((nrow(sample_treat_2)*m/(t-1))*0.4)))[,2]
    t6_borrar <- c(t5_borrar,t6_borrar)
    t6 <- as.data.frame(t6[-t6_borrar,]) 
    
    t7 <- cbind(ID = sample_treat_2$ID, t7 = sample_treat_2$V7, umbral = cut(x = sample_treat_2$V7, breaks = c(-Inf, quantile(x = sample_treat_2$V7, probs = 0.10), quantile(x = sample_treat_2$V7, probs = 0.90),Inf))) 
    t7_borrar <- sampling::strata(data = t7, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_2)*m/(t-1))*0.4), ceiling((nrow(sample_treat_2)*m/(t-1))*0.2), ceiling((nrow(sample_treat_2)*m/(t-1))*0.4)))[,2]
    t7_borrar <- c(t6_borrar,t7_borrar)
    t7 <- as.data.frame(t7[-t7_borrar,]) 
    
    t8 <- cbind(ID = sample_treat_2$ID, t8 = sample_treat_2$V8, umbral = cut(x = sample_treat_2$V8, breaks = c(-Inf, quantile(x = sample_treat_2$V8, probs = 0.10), quantile(x = sample_treat_2$V8, probs = 0.90),Inf))) 
    t8_borrar <- sampling::strata(data = t8, method = c("srswor"), stratanames = c("umbral"), size = c( ceiling((nrow(sample_treat_2)*m/(t-1))*0.4), ceiling((nrow(sample_treat_2)*m/(t-1))*0.2), ceiling((nrow(sample_treat_2)*m/(t-1))*0.4)))[,2]
    t8_borrar <- c(t7_borrar,t8_borrar)
    t8 <- as.data.frame(t8[-t8_borrar,]) 
    
    sample_treat_2_perdidas <- full_join(x = t1,    y = t2, by = "ID")
    sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t3, by = "ID")
    sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t4, by = "ID")
    sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t5, by = "ID")
    sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t6, by = "ID")
    sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t7, by = "ID")
    sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t8, by = "ID")
    sample_treat_2_perdidas <- mutate(sample_treat_2_perdidas, treat = rep(2, nrow(sample_treat_2_perdidas)))
    
    sample_treat_2_perdidas <- sample_treat_2_perdidas[,-c(4,6,8,10,12,14,16)]  
    
    missing_treat_2[(i/k)-(sample_min),j]    <- length(t8_borrar)/nrow(sample_treat_2_perdidas)
    
    # Uniendo las bases de los tratamientos con pérdidas de seguimiento
    
    sample <- bind_rows(sample_treat_1_perdidas,sample_treat_2_perdidas)
    
    # Pasando la base a formato ancho
    
    sample_long <- reshape(data = sample,varying = 2:(t+1), v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
    sample_long <- arrange(sample_long,ID,tiempo)
    sample_long$tiempo <- as.numeric(sample_long$tiempo)
    sample_long$tiempo <- (sample_long$tiempo-1)/(t-1)
    
    #Modelos
    intercanbiable <- gee(yij ~ treat + tiempo + treat * tiempo, id = ID, data = sample_long, family = gaussian, corstr = "exchangeable")
    #AR1            <- gee(yij ~ treat + tiempo + treat * tiempo, id = ID, data = sample_long, family = gaussian, corstr = "AR-M", Mv = 1)
    unstructured   <- gee(yij ~ treat + tiempo + treat * tiempo, id = ID, data = sample_long, family = gaussian, corstr = "unstructured")
    intercepto     <- lmer(yij ~ treat + tiempo + treat * tiempo + (1|ID), data = sample_long, REML = FALSE)
    pen_intercepto <- lmer(yij ~ treat + tiempo + treat * tiempo + (tiempo|ID), data = sample_long, REML = FALSE)
    
    #Completando las matrices con la decisión de la hipótesis
    Gee_intercanbiable[(i/k)-(sample_min),j] <-if(pnorm(as.matrix(intercanbiable$coefficients)[4,]/sqrt(intercanbiable$robust.variance[4,4]), 0, 1)        < 0.05) 1 else 0
    #Gee_AR1[(i/k)-(sample_min),j]            <-if(pnorm(as.matrix(AR1$coefficients)[4,]/sqrt(AR1$robust.variance[4,4]), 0, 1)                              < 0.05) 1 else 0
    Gee_unstructured[(i/k)-(sample_min),j]   <-if(pnorm(as.matrix(unstructured$coefficients)[4,]/sqrt(unstructured$robust.variance[4,4]), 0, 1)            < 0.05) 1 else 0
    Mixto_intercepto[(i/k)-(sample_min),j]   <-if(pt(as.matrix(intercepto@beta)[4,]/sqrt(intercepto@vcov_beta[4,4]), df = as.matrix(intercepto@Gp)[2,]-1)     < 0.05) 1 else 0
    Mixto_pen_inter[(i/k)-(sample_min),j]    <-if(pt(as.matrix(pen_intercepto@beta)[4,]/sqrt(pen_intercepto@vcov_beta[4,4]), df = as.matrix(pen_intercepto@Gp)[2,]-1)     < 0.05) 1 else 0
    
    
    
    
  }}
  
  #Base de datos
  Gee_inter  <- as.matrix(apply(X = Gee_intercanbiable, MARGIN = 1, FUN = mean))
  #Gee_AR     <- as.matrix(apply(X = Gee_AR1,            MARGIN = 1, FUN = mean))
  Gee_unst   <- as.matrix(apply(X = Gee_unstructured,   MARGIN = 1, FUN = mean))
  Mixto_inte <- as.matrix(apply(X = Mixto_intercepto,   MARGIN = 1, FUN = mean))
  Mixto_pen_ <- as.matrix(apply(X = Mixto_pen_inter,    MARGIN = 1, FUN = mean))
  
  missing_treat__1 <- as.matrix(apply(X = missing_treat_1,    MARGIN = 1, FUN = mean))
  missing_treat__2 <- as.matrix(apply(X = missing_treat_2,    MARGIN = 1, FUN = mean))
  
  Base <- as.data.frame(cbind(ID = (sample_min:(sample_max-1)*k)))
  Base <- mutate(Base,Gee_inter)
  #Base <- mutate(Base,Gee_AR)
  Base <- mutate(Base,Gee_unst)
  Base <- mutate(Base,Mixto_inte)
  Base <- mutate(Base,missing_treat__1)
  Base <- mutate(Base,missing_treat__2)
  
  
  return(Base = Base)
  
}
                                        
                                        

# Escenario 1
var_cov_v0iv1i(sd_v0i = 10.1, sd_v1i = 8, cor_v0iv1i = -0.75)

# Población
población <- Pob_ECA(n = 1000, 
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

yij_2_treat <- población$Treat_2$yij_2
yij_3_treat <- población$Treat_3$yij_3

#Escenarios --------------------------------------------------------------------

# Escenario
Escenario <- Comp_2_treat(yij_2_treat = yij_2_treat, sample_min = 10, sample_max = 20, repeticiones = 10, t = 4, k = 1)

# Escenario
Escenario <- Comp_3_treat(yij_2_treat = yij_2_treat, sample_min = 10, sample_max = 20, repeticiones = 10, t = 4, k = 1)
                                        
                                        
