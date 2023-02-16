{
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
  Comp_modelos     <- function(base, treat, n, repeticiones, t) {
    
    # Capturando la hora de inicio del la función
    Inicio <- DescTools::Now()
    
    # Semilla para fijar los resultados
    set.seed(16)
    
    # Librerías que utiliza la función
    library(dplyr)
    library(lmerTest)
    library(lme4)
    library(geepack)
    
    # Matrices que guardarán los resultados de las hipótesis planteadas para 2 y 3 tratamientos
    Gee_intercanbiable    <- matrix(0,length(n),repeticiones)
    Gee_AR1               <- matrix(0,length(n),repeticiones)
    Gee_unstructured      <- matrix(0,length(n),repeticiones) 
    Mixto_intercepto      <- matrix(0,length(n),repeticiones)
    Mixto_pen_inter       <- matrix(0,length(n),repeticiones)
    
    Gee_intercanbiable_3  <- matrix(0,length(n),repeticiones)
    Gee_AR1_3             <- matrix(0,length(n),repeticiones)
    Gee_unstructured_3    <- matrix(0,length(n),repeticiones) 
    Mixto_intercepto_3    <- matrix(0,length(n),repeticiones)
    Mixto_pen_inter_3     <- matrix(0,length(n),repeticiones)
    
    # Cálculos para 2 tratamientos
    if(treat == 2){ 
      
      # Bucle para 2 tratamientos
      for (i in 1:length(n)) {for(j in 1:repeticiones){ 
        
        sample_treat_1 <- base[sample(x = 1:(nrow(base)/2),            size = n[[i]], replace = FALSE),]
        sample_treat_2 <- base[sample(x = (nrow(base)/2+1):nrow(base), size = n[[i]], replace = FALSE),]
        sample         <- bind_rows(sample_treat_1,sample_treat_2)
        
        sample_long <- reshape(data = sample, varying = 1:t, v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
        sample_long <- arrange(sample_long,ID,tiempo)
        sample_long$treat <- as.factor(sample_long$treat)
        sample_long$tiempo <- as.numeric(sample_long$tiempo)
        sample_long$tiempo <- (sample_long$tiempo-1)/(t-1) # acá se estandariza el tiempo (de cero a uno) y garantiza que entre una medición y otra "t" tenga la misma distancia.
        
        #Modelos
        intercanbiable <- geeglm(formula = yij ~ treat + tiempo + treat * tiempo,      id = ID, data =  sample_long, family = gaussian, corstr = "exchangeable")
        AR1            <- geeglm(formula = yij ~ treat + tiempo + treat * tiempo,      id = ID, data =  sample_long, family = gaussian, corstr = "ar1")
        unstructured   <- geeglm(formula = yij ~ treat + tiempo + treat * tiempo,      id = ID, data =  sample_long, family = gaussian, corstr = "unstructured")
        intercepto     <- lmer  (formula = yij ~ treat + tiempo + treat * tiempo +  (1     |ID), data = sample_long, REML = FALSE)
        pen_intercepto <- lmer  (formula = yij ~ treat + tiempo + treat * tiempo +  (tiempo|ID), data = sample_long, REML = FALSE)
        
        #Completando las matrices con la decisión de la hipótesis
        
        Gee_intercanbiable [i,j]  <- if((1-(pnorm( abs(coef(summary(intercanbiable))[4,1] / coef(summary(intercanbiable))[4,2] ))))*2 < 0.05) 1 else 0
        Gee_AR1            [i,j]  <- if((1-(pnorm( abs(coef(summary(AR1))[4,1]            / coef(summary(AR1))[4,2]            ))))*2 < 0.05) 1 else 0
        Gee_unstructured   [i,j]  <- if((1-(pnorm( abs(coef(summary(unstructured))[4,1]   / coef(summary(unstructured))[4,2]   ))))*2 < 0.05) 1 else 0
        Mixto_intercepto   [i,j]  <- if(coef(summary(intercepto))    [4,5]                                                            < 0.05) 1 else 0
        Mixto_pen_inter    [i,j]  <- if(coef(summary(pen_intercepto))[4,5]                                                            < 0.05) 1 else 0
        
        print(c(n[[i]],j))
        
      }}
      
      #Base de datos para 2 tratamientos
      Gee_inter  <- as.matrix(apply(X = Gee_intercanbiable, MARGIN = 1, FUN = mean))
      Gee_AR     <- as.matrix(apply(X = Gee_AR1,            MARGIN = 1, FUN = mean))
      Gee_unst   <- as.matrix(apply(X = Gee_unstructured,   MARGIN = 1, FUN = mean))
      Mixto_inte <- as.matrix(apply(X = Mixto_intercepto,   MARGIN = 1, FUN = mean))
      Mixto_pen_ <- as.matrix(apply(X = Mixto_pen_inter,    MARGIN = 1, FUN = mean))
      
      Base <- as.data.frame(cbind(n = n))
      
      Base <- mutate(Base,Gee_inter)
      Base <- mutate(Base,Gee_AR)
      Base <- mutate(Base,Gee_unst)
      Base <- mutate(Base,Mixto_inte)
      Base <- mutate(Base,Mixto_pen_)
    }
    
    # Cálculos para 2 tratamientos
    if(treat == 3){    
      
      # Bucle para 3 tratamientos
      for (i in 1:length(n)) {for(j in 1:repeticiones){
        
        sample_treat_1 <- base[sample(x = 1:(nrow(base)/3),                    size = n[[i]], replace = FALSE),]
        sample_treat_2 <- base[sample(x = (nrow(base)/3+1):((nrow(base)/3)*2), size = n[[i]], replace = FALSE),]
        sample_treat_3 <- base[sample(x = ((nrow(base)/3)*2+1):nrow(base),     size = n[[i]], replace = FALSE),]
        
        sample <- bind_rows(sample_treat_1,sample_treat_2,sample_treat_3)
        
        sample_long <- reshape(data = sample, varying = 1:t, v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
        sample_long <- arrange(sample_long,ID,tiempo)
        sample_long$treat <- as.factor(sample_long$treat)
        sample_long$tiempo <- as.numeric(sample_long$tiempo)
        sample_long$tiempo <- (sample_long$tiempo-1)/(t-1) # acá se estandariza el tiempo (de cero a uno) y garantiza que entre una medición y otra "t" tenga la misma distancia.
        
        #Modelos
        intercanbiable <- geeglm(formula = yij ~ treat + tiempo + treat * tiempo,     id = ID,    data = sample_long, family = gaussian, corstr = "exchangeable")
        AR1            <- geeglm(formula = yij ~ treat + tiempo + treat * tiempo,     id = ID,    data = sample_long, family = gaussian, corstr = "ar1")
        unstructured   <- geeglm(formula = yij ~ treat + tiempo + treat * tiempo,     id = ID,    data = sample_long, family = gaussian, corstr = "unstructured")
        intercepto     <- lmer  (formula = yij ~ treat + tiempo + treat * tiempo + (1     |ID),   data = sample_long, REML = FALSE)
        pen_intercepto <- lmer  (formula = yij ~ treat + tiempo + treat * tiempo + (tiempo|ID),   data = sample_long, REML = FALSE)
        
        #Completando las matrices con la decisión de la hipótesis
        
        Gee_intercanbiable   [i,j]  <- if((1-(pnorm( abs( coef(summary(intercanbiable)) [5,1]  /   coef(summary(intercanbiable))[5,2] ))))*2 < 0.05) 1 else 0
        Gee_AR1              [i,j]  <- if((1-(pnorm( abs( coef(summary(AR1))            [5,1]  /   coef(summary(AR1))           [5,2] ))))*2 < 0.05) 1 else 0
        Gee_unstructured     [i,j]  <- if((1-(pnorm( abs( coef(summary(unstructured))   [5,1]  /   coef(summary(unstructured))  [5,2] ))))*2 < 0.05) 1 else 0
        Mixto_intercepto     [i,j]  <- if(coef(summary(intercepto))    [5,5]                                                                 < 0.05) 1 else 0
        Mixto_pen_inter      [i,j]  <- if(coef(summary(pen_intercepto))[5,5]                                                                 < 0.05) 1 else 0
        
        Gee_intercanbiable_3 [i,j]  <- if((1-(pnorm( abs( coef(summary(intercanbiable)) [6,1]  /   coef(summary(intercanbiable))[6,2] ))))*2 < 0.05) 1 else 0
        Gee_AR1_3            [i,j]  <- if((1-(pnorm( abs( coef(summary(AR1))            [6,1]  /   coef(summary(AR1))           [6,2] ))))*2 < 0.05) 1 else 0
        Gee_unstructured_3   [i,j]  <- if((1-(pnorm( abs( coef(summary(unstructured))   [6,1]  /   coef(summary(unstructured))  [6,2] ))))*2 < 0.05) 1 else 0
        Mixto_intercepto_3   [i,j]  <- if(coef(summary(intercepto))    [6,5]                                                                 < 0.05) 1 else 0
        Mixto_pen_inter_3    [i,j]  <- if(coef(summary(pen_intercepto))[6,5]                                                                 < 0.05) 1 else 0
        
        print(c(n[[i]],j))
        
      }}
      
      #Base de datos para 3 tratamientos
      Gee_inter     <- as.matrix(apply(X = Gee_intercanbiable,   MARGIN = 1, FUN = mean))
      Gee_AR        <- as.matrix(apply(X = Gee_AR1,              MARGIN = 1, FUN = mean))
      Gee_unst      <- as.matrix(apply(X = Gee_unstructured,     MARGIN = 1, FUN = mean))
      Mixto_inte    <- as.matrix(apply(X = Mixto_intercepto,     MARGIN = 1, FUN = mean))
      Mixto_pen_    <- as.matrix(apply(X = Mixto_pen_inter,      MARGIN = 1, FUN = mean))
      
      Gee_inter_3   <- as.matrix(apply(X = Gee_intercanbiable_3, MARGIN = 1, FUN = mean))
      Gee_AR_3      <- as.matrix(apply(X = Gee_AR1_3,            MARGIN = 1, FUN = mean))
      Gee_unst_3    <- as.matrix(apply(X = Gee_unstructured_3,   MARGIN = 1, FUN = mean))
      Mixto_inte_3  <- as.matrix(apply(X = Mixto_intercepto_3,   MARGIN = 1, FUN = mean))
      Mixto_pen_3   <- as.matrix(apply(X = Mixto_pen_inter_3,    MARGIN = 1, FUN = mean))
      
      Base <- as.data.frame(cbind(ID = n))
      
      Base <- mutate(Base,Gee_inter)
      Base <- mutate(Base,Gee_AR)
      Base <- mutate(Base,Gee_unst)
      Base <- mutate(Base,Mixto_inte)
      Base <- mutate(Base,Mixto_pen_)
      
      Base <- mutate(Base,Gee_inter_3)
      Base <- mutate(Base,Gee_AR_3)
      Base <- mutate(Base,Gee_unst_3)
      Base <- mutate(Base,Mixto_inte_3)
      Base <- mutate(Base,Mixto_pen_3)
      
    } else {"Favor indicar 2 o 3  tratamientos en 'treat'"}
    
    
    # Capturando la hora de término del la función
    Fin <- DescTools::Now()
    
    # Calculando la duración
    Duración <- Fin - Inicio; print(Duración)
    
    # Retornando los resultados
    return(Base = Base)
    
  }
  Comp_modelos_par <- function(base, treat, n, repeticiones, t) {
    
    # Capturando la hora de inicio del la función
    Inicio <- DescTools::Now()
    
    # Librerías que utiliza la función
    library(dplyr)
    library(lmerTest)
    library(lme4)
    library(geepack)
    library(parallel)
    library(foreach)
    
    # Matrices que guardarán los resultados de las hipótesis planteadas para 2 y 3 tratamientos
    Gee_intercanbiable    <- matrix(0,length(n),repeticiones)
    Gee_AR1               <- matrix(0,length(n),repeticiones)
    Gee_unstructured      <- matrix(0,length(n),repeticiones) 
    Mixto_intercepto      <- matrix(0,length(n),repeticiones)
    Mixto_pen_inter       <- matrix(0,length(n),repeticiones)
    
    Gee_intercanbiable_3  <- matrix(0,length(n),repeticiones)
    Gee_AR1_3             <- matrix(0,length(n),repeticiones)
    Gee_unstructured_3    <- matrix(0,length(n),repeticiones) 
    Mixto_intercepto_3    <- matrix(0,length(n),repeticiones)
    Mixto_pen_inter_3     <- matrix(0,length(n),repeticiones)
    
    # Ajuntando los parámetros para utilizar %dopar%
    doParallel::registerDoParallel( cl = parallel::makeCluster( parallel::detectCores() - 1, type = "PSOCK" ))
    
    # Cálculos para 2 tratamientos
    if(treat == 2){ 
      
      # Bucle para 2 tratamientos
      
      for (i in 1:length(n)) {
        i_result <- foreach (j = 1:repeticiones) %dopar%  {
          
          # Semilla para fijar los resultados
          set.seed(i * j)
          
          # Obteniendo la muestra
          sample_treat_1 <- base[sample(x = 1:(nrow(base)/2),            size = n[[i]], replace = FALSE),]
          sample_treat_2 <- base[sample(x = (nrow(base)/2+1):nrow(base), size = n[[i]], replace = FALSE),]
          sample         <- dplyr::bind_rows(sample_treat_1,sample_treat_2)
          
          # Pasando a formato largo la muestra obtenida
          sample_long <- reshape(data = sample, varying = 1:t, v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
          sample_long <- dplyr::arrange(sample_long,ID,tiempo)
          sample_long$treat <- as.factor(sample_long$treat)
          sample_long$tiempo <- as.numeric(sample_long$tiempo)
          sample_long$tiempo <- (sample_long$tiempo-1)/(t-1) # acá se estandariza el tiempo (de cero a uno) y garantiza que entre una medición y otra "t" tenga la misma distancia.
          
          #Ajustando los modelos
          intercanbiable <- geepack::geeglm(formula = yij ~ treat + tiempo + treat * tiempo,       id = ID, data = sample_long, family = gaussian, corstr = "exchangeable")
          AR1            <- geepack::geeglm(formula = yij ~ treat + tiempo + treat * tiempo,       id = ID, data = sample_long, family = gaussian, corstr = "ar1")
          unstructured   <- geepack::geeglm(formula = yij ~ treat + tiempo + treat * tiempo,       id = ID, data = sample_long, family = gaussian, corstr = "unstructured")
          intercepto     <- lmerTest::lmer     (formula = yij ~ treat + tiempo + treat * tiempo +  (1     |ID), data = sample_long, REML = FALSE)
          pen_intercepto <- lmerTest::lmer     (formula = yij ~ treat + tiempo + treat * tiempo +  (tiempo|ID), data = sample_long, REML = FALSE)
          
          #Guardando la decisión de la hipótesis planteada  
          gee_intercanbiable <- if((1-(pnorm( abs(coef(summary(intercanbiable))[4,1] / coef(summary(intercanbiable))[4,2] ))))*2 < 0.05) 1 else 0
          gee_AR1            <- if((1-(pnorm( abs(coef(summary(AR1))[4,1]            / coef(summary(AR1))[4,2]            ))))*2 < 0.05) 1 else 0
          gee_unstructured   <- if((1-(pnorm( abs(coef(summary(unstructured))[4,1]   / coef(summary(unstructured))[4,2]   ))))*2 < 0.05) 1 else 0
          mixto_intercepto   <- if(coef(summary(intercepto))    [4,5]                                                            < 0.05) 1 else 0
          mixto_pen_inter    <- if(coef(summary(pen_intercepto))[4,5]                                                            < 0.05) 1 else 0
          
          #Resultado que se guardará en i_result para cada i
          r <- list(
            j = j,
            gee_intercanbiable = as.integer(gee_intercanbiable),
            gee_AR1 = as.integer(gee_AR1),
            gee_unstructured = as.integer(gee_unstructured),
            mixto_intercepto = as.integer(mixto_intercepto),
            mixto_pen_inter = as.integer(mixto_pen_inter)
          )
          r
        }
        
        #Completando las matrices con la decisión de la hipótesis 
        for(r in i_result) {
          j <- r$j
          Gee_intercanbiable [i,j] <- r$gee_intercanbiable
          Gee_AR1            [i,j] <- r$gee_AR1
          Gee_unstructured   [i,j] <- r$gee_unstructured
          Mixto_intercepto   [i,j] <- r$mixto_intercepto
          Mixto_pen_inter    [i,j] <- r$mixto_pen_inter
        }
        print(n[[i]])
        print(DescTools::Now() - Inicio)
      }
      
      #Base de datos para 2 tratamientos
      Gee_inter  <- as.matrix(apply(X = Gee_intercanbiable, MARGIN = 1, FUN = mean))
      Gee_AR     <- as.matrix(apply(X = Gee_AR1,            MARGIN = 1, FUN = mean))
      Gee_unst   <- as.matrix(apply(X = Gee_unstructured,   MARGIN = 1, FUN = mean))
      Mixto_inte <- as.matrix(apply(X = Mixto_intercepto,   MARGIN = 1, FUN = mean))
      Mixto_pen_ <- as.matrix(apply(X = Mixto_pen_inter,    MARGIN = 1, FUN = mean))
      
      Base <- as.data.frame(cbind(n = n))
      
      Base <- mutate(Base,Gee_inter)
      Base <- mutate(Base,Gee_AR)
      Base <- mutate(Base,Gee_unst)
      Base <- mutate(Base,Mixto_inte)
      Base <- mutate(Base,Mixto_pen_)
    }
    
    # Cálculos para 3 tratamientos
    if(treat == 3){ 
      
      # Bucle para 3 tratamientos
      
      for (i in 1:length(n)) {
        i_result <- foreach (j = 1:repeticiones) %dopar%  {
          
          # Semilla para fijar los resultados
          set.seed(i * j)
          
          # Obteniendo la muestra
          sample_treat_1 <- base[sample(x = 1:(nrow(base)/3),                    size = n[[i]], replace = FALSE),]
          sample_treat_2 <- base[sample(x = (nrow(base)/3+1):((nrow(base)/3)*2), size = n[[i]], replace = FALSE),]
          sample_treat_3 <- base[sample(x = ((nrow(base)/3)*2+1):nrow(base),     size = n[[i]], replace = FALSE),]
          sample         <- dplyr::bind_rows(sample_treat_1,sample_treat_2,sample_treat_3)
          
          
          # Pasando a formato largo la muestra obtenida
          sample_long <- reshape(data = sample, varying = 1:t, v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
          sample_long <- dplyr::arrange(sample_long,ID,tiempo)
          sample_long$treat <- as.factor(sample_long$treat)
          sample_long$tiempo <- as.numeric(sample_long$tiempo)
          sample_long$tiempo <- (sample_long$tiempo-1)/(t-1) # acá se estandariza el tiempo (de cero a uno) y garantiza que entre una medición y otra "t" tenga la misma distancia.
          
          #Ajustando los modelos
          intercanbiable <- geepack::geeglm(formula = yij ~ treat + tiempo + treat * tiempo,       id = ID, data = sample_long, family = gaussian, corstr = "exchangeable")
          AR1            <- geepack::geeglm(formula = yij ~ treat + tiempo + treat * tiempo,       id = ID, data = sample_long, family = gaussian, corstr = "ar1")
          unstructured   <- geepack::geeglm(formula = yij ~ treat + tiempo + treat * tiempo,       id = ID, data = sample_long, family = gaussian, corstr = "unstructured")
          intercepto     <- lmerTest::lmer     (formula = yij ~ treat + tiempo + treat * tiempo +  (1     |ID), data = sample_long, REML = FALSE)
          pen_intercepto <- lmerTest::lmer     (formula = yij ~ treat + tiempo + treat * tiempo +  (tiempo|ID), data = sample_long, REML = FALSE)
          
          #Guardando la decisión de la hipótesis planteada  
          gee_intercanbiable <- if((1-(pnorm( abs(coef(summary(intercanbiable))[5,1] / coef(summary(intercanbiable))[5,2] ))))*2 < 0.05) 1 else 0
          gee_AR1            <- if((1-(pnorm( abs(coef(summary(AR1))           [5,1] / coef(summary(AR1))           [5,2] ))))*2 < 0.05) 1 else 0
          gee_unstructured   <- if((1-(pnorm( abs(coef(summary(unstructured))  [5,1] / coef(summary(unstructured))  [5,2] ))))*2 < 0.05) 1 else 0
          mixto_intercepto   <- if(coef(summary(intercepto))    [5,5]                                                            < 0.05) 1 else 0
          mixto_pen_inter    <- if(coef(summary(pen_intercepto))[5,5]                                                            < 0.05) 1 else 0
          
          gee_intercanbiable_3 <- if((1-(pnorm( abs(coef(summary(intercanbiable))[6,1] / coef(summary(intercanbiable))[6,2] ))))*2 < 0.05) 1 else 0
          gee_AR1_3            <- if((1-(pnorm( abs(coef(summary(AR1))           [6,1] / coef(summary(AR1))           [6,2] ))))*2 < 0.05) 1 else 0
          gee_unstructured_3   <- if((1-(pnorm( abs(coef(summary(unstructured))  [6,1]   / coef(summary(unstructured))[6,2] ))))*2 < 0.05) 1 else 0
          mixto_intercepto_3   <- if(coef(summary(intercepto))    [6,5]                                                            < 0.05) 1 else 0
          mixto_pen_inter_3    <- if(coef(summary(pen_intercepto))[6,5]                                                            < 0.05) 1 else 0
          
          #Resultado que se guardará en i_result para cada i
          r <- list(
            j = j,
            gee_intercanbiable =  as.integer(gee_intercanbiable),
            gee_AR1            =  as.integer(gee_AR1),
            gee_unstructured   =  as.integer(gee_unstructured),
            mixto_intercepto   =  as.integer(mixto_intercepto),
            mixto_pen_inter    =  as.integer(mixto_pen_inter),
            gee_intercanbiable_3  =  as.integer(gee_intercanbiable_3),
            gee_AR1_3             =  as.integer(gee_AR1_3),
            gee_unstructured_3    =  as.integer(gee_unstructured_3),
            mixto_intercepto_3    =  as.integer(mixto_intercepto_3),
            mixto_pen_inter_3     =  as.integer(mixto_pen_inter_3)
          )
          r
        }
        
        #Completando las matrices con la decisión de la hipótesis 
        for(r in i_result) {
          j <- r$j
          Gee_intercanbiable[i,j] <- r$gee_intercanbiable
          Gee_AR1[i,j] <- r$gee_AR1
          Gee_unstructured[i,j] <- r$gee_unstructured
          Mixto_intercepto[i,j] <- r$mixto_intercepto
          Mixto_pen_inter[i,j] <- r$mixto_pen_inter
          Gee_intercanbiable_3[i,j] <- r$gee_intercanbiable_3
          Gee_AR1_3[i,j] <- r$gee_AR1_3
          Gee_unstructured_3[i,j] <- r$gee_unstructured_3 
          Mixto_intercepto_3[i,j] <- r$mixto_intercepto_3
          Mixto_pen_inter_3[i,j] <- r$mixto_pen_inter_3
          
        }
        print(n[[i]])
        print(DescTools::Now() - Inicio)
      }
      
      #Base de datos para 3 tratamientos
      Gee_inter     <- as.matrix(apply(X = Gee_intercanbiable,   MARGIN = 1, FUN = mean))
      Gee_AR        <- as.matrix(apply(X = Gee_AR1,              MARGIN = 1, FUN = mean))
      Gee_unst      <- as.matrix(apply(X = Gee_unstructured,     MARGIN = 1, FUN = mean))
      Mixto_inte    <- as.matrix(apply(X = Mixto_intercepto,     MARGIN = 1, FUN = mean))
      Mixto_pen_    <- as.matrix(apply(X = Mixto_pen_inter,      MARGIN = 1, FUN = mean))
      
      Gee_inter_3   <- as.matrix(apply(X = Gee_intercanbiable_3, MARGIN = 1, FUN = mean))
      Gee_AR_3      <- as.matrix(apply(X = Gee_AR1_3,            MARGIN = 1, FUN = mean))
      Gee_unst_3    <- as.matrix(apply(X = Gee_unstructured_3,   MARGIN = 1, FUN = mean))
      Mixto_inte_3  <- as.matrix(apply(X = Mixto_intercepto_3,   MARGIN = 1, FUN = mean))
      Mixto_pen_3   <- as.matrix(apply(X = Mixto_pen_inter_3,    MARGIN = 1, FUN = mean))
      
      Base <- as.data.frame(cbind(ID = n))
      
      Base <- mutate(Base,Gee_inter)
      Base <- mutate(Base,Gee_AR)
      Base <- mutate(Base,Gee_unst)
      Base <- mutate(Base,Mixto_inte)
      Base <- mutate(Base,Mixto_pen_)
      
      Base <- mutate(Base,Gee_inter_3)
      Base <- mutate(Base,Gee_AR_3)
      Base <- mutate(Base,Gee_unst_3)
      Base <- mutate(Base,Mixto_inte_3)
      Base <- mutate(Base,Mixto_pen_3)
    }
    
    # Capturando la hora de término del la función
    Fin <- DescTools::Now()
    
    # Calculando la duración
    Duración <- Fin - Inicio; print(Duración)
    
    # Retornando los resultados
    return(Base = Base)
    
  }
  Comp_modelos_missing <- function(base, treat, n, repeticiones, t, q, mi, mc, ms) {
    
    # Capturando la hora de inicio del la función
    Inicio <- DescTools::Now()
    
    # Semilla para fijar los resultados
    set.seed(16)
    
    # Librerías que utiliza la función
    library(dplyr)
    library(lmerTest)
    library(lme4)
    library(geepack)
    
    # Matrices que guardarán los resultados de las hipótesis planteadas para 2 y 3 tratamientos
    Gee_intercanbiable    <- matrix(0,length(n),repeticiones)
    Gee_AR1               <- matrix(0,length(n),repeticiones)
    Gee_unstructured      <- matrix(0,length(n),repeticiones) 
    Mixto_intercepto      <- matrix(0,length(n),repeticiones)
    Mixto_pen_inter       <- matrix(0,length(n),repeticiones)
    
    Gee_intercanbiable_3  <- matrix(0,length(n),repeticiones)
    Gee_AR1_3             <- matrix(0,length(n),repeticiones)
    Gee_unstructured_3    <- matrix(0,length(n),repeticiones) 
    Mixto_intercepto_3    <- matrix(0,length(n),repeticiones)
    Mixto_pen_inter_3     <- matrix(0,length(n),repeticiones)
    
    # Bucle para 2 tratamientos y 8 mediciones
    if (treat == 2 & t == 4) {
      for (i in 1:length(n)) {for(j in 1:repeticiones){ 
        
        sample_treat_1 <- base[sample(x = 1:(nrow(base)/2),            size = n[[i]], replace = FALSE),]
        sample_treat_2 <- base[sample(x = (nrow(base)/2+1):nrow(base), size = n[[i]], replace = FALSE),]
        
        # Eliminando las observaciones para el tratamiento 1
        t1 <- as.data.frame(cbind(ID = sample_treat_1$ID, t1 = sample_treat_1$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_1$ID, t2 = sample_treat_1$V2, umbral = cut(x = sample_treat_1$V2, breaks = c(-Inf, quantile(x = sample_treat_1$V2, probs = q), quantile(x = sample_treat_1$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_1$ID, t3 = sample_treat_1$V3, umbral = cut(x = sample_treat_1$V3, breaks = c(-Inf, quantile(x = sample_treat_1$V3, probs = q), quantile(x = sample_treat_1$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_1$ID, t4 = sample_treat_1$V4, umbral = cut(x = sample_treat_1$V4, breaks = c(-Inf, quantile(x = sample_treat_1$V4, probs = q), quantile(x = sample_treat_1$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_1_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t3, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t4, by = "ID")
        sample_treat_1_perdidas <- mutate(sample_treat_1_perdidas, treat = rep(1, nrow(sample_treat_1)))
        
        # Eliminando las observaciones para el tratamiento 2
        
        t1 <- as.data.frame(cbind(ID = sample_treat_2$ID, t1 = sample_treat_2$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_2$ID, t2 = sample_treat_2$V2, umbral = cut(x = sample_treat_2$V2, breaks = c(-Inf, quantile(x = sample_treat_2$V2, probs = q), quantile(x = sample_treat_2$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_2$ID, t3 = sample_treat_2$V3, umbral = cut(x = sample_treat_2$V3, breaks = c(-Inf, quantile(x = sample_treat_2$V3, probs = q), quantile(x = sample_treat_2$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_2$ID, t4 = sample_treat_2$V4, umbral = cut(x = sample_treat_2$V4, breaks = c(-Inf, quantile(x = sample_treat_2$V4, probs = q), quantile(x = sample_treat_2$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_2_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t3, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t4, by = "ID")
        sample_treat_2_perdidas <- mutate(sample_treat_2_perdidas, treat = rep(2, nrow(sample_treat_2_perdidas)))
        
        # Muestra de ambos tratamiento con perdidas
        sample <- bind_rows(sample_treat_1_perdidas,sample_treat_2_perdidas)
        
        sample_long <- reshape(data = sample,varying = 2:(t+1), v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
        sample_long <- arrange(sample_long,ID,tiempo)
        sample_long$tiempo <- as.numeric(sample_long$tiempo)
        sample_long$treat <- as.factor(sample_long$treat)
        sample_long$tiempo <- (sample_long$tiempo-1)/(t-1)
        
        #Modelos
        intercanbiable <- geeglm (formula = yij ~ treat + tiempo + treat * tiempo, id = ID,      data = sample_long, family = gaussian, corstr = "exchangeable")
        AR1            <- geeglm (formula = yij ~ treat + tiempo + treat * tiempo, id = ID,      data = sample_long, family = gaussian, corstr = "ar1")        
        unstructured   <- geeglm (formula = yij ~ treat + tiempo + treat * tiempo, id = ID,      data = sample_long, family = gaussian, corstr = "unstructured")
        intercepto     <- lmer(formula = yij ~ treat + tiempo + treat * tiempo + (1     |ID), data = sample_long, REML = FALSE)
        pen_intercepto <- lmer(formula = yij ~ treat + tiempo + treat * tiempo + (tiempo|ID), data = sample_long, REML = FALSE)
        
        #Completando las matrices con la decisión de la hipótesis
        
        Gee_intercanbiable [i,j]  <- if((1-(pnorm(abs(as.matrix(intercanbiable$coefficients)[4,]/coef(summary(intercanbiable))[4,2]))))*2 < 0.05) 1 else 0
        Gee_AR1            [i,j]  <- if((1-(pnorm(abs(as.matrix(AR1$coefficients)[4,]           /coef(summary(AR1))[4,2]))))*2            < 0.05) 1 else 0
        Gee_unstructured   [i,j]  <- if((1-(pnorm(abs(as.matrix(unstructured$coefficients)[4,]  /coef(summary(unstructured))[4,2]))))*2   < 0.05) 1 else 0
        Mixto_intercepto   [i,j]  <- if(coef(summary(intercepto))[4,5]     < 0.05) 1 else 0
        Mixto_pen_inter    [i,j]  <- if(coef(summary(pen_intercepto))[4,5] < 0.05) 1 else 0
        
        print(c(n[[i]],j, round((sum(is.na(sample_long$yij))/nrow(sample_long))*100, 1) ))
      }}
    }
    if (treat == 2 & t == 8) {
      
      for (i in 1:length(n)) {for(j in 1:repeticiones){ 
        
        sample_treat_1 <- base[sample(x = 1:(nrow(base)/2),            size = n[[i]], replace = FALSE),]
        sample_treat_2 <- base[sample(x = (nrow(base)/2+1):nrow(base), size = n[[i]], replace = FALSE),]
        
        # Eliminando las observaciones para el tratamiento 1
        t1 <- as.data.frame(cbind(ID = sample_treat_1$ID, t1 = sample_treat_1$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_1$ID, t2 = sample_treat_1$V2, umbral = cut(x = sample_treat_1$V2, breaks = c(-Inf, quantile(x = sample_treat_1$V2, probs = q), quantile(x = sample_treat_1$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_1$ID, t3 = sample_treat_1$V3, umbral = cut(x = sample_treat_1$V3, breaks = c(-Inf, quantile(x = sample_treat_1$V3, probs = q), quantile(x = sample_treat_1$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_1$ID, t4 = sample_treat_1$V4, umbral = cut(x = sample_treat_1$V4, breaks = c(-Inf, quantile(x = sample_treat_1$V4, probs = q), quantile(x = sample_treat_1$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        t5 <- as.data.frame(cbind(ID = sample_treat_1$ID, t5 = sample_treat_1$V5, umbral = cut(x = sample_treat_1$V5, breaks = c(-Inf, quantile(x = sample_treat_1$V5, probs = q), quantile(x = sample_treat_1$V5, probs = (1-q)),Inf))))
        t5 <- t5 %>% filter(ID %in% t4$ID)
        n1 <- trunc((nrow(t5)*q)*mi); n2 <- trunc((nrow(t5)*(1-(q*2)))*mc); n3 <- trunc((nrow(t5)*q)*ms)
        borrar_1 <- if((nrow( t5 %>% filter(umbral == 1) )) > n1){true =  (t5 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t5 %>% filter(umbral == 2) )) > n2){true =  (t5 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t5 %>% filter(umbral == 3) )) > n3){true =  (t5 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_1$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_2$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_3$ID)
        
        t6 <- as.data.frame(cbind(ID = sample_treat_1$ID, t6 = sample_treat_1$V6, umbral = cut(x = sample_treat_1$V6, breaks = c(-Inf, quantile(x = sample_treat_1$V6, probs = q), quantile(x = sample_treat_1$V6, probs = (1-q)),Inf))))
        t6 <- t6 %>% filter(ID %in% t5$ID)
        n1 <- trunc((nrow(t6)*q)*mi); n2 <- trunc((nrow(t6)*(1-(q*2)))*mc); n3 <- trunc((nrow(t6)*q)*ms)
        borrar_1 <- if((nrow( t6 %>% filter(umbral == 1) )) > n1){true =  (t6 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t6 %>% filter(umbral == 2) )) > n2){true =  (t6 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t6 %>% filter(umbral == 3) )) > n3){true =  (t6 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_1$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_2$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_3$ID)
        
        t7 <- as.data.frame(cbind(ID = sample_treat_1$ID, t7 = sample_treat_1$V7, umbral = cut(x = sample_treat_1$V7, breaks = c(-Inf, quantile(x = sample_treat_1$V7, probs = q), quantile(x = sample_treat_1$V7, probs = (1-q)),Inf))))
        t7 <- t7 %>% filter(ID %in% t6$ID)
        n1 <- trunc((nrow(t7)*q)*mi); n2 <- trunc((nrow(t7)*(1-(q*2)))*mc); n3 <- trunc((nrow(t7)*q)*ms)
        borrar_1 <- if((nrow( t7 %>% filter(umbral == 1) )) > n1){true =  (t7 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t7 %>% filter(umbral == 2) )) > n2){true =  (t7 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t7 %>% filter(umbral == 3) )) > n3){true =  (t7 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_1$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_2$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_3$ID)
        
        t8 <- as.data.frame(cbind(ID = sample_treat_1$ID, t8 = sample_treat_1$V8, umbral = cut(x = sample_treat_1$V8, breaks = c(-Inf, quantile(x = sample_treat_1$V8, probs = q), quantile(x = sample_treat_1$V8, probs = (1-q)),Inf))))
        t8 <- t8 %>% filter(ID %in% t7$ID)
        n1 <- trunc((nrow(t8)*q)*mi); n2 <- trunc((nrow(t8)*(1-(q*2)))*mc); n3 <- trunc((nrow(t8)*q)*ms)
        borrar_1 <- if((nrow( t8 %>% filter(umbral == 1) )) > n1){true =  (t8 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t8 %>% filter(umbral == 2) )) > n2){true =  (t8 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t8 %>% filter(umbral == 3) )) > n3){true =  (t8 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_1$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_2$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_1_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t3, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t4, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t5, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t6, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t7, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t8, by = "ID")
        sample_treat_1_perdidas <- mutate(sample_treat_1_perdidas, treat = rep(1, nrow(sample_treat_1)))
        
        # Eliminando las observaciones para el tratamiento 2
        
        t1 <- as.data.frame(cbind(ID = sample_treat_2$ID, t1 = sample_treat_2$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_2$ID, t2 = sample_treat_2$V2, umbral = cut(x = sample_treat_2$V2, breaks = c(-Inf, quantile(x = sample_treat_2$V2, probs = q), quantile(x = sample_treat_2$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_2$ID, t3 = sample_treat_2$V3, umbral = cut(x = sample_treat_2$V3, breaks = c(-Inf, quantile(x = sample_treat_2$V3, probs = q), quantile(x = sample_treat_2$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_2$ID, t4 = sample_treat_2$V4, umbral = cut(x = sample_treat_2$V4, breaks = c(-Inf, quantile(x = sample_treat_2$V4, probs = q), quantile(x = sample_treat_2$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        t5 <- as.data.frame(cbind(ID = sample_treat_2$ID, t5 = sample_treat_2$V5, umbral = cut(x = sample_treat_2$V5, breaks = c(-Inf, quantile(x = sample_treat_2$V5, probs = q), quantile(x = sample_treat_2$V5, probs = (1-q)),Inf))))
        t5 <- t5 %>% filter(ID %in% t4$ID)
        n1 <- trunc((nrow(t5)*q)*mi); n2 <- trunc((nrow(t5)*(1-(q*2)))*mc); n3 <- trunc((nrow(t5)*q)*ms)
        borrar_1 <- if((nrow( t5 %>% filter(umbral == 1) )) > n1){true =  (t5 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t5 %>% filter(umbral == 2) )) > n2){true =  (t5 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t5 %>% filter(umbral == 3) )) > n3){true =  (t5 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_1$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_2$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_3$ID)
        
        t6 <- as.data.frame(cbind(ID = sample_treat_2$ID, t6 = sample_treat_2$V6, umbral = cut(x = sample_treat_2$V6, breaks = c(-Inf, quantile(x = sample_treat_2$V6, probs = q), quantile(x = sample_treat_2$V6, probs = (1-q)),Inf))))
        t6 <- t6 %>% filter(ID %in% t5$ID)
        n1 <- trunc((nrow(t6)*q)*mi); n2 <- trunc((nrow(t6)*(1-(q*2)))*mc); n3 <- trunc((nrow(t6)*q)*ms)
        borrar_1 <- if((nrow( t6 %>% filter(umbral == 1) )) > n1){true =  (t6 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t6 %>% filter(umbral == 2) )) > n2){true =  (t6 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t6 %>% filter(umbral == 3) )) > n3){true =  (t6 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_1$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_2$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_3$ID)
        
        t7 <- as.data.frame(cbind(ID = sample_treat_2$ID, t7 = sample_treat_2$V7, umbral = cut(x = sample_treat_2$V7, breaks = c(-Inf, quantile(x = sample_treat_2$V7, probs = q), quantile(x = sample_treat_2$V7, probs = (1-q)),Inf))))
        t7 <- t7 %>% filter(ID %in% t6$ID)
        n1 <- trunc((nrow(t7)*q)*mi); n2 <- trunc((nrow(t7)*(1-(q*2)))*mc); n3 <- trunc((nrow(t7)*q)*ms)
        borrar_1 <- if((nrow( t7 %>% filter(umbral == 1) )) > n1){true =  (t7 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t7 %>% filter(umbral == 2) )) > n2){true =  (t7 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t7 %>% filter(umbral == 3) )) > n3){true =  (t7 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_1$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_2$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_3$ID)
        
        t8 <- as.data.frame(cbind(ID = sample_treat_2$ID, t8 = sample_treat_2$V8, umbral = cut(x = sample_treat_2$V8, breaks = c(-Inf, quantile(x = sample_treat_2$V8, probs = q), quantile(x = sample_treat_2$V8, probs = (1-q)),Inf))))
        t8 <- t8 %>% filter(ID %in% t7$ID)
        n1 <- trunc((nrow(t8)*q)*mi); n2 <- trunc((nrow(t8)*(1-(q*2)))*mc); n3 <- trunc((nrow(t8)*q)*ms)
        borrar_1 <- if((nrow( t8 %>% filter(umbral == 1) )) > n1){true =  (t8 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t8 %>% filter(umbral == 2) )) > n2){true =  (t8 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t8 %>% filter(umbral == 3) )) > n3){true =  (t8 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_1$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_2$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_2_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t3, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t4, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t5, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t6, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t7, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t8, by = "ID")
        sample_treat_2_perdidas <- mutate(sample_treat_2_perdidas, treat = rep(2, nrow(sample_treat_2_perdidas)))
        
        # Muestra de ambos tratamiento con perdidas
        sample <- bind_rows(sample_treat_1_perdidas,sample_treat_2_perdidas)
        
        sample_long <- reshape(data = sample,varying = 2:(t+1), v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
        sample_long <- arrange(sample_long,ID,tiempo)
        sample_long$tiempo <- as.numeric(sample_long$tiempo)
        sample_long$treat <- as.factor(sample_long$treat)
        sample_long$tiempo <- (sample_long$tiempo-1)/(t-1)
        
        #Modelos
        intercanbiable <- geeglm (formula = yij ~ treat + tiempo + treat * tiempo, id = ID,      data = sample_long, family = gaussian, corstr = "exchangeable")
        AR1            <- geeglm (formula = yij ~ treat + tiempo + treat * tiempo, id = ID,      data = sample_long, family = gaussian, corstr = "ar1")        
        unstructured   <- geeglm (formula = yij ~ treat + tiempo + treat * tiempo, id = ID,      data = sample_long, family = gaussian, corstr = "unstructured")
        intercepto     <- lmer(formula = yij ~ treat + tiempo + treat * tiempo + (1     |ID), data = sample_long, REML = FALSE)
        pen_intercepto <- lmer(formula = yij ~ treat + tiempo + treat * tiempo + (tiempo|ID), data = sample_long, REML = FALSE)
        
        #Completando las matrices con la decisión de la hipótesis
        
        Gee_intercanbiable [i,j]  <- if((1-(pnorm(abs(as.matrix(intercanbiable$coefficients)[4,]/coef(summary(intercanbiable))[4,2]))))*2 < 0.05) 1 else 0
        Gee_AR1            [i,j]  <- if((1-(pnorm(abs(as.matrix(AR1$coefficients)[4,]           /coef(summary(AR1))[4,2]))))*2            < 0.05) 1 else 0
        Gee_unstructured   [i,j]  <- if((1-(pnorm(abs(as.matrix(unstructured$coefficients)[4,]  /coef(summary(unstructured))[4,2]))))*2   < 0.05) 1 else 0
        Mixto_intercepto   [i,j]  <- if(coef(summary(intercepto))[4,5]     < 0.05) 1 else 0
        Mixto_pen_inter    [i,j]  <- if(coef(summary(pen_intercepto))[4,5] < 0.05) 1 else 0
        
        print(c(n[[i]],j, round((sum(is.na(sample_long$yij))/nrow(sample_long))*100, 1) ))
      }}
    }
    if (treat == 3 & t == 4) {
      for (i in 1:length(n)) {for(j in 1:repeticiones){ 
        
        sample_treat_1 <- base[sample(x = 1:(nrow(base)/3),                    size = n[[i]], replace = FALSE),]
        sample_treat_2 <- base[sample(x = (nrow(base)/3+1):((nrow(base)/3)*2), size = n[[i]], replace = FALSE),]
        sample_treat_3 <- base[sample(x = ((nrow(base)/3)*2+1):nrow(base),     size = n[[i]], replace = FALSE),]
        
        # Eliminando las observaciones para el tratamiento 1
        t1 <- as.data.frame(cbind(ID = sample_treat_1$ID, t1 = sample_treat_1$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_1$ID, t2 = sample_treat_1$V2, umbral = cut(x = sample_treat_1$V2, breaks = c(-Inf, quantile(x = sample_treat_1$V2, probs = q), quantile(x = sample_treat_1$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_1$ID, t3 = sample_treat_1$V3, umbral = cut(x = sample_treat_1$V3, breaks = c(-Inf, quantile(x = sample_treat_1$V3, probs = q), quantile(x = sample_treat_1$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_1$ID, t4 = sample_treat_1$V4, umbral = cut(x = sample_treat_1$V4, breaks = c(-Inf, quantile(x = sample_treat_1$V4, probs = q), quantile(x = sample_treat_1$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_1_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t3, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t4, by = "ID")
        sample_treat_1_perdidas <- mutate(sample_treat_1_perdidas, treat = rep(1, nrow(sample_treat_1)))
        
        # Eliminando las observaciones para el tratamiento 2
        
        t1 <- as.data.frame(cbind(ID = sample_treat_2$ID, t1 = sample_treat_2$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_2$ID, t2 = sample_treat_2$V2, umbral = cut(x = sample_treat_2$V2, breaks = c(-Inf, quantile(x = sample_treat_2$V2, probs = q), quantile(x = sample_treat_2$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_2$ID, t3 = sample_treat_2$V3, umbral = cut(x = sample_treat_2$V3, breaks = c(-Inf, quantile(x = sample_treat_2$V3, probs = q), quantile(x = sample_treat_2$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_2$ID, t4 = sample_treat_2$V4, umbral = cut(x = sample_treat_2$V4, breaks = c(-Inf, quantile(x = sample_treat_2$V4, probs = q), quantile(x = sample_treat_2$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_2_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t3, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t4, by = "ID")
        sample_treat_2_perdidas <- mutate(sample_treat_2_perdidas, treat = rep(2, nrow(sample_treat_2_perdidas)))
        
        # Eliminando las observaciones para el tratamiento 3
        t1 <- as.data.frame(cbind(ID = sample_treat_3$ID, t1 = sample_treat_3$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_3$ID, t2 = sample_treat_3$V2, umbral = cut(x = sample_treat_3$V2, breaks = c(-Inf, quantile(x = sample_treat_3$V2, probs = q), quantile(x = sample_treat_3$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_3$ID, t3 = sample_treat_3$V3, umbral = cut(x = sample_treat_3$V3, breaks = c(-Inf, quantile(x = sample_treat_3$V3, probs = q), quantile(x = sample_treat_3$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_3$ID, t4 = sample_treat_3$V4, umbral = cut(x = sample_treat_3$V4, breaks = c(-Inf, quantile(x = sample_treat_3$V4, probs = q), quantile(x = sample_treat_3$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_3_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_3_perdidas <- full_join(x = sample_treat_3_perdidas, y = t3, by = "ID")
        sample_treat_3_perdidas <- full_join(x = sample_treat_3_perdidas, y = t4, by = "ID")
        sample_treat_3_perdidas <- mutate(sample_treat_3_perdidas, treat = rep(3, nrow(sample_treat_3)))
        
        # Muestra de ambos tratamiento con perdidas
        sample <- bind_rows(sample_treat_1_perdidas,sample_treat_2_perdidas,sample_treat_3_perdidas)
        
        sample_long <- reshape(data = sample,varying = 2:(t+1), v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
        sample_long <- arrange(sample_long,ID,tiempo)
        sample_long$tiempo <- as.numeric(sample_long$tiempo)
        sample_long$treat <- as.factor(sample_long$treat)
        sample_long$tiempo <- (sample_long$tiempo-1)/(t-1)
        
        #Modelos
        intercanbiable <- geeglm(formula = yij ~ treat + tiempo + treat * tiempo,     id = ID,    data = sample_long, family = gaussian, corstr = "exchangeable")
        AR1            <- geeglm(formula = yij ~ treat + tiempo + treat * tiempo,     id = ID,    data = sample_long, family = gaussian, corstr = "ar1")
        unstructured   <- geeglm(formula = yij ~ treat + tiempo + treat * tiempo,     id = ID,    data = sample_long, family = gaussian, corstr = "unstructured")
        intercepto     <- lmer  (formula = yij ~ treat + tiempo + treat * tiempo + (1     |ID),   data = sample_long, REML = FALSE)
        pen_intercepto <- lmer  (formula = yij ~ treat + tiempo + treat * tiempo + (tiempo|ID),   data = sample_long, REML = FALSE)
        
        #Completando las matrices con la decisión de la hipótesis
        Gee_intercanbiable   [i,j]  <- if((1-(pnorm( abs( coef(summary(intercanbiable)) [5,1]  /   coef(summary(intercanbiable))[5,2] ))))*2 < 0.05) 1 else 0
        Gee_AR1              [i,j]  <- if((1-(pnorm( abs( coef(summary(AR1))            [5,1]  /   coef(summary(AR1))           [5,2] ))))*2 < 0.05) 1 else 0
        Gee_unstructured     [i,j]  <- if((1-(pnorm( abs( coef(summary(unstructured))   [5,1]  /   coef(summary(unstructured))  [5,2] ))))*2 < 0.05) 1 else 0
        Mixto_intercepto     [i,j]  <- if(coef(summary(intercepto))    [5,5]                                                                 < 0.05) 1 else 0
        Mixto_pen_inter      [i,j]  <- if(coef(summary(pen_intercepto))[5,5]                                                                 < 0.05) 1 else 0
        
        Gee_intercanbiable_3 [i,j]  <- if((1-(pnorm( abs( coef(summary(intercanbiable)) [6,1]  /   coef(summary(intercanbiable))[6,2] ))))*2 < 0.05) 1 else 0
        Gee_AR1_3            [i,j]  <- if((1-(pnorm( abs( coef(summary(AR1))            [6,1]  /   coef(summary(AR1))           [6,2] ))))*2 < 0.05) 1 else 0
        Gee_unstructured_3   [i,j]  <- if((1-(pnorm( abs( coef(summary(unstructured))   [6,1]  /   coef(summary(unstructured))  [6,2] ))))*2 < 0.05) 1 else 0
        Mixto_intercepto_3   [i,j]  <- if(coef(summary(intercepto))    [6,5]                                                                 < 0.05) 1 else 0
        Mixto_pen_inter_3    [i,j]  <- if(coef(summary(pen_intercepto))[6,5]                                                                 < 0.05) 1 else 0
        
        print(c(n[[i]],j, round((sum(is.na(sample_long$yij))/nrow(sample_long))*100, 1) ))
        
      }}
    }
    if (treat == 3 & t == 8) {
      
      for (i in 1:length(n)) {for(j in 1:repeticiones){ 
        
        sample_treat_1 <- base[sample(x = 1:(nrow(base)/3),                    size = n[[i]], replace = FALSE),]
        sample_treat_2 <- base[sample(x = (nrow(base)/3+1):((nrow(base)/3)*2), size = n[[i]], replace = FALSE),]
        sample_treat_3 <- base[sample(x = ((nrow(base)/3)*2+1):nrow(base),     size = n[[i]], replace = FALSE),]
        
        # Eliminando las observaciones para el tratamiento 1
        t1 <- as.data.frame(cbind(ID = sample_treat_1$ID, t1 = sample_treat_1$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_1$ID, t2 = sample_treat_1$V2, umbral = cut(x = sample_treat_1$V2, breaks = c(-Inf, quantile(x = sample_treat_1$V2, probs = q), quantile(x = sample_treat_1$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_1$ID, t3 = sample_treat_1$V3, umbral = cut(x = sample_treat_1$V3, breaks = c(-Inf, quantile(x = sample_treat_1$V3, probs = q), quantile(x = sample_treat_1$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_1$ID, t4 = sample_treat_1$V4, umbral = cut(x = sample_treat_1$V4, breaks = c(-Inf, quantile(x = sample_treat_1$V4, probs = q), quantile(x = sample_treat_1$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        t5 <- as.data.frame(cbind(ID = sample_treat_1$ID, t5 = sample_treat_1$V5, umbral = cut(x = sample_treat_1$V5, breaks = c(-Inf, quantile(x = sample_treat_1$V5, probs = q), quantile(x = sample_treat_1$V5, probs = (1-q)),Inf))))
        t5 <- t5 %>% filter(ID %in% t4$ID)
        n1 <- trunc((nrow(t5)*q)*mi); n2 <- trunc((nrow(t5)*(1-(q*2)))*mc); n3 <- trunc((nrow(t5)*q)*ms)
        borrar_1 <- if((nrow( t5 %>% filter(umbral == 1) )) > n1){true =  (t5 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t5 %>% filter(umbral == 2) )) > n2){true =  (t5 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t5 %>% filter(umbral == 3) )) > n3){true =  (t5 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_1$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_2$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_3$ID)
        
        t6 <- as.data.frame(cbind(ID = sample_treat_1$ID, t6 = sample_treat_1$V6, umbral = cut(x = sample_treat_1$V6, breaks = c(-Inf, quantile(x = sample_treat_1$V6, probs = q), quantile(x = sample_treat_1$V6, probs = (1-q)),Inf))))
        t6 <- t6 %>% filter(ID %in% t5$ID)
        n1 <- trunc((nrow(t6)*q)*mi); n2 <- trunc((nrow(t6)*(1-(q*2)))*mc); n3 <- trunc((nrow(t6)*q)*ms)
        borrar_1 <- if((nrow( t6 %>% filter(umbral == 1) )) > n1){true =  (t6 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t6 %>% filter(umbral == 2) )) > n2){true =  (t6 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t6 %>% filter(umbral == 3) )) > n3){true =  (t6 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_1$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_2$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_3$ID)
        
        t7 <- as.data.frame(cbind(ID = sample_treat_1$ID, t7 = sample_treat_1$V7, umbral = cut(x = sample_treat_1$V7, breaks = c(-Inf, quantile(x = sample_treat_1$V7, probs = q), quantile(x = sample_treat_1$V7, probs = (1-q)),Inf))))
        t7 <- t7 %>% filter(ID %in% t6$ID)
        n1 <- trunc((nrow(t7)*q)*mi); n2 <- trunc((nrow(t7)*(1-(q*2)))*mc); n3 <- trunc((nrow(t7)*q)*ms)
        borrar_1 <- if((nrow( t7 %>% filter(umbral == 1) )) > n1){true =  (t7 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t7 %>% filter(umbral == 2) )) > n2){true =  (t7 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t7 %>% filter(umbral == 3) )) > n3){true =  (t7 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_1$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_2$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_3$ID)
        
        t8 <- as.data.frame(cbind(ID = sample_treat_1$ID, t8 = sample_treat_1$V8, umbral = cut(x = sample_treat_1$V8, breaks = c(-Inf, quantile(x = sample_treat_1$V8, probs = q), quantile(x = sample_treat_1$V8, probs = (1-q)),Inf))))
        t8 <- t8 %>% filter(ID %in% t7$ID)
        n1 <- trunc((nrow(t8)*q)*mi); n2 <- trunc((nrow(t8)*(1-(q*2)))*mc); n3 <- trunc((nrow(t8)*q)*ms)
        borrar_1 <- if((nrow( t8 %>% filter(umbral == 1) )) > n1){true =  (t8 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t8 %>% filter(umbral == 2) )) > n2){true =  (t8 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t8 %>% filter(umbral == 3) )) > n3){true =  (t8 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_1$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_2$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_1_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t3, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t4, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t5, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t6, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t7, by = "ID")
        sample_treat_1_perdidas <- full_join(x = sample_treat_1_perdidas, y = t8, by = "ID")
        sample_treat_1_perdidas <- mutate(sample_treat_1_perdidas, treat = rep(1, nrow(sample_treat_1)))
        
        # Eliminando las observaciones para el tratamiento 2
        
        t1 <- as.data.frame(cbind(ID = sample_treat_2$ID, t1 = sample_treat_2$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_2$ID, t2 = sample_treat_2$V2, umbral = cut(x = sample_treat_2$V2, breaks = c(-Inf, quantile(x = sample_treat_2$V2, probs = q), quantile(x = sample_treat_2$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_2$ID, t3 = sample_treat_2$V3, umbral = cut(x = sample_treat_2$V3, breaks = c(-Inf, quantile(x = sample_treat_2$V3, probs = q), quantile(x = sample_treat_2$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_2$ID, t4 = sample_treat_2$V4, umbral = cut(x = sample_treat_2$V4, breaks = c(-Inf, quantile(x = sample_treat_2$V4, probs = q), quantile(x = sample_treat_2$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        t5 <- as.data.frame(cbind(ID = sample_treat_2$ID, t5 = sample_treat_2$V5, umbral = cut(x = sample_treat_2$V5, breaks = c(-Inf, quantile(x = sample_treat_2$V5, probs = q), quantile(x = sample_treat_2$V5, probs = (1-q)),Inf))))
        t5 <- t5 %>% filter(ID %in% t4$ID)
        n1 <- trunc((nrow(t5)*q)*mi); n2 <- trunc((nrow(t5)*(1-(q*2)))*mc); n3 <- trunc((nrow(t5)*q)*ms)
        borrar_1 <- if((nrow( t5 %>% filter(umbral == 1) )) > n1){true =  (t5 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t5 %>% filter(umbral == 2) )) > n2){true =  (t5 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t5 %>% filter(umbral == 3) )) > n3){true =  (t5 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_1$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_2$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_3$ID)
        
        t6 <- as.data.frame(cbind(ID = sample_treat_2$ID, t6 = sample_treat_2$V6, umbral = cut(x = sample_treat_2$V6, breaks = c(-Inf, quantile(x = sample_treat_2$V6, probs = q), quantile(x = sample_treat_2$V6, probs = (1-q)),Inf))))
        t6 <- t6 %>% filter(ID %in% t5$ID)
        n1 <- trunc((nrow(t6)*q)*mi); n2 <- trunc((nrow(t6)*(1-(q*2)))*mc); n3 <- trunc((nrow(t6)*q)*ms)
        borrar_1 <- if((nrow( t6 %>% filter(umbral == 1) )) > n1){true =  (t6 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t6 %>% filter(umbral == 2) )) > n2){true =  (t6 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t6 %>% filter(umbral == 3) )) > n3){true =  (t6 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_1$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_2$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_3$ID)
        
        t7 <- as.data.frame(cbind(ID = sample_treat_2$ID, t7 = sample_treat_2$V7, umbral = cut(x = sample_treat_2$V7, breaks = c(-Inf, quantile(x = sample_treat_2$V7, probs = q), quantile(x = sample_treat_2$V7, probs = (1-q)),Inf))))
        t7 <- t7 %>% filter(ID %in% t6$ID)
        n1 <- trunc((nrow(t7)*q)*mi); n2 <- trunc((nrow(t7)*(1-(q*2)))*mc); n3 <- trunc((nrow(t7)*q)*ms)
        borrar_1 <- if((nrow( t7 %>% filter(umbral == 1) )) > n1){true =  (t7 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t7 %>% filter(umbral == 2) )) > n2){true =  (t7 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t7 %>% filter(umbral == 3) )) > n3){true =  (t7 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_1$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_2$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_3$ID)
        
        t8 <- as.data.frame(cbind(ID = sample_treat_2$ID, t8 = sample_treat_2$V8, umbral = cut(x = sample_treat_2$V8, breaks = c(-Inf, quantile(x = sample_treat_2$V8, probs = q), quantile(x = sample_treat_2$V8, probs = (1-q)),Inf))))
        t8 <- t8 %>% filter(ID %in% t7$ID)
        n1 <- trunc((nrow(t8)*q)*mi); n2 <- trunc((nrow(t8)*(1-(q*2)))*mc); n3 <- trunc((nrow(t8)*q)*ms)
        borrar_1 <- if((nrow( t8 %>% filter(umbral == 1) )) > n1){true =  (t8 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t8 %>% filter(umbral == 2) )) > n2){true =  (t8 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t8 %>% filter(umbral == 3) )) > n3){true =  (t8 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_1$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_2$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_2_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t3, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t4, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t5, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t6, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t7, by = "ID")
        sample_treat_2_perdidas <- full_join(x = sample_treat_2_perdidas, y = t8, by = "ID")
        sample_treat_2_perdidas <- mutate(sample_treat_2_perdidas, treat = rep(2, nrow(sample_treat_2_perdidas)))
        
        # Eliminando las observaciones para el tratamiento 3
        t1 <- as.data.frame(cbind(ID = sample_treat_3$ID, t1 = sample_treat_3$V1) )
        
        t2 <- as.data.frame(cbind(ID = sample_treat_3$ID, t2 = sample_treat_3$V2, umbral = cut(x = sample_treat_3$V2, breaks = c(-Inf, quantile(x = sample_treat_3$V2, probs = q), quantile(x = sample_treat_3$V2, probs = (1-q)),Inf))))
        n1 <- trunc((nrow(t2)*q)*mi); n2 <- trunc((nrow(t2)*(1-(q*2)))*mc); n3 <- trunc((nrow(t2)*q)*ms)
        borrar_1 <- if((nrow( t2 %>% filter(umbral == 1) )) > n1){true =  (t2 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t2 %>% filter(umbral == 2) )) > n2){true =  (t2 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t2 %>% filter(umbral == 3) )) > n3){true =  (t2 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_1$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_2$ID)
        t2 <- t2 %>% select(ID, t2) %>% filter(!ID %in% borrar_3$ID)
        
        t3 <- as.data.frame(cbind(ID = sample_treat_3$ID, t3 = sample_treat_3$V3, umbral = cut(x = sample_treat_3$V3, breaks = c(-Inf, quantile(x = sample_treat_3$V3, probs = q), quantile(x = sample_treat_3$V3, probs = (1-q)),Inf))))
        t3 <- t3 %>% filter(ID %in% t2$ID)
        n1 <- trunc((nrow(t3)*q)*mi); n2 <- trunc((nrow(t3)*(1-(q*2)))*mc); n3 <- trunc((nrow(t3)*q)*ms)
        borrar_1 <- if((nrow( t3 %>% filter(umbral == 1) )) > n1){true =  (t3 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t3 %>% filter(umbral == 2) )) > n2){true =  (t3 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t3 %>% filter(umbral == 3) )) > n3){true =  (t3 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_1$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_2$ID)
        t3 <- t3 %>% select(ID, t3) %>% filter(!ID %in% borrar_3$ID)
        
        t4 <- as.data.frame(cbind(ID = sample_treat_3$ID, t4 = sample_treat_3$V4, umbral = cut(x = sample_treat_3$V4, breaks = c(-Inf, quantile(x = sample_treat_3$V4, probs = q), quantile(x = sample_treat_3$V4, probs = (1-q)),Inf))))
        t4 <- t4 %>% filter(ID %in% t3$ID)
        n1 <- trunc((nrow(t4)*q)*mi); n2 <- trunc((nrow(t4)*(1-(q*2)))*mc); n3 <- trunc((nrow(t4)*q)*ms)
        borrar_1 <- if((nrow( t4 %>% filter(umbral == 1) )) > n1){true =  (t4 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t4 %>% filter(umbral == 2) )) > n2){true =  (t4 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t4 %>% filter(umbral == 3) )) > n3){true =  (t4 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_1$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_2$ID)
        t4 <- t4 %>% select(ID, t4) %>% filter(!ID %in% borrar_3$ID)
        
        t5 <- as.data.frame(cbind(ID = sample_treat_3$ID, t5 = sample_treat_3$V5, umbral = cut(x = sample_treat_3$V5, breaks = c(-Inf, quantile(x = sample_treat_3$V5, probs = q), quantile(x = sample_treat_3$V5, probs = (1-q)),Inf))))
        t5 <- t5 %>% filter(ID %in% t4$ID)
        n1 <- trunc((nrow(t5)*q)*mi); n2 <- trunc((nrow(t5)*(1-(q*2)))*mc); n3 <- trunc((nrow(t5)*q)*ms)
        borrar_1 <- if((nrow( t5 %>% filter(umbral == 1) )) > n1){true =  (t5 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t5 %>% filter(umbral == 2) )) > n2){true =  (t5 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t5 %>% filter(umbral == 3) )) > n3){true =  (t5 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_1$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_2$ID)
        t5 <- t5 %>% select(ID, t5) %>% filter(!ID %in% borrar_3$ID)
        
        t6 <- as.data.frame(cbind(ID = sample_treat_3$ID, t6 = sample_treat_3$V6, umbral = cut(x = sample_treat_3$V6, breaks = c(-Inf, quantile(x = sample_treat_3$V6, probs = q), quantile(x = sample_treat_3$V6, probs = (1-q)),Inf))))
        t6 <- t6 %>% filter(ID %in% t5$ID)
        n1 <- trunc((nrow(t6)*q)*mi); n2 <- trunc((nrow(t6)*(1-(q*2)))*mc); n3 <- trunc((nrow(t6)*q)*ms)
        borrar_1 <- if((nrow( t6 %>% filter(umbral == 1) )) > n1){true =  (t6 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t6 %>% filter(umbral == 2) )) > n2){true =  (t6 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t6 %>% filter(umbral == 3) )) > n3){true =  (t6 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_1$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_2$ID)
        t6 <- t6 %>% select(ID, t6) %>% filter(!ID %in% borrar_3$ID)
        
        t7 <- as.data.frame(cbind(ID = sample_treat_3$ID, t7 = sample_treat_3$V7, umbral = cut(x = sample_treat_3$V7, breaks = c(-Inf, quantile(x = sample_treat_3$V7, probs = q), quantile(x = sample_treat_3$V7, probs = (1-q)),Inf))))
        t7 <- t7 %>% filter(ID %in% t6$ID)
        n1 <- trunc((nrow(t7)*q)*mi); n2 <- trunc((nrow(t7)*(1-(q*2)))*mc); n3 <- trunc((nrow(t7)*q)*ms)
        borrar_1 <- if((nrow( t7 %>% filter(umbral == 1) )) > n1){true =  (t7 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t7 %>% filter(umbral == 2) )) > n2){true =  (t7 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t7 %>% filter(umbral == 3) )) > n3){true =  (t7 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_1$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_2$ID)
        t7 <- t7 %>% select(ID, t7) %>% filter(!ID %in% borrar_3$ID)
        
        t8 <- as.data.frame(cbind(ID = sample_treat_3$ID, t8 = sample_treat_3$V8, umbral = cut(x = sample_treat_3$V8, breaks = c(-Inf, quantile(x = sample_treat_3$V8, probs = q), quantile(x = sample_treat_3$V8, probs = (1-q)),Inf))))
        t8 <- t8 %>% filter(ID %in% t7$ID)
        n1 <- trunc((nrow(t8)*q)*mi); n2 <- trunc((nrow(t8)*(1-(q*2)))*mc); n3 <- trunc((nrow(t8)*q)*ms)
        borrar_1 <- if((nrow( t8 %>% filter(umbral == 1) )) > n1){true =  (t8 %>% filter(umbral == 1) %>% sample_n(size = n1, replace=FALSE))[1]}
        borrar_2 <- if((nrow( t8 %>% filter(umbral == 2) )) > n2){true =  (t8 %>% filter(umbral == 2) %>% sample_n(size = n2, replace=FALSE))[1]}
        borrar_3 <- if((nrow( t8 %>% filter(umbral == 3) )) > n3){true =  (t8 %>% filter(umbral == 3) %>% sample_n(size = n3, replace=FALSE))[1]}
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_1$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_2$ID)
        t8 <- t8 %>% select(ID, t8) %>% filter(!ID %in% borrar_3$ID)
        
        sample_treat_3_perdidas <- full_join(x = t1,    y = t2, by = "ID")
        sample_treat_3_perdidas <- full_join(x = sample_treat_3_perdidas, y = t3, by = "ID")
        sample_treat_3_perdidas <- full_join(x = sample_treat_3_perdidas, y = t4, by = "ID")
        sample_treat_3_perdidas <- full_join(x = sample_treat_3_perdidas, y = t5, by = "ID")
        sample_treat_3_perdidas <- full_join(x = sample_treat_3_perdidas, y = t6, by = "ID")
        sample_treat_3_perdidas <- full_join(x = sample_treat_3_perdidas, y = t7, by = "ID")
        sample_treat_3_perdidas <- full_join(x = sample_treat_3_perdidas, y = t8, by = "ID")
        sample_treat_3_perdidas <- mutate(sample_treat_3_perdidas, treat = rep(3, nrow(sample_treat_3)))
        
        # Muestra de ambos tratamiento con perdidas
        sample <- bind_rows(sample_treat_1_perdidas,sample_treat_2_perdidas,sample_treat_3_perdidas)
        
        sample_long <- reshape(data = sample,varying = 2:(t+1), v.names = "yij", timevar= "tiempo", idvar = "ID", direction = "long")
        sample_long <- arrange(sample_long,ID,tiempo)
        sample_long$tiempo <- as.numeric(sample_long$tiempo)
        sample_long$treat <- as.factor(sample_long$treat)
        sample_long$tiempo <- (sample_long$tiempo-1)/(t-1)
        
        #Modelos
        intercanbiable <- geeglm (formula = yij ~ treat + tiempo + treat * tiempo, id = ID,      data = sample_long, family = gaussian, corstr = "exchangeable")
        AR1            <- geeglm (formula = yij ~ treat + tiempo + treat * tiempo, id = ID,      data = sample_long, family = gaussian, corstr = "ar1")        
        unstructured   <- geeglm (formula = yij ~ treat + tiempo + treat * tiempo, id = ID,      data = sample_long, family = gaussian, corstr = "unstructured")
        intercepto     <- lmer(formula = yij ~ treat + tiempo + treat * tiempo + (1     |ID), data = sample_long, REML = FALSE)
        pen_intercepto <- lmer(formula = yij ~ treat + tiempo + treat * tiempo + (tiempo|ID), data = sample_long, REML = FALSE)
        
        #Completando las matrices con la decisión de la hipótesis
        
        Gee_intercanbiable   [i,j]  <- if((1-(pnorm( abs( coef(summary(intercanbiable)) [5,1]  /   coef(summary(intercanbiable))[5,2] ))))*2 < 0.05) 1 else 0
        Gee_AR1              [i,j]  <- if((1-(pnorm( abs( coef(summary(AR1))            [5,1]  /   coef(summary(AR1))           [5,2] ))))*2 < 0.05) 1 else 0
        Gee_unstructured     [i,j]  <- if((1-(pnorm( abs( coef(summary(unstructured))   [5,1]  /   coef(summary(unstructured))  [5,2] ))))*2 < 0.05) 1 else 0
        Mixto_intercepto     [i,j]  <- if(coef(summary(intercepto))    [5,5]                                                                 < 0.05) 1 else 0
        Mixto_pen_inter      [i,j]  <- if(coef(summary(pen_intercepto))[5,5]                                                                 < 0.05) 1 else 0
        
        Gee_intercanbiable_3 [i,j]  <- if((1-(pnorm( abs( coef(summary(intercanbiable)) [6,1]  /   coef(summary(intercanbiable))[6,2] ))))*2 < 0.05) 1 else 0
        Gee_AR1_3            [i,j]  <- if((1-(pnorm( abs( coef(summary(AR1))            [6,1]  /   coef(summary(AR1))           [6,2] ))))*2 < 0.05) 1 else 0
        Gee_unstructured_3   [i,j]  <- if((1-(pnorm( abs( coef(summary(unstructured))   [6,1]  /   coef(summary(unstructured))  [6,2] ))))*2 < 0.05) 1 else 0
        Mixto_intercepto_3   [i,j]  <- if(coef(summary(intercepto))    [6,5]                                                                 < 0.05) 1 else 0
        Mixto_pen_inter_3    [i,j]  <- if(coef(summary(pen_intercepto))[6,5]                                                                 < 0.05) 1 else 0
        
        print(c(n[[i]],j, round((sum(is.na(sample_long$yij))/nrow(sample_long))*100, 1) ))
        
      }}
    }
    
    #Completando las matrices con la decisión de la hipótesis para 2 y 3 brazos de tratamientos
    if (treat == 2) {
      #Base de datos para 2 tratamientos
      Gee_inter  <- as.matrix(apply(X = Gee_intercanbiable, MARGIN = 1, FUN = mean))
      Gee_AR     <- as.matrix(apply(X = Gee_AR1,            MARGIN = 1, FUN = mean))
      Gee_unst   <- as.matrix(apply(X = Gee_unstructured,   MARGIN = 1, FUN = mean))
      Mixto_inte <- as.matrix(apply(X = Mixto_intercepto,   MARGIN = 1, FUN = mean))
      Mixto_pen_ <- as.matrix(apply(X = Mixto_pen_inter,    MARGIN = 1, FUN = mean))
      
      Base <- as.data.frame(cbind(n = n ))
      Base <- mutate(Base,Gee_inter)
      Base <- mutate(Base,Gee_AR)
      Base <- mutate(Base,Gee_unst)
      Base <- mutate(Base,Mixto_inte)
      Base <- mutate(Base,Mixto_pen_)
    }
    if (treat == 3) {
      #Base de datos para 3 tratamientos
      Gee_inter     <- as.matrix(apply(X = Gee_intercanbiable,   MARGIN = 1, FUN = mean))
      Gee_AR        <- as.matrix(apply(X = Gee_AR1,              MARGIN = 1, FUN = mean))
      Gee_unst      <- as.matrix(apply(X = Gee_unstructured,     MARGIN = 1, FUN = mean))
      Mixto_inte    <- as.matrix(apply(X = Mixto_intercepto,     MARGIN = 1, FUN = mean))
      Mixto_pen_    <- as.matrix(apply(X = Mixto_pen_inter,      MARGIN = 1, FUN = mean))
      
      Gee_inter_3   <- as.matrix(apply(X = Gee_intercanbiable_3, MARGIN = 1, FUN = mean))
      Gee_AR_3      <- as.matrix(apply(X = Gee_AR1_3,            MARGIN = 1, FUN = mean))
      Gee_unst_3    <- as.matrix(apply(X = Gee_unstructured_3,   MARGIN = 1, FUN = mean))
      Mixto_inte_3  <- as.matrix(apply(X = Mixto_intercepto_3,   MARGIN = 1, FUN = mean))
      Mixto_pen_3   <- as.matrix(apply(X = Mixto_pen_inter_3,    MARGIN = 1, FUN = mean))
      
      Base <- as.data.frame(cbind(ID = n))
      
      Base <- mutate(Base,Gee_inter)
      Base <- mutate(Base,Gee_AR)
      Base <- mutate(Base,Gee_unst)
      Base <- mutate(Base,Mixto_inte)
      Base <- mutate(Base,Mixto_pen_)
      
      Base <- mutate(Base,Gee_inter_3)
      Base <- mutate(Base,Gee_AR_3)
      Base <- mutate(Base,Gee_unst_3)
      Base <- mutate(Base,Mixto_inte_3)
      Base <- mutate(Base,Mixto_pen_3)
    }
    
    # Capturando la hora de término del la función
    Fin <- DescTools::Now()
    
    # Calculando la duración
    Duración <- Fin - Inicio; print(Duración)
    
    # Retornando los resultados
    return(Base = Base)
  } 
}# Funciones necesarias para la simulación

  ##Escenario_31
  población <- Pob_ECA(n = 1000000, 
                       t = 4, 
                       b0 = 44.9, 
                       b1 = 0.11, 
                       b2 = -2.30, 
                       b3 = 0, 
                       b4 = -2.30, 
                       b5 = 0, 
                       var_v0i = 36, 
                       var_v1i = 100, 
                       cov_v0iv1i = 45, 
                       var_eij = 22.09)
  
  base <- población$Treat_3$yij_3
  print(DescTools::Now() - Inicio)
  Escenario_31 <- Comp_modelos_par(base = base, treat = 3, n = c(10,30,50,70,100,250,500,750,1000,1250,1500,1750,2000), repeticiones = 10000, t = 4);Escenario_31;"Escenario_31"
  print(DescTools::Now() - Inicio)
  write.table(Escenario_31, "Escenario_31.txt")  
  
      
    
