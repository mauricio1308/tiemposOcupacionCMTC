#' Cargue la matriz generadora, el t a evaluar (debe coincidir con las tasas de la matriz Q) y el alpha o vector de probabilidad inicial.
#'
#'Esta función permite calcular la aproximación de la matriz de ocupación dada una matriz generadora (matrizQ), un tiempo a evaluar (t) y el vector de probabilidad inicial (alpha)
#'
#'@param matrizQ Matriz generadora de la CMTC
#'@param t tiempo a evaluar (debe estar en las mismas unidades que están las tasas de la matriz generadora)
#'@param alpha Vector de probabilidad inicial
#'
#'@return El número de transiciones aproximadas que transcurren en el intervalo de tiempo dado por parámetro y la matriz de ocupación (M^t)
#'
#'@export

tiemposOcupacionCMTC <- function(matrizQ,t,alpha){
  library(markovchain)
  L <- -1/diag(matrizQ)
  matrizP <- generatorToTransitionMatrix(matrizQ)
  matrizMultiplicar <- matrizP
  contador <- 1
  valorEsperadoTiempo <- alpha%*%L
  
  while(valorEsperadoTiempo<t){
    contador=contador + 1
    valorEsperadoTiempo <- valorEsperadoTiempo +alpha%*%matrizMultiplicar%*%L
    #print(paste0(contador,"----",valorEsperadoTiempo))
    matrizMultiplicar <- matrizMultiplicar%*%matrizP
  }
  print(paste0(contador,"------",valorEsperadoTiempo))
  matrizOcupacion <- diag(nrow(matrizQ)) + matrizP
  resultado <- matrizP 
  
  for(i in 2:(contador-1)){
    resultado <- resultado%*%matrizP
    matrizOcupacion <- matrizOcupacion + resultado
  }
  
  M_t <- matrizOcupacion*L
  return(M_t)
}