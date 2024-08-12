
# here is the age-structured 3-strain model, including heterotypic immunity
# here is the age-structured 3-strain model, including heterotypic immunity

rm(list=ls())

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig5/sim-new/"))
#set wdset wd
#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4")


library(plyr)
library(dplyr)
library(matrixcalc)
library(Matrix)
library(ggplot2)
library(mvtnorm)
library(reshape2)


#helper functions
#need to modify pre
buildTMat_four_strain <- function(c, Npop, age.classes, surv.biwk, wane_hetero, age.brk, foi, age.mult.foi, recov, sigma, age.rate){
  
  s <- nage <- length(age.classes)
  
  
  
  if (length(foi)==1) foi <- rep(foi,nage) # 
  if (length(foi)==s) foi <- foi
  if (length(foi) > 1 & length(foi)<s){
    foi.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(foi)){
      foi.list[[i]] = rep(foi[i], age.dur[i])
    }
    foi = c(unlist(foi.list))
  } 
  
  
  if (length(wane_hetero)==1) wane_hetero <- rep(wane_hetero,nage) # 
  if (length(wane_hetero)==s) wane_hetero <- wane_hetero
  if (length(wane_hetero) > 1 & length(wane_hetero)<s){
    wane_hetero.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(wane_hetero)){
      wane_hetero.list[[i]] = rep(wane_hetero[i], age.dur[i])
    }
    wane_hetero = c(unlist(wane_hetero.list))
  }   
  
  #then, now that foi is appropriate length of all the age classes, 
  #multiply by the age.mult vector
  foi <- foi*age.mult.foi
  
  if (length(sigma)==1) sigma <- rep(sigma,nage) # 
  if (length(sigma)==s) sigma <- sigma
  if (length(sigma) > 1 & length(sigma)<s){
    sigma.list <- list()
    
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(sigma)){
      sigma.list[[i]] = rep(sigma[i], age.dur[i])
    }
    sigma = c(unlist(sigma.list))
  } 
  
  
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  
  
  
  mat1 <- matrix(0,c,c) 
  
  Tmat <- matrix(0,c*nage,c*nage) 
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek
    mat1[] <- 0
    
    #infections in naive susceptibles
    mat1[1,1] <- 1- (4*foi[j])
    mat1[2,1] <- foi[j]
    mat1[3,1] <- foi[j]
    mat1[4,1] <- foi[j]
    mat1[5,1] <- foi[j]
    
    #I1s recover
    mat1[2,2] <- 1-recov
    mat1[6,2] <- recov
    
    
    #I2s recover
    mat1[3,3] <- 1-recov
    mat1[7,3] <- recov
    
    #I3s recover
    mat1[4,4] <- 1-recov
    mat1[8,4] <- recov
    
    #I4s recover
    mat1[5,5] <- 1-recov
    mat1[9,5] <- recov
    
    
    #waning out of heterotypic immunity (P1 to S1) - those previously infected with I1
    mat1[6,6] <- 1- wane_hetero[j]
    mat1[10,6] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I2
    mat1[7,7] <- 1- wane_hetero[j]
    mat1[11,7] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I3
    mat1[8,8] <- 1- wane_hetero[j]
    mat1[12,8] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I4
    mat1[9,9] <- 1- wane_hetero[j]
    mat1[13,9] <- wane_hetero[j]
    
    #secondary infections - those previously infected with I1
    mat1[10,10] <- 1-(3*foi[j])
    mat1[14,10] <- foi[j] #I12
    mat1[15,10] <- foi[j] #I13
    mat1[16,10] <- foi[j] #I14
    
    #secondary infections - those previously infected with I2
    mat1[11,11] <- 1-(3*foi[j])
    mat1[17,11] <- foi[j] #I21
    mat1[18,11] <- foi[j] #I22
    mat1[19,11] <- foi[j] #I23
    
    #secondary infections - those previously infected with I3
    mat1[12,12] <-  1-(3*foi[j])
    mat1[20,12] <- foi[j] #I31
    mat1[21,12] <- foi[j] #I32
    mat1[22,12] <- foi[j] #I34
    
    #secondary infections - those previously infected with I4
    mat1[13,13] <-  1-(3*foi[j])
    mat1[23,13] <- foi[j] #I41
    mat1[24,13] <- foi[j] #I42
    mat1[25,13] <- foi[j] #I43
    
    #recovery from secondary infections - I12
    mat1[14,14] <- 1-recov
    mat1[26,14] <- recov
    
    #recovery from secondary infections - I13
    mat1[15,15] <- 1-recov
    mat1[27,15] <- recov
    
    #recovery from secondary infections - I14
    mat1[16,16] <- 1-recov
    mat1[28,16] <- recov
    
    #recovery from secondary infections - I21
    mat1[17,17] <- 1-recov
    mat1[29,17] <- recov
    
    #recovery from secondary infections - I23
    mat1[18,18] <- 1-recov
    mat1[30,18] <- recov
    
    #recovery from secondary infections - I24
    mat1[19,19] <- 1-recov
    mat1[31,19] <- recov
    
    #recovery from secondary infections - I31
    mat1[20,20] <- 1-recov
    mat1[32,30] <- recov
    
    #recovery from secondary infections - I32
    mat1[21,21] <- 1-recov
    mat1[33,21] <- recov
    
    #recovery from secondary infections - I34 to P34
    mat1[22,22] <- 1-recov
    mat1[34,22] <- recov
    
    #recovery from secondary infections - I41 to P41
    mat1[23,23] <- 1-recov
    mat1[35,23] <- recov
    
    #recovery from secondary infections - I42 to P42
    mat1[24,24] <- 1-recov
    mat1[36,24] <- recov
    
    #recovery from secondary infections - I43 to P43
    mat1[25,25] <- 1-recov
    mat1[37,25] <- recov
    
    #waning out of heterotypic immunity - those previously infected with I12 (P12 to S12)
    mat1[26,26] <- 1- wane_hetero[j]
    mat1[38,26] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I13
    mat1[27,27] <- 1- wane_hetero[j]
    mat1[39,27] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I14 (P14 to S14)
    mat1[28,28] <- 1- wane_hetero[j]
    mat1[40,28] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I21
    mat1[29,29] <- 1- wane_hetero[j]
    mat1[41,29] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I23
    mat1[30,30] <- 1- wane_hetero[j]
    mat1[42,30] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I24
    mat1[31,31] <- 1- wane_hetero[j]
    mat1[43,31] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I31
    mat1[32,32] <- 1- wane_hetero[j]
    mat1[44,32] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I32
    mat1[33,33] <- 1- wane_hetero[j]
    mat1[45,33] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I34
    mat1[34,34] <- 1- wane_hetero[j]
    mat1[46,34] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I41
    mat1[35,35] <- 1- wane_hetero[j]
    mat1[47,35] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I42
    mat1[36,36] <- 1- wane_hetero[j]
    mat1[48,36] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I43
    mat1[37,37] <- 1- wane_hetero[j]
    mat1[49,37] <- wane_hetero[j]
    
    
    #tertiary infections - those previously infected with I12 (S12 to I123 and I124)
    mat1[38,38] <- 1-2*foi[j]
    mat1[50,38] <- foi[j]
    mat1[51,38] <- foi[j]
    
    #tertiary infections - those previously infected with I13 (I132 and I34)
    mat1[39,39] <- 1-2*foi[j]
    mat1[52,39] <- foi[j]
    mat1[53,39] <- foi[j]
    
    
    #tertiary infections - those previously infected with I14 (I142 and I43)
    mat1[40,40] <- 1-2*foi[j]
    mat1[54,40] <- foi[j]
    mat1[55,40] <- foi[j]
    
    #tertiary infections - those previously infected with I21 (I213 + I214)
    mat1[41,41] <- 1-2*foi[j]
    mat1[56,41] <- foi[j]
    mat1[57,41] <- foi[j]
    
    #tertiary infections - those previously infected with I23 (I231 + I234)
    mat1[42,42] <- 1-foi[j]
    mat1[58,42] <- foi[j]
    mat1[59,42] <- foi[j]
    
    
    #tertiary infections - those previously infected with I24 (I241 + I243)
    mat1[43,43] <- 1-2*foi[j]
    mat1[60,43] <- foi[j]
    mat1[61,43] <- foi[j]
    
    
    #tertiary infections - those previously infected with I31 (I312+I314)
    mat1[44,44] <- 1-2*foi[j]
    mat1[62,44] <- foi[j]
    mat1[63,44] <- foi[j]
    
    #tertiary infections - those previously infected with I32 (I321+I324)
    mat1[45,45] <- 1-2*foi[j]
    mat1[64,45] <- foi[j]
    mat1[65,45] <- foi[j]
    
    #tertiary infections - those previously infected with I34 (I341+I342)
    mat1[46,46] <- 1-2*foi[j]
    mat1[66,46] <- foi[j]
    mat1[67,46] <- foi[j]
    
    
    #tertiary infections - those previously infected with I41 (I412+I413)
    mat1[47,47] <- 1-2*foi[j]
    mat1[68,47] <- foi[j]
    mat1[69,47] <- foi[j]
    
    #tertiary infections - those previously infected with I42 (I421+I423)
    mat1[48,48] <- 1-2*foi[j]
    mat1[70,48] <- foi[j]
    mat1[71,48] <- foi[j]
    
    #tertiary infections - those previously infected with I43 (I431+I432)
    mat1[49,49] <- 1-2*foi[j]
    mat1[72,49] <- foi[j]
    mat1[73,49] <- foi[j]
    
    
    #recovery from tertiary infections - I123
    mat1[50,50] <- 1-recov
    mat1[74,50] <- recov
    
    #recovery from tertiary infections - I124
    mat1[51,51] <- 1-recov
    mat1[75,51] <- recov
    
    #recovery from tertiary infections - I132
    mat1[52,52] <- 1-recov
    mat1[76,52] <- recov
    
    #recovery from tertiary infections - I134
    mat1[53,53] <- 1-recov
    mat1[77,53] <- recov
    
    #recovery from tertiary infections - I142
    mat1[54,54] <- 1-recov
    mat1[78,54] <- recov
    
    #recovery from tertiary infections - I143
    mat1[55,55] <- 1-recov
    mat1[79,55] <- recov
    
    
    #recovery from tertiary infections - I213
    mat1[56,56] <- 1-recov
    mat1[80,56] <- recov
    
    #recovery from tertiary infections - I214
    mat1[57,57] <- 1-recov
    mat1[81,57] <- recov
    
    #recovery from tertiary infections - I231
    mat1[58,58] <- 1-recov
    mat1[82,58] <- recov
    
    #recovery from tertiary infections - I234
    mat1[59,59] <- 1-recov
    mat1[83,59] <- recov
    
    
    #recovery from tertiary infections - I241
    mat1[60,60] <- 1-recov
    mat1[84,60] <- recov
    
    #recovery from tertiary infections - I243
    mat1[61,61] <- 1-recov
    mat1[85,61] <- recov
    
    
    #recovery from tertiary infections - I312
    mat1[62,62] <- 1-recov
    mat1[86,62] <- recov
    
    #recovery from tertiary infections - I314
    mat1[63,63] <- 1-recov
    mat1[87,63] <- recov
    
    #recovery from tertiary infections - I321
    mat1[64,64] <- 1-recov
    mat1[88,64] <- recov
    
    #recovery from tertiary infections - I324
    mat1[65,65] <- 1-recov
    mat1[89,65] <- recov
    
    #recovery from tertiary infections - I341
    mat1[66,66] <- 1-recov
    mat1[90,66] <- recov
    
    #recovery from tertiary infections - I342
    mat1[67,67] <- 1-recov
    mat1[91,67] <- recov
    
    
    
    
    #recovery from tertiary infections - I412
    mat1[68,68] <- 1-recov
    mat1[92,68] <- recov
    
    #recovery from tertiary infections - I413
    mat1[69,69] <- 1-recov
    mat1[93,69] <- recov
    
    #recovery from tertiary infections - I421
    mat1[70,70] <- 1-recov
    mat1[94,70] <- recov
    
    #recovery from tertiary infections - I423
    mat1[71,71] <- 1-recov
    mat1[95,71] <- recov
    
    #recovery from tertiary infections - I431
    mat1[72,72] <- 1-recov
    mat1[96,72] <- recov
    
    #recovery from tertiary infections - I432
    mat1[73,73] <- 1-recov
    mat1[97,73] <- recov
    
    
    #waning out of heterotypic immunity - P123 to S123
    mat1[74,74] <- 1- wane_hetero[j]
    mat1[98,74] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P124 to S124
    mat1[75,75] <- 1- wane_hetero[j]
    mat1[99,75] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P132 to S132
    mat1[76,76] <- 1- wane_hetero[j]
    mat1[100,76] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P134 to S134
    mat1[77,77] <- 1- wane_hetero[j]
    mat1[101,77] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P142
    mat1[78,78] <- 1- wane_hetero[j]
    mat1[102,78] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P143
    mat1[79,79] <- 1- wane_hetero[j]
    mat1[103,79] <- wane_hetero[j]
    
    
    
    #waning out of heterotypic immunity - P213
    mat1[80,80] <- 1- wane_hetero[j]
    mat1[104,80] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P214
    mat1[81,81] <- 1- wane_hetero[j]
    mat1[105,81] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P231
    mat1[82,82] <- 1- wane_hetero[j]
    mat1[106,82] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P234
    mat1[83,83] <- 1- wane_hetero[j]
    mat1[107,83] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P241
    mat1[84,84] <- 1- wane_hetero[j]
    mat1[108,84] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P243
    mat1[85,85] <- 1- wane_hetero[j]
    mat1[109,85] <- wane_hetero[j]
    
    
    
    #waning out of heterotypic immunity - P312
    mat1[86,86] <- 1- wane_hetero[j]
    mat1[110,86] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P314
    mat1[87,87] <- 1- wane_hetero[j]
    mat1[111,87] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P321
    mat1[88,88] <- 1- wane_hetero[j]
    mat1[112,88] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P324
    mat1[89,89] <- 1- wane_hetero[j]
    mat1[113,89] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P341
    mat1[90,90] <- 1- wane_hetero[j]
    mat1[114,90] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P342
    mat1[91,91] <- 1- wane_hetero[j]
    mat1[115,91] <- wane_hetero[j]
    
    
    
    
    #waning out of heterotypic immunity - P412
    mat1[92,92] <- 1- wane_hetero[j]
    mat1[116,92] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P413
    mat1[93,93] <- 1- wane_hetero[j]
    mat1[117,93] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P421
    mat1[94,94] <- 1- wane_hetero[j]
    mat1[118,94] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P423
    mat1[95,95] <- 1- wane_hetero[j]
    mat1[119,95] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P431
    mat1[96,96] <- 1- wane_hetero[j]
    mat1[120,96] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P432
    mat1[97,97] <- 1- wane_hetero[j]
    mat1[121,97] <- wane_hetero[j]
    
    
    
    #quaternary infections -  S123 to I1234
    mat1[98, 98] <- 1-foi[j]
    mat1[122,98] <- foi[j]
    
    #quaternary infections -  S124 to I1243
    mat1[99, 99] <- 1-foi[j]
    mat1[123,99] <- foi[j]
    
    #quaternary infections -  S132 to I1324
    mat1[100, 100] <- 1-foi[j]
    mat1[124,100] <- foi[j]
    
    #quaternary infections -  S134 to I1342
    mat1[101, 101] <- 1-foi[j]
    mat1[125,101] <- foi[j]
    
    #quaternary infections -  S142
    mat1[102, 102] <- 1-foi[j]
    mat1[126,102] <- foi[j]
    
    #quaternary infections -  S143
    mat1[103, 103] <- 1-foi[j]
    mat1[127,103] <- foi[j]
    
    
    
    
    #quaternary infections -  S213
    mat1[104, 104] <- 1-foi[j]
    mat1[128,104] <- foi[j]
    
    #quaternary infections -  S214
    mat1[105, 105] <- 1-foi[j]
    mat1[129,105] <- foi[j]
    
    #quaternary infections -  S231
    mat1[106, 106] <- 1-foi[j]
    mat1[130,106] <- foi[j]
    
    #quaternary infections -  S234
    mat1[107, 107] <- 1-foi[j]
    mat1[131,107] <- foi[j]
    
    #quaternary infections -  S241
    mat1[108, 108] <- 1-foi[j]
    mat1[132,108] <- foi[j]
    
    #quaternary infections -  S243
    mat1[109, 109] <- 1-foi[j]
    mat1[133,109] <- foi[j]
    
    
    
    
    #quaternary infections -  S312
    mat1[110, 110] <- 1-foi[j]
    mat1[134,110] <- foi[j]
    
    #quaternary infections -  S314
    mat1[111, 111] <- 1-foi[j]
    mat1[135,111] <- foi[j]
    
    #quaternary infections -  S321
    mat1[112, 112] <- 1-foi[j]
    mat1[136,112] <- foi[j]
    
    #quaternary infections -  S324
    mat1[113, 113] <- 1-foi[j]
    mat1[137, 113] <- foi[j]
    
    #quaternary infections -  S341
    mat1[114, 114] <- 1-foi[j]
    mat1[138, 114] <- foi[j]
    
    #quaternary infections -  S342
    mat1[115, 115] <- 1-foi[j]
    mat1[139,115] <- foi[j]
    
    
    
    
    #quaternary infections -  S412
    mat1[116, 116] <- 1-foi[j]
    mat1[140,116] <- foi[j]
    
    #quaternary infections -  S413
    mat1[117, 117] <- 1-foi[j]
    mat1[141,117] <- foi[j]
    
    #quaternary infections -  S421
    mat1[118, 118] <- 1-foi[j]
    mat1[142,118] <- foi[j]
    
    #quaternary infections -  S423
    mat1[119, 119] <- 1-foi[j]
    mat1[143, 119] <- foi[j]
    
    #quaternary infections -  S431
    mat1[120, 120] <- 1-foi[j]
    mat1[144, 120] <- foi[j]
    
    #quaternary infections -  S432
    mat1[121, 121] <- 1-foi[j]
    mat1[145, 121] <- foi[j]
    
    
    #recovery from quaternary infections into short or longterm heterotypic immunity to all serotypes
    #I1234 to Pms
    mat1[122,122] <- 1-recov
    mat1[146,122] <- recov
    
    mat1[123,123] <- 1-recov
    mat1[146,123] <- recov
    
    mat1[124,124] <- 1-recov
    mat1[146,124] <- recov
    
    mat1[125,125] <- 1-recov
    mat1[146,125] <- recov
    
    mat1[126,126] <- 1-recov
    mat1[146,126] <- recov
    
    mat1[127,127] <- 1-recov
    mat1[146,127] <- recov
    
    mat1[128,128] <- 1-recov
    mat1[146,128] <- recov
    
    mat1[129,129] <- 1-recov
    mat1[146,129] <- recov
    
    mat1[130,130] <- 1-recov
    mat1[146,130] <- recov
    
    mat1[131,131] <- 1-recov
    mat1[146,131] <- recov
    
    mat1[132,132] <- 1-recov
    mat1[146,132] <- recov
    
    mat1[133,133] <- 1-recov
    mat1[146,133] <- recov
    
    mat1[134,134] <- 1-recov
    mat1[146,134] <- recov
    
    mat1[135,135] <- 1-recov
    mat1[146,135] <- recov
    
    mat1[136,136] <- 1-recov
    mat1[146,136] <- recov
    
    mat1[137,137] <- 1-recov
    mat1[146,137] <- recov
    
    mat1[138,138] <- 1-recov
    mat1[146,138] <- recov
    
    mat1[139,139] <- 1-recov
    mat1[146,139] <- recov
    
    mat1[140,140] <- 1-recov
    mat1[146,140] <- recov
    
    mat1[141,141] <- 1-recov
    mat1[146,141] <- recov
    
    mat1[142,142] <- 1-recov
    mat1[146,142] <- recov
    
    mat1[143,143] <- 1-recov
    mat1[146,143] <- recov
    
    mat1[144,144] <- 1-recov
    mat1[146,144] <- recov
    
    mat1[145,145] <- 1-recov
    mat1[146,145] <- recov
    
    #then waning into longterm multivalent immunity
    mat1[146,146] <- 1-wane_hetero[j]
    mat1[147,146] <- wane_hetero[j]
    
    #once you reach immunity, you stay
    mat1[147,147] <- 1
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*c+1):((j-1)*c+c),((j-1)*c+1):((j-1)*c+c)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
  return(Tmat)
  
}
buildTMat_five_strain_intro <- function(c, Npop, age.classes, surv.biwk, wane_hetero, age.brk, foi, age.mult.foi, recov, sigma, age.rate){
  
  s <- nage <- length(age.classes)
  
  
  
  if (length(foi)==1) foi <- rep(foi,nage) # 
  if (length(foi)==s) foi <- foi
  if (length(foi) > 1 & length(foi)<s){
    foi.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(foi)){
      foi.list[[i]] = rep(foi[i], age.dur[i])
    }
    foi = c(unlist(foi.list))
  } 
  
  
  if (length(wane_hetero)==1) wane_hetero <- rep(wane_hetero,nage) # 
  if (length(wane_hetero)==s) wane_hetero <- wane_hetero
  if (length(wane_hetero) > 1 & length(wane_hetero)<s){
    wane_hetero.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(wane_hetero)){
      wane_hetero.list[[i]] = rep(wane_hetero[i], age.dur[i])
    }
    wane_hetero = c(unlist(wane_hetero.list))
  }   
  
  #then, now that foi is appropriate length of all the age classes, 
  #multiply by the age.mult vector
  foi <- foi*age.mult.foi
  
  if (length(sigma)==1) sigma <- rep(sigma,nage) # 
  if (length(sigma)==s) sigma <- sigma
  if (length(sigma) > 1 & length(sigma)<s){
    sigma.list <- list()
    
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(sigma)){
      sigma.list[[i]] = rep(sigma[i], age.dur[i])
    }
    sigma = c(unlist(sigma.list))
  } 
  
  
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  
  
  
  mat1 <- matrix(0,c,c) 
  
  Tmat <- matrix(0,c*nage,c*nage) 
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek
    mat1[] <- 0
    
    #infections in naive susceptibles
    mat1[1,1] <- 1- (5*foi[j])
    mat1[2,1] <- foi[j]
    mat1[3,1] <- foi[j]
    mat1[4,1] <- foi[j]
    mat1[5,1] <- foi[j]
    mat1[6,1] <- foi[j]
    
    #I1s recover
    mat1[2,2] <- 1-recov
    mat1[7,2] <- recov
    
    
    #I2s recover
    mat1[3,3] <- 1-recov
    mat1[8,3] <- recov
    
    #I3s recover
    mat1[4,4] <- 1-recov
    mat1[9,4] <- recov
    
    #I4s recover
    mat1[5,5] <- 1-recov
    mat1[10,5] <- recov
    
    #I5s recover
    mat1[6,6] <- 1-recov
    mat1[11,6] <- recov
    
    
    #waning out of heterotypic immunity (P1 to S1) - those previously infected with I1
    mat1[7,7] <- 1- wane_hetero[j]
    mat1[12,7] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I2
    mat1[8,8] <- 1- wane_hetero[j]
    mat1[13,8] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I3
    mat1[9,9] <- 1- wane_hetero[j]
    mat1[14,9] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I4
    mat1[10,10] <- 1- wane_hetero[j]
    mat1[15,10] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I5
    mat1[11,11] <- 1- wane_hetero[j]
    mat1[16,11] <- wane_hetero[j]
    
    #secondary infections - those previously infected with I1
    mat1[12,12] <- 1-(4*foi[j])
    mat1[17,12] <- foi[j] #I12
    mat1[18,12] <- foi[j] #I13
    mat1[19,12] <- foi[j] #I14
    mat1[20,12] <- foi[j] #I15
    
    #secondary infections - those previously infected with I2
    mat1[13,13] <- 1-(4*foi[j])
    mat1[21,13] <- foi[j] #I21
    mat1[22,13] <- foi[j] #I22
    mat1[23,13] <- foi[j] #I23
    mat1[24,13] <- foi[j] #I25
    
    #secondary infections - those previously infected with I3
    mat1[14,14] <-  1-(4*foi[j])
    mat1[25,14] <- foi[j] #I31
    mat1[26,14] <- foi[j] #I32
    mat1[27,14] <- foi[j] #I34
    mat1[28,14] <- foi[j] #I35
    
    #secondary infections - those previously infected with I4
    mat1[15,15] <-  1-(4*foi[j])
    mat1[29,15] <- foi[j] #I41
    mat1[30,15] <- foi[j] #I42
    mat1[31,15] <- foi[j] #I43
    mat1[32,15] <- foi[j] #I45
    
    
    #secondary infections - those previously infected with I5
    mat1[16,16] <-  1-(4*foi[j])
    mat1[33,16] <- foi[j] #I51
    mat1[34,16] <- foi[j] #I52
    mat1[35,16] <- foi[j] #I53
    mat1[36,16] <- foi[j] #I54
    
    #recovery from secondary infections - I12
    mat1[17,17] <- 1-recov
    mat1[37,17] <- recov
    
    #recovery from secondary infections - I13
    mat1[18,18] <- 1-recov
    mat1[38,18] <- recov
    
    #recovery from secondary infections - I14
    mat1[19,19] <- 1-recov
    mat1[39,19] <- recov
    
    #recovery from secondary infections - I15
    mat1[20,20] <- 1-recov
    mat1[40,20] <- recov
    
    #recovery from secondary infections - I21
    mat1[21,21] <- 1-recov
    mat1[41,21] <- recov
    
    #recovery from secondary infections - I23
    mat1[22,22] <- 1-recov
    mat1[42,22] <- recov
    
    #recovery from secondary infections - I24
    mat1[23,23] <- 1-recov
    mat1[43,23] <- recov
    
    #recovery from secondary infections - I25
    mat1[24,24] <- 1-recov
    mat1[44,24] <- recov
    
    #recovery from secondary infections - I31
    mat1[25,25] <- 1-recov
    mat1[45,25] <- recov
    
    #recovery from secondary infections - I32
    mat1[26,26] <- 1-recov
    mat1[46,26] <- recov
    
    #recovery from secondary infections - I34 to P34
    mat1[27,27] <- 1-recov
    mat1[47,27] <- recov
    
    #recovery from secondary infections - I35 to P35
    mat1[28,28] <- 1-recov
    mat1[48,28] <- recov
    
    #recovery from secondary infections - I41 to P41
    mat1[29,29] <- 1-recov
    mat1[49,29] <- recov
    
    #recovery from secondary infections - I42 to P42
    mat1[30,30] <- 1-recov
    mat1[50,30] <- recov
    
    #recovery from secondary infections - I43 to P43
    mat1[31,31] <- 1-recov
    mat1[51,31] <- recov
    
    #recovery from secondary infections - I45 to P45
    mat1[32,32] <- 1-recov
    mat1[52,32] <- recov
    
    
    #recovery from secondary infections - I51 to P51
    mat1[33,33] <- 1-recov
    mat1[53,33] <- recov
    
    #recovery from secondary infections - I52 to P52
    mat1[34,34] <- 1-recov
    mat1[54,34] <- recov
    
    #recovery from secondary infections - I53 to P53
    mat1[35,35] <- 1-recov
    mat1[55,35] <- recov
    
    #recovery from secondary infections - I54 to P54
    mat1[36,36] <- 1-recov
    mat1[56,36] <- recov
    
    #waning out of heterotypic immunity - those previously infected with I12 (P12 to S12)
    mat1[37,37] <- 1- wane_hetero[j]
    mat1[57,37] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I13
    mat1[38,38] <- 1- wane_hetero[j]
    mat1[58,38] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I14 (P14 to S14)
    mat1[39,39] <- 1- wane_hetero[j]
    mat1[59,39] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I15 (P15 to S15)
    mat1[40,40] <- 1- wane_hetero[j]
    mat1[60,40] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I21
    mat1[41,41] <- 1- wane_hetero[j]
    mat1[61,41] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I23
    mat1[42,42] <- 1- wane_hetero[j]
    mat1[62,42] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I24
    mat1[43,43] <- 1- wane_hetero[j]
    mat1[63,43] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I25
    mat1[44,44] <- 1- wane_hetero[j]
    mat1[64,44] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I31
    mat1[45,45] <- 1- wane_hetero[j]
    mat1[65,45] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I32
    mat1[46,46] <- 1- wane_hetero[j]
    mat1[66,46] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I34
    mat1[47,47] <- 1- wane_hetero[j]
    mat1[67,47] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I35
    mat1[48,48] <- 1- wane_hetero[j]
    mat1[68,48] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I41
    mat1[49,49] <- 1- wane_hetero[j]
    mat1[69,49] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I42
    mat1[50,50] <- 1- wane_hetero[j]
    mat1[70,50] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I43
    mat1[51,51] <- 1- wane_hetero[j]
    mat1[71,51] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I45
    mat1[52,52] <- 1- wane_hetero[j]
    mat1[72,52] <- wane_hetero[j]
    
    
    #waning out of heterotypic immunity - those previously infected with I51
    mat1[53,53] <- 1- wane_hetero[j]
    mat1[73,53] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I52
    mat1[54,54] <- 1- wane_hetero[j]
    mat1[74,54] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I53
    mat1[55,55] <- 1- wane_hetero[j]
    mat1[75,55] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I54
    mat1[56,56] <- 1- wane_hetero[j]
    mat1[76,56] <- wane_hetero[j]
    
    
    #tertiary infections - those previously infected with I12 (S12 to I123 and I124 and I125)
    mat1[57,57] <- 1-3*foi[j]
    mat1[77,57] <- foi[j]
    mat1[78,57] <- foi[j]
    mat1[79,57] <- foi[j]
    
    #tertiary infections - those previously infected with I13 (I132 and I134 and I135)
    mat1[58,58] <- 1-3*foi[j]
    mat1[80,58] <- foi[j]
    mat1[81,58] <- foi[j]
    mat1[82,58] <- foi[j]
    
    
    #tertiary infections - those previously infected with I14 (I142 and I143 and I145)
    mat1[59,59] <- 1-3*foi[j]
    mat1[83,59] <- foi[j]
    mat1[84,59] <- foi[j]
    mat1[85,59] <- foi[j]
    
    #tertiary infections - those previously infected with I15 (I152 and I153 and I154)
    mat1[60,60] <- 1-3*foi[j]
    mat1[86,60] <- foi[j]
    mat1[87,60] <- foi[j]
    mat1[88,60] <- foi[j]
    
    #tertiary infections - those previously infected with I21 (I213 + I214 + I215)
    mat1[61,61] <- 1-3*foi[j]
    mat1[89,61] <- foi[j]
    mat1[90,61] <- foi[j]
    mat1[91,61] <- foi[j]
    
    #tertiary infections - those previously infected with I23 (I231 + I234 + I235)
    mat1[62,62] <- 1-3*foi[j]
    mat1[92,62] <- foi[j]
    mat1[93,62] <- foi[j]
    mat1[94,62] <- foi[j]
    
    #tertiary infections - those previously infected with I24 (I241 + I243 + I245)
    mat1[63,63] <- 1-3*foi[j]
    mat1[95,63] <- foi[j]
    mat1[96,63] <- foi[j]
    mat1[97,63] <- foi[j]
    
    #tertiary infections - those previously infected with I25 (I251+I253+I254)
    mat1[64,64] <- 1-3*foi[j]
    mat1[98,64] <- foi[j]
    mat1[99,64] <- foi[j]
    mat1[100,64] <- foi[j]
    
    #tertiary infections - those previously infected with I31 (I312+I314+I315)
    mat1[65,65] <- 1-3*foi[j]
    mat1[101,65] <- foi[j]
    mat1[102,65] <- foi[j]
    mat1[103,65] <- foi[j]
    
    
    #tertiary infections - those previously infected with I32 (I321+I324+I325)
    mat1[66,66] <- 1-3*foi[j]
    mat1[104,66] <- foi[j]
    mat1[105,66] <- foi[j]
    mat1[106,66] <- foi[j]
    
    #tertiary infections - those previously infected with I34 (I341+I342+I345)
    mat1[67,67] <- 1-3*foi[j]
    mat1[107,67] <- foi[j]
    mat1[108,67] <- foi[j]
    mat1[109,67] <- foi[j]
    
    #tertiary infections - those previously infected with I35 (I351+I352+I354)
    mat1[68,68] <- 1-3*foi[j]
    mat1[110,68] <- foi[j]
    mat1[111,68] <- foi[j]
    mat1[112,68] <- foi[j]
    
    #tertiary infections - those previously infected with I41 (I412+I413+415)
    mat1[69,69] <- 1-3*foi[j]
    mat1[113,69] <- foi[j]
    mat1[114,69] <- foi[j]
    mat1[115,69] <- foi[j]
    
    #tertiary infections - those previously infected with I42 (I421+I423+I425)
    mat1[70,70] <- 1-3*foi[j]
    mat1[116,70] <- foi[j]
    mat1[117,70] <- foi[j]
    mat1[118,70] <- foi[j]
    
    #tertiary infections - those previously infected with I43 (I431+I432+I435)
    mat1[71,71] <- 1-3*foi[j]
    mat1[119,71] <- foi[j]
    mat1[120,71] <- foi[j]
    mat1[121,71] <- foi[j]
    
    #tertiary infections - those previously infected with I45 (I451+I452+I453)
    mat1[72,72] <- 1-3*foi[j]
    mat1[122,72] <- foi[j]
    mat1[123,72] <- foi[j]
    mat1[124,72] <- foi[j]
    
    
    #tertiary infections - those previously infected with I51 (I512+I513+514)
    mat1[73,73] <- 1-3*foi[j]
    mat1[125,73] <- foi[j]
    mat1[126,73] <- foi[j]
    mat1[127,73] <- foi[j]
    
    
    #tertiary infections - those previously infected with I52 (I521+I523+524)
    mat1[74,74] <- 1-3*foi[j]
    mat1[128,74] <- foi[j]
    mat1[129,74] <- foi[j]
    mat1[130,74] <- foi[j]
    
    #tertiary infections - those previously infected with I53 (I531+I532+534)
    mat1[75,75] <- 1-3*foi[j]
    mat1[131,75] <- foi[j]
    mat1[132,75] <- foi[j]
    mat1[133,75] <- foi[j]
    
    
    #tertiary infections - those previously infected with I54 (I541+I542+543)
    mat1[76,76] <- 1-3*foi[j]
    mat1[134,76] <- foi[j]
    mat1[135,76] <- foi[j]
    mat1[136,76] <- foi[j]
    
    
    #recovery from tertiary infections - I123
    mat1[77,77] <- 1-recov
    mat1[137,77] <- recov
    
    #recovery from tertiary infections - I124
    mat1[78,78] <- 1-recov
    mat1[138,78] <- recov
    
    #recovery from tertiary infections - I125
    mat1[79,79] <- 1-recov
    mat1[139,79] <- recov
    
    #recovery from tertiary infections - I132
    mat1[80,80] <- 1-recov
    mat1[140,80] <- recov
    
    #recovery from tertiary infections - I134
    mat1[81,81] <- 1-recov
    mat1[141,81] <- recov
    
    #recovery from tertiary infections - I135
    mat1[82,82] <- 1-recov
    mat1[142,82] <- recov
    
    #recovery from tertiary infections - I142
    mat1[83,83] <- 1-recov
    mat1[143,83] <- recov
    
    #recovery from tertiary infections - I143
    mat1[84,84] <- 1-recov
    mat1[144,84] <- recov
    
    #recovery from tertiary infections - I145
    mat1[85,85] <- 1-recov
    mat1[145,85] <- recov
    
    
    #recovery from tertiary infections - I152
    mat1[86,86] <- 1-recov
    mat1[146,86] <- recov
    
    #recovery from tertiary infections - I153
    mat1[87,87] <- 1-recov
    mat1[147,87] <- recov
    
    #recovery from tertiary infections - I154
    mat1[88,88] <- 1-recov
    mat1[148,88] <- recov
    
    #recovery from tertiary infections - I213
    mat1[89,89] <- 1-recov
    mat1[149,89] <- recov
    
    #recovery from tertiary infections - I214
    mat1[90,90] <- 1-recov
    mat1[150,90] <- recov
    
    #recovery from tertiary infections - I215
    mat1[91,91] <- 1-recov
    mat1[151,91] <- recov
    
    
    #recovery from tertiary infections - I231
    mat1[92,92] <- 1-recov
    mat1[152,92] <- recov
    
    #recovery from tertiary infections - I234
    mat1[93,93] <- 1-recov
    mat1[153,93] <- recov
    
    #recovery from tertiary infections - I235
    mat1[94,94] <- 1-recov
    mat1[154,94] <- recov
    
    #recovery from tertiary infections - I241
    mat1[95,95] <- 1-recov
    mat1[155,95] <- recov
    
    #recovery from tertiary infections - I243
    mat1[96,96] <- 1-recov
    mat1[156,96] <- recov
    
    #recovery from tertiary infections - I245
    mat1[97,97] <- 1-recov
    mat1[157,97] <- recov
    
    
    #recovery from tertiary infections - I251
    mat1[98,98] <- 1-recov
    mat1[158,98] <- recov
    
    #recovery from tertiary infections - I253
    mat1[99,99] <- 1-recov
    mat1[159,99] <- recov
    
    #recovery from tertiary infections - I254
    mat1[100,100] <- 1-recov
    mat1[160,100] <- recov
    
    #recovery from tertiary infections - I312
    mat1[101,101] <- 1-recov
    mat1[161,101] <- recov
    
    #recovery from tertiary infections - I314
    mat1[102,102] <- 1-recov
    mat1[162,102] <- recov
    
    #recovery from tertiary infections - I315
    mat1[103,103] <- 1-recov
    mat1[163,103] <- recov
    
    #recovery from tertiary infections - I321
    mat1[104,104] <- 1-recov
    mat1[164,104] <- recov
    
    #recovery from tertiary infections - I324
    mat1[105,105] <- 1-recov
    mat1[165,105] <- recov
    
    #recovery from tertiary infections - I325
    mat1[106,106] <- 1-recov
    mat1[166,106] <- recov
    
    #recovery from tertiary infections - I341
    mat1[107,107] <- 1-recov
    mat1[167,107] <- recov
    
    #recovery from tertiary infections - I342
    mat1[108,108] <- 1-recov
    mat1[168,108] <- recov
    
    #recovery from tertiary infections - I345
    mat1[109,109] <- 1-recov
    mat1[169,109] <- recov
    
    
    #recovery from tertiary infections - I351
    mat1[110,110] <- 1-recov
    mat1[170,110] <- recov
    
    #recovery from tertiary infections - I352
    mat1[111,111] <- 1-recov
    mat1[171,111] <- recov
    
    #recovery from tertiary infections - I354
    mat1[112,112] <- 1-recov
    mat1[172,112] <- recov
    
    
    #recovery from tertiary infections - I412
    mat1[113,113] <- 1-recov
    mat1[173,113] <- recov
    
    #recovery from tertiary infections - I413
    mat1[114,114] <- 1-recov
    mat1[174,114] <- recov
    
    #recovery from tertiary infections - I415
    mat1[115,115] <- 1-recov
    mat1[175,115] <- recov
    
    #recovery from tertiary infections - I421
    mat1[116,116] <- 1-recov
    mat1[176,116] <- recov
    
    #recovery from tertiary infections - I423
    mat1[117,117] <- 1-recov
    mat1[177,117] <- recov
    
    #recovery from tertiary infections - I425
    mat1[118,118] <- 1-recov
    mat1[178,118] <- recov
    
    #recovery from tertiary infections - I431
    mat1[119,119] <- 1-recov
    mat1[179,119] <- recov
    
    #recovery from tertiary infections - I432
    mat1[120,120] <- 1-recov
    mat1[180,120] <- recov
    
    #recovery from tertiary infections - I435
    mat1[121,121] <- 1-recov
    mat1[181,121] <- recov
    
    
    #recovery from tertiary infections - I451
    mat1[122,122] <- 1-recov
    mat1[182,122] <- recov
    
    #recovery from tertiary infections - I452
    mat1[123,123] <- 1-recov
    mat1[183,123] <- recov
    
    #recovery from tertiary infections - I453
    mat1[124,124] <- 1-recov
    mat1[184,124] <- recov
    
    #recovery from tertiary infections - I512
    mat1[125,125] <- 1-recov
    mat1[185,125] <- recov
    
    #recovery from tertiary infections - I513
    mat1[126,126] <- 1-recov
    mat1[186,126] <- recov
    
    #recovery from tertiary infections - I514
    mat1[127,127] <- 1-recov
    mat1[187,127] <- recov
    
    #recovery from tertiary infections - I521
    mat1[128,128] <- 1-recov
    mat1[188,128] <- recov
    
    #recovery from tertiary infections - I523
    mat1[129,129] <- 1-recov
    mat1[189,129] <- recov
    
    #recovery from tertiary infections - I524
    mat1[130,130] <- 1-recov
    mat1[190,130] <- recov
    
    #recovery from tertiary infections - I531
    mat1[131,131] <- 1-recov
    mat1[191,131] <- recov
    
    #recovery from tertiary infections - I532
    mat1[132,132] <- 1-recov
    mat1[192,132] <- recov
    
    #recovery from tertiary infections - I534
    mat1[133,133] <- 1-recov
    mat1[193,133] <- recov
    
    #recovery from tertiary infections - I541
    mat1[134,134] <- 1-recov
    mat1[194,134] <- recov
    
    #recovery from tertiary infections - I542
    mat1[135,135] <- 1-recov
    mat1[195,135] <- recov
    
    #recovery from tertiary infections - I543
    mat1[136,136] <- 1-recov
    mat1[196,136] <- recov
    
    
    #waning out of heterotypic immunity - P123 to S123
    mat1[137,137] <- 1- wane_hetero[j]
    mat1[197,137] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P124 to S124
    mat1[138,138] <- 1- wane_hetero[j]
    mat1[198,138] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P125 to S125
    mat1[139,139] <- 1- wane_hetero[j]
    mat1[199,139] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P132 to S132
    mat1[140,140] <- 1- wane_hetero[j]
    mat1[200,140] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P134 to S134
    mat1[141,141] <- 1- wane_hetero[j]
    mat1[201,141] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P135 to S135
    mat1[142,142] <- 1- wane_hetero[j]
    mat1[202,142] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P142
    mat1[143,143] <- 1- wane_hetero[j]
    mat1[203,143] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P143
    mat1[144,144] <- 1- wane_hetero[j]
    mat1[204,144] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P145
    mat1[145,145] <- 1- wane_hetero[j]
    mat1[205,145] <- wane_hetero[j]
    
    
    #waning out of heterotypic immunity - P152
    mat1[146,146] <- 1- wane_hetero[j]
    mat1[206,146] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P153
    mat1[147,147] <- 1- wane_hetero[j]
    mat1[207,147] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P154
    mat1[148,148] <- 1- wane_hetero[j]
    mat1[208,148] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P213
    mat1[149,149] <- 1- wane_hetero[j]
    mat1[209,149] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P214
    mat1[150,150] <- 1- wane_hetero[j]
    mat1[210,150] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P215
    mat1[151,151] <- 1- wane_hetero[j]
    mat1[211,151] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P231
    mat1[152,152] <- 1- wane_hetero[j]
    mat1[212,152] <- wane_hetero[j]
  
    
    #waning out of heterotypic immunity - P234
    mat1[153,153] <- 1- wane_hetero[j]
    mat1[213,153] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P235
    mat1[154,154] <- 1- wane_hetero[j]
    mat1[214,154] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P241
    mat1[155,155] <- 1- wane_hetero[j]
    mat1[215,155] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P243
    mat1[156,156] <- 1- wane_hetero[j]
    mat1[216,156] <- wane_hetero[j]

    #waning out of heterotypic immunity - P245
    mat1[157,157] <- 1- wane_hetero[j]
    mat1[217,157] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P251
    mat1[158,158] <- 1- wane_hetero[j]
    mat1[218,158] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P253
    mat1[159,159] <- 1- wane_hetero[j]
    mat1[219,159] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P254
    mat1[160,160] <- 1- wane_hetero[j]
    mat1[220,160] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P312
    mat1[161,161] <- 1- wane_hetero[j]
    mat1[221,161] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P314
    mat1[162,162] <- 1- wane_hetero[j]
    mat1[222,162] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P315
    mat1[163,163] <- 1- wane_hetero[j]
    mat1[223,163] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P321
    mat1[164,164] <- 1- wane_hetero[j]
    mat1[224,164] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P324
    mat1[165,165] <- 1- wane_hetero[j]
    mat1[225,165] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P325
    mat1[166,166] <- 1- wane_hetero[j]
    mat1[226,166] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P341
    mat1[167,167] <- 1- wane_hetero[j]
    mat1[227,167] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P342
    mat1[168,168] <- 1- wane_hetero[j]
    mat1[228,168] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P345
    mat1[169,169] <- 1- wane_hetero[j]
    mat1[229,169] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P351
    mat1[170,170] <- 1- wane_hetero[j]
    mat1[230,170] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P352
    mat1[171,171] <- 1- wane_hetero[j]
    mat1[231,171] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P354
    mat1[172,172] <- 1- wane_hetero[j]
    mat1[232,172] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P412
    mat1[173,173] <- 1- wane_hetero[j]
    mat1[233,173] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P413
    mat1[174,174] <- 1- wane_hetero[j]
    mat1[234,174] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P415
    mat1[175,175] <- 1- wane_hetero[j]
    mat1[235,175] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P421
    mat1[176,176] <- 1- wane_hetero[j]
    mat1[236,176] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P423
    mat1[177,177] <- 1- wane_hetero[j]
    mat1[237,177] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P425
    mat1[178,178] <- 1- wane_hetero[j]
    mat1[238,178] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P431
    mat1[179,179] <- 1- wane_hetero[j]
    mat1[239,179] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P432
    mat1[180,180] <- 1- wane_hetero[j]
    mat1[240,180] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P435
    mat1[181,181] <- 1- wane_hetero[j]
    mat1[241,181] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P451
    mat1[182,182] <- 1- wane_hetero[j]
    mat1[242,182] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P452
    mat1[183,183] <- 1- wane_hetero[j]
    mat1[243,183] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P453
    mat1[184,184] <- 1- wane_hetero[j]
    mat1[244,184] <- wane_hetero[j]
    
    
      
    #waning out of heterotypic immunity - P512
    mat1[185,185] <- 1- wane_hetero[j]
    mat1[245,185] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P513
    mat1[186,186] <- 1- wane_hetero[j]
    mat1[246,186] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P514
    mat1[187,187] <- 1- wane_hetero[j]
    mat1[247,187] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P521
    mat1[188,188] <- 1- wane_hetero[j]
    mat1[248,188] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P523
    mat1[189,189] <- 1- wane_hetero[j]
    mat1[249,189] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P524
    mat1[190,190] <- 1- wane_hetero[j]
    mat1[250,190] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P531
    mat1[191,191] <- 1- wane_hetero[j]
    mat1[251,191] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P532
    mat1[192,192] <- 1- wane_hetero[j]
    mat1[252,192] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P534
    mat1[193,193] <- 1- wane_hetero[j]
    mat1[253,193] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P541
    mat1[194,194] <- 1- wane_hetero[j]
    mat1[254,194] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P542
    mat1[195,195] <- 1- wane_hetero[j]
    mat1[255,195] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P543
    mat1[196,196] <- 1- wane_hetero[j]
    mat1[256,196] <- wane_hetero[j]
    
    
    #quaternary infections -  S123 to I1234 or I1235
    mat1[197, 197] <- 1-2*foi[j]
    mat1[257,197] <- foi[j]
    mat1[258,197] <- foi[j]
    
    #quaternary infections -  S124 to I1243 or I1245
    mat1[198, 198] <- 1-2*foi[j]
    mat1[259,198] <- foi[j]
    mat1[260,198] <- foi[j]
    
    #quaternary infections -  S125
    mat1[199, 199] <- 1-2*foi[j]
    mat1[261,199] <- foi[j]
    mat1[262,199] <- foi[j]
    
    #quaternary infections -  S132 to I1324 or I1325
    mat1[200, 200] <- 1-2*foi[j]
    mat1[263,200] <- foi[j]
    mat1[264,200] <- foi[j]
    
    #quaternary infections -  S134 to I1342 or I1345
    mat1[201, 201] <- 1-2*foi[j]
    mat1[265,201] <- foi[j]
    mat1[266,201] <- foi[j]
    
    #quaternary infections -  S135
    mat1[202, 202] <- 1-2*foi[j]
    mat1[267,202] <- foi[j]
    mat1[268,202] <- foi[j]
    
    #quaternary infections -  S142
    mat1[203, 203] <- 1-2*foi[j]
    mat1[269,203] <- foi[j]
    mat1[270,203] <- foi[j]
    
    #quaternary infections -  S143
    mat1[204, 204] <- 1-2*foi[j]
    mat1[271,204] <- foi[j]
    mat1[272,204] <- foi[j]
    
    #quaternary infections -  S145
    mat1[205, 205] <- 1-2*foi[j]
    mat1[273,205] <- foi[j]
    mat1[274,205] <- foi[j]
    
    #quaternary infections -  S152
    mat1[206, 206] <- 1-2*foi[j]
    mat1[275,206] <- foi[j]
    mat1[276,206] <- foi[j]
    
    #quaternary infections -  S153
    mat1[207, 207] <- 1-2*foi[j]
    mat1[277,207] <- foi[j]
    mat1[278,207] <- foi[j]
    
    #quaternary infections -  S154
    mat1[208, 208] <- 1-2*foi[j]
    mat1[279,208] <- foi[j]
    mat1[280,208] <- foi[j]
    
    #quaternary infections -  S213
    mat1[209, 209] <- 1-2*foi[j]
    mat1[281,209] <- foi[j]
    mat1[282,209] <- foi[j]
    
    #quaternary infections -  S214
    mat1[210, 210] <- 1-2*foi[j]
    mat1[283,210] <- foi[j]
    mat1[284,210] <- foi[j]
    
    #quaternary infections -  S215
    mat1[211, 211] <- 1-2*foi[j]
    mat1[285,211] <- foi[j]
    mat1[286,211] <- foi[j]
    
    #quaternary infections -  S231
    mat1[212, 212] <- 1-2*foi[j]
    mat1[287,212] <- foi[j]
    mat1[288,212] <- foi[j]
    
    #quaternary infections -  S234
    mat1[213, 213] <- 1-2*foi[j]
    mat1[289,213] <- foi[j]
    mat1[290,213] <- foi[j]
    
    #quaternary infections -  S235
    mat1[214, 214] <- 1-2*foi[j]
    mat1[291,214] <- foi[j]
    mat1[292,214] <- foi[j]
    
    #quaternary infections -  S241
    mat1[215, 215] <- 1-2*foi[j]
    mat1[293,215] <- foi[j]
    mat1[294,215] <- foi[j]
    
    #quaternary infections -  S243
    mat1[216, 216] <- 1-2*foi[j]
    mat1[295,216] <- foi[j]
    mat1[296,216] <- foi[j]
    
    #quaternary infections -  S245
    mat1[217, 217] <- 1-2*foi[j]
    mat1[297,217] <- foi[j]
    mat1[298,217] <- foi[j]
    
    
    #quaternary infections -  S251
    mat1[218, 218] <- 1-2*foi[j]
    mat1[299,218] <- foi[j]
    mat1[300,218] <- foi[j]
    
    
    #quaternary infections -  S253
    mat1[219, 219] <- 1-2*foi[j]
    mat1[301,219] <- foi[j]
    mat1[302,219] <- foi[j]
    
    #quaternary infections -  S254
    mat1[220, 220] <- 1-2*foi[j]
    mat1[303,220] <- foi[j]
    mat1[304,220] <- foi[j]
    
    #quaternary infections -  S312
    mat1[221, 221] <- 1-2*foi[j]
    mat1[305,221] <- foi[j]
    mat1[306,221] <- foi[j]
    
    #quaternary infections -  S314
    mat1[222, 222] <- 1-2*foi[j]
    mat1[307,222] <- foi[j]
    mat1[308,222] <- foi[j]
    
    #quaternary infections -  S315
    mat1[223, 223] <- 1-2*foi[j]
    mat1[309,223] <- foi[j]
    mat1[310,223] <- foi[j]
    
    #quaternary infections -  S321
    mat1[224, 224] <- 1-2*foi[j]
    mat1[311,224] <- foi[j]
    mat1[312,224] <- foi[j]
    
    #quaternary infections -  S324
    mat1[225, 225] <- 1-2*foi[j]
    mat1[313,225] <- foi[j]
    mat1[314,225] <- foi[j]
    
    #quaternary infections -  S325
    mat1[226, 226] <- 1-2*foi[j]
    mat1[315,226] <- foi[j]
    mat1[316,226] <- foi[j]
    
    #quaternary infections -  S341
    mat1[227, 227] <- 1-2*foi[j]
    mat1[317,227] <- foi[j]
    mat1[318,227] <- foi[j]
    
    #quaternary infections -  S342
    mat1[228, 228] <- 1-2*foi[j]
    mat1[319,228] <- foi[j]
    mat1[320,228] <- foi[j]
    
    #quaternary infections -  S345
    mat1[229, 229] <- 1-2*foi[j]
    mat1[321,229] <- foi[j]
    mat1[322,229] <- foi[j]
    
    
    #quaternary infections -  S351
    mat1[230, 230] <- 1-2*foi[j]
    mat1[323,230] <- foi[j]
    mat1[324,230] <- foi[j]
    
    #quaternary infections -  S352
    mat1[231, 231] <- 1-2*foi[j]
    mat1[325,231] <- foi[j]
    mat1[326,231] <- foi[j]
    
    #quaternary infections -  S354
    mat1[232, 232] <- 1-2*foi[j]
    mat1[327,232] <- foi[j]
    mat1[328,232] <- foi[j]
    
    #quaternary infections -  S412
    mat1[233, 233] <- 1-2*foi[j]
    mat1[329,233] <- foi[j]
    mat1[330,233] <- foi[j]
    
    #quaternary infections -  S413
    mat1[234, 234] <- 1-2*foi[j]
    mat1[331,234] <- foi[j]
    mat1[332,234] <- foi[j]
    
    #quaternary infections -  S415
    mat1[235, 235] <- 1-2*foi[j]
    mat1[333,235] <- foi[j]
    mat1[334,235] <- foi[j]
    
    #quaternary infections -  S421
    mat1[236, 236] <- 1-2*foi[j]
    mat1[335,236] <- foi[j]
    mat1[336,236] <- foi[j]
    
    #quaternary infections -  S423
    mat1[237, 237] <- 1-2*foi[j]
    mat1[337,237] <- foi[j]
    mat1[338,237] <- foi[j]
    
    
    #quaternary infections -  S425
    mat1[238, 238] <- 1-2*foi[j]
    mat1[339,238] <- foi[j]
    mat1[340,238] <- foi[j]
    
    
    #quaternary infections -  S431
    mat1[239, 239] <- 1-2*foi[j]
    mat1[341,239] <- foi[j]
    mat1[342,239] <- foi[j]
    
    #quaternary infections -  S432
    mat1[240, 240] <- 1-2*foi[j]
    mat1[343,240] <- foi[j]
    mat1[344,240] <- foi[j]
    
    #quaternary infections -  S435
    mat1[241, 241] <- 1-2*foi[j]
    mat1[345,241] <- foi[j]
    mat1[346,241] <- foi[j]
    
    #quaternary infections -  S451
    mat1[242, 242] <- 1-2*foi[j]
    mat1[347,242] <- foi[j]
    mat1[348,242] <- foi[j]
    
    #quaternary infections -  S452
    mat1[243, 243] <- 1-2*foi[j]
    mat1[349,243] <- foi[j]
    mat1[350,243] <- foi[j]
    
    #quaternary infections -  S453
    mat1[244, 244] <- 1-2*foi[j]
    mat1[351,244] <- foi[j]
    mat1[352,244] <- foi[j]
    
    #quaternary infections -  S512
    mat1[245, 245] <- 1-2*foi[j]
    mat1[353,245] <- foi[j]
    mat1[354,245] <- foi[j]
    
    #quaternary infections -  S513
    mat1[246, 246] <- 1-2*foi[j]
    mat1[355,246] <- foi[j]
    mat1[356,246] <- foi[j]
    
    #quaternary infections -  S514
    mat1[247, 247] <- 1-2*foi[j]
    mat1[357,247] <- foi[j]
    mat1[358,247] <- foi[j]
    
    #quaternary infections -  S521
    mat1[248, 248] <- 1-2*foi[j]
    mat1[359,248] <- foi[j]
    mat1[360,248] <- foi[j]
    
    #quaternary infections -  S523
    mat1[249, 249] <- 1-2*foi[j]
    mat1[361,249] <- foi[j]
    mat1[362,249] <- foi[j]
    
    #quaternary infections -  S524
    mat1[250, 250] <- 1-2*foi[j]
    mat1[363,250] <- foi[j]
    mat1[364,250] <- foi[j]
    
    #quaternary infections -  S531
    mat1[251, 251] <- 1-2*foi[j]
    mat1[365,251] <- foi[j]
    mat1[366,251] <- foi[j]
    
    #quaternary infections -  S532
    mat1[252, 252] <- 1-2*foi[j]
    mat1[367,252] <- foi[j]
    mat1[368,252] <- foi[j]
    
    #quaternary infections -  S534
    mat1[253, 253] <- 1-2*foi[j]
    mat1[369,253] <- foi[j]
    mat1[370,253] <- foi[j]
    
    #quaternary infections -  S541
    mat1[254, 254] <- 1-2*foi[j]
    mat1[371,254] <- foi[j]
    mat1[372,254] <- foi[j]
    
    #quaternary infections -  S542
    mat1[255, 255] <- 1-2*foi[j]
    mat1[373,255] <- foi[j]
    mat1[374,255] <- foi[j]
    
    #quaternary infections -  S543
    mat1[256, 256] <- 1-2*foi[j]
    mat1[375,256] <- foi[j]
    mat1[376,256] <- foi[j]
    
    
    #recovery from quaternary infections into shorttrm heterotypic immunity
    
    for(i in 0:119){
      #recovery from quaternary infections into shorttrm heterotypic immunity
      mat1[(257+i), (257+i)] <- 1-recov
      mat1[(377+i), (257+i)] <- recov 
      
      #waning out of heterotypic immunity
      mat1[(377+i), (377+i)] <- 1- wane_hetero[j]  
      mat1[(497+i), (377+i)] <- wane_hetero[j]  
      
      #5th order infections - only one option for each
      mat1[(497+i), (497+i)] <- 1- foi[j]  
      mat1[(617+i), (497+i)] <- foi[j]  
      
      #recovery from 5th-order infections into longterm heterotypic immunity to all serotypes (Pm = 736) 
      mat1[(617+i), (617+i)] <-  1-recov
      mat1[(737), (617+i)] <- recov
      
    }
   
    #once you reach immunity, you stay
    mat1[737, 737] <- 1
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*c+1):((j-1)*c+c),((j-1)*c+1):((j-1)*c+c)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
  return(Tmat)
  
}
buildTMat_five_strain_after <- function(c, Npop, age.classes, surv.biwk, wane_hetero, age.brk, foi, age.mult.foi, recov, sigma, age.rate){
  
  
  #remove foi for serotype 2, as it is eliminated
  
  s <- nage <- length(age.classes)
  
  
  
  if (length(foi)==1) foi <- rep(foi,nage) # 
  if (length(foi)==s) foi <- foi
  if (length(foi) > 1 & length(foi)<s){
    foi.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(foi)){
      foi.list[[i]] = rep(foi[i], age.dur[i])
    }
    foi = c(unlist(foi.list))
  } 
  
  
  if (length(wane_hetero)==1) wane_hetero <- rep(wane_hetero,nage) # 
  if (length(wane_hetero)==s) wane_hetero <- wane_hetero
  if (length(wane_hetero) > 1 & length(wane_hetero)<s){
    wane_hetero.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(wane_hetero)){
      wane_hetero.list[[i]] = rep(wane_hetero[i], age.dur[i])
    }
    wane_hetero = c(unlist(wane_hetero.list))
  }   
  
  #then, now that foi is appropriate length of all the age classes, 
  #multiply by the age.mult vector
  foi <- foi*age.mult.foi
  
  if (length(sigma)==1) sigma <- rep(sigma,nage) # 
  if (length(sigma)==s) sigma <- sigma
  if (length(sigma) > 1 & length(sigma)<s){
    sigma.list <- list()
    
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(sigma)){
      sigma.list[[i]] = rep(sigma[i], age.dur[i])
    }
    sigma = c(unlist(sigma.list))
  } 
  
  
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  
  
  
  mat1 <- matrix(0,c,c) 
  
  Tmat <- matrix(0,c*nage,c*nage) 
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek
    mat1[] <- 0
    
    #infections in naive susceptibles
    mat1[1,1] <- 1- (4*foi[j])
    #mat1[2,1] <- foi[j]
    mat1[3,1] <- foi[j]
    mat1[4,1] <- foi[j]
    mat1[5,1] <- foi[j]
    mat1[6,1] <- foi[j]
    
    #I1s recover
    mat1[2,2] <- 1-recov
    mat1[7,2] <- recov
    
    
    #I2s recover
    mat1[3,3] <- 1-recov
    mat1[8,3] <- recov
    
    #I3s recover
    mat1[4,4] <- 1-recov
    mat1[9,4] <- recov
    
    #I4s recover
    mat1[5,5] <- 1-recov
    mat1[10,5] <- recov
    
    #I5s recover
    mat1[6,6] <- 1-recov
    mat1[11,6] <- recov
    
    
    #waning out of heterotypic immunity (P1 to S1) - those previously infected with I1
    mat1[7,7] <- 1- wane_hetero[j]
    mat1[12,7] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I2
    mat1[8,8] <- 1- wane_hetero[j]
    mat1[13,8] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I3
    mat1[9,9] <- 1- wane_hetero[j]
    mat1[14,9] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I4
    mat1[10,10] <- 1- wane_hetero[j]
    mat1[15,10] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I5
    mat1[11,11] <- 1- wane_hetero[j]
    mat1[16,11] <- wane_hetero[j]
    
    #secondary infections - those previously infected with I1
    mat1[12,12] <- 1-(3*foi[j])
    #mat1[17,12] <- foi[j] #I12
    mat1[18,12] <- foi[j] #I13
    mat1[19,12] <- foi[j] #I14
    mat1[20,12] <- foi[j] #I15
    
    #secondary infections - those previously infected with I2
    mat1[13,13] <- 1-(4*foi[j])
    mat1[21,13] <- foi[j] #I21
    mat1[22,13] <- foi[j] #I22
    mat1[23,13] <- foi[j] #I23
    mat1[24,13] <- foi[j] #I25
    
    #secondary infections - those previously infected with I3
    mat1[14,14] <-  1-(3*foi[j])
    mat1[25,14] <- foi[j] #I31
    #mat1[26,14] <- foi[j] #I32
    mat1[27,14] <- foi[j] #I34
    mat1[28,14] <- foi[j] #I35
    
    #secondary infections - those previously infected with I4
    mat1[15,15] <-  1-(3*foi[j])
    mat1[29,15] <- foi[j] #I41
    #mat1[30,15] <- foi[j] #I42
    mat1[31,15] <- foi[j] #I43
    mat1[32,15] <- foi[j] #I45
    
    
    #secondary infections - those previously infected with I5
    mat1[16,16] <-  1-(3*foi[j])
    mat1[33,16] <- foi[j] #I51
    #mat1[34,16] <- foi[j] #I52
    mat1[35,16] <- foi[j] #I53
    mat1[36,16] <- foi[j] #I54
    
    #recovery from secondary infections - I12
    mat1[17,17] <- 1-recov
    mat1[37,17] <- recov
    
    #recovery from secondary infections - I13
    mat1[18,18] <- 1-recov
    mat1[38,18] <- recov
    
    #recovery from secondary infections - I14
    mat1[19,19] <- 1-recov
    mat1[39,19] <- recov
    
    #recovery from secondary infections - I15
    mat1[20,20] <- 1-recov
    mat1[40,20] <- recov
    
    #recovery from secondary infections - I21
    mat1[21,21] <- 1-recov
    mat1[41,21] <- recov
    
    #recovery from secondary infections - I23
    mat1[22,22] <- 1-recov
    mat1[42,22] <- recov
    
    #recovery from secondary infections - I24
    mat1[23,23] <- 1-recov
    mat1[43,23] <- recov
    
    #recovery from secondary infections - I25
    mat1[24,24] <- 1-recov
    mat1[44,24] <- recov
    
    #recovery from secondary infections - I31
    mat1[25,25] <- 1-recov
    mat1[45,25] <- recov
    
    #recovery from secondary infections - I32
    mat1[26,26] <- 1-recov
    mat1[46,26] <- recov
    
    #recovery from secondary infections - I34 to P34
    mat1[27,27] <- 1-recov
    mat1[47,27] <- recov
    
    #recovery from secondary infections - I35 to P35
    mat1[28,28] <- 1-recov
    mat1[48,28] <- recov
    
    #recovery from secondary infections - I41 to P41
    mat1[29,29] <- 1-recov
    mat1[49,29] <- recov
    
    #recovery from secondary infections - I42 to P42
    mat1[30,30] <- 1-recov
    mat1[50,30] <- recov
    
    #recovery from secondary infections - I43 to P43
    mat1[31,31] <- 1-recov
    mat1[51,31] <- recov
    
    #recovery from secondary infections - I45 to P45
    mat1[32,32] <- 1-recov
    mat1[52,32] <- recov
    
    
    #recovery from secondary infections - I51 to P51
    mat1[33,33] <- 1-recov
    mat1[53,33] <- recov
    
    #recovery from secondary infections - I52 to P52
    mat1[34,34] <- 1-recov
    mat1[54,34] <- recov
    
    #recovery from secondary infections - I53 to P53
    mat1[35,35] <- 1-recov
    mat1[55,35] <- recov
    
    #recovery from secondary infections - I54 to P54
    mat1[36,36] <- 1-recov
    mat1[56,36] <- recov
    
    #waning out of heterotypic immunity - those previously infected with I12 (P12 to S12)
    mat1[37,37] <- 1- wane_hetero[j]
    mat1[57,37] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I13
    mat1[38,38] <- 1- wane_hetero[j]
    mat1[58,38] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I14 (P14 to S14)
    mat1[39,39] <- 1- wane_hetero[j]
    mat1[59,39] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I15 (P15 to S15)
    mat1[40,40] <- 1- wane_hetero[j]
    mat1[60,40] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I21
    mat1[41,41] <- 1- wane_hetero[j]
    mat1[61,41] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I23
    mat1[42,42] <- 1- wane_hetero[j]
    mat1[62,42] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I24
    mat1[43,43] <- 1- wane_hetero[j]
    mat1[63,43] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I25
    mat1[44,44] <- 1- wane_hetero[j]
    mat1[64,44] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I31
    mat1[45,45] <- 1- wane_hetero[j]
    mat1[65,45] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I32
    mat1[46,46] <- 1- wane_hetero[j]
    mat1[66,46] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I34
    mat1[47,47] <- 1- wane_hetero[j]
    mat1[67,47] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I35
    mat1[48,48] <- 1- wane_hetero[j]
    mat1[68,48] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I41
    mat1[49,49] <- 1- wane_hetero[j]
    mat1[69,49] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I42
    mat1[50,50] <- 1- wane_hetero[j]
    mat1[70,50] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I43
    mat1[51,51] <- 1- wane_hetero[j]
    mat1[71,51] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I45
    mat1[52,52] <- 1- wane_hetero[j]
    mat1[72,52] <- wane_hetero[j]
    
    
    #waning out of heterotypic immunity - those previously infected with I51
    mat1[53,53] <- 1- wane_hetero[j]
    mat1[73,53] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I52
    mat1[54,54] <- 1- wane_hetero[j]
    mat1[74,54] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I53
    mat1[55,55] <- 1- wane_hetero[j]
    mat1[75,55] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I54
    mat1[56,56] <- 1- wane_hetero[j]
    mat1[76,56] <- wane_hetero[j]
    
    
    #tertiary infections - those previously infected with I12 (S12 to I123 and I124 and I125)
    mat1[57,57] <- 1-3*foi[j]
    mat1[77,57] <- foi[j]
    mat1[78,57] <- foi[j]
    mat1[79,57] <- foi[j]
    
    #tertiary infections - those previously infected with I13 (I132 and I134 and I135)
    mat1[58,58] <- 1-2*foi[j]
    #mat1[80,58] <- foi[j]
    mat1[81,58] <- foi[j]
    mat1[82,58] <- foi[j]
    
    
    #tertiary infections - those previously infected with I14 (I142 and I143 and I145)
    mat1[59,59] <- 1-2*foi[j]
    #mat1[83,59] <- foi[j]
    mat1[84,59] <- foi[j]
    mat1[85,59] <- foi[j]
    
    #tertiary infections - those previously infected with I15 (I152 and I153 and I154)
    mat1[60,60] <- 1-2*foi[j]
    #mat1[86,60] <- foi[j]
    mat1[87,60] <- foi[j]
    mat1[88,60] <- foi[j]
    
    #tertiary infections - those previously infected with I21 (I213 + I214 + I215)
    mat1[61,61] <- 1-3*foi[j]
    mat1[89,61] <- foi[j]
    mat1[90,61] <- foi[j]
    mat1[91,61] <- foi[j]
    
    #tertiary infections - those previously infected with I23 (I231 + I234 + I235)
    mat1[62,62] <- 1-3*foi[j]
    mat1[92,62] <- foi[j]
    mat1[93,62] <- foi[j]
    mat1[94,62] <- foi[j]
    
    #tertiary infections - those previously infected with I24 (I241 + I243 + I245)
    mat1[63,63] <- 1-3*foi[j]
    mat1[95,63] <- foi[j]
    mat1[96,63] <- foi[j]
    mat1[97,63] <- foi[j]
    
    #tertiary infections - those previously infected with I25 (I251+I253+I254)
    mat1[64,64] <- 1-3*foi[j]
    mat1[98,64] <- foi[j]
    mat1[99,64] <- foi[j]
    mat1[100,64] <- foi[j]
    
    #tertiary infections - those previously infected with I31 (I312+I314+I315)
    mat1[65,65] <- 1-2*foi[j]
    #mat1[101,65] <- foi[j]
    mat1[102,65] <- foi[j]
    mat1[103,65] <- foi[j]
    
    
    #tertiary infections - those previously infected with I32 (I321+I324+I325)
    mat1[66,66] <- 1-3*foi[j]
    mat1[104,66] <- foi[j]
    mat1[105,66] <- foi[j]
    mat1[106,66] <- foi[j]
    
    #tertiary infections - those previously infected with I34 (I341+I342+I345)
    mat1[67,67] <- 1-2*foi[j]
    mat1[107,67] <- foi[j]
    #mat1[108,67] <- foi[j]
    mat1[109,67] <- foi[j]
    
    #tertiary infections - those previously infected with I35 (I351+I352+I354)
    mat1[68,68] <- 1-2*foi[j]
    mat1[110,68] <- foi[j]
    #mat1[111,68] <- foi[j]
    mat1[112,68] <- foi[j]
    
    #tertiary infections - those previously infected with I41 (I412+I413+415)
    mat1[69,69] <- 1-2*foi[j]
    #mat1[113,69] <- foi[j]
    mat1[114,69] <- foi[j]
    mat1[115,69] <- foi[j]
    
    #tertiary infections - those previously infected with I42 (I421+I423+I425)
    mat1[70,70] <- 1-3*foi[j]
    mat1[116,70] <- foi[j]
    mat1[117,70] <- foi[j]
    mat1[118,70] <- foi[j]
    
    #tertiary infections - those previously infected with I43 (I431+I432+I435)
    mat1[71,71] <- 1-2*foi[j]
    mat1[119,71] <- foi[j]
    #mat1[120,71] <- foi[j]
    mat1[121,71] <- foi[j]
    
    #tertiary infections - those previously infected with I45 (I451+I452+I453)
    mat1[72,72] <- 1-2*foi[j]
    mat1[122,72] <- foi[j]
    #mat1[123,72] <- foi[j]
    mat1[124,72] <- foi[j]
    
    
    #tertiary infections - those previously infected with I51 (I512+I513+514)
    mat1[73,73] <- 1-2*foi[j]
    #mat1[125,73] <- foi[j]
    mat1[126,73] <- foi[j]
    mat1[127,73] <- foi[j]
    
    
    #tertiary infections - those previously infected with I52 (I521+I523+524)
    mat1[74,74] <- 1-3*foi[j]
    mat1[128,74] <- foi[j]
    mat1[129,74] <- foi[j]
    mat1[130,74] <- foi[j]
    
    #tertiary infections - those previously infected with I53 (I531+I532+534)
    mat1[75,75] <- 1-2*foi[j]
    mat1[131,75] <- foi[j]
    #mat1[132,75] <- foi[j]
    mat1[133,75] <- foi[j]
    
    
    #tertiary infections - those previously infected with I54 (I541+I542+543)
    mat1[76,76] <- 1-2*foi[j]
    mat1[134,76] <- foi[j]
    #mat1[135,76] <- foi[j]
    mat1[136,76] <- foi[j]
    
    
    #recovery from tertiary infections - I123
    mat1[77,77] <- 1-recov
    mat1[137,77] <- recov
    
    #recovery from tertiary infections - I124
    mat1[78,78] <- 1-recov
    mat1[138,78] <- recov
    
    #recovery from tertiary infections - I125
    mat1[79,79] <- 1-recov
    mat1[139,79] <- recov
    
    #recovery from tertiary infections - I132
    mat1[80,80] <- 1-recov
    mat1[140,80] <- recov
    
    #recovery from tertiary infections - I134
    mat1[81,81] <- 1-recov
    mat1[141,81] <- recov
    
    #recovery from tertiary infections - I135
    mat1[82,82] <- 1-recov
    mat1[142,82] <- recov
    
    #recovery from tertiary infections - I142
    mat1[83,83] <- 1-recov
    mat1[143,83] <- recov
    
    #recovery from tertiary infections - I143
    mat1[84,84] <- 1-recov
    mat1[144,84] <- recov
    
    #recovery from tertiary infections - I145
    mat1[85,85] <- 1-recov
    mat1[145,85] <- recov
    
    
    #recovery from tertiary infections - I152
    mat1[86,86] <- 1-recov
    mat1[146,86] <- recov
    
    #recovery from tertiary infections - I153
    mat1[87,87] <- 1-recov
    mat1[147,87] <- recov
    
    #recovery from tertiary infections - I154
    mat1[88,88] <- 1-recov
    mat1[148,88] <- recov
    
    #recovery from tertiary infections - I213
    mat1[89,89] <- 1-recov
    mat1[149,89] <- recov
    
    #recovery from tertiary infections - I214
    mat1[90,90] <- 1-recov
    mat1[150,90] <- recov
    
    #recovery from tertiary infections - I215
    mat1[91,91] <- 1-recov
    mat1[151,91] <- recov
    
    
    #recovery from tertiary infections - I231
    mat1[92,92] <- 1-recov
    mat1[152,92] <- recov
    
    #recovery from tertiary infections - I234
    mat1[93,93] <- 1-recov
    mat1[153,93] <- recov
    
    #recovery from tertiary infections - I235
    mat1[94,94] <- 1-recov
    mat1[154,94] <- recov
    
    #recovery from tertiary infections - I241
    mat1[95,95] <- 1-recov
    mat1[155,95] <- recov
    
    #recovery from tertiary infections - I243
    mat1[96,96] <- 1-recov
    mat1[156,96] <- recov
    
    #recovery from tertiary infections - I245
    mat1[97,97] <- 1-recov
    mat1[157,97] <- recov
    
    
    #recovery from tertiary infections - I251
    mat1[98,98] <- 1-recov
    mat1[158,98] <- recov
    
    #recovery from tertiary infections - I253
    mat1[99,99] <- 1-recov
    mat1[159,99] <- recov
    
    #recovery from tertiary infections - I254
    mat1[100,100] <- 1-recov
    mat1[160,100] <- recov
    
    #recovery from tertiary infections - I312
    mat1[101,101] <- 1-recov
    mat1[161,101] <- recov
    
    #recovery from tertiary infections - I314
    mat1[102,102] <- 1-recov
    mat1[162,102] <- recov
    
    #recovery from tertiary infections - I315
    mat1[103,103] <- 1-recov
    mat1[163,103] <- recov
    
    #recovery from tertiary infections - I321
    mat1[104,104] <- 1-recov
    mat1[164,104] <- recov
    
    #recovery from tertiary infections - I324
    mat1[105,105] <- 1-recov
    mat1[165,105] <- recov
    
    #recovery from tertiary infections - I325
    mat1[106,106] <- 1-recov
    mat1[166,106] <- recov
    
    #recovery from tertiary infections - I341
    mat1[107,107] <- 1-recov
    mat1[167,107] <- recov
    
    #recovery from tertiary infections - I342
    mat1[108,108] <- 1-recov
    mat1[168,108] <- recov
    
    #recovery from tertiary infections - I345
    mat1[109,109] <- 1-recov
    mat1[169,109] <- recov
    
    
    #recovery from tertiary infections - I351
    mat1[110,110] <- 1-recov
    mat1[170,110] <- recov
    
    #recovery from tertiary infections - I352
    mat1[111,111] <- 1-recov
    mat1[171,111] <- recov
    
    #recovery from tertiary infections - I354
    mat1[112,112] <- 1-recov
    mat1[172,112] <- recov
    
    
    #recovery from tertiary infections - I412
    mat1[113,113] <- 1-recov
    mat1[173,113] <- recov
    
    #recovery from tertiary infections - I413
    mat1[114,114] <- 1-recov
    mat1[174,114] <- recov
    
    #recovery from tertiary infections - I415
    mat1[115,115] <- 1-recov
    mat1[175,115] <- recov
    
    #recovery from tertiary infections - I421
    mat1[116,116] <- 1-recov
    mat1[176,116] <- recov
    
    #recovery from tertiary infections - I423
    mat1[117,117] <- 1-recov
    mat1[177,117] <- recov
    
    #recovery from tertiary infections - I425
    mat1[118,118] <- 1-recov
    mat1[178,118] <- recov
    
    #recovery from tertiary infections - I431
    mat1[119,119] <- 1-recov
    mat1[179,119] <- recov
    
    #recovery from tertiary infections - I432
    mat1[120,120] <- 1-recov
    mat1[180,120] <- recov
    
    #recovery from tertiary infections - I435
    mat1[121,121] <- 1-recov
    mat1[181,121] <- recov
    
    
    #recovery from tertiary infections - I451
    mat1[122,122] <- 1-recov
    mat1[182,122] <- recov
    
    #recovery from tertiary infections - I452
    mat1[123,123] <- 1-recov
    mat1[183,123] <- recov
    
    #recovery from tertiary infections - I453
    mat1[124,124] <- 1-recov
    mat1[184,124] <- recov
    
    #recovery from tertiary infections - I512
    mat1[125,125] <- 1-recov
    mat1[185,125] <- recov
    
    #recovery from tertiary infections - I513
    mat1[126,126] <- 1-recov
    mat1[186,126] <- recov
    
    #recovery from tertiary infections - I514
    mat1[127,127] <- 1-recov
    mat1[187,127] <- recov
    
    #recovery from tertiary infections - I521
    mat1[128,128] <- 1-recov
    mat1[188,128] <- recov
    
    #recovery from tertiary infections - I523
    mat1[129,129] <- 1-recov
    mat1[189,129] <- recov
    
    #recovery from tertiary infections - I524
    mat1[130,130] <- 1-recov
    mat1[190,130] <- recov
    
    #recovery from tertiary infections - I531
    mat1[131,131] <- 1-recov
    mat1[191,131] <- recov
    
    #recovery from tertiary infections - I532
    mat1[132,132] <- 1-recov
    mat1[192,132] <- recov
    
    #recovery from tertiary infections - I534
    mat1[133,133] <- 1-recov
    mat1[193,133] <- recov
    
    #recovery from tertiary infections - I541
    mat1[134,134] <- 1-recov
    mat1[194,134] <- recov
    
    #recovery from tertiary infections - I542
    mat1[135,135] <- 1-recov
    mat1[195,135] <- recov
    
    #recovery from tertiary infections - I543
    mat1[136,136] <- 1-recov
    mat1[196,136] <- recov
    
    
    #waning out of heterotypic immunity - P123 to S123
    mat1[137,137] <- 1- wane_hetero[j]
    mat1[197,137] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P124 to S124
    mat1[138,138] <- 1- wane_hetero[j]
    mat1[198,138] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P125 to S125
    mat1[139,139] <- 1- wane_hetero[j]
    mat1[199,139] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P132 to S132
    mat1[140,140] <- 1- wane_hetero[j]
    mat1[200,140] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P134 to S134
    mat1[141,141] <- 1- wane_hetero[j]
    mat1[201,141] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P135 to S135
    mat1[142,142] <- 1- wane_hetero[j]
    mat1[202,142] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P142
    mat1[143,143] <- 1- wane_hetero[j]
    mat1[203,143] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P143
    mat1[144,144] <- 1- wane_hetero[j]
    mat1[204,144] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P145
    mat1[145,145] <- 1- wane_hetero[j]
    mat1[205,145] <- wane_hetero[j]
    
    
    #waning out of heterotypic immunity - P152
    mat1[146,146] <- 1- wane_hetero[j]
    mat1[206,146] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P153
    mat1[147,147] <- 1- wane_hetero[j]
    mat1[207,147] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P154
    mat1[148,148] <- 1- wane_hetero[j]
    mat1[208,148] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P213
    mat1[149,149] <- 1- wane_hetero[j]
    mat1[209,149] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P214
    mat1[150,150] <- 1- wane_hetero[j]
    mat1[210,150] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P215
    mat1[151,151] <- 1- wane_hetero[j]
    mat1[211,151] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P231
    mat1[152,152] <- 1- wane_hetero[j]
    mat1[212,152] <- wane_hetero[j]
    
    
    #waning out of heterotypic immunity - P234
    mat1[153,153] <- 1- wane_hetero[j]
    mat1[213,153] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P235
    mat1[154,154] <- 1- wane_hetero[j]
    mat1[214,154] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P241
    mat1[155,155] <- 1- wane_hetero[j]
    mat1[215,155] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P243
    mat1[156,156] <- 1- wane_hetero[j]
    mat1[216,156] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P245
    mat1[157,157] <- 1- wane_hetero[j]
    mat1[217,157] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P251
    mat1[158,158] <- 1- wane_hetero[j]
    mat1[218,158] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P253
    mat1[159,159] <- 1- wane_hetero[j]
    mat1[219,159] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P254
    mat1[160,160] <- 1- wane_hetero[j]
    mat1[220,160] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P312
    mat1[161,161] <- 1- wane_hetero[j]
    mat1[221,161] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P314
    mat1[162,162] <- 1- wane_hetero[j]
    mat1[222,162] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P315
    mat1[163,163] <- 1- wane_hetero[j]
    mat1[223,163] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P321
    mat1[164,164] <- 1- wane_hetero[j]
    mat1[224,164] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P324
    mat1[165,165] <- 1- wane_hetero[j]
    mat1[225,165] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P325
    mat1[166,166] <- 1- wane_hetero[j]
    mat1[226,166] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P341
    mat1[167,167] <- 1- wane_hetero[j]
    mat1[227,167] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P342
    mat1[168,168] <- 1- wane_hetero[j]
    mat1[228,168] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P345
    mat1[169,169] <- 1- wane_hetero[j]
    mat1[229,169] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P351
    mat1[170,170] <- 1- wane_hetero[j]
    mat1[230,170] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P352
    mat1[171,171] <- 1- wane_hetero[j]
    mat1[231,171] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P354
    mat1[172,172] <- 1- wane_hetero[j]
    mat1[232,172] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P412
    mat1[173,173] <- 1- wane_hetero[j]
    mat1[233,173] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P413
    mat1[174,174] <- 1- wane_hetero[j]
    mat1[234,174] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P415
    mat1[175,175] <- 1- wane_hetero[j]
    mat1[235,175] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P421
    mat1[176,176] <- 1- wane_hetero[j]
    mat1[236,176] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P423
    mat1[177,177] <- 1- wane_hetero[j]
    mat1[237,177] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P425
    mat1[178,178] <- 1- wane_hetero[j]
    mat1[238,178] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P431
    mat1[179,179] <- 1- wane_hetero[j]
    mat1[239,179] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P432
    mat1[180,180] <- 1- wane_hetero[j]
    mat1[240,180] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P435
    mat1[181,181] <- 1- wane_hetero[j]
    mat1[241,181] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P451
    mat1[182,182] <- 1- wane_hetero[j]
    mat1[242,182] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P452
    mat1[183,183] <- 1- wane_hetero[j]
    mat1[243,183] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P453
    mat1[184,184] <- 1- wane_hetero[j]
    mat1[244,184] <- wane_hetero[j]
    
    
    
    #waning out of heterotypic immunity - P512
    mat1[185,185] <- 1- wane_hetero[j]
    mat1[245,185] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P513
    mat1[186,186] <- 1- wane_hetero[j]
    mat1[246,186] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P514
    mat1[187,187] <- 1- wane_hetero[j]
    mat1[247,187] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P521
    mat1[188,188] <- 1- wane_hetero[j]
    mat1[248,188] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P523
    mat1[189,189] <- 1- wane_hetero[j]
    mat1[249,189] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P524
    mat1[190,190] <- 1- wane_hetero[j]
    mat1[250,190] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P531
    mat1[191,191] <- 1- wane_hetero[j]
    mat1[251,191] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P532
    mat1[192,192] <- 1- wane_hetero[j]
    mat1[252,192] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P534
    mat1[193,193] <- 1- wane_hetero[j]
    mat1[253,193] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P541
    mat1[194,194] <- 1- wane_hetero[j]
    mat1[254,194] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P542
    mat1[195,195] <- 1- wane_hetero[j]
    mat1[255,195] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - P543
    mat1[196,196] <- 1- wane_hetero[j]
    mat1[256,196] <- wane_hetero[j]
    
    
    #quaternary infections -  S123 to I1234 or I1235
    mat1[197, 197] <- 1-2*foi[j]
    mat1[257,197] <- foi[j]
    mat1[258,197] <- foi[j]
    
    #quaternary infections -  S124 to I1243 or I1245
    mat1[198, 198] <- 1-2*foi[j]
    mat1[259,198] <- foi[j]
    mat1[260,198] <- foi[j]
    
    #quaternary infections -  S125
    mat1[199, 199] <- 1-2*foi[j]
    mat1[261,199] <- foi[j]
    mat1[262,199] <- foi[j]
    
    #quaternary infections -  S132 to I1324 or I1325
    mat1[200, 200] <- 1-2*foi[j]
    mat1[263,200] <- foi[j]
    mat1[264,200] <- foi[j]
    
    #quaternary infections -  S134 to I1342 or I1345
    mat1[201, 201] <- 1-foi[j]
    #mat1[265,201] <- foi[j]
    mat1[266,201] <- foi[j]
    
    #quaternary infections -  S135
    mat1[202, 202] <- 1-foi[j]
    #mat1[267,202] <- foi[j]
    mat1[268,202] <- foi[j]
    
    #quaternary infections -  S142
    mat1[203, 203] <- 1-2*foi[j]
    mat1[269,203] <- foi[j]
    mat1[270,203] <- foi[j]
    
    #quaternary infections -  S143
    mat1[204, 204] <- 1-foi[j]
    #mat1[271,204] <- foi[j]
    mat1[272,204] <- foi[j]
    
    #quaternary infections -  S145
    mat1[205, 205] <- 1-foi[j]
    #mat1[273,205] <- foi[j]
    mat1[274,205] <- foi[j]
    
    #quaternary infections -  S152
    mat1[206, 206] <- 1-2*foi[j]
    mat1[275,206] <- foi[j]
    mat1[276,206] <- foi[j]
    
    #quaternary infections -  S153
    mat1[207, 207] <- 1-foi[j]
    #mat1[277,207] <- foi[j]
    mat1[278,207] <- foi[j]
    
    #quaternary infections -  S154
    mat1[208, 208] <- 1-foi[j]
    #mat1[279,208] <- foi[j]
    mat1[280,208] <- foi[j]
    
    #quaternary infections -  S213
    mat1[209, 209] <- 1-2*foi[j]
    mat1[281,209] <- foi[j]
    mat1[282,209] <- foi[j]
    
    #quaternary infections -  S214
    mat1[210, 210] <- 1-2*foi[j]
    mat1[283,210] <- foi[j]
    mat1[284,210] <- foi[j]
    
    #quaternary infections -  S215
    mat1[211, 211] <- 1-2*foi[j]
    mat1[285,211] <- foi[j]
    mat1[286,211] <- foi[j]
    
    #quaternary infections -  S231
    mat1[212, 212] <- 1-2*foi[j]
    mat1[287,212] <- foi[j]
    mat1[288,212] <- foi[j]
    
    #quaternary infections -  S234
    mat1[213, 213] <- 1-2*foi[j]
    mat1[289,213] <- foi[j]
    mat1[290,213] <- foi[j]
    
    #quaternary infections -  S235
    mat1[214, 214] <- 1-2*foi[j]
    mat1[291,214] <- foi[j]
    mat1[292,214] <- foi[j]
    
    #quaternary infections -  S241
    mat1[215, 215] <- 1-2*foi[j]
    mat1[293,215] <- foi[j]
    mat1[294,215] <- foi[j]
    
    #quaternary infections -  S243
    mat1[216, 216] <- 1-2*foi[j]
    mat1[295,216] <- foi[j]
    mat1[296,216] <- foi[j]
    
    #quaternary infections -  S245
    mat1[217, 217] <- 1-2*foi[j]
    mat1[297,217] <- foi[j]
    mat1[298,217] <- foi[j]
    
    
    #quaternary infections -  S251
    mat1[218, 218] <- 1-2*foi[j]
    mat1[299,218] <- foi[j]
    mat1[300,218] <- foi[j]
    
    
    #quaternary infections -  S253
    mat1[219, 219] <- 1-2*foi[j]
    mat1[301,219] <- foi[j]
    mat1[302,219] <- foi[j]
    
    #quaternary infections -  S254
    mat1[220, 220] <- 1-2*foi[j]
    mat1[303,220] <- foi[j]
    mat1[304,220] <- foi[j]
    
    #quaternary infections -  S312
    mat1[221, 221] <- 1-2*foi[j]
    mat1[305,221] <- foi[j]
    mat1[306,221] <- foi[j]
    
    #quaternary infections -  S314
    mat1[222, 222] <- 1-foi[j]
    #mat1[307,222] <- foi[j]
    mat1[308,222] <- foi[j]
    
    #quaternary infections -  S315
    mat1[223, 223] <- 1-foi[j]
    #mat1[309,223] <- foi[j]
    mat1[310,223] <- foi[j]
    
    #quaternary infections -  S321
    mat1[224, 224] <- 1-2*foi[j]
    mat1[311,224] <- foi[j]
    mat1[312,224] <- foi[j]
    
    #quaternary infections -  S324
    mat1[225, 225] <- 1-2*foi[j]
    mat1[313,225] <- foi[j]
    mat1[314,225] <- foi[j]
    
    #quaternary infections -  S325
    mat1[226, 226] <- 1-2*foi[j]
    mat1[315,226] <- foi[j]
    mat1[316,226] <- foi[j]
    
    #quaternary infections -  S341
    mat1[227, 227] <- 1-foi[j]
    #mat1[317,227] <- foi[j]
    mat1[318,227] <- foi[j]
    
    #quaternary infections -  S342
    mat1[228, 228] <- 1-2*foi[j]
    mat1[319,228] <- foi[j]
    mat1[320,228] <- foi[j]
    
    #quaternary infections -  S345
    mat1[229, 229] <- 1-foi[j]
    mat1[321,229] <- foi[j]
    #mat1[322,229] <- foi[j]
    
    
    #quaternary infections -  S351
    mat1[230, 230] <- 1-foi[j]
    #mat1[323,230] <- foi[j]
    mat1[324,230] <- foi[j]
    
    #quaternary infections -  S352
    mat1[231, 231] <- 1-2*foi[j]
    mat1[325,231] <- foi[j]
    mat1[326,231] <- foi[j]
    
    #quaternary infections -  S354
    mat1[232, 232] <- 1-foi[j]
    mat1[327,232] <- foi[j]
    #mat1[328,232] <- foi[j]
    
    #quaternary infections -  S412
    mat1[233, 233] <- 1-2*foi[j]
    mat1[329,233] <- foi[j]
    mat1[330,233] <- foi[j]
    
    #quaternary infections -  S413
    mat1[234, 234] <- 1-foi[j]
    #mat1[331,234] <- foi[j]
    mat1[332,234] <- foi[j]
    
    #quaternary infections -  S415
    mat1[235, 235] <- 1-foi[j]
    #mat1[333,235] <- foi[j]
    mat1[334,235] <- foi[j]
    
    #quaternary infections -  S421
    mat1[236, 236] <- 1-2*foi[j]
    mat1[335,236] <- foi[j]
    mat1[336,236] <- foi[j]
    
    #quaternary infections -  S423
    mat1[237, 237] <- 1-2*foi[j]
    mat1[337,237] <- foi[j]
    mat1[338,237] <- foi[j]
    
    
    #quaternary infections -  S425
    mat1[238, 238] <- 1-2*foi[j]
    mat1[339,238] <- foi[j]
    mat1[340,238] <- foi[j]
    
    
    #quaternary infections -  S431
    mat1[239, 239] <- 1-foi[j]
    #mat1[341,239] <- foi[j]
    mat1[342,239] <- foi[j]
    
    #quaternary infections -  S432
    mat1[240, 240] <- 1-2*foi[j]
    mat1[343,240] <- foi[j]
    mat1[344,240] <- foi[j]
    
    #quaternary infections -  S435
    mat1[241, 241] <- 1-foi[j]
    mat1[345,241] <- foi[j]
    #mat1[346,241] <- foi[j]
    
    #quaternary infections -  S451
    mat1[242, 242] <- 1-foi[j]
    #mat1[347,242] <- foi[j]
    mat1[348,242] <- foi[j]
    
    #quaternary infections -  S452
    mat1[243, 243] <- 1-2*foi[j]
    mat1[349,243] <- foi[j]
    mat1[350,243] <- foi[j]
    
    #quaternary infections -  S453
    mat1[244, 244] <- 1-foi[j]
    mat1[351,244] <- foi[j]
    #mat1[352,244] <- foi[j]
    
    #quaternary infections -  S512
    mat1[245, 245] <- 1-2*foi[j]
    mat1[353,245] <- foi[j]
    mat1[354,245] <- foi[j]
    
    #quaternary infections -  S513
    mat1[246, 246] <- 1-foi[j]
    #mat1[355,246] <- foi[j]
    mat1[356,246] <- foi[j]
    
    #quaternary infections -  S514
    mat1[247, 247] <- 1-foi[j]
    #mat1[357,247] <- foi[j]
    mat1[358,247] <- foi[j]
    
    #quaternary infections -  S521
    mat1[248, 248] <- 1-2*foi[j]
    mat1[359,248] <- foi[j]
    mat1[360,248] <- foi[j]
    
    #quaternary infections -  S523
    mat1[249, 249] <- 1-2*foi[j]
    mat1[361,249] <- foi[j]
    mat1[362,249] <- foi[j]
    
    #quaternary infections -  S524
    mat1[250, 250] <- 1-2*foi[j]
    mat1[363,250] <- foi[j]
    mat1[364,250] <- foi[j]
    
    #quaternary infections -  S531
    mat1[251, 251] <- 1-foi[j]
    #mat1[365,251] <- foi[j]
    mat1[366,251] <- foi[j]
    
    #quaternary infections -  S532
    mat1[252, 252] <- 1-2*foi[j]
    mat1[367,252] <- foi[j]
    mat1[368,252] <- foi[j]
    
    #quaternary infections -  S534
    mat1[253, 253] <- 1-foi[j]
    mat1[369,253] <- foi[j]
    #mat1[370,253] <- foi[j]
    
    #quaternary infections -  S541
    mat1[254, 254] <- 1-foi[j]
    #mat1[371,254] <- foi[j]
    mat1[372,254] <- foi[j]
    
    #quaternary infections -  S542
    mat1[255, 255] <- 1-2*foi[j]
    mat1[373,255] <- foi[j]
    mat1[374,255] <- foi[j]
    
    #quaternary infections -  S543
    mat1[256, 256] <- 1-foi[j]
    mat1[375,256] <- foi[j]
    #mat1[376,256] <- foi[j]
    
    
    #recovery from quaternary infections into shorttrm heterotypic immunity
    
    for(i in 0:119){
      #recovery from quaternary infections into shorttrm heterotypic immunity
      mat1[(257+i), (257+i)] <- 1-recov
      mat1[(377+i), (257+i)] <- recov 
      
      #waning out of heterotypic immunity
      mat1[(377+i), (377+i)] <- 1- wane_hetero[j]  
      mat1[(497+i), (377+i)] <- wane_hetero[j]  
      
      #5th order infections - only one option for each - need to eliminate the ones getting infected by serotype 2 here
      
      mat1[(497+i), (497+i)] <- 1- foi[j]  
      mat1[(617+i), (497+i)] <- foi[j]  
      
      #recovery from 5th-order infections into longterm heterotypic immunity to all serotypes (Pm = 736) 
      mat1[(617+i), (617+i)] <-  1-recov
      mat1[(737), (617+i)] <- recov
      
    }
    
    #then, for those 5th order infections ending in serotype 2, replace with 0
    mat1[506,506] <- 1
    mat1[626,506] <- 0
    
    mat1[508,508] <- 1
    mat1[628,508] <- 0
    
    
    mat1[512,512] <- 1
    mat1[632,512] <- 0
    
    
    mat1[514,514] <- 1
    mat1[634,514] <- 0
    
    
    mat1[518,518] <- 1
    mat1[638,518] <- 0
    
    mat1[520,520] <- 1
    mat1[640,520] <- 0
    
    mat1[550,550] <- 1
    mat1[670,550] <- 0
    
    mat1[558,558] <- 1
    mat1[678,558] <- 0
    
    mat1[561,561] <- 1
    mat1[681,561] <- 0
    
    mat1[564,564] <- 1
    mat1[684,564] <- 0
    
    mat1[567,567] <- 1
    mat1[687,567] <- 0
    
    mat1[572,572] <- 1
    mat1[692,572] <- 0
    
    mat1[574,574] <- 1
    mat1[694,574] <- 0
    
    mat1[583,583] <- 1
    mat1[703,583] <- 0
    
    mat1[585,585] <- 1
    mat1[705,585] <- 0
    
    mat1[588,588] <- 1
    mat1[708,588] <- 0
    
    mat1[591,591] <- 1
    mat1[711,581] <- 0
    
    mat1[596,596] <- 1
    mat1[716,596] <- 0
    
    mat1[598,598] <- 1
    mat1[718,598] <- 0
    
    mat1[606,606] <- 1
    mat1[726,606] <- 0
    
    mat1[609,609] <- 1
    mat1[729,609] <- 0
    
    mat1[612,612] <- 1
    mat1[732,612] <- 0
    
    mat1[615,615] <- 1
    mat1[735,615] <- 0
    
    
    
    #once you reach immunity, you stay
    mat1[737, 737] <- 1
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*c+1):((j-1)*c+c),((j-1)*c+1):((j-1)*c+c)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
  return(Tmat)
  
}
transform.vect <- function(vec, s, c){
  vec2 <- t(commutation.matrix(r=s,c=c))%*%vec
  return(vec2)
}
find.biweek = function(t, times){
  
  biwks <- sort(unique(round(revtrunc(times),4)))
  biwks <- c(biwks[2:length(biwks)], biwks[1])
  this.wk = round(revtrunc(times[t]), 4)
  this.biwk <- which(biwks==this.wk)
  
  
  return(this.biwk)
  
  
}
get.age.struct = function(pop, s){ #function that collapses population vector in disease state form to age structure
  age.mat <- mat_split(M=matrix(pop, ncol=1), r=s, c=1)
  age.mat.list <- c()
  for (i in 1:dim(age.mat)[3]){
    age.mat.list[[i]] <- age.mat[,,i]
  }
  age.mat <- Reduce('+', age.mat.list)
  return(age.mat)
}
get.age.struct.M = function(pop, c){ #function that collapses population vector in disease state form to age structure
  age.mat <- mat_split(M=matrix(pop, ncol=1), r=c, c=1)
  age.mat.list <- c()
  for (i in 1:dim(age.mat)[3]){
    age.mat.list[[i]] <- age.mat[,,i]
  }
  #zage.mat <- Reduce('+', age.mat.list)
  age.mat <- lapply(age.mat.list, sum)
  age.mat.dat <- c( c(1:length(age.mat)), unlist(age.mat))
  names(age.mat.dat) <- c("age", "pop")
  return(age.mat)
}
prev.by.age <- function(dat){
  dat$n_age = sum(dat$count)
  dat$prev = dat$count/dat$n_age
  #and seroprev
  dat$seropos = sum(dat$count[dat$class=="M"], dat$count[dat$class=="R"])
  dat$seroprev = dat$seropos/dat$n_age
  return(dat)
}
revtrunc = function(x){
  newx = x - floor(x)
  return(newx)
}
spline.fit = function(data){
  #replace any Inf vals only the okay values
  #data <- data[complete.cases(data),]
  data$seroprevalence[data$seroprevalence==Inf] <- 1
  data$seroprevalence[data$seroprevalence<0] <- 0
  # if(length(data$age)>1){
  spline1 = with(data, smooth.spline(x=age, y=seroprevalence)) 
  #  return(spline1)
  # }else{
  #  return(NA)
  #  }
}
pred.fit <- function(spline1, ages.new){
  prev.comp = predict(spline1, x=ages.new)
  return(prev.comp$y)
}
get.seas.seroprev = function(dat){
  #get doy in bat calendar
  #convert to biweek
  #calc seroprev by biweek
  
  dat$doy <- yday(dat$date)
  #now correct for birthday of bat in questiob
  if(unique(dat$species=="Pteropus rufus")){
    bday <- yday("2015-10-01") #320
  } else if(unique(dat$species=="Eidolon dupreanum")){
    bday = yday("2015-11-01")
  }
  #write over for doy
  dat$new_doy = NA
  for (i in 1:length(dat$doy)){
    if(dat$doy[i] > bday){
      dat$new_doy[i] <- dat$doy[i] - bday
    } else if (dat$doy[i] <= bday){
      dat$new_doy[i] <- dat$doy[i] + (365-bday)
    }
  }
  
  dat$doy <- dat$new_doy
  dat <- dplyr::select(dat, -(new_doy))
  
  #now convert to biweek
  dat$biwk <- NA
  brk <- seq(0,365, by=14)
  for (i in 1:length(dat$biwk)){
    tmp <- brk[dat$doy[i] >=  brk]
    tmp2 <- tmp[length(tmp)] 
    dat$biwk[i] <- which(brk==tmp2)
  }
  
  #now calc seroprev by biweek
  biwk.sero <- ddply(dat, .(biwk), summarize, seropos = sum(prev), sero_lci = sum(prev_lci), sero_uci = sum(prev_uci), n=length(prev))  
  biwk.sero$seroprev = biwk.sero$seropos/biwk.sero$n
  biwk.sero$seroprev_lci = biwk.sero$sero_lci/biwk.sero$n
  biwk.sero$seroprev_uci = biwk.sero$sero_uci/biwk.sero$n
  biwk.sero$biwk = as.numeric(biwk.sero$biwk)
  
  return(biwk.sero)
}
get.mod.seas = function(mod.out){
  #convert time to doy
  mod.out$doy = mod.out$time*365
  #then to biwk
  brk <- seq(0,365, by=14) #doy breaks by biweek
  brk <- brk[-length(brk)]
  #brk.seq = 1:length(brk)
  mod.out$biwk <- NA
  for (i in 1:length(mod.out$biwk)){
    tmp <- brk[mod.out$doy[i] >=  brk]
    tmp2 <- tmp[length(tmp)] 
    mod.out$biwk[i] <- which(brk==tmp2)
  }
  #then return
  return(mod.out)
}
cum.cases = function(df){
  df$cumcases = cumsum(df$count)
  return(df)
}
age.sero.plot = function(dat){
  with(dat, plot(age, seroprevalence, type = "b", ylim=c(0,1)))
}
get.seroprev.dat = function(data, vis_split, cutoff){
  
  visbin = seq(3,floor(max(data$age)), vis_split)
  #but breakdown the early years into more
  visbin = c(c(0,.5, 1, 2),   visbin)
  data$age_year <- NA
  
  for (i in 1:length(data$age_year)){
    tmp = visbin[data$age[i] > visbin]
    data$age_year[i] <- tmp[length(tmp)] 
  }
  
  if(cutoff=="mean"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev)/length(prev), count=length(prev))  
  }else if (cutoff=="uci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_uci)/length(prev_uci), count=length(prev_uci)) 
  }else if (cutoff=="lci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_lci)/length(prev_lci), count=length(prev_lci)) 
  }
  
  #you want the midpoint between either end of the age class now
  vect.age.yr = sort(unique(data$age_year))
  dat.sum2$age_plot = NA
  
  for (i in 1:(length(dat.sum2$age_year)-1)){
    dat.sum2$age_plot [i] = ((dat.sum2$age_year[i + 1] - dat.sum2$age_year[i])/2) + dat.sum2$age_year[i]
  }
  dat.sum2$age_plot[length(dat.sum2$age_plot)] =  (ceiling(max(data$age)) - dat.sum2$age_year[length(dat.sum2$age_year)])/2  + dat.sum2$age_year[length(dat.sum2$age_year)]
  
  names(dat.sum2)  <- c("real_age", "prevalence", "count", "age_year") #<- names(dat.sum2.tmp)
  
  # dat.sum2 <-  rbind(dat.sum2.tmp, dat.sum2)
  
  dat.sum3 = dat.sum2
  dat.sum3$class = "seropositive"
  rownames(dat.sum3) <- c()
  
  return(dat.sum3)
  
}
build.pop.mat = function(surv, surv_juv, s, adult_fec_vector){
  pop.mat = matrix(0,  nrow=(s-1), ncol = (s-1))
  diag(pop.mat) = surv
  diag(pop.mat)[1] = surv_juv
  col_s = c(rep(0, s-2), surv)
  pop.mat = cbind(pop.mat, col_s)
  row1=adult_fec_vector*surv #birth rates* survival rates for the entire age distribution
  #row1 = c(0, rep((adult_fec*surv), (s-1))) #bats reproduce for the first time at the end of the second year of life. good for E. dup and P. ruf
  pop.mat = rbind(row1,pop.mat)
  return(pop.mat)
  
}#for stable age distribution
mat_split <- function(M, r, c){
  nr <- ceiling(nrow(M)/r)
  nc <- ceiling(ncol(M)/c)
  newM <- matrix(NA, nr*r, nc*c)
  newM[1:nrow(M), 1:ncol(M)] <- M
  
  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats)<-c(r, c, N)
  return(mats)
}#for matrix slicing. give it the number of rows and columns you want in the resulting matrices
stack.age = function(dat, s){
  dat.new= list()
  for (i in 1:s){
    dat.new[[i]] = dat[i,]
  }
  dat.new = c(unlist(dat.new))
  return(dat.new)
}
stack.class = function(mat.dat, c){
  
  new.dat = list()
  for (i in 1:c){
    new.dat[[i]] = mat.dat[i,]
  }
  new.dat = c(unlist(new.dat))
  return(new.dat)
}
sum.yr.all <- function(df){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(df$age))
  
  df.out <- merge(df.out, df.sum, by="age", all.x = T, sort = F)
  df.out$Nage[is.na(df.out$Nage)] <- 0
  df.out$year[is.na(df.out$year)] <- unique(df$year)
  #df.out <- rbind(c(0,0), df.out)
  #bind
  df.add <- cbind.data.frame(age=0, year = unique(df$year), Nage=0)
  df.out <- rbind(df.add, df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.sum$Nage)
  df.out$cum_prop_cases <- df.out$cum_cases/df.out$n
  
  return(df.out)
  
}
divide <- function(df, ntyr1){
  dfnew <- (df/ntyr1)
  return(dfnew)
  
}
divide.rate <- function(df, ntyr1){
  dfnew = (1-exp(-(df/ntyr1)))
  
  return(dfnew)
  
}
replicate.data <- function(df, slim.quant){
  #print(unique(df$year))
  if(df$count>0){
    new.dat = cbind.data.frame(age=rep(df$age,(df$count)), case=rep(1, (df$count)))  
    new.dat$year <- unique(df$year)
    
    #then, cut to 5% of the cases:
    #all should be the same, so just take the top 5% of rows
    n.row= round(nrow(new.dat)*(slim.quant),0)
    
    new.dat <- new.dat[1:n.row,]
    
    
    return(new.dat)
  }
  
  
}
replicate.data.type <- function(df, slim.quant){
  #print(unique(df$year))
  if(df$count>0){
    new.dat = cbind.data.frame(age=rep(df$age,(df$count)), case=rep(1, (df$count)))  
    new.dat$year <- unique(df$year)
    new.dat$state <- unique(df$state)
    
    #then, cut to 5% of the cases:
    #all should be the same, so just take the top 5% of rows
    n.row= round(nrow(new.dat)*(slim.quant),0)
    
    new.dat <- new.dat[1:n.row,]
    
    
    return(new.dat)
  }
  
  
}
mean.age <- function(df){
  df$mult <- df$age*df$count
  mean.age <- sum(df$mult)/sum(df$count)
  df2 <- cbind.data.frame(year=unique(df$year), mean_age=mean.age, state = unique(df$state))
  return(df2)
}
plot.ages.all <- function(dat, save.plot, view.plot, filename, slim.quant, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #select only those viewed as "cases"
  
  if(length(unique(dat$age))>1){
    dat2 = subset(dat1, age<max(dat$age))  
  }
  
  # 
  #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
  df.sum = ddply(dat2,.(year,age),summarise, count=sum(count))
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  
  
  
  
  #split by a year
  df.year <- dlply(df.sum,.(year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year
  df.age <- dlply(df.sum,.(year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data, slim.quant=slim.quant))
  
  
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  
  #and plot
  p1 <- ggplot(dat.age) + 
    geom_jitter(aes(x=year, y=age), width=.09, height=.09, size=.09, alpha=.8, show.legend = F) +
    geom_violin(aes(x=year,y=age, group=year),  color="gray55", draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") +
    geom_hline(aes(yintercept=1), color="red") + coord_cartesian(ylim=c(0,max(dat$age)))#,xlim=c(2015,2020))
  
  
  if(save.plot==TRUE){
    ggsave(file = filename,
           plot= p1,
           units="mm",  
           width=80, 
           height=55, 
           scale=3, 
           dpi=300)
    
  }
  if(view.plot==TRUE){
    print(p1)
  }
  
  
  
  return(mean.df)
}
process.plot.age.cum.triple.type <-function(dat, count.type, perc.obs, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #select only those viewed as "cases"
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32" | class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  denv.case$type = "secondary"
  
  if(count.type=="tertiary"){
    
    denv.case.tert = subset(dat1, class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
    denv.case.tert$type <- "tertiary"  
    denv.case.tert$count <- denv.case.tert$count*perc.obs
    denv.case <- rbind( denv.case,  denv.case.tert)
    
  }
  
  
  
  denv.split <- dlply(denv.case, .(year))
  
  
  denv.cum <- lapply(denv.split, cum.sum.year)
  denv.dat <- data.table::rbindlist(denv.cum)
  denv.dat$year <- as.factor(denv.dat$year)
  
  #head(denv.dat)
  
  p1 <- ggplot(data=denv.dat) + geom_line(aes(x=age, y=cum_prop_cases, color=year)) +
    scale_color_viridis_d(direction=-1,option="turbo")# + facet_wrap(~year) 
  print(p1)
  
  return(denv.dat)
}
check.equil <- function(dat){
  dat.ts <- ddply(dat, .(time, class), summarise, count=sum(count))
  #and get total by time
  dat.N  <- ddply(dat, .(time), summarise, Ntot=sum(count))
  
  dat.ts <- merge(dat.ts, dat.N, by="time")
  dat.ts$proportion <- dat.ts$count/dat.ts$Ntot
  
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=proportion, color=class)) #+ coord_cartesian(ylim=c(0,.1))
  print(p1)
  return(dat.ts)
}
check.equil.cases.only <- function(dat){
  #and get total by time
  denv.case = subset(dat,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32" | class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
  denv.case$serotype <- NA
  denv.case$serotype[denv.case$class=="I12" | denv.case$class=="I32" | denv.case$class=="I132" | denv.case$class=="I312" ] <- "I2"
  denv.case$serotype[denv.case$class=="I13" | denv.case$class=="I23" | denv.case$class=="I123" | denv.case$class=="I213"] <- "I3"
  denv.case$serotype[denv.case$class=="I31" | denv.case$class=="I21"| denv.case$class=="I321" | denv.case$class=="I231"] <- "I1"
  
  #dat.ts <- ddply(denv.case, .(time, class), summarise, count=sum(count))
  dat.ser <- ddply(denv.case, .(time, serotype), summarise, count=sum(count))
  
  
  dat.N  <- ddply(denv.case, .(time), summarise, Ntot=sum(count))
  
  #dat.ts <- merge(dat.ts, dat.N, by="time")
  dat.ser <- merge(dat.ser, dat.N, by="time")
  
  #dat.ts$proportion <- dat.ts$count/dat.ts$Ntot
  dat.ser$proportion <- dat.ser$count/dat.ser$Ntot
  #dat.ts$proportion[is.na(dat.ts$proportion)] <- 0
  dat.ser$proportion[is.na(dat.ser$proportion)] <- 0
  
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=proportion, color=class)) #+ coord_cartesian(ylim=c(0,.1))
  p1 <- ggplot(dat.ser) + geom_line(aes(x=time, y=proportion, color=serotype)) +facet_grid(serotype~.) + coord_cartesian(ylim=c(0,1))
  print(p1)
  return(dat.ser)
}
cum.sum.year <- function(df){
  df.sum <- ddply(df,.(age), summarise, cases=sum(count))
  
  df.sum$cum_cases = cumsum(df.sum$cases)
  df.sum$cum_prop_cases <- df.sum$cum_cases/sum(df.sum$cases)
  df.sum$year <- unique(df$year)
  return(df.sum)
}
process.all <- function(df){
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "S", "I1","P1","I2", "P2","I3","I4", "PM")
  dat.tot$N = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  #par(mfrow = c(1,1))
  #with(dat.tot, plot(time, N, type="l")) #ylim =c(0,1.2*max(N))))
  #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  
  # prop.tot = dat.tot
  # 
  # prop.tot$S = prop.tot$S/prop.tot$N
  # prop.tot$I = prop.tot$Ii/prop.tot$N
  # prop.tot$R = prop.tot$R/prop.tot$N
  # 
  
  
  
  
  dat.melt = melt(dat.tot,id.vars = "time")
  
  #head(dat.melt)
  
  #take after burnin if we assume this is a virus at equilibrium
  dat.melt <- subset(dat.melt, time>length(burnin_lambda))
  dat.melt <- subset(dat.melt, !is.na(value))
  
  #dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5])))
  names(dat.melt) = c("time", "class", "count")
  dat.N = subset(dat.melt, class =="N")
  dat.melt = subset(dat.melt, class !="N")
  dat.N <- dplyr::select(dat.N, -(class))
  names(dat.N)[names(dat.N)=="count"] <- "Ntot"
  
  dat.melt <- merge(dat.melt, dat.N, by="time")
  #head(dat.melt)
  #dat.tot$class = rep(c("M", "S", "I", "R"), each= length(times))
  dat.melt$class = factor(dat.melt$class, levels=c("S", "I1","P1","I2", "P2","I3","I4", "PM"))
  dat.melt$proportion <- dat.melt$count/dat.melt$Ntot
  
  dat.melt$time <- dat.melt$time - length(burnin_lambda)
  # colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
  #ggplot(data=dat.melt) + geom_line(aes(x=time, y=count, color=class)) #+ scale_color_manual(values=colz)
  
  #dat.melt$proportion = dat.melt$count/sim_pop
  #ggplot(data=dat.melt) + geom_line(aes(x=time, y=proportion, color=class)) #+ scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  #now, get the age-structured incidence - incidence is "I3 + I4"
  
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("S", "I1","P1","I2", "P2","I3","I4", "PM"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  
  #plot age structured cumulative incidence
  age.incidence = subset(age.dat.tot, class=="I3" | class == "I4")
  
  #and collect everyone else
  age.non = subset(age.dat.tot, class!="I3" & class!="I4")
  #head(age.incidence )
  
  age.incidence$year = trunc(age.incidence$times)
  age.non$year = trunc(age.non$times)
  age.incidence <- arrange(age.incidence, year, age)
  
  age.non <- arrange(age.non, year, age)
  
  names(age.non) <- c("times", "count_non_infectious", "class", "age", "year")
  
  age.non <- dplyr::select(age.non, -(class), -(times))
  
  age.non <- ddply(age.non, .(year, age), summarise, count_non_infectious = sum(count_non_infectious))
  
  
  age.sum <- ddply(age.incidence, .(year, age), summarise, count = sum(count))
  
  age.sum <- merge(age.sum, age.non, by=c("year", "age"))
  
  #and just return this counts of years, ages, and those infected and not
  
  
  return(age.sum)#, seas.prev))
  
}
add.year <- function(dat, add.on){
  dat$year <- trunc(dat$time) + add.on
  return(dat)
}
sim.SIR.age.four.clim.maintain <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, rate_hetero, clim.vect, births, pop_vector, sigma){
  
  #check on climate vector - should be length of steps within a year (often 26)
  #(this could be modulated to be the length of the time series if you wanted)
  if(length(clim.vect)==1){
    clim.vect = rep(clim.vect, ntyr)
  }
  
  
  #first, if mort, births, sigma, or foi are shorter than the time series,
  #make them match in length here.
  
  #this will also repeat if they come in as a vector (age-structured)
  
  #make as a list
  if(length(foi)<yrs){
    foi <- rep(list(foi), yrs)
  }else{
    foi <- as.list(foi)
  }
  
  if(length(sigma)< yrs){
    sigma = rep(list(sigma), yrs)
  }else{
    sigma <- as.list(sigma)
  }
  
  # if(length(mort)<yrs){
  #   mort = rep(list(mort), yrs)
  # }else{
  #   mort = as.list(mort)
  # }
  # 
  if(length(births)<yrs){
    births = rep(list(births), yrs)
  }else{
    births = as.list(births)
  }
  
  # 
  if(length(rate_hetero)<yrs){
  rate_hetero = rep(list(rate_hetero), yrs)
  }else{
    rate_hetero = as.list(rate_hetero)
  }
  
  
  #
  #number of epidemic classes
  c=147
  
  
  #length of time series
  times <-   seq(0, yrs, by =1/ntyr)
  times <-   times + (year.end-yrs)
  
  #split the births up by biweek too
  birth_vector_biwk = as.list(c(unlist(lapply(lapply(births, divide, ntyr1 =ntyr), rep, ntyr))))
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  #import the stable age structure from the literature
  #mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  #stab.struct = Re(eigen(mat1)$vector[,1])
  #stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  #lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  pop.mat = pop_vector #number of people per age class
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  #structure the population for 3 strains even though we only have 2 to start. they will get added later
  
  PM_init = rep(0, s)
  
  PMs_init = rep(0, s)
  
  I4321_init = rep(0, s)
  I4312_init = rep(0, s)
  I4231_init = rep(0, s)
  I4213_init = rep(0, s)
  I4132_init = rep(0, s)
  I4123_init = rep(0, s)
  
  I3421_init = rep(0, s)
  I3412_init = rep(0, s)
  I3241_init = rep(0, s)
  I3214_init = rep(0, s)
  I3142_init = rep(0, s)
  I3124_init = rep(0, s)
  
  I2431_init = rep(0, s)
  I2413_init = rep(0, s)
  I2341_init = rep(0, s)
  I2314_init = rep(0, s)
  I2143_init = rep(0, s)
  I2134_init = rep(0, s)
  
  I1432_init = rep(0, s)
  I1423_init = rep(0, s)
  I1342_init = rep(0, s)
  I1324_init = rep(0, s)
  I1243_init = rep(0, s)
  I1234_init = rep(0, s)
  
  
  S432_init = rep(0, s)
  S431_init = rep(0, s)
  S423_init = rep(0, s)
  S421_init = rep(0, s)
  S413_init = rep(0, s)
  S412_init = rep(0, s)
  
  S342_init = rep(0, s)
  S341_init = rep(0, s)
  S324_init = rep(0, s)
  S321_init = rep(0, s)
  S314_init = rep(0, s)
  S312_init = rep(0, s)
  
  S243_init = rep(0, s)
  S241_init = rep(0, s)
  S234_init = rep(0, s)
  S231_init = rep(0, s)
  S214_init = rep(0, s)
  S213_init = rep(0, s)
  
  S143_init = rep(0, s)
  S142_init = rep(0, s)
  S134_init = rep(0, s)
  S132_init = rep(0, s)
  S124_init = rep(0, s)
  S123_init = rep(0, s)
  
  
  P432_init = rep(0, s)
  P431_init = rep(0, s)
  P423_init = rep(0, s)
  P421_init = rep(0, s)
  P413_init = rep(0, s)
  P412_init = rep(0, s)
  
  P342_init = rep(0, s)
  P341_init = rep(0, s)
  P324_init = rep(0, s)
  P321_init = rep(0, s)
  P314_init = rep(0, s)
  P312_init = rep(0, s)
  
  P243_init = rep(0, s)
  P241_init = rep(0, s)
  P234_init = rep(0, s)
  P231_init = rep(0, s)
  P214_init = rep(0, s)
  P213_init = rep(0, s)
  
  P143_init = rep(0, s)
  P142_init = rep(0, s)
  P134_init = rep(0, s)
  P132_init = rep(0, s)
  P124_init = rep(0, s)
  P123_init = rep(0, s)
  
  
  I432_init = rep(0, s)
  I431_init = rep(0, s)
  I423_init = rep(0, s)
  I421_init = rep(0, s)
  I413_init = rep(0, s)
  I412_init = rep(0, s)
  
  I342_init = rep(0, s)
  I341_init = rep(0, s)
  I324_init = rep(0, s)
  I321_init = rep(0, s)
  I314_init = rep(0, s)
  I312_init = rep(0, s)
  
  I243_init = rep(0, s)
  I241_init = rep(0, s)
  I234_init = rep(0, s)
  I231_init = rep(0, s)
  I214_init = rep(0, s)
  I213_init = rep(0, s)
  
  I143_init = rep(0, s)
  I142_init = rep(0, s)
  I134_init = rep(0, s)
  I132_init = rep(0, s)
  I124_init = rep(0, s)
  I123_init = rep(0, s)
  
  S43_init = rep(0, s)
  S42_init = rep(0, s)
  S41_init = rep(0, s)
  S34_init = rep(0, s)
  S32_init = rep(0, s)
  S31_init = rep(0, s)
  S24_init = rep(0, s)
  S23_init = rep(0, s)
  S21_init = rep(0, s)
  S14_init = rep(0, s)
  S13_init = rep(0, s)
  S12_init = rep(0, s)
  
  P43_init = rep(0, s)
  P42_init = rep(0, s)
  P41_init = rep(0, s)
  P34_init = rep(0, s)
  P32_init = rep(0, s)
  P31_init = rep(0, s)
  P24_init = rep(0, s)
  P23_init = rep(0, s)
  P21_init = rep(0, s)
  P14_init = rep(0, s)
  P13_init = rep(0, s)
  P12_init = rep(0, s)
  
  I43_init = rep(0, s)
  I42_init = rep(0, s)
  I41_init = rep(0, s)
  I34_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I24_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I14_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  
  S4_init = rep(0, s)
  S3_init = rep(0, s)
  S2_init = rep(0, s)
  S1_init = rep(0, s)
  
  P4_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  
  I4_init = rep(0, s); I4_init[1] = 5 #comment out if you just want to check demography
  I3_init = rep(0, s); I3_init[1] = 5 #comment out if you just want to check demography
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  
  S_init = pop.mat - I1_init - I2_init - I3_init - I4_init
  
  
  N_tot = cbind(S_init, 
                I1_init, I2_init,I3_init,I4_init,
                P1_init, P2_init,P3_init,P4_init, 
                S1_init,S2_init,S3_init,S4_init,
                I12_init, I13_init, I14_init, 
                I21_init, I23_init, I34_init, 
                I31_init, I32_init, I34_init,
                I41_init, I42_init, I43_init,
                P12_init, P13_init, P14_init, 
                P21_init, P23_init, P34_init, 
                P31_init, P32_init, P34_init,
                P41_init, P42_init, P43_init,
                S12_init, S13_init, S14_init, 
                S21_init, S23_init, S34_init, 
                S31_init, S32_init, S34_init,
                S41_init, S42_init, S43_init,
                I143_init,I142_init,I134_init,I132_init,I124_init,I123_init,
                I213_init,I214_init, I231_init, I234_init, I241_init, I243_init,
                I312_init, I314_init, I321_init, I324_init, I341_init, I342_init,
                I412_init, I413_init, I421_init, I423_init,I431_init, I432_init, 
                P143_init,P142_init,P134_init,P132_init,P124_init,P123_init,
                P213_init,P214_init, P231_init, P234_init, P241_init, P243_init,
                P312_init, P314_init, P321_init, P324_init, P341_init, P342_init,
                P412_init, P413_init, P421_init, P423_init,P431_init, P432_init, 
                S143_init,S142_init,S134_init,S132_init,S124_init,S123_init,
                S213_init,S214_init, S231_init, S234_init, S241_init, S243_init,
                S312_init, S314_init, S321_init, S324_init, S341_init, S342_init,
                S412_init, S413_init, S421_init, S423_init,S431_init, S432_init,  
                I1234_init, I1243_init, I1324_init, I1342_init,I1423_init,I1432_init,
                I2134_init, I2143_init, I2314_init, I2341_init,I2413_init, I2431_init, 
                I3124_init, I3142_init, I3214_init, I3241_init, I3412_init, I3421_init, 
                I4123_init, I4132_init, I4213_init, I4231_init, I4312_init, I4321_init,
                PMs_init,
                PM_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  
  N_pop_ts[,1] <- M_pop # by age class
  
  #and fill in for births <- every epidemic class can give birth
  #transform birth vector into one distributed by age and epidemic class
  
  
  # stab.struct <- get.age.struct.M(M_pop, c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  
  # foi gets spread across the biweeks as well
  # (it comes in as 1 value per year)
  # take your Muench-estimated rates of annual infection and convert
  # to biweekly probabilities of infection
  foi.biweek = as.list(c(unlist(lapply(lapply(foi, divide.rate, ntyr1 =ntyr), rep, ntyr)))) #here, it could easily be modulated to be seasonal - (including forced by precip/temp)
  
  
  #and make annual waning immunity rate into probability of waning immunity per biweek
  sigma.biweek = as.list(c(unlist(lapply(lapply(sigma, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #and the probability of waning heterotypic immunity per biweek
  hetero.biweek = as.list(c(unlist(lapply(lapply(rate_hetero, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(hetero.biweek)!=(length(times)-1) |length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
    print("input vectors of unequal length")
  }
  #and your aging probability to biweeks (aging probability per year is 100%)
  #age.biwk = rep(1-exp(-(rep(1, yrs)/ntyr)), each=26)
  
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 2:(length(times))){
    
    biwk1 <- find.biweek(t=i, times=times)
    clim.mod <- clim.vect[biwk1]
    
    #and select the appropriate age vector, based on the year
    age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
    
    
    
    if(nrow(age.mult.sub)>11){
      print("error with age multiplier sub-selection")
    }
    
    age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
    age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
    age.mult.sub$dur[age.mult.sub$dur<0] <- 0
    age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
    
    
    #and the mortality rates
    if(times[i]<min(mort$year)){
      mort.df = subset(mort, year==min(year))
    }else{
      mort.df = subset(mort, year==trunc(times[i]))
    }
    #make annual mortality rate by age into probability of mortality per biweek
    #these are deaths per 1000 people.
    #we multiply by the number in the age classes, so need to divide by 1000 here
    mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
    mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
    
    
      
      #then, let her go from there!
      
      Tmat <- buildTMat_four_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }
    
  
  
  
  #the first column is the initial conditions. everything after that represents transitions within the year
  #should take the ceiling of each biweek to total within the year
  
  
  # stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  if(s>1){
    N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  }
  
  N.sim.dat <- as.data.frame(N_pop_ts)
  
  names(N.sim.dat) <- times
  
  N.sim.dat$class <- rep(c("S", "I1","I2", "I3", "I4",  "P1", "P2", "P3", "P4", "S1", "S2",  "S3", "S4",
                           "I12", "I13", "I14","I21", "I23", "I24", "I31", "I32", "I34", "I41", "I42", "I43",
                           "P12", "P13", "P14","P21", "P23", "P24", "P31", "P32", "P34",  "P41", "P42", "P43",
                           "S12", "S13", "S14","S21", "S23", "S24", "S31", "S32", "S34", "S41", "S42", "S43",
                           "I123","I124", "I132", "I134", "I142", "I143",
                           "I213", "I214","I231","I234", "I241", "I243",
                           "I312","I314", "I321", "I324", "I341", "I342",
                           "I412", "I413", "I421", "I423", "I431", "I432",
                           "P123","P124", "P132", "P134", "P142", "P143",
                           "P213", "P214","P231","P234", "P241", "P243",
                           "P312","P314", "P321", "P324", "P341", "P342",
                           "P412", "P413", "P421", "P423", "P431", "P432",
                           "S123","S124", "S132", "S134", "S142", "S143",
                           "S213", "S214","S231","S234", "S241", "S243",
                           "S312","S314", "S321", "S324", "S341", "S342",
                           "S412", "S413", "S421", "S423", "S431", "S432",
                           "I1234","I1243", "I1324", "I1342", "I1423", "I1432",
                           "I2134", "I2143","I2314","I2341", "I2413", "I2431",
                           "I3124","I3142", "I3214", "I3241", "I3412", "I3421",
                           "I4123", "I4132", "I4213", "I4231", "I4312", "I4321",
                           "Pms", "Pm"), 
                         each = s)
  N.sim.dat$age <- rep(1:s, length(unique(N.sim.dat$class)))
  
  #and melt
  N.sim.df <- melt(N.sim.dat, id=c("class", "age"), value.name = "count", variable.name = "time")
  
  N.sim.df$time <- as.numeric(as.character(N.sim.df$time))
  N.sim.df$count <- as.numeric(as.character(N.sim.df$count))
  
  #and add year
  N.sim.df$year <- trunc(N.sim.df$time)
  #N.sim.df$time <- N.sim.df$time-1
  #N.sim.df$year <- N.sim.df$year-1
  
  
  
  return( N.sim.df)
}
sim.SIR.age.four.clim.wane <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, rate_hetero, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro, prop_imm){
  
  
  
  
  
  #check on climate vector - should be length of steps within a year (often 26)
  #(this could be modulated to be the length of the time series if you wanted)
  if(length(clim.vect)==1){
    clim.vect = rep(clim.vect, ntyr)
  }
  
  
  #first, if mort, births, sigma, or foi are shorter than the time series,
  #make them match in length here.
  
  #this will also repeat if they come in as a vector (age-structured)
  
  #make as a list
  if(length(foi)<yrs){
    foi <- rep(list(foi), yrs)
  }else{
    foi <- as.list(foi)
  }
  
  if(length(sigma)< yrs){
    sigma = rep(list(sigma), yrs)
  }else{
    sigma <- as.list(sigma)
  }
  
  # if(length(mort)<yrs){
  #   mort = rep(list(mort), yrs)
  # }else{
  #   mort = as.list(mort)
  # }
  # 
  if(length(births)<yrs){
    births = rep(list(births), yrs)
  }else{
    births = as.list(births)
  }
  
  # 
  if(length(rate_hetero)<yrs){
    rate_hetero = rep(list(rate_hetero), yrs)
  }else{
    rate_hetero = as.list(rate_hetero)
  }
  
  
  #
  #number of epidemic classes
  c=147
  
  
  #length of time series
  times <-   seq(0, yrs, by =1/ntyr)
  times <-   times + (year.end-yrs)
  
  #split the births up by biweek too
  birth_vector_biwk = as.list(c(unlist(lapply(lapply(births, divide, ntyr1 =ntyr), rep, ntyr))))
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  #import the stable age structure from the literature
  #mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  #stab.struct = Re(eigen(mat1)$vector[,1])
  #stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  #lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  pop.mat = pop_vector #number of people per age class
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  #structure the population for 3 strains even though we only have 2 to start. they will get added later
  
  PM_init = rep(0, s)
  
  PMs_init = rep(0, s)
  
  I4321_init = rep(0, s)
  I4312_init = rep(0, s)
  I4231_init = rep(0, s)
  I4213_init = rep(0, s)
  I4132_init = rep(0, s)
  I4123_init = rep(0, s)
  
  I3421_init = rep(0, s)
  I3412_init = rep(0, s)
  I3241_init = rep(0, s)
  I3214_init = rep(0, s)
  I3142_init = rep(0, s)
  I3124_init = rep(0, s)
  
  I2431_init = rep(0, s)
  I2413_init = rep(0, s)
  I2341_init = rep(0, s)
  I2314_init = rep(0, s)
  I2143_init = rep(0, s)
  I2134_init = rep(0, s)
  
  I1432_init = rep(0, s)
  I1423_init = rep(0, s)
  I1342_init = rep(0, s)
  I1324_init = rep(0, s)
  I1243_init = rep(0, s)
  I1234_init = rep(0, s)
  
  
  S432_init = rep(0, s)
  S431_init = rep(0, s)
  S423_init = rep(0, s)
  S421_init = rep(0, s)
  S413_init = rep(0, s)
  S412_init = rep(0, s)
  
  S342_init = rep(0, s)
  S341_init = rep(0, s)
  S324_init = rep(0, s)
  S321_init = rep(0, s)
  S314_init = rep(0, s)
  S312_init = rep(0, s)
  
  S243_init = rep(0, s)
  S241_init = rep(0, s)
  S234_init = rep(0, s)
  S231_init = rep(0, s)
  S214_init = rep(0, s)
  S213_init = rep(0, s)
  
  S143_init = rep(0, s)
  S142_init = rep(0, s)
  S134_init = rep(0, s)
  S132_init = rep(0, s)
  S124_init = rep(0, s)
  S123_init = rep(0, s)
  
  
  P432_init = rep(0, s)
  P431_init = rep(0, s)
  P423_init = rep(0, s)
  P421_init = rep(0, s)
  P413_init = rep(0, s)
  P412_init = rep(0, s)
  
  P342_init = rep(0, s)
  P341_init = rep(0, s)
  P324_init = rep(0, s)
  P321_init = rep(0, s)
  P314_init = rep(0, s)
  P312_init = rep(0, s)
  
  P243_init = rep(0, s)
  P241_init = rep(0, s)
  P234_init = rep(0, s)
  P231_init = rep(0, s)
  P214_init = rep(0, s)
  P213_init = rep(0, s)
  
  P143_init = rep(0, s)
  P142_init = rep(0, s)
  P134_init = rep(0, s)
  P132_init = rep(0, s)
  P124_init = rep(0, s)
  P123_init = rep(0, s)
  
  
  I432_init = rep(0, s)
  I431_init = rep(0, s)
  I423_init = rep(0, s)
  I421_init = rep(0, s)
  I413_init = rep(0, s)
  I412_init = rep(0, s)
  
  I342_init = rep(0, s)
  I341_init = rep(0, s)
  I324_init = rep(0, s)
  I321_init = rep(0, s)
  I314_init = rep(0, s)
  I312_init = rep(0, s)
  
  I243_init = rep(0, s)
  I241_init = rep(0, s)
  I234_init = rep(0, s)
  I231_init = rep(0, s)
  I214_init = rep(0, s)
  I213_init = rep(0, s)
  
  I143_init = rep(0, s)
  I142_init = rep(0, s)
  I134_init = rep(0, s)
  I132_init = rep(0, s)
  I124_init = rep(0, s)
  I123_init = rep(0, s)
  
  S43_init = rep(0, s)
  S42_init = rep(0, s)
  S41_init = rep(0, s)
  S34_init = rep(0, s)
  S32_init = rep(0, s)
  S31_init = rep(0, s)
  S24_init = rep(0, s)
  S23_init = rep(0, s)
  S21_init = rep(0, s)
  S14_init = rep(0, s)
  S13_init = rep(0, s)
  S12_init = rep(0, s)
  
  P43_init = rep(0, s)
  P42_init = rep(0, s)
  P41_init = rep(0, s)
  P34_init = rep(0, s)
  P32_init = rep(0, s)
  P31_init = rep(0, s)
  P24_init = rep(0, s)
  P23_init = rep(0, s)
  P21_init = rep(0, s)
  P14_init = rep(0, s)
  P13_init = rep(0, s)
  P12_init = rep(0, s)
  
  I43_init = rep(0, s)
  I42_init = rep(0, s)
  I41_init = rep(0, s)
  I34_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I24_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I14_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  
  S4_init = rep(0, s)
  S3_init = rep(0, s)
  S2_init = rep(0, s)
  S1_init = rep(0, s)
  
  P4_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  
  I4_init = rep(0, s); I4_init[1] = 5 #comment out if you just want to check demography
  I3_init = rep(0, s); I3_init[1] = 5 #comment out if you just want to check demography
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  
  S_init = pop.mat - I1_init - I2_init - I3_init - I4_init
  
  
  N_tot = cbind(S_init, 
                I1_init, I2_init,I3_init,I4_init,
                P1_init, P2_init,P3_init,P4_init, 
                S1_init,S2_init,S3_init,S4_init,
                I12_init, I13_init, I14_init, 
                I21_init, I23_init, I34_init, 
                I31_init, I32_init, I34_init,
                I41_init, I42_init, I43_init,
                P12_init, P13_init, P14_init, 
                P21_init, P23_init, P34_init, 
                P31_init, P32_init, P34_init,
                P41_init, P42_init, P43_init,
                S12_init, S13_init, S14_init, 
                S21_init, S23_init, S34_init, 
                S31_init, S32_init, S34_init,
                S41_init, S42_init, S43_init,
                I143_init,I142_init,I134_init,I132_init,I124_init,I123_init,
                I213_init,I214_init, I231_init, I234_init, I241_init, I243_init,
                I312_init, I314_init, I321_init, I324_init, I341_init, I342_init,
                I412_init, I413_init, I421_init, I423_init,I431_init, I432_init, 
                P143_init,P142_init,P134_init,P132_init,P124_init,P123_init,
                P213_init,P214_init, P231_init, P234_init, P241_init, P243_init,
                P312_init, P314_init, P321_init, P324_init, P341_init, P342_init,
                P412_init, P413_init, P421_init, P423_init,P431_init, P432_init, 
                S143_init,S142_init,S134_init,S132_init,S124_init,S123_init,
                S213_init,S214_init, S231_init, S234_init, S241_init, S243_init,
                S312_init, S314_init, S321_init, S324_init, S341_init, S342_init,
                S412_init, S413_init, S421_init, S423_init,S431_init, S432_init,  
                I1234_init, I1243_init, I1324_init, I1342_init,I1423_init,I1432_init,
                I2134_init, I2143_init, I2314_init, I2341_init,I2413_init, I2431_init, 
                I3124_init, I3142_init, I3214_init, I3241_init, I3412_init, I3421_init, 
                I4123_init, I4132_init, I4213_init, I4231_init, I4312_init, I4321_init,
                PMs_init,
                PM_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  
  N_pop_ts[,1] <- M_pop # by age class
  
  #and fill in for births <- every epidemic class can give birth
  #transform birth vector into one distributed by age and epidemic class
  
  
  # stab.struct <- get.age.struct.M(M_pop, c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  
  # foi gets spread across the biweeks as well
  # (it comes in as 1 value per year)
  # take your Muench-estimated rates of annual infection and convert
  # to biweekly probabilities of infection
  foi.biweek = as.list(c(unlist(lapply(lapply(foi, divide.rate, ntyr1 =ntyr), rep, ntyr)))) #here, it could easily be modulated to be seasonal - (including forced by precip/temp)
  
  
  #and make annual waning immunity rate into probability of waning immunity per biweek
  sigma.biweek = as.list(c(unlist(lapply(lapply(sigma, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #and the probability of waning heterotypic immunity per biweek
  hetero.biweek = as.list(c(unlist(lapply(lapply(rate_hetero, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(hetero.biweek)!=(length(times)-1) |length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
    print("input vectors of unequal length")
  }
  #and your aging probability to biweeks (aging probability per year is 100%)
  #age.biwk = rep(1-exp(-(rep(1, yrs)/ntyr)), each=26)
  
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 2:(length(times))){
    
    
    
    #here are the dynamics before and after the immune waning event
    if((times[i]<((yr.intro) + biwk.intro/26)) | (times[i]>((yr.intro) + biwk.intro/26))){ 
      
      
      
      biwk1 <- find.biweek(t=i, times=times)
      clim.mod <- clim.vect[biwk1]
      
      #and select the appropriate age vector, based on the year
      age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
      
      
      
      if(nrow(age.mult.sub)>11){
        print("error with age multiplier sub-selection")
      }
      
      age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
      age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
      age.mult.sub$dur[age.mult.sub$dur<0] <- 0
      age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
      
      
      #and the mortality rates
      if(times[i]<min(mort$year)){
        mort.df = subset(mort, year==min(year))
      }else{
        mort.df = subset(mort, year==trunc(times[i]))
      }
      #make annual mortality rate by age into probability of mortality per biweek
      #these are deaths per 1000 people.
      #we multiply by the number in the age classes, so need to divide by 1000 here
      mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
      mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
      
      
      #then, let her go from there!
      
      Tmat <- buildTMat_four_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }else if(times[i]==((yr.intro) + biwk.intro/26)){
      
      
      
      biwk1 <- find.biweek(t=i, times=times)
      clim.mod <- clim.vect[biwk1]
      
      #and select the appropriate age vector, based on the year
      age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
      
      
      
      if(nrow(age.mult.sub)>11){
        print("error with age multiplier sub-selection")
      }
      
      age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
      age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
      age.mult.sub$dur[age.mult.sub$dur<0] <- 0
      age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
      
      
      #and the mortality rates
      if(times[i]<min(mort$year)){
        mort.df = subset(mort, year==min(year))
      }else{
        mort.df = subset(mort, year==trunc(times[i]))
      }
      #make annual mortality rate by age into probability of mortality per biweek
      #these are deaths per 1000 people.
      #we multiply by the number in the age classes, so need to divide by 1000 here
      mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
      mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
      
      
      # in this iteration, it is a 4-strain matrix, but we arbitrarily move some proportion of those with previous 
      # exposure to DENV-2 to a state of naivity for that seroptype at the beginning of the year of introduction
      # so any previous "S" class that had a 2 now gets moved to an S class without a 2
      
      #need to find the index that corresponds to the disease state of interest
      #vertical indices are the length of the number of epidemic classes (c) x the number of age classes (s)
      #matrix is stacked with all the ages of a given disease status vertically on top of one another
      #so [1,1] is S age 1, and [2,1] is S age 2
      
      #list the first index in each disease class
      
      
      
      # S2 to S
      tmp <- prop_imm*(N_pop_ts[(10*s +1):(11*s), (i-1)]) #S2
      N_pop_ts[(10*s +1):(11*s), (i-1)] <- N_pop_ts[(10*s +1):(11*s), (i-1)]-tmp
      N_pop_ts[(1):(1*s), (i-1)] <-  N_pop_ts[(1):(1*s), (i-1)] + tmp #S
      
      
      # S12 to S1
      tmp <- prop_imm*(N_pop_ts[(37*s +1):(38*s), (i-1)]) #S12
      N_pop_ts[(37*s +1):(38*s), (i-1)] <- N_pop_ts[(37*s +1):(38*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + tmp #S1
      
      
      # S32 to S3
      tmp <- prop_imm*(N_pop_ts[(44*s +1):(45*s), (i-1)]) #S32
      N_pop_ts[(44*s +1):(45*s), (i-1)] <- N_pop_ts[(44*s +1):(45*s), (i-1)]-tmp
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + tmp #S3
      
      # S42 to S4
      tmp <- prop_imm*(N_pop_ts[(47*s +1):(48*s), (i-1)]) #S42
      N_pop_ts[(47*s +1):(48*s), (i-1)] <- N_pop_ts[(47*s +1):(48*s), (i-1)]-tmp
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + tmp #S4
      
      
      
      # S123 to S13
      tmp <- prop_imm*(N_pop_ts[(97*s +1):(98*s), (i-1)]) #S123
      N_pop_ts[(97*s +1):(98*s), (i-1)] <- N_pop_ts[(97*s +1):(98*s), (i-1)]-tmp
      N_pop_ts[(38*s +1):(39*s), (i-1)]   <-  N_pop_ts[(38*s +1):(39*s), (i-1)]  + tmp #S13
      
      
      # S124 to S14
      tmp <- prop_imm*(N_pop_ts[(98*s +1):(99*s), (i-1)]) #S124
      N_pop_ts[(98*s +1):(99*s), (i-1)] <- N_pop_ts[(98*s +1):(99*s), (i-1)]-tmp
      N_pop_ts[(39*s +1):(40*s), (i-1)]   <-  N_pop_ts[(39*s +1):(40*s), (i-1)]  + tmp #S14
      
      
      # S132 to S13
      tmp <- prop_imm*(N_pop_ts[(99*s +1):(100*s), (i-1)]) #S132
      N_pop_ts[(99*s +1):(100*s), (i-1)] <- N_pop_ts[(99*s +1):(100*s), (i-1)]-tmp
      N_pop_ts[(38*s +1):(39*s), (i-1)]   <-  N_pop_ts[(38*s +1):(39*s), (i-1)]  + tmp #S13
      
      
      # S142 to S14
      tmp <- prop_imm*(N_pop_ts[(101*s +1):(102*s), (i-1)]) #S142
      N_pop_ts[(101*s +1):(102*s), (i-1)] <- N_pop_ts[(101*s +1):(102*s), (i-1)]-tmp
      N_pop_ts[(39*s +1):(40*s), (i-1)]   <-  N_pop_ts[(39*s +1):(40*s), (i-1)]  + tmp #S14
      
      
      
      
      # S213 to S13
      tmp <- prop_imm*(N_pop_ts[(103*s +1):(104*s), (i-1)]) #S213
      N_pop_ts[(103*s +1):(104*s), (i-1)] <- N_pop_ts[(103*s +1):(104*s), (i-1)]-tmp
      N_pop_ts[(38*s +1):(39*s), (i-1)]   <-  N_pop_ts[(38*s +1):(39*s), (i-1)]  + tmp #S13
      
      # S214 to S14
      tmp <- prop_imm*(N_pop_ts[(104*s +1):(105*s), (i-1)]) #S214
      N_pop_ts[(104*s +1):(105*s), (i-1)] <- N_pop_ts[(104*s +1):(105*s), (i-1)]-tmp
      N_pop_ts[(39*s +1):(40*s), (i-1)]   <-  N_pop_ts[(39*s +1):(40*s), (i-1)]  + tmp #S14
      
      # S231 to S31
      tmp <- prop_imm*(N_pop_ts[(105*s +1):(106*s), (i-1)]) #S231
      N_pop_ts[(105*s +1):(106*s), (i-1)] <- N_pop_ts[(105*s +1):(106*s), (i-1)]-tmp
      N_pop_ts[(43*s +1):(44*s), (i-1)]   <-  N_pop_ts[(43*s +1):(44*s), (i-1)]  + tmp #S31
      
      # S234 to S34
      tmp <- prop_imm*(N_pop_ts[(106*s +1):(107*s), (i-1)]) #S234
      N_pop_ts[(106*s +1):(107*s), (i-1)] <- N_pop_ts[(106*s +1):(107*s), (i-1)]-tmp
      N_pop_ts[(45*s +1):(46*s), (i-1)]   <-  N_pop_ts[(45*s +1):(46*s), (i-1)]  + tmp #S34
      
      # S241 to S41
      tmp <- prop_imm*(N_pop_ts[(107*s +1):(108*s), (i-1)]) #S241
      N_pop_ts[(107*s +1):(108*s), (i-1)] <- N_pop_ts[(107*s +1):(108*s), (i-1)]-tmp
      N_pop_ts[(46*s +1):(47*s), (i-1)]   <-  N_pop_ts[(46*s +1):(47*s), (i-1)]  + tmp #S41
      
      # S243 to S43
      tmp <- prop_imm*(N_pop_ts[(108*s +1):(109*s), (i-1)]) #S243
      N_pop_ts[(108*s +1):(109*s), (i-1)] <- N_pop_ts[(108*s +1):(109*s), (i-1)]-tmp
      N_pop_ts[(48*s +1):(49*s), (i-1)]   <-  N_pop_ts[(48*s +1):(49*s), (i-1)]  + tmp #S43
      
      
      
      
      # S312 to S31
      tmp <- prop_imm*(N_pop_ts[(109*s +1):(110*s), (i-1)]) #S312
      N_pop_ts[(109*s +1):(110*s), (i-1)] <- N_pop_ts[(109*s +1):(110*s), (i-1)]-tmp
      N_pop_ts[(43*s +1):(44*s), (i-1)]   <-  N_pop_ts[(43*s +1):(44*s), (i-1)]  + tmp #S31
      
      # S342 to S34
      tmp <- prop_imm*(N_pop_ts[(114*s +1):(115*s), (i-1)]) #S342
      N_pop_ts[(114*s +1):(115*s), (i-1)] <- N_pop_ts[(114*s +1):(115*s), (i-1)]-tmp
      N_pop_ts[(45*s +1):(46*s), (i-1)]   <-  N_pop_ts[(45*s +1):(46*s), (i-1)]  + tmp #S34
      
      # S321 to S31
      tmp <- prop_imm*(N_pop_ts[(111*s +1):(112*s), (i-1)]) #S321
      N_pop_ts[(111*s +1):(112*s), (i-1)] <- N_pop_ts[(111*s +1):(112*s), (i-1)]-tmp
      N_pop_ts[(43*s +1):(44*s), (i-1)]   <-  N_pop_ts[(43*s +1):(44*s), (i-1)]  + tmp #S31
      
      # S324 to S34
      tmp <- prop_imm*(N_pop_ts[(112*s +1):(113*s), (i-1)]) #S324
      N_pop_ts[(112*s +1):(113*s), (i-1)] <- N_pop_ts[(112*s +1):(113*s), (i-1)]-tmp
      N_pop_ts[(45*s +1):(46*s), (i-1)]   <-  N_pop_ts[(45*s +1):(46*s), (i-1)]  + tmp #S34
      
      
      
      
      # S412 to S41
      tmp <- prop_imm*(N_pop_ts[(115*s +1):(116*s), (i-1)]) #S412
      N_pop_ts[(115*s +1):(116*s), (i-1)] <- N_pop_ts[(115*s +1):(116*s), (i-1)]-tmp
      N_pop_ts[(46*s +1):(47*s), (i-1)]   <-  N_pop_ts[(46*s +1):(47*s), (i-1)]  + tmp #S41
      
      
      # S423 to S43
      tmp <- prop_imm*(N_pop_ts[(118*s +1):(119*s), (i-1)]) #S423
      N_pop_ts[(118*s +1):(119*s), (i-1)] <- N_pop_ts[(118*s +1):(119*s), (i-1)]-tmp
      N_pop_ts[(48*s +1):(49*s), (i-1)]   <-  N_pop_ts[(48*s +1):(49*s), (i-1)]  + tmp #S43
      
      # S421 to S41
      tmp <- prop_imm*(N_pop_ts[(117*s +1):(118*s), (i-1)]) #S421 
      N_pop_ts[(117*s +1):(118*s), (i-1)] <- N_pop_ts[(117*s +1):(118*s), (i-1)]-tmp
      N_pop_ts[(46*s +1):(47*s), (i-1)]   <-  N_pop_ts[(46*s +1):(47*s), (i-1)]  + tmp #S41
      
      # S432 to S43
      tmp <- prop_imm*(N_pop_ts[(120*s +1):(121*s), (i-1)]) #S432
      N_pop_ts[(120*s +1):(121*s), (i-1)] <- N_pop_ts[(120*s +1):(121*s), (i-1)]-tmp
      N_pop_ts[(48*s +1):(49*s), (i-1)]   <-  N_pop_ts[(48*s +1):(49*s), (i-1)]  + tmp #S43
      
      
      # Sm to S134, S143, S314, S341, S413, S431 (evenly)
      tmp <- prop_imm*(N_pop_ts[(146*s +1):(147*s), (i-1)]) #Sm 
      N_pop_ts[(146*s +1):(147*s), (i-1)] <- N_pop_ts[(146*s +1):(147*s), (i-1)] - tmp
      N_pop_ts[(100*s +1):(101*s), (i-1)]   <-  N_pop_ts[(100*s +1):(101*s), (i-1)]  + (tmp/6) #S134
      N_pop_ts[(102*s +1):(103*s), (i-1)]   <-  N_pop_ts[(102*s +1):(103*s), (i-1)]  + (tmp/6) #S143
      
      N_pop_ts[(110*s +1):(111*s), (i-1)]   <-  N_pop_ts[(110*s +1):(111*s), (i-1)]  + (tmp/6) #S314
      N_pop_ts[(113*s +1):(114*s), (i-1)]   <-  N_pop_ts[(113*s +1):(114*s), (i-1)]  + (tmp/6) #S341
      
      N_pop_ts[(116*s +1):(117*s), (i-1)]   <-  N_pop_ts[(116*s +1):(117*s), (i-1)]  + (tmp/6) #S413
      N_pop_ts[(119*s +1):(120*s), (i-1)]   <-  N_pop_ts[(119*s +1):(120*s), (i-1)]  + (tmp/6) #S431
      
      
      
      
      #then, let her go from there!
      
      Tmat <- buildTMat_four_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }
    
    
    
  }
  
  
  
  
  #the first column is the initial conditions. everything after that represents transitions within the year
  #should take the ceiling of each biweek to total within the year
  
  
  # stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  if(s>1){
    N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  }
  
  N.sim.dat <- as.data.frame(N_pop_ts)
  
  names(N.sim.dat) <- times
  
  N.sim.dat$class <- rep(c("S", "I1","I2", "I3", "I4",  "P1", "P2", "P3", "P4", "S1", "S2",  "S3", "S4",
                           "I12", "I13", "I14","I21", "I23", "I24", "I31", "I32", "I34", "I41", "I42", "I43",
                           "P12", "P13", "P14","P21", "P23", "P24", "P31", "P32", "P34",  "P41", "P42", "P43",
                           "S12", "S13", "S14","S21", "S23", "S24", "S31", "S32", "S34", "S41", "S42", "S43",
                           "I123","I124", "I132", "I134", "I142", "I143",
                           "I213", "I214","I231","I234", "I241", "I243",
                           "I312","I314", "I321", "I324", "I341", "I342",
                           "I412", "I413", "I421", "I423", "I431", "I432",
                           "P123","P124", "P132", "P134", "P142", "P143",
                           "P213", "P214","P231","P234", "P241", "P243",
                           "P312","P314", "P321", "P324", "P341", "P342",
                           "P412", "P413", "P421", "P423", "P431", "P432",
                           "S123","S124", "S132", "S134", "S142", "S143",
                           "S213", "S214","S231","S234", "S241", "S243",
                           "S312","S314", "S321", "S324", "S341", "S342",
                           "S412", "S413", "S421", "S423", "S431", "S432",
                           "I1234","I1243", "I1324", "I1342", "I1423", "I1432",
                           "I2134", "I2143","I2314","I2341", "I2413", "I2431",
                           "I3124","I3142", "I3214", "I3241", "I3412", "I3421",
                           "I4123", "I4132", "I4213", "I4231", "I4312", "I4321",
                           "Pms", "Pm"), 
                         each = s)
  N.sim.dat$age <- rep(1:s, length(unique(N.sim.dat$class)))
  
  #and melt
  N.sim.df <- melt(N.sim.dat, id=c("class", "age"), value.name = "count", variable.name = "time")
  
  N.sim.df$time <- as.numeric(as.character(N.sim.df$time))
  N.sim.df$count <- as.numeric(as.character(N.sim.df$count))
  
  #and add year
  #N.sim.df$year <- ceiling(N.sim.df$time)
  N.sim.df$year <- trunc(N.sim.df$time)
  #N.sim.df$time <- N.sim.df$time-1
  #N.sim.df$year <- N.sim.df$year-1
  
  
  
  return( N.sim.df)
}
sim.SIR.age.four.clim.wane.intro <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, rate_hetero, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro){
  
  
  
  
  
  #check on climate vector - should be length of steps within a year (often 26)
  #(this could be modulated to be the length of the time series if you wanted)
  if(length(clim.vect)==1){
    clim.vect = rep(clim.vect, ntyr)
  }
  
  
  #first, if mort, births, sigma, or foi are shorter than the time series,
  #make them match in length here.
  
  #this will also repeat if they come in as a vector (age-structured)
  
  #make as a list
  if(length(foi)<yrs){
    foi <- rep(list(foi), yrs)
  }else{
    foi <- as.list(foi)
  }
  
  if(length(sigma)< yrs){
    sigma = rep(list(sigma), yrs)
  }else{
    sigma <- as.list(sigma)
  }
  
  # if(length(mort)<yrs){
  #   mort = rep(list(mort), yrs)
  # }else{
  #   mort = as.list(mort)
  # }
  # 
  if(length(births)<yrs){
    births = rep(list(births), yrs)
  }else{
    births = as.list(births)
  }
  
  # 
  if(length(rate_hetero)<yrs){
    rate_hetero = rep(list(rate_hetero), yrs)
  }else{
    rate_hetero = as.list(rate_hetero)
  }
  
  
  #
  #number of epidemic classes
  c=737
  
  
  #length of time series
  times <-   seq(0, yrs, by =1/ntyr)
  times <-   times + (year.end-yrs)
  
  #split the births up by biweek too
  birth_vector_biwk = as.list(c(unlist(lapply(lapply(births, divide, ntyr1 =ntyr), rep, ntyr))))
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  #import the stable age structure from the literature
  #mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  #stab.struct = Re(eigen(mat1)$vector[,1])
  #stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  #lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  pop.mat = pop_vector #number of people per age class
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  #structure the population for 3 strains even though we only have 2 to start. they will get added later
  
  
  
  I5_init = rep(0, s);  #none to start
  I4_init = rep(0, s); I4_init[1] = 5 #comment out if you just want to check demography
  I3_init = rep(0, s); I3_init[1] = 5 #comment out if you just want to check demography
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  S_init = pop.mat - I1_init - I2_init - I3_init - I4_init
  
  #each column is a state variable
  #each row is an age class
  N_tot = matrix(0, nrow=s, ncol=c)
  
  #and seed the infections 
  N_tot[,1:6] <- cbind(S_init,I1_init, I2_init,I3_init,I4_init,I5_init)
  # 
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  
  N_pop_ts[,1] <- M_pop # by age class
  
  #and fill in for births <- every epidemic class can give birth
  #transform birth vector into one distributed by age and epidemic class
  
  
  # stab.struct <- get.age.struct.M(M_pop, c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  
  # foi gets spread across the biweeks as well
  # (it comes in as 1 value per year)
  # take your Muench-estimated rates of annual infection and convert
  # to biweekly probabilities of infection
  foi.biweek = as.list(c(unlist(lapply(lapply(foi, divide.rate, ntyr1 =ntyr), rep, ntyr)))) #here, it could easily be modulated to be seasonal - (including forced by precip/temp)
  
  
  #and make annual waning immunity rate into probability of waning immunity per biweek
  sigma.biweek = as.list(c(unlist(lapply(lapply(sigma, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #and the probability of waning heterotypic immunity per biweek
  hetero.biweek = as.list(c(unlist(lapply(lapply(rate_hetero, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(hetero.biweek)!=(length(times)-1) |length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
    print("input vectors of unequal length")
  }
  #and your aging probability to biweeks (aging probability per year is 100%)
  #age.biwk = rep(1-exp(-(rep(1, yrs)/ntyr)), each=26)
  
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 2:(length(times))){
    
    
    
    #here are the dynamics before the immune waning event
    if((times[i]<((yr.intro) + biwk.intro/26))){ 
      
      
      
      biwk1 <- find.biweek(t=i, times=times)
      clim.mod <- clim.vect[biwk1]
      
      #and select the appropriate age vector, based on the year
      age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
      
      
      
      if(nrow(age.mult.sub)>11){
        print("error with age multiplier sub-selection")
      }
      
      age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
      age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
      age.mult.sub$dur[age.mult.sub$dur<0] <- 0
      age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
      
      
      #and the mortality rates
      if(times[i]<min(mort$year)){
        mort.df = subset(mort, year==min(year))
      }else{
        mort.df = subset(mort, year==trunc(times[i]))
      }
      #make annual mortality rate by age into probability of mortality per biweek
      #these are deaths per 1000 people.
      #we multiply by the number in the age classes, so need to divide by 1000 here
      mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
      mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
      
      
      #then, let her go from there!
      
      Tmat <- buildTMat_five_strain_intro(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }else if(times[i]==((yr.intro) + biwk.intro/26)){
      #here are the dynamics at the immune waning event : 
      #introduce new infecteds into infected class 5
      
      N_pop_ts[6, (i-1)] <- 5 #I5 -- here at the lowest age class
      N_pop_ts[1, (i-1)] <-(N_pop_ts[1, (i-1)]) - 5 #remove from susceptibles
      
      
      
      biwk1 <- find.biweek(t=i, times=times)
      clim.mod <- clim.vect[biwk1]
      
      #and select the appropriate age vector, based on the year
      age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
      
      
      
      if(nrow(age.mult.sub)>11){
        print("error with age multiplier sub-selection")
      }
      
      age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
      age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
      age.mult.sub$dur[age.mult.sub$dur<0] <- 0
      age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
      
      
      #and the mortality rates
      if(times[i]<min(mort$year)){
        mort.df = subset(mort, year==min(year))
      }else{
        mort.df = subset(mort, year==trunc(times[i]))
      }
      #make annual mortality rate by age into probability of mortality per biweek
      #these are deaths per 1000 people.
      #we multiply by the number in the age classes, so need to divide by 1000 here
      mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
      mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
      
      
      
      #then, let her go from there!
      
      Tmat <- buildTMat_five_strain_intro(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }else if ((times[i]>((yr.intro) + biwk.intro/26))){
      #here are the dynamics after the immune waning event : 
      #use the new Tmat with no foi for serotype 2
      
      #here are the dynamics at the immune waning event : 
      #introduce new infecteds into infected class 5
      
      
      biwk1 <- find.biweek(t=i, times=times)
      clim.mod <- clim.vect[biwk1]
      
      #and select the appropriate age vector, based on the year
      age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
      
      
      
      if(nrow(age.mult.sub)>11){
        print("error with age multiplier sub-selection")
      }
      
      age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
      age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
      age.mult.sub$dur[age.mult.sub$dur<0] <- 0
      age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
      
      
      #and the mortality rates
      if(times[i]<min(mort$year)){
        mort.df = subset(mort, year==min(year))
      }else{
        mort.df = subset(mort, year==trunc(times[i]))
      }
      #make annual mortality rate by age into probability of mortality per biweek
      #these are deaths per 1000 people.
      #we multiply by the number in the age classes, so need to divide by 1000 here
      mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
      mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
      
      
      
      #then, let her go from there!
      
      Tmat <- buildTMat_five_strain_after(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
      
    }
    
    
    
  }
  
  
  
  
  #the first column is the initial conditions. everything after that represents transitions within the year
  #should take the ceiling of each biweek to total within the year
  
  
  # stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  if(s>1){
    N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  }
  
  N.sim.dat <- as.data.frame(N_pop_ts)
  
  names(N.sim.dat) <- times
  
  N.sim.dat$class <- rep( rev(c("PM",
                                "I54321",
                                "I54312",
                                "I54231",
                                "I54213",
                                "I54132",
                                "I54123",
                                "I53421",
                                "I53412",
                                "I53241",
                                "I53214",
                                "I53142",
                                "I53124",
                                "I52431",
                                "I52413",
                                "I52341",
                                "I52314",
                                "I52143",
                                "I52134",
                                "I51432",
                                "I51423",
                                "I51342",
                                "I51324",
                                "I51243",
                                "I51234",
                                "I45321",
                                "I45312",
                                "I45231",
                                "I45213",
                                "I45132",
                                "I45123",
                                "I43521",
                                "I43512",
                                "I43251",
                                "I43215",
                                "I43152",
                                "I43125",
                                "I42531",
                                "I42513",
                                "I42351",
                                "I42315",
                                "I42153",
                                "I42135",
                                "I41532",
                                "I41523",
                                "I41352",
                                "I41325",
                                "I41253",
                                "I41235",
                                "I35421",
                                "I35412",
                                "I35241",
                                "I35214",
                                "I35142",
                                "I35124",
                                "I34521",
                                "I34512",
                                "I34251",
                                "I34215",
                                "I34152",
                                "I34125",
                                "I32541",
                                "I32514",
                                "I32451",
                                "I32415",
                                "I32154",
                                "I32145",
                                "I31542",
                                "I31524",
                                "I31452",
                                "I31425",
                                "I31254",
                                "I31245",
                                "I25431",
                                "I25413",
                                "I25341",
                                "I25314",
                                "I25143",
                                "I25134",
                                "I24531",
                                "I24513",
                                "I24351",
                                "I24315",
                                "I24153",
                                "I24135",
                                "I23541",
                                "I23514",
                                "I23451",
                                "I23415",
                                "I23154",
                                "I23145",
                                "I21543",
                                "I21534",
                                "I21453",
                                "I21435",
                                "I21354",
                                "I21345",
                                "I15432",
                                "I15423",
                                "I15342",
                                "I15324",
                                "I15243",
                                "I15234",
                                "I14532",
                                "I14523",
                                "I14352",
                                "I14325",
                                "I14253",
                                "I14235",
                                "I13542",
                                "I13524",
                                "I13452",
                                "I13425",
                                "I13254",
                                "I13245",
                                "I12543",
                                "I12534",
                                "I12453",
                                "I12435",
                                "I12354",
                                "I12345",
                                "S5432",
                                "S5431",
                                "S5423",
                                "S5421",
                                "S5413",
                                "S5412",
                                "S5342",
                                "S5341",
                                "S5324",
                                "S5321",
                                "S5314",
                                "S5312",
                                "S5243",
                                "S5241",
                                "S5234",
                                "S5231",
                                "S5214",
                                "S5213",
                                "S5143",
                                "S5142",
                                "S5134",
                                "S5132",
                                "S5124",
                                "S5123",
                                "S4532",
                                "S4531",
                                "S4523",
                                "S4521",
                                "S4513",
                                "S4512",
                                "S4352",
                                "S4351",
                                "S4325",
                                "S4321",
                                "S4315",
                                "S4312",
                                "S4253",
                                "S4251",
                                "S4235",
                                "S4231",
                                "S4215",
                                "S4213",
                                "S4153",
                                "S4152",
                                "S4135",
                                "S4132",
                                "S4125",
                                "S4123",
                                "S3542",
                                "S3541",
                                "S3524",
                                "S3521",
                                "S3514",
                                "S3512",
                                "S3452",
                                "S3451",
                                "S3425",
                                "S3421",
                                "S3415",
                                "S3412",
                                "S3254",
                                "S3251",
                                "S3245",
                                "S3241",
                                "S3215",
                                "S3214",
                                "S3154",
                                "S3152",
                                "S3145",
                                "S3142",
                                "S3125",
                                "S3124",
                                "S2543",
                                "S2541",
                                "S2534",
                                "S2531",
                                "S2514",
                                "S2513",
                                "S2453",
                                "S2451",
                                "S2435",
                                "S2431",
                                "S2415",
                                "S2413",
                                "S2354",
                                "S2351",
                                "S2345",
                                "S2341",
                                "S2315",
                                "S2314",
                                "S2154",
                                "S2153",
                                "S2145",
                                "S2143",
                                "S2135",
                                "S2134",
                                "S1543",
                                "S1542",
                                "S1534",
                                "S1532",
                                "S1524",
                                "S1523",
                                "S1453",
                                "S1452",
                                "S1435",
                                "S1432",
                                "S1425",
                                "S1423",
                                "S1354",
                                "S1352",
                                "S1345",
                                "S1342",
                                "S1325",
                                "S1324",
                                "S1254",
                                "S1253",
                                "S1245",
                                "S1243",
                                "S1235",
                                "S1234",
                                "P5432",
                                "P5431",
                                "P5423",
                                "P5421",
                                "P5413",
                                "P5412",
                                "P5342",
                                "P5341",
                                "P5324",
                                "P5321",
                                "P5314",
                                "P5312",
                                "P5243",
                                "P5241",
                                "P5234",
                                "P5231",
                                "P5214",
                                "P5213",
                                "P5143",
                                "P5142",
                                "P5134",
                                "P5132",
                                "P5124",
                                "P5123",
                                "P4532",
                                "P4531",
                                "P4523",
                                "P4521",
                                "P4513",
                                "P4512",
                                "P4352",
                                "P4351",
                                "P4325",
                                "P4321",
                                "P4315",
                                "P4312",
                                "P4253",
                                "P4251",
                                "P4235",
                                "P4231",
                                "P4215",
                                "P4213",
                                "P4153",
                                "P4152",
                                "P4135",
                                "P4132",
                                "P4125",
                                "P4123",
                                "P3542",
                                "P3541",
                                "P3524",
                                "P3521",
                                "P3514",
                                "P3512",
                                "P3452",
                                "P3451",
                                "P3425",
                                "P3421",
                                "P3415",
                                "P3412",
                                "P3254",
                                "P3251",
                                "P3245",
                                "P3241",
                                "P3215",
                                "P3214",
                                "P3154",
                                "P3152",
                                "P3145",
                                "P3142",
                                "P3125",
                                "P3124",
                                "P2543",
                                "P2541",
                                "P2534",
                                "P2531",
                                "P2514",
                                "P2513",
                                "P2453",
                                "P2451",
                                "P2435",
                                "P2431",
                                "P2415",
                                "P2413",
                                "P2354",
                                "P2351",
                                "P2345",
                                "P2341",
                                "P2315",
                                "P2314",
                                "P2154",
                                "P2153",
                                "P2145",
                                "P2143",
                                "P2135",
                                "P2134",
                                "P1543",
                                "P1542",
                                "P1534",
                                "P1532",
                                "P1524",
                                "P1523",
                                "P1453",
                                "P1452",
                                "P1435",
                                "P1432",
                                "P1425",
                                "P1423",
                                "P1354",
                                "P1352",
                                "P1345",
                                "P1342",
                                "P1325",
                                "P1324",
                                "P1254",
                                "P1253",
                                "P1245",
                                "P1243",
                                "P1235",
                                "P1234",
                                "I5432",
                                "I5431",
                                "I5423",
                                "I5421",
                                "I5413",
                                "I5412",
                                "I5342",
                                "I5341",
                                "I5324",
                                "I5321",
                                "I5314",
                                "I5312",
                                "I5243",
                                "I5241",
                                "I5234",
                                "I5231",
                                "I5214",
                                "I5213",
                                "I5143",
                                "I5142",
                                "I5134",
                                "I5132",
                                "I5124",
                                "I5123",
                                "I4532",
                                "I4531",
                                "I4523",
                                "I4521",
                                "I4513",
                                "I4512",
                                "I4352",
                                "I4351",
                                "I4325",
                                "I4321",
                                "I4315",
                                "I4312",
                                "I4253",
                                "I4251",
                                "I4235",
                                "I4231",
                                "I4215",
                                "I4213",
                                "I4153",
                                "I4152",
                                "I4135",
                                "I4132",
                                "I4125",
                                "I4123",
                                "I3542",
                                "I3541",
                                "I3524",
                                "I3521",
                                "I3514",
                                "I3512",
                                "I3452",
                                "I3451",
                                "I3425",
                                "I3421",
                                "I3415",
                                "I3412",
                                "I3254",
                                "I3251",
                                "I3245",
                                "I3241",
                                "I3215",
                                "I3214",
                                "I3154",
                                "I3152",
                                "I3145",
                                "I3142",
                                "I3125",
                                "I3124",
                                "I2543",
                                "I2541",
                                "I2534",
                                "I2531",
                                "I2514",
                                "I2513",
                                "I2453",
                                "I2451",
                                "I2435",
                                "I2431",
                                "I2415",
                                "I2413",
                                "I2354",
                                "I2351",
                                "I2345",
                                "I2341",
                                "I2315",
                                "I2314",
                                "I2154",
                                "I2153",
                                "I2145",
                                "I2143",
                                "I2135",
                                "I2134",
                                "I1543",
                                "I1542",
                                "I1534",
                                "I1532",
                                "I1524",
                                "I1523",
                                "I1453",
                                "I1452",
                                "I1435",
                                "I1432",
                                "I1425",
                                "I1423",
                                "I1354",
                                "I1352",
                                "I1345",
                                "I1342",
                                "I1325",
                                "I1324",
                                "I1254",
                                "I1253",
                                "I1245",
                                "I1243",
                                "I1235",
                                "I1234",
                                "S543",
                                "S542",
                                "S541",
                                "S534",
                                "S532",
                                "S531",
                                "S524",
                                "S523",
                                "S521",
                                "S514",
                                "S513",
                                "S512",
                                "S453",
                                "S452",
                                "S451",
                                "S435",
                                "S432",
                                "S431",
                                "S425",
                                "S423",
                                "S421",
                                "S415",
                                "S413",
                                "S412",
                                "S354",
                                "S352",
                                "S351",
                                "S345",
                                "S342",
                                "S341",
                                "S325",
                                "S324",
                                "S321",
                                "S315",
                                "S314",
                                "S312",
                                "S254",
                                "S253",
                                "S251",
                                "S245",
                                "S243",
                                "S241",
                                "S235",
                                "S234",
                                "S231",
                                "S215",
                                "S214",
                                "S213",
                                "S154",
                                "S153",
                                "S152",
                                "S145",
                                "S143",
                                "S142",
                                "S135",
                                "S134",
                                "S132",
                                "S125",
                                "S124",
                                "S123",
                                "P543",
                                "P542",
                                "P541",
                                "P534",
                                "P532",
                                "P531",
                                "P524",
                                "P523",
                                "P521",
                                "P514",
                                "P513",
                                "P512",
                                "P453",
                                "P452",
                                "P451",
                                "P435",
                                "P432",
                                "P431",
                                "P425",
                                "P423",
                                "P421",
                                "P415",
                                "P413",
                                "P412",
                                "P354",
                                "P352",
                                "P351",
                                "P345",
                                "P342",
                                "P341",
                                "P325",
                                "P324",
                                "P321",
                                "P315",
                                "P314",
                                "P312",
                                "P254",
                                "P253",
                                "P251",
                                "P245",
                                "P243",
                                "P241",
                                "P235",
                                "P234",
                                "P231",
                                "P215",
                                "P214",
                                "P213",
                                "P154",
                                "P153",
                                "P152",
                                "P145",
                                "P143",
                                "P142",
                                "P135",
                                "P134",
                                "P132",
                                "P125",
                                "P124",
                                "P123",
                                "I543",
                                "I542",
                                "I541",
                                "I534",
                                "I532",
                                "I531",
                                "I524",
                                "I523",
                                "I521",
                                "I514",
                                "I513",
                                "I512",
                                "I453",
                                "I452",
                                "I451",
                                "I435",
                                "I432",
                                "I431",
                                "I425",
                                "I423",
                                "I421",
                                "I415",
                                "I413",
                                "I412",
                                "I354",
                                "I352",
                                "I351",
                                "I345",
                                "I342",
                                "I341",
                                "I325",
                                "I324",
                                "I321",
                                "I315",
                                "I314",
                                "I312",
                                "I254",
                                "I253",
                                "I251",
                                "I245",
                                "I243",
                                "I241",
                                "I235",
                                "I234",
                                "I231",
                                "I215",
                                "I214",
                                "I213",
                                "I154",
                                "I153",
                                "I152",
                                "I145",
                                "I143",
                                "I142",
                                "I135",
                                "I134",
                                "I132",
                                "I125",
                                "I124",
                                "I123",
                                "S54",
                                "S53",
                                "S52",
                                "S51",
                                "S45",
                                "S43",
                                "S42",
                                "S41",
                                "S35",
                                "S34",
                                "S32",
                                "S31",
                                "S25",
                                "S24",
                                "S23",
                                "S21",
                                "S15",
                                "S14",
                                "S13",
                                "S12",
                                "P54",
                                "P53",
                                "P52",
                                "P51",
                                "P45",
                                "P43",
                                "P42",
                                "P41",
                                "P35",
                                "P34",
                                "P32",
                                "P31",
                                "P25",
                                "P24",
                                "P23",
                                "P21",
                                "P15",
                                "P14",
                                "P13",
                                "P12",
                                "I54",
                                "I53",
                                "I52",
                                "I51",
                                "I45",
                                "I43",
                                "I42",
                                "I41",
                                "I35",
                                "I34",
                                "I32",
                                "I31",
                                "I25",
                                "I24",
                                "I23",
                                "I21",
                                "I15",
                                "I14",
                                "I13",
                                "I12",
                                "S5",
                                "S4",
                                "S3",
                                "S2",
                                "S1",
                                "P5",
                                "P4",
                                "P3",
                                "P2",
                                "P1",
                                "I5",
                                "I4",
                                "I3",
                                "I2",
                                "I1",
                                "S")), each = s)
  N.sim.dat$age <- rep(1:s, length(unique(N.sim.dat$class)))
  
  #and melt
  N.sim.df <- melt(N.sim.dat, id=c("class", "age"), value.name = "count", variable.name = "time")
  
  N.sim.df$time <- as.numeric(as.character(N.sim.df$time))
  N.sim.df$count <- as.numeric(as.character(N.sim.df$count))
  
  #and add year
  #N.sim.df$year <- ceiling(N.sim.df$time)
  N.sim.df$year <- trunc(N.sim.df$time)
  #N.sim.df$time <- N.sim.df$time-1
  #N.sim.df$year <- N.sim.df$year-1
  
  
  
  return( N.sim.df)
}
sim.SIR.age.four.clim.wane.sec <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, rate_hetero, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro, prop_imm){
  
  
  
  
  
  #check on climate vector - should be length of steps within a year (often 26)
  #(this could be modulated to be the length of the time series if you wanted)
  if(length(clim.vect)==1){
    clim.vect = rep(clim.vect, ntyr)
  }
  
  
  #first, if mort, births, sigma, or foi are shorter than the time series,
  #make them match in length here.
  
  #this will also repeat if they come in as a vector (age-structured)
  
  #make as a list
  if(length(foi)<yrs){
    foi <- rep(list(foi), yrs)
  }else{
    foi <- as.list(foi)
  }
  
  if(length(sigma)< yrs){
    sigma = rep(list(sigma), yrs)
  }else{
    sigma <- as.list(sigma)
  }
  
  # if(length(mort)<yrs){
  #   mort = rep(list(mort), yrs)
  # }else{
  #   mort = as.list(mort)
  # }
  # 
  if(length(births)<yrs){
    births = rep(list(births), yrs)
  }else{
    births = as.list(births)
  }
  
  # 
  if(length(rate_hetero)<yrs){
    rate_hetero = rep(list(rate_hetero), yrs)
  }else{
    rate_hetero = as.list(rate_hetero)
  }
  
  
  #
  #number of epidemic classes
  c=147
  
  
  #length of time series
  times <-   seq(0, yrs, by =1/ntyr)
  times <-   times + (year.end-yrs)
  
  #split the births up by biweek too
  birth_vector_biwk = as.list(c(unlist(lapply(lapply(births, divide, ntyr1 =ntyr), rep, ntyr))))
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  #import the stable age structure from the literature
  #mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  #stab.struct = Re(eigen(mat1)$vector[,1])
  #stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  #lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  pop.mat = pop_vector #number of people per age class
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  #structure the population for 3 strains even though we only have 2 to start. they will get added later
  
  PM_init = rep(0, s)
  
  PMs_init = rep(0, s)
  
  I4321_init = rep(0, s)
  I4312_init = rep(0, s)
  I4231_init = rep(0, s)
  I4213_init = rep(0, s)
  I4132_init = rep(0, s)
  I4123_init = rep(0, s)
  
  I3421_init = rep(0, s)
  I3412_init = rep(0, s)
  I3241_init = rep(0, s)
  I3214_init = rep(0, s)
  I3142_init = rep(0, s)
  I3124_init = rep(0, s)
  
  I2431_init = rep(0, s)
  I2413_init = rep(0, s)
  I2341_init = rep(0, s)
  I2314_init = rep(0, s)
  I2143_init = rep(0, s)
  I2134_init = rep(0, s)
  
  I1432_init = rep(0, s)
  I1423_init = rep(0, s)
  I1342_init = rep(0, s)
  I1324_init = rep(0, s)
  I1243_init = rep(0, s)
  I1234_init = rep(0, s)
  
  
  S432_init = rep(0, s)
  S431_init = rep(0, s)
  S423_init = rep(0, s)
  S421_init = rep(0, s)
  S413_init = rep(0, s)
  S412_init = rep(0, s)
  
  S342_init = rep(0, s)
  S341_init = rep(0, s)
  S324_init = rep(0, s)
  S321_init = rep(0, s)
  S314_init = rep(0, s)
  S312_init = rep(0, s)
  
  S243_init = rep(0, s)
  S241_init = rep(0, s)
  S234_init = rep(0, s)
  S231_init = rep(0, s)
  S214_init = rep(0, s)
  S213_init = rep(0, s)
  
  S143_init = rep(0, s)
  S142_init = rep(0, s)
  S134_init = rep(0, s)
  S132_init = rep(0, s)
  S124_init = rep(0, s)
  S123_init = rep(0, s)
  
  
  P432_init = rep(0, s)
  P431_init = rep(0, s)
  P423_init = rep(0, s)
  P421_init = rep(0, s)
  P413_init = rep(0, s)
  P412_init = rep(0, s)
  
  P342_init = rep(0, s)
  P341_init = rep(0, s)
  P324_init = rep(0, s)
  P321_init = rep(0, s)
  P314_init = rep(0, s)
  P312_init = rep(0, s)
  
  P243_init = rep(0, s)
  P241_init = rep(0, s)
  P234_init = rep(0, s)
  P231_init = rep(0, s)
  P214_init = rep(0, s)
  P213_init = rep(0, s)
  
  P143_init = rep(0, s)
  P142_init = rep(0, s)
  P134_init = rep(0, s)
  P132_init = rep(0, s)
  P124_init = rep(0, s)
  P123_init = rep(0, s)
  
  
  I432_init = rep(0, s)
  I431_init = rep(0, s)
  I423_init = rep(0, s)
  I421_init = rep(0, s)
  I413_init = rep(0, s)
  I412_init = rep(0, s)
  
  I342_init = rep(0, s)
  I341_init = rep(0, s)
  I324_init = rep(0, s)
  I321_init = rep(0, s)
  I314_init = rep(0, s)
  I312_init = rep(0, s)
  
  I243_init = rep(0, s)
  I241_init = rep(0, s)
  I234_init = rep(0, s)
  I231_init = rep(0, s)
  I214_init = rep(0, s)
  I213_init = rep(0, s)
  
  I143_init = rep(0, s)
  I142_init = rep(0, s)
  I134_init = rep(0, s)
  I132_init = rep(0, s)
  I124_init = rep(0, s)
  I123_init = rep(0, s)
  
  S43_init = rep(0, s)
  S42_init = rep(0, s)
  S41_init = rep(0, s)
  S34_init = rep(0, s)
  S32_init = rep(0, s)
  S31_init = rep(0, s)
  S24_init = rep(0, s)
  S23_init = rep(0, s)
  S21_init = rep(0, s)
  S14_init = rep(0, s)
  S13_init = rep(0, s)
  S12_init = rep(0, s)
  
  P43_init = rep(0, s)
  P42_init = rep(0, s)
  P41_init = rep(0, s)
  P34_init = rep(0, s)
  P32_init = rep(0, s)
  P31_init = rep(0, s)
  P24_init = rep(0, s)
  P23_init = rep(0, s)
  P21_init = rep(0, s)
  P14_init = rep(0, s)
  P13_init = rep(0, s)
  P12_init = rep(0, s)
  
  I43_init = rep(0, s)
  I42_init = rep(0, s)
  I41_init = rep(0, s)
  I34_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I24_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I14_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  
  S4_init = rep(0, s)
  S3_init = rep(0, s)
  S2_init = rep(0, s)
  S1_init = rep(0, s)
  
  P4_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  
  I4_init = rep(0, s); I4_init[1] = 5 #comment out if you just want to check demography
  I3_init = rep(0, s); I3_init[1] = 5 #comment out if you just want to check demography
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  
  S_init = pop.mat - I1_init - I2_init - I3_init - I4_init
  
  
  N_tot = cbind(S_init, 
                I1_init, I2_init,I3_init,I4_init,
                P1_init, P2_init,P3_init,P4_init, 
                S1_init,S2_init,S3_init,S4_init,
                I12_init, I13_init, I14_init, 
                I21_init, I23_init, I34_init, 
                I31_init, I32_init, I34_init,
                I41_init, I42_init, I43_init,
                P12_init, P13_init, P14_init, 
                P21_init, P23_init, P34_init, 
                P31_init, P32_init, P34_init,
                P41_init, P42_init, P43_init,
                S12_init, S13_init, S14_init, 
                S21_init, S23_init, S34_init, 
                S31_init, S32_init, S34_init,
                S41_init, S42_init, S43_init,
                I143_init,I142_init,I134_init,I132_init,I124_init,I123_init,
                I213_init,I214_init, I231_init, I234_init, I241_init, I243_init,
                I312_init, I314_init, I321_init, I324_init, I341_init, I342_init,
                I412_init, I413_init, I421_init, I423_init,I431_init, I432_init, 
                P143_init,P142_init,P134_init,P132_init,P124_init,P123_init,
                P213_init,P214_init, P231_init, P234_init, P241_init, P243_init,
                P312_init, P314_init, P321_init, P324_init, P341_init, P342_init,
                P412_init, P413_init, P421_init, P423_init,P431_init, P432_init, 
                S143_init,S142_init,S134_init,S132_init,S124_init,S123_init,
                S213_init,S214_init, S231_init, S234_init, S241_init, S243_init,
                S312_init, S314_init, S321_init, S324_init, S341_init, S342_init,
                S412_init, S413_init, S421_init, S423_init,S431_init, S432_init,  
                I1234_init, I1243_init, I1324_init, I1342_init,I1423_init,I1432_init,
                I2134_init, I2143_init, I2314_init, I2341_init,I2413_init, I2431_init, 
                I3124_init, I3142_init, I3214_init, I3241_init, I3412_init, I3421_init, 
                I4123_init, I4132_init, I4213_init, I4231_init, I4312_init, I4321_init,
                PMs_init,
                PM_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  
  N_pop_ts[,1] <- M_pop # by age class
  
  #and fill in for births <- every epidemic class can give birth
  #transform birth vector into one distributed by age and epidemic class
  
  
  # stab.struct <- get.age.struct.M(M_pop, c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  
  # foi gets spread across the biweeks as well
  # (it comes in as 1 value per year)
  # take your Muench-estimated rates of annual infection and convert
  # to biweekly probabilities of infection
  foi.biweek = as.list(c(unlist(lapply(lapply(foi, divide.rate, ntyr1 =ntyr), rep, ntyr)))) #here, it could easily be modulated to be seasonal - (including forced by precip/temp)
  
  
  #and make annual waning immunity rate into probability of waning immunity per biweek
  sigma.biweek = as.list(c(unlist(lapply(lapply(sigma, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #and the probability of waning heterotypic immunity per biweek
  hetero.biweek = as.list(c(unlist(lapply(lapply(rate_hetero, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(hetero.biweek)!=(length(times)-1) |length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
    print("input vectors of unequal length")
  }
  #and your aging probability to biweeks (aging probability per year is 100%)
  #age.biwk = rep(1-exp(-(rep(1, yrs)/ntyr)), each=26)
  
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 2:(length(times))){
    
    
    
    #here are the dynamics before and after the immune waning event
    if((times[i]<((yr.intro) + biwk.intro/26)) | (times[i]>((yr.intro) + biwk.intro/26))){ 
      
      
      
      biwk1 <- find.biweek(t=i, times=times)
      clim.mod <- clim.vect[biwk1]
      
      #and select the appropriate age vector, based on the year
      age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
      
      
      
      if(nrow(age.mult.sub)>11){
        print("error with age multiplier sub-selection")
      }
      
      age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
      age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
      age.mult.sub$dur[age.mult.sub$dur<0] <- 0
      age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
      
      
      #and the mortality rates
      if(times[i]<min(mort$year)){
        mort.df = subset(mort, year==min(year))
      }else{
        mort.df = subset(mort, year==trunc(times[i]))
      }
      #make annual mortality rate by age into probability of mortality per biweek
      #these are deaths per 1000 people.
      #we multiply by the number in the age classes, so need to divide by 1000 here
      mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
      mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
      
      
      #then, let her go from there!
      
      Tmat <- buildTMat_four_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }else if(times[i]==((yr.intro) + biwk.intro/26)){
      
      
      
      biwk1 <- find.biweek(t=i, times=times)
      clim.mod <- clim.vect[biwk1]
      
      #and select the appropriate age vector, based on the year
      age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
      
      
      
      if(nrow(age.mult.sub)>11){
        print("error with age multiplier sub-selection")
      }
      
      age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
      age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
      age.mult.sub$dur[age.mult.sub$dur<0] <- 0
      age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
      
      
      #and the mortality rates
      if(times[i]<min(mort$year)){
        mort.df = subset(mort, year==min(year))
      }else{
        mort.df = subset(mort, year==trunc(times[i]))
      }
      #make annual mortality rate by age into probability of mortality per biweek
      #these are deaths per 1000 people.
      #we multiply by the number in the age classes, so need to divide by 1000 here
      mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
      mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
      
      
      # in this iteration, it is a 4-strain matrix, but we arbitrarily move some proportion of those with previous 
      # exposure to DENV-2 to a state of naivity for that seroptype at the beginning of the year of introduction
      # so any previous "S" class that had a 2 now gets moved to an S class without a 2
      
      #need to find the index that corresponds to the disease state of interest
      #vertical indices are the length of the number of epidemic classes (c) x the number of age classes (s)
      #matrix is stacked with all the ages of a given disease status vertically on top of one another
      #so [1,1] is S age 1, and [2,1] is S age 2
      
      #list the first index in each disease class
      
      
      
      # S2 to S
      #tmp <- prop_imm*(N_pop_ts[(10*s +1):(11*s), (i-1)]) #S2
      #N_pop_ts[(10*s +1):(11*s), (i-1)] <- N_pop_ts[(10*s +1):(11*s), (i-1)]-tmp
      #N_pop_ts[(1):(1*s), (i-1)] <-  N_pop_ts[(1):(1*s), (i-1)] + tmp #S
      
      
      # S12 to S1
      tmp <- prop_imm*(N_pop_ts[(37*s +1):(38*s), (i-1)]) #S12
      N_pop_ts[(37*s +1):(38*s), (i-1)] <- N_pop_ts[(37*s +1):(38*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + tmp #S1
      
      
      # S32 to S3
      tmp <- prop_imm*(N_pop_ts[(44*s +1):(45*s), (i-1)]) #S32
      N_pop_ts[(44*s +1):(45*s), (i-1)] <- N_pop_ts[(44*s +1):(45*s), (i-1)]-tmp
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + tmp #S3
      
      # S42 to S4
      tmp <- prop_imm*(N_pop_ts[(47*s +1):(48*s), (i-1)]) #S42
      N_pop_ts[(47*s +1):(48*s), (i-1)] <- N_pop_ts[(47*s +1):(48*s), (i-1)]-tmp
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + tmp #S4
      
      
      
      # S123 to S13
      tmp <- prop_imm*(N_pop_ts[(97*s +1):(98*s), (i-1)]) #S123
      N_pop_ts[(97*s +1):(98*s), (i-1)] <- N_pop_ts[(97*s +1):(98*s), (i-1)]-tmp
      N_pop_ts[(38*s +1):(39*s), (i-1)]   <-  N_pop_ts[(38*s +1):(39*s), (i-1)]  + tmp #S13
      
      
      # S124 to S14
      tmp <- prop_imm*(N_pop_ts[(98*s +1):(99*s), (i-1)]) #S124
      N_pop_ts[(98*s +1):(99*s), (i-1)] <- N_pop_ts[(98*s +1):(99*s), (i-1)]-tmp
      N_pop_ts[(39*s +1):(40*s), (i-1)]   <-  N_pop_ts[(39*s +1):(40*s), (i-1)]  + tmp #S14
      
      
      # S132 to S13
      tmp <- prop_imm*(N_pop_ts[(99*s +1):(100*s), (i-1)]) #S132
      N_pop_ts[(99*s +1):(100*s), (i-1)] <- N_pop_ts[(99*s +1):(100*s), (i-1)]-tmp
      N_pop_ts[(38*s +1):(39*s), (i-1)]   <-  N_pop_ts[(38*s +1):(39*s), (i-1)]  + tmp #S13
      
      
      # S142 to S14
      tmp <- prop_imm*(N_pop_ts[(101*s +1):(102*s), (i-1)]) #S142
      N_pop_ts[(101*s +1):(102*s), (i-1)] <- N_pop_ts[(101*s +1):(102*s), (i-1)]-tmp
      N_pop_ts[(39*s +1):(40*s), (i-1)]   <-  N_pop_ts[(39*s +1):(40*s), (i-1)]  + tmp #S14
      
      
      
      
      # S213 to S13
      tmp <- prop_imm*(N_pop_ts[(103*s +1):(104*s), (i-1)]) #S213
      N_pop_ts[(103*s +1):(104*s), (i-1)] <- N_pop_ts[(103*s +1):(104*s), (i-1)]-tmp
      N_pop_ts[(38*s +1):(39*s), (i-1)]   <-  N_pop_ts[(38*s +1):(39*s), (i-1)]  + tmp #S13
      
      # S214 to S14
      tmp <- prop_imm*(N_pop_ts[(104*s +1):(105*s), (i-1)]) #S214
      N_pop_ts[(104*s +1):(105*s), (i-1)] <- N_pop_ts[(104*s +1):(105*s), (i-1)]-tmp
      N_pop_ts[(39*s +1):(40*s), (i-1)]   <-  N_pop_ts[(39*s +1):(40*s), (i-1)]  + tmp #S14
      
      # S231 to S31
      tmp <- prop_imm*(N_pop_ts[(105*s +1):(106*s), (i-1)]) #S231
      N_pop_ts[(105*s +1):(106*s), (i-1)] <- N_pop_ts[(105*s +1):(106*s), (i-1)]-tmp
      N_pop_ts[(43*s +1):(44*s), (i-1)]   <-  N_pop_ts[(43*s +1):(44*s), (i-1)]  + tmp #S31
      
      # S234 to S34
      tmp <- prop_imm*(N_pop_ts[(106*s +1):(107*s), (i-1)]) #S234
      N_pop_ts[(106*s +1):(107*s), (i-1)] <- N_pop_ts[(106*s +1):(107*s), (i-1)]-tmp
      N_pop_ts[(45*s +1):(46*s), (i-1)]   <-  N_pop_ts[(45*s +1):(46*s), (i-1)]  + tmp #S34
      
      # S241 to S41
      tmp <- prop_imm*(N_pop_ts[(107*s +1):(108*s), (i-1)]) #S241
      N_pop_ts[(107*s +1):(108*s), (i-1)] <- N_pop_ts[(107*s +1):(108*s), (i-1)]-tmp
      N_pop_ts[(46*s +1):(47*s), (i-1)]   <-  N_pop_ts[(46*s +1):(47*s), (i-1)]  + tmp #S41
      
      # S243 to S43
      tmp <- prop_imm*(N_pop_ts[(108*s +1):(109*s), (i-1)]) #S243
      N_pop_ts[(108*s +1):(109*s), (i-1)] <- N_pop_ts[(108*s +1):(109*s), (i-1)]-tmp
      N_pop_ts[(48*s +1):(49*s), (i-1)]   <-  N_pop_ts[(48*s +1):(49*s), (i-1)]  + tmp #S43
      
      
      
      
      # S312 to S31
      tmp <- prop_imm*(N_pop_ts[(109*s +1):(110*s), (i-1)]) #S312
      N_pop_ts[(109*s +1):(110*s), (i-1)] <- N_pop_ts[(109*s +1):(110*s), (i-1)]-tmp
      N_pop_ts[(43*s +1):(44*s), (i-1)]   <-  N_pop_ts[(43*s +1):(44*s), (i-1)]  + tmp #S31
      
      # S342 to S34
      tmp <- prop_imm*(N_pop_ts[(114*s +1):(115*s), (i-1)]) #S342
      N_pop_ts[(114*s +1):(115*s), (i-1)] <- N_pop_ts[(114*s +1):(115*s), (i-1)]-tmp
      N_pop_ts[(45*s +1):(46*s), (i-1)]   <-  N_pop_ts[(45*s +1):(46*s), (i-1)]  + tmp #S34
      
      # S321 to S31
      tmp <- prop_imm*(N_pop_ts[(111*s +1):(112*s), (i-1)]) #S321
      N_pop_ts[(111*s +1):(112*s), (i-1)] <- N_pop_ts[(111*s +1):(112*s), (i-1)]-tmp
      N_pop_ts[(43*s +1):(44*s), (i-1)]   <-  N_pop_ts[(43*s +1):(44*s), (i-1)]  + tmp #S31
      
      # S324 to S34
      tmp <- prop_imm*(N_pop_ts[(112*s +1):(113*s), (i-1)]) #S324
      N_pop_ts[(112*s +1):(113*s), (i-1)] <- N_pop_ts[(112*s +1):(113*s), (i-1)]-tmp
      N_pop_ts[(45*s +1):(46*s), (i-1)]   <-  N_pop_ts[(45*s +1):(46*s), (i-1)]  + tmp #S34
      
      
      
      
      # S412 to S41
      tmp <- prop_imm*(N_pop_ts[(115*s +1):(116*s), (i-1)]) #S412
      N_pop_ts[(115*s +1):(116*s), (i-1)] <- N_pop_ts[(115*s +1):(116*s), (i-1)]-tmp
      N_pop_ts[(46*s +1):(47*s), (i-1)]   <-  N_pop_ts[(46*s +1):(47*s), (i-1)]  + tmp #S41
      
      
      # S423 to S43
      tmp <- prop_imm*(N_pop_ts[(118*s +1):(119*s), (i-1)]) #S423
      N_pop_ts[(118*s +1):(119*s), (i-1)] <- N_pop_ts[(118*s +1):(119*s), (i-1)]-tmp
      N_pop_ts[(48*s +1):(49*s), (i-1)]   <-  N_pop_ts[(48*s +1):(49*s), (i-1)]  + tmp #S43
      
      # S421 to S41
      tmp <- prop_imm*(N_pop_ts[(117*s +1):(118*s), (i-1)]) #S421 
      N_pop_ts[(117*s +1):(118*s), (i-1)] <- N_pop_ts[(117*s +1):(118*s), (i-1)]-tmp
      N_pop_ts[(46*s +1):(47*s), (i-1)]   <-  N_pop_ts[(46*s +1):(47*s), (i-1)]  + tmp #S41
      
      # S432 to S43
      tmp <- prop_imm*(N_pop_ts[(120*s +1):(121*s), (i-1)]) #S432
      N_pop_ts[(120*s +1):(121*s), (i-1)] <- N_pop_ts[(120*s +1):(121*s), (i-1)]-tmp
      N_pop_ts[(48*s +1):(49*s), (i-1)]   <-  N_pop_ts[(48*s +1):(49*s), (i-1)]  + tmp #S43
      
      
      # Sm to S134, S143, S314, S341, S413, S431 (evenly)
      tmp <- prop_imm*(N_pop_ts[(146*s +1):(147*s), (i-1)]) #Sm 
      N_pop_ts[(146*s +1):(147*s), (i-1)] <- N_pop_ts[(146*s +1):(147*s), (i-1)] - tmp
      N_pop_ts[(100*s +1):(101*s), (i-1)]   <-  N_pop_ts[(100*s +1):(101*s), (i-1)]  + (tmp/6) #S134
      N_pop_ts[(102*s +1):(103*s), (i-1)]   <-  N_pop_ts[(102*s +1):(103*s), (i-1)]  + (tmp/6) #S143
      
      N_pop_ts[(110*s +1):(111*s), (i-1)]   <-  N_pop_ts[(110*s +1):(111*s), (i-1)]  + (tmp/6) #S314
      N_pop_ts[(113*s +1):(114*s), (i-1)]   <-  N_pop_ts[(113*s +1):(114*s), (i-1)]  + (tmp/6) #S341
      
      N_pop_ts[(116*s +1):(117*s), (i-1)]   <-  N_pop_ts[(116*s +1):(117*s), (i-1)]  + (tmp/6) #S413
      N_pop_ts[(119*s +1):(120*s), (i-1)]   <-  N_pop_ts[(119*s +1):(120*s), (i-1)]  + (tmp/6) #S431
      
      
      
      
      #then, let her go from there!
      
      Tmat <- buildTMat_four_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }
    
    
    
  }
  
  
  
  
  #the first column is the initial conditions. everything after that represents transitions within the year
  #should take the ceiling of each biweek to total within the year
  
  
  # stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  if(s>1){
    N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  }
  
  N.sim.dat <- as.data.frame(N_pop_ts)
  
  names(N.sim.dat) <- times
  
  N.sim.dat$class <- rep(c("S", "I1","I2", "I3", "I4",  "P1", "P2", "P3", "P4", "S1", "S2",  "S3", "S4",
                           "I12", "I13", "I14","I21", "I23", "I24", "I31", "I32", "I34", "I41", "I42", "I43",
                           "P12", "P13", "P14","P21", "P23", "P24", "P31", "P32", "P34",  "P41", "P42", "P43",
                           "S12", "S13", "S14","S21", "S23", "S24", "S31", "S32", "S34", "S41", "S42", "S43",
                           "I123","I124", "I132", "I134", "I142", "I143",
                           "I213", "I214","I231","I234", "I241", "I243",
                           "I312","I314", "I321", "I324", "I341", "I342",
                           "I412", "I413", "I421", "I423", "I431", "I432",
                           "P123","P124", "P132", "P134", "P142", "P143",
                           "P213", "P214","P231","P234", "P241", "P243",
                           "P312","P314", "P321", "P324", "P341", "P342",
                           "P412", "P413", "P421", "P423", "P431", "P432",
                           "S123","S124", "S132", "S134", "S142", "S143",
                           "S213", "S214","S231","S234", "S241", "S243",
                           "S312","S314", "S321", "S324", "S341", "S342",
                           "S412", "S413", "S421", "S423", "S431", "S432",
                           "I1234","I1243", "I1324", "I1342", "I1423", "I1432",
                           "I2134", "I2143","I2314","I2341", "I2413", "I2431",
                           "I3124","I3142", "I3214", "I3241", "I3412", "I3421",
                           "I4123", "I4132", "I4213", "I4231", "I4312", "I4321",
                           "Pms", "Pm"), 
                         each = s)
  N.sim.dat$age <- rep(1:s, length(unique(N.sim.dat$class)))
  
  #and melt
  N.sim.df <- melt(N.sim.dat, id=c("class", "age"), value.name = "count", variable.name = "time")
  
  N.sim.df$time <- as.numeric(as.character(N.sim.df$time))
  N.sim.df$count <- as.numeric(as.character(N.sim.df$count))
  
  #and add year
  #N.sim.df$year <- ceiling(N.sim.df$time)
  N.sim.df$year <- trunc(N.sim.df$time)
  #N.sim.df$time <- N.sim.df$time-1
  #N.sim.df$year <- N.sim.df$year-1
  
  
  
  return( N.sim.df)
}
sim.SIR.age.four.clim.wane.allsec <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, rate_hetero, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro, prop_imm){
  
  
  
  
  
  #check on climate vector - should be length of steps within a year (often 26)
  #(this could be modulated to be the length of the time series if you wanted)
  if(length(clim.vect)==1){
    clim.vect = rep(clim.vect, ntyr)
  }
  
  
  #first, if mort, births, sigma, or foi are shorter than the time series,
  #make them match in length here.
  
  #this will also repeat if they come in as a vector (age-structured)
  
  #make as a list
  if(length(foi)<yrs){
    foi <- rep(list(foi), yrs)
  }else{
    foi <- as.list(foi)
  }
  
  if(length(sigma)< yrs){
    sigma = rep(list(sigma), yrs)
  }else{
    sigma <- as.list(sigma)
  }
  
  # if(length(mort)<yrs){
  #   mort = rep(list(mort), yrs)
  # }else{
  #   mort = as.list(mort)
  # }
  # 
  if(length(births)<yrs){
    births = rep(list(births), yrs)
  }else{
    births = as.list(births)
  }
  
  # 
  if(length(rate_hetero)<yrs){
    rate_hetero = rep(list(rate_hetero), yrs)
  }else{
    rate_hetero = as.list(rate_hetero)
  }
  
  
  #
  #number of epidemic classes
  c=147
  
  
  #length of time series
  times <-   seq(0, yrs, by =1/ntyr)
  times <-   times + (year.end-yrs)
  
  #split the births up by biweek too
  birth_vector_biwk = as.list(c(unlist(lapply(lapply(births, divide, ntyr1 =ntyr), rep, ntyr))))
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  #import the stable age structure from the literature
  #mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  #stab.struct = Re(eigen(mat1)$vector[,1])
  #stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  #lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  pop.mat = pop_vector #number of people per age class
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  #structure the population for 3 strains even though we only have 2 to start. they will get added later
  
  PM_init = rep(0, s)
  
  PMs_init = rep(0, s)
  
  I4321_init = rep(0, s)
  I4312_init = rep(0, s)
  I4231_init = rep(0, s)
  I4213_init = rep(0, s)
  I4132_init = rep(0, s)
  I4123_init = rep(0, s)
  
  I3421_init = rep(0, s)
  I3412_init = rep(0, s)
  I3241_init = rep(0, s)
  I3214_init = rep(0, s)
  I3142_init = rep(0, s)
  I3124_init = rep(0, s)
  
  I2431_init = rep(0, s)
  I2413_init = rep(0, s)
  I2341_init = rep(0, s)
  I2314_init = rep(0, s)
  I2143_init = rep(0, s)
  I2134_init = rep(0, s)
  
  I1432_init = rep(0, s)
  I1423_init = rep(0, s)
  I1342_init = rep(0, s)
  I1324_init = rep(0, s)
  I1243_init = rep(0, s)
  I1234_init = rep(0, s)
  
  
  S432_init = rep(0, s)
  S431_init = rep(0, s)
  S423_init = rep(0, s)
  S421_init = rep(0, s)
  S413_init = rep(0, s)
  S412_init = rep(0, s)
  
  S342_init = rep(0, s)
  S341_init = rep(0, s)
  S324_init = rep(0, s)
  S321_init = rep(0, s)
  S314_init = rep(0, s)
  S312_init = rep(0, s)
  
  S243_init = rep(0, s)
  S241_init = rep(0, s)
  S234_init = rep(0, s)
  S231_init = rep(0, s)
  S214_init = rep(0, s)
  S213_init = rep(0, s)
  
  S143_init = rep(0, s)
  S142_init = rep(0, s)
  S134_init = rep(0, s)
  S132_init = rep(0, s)
  S124_init = rep(0, s)
  S123_init = rep(0, s)
  
  
  P432_init = rep(0, s)
  P431_init = rep(0, s)
  P423_init = rep(0, s)
  P421_init = rep(0, s)
  P413_init = rep(0, s)
  P412_init = rep(0, s)
  
  P342_init = rep(0, s)
  P341_init = rep(0, s)
  P324_init = rep(0, s)
  P321_init = rep(0, s)
  P314_init = rep(0, s)
  P312_init = rep(0, s)
  
  P243_init = rep(0, s)
  P241_init = rep(0, s)
  P234_init = rep(0, s)
  P231_init = rep(0, s)
  P214_init = rep(0, s)
  P213_init = rep(0, s)
  
  P143_init = rep(0, s)
  P142_init = rep(0, s)
  P134_init = rep(0, s)
  P132_init = rep(0, s)
  P124_init = rep(0, s)
  P123_init = rep(0, s)
  
  
  I432_init = rep(0, s)
  I431_init = rep(0, s)
  I423_init = rep(0, s)
  I421_init = rep(0, s)
  I413_init = rep(0, s)
  I412_init = rep(0, s)
  
  I342_init = rep(0, s)
  I341_init = rep(0, s)
  I324_init = rep(0, s)
  I321_init = rep(0, s)
  I314_init = rep(0, s)
  I312_init = rep(0, s)
  
  I243_init = rep(0, s)
  I241_init = rep(0, s)
  I234_init = rep(0, s)
  I231_init = rep(0, s)
  I214_init = rep(0, s)
  I213_init = rep(0, s)
  
  I143_init = rep(0, s)
  I142_init = rep(0, s)
  I134_init = rep(0, s)
  I132_init = rep(0, s)
  I124_init = rep(0, s)
  I123_init = rep(0, s)
  
  S43_init = rep(0, s)
  S42_init = rep(0, s)
  S41_init = rep(0, s)
  S34_init = rep(0, s)
  S32_init = rep(0, s)
  S31_init = rep(0, s)
  S24_init = rep(0, s)
  S23_init = rep(0, s)
  S21_init = rep(0, s)
  S14_init = rep(0, s)
  S13_init = rep(0, s)
  S12_init = rep(0, s)
  
  P43_init = rep(0, s)
  P42_init = rep(0, s)
  P41_init = rep(0, s)
  P34_init = rep(0, s)
  P32_init = rep(0, s)
  P31_init = rep(0, s)
  P24_init = rep(0, s)
  P23_init = rep(0, s)
  P21_init = rep(0, s)
  P14_init = rep(0, s)
  P13_init = rep(0, s)
  P12_init = rep(0, s)
  
  I43_init = rep(0, s)
  I42_init = rep(0, s)
  I41_init = rep(0, s)
  I34_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I24_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I14_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  
  S4_init = rep(0, s)
  S3_init = rep(0, s)
  S2_init = rep(0, s)
  S1_init = rep(0, s)
  
  P4_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  
  I4_init = rep(0, s); I4_init[1] = 5 #comment out if you just want to check demography
  I3_init = rep(0, s); I3_init[1] = 5 #comment out if you just want to check demography
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  
  S_init = pop.mat - I1_init - I2_init - I3_init - I4_init
  
  
  N_tot = cbind(S_init, 
                I1_init, I2_init,I3_init,I4_init,
                P1_init, P2_init,P3_init,P4_init, 
                S1_init,S2_init,S3_init,S4_init,
                I12_init, I13_init, I14_init, 
                I21_init, I23_init, I34_init, 
                I31_init, I32_init, I34_init,
                I41_init, I42_init, I43_init,
                P12_init, P13_init, P14_init, 
                P21_init, P23_init, P34_init, 
                P31_init, P32_init, P34_init,
                P41_init, P42_init, P43_init,
                S12_init, S13_init, S14_init, 
                S21_init, S23_init, S34_init, 
                S31_init, S32_init, S34_init,
                S41_init, S42_init, S43_init,
                I143_init,I142_init,I134_init,I132_init,I124_init,I123_init,
                I213_init,I214_init, I231_init, I234_init, I241_init, I243_init,
                I312_init, I314_init, I321_init, I324_init, I341_init, I342_init,
                I412_init, I413_init, I421_init, I423_init,I431_init, I432_init, 
                P143_init,P142_init,P134_init,P132_init,P124_init,P123_init,
                P213_init,P214_init, P231_init, P234_init, P241_init, P243_init,
                P312_init, P314_init, P321_init, P324_init, P341_init, P342_init,
                P412_init, P413_init, P421_init, P423_init,P431_init, P432_init, 
                S143_init,S142_init,S134_init,S132_init,S124_init,S123_init,
                S213_init,S214_init, S231_init, S234_init, S241_init, S243_init,
                S312_init, S314_init, S321_init, S324_init, S341_init, S342_init,
                S412_init, S413_init, S421_init, S423_init,S431_init, S432_init,  
                I1234_init, I1243_init, I1324_init, I1342_init,I1423_init,I1432_init,
                I2134_init, I2143_init, I2314_init, I2341_init,I2413_init, I2431_init, 
                I3124_init, I3142_init, I3214_init, I3241_init, I3412_init, I3421_init, 
                I4123_init, I4132_init, I4213_init, I4231_init, I4312_init, I4321_init,
                PMs_init,
                PM_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  
  N_pop_ts[,1] <- M_pop # by age class
  
  #and fill in for births <- every epidemic class can give birth
  #transform birth vector into one distributed by age and epidemic class
  
  
  # stab.struct <- get.age.struct.M(M_pop, c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  
  # foi gets spread across the biweeks as well
  # (it comes in as 1 value per year)
  # take your Muench-estimated rates of annual infection and convert
  # to biweekly probabilities of infection
  foi.biweek = as.list(c(unlist(lapply(lapply(foi, divide.rate, ntyr1 =ntyr), rep, ntyr)))) #here, it could easily be modulated to be seasonal - (including forced by precip/temp)
  
  
  #and make annual waning immunity rate into probability of waning immunity per biweek
  sigma.biweek = as.list(c(unlist(lapply(lapply(sigma, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #and the probability of waning heterotypic immunity per biweek
  hetero.biweek = as.list(c(unlist(lapply(lapply(rate_hetero, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(hetero.biweek)!=(length(times)-1) |length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
    print("input vectors of unequal length")
  }
  #and your aging probability to biweeks (aging probability per year is 100%)
  #age.biwk = rep(1-exp(-(rep(1, yrs)/ntyr)), each=26)
  
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 2:(length(times))){
    #print(times[i])
    
    
    
    #here are the dynamics before and after the immune waning event
    if((times[i]<((yr.intro) + biwk.intro/26)) | (times[i]>((yr.intro) + biwk.intro/26))){ 
      
      
      
      biwk1 <- find.biweek(t=i, times=times)
      clim.mod <- clim.vect[biwk1]
      
      #and select the appropriate age vector, based on the year
      age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
      
      
      
      if(nrow(age.mult.sub)>11){
        print("error with age multiplier sub-selection")
      }
      
      age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
      age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
      age.mult.sub$dur[age.mult.sub$dur<0] <- 0
      age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
      
      
      #and the mortality rates
      if(times[i]<min(mort$year)){
        mort.df = subset(mort, year==min(year))
      }else{
        mort.df = subset(mort, year==trunc(times[i]))
      }
      #make annual mortality rate by age into probability of mortality per biweek
      #these are deaths per 1000 people.
      #we multiply by the number in the age classes, so need to divide by 1000 here
      mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
      mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
      
      
      #then, let her go from there!
      
      Tmat <- buildTMat_four_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }else if(times[i]==((yr.intro) + biwk.intro/26)){
      
      
      
      biwk1 <- find.biweek(t=i, times=times)
      clim.mod <- clim.vect[biwk1]
      
      #and select the appropriate age vector, based on the year
      age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
      
      
      
      if(nrow(age.mult.sub)>11){
        print("error with age multiplier sub-selection")
      }
      
      age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
      age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
      age.mult.sub$dur[age.mult.sub$dur<0] <- 0
      age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
      
      
      #and the mortality rates
      if(times[i]<min(mort$year)){
        mort.df = subset(mort, year==min(year))
      }else{
        mort.df = subset(mort, year==trunc(times[i]))
      }
      #make annual mortality rate by age into probability of mortality per biweek
      #these are deaths per 1000 people.
      #we multiply by the number in the age classes, so need to divide by 1000 here
      mort.df$deaths_per_cap <- mort.df$deaths_per_1000_per_age/1000
      mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
      
      
      # in this iteration, it is a 4-strain matrix, but we arbitrarily move some proportion of those with previous 
      # exposure to DENV-2 to a state of naivity for that seroptype at the beginning of the year of introduction
      # so any previous "S" class that had a 2 now gets moved to an S class without a 2
      
      #need to find the index that corresponds to the disease state of interest
      #vertical indices are the length of the number of epidemic classes (c) x the number of age classes (s)
      #matrix is stacked with all the ages of a given disease status vertically on top of one another
      #so [1,1] is S age 1, and [2,1] is S age 2
      
      #list the first index in each disease class
      
      
      
      # S2 to S
      #tmp <- prop_imm*(N_pop_ts[(10*s +1):(11*s), (i-1)]) #S2
      #N_pop_ts[(10*s +1):(11*s), (i-1)] <- N_pop_ts[(10*s +1):(11*s), (i-1)]-tmp
      #N_pop_ts[(1):(1*s), (i-1)] <-  N_pop_ts[(1):(1*s), (i-1)] + tmp #S
      
      
      # S12 to S1
      tmp <- prop_imm*(N_pop_ts[(37*s +1):(38*s), (i-1)]) #S12
      N_pop_ts[(37*s +1):(38*s), (i-1)] <- N_pop_ts[(37*s +1):(38*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + tmp #S1
      
      
      # S32 to S3
      tmp <- prop_imm*(N_pop_ts[(44*s +1):(45*s), (i-1)]) #S32
      N_pop_ts[(44*s +1):(45*s), (i-1)] <- N_pop_ts[(44*s +1):(45*s), (i-1)]-tmp
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + tmp #S3
      
      # S42 to S4
      tmp <- prop_imm*(N_pop_ts[(47*s +1):(48*s), (i-1)]) #S42
      N_pop_ts[(47*s +1):(48*s), (i-1)] <- N_pop_ts[(47*s +1):(48*s), (i-1)]-tmp
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + tmp #S4
      
      
      
      # S123 to half S1 and half S3
      tmp <- prop_imm*(N_pop_ts[(97*s +1):(98*s), (i-1)]) #S123
      N_pop_ts[(97*s +1):(98*s), (i-1)] <- N_pop_ts[(97*s +1):(98*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      
      
      
      # S124 to half S1 and half S4
      tmp <- prop_imm*(N_pop_ts[(98*s +1):(99*s), (i-1)]) #S124
      N_pop_ts[(98*s +1):(99*s), (i-1)] <- N_pop_ts[(98*s +1):(99*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      
      # S132 to half S1 and  and half S3
      tmp <- prop_imm*(N_pop_ts[(99*s +1):(100*s), (i-1)]) #S132
      N_pop_ts[(99*s +1):(100*s), (i-1)] <- N_pop_ts[(99*s +1):(100*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      
      
      # S142 to half S1 and half S4
      tmp <- prop_imm*(N_pop_ts[(101*s +1):(102*s), (i-1)]) #S142
      N_pop_ts[(101*s +1):(102*s), (i-1)] <- N_pop_ts[(101*s +1):(102*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      # S213 to half S1 and half S3
      tmp <- prop_imm*(N_pop_ts[(103*s +1):(104*s), (i-1)]) #S213
      N_pop_ts[(103*s +1):(104*s), (i-1)] <- N_pop_ts[(103*s +1):(104*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      
      
      # S214 to half S1  and half S4
      tmp <- prop_imm*(N_pop_ts[(104*s +1):(105*s), (i-1)]) #S214
      N_pop_ts[(104*s +1):(105*s), (i-1)] <- N_pop_ts[(104*s +1):(105*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      # S231 to half S3  and half S1
      tmp <- prop_imm*(N_pop_ts[(105*s +1):(106*s), (i-1)]) #S231
      N_pop_ts[(105*s +1):(106*s), (i-1)] <- N_pop_ts[(105*s +1):(106*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      
      
      # S234 to half S3 and half S4
      tmp <- prop_imm*(N_pop_ts[(106*s +1):(107*s), (i-1)]) #S234
      N_pop_ts[(106*s +1):(107*s), (i-1)] <- N_pop_ts[(106*s +1):(107*s), (i-1)]-tmp
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      
      
      # S241 to half S4 and half S1
      tmp <- prop_imm*(N_pop_ts[(107*s +1):(108*s), (i-1)]) #S241
      N_pop_ts[(107*s +1):(108*s), (i-1)] <- N_pop_ts[(107*s +1):(108*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      # S243 to half S4 and half S3
      tmp <- prop_imm*(N_pop_ts[(108*s +1):(109*s), (i-1)]) #S243
      N_pop_ts[(108*s +1):(109*s), (i-1)] <- N_pop_ts[(108*s +1):(109*s), (i-1)]-tmp
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      # S312 to half S3 and half S1
      tmp <- prop_imm*(N_pop_ts[(109*s +1):(110*s), (i-1)]) #S312
      N_pop_ts[(109*s +1):(110*s), (i-1)] <- N_pop_ts[(109*s +1):(110*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      
      
      # S342 to half S3 and half S4
      tmp <- prop_imm*(N_pop_ts[(114*s +1):(115*s), (i-1)]) #S342
      N_pop_ts[(114*s +1):(115*s), (i-1)] <- N_pop_ts[(114*s +1):(115*s), (i-1)]-tmp
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      # S321 to half S3 and half S1
      tmp <- prop_imm*(N_pop_ts[(111*s +1):(112*s), (i-1)]) #S321
      N_pop_ts[(111*s +1):(112*s), (i-1)] <- N_pop_ts[(111*s +1):(112*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      
      
      # S324 to half S3 and half S4
      tmp <- prop_imm*(N_pop_ts[(112*s +1):(113*s), (i-1)]) #S324
      N_pop_ts[(112*s +1):(113*s), (i-1)] <- N_pop_ts[(112*s +1):(113*s), (i-1)]-tmp
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      # S412 to half S4 and half S1
      tmp <- prop_imm*(N_pop_ts[(115*s +1):(116*s), (i-1)]) #S412
      N_pop_ts[(115*s +1):(116*s), (i-1)] <- N_pop_ts[(115*s +1):(116*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      # S423 to half S4 and half S3
      tmp <- prop_imm*(N_pop_ts[(118*s +1):(119*s), (i-1)]) #S423
      N_pop_ts[(118*s +1):(119*s), (i-1)] <- N_pop_ts[(118*s +1):(119*s), (i-1)]-tmp
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      # S421 to half S4 and half S1
      tmp <- prop_imm*(N_pop_ts[(117*s +1):(118*s), (i-1)]) #S421 
      N_pop_ts[(117*s +1):(118*s), (i-1)] <- N_pop_ts[(117*s +1):(118*s), (i-1)]-tmp
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/2) #S1
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      # S432 to half S4 and half S3
      tmp <- prop_imm*(N_pop_ts[(120*s +1):(121*s), (i-1)]) #S432
      N_pop_ts[(120*s +1):(121*s), (i-1)] <- N_pop_ts[(120*s +1):(121*s), (i-1)]-tmp
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/2) #S3
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/2) #S4
      
      
      
      # Sm to S1, S3, S4 (evenly)
      tmp <- prop_imm*(N_pop_ts[(146*s +1):(147*s), (i-1)]) #Sm 
      N_pop_ts[(9*s+1):(10*s), (i-1)] <-  N_pop_ts[(9*s+1):(10*s), (i-1)] + (tmp/3) #S1
      N_pop_ts[(11*s +1):(12*s), (i-1)]  <-  N_pop_ts[(11*s +1):(12*s), (i-1)]  + (tmp/3) #S3
      N_pop_ts[(12*s +1):(13*s), (i-1)]  <-  N_pop_ts[(12*s +1):(13*s), (i-1)]  + (tmp/3) #S4
      
      
      
      #then, let her go from there!
      
      Tmat <- buildTMat_four_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }
    
    
    
  }
  
  
  
  
  #the first column is the initial conditions. everything after that represents transitions within the year
  #should take the ceiling of each biweek to total within the year
  
  
  # stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  if(s>1){
    N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  }
  
  N.sim.dat <- as.data.frame(N_pop_ts)
  
  names(N.sim.dat) <- times
  
  N.sim.dat$class <- rep(c("S", "I1","I2", "I3", "I4",  "P1", "P2", "P3", "P4", "S1", "S2",  "S3", "S4",
                           "I12", "I13", "I14","I21", "I23", "I24", "I31", "I32", "I34", "I41", "I42", "I43",
                           "P12", "P13", "P14","P21", "P23", "P24", "P31", "P32", "P34",  "P41", "P42", "P43",
                           "S12", "S13", "S14","S21", "S23", "S24", "S31", "S32", "S34", "S41", "S42", "S43",
                           "I123","I124", "I132", "I134", "I142", "I143",
                           "I213", "I214","I231","I234", "I241", "I243",
                           "I312","I314", "I321", "I324", "I341", "I342",
                           "I412", "I413", "I421", "I423", "I431", "I432",
                           "P123","P124", "P132", "P134", "P142", "P143",
                           "P213", "P214","P231","P234", "P241", "P243",
                           "P312","P314", "P321", "P324", "P341", "P342",
                           "P412", "P413", "P421", "P423", "P431", "P432",
                           "S123","S124", "S132", "S134", "S142", "S143",
                           "S213", "S214","S231","S234", "S241", "S243",
                           "S312","S314", "S321", "S324", "S341", "S342",
                           "S412", "S413", "S421", "S423", "S431", "S432",
                           "I1234","I1243", "I1324", "I1342", "I1423", "I1432",
                           "I2134", "I2143","I2314","I2341", "I2413", "I2431",
                           "I3124","I3142", "I3214", "I3241", "I3412", "I3421",
                           "I4123", "I4132", "I4213", "I4231", "I4312", "I4321",
                           "Pms", "Pm"), 
                         each = s)
  N.sim.dat$age <- rep(1:s, length(unique(N.sim.dat$class)))
  
  #and melt
  N.sim.df <- melt(N.sim.dat, id=c("class", "age"), value.name = "count", variable.name = "time")
  
  N.sim.df$time <- as.numeric(as.character(N.sim.df$time))
  N.sim.df$count <- as.numeric(as.character(N.sim.df$count))
  
  #and add year
  #N.sim.df$year <- ceiling(N.sim.df$time)
  N.sim.df$year <- trunc(N.sim.df$time)
  #N.sim.df$time <- N.sim.df$time-1
  #N.sim.df$year <- N.sim.df$year-1
  
  
  
  return( N.sim.df)
}

check.equil <- function(dat){
  dat.ts <- ddply(dat, .(time, class), summarise, count=sum(count))
  #and get total by time
  dat.N  <- ddply(dat, .(time), summarise, Ntot=sum(count))
  
  dat.ts <- merge(dat.ts, dat.N, by="time")
  dat.ts$proportion <- dat.ts$count/dat.ts$Ntot
  
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=proportion, color=class)) #+ coord_cartesian(ylim=c(0,.1))
  print(p1)
  return(dat.ts)
}
plot.cases.annual <- function(dat, year.start){
  dat1 = subset(dat, year>=year.start)
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  #dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts <- ddply(denv.case, .(year), summarise, count=sum(count))
  dat.ts$tot_count <- dat.ts$count#/.001
  #dat.ts$tot_count <- dat.ts$count/.001
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count)) #+ coord_cartesian(ylim=c(0,.1))
  p1 <- ggplot(dat.ts) + theme_bw() +
        geom_vline(aes(xintercept=2007), linetype=2)+
        geom_vline(aes(xintercept=1990), linetype=2)+
        geom_vline(aes(xintercept=2012), linetype=2)+
        geom_vline(aes(xintercept=2019), linetype=2) +
        geom_line(aes(x=year, y=tot_count)) + 
          theme(panel.grid = element_blank())
  print(p1)
  return(dat.ts)
}
plot.cases.annual.triple.type <- function(dat, perc.obs, count.type, year.start){
  dat1 = subset(dat, year>=year.start)

  
  #select only those viewed as "cases"
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32" | class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  denv.case$type = "secondary"
  
  if(count.type=="tertiary"){
    
    denv.case.tert = subset(dat1, class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
    if(length(unique(denv.case.tert$age))>1){
      denv.case.tert = subset(denv.case.tert, age<max(denv.case$age))  
    }
    
    denv.case.tert$type <- "tertiary"  
    denv.case.tert$count <- denv.case.tert$count*perc.obs
    denv.case <- rbind( denv.case,  denv.case.tert)
    
  }
  
  #dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts <- ddply(denv.case, .(year), summarise, count=sum(count))
  dat.ts$tot_count <- dat.ts$count#/.001
  #dat.ts$tot_count <- dat.ts$count/.001
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count)) #+ coord_cartesian(ylim=c(0,.1))
  p1 <- ggplot(dat.ts) + theme_bw() +
    geom_vline(aes(xintercept=2007), linetype=2)+
    geom_vline(aes(xintercept=1990), linetype=2)+
    geom_vline(aes(xintercept=2012), linetype=2)+
    geom_vline(aes(xintercept=2019), linetype=2) +
    geom_line(aes(x=year, y=tot_count)) + 
    theme(panel.grid = element_blank())
  print(p1)
  return(dat.ts)
}
cum.sum.year <- function(df){
  df.sum <- ddply(df,.(age), summarise, cases=sum(count))
  
  df.sum$cum_cases = cumsum(df.sum$cases)
  df.sum$cum_prop_cases <- df.sum$cum_cases/sum(df.sum$cases)
  df.sum$year <- unique(df$year)
  return(df.sum)
}
process.all <- function(df){
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "S", "I1","P1","I2", "P2","I3","I4", "PM")
  dat.tot$N = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  #par(mfrow = c(1,1))
  #with(dat.tot, plot(time, N, type="l")) #ylim =c(0,1.2*max(N))))
  #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  
  # prop.tot = dat.tot
  # 
  # prop.tot$S = prop.tot$S/prop.tot$N
  # prop.tot$I = prop.tot$Ii/prop.tot$N
  # prop.tot$R = prop.tot$R/prop.tot$N
  # 
  
  
  
  
  dat.melt = melt(dat.tot,id.vars = "time")
  
  #head(dat.melt)
  
  #take after burnin if we assume this is a virus at equilibrium
  dat.melt <- subset(dat.melt, time>length(burnin_lambda))
  dat.melt <- subset(dat.melt, !is.na(value))
  
  #dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5])))
  names(dat.melt) = c("time", "class", "count")
  dat.N = subset(dat.melt, class =="N")
  dat.melt = subset(dat.melt, class !="N")
  dat.N <- dplyr::select(dat.N, -(class))
  names(dat.N)[names(dat.N)=="count"] <- "Ntot"
  
  dat.melt <- merge(dat.melt, dat.N, by="time")
  #head(dat.melt)
  #dat.tot$class = rep(c("M", "S", "I", "R"), each= length(times))
  dat.melt$class = factor(dat.melt$class, levels=c("S", "I1","P1","I2", "P2","I3","I4", "PM"))
  dat.melt$proportion <- dat.melt$count/dat.melt$Ntot
  
  dat.melt$time <- dat.melt$time - length(burnin_lambda)
  # colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
  #ggplot(data=dat.melt) + geom_line(aes(x=time, y=count, color=class)) #+ scale_color_manual(values=colz)
  
  #dat.melt$proportion = dat.melt$count/sim_pop
  #ggplot(data=dat.melt) + geom_line(aes(x=time, y=proportion, color=class)) #+ scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  #now, get the age-structured incidence - incidence is "I3 + I4"
  
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("S", "I1","P1","I2", "P2","I3","I4", "PM"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  
  #plot age structured cumulative incidence
  age.incidence = subset(age.dat.tot, class=="I3" | class == "I4")
  
  #and collect everyone else
  age.non = subset(age.dat.tot, class!="I3" & class!="I4")
  #head(age.incidence )
  
  age.incidence$year = trunc(age.incidence$times)
  age.non$year = trunc(age.non$times)
  age.incidence <- arrange(age.incidence, year, age)
  
  age.non <- arrange(age.non, year, age)
  
  names(age.non) <- c("times", "count_non_infectious", "class", "age", "year")
  
  age.non <- dplyr::select(age.non, -(class), -(times))
  
  age.non <- ddply(age.non, .(year, age), summarise, count_non_infectious = sum(count_non_infectious))
  
  
  age.sum <- ddply(age.incidence, .(year, age), summarise, count = sum(count))
  
  age.sum <- merge(age.sum, age.non, by=c("year", "age"))
  
  #and just return this counts of years, ages, and those infected and not
  
  
  return(age.sum)#, seas.prev))
  
}

# first, sim Cambodia 2 strain to endemic equilibrium for 50 yrs prior to 1960,
# then sim with birth and death rates (and starting at the correct population distribution) 
# through to 2020. Assess how well this recapitulates what we find in the data


#load climate from TSIR first:
clim.dat <- read.csv(file = "beta_TSIR_fit_province.csv", header = T, stringsAsFactors = F)
head(clim.dat)

#then, get the median beta across epidemic periods and provinces
#we just want this to be an input of seasonality to our transmission rate
beta.med <- ddply(clim.dat,.(biweek), summarise, beta = median(beta))
beta.med$clim_vect <- scales::rescale(beta.med$beta, to=c(.5,1.5), from=c(range(beta.med$beta)))



#load birth and death rates for Cambodia
pop.dat <- read.csv(file = "pop_data_full.csv", header = T, stringsAsFactors = F)


birth.dat = subset(pop.dat, metric=="births per\n1000 ppl")
death.dat = subset(pop.dat, metric == "deaths per\n1000 ppl")


mort.dat <-  read.csv(file = "cambodia_age_specific_mort_through_time.csv", header = T, stringsAsFactors = F)
names(mort.dat) <- c("year", seq(0,100,1))
mort.melt <- melt(mort.dat, id.vars = "year", variable.name = "age", value.name = "count")
mort.melt = arrange(mort.melt, year, age)
mort.melt$age <- as.numeric(as.character(mort.melt$age))
mort.melt$year <- as.numeric(as.character(mort.melt$year))

mort.sum <- ddply(mort.melt,.(year), summarise, total_deaths = sum(count))
mort.melt <- merge(mort.melt, mort.sum, by="year", all.x = T, sort = F)
mort.melt$perc_deaths_per_age <- mort.melt$count/mort.melt$total_deaths

#and all deaths
mort.melt <- merge(mort.melt, death.dat, by="year", all.x = T)
head(mort.melt)
mort.melt$value[is.na(mort.melt$value)] <- mort.melt$value[mort.melt$year==1960]
mort.melt$metric[is.na(mort.melt$metric)] <- mort.melt$metric[mort.melt$year==1960]

mort.melt$deaths_per_1000_per_age <- mort.melt$perc_deaths_per_age*mort.melt$value

mort.melt <- dplyr::select(mort.melt, year, age, deaths_per_1000_per_age)


pop.dist <- read.csv(file = "cambodia_pop_dist_through_time.csv", header = T, stringsAsFactors = F)
head(pop.dist)
names(pop.dist) <- c("year", seq(0,100,1))
pop.melt <- melt(pop.dist, id.vars = "year", variable.name = "age", value.name = "count")
pop.melt = arrange(pop.melt, year, age)
pop.melt$age <- as.numeric(as.character(pop.melt$age))
pop.melt$year <- as.numeric(as.character(pop.melt$year))

#and replicate mortality to match the length?
# year.sub = subset(mort.melt, year == min(mort.melt$year))
# year.rep = 

#load the foi fits for cambodia
fit.dat <- read.csv(file = "prov-fits-FOI.csv", stringsAsFactors = F, header = T)
nrow(fit.dat[fit.dat$provname=="National",]) #40 years - run for 40 before this as burnin

#load age structure
age.fit <- read.csv(file = "age-mult-profile.csv", stringsAsFactors = F, header = T)
age.fit$year_min <-   0 #for simulation
age.fit$year_max <-   100 #for simulation
#age.fit$year_max[age.fit$year_range=="<=2010"] <- 2010
age.fit$year_min[age.fit$year_range=="2002-2010"] <- 0
age.fit$year_max[age.fit$year_range=="2002-2010"] <- 2010.999999
age.fit$year_min[age.fit$year_range=="2011-2020"] <- 2011
age.fit$year_max[age.fit$year_range=="2011-2020"] <- 2025

age.fit = subset(age.fit, age_max <=4)
pop.melt = subset(pop.melt, age<=3)
#first, just sim normal
#sim here, using foi at the National level, but replacing the too-low values
fit.dat$lambda[fit.dat$provname=="National" &fit.dat$year<1999] <- .9
mort.melt = subset(mort.melt, age<=3)

#run different scenarios
# you wane one step back in Abs to just one target serotype BUT you can't go back to naive immunity
# but everything wanes back to only to monotypic immunity so you re-experience a secondary infection
out.cam.wane.test.intro = sim.SIR.age.four.clim.wane.intro(yrs=10,
                                              ntyr=26,
                                              s=4, 
                                              foi=fit.dat$lambda[fit.dat$provname=="National"][1:10], 
                                              births =  birth.dat$value[52:61], # these are per 1000
                                              pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                              recov=1,
                                              age.mult.df=age.fit, 
                                              clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                              age.brk=NA,#this is a placeholder for additional age structure
                                              mort=mort.melt, 
                                              rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                              sigma=0, #nothing at first
                                              yr.intro = 2015,
                                              year.end = 2021,
                                              biwk.intro = 1)

#you wane one step back in Abs to just one target serotype BUT you can't go back to naive immunity, only to monotypic or duotypic, tritypic
out.cam.wane.test.sec = sim.SIR.age.four.clim.wane.sec(yrs=40,
                                                      ntyr=26,
                                                      s=6, 
                                                      foi=fit.dat$lambda[fit.dat$provname=="National"][1:40], 
                                                      births =  birth.dat$value[22:61], # these are per 1000
                                                      pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                      recov=1,
                                                      age.mult.df=age.fit, 
                                                      clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                              age.brk=NA,#this is a placeholder for additional age structure
                                                              mort=mort.melt, 
                                                              rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                              sigma=0, #nothing at first
                                                              yr.intro = 2015,
                                                              year.end = 2021,
                                                              biwk.intro = 1,
                                                              prop_imm = .8)


out.cam.wane.test.normal = sim.SIR.age.four.clim.wane(yrs=40,
                                                       ntyr=26,
                                                       s=6, 
                                                       foi=fit.dat$lambda[fit.dat$provname=="National"][1:40], 
                                                       births =  birth.dat$value[22:61], # these are per 1000
                                                       pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                       recov=1,
                                                       age.mult.df=age.fit, 
                                                       clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                       age.brk=NA,#this is a placeholder for additional age structure
                                                       mort=mort.melt, 
                                                       rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                       sigma=0, #nothing at first
                                                       yr.intro = 2015,
                                                       year.end = 2021,
                                                       biwk.intro = 1,
                                                       prop_imm = .8)



out.cam.maintain = sim.SIR.age.four.clim.maintain(yrs=40,
                                                       ntyr=26,
                                                       s=6, 
                                                       foi=fit.dat$lambda[fit.dat$provname=="National"][1:40], 
                                                       births =  birth.dat$value[22:61], # these are per 1000
                                                       pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                       recov=1,
                                                       age.mult.df=age.fit, 
                                                       clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                       age.brk=NA,#this is a placeholder for additional age structure
                                                       mort=mort.melt, 
                                                       rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                       sigma=0, #nothing at first
                                                       #yr.intro = 2015,
                                                       year.end = 2021)
                                                       #biwk.intro = 1,
                                                       #prop_imm = .8)

#save(out.cam.wane.test.sec.all, file = "cam-sim-4-sero-wane-test-sec-all.Rdata")

#load("cam-sim-4-sero-wane-test-sec-all.Rdata")

library(stringr)

summarise.age.dist <- function(dat, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #summarise into primary, secondary, tertiary, quaternary cases
  
  dat1$case_type <- (str_extract(dat1$class, "[aA-zZ]+"))
  dat1$case_sum <- nchar(dat1$class)
  dat1$state <- NA  
  dat1$state[dat1$case_type=="S"] <- "S"
  dat1$state[dat1$case_type=="P" | dat1$case_type=="Pms"] <- "Temporary-Heterotypic-Immunity"
  dat1$state[dat1$case_type=="Pm" ] <- "Pm"
  dat1$state[dat1$case_type=="I" & dat1$case_sum==2] <- "Primary-Infection"
  dat1$state[dat1$case_type=="I" & dat1$case_sum==3] <- "Secondary-Infection"
  dat1$state[dat1$case_type=="I" & dat1$case_sum==4] <- "Tertiary-Infection"
  dat1$state[dat1$case_type=="I" & dat1$case_sum==5] <- "Quaternary-Infection"
  dat1$state[dat1$case_type=="I" & dat1$case_sum==6] <- "Fifth-Infection"
  
  #tracking reinfections with serotype 2
  
  #and sum by year
  df.sum <- ddply(dat1,.(year, age, state), summarise, count = sum(count))
  
  
  
  df.sum$count<- round(df.sum$count,0)
  
  
  df.sum <- df.sum[complete.cases(df.sum),]
  
  #and just focus on infections
  #df.sum.1 <- df.sum
  
  #and return these
  df.sum.I = subset(df.sum, state!="Pm" & state!="Temporary-Heterotypic-Immunity" & state!="S")
  return(df.sum.I)
}
out.age.wane.sec.all <- summarise.age.dist(dat = out.cam.wane.test.sec.all, year.start = 1981)
out.age.wane.sec <- summarise.age.dist(dat = out.cam.wane.test.sec, year.start = 1981)
out.age.wane <- summarise.age.dist(dat = out.cam.wane.test.normal, year.start = 1981)
out.age.maintain <- summarise.age.dist(dat = out.cam.maintain, year.start = 1981)

#save(out.age.wane, file = "cam-age-sum-wane-4-sero.Rdata")  


#and plot cases
plot.age.dist <- function(df, slim.quant, view.plot, save.plot, filename){
  #split by a year
  df.year <- dlply(df,.(year, state))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year and type
  df.age <- dlply(df,.(state, year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data.type, slim.quant=slim.quant))
  
  dat.age$state <- factor(dat.age$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection"))
  mean.df$state <- factor(mean.df$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection"))
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  #colz = c("secondary"="black", "tertiary"="royalblue3")
  
  p1 <- ggplot(dat.age) + facet_grid(~state) +
    geom_jitter(aes(x=year, y=age, color=state), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    geom_violin(aes(x=year,y=age, group=year), color="gray60",  draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") + #scale_color_manual(values=colz) +
    geom_hline(aes(yintercept=1), color="red") #+ coord_cartesian(ylim=c(0,80))#,xlim=c(2015,2020))
  
  
  if(save.plot==TRUE){
    ggsave(file = filename,
           plot= p1,
           units="mm",  
           width=100, 
           height=55, 
           scale=3, 
           dpi=300)
    
  }
  if(view.plot==TRUE){
    print(p1)
  }
  
  
  
  return(mean.df)
}
plot.age.dist.five <- function(df, slim.quant, view.plot, save.plot, filename){
  #split by a year
  df.year <- dlply(df,.(year, state))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year and type
  df.age <- dlply(df,.(state, year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data.type, slim.quant=slim.quant))
  
  dat.age$state <- factor(dat.age$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection", "Fifth-Infection"))
  mean.df$state <- factor(mean.df$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection", "Fifth-Infection"))
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  #colz = c("secondary"="black", "tertiary"="royalblue3")
  
  p1 <- ggplot(dat.age) + facet_grid(~state) +
    geom_jitter(aes(x=year, y=age, color=state), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    geom_violin(aes(x=year,y=age, group=year), color="gray60",  draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") + #scale_color_manual(values=colz) +
    geom_hline(aes(yintercept=1), color="red") #+ coord_cartesian(ylim=c(0,80))#,xlim=c(2015,2020))
  
  
  if(save.plot==TRUE){
    ggsave(file = filename,
           plot= p1,
           units="mm",  
           width=100, 
           height=55, 
           scale=3, 
           dpi=300)
    
  }
  if(view.plot==TRUE){
    print(p1)
  }
  
  
  
  return(mean.df)
}

plot.age.dist.five(df=out.age,slim.quant = 0.05, view.plot = T, save.plot=F, filename = NA) 
#this is when waning sends all cases back to monotypic immunity so they can re-experience secondary exposures
#this is the closest mechanistically to our ferguson model
#problem is we are allowing for re-exposure with other cases too...
#and we should get elevated mortality here
plot.age.dist(df=out.age.wane.sec.all,slim.quant = 0.05, view.plot = T, save.plot=F, filename = NA) 
# pronounced spike in age in secondary cases, then longer term elevated mean age of secondary cases

#this takes waning immunity back just one level for the target antigen, but you can't re-experience primary infection
# we see elevated age in the shorterm in secondary cases, but it declines downstream (because those primary infections are not happening in the older)
plot.age.dist(df=out.age.wane.sec,slim.quant = 0.05, view.plot = T, save.plot=F, filename = NA) 

#this is similar to above: you get spikes in age in the secondary cases in the one year,
#but then decreased age downstream. you also see reductions after the introduction in the 
#average age of primary infection because people are going through the cycle a second time
#maybe this is right-- the spike tends to happen in the year of introduction, then decrease 
#after in the data...
plot.age.dist(df=out.age.wane,slim.quant = 0.05, view.plot = T, save.plot=F, filename = NA) 


#generally increasing average age of infection across thw whole time series
plot.age.dist(df=out.age.maintain,slim.quant = 0.05, view.plot = T, save.plot=F, filename = NA) 

#and what would we expect for the age and time trend in pathology if it is linked to secondary cases?
# you;ll need to be able to differentiate between a re-exposure and a true secondary case 


#what do pathogenic cases look 


