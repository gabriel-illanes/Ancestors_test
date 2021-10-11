library(data.table)
library(magrittr)
#setwd()

lindeley=function(x, y = finalesChr[1:21]+1){
  lin=c(0,x)
  L=length(lin)
  lenghts=vector(mode='numeric',length=(L-1))
  i=1
  while(i < L){
    while(x[i]==1 & i < L){
      if(!i %in% y){
        lin[i+1]=lin[i]+x[i]
        i=i+1  
      }
      else{
        lin[i] = 1
        lengths[i-1] = lin[i-1]
        i = i+1
      }
    }
    lengths[i-1]=lin[i]
    i=i+1
  }  
  lin=lin[-1]
  return(list('chunks'=lin,'lenghts'=lenghts))
}

findChunk=function(lgChunks, valMax=max(lgChunks)){
  end<-which(lgChunks==valMax)
  beg<-max(which(lgChunks[1:end]==0))+1
  return(seq(beg,end))
}

ped=fread('data/newer_nat_corrected.txt', data.table=FALSE, sep=' ')
pedDT=fread('data/newer_nat_corrected.txt', sep=' ') #en formato data. table
fam=fread('data/eur_nat_afro_uru.fam', data.table = FALSE, sep=' ')
map=read.table('data/newer_native_posta_corrected.map')

pops<- read.table('data/unmasked.ind')[,c(1,3)] %>% merge(fam[,1:2],.,by='V1', sort=FALSE) 
hapsFam=rep(fam[,1], each=2)
hapsPopSinInfo=c(as.character(rep(pops[,3], each=2)),rep('UrA',20),rep('UrN',20))
hapsPopSinInfo2=rep(c(as.character(pops[,3]),fam[526:545,1]), each=2)

centromeres=read.csv('data/posCentromeros.csv',sep=';')
finChr=rep(0,22)
centers=rep(0,22)
for (i in 1:22){finChr[i+1]=finChr[i]+which.max(map[map$V1==i,4])
centers[i]=finChr[i]+which.min(abs(map[map$V1==i,4] - centromeres[i,2]*1e6))}
finChr=finChr[-1]

newPed=ped

matInfoAux=!is.na(ped)+0
PropInfo=apply(matInfoAux,1,sum)/ncol(matInfoAux)
zeroInfo=which(PropInfo<0.05)
tabla=table(hapsPopSinInfo[zeroInfo])
tabla
matInfo=matInfoAux[-zeroInfo,]
hapsPop=hapsPopSinInfo[-zeroInfo]
hapsPop2=hapsPopSinInfo2[-zeroInfo]
summary(PropInfo[-zeroInfo])
pedInfo=ped[-zeroInfo,]

namesPop = names(table(hapsPop))
pedInfo = as.matrix(pedInfo)
pedInfoD = matrix(NA, ncol = 363578, nrow = 441)
for(i in 1:length(table(hapsPop))){
  m = pedInfo[hapsPop == namesPop[i], ]
  m[is.na(m)] = mean(m, na.rm = TRUE)
  pedInfoD[hapsPop == namesPop[i], ] = m
}
pedInfoD = apply(pedInfoD, 1, cumsum)

ChunksNatives=t(apply(matInfoAux[1051:1090,],1,lindeley))
lenghtsChunksNatives=t(sapply(seq(1:40),function(x) ChunksNatives[[x]]$chunks))
lengthsChunks=t(sapply(seq(1:40),function(x) ChunksNatives[[x]]$lenghts))

max_statistic = function(lenghtsChunks, map, finChr){
  est = matrix(NA, 40, 22)
  for(i in 1:22){
    for(j in 1:40){
      cr = lenghtsChunks[j, map[,1] == i]
      l_cr = length(cr)
      fin_pos = c(1:l_cr)[!cr == 0]
      max_cm = 0;
      if(!length(fin_pos) == 0){
        pos_dif = cr[fin_pos]
        for(k in 1:length(fin_pos)){
          if(i == 1){
            max_cm = max(max_cm, map[fin_pos[k],3] - map[fin_pos[k]-pos_dif[k]+1,3])  
          } else{
            max_cm = max(max_cm, map[finChr[i-1] + fin_pos[k],3] - map[finChr[i-1]+fin_pos[k]-pos_dif[k]+1,3])
          }
          
        }
      }
      est[j,i] = max_cm/100
    }
  }
  return(est)
}

max_stat = max_statistic(lenghtsChunks, map, finChr)

sum_statistic = function(lenghtsChunks, map, finChr){
  est = matrix(NA, 40, 22)
  for(i in 1:22){
    for(j in 1:40){
      cr = lenghtsChunks[j, map[,1] == i]
      l_cr = length(cr)
      fin_pos = c(1:l_cr)[!cr == 0]
      sum_cm = 0;
      if(!length(pos_fin) == 0){
        pos_dif = cr[fin_pos]
        for(k in 1:length(fin_pos)){
          if(i == 1){
            sum_cm = sum_cm + map[fin_pos[k],3] - map[fin_pos[k]- pos_dif[k] + 1,3]
          } else{
            sum_cm = sum_cm + map[finChr[i-1] + fin_pos[k],3] - map[finChr[i-1] + fin_pos[k]- pos_dif[k] + 1,3]
          }
        }
      }
      est[j,i] = sum_cm/100
    }
  }
  return(est)
}

sum_stat = sum_statistic(lenghtsChunks, map, finChr)

write.table(max_stat, file = "max_stat.csv", row.names = FALSE, col.names = FALSE)
write.table(sum_stat, file = "sum_stat.csv", row.names = FALSE, col.names = FALSE)
write.table(map[finChr, 3], file = "lenghts_chr.csv", row.names = FALSE, col.names = FALSE)


