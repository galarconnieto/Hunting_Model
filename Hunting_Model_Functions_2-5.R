




status <- function(x){
  p <- plot_ly(z = x, colors = colorRamp(c("red", "blue", "dark green")), type = "heatmap")
  print(p)
}

## Disperse funtion 

disperse <- function(Curr.abundance, Carr.capacity) {
  
  ### Check for cells above K
  above <- which(Curr.abundance > Carr.capacity, arr.ind = TRUE)
  
  if (nrow(above) > 0) {
    above <- above[sample(1:nrow(above)),,drop=FALSE]
    
    for (i in 1:nrow(above)) {
      cur <- above[i,]
      neighbours <- expand.grid(c(cur[1]-1,cur[1],cur[1]+1),c(cur[2]-1,cur[2],cur[2]+1))
      neighbours <- neighbours[which(neighbours[,1]>0 & neighbours[,1]<=nrow(Curr.abundance) & neighbours[,2]>0 & neighbours[,2]<=ncol(Curr.abundance)),]
      #neighbours <- neighbours[which(neighbours[1,]>0 & neighbours[1,]<=ncol(Curr.abundance) & neighbours[2,]>0 & neighbours[2,]<=ncol(Curr.abundance)),]
      #Curr.ab.local <- Curr.abundance[cbind(neighbours[,1],neighbours[,2])]
      Curr.car.local <- Carr.capacity[cbind(neighbours[,1],neighbours[,2])]
      real.neighbours <- cbind(neighbours,Curr.car.local)
      real.neighbours <- real.neighbours[which(real.neighbours[,3]/real.neighbours[,3]== 1),]
      Curr.ab.local <- Curr.abundance[cbind(real.neighbours[,1],real.neighbours[,2])]
      Curr.car.neig <- Carr.capacity[cbind(real.neighbours[,1],real.neighbours[,2])]
      probs <- 1-1/(1+exp(Curr.car.neig-Curr.ab.local))
      
      ##print(sum(probs))
      ##print(((1*(as.integer(length(probs))))/2))
      
      if (sum(probs, na.rm = T)>0){
        new.cells <- sample(1:length(probs), Curr.abundance[cur[1],cur[2]],replace=TRUE,prob=probs)
        upd <- table(factor(new.cells,levels=c(1:length(probs))))
        upd.cur <- upd[which(as.numeric(names(upd))==which(real.neighbours[,1]==cur[1] & real.neighbours[,2]==cur[2]))]
        upd.others <- upd[which(as.numeric(names(upd))!=which(real.neighbours[,1]==cur[1] & real.neighbours[,2]==cur[2]))]
        Curr.abundance[cbind(real.neighbours[as.numeric(names(upd.cur)),1],real.neighbours[as.numeric(names(upd.cur)),2])] <- 
          as.numeric(upd.cur)
        Curr.abundance[cbind(real.neighbours[as.numeric(names(upd.others)),1],real.neighbours[as.numeric(names(upd.others)),2])] <- 
          Curr.abundance[cbind(real.neighbours[as.numeric(names(upd.others)),1],real.neighbours[as.numeric(names(upd.others)),2])]+as.numeric(upd.others)
        ##print(probs)
        ##print(Curr.car.local)
        ##print(Curr.ab.local)
      }
      else {
        
        #Curr.ab.local[cur] <- (Curr.ab.local[cur]-(Curr.ab.local[cur]*0.01))
        #print(c("nowhere to go! at: ", as.character(cur)))
      }
    }
  }
  return(Curr.abundance)
  
}


##Population dynamics

popdynamics <- function(Curr.abundance, Carr.capacity, R, D){
  # new cell size as funciotn of current grid and dispersal distance in meters
  #print(c(sum(Curr.abundance),"sex!!"))
  num.cell <- D %/% 250
  curr.grid <- Curr.abundance
  num.row <- (nrow(curr.grid) %/% num.cell)
  if ((nrow(Curr.abundance)/num.cell) > num.row){
    num.row <- num.row+1
  }
  num.col <- ncol(curr.grid) %/% num.cell
  if ((ncol(Curr.abundance)/num.cell) > num.col){
    num.col <- num.col+1
  }
  new.pop.matrix <- matrix(1:(num.col*num.row), nrow = num.row, ncol = num.col, byrow = T)
  
  # new area
  for (i in (1:ncol(new.pop.matrix))){
    a <- 1 + (num.cell*(i-1))
    b <- num.cell*i
    if(b>ncol(Curr.abundance)){
      b=ncol(Curr.abundance)
    }
    for (j in (1:nrow(new.pop.matrix))){
      c <- 1 + (num.cell*(j-1))
      d <- num.cell*j
      if(d>nrow(Curr.abundance)){
        d=nrow(Curr.abundance)
      }
      ##print(c(a,b,c,d))
      # select area
      new.abu <- Curr.abundance[c:d, a:b]
      #big.abu <- as.matrix(expand.grid(c(c:d) , c(a:b)))
      local.abundance <- sum(Curr.abundance[c:d , a:b], na.rm = TRUE)
      ##print(local.abundance)
      local.carr.cap <- sum(Carr.capacity[c:d , a:b], na.rm = TRUE)
      ##print(c("K=", local.carr.cap))
      if(local.carr.cap==0){
        break
      }
      new.pop <- round(local.abundance*R*(1-(local.abundance/local.carr.cap)))
      if(new.pop==0){
        break
      }
      if((local.abundance+new.pop)<0){
        break #new.pop <- local.abundance
      }
      ##print(c("np", new.pop))
      #where <- which(x[c:d , a:b] > 0, arr.ind = TRUE)
      where <- which(Curr.abundance[c:d, a:b] > 0, arr.ind = TRUE)
      if(length(where)== 0){
        break #new.pop <- local.abundance
      }
      #print(where)
      probs <- new.abu[cbind(where[,1], where[,2])]
      probs <- 1-1/(1+exp(probs))
      ##print(probs)
      ##print(c("i is",i))
      ##print(j)
      ##print(new.pop)
      
      if ((abs(new.pop)>0)){ #((sum(probs, na.rm = TRUE)>0) & 
        ##print("son mas que 0")
        
        change <- sample(1:length(where[,2]), size =abs(new.pop), replace = T, prob = probs)
        #change <- sample(1:length(where[,2]), size = 3, replace = T) #prob = probs)
        #print(c("ch", change))
        upd <- matrix(where[change,], ncol = 2)
        #print("upd")
        #print(upd)
        #if(abs(new.pop)==0){
        # print("es cero")
        #}else{
        ddd <- c(1:abs(new.pop))
        #print(ddd)
        for (zz in seq_along(ddd)){
          #print(zz)
          cur <- upd[zz,]
          #old <- new.abu[cur[1], cur[2]]
          #print(c("old was", old)) 
          #print(" cur")
          #print(cur)
          if(new.pop>0){
            new.abu[cur[1], cur[2]] <- new.abu[cur[1], cur[2]] + 1
          }else{
            if((new.abu[cur[1], cur[2]])==0){
              print("no animals left")
              break #new.abu[cur[1], cur[2]] <- 0
            }else{new.abu[cur[1], cur[2]] <- new.abu[cur[1], cur[2]] - 1}
          }
          
          #if(new>0){
          # break
          #print("aqui es")
          #}
          #print(c("new is", new))
        }
        #}
      }
      else{
        #print("negative probs")
      }
      #print(x[c:d, a:b])
      #print(new.abu)
      Curr.abundance[c:d, a:b] <- new.abu
    }
    #print(x)
    
  }
  #print(c(sum(Curr.abundance), "new abundance"))
  return(Curr.abundance)
}


## Movement of animals per month
movement.month <- function(Curr.abundance.All, Carr.capacity.small, Carr.capacity.medium,
                           Carr.capacity.large, Carr.capacity.monkeys){
  #print('disperse!')
  Curr.abundance.small <- Curr.abundance.All$small
  for (mov.small in 1:3){
  Curr.abundance.small <- disperse(Curr.abundance.small, Carr.capacity.small)
  }
  
  Curr.abundance.monkeys <- Curr.abundance.All$monkeys
  for (mov.monkeys in 1:3){
    Curr.abundance.monkeys <- disperse(Curr.abundance.monkeys, Carr.capacity.monkeys)
  }
  Curr.abundance.medium <- Curr.abundance.All$medium
  for (mov.medium in 8){
    Curr.abundance.medium <- disperse(Curr.abundance.medium, Carr.capacity.medium)
  }
  Curr.abundance.large <- Curr.abundance.All$large
  for (mov.large in 1:8){
    Curr.abundance.large <- disperse(Curr.abundance.large, Carr.capacity.large)
  }
  curr.abund.after.movement <- list("small" = Curr.abundance.small, "medium"= Curr.abundance.medium,
                                    "large"= Curr.abundance.large, "monkeys"= Curr.abundance.monkeys)
  return(curr.abund.after.movement )
}



# calculate encounter rate of a population, using the local abundance per km-sq 
enc.rate <- function(x){(-1.052e-04)*x^2 + (1.322e-02)*x + (7.920e-03)}

walk <- function(start.cell, distance){
  tot.cells <- (distance/250) # number of cells to complete hunting trip
  xx <- distances(start.cell)
  cinco.mil <- which(xx==distance, arr.ind=T) ## Which cells are at max distance
  cinco.mil <- cinco.mil[sample(1:nrow(cinco.mil)),,drop= FALSE] # Randomize cell numbers
  #print(cinco.mil)
  a.donde <- sample(1:length(cinco.mil[,2]), 1) # Selection of final cell 
  #print(a.donde)
  end.cell <- cinco.mil[a.donde,] # defines cell coordinates of final cell
  #print(c("end cell =", end.cell))
  #end.cell <- walk(start.cell, 5000)
  enx <- end.cell[2] - start.cell[2] ## calculates the movement (in cells) on x axix
  
  #print(enx)
  eny <- end.cell[1] - start.cell[1] ## calculates the movement (in cells) on y axix
  #print(eny)
  movement <- matrix(0, ncol = 2, nrow = tot.cells) ## This part creates a path for the trip, generating movements (steps) 
  if(enx!=0){
    movement[,2][(tot.cells+1-(abs(enx))):tot.cells] <- enx/abs(enx)
  }else{
    #movement[,2][(tot.cells+1-(abs(enx))):tot.cells] <- 0
  }
  
  if(eny!=0){
    movement[,1][1:abs(eny)] <- eny/abs(eny)
  }else{
    #movement[,1][1:abs(eny)] <-0
  }
  
  
  rrr <- movement[sample(1:nrow(movement)),,drop=FALSE]  ## randomization of the steps to get to the end.cell
  #print(rrr)
  go.to <- matrix(NA, ncol = 2, nrow = tot.cells) # creates a matrix with the cell numbers for the path
  xxx <- start.cell[2]
  yyy <- start.cell[1]
  
  for (move in 1:tot.cells){
    
    xxx <- xxx + rrr[move,2]
    yyy <- yyy + rrr[move,1]
    go.to[move,] <- c(yyy,xxx)
  }
  return(go.to)
}


distances <- function(start.cell){
  dist.to.cell <- matrix(NA, nrow(template.mtx), ncol(template.mtx),  byrow=T)
  #start.cell <- c(1,1)
  x1<- start.cell[2]
  #print(x1)
  y1<- start.cell[1]
  #print(y1)
  
  for (dd in 1:length(cell.nums)){
    end.cell <- c(which(cell.nums==dd, arr.ind = TRUE))
    x2 <- end.cell[2]
    #print(x2)
    y2 <- end.cell[1]
    #print(y2)
    #print
    dist <- round(sqrt((abs(x1 - x2))^2+(abs(y1 - y2))^2)) * 250
    #print(dist)
    dist.to.cell[y2,x2] <- dist
  }
  return(dist.to.cell)
}



### Simulating a single hunting event. Hunter walks a "distance" in meters,
###searching for large and medium on the outbond walk and small and monkeys in the way back.
### This function returns a hunting event report including time (current time period), location(y),
###ocation(x), animal hunted (size category), and distance walked
hunt.event <- function(start.cell, distance){
  event <- cbind(month, NA, NA, "none", 0 , (distance*2), people, budget)
  one.walk <- walk(start.cell, distance)
  return.walk <- matrix(NA, ncol = 2, nrow = length(one.walk[,1]))
  return.walk[,1] <- rev(one.walk[,1])
  return.walk[,2] <- rev(one.walk[,2])
  for (eee in 1:length(one.walk[,2])){
    #print(eee)
    cell <- one.walk[eee,]
    ########Calculates encounter rates for each size group. Since hunter are capable of tracking animals,
    #for large and medium sized animals, the ER is 1 if animals are present,
    if(as.integer(Curr.abundance.All$large[cell[1],cell[2]]) > 0){
      event <- cbind(month, cell[1], cell[2], "large", 120, (eee*250), people, budget)
      break}
    if(as.integer(Curr.abundance.All$medium[cell[1],cell[2]]) > 0){
      event <- cbind(month, cell[1], cell[2], "medium", 26.25, (eee*250), people, budget)
      break}
  }

### Return walk. SEarching for small and monkeys  
  if (event[5]==0){ 
    for (rrr in 1:length(return.walk[,2])){
      #print(rrr+20)
      cell <- return.walk[rrr,]
      #print(cell)
      if(Curr.abundance.All$small[cell[1],cell[2]]>0){
        small <- enc.rate(as.integer(Curr.abundance.All$small[cell[1],cell[2]])*16)  
      } else{small=0}
      
      #print("small")
      #print(small)
      
      if (small>0){
        
        #(sum(Curr.abundance.All$small))
        encounter <- sample(c(0,1), 1, prob = c((1-small), small))
        #print(c("encs", encounter))
        if(encounter==1){
          event <- cbind(month, cell[1], cell[2], "small", 3.5, (distance + (rrr*250)), people, budget)
          break
        }
      }
      
      if(Curr.abundance.All$monkeys[cell[1],cell[2]]>0){
        monkeys <- enc.rate(as.integer(Curr.abundance.All$monkeys[cell[1],cell[2]])*16)  
      } else{monkeys=0}
      
      #print("monkeys")
      #print(monkeys)
      
      if (monkeys>0){
        #(sum(Curr.abundance.All$monkeys))
        encounter <- sample(c(0,1), 1, prob = c((1-monkeys), monkeys))
        #print(c("encm", encounter))
        if(encounter==1){
          event <- cbind(month, cell[1], cell[2], "monkeys", 4, (distance + (rrr*250)), people, budget)
          break
        }
      }
    }
  }
  
  
  #if (rrr == length(return.walk[,2]){
  #print(event)
  return(event)

}


hunt.event.close <- function(start.cell, distance){
  event <- cbind(month, NA, NA, "none", 0 , (distance*2), people, budget)
  one.walk <- walk(start.cell, distance)
  return.walk <- matrix(NA, ncol = 2, nrow = length(one.walk[,1]))
  return.walk[,1] <- rev(one.walk[,1])
  return.walk[,2] <- rev(one.walk[,2])
  for (eee in 1:length(one.walk[,2])){
    #print(eee)
    cell <- one.walk[eee,]
    ########Calculates encounter rates for each size group. Since hunter are capable of tracking animals,
    #for large and medium sized animals, the ER is 1 if animals are present,
    if(as.integer(Curr.abundance.All$large[cell[1],cell[2]]) > 0){
      event <- cbind(month, cell[1], cell[2], "large", 120, (eee*250), people, budget)
      break}
    if(as.integer(Curr.abundance.All$medium[cell[1],cell[2]]) > 0){
      event <- cbind(month, cell[1], cell[2], "medium", 26.25, (eee*250), people, budget)
      break}
    
    if(Curr.abundance.All$small[cell[1],cell[2]]>0){
      small <- enc.rate(as.integer(Curr.abundance.All$small[cell[1],cell[2]])*16)  
    } else{small=0}
    
    #print("small")
    #print(small)
    
    if (small>0){
      
      #(sum(Curr.abundance.All$small))
      encounter <- sample(c(0,1), 1, prob = c((1-small), small))
      #print(c("encs", encounter))
      if(encounter==1){
        event <- cbind(month, cell[1], cell[2], "small", 3.5, (distance + (eee*250)), people, budget)
        break
      }
    }
    
    if(Curr.abundance.All$monkeys[cell[1],cell[2]]>0){
      monkeys <- enc.rate(as.integer(Curr.abundance.All$monkeys[cell[1],cell[2]])*16)  
    } else{monkeys=0}
    
    #print("monkeys")
    #print(monkeys)
    
    if (monkeys>0){
      #(sum(Curr.abundance.All$monkeys))
      encounter <- sample(c(0,1), 1, prob = c((1-monkeys), monkeys))
      #print(c("encm", encounter))
      if(encounter==1){
        event <- cbind(month, cell[1], cell[2], "monkeys", 4, (distance + (eee*250)), people, budget)
        break
      }
    }
    
  }
  
  
  ### Return walk. SEarching for small and monkeys  
  if (event[5]==0){ 
    for (rrr in 1:length(return.walk[,2])){
      #print(rrr+20)
      cell <- return.walk[rrr,]
      #print(cell)
      (Curr.abundance.All$small[cell[1],cell[2]])
      if(Curr.abundance.All$small[cell[1],cell[2]]>0){
        small <- enc.rate(as.integer(Curr.abundance.All$small[cell[1],cell[2]])*16)  
      } else{small=0}
      
      #print("small")
      #print(small)
      
      if (small>0){
        
        #(sum(Curr.abundance.All$small))
        encounter <- sample(c(0,1), 1, prob = c((1-small), small))
        #print(c("encs", encounter))
        if(encounter==1){
          event <- cbind(month, cell[1], cell[2], "small", 3.5, (distance + (rrr*250)), people, budget)
          break
        }
      }
      
      if(Curr.abundance.All$monkeys[cell[1],cell[2]]>0){
        monkeys <- enc.rate(as.integer(Curr.abundance.All$monkeys[cell[1],cell[2]])*16)  
      } else{monkeys=0}
      
      #print("monkeys")
      #print(monkeys)
      
      if (monkeys>0){
        #(sum(Curr.abundance.All$monkeys))
        encounter <- sample(c(0,1), 1, prob = c((1-monkeys), monkeys))
        #print(c("encm", encounter))
        if(encounter==1){
          event <- cbind(month, cell[1], cell[2], "monkeys", 4, (distance + (rrr*250)), people, budget)
          break
        }
      }
    }
  }
  
  
  #if (rrr == length(return.walk[,2]){
  #print(event)
  return(event)
  
}


#Calculate near cells
distance.from.village <- distances(c(139,122))
near <- which(distance.from.village <= 1500, arr.ind = T)

far.away <- which(distance.from.village == 1500, arr.ind = T)


## HUNTING

# calculate encounter rate of a population, using the local abundance per km-sq 
enc.rate <- function(x){(-1.052e-04)*x^2 + (1.322e-02)*x + (7.920e-03)}




#Create matrix with cell numbers
#cell.nums <- matrix(1:length(template.mtx), nrow(template.mtx), ncol(template.mtx),  byrow=T)



