SpecificDigit <- function(digit){
  max <- as.data.frame(do.call(pmax, digit))
  all <- cbind(max, digit)
  all <- all[all$`do.call(pmax, digit)` >= 3, ]
  test <- all[, c(2:ncol(all))]
  
  tau_max <- as.data.frame((rowSums(1 - (test/(do.call(pmax, test)))))/(length(test) - 1))
  names(tau_max)[1] <- "Tau"
  x_max <- merge(tau_max, test, by = "row.names", all = TRUE)
  x_max$Row.names <- as.character(x_max$Row.names)
  row.names(x_max) <- x_max$Row.names
  
  filter_max <- as.data.frame(x_max[which(x_max$Tau > 0.5), ])
  
  tissue_name <- paste("MAX_", colnames(filter_max[, 3:ncol(filter_max)]))
  rownames(filter_max) <- filter_max$Row.names
  filter_max1 <- filter_max[, -1]
  filter_max2 <- filter_max1
  for(i in seq(nrow(filter_max1))){
    a <- as.matrix(filter_max1[i, 2:ncol(filter_max1)])
    filter_max2[i, ncol(filter_max)] <- tissue_name[which(a == a[which.max(a)], arr.ind = T)[2]]
  }
  filter_max2$gene <- rownames(filter_max2)
  
  tau_min <- as.data.frame(1 - (1/((rowSums((test/(do.call(pmin, test)))))/(length(test) - 1))))
  names(tau_min)[1] <- "Tau"
  x_min <- merge(tau_min, test, by = "row.names", all = TRUE)
  x_min$Row.names <- as.character(x_min$Row.names)
  row.names(x_min) <- x_min$Row.names
  filter_min <- as.data.frame(x_min[which(x_min$Tau > 0.7), ])
  tissue_name <- paste("MIN_", colnames(filter_min[, 3:ncol(filter_min)]))
  rownames(filter_min) <- filter_min$Row.names
  filter_min1 <- filter_min[, -1]
  filter_min2 <- filter_min1
  for(i in seq(nrow(filter_min1))){
    a <- as.matrix(filter_min1[i, 2:ncol(filter_min1)])
    filter_min2[i, ncol(filter_min)] <- tissue_name[which(a == a[which.min(a)], arr.ind = T)[2]]
  }
  filter_min2$gene <- rownames(filter_min2)
  
  mergeAll <- rbind(filter_max2, filter_min2)
  names(mergeAll)[length(mergeAll) - 1] <- "Tau_info"
  as.data.frame(mergeAll)
}
