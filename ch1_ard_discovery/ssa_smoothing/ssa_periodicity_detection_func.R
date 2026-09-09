
' Function to detect periodicity by monitoring the number of gradient switches. '

# gradient switches >=10 times in the series. Reversal has to cover >0.5 of the preceding change 
# (e.g. previous segment was +0.1, this segment has to be < -0.025 to count as a gradient switch)

gradient_switch_counter <- function(test_vec){
  storage_vec <- rep(0, 3) # (1st extreme, 2nd extreme, 3rd extreme), captures the bounds of two consecutive segments
  gradient_id_12 <- 0 # +1 if gradient is currently (+)ve; -1 if currently (-)ve
  gradient_id_23 <- 0 
  gradient_switch_count <- 0
  
  for (i in 1:length(test_vec)){
    if (storage_vec[1] == 0){
      storage_vec[1] <- test_vec[i] # store first value in the series vector
    }
    else { # if we are beyond the 1st value in the series vector
      if (storage_vec[2] == 0){ # if we are at the second value in the series vector
        storage_vec[2] <- test_vec[i] # store second value in the series vector
        if (storage_vec[2] > storage_vec[1]) gradient_id_12 <- 1 # initialise gradient_id_12
        else if (storage_vec[2] < storage_vec[1]) gradient_id_12 <- -1
        # gradient_id_12 maintained at 0 if the two values are the same
      }
      else { # if we are beyond the 1st two values in the series vector
        # no third value in storage_vec => first time this gradient switch is being reported
        # only if there's no third value in storage_vec do you compare the new value to storage_vec[2]
        if (storage_vec[3] == 0){
          if (test_vec[i] > storage_vec[2]){ # if new->2 gradient is positive
            if (gradient_id_12 == 1){ # 1st segment continues if 1st segment's gradient is positive
              storage_vec[2] <- test_vec[i]
            }
            else if (gradient_id_12 == -1 | gradient_id_12 == 0){ # 1st gradient switch if 1st segment's gradient is negative
              storage_vec[3] <- test_vec[i]
              gradient_id_23 <- 1
            }
          }
          else if (test_vec[i] < storage_vec[2]){ # if new->2 gradient is negative
            if (gradient_id_12 == -1 ){ # 1st segment continues if 1st segment's gradient is negative
              storage_vec[2] <- test_vec[i]
            }
            else if (gradient_id_12 == 1 | gradient_id_12 == 0){ # 1st gradient switch if 1st segment's gradient is positive
              storage_vec[3] <- test_vec[i]
              gradient_id_23 <- -1
            }
          }
        }
        # if there's a third value in storage_vec (i.e. 2nd segment already recorded) -> do the switch logic
        else{
          if (test_vec[i] > storage_vec[3]){ # if new->3 is positive
            if (gradient_id_23 == 1){ # 2nd segment continues if 2nd segment's gradient is positive
              storage_vec[3] <- test_vec[i]
            }
            else if (gradient_id_23 == -1){ # if 2nd segment's gradient is negative -> the 2nd gradient switch: the end of the 2nd segment
              if (abs(storage_vec[3] - storage_vec[2]) > 0.5 * abs(storage_vec[2] - storage_vec[1])){ # confirm magnitude of 1st switch
                gradient_switch_count <- gradient_switch_count + 1 # count the 1st switch
                storage_vec[1] <- storage_vec[2] # translate recorded points s.t. the 2nd segment is now the 1st segment
                storage_vec[2] <- storage_vec[3]
                storage_vec[3] <- test_vec[i]
                gradient_id_12 <- gradient_id_23
                if (storage_vec[3] > storage_vec[2]) gradient_id_23 <- 1 # re-initialise gradient_id_23
                else if (storage_vec[3] < storage_vec[2]) gradient_id_23 <- -1
              }
            }
          }
          else if (test_vec[i] < storage_vec[3]){ # if new->3 is negative
            if (gradient_id_23 == -1){ # 2nd segment continues if 2nd segment's gradient is negative
              storage_vec[3] <- test_vec[i]
            }
            else if (gradient_id_23 == 1){ # if 2nd segment's gradient is positive -> the 2nd gradient switch: the end of the 2nd segment
              if (abs(storage_vec[3] - storage_vec[2]) > 0.5 * abs(storage_vec[2] - storage_vec[1])){ # confirm magnitude of 1st switch
                gradient_switch_count <- gradient_switch_count + 1 # count the 1st switch
                storage_vec[1] <- storage_vec[2] # translate recorded points s.t. the 2nd segment is now the 1st segment
                storage_vec[2] <- storage_vec[3]
                storage_vec[3] <- test_vec[i]
                gradient_id_12 <- gradient_id_23
                if (storage_vec[3] > storage_vec[2]) gradient_id_23 <- 1 # re-initialise gradient_id_23
                else if (storage_vec[3] < storage_vec[2]) gradient_id_23 <- -1
              }
            }
          }
        }
      }
    }
  }
  gradient_switch_count
}