### Monty Hall
#by Levi Raskin
#https://github.com/Levi-Raskin

#See readme for longer explanation of problem.
#Short version: you have three doors to choose from. 
#Behind one of them is a car and behind the other two are goats.
#After you select one door, the host (Monty Hall) reveals one of the goat-doors
#Do you swap?

#Made for Bryn Mawr Biology class B215: Biostatistics with R fall 2023

# Not swapping ------------------------------------------------------------

#To test no swapping, all we have to do is create a vector containing the options
#Then, sample twice-- once with the original choice and once with the correct choice
#Finally, we compare the two and same the T/F output to the results vector
#Comparing the proportion of T to F gives us the % the original choice was correct
#Should be 1/3rd and that's what we observe

#vector containing possible choices
doors <- 1:3

#an empty vector to store our results in
results <- rep(NA, 1000)

for(i in 1:1000){
  #initial choice of which door (1, 2, 3)
  choice <- sample(doors, 1)
  #which door is correct? 1, 2, or 3
  correct <- sample(doors, 1)
  #is our choice the sample as the correct door?
  results[i] <- choice==correct
  print(i)
}

table(results)

print(paste("Proportion correct:", table(results)[2]/sum(table(results))))

# Swapping ----------------------------------------------------------------

#Swapping is much trickier to code
#We still create a vector containing the options of doors (here, 1:3)
#Then, we sample our first choice; this will be discarded later
#Again, we sample the correct choice
#Then, we have two cases. The first case is if we chose correctly and the second is if we didn't
#This code could be condensed to one for loop here, by comparing these two samples and storing as a sep. results vector
#However, I split this code to two for loops to make it easier for students to understand what I was doing
#We need the host to reveal a door and, critically, it can *neither* be the door you picked originally *nor* can it be the door contianing a car
#Thus, if we chose the correct door originally, then the host reveals one of the other doors-- e.g. if correct+choice is door 1, host can reveal either 2 or 3
#The choice to reveal one of the other doors is done randomly, but it doesn't really matter- you're gonna swap to an incorrect door regardless
#But, if we chose an incorrect door originally, then the host needs to reveal to us a door that is neither our door, nor the car door
#Thus, the wrong choice revealed to us is just the doors vector, minus the correct choice and your original choice (making use of the fact that doors 1, 2, and 3 are in elements 1, 2 and 3)
#For example, this would be equivalent to correct door 1, choice door 2, revealed door *must* be door 3
#Finally, our new choice (i.e., the door we swapped to) is the doors vector minus our original choice and the wrong door
#Again, here, the code could be more efficient, but I broke it to two steps to both show students the "which" function and help break it down'
#finally, we compare our final choice (after swapping) to the correct door and save the T/F as a vector

#door 1, 2, or 3
doors <- 1:3

# we'll be running the experiment 10000 times and saving our correct/incorrect to this vector
results <- rep(NA, 1000)

# using a for loop to repeat the choosing of door 1, 2, or 3 1k times
for(i in 1:1000){#initial random choice of which door (1, 2, 3)
  choice <- sample(doors, 1)
  #which door is correct; random and independent from the choice 
  correct <- sample(doors, 1)
  #which door host reveals (which has to have a goat)
  if(choice== correct){
  
    #can condense to a single for loop for both swapping and not swapping here:
    #results_noswap[i] <- TRUE
    
    #from the doors, remove both the choice and correct positions
    wrongdoor <- sample(doors[-choice], size = 1)
    #to condense:
    #results_noswap[i] <- FALSE
    }
  else{wrongdoor <- doors[-c(choice,correct)]}
  #new choices not including the wrong door option
  newchoices <- doors[-wrongdoor]
  #here, we swap by removing our original choice from the new choices vector
  finalchoice <- newchoices[-which(newchoices==choice)]
  results[i] <- finalchoice==correct
  print(i)
}

table(results)

print(paste("Proportion correct:", table(results)[2]/sum(table(results))))

