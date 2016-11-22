# FishSimulation
This is a simple cellular automata/random walk simulation written in matlab which simulates a fish population with an 
external food source (plankton) where fish have a chance to migrate to  an adjacent cell with the greatest increae in 
plankton. 

Rules
Each cell or school can support no more than "MaxfishPerCell" fish.
If there is no plankton present in a cell, one fish will die in that cell.
For every plankton present in the cell, the number of fish in the cell will increae by 1, but not above "MaxFishPerCell".
Fish cannot eat plankton and reproduce if there is less than 2 fish in the cell.
Fish have a chance to migrate to the adjacent cell with the largest increase in planton.
If there is greater than "MaxFishPerCell" fish in a cell, all fish die of overpopulation.
If there are no live fish in the cell but between 2 to "MaxFishPerCell" * 4 fish in adjacent cells, "FishBornInNewCell" new fish are born into that cell.
