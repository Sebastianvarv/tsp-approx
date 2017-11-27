# tsp-approx
Quite efficient and good approximation to TSP

All the code is in the tsp.py file. I did also include test data.

I used multistart local search, which has proven itself very well. It creates iteration number of times a random
solution and begins to optimize it until it finds local minimum. To get actually good result we have to run it enough
times to randomly hit the global mininimum or something close to that. For optimization it uses variation
of 2-opt where edges are swapped and then it checks if that did help to improve the current solution, if it does then
new solution is saved as the best. I achieved really good performance since I created distances dictionary (matrix) at
the beginning and used it afterwards, it is much faster than computing the same distances over and over again. Also if
I swapped some edges then I computed distance locally using the previously created distances dictionary. Also I did
make list of all the closest nodes, which has elements sorted in distance order which means if I want to find closest
node fot the 5-th node then I just fetch 5-th element from that list and get first (the closest) from there.

There are data sets with different sizes in data directory, each data set consists of 2 files. Usual text files contain
coordinates of the points and .swog files are for visualizing them. SWOG files can be used at
http://biit.cs.ut.ee/SWOG/index.cgi to visualize your data.

Current best result on 1000 elements data set is 15.1 seconds on 2 core i5 laptop with distance 24516. Feel free to
clone this repo and try to improve the result.