# SBF (Sparse Binary Factorization)
SBF is a solution to the community discovery problem from undirected graphs that uses Metropolis-Hastings sampling (you don't need to know anything about that in order to use this package). 
The discovered communities can be overlapping or disjoint. 
The implementation runs faster than and scales to bigger graphs than any other community discovery packages we have tried - for example, it has been used to discover 500,000 communities from a graph with 100 Million nodes and 5 Billion edges in less than 2 hours using 16 threads. 

# Quick Start 

After cloning the repo, build as follows (for the impatient, use option `-DskipTests` to skip tests):

```
$ mvn package 
```
If you don't have maven on your computer, you'll need to first install that.


You can see two examples in the `examples` folder. 
For example, unweighted8node.config has the config needed to run on the simple 8-node (unweighted) graph specified in unweighted8node.txt. 
You can run it as follows, assuming `mvn package` succeeded.

```
$ cd examples
$ java -jar ../target/sbf-1.0.0.jar unweighted8node.config
```

Once the application finishes, you can see the output in unweighted8node.assignments.

Similarly, you can also run it using `weighted8node.config`.

`examples/runExamples.sh` does all of this. 

