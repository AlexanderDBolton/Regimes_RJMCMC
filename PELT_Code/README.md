### The PELT Algorithm

The Pruned Exact Linear Time (PELT) algorithm was created by R. Killick, P. Fearnhead and I.A. Eckley. Their paper on the algorithm can be found at https://arxiv.org/pdf/1101.1438.pdf. Given data and a model to fit to the data (e.g. Bernoulli process, Markov chain, Poisson process, Gaussian) PELT finds the set of change points that minimises a given cost function. In this implementation the available cost functions are AIC and BIC. In the worst case the algorithm takes O(n^2) time (where n is the number of data points) but the algorithm will take O(n) time if the change points occur linearly with n.

### Input parameters

```
yourPrefix_everything.txt
yourCostFunction
dataLength
diff
```

