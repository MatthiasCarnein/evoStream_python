# evoStream - Evolutionary Stream Clustering Utilizing Idle Times


This is the implementation of an evolutionary stream clustering algorithm as proposed in our article in the Journal of Big Data Research.

> Carnein M. and Trautmann H. (2018), "evoStream - Evolutionary Stream Clustering Utilizing Idle Times", Big Data Research., May, 2018. Vol. 14, pp. 101 - 111. 

The online component uses a simplified version of `DBSTREAM` to generate micro-clusters.
The micro-clusters are then incrementally reclustered using an evloutionary algorithm.
Evolutionary algorithms create slight variations by combining and randomly modifying existing solutions.
By iteratively selecting better solutions, an evolutionary pressure is created which improves the clustering over time.
Since the evolutionary algorithm is incremental, it is possible to apply it between observations, e.g. in the idle time of the stream.
Whenever there is idle time, we can call the `recluster` function of the reference class to improve the macro-clusters (see example).
The evolutionary algorithm can also be applied as a traditional reclustering step, or a combination of both.
In addition, this implementation also allows to evaluate a fixed number of generations after each observation.

## Installation

This is the Python port of evoStream. It is based on the C++ implementation with wrappers for Python.

In order to install the module, run the following command in the modules main directory:

```
python setup.py install --force
```

For convenience, the command can be issued using the `install.bat` or `install.sh` files.


## Usage

Once installed, the interfaces are the same as in the C++ and R implementations:


```Python
import evoStream

evo = evoStream.EvoStream(0.05, 0.001, 100, 4, .8, .001, 100, 2*4, 1000) ## init
evo.cluster([10.0, 20.0, 30.0]) ## read observation
evo.get_microweights()
evo.get_microclusters()
evo.get_macroclusters()
evo.get_macroweights()
evo.recluster(100) ## evaluate 100 more macro solutions
evo.microToMacro()



## Full Example: Read CSV file (here: comma-separated, numeric values)
import csv

evo = evoStream.EvoStream(0.05, 0.001, 100, 4, .8, .001, 100, 2*4, 1000);

with open('data.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        evo.cluster(row)
        evo.recluster(1) # evaluate 1 generation after every observation. This can be adapted to the available time


print("Micro Clusters:")
x = evo.get_microclusters()
print(x)

print("\nMicro Weights:")
x = evo.get_microweights()
print(x)

print("\nMacro Clusters (here: performs an additional 1000 reclustering steps, see parameter)")
x = evo.get_macroclusters()
print(x)

print("\nMacro Weights")
x = evo.get_macroweights()
print(x)

print("\nAssignment of Micro Clusters to Macro Clusters")
x = evo.microToMacro()
print(x)
```


## Related Implementations

The original implementation is available as an R-package here: [https://wiwi-gitlab.uni-muenster.de/m_carn01/evoStream](https://wiwi-gitlab.uni-muenster.de/m_carn01/evoStream)

There is also a C++ port of evoStream. It is available here: [https://wiwi-gitlab.uni-muenster.de/m_carn01/evoStream_C](https://wiwi-gitlab.uni-muenster.de/m_carn01/evoStream_C)

