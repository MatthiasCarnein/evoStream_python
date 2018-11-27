import evoStream

evo = evoStream.EvoStream(0.05, 0.001, 100, 4, .8, .001, 100, 2*4, 250) ## init
evo.cluster([10.0, 20.0, 30.0]) ## read observation
evo.get_microweights()
evo.get_microclusters()
evo.get_macroclusters()
evo.get_macroweights()
evo.recluster(100) ## evaluate 100 more macro solutions
evo.microToMacro()



## Full Example: Read CSV file (here: comma-separated, numeric values)
import csv

evo = evoStream.EvoStream(0.05, 0.001, 100, 4, .8, .001, 100, 2*4, 250);

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

print("\nMacro Clusters (here: performs an additional 250 reclustering steps, see parameter)")
x = evo.get_macroclusters()
print(x)

print("\nMacro Weights")
x = evo.get_macroweights()
print(x)

print("\nAssignment of Micro Clusters to Macro Clusters")
x = evo.microToMacro()
print(x)
