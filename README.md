Basic code to explore jet substructure at the MCEG level. Workflow is divided into two parts: simulation event tree generation and analysis on resulting trees.

## Step 1 - Generate events

The dir ```genSimuTrees/``` contains the code to run simulation and store output in a flat tree structure. Currently, ```genPythiaTrees.cxx``` runs the Pythia8 MCEG according to the settings stored in steering files (in ```steerFiles/```) and saves the relevant output in a root tree.

```
cd genSimuTrees
make 
./runPyTree.exe steerFiles/dis_18x275 1000 outfile.root
```

## Step 2 - Analyze generated events

The dir ```analysis/``` contains the code to read in a previously generated simulation tree and run analyses on the events within. The program ```analyzeSimuJetSubstructure.cxx``` clusters jets from the particles saved in the simulation tree and runs various substructure analyses. The file ```histogramUtilities.h``` defines and handles filling many of the relevant histograms.

Finally, the ```treeUtilities.h``` file in the top-level directory defines the common tree structure and provides helper functions that both the tree generation and analysis codes use to read/write the trees.

## Requirements

- ROOT
- [PYTHIA 8](https://pythia.org)
- [FastJet](http://fastjet.fr/quickstart.html)
- [FastJet Contrib](https://fastjet.hepforge.org/contrib/)
