# Using the TUSCAN predictor

**Building a model**  
To build classifier, run
```
python ModelBuilder.py -c ClassificationMatrix.txt ClassificationFeatures.txt
```

To build regressor, run
```
python ModelBuilder.py -r RegressionMatrix.txt RegressionFeatures.txt
```

This builds models and stores them in joblib binaries

**Predicting using a model**  
Run ModelPredictor with FeatureMatrix.py, your joblib files and the relevant list of features in the same directory

TUSCAN requires the following Python packages to work:

- sys
- cPickle
- sklearn (v0.18)
- numpy
- pybedtools
- argparse
- collections


Usage:
```
python ModelPredictor.py -i Infile (.fa/.bed/.txt) -o Outfile.txt -m Regression/Classification -t InfileType (fa/bed/txt) -g genome (if supplying bed file)
```

Use -h for Help
