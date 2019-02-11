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
Run ModelPredictor with your joblib files in the same directory

TUSCAN requires the following Python packages to work:
Native:
- sys
- argparse
- collections
- os
- re
- string
- multiprocessing

Additional:
- sklearn
- numpy
- pybedtools 

Usage:
```
python TUSCAN.py -o Outfile.txt -m Regression/Classification -g Infile.fa [-s startPosition] [-f endPosition] [-c Chromosome] [-e (Include if extracting region from genome)] [-t Number of threads]
```

Please see the head of TUSCAN.py for more details

Use -h for Help
