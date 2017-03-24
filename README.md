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

Usage:
```
python ModelPredictor.py -i Infile -o Outfile -m Regression/Classification -t InfileType (fasta/bed/txt) -g genome (if supplying bed file)
```

Use -h for Help
