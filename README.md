To install dependencies:

1) Install Anaconda or miniconda (preferably).

2) Go to folder containing the adcircpy.yml file and run the command "conda env create -f adcpy.yml"

3) To activate the environment run "source activate adcpy"

4) To use the package on your project, symlink the source code folder to your project's directory, and it will allow you to import.  Alternatively use:
```python
    import sys
    sys.path.append("/path/to/adcpy")
```
For questions: <br />
jreniel@gmail.com


