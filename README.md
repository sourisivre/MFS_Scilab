# MFS_Scilab
### Method of Fundamental Solutions implemented in SciLAB

In this repository I make avaiable the basic code needed to apply the Method of Fundamental Solutions (MFS) to the Laplacian Operator in a 2D space.
Please check the following **[publication](https://www.scielo.br/j/rbrh/a/QCbHXbtzBDgBY5TzN9XTFnP/?lang=en)** for further details about method and its applications.

## Files

Files are organized in this way.

* `MFS_Core.sce` - contains the main function to perform the MFS numerical computing
* `Example.sce` - contains an example of the MFS application
* `Database.xls` - contains data to perform `Example.sce`

## Using in SciLAB

Files are scripts. Please load them in script editor before running.

1. After downloading all files above, you will need to replace the path to access the `Database.xls` in your system where is indicated in code
2. Now you can run `Example.sce`

### Important observations
* `MFS_Core.sce` is not able to perform alone a MFS test. You will need to create a bunch of matrixes to feed main functions.
* `Example.sce` provides enough details for you to replicate and adapt the code to your projects.
