# MLNeams

To build it just type "make". That will create a set of executables full_p.a, RS1_p.a, etc
Typing ./full_p.a will print the usage information.

But for reference the usage is 

```console
./full_p < matrix file > < output file > < percentage as decimal > <maxrowstosample> <edge samples> <solve> <matvecs >
```
<percentage as decimal> is a decimal between zero and one that determines the percentage of the total columns to sample as interior points.

< maxrowstosample> is a integer you can set to put a cap on the number of rows to sample ( i.e sample 20% of rows up to 500 columns ). If you do
not want a cap, set this to -1.

<edge samples> this is a integer to define the number of edge samples  (same as before) . The total number of samples is  ~=  edgesamples + (#columns)*percentage

< sovle > is either 0 or 1, depending on if you want to solve the matrix as well

<matvecs> is 1 if you want to use matrix vector multiplications to extract the values. If you set it to zero, the code will use MatGetColumn(...) to extract the samples.

e.g For a pretend matrix stored at ./pretend_matrix.pbin with 1000 columns

./full.a     ./pretend_matrix.pbin     ./outputfile    0.6      -1    10   1    1     -- > this will do feature extraction with 600 interior points and 10 edge points (610 samples) , then solve the matrix

./full.a     ./pretend_matrix.pbin     ./outputfile    0.6      250    10   1    1     -- > this will do feature extraction with 250 interior points and 10 edge points (260 samples because the max was met

./full.a     ./pretend_matrix.pbin     ./outputfile    0.3      -1   5  0    1     -- > this will do feature extraction with 300 interior points and 5 edge points (305 samples) , but not solve the matrix

./full.a     ./pretend_matrix.pbin     ./outputfile    0.3      -1   5  0    0     -- > this will do feature extraction with 300 interior points and 5 edge points (305 samples) using MatGetColumn to read the sample columns from file rather than by using matrix vector multiplications.
