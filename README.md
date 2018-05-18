# Demo decrypt md5 digest using MPI library run on Rocks Cluster.

## How to run?

First compile file:

```
mpicc md5.c -o md5
```

And then run with numbers of rank (process):

```
mpirun md5 -np 4 md5
```
