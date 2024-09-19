# Thesis-Project
## General Notes

In order to an implementative choice, in the code, it is used the camelCase notation for coherence with
the CPLEX library although it is not the native notation of C++. 


## How to compile

Command to compile a file in Cluster, after login in `login.dei.unipd.it` log in `arrow-16.dei.unipd.it` for compile for getting the binary:

First command :
```bash
g++ -I/nfsd/rop/sw/ibm/cos221/concert/include -I/nfsd/rop/sw/ibm/cos221/cplex/include -Wno-deprecated-declarations -std=c++17 -c Utility.cpp -o Utility.o

g++ -I/nfsd/rop/sw/ibm/cos221/concert/include -I/nfsd/rop/sw/ibm/cos221/cplex/include -Wno-deprecated-declarations -std=c++17 -c ACS.cpp -o acs.o
```

Second command: 
```bash
g++ Utility.o acs.o  -o acs -L/nfsd/rop/sw/ibm/cos221/cplex/lib/x86-64_linux/static_pic/ -L/nfsd/rop/sw/ibm/cos221/concert/lib/x86-64_linux/static_pic/  -lilocplex -lconcert -lcplex -lm -lpthread -Wno-deprecated-declarations
```

Then after compile all jobs thato you need, logout from arrow and return to login where you have to subbimt the jobs with `launch.py` in `/src` path so put this command

```bash
python3 launch.py
```
***

For cleaning the directory of log like `/out_csv`, `/out_problem` or `/out_terminal`, is possible to use the following command:
 
```bash
find out_csv out_problem out_terminal/*/* -type f -delete
```
***
In file Test.numbers you can find the general results.
