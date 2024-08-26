#include <iostream>
#include <string>
#include <vector>
#include <ilcplex/ilocplex.h>

#include "Utility.hpp"


// Function FMIP, resolve a sub-MIP problem to find a vector of feasible solutions


// Function OMIP, resolve another sub-MIP problem to improve the vector returned from FMIP with respect to the original objective form the initial MIP



ILOSTLBEGIN



int main(int argc, char **argv){

    // Costruction of a CPLEX enviroment
    IloEnv env;

    // Declaration of a vector of MIP problems to resolve
    std::vector<std::string> problems;

    std::string cmdline = readCommandline(argc, argv);
    problems.push_back(cmdline);
    const char *path = problems[0].c_str();


    try
    {
        IloModel mip(env);
        IloCplex cplexMip(env);

        // Declaration of some Extractables object for the general MIP model
        IloObjective obj(env); 
        IloNumVarArray var(env);
        IloRangeArray rng(env);

        // (Decmpression and) Reading from a file .mps or .gz and saving the data from the problem (f.obj, variables, constraints, ...)
        cplexMip.importModel(mip, path, obj, var, rng);






        //#### Definizione delta 1 e 2 per usarli nel FMIP
        IloInt m = rng.getSize();

        IloNumVarArray delta1(env, m + 1, 0.0 , IloInfinity, IloNumVar::Type::Float);
        IloNumVarArray delta2(env, m + 1, 0.0 , IloInfinity, IloNumVar::Type::Float);

        

        std::cout << "Objective function of General MIP " << obj.getExpr() << std::endl;
        // Creation of the FMIP sub-problem and his parameters
        try
        {
            IloModel fmip(env);
            IloCplex cplexFmip(env);

            // Declaration of some Extractables object for the FMIP sub-problem
            IloObjective objF(env);
            IloNumVarArray varF(env);
            IloRangeArray rngF(env);

            IloExpr objFExpr(env);
            for (int i = 0; i <= m; i++)
            {
                objFExpr = delta1[i] + delta2[i];
            }

            fmip.add(IloMinimize(env, objFExpr));
        }
        catch (IloException &e)
        {
            std::cerr << "Concert exception caught: " << e << endl;
            env.end();
        }
        catch (...)
        {
            std::cerr << "Unknown exception caught" << endl;
            env.end();
        }

        // Creation of the OMIP sub-problem and his parameters
        try
        {
            IloModel omip(env);
            IloCplex cplexOmip(env);

            // Declaration of some Extractables object for the general MIP model
            IloObjective objO(env);
            IloNumVarArray varO(env);
            IloRangeArray rngO(env);

            objO = obj;
            std::cout << "Objective function of OMIP "<< objO.getExpr() << std::endl;
        }
        catch (IloException &e)
        {
            std::cerr << "Concert exception caught: " << e << endl;
            env.end();
        }
        catch (...)
        {
            std::cerr << "Unknown exception caught" << endl;
            env.end();
        }
        /*
                while ()
                {

                }
        */
    }
    catch (IloException &e)
    {
        std::cerr << "Concert exception caught: " << e << endl;
        env.end();
    }
    catch (...)
    {
        std::cerr << "Unknown exception caught" << endl;
        env.end(); 
    }

    // Frees the resources of the environment
    env.end();

    return 0;
}




