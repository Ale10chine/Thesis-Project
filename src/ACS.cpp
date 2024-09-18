#include <iostream>
#include <string>
#include <vector>
#include <random> 
#include <unordered_set>
#include <cstdlib> 
#include <iomanip> 
#include <cmath>
#include <chrono>
#include <fstream>
#include <thread> 

#include <ilcplex/ilocplex.h> // Cplex library
#include "Utility.hpp" // Utility library with some useful function 

// Only for test in cluster
//#include "/home/chinelloal/Thesis Project/include/Utility.hpp" 
//#include </nfsd/rop/sw/ibm/cos221/cplex/include/ilcplex/ilocplex.h>

// > FMIP, resolve a sub-MIP problem to find a vector of feasible solutions
// > OMIP, resolve another sub-MIP problem to improve the vector returned from FMIP with respect to 
//   the original objective form the initial MIP

ILOSTLBEGIN

// --- Foundamental functions useful for the building of the ACS algorithm ---
std::vector<IloNumArray> startV1(IloEnv &env, IloNumVarArray &v, IloIntArray &setI,
                                 const IloInt &n, const IloInt &m);

void printVector(IloEnv &env, const std::vector<IloNumArray> &v);

void printVector(IloEnv &env, const std::vector<IloInt> &v);

std::vector<IloInt> randomVecGenerator(int indX, int paramPercentage);

void printBounds(IloEnv &env, const std::vector<std::pair<IloNum, IloNum>> &vec);

void variableFixing(IloNumVarArray &v, const IloIntArray &setI, IloIntArray &setF,
                    const std::vector<IloInt> &rVec, const std::vector<IloNumArray> &sVec,
                    std::vector<std::pair<IloNum, IloNum>> &bounds, int paramPercentage);

void setFGenerator(const IloIntArray &setI, IloIntArray &setF,
                   const std::vector<IloInt> &rVec, int paramPercentage);

void updateVector(IloEnv &env, const IloIntArray &setI, std::vector<IloNumArray> &sVec,
                  const IloNumArray &vals, const IloInt &n, const IloInt &m);

void restoreBounds(IloNumVarArray &v, const IloIntArray &setF,
                   const std::vector<std::pair<IloNum, IloNum>> &bounds, int paramPercentage);

IloNum deltaCalculator(IloNum &delta, const std::vector<IloNumArray> &sVec, const IloInt &m);


IloNum deltaCalculator(const std::vector<IloNumArray> &sVec, const IloInt &m);

std::string primalGapCalculator(std::string objVal1, std::string objVal2);


// --- Useful for the conversion from numeric status to string ---
std::string getStatusDescription(IloAlgorithm::Status status);

// --- Function only useful for debugging process ---
void printModel(IloEnv &env, const IloNumVarArray &v, const IloRangeArray &r, const IloObjective &o); 
void printObj(IloEnv &env, const IloObjective &o);                                                    
void printVar(IloEnv &env, const IloNumVarArray &v);
void printRng(IloEnv &env, const IloRangeArray &r);

void printVector(IloEnv &env, const std::vector<IloNumArray> &v, IloInt n);


// Define de T parameter and C parameter, seed and percentage parameter 
// will be passed from the command line
constexpr int T_LIM = 30;
constexpr int C_ITER = 10;

// Main function
int main(int argc, char **argv)
{

    // Variable necessary for the printing on the log file (that could be .csv , .txt or .lp)
    std::string s;
    const char *msg;

    // Counters to keep track of the type of results to which ACS converges
    int resolved = 0;
    int timeLimit = 0;
    int iterLimit = 0;
    int error = 0; // (*)
    int errorF = 0; // (*)
    int except = 0; // (*)

    // Number of problem to resolve (cycle / 240 max) and p is the current indices of the problem
    int cycle; 
    int p = 0; 

    // Declaration of a vector of MIP problems to resolve
    std::vector<std::string> problems;
    
    // Setting parameters from command line input: the user can either select a specific instance by name
    // or choose a total number of problems based on their size
    int seed;    
    float percentage;
    if (argc == 3)
    {
        problems = readDirectory("../benchmark");
        cycle = 240;   
        seed = std::stoi(argv[1]);
        percentage = std::stoi(argv[2]) * 0.01; 
    }
    else if (argc == 4)
    {  
        std::string cmdline = readCommandline(argc, argv); 
        problems.push_back(cmdline);
        cycle = 1;
        seed = std::stoi(argv[2]);
        percentage = std::stoi(argv[3]) * 0.01;
    }
    else
    {
        std::cerr << "\nError entering the input data, try again in the following format: \n"
                  << "/executable \'seed\' \'percentage\' \n"
                  << "/executable \'problem.mps.gz\' \'seed\' \'percentage\' \n" 
                  << std::endl;
        return 1;
    }

    // Initialization of the random number generator with the specified seed
    std::srand(seed); 

    // Static decision for match a correct subpath for the directory of log like
    // out_csv, out_problem, out_terminal
    std::string subPath;

    if(seed == 12345){ // seed1
        if(percentage < 0.3){
            subPath = "seed1_p1/";
        }else if(percentage > 0.3 and percentage < 0.6){
            subPath = "seed1_p2/";
        }else if(percentage > 0.6){
            subPath = "seed1_p3/";
        }
    }

    if(seed == 54321){ // seed2
        if(percentage < 0.3){
            subPath = "seed2_p1/";
        }else if(percentage > 0.3 and percentage < 0.6){
            subPath = "seed2_p2/";
        }else if(percentage > 0.6){
            subPath = "seed2_p3/";
        }
    }

    if(seed == 1){ // seed3
        if(percentage < 0.3){
            subPath = "seed3_p1/";
        }else if(percentage > 0.3 and percentage < 0.6){
            subPath = "seed3_p2/";
        }else if(percentage > 0.6){
            subPath = "seed3_p3/";
        }
    }

    // Declaration of time variables to measure the computation time of the algorithm for a given problem
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> elapsed;


    // Open .csv file in append mode (only add new row without overwriting)
    std::ofstream csv("../out_csv/"+subPath+"Result.csv", std::ios::app);

    // If file doesn't yet exist or is void, write the header
    std::ifstream checkFile("../out_csv/"+subPath+"Result.csv");
    if (checkFile.peek() == std::ifstream::traits_type::eof())
    {
        csv << "\"ProblemName\",\"Type\",\"Status\",\"Time (sec)\",\"Iteration (n)\","
            << "\"Objective value (ACS)\",\"Objective value (optimal)\","
            << "\"PrimalGap\"\n";
    }
    checkFile.close();


    // --- Start of the reading the problem to compute with ACS ---
    for (; p < cycle; p++)
    {
        // Saving the path of the problem in a temporaneos variable
        const char *path = problems[p].c_str();
        std::cout<< problems[p]<< endl;
        
        //  Costruction of a CPLEX enviroment
        IloEnv env;

        // Open the log file .txt to write the terminal output of each problem 
        s = parserLog(path,("../out_terminal/"+subPath).c_str()) + ".txt";
        std::ofstream logFile(s);

        if (!logFile.is_open())
        {
            std::cerr << "Error opening ACS execution log file"<< std::endl;
            return 1;
        }

        // Initialization of string variable for the print in the .csv log file for data extraction
        std::string statusCSV = "";
        std::string timeCSV = "";
        std::string iteractionCSV = "";
        std::string objACSCSV = ""; 
        std::string objOpCSV = "";
        std::string typeCSV = "";
        std::string primalGapCSV ="";
        std::vector<std::string> csvRead; // For info about "BenchmarkSet.csv" file

        // --- Creation of the FMIP sub-problem and his parameters ---
        try
        {
            ///////++++++ FMIP ++++++////////
            IloModel fmip(env);
            IloCplex cplexFmip(env); // In teoria questo serve solo per riempire il modello

            // Declaration of some Extractables object for the FMIP sub-problem
            IloObjective objF(env);
            IloNumVarArray varF(env);
            IloRangeArray rngF(env);

            // Support variables for the print of the output (optional, useful for debugging)
            IloNumArray valsF(env);

            // Initialization of variables of FMIP by reading from file
            cplexFmip.importModel(fmip, path, objF, varF, rngF); 

            // -----------------------------------------------------------------------------------//

            ///////++++++ OMIP +++++++////////
            IloModel omip(env);
            IloCplex cplexOmip(env); // In teoria questo serve solo per riempire il modello

            // Declaration of some Extractables object for the OMIP model
            IloObjective objO(env);
            IloNumVarArray varO(env);
            IloRangeArray rngO(env);

            // Support variables for the print of the output (optional, useful for debugging)
            IloNumArray valsO(env);

            // Initialization of variables of OMIP by reading from file
            cplexOmip.importModel(omip, path, objO, varO, rngO);

            // ---------------------------------------------------------------------------------- //

            // # of costraints, note that this size is + 1 for the OMIP (m + 1)
            IloInt m = rngF.getSize(); 

            // Size of variable without slacks (Delta+ and Delta-) 
            IloInt n = varF.getSize();

            IloIntArray setI(env); // setI 
            IloIntArray setF(env); // setF (F c I) 

            // Random vector for the choice of what indices of I is going to be fixed in F
            std::vector<IloInt> randomVec;                      

            // Data structor for memorizing the starting bounds of variable before the fixing 
            std::vector<std::pair<IloNum, IloNum>> startBounds; 

            // DeltaUB costant for the m + 1 costraint in OMIP 
            IloNum deltaUB; 

            // Inizialization of Delta variables for FMIP 
            IloNumVarArray dpF(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);
            IloNumVarArray dnF(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);
            // and fot OMIP
            IloNumVarArray dpO(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);
            IloNumVarArray dnO(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);

            for (int i = 0; i < m; i++)
            {
                string name1 = "dplus" + to_string(i + 1);
                dpF[i].setName(name1.c_str());
                dpO[i].setName(name1.c_str());

                string name2 = "dmin" + to_string(i + 1);
                dnF[i].setName(name2.c_str());
                dnO[i].setName(name2.c_str());
            }

            // ---------------------- CREATION OF FMIP ---------------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            // Variables
            varF.add(dpF);
            varF.add(dnF);
            fmip.add(dpF);
            fmip.add(dnF);

            // -------------------------------------------------------------------
            
            // Obj. function of the FMIP : sum_{i = 0}^{m} (delta1[i] + delta2[i])
            fmip.remove(objF); // Remove the current obj.function

            IloExpr objFExpr(env);
            for (int i = 0; i < m; i++)
            {
                objFExpr += dpF[i] + dnF[i];
            }
            objF = IloMinimize(env, objFExpr, "MINIMIZE");
            fmip.add(objF);
           

            // -------------------------------------------------------------------
            
            // Constraints of the FMIP : A * x[i] + I_m * dpF[i] - I_m * dnF[i] = b
            IloExpr expr(env);
            for (int i = 0; i < m; i++)
            {
                expr = rngF[i].getExpr();

                if (rngF[i].getLb() == rngF[i].getUb()) // Caso '='
                {
                    expr += dpF[i] - dnF[i];
                }
                else if (rngF[i].getLb() > -IloInfinity)
                {
                    expr += dpF[i]; // Caso '>='
                }
                else if (rngF[i].getUb() < IloInfinity) // Caso '<='
                {
                    expr += -dnF[i];
                }
                else
                {
                    env.out() << "Error in the creation of FMIP" << endl;
                }

                rngF[i].setExpr(expr);
            }

            // ---------------------- CREATION OF OMIP ---------------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            // Variables
            omip.add(dpO);
            omip.add(dnO);
            varO.add(dpO);
            varO.add(dnO);

            // -------------------------------------------------------------------

            // Obj. function of the OMIP is the same of the initial MIP

            // -------------------------------------------------------------------

            // Constraints of the OMIP : A * x[i] + I_m * dpO[i] - I_m * dnO[i]
            //                           sum_{i = 0}^{m}(dpO[i] + dnO[i]) <= deltaUB
            
            IloExpr expr2(env); // For the additional costraint
            for (int i = 0; i < m; i++)
            {
                expr = rngO[i].getExpr();

                if (rngO[i].getLb() == rngO[i].getUb()) // Case '='
                {
                    expr += dpO[i] - dnO[i];
                    expr2 += dpO[i] + dnO[i];
                }
                else if (rngO[i].getLb() > -IloInfinity)
                {
                    expr += dpO[i]; // Case '>='
                    expr2 += dpO[i];
                }
                else if (rngO[i].getUb() < IloInfinity) // Case '<='
                {
                    expr += -dnO[i];
                    expr2 += dnO[i];
                }
                else
                {
                    env.out() << "Error in the creation of OMIP" << endl;
                }
                rngO[i].setExpr(expr);
                //expr2 += dpO[i] + dnO[i];
            }

            IloRange lastRange = (expr2 <= deltaUB);

            rngO.add(lastRange);
            omip.add(lastRange);

            // -------------------------------------------------------------------------------------

//---------------------------------------------- ACS --------------------------------------------------------//


            // ---------------- INITIALITATION OF STARTING VECTOR ----------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            // startV1 returns a starting vector and initialize the I set with olny integer variables
            // beacause they are the only possible variables to fix in this algorithm
            std::vector<IloNumArray> startingVector = startV1(env, varF, setI, n, m); 
            // printVector(env, startingVector);

            // Initializing the percentual parameter of how many variables ACS has to fix in each sub-Mip
            const int paramPercentage = setI.getSize() * percentage; 
            //--------------------------------------------------------------------------------------
            
            // Boolean variable for the first entrance in FMIP
            bool init = true; 
            // Initialising the DeltaUb variable for enter in the algorithm
            deltaUB = 0.1;
            // Initialising the iteration counter of the ACS
            int counter = 0;

            // Initialising variable for the management of the sub-Mip status
            IloCplex::Status statusF;
            IloCplex::Status statusO;
            IloAlgorithm::Status AstatusF;
            IloAlgorithm::Status AstatusO;
            bool solvOmip = true;
            bool solvFmip = true;

            // Variable for memorize the value of obj of the FMIP and OMIP 
            IloNum objFval;
            IloNum objOval;
            
            // Starting to mesure the time of computation of the ACS on problem i
            start = std::chrono::high_resolution_clock::now();

            while (deltaUB != 0 and counter < C_ITER) 
            {
                
                // FMIP condition, if deltaUB = 0 the FMIP is not processed, only OMIP
                if (deltaCalculator(startingVector, m) > 0 or init)
                {
                    // -------------------- RESOLUTION OF THE FMIP -----------------------
                    // -------------------------------------------------------------------
                    // -------------------------------------------------------------------

                    init = false; 

                    // Celaning the F set of the previous sub-Mip
                    setF.clear(); 

                    // In randomVec is stored the index of what variable we are going to fix from set I
                    // Note: startingVector[0].getSize() - 1 = setI.getSize()
                    randomVec = randomVecGenerator(startingVector[0].getSize() - 1, paramPercentage); 
                    
                    // Generation of set F with randomic values
                    setFGenerator(setI, setF, randomVec, paramPercentage);

                    // Variable fixing, the initial bounds are stored and then fixed with the initial
                    // values ​​of the starting vector or with the values ​​of the solutions after the first iteration
                    variableFixing(varF, setI, setF, randomVec, startingVector, startBounds, paramPercentage);

                    IloCplex cplexFmip(env);
                    cplexFmip.setParam(IloCplex::TiLim, T_LIM);

                    // Extraction of the model with C++ API, for the creation of the matrix that
                    // are needet from Cplex
                    cplexFmip.extract(fmip);


                    // Print the FMIP model in log file .lp
/*
                    s = parserLog(path, ("../out_problem/"+subPath+"FMIP").c_str(), true);
                    msg = s.c_str();
                    cplexFmip.exportModel(msg); 
*/
                    logAndPrint(logFile,"\n\n\n ---------------- RESOLUTION FMIP MODEL : " 
                                 + std::to_string(counter +1) + " --------------------\n\n\n");
                    
                    // Call to the Cplex resolutor
                    cplexFmip.solve();

                    // Save the status
                    statusF = cplexFmip.getCplexStatus(); // For timeLimit
                    AstatusF = cplexFmip.getStatus(); // For obj status
                    objFval = cplexFmip.getObjValue();
                    cplexFmip.getValues(valsF, varF);


                    // Managing of the FMIP status
                    if(AstatusF == IloAlgorithm::Infeasible){ 
                        solvFmip = false;
                        break;
                    }

                    if (statusF == IloCplex::AbortTimeLim)
                    {
                        if (AstatusF == IloAlgorithm::Feasible)
                        {
                            env.out() << " FMIP feasible" << endl;
                        }
                        else if (AstatusF == IloAlgorithm::Optimal)
                        {
                            env.out() << "FMIP optimal" << endl;
                        }
                        else
                        {
                            env.out() << "FMIP infeasible" << endl;
                            solvFmip = false;
                            break;
                        }
                    }
                    // Writing the statues and the objective value in the log .txt and terminal
                    logAndPrint(logFile,"\n\n\nSolution status = "+ getStatusDescription(AstatusF) + "\n");
                    logAndPrint(logFile,"Solution value = "+ std::to_string(objFval) + "\n");
                    

                    // Call to updateVector function to get the new update of the initial solution vector
                    // (the variable was called  startingVector), note that the solution vector,
                    // in x^, is composed only with integer values
                    updateVector(env, setI, startingVector, valsF, n, m); 

                    // Call to the restoreBounds function for restore the initial bounds of the variable
                    restoreBounds(varF, setF, startBounds, paramPercentage);
                    
                    // Cleaning of previously saved initial bounds 
                    startBounds.clear(); 

                    // Dealloco le matrici create usando extract()
                    cplexFmip.end();
                }

                //----------------------------------------------------------------------------------
                
                // ------------------------ UPDATE OF DELTAUB ------------------------
                // -------------------------------------------------------------------
                // -------------------------------------------------------------------

                
                deltaCalculator(deltaUB, startingVector, m);
                
                logAndPrint(logFile, "\n\n\n DeltaUB POST FMIP: " 
                            + std::to_string(deltaUB) + "\n\n");

                // Update of the rispective costants DeltaUB in the last costrint of OMIP 
                rngO[m].setBounds(rngO[m].getLb(), deltaUB); 

                // ---------------------------------------------------------------------------------
                
                // -------------------- RESOLUTION OF THE OMIP -----------------------
                // -------------------------------------------------------------------
                // -------------------------------------------------------------------

                // Celaning the F set of the previous sub-Mip
                setF.clear(); 

                // In randomVec is stored the index of what variable we are going to fix from set I
                // Note: startingVector[0].getSize() - 1 = setI.getSize()
                randomVec = randomVecGenerator(startingVector[0].getSize() - 1, paramPercentage); 
                
                // Generation of set F with randomic values
                setFGenerator(setI, setF, randomVec, paramPercentage);

                // Variable fixing, the initial bounds are stored and then fixed with the initial
                // values ​​of the starting vector or with the values ​​of the solutions after the first iteration 
                variableFixing(varO, setI, setF, randomVec, startingVector, startBounds, paramPercentage);


                IloCplex cplexOmip(env);
                cplexOmip.setParam(IloCplex::TiLim, T_LIM);
                // Extraction of the model with C++ API, for the creation of the matrix that
                // are needet from Cplex
                cplexOmip.extract(omip);

                // Print the OMIP model in log file .lp
/*
                s = parserLog(path, ("../out_problem/"+subPath+"OMIP").c_str(), true);
                msg = s.c_str();
                cplexOmip.exportModel(msg); 
*/

                logAndPrint(logFile,"\n\n\n ---------------- RESOLUTION OMIP MODEL: "
                             + std::to_string(counter +1) + " --------------------\n\n\n");

                    
                // Call to the Cplex resolutor and managing the status 
                if (cplexOmip.solve())
                {
                    statusO = cplexOmip.getCplexStatus(); // For timeLimit
                    AstatusO = cplexOmip.getStatus();     // For obj status
                    objOval = cplexOmip.getObjValue();    
                    cplexOmip.getValues(valsO, varO);     

                    // Writing the statues and the objective value in the log .txt and terminal
                    logAndPrint(logFile,"\n\n\nSolution status = "+ getStatusDescription(AstatusO) + "\n");
                    logAndPrint(logFile,"Solution value = "+ std::to_string(objOval) + "\n");
                    
                }
                else
                {
                    statusO = cplexOmip.getCplexStatus(); // For timeLimit
                    AstatusO = cplexOmip.getStatus();     // For obj status

                    // Writing the statues and the objective value in the log .txt and terminal
                    logAndPrint(logFile,"\n\n\nSolution status = "+ getStatusDescription(AstatusO) + "\n");
                    
                    if(AstatusO == IloAlgorithm::Infeasible){ // theoretically it should never happen
                        solvOmip = false;
                        break;

                    }
                }

                if (statusO == IloCplex::AbortTimeLim)
                {
                    env.out() << "I aborted in OMIP because I exceeded the time limit" << endl;
                    if (AstatusO != IloAlgorithm::Feasible and AstatusO != IloAlgorithm::Optimal)
                    {
                        solvOmip = false; // Time limit exceeded for this OMIP
                        break;
                    }

                } 

                // Call to updateVector function to get the new update of the initial solution vector
                // (the variable was called  startingVector), note that the solution vector,
                // in x^, is composed only with integer values
                updateVector(env, setI, startingVector, valsO, n, m); 


                // Call to the restoreBounds function for restore the initial bounds of the variable
                restoreBounds(varO, setF, startBounds, paramPercentage);

                // Cleaning of previously saved initial bounds 
                startBounds.clear(); 

                // Iteration limit counter grows
                counter++;

                // Print of the value at the before another itercation of the new FMIP
                logAndPrint(logFile, "\n\n\n DeltaUB post OMIP: " 
                            + std::to_string(deltaCalculator(startingVector, m)) +"\n\n");

                cplexOmip.end();
            }

            // End of mesuration of the time of computation of the ACS on problem i
            end = std::chrono::high_resolution_clock::now();
            // Calulating the time eplased from the start
            elapsed = end - start;

            // -------------------------------------------------------------------------------------
                
            // -------------------- MANAGING OF THE STATUS -----------------------
            // -------------------------------------------------------------------
            // ------------------------------------------------------------------- 
            if (!solvFmip) // theoretically it should never happen (*)
            {
                logAndPrint(logFile, "FMIP bug, infeasible\n");
                statusCSV = "Bug FMIP";
                errorF++;
            }
            else if (!solvOmip)
            {
                if (AstatusO == IloAlgorithm::Infeasible) // theoretically it should never happen (*)
                {
                    logAndPrint(logFile, "OMIP bug, infeasible\n");
                    statusCSV = "Bug OMIP";
                    error++;
                }
                else
                {
                    logAndPrint(logFile, "TimeLimit Exceeded for the OMIP --> Mip not resolved \n ");
                    logAndPrint(logFile, "Status OMIP :" + getStatusDescription(AstatusO) + "\n");
                    statusCSV = "Solution not found (time limit exceded)";
                    timeLimit++;
                }
            }
            else if (deltaUB == 0)
            {
                logAndPrint(logFile, "OMIP resolved and slack are equal to 0 --> Found feasible solution for Mip \n");
                logAndPrint(logFile, "Status OMIP :" + getStatusDescription(AstatusO) + "\n");
                statusCSV = "Feasible";
                objACSCSV = std::to_string(objOval);
                resolved++;
            }
            else if (counter == C_ITER)
            {
                logAndPrint(logFile, "Iteration limit exceded --> Mip not resolved \n");
                statusCSV = "Solution not found (iteration limit exceded)";
                iterLimit++;
            }
            else
            {
                logAndPrint(logFile, "Error \n");
            }

            // Managing and Print data for the extraction in .csv file
            s = parserLog(path,""); // In this case parser log return the name of the problem (es: markshare_4_0)
            timeCSV = std::to_string(elapsed.count());
            iteractionCSV = std::to_string(counter);
            
            csvRead = searchInstanceInCSV(s); // Use this funtion for saving two information useful for the extraction from "BenchmarkSet.csv"
            objOpCSV = csvRead[0]; // Mip optimal value already known and posted in MIPLIB 2017 benchmark set ("BenchmarkSet.csv")
            typeCSV = csvRead[1]; // Mip type like hard, open or easy

            // Computation of the primal gap to insert in .csv
            if(statusCSV == "Feasible"){
                primalGapCSV = primalGapCalculator(objACSCSV, objOpCSV);
            }

            // Call to function that actually print data on .csv
            csvPrintLine(csv,s,statusCSV,timeCSV,iteractionCSV,objACSCSV,objOpCSV,typeCSV,primalGapCSV);

            logAndPrint(logFile, "ACS ended for" + string(path) + "...\n\n" );
        }
        catch (IloException &e)
        {
            std::cerr << "Concert exception caught: " << e << endl;
            except++; // (*)
        }
        catch (...)
        {
            std::cerr << "Unknown exception caught" << endl;
        }
        
        // Close the .txt log
        logFile.close();
        // Frees the resources of the environment for the i problem
        env.end(); 
        
    }
    // Close the .csv file after the printation of all data from all problem
    csv.close();

    // The case signed with this symbol are (*) impossible, so the cheek is only for debugging
    std::cout<<"Resolved : "<<resolved<<endl;
    std::cout<<"Not resolved : "<<iterLimit <<" for Iteration limit exceeded"<<endl;
    std::cout<<"Not resolved : "<<timeLimit<<" for Time limit exceded"<<endl;
    //std::cout<<"Bug  "<<error<<" Infeasible OMIP"<<endl; 
    //std::cout<<"Bug  "<<errorF<<" Infeasible FMIP "<<endl; 
    //std::cout<<"Bug  "<<except<<" Magnitude exeeded "<<endl; 

    return 0;
}



// ------------------------------------------------------------------------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------------//
// ----------------------------------------- FUNCTIONS ------------------------------------------ // 
// -----------------------------------------------------------------------------------------------//

/*
    Starting Vector algorithm, sets the initial value of the x^ vector's value with randomic values
    taken from the range of the variable considerated in the iteration.
    It is foundamental to know n that is the size of the IlNumVar without Delta (slacks), it is foundamental
    to to pass for reference setI, that will be updated in this function
*/
std::vector<IloNumArray> startV1(IloEnv &env, IloNumVarArray &v, IloIntArray &setI, const IloInt &n, const IloInt &m)
{
    IloNumArray x(env);       // x^ vector, his size is the dimension of I 
    IloNumArray Delta1(env,m); // delta+ vector, his size is m
    IloNumArray Delta2(env,m); // delta- vector, his size is m

    // Choosing the randomic value for the integer variable that can be fixed
    for (IloInt i = 0; i < n; i++) // I scroll until the slack variable, Note that I is a subset of the first n IloNumVar
    {
        if (v[i].getType() == ILOBOOL || v[i].getType() == ILOINT)
        {
            IloInt lb = static_cast<IloInt>(v[i].getLb());
            IloInt ub = static_cast<IloInt>(v[i].getUb());
            
            // For limit big scale to the solver that could cause problem in during resolution
            if(ub > 100000) ub = 100000; // 100000 

            IloInt randomNumber = lb + (std::rand() % (ub - lb + 1)); // Generation of a pseudorando value between the bounds of the i variable
            x.add(static_cast<double>(randomNumber)); // Casting to double because x is an array of double (IloNumArray)
            setI.add(i);
     
        }
    }
    // Creation of the starting vector[x^,delta+, delta-] with the x's still to be fixed and empty Delta (slacks) 
    std::vector<IloNumArray> sV;
    sV.push_back(x);
    sV.push_back(Delta1);
    sV.push_back(Delta2);

    return sV;
}

/*
    Function to print std::vector<IloNumArray> vector like startingVector :
        
        +--------------------------------------------------------------+
        |  +------------+  +------------------+  +------------------+  |
        |  |  x_1       |  |  dplus_1         |  |  dmin_1          |  |
        |  |  x_2       |  |  dplus_2         |  |  dmin_2          |  |
        |  |  ...       |  |  ...             |  |  ...             |  |
        |  |  x_i       |  |  dplus_m         |  |  dmin_m          |  |
        |  +------------+  +------------------+  +------------------+  |
        +--------------------------------------------------------------+
*/
void printVector(IloEnv &env, const std::vector<IloNumArray> &v)
{
    env.out() << " [ ";

    for (int i = 0; i < v.size(); i++)
    {
        env.out() << "{ ";
        for (int j = 0; j < v[i].getSize(); j++)
        {
            env.out() << v[i][j] << " ";
        }
        env.out() << "}";
    }
    env.out() << " ]" << endl;
}
void printVector(IloEnv &env, const std::vector<IloNumArray> &v, IloInt n)
{
    env.out() << " [ ";

    for (int i = 1; i < v.size(); i++)
    {
        env.out() << "{ ";
        for (int j = 0; j < v[i].getSize(); j++)
        {
            env.out() << v[i][j] << " ";
        }
        env.out() << "}";
    }
    env.out() << " ]" << endl;
}
/*
    Another ovveride of the previous function, is useful for print vector like randomVec, that is
    a normal vector of IloInt (akaLong)
*/
void printVector(IloEnv &env, const std::vector<IloInt> &v)
{
    env.out() << " [ ";

    for (int i = 0; i < v.size(); i++)
    {
       env.out() << v[i] << " ";
    }
    env.out() << " ]" << endl;
}

/*
    randomVecGenerator returns random indices from the interval [ 0 - setI.getSize()], this is useful
    after to fix indices of I that will be stored in F
*/
std::vector<IloInt> randomVecGenerator(int indX, int paramPercentage)
{
    std::unordered_set<IloInt> randomSet;

    IloInt randomNumber;
    
    // Generation of a set of pseudorandom numbers, for this type of data structer they can't happen again
    while (randomSet.size() < paramPercentage)
    {
        randomNumber = std::rand() % (indX + 1); // indX = startingVec[0].getSize()-1 = setI.getSize() - 1
        randomSet.insert(randomNumber);
    }

    // Conversion of the set in a std::vector
    std::vector<IloInt> randomVec(randomSet.begin(), randomSet.end());
    return randomVec;
}

/*
    This function will return in output the print of the vector std::vector<std::pair<IloNum, IloNum>>
    like startBounds:

               x_1          x_2        x_3 ...........     x_i
        +-----------------------------------------------------------+
        |  +---------+  +---------+  +---------+       +---------+  |
        |  | (0, 1)  |  | (0, 1)  |  | (50, 60)| ....  | (3, 10) |  |
        |  +---------+  +---------+  +---------+       +---------+  |
        +-----------------------------------------------------------+
*/
void printBounds(IloEnv &env, const std::vector<std::pair<IloNum, IloNum>> &vec)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        std::cout << "Pair " << i << ": LB = " << vec[i].first << ", UB = " << vec[i].second << std::endl;
    }
}


/*
    Variable Fixing Algorithm that memorized the previously bounds in the startBounds vector
    and fix the new bounds to the variable scrolling the indices saved in setF
*/
void variableFixing(IloNumVarArray &v, const IloIntArray &setI, IloIntArray &setF,
                    const std::vector<IloInt> &rVec, const std::vector<IloNumArray> &sVec,
                    std::vector<std::pair<IloNum, IloNum>> &bounds, int paramPercentage)
{
    for (int i = 0; i < paramPercentage; i++)
    {
        bounds.push_back({v[setF[i]].getLb(), v[setF[i]].getUb()}); // Memorize Lb and Ub from the variable with indices £ F
        IloInt a = IloRound(sVec[0][rVec[i]]); // IloInt <- IloNum, fixing variable are surely integer

        // (Optional prevention) security control in case of particular error from cplex
        if(a < v[setF[i]].getLb()){
            a = v[setF[i]].getLb();
        }else if(a > v[setF[i]].getUb()){
            a = v[setF[i]].getUb();
        }

        // Setting the new bounds, the varible with indices £ F is now fixed
        v[setF[i]].setBounds(a,a); 
       
    }
}

/*
    Function that modifies for refercence the content of setF, in other words, this function
    fills the setF
*/
void setFGenerator(const IloIntArray &setI, IloIntArray &setF,
                   const std::vector<IloInt> &rVec, int paramPercentage)
{
    for (int i = 0; i < paramPercentage; i++)
    {
        setF.add(setI[rVec[i]]); // Memorize of the index of I in F (F c I)
    }
}

/*
    UpdateVector algorithm, is used for update the initial vector , in main called "startingVector",
    for understand his procedure we have to print it. The fact is that it's a "game" of indexes so it's
    difficoult to explain with only words
*/
void updateVector(IloEnv &env, const IloIntArray &setI, std::vector<IloNumArray> &sVec,
                  const IloNumArray &vals, const IloInt &n, const IloInt &m)
{

    // Update x^ vector
    // env.out() << "n :" << n << " m : " << m << endl;
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        if (setI[j] == i)
        {
            sVec[0][j] = vals[i];
            j++;
        }
        if (j == setI.getSize())
        {
            break;
        }
    }

    // Update slacks vector Delta+ and Delta-
    j = 0;
    int c = 0;
    for (int i = n; i < 2 * m + n; i++, j++)
    {
        if (j < m)
        {
            sVec[1][j] = vals[i];
        }
        else
        {
            sVec[2][c] = vals[i];
            c++;
        }
    }
}

/*
    Function who restore initial bounds usign reference with the variables of IloNumVarArray vars of
    FMIP or OMIP problem
*/
void restoreBounds(IloNumVarArray &v, const IloIntArray &setF,
                   const std::vector<std::pair<IloNum, IloNum>> &bounds, int paramPercentage)
{
    for (int i = 0; i < paramPercentage; i++)
    {
        if (bounds[i].first < bounds[i].second)
        {   
            v[setF[i]].setBounds(bounds[i].first, bounds[i].second);
        }
        else
        {
            v[setF[i]].setBounds(bounds[i].second, bounds[i].first);
        }
    }
}
/*
    Function for deltaUB computation
*/
IloNum deltaCalculator(const std::vector<IloNumArray> &sVec, const IloInt &m){

    IloNum delta = 0; // DeltaUB wil be restored for the new computation
    for (int i = 0; i < m; i++)
    {
        IloNum tmp = sVec[1][i] + sVec[2][i];
        delta = delta + tmp;
    }

    return delta;
}
/*
    Function for deltaUB computation (version with modify for reference)
*/
IloNum deltaCalculator(IloNum &delta, const std::vector<IloNumArray> &sVec,const IloInt &m)
{
    delta = 0; // DeltaUB wil be restored for the new computation
    for (int i = 0; i < m; i++)
    {
        IloNum tmp = sVec[1][i] + sVec[2][i];
        delta = delta + tmp;
    }

    return delta;
}

/*
    Simple function for convert the status from numeric to string
*/
std::string getStatusDescription(IloAlgorithm::Status status) {
    switch (status) {
        case IloAlgorithm::Unknown:
            return "Unknown";
        case IloAlgorithm::Feasible:
            return "Feasible";
        case IloAlgorithm::Optimal:
            return "Optimal";
        case IloAlgorithm::Infeasible:
            return "Infeasible";
        case IloAlgorithm::Unbounded:
            return "Unbounded";
        case IloAlgorithm::InfeasibleOrUnbounded:
            return "Infeasible or Unbounded";
        case IloAlgorithm::Error:
            return "Error";
        default:
            return "Other Unkonwn problems";
    }
}

/*
    Function for computation of the primal gap between feasible solution of the ACS
    algorithm and optimal MIP solution of the same problem
*/
std::string primalGapCalculator(std::string objVal1, std::string objVal2){

    double primalGap;

    double oACS = stod(objVal1); // c^T*x^
    double oOp = stod(objVal2);  // c^T*x

    if(std::fabs(oACS) == std::fabs(oOp)){  //case |c^T*x^| = |c^T*x^|
        primalGap = 0;
    }else if(oACS * oOp < 0){  // case c^T*x^ * c^T*x < 0
        primalGap = 1;
    }else {  // otherwise
        primalGap = (std::fabs(oACS - oOp))/(std::fmax(std::fabs(oACS),std::fabs(oOp)));
    }
    
    return std::to_string(primalGap);

}

// -----------------------------------------------------------------------------------------------//
// ---------------------------- FUNCTION FOR DEBUGGING -------------------------------------------//
// -----------------------------------------------------------------------------------------------//

/*
    Function to print variable of the model 
*/
void printModel(IloEnv &env, const IloNumVarArray &v, const IloRangeArray &r, const IloObjective &o){
    // Obj
    env.out() <<"Name: "<<o.getName() << " Sense: " << o.getSense() << " Expr: " << o.getExpr() << endl << endl;
    // Vars
    for (int i = 0; i < v.getSize(); i++)
    {
        env.out() << v[i].getName() << " " << v[i].getType() << "  " << v[i].getId() <<"   " 
                  << v[i].getLb() << " <= " << v[i].getName() << " <= " <<v[i].getUb()<<endl;
        
    }
    env.out() << endl;
    // Constraints
    for (int i = 0; i < r.getSize(); i++)
    {
        env.out() << r[i].getExpr(); 
        if(r[i].getLb() == r[i].getUb()){
            env.out() << " = " << r[i].getLb();
        } else if(r[i].getLb() > -IloInfinity){
            env.out() << " >= " << r[i].getLb();
        } else if (r[i].getUb() < IloInfinity){
            env.out() << " <= "<< r[i].getUb();
        } else {
            env.out() << "No bounds for this constraint";
        }
        env.out() << endl << endl;
    }
    env.out() << endl;

}
/*
    Function to print only function objective of the model
*/
void printObj(IloEnv &env, const IloObjective &o)
{
    env.out() << "Name: " << o.getName() << " Sense: " << o.getSense() << " Expr: " << o.getExpr() << endl
              << endl;
}

/*
    Function to print only variable of the model
*/
void printVar(IloEnv &env, const IloNumVarArray &v)
{
    for (int i = 0; i < v.getSize(); i++)
    {
        env.out() << v[i].getName() << " " << v[i].getType() << "  " << v[i].getId() << "   " 
                  << v[i].getLb() << " <= " << v[i].getName() << " <= " << v[i].getUb() << endl;
       
    }
    env.out() << endl;
}

/*
    Function to print only costraints of the model
*/
void printRng(IloEnv &env, const IloRangeArray &r)
{
    // Scrool every range 
    for (int i = 0; i < r.getSize(); i++)
    {
        // Extraction of the expression from the current costraint
        IloExpr expr = r[i].getExpr();
        bool first = true;

        // Iteration on terms of the i-range
        for (IloExpr::LinearIterator it = expr.getLinearIterator(); it.ok(); ++it)
        {
            // Get current variable
            IloNumVar v = it.getVar();

            // And get her coefficent 
            double coeff = it.getCoef();

            // Print the variable with respective id (useful for duplicate variable debugging) 

            if (!first)
            {
                env.out() << " + ";
            }
            first = false;

            env.out() << coeff << " * " << v.getName() << "(" << v.getId() << ")";
        }

        //  Print the right side of the costraint
        if (r[i].getLb() == r[i].getUb())
        {
            env.out() << " = " << r[i].getLb();
        }
        else if (r[i].getLb() > -IloInfinity)
        {
            env.out() << " >= " << r[i].getLb();
        }
        else if (r[i].getUb() < IloInfinity)
        {
            env.out() << " <= " << r[i].getUb();
        }
        else
        {
            env.out() << "No bounds for this constraint";
        }

        env.out() << endl
                  << endl;

        //  Release of the expression for a correct management of the memory
        expr.end();
    }
}

/*
g++ -std=c++17 -c Utility.cpp -o Utility.o
g++ -std=c++17 -c ACS.cpp -o acs.o


g++ Utility.o ACS.o  -o acs -lilocplex -lconcert -lcplex -lm -lpthread -Wno-deprecated-declarations


./acs markshare_4_0.mps.gz
*/

