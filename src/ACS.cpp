#include <iostream>
#include <string>
#include <vector>
#include <random> // Forse se non serve più toglila, usiamo cstdlib
#include <unordered_set>
#include <cstdlib> // Per rand() e srand()
#include <iomanip> // per il std::fixed

#include <ilcplex/ilocplex.h>

#include "Utility.hpp"

// FMIP, resolve a sub-MIP problem to find a vector of feasible solutions
// OMIP, resolve another sub-MIP problem to improve the vector returned from FMIP with respect to the original objective form the initial MIP

/* Commenti e note:
 ####Come ottenere la stampa del vincolo intero con anche il GE, LE , EQ per ora non ci interessa
 Nota: non è automaticamente implementata la deep copy tra oggetti di cplex
 Nota: Con addMIPStart imponi solo un suggerimento al risolutore per dove partire a risolvere il problema, ma di fatto non fissi le variabili
 Teoricamente CPLEX riconosce automaticamente il fatto che certe variabili appartengono ad un determinato modello vedendo implicitamente come vengono
 definiti vincoli ecc gestendo però tre modelli in un unico enviroment è buona prassi agiungere ogni ogetto Extractables dentro il modello rispettivo
 Nota: Se non setti il name nella funzione obbiettivo, e provi a stamparla ti torna un seg fouult a causa di un puntatore nullptr
 Nota: In cplex non esisteun modo per recuperare variabili, vincoli ecc da un modello, bisogna tenerne traccia separatamente
*/
ILOSTLBEGIN

    //void cloneVar(IloEnv &env, const IloNumVarArray &a, IloNumVarArray &b);                                                 // copy from a to b (variables)
    //void cloneRng(IloEnv &env, const IloRangeArray &a, IloRangeArray &b, const IloNumVarArray &varA, IloNumVarArray &varB); // copy froma a to b (constraints)
    //void cloneObj(IloEnv &env, const IloObjective &a, IloObjective &b);                                                     // copy from a to b (function object)
void printModel(IloEnv &env, const IloNumVarArray &v, const IloRangeArray &r, const IloObjective &o);                   // Get a view of the model
void printObj(IloEnv &env, const IloObjective &o);                                                                      // View more detailed for each component of the model
void printVar(IloEnv &env, const IloNumVarArray &v, const IloCplex &cplex);
void printRng(IloEnv &env, const IloRangeArray &r);

// Fondamentale questa n perchè ci dice quanto è lungo il vettore di var senza l'aggiunta delle slack
// Fondamentale anche il passaggio per reference di setI, che verrà aggiornato durante il programma
std::vector<IloNumArray> startV1(IloEnv &env, IloNumVarArray &v, IloIntArray &setI, const IloInt &n, const IloInt &m);
void printVector(IloEnv &env, const std::vector<IloNumArray> &v);
void printVector(IloEnv &env, const std::vector<IloInt> &v);

std::vector<IloInt> randomVecGenerator(int indX, int paramPercentage);                                 // Passo il range, ovvero tra i possibili estremi di x^ cioè tra i suoi indici [ 0 - n-1 ]
void printBounds(IloEnv &env, const std::vector<std::pair<IloNum, IloNum>> &vec); // Stampa i range che ogni volta vengono salvati delle variabili da fissare
// Funzione per il fissaggio delle variabili, oltre a queste ultime modifica F e i bounds di startingBounds
void variableFixing(IloNumVarArray &v, const IloIntArray &setI, IloIntArray &setF, const std::vector<IloInt> &rVec, const std::vector<IloNumArray> &sVec, std::vector<std::pair<IloNum, IloNum>> &bounds, IloCplex &cplex, int paramPercentage);
void setFGenerator(const IloIntArray &setI, IloIntArray &setF, const std::vector<IloInt> &rVec, int paramPercentage);

void updateVector(IloEnv &env, const IloIntArray &setI, std::vector<IloNumArray> &sVec, const IloNumArray &vals, const IloInt &n, const IloInt &m);
void restoreBounds(IloNumVarArray &v, const IloIntArray &setF, const std::vector<std::pair<IloNum, IloNum>> &bounds,IloCplex &cplex, int paramPercentage);
IloInt deltaCalculator(IloInt delta, const std::vector<IloNumArray> &sVec,const IloInt &m);


double getDirectModelLB(int index, CPXENVptr c_env, CPXLPptr c_lp); //NOTA: Ho impelmentato queste due funzini in quanto per bypassare il set bounds
double getDirectModelUB(int index, CPXENVptr c_env, CPXLPptr c_lp);
void setDirectModelBounds(int index, CPXENVptr c_env, CPXLPptr c_lp, char what, IloNum val);


int main(int argc, char **argv)
{
    const int seed = 12345;
    std::srand(seed); // Inizializza il generatore di numeri casuali

    // IloEnv envI;
    //  Costruction of a CPLEX enviroment
    IloEnv env;

    // Declaration of a vector of MIP problems to resolve
    std::vector<std::string> problems;

    // Variabili necessari per la stampa in un file di log
    std::string s;
    const char *msg;

    std::string cmdline = readCommandline(argc, argv);
    problems.push_back(cmdline);
    const char *path = problems[0].c_str();

    // Creation of the FMIP sub-problem and his parameters
    try
    {
        ///////++++++ FMIP ++++++////////
        IloModel fmip(env);
        IloCplex cplexFmip(env);

        // Declaration of some Extractables object for the FMIP sub-problem
        IloObjective objF(env);
        IloNumVarArray varF(env);
        IloRangeArray rngF(env);
        
        //Variabili di appoggio per la stampa dell'output
        IloNumArray valsF(env);

        // Riempimento variabili FMIP
        cplexFmip.importModel(fmip, path, objF, varF, rngF); // Importando il modello non serve riaggiungere le variabili

// -------------------------------------------------------------------------------------------------------------- //

        ///////++++++ OMIP +++++++////////
        IloModel omip(env);
        IloCplex cplexOmip(env);

        // Declaration of some Extractables object for the OMIP model
        IloObjective objO(env);
        IloNumVarArray varO(env);
        IloRangeArray rngO(env);
        
        // variabili di appoggio per la stampa dell'output
        IloNumArray valsO(env);

        // riempimento variabili OMIP
        cplexOmip.importModel(omip, path, objO, varO, rngO);


// ---------------------------------------------------------------------------------------------------------------- //

        // #### Definizione delta 1 e 2 (variabili) per usarli nel FMIP e nel OMIP
        IloInt m = rngF.getSize(); // # dei vincoli è + 1 nel OMIP (NOTA)

        // Mi salvo qui la grandezza delle variabili cosi posso fissare solo queste, i delta non li devo fissare nel FMIP!!
        IloInt n = varF.getSize();

        IloIntArray setI(env); // setI (F c I) :  array per memorizzare dove si trovano le variabili inzizialemnte (solo quelle fissabili, cioè solo var intere)
        IloIntArray setF(env); // setF (F c I) : array per memorizzare solo la percentuale di incidi selezionati in modo randomico che fanno parte di I

        std::vector<IloInt> randomVec; // Vettore randomico per la scelta di quali indici di I fissare in F
        std::vector<std::pair<IloNum, IloNum>> startBounds; // Struttura dati che mi permette di memeorizzare i bound inizziali durante il fissaggio delle variabili

        IloNum deltaUB; // Vettore che fornirà il bound per il vincolo aggiuntivo dell' OMIP

        IloNumVarArray dpF(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);
        IloNumVarArray dnF(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);

        IloNumVarArray dpO(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);
        IloNumVarArray dnO(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);

        env.out() << varF.getSize() << endl; // Per vedere quanto grande è il modello fino a qua

        // Assegnazione nome alle variabili
        for (int i = 0; i < m; i++)
        {
            string name1 = "dplus" + to_string(i + 1);
            dpF[i].setName(name1.c_str());
            dpO[i].setName(name1.c_str());

            string name2 = "dmin" + to_string(i + 1);
            dnF[i].setName(name2.c_str());
            dnO[i].setName(name2.c_str());
        }

        // ------------------------- CREAZIONE VARIABILI FMIP-----------------
        // -------------------------------------------------------------------
        // -------------------------------------------------------------------

        varF.add(dpF);
        varF.add(dnF);
        fmip.add(dpF);
        fmip.add(dnF);
        //fmip.add(varF);

        //-------------------- CREAZIONE VINCOLI FMIP --------------------
        // -------------------------------------------------------------------
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
               expr += dpF[i];   // Caso '>='
            }
            else if (rngF[i].getUb() < IloInfinity) // Caso '<='
            {
               expr += - dnF[i];
            }
            else
            {
                env.out() << "Errore nella creazione di FMIP" << endl;
            }

            rngF[i].setExpr(expr);
        }
        // fmip.add(rngF);

        //-------------------- CREAZIONE FUNZIONE OBBIETTIVO FMIP--------------
        // -------------------------------------------------------------------
        // -------------------------------------------------------------------

        // Obj. function of the FMIP : sum_{i = 0}^{m} (delta1[i] + delta2[i])
        fmip.remove(objF); // Tolgo la funzione obbiettivo già presente

        IloExpr objFExpr(env);
        for (int i = 0; i < m; i++)
        {
            objFExpr += dpF[i] + dnF[i];
        }
        objF = IloMinimize(env, objFExpr, "MINIMIZE");
        fmip.add(objF);

        // Stampa modello con modifiche su costraints e function object
//        env.out() << "Modello FMIP dopo la creazione: " << endl;
//        printModel(env, varF, rngF, objF);

// -------------------------------------------------------------------------------------------------------------------//

        //-------------------- CREAZIONE VARIABILI OMIP--------------------
        // -------------------------------------------------------------------
        // -------------------------------------------------------------------
        omip.add(dpO);
        omip.add(dnO);
        varO.add(dpO);
        varO.add(dnO);

        //-------------------- CREAZIONE VINCOLI OMIP --------------------
        // -------------------------------------------------------------------
        // -------------------------------------------------------------------

        // Constraints of the OMIP : A * x[i] + I_m * dpO[i] - I_m * dnO[i]
        //                           sum_{i = 0}^{m}(dpO[i] + dnO[i]) <= deltaUB
        // nota expr già definita in precedenza riporta su nelle dichiarazioni anche expr2
        IloExpr expr2(env);
        for (int i = 0; i < m; i++)
        {
            expr = rngO[i].getExpr();

            if (rngO[i].getLb() == rngO[i].getUb()) // Caso '='
            {
               expr += dpO[i] - dnO[i]; 
            }
            else if (rngO[i].getLb() > -IloInfinity)
            {
               expr += dpO[i];   // Caso '>='
            }
            else if (rngO[i].getUb() < IloInfinity) // Caso '<='
            {
               expr += - dnO[i];
            }
            else
            {
                env.out() << "Errore nella creazione di FMIP" << endl;
            }

            rngO[i].setExpr(expr);

            expr2 += dpO[i] + dnO[i];
        }

        IloRange lastRange = (expr2 <= deltaUB);

        rngO.add(lastRange);
        omip.add(lastRange);

        //-------------------- CREAZIONE FUNZIONE OBBIETTIVO OMIP-------------
        // -------------------------------------------------------------------
        // -------------------------------------------------------------------

        // Non serve fare modifiche in quanto rimane quela del problema iniziale


        // Stampa modello con modifiche su costraints e aggiunta costaint extra 
//        env.out() << "Modello OMIP dopo la creazione: " << endl;
        //printModel(env,varO,rngO,objO);
  
  

//--------------------------------------------------------------------------------------------------------------------//
        // ######### Starting vector for the FMIP ########### Poi qua si può implementare ALG2
        // startV1 modifica sia le varF con reference che il setI
        std::vector<IloNumArray> startingVector = startV1(env, varF, setI, n, m); // Ritorna uno starting vector e modifica il setI aggiungendo gli indici relativi alle variabili intere tale per cui poter essere fissate
        env.out() << "Starting Vector: " << endl;
        printVector(env, startingVector);

        const int paramPercentage = setI.getSize()*0.2; //  PARAM : fisso ogni volta n variabili
        std::cout<<"paramPercentage :"<<paramPercentage<<endl;

        bool start = true; // Boleana che serve solo per entrare nel primo FMIP
//---------------------------------------------------------------------------------------------------------------------//
        deltaUB = 1;
        int counter = 1;
        
        std::cout<<"Dimensione starting vector: "<<startingVector[0].getSize()<<endl;
        std::cout<<"Dimensione di I "<<setI.getSize()<<endl;

        IloAlgorithm::Status status;

        while (status != IloAlgorithm::Optimal or deltaUB > 0)
        {
            
            if (deltaCalculator(deltaUB,startingVector,m) > 0 || start)
            {
                start = false; // Variabile che rimane accesa solo all'inizio

                // ######## Fissaggio variabili ##################
                setF.clear(); // Svuoto il precedente setF
                // Si crea in F (setF) una percentuale di indici randomici, ai quali verranno fissate le variabili scelte dallo starting vector
                randomVec = randomVecGenerator(startingVector[0].getSize() - 1, paramPercentage); // Questa grandezza coincide con la dimensione di I
                // Genero un F randomico
                setFGenerator(setI, setF, randomVec, paramPercentage);
                //cout<<"\n\nDimensione di F "<<setF.getSize()<<endl; 
                //cout<<"F pre FMIP= "<<setF<<endl;

                /*
                s = parserLog(path, "../out/BFixFMIP", true);
                msg = s.c_str();
                cplexFmip.exportModel(msg); // Esportazione output in un file di log dedicato
                */

                cplexFmip.extract(fmip);
                // Funzione per il fissaggio delle variabili, mantiene in memoria i bounds inizziali e li modifica, modifica F
                variableFixing(varF, setI, setF, randomVec, startingVector, startBounds, cplexFmip, paramPercentage);

                env.out()<<" (1) Variabili dopo il fissaggio FMIP\n"<<endl;
                printVar(env, varF, cplexFmip);

  //            env.out() << " \nF = " << setF << endl;
  //              env.out() << "\nRandVec = ";
  //              printVector(env, randomVec);

//                env.out() << "\nstartBounds : " << endl;
//                printBounds(env, startBounds);
//                std::cout << "\n\n\n";

//                env.out() << "Modello FMIP dopo il fissaggio delle variabili e la creazione dello starting vector: " << endl;
//                printModel(env, varF, rngF, objF);

                // -------------------- RISOLUZIONE + FILE LOG FMIP--------------------
                // -------------------------------------------------------------------
                // -------------------------------------------------------------------
        //        cout<<"Starting vector Prima della risoluzione del FMIP"<<endl;
        //        printVector(env, startingVector);
        //        cout<<"Variabili Modello FMIP: "<<endl;
        //        printVar(env,varF);

                
                // Risoluzione del modello
                //cplexFmip.extract(fmip);
    /*
            env.out() << "\n\n (A3) Variabile dopo EXTRAXT FMIP + altra variabile" << endl;

                     env.out() << varF[8763].getName() << " " << varF[8763].getType() << "  " << varF[8763].getId() << "   " << varF[8763].getLb() << " <= " << varF[8763].getName() << " <= " << varF[8763].getUb();
            env.out() << endl;
            env.out() << varF[8759].getName() << " " << varF[8759].getType() << "  " << varF[8759].getId() << "   " << varF[8759].getLb() << " <= " << varF[8759].getName() << " <= " << varF[8759].getUb();
            env.out() << endl;  
    */
                // Stampa output in file di log

                s = parserLog(path, "../out/FMIP", true);
                msg = s.c_str();
                cplexFmip.exportModel(msg); // Esportazione output in un file di log dedicato

                if (counter == 1)
                {
                // Stampa output in file di log
               
                s = parserLog(path, "../out/1FMIP", true);
                msg = s.c_str();
                cplexFmip.exportModel(msg); // Esportazione output in un file di log dedicato
                }

                std::cout << "\n\n\n ---------------- RISOLUZIONE MODELLO FMIP : " <<counter<<" --------------------\n\n\n";
                

                //cplexFmip.extract(fmip);
                cplexFmip.solve(); // Risoluzione di FMIP
                cplexFmip.getValues(valsF, varF);

                env.out() << "\n\n\n\nSolution status = " << cplexFmip.getStatus() << std::endl;
                env.out() << "Solution value  = " << cplexFmip.getObjValue() << "\n\n\n"
                          << endl;
//                env.out() << "Values = " << valsF << std::endl;
                env.out() << "Dimensione varF: " << varF.getSize() << endl;
                env.out() << "Dimesnione valsF: " << valsF.getSize() << endl;

                //env.out()<<valsF<<endl;
//                cout << "Lo starting vector era cosi' inzizialmente" << endl;
//                printVector(env, startingVector);

//                env.out() << "Dimesnione di I : " << setI.getSize() << endl;
                updateVector(env, setI, startingVector, valsF, n, m); // Riempio lo starting Vector filtrando i risultati di vals dati da cplex

        //        env.out() << "Starting Vector (post risoluzione FMIP): " << endl;
        //        printVector(env, startingVector);

                // Finito di risolvere, bisogna risettare le variabili come all'inziio, sfruttando startBounds
                // Dopodichè si restituirà il vettore inziziale come output [x^,Delta+, delta-] dove
                // Averrà una Ricerca nell'output, filtrando solo le variabili intere.
                // Inoltre bisognerà salvarsi i valori delle soluzioni delta che serviranno per costruire DeltaUB
                cplexFmip.extract(fmip);

                // Ripristino i bounds iniziali delle variabili del problema
                restoreBounds(varF, setF, startBounds,cplexFmip ,paramPercentage);
                env.out()<<"(2) Variabili dopo il restore bounds FMIP\n"<<endl;
                printVar(env, varF, cplexFmip);
//                env.out() << "FMIP riprristinato ai bounds originari" << endl;
//                printModel(env, varF, rngF, objF);

                startBounds.clear(); // Pulisco i bounds memorizzati in favore del riempimento ch si farà in OMIP
//                env.out() << "Bounds dopo la pulizia: " << endl;
//                printBounds(env, startBounds);
            }

            //--------------------------------------------------------------------------------------------------------------------------------------------//

            // Creazione di DeltaUB da passare come parametro poi per il vincolo del OMIP, sarà anche la condizione necessaria per accedere al FMIP

            deltaUB = deltaCalculator(deltaUB,startingVector,m);
            env.out()<<"-------DeltaUB POST FMIP-----------: "<<deltaUB<<endl<<endl;
            rngO[m].setBounds(rngO[m].getLb(), deltaUB); // Aggiorno il DeltaUB nell'OMIP, NOTA questo setBOUNDS non è il metodo che da problemi, è su IloRange

            // (IMPORTANTE, ragiona su queli elementi farla) Pulizia variabili o strutture di memoria che devono essere svuotate (NOTA)

            // -------------------------------------------------------------------------------------------------------------------------------------------//

            // ######## Fissaggio variabili ##################
            setF.clear(); // Svuoto il precedente setF

            randomVec = randomVecGenerator(startingVector[0].getSize() - 1, paramPercentage); // Questa grandezza coincide con la dimensione di I
            // Randomizzazione di F
            setFGenerator(setI, setF, randomVec, paramPercentage);
            //cout << "\n\nDimensione di F " << setF.getSize() << endl;
            //cout<<"F pre OMIP = "<<setF<<endl;


            cplexOmip.extract(omip);


            s = parserLog(path, "../out/BFixOMIP", true);
            msg = s.c_str();
            cplexFmip.exportModel(msg); // Esportazione output in un file di log dedicato

            variableFixing(varO, setI, setF, randomVec, startingVector, startBounds, cplexOmip, paramPercentage);

            env.out()<<"(3) Variabili dopo il fissaggio OMIP\n"<<endl;
            printVar(env, varO, cplexOmip);
 //           cout << "Set F in OMIP = " << setF << endl;

 //           env.out() << "Modello OMIP dopo il fissaggio delle variabili : " << endl;
 //           printModel(env, varO, rngO, objO);

            // -------------------- RISOLUZIONE + FILE LOG OMIP--------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------
/*
            env.out() << "\n\n(2) Variabile Bastarda DOPO DEL FISSAGIO OMIP e prima dell'estrazione del modello" << endl;

            env.out() << varO[8763].getName() << " " << varO[8763].getType() << "  " << varO[8763].getId() << "   " << varO[8763].getLb() << " <= " << varO[8763].getName() << " <= " << varO[8763].getUb();
            env.out() << endl;
*/            
            //printModel(env, varO, rngO, objO);

     //       env.out() << "Starting Vector (pre risoluzione OMIP): " << endl;
     //       printVector(env, startingVector);
     //      cout<<"Variabili Modello OMIP: "<<endl;
     //       printVar(env,varO);

            // Risoluzione del modello
            //cplexOmip.extract(omip);
/*
            env.out() << "\n\n (3) Variabile Bastarda dopo l'estrazione del modello" << endl;

            env.out() << varO[8763].getName() << " " << varO[8763].getType() << "  " << varO[8763].getId() << "   " << varO[8763].getLb() << " <= " << varO[8763].getName() << " <= " << varO[8763].getUb();
            env.out() << endl;
*/
            // Stampa output in file di log

            s = parserLog(path, "../out/OMIP", true);
            msg = s.c_str();
            cplexOmip.exportModel(msg); // Esportazione output in un file di log dedicato

            if (counter == 1)
            {
                // Stampa output in file di log

                s = parserLog(path, "../out/1OMIP", true);
                msg = s.c_str();
                cplexOmip.exportModel(msg); // Esportazione output in un file di log dedicato
            }

            std::cout << "\n\n\n ---------------- RISOLUZIONE MODELLO OMIP: " <<counter<<" --------------------\n\n\n";

            
            /*
            CPXENVptr c_env = (CPXENVptr)cplexOmip.getImpl()->getCplexEnv();  //CPXENVptr c_env = (CPXENVptr)env.getImpl()->getEnv();
            CPXLPptr c_lp = (CPXLPptr)cplexOmip.getImpl()->getCplexLp();

            int status = CPXmipopt(c_env, c_lp);

            if (status != 0)
            {
                std::cerr << "Errore nella risoluzione del problema MIP: " << status << std::endl;
            }
            else
            {
                std::cout << "Problema risolto con successo!" << std::endl;
            }
            */
            //cplexOmip.extract(omip);
            cplexOmip.solve(); // Risoluzione di OMIP
            cplexOmip.getValues(valsO, varO);

            env.out() << "\n\n\n\nSolution status = " << cplexOmip.getStatus() << std::endl;

            status = cplexOmip.getStatus();

            env.out() << "Solution value  = " << cplexOmip.getObjValue() << "\n\n\n"
                      << endl;
//            env.out() << "Values = " << valsO << std::endl;
            env.out() << "Dimensione varO: " << varO.getSize() << endl;
            env.out() << "Dimesnione valsO: " << valsO.getSize() << endl;

            //env.out()<<valsO<<endl;
//            cout << "Lo starting vector in uscita dall' FMIP e in entrata nell' OMIP era:" << endl;
//            printVector(env, startingVector);

/*
            env.out() << "\n\n (4) Variabile Bastarda dopo la risoluzione del modello" << endl;

            env.out() << varO[8763].getName() << " " << varO[8763].getType() << "  " << varO[8763].getId() << "   " << varO[8763].getLb() << " <= " << varO[8763].getName() << " <= " << varO[8763].getUb();
            env.out() << endl;
*/
            //            env.out() << "Dimesnione di I : " << setI.getSize() << endl;
            updateVector(env, setI, startingVector, valsO, n, m); // Riempio lo starting Vector filtrando i risultati di vals dati da cplex, occhio che qua m = m+1
/*
            env.out() << "\n\n (5) Variabile Bastarda dopo l'update del vector'" << endl;

            env.out() << varO[8763].getName() << " " << varO[8763].getType() << "  " << varO[8763].getId() << "   " << varO[8763].getLb() << " <= " << varO[8763].getName() << " <= " << varO[8763].getUb();
            env.out() << endl;
*/            
            //      env.out() << "Starting Vector (post risoluzione OMIP): " << endl;
            //       printVector(env, startingVector);

            // Finito di risolvere, bisogna risettare le variabili come all'inziio, sfruttando startBounds
            // Dopodichè si restituirà il vettore inziziale come output [x^,Delta+, delta-] dove
            // Averrà una Ricerca nell'output, filtrando solo le variabili intere.
            // Inoltre bisognerà salvarsi i valori delle soluzioni delta che serviranno per costruire DeltaUB
      
      //      env.out() <<"In OMIP Bounds in cui fare il restore"<<endl;
      //      printBounds(env, startBounds);
            // Ripristino i bounds iniziali delle variabili del problema
            
            cplexOmip.extract(omip);

            restoreBounds(varO, setF, startBounds,cplexOmip,paramPercentage);
            env.out()<<"(4) Variabili dopo il restore bounds OMIP\n"<<endl;
            printVar(env, varO, cplexOmip);
 /*           env.out() << "\n\n (6) Variabile Bastarda dopo il restore dei bounds" << endl;

            env.out() << varO[8763].getName() << " " << varO[8763].getType() << "  " << varO[8763].getId() << "   " << varO[8763].getLb() << " <= " << varO[8763].getName() << " <= " << varO[8763].getUb();
            env.out() << endl;
 */           
//            env.out() << "FMIP riprristinato ai bounds originari" << endl;
//            printModel(env, varO, rngO, objO);

            startBounds.clear(); // Pulisco i bounds memorizzati in favore del riempimento ch si farà in OMIP
//            env.out() << "Bounds dopo la pulizia: " << endl;
//            printBounds(env, startBounds);
            counter++;

            std::cout<<"\n\n\nDeltaUB post OMIP (quella controlla per entrare nel if) = "<<deltaCalculator(deltaUB,startingVector,m)<<endl;

            /*
            s = parserLog(path, "../sol/", false);
            msg = s.c_str();
            cplexOmip.writeSolutions(msg);
            */
        }

        std::cout<<"ACS ended for"<< string(path)<<"..."<<endl;
    }
    catch (IloException &e)
    {
        std::cerr << "Concert exception caught: " << e << endl;
        //env.end();
    }
    catch (...)
    {
        std::cerr << "Unknown exception caught" << endl;
        //env.end();
    }

    env.end(); // Frees the resources of the environment

    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------//

// ------------------------------------------------------------------------- //
// ------------------------------ FUNCTIONS -------------------------------- // 
// ------------------------------------------------------------------------- //


void printModel(IloEnv &env, const IloNumVarArray &v, const IloRangeArray &r, const IloObjective &o){
    // Obj
    env.out() <<"Nome: "<<o.getName() << " Sense: " << o.getSense() << " Expr: " << o.getExpr() << endl << endl;
    // Vars
    for (int i = 0; i < v.getSize(); i++)
    {
        env.out() << v[i].getName() << " " << v[i].getType() << "  " << v[i].getId() <<"   " <<v[i].getLb() << " <= " << v[i].getName() << " <= " <<v[i].getUb();
        env.out() << endl;
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

void printObj(IloEnv &env, const IloObjective &o)
{
    env.out() << "Nome: " << o.getName() << " Sense: " << o.getSense() << " Expr: " << o.getExpr() << endl
              << endl;
}

void printVar(IloEnv &env, const IloNumVarArray &v , const IloCplex &cplex)
{
    CPXENVptr c_env = (CPXENVptr)cplex.getImpl()->getCplexEnv();
    CPXLPptr c_lp = (CPXLPptr)cplex.getImpl()->getCplexLp();

    for (int i = 0; i < v.getSize(); i++)
    {
        env.out() <<"Stampa con C++ API :"<< v[i].getName() << " " << v[i].getType() << " (id : " << v[i].getId() << " )  " << v[i].getLb() << " <= " << v[i].getName() << " <= " << v[i].getUb();
        env.out() << endl;

        env.out()<<"(Stampa con funzioni inedite ) "<<v[i].getName()<<" (Id : "<<v[i].getId()<<" ) è fissata a LB : "<<getDirectModelLB(i,c_env,c_lp)<<" e UB: "<<getDirectModelUB(i,c_env,c_lp)<<endl<<endl;
    }
    env.out() << endl;
}

void printRng(IloEnv &env, const IloRangeArray &r)
{
    // Stampa Range
    for (int i = 0; i < r.getSize(); i++)
    {
        // Estrai l'espressione dal vincolo corrente
        IloExpr expr = r[i].getExpr();
        bool first = true;

        // Itera sui termini nell'espressione
        for (IloExpr::LinearIterator it = expr.getLinearIterator(); it.ok(); ++it)
        {
            // Ottieni la variabile corrente
            IloNumVar v = it.getVar();

            // Ottieni il coefficiente della variabile
            double coeff = it.getCoef();

            // Stampa il coefficiente, il nome della variabile e il suo ID

            if (!first)
            {
                env.out() << " + ";
            }
            first = false;

            env.out() << coeff << " * " << v.getName() << "(" << v.getId() << ")";
        }

        // Stampa la parte destra del vincolo
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

        // Rilascia l'espressione per evitare perdite di memoria
        expr.end();
    }
}

// Fondamentale questa dim perche ci dice quanto è lungo il vettore di var senza l'aggiunta delle slack
std::vector<IloNumArray> startV1(IloEnv &env, IloNumVarArray &v, IloIntArray &setI, const IloInt &n, const IloInt &m)
{
    IloNumArray x(env);       // x^ non sappiamo ancora la dimensione, sarà grande tanto quanto il numero di variabili intere da fissare
    IloNumArray Delta1(env,m); // delta+ la sua dimensione è sicuramente m
    IloNumArray Delta2(env,m); // delta- la sua dimensione è sicuramente m

    // Inizializza il generatore di numeri casuali con un seme specifico
        

    // Prendo solo le variabili intere e metto dentro x^ un valore compreso nei bound
    for (IloInt i = 0; i < n; i++) // scorro fino a tutte le variabili prima di delta1_1, ...
    {
        if (v[i].getType() == ILOBOOL || v[i].getType() == ILOINT)
        {
            IloInt lb = static_cast<IloInt>(v[i].getLb());
            IloInt ub = static_cast<IloInt>(v[i].getUb());
            IloInt randomNumber = lb + (std::rand() % (ub - lb + 1)); // Generazione di un valore pseudocasuale nei bound
            x.add(static_cast<double>(randomNumber)); // Aggiungi come double perché x è un IloNumArray
            setI.add(i);
     
        }
    }
    env.out() << "I = " << setI << std::endl;
    env.out() << "x^ = " << x << "\t  dim x^ : " << x.getSize() << std::endl;

    // Creazione di starting vector[x^,delta+, delta-] con x fissati e i delta per ora vuoti
    std::vector<IloNumArray> sV;
    sV.push_back(x);
    sV.push_back(Delta1);
    sV.push_back(Delta2);

    return sV;
}

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

void printVector(IloEnv &env, const std::vector<IloInt> &v)
{
    env.out() << " [ ";

    for (int i = 0; i < v.size(); i++)
    {
       env.out() << v[i] << " ";
    }
    env.out() << " ]" << endl;
}

std::vector<IloInt> randomVecGenerator(int indX, int paramPercentage)
{
    std::unordered_set<IloInt> randomSet;

    IloInt randomNumber;
    
    // Genero un set di numeri casuali quindi non ripetibili
    while (randomSet.size() < paramPercentage)
    {
        randomNumber = std::rand() % (indX + 1); // Genera un numero casuale nell'intervallo [0, indX]
        randomSet.insert(randomNumber);
    }

    // Converto in un std::vector
    std::vector<IloInt> randomVec(randomSet.begin(), randomSet.end());
    return randomVec;
}

void printBounds(IloEnv &env, const std::vector<std::pair<IloNum, IloNum>> &vec)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        std::cout << "Pair " << i << ": LB = " << vec[i].first << ", UB = " << vec[i].second << std::endl;
    }
}

void variableFixing(IloNumVarArray &v, const IloIntArray &setI, IloIntArray &setF, const std::vector<IloInt> &rVec, const std::vector <IloNumArray> &sVec ,std::vector<std::pair<IloNum, IloNum>> &bounds, IloCplex &cplex, int paramPercentage)
{

    CPXENVptr c_env = (CPXENVptr)cplex.getImpl()->getCplexEnv();
    CPXLPptr c_lp = (CPXLPptr)cplex.getImpl()->getCplexLp();

    // Fissaggio delle variabili, Genero gli indici del setF, e memorizzazione dei bound iniziali
    for (int i = 0; i < paramPercentage; i++)
    {
        
    
         std::cout << std::fixed << std::setprecision(5) << "Variabile che modifico in questo variable Fixing: " << v[setF[i]].getName() << " (Id " << v[setF[i]].getId() << " ) " << "con LB : " << sVec[0][rVec[i]] << " e UB :" << sVec[0][rVec[i]] << std::endl;
    

        bounds.push_back({v[setF[i]].getLb(), v[setF[i]].getUb()}); // Memorizzo Lb e Ub della varibile con indice £ F
        // env.out()<< "varF[ "<< setF[i] << " ] "<<"  Lb : "<<varF[setF[i]].getLb()<< " Ub : "<< varF[setF[i]].getUb()<< endl;

        //IloInt a = sVec[0][rVec[i]]; // IloInt <- IloNum
        IloNum a = sVec[0][rVec[i]]; // IloNum <- IloNum

        //if(std::string(v[setF[i]].getName()) == "S_15"){cout<<"Il valore di a che sto per settare vale : "<<a<<endl; }

/*
        if(a < v[setF[i]].getLb()){
            a = v[setF[i]].getLb();
            cout<<"\n\nHo aggiornato a con LB , qualcosa non andava\n\n";
        }else if(a > v[setF[i]].getUb()){
            a = v[setF[i]].getUb();
            cout<<"\n\n Ho aggiornato a con UB, qualcosa non andava\n\n";
        }
*/

        //v[setF[i]].setBounds(a,a); 
        //v[setF[i]].setLB(a);
        //v[setF[i]].setUB(a);

        setDirectModelBounds((int)setF[i],c_env,c_lp,'B',a);  //setto entrambi i bounds allo stesso valore 
    /*
        char lu = 'B';
        //double bd = 
        int ind = ((int)setF[i]);
        cout<<"Indice : "<<ind<<endl;
        cout<<"N colonne :"<<CPXgetnumcols(c_env, c_lp)<<endl;


        if(CPXchgbds( c_env, c_lp, 1, &ind, &lu, &a ) > 0){
                cout<<"Errorre nella funzione di C (setBounds variabile fixing)"<<endl;
        }
    */    
        int index = (int) setF[i];
        //v[setF[i]].setBounds(sVec[0][rVec[i]],sVec[0][rVec[i]]); // Fisso la variabile, ovvero setto entrambi i bound con il valore salvato in x^
        //v[setF[i]].setBounds(std::ceil(sVec[0][rVec[i]]),std::ceil(sVec[0][rVec[i]])); // Fisso la variabile, ovvero setto entrambi i bound con il valore salvato in x^
        
        std::cout<<"Dopo la riga di codice che fissa , "<<v[setF[i]].getName()<<" (Id "<<v[setF[i]].getId()<<" è fissata a LB : "<<getDirectModelLB(index,c_env,c_lp)<<" e UB: "<<getDirectModelUB(index,c_env,c_lp)<<endl<<endl;
        
        
    }
}

void setFGenerator(const IloIntArray &setI, IloIntArray &setF, const std::vector<IloInt> &rVec, int paramPercentage){
      for(int i = 0; i < paramPercentage; i++){
            setF.add(setI[rVec[i]]); // Memorizzo gli indici di I in F (F c I)
        }
}
void updateVector(IloEnv &env, const IloIntArray &setI, std::vector <IloNumArray> &sVec,const IloNumArray &vals, const IloInt &n, const IloInt &m)
{

//    env.out() << "n :" << n << " m : " << m << endl;
    int j = 0;
    for (int i = 0; i < n; i++)
    {
//        cout<<"flag = SetI[ "<<j<<" ] == "<<i<<endl;

        if (setI[j] == i)
        {
//            env.out() << i << " ";
            sVec[0][j] = vals[i];
//            env.out() << " x :" << sVec[0][j] << endl;
            j++;
        }
        if (j == setI.getSize())
        {
//            cout << "Brekko" << endl;
            break;
        }
    }

    // Al primo ciclo sono vuoti i delta quindi bisognerà aggiungerli
    j = 0;
    int c = 0;
    for (int i = n; i < 2 * m + n; i++, j++)
    {
        
        if (j < m)
        {
            sVec[1][j] = vals[i];
//            env.out() << " Delta+ :" << sVec[1][j] << " j :" << j << " i : " << i << endl;
        }
        else
        {
            sVec[2][c] = vals[i];
//            env.out() << " Delta- :" << sVec[2][c] << " c :" << c << " j :" << j << " i : " << i << endl;
            c++;
        }
    }
}

void restoreBounds(IloNumVarArray &v, const IloIntArray &setF, const std::vector<std::pair<IloNum,IloNum>> &bounds, IloCplex &cplex, int paramPercentage)
{
    CPXENVptr c_env = (CPXENVptr)cplex.getImpl()->getCplexEnv();
    CPXLPptr c_lp = (CPXLPptr)cplex.getImpl()->getCplexLp();

    for (int i = 0; i < paramPercentage; i++)
    {
        if (bounds[i].first < bounds[i].second)
        {   
            int index = ((int) setF[i]);
            std::cout<<"Restore dei bounds della variabile "<<v[setF[i]].getName()<<"Con LB :"<<bounds[i].first<<" e UB :"<<bounds[i].second<<endl;
            setDirectModelBounds(index,c_env,c_lp,bounds[i].first,'L');
            setDirectModelBounds(index,c_env,c_lp,bounds[i].second,'U');

            //v[setF[i]].setBounds(bounds[i].first, bounds[i].second);
        }
        else
        {
            std::cout<<"Restore dei bounds della variabile "<<v[setF[i]].getName()<<"Con LB :"<<bounds[i].second<<" e UB :"<<bounds[i].first<<endl;
            int index = ((int) setF[i]);
            setDirectModelBounds(index,c_env,c_lp,bounds[i].second,'L');
            setDirectModelBounds(index,c_env,c_lp,bounds[i].first,'U');

            //v[setF[i]].setBounds(bounds[i].second, bounds[i].first);
        }
    }
}

IloInt deltaCalculator(IloInt delta, const std::vector<IloNumArray> &sVec,const IloInt &m)
{
    delta = 0; // Azzero ogni volta la somma, cosi da poter ricalcolare quella nuova
    for (int i = 0; i < m; i++)
    {
        delta += sVec[1][i] + sVec[2][i];
    }
    return delta;
}

double getDirectModelLB(int index, CPXENVptr c_env, CPXLPptr c_lp) {
    //int index = var.getId();
    double lb;
    CPXgetlb(c_env, c_lp, &lb, index, index);
    return lb;
}

double getDirectModelUB(int index, CPXENVptr c_env, CPXLPptr c_lp) {
    //int index = var.getId();
    double ub;
    CPXgetub(c_env, c_lp, &ub, index, index);
    return ub;
}

void setDirectModelBounds(int index, CPXENVptr c_env, CPXLPptr c_lp, char what, IloNum val)
{
    //char lu = 'B';
    //cout << "Indice : " << index << endl;
    //cout << "N colonne :" << CPXgetnumcols(c_env, c_lp) << endl;

    if (CPXchgbds(c_env, c_lp, 1, &index, &what, &val) > 0)
    {
        std::cout << "Errorre nella funzione di C (setBounds variabile fixing)" << endl;
    }
}
/*
g++ -std=c++17 -c Utility.cpp -o Utility.o
g++ -std=c++17 -c ACS.cpp -o acs.o


g++ Utility.o ACS.o  -o acs -lilocplex -lconcert -lcplex -lm -lpthread -Wno-deprecated-declarations


./acs markshare_4_0.mps.gz
*/


/*
// 3 Funzioni di clonazione che funzionano ma hanno il problema che la clonazione dei range soffre di una complessità
// quadratica che diventa molto pesante per problemi con molti vincoli e variabili in cui conviene sfruttare
// L'ottimizzazione della funzione .importModel()

void cloneVar(IloEnv &env,const IloNumVarArray &a, IloNumVarArray &b){

    for (int i = 0; i < a.getSize(); i++)
    {
        b.add(IloNumVar(env, a[i].getLb(), a[i].getUb(), a[i].getType(), a[i].getName()));
    }
}

// La copia di un range è molto più complessa: pioichè il costruttore di default mi crea delle variabili generando duplicazione
// faccio in modo di andare a copiare in un espressione temporanea il vincolo dal modello vecchio, dopodichè sostituisco variabile per variabile
// con quelle del nuovo modello per evitare la duplicazione del costruttore
void cloneRng(IloEnv &env, const IloRangeArray &a, IloRangeArray &b, const IloNumVarArray &varA, IloNumVarArray &varB) 
{
    for (int i = 0; i < a.getSize(); i++)
    {
        IloExpr expr(env);
        IloExpr originalExpr = a[i].getExpr();

        // Itera attraverso i termini nell'espressione del vincolo originale
        for (IloExpr::LinearIterator it = originalExpr.getLinearIterator(); it.ok(); ++it)
        {
            int varIndex = varA.find(it.getVar()); // Trova l'indice della variabile nel modello originale
            // IloNumVar correspondingVar = varF[varIndex]; // Usa la variabile corrispondente in varF
            expr += it.getCoef() * varB[varIndex]; // Ricostruisci l'espressione usando le variabili di fmip
        }

        // Crea un nuovo vincolo nel modello FMIP usando l'espressione ricostruita
        IloRange newRange(env, a[i].getLb(), expr, a[i].getUb(), a[i].getName());
        b.add(newRange); // Aggiungi il nuovo vincolo all'array rngF

        expr.end(); // Libera la memoria dell'espressione
    }
}
void cloneObj(IloEnv &env, const IloObjective &a, IloObjective &b)
{
    //std::cout << "Stampa f obb copiata:" << endl;

    b = IloObjective(env, a.getExpr(), a.getSense(), a.getName());
    //env.out() << b.getExpr() << endl;
}

*/




/* . Facendo cosi funzionerebbe ci sarebbe il problema di capire come gestire quello che ritorna
     nel caso dei vincoli dato che non esistono cast o metodi che aiutano a convertire i costraints in 
     Ilorange
            //clono il modello e aggiungo quello che mi serve:
            omip = IloGetClone(env,mip);

            
            // Stampa modello
            env.out() << "\n\nModello MIP dopo l'estrazione:" << endl
                      << endl;
            printModel(env, var, rng, obj);

            IloConstraintArray cosO(env);
            int a, b, z = 0;
            // Itera sugli estratti del modello clonato
            for (IloModel::Iterator it(omip); it.ok(); ++it)
            {
                IloExtractable ext = *it; // Ottieni l'elemento corrente

                if (ext.isVariable())
                {
                    varO.add(ext.asVariable());
                    a++;
                }
                else if (ext.isNumExpr())
                {
                    rngO.add((IloRange)ext);
                    b++;
                }
                else if (ext.isObjective())
                {
                    objO = ext.asObjective();
                    z++;
                }
            }

            //objO = IloGetClone(env, obj);
            objO.setName("MINIMIZE");
            //cloneObj(env,obj,objO);

            //IloNumVar(env, a[i].getLb(), a[i].getUb(), a[i].getType(), a[i].getName());

            cout << " n var: " << a << endl;
            cout << " n costraints: " << b << endl;
            cout << " n f obb: " << z << endl;

*/
/*
CON IL GENERATORE COMPLETAMENTE RANDOMICO


// Fondamentale questa dim perche ci dice quanto è lungo il vettore di var senza l'aggiunta delle slack
std::vector<IloNumArray> startV1(IloEnv &env, IloNumVarArray &v, IloIntArray &setI, const IloInt &n, const IloInt &m)
{

    IloNumArray x(env);      // x^ non sappiamo ancora la dimensione, sarà grande tanto quanto il numero di variabili intere da fissare
    IloNumArray Delta1(env,m); // delta+ la sua dimensione è sicuramente m
    IloNumArray Delta2(env,m); // delta- la sua dimensione è sicuramente m

    // IloInt n = var.getSize(); // Occhio da qua prendi la dimensione dopo l'aggiunta delle nuove variabili

    std::random_device rd;
    std::mt19937 gen(rd());
    IloInt randomNumber;
    std::uniform_int_distribution<> distr;

    int rel = 0;
    // Prendo solo le variabili intere e metto dentro x^ un valore compreso nei bound
    for (int i = 0; i < n; i++) // scorro fino a tutte le variabili prima di delta1_1, ...
    {
        // env.out()<<"flag"<<endl;
        if (v[i].getType() == ILOBOOL || v[i].getType() == ILOINT)
        {                                                                        // env.out()<<"flag2"<<endl;
            distr = std::uniform_int_distribution<>(v[i].getLb(), v[i].getUb()); // Fissaggio tramite criterio casuale di un valore
            randomNumber = distr(gen);
            x.add(randomNumber);
//            std::cout << i << ") " << rel << ") " << x[rel++] << endl;
            setI.add(i);
        }
    }
//    env.out() << "I = " << setI << endl;
//    env.out() << "x^ = " << x << "\t  dim x^ : " << x.getSize() << endl;
//    env.out() << "Rel = " << rel << endl;

    // Creazione di starting vector[x^,delta+, delta-] con x fissati e i delta per ora vuoti
    std::vector<IloNumArray> sV;
    sV.push_back(x);
    sV.push_back(Delta1);
    sV.push_back(Delta2);

    return sV;
}

// CON IL GENERATORE COMPLETAMENTE RANDOMICO

std::vector<IloInt> randomVecGenerator(int indX)
{
    std::unordered_set<IloInt> randomSet;
    std::random_device rd;
    std::mt19937 gen(rd());
    IloInt randomNumber;
    std::uniform_int_distribution<> distr(0, indX); // startingVector[0] = x^

    // Genero un set di random number quindi non ripetibili
    while (randomSet.size() < paramPercentage)
    {
        randomNumber = distr(gen);
        randomSet.insert(randomNumber);
    }

    // Converto in un std::vector
    std::vector<IloInt> randomVec(randomSet.begin(), randomSet.end());
    return randomVec;
}

*/
