#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <unordered_set>
#include <ilcplex/ilocplex.h>

#include "Utility.hpp"

// Function FMIP, resolve a sub-MIP problem to find a vector of feasible solutions

// Function OMIP, resolve another sub-MIP problem to improve the vector returned from FMIP with respect to the original objective form the initial MIP

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

void cloneVar(IloEnv &env, const IloNumVarArray &a, IloNumVarArray &b); // copy from a to b (variables)
void cloneRng(IloEnv &env, const IloRangeArray &a, IloRangeArray &b, const IloNumVarArray &varA, IloNumVarArray &varB); // copy froma a to b (constraints)
void cloneObj(IloEnv &env, const IloObjective &a, IloObjective &b); // copy from a to b (function object)
void printModel(IloEnv &env, const IloNumVarArray &v, const IloRangeArray &r, const IloObjective &o); // Get a view of the model
void printObj(IloEnv &env, const IloObjective &o); // View more detailed for each component of the model
void printVar(IloEnv &env, const IloNumVarArray &v);
void printRng(IloEnv &env, const IloRangeArray &r);

// Fondamentale questa n perchè ci dice quanto è lungo il vettore di var senza l'aggiunta delle slack
// Fondamentale anche il passaggio per reference di setI, che verrà aggiornato durante il programma
std::vector<IloNumArray> startV1(IloEnv &env, IloNumVarArray &v, IloNumArray &setI ,const IloInt &n); 
void printVector(IloEnv &env, const std::vector<IloNumArray> & v);
std::vector<IloInt> randomVecGenerator(int indX); // Passo il range, ovvero tra i possibili estremi di x^ cioè tra i suoi indici [ 0 - n-1 ]


const int paramPercentage = 10; //  PARAM : fisso ogni volta 10 variabili

int main(int argc, char **argv)
{

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

        env.out() << var.getSize() << endl;

        mip.add(obj);
        mip.add(var);
        mip.add(rng);

        // Stampa modello 
        env.out() << "\n\nModello MIP dopo l'estrazione:" << endl <<endl;
        printModel(env, var, rng,obj);

        // Creation of the FMIP sub-problem and his parameters
        try
        {
            IloModel fmip(env);
            IloCplex cplexFmip(env);

            // Declaration of some Extractables object for the FMIP sub-problem
            IloObjective objF(env);
            IloNumVarArray varF(env);
            IloRangeArray rngF(env);

            // #### Definizione delta 1 e 2 (variabili) per usarli nel FMIP e nel OMIP
            IloInt m = rng.getSize(); // # dei vincoli
            // Mi salvo qui la grandezza delle variabili cosi posso fissare solo queste, i delta non li devo fissare nel FMIP!!
            IloInt n = var.getSize();

            IloNumArray setI(env); // setI (F c I) :  array per memorizzare dove si trovano le variabili inzizialemnte (solo quelle fissabili, cioè solo var intere)
            IloIntArray setF(env); // setF (F c I) : array per memorizzare solo la percentuale di incidi selezionati in modo randomico



            //-------------------- CREAZIONE VARIABILI FMIP--------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            IloNumVarArray dp(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);
            IloNumVarArray dn(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);

            env.out() << var.getSize() << endl; // Per vedere quanto grande è il modello fino a qua

            // Assegnazione nome alle variabili
            for (int i = 0; i < m; i++)
            {
                string name1 = "dplus" + to_string(i + 1);
                dp[i].setName(name1.c_str());

                string name2 = "dmin" + to_string(i + 1);
                dn[i].setName(name2.c_str());
            }

            // Deep copy dal MIP a FMIP
            // Copia variabili:
            std::cout << "FUNCTION cloneVar... \n";
            cloneVar(env, var, varF);
            std::cout<<endl;

            env.out() << var.getSize() << endl;
            

           
            // ######### Starting vector for the FMIP ########### Poi qua si può implementare ALG2
            // startV1 modifica sia le varF con reference che il setI
            std::vector <IloNumArray> startingVector = startV1(env,varF,setI,n); // Ritorna uno starting vector e modifica il setI aggiungendo gli indici relativi alle variabili intere tale per cui poter essere fissate
            env.out()<<"Starting Vector: "<< endl;
            printVector(env,startingVector);


            // ######## Fissaggio variabili ##################
            // Si crea in F (setF) una percentuale di indici randomici, ai quali verranno fissate le variabili scelte dallo starting vector
            std::vector<IloInt> randomVec = randomVecGenerator(startingVector[0].getSize()-1);

            std::vector<std::pair<IloInt, IloInt>> startBounds; // Struttura dati che mi permette di memeorizzare i bound inizziali durante il fissaggio delle variabili

            // Fissaggio delle variabili, Genero gli indici del setF, e memorizzazione dei bound iniziali
            for (int i = 0; i < paramPercentage; i++)
            {
               // env.out() << "Random generato : "<< randomVec[i] <<endl;
                setF.add(setI[randomVec[i]]); // Memorizzo gli indici di I in F (F c I)

                startBounds.push_back({varF[setF[i]].getLb(), varF[setF[i]].getUb()}); // Memorizzo Lb e Ub della varibile con indice £ F
                //env.out()<< "varF[ "<< setF[i] << " ] "<<"  Lb : "<<varF[setF[i]].getLb()<< " Ub : "<< varF[setF[i]].getUb()<< endl;

                varF[setF[i]].setBounds(startingVector[0][randomVec[i]],startingVector[0][randomVec[i]]); // Fisso la variabile, ovvero setto entrambi i bound con il valore salvato in x^

            }

            env.out() << " \nF = " << setF << endl;
            env.out() << "\nstartBounds : "<<endl;
            for (size_t i = 0; i < startBounds.size(); i++)
            {
                std::cout << "Pair " << i << ": LB = " << startBounds[i].first << ", UB = " << startBounds[i].second << std::endl;
            }

            std::cout<<"\n\n\n";
            varF.add(dp);
            varF.add(dn);
            fmip.add(varF);
       

            //-------------------- CREAZIONE VINCOLI FMIP --------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            // Copia range
            std::cout << "FUNCTION cloneRng... \n";
            cloneRng(env, rng, rngF, var, varF);
            std::cout<<endl;            

            // Constraints of the FMIP : A * x[i] + I_m * delta1[i] + I_m * delta2[i]
            IloExpr expr (env);
            for (int i = 0; i < m; i++)
            {
                expr = rngF[i].getExpr();
                expr += dp[i] - dn[i];   
                rngF[i].setExpr(expr);
                
            }
            fmip.add(rngF);

            //-------------------- CREAZIONE FUNZIONE OBBIETTIVO FMIP--------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            // Obj. function of the FMIP : sum_{i = 0}^{m} (delta1[i] + delta2[i])
            IloExpr objFExpr(env);
            for (int i = 0; i < m; i++)
            {
                objFExpr += dp[i] + dn[i];
            }
            objF = IloMinimize(env, objFExpr, "MINIMIZE");
            fmip.add(objF);

        


            // Stampa modello con modifiche su costraints e function object
            env.out() << "Modello FMIP aggiornato: " << endl;
            printModel(env,varF,rngF,objF);

            env.out() << "\n\nModello MIP: " << endl;
            printModel(env,var,rng,obj);


            // -------------------- RISOLUZIONE + FILE LOG FMIP--------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------
            //Risoluzione del modello
            cplexFmip.extract(fmip);

            // Stampa output in file di log
            std::string s;
            const char *msg;

            s = parserLog(path, "../out/", true); 
            msg = s.c_str();                      

            cplexFmip.exportModel(msg); // Esportazione output in un file di log dedicato


            std::cout << "\n\n\n ---------------- RISOLUZIONE MODELLO --------------------\n\n\n";
            cplexFmip.solve(); // Risoluzione di FMIP

            std::cout << "\n\n\n\n";
            env.out() << "Solution status = " << cplexFmip.getStatus() << std::endl;
            env.out() << "Solution value  = " << cplexFmip.getObjValue() << "\n\n\n"
                      << std::endl;

            
            IloNumArray vals(env); 
            cplexFmip.getValues(vals, varF);        
            env.out() << "Values = " << vals << std::endl; 

            env.out() << varF.getSize() << endl;
            env.out() << vals.getSize() <<endl;

            // Sistemare ::::::::::::::::: ripartire da qui
            /*
            int j = 0;
            for (int i = vals.getSize() - (2*m); i < vals.getSize(); i++)
            {
                if(j < m) {

                    startingVector[1][j++] = vals[i];
                } else {
                    startingVector[2][j++] = vals[i];
                }
            }
            */    
            env.out() << "Starting Vector (post risoluzione FMIP): " << endl;
            printVector(env,startingVector);

            // Finito di risolvere, bisogna risettare le variabili come all'inziio, sfruttando startBounds
            // Dopodichè si restituirà il vettore inziziale come output [x^,Delta+, delta-] dove
            // Averrà una Ricerca nell'output, filtrando solo le variabili intere.
            // Inoltre bisognerà salvarsi i valori delle soluzioni delta che serviranno per costruire DeltaUB
            //...
            //..
            //.
    
            
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


        // Creazione di DeltaUB da passare come parametro poi per il vincolo del OMIP
        //...
        //..
        //.

        // Creation of the OMIP sub-problem and his parameters
        try
        {
            IloModel omip(env);
            IloCplex cplexOmip(env);

            // Declaration of some Extractables object for the general MIP model
            IloObjective objO(env);
            IloNumVarArray varO(env);
            IloRangeArray rngO(env);

            //objO = obj;
            //std::cout << "Objective function of OMIP " << objO.getExpr() << std::endl;
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
        } else if(r[i].getLb() > IloInfinity){
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

void printVar(IloEnv &env, const IloNumVarArray &v)
{
    for (int i = 0; i < v.getSize(); i++)
    {
        env.out() << v[i].getName() << " " << v[i].getType() << "  " << v[i].getId() << "   " << v[i].getLb() << " <= " << v[i].getName() << " <= " << v[i].getUb();
        env.out() << endl;
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
std::vector<IloNumArray> startV1(IloEnv &env, IloNumVarArray &v, IloNumArray &setI, const IloInt &n)
{

    IloNumArray x(env);      // x^
   
    IloNumArray Delta1(env); // delta+
    IloNumArray Delta2(env); // delta-

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
            distr = std::uniform_int_distribution<>(v[i].getLB(), v[i].getUB()); // Fissaggio tramite criterio casuale di un valore
            randomNumber = distr(gen);
            x.add(randomNumber);
            //std::cout << i << ") " << rel << ") " << x[rel++] << endl;
            setI.add(i);
        }
    }
    //env.out() << "I = " << setI << endl;
    //env.out() << "x^ = " << x << "\n " << x.getSize() << endl;
    //env.out() << "Rel = " << rel << endl;

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

/*
g++ -std=c++17 -c Utility.cpp -o Utility.o
g++ -std=c++17 -c ACS.cpp -o acs.o


g++ Utility.o ACS.o  -o acs -lilocplex -lconcert -lcplex -lm -lpthread -Wno-deprecated-declarations


./acs markshare_4_0.mps.gz
*/


