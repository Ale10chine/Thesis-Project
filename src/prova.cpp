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

void cloneVar(IloEnv &env, const IloNumVarArray &a, IloNumVarArray &b);                                                 // copy from a to b (variables)
void cloneRng(IloEnv &env, const IloRangeArray &a, IloRangeArray &b, const IloNumVarArray &varA, IloNumVarArray &varB); // copy froma a to b (constraints)
void cloneObj(IloEnv &env, const IloObjective &a, IloObjective &b);                                                     // copy from a to b (function object)
void printModel(IloEnv &env, const IloNumVarArray &v, const IloRangeArray &r, const IloObjective &o);                   // Get a view of the model
void printObj(IloEnv &env, const IloObjective &o);                                                                      // View more detailed for each component of the model
void printVar(IloEnv &env, const IloNumVarArray &v);
void printRng(IloEnv &env, const IloRangeArray &r);

// Fondamentale questa n perchè ci dice quanto è lungo il vettore di var senza l'aggiunta delle slack
// Fondamentale anche il passaggio per reference di setI, che verrà aggiornato durante il programma
std::vector<IloNumArray> startV1(IloEnv &env, IloNumVarArray &v, IloIntArray &setI, const IloInt &n, const IloInt &m);
void printVector(IloEnv &env, const std::vector<IloNumArray> &v);

std::vector<IloInt> randomVecGenerator(int indX);                                 // Passo il range, ovvero tra i possibili estremi di x^ cioè tra i suoi indici [ 0 - n-1 ]
void printBounds(IloEnv &env, const std::vector<std::pair<IloNum, IloNum>> &vec); // Stampa i range che ogni volta vengono salvati delle variabili da fissare
// Funzione per il fissaggio delle variabili, oltre a queste ultime modifica F e i bounds di startingBounds
void variableFixing(IloNumVarArray &v, const IloIntArray &setI, IloIntArray &setF, const std::vector<IloInt> &rVec, const std::vector<IloNumArray> &sVec, std::vector<std::pair<IloNum, IloNum>> &bounds);

void updateVector(IloEnv &env, const IloIntArray &setI, std::vector<IloNumArray> &sVec, const IloNumArray &vals, const IloInt &n, const IloInt &m);
void restoreBounds(IloNumVarArray &v, const IloIntArray &setF, const std::vector<std::pair<IloNum, IloNum>> &bounds);

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
        env.out() << "\n\nModello MIP dopo l'estrazione:" << endl
                  << endl;
        printModel(env, var, rng, obj);

        // Creation of the FMIP sub-problem and his parameters
        try
        {
            IloModel fmip(env);
            IloCplex cplexFmip(env);

            // Declaration of some Extractables object for the FMIP sub-problem
            IloObjective objF(env);
            IloNumVarArray varF(env);
            IloRangeArray rngF(env);

            IloModel omip(env);
            IloCplex cplexOmip(env);

            // Declaration of some Extractables object for the OMIP model
            IloObjective objO(env);
            IloNumVarArray varO(env);
            IloRangeArray rngO(env);

            // #### Definizione delta 1 e 2 (variabili) per usarli nel FMIP e nel OMIP
            IloInt m = rng.getSize(); // # dei vincoli è + 1 nel OMIP (NOTA)

            // Mi salvo qui la grandezza delle variabili cosi posso fissare solo queste, i delta non li devo fissare nel FMIP!!
            IloInt n = var.getSize();

            IloIntArray setI(env); // setI (F c I) :  array per memorizzare dove si trovano le variabili inzizialemnte (solo quelle fissabili, cioè solo var intere)
            IloIntArray setF(env); // setF (F c I) : array per memorizzare solo la percentuale di incidi selezionati in modo randomico

            std::vector<std::pair<IloNum, IloNum>> startBounds; // Struttura dati che mi permette di memeorizzare i bound inizziali durante il fissaggio delle variabili

            IloNum deltaUB; // Vettore che fornirà il bound per il vincolo aggiuntivo dell' OMIP

            IloNumVarArray dpF(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);
            IloNumVarArray dnF(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);

            IloNumVarArray dpO(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);
            IloNumVarArray dnO(env, m, 0.0, IloInfinity, IloNumVar::Type::Float);

            env.out() << var.getSize() << endl; // Per vedere quanto grande è il modello fino a qua

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

            //-------------------- CREAZIONE VARIABILI FMIP--------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------


            // Deep copy dal MIP a FMIP
            // Copia variabili:
            std::cout << "FUNCTION cloneVar... \n";
            cloneVar(env, var, varF);
            std::cout << endl;

            env.out() << var.getSize() << endl;

            // ######### Starting vector for the FMIP ########### Poi qua si può implementare ALG2
            // startV1 modifica sia le varF con reference che il setI
            std::vector<IloNumArray> startingVector = startV1(env, varF, setI, n, m); // Ritorna uno starting vector e modifica il setI aggiungendo gli indici relativi alle variabili intere tale per cui poter essere fissate
            env.out() << "Starting Vector: " << endl;
            printVector(env, startingVector);

            // ######## Fissaggio variabili ##################
            // Si crea in F (setF) una percentuale di indici randomici, ai quali verranno fissate le variabili scelte dallo starting vector
            std::vector<IloInt> randomVec = randomVecGenerator(startingVector[0].getSize() - 1);
            // Funzione per il fissaggio delle variabili, mantiene in memoria i bounds inizziali e li modifica, modifica F
            variableFixing(varF, setI, setF, randomVec, startingVector, startBounds);

            env.out() << " \nF = " << setF << endl;
            env.out() << "\nstartBounds : " << endl;
            printBounds(env, startBounds);
            std::cout << "\n\n\n";

            varF.add(dpF);
            varF.add(dnF);
            fmip.add(varF);

            //-------------------- CREAZIONE VINCOLI FMIP --------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            // Copia range
            std::cout << "FUNCTION cloneRng... \n";
            cloneRng(env, rng, rngF, var, varF);
            std::cout << endl;

            // Constraints of the FMIP : A * x[i] + I_m * dpF[i] - I_m * dnF[i]
            IloExpr expr(env);
            for (int i = 0; i < m; i++)
            {
                expr = rngF[i].getExpr();
                expr += dpF[i] - dnF[i];
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
                objFExpr += dpF[i] + dnF[i];
            }
            objF = IloMinimize(env, objFExpr, "MINIMIZE");
            fmip.add(objF);

            // Stampa modello con modifiche su costraints e function object
            env.out() << "Modello FMIP aggiornato: " << endl;
            printModel(env, varF, rngF, objF);

            env.out() << "\n\nModello MIP: " << endl;
            printModel(env, var, rng, obj);

            // -------------------- RISOLUZIONE + FILE LOG FMIP--------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------
            // Risoluzione del modello
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
            env.out() << vals.getSize() << endl;

            printVector(env, startingVector);

            updateVector(env, setI, startingVector, vals, n, m); // Riempio lo starting Vector filtrando i risultati di vals dati da cplex

            env.out() << "Starting Vector (post risoluzione FMIP): " << endl;
            printVector(env, startingVector);

            // Finito di risolvere, bisogna risettare le variabili come all'inziio, sfruttando startBounds
            // Dopodichè si restituirà il vettore inziziale come output [x^,Delta+, delta-] dove
            // Averrà una Ricerca nell'output, filtrando solo le variabili intere.
            // Inoltre bisognerà salvarsi i valori delle soluzioni delta che serviranno per costruire DeltaUB
            //...
            //..
            //.

            // Ripristino i bounds iniziali delle variabili del problema
            restoreBounds(varF, setF, startBounds);
            env.out() << "FMIP riprristinato ai bounds originari" << endl;
            printModel(env, varF, rngF, objF);

            // Creazione di DeltaUB da passare come parametro poi per il vincolo del OMIP, sarà anche la condizione necessaria per accedere al FMIP
            for (int i = 0; i < m; i++)
            {
                deltaUB += startingVector[1][i] + startingVector[2][i];
            }
            env.out() << "DeltaUB: " << deltaUB << endl;

            // (IMPORTANTE, ragiona su queli elementi farla) Pulizia variabili o strutture di memoria che devono essere svuotate (NOTA)




            // Creation of the OMIP sub-problem and his parameters

            //-------------------- CREAZIONE VARIABILI OMIP--------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            // Deep copy dal MIP a OMIP
            // Copia variabili:
            std::cout << "FUNCTION cloneVar... \n";
            cloneVar(env, var, varO);
            std::cout << endl;

            // Fissaggio delle variabili leggendo il vettore x^ ?
            // ...
            // ..
            // .

            varO.add(dpO);
            varO.add(dnF);
            omip.add(varO);

            //-------------------- CREAZIONE VINCOLI OMIP --------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            // Copia range
            std::cout << "FUNCTION cloneRng... \n";
            cloneRng(env, rng, rngO, var, varO);
            std::cout << endl;

            // Constraints of the OMIP : A * x[i] + I_m * dpO[i] - I_m * dnO[i]
            //                           sum_{i = 0}^{m}(dpO[i] + dnO[i]) <= deltaUB
            // nota expr già definita in precedenza riporta su nelle dichiarazioni anche expr2
            IloExpr expr2(env);
            for (int i = 0; i < m; i++)
            {
                expr = rngO[i].getExpr();
                expr += dpO[i] - dnO[i];
                rngO[i].setExpr(expr);

                expr2 += dpO[i] + dnO[i];   
            }

            IloRange lastRange = (expr2 <= deltaUB);
            
            rngO.add(lastRange);
            omip.add(rngO);
            

            //-------------------- CREAZIONE FUNZIONE OBBIETTIVO OMIP--------------------
            // -------------------------------------------------------------------
            // -------------------------------------------------------------------

            // Copia range
            std::cout << "FUNCTION fobb (omip)... \n";
            cloneObj(env, obj, objF);
            std::cout << endl;

            
            //objO.setName("MINIMIZE");

            omip.add(objO);

            printObj(env,objF);            
            env.out() <<" OMIP" <<endl;
            env.out() << rngO.getSize()<<endl;
            //printModel(env, varO, rngO, objO); INDAGARE SUL PERCHE QUA DA SEG FOULT


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

// ------------------------------------------------------------------------------------------------------------------------------------------------------//

// ------------------------------------------------------------------------- //
// ------------------------------ FUNCTIONS -------------------------------- // 
// ------------------------------------------------------------------------- //

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

void printBounds(IloEnv &env, const std::vector<std::pair<IloNum, IloNum>> &vec)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        std::cout << "Pair " << i << ": LB = " << vec[i].first << ", UB = " << vec[i].second << std::endl;
    }
}

void variableFixing(IloNumVarArray &v, const IloIntArray &setI, IloIntArray &setF, const std::vector<IloInt> &rVec, const std::vector <IloNumArray> &sVec ,std::vector<std::pair<IloNum, IloNum>> &bounds)
{

    // Fissaggio delle variabili, Genero gli indici del setF, e memorizzazione dei bound iniziali
    for (int i = 0; i < paramPercentage; i++)
    {
        // env.out() << "Random generato : "<< randomVec[i] <<endl;
        setF.add(setI[rVec[i]]); // Memorizzo gli indici di I in F (F c I)

        bounds.push_back({v[setF[i]].getLb(), v[setF[i]].getUb()}); // Memorizzo Lb e Ub della varibile con indice £ F
        // env.out()<< "varF[ "<< setF[i] << " ] "<<"  Lb : "<<varF[setF[i]].getLb()<< " Ub : "<< varF[setF[i]].getUb()<< endl;

        v[setF[i]].setBounds(sVec[0][rVec[i]], sVec[0][rVec[i]]); // Fisso la variabile, ovvero setto entrambi i bound con il valore salvato in x^
    }
}

void updateVector(IloEnv &env, const IloIntArray &setI, std::vector <IloNumArray> &sVec,const IloNumArray &vals, const IloInt &n, const IloInt &m )
{

    env.out() << "n :" << n << " m : " << m << endl;
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        if (setI[j] == i)
        {
            env.out() << i << " ";
            sVec[0][j] = vals[i];
            env.out() << " x :" << sVec[0][j] << endl;
            j++;
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
            env.out() << " Delta+ :" << sVec[1][j] << " j :" << j << " i : " << i << endl;
        }
        else
        {
            sVec[2][c] = vals[i];
            env.out() << " Delta- :" << sVec[2][c] << " c :" << c << " j :" << j << " i : " << i << endl;
            c++;
        }
    }
}

void restoreBounds(IloNumVarArray &v, const IloIntArray &setF, const std::vector<std::pair<IloNum,IloNum>> &bounds)
{
    for (int i = 0; i < paramPercentage; i++)
    {
        v[setF[i]].setBounds(bounds[i].first, bounds[i].second);
    }
}
/*
g++ -std=c++17 -c Utility.cpp -o Utility.o
g++ -std=c++17 -c ACS.cpp -o acs.o


g++ Utility.o ACS.o  -o acs -lilocplex -lconcert -lcplex -lm -lpthread -Wno-deprecated-declarations


./acs markshare_4_0.mps.gz
*/


