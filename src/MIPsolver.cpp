#include <string>
#include <vector>
#include <filesystem> // Libreria della STL per gestione di funzioni sul filesystem
#include <ilcplex/ilocplex.h>

#include "Utility.hpp"

/*
ILOSTLBEGIN aiuta a garantire che le definizioni STL (come std::vector) siano gestite correttamente accanto ai componenti CPLEX.

Differenze con Macro Classiche:
 > Macro Classiche: Usano #define per definire costanti,
  espressioni, o blocchi di codice. Esempi sono #define PI 3.14 o #define MAX(x, y) ((x) > (y) ? (x) : (y)).
 > ILOSTLBEGIN: Non è una macro classica del preprocessore, ma una macro di convenzione definita specificamente per
  garantire la compatibilità tra CPLEX e STL. Serve a configurare l'ambiente di sviluppo in modo che le librerie possano coesistere senza conflitti.
*/
ILOSTLBEGIN
namespace fs = std::filesystem; // Creo un alias per semplicità


int main(int argc, char **argv) // Mi prendo il numero di parametri passati da command line (argc) e i file stessi con un puntatore di puntatore a caratteri
{

    IloEnv env; // Costruisce un ambiente CPLEX (CPLEX enviroment), solitamente è sempre il primo oggetto creato in qualsiasi appl.
                // E' importante in quanto solitmaente è richiesto da molti costruttori di classi della Concert technology
                // Si noti che IloEnv è una Handle class ovvero l'oggetto env è un puntatore ad un implementation object
                // Usando il copy constructor IloEnv env 2 = env; si finirebbe per creare un oggetto che punta a quello che punta env
                // Alla fine dell'utilizzo di un oggetto IloEnv  deve sempre essere distrutto tramite .end();

    bool choice = false; // rimarrà falsa se si desidera leggere dalla directory benchmark una serie di porblemi,
                         // diventerà vera nel caso in cui si decida di passare un nome di un problema da terminale

    int num; // Numero di cicli da ripetere la risoluzione del modello, in caso ci fossero più problemi da risolvere

    // Lettura di un unico problema MIP da terminale
    std::string cmdline = readCommandline(argc, argv);

    // Lettura di 4 pronblemi MIP dalla cartella benchmark
    int p = 0; // Inizizalizzo il contatore del problema p-esimo (MIP) da risolvere

    if (!cmdline.empty())
        choice = true;

    // const char* path;
    // path = cmdline.c_str(); // Per convertire una std::string in const char *
    vector<string> problems;

    if (choice)
    {
        problems.push_back(cmdline);
        num = 1;
    }
    else
    {
        
        std::cout<< "FLAG \n\n\n";
        problems = readDirectory("../benchmark");
        num = 4;
    }

    for (; p < num; p++)
    {
        const char *path = problems[p].c_str(); // Per convertire una std::string in const char *
        // cout << path << endl;
        // cout << p << endl;

        if (problems[p] != "../benchmark/markshare2.mps.gz") // Tolgo le istanze HARD per il momento (aggiungine per test piu grandi)
        {
            try
            {
                IloModel model(env); // Per definire il problema di ottimizzazione
                IloCplex cplex(env); // Classe per risolvere i modelli di ottimizzazione

                // Oggetti Extractables che derivano dalla libreria IloExtractable, gl iextractable sono oggetti che possono essere aggiunti o esratti
                // da un modello di ottimizzazione.
                IloObjective obj(env); // Puoi anche non passare env, e usare il costruttore di default in quanto viene passato env in importModel
                IloNumVarArray var(env);
                IloRangeArray rng(env);

                std::cout << "\n\n\n ---------------- IMOPORTAZIONE MODELLO: " << path + 13 << " --------------------\n\n\n"; // offset di 13 per cancellare ../benchmark
                cplex.importModel(model, path, obj, var, rng);                                                           // Lettura dai file gz o mps

                std::cout << "\n\n\n ESTRAZIONE MODELLO...";
                cplex.extract(model); // Estrazione del modello

                std::cout << "\n\n\n ---------------- RISOLUZIONE MODELLO --------------------\n\n\n";
                cplex.solve(); // Risoluzione del modello
                               // Grazie a questa funzione viene riconosciuto automaticamente se si tratta di un modello MIP o LP da risolvere
                               // In base al tipo di variabili e ai loro domini forniti
                               // Si possono comunque settare diversi parametri:
                               // cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.01); // Imposta il gap tollerato per la soluzione MIP
                               // cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1000);      // Limita il numero di nodi nel Branch-and-Bound
                               // cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Auto); // CPLEX sceglie automaticamente

                std::cout << "\n\n\n\n";
                env.out() << "Solution status = " << cplex.getStatus() << std::endl;
                env.out() << "Solution value  = " << cplex.getObjValue() << "\n\n\n" << std::endl;

                /*This method puts the solution values of the integer variables specified by the array var into the array val.
                 Array val is resized to the same size as array var, and val[i] will contain the solution val */
                IloNumArray vals(env); // Simile a std::vector<double> ma integrato a cplex, utile in LP, non è un extractables
                // IloIntArray vals(env); // Simile a std::vector<int> ma integrato a cplex, utile in MIP, non è un extractables si passa env sempre per motivi di integrazione con l'ecosistema cplex
                cplex.getValues(vals, var);               // Non esiste un ovveride di getValues che utilizzi IloIntArray, bisognerebbe castare una a una i valori delle variabili
                env.out() << "Values = " << vals << std::endl; // oggetto di tipo IloOstream che rappresenta uno stream di output associatto
                                                          // all'ambiente IloEnv, è supportata la formattazione come per std::cout (operatori <<)
                                                          // E' possibile anche settare la verbosità ovvero la quantità di informazioni stampate nel seguente modo:
                                                          // env.setParam(IloCplex::Param::MIP::Display, 2);  // Imposta il livello di verbosità del MIP

                // Stampa output in file di log
                std::string s;
                const char *msg;

                s = parserLog(path, "../out/", true);  // Chiamo il parser che mi permette di stampare nella directory corretta
                msg = s.c_str(); // Per convertire una std::string in const char *

                cplex.exportModel(msg); // Esportazione output in un file di log dedicato

                s = parserLog(path,"../sol/", false);
                msg = s.c_str();
                cplex.writeSolutions(msg);

                model.end();
            }

            catch (IloException &e)
            {
                std::cerr << "Concert exception caught: " << e << endl;
                env.end(); // Chiusura dell'enviroment anche nel caso di eccezione per prevenzione
            }
            catch (...)
            {
                std::cerr << "Unknown exception caught" << endl;
                env.end(); // Chiusura dell'enviroment anche nel caso di eccezione per prevenzione
            }
        }
        else
        {
            num++; // Nel caso incontri dei porblemi hard li skippa al momento (per facilità nei test) 
        }
    }

    env.end(); // Questa funzione serve a liberare le risorse associate all'ambiente
               // e a garantire che tutte le risorse allocate siano correttamente deallocate
    return 0;
}
