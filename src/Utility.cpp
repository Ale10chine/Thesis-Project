#include "Utility.hpp"
#include <iostream>
#include <filesystem>
#include <sstream>

namespace fs = std::filesystem; // Creo un alias per semplicità

std::string readCommandline(int argc, char **argv)
{
    if (argc != 4)
    {
        std::cout << "Inserire un unico file: \"./fileoggetto filedata.msp.gz\"" << std::endl;
        return "";
    }
    std::string cmdline = "../benchmark/" + std::string(argv[1]); // Conversione argv[1] che è un char * in una std::string
    return cmdline;
}

std::vector<std::string> readDirectory(const std::string &path)
{

    std::vector<std::string> problems;
    int i = 0;
    for (const auto &file : fs::directory_iterator(path))
    {
        if (fs::is_regular_file(file.status()))
        {                                                                          // Controllo che il file sia regolare
            problems.push_back("../benchmark/" + file.path().filename().string()); // Aggiunge il nome del path del file sottoformato di stringa al vettore di stringhe
        }
    }
    // Aggiungere ordinamento per grandezza dei file
    sort(problems.begin(), problems.end(), [](const std::string &s1, const std::string &s2)
         { return (fs::file_size(s1) < fs::file_size(s2)); });

    for (int i = 0; i < problems.size(); i++)
    {
        std::cout << problems[i] << std::endl;
    }
    return problems;
}

// Partendo da 13 sto puntando a questo elemento ../benchmark/problema1.mps
//                                                            ^
std::string parserLog(const char *path, const char *dir, bool type)
{
    std::string tmp;
    for (int i = 13; i < strlen(path); i++)
    {

        // cout << strlen(path) << endl;
        const char *c = path + i;
        if (*c == '.')
        {
            tmp = std::string(path).substr(13, i - 13); // ... = (inizio , lunghhezza sottostringa da estrarre)
            break;
        }
    }

    // Selezione formato in base a se si tratta di soluzioni o di importazione modello
    if (type)
    {
        tmp = dir + tmp + "_log.lp";
    }
    else
    {
        tmp = dir + tmp + "_log.sol";
    }

    return tmp;
}

std::string parserLog(const char *path, const char *dir) // dato un percorso "../benchmark/markshare_4_0.mps.gz" e il nome di una directory,
                                                         // Torna dir + markshare_4_0 (senza estensione)
{
    std::string tmp;
    for (int i = 13; i < strlen(path); i++)
    {

        // cout << strlen(path) << endl;
        const char *c = path + i;
        if (*c == '.')
        {
            tmp = std::string(path).substr(13, i - 13); // ... = (inizio , lunghhezza sottostringa da estrarre)
            break;
        }
    }
    
    return dir + tmp;
}

void logAndPrint(std::ofstream &logFile, const std::string &message) {
    std::cout << message;  // Stampa su terminale
    logFile << message;    // Stampa su file
}

// Funzione che cerca il nome dell'istanza e restituisce StatusStat e ObjectiveObje
std::vector<std::string> searchInstanceInCSV(const std::string &instanceName) {
    std::ifstream file("../BenchmarkSet.csv");
    std::string line, instance, status, objective;
    std::vector<std::string> data;

    if (file.is_open()) {
        // Leggi e ignora la prima riga di intestazioni
        std::getline(file, line);

        //int i = 1;
        //std::cout << "FALG 1 \n";

        // Leggi il CSV riga per riga
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            std::string temp;

            //std::cout << "FALG 2: " << i++ << "  ";

            // Leggi la prima colonna (InstanceInst.)
            instance = readQuotedCSVValue(ss);
            //std::cout << "Nella colonna leggo: " << instance << "\n";

            // Se l'istanza corrisponde, cattura i valori di StatusStat e ObjectiveObje
            if (instance == instanceName) {
                //std::cout << "FALG 3 \n";

                // Leggi la seconda colonna (StatusStat.)
                status = readQuotedCSVValue(ss);

                // Salta le colonne intermedie
                for (int i = 0; i < 8; ++i) {
                    readQuotedCSVValue(ss);
                }

                // Leggi l'undicesima colonna (ObjectiveObje.)
                objective = readQuotedCSVValue(ss);

                // Stampa i risultati
                //std::cout << "Status: " << status << "\nObjective: " << objective << std::endl;

                data.push_back(objective);
                data.push_back(status);

                return data;
            }
        }

        file.close();
    } else {
        std::cerr << "Errore nell'apertura del file!" << std::endl;
    }

    std::cout << "Instance: " << instanceName << " non trovato nel CSV!" << std::endl;
    return data;
}

// Funzione per leggere un valore da un CSV che è racchiuso tra virgolette
std::string readQuotedCSVValue(std::istringstream &ss) {
    std::string value;
    char ch;

    // Assumi che il valore inizi con una virgolette
    if (ss.get(ch) && ch == '"') {
        while (ss.get(ch) && ch != '"') {
            value += ch;
        }
    }
    
    // Controlla il carattere successivo dopo le virgolette
    // (può essere una virgola, fine riga o altro)
    if (ss.peek() == ',') {
        ss.get(); // Consuma la virgola
    }

    return value;
}


void csvPrintLine(std::ofstream &csv, std::string probName, std::string status, std::string time, std::string objACS, std::string objOp, std::string type, std::string primalGap)
{
    if(status != "Feasible"){
        //time = " / ";
        objACS = " / ";
        primalGap = " / ";
    }

    csv << "\"" << probName << "\","
        << "\"" << status << "\","
        << "\"" << time << "\","
        << "\"" << objACS << "\","
        << "\"" << objOp << "\","
        << "\"" << type << "\","
        << "\"" << primalGap << "\"\n";
}