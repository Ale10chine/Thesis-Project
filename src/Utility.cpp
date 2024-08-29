#include "Utility.hpp"
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem; // Creo un alias per semplicità

std::string readCommandline(int argc, char **argv)
{
    if (argc != 2)
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
