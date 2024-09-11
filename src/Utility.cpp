#include "Utility.hpp"
#include <iostream>
#include <filesystem>
#include <sstream>
#include <algorithm> 
#include <cstring>   

//Only for test in cluster
//#include "/home/chinelloal/Thesis Project/include/Utility.hpp"

// I create an alias for semplicity
namespace fs = std::filesystem; 

/*
    Function to read a problem from input by terminal
*/
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

/*
    Function to read in order of size (in byte) the problem of the benchmark set
*/
std::vector<std::string> readDirectory(const std::string &path)
{

    std::vector<std::string> problems;
    int i = 0;
    for (const auto &file : fs::directory_iterator(path))
    {
        // Check if the file is regular
        if (fs::is_regular_file(file.status()))
        {                                                                         
            // Add the name of the file path on string format in the vector of string
            problems.push_back("../benchmark/" + file.path().filename().string()); 
        }
    }
    //  Sorting for size (in byte) using a lambda function
    sort(problems.begin(), problems.end(), [](const std::string &s1, const std::string &s2)
         { return (fs::file_size(s1) < fs::file_size(s2)); });

    for (int i = 0; i < problems.size(); i++)
    {
        std::cout << problems[i] << std::endl;
    }
    return problems;
}

/*
    Function to parse some path, putting the path name of the benchmark problem and a destination
    this function will return a parser (the problem name without extension like .mps.gz) combined
    with two type of extencion, that are .lp or .sol useful for print the model with the C++ API
    function that is called "exportModel"

    Starting from 13 i'm pointing at this element ../benchmark/problema1.mps
                                                               ^
*/
std::string parserLog(const char *path, const char *dir, bool type)
{
    std::string tmp;
    for (int i = 13; i < strlen(path); i++)
    {

        // cout << strlen(path) << endl;
        const char *c = path + i;
        if (*c == '.')
        {   
            // ... = (start, dim of the substring thath will be extracted
            tmp = std::string(path).substr(13, i - 13); 
            break;
        }
    }
    // Choice of the extenction  
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

/*
    Overriding of the previous parserLog function, this will return only the name of the problem
    without any extenton
    
    Input: "../benchmark/markshare_4_0.mps.gz" and a name of directory, will return
    dir + "markshare_4_0" (without any extencion), where if dir is "" this function will return
    only the name of the problem like "markashare_4_0"
*/
std::string parserLog(const char *path, const char *dir) 
{
    std::string tmp;
    for (int i = 13; i < strlen(path); i++)
    {

        const char *c = path + i;
        if (*c == '.')
        {

            // ... = (start, dim of the substring thath will be extracted
            tmp = std::string(path).substr(13, i - 13); 
            break;
        }
    }
    
    return dir + tmp;
}

void logAndPrint(std::ofstream &logFile, const std::string &message)
{
    // Print on terminal
    std::cout << message;

    // Print on file log
    logFile << message;
}

// Funzione che cerca il nome dell'istanza e restituisce StatusStat e ObjectiveObje
/*
    Function thath given a porblem name like "markshare_4_0" search in the "BencharkSet.csv"
    the corrispondent row for memorize two useful parameter, that are the optimal value of the 
    objective function, and the type of the problem like hard, open or easy. For a better 
    understanding it is advisable to open and watch the structore of this .csv
*/
std::vector<std::string> searchInstanceInCSV(const std::string &instanceName) {
    std::ifstream file("../BenchmarkSet.csv");
    std::string line, instance, status, objective;
    std::vector<std::string> data;

    if (file.is_open()) {
        // Open and readonly the first line, ignoring that line
        std::getline(file, line);

        // Read the csv line by line
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            std::string temp;

            // Read the first column (InstanceInst.) 
            instance = readQuotedCSVValue(ss);

            // If we found the row we capture the ObjValue and the type of the problem 
            if (instance == instanceName) {

                // Read the second column that is (StatusStat.), so the type of the problem
                status = readQuotedCSVValue(ss);

                // Skip the intermediate columns 
                for (int i = 0; i < 8; ++i) {
                    readQuotedCSVValue(ss);
                }

                // Read the eleven column that is (ObjectiveObje.), so the optimal value of obj.f.
                objective = readQuotedCSVValue(ss);

                // Storing data in a vector for returing in main function
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

/*
    Function to read a value from a .csv file that is enclosed in quotation marks (" "),
    this function is useful for the parsing and for the read of data in the previously function 
    "searchInstanceInCSV"
*/
std::string readQuotedCSVValue(std::istringstream &ss) {
    std::string value;
    char ch;

    // For the structor of BenchmarkSet.csv we already kown that every istance start whit "" 
    if (ss.get(ch) && ch == '"') {
        while (ss.get(ch) && ch != '"') {
            value += ch;
        }
    }
    
    // Check the char after the "" that could be a , 
    if (ss.peek() == ',') {
        // Skip the ,
        ss.get(); 
    }

    return value;
}

/*
    Function to print the data row of each problem in the csv file with the stats catched from
    ACS computation
*/
void csvPrintLine(std::ofstream &csv, std::string probName, std::string status, std::string time,
                  std::string iteraction, std::string objACS, std::string objOp, std::string type, std::string primalGap)
{
    if(status != "Feasible"){
        //time = " / ";
        objACS = " / ";
        primalGap = " / ";
    }

    csv << "\"" << probName << "\","
        << "\"" << type << "\","
        << "\"" << status << "\","
        << "\"" << time << "\","
        << "\"" << iteraction << "\","
        << "\"" << objACS << "\","
        << "\"" << objOp << "\","
        << "\"" << primalGap << "\"\n";
}