#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <vector>
#include <string>
#include <fstream>


std::vector<std::string> readDirectory(const std::string &path);

std::string readCommandline(int argc, char **argv); 

std::string parserLog(const char *path, const char *dir, bool type); 

std::string parserLog(const char *path, const char *dir);

void logAndPrint(std::ofstream &logFile, const std::string &message); 

std::vector<std::string> searchInstanceInCSV(const std::string &instanceName);

std::string readQuotedCSVValue(std::istringstream &ss);

void csvPrintLine(std::ofstream &csv, std::string probName, std::string status, std::string time,
                  std::string iteraction, std::string objACS, std::string objOp, std::string type, std::string primalGap);

#endif