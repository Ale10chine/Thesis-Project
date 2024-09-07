#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <vector>
#include <string>
#include <fstream>

// ***togli o aggiungi static in caso dovessi usare queste due funzioni in altri file .cpp/ per questioni di leggibilità (indaga (?))
std::vector<std::string> readDirectory(const std::string &path);
std::string readCommandline(int argc, char **argv); // Restituisce una stringa vuota in caso di inserimento scorretto da terminale dell'utente
std::string parserLog(const char *path, const char *dir, bool type); // Per scrivere in file .lp o .sol
std::string parserLog(const char *path, const char *dir);
void logAndPrint(std::ofstream &logFile, const std::string &message); // Funzione per stampare contemporaneamente da terminale e su un file di log
std::vector<std::string> searchInstanceInCSV(const std::string &instanceName);
std::string readQuotedCSVValue(std::istringstream &ss);// Funzione per leggere un valore da un CSV che è racchiuso tra virgolette
void csvPrintLine(std::ofstream &csv, std::string probName, std::string status, std::string time, std::string objACS, std::string objOp, std::string type, std::string primalGap);


#endif