#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <vector>
#include <string>

// ***togli o aggiungi static in caso dovessi usare queste due funzioni in altri file .cpp/ per questioni di leggibilità (indaga (?))
std::vector<std::string> readDirectory(const std::string &path);
std::string readCommandline(int argc, char **argv); // Restituisce una stringa vuota in caso di inserimento scorretto da terminale dell'utente
std::string parserLog(const char *path, const char *dir, bool type);

#endif