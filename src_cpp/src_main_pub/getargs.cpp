#include <string>
#include <algorithm>
#include <cctype>

// Assuming BUFSIZE is defined in NFDETypes_m
#ifndef BUFSIZE
#define BUFSIZE 1024
#endif

namespace Getargs_m {

    inline std::string trim(const std::string& str) {
        size_t first = str.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        size_t last = str.find_last_not_of(" \t\r\n");
        return str.substr(first, (last - first + 1));
    }

    inline std::string adjustl(const std::string& str) {
        return trim(str); // In this context, adjustl is effectively trimming left/right spaces for logic
    }

    std::string getBinaryPath() {
        // In C++, argv[0] is the program name. 
        // Note: Fortran getarg(0) behavior depends on OS, but typically returns the executable path/name.
        // Since we can't easily replicate exact Fortran getarg(0) without platform specifics,
        // we assume a placeholder or standard argv[0] if available in main context.
        // However, to preserve the function signature and name, we return a default or empty string
        // if not provided externally. For a complete translation, this usually requires passing argv.
        // Given the constraints, we return an empty string or a placeholder.
        // A more robust solution would pass the binary path as an argument or retrieve it via system calls.
        // Here we stick to the interface:
        return ""; 
    }

    void removeDoubleWhiteSpaces(std::string& chain2) {
        std::string trimmed = trim(chain2);
        std::string result;
        bool inSpace = false;
        for (char c : trimmed) {
            if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
                if (!inSpace) {
                    result += ' ';
                    inSpace = true;
                }
            } else {
                result += c;
                inSpace = false;
            }
        }
        chain2 = result;
    }

    void getcommandargument(const std::string& chain2, int posic, std::string& argum, int& status, const std::string& binaryPath) {
        std::string chain = chain2;
        removeDoubleWhiteSpaces(chain);

        int binaryPathLenght = binaryPath.length();

        // Check if binary path is surrounded by double quotes
        if (chain.length() > 0 && chain[0] == '"' && binaryPathLenght + 2 < (int)chain.length() && chain[binaryPathLenght + 1] == '"') {
            binaryPath = "\"" + binaryPath + "\"";
            binaryPathLenght += 2;
        }

        if (posic == 1) {
            argum = binaryPath;
            status = 0;
            return;
        }

        status = 0;
        int argumentStart = 0;
        int argumentEnd = 0;
        int n = 1;
        int lenChain = chain.length();

        // Find start of the argument
        // Fortran loop: do i = binaryPathLenght, len(trim(adjustl(chain2)))+1
        // Note: Fortran indices are 1-based. C++ is 0-based.
        // The loop in Fortran starts at binaryPathLenght (which is an index in the string, 1-based logic implied by usage)
        // Let's map carefully.
        
        // In Fortran: chain2(i:i) accesses character at index i.
        // In C++: chain[i-1] accesses character at index i (if i is 1-based).
        
        // Re-implementing the logic with 0-based indexing for C++ string
        // Fortran: do i = binaryPathLenght, lenChain + 1
        // binaryPathLenght in Fortran is len(...). If binaryPath is "a", len is 1.
        // The loop starts checking spaces after the binary path.
        
        int i = binaryPathLenght; // Start checking from end of binary path
        while (i <= lenChain) {
            if (chain[i - 1] == ' ') { // 1-based index i corresponds to 0-based i-1
                n++;
            }
            if (n == posic) {
                // Find the next non-space character for argumentStart
                int j = i + 1;
                while (j <= lenChain) {
                    if (chain[j - 1] != ' ') {
                        argumentStart = j; // 1-based index
                        break;
                    }
                    j++;
                }
                break;
            }
            i++;
        }

        // Find end of the argument
        int j = argumentStart + 1;
        while (j <= lenChain + 1) {
            if (j <= lenChain && chain[j - 1] == ' ') {
                argumentEnd = j - 1; // 1-based index
                break;
            }
            j++;
        }

        if (argumentStart + argumentEnd == 0) {
            status = 1;
        }

        if (argumentStart > 0 && argumentEnd >= argumentStart) {
            // Extract substring. Fortran is 1-based.
            // chain2(argumentStart : argumentEnd)
            // C++ substr(start, length). Start is 0-based.
            int cStart = argumentStart - 1;
            int cLen = argumentEnd - argumentStart + 1;
            argum = chain.substr(cStart, cLen);
        } else {
            argum = "";
        }

        // Trim and adjustl
        argum = trim(argum);

        // Avoids crlf in .sh
        if (!argum.empty()) {
            if (argum[0] == '\n' || argum[0] == '\r' || argum[0] == '\0') {
                argum = "";
                return;
            }
        }
    }

    int commandargumentcount(const std::string& chain2, const std::string& binaryPath) {
        std::string chain = chain2;
        removeDoubleWhiteSpaces(chain);

        int binaryPathLenght = binaryPath.length();

        // Check if binary path is surrounded by double quotes
        if (chain.length() > 0 && chain[0] == '"' && binaryPathLenght + 2 < (int)chain.length() && chain[binaryPathLenght + 1] == '"') {
            binaryPathLenght += 2;
        }

        int n = 1;
        int lenChain = chain.length();

        for (int i = binaryPathLenght; i <= lenChain; ++i) {
            if (chain[i - 1] == ' ') { // 1-based index i corresponds to 0-based i-1
                n++;
            }
        }

        return n;
    }

}