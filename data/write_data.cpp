//write_data.cpp
#include "write_data.h"

// Function to convert a vector of vectors to a string
std::string vectorToString(const Lattice& M) {
    std::ostringstream oss;
    for (const Vector& innerVec : M) {
        for (size_t i = 0; i < innerVec.size(); ++i) {
            oss << innerVec[i];
            if (i != innerVec.size() - 1) {
                oss << " "; // Space-separated values within inner vectors
            }
        }
        oss << ";"; // Semicolon-separated inner vectors
    }
    return oss.str();
}

void write_data(const Result *results, size_t size){
    // Print the current working directory
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::cout << "Current working directory: " << cwd << std::endl;
    } else {
        perror("getcwd() error");
    }
    // Create an ofstream object to write to a file
    std::ofstream outFile("output_data.csv");

    // Check if the file is open
    if (outFile.is_open()) {
    // Write the CSV header
    outFile << "ID,back_test,weights\n";

    // Write each data entry to the CSV file
    for (size_t i = 0; i < size; ++i) {
        const Result& data = results[i];
        outFile << data.id << "," << vectorToString(data.back_test.getMatrix()) << "," << vectorToString(data.weights.getMatrix()) << "\n";
    }

    // Close the file
    outFile.close();
    std::cout << "Writing to CSV file completed successfully." << std::endl;
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}