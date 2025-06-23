#include "dnadb.h" 
#include <math.h> 
#include <algorithm> 
#include <random> 
#include <vector>
using namespace std;

// Enum to define different types of random number distributions
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL, SHUFFLE};

// Random class for generating various types of random numbers and sequences
class Random {
public:
    Random(){} // Default constructor

    // Constructor to initialize the random number generator with specific parameters
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            // Initialize for normal distribution 
            m_generator = std::mt19937(m_device()); // Uses a non-deterministic seed from hardware
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            // Initialize for uniform integer distribution with a fixed seed
            m_generator = std::mt19937(10);
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else if (type == UNIFORMREAL) {
            // Initialize for uniform real number distribution with a fixed seed
            m_generator = std::mt19937(10);
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
        else {
            // Initialize for shuffling, using a non-deterministic seed
            m_generator = std::mt19937(m_device());
        }
    }

    // Allows setting a custom seed for the random number generator
    void setSeed(int seedNum){
        m_generator = std::mt19937(seedNum);
    }

    // Initializes the random generator for uniform integer distribution
    void init(int min, int max){
        m_min = min;
        m_max = max;
        m_type = UNIFORMINT;
        m_generator = std::mt19937(10); // Using a fixed seed for reproducibility
        m_unidist = std::uniform_int_distribution<>(min,max);
    }

    // Populates a vector with numbers from min to max and shuffles them
    void getShuffle(vector<int> & array){
        for (int i = m_min; i<=m_max; i++){
            array.push_back(i);
        }
        shuffle(array.begin(),array.end(),m_generator);
    }

    // Populates an array with numbers from min to max and shuffles them
    void getShuffle(int array[]){
        vector<int> temp;
        for (int i = m_min; i<=m_max; i++){
            temp.push_back(i);
        }
        std::shuffle(temp.begin(), temp.end(), m_generator); // Use std::shuffle for randomizing
        vector<int>::iterator it;
        int i = 0;
        for (it=temp.begin(); it != temp.end(); it++){
            array[i] = *it;
            i++;
        }
    }

    // Generates a random integer based on the selected distribution type
    int getRandNum(){
        int result = 0;
        if(m_type == NORMAL){
            // For normal distribution, ensure result is within min and max
            result = m_min - 1; // Initialize to an invalid value to enter loop
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            // For uniform integer distribution
            result = m_unidist(m_generator);
        }
        return result;
    }

    // Generates a random real number, rounded to two decimal places
    double getRealRandNum(){
        double result = m_uniReal(m_generator);
        result = std::floor(result*100.0)/100.0; // Rounds down to two decimal places
        return result;
    }

    // Generates a random string of specified size, using ASCII characters
    string getRandString(int size){
        string output = "";
        for (int i=0;i<size;i++){
            output = output + (char)getRandNum(); // Appends random character
        }
        return output;
    }
    
    // Getter for minimum value
    int getMin(){return m_min;}
    // Getter for maximum value
    int getMax(){return m_max;}

private:
    int m_min; // Minimum value for random generation
    int m_max; // Maximum value for random generation
    RANDOM m_type; // Type of random distribution
    std::random_device m_device; // Non-deterministic random number generator source
    std::mt19937 m_generator; // Mersenne Twister pseudo-random number generator
    std::normal_distribution<> m_normdist; // Normal distribution object
    std::uniform_int_distribution<> m_unidist; // Uniform integer distribution object
    std::uniform_real_distribution<double> m_uniReal; // Uniform real distribution object
};


class Tester{
    
};

// Function declaration for hashing DNA sequences
unsigned int hashCode(const string str);
// Function declaration for generating DNA sequences
string sequencer(int size, int seedNum);

// Main function where the DNA database operations are tested
int main(){
    vector<DNA> dataList; // Stores DNA objects for later verification
    Random RndLocation(MINLOCID,MAXLOCID); // Random generator for location IDs
    // Initialize DnaDb with a minimum prime size, the hashCode function, and a collision policy
    DnaDb dnadb(MINPRIME, hashCode, DOUBLEHASH); 
    bool result = true; // Flag for checking data integrity
    
    cout << "Inserting 49 data nodes!" << endl; 
    for (int i=0;i<49;i++){
        // Generate a random DNA object with sequence and location ID
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        // Store the generated DNA object in a list
        dataList.push_back(dataObj);
        // Insert the DNA object into our custom DNA database
        if (!dnadb.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
    }

    // Demonstrate removing specific data nodes from the database
    cout << "Removing data node " << &dataList[5] << endl;
    dnadb.remove(dataList[5]);
    cout << "Removing data node " << &dataList[15] << endl;
    dnadb.remove(dataList[15]);
    dnadb.dump(); // Display the current state of the database after removals
    
    // Verify if all (non-removed) data points still exist in the database
    cout << endl << "Checking whether all data exist in the DB:" << endl;
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++){
        // Attempt to retrieve DNA object from the database using its sequence and location ID
        DNA anObj = dnadb.getDNA((*it).getSequence(), (*it).getLocId());
        // Compare the retrieved object with the original
        bool foundIt = (*it == anObj);
        result = result && foundIt; // Update overall result
        if (!foundIt){
            // Report if any data point is missing
            cout << "Data point " << 
            (*it).getSequence() << "(" <<
            (*it).getLocId() << ")" <<
            " is missing!" << endl;
        }
    }
    // Final check for overall data integrity
    if (result)
        cout << "\tAll data points exist in the DnaDb object!\n";
    
    return 0; // Indicate successful execution
}

// Hash function implementation 
unsigned int hashCode(const string str) {
    unsigned int val = 0 ;
    const unsigned int thirtyThree = 33 ;   // A common multiplier in hash functions
    for ( int i = 0 ; i < str.length(); i++)
        val = val * thirtyThree + str[i] ; // Accumulates hash value
    return val ;
}

// Function to generate a random DNA sequence
string sequencer(int size, int seedNum){
    string sequence = "";
    // Random object to pick characters (A, C, G, T)
    Random rndObject(0,3); // Range 0-3 for indexing ALPHA array
    rndObject.setSeed(seedNum); // Set seed for sequence generation
    // Loop to build the DNA sequence character by character
    for (int i=0;i<size;i++){
        // ALPHA is assumed to be a global constant array/string like "ACGT"
        sequence = sequence + ALPHA[rndObject.getRandNum()]; 
    }
    return sequence;
}