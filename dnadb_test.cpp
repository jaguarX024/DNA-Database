#include "dnadb.h" 
#include <math.h> 
#include <algorithm> 
#include <random> 
#include <vector> 
using namespace std;

// This class will contain methods to test the DnaDb functionality
class Tester{
    public:
    // Declaration of test methods for various DnaDb operations
    bool testInsertion();
    bool testFindDnaErrorCase();
    bool testFindDnaNormalCase();
    bool testRemoveDnaNormalCase();
    bool testRemoveDnaErrorCase();
    
};

// Hash code function to generate a hash value for a given string (DNA sequence)
unsigned int hashCode(const string str) {
    unsigned int val = 0 ;
    const unsigned int thirtyThree = 33 ;   
    for ( int i = 0 ; i < str.length(); i++)
        val = val * thirtyThree + str[i] ; // Applies a common string hashing algorithm
    return val ;
}

// Implements the test for inserting DNA objects into the database
bool Tester::testInsertion(){
    // Define sample DNA sequences for insertion testing
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DNA gene3("CAGTA", 100002, false);
    DNA gene4("CATCT", 100003, false);
    DNA gene5("GAGCT", 100004, false);

    // Create a DnaDb instance with a minimum prime size, our hash function, and linear probing
    DnaDb database(MINPRIME, hashCode, LINEAR);
    // Insert two DNA objects into the database
    database.insert(gene1);
    database.insert(gene2);

    // Verify if the inserted elements are at their expected positions 
    bool bgene1= *database.m_currentTable[ database.m_hash(gene1.getSequence()) % database.m_currentCap ] == gene1;
    bool bgene2= *database.m_currentTable[database.m_hash(gene2.getSequence())% database.m_currentCap] == gene2;
    
    // Return true if both insertions are verified, false otherwise
    return (bgene1 && bgene2);
}

// Implements a test for finding a DNA object that *does not* exist in the database
bool Tester::testFindDnaErrorCase(){
    // Insert dna in the database
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DnaDb database(MINPRIME, hashCode, LINEAR);
    database.insert(gene1);
    database.insert(gene2);

    // Try to find a dna that does not exist in the database
    DNA toFind= database.getDNA("AGCGC", 100003); 
    
    // An empty DNA object (default constructor) signifies not found
    return toFind == DNA();
}

// Implements the test for successfully finding an existing DNA object in the database
bool Tester::testFindDnaNormalCase(){
    // Insert dna in the database
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DnaDb database(MINPRIME, hashCode, LINEAR);
    database.insert(gene1);
    database.insert(gene2);

    // Try to find a dna that does not exist in the database
    DNA toFind= database.getDNA("GTTTT", 100000); 
    
    // Return true if the found object matches the original, false otherwise
    return toFind == gene1;
}

// Implements a test for successfully removing an existing DNA object from the database
bool Tester::testRemoveDnaNormalCase(){
    // Insert dna in the database
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DnaDb database(MINPRIME, hashCode, LINEAR);
    database.insert(gene1);
    database.insert(gene2);

    // Attempt to remove an existing DNA object
    return database.remove(gene1);
}

// Implements a test for attempting to remove a DNA object that *does not* exist
bool Tester::testRemoveDnaErrorCase(){
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DnaDb database(MINPRIME, hashCode, LINEAR);
    database.insert(gene1);
    database.insert(gene2);

    // DNA that does not exist in the database
    DNA gene3("TCACG", 100000, false);
    // Attempt to remove the non-existent DNA object. The remove function should return false in this error case
    return (database.remove(gene3)==false);
}

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
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else if (type == UNIFORMREAL) { //the case of UNIFORMREAL to generate real numbers
            // Initialize for uniform real number distribution with a fixed seed
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
        else { //the case of SHUFFLE to generate every number only once
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
        m_generator = std::mt19937(10);// 10 is the fixed seed value
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
        std::shuffle(temp.begin(), temp.end(), m_generator);
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
        result = std::floor(result*100.0)/100.0;
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

// sequencer function to generate random DNA sequences
string sequencer(int size, int seedNum){
    //this function returns a random DNA sequence
    // size param specifies the size of string
    string sequence = "";
    Random rndObject(0,3);
    rndObject.setSeed(seedNum);
    for (int i=0;i<size;i++){
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}

// Main function to run all the tests for the DnaDb
int main(){
    Tester tester; // Create an instance of the Tester class

    // Execute and print results for each test case
    cout<<"Test insertion operation in the database: "<<(tester.testInsertion()? "Passed": "Failed")<<endl;
    cout<<"Test a normal case for finding a dna in the database : "<<(tester.testFindDnaNormalCase()? "Passed": "Failed")<<endl;
    cout<<"Test a normal case of removing a dna from the database : "<<(tester.testRemoveDnaNormalCase()? "Passed": "Failed")<<endl;
    cout<<"Test an error case of removing a dna from the database : "<<(tester.testRemoveDnaErrorCase()? "Passed": "Failed")<<endl;
    cout<<"Test finding a dna that does not exist in the database : "<<(tester.testFindDnaErrorCase()? "Passed": "Failed")<<endl;
    
    return 0; // Indicate successful execution of tests
}