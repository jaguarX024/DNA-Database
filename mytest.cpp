#include "dnadb.h"
#include <math.h>
#include <algorithm>
#include <random>
#include <vector>
using namespace std;

class Tester{
    public:
    bool testInsertion();
    bool testFindDnaErrorCase();
    bool testFindDnaNormalCase();
    bool testRemoveDnaNormalCase();
    bool testRemoveDnaErrorCase();
    bool testRehash();
};
//hashcode function
unsigned int hashCode(const string str) {
    unsigned int val = 0 ;
    const unsigned int thirtyThree = 33 ;  
    for ( int i = 0 ; i < str.length(); i++)
       val = val * thirtyThree + str[i] ;
    return val ;
 }
//test insertion operation in the hash table
bool Tester::testInsertion(){
    //DNA to insert
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DNA gene3("CAGTA", 100002, false);
    DNA gene4("CATCT", 100003, false);
    DNA gene5("GAGCT", 100004, false);

    //insert DNAs into a database
    DnaDb database(MINPRIME, hashCode, LINEAR);
    database.insert(gene1);
    database.insert(gene2);

    //check if elements are inserted at the right index
    bool bgene1= *database.m_currentTable[ database.m_hash(gene1.getSequence()) % database.m_currentCap ] == gene1;
    bool bgene2= *database.m_currentTable[database.m_hash(gene2.getSequence())% database.m_currentCap] == gene2;
   
    return (bgene1 && bgene2);
    
}

//test an error case for finding a dna in the database
bool Tester::testFindDnaErrorCase(){
    //insert dna in the database
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DnaDb database(MINPRIME, hashCode, LINEAR);
    database.insert(gene1);
    database.insert(gene2);

    //try to find a dna that does not exist in the database
    DNA toFind= database.getDNA("AGCGC", 100003);
    
    return toFind== DNA();
}

//test the find operation normal case
bool Tester::testFindDnaNormalCase(){
    //insert dna in the database
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DnaDb database(MINPRIME, hashCode, LINEAR);
    database.insert(gene1);
    database.insert(gene2);

    //try to find a dna that does not exist in the database
    DNA toFind= database.getDNA("GTTTT", 100000);
    
    return toFind == gene1;
}

//test a normal case of removing an element from a database
bool Tester::testRemoveDnaNormalCase(){
    //insert dna in the database
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DnaDb database(MINPRIME, hashCode, LINEAR);
    database.insert(gene1);
    database.insert(gene2);

    return database.remove(gene1);
}

//test removing an element that does not exist in the database
bool Tester::testRemoveDnaErrorCase(){
    DNA gene1("GTTTT", 100000, false);
    DNA gene2("AGCGC", 100001, false);
    DnaDb database(MINPRIME, hashCode, LINEAR);
    database.insert(gene1);
    database.insert(gene2);

    //dna that does not exist in the database
    DNA gene3("TCACG", 100000, false);
    return (database.remove(gene3)==false);
}

//Random class
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL, SHUFFLE};
class Random {
public:
    Random(){}
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else if (type == UNIFORMREAL) { //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
        else { //the case of SHUFFLE to generate every number only once
            m_generator = std::mt19937(m_device());
        }
    }
    void setSeed(int seedNum){
        m_generator = std::mt19937(seedNum);
    }
    void init(int min, int max){
        m_min = min;
        m_max = max;
        m_type = UNIFORMINT;
        m_generator = std::mt19937(10);// 10 is the fixed seed value
        m_unidist = std::uniform_int_distribution<>(min,max);
    }
    void getShuffle(vector<int> & array){
        // this function provides a list of all values between min and max
        for (int i = m_min; i<=m_max; i++){
            array.push_back(i);
        }
        shuffle(array.begin(),array.end(),m_generator);
    }

    void getShuffle(int array[]){
        // this function provides a list of all values between min and max
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

    int getRandNum(){
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        double result = m_uniReal(m_generator);
        result = std::floor(result*100.0)/100.0;
        return result;
    }

    string getRandString(int size){
        // the parameter size specifies the length of string we ask for
        // to use ASCII char the number range in constructor must be set to 97 - 122
        string output = "";
        for (int i=0;i<size;i++){
            output = output + (char)getRandNum();
        }
        return output;
    }
    
    int getMin(){return m_min;}
    int getMax(){return m_max;}
    private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;
    std::uniform_int_distribution<> m_unidist;
    std::uniform_real_distribution<double> m_uniReal;
};

//sequencer function
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




int main(){
    Tester tester;
    //test insertion operation in the hash table
    cout<<"Test insertion operation in the database: "<<(tester.testInsertion()? "Passed": "Failed")<<endl;
    
    //test a normal case for finding a dna in the database
    cout<<"Test a normal case for finding a dna in the database : "<<(tester.testFindDnaNormalCase()? "Passed": "Failed")<<endl;

    //test a normal case of removing an element from a database
    cout<<"Test a normal case of removing a dna from the database : "<<(tester.testRemoveDnaNormalCase()? "Passed": "Failed")<<endl;

     //test an error case of removing an element from a database
     cout<<"Test an error case of removing a dna from the database : "<<(tester.testRemoveDnaErrorCase()? "Passed": "Failed")<<endl;

    //test an error case for finding a dna in the database
    cout<<"Test finding a dna that does not exist in the database : "<<(tester.testFindDnaErrorCase()? "Passed": "Failed")<<endl;

    
    
    return 0;
}