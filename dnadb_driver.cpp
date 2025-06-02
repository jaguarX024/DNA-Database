#include "dnadb.h"
#include <math.h>
#include <algorithm>
#include <random>
#include <vector>
using namespace std;
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

class Tester{
    
};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

int main(){
    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode, DOUBLEHASH);
    bool result = true;
    
    cout << "Inserting 49 data nodes!" << endl; 
    for (int i=0;i<49;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum(), true);
        // saving data for later use
        dataList.push_back(dataObj);
        // inserting data in to the DB object
        if (!dnadb.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
    }
    // try to delete some data node
    cout << "Removing data node " << &dataList[5] << endl;
    dnadb.remove(dataList[5]);
    cout << "Removing data node " << &dataList[15] << endl;
    dnadb.remove(dataList[15]);
    dnadb.dump();
    
    // checking whether all data points are there
    cout << endl << "Checking whether all data exist in the DB:" << endl;
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++){
        DNA anObj = dnadb.getDNA((*it).getSequence(), (*it).getLocId());
        bool foundIt = (*it == anObj);
        result = result && foundIt;
        if (!foundIt){
            cout << "Data point " << 
            (*it).getSequence() << "(" <<
            (*it).getLocId() << ")" <<
            " is missing!" << endl;
        }
    }
    if (result)
        cout << "\tAll data points exist in the DnaDb object!\n";
    
    return 0;
}

unsigned int hashCode(const string str) {
    unsigned int val = 0 ;
    const unsigned int thirtyThree = 33 ;  
    for ( int i = 0 ; i < str.length(); i++)
       val = val * thirtyThree + str[i] ;
    return val ;
 }
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
