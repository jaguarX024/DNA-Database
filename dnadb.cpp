#include "dnadb.h"
DnaDb::DnaDb(int size, hash_fn hash, prob_t probing = DEFPOLCY){
    m_currentCap= size;

    //check that size is within the appropriate range and is prime
    if(m_currentCap <MINPRIME)
        m_currentCap= MINPRIME;
    else if(m_currentCap> MAXPRIME)
        m_currentCap= MAXPRIME;
    else{
        if(!isPrime(m_currentCap)){
            m_currentCap= findNextPrime(m_currentCap);
        }
    }

    //hash function
    m_hash= hash;
    //current probing policy
    m_currProbing= probing;
    //memory creation for the current table
    m_currentTable= new DNA*[m_currentCap];
    for(int i=0; i< m_currentCap; i++){
        m_currentTable[i]= nullptr;
    }

    //member variables initialization
    m_currentSize=0;
    m_currNumDeleted=0;
    m_transferIndex=0;
    m_newPolicy= probing;
    m_oldTable= nullptr;
    m_oldCap=0;
    m_oldSize=0;
    m_oldNumDeleted=0;
    m_oldProbing= probing;

    
}

DnaDb::~DnaDb(){
    //delete all pointers in the current table and deallocate the memory for the table
    if (m_currentTable){
        for(int i=0; i<m_currentCap; i++){
            delete m_currentTable[i];
        }
        delete[] m_currentTable;
    
    }
    
    //delete all pointers in the old table and free memory allocated for it
    if (m_oldTable){
        for(int i=0; i<m_oldCap; i++){
            delete m_oldTable[i];
        }
        delete[] m_oldTable;
    }
    
}

void DnaDb::changeProbPolicy(prob_t policy){
    m_newPolicy= policy;
}

bool DnaDb::insert(DNA dna){
    //check location
    if (dna.getLocId() < MINLOCID || dna.getLocId()> MAXLOCID){
        return false;
    }
    //compute index
    unsigned int hashValue= m_hash(dna.getSequence());
    int originalIndex= hashValue % m_currentCap;
    int index= originalIndex;

    //collision handling
    int i=1;
    while(m_currentTable[index] != nullptr && m_currentTable[index]->m_used){
        if (*m_currentTable[index]== dna)
            return false;
        
        //apply probing policy to get the next index
        switch(m_currProbing){
            case QUADRATIC: 
                index= (originalIndex+ i*i)% m_currentCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2= 11 - (hashValue % 11);
                    index = (originalIndex + i * hash2) % m_currentCap;
                    break;
                }
            case LINEAR:
                index= (originalIndex + i)% m_currentCap;
                break;
        }
        i++;
    }
    //insert DNA object
    if (m_currentTable[index]== nullptr){
        m_currentTable[index]= new DNA(dna);

    }
    else{
        *m_currentTable[index]=dna;
    }

    m_currentTable[index]->m_used= true;
    m_currentSize++;

    //check loadfactor for possible rehash
    if (lambda()> 0.5){
        rehash();
        incrementalRehash();
    }
    else if (m_oldTable != nullptr){
        incrementalRehash();
    }

    return true;
}

bool DnaDb::remove(DNA dna){

    //determine the index of the element to be deleted in the current table
    unsigned int hashValue= m_hash(dna.getSequence());
    int originalIndex_curr= hashValue% m_currentCap;
    int index= originalIndex_curr;

    //scan the current table to find and delete dna
    int i=0;
    while(i< m_currentCap){
        if (m_currentTable[index]!= nullptr && m_currentTable[index]->m_used){
            if (*m_currentTable[index] == dna){
                m_currentTable[index]->m_used= false;
                m_currNumDeleted++;
                //check if the current table needs to be rehashed
                if((float)m_currNumDeleted > 0.8 * m_currentSize){
                    rehash();
                }
                    
                incrementalRehash(); // Call incrementalRehash()
                return true;
            }

            
        }

        //find the next index to be checked according to the probing policy
        switch(m_currProbing){
            case QUADRATIC: 
                index= (originalIndex_curr+ i*i)% m_currentCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2= 11 - (hashValue % 11);
                    index = (originalIndex_curr + i * hash2) % m_currentCap;
                    break;
                }
            case LINEAR:
                index= (originalIndex_curr + i)% m_currentCap;
                break;
        }
        i++;
        
    }


    //if dna is not deleted, check if an old table exists and remove dna from it
    if(m_oldTable){
        //determine the index of the element to be deleted in the old table
    
    int originalIndex_old= hashValue% m_oldCap;
    index= originalIndex_old;

    //scan the old table to find and delete dna
    int j=0;
    while(j< m_oldCap){
        if (m_oldTable[index]!= nullptr && m_oldTable[index]->m_used){
            if (*m_oldTable[index] == dna){
                m_oldTable[index]->m_used= false;
                
                m_oldNumDeleted++;
                incrementalRehash();
                return true;
            }

            
        }

        //find the next index to be checked according to the probing policy
        switch(m_oldProbing){
            case QUADRATIC: 
                index= (originalIndex_old+ j*j)% m_oldCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2= 11 - (hashValue % 11);
                    index = (originalIndex_old + j * hash2) % m_oldCap;
                    break;
                }
            case LINEAR:
                index= (originalIndex_old + j)% m_oldCap;
                break;
        }
        j++;
        
    }


    }
    
    return false;

}

const DNA DnaDb::getDNA(string sequence, int location) const{
    //get index using sequence
    unsigned int hashValue= m_hash(sequence);
    int originalIndex_curr= hashValue % m_currentCap;
    int index= originalIndex_curr;
    //scan the current table to find and delete dna
    int i=0;
    while(i< m_currentCap){
        if (m_currentTable[index]!= nullptr && m_currentTable[index]->m_used){
            if (m_currentTable[index]->m_location == location && m_currentTable[index]->m_sequence == sequence){
                return *m_currentTable[index];
            }

            
        }

        //find the next index to be checked according to the probing policy
        switch(m_currProbing){
            case QUADRATIC: 
                index= (originalIndex_curr+ i*i)% m_currentCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2= 11 - (hashValue % 11);
                    index = (originalIndex_curr + i * hash2) % m_currentCap;
                    break;
                }
            case LINEAR:
                index= (originalIndex_curr + i)% m_currentCap;
                break;
        }
        i++;
        
    }
    if (m_oldTable){
        int originalIndex_old= hashValue % m_oldCap;
        index= originalIndex_old;

        //find element in the old table 
        int j=0;
    while(j< m_oldCap){
        if (m_oldTable[index]!= nullptr && m_oldTable[index]->m_used){
            if (m_oldTable[index]->m_location == location && m_oldTable[index]->m_sequence == sequence ){
                return *m_oldTable[index];
            }

            
        }

        //find the next index to be checked according to the probing policy
        switch(m_oldProbing){
            case QUADRATIC: 
                index= (originalIndex_old+ j*j)% m_oldCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2= 11 - (hashValue % 11);
                    index = (originalIndex_old + j * hash2) % m_oldCap;
                    break;
                }
            case LINEAR:
                index= (originalIndex_old + j)% m_oldCap;
                break;
        }
        j++;
        
    }

    

    }
    return DNA();
}

bool DnaDb::updateLocId(DNA dna, int location){
    //get index using the dna sequence
    unsigned int hashValue= m_hash(dna.getSequence());
    int originalIndex_curr= hashValue % m_currentCap;
    int index= originalIndex_curr;
    //scan the current table to find and delete dna
    int i=0;
    while(i< m_currentCap){
        if (m_currentTable[index]!= nullptr && m_currentTable[index]->m_used){
            if ( *m_currentTable[index] == dna){
                m_currentTable[index]->m_location= location;
                return true;
            }

            
        }

        //find the next index to be checked according to the probing policy
        switch(m_currProbing){
            case QUADRATIC: 
                index= (originalIndex_curr+ i*i)% m_currentCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2= 11 - (hashValue % 11);
                    index = (originalIndex_curr + i * hash2) % m_currentCap;
                    break;
                }
            case LINEAR:
                index= (originalIndex_curr + i)% m_currentCap;
                break;
        }
        i++;
        
    }
    //if dna not found in the current table, check in the old table if it exists
    if (m_oldTable){
        int originalIndex_old= hashValue % m_oldCap;
        index= originalIndex_old;

        //find element in the old table 
        int j=0;
    while(j< m_oldCap){
        if (m_oldTable[index]!= nullptr && m_oldTable[index]->m_used){
            if ( *m_oldTable[index]== dna ){
               m_oldTable[index]->m_location= location;
               return true;
            }

            
        }

        //find the next index to be checked according to the probing policy
        switch(m_oldProbing){
            case QUADRATIC: 
                index= (originalIndex_old+ j*j)% m_oldCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2= 11 - (hashValue % 11);
                    index = (originalIndex_old + j * hash2) % m_oldCap;
                    break;
                }
            case LINEAR:
                index= (originalIndex_old + j)% m_oldCap;
                break;
        }
        j++;
        
    }

    }

    //dna is not found
    return false;
}

float DnaDb::lambda() const {
      return (float)m_currentSize/ (float)m_currentCap;
}

float DnaDb::deletedRatio() const {
    return (float)m_currNumDeleted / (float)m_currentSize;
}

void DnaDb::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool DnaDb::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) { 
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0) 
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

//function to transfer elements from old to new table when the load factor is >0.5
void DnaDb::rehash(){
    int newCap = findNextPrime(4 * (m_currentSize - m_currNumDeleted));
    DNA** newTable = new DNA*[newCap];
    for (int i = 0; i < newCap; ++i) {
        newTable[i] = nullptr;
    }

    //set new and old tables
    m_oldTable = m_currentTable;
    m_oldCap = m_currentCap;
    m_oldSize = m_currentSize;
    m_oldNumDeleted = m_currNumDeleted;
    m_oldProbing = m_currProbing;

    m_currentTable = newTable;
    m_currentCap = newCap;
    m_currentSize = 0;
    m_currNumDeleted = 0;
    m_currProbing = m_newPolicy; 

    m_transferIndex = 0;

}

void DnaDb::incrementalRehash() {
    if (m_oldTable == nullptr) {
        return; // Nothing to rehash
    }

    int transferCount = std::floor(m_oldCap * 0.25);
    for (int i = 0; i < transferCount; ++i) {
        int index = (m_transferIndex + i) % m_oldCap;
        if (m_oldTable[index] != nullptr && m_oldTable[index]->m_used) {
            // Rehash into new table
            unsigned int hashValue = m_hash(m_oldTable[index]->m_sequence);
            int newIndex = hashValue % m_currentCap;
            int j = 0;
            while (m_currentTable[newIndex] != nullptr && m_currentTable[newIndex]->m_used) {
                // Apply probing policy to find the next available slot
                  switch (m_currProbing) {
                    case QUADRATIC:
                        newIndex = (hashValue % m_currentCap + j * j) % m_currentCap;
                        break;
                    case DOUBLEHASH:
                    {
                        unsigned int hash2 = 11 - (hashValue % 11); 
                        newIndex = (hashValue % m_currentCap + j * hash2) % m_currentCap;
                        break;
                    }
                    case LINEAR:
                        newIndex = (hashValue % m_currentCap + j) % m_currentCap;
                        break;
                }
                j++;
            }
            m_currentTable[newIndex] = m_oldTable[index];
            m_currentTable[newIndex]->m_used = true;
            m_currentSize++;
            m_oldSize--;
        }
        delete m_oldTable[index];
        m_oldTable[index] = nullptr;
    }

    m_transferIndex = (m_transferIndex + transferCount) % m_oldCap;

    if (m_oldSize == 0) {
        delete[] m_oldTable;
        m_oldTable = nullptr;
        m_oldCap = 0;
    }
}