#include "dnadb.h" 

// DnaDb constructor to initialize our hash table
DnaDb::DnaDb(int size, hash_fn hash, prob_t probing = DEFPOLCY){
    m_currentCap = size;

    // Ensure the initial capacity is within the valid prime range
    if(m_currentCap < MINPRIME)
        m_currentCap = MINPRIME;
    else if(m_currentCap > MAXPRIME)
        m_currentCap = MAXPRIME;
    else{
        // If not prime, find the next prime number for capacity
        if(!isPrime(m_currentCap)){
            m_currentCap = findNextPrime(m_currentCap);
        }
    }

    
    m_hash = hash;
    // Set the initial probing policy for collision resolution
    m_currProbing = probing;
    // Dynamically allocate memory for the current hash table array
    m_currentTable = new DNA*[m_currentCap];
    // Initialize all pointers in the new table to nullptr
    for(int i=0; i< m_currentCap; i++){
        m_currentTable[i] = nullptr;
    }

    // Initialize all other member variables related to table state and rehashing
    m_currentSize = 0; // Number of active elements in the current table
    m_currNumDeleted = 0; // Number of deleted elements in the current table
    m_transferIndex = 0; // Tracks progress of incremental rehash
    m_newPolicy = probing; // Stores the policy for the next rehash
    m_oldTable = nullptr; // Pointer to the old table during rehashing
    m_oldCap = 0; // Capacity of the old table
    m_oldSize = 0; // Size of the old table
    m_oldNumDeleted = 0; // Number of deleted items in the old table
    m_oldProbing = probing; // Probing policy of the old table
}

// Destructor to properly deallocate all dynamically allocated memory
DnaDb::~DnaDb(){
    // Delete all DNA objects and then the array for the current table
    if (m_currentTable){
        for(int i=0; i<m_currentCap; i++){
            delete m_currentTable[i];
        }
        delete[] m_currentTable;
    }
    
    // Delete all DNA objects and then the array for the old table, if it exists
    if (m_oldTable){
        for(int i=0; i<m_oldCap; i++){
            delete m_oldTable[i];
        }
        delete[] m_oldTable;
    }
}

// Allows changing the probing policy for future rehashes
void DnaDb::changeProbPolicy(prob_t policy){
    m_newPolicy = policy;
}

// Inserts a DNA object into the hash table
bool DnaDb::insert(DNA dna){
    // Return false if the location ID is out of bounds
    if (dna.getLocId() < MINLOCID || dna.getLocId() > MAXLOCID){
        return false;
    }

    // Compute the hash value and initial index for the DNA object
    unsigned int hashValue = m_hash(dna.getSequence());
    int originalIndex = hashValue % m_currentCap;
    int index = originalIndex;

    // Handle collisions by probing until an empty or deleted slot is found
    int i = 1;
    while(m_currentTable[index] != nullptr && m_currentTable[index]->m_used){
        // If the DNA already exists, return false
        if (*m_currentTable[index] == dna)
            return false;
        
        // Apply the current probing policy to find the next index
        switch(m_currProbing){
            case QUADRATIC: 
                index = (originalIndex + i*i) % m_currentCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2 = 11 - (hashValue % 11); // Second hash function for double hashing
                    index = (originalIndex + i * hash2) % m_currentCap;
                    break;
                }
            case LINEAR:
                index = (originalIndex + i) % m_currentCap;
                break;
        }
        i++;
    }
    // Insert the new DNA object or overwrite a previously deleted one
    if (m_currentTable[index] == nullptr){
        m_currentTable[index] = new DNA(dna);
    }
    else{
        *m_currentTable[index] = dna;
    }

    m_currentTable[index]->m_used = true; // Mark the slot as used
    m_currentSize++; // Increment the count of active elements

    // Check load factor and trigger rehash if necessary
    if (lambda() > 0.5){
        rehash(); // Perform a full rehash to a larger table
        incrementalRehash(); // Start incremental transfer if rehashing is in progress
    }
    // If a rehash is already in progress, continue incremental transfer
    else if (m_oldTable != nullptr){
        incrementalRehash();
    }

    return true; 
}

// Removes a DNA object from the hash table
bool DnaDb::remove(DNA dna){

    // Determine the initial index for the element in the current table
    unsigned int hashValue = m_hash(dna.getSequence());
    int originalIndex_curr = hashValue % m_currentCap;
    int index = originalIndex_curr;

    // Scan the current table to find and mark the DNA as deleted
    int i = 0;
    while(i < m_currentCap){ // Iterate through potential slots
        if (m_currentTable[index] != nullptr && m_currentTable[index]->m_used){
            if (*m_currentTable[index] == dna){
                m_currentTable[index]->m_used = false; // Mark as logically deleted
                m_currNumDeleted++; // Increment deleted count
                // Trigger rehash if the ratio of deleted elements is too high
                if((float)m_currNumDeleted > 0.8 * m_currentSize){
                    rehash();
                }
                incrementalRehash(); // Continue incremental rehash if active
                return true; // DNA object successfully marked for deletion
            }
        }
        // Move to the next index based on current probing policy
        switch(m_currProbing){
            case QUADRATIC: 
                index = (originalIndex_curr + i*i) % m_currentCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2 = 11 - (hashValue % 11);
                    index = (originalIndex_curr + i * hash2) % m_currentCap;
                    break;
                }
            case LINEAR:
                index = (originalIndex_curr + i) % m_currentCap;
                break;
        }
        i++;
    }

    // If not found in the current table, check the old table if a rehash is in progress
    if(m_oldTable){
        int originalIndex_old = hashValue % m_oldCap;
        index = originalIndex_old;

        // Scan the old table to find and mark for deletion
        int j = 0;
        while(j < m_oldCap){ // Iterate through potential slots in the old table
            if (m_oldTable[index] != nullptr && m_oldTable[index]->m_used){
                if (*m_oldTable[index] == dna){
                    m_oldTable[index]->m_used = false; // Mark as logically deleted in old table
                    m_oldNumDeleted++; // Increment old table's deleted count
                    incrementalRehash(); // Continue incremental rehash
                    return true; // DNA object successfully marked for deletion
                }
            }
            // Move to the next index based on old probing policy
            switch(m_oldProbing){
                case QUADRATIC: 
                    index = (originalIndex_old + j*j) % m_oldCap;
                    break;
                case DOUBLEHASH:
                    {
                        unsigned int hash2 = 11 - (hashValue % 11);
                        index = (originalIndex_old + j * hash2) % m_oldCap;
                        break;
                    }
                case LINEAR:
                    index = (originalIndex_old + j) % m_oldCap;
                    break;
            }
            j++;
        }
    }
    
    return false; // DNA object not found or not deleted
}

// Retrieves a DNA object based on its sequence and location ID
const DNA DnaDb::getDNA(string sequence, int location) const{
    // Compute the initial index for the sequence in the current table
    unsigned int hashValue = m_hash(sequence);
    int originalIndex_curr = hashValue % m_currentCap;
    int index = originalIndex_curr;

    // Scan the current table for the DNA object
    int i = 0;
    while(i < m_currentCap){
        if (m_currentTable[index] != nullptr && m_currentTable[index]->m_used){
            // Check if both sequence and location ID match
            if (m_currentTable[index]->m_location == location && m_currentTable[index]->m_sequence == sequence){
                return *m_currentTable[index]; // Return the found DNA object
            }
        }
        // Move to the next index based on current probing policy
        switch(m_currProbing){
            case QUADRATIC: 
                index = (originalIndex_curr + i*i) % m_currentCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2 = 11 - (hashValue % 11);
                    index = (originalIndex_curr + i * hash2) % m_currentCap;
                    break;
                }
            case LINEAR:
                index = (originalIndex_curr + i) % m_currentCap;
                break;
        }
        i++;
    }

    // If not found in the current table, check the old table if a rehash is in progress
    if (m_oldTable){
        int originalIndex_old = hashValue % m_oldCap;
        index = originalIndex_old;

        // Scan the old table for the DNA object
        int j = 0;
        while(j < m_oldCap){
            if (m_oldTable[index] != nullptr && m_oldTable[index]->m_used){
                // Check if both sequence and location ID match
                if (m_oldTable[index]->m_location == location && m_oldTable[index]->m_sequence == sequence ){
                    return *m_oldTable[index]; // Return the found DNA object
                }
            }
            // Move to the next index based on old probing policy
            switch(m_oldProbing){
                case QUADRATIC: 
                    index = (originalIndex_old + j*j) % m_oldCap;
                    break;
                case DOUBLEHASH:
                    {
                        unsigned int hash2 = 11 - (hashValue % 11);
                        index = (originalIndex_old + j * hash2) % m_oldCap;
                        break;
                    }
                case LINEAR:
                    index = (originalIndex_old + j) % m_oldCap;
                    break;
            }
            j++;
        }
    }
    // If DNA object is not found in either table, return a default-constructed (empty) DNA object
    return DNA();
}

// Updates the location ID of an existing DNA object
bool DnaDb::updateLocId(DNA dna, int location){
    // Get initial index for the current table
    unsigned int hashValue = m_hash(dna.getSequence());
    int originalIndex_curr = hashValue % m_currentCap;
    int index = originalIndex_curr;

    // Scan the current table to find and update the DNA object's location
    int i = 0;
    while(i < m_currentCap){
        if (m_currentTable[index] != nullptr && m_currentTable[index]->m_used){
            if ( *m_currentTable[index] == dna){ // Check for equality of DNA objects
                m_currentTable[index]->m_location = location; // Update the location ID
                return true; // Update successful
            }
        }
        // Move to the next index based on current probing policy
        switch(m_currProbing){
            case QUADRATIC: 
                index = (originalIndex_curr + i*i) % m_currentCap;
                break;
            case DOUBLEHASH:
                {
                    unsigned int hash2 = 11 - (hashValue % 11);
                    index = (originalIndex_curr + i * hash2) % m_currentCap;
                    break;
                }
            case LINEAR:
                index = (originalIndex_curr + i) % m_currentCap;
                break;
        }
        i++;
    }
    // If not found in the current table, check the old table if a rehash is in progress
    if (m_oldTable){
        int originalIndex_old = hashValue % m_oldCap;
        index = originalIndex_old;

        // Scan the old table to find and update the DNA object's location
        int j = 0;
        while(j < m_oldCap){
            if (m_oldTable[index] != nullptr && m_oldTable[index]->m_used){
                if ( *m_oldTable[index] == dna ){ // Check for equality of DNA objects
                    m_oldTable[index]->m_location = location; // Update location ID
                    return true; // Update successful
                }
            }
            // Move to the next index based on old probing policy
            switch(m_oldProbing){
                case QUADRATIC: 
                    index = (originalIndex_old + j*j) % m_oldCap;
                    break;
                case DOUBLEHASH:
                    {
                        unsigned int hash2 = 11 - (hashValue % 11);
                        index = (originalIndex_old + j * hash2) % m_oldCap;
                        break;
                    }
                case LINEAR:
                    index = (originalIndex_old + j) % m_oldCap;
                    break;
            }
            j++;
        }
    }
    // DNA object not found in either table
    return false;
}

// Calculates the load factor of the current hash table
float DnaDb::lambda() const {
    return (float)m_currentSize / (float)m_currentCap;
}

// Calculates the ratio of deleted items to active items in the current hash table
float DnaDb::deletedRatio() const {
    return (float)m_currNumDeleted / (float)m_currentSize;
}

// Prints the contents of both the current and old hash tables
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

// Checks if a given number is prime
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

// Finds the next prime number greater than or equal to current, within defined bounds
int DnaDb::findNextPrime(int current){
    // We always stay within the range [MINPRIME-MAXPRIME]
    // The smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME - 1; // Adjust to start searching from just before MINPRIME
    for (int i = current; i < MAXPRIME; i++) { 
        for (int j = 2; j * j <= i; j++) { // Optimized primality test
            if (i % j == 0) 
                break; // Not prime
            else if (j + 1 > sqrt(i) && i != current) { 
                return i;
            }
        }
    }
    // If a user tries to go over MAXPRIME, return MAXPRIME
    return MAXPRIME;
}

// Initiates a rehash operation, creating a new, larger table
void DnaDb::rehash(){
    // Determine the new capacity (next prime after 4 times the number of active elements)
    int newCap = findNextPrime(4 * (m_currentSize - m_currNumDeleted));
    DNA** newTable = new DNA*[newCap]; // Allocate memory for the new table
    // Initialize pointers in the new table to nullptr
    for (int i = 0; i < newCap; ++i) {
        newTable[i] = nullptr;
    }

    // Move the current table to the 'old' table state for incremental rehashing
    m_oldTable = m_currentTable;
    m_oldCap = m_currentCap;
    m_oldSize = m_currentSize;
    m_oldNumDeleted = m_currNumDeleted;
    m_oldProbing = m_currProbing;

    // Set up the new table as the current table
    m_currentTable = newTable;
    m_currentCap = newCap;
    m_currentSize = 0; // Reset current size as items will be transferred
    m_currNumDeleted = 0; // Reset deleted count for the new table
    m_currProbing = m_newPolicy; // Apply the new probing policy

    m_transferIndex = 0; // Reset the transfer index for incremental rehash
}

// Performs incremental rehash, moving a portion of elements from the old to the new table
void DnaDb::incrementalRehash() {
    // Return if there's no old table to rehash from
    if (m_oldTable == nullptr) {
        return; 
    }

    // Determine how many elements to transfer in this increment (25% of old capacity)
    int transferCount = std::floor(m_oldCap * 0.25);
    // Iterate through a segment of the old table
    for (int i = 0; i < transferCount; ++i) {
        int index = (m_transferIndex + i) % m_oldCap; // Calculate index in old table

        // If the slot in the old table is used (not null and not deleted)
        if (m_oldTable[index] != nullptr && m_oldTable[index]->m_used) {
            // Rehash the element into the new table
            unsigned int hashValue = m_hash(m_oldTable[index]->m_sequence);
            int newIndex = hashValue % m_currentCap;
            int j = 0;
            // Probe for an empty slot in the new table using the current probing policy
            while (m_currentTable[newIndex] != nullptr && m_currentTable[newIndex]->m_used) {
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
            // Transfer the pointer from the old table to the new table
            m_currentTable[newIndex] = m_oldTable[index];
            m_currentTable[newIndex]->m_used = true; // Ensure it's marked as used
            m_currentSize++; // Increment current table's size
            m_oldSize--; // Decrement old table's size
        }
        // Set the old table's slot to nullptr as the element has been moved or was empty/deleted
        m_oldTable[index] = nullptr; // Ensure memory is managed (moved, not deleted here)
    }

    // Update the starting index for the next incremental transfer
    m_transferIndex = (m_transferIndex + transferCount) % m_oldCap;

    // If all elements from the old table have been transferred, deallocate the old table
    if (m_oldSize == 0) {
        delete[] m_oldTable;
        m_oldTable = nullptr;
        m_oldCap = 0;
    }
}