# DNA-Database
This project develops a database application to efficiently manage DNA sample information. It prioritizes speed for insert, find, and remove operations using a hash table.
The application uses DNA sequences (combinations of Adenine (A), Cytosine (C), Guanine (G), and Thymine (T)) as keys. To handle potential collisions and clustering in the hash table, it implements multiple probing forms for collision resolution.

**Key Features:**

* **Efficient Data Management:** Utilizes a hash table for fast insert, find, and remove operations on DNA sample information.

* **Dynamic Resizing (Rehashing):** The hash table automatically rehashes to a new, larger table when certain load factors are exceeded or when too many deletions occur.

* **Incremental Rehashing:** Rehashing is performed incrementally during regular operations (insert/remove), transferring 25% of live nodes at a time to maintain performance. Deleted buckets are not transferred.

* **Collision Handling Policies:** Supports various probing methods for collision resolution.

* **Adaptable Collision Policy:** Allows users to change the collision handling policy, which will be applied to the new table during the next rehash.

**Classes:**

* **DnaDb**: This class implements the core database functionality, managing the hash table, handling insertions, deletions, finds, and overseeing the rehashing process.

* **DNA**: This class represents a DNA sample, with its key attribute being the DNA sequence.

**Rehashing Logic:**

* **Insertion Trigger:** If the load factor exceeds 0.5 after an insertion, the table rehashes to a new prime-sized table.

* **Deletion Trigger:** If the number of deleted buckets exceeds 80% of the total occupied buckets after a deletion, the table rehashes to a new prime-sized table.

* **Deleted Buckets:** During rehashing, deleted buckets are permanently removed and not transferred to the new table.
