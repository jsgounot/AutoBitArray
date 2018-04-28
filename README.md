
AutoBitArray
===============

`AutoBitArray` (ABA) is a wrapper around the `BitArray` module allowing you to store any sequence containing hashable elements using a constant bit size. Originally designed for genomic sequences, this module provides a convenient way to store encode and decode sequence in a fast and reliable way using a minimum memory size.  Beside using low memory, bitarray are also mutable, allowing you to mutate your genomic sequence easily. However, if you look to a more convenient way to use mutable sequences, look at the [MutableSeq](http://biopython.org/DIST/docs/api/Bio.Seq.MutableSeq-class.html) array from BioPython which is not memory efficient but contains lot of usefull methods.

`AutoBitArray` takes two arguments, the `seq` argument which can be an iterator or a string, and a `char_size` which will indicates with how many bites a character should be encoded. Note that the limit of unique characters and the size of the bitarray is defined as 2 ^ `char_size`. For example, a `char_size` of 2 will allows 2 ^ 2 possibilities (00, 01, 10, 00) while a `char_size` of 3 will allows 2 ^ 3 possibilities. The default value of `char_size` is set to 3.

Usage
-----

```python
from aba import AutoBitArray

seq = AutoBitArray("ATGC", 2)
seq.extend("TA")
seq[0] = "GA"
print (seq) # GATGCTA
```

To retrieve all the elements, use the `decode` function.

```python
seq = AutoBitArray("ATGC" * 10, 2)
print ( "".join(seq.decode()) )
```

`MutableBitArray` can be sliced, allowing to decode only a part of the array.

```python
seq = AutoBitArray("ATGCA", 2)
for i in range(0, len(seq), 2) :
    print ("".join(seq[i:i+2].decode()))
```

If you're not sure of the `char_size` value to use :

```python
from aba import BitDict

# If you know the number of unique elements
print (BitDict.char_size_from_nunique(2)) # -> 1

# Or if you have the sequence you want to encode
print (BitDict.char_size_from_nunique("ATGC")) # -> 2
```

Memory size
-----------

Test using python 3.5.2 on a 64 bits system and a sequence of 10k elements with 4 possible characters.

| Method   | Bytes size |
|----------|-----------:|
| String   |     10.049 |
| List     |     90.312 |
| Array    |     40.264 |
| BitArray |     2.556  |

BitDict object
--------------

To encode and decode bits sequences, `MutableBitArray` uses an internal dictionary named `BitDict`. This dictionary is filled whenever a new character is observed in the input sequence. This is a requirement for the proper functioning of the instance, however it blocks the possibility to concatenate two bit arrays for which their dictionaries differ. To bypass this constraint, a `BitDict` can be initialized prior to the bit array like this :

```python
bd = BitDict()
dna = AutoBitArray("TGG", bit_dict=bd)
dna2 = AutoBitArray("ATT", bit_dict=bd)
dna3 = dna + dna2
print ( "".join(dna3.decode()) )
```

BioMba
------

BioMba contains two different class similar to biopython `SeqRecord` and `MultipleSeqAlignement`. These class are designed to produce multiple sequences alignment without string objects. The current object contains basic methods to slice the alignment and to save it into several formats (for now only fasta, phylip, structure and nexus are available).

```python
from random import choice
from bioAba import Record, MultipleSeqAlignement

bases = "ATGC"
data = {"name-%i" %j : "".join(choice(bases) for i in range(100)) for j in range(10)}
records = [Record(name, sequence) for name, sequence in data.items()]
msa = MultipleSeqAlignement(records)

# if you want to switch to biopython :
msa_biopython = msa.to_biopython()

outfile = "./myalignment.fasta"
msa.save(outfile, "fasta")
```
